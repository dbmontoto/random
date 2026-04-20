/**
 * @file gibbs_tp_optimizer.cpp
 * @brief Implementation of EOS-agnostic multi-phase Gibbs minimizer at constant T,P.
 */

#include "thermo/equilibrium/optimization/gibbs_tp_optimizer.h"

#include "thermo/core/constants.h"
#include "thermo/equilibrium/density_selection.h"
#include "thermo/equilibrium/internal/k_values.h"
#include "thermo/equilibrium/internal/rachford_rice.h"
#include "thermo/equilibrium/optimization/phase_allocations.h"
#include "thermo/equilibrium/optimization/internal/gibbs_tp_outer_loop_problem.h"
#include "thermo/equilibrium/optimization/outer_loop_helpers.h"
#include "thermo/equilibrium/optimization/phase_postprocess.h"
#include "thermo/equilibrium/optimization/phase_stability_outer_loop.h"
#include "thermo/equilibrium/optimization/tp_phase_split.h"
#include "thermo/factory/numerics_factory.h"
#include "thermo/factory/stability_factory.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

GibbsTPOptimizer::GibbsTPOptimizer(EOSPtr eos, Stability::StabilityAnalyzerPtr stability)
    : eos_(std::move(eos)), stability_(std::move(stability))
{
    if (!eos_) {
        throw std::invalid_argument("GibbsTPOptimizer: EOS pointer is null");
    }
    if (!stability_) {
        stability_ = Factory::StabilityFactory::create(eos_);
    }
}

EquilibriumState GibbsTPOptimizer::minimizeTP(
    double T,
    double P,
    const std::vector<double>& z_in,
    const Config::EquilibriumOptimizerConfig& cfg) const
{
    EquilibriumState out;
    out.temperature = T;
    out.pressure = P;
    out.z = z_in;
    out.method_used = "GibbsOptimizer(T,P)";

    const int nc = eos_->numComponents();
    if (static_cast<int>(z_in.size()) != nc) {
        out.converged = false;
        out.message = "Composition size mismatch";
        return out;
    }
    if (!(T > 0.0 && P > 0.0)) {
        out.converged = false;
        out.message = "Invalid T/P";
        return out;
    }

    const std::vector<double> z = clampAndNormalize(z_in, Constants::MIN_MOLE_FRACTION);

    // Early stability check: single-phase return.
    if (cfg.use_phase_stability_outer_loop) {
        const auto stab_cfg = stabilityConfigFromEquilibriumConfig(cfg);
        const auto stab = stability_->analyze(T, P, z, stab_cfg);
        if (stab.is_stable) {
            const Core::Units::PropertyVariable pressure_tolerance(cfg.density_newton_tol, Core::Units::Unit::PA);
            const auto ev = Detail::evalAtTP(
                *eos_, *stability_, T, P, z, PhaseType::Unknown,
                cfg.phase_detection_space, 0.0, cfg.max_density_newton_iters, pressure_tolerance
            );
            PhaseState ps;
            ps.type = ev.phase;
            ps.density = ev.rho;
            ps.compressibility = eos_->compressibility(T, ps.density, z);
            ps.x = z;
            ps.phi = ev.phi;
            ps.fraction = 1.0;
            out.phases = {ps};
            out.converged = true;
            out.objective_value = phaseGibbsRT(ps.x, ev.lnphi, cfg.composition_min);
            out.infeasibility = 0.0;
            out.score_value = out.objective_value;
            out.message = "Single phase (stable)";
            return out;
        }
    }

    const int max_phases = std::clamp(cfg.max_phases, 1, 3);
    int M = std::min(2, max_phases);

    // Initial allocations from a Wilson K-based 2-phase guess.
    AllocState alloc;
    alloc.nc = nc;
    alloc.M = M;
    alloc.w.assign(static_cast<size_t>(nc * M), 0.0);
    if (M == 1) {
        for (int i = 0; i < nc; ++i) alloc.w[static_cast<size_t>(i)] = 1.0;
    } else {
        const auto K = Equilibrium::Internal::wilsonK(*eos_, T, P);
        double beta_v = 0.5;
        (void)Equilibrium::Internal::solveRachfordRiceBeta(z, K, beta_v);
        beta_v = std::clamp(beta_v, 1e-6, 1.0 - 1e-6);

        std::vector<double> x_liq(static_cast<size_t>(nc), 0.0), y_vap(static_cast<size_t>(nc), 0.0);
        for (int i = 0; i < nc; ++i) {
            const double denom = 1.0 + beta_v * (K[i] - 1.0);
            x_liq[i] = z[i] / std::max(denom, Constants::MIN_MOLE_FRACTION);
            y_vap[i] = K[i] * x_liq[i];
        }
        x_liq = clampAndNormalize(x_liq, Constants::MIN_MOLE_FRACTION);
        y_vap = clampAndNormalize(y_vap, Constants::MIN_MOLE_FRACTION);

        for (int i = 0; i < nc; ++i) {
            const double zi = std::max(z[i], Constants::MIN_MOLE_FRACTION);
            const double wvap = std::clamp((beta_v * y_vap[i]) / zi, 0.0, 1.0);
            alloc.w[static_cast<size_t>(i * M + 0)] = wvap;
            alloc.w[static_cast<size_t>(i * M + 1)] = 1.0 - wvap;
        }
    }

    std::vector<double> u = packLogitsFromAllocations(alloc, cfg.composition_min);

    auto optimizer = Factory::NumericsFactory::createBFGS();
    Internal::GibbsTPOuterLoopProblem problem{*eos_, *stability_, T, P, z, nc, cfg, out.method_used};
    try {
        const auto res = Numerics::Optimization::runPhaseStabilityOuterLoop(
            *optimizer, problem, M, max_phases, u, outerLoopConfigFrom(cfg));
        out = res.best_state;
    } catch (const std::exception& e) {
        out.converged = false;
        out.message = e.what();
        return out;
    }
    if (!out.converged && out.message.empty()) {
        out.message = "Did not converge";
    }
    return out;
}

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo
