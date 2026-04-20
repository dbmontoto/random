/**
 * @file reactive_gibbs_tp_optimizer.cpp
 * @brief Implementation of EOS-agnostic TP Gibbs minimizer with reaction extents.
 */

#include "thermo/equilibrium/optimization/reactive_gibbs_tp_optimizer.h"

#include "thermo/core/constants.h"
#include "thermo/equilibrium/density_selection.h"
#include "thermo/equilibrium/optimization/internal/reactive_tp_outer_loop_problem.h"
#include "thermo/equilibrium/optimization/phase_allocations.h"
#include "thermo/equilibrium/optimization/variable_blocks.h"
#include "thermo/equilibrium/optimization/outer_loop_helpers.h"
#include "thermo/equilibrium/optimization/phase_postprocess.h"
#include "thermo/equilibrium/optimization/phase_stability_outer_loop.h"
#include "thermo/equilibrium/optimization/tp_phase_split.h"
#include "thermo/equilibrium/optimization/variable_layout.h"
#include "thermo/factory/numerics_factory.h"
#include "thermo/factory/stability_factory.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

ReactiveGibbsTPOptimizer::ReactiveGibbsTPOptimizer(EOSPtr eos, Stability::StabilityAnalyzerPtr stability)
    : eos_(std::move(eos)), stability_(std::move(stability))
{
    if (!eos_) {
        throw std::invalid_argument("ReactiveGibbsTPOptimizer: EOS pointer is null");
    }
    if (!stability_) {
        stability_ = Factory::StabilityFactory::create(eos_);
    }
}

EquilibriumState ReactiveGibbsTPOptimizer::minimizeTPReactive(
    double T,
    double P,
    const std::vector<double>& n0,
    const Reactions::ReactionSystem& system,
    const std::vector<double>& mu0,
    const ReactiveGibbsTPConfig& cfg) const
{
    EquilibriumState out;
    out.temperature = T;
    out.pressure = P;
    out.method_used = "ReactiveGibbsOptimizer(T,P)";
    out.z.clear();

    if (!(T > 0.0 && P > 0.0) || !std::isfinite(T) || !std::isfinite(P)) {
        out.converged = false;
        out.message = "Invalid T/P";
        return out;
    }

    const int nc = eos_->numComponents();
    if (static_cast<int>(n0.size()) != nc) {
        out.converged = false;
        out.message = "n0 size mismatch";
        return out;
    }
    if (static_cast<int>(mu0.size()) != nc) {
        out.converged = false;
        out.message = "mu0 size mismatch";
        return out;
    }
    if (system.numComponents() != nc) {
        out.converged = false;
        out.message = "Reaction system component mismatch";
        return out;
    }
    if (!(std::isfinite(cfg.standard_pressure) && cfg.standard_pressure > 0.0)) {
        out.converged = false;
        out.message = "Invalid standard_pressure";
        return out;
    }
    if (!(std::isfinite(cfg.min_moles) && cfg.min_moles >= 0.0)) {
        out.converged = false;
        out.message = "Invalid min_moles";
        return out;
    }

    const int nr = system.numReactions();

    // Initial composition from n0.
    double N0 = 0.0;
    for (double ni : n0) N0 += std::max(ni, cfg.min_moles);
    if (!(std::isfinite(N0) && N0 > 0.0)) {
        out.converged = false;
        out.message = "Invalid initial moles";
        return out;
    }
    std::vector<double> z0(static_cast<size_t>(nc), 0.0);
    for (int i = 0; i < nc; ++i) z0[static_cast<size_t>(i)] = std::max(n0[static_cast<size_t>(i)], cfg.min_moles) / N0;
    z0 = clampAndNormalize(z0, Constants::MIN_MOLE_FRACTION);

    // Start single-phase; outer loop may add phases via stability analysis.
    const int max_phases = std::clamp(cfg.equilibrium.max_phases, 1, 3);
    int M = 1;

    // Extent bounds + initial u_xi.
    std::vector<Reactions::ExtentBounds> xi_bounds;
    xi_bounds.reserve(static_cast<size_t>(nr));
    std::vector<double> xi0(static_cast<size_t>(nr), 0.0);
    for (int r = 0; r < nr; ++r) {
        xi_bounds.push_back(system.extentBounds(r, n0, cfg.min_moles));
        double x = 0.0;
        if (x < xi_bounds.back().min || x > xi_bounds.back().max) {
            x = 0.5 * (xi_bounds.back().min + xi_bounds.back().max);
        }
        xi0[static_cast<size_t>(r)] = std::clamp(x, xi_bounds.back().min, xi_bounds.back().max);
    }
    const std::vector<double> u_xi = encodeExtents(xi0, xi_bounds);

    // Initial allocations (M=1).
    AllocState alloc;
    alloc.nc = nc;
    alloc.M = M;
    alloc.w.assign(static_cast<size_t>(nc * M), 1.0);

    std::vector<double> u_alloc = packLogitsFromAllocations(alloc, cfg.equilibrium.composition_min);

    auto optimizer = Factory::NumericsFactory::createBFGS();
    std::vector<double> vars;
    vars.reserve(u_xi.size() + u_alloc.size());
    vars.insert(vars.end(), u_xi.begin(), u_xi.end());
    vars.insert(vars.end(), u_alloc.begin(), u_alloc.end());

    Internal::ReactiveTPOuterLoopProblem problem{
        *eos_,
        *stability_,
        T,
        P,
        n0,
        system,
        mu0,
        cfg,
        xi_bounds,
        nc,
        nr,
        out.method_used
    };

    const auto res = Numerics::Optimization::runPhaseStabilityOuterLoop(
        *optimizer, problem, M, max_phases, vars, outerLoopConfigFrom(cfg.equilibrium));
    out = res.best_state;
    if (!out.converged && out.message.empty()) {
        out.message = "Did not converge";
    }
    return out;
}

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo
