/**
 * @file tp_phase_split.h
 * @brief Shared evaluation utilities for TP phase-split optimization problems.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_TP_PHASE_SPLIT_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_TP_PHASE_SPLIT_H

#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/density_selection.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/equilibrium/optimization/phase_allocations.h"

#include <cmath>
#include <limits>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {
namespace TP {

struct PhaseSplitOptions {
    bool build_phases = false;
    // If > 0, adds this amount to infeasibility for any phase that evaluates to non-finite or non-positive rho.
    // When this triggers, the phase's Gibbs contribution can optionally be skipped.
    double invalid_rho_penalty = 0.0;
    bool skip_objective_on_invalid_rho = false;
};

struct PhaseSplitEvaluation {
    std::vector<double> betas;
    std::vector<std::vector<double>> x;
    std::vector<PhaseState> phases;
    double phase_objective = 0.0; // sum_k beta_k * g_k/(RT)
    double infeasibility = 0.0;   // beta-min penalty + optional rho-invalid penalties
};

inline PhaseSplitEvaluation evaluatePhaseSplit(
    const EOS& eos,
    const Stability::IStabilityAnalyzer& stability,
    double T,
    double P,
    const std::vector<double>& z,
    int M,
    const double* u_alloc,
    int u_alloc_size,
    const Config::EquilibriumOptimizerConfig& cfg,
    const PhaseSplitOptions& opts = PhaseSplitOptions{})
{
    const int nc = eos.numComponents();
    PhaseSplitEvaluation out;
    if (opts.build_phases) out.phases.reserve(static_cast<size_t>(M));

    const auto st = unpackAllocations(u_alloc, u_alloc_size, nc, M, cfg.composition_min);
    computeBetasAndCompositions(z, st, cfg.beta_min, out.betas, out.x, out.infeasibility);
    const Core::Units::PropertyVariable pressure_tolerance(cfg.density_newton_tol, Core::Units::Unit::PA);

    for (int k = 0; k < M; ++k) {
        const auto ev = DMThermo::Equilibrium::Detail::evalAtTP(
            eos,
            stability,
            T,
            P,
            out.x[static_cast<size_t>(k)],
            hintForPhaseIndex(k),
            cfg.phase_detection_space,
            0.0,
            cfg.max_density_newton_iters,
            pressure_tolerance
        );

        const bool rho_ok = std::isfinite(ev.rho) && (ev.rho > 0.0);
        if (!rho_ok && opts.invalid_rho_penalty > 0.0) {
            out.infeasibility += opts.invalid_rho_penalty;
            if (opts.skip_objective_on_invalid_rho) {
                if (opts.build_phases) {
                    PhaseState ps;
                    ps.type = ev.phase;
                    ps.density = ev.rho;
                    ps.compressibility = 0.0;
                    ps.x = out.x[static_cast<size_t>(k)];
                    ps.phi = ev.phi;
                    ps.fraction = out.betas[static_cast<size_t>(k)];
                    out.phases.push_back(std::move(ps));
                }
                continue;
            }
        }

        out.phase_objective += out.betas[static_cast<size_t>(k)] *
            phaseGibbsRT(out.x[static_cast<size_t>(k)], ev.lnphi, cfg.composition_min);

        if (opts.build_phases) {
            PhaseState ps;
            ps.type = ev.phase;
            ps.density = ev.rho;
            ps.compressibility = rho_ok ? eos.compressibility(T, ps.density, out.x[static_cast<size_t>(k)]) : 0.0;
            ps.x = out.x[static_cast<size_t>(k)];
            ps.phi = ev.phi;
            ps.fraction = out.betas[static_cast<size_t>(k)];
            out.phases.push_back(std::move(ps));
        }
    }

    return out;
}

inline PhaseSplitEvaluation evaluatePhaseSplit(
    const EOS& eos,
    const Stability::IStabilityAnalyzer& stability,
    double T,
    double P,
    const std::vector<double>& z,
    int M,
    const std::vector<double>& u_alloc,
    const Config::EquilibriumOptimizerConfig& cfg,
    const PhaseSplitOptions& opts = PhaseSplitOptions{})
{
    const double* ptr = u_alloc.empty() ? nullptr : u_alloc.data();
    return evaluatePhaseSplit(eos, stability, T, P, z, M, ptr, static_cast<int>(u_alloc.size()), cfg, opts);
}

} // namespace TP
} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_TP_PHASE_SPLIT_H
