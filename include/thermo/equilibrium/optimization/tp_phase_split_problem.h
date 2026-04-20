/**
 * @file tp_phase_split_problem.h
 * @brief TP phase-split optimization problem wrapper (thermo-specific).
 *
 * Provides a reusable "problem" object for a fixed (T, P, z, M) that:
 *   - exposes objective terms (phase objective + infeasibility)
 *   - decodes optimizer results into an EquilibriumState with phase post-processing
 *
 * This is intended to be used by multiple workflows (flash, stability outer loops,
 * saturation, reactions) while keeping Numerics thermo-agnostic.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_TP_PHASE_SPLIT_PROBLEM_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_TP_PHASE_SPLIT_PROBLEM_H

#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/equilibrium/optimization/equilibrium_state.h"
#include "thermo/equilibrium/optimization/phase_postprocess.h"
#include "thermo/equilibrium/optimization/tp_phase_split.h"
#include "thermo/numerics/optimization_problem.h"

#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {
namespace TP {

struct TPPhaseSplitProblem {
    using State = EquilibriumState;

    const EOS& eos;
    const Stability::IStabilityAnalyzer& stability;
    double T = 0.0;
    double P = 0.0;
    const std::vector<double>& z;
    int nc = 0;
    int M = 0;
    const Config::EquilibriumOptimizerConfig& cfg;
    const std::string& method_used;

    Numerics::Optimization::ProblemEvaluation evaluate(const std::vector<double>& vars) const {
        const auto e = TP::evaluatePhaseSplit(
            eos, stability, T, P, z, M, vars, cfg, TP::PhaseSplitOptions{});
        Numerics::Optimization::ProblemEvaluation out;
        out.objective = e.phase_objective;
        out.infeasibility = e.infeasibility;
        return out;
    }

    EquilibriumState decode(const Numerics::Optimization::OptimizationResult& opt) const {
        const auto e = TP::evaluatePhaseSplit(
            eos, stability, T, P, z, M, opt.x_optimal, cfg, TP::PhaseSplitOptions{true});
        std::vector<PhaseState> phases = e.phases;

        Numerics::Optimization::ProblemEvaluation eval;
        eval.objective = e.phase_objective;
        eval.infeasibility = e.infeasibility;

        PhasePostProcessOptions pp;
        pp.prune_fraction = cfg.beta_min * 10.0;
        pp.allow_empty = false;
        pp.merge_duplicates = false;
        pp.sort_by_type = true;
        pp.canonicalize_liquid = true;
        const bool ok = postProcessPhases(phases, pp);
        if (!ok) {
            // Avoid throwing from a decode path that is frequently used inside outer solves.
            // If the candidate collapses, fall back to a best-effort single-phase state at z.
            const Core::Units::PropertyVariable pressure_tolerance(cfg.density_newton_tol, Core::Units::Unit::PA);
            const auto ev = Detail::evalAtTP(
                eos, stability, T, P, z, PhaseType::Unknown,
                cfg.phase_detection_space, 0.0, cfg.max_density_newton_iters, pressure_tolerance
            );
            PhaseState ps;
            ps.type = ev.phase;
            ps.density = ev.rho;
            ps.compressibility = eos.compressibility(T, ps.density, z);
            ps.x = z;
            ps.phi = ev.phi;
            ps.fraction = 1.0;
            phases = {ps};

            eval.objective = phaseGibbsRT(ps.x, ev.lnphi, cfg.composition_min);
            eval.infeasibility = std::numeric_limits<double>::infinity();
        }

        EquilibriumState candidate;
        candidate.temperature = T;
        candidate.pressure = P;
        candidate.z = z;
        candidate.phases = std::move(phases);
        candidate.iterations = opt.iterations;
        candidate.objective_value = eval.objective;
        candidate.infeasibility = eval.infeasibility;
        candidate.score_value = eval.penalized(cfg.infeasibility_penalty);
        candidate.converged = ok && opt.converged;
        candidate.method_used = method_used;
        candidate.message = ok ? opt.message : (opt.message + " (phase collapse -> single-phase fallback)");
        return candidate;
    }
};

} // namespace TP
} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_TP_PHASE_SPLIT_PROBLEM_H
