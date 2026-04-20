/**
 * @file gibbs_tp_outer_loop_problem.h
 * @brief Gibbs TP outer-loop problem definition (internal).
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_INTERNAL_GIBBS_TP_OUTER_LOOP_PROBLEM_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_INTERNAL_GIBBS_TP_OUTER_LOOP_PROBLEM_H

#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/equilibrium/optimization/phase_allocations.h"
#include "thermo/equilibrium/optimization/phase_postprocess.h"
#include "thermo/equilibrium/optimization/phase_stability_outer_loop.h"
#include "thermo/equilibrium/optimization/tp_phase_split_problem.h"
#include "thermo/equilibrium/optimization/tp_phase_split.h"
#include "thermo/numerics/ioptimizer.h"
#include "thermo/numerics/optimization_problem.h"

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {
namespace Internal {

struct GibbsTPOuterLoopProblem {
    using State = EquilibriumState;

    const EOS& eos;
    const Stability::IStabilityAnalyzer& stability;
    double T = 0.0;
    double P = 0.0;
    const std::vector<double>& z;
    int nc = 0;
    const Config::EquilibriumOptimizerConfig& cfg;
    const std::string& method_used;

    Numerics::Optimization::ProblemEvaluation evaluate(int M, const std::vector<double>& vars) const {
        const TP::TPPhaseSplitProblem prob{eos, stability, T, P, z, nc, M, cfg, method_used};
        return prob.evaluate(vars);
    }

    double penaltyWeight() const { return cfg.infeasibility_penalty; }

    double score(const EquilibriumState& st) const { return st.score_value; }

    EquilibriumState decode(int M, const Numerics::Optimization::OptimizationResult& opt) const {
        const TP::TPPhaseSplitProblem prob{eos, stability, T, P, z, nc, M, cfg, method_used};
        return prob.decode(opt);
    }

    std::optional<std::vector<double>> unstableComposition(int, const EquilibriumState& candidate) const {
        return findUnstableComposition(stability, T, P, candidate.phases, cfg);
    }

    std::vector<double> expandInitialGuess(
        int Mold,
        int,
        const Numerics::Optimization::OptimizationResult& opt,
        const EquilibriumState& candidate,
        const std::vector<double>& unstable_w) const
    {
        const AllocState st_old = unpackAllocations(opt.x_optimal, nc, Mold, cfg.composition_min);
        const AllocState st_new = appendPhaseFromComposition(st_old, candidate.z, unstable_w, 1e-3);
        return packLogitsFromAllocations(st_new, cfg.composition_min);
    }
};

} // namespace Internal
} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_INTERNAL_GIBBS_TP_OUTER_LOOP_PROBLEM_H
