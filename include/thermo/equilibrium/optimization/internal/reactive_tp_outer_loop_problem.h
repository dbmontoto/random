/**
 * @file reactive_tp_outer_loop_problem.h
 * @brief Reactive Gibbs TP outer-loop problem definition (internal).
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_INTERNAL_REACTIVE_TP_OUTER_LOOP_PROBLEM_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_INTERNAL_REACTIVE_TP_OUTER_LOOP_PROBLEM_H

#include "thermo/core/constants.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/equilibrium/optimization/phase_allocations.h"
#include "thermo/equilibrium/optimization/phase_postprocess.h"
#include "thermo/equilibrium/optimization/phase_stability_outer_loop.h"
#include "thermo/equilibrium/optimization/reactive_gibbs_tp_optimizer.h"
#include "thermo/equilibrium/optimization/reactive_tp_problem.h"
#include "thermo/equilibrium/optimization/tp_phase_split.h"
#include "thermo/equilibrium/optimization/variable_blocks.h"
#include "thermo/equilibrium/optimization/variable_layout.h"
#include "thermo/numerics/ioptimizer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {
namespace Internal {

struct ReactiveTPOuterLoopProblem {
    using State = EquilibriumState;

    const EOS& eos;
    const Stability::IStabilityAnalyzer& stability;
    double T = 0.0;
    double P = 0.0;
    const std::vector<double>& n0;
    const Reactions::ReactionSystem& system;
    const std::vector<double>& mu0;
    const ReactiveGibbsTPConfig& cfg;
    const std::vector<Reactions::ExtentBounds>& xi_bounds;
    int nc = 0;
    int nr = 0;
    const std::string& method_used;

    Numerics::Optimization::ProblemEvaluation evaluate(int M, const std::vector<double>& vars) const {
        const ReactiveTPProblem prob{
            eos, stability, T, P, n0, system, mu0, cfg, xi_bounds, nc, nr, M, method_used
        };
        return prob.evaluate(vars);
    }

    double penaltyWeight() const { return cfg.equilibrium.infeasibility_penalty; }

    double score(const EquilibriumState& st) const { return st.score_value; }

    EquilibriumState decode(int M, const Numerics::Optimization::OptimizationResult& opt) const {
        const ReactiveTPProblem prob{
            eos, stability, T, P, n0, system, mu0, cfg, xi_bounds, nc, nr, M, method_used
        };
        return prob.decode(opt);
    }

    std::optional<std::vector<double>> unstableComposition(int, const EquilibriumState& candidate) const {
        return findUnstableComposition(stability, T, P, candidate.phases, cfg.equilibrium);
    }

    std::vector<double> expandInitialGuess(
        int Mold,
        int,
        const Numerics::Optimization::OptimizationResult& opt,
        const EquilibriumState& candidate,
        const std::vector<double>& unstable_w) const
    {
        const auto layout_old = ExtentAllocLayout::make(nr, nc, Mold);
        if (!layout_old.matches(opt.x_optimal)) {
            throw std::runtime_error("ReactiveTPOuterLoopProblem::expandInitialGuess: optimizer variable size mismatch");
        }

        const auto uxi = layout_old.extentSlice(opt.x_optimal);
        const auto ualloc_u = layout_old.allocSlice(opt.x_optimal);

        const AllocState st_old = unpackAllocations(ualloc_u.data, ualloc_u.size, nc, Mold, cfg.equilibrium.composition_min);
        const AllocState st_new = appendPhaseFromComposition(st_old, candidate.z, unstable_w, 1e-3);
        const std::vector<double> u_alloc_new = packLogitsFromAllocations(st_new, cfg.equilibrium.composition_min);

        std::vector<double> vars;
        vars.reserve(static_cast<size_t>(nr) + u_alloc_new.size());
        vars.insert(vars.end(), uxi.data, uxi.data + static_cast<ptrdiff_t>(uxi.size));
        vars.insert(vars.end(), u_alloc_new.begin(), u_alloc_new.end());
        return vars;
    }
};

} // namespace Internal
} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_INTERNAL_REACTIVE_TP_OUTER_LOOP_PROBLEM_H
