/**
 * @file reactive_tp_problem.h
 * @brief Reactive TP optimization problem wrapper (thermo-specific).
 *
 * Encapsulates the "reactive Gibbs at (T,P) + phase split" objective for a fixed phase count M.
 * This wrapper exposes:
 *   - objective terms (physical objective + infeasibility) for penalized optimization, and
 *   - decoding into an EquilibriumState (including reaction extents).
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_REACTIVE_TP_PROBLEM_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_REACTIVE_TP_PROBLEM_H

#include "thermo/core/constants.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/equilibrium/optimization/equilibrium_state.h"
#include "thermo/equilibrium/optimization/phase_postprocess.h"
#include "thermo/equilibrium/optimization/reactive_gibbs_tp_optimizer.h"
#include "thermo/equilibrium/optimization/tp_phase_split.h"
#include "thermo/equilibrium/optimization/variable_blocks.h"
#include "thermo/equilibrium/optimization/variable_layout.h"
#include "thermo/numerics/optimization_problem.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

struct ReactiveTPProblem {
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
    int M = 0;
    const std::string& method_used;

    Numerics::Optimization::ProblemEvaluation evaluate(const std::vector<double>& vars) const {
        const auto layout = ExtentAllocLayout::make(nr, nc, M);
        if (!layout.matches(vars)) {
            return Numerics::Optimization::ProblemEvaluation{
                std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        }

        const auto xi_u = layout.extentSlice(vars);
        const auto ualloc_u = layout.allocSlice(vars);

        const auto xi = decodeExtents(xi_u, xi_bounds);

        const auto n = system.molesFromExtents(n0, xi);
        double N = 0.0;
        double infeas = 0.0;
        for (double ni : n) {
            if (!std::isfinite(ni)) {
                return Numerics::Optimization::ProblemEvaluation{
                    std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
            }
            if (ni < cfg.min_moles) {
                const double d = cfg.min_moles - ni;
                infeas += d * d;
            }
            N += std::max(ni, cfg.min_moles);
        }
        if (!(std::isfinite(N) && N > 0.0)) {
            return Numerics::Optimization::ProblemEvaluation{
                std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        }

        std::vector<double> z(static_cast<size_t>(nc), 0.0);
        for (int i = 0; i < nc; ++i) z[static_cast<size_t>(i)] = std::max(n[static_cast<size_t>(i)], cfg.min_moles) / N;
        z = clampAndNormalize(z, Constants::MIN_MOLE_FRACTION);

        const double invRT = 1.0 / (Constants::GAS_CONSTANT * T);
        const double lnP = safeLog(P / cfg.standard_pressure, cfg.equilibrium.composition_min);
        double mu_term = 0.0;
        for (int i = 0; i < nc; ++i) {
            mu_term += z[static_cast<size_t>(i)] * (mu0[static_cast<size_t>(i)] * invRT);
        }
        double obj = N * (lnP + mu_term);

        const auto pe = TP::evaluatePhaseSplit(
            eos,
            stability,
            T,
            P,
            z,
            M,
            ualloc_u.data,
            ualloc_u.size,
            cfg.equilibrium,
            TP::PhaseSplitOptions{false, 1.0, true}
        );
        obj += N * pe.phase_objective;
        infeas += pe.infeasibility;

        if (!std::isfinite(obj)) {
            return Numerics::Optimization::ProblemEvaluation{
                std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
        }
        return Numerics::Optimization::ProblemEvaluation{obj, infeas};
    }

    EquilibriumState decode(const Numerics::Optimization::OptimizationResult& opt) const {
        const auto layout = ExtentAllocLayout::make(nr, nc, M);
        if (!layout.matches(opt.x_optimal)) {
            throw std::runtime_error("ReactiveTPProblem::decode: optimizer variable size mismatch");
        }
        const auto uxi_u = layout.extentSlice(opt.x_optimal);
        const auto ualloc_u = layout.allocSlice(opt.x_optimal);

        auto xi = decodeExtents(uxi_u, xi_bounds);
        const auto n = system.molesFromExtents(n0, xi);
        double N = 0.0;
        double infeas = 0.0;
        for (double ni : n) {
            if (!std::isfinite(ni)) {
                throw std::runtime_error("ReactiveTPProblem::decode: non-finite moles from extents");
            }
            if (ni < cfg.min_moles) {
                const double d = cfg.min_moles - ni;
                infeas += d * d;
            }
            N += std::max(ni, cfg.min_moles);
        }
        if (!(std::isfinite(N) && N > 0.0)) {
            throw std::runtime_error("ReactiveTPProblem::decode: invalid total moles from extents");
        }

        std::vector<double> z(static_cast<size_t>(nc), 0.0);
        for (int i = 0; i < nc; ++i) z[static_cast<size_t>(i)] = std::max(n[static_cast<size_t>(i)], cfg.min_moles) / N;
        z = clampAndNormalize(z, Constants::MIN_MOLE_FRACTION);

        const double invRT = 1.0 / (Constants::GAS_CONSTANT * T);
        const double lnP = safeLog(P / cfg.standard_pressure, cfg.equilibrium.composition_min);
        double mu_term = 0.0;
        for (int i = 0; i < nc; ++i) {
            mu_term += z[static_cast<size_t>(i)] * (mu0[static_cast<size_t>(i)] * invRT);
        }

        const auto pe = TP::evaluatePhaseSplit(
            eos,
            stability,
            T,
            P,
            z,
            M,
            ualloc_u.data,
            ualloc_u.size,
            cfg.equilibrium,
            TP::PhaseSplitOptions{true}
        );
        std::vector<PhaseState> phases = pe.phases;
        const double obj = N * (lnP + mu_term) + N * pe.phase_objective;
        infeas += pe.infeasibility;

        PhasePostProcessOptions pp;
        pp.prune_fraction = cfg.equilibrium.beta_min * 10.0;
        pp.allow_empty = true;
        pp.merge_duplicates = false;
        pp.sort_by_type = false;
        pp.canonicalize_liquid = false;
        (void)postProcessPhases(phases, pp);

        EquilibriumState st;
        st.temperature = T;
        st.pressure = P;
        st.z = std::move(z);
        st.phases = std::move(phases);
        st.reaction_extents = std::move(xi);
        st.iterations = opt.iterations;

        Numerics::Optimization::ProblemEvaluation eval;
        eval.objective = obj;
        eval.infeasibility = infeas;
        st.objective_value = eval.objective;
        st.infeasibility = eval.infeasibility;
        st.score_value = eval.penalized(cfg.equilibrium.infeasibility_penalty);
        st.converged = opt.converged;
        st.method_used = method_used;
        st.message = opt.message;
        return st;
    }
};

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_REACTIVE_TP_PROBLEM_H
