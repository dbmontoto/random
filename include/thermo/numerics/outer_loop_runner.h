/**
 * @file outer_loop_runner.h
 * @brief Generic outer-loop runner for "optimize + adapt + repeat" workflows.
 *
 * This is intentionally thermo-agnostic. A caller supplies a Problem that can:
 *   - evaluate objective terms (objective + infeasibility) for a given "phase count" M and variable vector
 *   - provide a penalty weight used to form a penalized objective
 *   - decode an optimizer result into a domain-specific State
 *   - score a decoded State to select the best candidate
 *   - decide whether to expand the problem (e.g., add a phase) and produce a new initial guess
 */

#ifndef THERMO_NUMERICS_OUTER_LOOP_RUNNER_H
#define THERMO_NUMERICS_OUTER_LOOP_RUNNER_H

#include "thermo/numerics/ioptimizer.h"
#include "thermo/numerics/optimization_problem.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace DMThermo {
namespace Numerics {
namespace Optimization {

struct OuterLoopConfig {
    bool enabled = true;
    int max_outer_iterations = 6;
    Config::OptimizerConfig optimizer = Config::OptimizerConfig::defaults();
};

/**
 * @brief Explicit "outer-loop problem" interface for type-erased use cases.
 *
 * The existing template runner remains available and is typically preferred when
 * concrete problem types are known at compile time. This interface enables
 * runtime selection/composition of problems without forcing the runner itself
 * to be virtual.
 */
template <class StateT>
class IOuterLoopProblem {
public:
    using State = StateT;
    virtual ~IOuterLoopProblem() = default;

    virtual ProblemEvaluation evaluate(int phases, const std::vector<double>& vars) const = 0;
    virtual double penaltyWeight() const = 0;
    virtual State decode(int phases, const OptimizationResult& opt) const = 0;
    virtual double score(const State& st) const = 0;
    virtual std::optional<std::vector<double>> unstableComposition(int phases, const State& st) const = 0;
    virtual std::vector<double> expandInitialGuess(
        int phases_old,
        int phases_new,
        const OptimizationResult& opt,
        const State& st,
        const std::vector<double>& unstable_w) const = 0;
};

/**
 * @brief Adapts a concrete Problem type (with the template runner interface) to IOuterLoopProblem.
 */
template <class Problem>
class OuterLoopProblemAdapter final : public IOuterLoopProblem<typename Problem::State> {
public:
    using State = typename Problem::State;

    explicit OuterLoopProblemAdapter(const Problem& problem) : problem_(problem) {}

    ProblemEvaluation evaluate(int phases, const std::vector<double>& vars) const override { return problem_.evaluate(phases, vars); }
    double penaltyWeight() const override { return problem_.penaltyWeight(); }
    State decode(int phases, const OptimizationResult& opt) const override { return problem_.decode(phases, opt); }
    double score(const State& st) const override { return problem_.score(st); }
    std::optional<std::vector<double>> unstableComposition(int phases, const State& st) const override {
        return problem_.unstableComposition(phases, st);
    }
    std::vector<double> expandInitialGuess(
        int phases_old,
        int phases_new,
        const OptimizationResult& opt,
        const State& st,
        const std::vector<double>& unstable_w) const override
    {
        return problem_.expandInitialGuess(phases_old, phases_new, opt, st, unstable_w);
    }

private:
    const Problem& problem_;
};

template <class State>
struct OuterLoopRunResult {
    State best_state{};
    double best_score = std::numeric_limits<double>::infinity();
    int phases_final = 0;
    int outer_iterations = 0;
};

namespace Detail {

template <class Problem>
inline OuterLoopRunResult<typename Problem::State> runPhaseStabilityOuterLoopImpl(
    const IOptimizer& optimizer,
    Problem& problem,
    int phases_initial,
    int phases_max,
    std::vector<double> vars_initial,
    const OuterLoopConfig& cfg)
{
    using State = typename Problem::State;

    OuterLoopRunResult<State> out;

    const int max_phases_clamped = std::max(1, phases_max);
    int M = std::clamp(phases_initial, 1, max_phases_clamped);
    std::vector<double> vars = std::move(vars_initial);

    const int max_outer = std::max(1, cfg.max_outer_iterations);
    bool have_best = false;

    for (int outer = 0; outer < max_outer; ++outer) {
        out.outer_iterations = outer + 1;

        auto objective = [&](const std::vector<double>& v) -> double {
            try {
                return problem.evaluate(M, v).penalized(problem.penaltyWeight());
            } catch (...) {
                return std::numeric_limits<double>::infinity();
            }
        };

        const auto opt = optimizer.minimize(objective, vars, cfg.optimizer);
        const State candidate = problem.decode(M, opt);

        const double score = problem.score(candidate);
        if (!have_best) {
            out.best_state = candidate;
            out.best_score = score;
            have_best = true;
        } else {
            const bool score_ok = std::isfinite(score);
            const bool best_ok = std::isfinite(out.best_score);
            if (score_ok && (!best_ok || score < out.best_score)) {
                out.best_state = candidate;
                out.best_score = score;
            }
        }

        if (!cfg.enabled || M >= phases_max) {
            break;
        }

        const auto unstable = problem.unstableComposition(M, candidate);
        if (!unstable.has_value()) {
            break;
        }

        const int Mold = M;
        M = std::min(M + 1, phases_max);
        vars = problem.expandInitialGuess(Mold, M, opt, candidate, unstable.value());
    }

    out.phases_final = M;
    return out;
}

} // namespace Detail

template <class Problem>
inline OuterLoopRunResult<typename Problem::State> runPhaseStabilityOuterLoop(
    const IOptimizer& optimizer,
    Problem& problem,
    int phases_initial,
    int phases_max,
    std::vector<double> vars_initial,
    const OuterLoopConfig& cfg)
{
    return Detail::runPhaseStabilityOuterLoopImpl(
        optimizer, problem, phases_initial, phases_max, std::move(vars_initial), cfg);
}

/**
 * @brief Type-erased runner overload using IOuterLoopProblem.
 */
template <class StateT>
inline OuterLoopRunResult<StateT> runPhaseStabilityOuterLoopErased(
    const IOptimizer& optimizer,
    const IOuterLoopProblem<StateT>& problem,
    int phases_initial,
    int phases_max,
    std::vector<double> vars_initial,
    const OuterLoopConfig& cfg)
{
    struct Wrapper {
        using State = StateT;
        const IOuterLoopProblem<State>& p;

        ProblemEvaluation evaluate(int phases, const std::vector<double>& vars) const { return p.evaluate(phases, vars); }
        double penaltyWeight() const { return p.penaltyWeight(); }
        State decode(int phases, const OptimizationResult& opt) const { return p.decode(phases, opt); }
        double score(const State& st) const { return p.score(st); }
        std::optional<std::vector<double>> unstableComposition(int phases, const State& st) const {
            return p.unstableComposition(phases, st);
        }
        std::vector<double> expandInitialGuess(
            int phases_old,
            int phases_new,
            const OptimizationResult& opt,
            const State& st,
            const std::vector<double>& unstable_w) const
        {
            return p.expandInitialGuess(phases_old, phases_new, opt, st, unstable_w);
        }
    };

    Wrapper wrapper{problem};
    return Detail::runPhaseStabilityOuterLoopImpl(
        optimizer, wrapper, phases_initial, phases_max, std::move(vars_initial), cfg);
}

} // namespace Optimization
} // namespace Numerics
} // namespace DMThermo

#endif // THERMO_NUMERICS_OUTER_LOOP_RUNNER_H
