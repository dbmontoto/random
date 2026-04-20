/**
 * @file optimization_problem.h
 * @brief Thermo-agnostic interfaces for penalized optimization problems.
 *
 * Many thermo workflows use an unconstrained optimizer (e.g., BFGS) together with:
 *  - variable transforms (to enforce bounds/simplex constraints), and/or
 *  - explicit penalties for infeasibility (constraint violations).
 *
 * This header provides a small, reusable interface that makes those patterns explicit
 * while keeping Numerics independent of any thermo domain types.
 */

#ifndef THERMO_NUMERICS_OPTIMIZATION_PROBLEM_H
#define THERMO_NUMERICS_OPTIMIZATION_PROBLEM_H

#include "thermo/numerics/ioptimizer.h"

#include <cmath>
#include <limits>
#include <vector>

namespace DMThermo {
namespace Numerics {
namespace Optimization {

struct ProblemEvaluation {
    double objective = std::numeric_limits<double>::infinity();
    double infeasibility = 0.0;

    bool isFinite() const { return std::isfinite(objective) && std::isfinite(infeasibility); }

    double penalized(double penalty_weight) const {
        if (!std::isfinite(penalty_weight)) return std::numeric_limits<double>::infinity();
        if (!isFinite()) return std::numeric_limits<double>::infinity();
        return objective + penalty_weight * infeasibility;
    }
};

template <class StateT>
class IOptimizationProblem {
public:
    using State = StateT;
    virtual ~IOptimizationProblem() = default;

    // Evaluate objective terms in the unconstrained optimizer variable space.
    // The problem owns any transforms required to map vars -> physical quantities.
    virtual ProblemEvaluation evaluate(const std::vector<double>& vars) const = 0;

    // Decode the final optimization result into a domain-specific state.
    virtual State decode(const OptimizationResult& opt) const = 0;
};

template <class Problem>
class OptimizationProblemAdapter final : public IOptimizationProblem<typename Problem::State> {
public:
    using State = typename Problem::State;

    explicit OptimizationProblemAdapter(const Problem& problem) : problem_(problem) {}

    ProblemEvaluation evaluate(const std::vector<double>& vars) const override { return problem_.evaluate(vars); }
    State decode(const OptimizationResult& opt) const override { return problem_.decode(opt); }

private:
    const Problem& problem_;
};

template <class State>
struct ProblemRunResult {
    OptimizationResult opt;
    State state{};
    ProblemEvaluation evaluation{};
};

template <class Problem>
inline ProblemRunResult<typename Problem::State> runPenalizedOptimization(
    const IOptimizer& optimizer,
    const Problem& problem,
    const std::vector<double>& vars0,
    double penalty_weight,
    const Config::OptimizerConfig& cfg = Config::OptimizerConfig::defaults())
{
    auto objective = [&](const std::vector<double>& v) -> double {
        try {
            return problem.evaluate(v).penalized(penalty_weight);
        } catch (...) {
            return std::numeric_limits<double>::infinity();
        }
    };

    ProblemRunResult<typename Problem::State> out;
    out.opt = optimizer.minimize(objective, vars0, cfg);
    out.state = problem.decode(out.opt);
    out.evaluation = problem.evaluate(out.opt.x_optimal);
    return out;
}

template <class StateT>
inline ProblemRunResult<StateT> runPenalizedOptimizationErased(
    const IOptimizer& optimizer,
    const IOptimizationProblem<StateT>& problem,
    const std::vector<double>& vars0,
    double penalty_weight,
    const Config::OptimizerConfig& cfg = Config::OptimizerConfig::defaults())
{
    struct Wrapper {
        using State = StateT;
        const IOptimizationProblem<State>& p;

        ProblemEvaluation evaluate(const std::vector<double>& vars) const { return p.evaluate(vars); }
        State decode(const OptimizationResult& opt) const { return p.decode(opt); }
    };

    const Wrapper wrapper{problem};
    return runPenalizedOptimization(optimizer, wrapper, vars0, penalty_weight, cfg);
}

} // namespace Optimization
} // namespace Numerics
} // namespace DMThermo

#endif // THERMO_NUMERICS_OPTIMIZATION_PROBLEM_H

