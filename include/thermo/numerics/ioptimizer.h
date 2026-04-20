/**
 * @file ioptimizer.h
 * @brief Abstract interfaces for optimization algorithms
 */

#ifndef THERMO_NUMERICS_IOPTIMIZER_H
#define THERMO_NUMERICS_IOPTIMIZER_H

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <optional>
#include "thermo/config/solver_config.h"

namespace DMThermo {
namespace Numerics {
namespace Optimization {

/**
 * @brief Result from optimization
 */
struct OptimizationResult {
    /// Optimal point
    std::vector<double> x_optimal;

    /// Optimal function value
    double f_optimal = 0.0;

    /// Did optimization converge?
    bool converged = false;

    /// Number of iterations
    int iterations = 0;

    /// Number of function evaluations
    int function_evaluations = 0;

    /// Final gradient norm (for unconstrained)
    double gradient_norm = 0.0;

    /// Status message
    std::string message;
};

/**
 * @brief Abstract interface for unconstrained optimizers
 *
 * Minimizes f(x) over R^n
 */
class IOptimizer {
public:
    virtual ~IOptimizer() = default;

    /**
     * @brief Minimize function starting from initial point
     *
     * @param f Objective function
     * @param x0 Initial point
     * @param config Solver configuration
     * @return Optimization result
     */
    virtual OptimizationResult minimize(
        std::function<double(const std::vector<double>&)> f,
        const std::vector<double>& x0,
        const Config::OptimizerConfig& config = Config::OptimizerConfig::defaults()
    ) const = 0;

    /**
     * @brief Minimize with user-provided gradient
     *
     * @param f Objective function
     * @param grad Gradient function
     * @param x0 Initial point
     * @param config Solver configuration
     * @return Optimization result
     */
    virtual OptimizationResult minimizeWithGradient(
        std::function<double(const std::vector<double>&)> f,
        std::function<std::vector<double>(const std::vector<double>&)> grad,
        const std::vector<double>& x0,
        const Config::OptimizerConfig& config = Config::OptimizerConfig::defaults()
    ) const = 0;

    /**
     * @brief Get algorithm name
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Bounds for box-constrained optimization
 */
struct Bounds {
    std::vector<double> lower;
    std::vector<double> upper;

    static Bounds none(int n) {
        return Bounds{
            std::vector<double>(n, -1e30),
            std::vector<double>(n, 1e30)
        };
    }

    static Bounds positive(int n) {
        return Bounds{
            std::vector<double>(n, 0.0),
            std::vector<double>(n, 1e30)
        };
    }

    static Bounds unit(int n) {
        return Bounds{
            std::vector<double>(n, 0.0),
            std::vector<double>(n, 1.0)
        };
    }
};

/**
 * @brief Abstract interface for box-constrained optimizers
 *
 * Minimizes f(x) subject to lower <= x <= upper
 */
class IBoxConstrainedOptimizer {
public:
    virtual ~IBoxConstrainedOptimizer() = default;

    /**
     * @brief Minimize function with box constraints
     *
     * @param f Objective function
     * @param x0 Initial point (must be feasible)
     * @param bounds Lower and upper bounds
     * @param config Solver configuration
     * @return Optimization result
     */
    virtual OptimizationResult minimize(
        std::function<double(const std::vector<double>&)> f,
        const std::vector<double>& x0,
        const Bounds& bounds,
        const Config::OptimizerConfig& config = Config::OptimizerConfig::defaults()
    ) const = 0;

    /**
     * @brief Get algorithm name
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Shared pointer type aliases
 */
using OptimizerPtr = std::shared_ptr<IOptimizer>;
using BoxConstrainedOptimizerPtr = std::shared_ptr<IBoxConstrainedOptimizer>;

} // namespace Optimization
} // namespace Numerics
} // namespace DMThermo

#endif // THERMO_NUMERICS_IOPTIMIZER_H
