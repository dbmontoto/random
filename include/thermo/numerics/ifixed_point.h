/**
 * @file ifixed_point.h
 * @brief Abstract interfaces for fixed-point iteration solvers
 */

#ifndef THERMO_NUMERICS_IFIXED_POINT_H
#define THERMO_NUMERICS_IFIXED_POINT_H

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include "thermo/config/solver_config.h"

namespace DMThermo {
namespace Numerics {
namespace FixedPoint {

/**
 * @brief Result from fixed-point iteration
 */
struct FixedPointResult {
    /// Converged fixed point
    std::vector<double> x;

    /// Did iteration converge?
    bool converged = false;

    /// Number of iterations
    int iterations = 0;

    /// Final residual (||x_new - x_old||)
    double residual = 0.0;

    /// Status message
    std::string message;

    /// Iteration history (if verbose)
    std::vector<double> residual_history;
};

/**
 * @brief Abstract interface for fixed-point iteration solvers
 *
 * Finds x such that x = g(x) (fixed point of g)
 */
class IFixedPointSolver {
public:
    virtual ~IFixedPointSolver() = default;

    /**
     * @brief Solve fixed-point problem x = g(x)
     *
     * @param g Fixed-point function
     * @param x0 Initial guess
     * @param config Solver configuration
     * @return Fixed-point result
     */
    virtual FixedPointResult solve(
        std::function<std::vector<double>(const std::vector<double>&)> g,
        const std::vector<double>& x0,
        const Config::FixedPointConfig& config = Config::FixedPointConfig::defaults()
    ) const = 0;

    /**
     * @brief Get algorithm name
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Direct substitution solver (no acceleration)
 */
class DirectSubstitution : public IFixedPointSolver {
public:
    FixedPointResult solve(
        std::function<std::vector<double>(const std::vector<double>&)> g,
        const std::vector<double>& x0,
        const Config::FixedPointConfig& config = Config::FixedPointConfig::defaults()
    ) const override;

    std::string name() const override { return "DirectSubstitution"; }
};

/**
 * @brief GDEM (General Dominant Eigenvalue Method) accelerated solver
 */
class GDEMAccelerated : public IFixedPointSolver {
public:
    FixedPointResult solve(
        std::function<std::vector<double>(const std::vector<double>&)> g,
        const std::vector<double>& x0,
        const Config::FixedPointConfig& config = Config::FixedPointConfig::defaults()
    ) const override;

    std::string name() const override { return "GDEM"; }
};

/**
 * @brief Wegstein accelerated solver
 */
class WegsteinAccelerated : public IFixedPointSolver {
public:
    FixedPointResult solve(
        std::function<std::vector<double>(const std::vector<double>&)> g,
        const std::vector<double>& x0,
        const Config::FixedPointConfig& config = Config::FixedPointConfig::defaults()
    ) const override;

    std::string name() const override { return "Wegstein"; }
};

/**
 * @brief Shared pointer type alias
 */
using FixedPointSolverPtr = std::shared_ptr<IFixedPointSolver>;

} // namespace FixedPoint
} // namespace Numerics
} // namespace DMThermo

#endif // THERMO_NUMERICS_IFIXED_POINT_H
