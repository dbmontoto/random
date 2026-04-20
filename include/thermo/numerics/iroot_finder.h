/**
 * @file iroot_finder.h
 * @brief Abstract interfaces for root finding algorithms
 */

#ifndef THERMO_NUMERICS_IROOT_FINDER_H
#define THERMO_NUMERICS_IROOT_FINDER_H

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include "thermo/config/solver_config.h"

namespace DMThermo {
namespace Numerics {
namespace RootFinding {

/**
 * @brief Result from scalar root finding
 */
struct ScalarRootResult {
    double root = 0.0;
    bool converged = false;
    int iterations = 0;
    double residual = 0.0;
    std::string message;
};

/**
 * @brief Result from multi-dimensional root finding
 */
struct VectorRootResult {
    std::vector<double> root;
    bool converged = false;
    int iterations = 0;
    double residual_norm = 0.0;
    std::string message;
};

/**
 * @brief Abstract interface for scalar root finders
 *
 * Finds x such that f(x) = 0
 */
class IScalarRootFinder {
public:
    virtual ~IScalarRootFinder() = default;

    /**
     * @brief Find root starting from initial guess
     *
     * @param f Function to find root of
     * @param x0 Initial guess
     * @param config Solver configuration
     * @return Root result
     */
    virtual ScalarRootResult solve(
        std::function<double(double)> f,
        double x0,
        const Config::RootFinderConfig& config = Config::RootFinderConfig::defaults()
    ) const = 0;

    /**
     * @brief Find root within brackets [a, b]
     *
     * Requires f(a) and f(b) have opposite signs.
     *
     * @param f Function to find root of
     * @param a Lower bracket
     * @param b Upper bracket
     * @param config Solver configuration
     * @return Root result
     */
    virtual ScalarRootResult solveBracketed(
        std::function<double(double)> f,
        double a,
        double b,
        const Config::RootFinderConfig& config = Config::RootFinderConfig::defaults()
    ) const = 0;

    /**
     * @brief Get algorithm name
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Abstract interface for vector root finders
 *
 * Finds x such that F(x) = 0 (vector equation)
 */
class IVectorRootFinder {
public:
    virtual ~IVectorRootFinder() = default;

    /**
     * @brief Find root of vector function
     *
     * @param F Function from R^n to R^n
     * @param x0 Initial guess
     * @param config Solver configuration
     * @return Root result
     */
    virtual VectorRootResult solve(
        std::function<std::vector<double>(const std::vector<double>&)> F,
        const std::vector<double>& x0,
        const Config::OptimizerConfig& config = Config::OptimizerConfig::defaults()
    ) const = 0;

    /**
     * @brief Find root with user-provided Jacobian
     *
     * @param F Function from R^n to R^n
     * @param J Jacobian function (returns n x n matrix as flat vector)
     * @param x0 Initial guess
     * @param config Solver configuration
     * @return Root result
     */
    virtual VectorRootResult solveWithJacobian(
        std::function<std::vector<double>(const std::vector<double>&)> F,
        std::function<std::vector<double>(const std::vector<double>&)> J,
        const std::vector<double>& x0,
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
using ScalarRootFinderPtr = std::shared_ptr<IScalarRootFinder>;
using VectorRootFinderPtr = std::shared_ptr<IVectorRootFinder>;

} // namespace RootFinding
} // namespace Numerics
} // namespace DMThermo

#endif // THERMO_NUMERICS_IROOT_FINDER_H
