/**
 * @file numerics_factory.h
 * @brief Factory for creating numerical solver objects
 */

#ifndef THERMO_FACTORY_NUMERICS_FACTORY_H
#define THERMO_FACTORY_NUMERICS_FACTORY_H

#include <memory>
#include <string>
#include "thermo/numerics/iroot_finder.h"
#include "thermo/numerics/ioptimizer.h"
#include "thermo/numerics/ifixed_point.h"

namespace DMThermo {
namespace Factory {

/**
 * @brief Factory for creating numerical solver implementations
 */
class NumericsFactory {
public:
    // =========================================================================
    // Root Finders
    // =========================================================================

    /**
     * @brief Create default scalar root finder (Newton-Raphson)
     */
    static Numerics::RootFinding::ScalarRootFinderPtr createScalarRootFinder();

    /**
     * @brief Create Newton-Raphson scalar root finder
     */
    static Numerics::RootFinding::ScalarRootFinderPtr createNewtonRaphson();

    /**
     * @brief Create Brent's method scalar root finder
     */
    static Numerics::RootFinding::ScalarRootFinderPtr createBrent();

    /**
     * @brief Create secant method scalar root finder
     */
    static Numerics::RootFinding::ScalarRootFinderPtr createSecant();

    /**
     * @brief Create default vector root finder (Newton multi-D)
     */
    static Numerics::RootFinding::VectorRootFinderPtr createVectorRootFinder();

    /**
     * @brief Create Newton multi-dimensional root finder
     */
    static Numerics::RootFinding::VectorRootFinderPtr createNewtonMultiD();

    // =========================================================================
    // Optimizers
    // =========================================================================

    /**
     * @brief Create default unconstrained optimizer
     */
    static Numerics::Optimization::OptimizerPtr createOptimizer();

    /**
     * @brief Create BFGS optimizer
     */
    static Numerics::Optimization::OptimizerPtr createBFGS();

    /**
     * @brief Create default box-constrained optimizer
     */
    static Numerics::Optimization::BoxConstrainedOptimizerPtr createBoxConstrainedOptimizer();

    // =========================================================================
    // Fixed-Point Solvers
    // =========================================================================

    /**
     * @brief Create default fixed-point solver (GDEM accelerated)
     */
    static Numerics::FixedPoint::FixedPointSolverPtr createFixedPointSolver();

    /**
     * @brief Create direct substitution solver (no acceleration)
     */
    static Numerics::FixedPoint::FixedPointSolverPtr createDirectSubstitution();

    /**
     * @brief Create GDEM accelerated solver
     */
    static Numerics::FixedPoint::FixedPointSolverPtr createGDEM();

    /**
     * @brief Create Wegstein accelerated solver
     */
    static Numerics::FixedPoint::FixedPointSolverPtr createWegstein();
};

} // namespace Factory
} // namespace DMThermo

#endif // THERMO_FACTORY_NUMERICS_FACTORY_H
