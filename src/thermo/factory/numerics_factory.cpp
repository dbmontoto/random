/**
 * @file numerics_factory.cpp
 * @brief Implementation of numerical solver factory.
 */

#include "thermo/factory/numerics_factory.h"
#include "thermo/numerics/iroot_finder.h"
#include "thermo/numerics/ioptimizer.h"
#include "thermo/numerics/ifixed_point.h"

namespace DMThermo {
namespace Numerics {
namespace RootFinding {
std::shared_ptr<IScalarRootFinder> makeNewtonRaphson();
std::shared_ptr<IScalarRootFinder> makeBrent();
std::shared_ptr<IScalarRootFinder> makeSecant();
} // namespace RootFinding
namespace Optimization {
OptimizerPtr makeBFGSOptimizer();
} // namespace Optimization
} // namespace Numerics

namespace Factory {

Numerics::RootFinding::ScalarRootFinderPtr NumericsFactory::createScalarRootFinder() {
    return createNewtonRaphson();
}

Numerics::RootFinding::ScalarRootFinderPtr NumericsFactory::createNewtonRaphson() {
    return Numerics::RootFinding::makeNewtonRaphson();
}

Numerics::RootFinding::ScalarRootFinderPtr NumericsFactory::createBrent() {
    return Numerics::RootFinding::makeBrent();
}

Numerics::RootFinding::ScalarRootFinderPtr NumericsFactory::createSecant() {
    return Numerics::RootFinding::makeSecant();
}

Numerics::RootFinding::VectorRootFinderPtr NumericsFactory::createVectorRootFinder() {
    // Not implemented yet; most equilibrium solves use custom structure.
    return nullptr;
}

Numerics::RootFinding::VectorRootFinderPtr NumericsFactory::createNewtonMultiD() {
    return nullptr;
}

Numerics::Optimization::OptimizerPtr NumericsFactory::createOptimizer() {
    return createBFGS();
}

Numerics::Optimization::OptimizerPtr NumericsFactory::createBFGS() {
    return Numerics::Optimization::makeBFGSOptimizer();
}

Numerics::Optimization::BoxConstrainedOptimizerPtr NumericsFactory::createBoxConstrainedOptimizer() {
    return nullptr;
}

Numerics::FixedPoint::FixedPointSolverPtr NumericsFactory::createFixedPointSolver() {
    return createGDEM();
}

Numerics::FixedPoint::FixedPointSolverPtr NumericsFactory::createDirectSubstitution() {
    return std::make_shared<Numerics::FixedPoint::DirectSubstitution>();
}

Numerics::FixedPoint::FixedPointSolverPtr NumericsFactory::createGDEM() {
    return std::make_shared<Numerics::FixedPoint::GDEMAccelerated>();
}

Numerics::FixedPoint::FixedPointSolverPtr NumericsFactory::createWegstein() {
    return std::make_shared<Numerics::FixedPoint::WegsteinAccelerated>();
}

} // namespace Factory
} // namespace DMThermo

