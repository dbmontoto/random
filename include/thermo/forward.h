/**
 * @file forward.h
 * @brief Forward declarations for the Thermo library
 *
 * Use this header when you only need type declarations without
 * full definitions, to reduce compilation dependencies.
 */

#ifndef THERMO_FORWARD_H
#define THERMO_FORWARD_H

#include <memory>
#include <vector>

namespace DMThermo {

// Core EOS interface (defined directly in DMThermo::)
class EOS;
using EOSPtr = std::shared_ptr<EOS>;
using EOSConstPtr = std::shared_ptr<const EOS>;

// Core types
namespace Core {
    class Component;
    class Mixture;
    class State;
    struct AssociationParams;
    struct PolymerParams;
    struct BinaryParameters;

    using ComponentPtr = std::shared_ptr<const Component>;
    using MixturePtr = std::shared_ptr<const Mixture>;
}

// PC-SAFT implementation
namespace PCSaft {
    class PCSaftEOS;

    using PCSaftEOSPtr = std::shared_ptr<PCSaftEOS>;
}

// Equilibrium
namespace Equilibrium {
    // Generic implementations
    class GibbsFlash;
    class HelmholtzFlash;
    class TPDStability;

    using GibbsFlashPtr = std::shared_ptr<GibbsFlash>;
    using HelmholtzFlashPtr = std::shared_ptr<HelmholtzFlash>;
    using TPDStabilityPtr = std::shared_ptr<TPDStability>;

    // Interfaces (for custom implementations)
    namespace Flash {
        class IFlashSolver;
        struct FlashResult;
        struct PhaseState;

        using FlashSolverPtr = std::shared_ptr<IFlashSolver>;
    }

    namespace Stability {
        class IStabilityAnalyzer;
        struct StabilityResult;
        struct TPDTrialResult;
        struct DensityRootResult;
        struct SpinodalResult;

        using StabilityAnalyzerPtr = std::shared_ptr<IStabilityAnalyzer>;
    }
}

// Numerics
namespace Numerics {
    class DensitySolver;

    using DensitySolverPtr = std::shared_ptr<DensitySolver>;

    namespace RootFinding {
        class IScalarRootFinder;
        class IVectorRootFinder;
        struct ScalarRootResult;
        struct VectorRootResult;

        using ScalarRootFinderPtr = std::shared_ptr<IScalarRootFinder>;
        using VectorRootFinderPtr = std::shared_ptr<IVectorRootFinder>;
    }

    namespace Optimization {
        class IOptimizer;
        class IBoxConstrainedOptimizer;
        struct OptimizationResult;
        struct Bounds;

        using OptimizerPtr = std::shared_ptr<IOptimizer>;
        using BoxConstrainedOptimizerPtr = std::shared_ptr<IBoxConstrainedOptimizer>;
    }

    namespace FixedPoint {
        class IFixedPointSolver;
        struct FixedPointResult;

        using FixedPointSolverPtr = std::shared_ptr<IFixedPointSolver>;
    }
}

// Config
namespace Config {
    struct FlashConfig;
    struct PHFlashConfig;
    struct PSFlashConfig;
    struct TVFlashConfig;
    struct StabilityConfig;
    struct DensityRootConfig;
    struct SpinodalConfig;
    struct RootFinderConfig;
    struct OptimizerConfig;
    struct FixedPointConfig;
}

// Factory
namespace Factory {
    class EOSFactory;
    class FlashFactory;
    class StabilityFactory;
    class NumericsFactory;
}

// Enums (need full definition for use)
enum class PhaseType;
enum class AssociationScheme;

} // namespace DMThermo

#endif // THERMO_FORWARD_H
