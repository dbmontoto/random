/**
 * @file istability.h
 * @brief Abstract interface for phase stability analyzers
 *
 * This interface defines the contract for stability analysis implementations.
 * Stability analysis determines if a phase will split at given conditions.
 */

#ifndef THERMO_EQUILIBRIUM_ISTABILITY_H
#define THERMO_EQUILIBRIUM_ISTABILITY_H

#include <memory>
#include <string>
#include <vector>
#include "stability_result.h"
#include "thermo/config/stability_config.h"
#include "thermo/eos.h"

namespace DMThermo {
namespace Equilibrium {
namespace Stability {

/**
 * @brief Abstract interface for phase stability analyzers
 *
 * All methods are const to ensure thread safety.
 *
 * Example usage:
 * @code
 * auto stability = Factory::StabilityFactory::create(eos);
 * StabilityResult result = stability->analyze(300.0, 1e5, {0.5, 0.5});
 * if (!result.is_stable) {
 *     // Phase will split - use incipient composition as initial guess
 *     auto w = result.unstable_composition.value();
 * }
 * @endcode
 */
class IStabilityAnalyzer {
public:
    virtual ~IStabilityAnalyzer() = default;

    // =========================================================================
    // Core Stability Analysis (Thread-Safe)
    // =========================================================================

    /**
     * @brief Complete stability analysis at given conditions
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param z Feed composition (mole fractions)
     * @param config Analysis configuration
     * @return StabilityResult with stability determination
     */
    virtual StabilityResult analyze(
        double T,
        double P,
        const std::vector<double>& z,
        const Config::StabilityConfig& config = Config::StabilityConfig::defaults()
    ) const = 0;

    /**
     * @brief Quick stability check (faster, less thorough)
     *
     * Uses fewer trials and simpler algorithm for quick screening.
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param z Feed composition
     * @return true if likely stable, false if definitely unstable
     */
    virtual bool isStable(
        double T,
        double P,
        const std::vector<double>& z
    ) const = 0;

    // =========================================================================
    // TPD (Tangent Plane Distance) Methods
    // =========================================================================

    /**
     * @brief Calculate TPD for a specific trial composition
     *
     * TPD = sum_i w_i * [mu_i(w) - mu_i(z)]
     * Negative TPD indicates instability.
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param z Feed composition
     * @param w Trial composition
     * @return TPD value (negative = unstable)
     */
    virtual double tangentPlaneDistance(
        double T,
        double P,
        const std::vector<double>& z,
        const std::vector<double>& w
    ) const = 0;

    /**
     * @brief Minimize TPD starting from given trial composition
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param z Feed composition
     * @param w_initial Initial trial composition
     * @param config Configuration
     * @return TPDTrialResult with converged composition and TPD
     */
    virtual TPDTrialResult minimizeTPD(
        double T,
        double P,
        const std::vector<double>& z,
        const std::vector<double>& w_initial,
        const Config::StabilityConfig& config = Config::StabilityConfig::defaults()
    ) const = 0;

    // =========================================================================
    // Density Root Methods
    // =========================================================================

    /**
     * @brief Find all density roots at given T, P, composition
     *
     * Finds vapor, metastable, and liquid roots if they exist.
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param x Composition
     * @param config Configuration
     * @return DensityRootResult with all roots
     */
    virtual DensityRootResult findDensityRoots(
        double T,
        double P,
        const std::vector<double>& x,
        const Config::DensityRootConfig& config = Config::DensityRootConfig{}
    ) const = 0;

    /**
     * @brief Find vapor density root
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param x Composition
     * @return Vapor density [mol/m³], or 0 if not found
     */
    virtual double findVaporDensity(
        double T,
        double P,
        const std::vector<double>& x
    ) const = 0;

    /**
     * @brief Find liquid density root
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param x Composition
     * @return Liquid density [mol/m³], or 0 if not found
     */
    virtual double findLiquidDensity(
        double T,
        double P,
        const std::vector<double>& x
    ) const = 0;

    // =========================================================================
    // Spinodal Methods
    // =========================================================================

    /**
     * @brief Detect spinodal points at given temperature
     *
     * Spinodal points are where dP/drho = 0.
     *
     * @param T Temperature [K]
     * @param x Composition
     * @param config Configuration
     * @return SpinodalResult with spinodal densities
     */
    virtual SpinodalResult findSpinodalPoints(
        double T,
        const std::vector<double>& x,
        const Config::SpinodalConfig& config = Config::SpinodalConfig{}
    ) const = 0;

    /**
     * @brief Check if density is within spinodal region
     *
     * @param T Temperature [K]
     * @param rho Density [mol/m³]
     * @param x Composition
     * @return true if in spinodal (absolutely unstable) region
     */
    virtual bool isInSpinodalRegion(
        double T,
        double rho,
        const std::vector<double>& x
    ) const = 0;

    // =========================================================================
    // Phase Classification
    // =========================================================================

    /**
     * @brief Classify phase type based on density
     *
     * @param T Temperature [K]
     * @param rho Density [mol/m³]
     * @param x Composition
     * @return Phase type classification
     */
    virtual PhaseType classifyPhase(
        double T,
        double rho,
        const std::vector<double>& x
    ) const = 0;

    // =========================================================================
    // Trial Composition Generation
    // =========================================================================

    /**
     * @brief Generate trial compositions for stability analysis
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param z Feed composition
     * @param num_trials Number of trials to generate
     * @return Vector of trial compositions
     */
    virtual std::vector<std::vector<double>> generateTrialCompositions(
        double T,
        double P,
        const std::vector<double>& z,
        int num_trials
    ) const = 0;

    // =========================================================================
    // Metadata
    // =========================================================================

    /**
     * @brief Get underlying EOS
     * @return Shared pointer to the EOS
     */
    virtual EOSPtr eos() const = 0;
};

/**
 * @brief Shared pointer type alias for IStabilityAnalyzer
 */
using StabilityAnalyzerPtr = std::shared_ptr<IStabilityAnalyzer>;

} // namespace Stability
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_ISTABILITY_H
