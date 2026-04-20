/**
 * @file iflash.h
 * @brief Abstract interface for flash solvers
 *
 * This interface defines the contract for all flash calculation implementations.
 * Flash calculations determine phase equilibrium at specified conditions.
 */

#ifndef THERMO_EQUILIBRIUM_IFLASH_H
#define THERMO_EQUILIBRIUM_IFLASH_H

#include <memory>
#include <string>
#include <vector>
#include "flash_result.h"
#include "thermo/config/flash_config.h"
#include "thermo/eos.h"

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

/**
 * @brief Abstract interface for flash solvers
 *
 * All methods are const to ensure thread safety. Implementations
 * should be stateless except for the EOS reference.
 *
 * Example usage:
 * @code
 * auto flash = Factory::FlashFactory::createGibbs(eos);
 * FlashResult result = flash->solvePT(300.0, 1e5, {0.5, 0.5});
 * if (result.converged && result.is_two_phase) {
 *     double beta = result.vapor_fraction;
 * }
 * @endcode
 */
class IFlashSolver {
public:
    virtual ~IFlashSolver() = default;

    // =========================================================================
    // Core Flash Methods (Thread-Safe)
    // =========================================================================

    /**
     * @brief PT Flash: Given T and P, find equilibrium phases
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param z Overall mole fractions
     * @param config Solver configuration
     * @return FlashResult with phase compositions and properties
     */
    virtual FlashResult solvePT(
        double T,
        double P,
        const std::vector<double>& z,
        const Config::FlashConfig& config = Config::FlashConfig::defaults()
    ) const = 0;

    /**
     * @brief TV Flash: Given T and V, find equilibrium phases
     *
     * Uses Helmholtz minimization directly since V is natural variable.
     *
     * @param T Temperature [K]
     * @param V Molar volume [m³/mol]
     * @param z Overall mole fractions
     * @param config Solver configuration
     * @return FlashResult with phase compositions and properties
     */
    virtual FlashResult solveTV(
        double T,
        double V,
        const std::vector<double>& z,
        const Config::TVFlashConfig& config = Config::TVFlashConfig{}
    ) const = 0;

    /**
     * @brief PH Flash: Given P and H, find T and equilibrium phases
     *
     * @param P Pressure [Pa]
     * @param H Molar enthalpy [J/mol]
     * @param z Overall mole fractions
     * @param config Solver configuration
     * @return FlashResult with phase compositions and properties
     */
    virtual FlashResult solvePH(
        double P,
        double H,
        const std::vector<double>& z,
        const Config::PHFlashConfig& config = Config::PHFlashConfig{}
    ) const = 0;

    /**
     * @brief PS Flash: Given P and S, find T and equilibrium phases
     *
     * @param P Pressure [Pa]
     * @param S Molar entropy [J/(mol·K)]
     * @param z Overall mole fractions
     * @param config Solver configuration
     * @return FlashResult with phase compositions and properties
     */
    virtual FlashResult solvePS(
        double P,
        double S,
        const std::vector<double>& z,
        const Config::PSFlashConfig& config = Config::PSFlashConfig{}
    ) const = 0;

    // =========================================================================
    // Multi-Phase Methods
    // =========================================================================

    /**
     * @brief Three-phase (VLLE) flash calculation
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param z Overall mole fractions
     * @param config Solver configuration
     * @return FlashResult with up to 3 phases
     */
    virtual FlashResult solveVLLE(
        double T,
        double P,
        const std::vector<double>& z,
        const Config::FlashConfig& config = Config::FlashConfig::threePhase()
    ) const = 0;

    // =========================================================================
    // Utility Methods
    // =========================================================================

    /**
     * @brief Get initial K-value estimates
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @return Vector of K-values from Wilson correlation
     */
    virtual std::vector<double> estimateKValues(
        double T,
        double P
    ) const = 0;

    /**
     * @brief Check if flash converged to valid result
     *
     * @param result Flash result to validate
     * @return true if result is physically valid
     */
    virtual bool isValidResult(const FlashResult& result) const = 0;

    // =========================================================================
    // Metadata
    // =========================================================================

    /**
     * @brief Get algorithm name
     * @return Name of the flash algorithm (e.g., "Gibbs-SS", "Helmholtz-Newton")
     */
    virtual std::string algorithmName() const = 0;

    /**
     * @brief Get underlying EOS
     * @return Shared pointer to the EOS
     */
    virtual EOSPtr eos() const = 0;
};

/**
 * @brief Shared pointer type alias for IFlashSolver
 */
using FlashSolverPtr = std::shared_ptr<IFlashSolver>;

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_IFLASH_H
