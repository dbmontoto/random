/**
 * @file flash_result.h
 * @brief Result structures for flash calculations
 */

#ifndef THERMO_EQUILIBRIUM_FLASH_RESULT_H
#define THERMO_EQUILIBRIUM_FLASH_RESULT_H

#include <vector>
#include <string>
#include <optional>
#include "thermo/core/types.h"

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

/**
 * @brief Detailed result from flash calculations
 */
struct FlashResult {
    // =========================================================================
    // Convergence Information
    // =========================================================================

    /// Did the calculation converge?
    bool converged = false;

    /// Number of iterations to converge
    int iterations = 0;

    /// Final residual norm
    double residual = 0.0;

    /// Convergence message or error description
    std::string message;

    // =========================================================================
    // State Variables
    // =========================================================================

    /// Temperature [K]
    double temperature = 0.0;

    /// Pressure [Pa]
    double pressure = 0.0;

    /// Feed composition (input)
    std::vector<double> z;

    // =========================================================================
    // Phase Information
    // =========================================================================

    /// Number of equilibrium phases (1, 2, or 3)
    int num_phases = 1;

    /// All phase states (ordered: vapor, liquid1, liquid2 if present)
    std::vector<PhaseState> phases;

    /// Vapor phase fraction (beta) [-], 0 if no vapor
    double vapor_fraction = 0.0;

    /// Is the system in two-phase region?
    bool is_two_phase = false;

    /// Is the system in three-phase (VLLE) region?
    bool is_three_phase = false;

    // =========================================================================
    // Convenience Accessors
    // =========================================================================

    /**
     * @brief Get vapor phase state (if present)
     * @return Pointer to vapor phase, nullptr if not present
     */
    const PhaseState* vapor() const {
        for (const auto& phase : phases) {
            if (phase.type == PhaseType::Vapor) {
                return &phase;
            }
        }
        return nullptr;
    }

    /**
     * @brief Get liquid phase state (if present)
     * @return Pointer to liquid phase, nullptr if not present
     */
    const PhaseState* liquid() const {
        for (const auto& phase : phases) {
            if (phase.type == PhaseType::Liquid ||
                phase.type == PhaseType::Liquid1) {
                return &phase;
            }
        }
        return nullptr;
    }

    /**
     * @brief Get second liquid phase (in VLLE)
     * @return Pointer to second liquid, nullptr if not present
     */
    const PhaseState* liquid2() const {
        for (const auto& phase : phases) {
            if (phase.type == PhaseType::Liquid2) {
                return &phase;
            }
        }
        return nullptr;
    }

    /**
     * @brief Get K-values (vapor/liquid equilibrium ratios)
     * @return K_i = y_i / x_i for each component
     */
    std::vector<double> kValues() const {
        const PhaseState* vap = vapor();
        const PhaseState* liq = liquid();
        if (!vap || !liq) return {};

        std::vector<double> K(vap->x.size());
        for (size_t i = 0; i < K.size(); ++i) {
            K[i] = (liq->x[i] > 1e-15) ? vap->x[i] / liq->x[i] : 1e6;
        }
        return K;
    }

    // =========================================================================
    // Optional Thermodynamic Properties
    // =========================================================================

    /// Mixture enthalpy [J/mol]
    std::optional<double> enthalpy;

    /// Mixture entropy [J/(mol·K)]
    std::optional<double> entropy;

    /// Mixture Gibbs energy [J/mol]
    std::optional<double> gibbs_energy;

    /// Heat capacity at constant pressure [J/(mol·K)]
    std::optional<double> Cp;

    /// Heat capacity at constant volume [J/(mol·K)]
    std::optional<double> Cv;

    // =========================================================================
    // Debug Information
    // =========================================================================

    /// Iteration history (if verbose enabled)
    std::vector<double> residual_history;

    /// Method used for convergence
    std::string method_used;

    /// Stability test result (if performed)
    bool stability_test_performed = false;
    bool feed_was_stable = true;
};

/**
 * @brief Simplified result for batch calculations
 */
struct FlashResultSimple {
    bool converged = false;
    double temperature = 0.0;
    double pressure = 0.0;
    double vapor_fraction = 0.0;
    std::vector<double> x_liquid;
    std::vector<double> y_vapor;
};

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_FLASH_RESULT_H
