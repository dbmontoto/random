/**
 * @file flash_factory.h
 * @brief Factory for creating flash solver objects
 */

#ifndef THERMO_FACTORY_FLASH_FACTORY_H
#define THERMO_FACTORY_FLASH_FACTORY_H

#include <memory>
#include <string>
#include <functional>
#include <unordered_map>
#include "thermo/equilibrium/iflash.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/eos.h"

namespace DMThermo {
namespace Factory {

/**
 * @brief Factory for creating flash solver implementations
 *
 * Example usage:
 * @code
 * // Create default flash solver (auto-selects method)
 * auto flash = FlashFactory::create(eos);
 *
 * // Create specific solver type
 * auto gibbs = FlashFactory::createGibbs(eos);
 * auto helmholtz = FlashFactory::createHelmholtz(eos);
 * @endcode
 */
class FlashFactory {
public:
    /// Creator function type
    using FlashCreator = std::function<Equilibrium::Flash::FlashSolverPtr(
        EOSPtr,
        Equilibrium::Stability::StabilityAnalyzerPtr
    )>;

    // =========================================================================
    // Generic Creation
    // =========================================================================

    /**
     * @brief Create flash solver by method name
     *
     * @param method Method name (e.g., "Gibbs", "Helmholtz", "Auto")
     * @param eos Equation of state
     * @param stability Stability analyzer (optional, will create if null)
     * @return Flash solver instance
     */
    static Equilibrium::Flash::FlashSolverPtr create(
        const std::string& method,
        EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability = nullptr
    );

    /**
     * @brief Create default flash solver (auto-selects best method)
     *
     * @param eos Equation of state
     * @return Flash solver instance
     */
    static Equilibrium::Flash::FlashSolverPtr create(EOSPtr eos);

    // =========================================================================
    // Specific Solver Creation
    // =========================================================================

    /**
     * @brief Create Gibbs energy minimization solver
     *
     * Uses successive substitution with K-values.
     *
     * @param eos Equation of state
     * @param stability Stability analyzer (optional)
     * @return Gibbs flash solver
     */
    static Equilibrium::Flash::FlashSolverPtr createGibbs(
        EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability = nullptr
    );

    /**
     * @brief Create Helmholtz energy minimization solver
     *
     * Uses Newton-based direct minimization.
     *
     * @param eos Equation of state
     * @param stability Stability analyzer (optional)
     * @return Helmholtz flash solver
     */
    static Equilibrium::Flash::FlashSolverPtr createHelmholtz(
        EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability = nullptr
    );

    /**
     * @brief Create combined solver (Gibbs start, Helmholtz refinement)
     *
     * @param eos Equation of state
     * @param stability Stability analyzer (optional)
     * @return Combined flash solver
     */
    static Equilibrium::Flash::FlashSolverPtr createCombined(
        EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability = nullptr
    );

    // =========================================================================
    // Registration
    // =========================================================================

    /**
     * @brief Register custom flash solver creator
     *
     * @param name Solver type name
     * @param creator Function to create solver
     */
    static void registerFlashSolver(const std::string& name, FlashCreator creator);

    /**
     * @brief Get list of registered solver types
     */
    static std::vector<std::string> registeredTypes();

private:
    static std::unordered_map<std::string, FlashCreator>& registry();
    static void initializeBuiltinSolvers();
    static bool initialized_;
};

} // namespace Factory
} // namespace DMThermo

#endif // THERMO_FACTORY_FLASH_FACTORY_H
