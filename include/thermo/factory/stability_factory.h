/**
 * @file stability_factory.h
 * @brief Factory for creating stability analyzer objects
 */

#ifndef THERMO_FACTORY_STABILITY_FACTORY_H
#define THERMO_FACTORY_STABILITY_FACTORY_H

#include <memory>
#include <string>
#include "thermo/equilibrium/istability.h"
#include "thermo/eos.h"

namespace DMThermo {
namespace Factory {

/**
 * @brief Factory for creating stability analyzer implementations
 *
 * Example usage:
 * @code
 * auto stability = StabilityFactory::create(eos);
 * auto result = stability->analyze(T, P, z);
 * @endcode
 */
class StabilityFactory {
public:
    // =========================================================================
    // Creation Methods
    // =========================================================================

    /**
     * @brief Create default stability analyzer
     *
     * @param eos Equation of state
     * @return Stability analyzer instance
     */
    static Equilibrium::Stability::StabilityAnalyzerPtr create(EOSPtr eos);

    /**
     * @brief Create stability analyzer with TPD method only
     *
     * @param eos Equation of state
     * @return TPD-based stability analyzer
     */
    static Equilibrium::Stability::StabilityAnalyzerPtr createTPD(EOSPtr eos);

    /**
     * @brief Create stability analyzer with spinodal detection
     *
     * @param eos Equation of state
     * @return Spinodal-based stability analyzer
     */
    static Equilibrium::Stability::StabilityAnalyzerPtr createWithSpinodal(EOSPtr eos);

    /**
     * @brief Create stability analyzer optimized for LLE systems
     *
     * @param eos Equation of state
     * @return LLE-optimized stability analyzer
     */
    static Equilibrium::Stability::StabilityAnalyzerPtr createForLLE(EOSPtr eos);

private:
    static void initializeBuiltinAnalyzers();
    static bool initialized_;
};

} // namespace Factory
} // namespace DMThermo

#endif // THERMO_FACTORY_STABILITY_FACTORY_H
