/**
 * @file thermo.h
 * @brief Main header for the Thermo library
 *
 * This header includes all public interfaces for the thermodynamic library.
 * For most use cases, including this single header is sufficient.
 *
 * Example usage:
 * @code
 * #include <thermo/thermo.h>
 *
 * using namespace DMThermo;
 *
 * // Create components and mixture
 * Core::Component methane("methane", 1.0, 3.7039, 150.03);
 * Core::Component ethane("ethane", 1.6069, 3.5206, 191.42);
 * Core::Mixture mixture({methane, ethane});
 *
 * // Create EOS and flash solver
 * auto eos = Factory::EOSFactory::createPCSaft(mixture);
 * auto flash = Factory::FlashFactory::create(eos);
 *
 * // Perform flash calculation
 * auto result = flash->solvePT(300.0, 1e6, {0.5, 0.5});
 * if (result.converged && result.is_two_phase) {
 *     std::cout << "Vapor fraction: " << result.vapor_fraction << std::endl;
 * }
 * @endcode
 */

#ifndef THERMO_THERMO_H
#define THERMO_THERMO_H

// Core types and constants
#include "core/constants.h"
#include "core/types.h"
#include "core/component.h"
#include "core/mixture.h"
#include "core/state.h"

// Configuration
#include "config/flash_config.h"
#include "config/stability_config.h"
#include "config/solver_config.h"
#include "config/error_handling.h"

// EOS interface and implementations
#include "eos.h"                          // DMThermo::EOS abstract interface
#include "eos/PCSAFT.h"                   // DMThermo::PCSaft::PCSaftEOS
#include "eos/PengRobinson.h"             // DMThermo::Cubic::PengRobinsonEOS
#include "eos/VTPR.h"                     // DMThermo::Cubic::VTPREOS
#include "eos/PRVT.h"                     // DMThermo::Cubic::PRVTEOS
#include "eos/SRK.h"                      // DMThermo::Cubic::SRKEOS

// Equilibrium interfaces
#include "equilibrium/iflash.h"
#include "equilibrium/istability.h"
#include "equilibrium/flash_result.h"
#include "equilibrium/stability_result.h"
#include "equilibrium/tpd_stability.h"     // EOS-based stability analyzer
#include "equilibrium/gibbs_multi_flash.h" // EOS-agnostic multi-phase PT flash solver
#include "equilibrium/pure_saturation.h"   // EOS-based pure saturation solver

// Numerical methods interfaces
#include "numerics/iroot_finder.h"
#include "numerics/ioptimizer.h"
#include "numerics/ifixed_point.h"

// Properties + high-level engine
#include "properties/thermo_properties.h"
#include "engine/property_engine.h"

// Factories
#include "factory/eos_factory.h"
#include "factory/flash_factory.h"
#include "factory/stability_factory.h"
#include "factory/numerics_factory.h"
#include "factory/mixture_factory.h"
#include "factory/reaction_factory.h"

// Databanks (CSV-backed)
#include "data/databanks.h"

// Reactions (work-in-progress)
#include "reactions/reaction_system.h"
#include "reactions/gibbs_reaction_equilibrium.h"
#include "reactions/thermochemistry.h"
#include "reactions/ireaction_equilibrium.h"
#include "reactions/auto_reaction_equilibrium.h"
#include "reactions/reactive_flash.h"

/**
 * @namespace DMThermo
 * @brief Root namespace for the thermodynamic library
 *
 * The DMThermo namespace contains all thermodynamic modeling functionality
 * organized into sub-namespaces:
 *
 * - **Core**: Fundamental types (Component, Mixture, State)
 * - **EOS**: Equation of state interfaces and implementations
 * - **Equilibrium**: Phase equilibrium calculations (Flash, Stability)
 * - **Properties**: Thermodynamic property calculations
 * - **Numerics**: Numerical methods (root finding, optimization)
 * - **Config**: Configuration structures
 * - **Factory**: Object creation factories
 * - **Constants**: Physical and numerical constants
 */
namespace DMThermo {

/**
 * @brief Library version information
 */
struct Version {
    static constexpr int major = 1;
    static constexpr int minor = 0;
    static constexpr int patch = 0;
    static constexpr const char* string = "1.0.0";
    static constexpr const char* name = "Enterprise Edition";
};

} // namespace DMThermo

#endif // THERMO_THERMO_H
