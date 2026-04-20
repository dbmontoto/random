/**
 * @file constants.h
 * @brief Physical and mathematical constants for thermodynamic calculations
 *
 * This file is part of the DMThermo namespace architecture.
 * All constants use SI units unless otherwise specified.
 */

#ifndef THERMO_CORE_CONSTANTS_H
#define THERMO_CORE_CONSTANTS_H

namespace DMThermo {
namespace Constants {

/// Universal gas constant [J/(mol·K)]
constexpr double GAS_CONSTANT = 8.314462618;

/// Avogadro's number [1/mol]
constexpr double AVOGADRO = 6.02214076e23;

/// Boltzmann constant [J/K]
constexpr double BOLTZMANN = 1.380649e-23;

/// Pi constant
constexpr double PI = 3.14159265358979323846;

/// Default absolute tolerance for convergence
constexpr double DEFAULT_TOLERANCE = 1e-10;

/// Default relative tolerance for convergence
constexpr double DEFAULT_REL_TOLERANCE = 1e-8;

/// Default maximum iterations
constexpr int DEFAULT_MAX_ITERATIONS = 100;

/// Minimum density for calculations [mol/m³]
constexpr double MIN_DENSITY = 1e-10;

/// Maximum density (close-packing limit) [mol/m³]
constexpr double MAX_DENSITY = 100000.0;

/// Minimum temperature [K]
constexpr double MIN_TEMPERATURE = 1.0;

/// Maximum temperature [K]
constexpr double MAX_TEMPERATURE = 10000.0;

/// Minimum pressure [Pa]
constexpr double MIN_PRESSURE = 1.0;

/// Maximum pressure [Pa]
constexpr double MAX_PRESSURE = 1e10;

/// Minimum mole fraction considered non-zero
constexpr double MIN_MOLE_FRACTION = 1e-15;

} // namespace Constants
} // namespace DMThermo

#endif // THERMO_CORE_CONSTANTS_H
