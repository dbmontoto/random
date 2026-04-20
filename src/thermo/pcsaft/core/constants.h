/**
 * @file constants.h
 * @brief Central repository for physical and mathematical constants
 *
 * This header file centralizes all constants used throughout the PC-SAFT
 * implementation to avoid duplication and ensure consistency.
 */

#ifndef PCSAFT_CONSTANTS_H
#define PCSAFT_CONSTANTS_H

#include <cmath>

namespace pcsaft {

// ============================================================================
// Mathematical Constants
// ============================================================================

/// Pi constant with high precision
constexpr double PI = 3.14159265358979323846;

/// Euler's number
constexpr double E = 2.71828182845904523536;

/// Square root of 2
constexpr double SQRT2 = 1.41421356237309504880;

// ============================================================================
// Physical Constants
// ============================================================================

/// Universal gas constant [J/(mol·K)]
constexpr double R = 8.314462618;

/// Avogadro's number [1/mol]
constexpr double N_A = 6.02214076e23;

/// Boltzmann constant [J/K]
constexpr double k_B = R / N_A;

/// Standard pressure [Pa]
constexpr double P_STD = 101325.0;

/// Standard temperature [K]
constexpr double T_STD = 298.15;

/// Zero Celsius in Kelvin
constexpr double T_ZERO_C = 273.15;

// ============================================================================
// Numerical Constants
// ============================================================================

/// Machine epsilon for double precision
constexpr double EPSILON = 1e-15;

/// Default numerical derivative step size (relative)
constexpr double DERIV_STEP = 1e-6;

/// Default convergence tolerance
constexpr double DEFAULT_TOL = 1e-8;

/// Maximum iterations for iterative solvers
constexpr int MAX_ITERATIONS = 100;

/// Value representing numerical divergence
constexpr double DIVERGENCE_VALUE = 1e10;

// ============================================================================
// PC-SAFT Specific Constants
// ============================================================================

/// Maximum packing fraction before numerical issues
constexpr double ETA_MAX = 0.74;

/// Packing fraction where liquid behavior becomes problematic
constexpr double ETA_LIQUID_THRESHOLD = 0.35;

/// Default temperature range bounds [K]
constexpr double T_MIN = 100.0;
constexpr double T_MAX = 1000.0;

/// Default pressure range bounds [Pa]
constexpr double P_MIN = 1e3;   // 0.01 bar
constexpr double P_MAX = 1e8;   // 1000 bar

/// Default density range bounds [mol/m³]
constexpr double RHO_MIN = 1e-3;
constexpr double RHO_MAX = 1e6;

/// Phase classification threshold [mol/m³]
/// (Note: This should be replaced with proper phase classification)
constexpr double RHO_PHASE_THRESHOLD = 10000.0;

// ============================================================================
// Unit Conversion Factors
// ============================================================================

/// Bar to Pascal conversion
constexpr double BAR_TO_PA = 1e5;

/// Pascal to bar conversion
constexpr double PA_TO_BAR = 1e-5;

/// Atmosphere to Pascal conversion
constexpr double ATM_TO_PA = 101325.0;

/// mmHg to Pascal conversion
constexpr double MMHG_TO_PA = 133.322;

/// Angstrom to meter conversion
constexpr double ANGSTROM_TO_M = 1e-10;

/// kcal/mol to J/mol conversion
constexpr double KCAL_TO_J = 4184.0;

// ============================================================================
// Utility Functions for Constants
// ============================================================================

/**
 * @brief Check if a value is effectively zero
 * @param value The value to check
 * @return true if |value| < EPSILON
 */
inline bool isZero(double value) {
    return std::abs(value) < EPSILON;
}

/**
 * @brief Check if a value represents divergence
 * @param value The value to check
 * @return true if |value| > DIVERGENCE_VALUE
 */
inline bool isDiverged(double value) {
    return std::abs(value) > DIVERGENCE_VALUE;
}

/**
 * @brief Smooth lower bound clamping (autodiff-safe)
 * @param value The value to clamp
 * @param lower_bound Minimum allowed value
 * @param epsilon Smoothing parameter
 * @return Smoothly clamped value
 *
 * Uses: value + epsilon + sqrt((value - lower_bound)^2 + epsilon^2)
 * This is differentiable everywhere and smoothly enforces value >= lower_bound
 */
inline double smoothClampLower(double value, double lower_bound, double epsilon = 1e-10) {
    double delta = value - lower_bound;
    return lower_bound + 0.5 * (delta + std::sqrt(delta * delta + epsilon * epsilon));
}

/**
 * @brief Smooth upper bound clamping (autodiff-safe)
 * @param value The value to clamp
 * @param upper_bound Maximum allowed value
 * @param epsilon Smoothing parameter
 * @return Smoothly clamped value
 */
inline double smoothClampUpper(double value, double upper_bound, double epsilon = 1e-10) {
    double delta = value - upper_bound;
    return upper_bound + 0.5 * (delta - std::sqrt(delta * delta + epsilon * epsilon));
}

/**
 * @brief Smooth clamping to range [lower, upper] (autodiff-safe)
 * @param value The value to clamp
 * @param lower_bound Minimum allowed value
 * @param upper_bound Maximum allowed value
 * @param epsilon Smoothing parameter
 * @return Smoothly clamped value
 */
inline double smoothClamp(double value, double lower_bound, double upper_bound, double epsilon = 1e-10) {
    return smoothClampUpper(smoothClampLower(value, lower_bound, epsilon), upper_bound, epsilon);
}

// ============================================================================
// Templated Smooth Clamping (for autodiff compatibility)
// ============================================================================

/**
 * @brief Templated smooth lower bound clamping (autodiff-safe)
 * Works with autodiff::dual and autodiff::dual2nd types
 */
template<typename T>
T smoothClampLowerT(T value, double lower_bound, double epsilon = 1e-10) {
    T delta = value - lower_bound;
    using std::sqrt;
    return lower_bound + 0.5 * (delta + sqrt(delta * delta + epsilon * epsilon));
}

/**
 * @brief Templated smooth upper bound clamping (autodiff-safe)
 */
template<typename T>
T smoothClampUpperT(T value, double upper_bound, double epsilon = 1e-10) {
    T delta = value - upper_bound;
    using std::sqrt;
    return upper_bound + 0.5 * (delta - sqrt(delta * delta + epsilon * epsilon));
}

/**
 * @brief Templated smooth clamping to range [lower, upper] (autodiff-safe)
 */
template<typename T>
T smoothClampT(T value, double lower_bound, double upper_bound, double epsilon = 1e-10) {
    return smoothClampUpperT(smoothClampLowerT(value, lower_bound, epsilon), upper_bound, epsilon);
}

/**
 * @brief Safe division with small epsilon to avoid division by zero (autodiff-safe)
 * @param numerator Numerator
 * @param denominator Denominator
 * @param epsilon Small value to add to denominator
 * @return numerator / (denominator + sign(denominator) * epsilon)
 */
inline double safeDivide(double numerator, double denominator, double epsilon = 1e-100) {
    double sign = (denominator >= 0.0) ? 1.0 : -1.0;
    return numerator / (denominator + sign * epsilon);
}

} // namespace pcsaft

#endif // PCSAFT_CONSTANTS_H
