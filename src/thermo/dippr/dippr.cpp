/**
 * @file dippr.cpp
 * @brief Implementation of DIPPR Calculator and Coefficients
 */

#include "thermo/dippr/dippr.h"
#include "thermo/dippr/equations.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace DIPPR {

namespace {

std::string formatRange(double T, const Coefficients& c) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << "(T=" << std::setprecision(6) << T << " K, valid [" << c.T_min << ", " << c.T_max << "] K)";
    return oss.str();
}

} // namespace

// ============================================================================
// Coefficients Implementation
// ============================================================================

double Coefficients::evaluate(double T, double Tc) const {
    double result = 0.0;

    switch (equation) {
        case EquationType::EQ_100:
            // Polynomial: Y = A + B*T + C*T^2 + D*T^3 + E*T^4
            result = A + B * T + C * T * T
                   + D * T * T * T + E * T * T * T * T;
            break;

        case EquationType::EQ_101:
            // Extended Antoine: Y = exp(A + B/T + C*ln(T) + D*T^E)
            if (T <= 0.0) {
                throw std::runtime_error("DIPPR: EQ_101 requires T > 0");
            }
            result = std::exp(A + B / T + C * std::log(T)
                            + D * std::pow(T, E));
            break;

        case EquationType::EQ_102:
            // Power law with rational correction:
            // Y = A*T^B / (1 + C/T + D/T^2)
            if (T <= 0.0) {
                throw std::runtime_error("DIPPR: EQ_102 requires T > 0");
            }
            {
                const double invT = 1.0 / T;
                const double denom = 1.0 + C * invT + D * invT * invT;
                if (denom == 0.0) {
                    throw std::runtime_error("DIPPR: EQ_102 invalid denominator");
                }
                result = (A * std::pow(T, B)) / denom;
            }
            break;

        case EquationType::EQ_105:
            // Rackett: Y = A / B^(1 + (1-T/C)^D)
            // Note: C stores Tc for this equation
            {
                const double Tc_use = C;
                if (Tc_use <= 0.0) {
                    throw std::runtime_error("DIPPR: EQ_105 requires a valid critical temperature in coefficient C");
                }
                const double Tr = T / Tc_use;
                if (Tr >= 1.0) {
                    throw std::runtime_error("DIPPR: EQ_105 is not valid at/above the critical temperature");
                }
                const double exponent = 1.0 + std::pow(1.0 - Tr, D);
                result = A / std::pow(B, exponent);
            }
            break;

        case EquationType::EQ_106:
            // Watson: Y = A * (1 - Tr)^(B + C*Tr + D*Tr^2 + E*Tr^3)
            // Tc must be provided.
            {
                if (Tc <= 0.0) {
                    throw std::runtime_error("DIPPR: EQ_106 requires a valid critical temperature");
                }
                const double Tc_use = Tc;

                const double Tr = T / Tc_use;
                if (Tr >= 1.0) {
                    throw std::runtime_error("DIPPR: EQ_106 is not valid at/above the critical temperature");
                }

                const double tau = 1.0 - Tr;

                const double exponent = (B + C * Tr + D * Tr * Tr + E * Tr * Tr * Tr);

                result = A * std::pow(tau, exponent);
            }
            break;

        case EquationType::EQ_107:
            // Aly-Lee: Y = A + B*((C/T)/sinh(C/T))^2 + D*((E/T)/cosh(E/T))^2
            // Note: Coefficients F and G are accepted for CSV schema compatibility but are unused here.
            {
                const double CT = C / T;
                const double ET = E / T;
                const double sinh_CT = std::sinh(CT);
                const double cosh_ET = std::cosh(ET);

                // Handle the removable singularity at C==0 (CT/sinh(CT) -> 1 as CT -> 0).
                const double ratio1 = (std::abs(CT) < 1e-12) ? 1.0 : (CT / sinh_CT);
                const double ratio2 = (ET == 0.0) ? 0.0 : (ET / cosh_ET);
                result = A
                       + B * (ratio1 * ratio1)
                       + D * (ratio2 * ratio2);
            }
            break;

        case EquationType::EQ_114:
            // Reduced-temperature Cp form (rare; alternate liquid heat capacity):
            // Y = A^2/tau + B - 2*A*C*tau - A*D*tau^2 - (1/3)*C^2*tau^3 - (1/2)*C*D*tau^4 - (1/5)*D^2*tau^5
            // tau = 1 - T/Tc
            {
                if (Tc <= 0.0) {
                    throw std::runtime_error("DIPPR: EQ_114 requires a valid critical temperature");
                }
                const double Tr = T / Tc;
                const double tau = 1.0 - Tr;
                if (!(tau > 0.0)) {
                    throw std::runtime_error("DIPPR: EQ_114 is not valid at/above the critical temperature");
                }
                const double A2_over_tau = (A * A) / tau;
                const double tau2 = tau * tau;
                const double tau3 = tau2 * tau;
                const double tau4 = tau3 * tau;
                const double tau5 = tau4 * tau;
                result =
                    A2_over_tau +
                    B -
                    2.0 * A * C * tau -
                    A * D * tau2 -
                    (1.0 / 3.0) * (C * C) * tau3 -
                    0.5 * C * D * tau4 -
                    0.2 * (D * D) * tau5;
            }
            break;

        case EquationType::EQ_123:
            // Reduced-temperature liquid thermal conductivity form (requires Tc):
            // Y = A * [1 + B*tau^0.35 + C*tau^(2/3) + D*tau + E*tau^(4/3)]
            // tau = 1 - T/Tc
            {
                if (Tc <= 0.0) {
                    throw std::runtime_error("DIPPR: EQ_123 requires a valid critical temperature");
                }
                const double Tr = T / Tc;
                const double tau = 1.0 - Tr;
                if (!(tau > 0.0)) {
                    throw std::runtime_error("DIPPR: EQ_123 is not valid at/above the critical temperature");
                }
                const double tau_035 = std::pow(tau, 0.35);
                const double tau_23 = std::pow(tau, 2.0 / 3.0);
                const double tau_43 = std::pow(tau, 4.0 / 3.0);
                result = A * (1.0 + B * tau_035 + C * tau_23 + D * tau + E * tau_43);
            }
            break;

        default:
            throw std::runtime_error("DIPPR: Unknown equation type");
    }

    return result;
}

// ============================================================================
// Calculator Implementation
// ============================================================================

Calculator::Calculator(const Parameters& params)
    : params_(params), initialized_(true) {}

void Calculator::setParameters(const Parameters& params) {
    params_ = params;
    initialized_ = true;
}

void Calculator::checkInitialized() const {
    if (!initialized_) {
        throw std::runtime_error("DIPPR::Calculator: Parameters not initialized");
    }
}

bool Calculator::isValidTemperature(double T, const Coefficients& coeffs) const {
    return coeffs.inRange(T);
}

double Calculator::evaluate(double T, const Coefficients& coeffs, double Tc) {
    return coeffs.evaluate(T, Tc);
}

// ============================================================================
// Primary Property Calculations
// ============================================================================

double Calculator::liquidDensity(double T) const {
    checkInitialized();
    if (!params_.has_liquid_density) {
        throw std::runtime_error("DIPPR: Liquid density correlation not available in databank");
    }
    if (!params_.liquid_density.isValid()) {
        throw std::runtime_error("DIPPR: Liquid density correlation is missing required parameters in databank");
    }
    if (!params_.liquid_density.inRange(T)) {
        throw std::runtime_error("DIPPR: Temperature out of range for liquid density correlation " + formatRange(T, params_.liquid_density));
    }

    // DIPPR 105 gives density in kmol/m^3, convert to mol/m^3
    double rho_kmol = params_.liquid_density.evaluate(T, params_.Tc);
    return rho_kmol * 1000.0;  // mol/m^3
}

double Calculator::vaporPressure(double T) const {
    checkInitialized();
    if (!params_.has_vapor_pressure) {
        throw std::runtime_error("DIPPR: Vapor pressure correlation not available in databank");
    }
    if (!params_.vapor_pressure.isValid()) {
        throw std::runtime_error("DIPPR: Vapor pressure correlation is missing required parameters in databank");
    }
    if (!params_.vapor_pressure.inRange(T)) {
        throw std::runtime_error("DIPPR: Temperature out of range for vapor pressure correlation " + formatRange(T, params_.vapor_pressure));
    }

    // DIPPR 101 gives pressure in Pa
    return params_.vapor_pressure.evaluate(T, params_.Tc);
}

double Calculator::heatOfVaporization(double T) const {
    checkInitialized();
    if (!params_.has_heat_of_vaporization) {
        throw std::runtime_error("DIPPR: Heat of vaporization correlation not available in databank");
    }
    if (!params_.heat_of_vaporization.isValid()) {
        throw std::runtime_error("DIPPR: Heat of vaporization correlation is missing required parameters in databank");
    }
    if (!params_.heat_of_vaporization.inRange(T)) {
        throw std::runtime_error("DIPPR: Temperature out of range for heat of vaporization correlation " + formatRange(T, params_.heat_of_vaporization));
    }

    // DIPPR 106 gives Hvap in J/kmol, convert to J/mol
    double Hvap_kmol = params_.heat_of_vaporization.evaluate(T, params_.Tc);
    return Hvap_kmol / 1000.0;  // J/mol
}

double Calculator::liquidCp(double T) const {
    checkInitialized();
    if (!params_.has_liquid_heat_capacity) {
        throw std::runtime_error("DIPPR: Liquid Cp correlation not available in databank");
    }
    if (!params_.liquid_heat_capacity.isValid()) {
        throw std::runtime_error("DIPPR: Liquid Cp correlation is missing required parameters in databank");
    }
    if (!params_.liquid_heat_capacity.inRange(T)) {
        throw std::runtime_error("DIPPR: Temperature out of range for liquid Cp correlation " + formatRange(T, params_.liquid_heat_capacity));
    }

    // DIPPR 100 gives Cp in J/(kmol*K), convert to J/(mol*K)
    double Cp_kmol = params_.liquid_heat_capacity.evaluate(T, params_.Tc);
    return Cp_kmol / 1000.0;  // J/(mol*K)
}

double Calculator::liquidViscosity(double T) const {
    checkInitialized();
    if (!params_.has_liquid_viscosity) {
        throw std::runtime_error("DIPPR: Liquid viscosity correlation not available in databank");
    }
    if (!params_.liquid_viscosity.isValid()) {
        throw std::runtime_error("DIPPR: Liquid viscosity correlation is missing required parameters in databank");
    }
    if (!params_.liquid_viscosity.inRange(T)) {
        throw std::runtime_error("DIPPR: Temperature out of range for liquid viscosity correlation " + formatRange(T, params_.liquid_viscosity));
    }

    // DIPPR 101 gives viscosity in Pa*s
    return params_.liquid_viscosity.evaluate(T, params_.Tc);
}

double Calculator::vaporViscosity(double T) const {
    checkInitialized();
    if (!params_.has_vapor_viscosity) {
        throw std::runtime_error("DIPPR: Vapor viscosity correlation not available in databank");
    }
    if (!params_.vapor_viscosity.isValid()) {
        throw std::runtime_error("DIPPR: Vapor viscosity correlation is missing required parameters in databank");
    }
    if (!params_.vapor_viscosity.inRange(T)) {
        throw std::runtime_error("DIPPR: Temperature out of range for vapor viscosity correlation " + formatRange(T, params_.vapor_viscosity));
    }

    // DIPPR 102 gives viscosity in Pa*s
    return params_.vapor_viscosity.evaluate(T, params_.Tc);
}

double Calculator::liquidThermalConductivity(double T) const {
    checkInitialized();
    if (!params_.has_liquid_thermal_conductivity) {
        throw std::runtime_error("DIPPR: Liquid thermal conductivity correlation not available in databank");
    }
    if (!params_.liquid_thermal_conductivity.isValid()) {
        throw std::runtime_error("DIPPR: Liquid thermal conductivity correlation is missing required parameters in databank");
    }
    if (!params_.liquid_thermal_conductivity.inRange(T)) {
        throw std::runtime_error("DIPPR: Temperature out of range for liquid thermal conductivity correlation " + formatRange(T, params_.liquid_thermal_conductivity));
    }

    // Units in csv schema are W/(m*K); assume coefficients are already for that unit system.
    return params_.liquid_thermal_conductivity.evaluate(T, params_.Tc);
}

double Calculator::surfaceTension(double T) const {
    checkInitialized();
    if (!params_.has_surface_tension) {
        throw std::runtime_error("DIPPR: Surface tension correlation not available in databank");
    }
    if (!params_.surface_tension.isValid()) {
        throw std::runtime_error("DIPPR: Surface tension correlation is missing required parameters in databank");
    }
    if (!params_.surface_tension.inRange(T)) {
        throw std::runtime_error("DIPPR: Temperature out of range for surface tension correlation " + formatRange(T, params_.surface_tension));
    }

    // DIPPR 106 gives surface tension in N/m
    return params_.surface_tension.evaluate(T, params_.Tc);
}

double Calculator::idealGasCp(double T) const {
    checkInitialized();
    if (!params_.has_ideal_gas_cp) {
        throw std::runtime_error("DIPPR: Ideal gas Cp correlation not available in databank");
    }
    if (!params_.ideal_gas_cp.isValid()) {
        throw std::runtime_error("DIPPR: Ideal gas Cp correlation is missing required parameters in databank");
    }
    if (!params_.ideal_gas_cp.inRange(T)) {
        throw std::runtime_error("DIPPR: Temperature out of range for ideal gas Cp correlation " + formatRange(T, params_.ideal_gas_cp));
    }

    // DIPPR 107 gives Cp in J/(kmol*K), convert to J/(mol*K)
    double Cp_kmol = params_.ideal_gas_cp.evaluate(T, params_.Tc);
    return Cp_kmol / 1000.0;  // J/(mol*K)
}

// ============================================================================
// Derived Properties
// ============================================================================

double Calculator::liquidEnthalpy(double T, double T_ref) const {
    checkInitialized();
    if (!params_.has_liquid_heat_capacity) {
        throw std::runtime_error("DIPPR: Cannot calculate enthalpy without Cp");
    }
    if (!params_.liquid_heat_capacity.isValid()) {
        throw std::runtime_error("DIPPR: Cannot calculate enthalpy with missing/invalid Cp correlation");
    }

    // The default reference temperature (298.15 K) may fall outside a correlation's validity
    // window (e.g., cryogenic liquids). In that case, pick the nearest in-range reference so
    // we can still compute *relative* enthalpy consistently within the valid range.
    if (!params_.liquid_heat_capacity.inRange(T_ref)) {
        T_ref = std::clamp(T_ref, params_.liquid_heat_capacity.T_min, params_.liquid_heat_capacity.T_max);
    }

    // Coefficients are in J/(kmol*K) -> integral is J/kmol -> convert to J/mol.
    return params_.liquid_heat_capacity.integrate(T, T_ref, params_.Tc) / 1000.0;
}

double Calculator::liquidEntropy(double T, double T_ref) const {
    checkInitialized();
    if (!params_.has_liquid_heat_capacity) {
        throw std::runtime_error("DIPPR: Cannot calculate entropy without Cp");
    }
    if (!params_.liquid_heat_capacity.isValid()) {
        throw std::runtime_error("DIPPR: Cannot calculate entropy with missing/invalid Cp correlation");
    }

    // See note in liquidEnthalpy regarding reference temperatures.
    if (!params_.liquid_heat_capacity.inRange(T_ref)) {
        T_ref = std::clamp(T_ref, params_.liquid_heat_capacity.T_min, params_.liquid_heat_capacity.T_max);
    }

    // Coefficients are in J/(kmol*K) -> integral is J/(kmol*K) -> convert to J/(mol*K).
    return params_.liquid_heat_capacity.integrateOverT(T, T_ref, params_.Tc) / 1000.0;
}

double Calculator::idealGasEnthalpy(double T, double T_ref) const {
    checkInitialized();
    if (!params_.has_ideal_gas_cp) {
        throw std::runtime_error("DIPPR: Cannot calculate enthalpy without Cp");
    }
    if (!params_.ideal_gas_cp.isValid()) {
        throw std::runtime_error("DIPPR: Cannot calculate enthalpy with missing/invalid ideal-gas Cp correlation");
    }

    // Coefficients are in J/(kmol*K) -> integral is J/kmol -> convert to J/mol.
    return params_.ideal_gas_cp.integrate(T, T_ref, params_.Tc) / 1000.0;
}

double Calculator::idealGasEntropy(double T, double T_ref) const {
    checkInitialized();
    if (!params_.has_ideal_gas_cp) {
        throw std::runtime_error("DIPPR: Cannot calculate entropy without Cp");
    }
    if (!params_.ideal_gas_cp.isValid()) {
        throw std::runtime_error("DIPPR: Cannot calculate entropy with missing/invalid ideal-gas Cp correlation");
    }

    // Coefficients are in J/(kmol*K) -> integral is J/(kmol*K) -> convert to J/(mol*K).
    return params_.ideal_gas_cp.integrateOverT(T, T_ref, params_.Tc) / 1000.0;
}

} // namespace DIPPR
