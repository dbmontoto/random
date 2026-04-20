/**
 * @file dippr.h
 * @brief DIPPR correlation calculator
 *
 * Main calculator class for evaluating DIPPR correlations.
 * This class is EOS-agnostic and can be used with any thermodynamic model.
 *
 * Reference: DIPPR 801 Database, AIChE
 */

#ifndef DIPPR_DIPPR_H
#define DIPPR_DIPPR_H

#include "equations.h"

namespace DIPPR {

/**
 * @brief DIPPR correlation calculator
 *
 * Provides methods to calculate thermophysical properties using
 * DIPPR correlations. Designed to work alongside any equation of state
 * for hybrid property methods.
 *
 * Usage:
 * @code
 * DIPPR::Parameters params = DIPPR::Database::getInstance().getByName("n-hexane");
 * DIPPR::Calculator calc(params);
 * double Psat = calc.vaporPressure(300.0);  // Pa
 * double rho_L = calc.liquidDensity(300.0); // mol/m^3
 * @endcode
 */
class Calculator {
public:
    /**
     * @brief Constructor with parameters
     * @param params DIPPR parameters for the component
     */
    explicit Calculator(const Parameters& params);

    /**
     * @brief Default constructor (must call setParameters before use)
     */
    Calculator() = default;

    /**
     * @brief Set DIPPR parameters
     */
    void setParameters(const Parameters& params);

    /**
     * @brief Get the parameters
     */
    const Parameters& getParameters() const { return params_; }

    // ========================================================================
    // Primary Property Calculations
    // ========================================================================

    /**
     * @brief Calculate saturated liquid density
     * @param T Temperature [K]
     * @return Liquid molar density [mol/m^3]
     */
    double liquidDensity(double T) const;

    /**
     * @brief Calculate vapor pressure
     * @param T Temperature [K]
     * @return Vapor pressure [Pa]
     */
    double vaporPressure(double T) const;

    /**
     * @brief Calculate heat of vaporization
     * @param T Temperature [K]
     * @return Heat of vaporization [J/mol]
     */
    double heatOfVaporization(double T) const;

    /**
     * @brief Calculate liquid heat capacity
     * @param T Temperature [K]
     * @return Cp [J/(mol*K)]
     */
    double liquidCp(double T) const;

    /**
     * @brief Calculate liquid viscosity
     * @param T Temperature [K]
     * @return Viscosity [Pa*s]
     */
    double liquidViscosity(double T) const;

    /**
     * @brief Calculate vapor viscosity
     * @param T Temperature [K]
     * @return Viscosity [Pa*s]
     */
    double vaporViscosity(double T) const;

    /**
     * @brief Calculate liquid thermal conductivity
     * @param T Temperature [K]
     * @return Thermal conductivity [W/(m*K)]
     */
    double liquidThermalConductivity(double T) const;

    /**
     * @brief Calculate surface tension
     * @param T Temperature [K]
     * @return Surface tension [N/m]
     */
    double surfaceTension(double T) const;

    /**
     * @brief Calculate ideal gas heat capacity
     * @param T Temperature [K]
     * @return Cp [J/(mol*K)]
     */
    double idealGasCp(double T) const;

    // ========================================================================
    // Derived Properties
    // ========================================================================

    /**
     * @brief Calculate liquid enthalpy relative to reference state
     * @param T Temperature [K]
     * @param T_ref Reference temperature [K] (default: 298.15 K)
     * @return Enthalpy [J/mol]
     */
    double liquidEnthalpy(double T, double T_ref = 298.15) const;

    /**
     * @brief Calculate liquid entropy relative to reference state
     * @param T Temperature [K]
     * @param T_ref Reference temperature [K]
     * @return Entropy [J/(mol*K)]
     */
    double liquidEntropy(double T, double T_ref = 298.15) const;

    /**
     * @brief Calculate ideal gas enthalpy relative to reference state
     * @param T Temperature [K]
     * @param T_ref Reference temperature [K] (default: 298.15 K)
     * @return Enthalpy [J/mol]
     */
    double idealGasEnthalpy(double T, double T_ref = 298.15) const;

    /**
     * @brief Calculate ideal gas entropy relative to reference state
     * @param T Temperature [K]
     * @param T_ref Reference temperature [K]
     * @return Entropy [J/(mol*K)]
     */
    double idealGasEntropy(double T, double T_ref = 298.15) const;

    // ========================================================================
    // Critical Properties Access
    // ========================================================================

    double Tc() const { return params_.Tc; }
    double Pc() const { return params_.Pc; }
    double Vc() const { return params_.Vc; }
    double Zc() const { return params_.Zc; }
    double omega() const { return params_.omega; }

    /**
     * @brief Calculate reduced temperature
     */
    double reducedTemperature(double T) const { return T / params_.Tc; }

    /**
     * @brief Check if temperature is in valid range for a property
     */
    bool isValidTemperature(double T, const Coefficients& coeffs) const;

    // ========================================================================
    // Static Equation Evaluation
    // ========================================================================

    /**
     * @brief Evaluate a DIPPR correlation (static method)
     * @param T Temperature [K]
     * @param coeffs Correlation coefficients
     * @param Tc Critical temperature [K] (for EQ_106)
     * @return Property value
     */
    static double evaluate(double T, const Coefficients& coeffs, double Tc = 0.0);

private:
    Parameters params_;
    bool initialized_ = false;

    void checkInitialized() const;
};

} // namespace DIPPR

#endif // DIPPR_DIPPR_H
