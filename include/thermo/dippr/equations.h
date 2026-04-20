/**
 * @file equations.h
 * @brief DIPPR equation types and coefficient structures
 *
 * DIPPR (Design Institute for Physical Properties) provides standardized
 * correlations for thermophysical properties. This header defines the
 * equation types and coefficient structures used by all DIPPR correlations.
 *
 * Reference: DIPPR 801 Database, AIChE
 */

#ifndef DIPPR_EQUATIONS_H
#define DIPPR_EQUATIONS_H

#include <string>
#include <optional>
#include <string_view>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace DIPPR {

/**
 * @brief DIPPR equation types
 *
 * Standard equation forms used by DIPPR correlations:
 * - EQ_100: Polynomial (liquid Cp, some viscosities)
 * - EQ_101: Extended Antoine/Wagner (vapor pressure, some viscosities)
 * - EQ_102: Power law with rational correction (vapor viscosity, etc.)
 * - EQ_105: Rackett equation (liquid density)
 * - EQ_106: Watson/corresponding states (heat of vaporization, surface tension)
 * - EQ_107: Aly-Lee (ideal gas heat capacity)
 * - EQ_123: Reduced-temperature form for saturated liquid thermal conductivity (requires Tc)
 * - EQ_114: Alternate Cp form using reduced temperature (rare; liquid Cp)
 */
enum class EquationType {
    EQ_100,  ///< Polynomial: Y = A + B*T + C*T^2 + D*T^3 + E*T^4
    EQ_101,  ///< Extended Antoine: Y = exp(A + B/T + C*ln(T) + D*T^E)
    EQ_102,  ///< Power law: Y = A*T^B / (1 + C/T + D/T^2)
    EQ_105,  ///< Rackett: Y = A / B^(1 + (1-T/C)^D)
    EQ_106,  ///< Watson: Y = A * (1 - Tr)^(B + C*Tr + D*Tr^2 + E*Tr^3)
    EQ_107,  ///< Aly-Lee: Y = A + B*((C/T)/sinh(C/T))^2 + D*((E/T)/cosh(E/T))^2
    EQ_114,  ///< Reduced-temperature Cp form (requires Tc)
    EQ_123   ///< Reduced-temperature liquid thermal conductivity form (requires Tc)
};

/**
 * @brief Map a DIPPR equation form number (e.g., 101) to an EquationType.
 *
 * Keep equation-form mappings centralized here so CSV adapters and legacy loaders
 * don't need to duplicate switch statements.
 */
inline std::optional<EquationType> equationTypeFromForm(int eq_form) {
    switch (eq_form) {
        case 100: return EquationType::EQ_100;
        case 101: return EquationType::EQ_101;
        case 102: return EquationType::EQ_102;
        case 105: return EquationType::EQ_105;
        case 106: return EquationType::EQ_106;
        case 107: return EquationType::EQ_107;
        case 114: return EquationType::EQ_114;
        case 123: return EquationType::EQ_123;
        default: return std::nullopt;
    }
}

inline std::string_view equationTypeName(EquationType eq) {
    switch (eq) {
        case EquationType::EQ_100: return "EQ_100";
        case EquationType::EQ_101: return "EQ_101";
        case EquationType::EQ_102: return "EQ_102";
        case EquationType::EQ_105: return "EQ_105";
        case EquationType::EQ_106: return "EQ_106";
        case EquationType::EQ_107: return "EQ_107";
        case EquationType::EQ_114: return "EQ_114";
        case EquationType::EQ_123: return "EQ_123";
        default: return "Unknown";
    }
}

/**
 * @brief DIPPR correlation coefficients
 *
 * Generic structure to hold coefficients for any DIPPR equation type.
 * The interpretation of coefficients depends on the equation type.
 */
struct Coefficients {
    EquationType equation = EquationType::EQ_100;
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double D = 0.0;
    double E = 0.0;
    double F = 0.0;
    double G = 0.0;
    double T_min = 0.0;    ///< Minimum valid temperature [K]
    double T_max = 1000.0; ///< Maximum valid temperature [K]

    Coefficients() = default;

    /**
     * @brief Constructor with explicit equation type
     */
    Coefficients(EquationType eq, double a, double b, double c, double d, double e,
                 double tmin = 0.0, double tmax = 1000.0)
        : equation(eq), A(a), B(b), C(c), D(d), E(e), F(0.0), G(0.0), T_min(tmin), T_max(tmax) {}

    /**
     * @brief Constructor with extended coefficients (for equation forms that store extra parameters)
     */
    Coefficients(EquationType eq, double a, double b, double c, double d, double e, double f, double g,
                 double tmin, double tmax)
        : equation(eq), A(a), B(b), C(c), D(d), E(e), F(f), G(g), T_min(tmin), T_max(tmax) {}

    /**
     * @brief Backward-compatible constructor (defaults to polynomial EQ_100)
     *
     * For polynomial Cp: Cp = A + B*T + C*T^2 + D*T^3 + E*T^4
     */
    Coefficients(double a, double b, double c, double d, double e,
                 double tmin = 200.0, double tmax = 1500.0)
        : equation(EquationType::EQ_100), A(a), B(b), C(c), D(d), E(e), F(0.0), G(0.0), T_min(tmin), T_max(tmax) {}

    /**
     * @brief Check if coefficients are valid (non-default)
     */
    bool isValid() const {
        if (!(std::isfinite(T_min) && std::isfinite(T_max) && T_max > T_min && T_max > 0.0)) {
            return false;
        }

        auto finite = [](double v) { return std::isfinite(v); };
        auto nonzero = [&](double v) { return finite(v) && (v != 0.0); };
        const bool any_nonzero =
            nonzero(A) || nonzero(B) || nonzero(C) || nonzero(D) ||
            nonzero(E) || nonzero(F) || nonzero(G);
        if (!any_nonzero) {
            return false;
        }

        switch (equation) {
            case EquationType::EQ_105:
                // Rackett: requires A,B,C(=Tc),D and B,C positive.
                return finite(A) && finite(B) && finite(C) && finite(D)
                    && (B > 0.0) && (C > 0.0);
            case EquationType::EQ_102:
                // Power law with rational correction: requires A,B,C,D. E/F/G are accepted for schema compatibility.
                return finite(A) && finite(B) && finite(C) && finite(D);
            case EquationType::EQ_107:
                // Aly-Lee: requires A,B,C,D,E. F/G are accepted for schema compatibility.
                return finite(A) && finite(B) && finite(C) && finite(D) && finite(E) && finite(F) && finite(G);
            case EquationType::EQ_114:
                // Reduced-temperature Cp form: requires A,B,C,D; E/F/G are accepted for schema compatibility.
                return finite(A) && finite(B) && finite(C) && finite(D);
            case EquationType::EQ_123:
                // Reduced-temperature liquid thermal conductivity form: requires A,B,C,D,E (Tc is supplied at evaluation time).
                return finite(A) && finite(B) && finite(C) && finite(D) && finite(E);
            case EquationType::EQ_100:
            case EquationType::EQ_101:
            case EquationType::EQ_106:
            default:
                return finite(A) && finite(B) && finite(C) && finite(D) && finite(E) && finite(F) && finite(G);
        }
    }

    /**
     * @brief Check if temperature is in valid range
     */
    bool inRange(double T) const {
        return T >= T_min && T <= T_max;
    }

    /**
     * @brief Evaluate the DIPPR equation at given temperature
     * @param T Temperature [K]
     * @param Tc Critical temperature [K] (required for EQ_106)
     * @return Calculated property value
     */
    double evaluate(double T, double Tc = 0.0) const;

    /**
     * @brief Calculate property value (alias for evaluate, for backward compatibility)
     * @param T Temperature [K]
     * @return Calculated property value
     *
     * For polynomial equations (EQ_100), this computes:
     *   Y = A + B*T + C*T^2 + D*T^3 + E*T^4
     */
    double calculate(double T) const {
        if (!isValid()) {
            throw std::runtime_error("DIPPR: Correlation coefficients are missing/invalid");
        }
        if (!inRange(T)) {
            throw std::runtime_error("DIPPR: Temperature out of correlation validity range");
        }
        return evaluate(T);
    }

    /**
     * @brief Calculate integral of property/T from T_ref to T (for entropy)
     *
     * For polynomial Cp (EQ_100):
     *   integral(Cp/T)dT = A*ln(T/T_ref) + B*(T-T_ref) + C/2*(T^2-T_ref^2)
     *                    + D/3*(T^3-T_ref^3) + E/4*(T^4-T_ref^4)
     *
     * @param T Temperature [K]
     * @param T_ref Reference temperature [K]
     * @return Integral value [same units as property / K]
     */
    double integrateOverT(double T, double T_ref, double Tc = 0.0) const {
        if (T <= 0.0 || T_ref <= 0.0) {
            throw std::runtime_error("DIPPR: Temperature must be positive for integration");
        }
        if (!isValid()) {
            throw std::runtime_error("DIPPR: Cannot integrate with missing/invalid correlation coefficients");
        }
        if (!inRange(T) || !inRange(T_ref)) {
            throw std::runtime_error("DIPPR: Temperature out of correlation validity range for integration");
        }

        const double T_calc = T;
        const double T_ref_calc = T_ref;
        if (T_calc == T_ref_calc) {
            return 0.0;
        }

        if (equation == EquationType::EQ_100) {
            return A * std::log(T_calc / T_ref_calc)
                 + B * (T_calc - T_ref_calc)
                 + C / 2.0 * (T_calc * T_calc - T_ref_calc * T_ref_calc)
                 + D / 3.0 * (T_calc * T_calc * T_calc - T_ref_calc * T_ref_calc * T_ref_calc)
                 + E / 4.0 * (T_calc * T_calc * T_calc * T_calc - T_ref_calc * T_ref_calc * T_ref_calc * T_ref_calc);
        }

        // Numerical integration (Simpson's rule) for non-polynomial equations.
        double a = T_ref_calc;
        double b = T_calc;
        double sign = 1.0;
        if (b < a) {
            std::swap(a, b);
            sign = -1.0;
        }

        constexpr int n_intervals = 256; // Must be even
        const double h = (b - a) / n_intervals;
        auto f = [this, Tc](double temp) { return this->evaluate(temp, Tc) / temp; };

        double sum = f(a) + f(b);
        for (int i = 1; i < n_intervals; ++i) {
            const double x = a + h * i;
            sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
        }
        return sign * (h / 3.0) * sum;
    }

    /**
     * @brief Calculate integral of property from T_ref to T (for enthalpy)
     *
     * For polynomial Cp (EQ_100):
     *   integral(Cp)dT = A*(T-T_ref) + B/2*(T^2-T_ref^2) + C/3*(T^3-T_ref^3)
     *                  + D/4*(T^4-T_ref^4) + E/5*(T^5-T_ref^5)
     *
     * @param T Temperature [K]
     * @param T_ref Reference temperature [K]
     * @return Integral value [same units as property * K]
     */
    double integrate(double T, double T_ref, double Tc = 0.0) const {
        if (!isValid()) {
            throw std::runtime_error("DIPPR: Cannot integrate with missing/invalid correlation coefficients");
        }
        if (!inRange(T) || !inRange(T_ref)) {
            throw std::runtime_error("DIPPR: Temperature out of correlation validity range for integration");
        }

        const double T_calc = T;
        const double T_ref_calc = T_ref;
        if (T_calc == T_ref_calc) {
            return 0.0;
        }

        if (equation == EquationType::EQ_100) {
            double T2 = T_calc * T_calc;
            double T3 = T2 * T_calc;
            double T4 = T3 * T_calc;
            double T5 = T4 * T_calc;
            double Tr2 = T_ref_calc * T_ref_calc;
            double Tr3 = Tr2 * T_ref_calc;
            double Tr4 = Tr3 * T_ref_calc;
            double Tr5 = Tr4 * T_ref_calc;
            return A * (T_calc - T_ref_calc)
                 + B / 2.0 * (T2 - Tr2)
                 + C / 3.0 * (T3 - Tr3)
                 + D / 4.0 * (T4 - Tr4)
                 + E / 5.0 * (T5 - Tr5);
        }

        // Numerical integration (Simpson's rule) for non-polynomial equations.
        double a = T_ref_calc;
        double b = T_calc;
        double sign = 1.0;
        if (b < a) {
            std::swap(a, b);
            sign = -1.0;
        }

        constexpr int n_intervals = 256; // Must be even
        const double h = (b - a) / n_intervals;
        auto f = [this, Tc](double temp) { return this->evaluate(temp, Tc); };

        double sum = f(a) + f(b);
        for (int i = 1; i < n_intervals; ++i) {
            const double x = a + h * i;
            sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
        }
        return sign * (h / 3.0) * sum;
    }
};

/**
 * @brief Complete set of DIPPR parameters for a component
 *
 * Contains all commonly used DIPPR correlations for a single component.
 * Each correlation is optional - check the has_* flags before use.
 */
struct Parameters {
    std::string name;
    std::string cas;

    // Critical properties (required for reduced property correlations)
    double Tc = 0.0;     ///< Critical temperature [K]
    double Pc = 0.0;     ///< Critical pressure [Pa]
    double Vc = 0.0;     ///< Critical molar volume [m^3/mol]
    double Zc = 0.0;     ///< Critical compressibility factor
    double omega = 0.0;  ///< Acentric factor

    // Correlation coefficients for various properties
    Coefficients liquid_density;       ///< Saturated liquid density [kmol/m^3]
    Coefficients vapor_pressure;       ///< Vapor pressure [Pa]
    Coefficients heat_of_vaporization; ///< Heat of vaporization [J/kmol]
    Coefficients liquid_heat_capacity; ///< Liquid heat capacity [J/(kmol*K)]
    Coefficients liquid_viscosity;     ///< Liquid viscosity [Pa*s]
    Coefficients vapor_viscosity;      ///< Vapor viscosity [Pa*s]
    Coefficients liquid_thermal_conductivity; ///< Liquid thermal conductivity [W/(m*K)]
    Coefficients surface_tension;      ///< Surface tension [N/m]
    Coefficients ideal_gas_cp;         ///< Ideal gas heat capacity [J/(kmol*K)]

    bool has_liquid_density = false;
    bool has_vapor_pressure = false;
    bool has_heat_of_vaporization = false;
    bool has_liquid_heat_capacity = false;
    bool has_liquid_viscosity = false;
    bool has_vapor_viscosity = false;
    bool has_liquid_thermal_conductivity = false;
    bool has_surface_tension = false;
    bool has_ideal_gas_cp = false;
};

} // namespace DIPPR

#endif // DIPPR_EQUATIONS_H
