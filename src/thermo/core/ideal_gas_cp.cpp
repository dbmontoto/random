/**
 * @file ideal_gas_cp.cpp
 * @brief Implementation of Core ideal-gas Cp correlations
 */

#include "thermo/core/ideal_gas_cp.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace DMThermo {
namespace Core {

std::optional<CorrelationForm> correlationFormFromInt(int form) {
    switch (form) {
        case 100: return CorrelationForm::EQ_100;
        case 101: return CorrelationForm::EQ_101;
        case 102: return CorrelationForm::EQ_102;
        case 105: return CorrelationForm::EQ_105;
        case 106: return CorrelationForm::EQ_106;
        case 107: return CorrelationForm::EQ_107;
        case 114: return CorrelationForm::EQ_114;
        case 123: return CorrelationForm::EQ_123;
        default: return std::nullopt;
    }
}

namespace {

inline bool finite(double v) {
    return std::isfinite(v);
}

} // namespace

bool IdealGasCpCorrelation::isValid() const {
    // Form==0 is the sentinel "not set".
    const int form_i = static_cast<int>(form);
    if (!correlationFormFromInt(form_i).has_value()) return false;

    if (!(finite(A) && finite(B) && finite(C) && finite(D) && finite(E) && finite(F) && finite(G))) return false;
    if (!(finite(T_min) && finite(T_max))) return false;
    if (T_max < T_min) return false;

    // Basic per-form requirements (mirrors legacy DIPPR::Coefficients::isValid()).
    switch (form) {
        case CorrelationForm::EQ_102:
        case CorrelationForm::EQ_105:
            // EQ_102 and EQ_105 typically use A-D (E/F/G are accepted for schema compatibility).
            return finite(A) && finite(B) && finite(C) && finite(D);
        case CorrelationForm::EQ_107:
            // Aly-Lee uses A,B,C,D,E.
            return finite(A) && finite(B) && finite(C) && finite(D) && finite(E);
        case CorrelationForm::EQ_114:
            // Reduced-temperature Cp form: requires A,B,C,D.
            return finite(A) && finite(B) && finite(C) && finite(D);
        case CorrelationForm::EQ_123:
            // Reduced-temperature form: requires A,B,C,D,E.
            return finite(A) && finite(B) && finite(C) && finite(D) && finite(E);
        case CorrelationForm::EQ_100:
        case CorrelationForm::EQ_101:
        case CorrelationForm::EQ_106:
        default:
            return finite(A) && finite(B) && finite(C) && finite(D) && finite(E) && finite(F) && finite(G);
    }
}

bool IdealGasCpCorrelation::inRange(double T) const {
    return T >= T_min && T <= T_max;
}

double IdealGasCpCorrelation::evaluate(double T, double Tc) const {
    switch (form) {
        case CorrelationForm::EQ_100:
            return A + B * T + C * T * T + D * T * T * T + E * T * T * T * T;

        case CorrelationForm::EQ_101:
            if (T <= 0.0) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_101 requires T > 0");
            }
            return std::exp(A + B / T + C * std::log(T) + D * std::pow(T, E));

        case CorrelationForm::EQ_102: {
            if (T <= 0.0) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_102 requires T > 0");
            }
            const double invT = 1.0 / T;
            const double denom = 1.0 + C * invT + D * invT * invT;
            if (denom == 0.0) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_102 invalid denominator");
            }
            return (A * std::pow(T, B)) / denom;
        }

        case CorrelationForm::EQ_105: {
            // Rackett: C stores Tc for this equation form.
            const double Tc_use = C;
            if (Tc_use <= 0.0) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_105 requires Tc in coefficient C");
            }
            const double Tr = T / Tc_use;
            if (Tr >= 1.0) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_105 not valid at/above Tc");
            }
            const double exponent = 1.0 + std::pow(1.0 - Tr, D);
            return A / std::pow(B, exponent);
        }

        case CorrelationForm::EQ_106: {
            if (Tc <= 0.0) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_106 requires Tc");
            }
            const double Tr = T / Tc;
            if (Tr >= 1.0) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_106 not valid at/above Tc");
            }
            const double tau = 1.0 - Tr;
            const double exponent = (B + C * Tr + D * Tr * Tr + E * Tr * Tr * Tr);
            return A * std::pow(tau, exponent);
        }

        case CorrelationForm::EQ_107: {
            const double CT = C / T;
            const double ET = E / T;
            const double sinh_CT = std::sinh(CT);
            const double cosh_ET = std::cosh(ET);

            const double ratio1 = (std::abs(CT) < 1e-12) ? 1.0 : (CT / sinh_CT);
            const double ratio2 = (ET == 0.0) ? 0.0 : (ET / cosh_ET);
            return A + B * (ratio1 * ratio1) + D * (ratio2 * ratio2);
        }

        case CorrelationForm::EQ_114: {
            if (Tc <= 0.0) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_114 requires Tc");
            }
            const double Tr = T / Tc;
            const double tau = 1.0 - Tr;
            if (!(tau > 0.0)) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_114 not valid at/above Tc");
            }
            const double A2_over_tau = (A * A) / tau;
            const double tau2 = tau * tau;
            const double tau3 = tau2 * tau;
            const double tau4 = tau3 * tau;
            const double tau5 = tau4 * tau;
            return A2_over_tau + B - 2.0 * A * C * tau - A * D * tau2 - (1.0 / 3.0) * (C * C) * tau3 -
                   0.5 * C * D * tau4 - 0.2 * (D * D) * tau5;
        }

        case CorrelationForm::EQ_123: {
            if (Tc <= 0.0) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_123 requires Tc");
            }
            const double Tr = T / Tc;
            const double tau = 1.0 - Tr;
            if (!(tau > 0.0)) {
                throw std::runtime_error("Core::IdealGasCpCorrelation: EQ_123 not valid at/above Tc");
            }
            const double tau_035 = std::pow(tau, 0.35);
            const double tau_23 = std::pow(tau, 2.0 / 3.0);
            const double tau_43 = std::pow(tau, 4.0 / 3.0);
            return A * (1.0 + B * tau_035 + C * tau_23 + D * tau + E * tau_43);
        }

        default:
            throw std::runtime_error("Core::IdealGasCpCorrelation: unknown correlation form");
    }
}

double IdealGasCpCorrelation::calculate(double T, double Tc) const {
    if (!isValid()) {
        throw std::runtime_error("Core::IdealGasCpCorrelation: missing/invalid coefficients");
    }
    if (!inRange(T)) {
        throw std::runtime_error("Core::IdealGasCpCorrelation: temperature out of validity range");
    }
    return evaluate(T, Tc);
}

double IdealGasCpCorrelation::integrateOverT(double T, double T_ref, double Tc) const {
    if (T <= 0.0 || T_ref <= 0.0) {
        throw std::runtime_error("Core::IdealGasCpCorrelation: temperature must be positive for integration");
    }
    if (!isValid()) {
        throw std::runtime_error("Core::IdealGasCpCorrelation: cannot integrate with missing/invalid coefficients");
    }
    if (!inRange(T) || !inRange(T_ref)) {
        throw std::runtime_error("Core::IdealGasCpCorrelation: temperature out of validity range for integration");
    }
    if (T == T_ref) return 0.0;

    if (form == CorrelationForm::EQ_100) {
        return A * std::log(T / T_ref) + B * (T - T_ref) + C / 2.0 * (T * T - T_ref * T_ref) +
               D / 3.0 * (T * T * T - T_ref * T_ref * T_ref) +
               E / 4.0 * (T * T * T * T - T_ref * T_ref * T_ref * T_ref);
    }

    double a = T_ref;
    double b = T;
    double sign = 1.0;
    if (b < a) {
        std::swap(a, b);
        sign = -1.0;
    }
    constexpr int n_intervals = 256;
    const double h = (b - a) / n_intervals;
    auto f = [this, Tc](double temp) { return this->evaluate(temp, Tc) / temp; };

    double sum = f(a) + f(b);
    for (int i = 1; i < n_intervals; ++i) {
        const double x = a + h * i;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
    }
    return sign * (h / 3.0) * sum;
}

double IdealGasCpCorrelation::integrate(double T, double T_ref, double Tc) const {
    if (!isValid()) {
        throw std::runtime_error("Core::IdealGasCpCorrelation: cannot integrate with missing/invalid coefficients");
    }
    if (!inRange(T) || !inRange(T_ref)) {
        throw std::runtime_error("Core::IdealGasCpCorrelation: temperature out of validity range for integration");
    }
    if (T == T_ref) return 0.0;

    if (form == CorrelationForm::EQ_100) {
        const double T2 = T * T;
        const double T3 = T2 * T;
        const double T4 = T3 * T;
        const double T5 = T4 * T;
        const double Tr2 = T_ref * T_ref;
        const double Tr3 = Tr2 * T_ref;
        const double Tr4 = Tr3 * T_ref;
        const double Tr5 = Tr4 * T_ref;
        return A * (T - T_ref) + B / 2.0 * (T2 - Tr2) + C / 3.0 * (T3 - Tr3) + D / 4.0 * (T4 - Tr4) +
               E / 5.0 * (T5 - Tr5);
    }

    double a = T_ref;
    double b = T;
    double sign = 1.0;
    if (b < a) {
        std::swap(a, b);
        sign = -1.0;
    }
    constexpr int n_intervals = 256;
    const double h = (b - a) / n_intervals;
    auto f = [this, Tc](double temp) { return this->evaluate(temp, Tc); };

    double sum = f(a) + f(b);
    for (int i = 1; i < n_intervals; ++i) {
        const double x = a + h * i;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
    }
    return sign * (h / 3.0) * sum;
}

} // namespace Core
} // namespace DMThermo

