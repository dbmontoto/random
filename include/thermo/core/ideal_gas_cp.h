/**
 * @file ideal_gas_cp.h
 * @brief Core ideal-gas heat capacity correlations (legacy-independent)
 */

#ifndef THERMO_CORE_IDEAL_GAS_CP_H
#define THERMO_CORE_IDEAL_GAS_CP_H

#include <optional>

namespace DMThermo {
namespace Core {

/**
 * @brief Supported correlation equation forms (DIPPR-style form numbers).
 *
 * We store equation forms as their traditional numeric identifiers (e.g. 107),
 * but keep the implementation in Core so public headers do not depend on
 * `thermo/legacy/dippr/...`.
 */
enum class CorrelationForm : int {
    EQ_100 = 100,
    EQ_101 = 101,
    EQ_102 = 102,
    EQ_105 = 105,
    EQ_106 = 106,
    EQ_107 = 107,
    EQ_114 = 114,
    EQ_123 = 123
};

std::optional<CorrelationForm> correlationFormFromInt(int form);

/**
 * @brief Ideal-gas heat capacity correlation (Cp(T)).
 *
 * Units are whatever the stored coefficients imply (commonly J/(mol*K) in the new
 * databank pipeline). Validity checks and evaluation mirror the prior DIPPR
 * coefficient behavior but without exposing DIPPR types in the public API.
 */
struct IdealGasCpCorrelation {
    CorrelationForm form = static_cast<CorrelationForm>(0);
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double D = 0.0;
    double E = 0.0;
    double F = 0.0;
    double G = 0.0;
    double T_min = 0.0;     ///< Minimum valid temperature [K]
    double T_max = 1000.0;  ///< Maximum valid temperature [K]

    bool isValid() const;
    bool inRange(double T) const;

    double evaluate(double T, double Tc = 0.0) const;
    double calculate(double T, double Tc = 0.0) const;
    double integrateOverT(double T, double T_ref, double Tc = 0.0) const;
    double integrate(double T, double T_ref, double Tc = 0.0) const;
};

} // namespace Core
} // namespace DMThermo

#endif // THERMO_CORE_IDEAL_GAS_CP_H
