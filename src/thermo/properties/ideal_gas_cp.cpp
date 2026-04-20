#include "thermo/properties/ideal_gas_cp.h"

#include "thermo/legacy/dippr/equations.h"

#include <string>
#include <stdexcept>

namespace DMThermo {
namespace Properties {

namespace {

DIPPR::Coefficients toDipprCoefficients(const Core::IdealGasCpCorrelation& corr) {
    const auto eq = DIPPR::equationTypeFromForm(corr.eq_form);
    if (!eq) {
        throw std::runtime_error("Unsupported ideal-gas Cp equation form " + std::to_string(corr.eq_form));
    }
    DIPPR::Coefficients coeffs(
        *eq,
        corr.A, corr.B, corr.C, corr.D, corr.E, corr.F, corr.G,
        corr.Tmin, corr.Tmax);
    if (!coeffs.isValid()) {
        throw std::runtime_error("Invalid ideal-gas Cp correlation coefficients for eq_form " + std::to_string(corr.eq_form));
    }
    return coeffs;
}

} // namespace

double evaluateIdealGasCp(const Core::IdealGasCpCorrelation& corr, double T) {
    return toDipprCoefficients(corr).calculate(T);
}

double integrateIdealGasCp(const Core::IdealGasCpCorrelation& corr, double T, double T_ref) {
    return toDipprCoefficients(corr).integrate(T, T_ref);
}

double integrateIdealGasCpOverT(const Core::IdealGasCpCorrelation& corr, double T, double T_ref) {
    return toDipprCoefficients(corr).integrateOverT(T, T_ref);
}

} // namespace Properties
} // namespace DMThermo
