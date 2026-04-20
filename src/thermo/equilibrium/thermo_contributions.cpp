/**
 * @file thermo_contributions.cpp
 * @brief Shared ideal-gas and residual contributions.
 */

#include "thermo/equilibrium/thermo_contributions.h"

#include "thermo/core/constants.h"
#include "thermo/properties/ideal_gas_cp.h"

#include <cmath>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Contributions {

namespace {

constexpr double T_REF = 298.15;
constexpr double P_REF = 1e5;

} // namespace

IdealProps idealGasProps(const Core::Mixture& mix, double T, double P, const std::vector<double>& x) {
    IdealProps out;
    const int nc = mix.numComponents();
    if (static_cast<int>(x.size()) != nc) {
        throw std::invalid_argument("idealGasProps: composition size mismatch");
    }

    const double R = Constants::GAS_CONSTANT;
    const double cp_fallback = 2.5 * R;

    double h = 0.0;
    double sT = 0.0;
    double smix = 0.0;
    for (int i = 0; i < nc; ++i) {
        const auto& comp = mix.component(i);
        if (!comp.hasIdealGasCp()) {
            throw std::runtime_error("idealGasProps: missing ideal-gas Cp for component '" + comp.name() + "'");
        }
        const auto& cp = comp.idealGasCpCorrelation();
        h += x[static_cast<size_t>(i)] * Properties::integrateIdealGasCp(cp, T, T_REF);
        sT += x[static_cast<size_t>(i)] * Properties::integrateIdealGasCpOverT(cp, T, T_REF);
        smix += (x[static_cast<size_t>(i)] > Constants::MIN_MOLE_FRACTION) ? x[static_cast<size_t>(i)] * std::log(x[static_cast<size_t>(i)]) : 0.0;
    }

    out.h_ig = h;
    out.s_ig = sT - R * std::log(std::max(P, 1.0) / P_REF) - R * smix;
    return out;
}

ResidualProps residualProps(const EOS& eos, double T, double rho, const std::vector<double>& x, double P_target) {
    ResidualProps out;
    const double R = Constants::GAS_CONSTANT;
    const int nc = eos.numComponents();
    if (static_cast<int>(x.size()) != nc) {
        throw std::invalid_argument("residualProps: composition size mismatch");
    }

    const auto res = eos.calculate(T, rho, x, DerivativeSpec{true, false, false, false});
    if (!res.success || !res.da_dT.has_value()) {
        throw std::runtime_error("residualProps: EOS derivatives unavailable");
    }

    const double a = res.a_residual;
    const double da_dT = res.da_dT.value();
    const double s_res_tv = -R * (a + T * da_dT);

    const double P = (P_target > 0.0) ? P_target : eos.pressure(T, rho, x);
    out.h_res = -R * T * T * da_dT + (P - rho * R * T) / rho;

    const double Z = P / (rho * R * T);
    out.s_res = s_res_tv + R * std::log(std::max(Z, 1e-300));

    if (res.ln_phi.size() != static_cast<size_t>(nc)) {
        throw std::runtime_error("residualProps: EOS did not provide ln_phi");
    }
    double g_res_over_rt = 0.0;
    for (int i = 0; i < nc; ++i) {
        g_res_over_rt += x[static_cast<size_t>(i)] * res.ln_phi[static_cast<size_t>(i)];
    }
    out.g_res = R * T * g_res_over_rt;
    return out;
}

} // namespace Contributions
} // namespace Equilibrium
} // namespace DMThermo
