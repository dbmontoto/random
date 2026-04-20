/**
 * @file hs_support.h
 * @brief Shared helpers for PH/PS flashes (ideal-gas Cp requirements + residual H/S).
 */

#ifndef THERMO_EQUILIBRIUM_INTERNAL_HS_SUPPORT_H
#define THERMO_EQUILIBRIUM_INTERNAL_HS_SUPPORT_H

#include "thermo/core/constants.h"
#include "thermo/core/mixture.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/thermo_contributions.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Internal {

inline const Core::Mixture& requireCoreMixture(const EOS& eos, const char* label, const char* hint = nullptr) {
    const auto* mix = eos.coreMixture();
    if (!mix) {
        std::string msg = std::string(label) + ": requires EOS exposing Core::Mixture";
        if (hint && *hint) msg += " " + std::string(hint);
        throw std::runtime_error(msg);
    }
    return *mix;
}

inline void requireIdealGasCp(const Core::Mixture& mix, const char* label) {
    std::string missing;
    for (int i = 0; i < mix.numComponents(); ++i) {
        if (mix.component(i).hasIdealGasCp()) continue;
        if (!missing.empty()) missing += ", ";
        missing += mix.component(i).name();
        const auto& cas = mix.component(i).cas();
        if (!cas.empty()) missing += " (" + cas + ")";
    }
    if (!missing.empty()) {
        throw std::runtime_error(std::string(label) + ": requires DIPPR ideal-gas Cp (ICP_*) for: " + missing);
    }
}

inline Contributions::ResidualProps residualPropsForHS(
    const EOS& eos,
    double T,
    double rho,
    const std::vector<double>& x,
    double P_target,
    const char* label)
{
    const int nc = eos.numComponents();
    if (static_cast<int>(x.size()) != nc) {
        throw std::runtime_error(std::string(label) + ": composition size mismatch");
    }
    if (!(std::isfinite(rho) && rho > 0.0)) {
        throw std::runtime_error(std::string(label) + ": invalid phase density (rho)");
    }
    const auto res = eos.calculate(T, rho, x, DerivativeSpec{true, false, false, false});
    if (!res.success || !res.da_dT.has_value()) {
        throw std::runtime_error(std::string(label) + ": requires EOS temperature derivative da_dT to compute H/S");
    }
    const double R = Constants::GAS_CONSTANT;
    Contributions::ResidualProps out;
    const double a = res.a_residual;
    const double da_dT = res.da_dT.value();
    const double P = (P_target > 0.0) ? P_target : eos.pressure(T, rho, x);
    out.h_res = -R * T * T * da_dT + (P - rho * R * T) / rho;

    const double Z = P / (rho * R * T);
    out.s_res = -R * (a + T * da_dT) + R * std::log(std::max(Z, 1e-300));

    if (res.ln_phi.size() != static_cast<size_t>(nc)) {
        throw std::runtime_error(std::string(label) + ": EOS did not provide ln_phi");
    }
    double g_res_over_rt = 0.0;
    for (int i = 0; i < nc; ++i) {
        g_res_over_rt += x[static_cast<size_t>(i)] * res.ln_phi[static_cast<size_t>(i)];
    }
    out.g_res = R * T * g_res_over_rt;
    return out;
}

} // namespace Internal
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_INTERNAL_HS_SUPPORT_H
