#include "thermo/properties/thermo_properties.h"
#include "thermo/properties/ideal_gas_cp.h"

#include "thermo/core/constants.h"

#include <cmath>
#include <stdexcept>

namespace DMThermo {
namespace Properties {

namespace {

double idealCpMix(const Core::Mixture& mix, double T, const std::vector<double>& x) {
    const int n = mix.numComponents();
    const double cp_fallback = 2.5 * Constants::GAS_CONSTANT;
    double cp = 0.0;
    for (int i = 0; i < n; ++i) {
        const auto& comp = mix.component(i);
        if (!comp.hasIdealGasCp()) {
            throw std::runtime_error("Missing ideal-gas Cp correlation for component '" + comp.name() + "'");
        }
        cp += x[i] * evaluateIdealGasCp(comp.idealGasCpCorrelation(), T);
    }
    return cp;
}

double idealH(const Core::Mixture& mix, double T, const std::vector<double>& x, double T_ref) {
    const int n = mix.numComponents();
    const double cp_fallback = 2.5 * Constants::GAS_CONSTANT;
    double h = 0.0;
    for (int i = 0; i < n; ++i) {
        const auto& comp = mix.component(i);
        if (!comp.hasIdealGasCp()) {
            throw std::runtime_error("Missing ideal-gas Cp correlation for component '" + comp.name() + "'");
        }
        h += x[i] * integrateIdealGasCp(comp.idealGasCpCorrelation(), T, T_ref);
    }
    return h;
}

double idealS_thermal(const Core::Mixture& mix, double T, const std::vector<double>& x, double T_ref) {
    const int n = mix.numComponents();
    const double cp_fallback = 2.5 * Constants::GAS_CONSTANT;
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        const auto& comp = mix.component(i);
        if (!comp.hasIdealGasCp()) {
            throw std::runtime_error("Missing ideal-gas Cp correlation for component '" + comp.name() + "'");
        }
        s += x[i] * integrateIdealGasCpOverT(comp.idealGasCpCorrelation(), T, T_ref);
    }
    return s;
}

double idealMixingEntropy(const std::vector<double>& x) {
    const double R = Constants::GAS_CONSTANT;
    double smix = 0.0;
    for (double xi : x) {
        if (xi > 0.0) {
            smix += -R * xi * std::log(xi);
        }
    }
    return smix;
}

} // namespace

ThermoProperties calculateTV(
    const EOS& eos,
    const Core::Mixture& mixture,
    double T,
    double rho,
    const std::vector<double>& x,
    const ReferenceState& ref)
{
    ThermoProperties out;
    out.T = T;
    out.rho = rho;
    out.v = (rho > 0.0) ? (1.0 / rho) : 0.0;
    out.model = eos.name();

    try {
        if (rho <= 0.0) {
            throw std::invalid_argument("rho must be positive");
        }
        if (!mixture.isValidComposition(x)) {
            throw std::invalid_argument("invalid composition vector");
        }

        // EOS residuals (dimensionless a_res = A_res/(RT)).
        //
        // Cp/Cv require second derivatives of a_res with respect to (T, rho). Not all EOS models
        // implement second derivatives; when unavailable, fall back to ideal-gas Cp/Cv while still
        // computing H/S/U/A/G from first derivatives.
        EOSResult r = eos.calculate(T, rho, x, DerivativeSpec{true, true, false, true});
        if (!r.success) {
            r = eos.calculate(T, rho, x, DerivativeSpec{true, true, false, false});
        }
        if (!r.success) {
            throw std::runtime_error("EOS calculation failed");
        }
        if (!r.da_dT.has_value()) {
            throw std::runtime_error("EOS did not provide da_dT (required for S/U/H residuals)");
        }

        out.P = r.pressure;
        out.Z = r.compressibility;

        const double R = Constants::GAS_CONSTANT;

        // Ideal-gas contributions:
        // - Use Cp(T) correlations per component; reference at (T_ref, P_ref).
        // - s_ig(T,P,x) = integral(Cp/T)dT - R ln(P/P_ref) + s_mix
        const double H_ig = idealH(mixture, T, x, ref.T_ref);
        const double S_ig = idealS_thermal(mixture, T, x, ref.T_ref)
                          - R * std::log(std::max(out.P, 1.0) / ref.P_ref)
                          + idealMixingEntropy(x);
        const double Cp_ig = idealCpMix(mixture, T, x);
        const double Cv_ig = Cp_ig - R;

        // Residual contributions from reduced Helmholtz energy:
        // A_res = a_res * R * T
        // S_res(T,P) = S - S_ig(T,P) = S_res(T,V) + R ln(Z)
        //           = -R*(a_res + T*da/dT)|rho + R ln(Z)
        // U_res = A_res + T*S_res = -R*T^2*da/dT
        // H_res = U_res + (P - rho*R*T)/rho = U_res + R*T*(Z-1)
        const double a_res = r.a_residual;
        const double da_dT = r.da_dT.value();
        const double S_res = -R * (a_res + T * da_dT) + R * std::log(std::max(out.Z, 1e-300));
        const double U_res = -R * T * T * da_dT;
        const double H_res = U_res + R * T * (out.Z - 1.0);

        out.H = H_ig + H_res;
        out.S = S_ig + S_res;
        out.U = (H_ig - R * T) + U_res; // U_ig = H_ig - R*T
        out.G = out.H - T * out.S;
        out.A = out.U - T * out.S;

        // Heat capacities.
        //
        // Cv_res from U_res = -R*T^2*da/dT:
        //   Cv_res = (dU_res/dT)|rho = -R*(2T*da/dT + T^2*d2a/dT2)
        double Cv_total = Cv_ig;
        if (r.d2a_dT2.has_value()) {
            Cv_total += -R * (2.0 * T * da_dT + T * T * r.d2a_dT2.value());
        }
        out.Cv = Cv_total;

        // Cp from the thermodynamic identity:
        //   Cp = Cv + T*(dP/dT|rho)^2 / (rho^2 * dP/drho|T)
        //
        // with:
        //   P = rho*R*T + rho^2*R*T*da/drho
        //   dP/dT|rho = rho*R*Z + rho^2*R*T*d2a/(dT drho)
        //   dP/drho|T = R*T + 2*rho*R*T*da/drho + rho^2*R*T*d2a/drho2
        double Cp_total = Cp_ig;
        if (r.da_drho.has_value() && r.d2a_dTdrho.has_value() && r.d2a_drho2.has_value() && rho > 0.0) {
            const double da_drho = r.da_drho.value();
            const double d2a_dTdrho = r.d2a_dTdrho.value();
            const double d2a_drho2 = r.d2a_drho2.value();

            const double dP_dT_rho = rho * R * out.Z + rho * rho * R * T * d2a_dTdrho;
            const double dP_drho_T = R * T + 2.0 * rho * R * T * da_drho + rho * rho * R * T * d2a_drho2;

            const double denom = rho * rho * dP_drho_T;
            if (std::isfinite(dP_dT_rho) && std::isfinite(denom) && denom > 0.0) {
                const double delta = T * (dP_dT_rho * dP_dT_rho) / denom;
                const double candidate = Cv_total + delta;
                if (std::isfinite(candidate)) {
                    Cp_total = candidate;
                }
            }
        }
        out.Cp = Cp_total;

        out.success = true;
        return out;

    } catch (const std::exception& e) {
        out.success = false;
        out.message = e.what();
        return out;
    }
}

} // namespace Properties
} // namespace DMThermo
