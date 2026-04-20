/**
 * @file PRVT.cpp
 * @brief Implementation of PR-VT EOS (volume-translated Peng-Robinson)
 */

#include "thermo/eos/PRVT.h"
#include "thermo/core/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DMThermo {
namespace Cubic {

namespace {
constexpr double ZC_PR = 0.307;
} // namespace

PRVTEOS::PRVTEOS(const Core::Mixture& mixture)
    : mixture_(mixture)
    , nc_(mixture_.numComponents())
    , pr_(mixture_)
{
    c_i_.assign(nc_, 0.0);

    const double R = Constants::GAS_CONSTANT;
    for (int i = 0; i < nc_; ++i) {
        const auto& comp = mixture_.component(i);
        const double Tc = comp.Tc();
        const double Pc = comp.Pc();
        const double Zc = comp.Zc();

        if (std::isfinite(Zc) && std::isfinite(Tc) && std::isfinite(Pc) && Tc > 0.0 && Pc > 0.0) {
            c_i_[i] = (Zc - ZC_PR) * (R * Tc / Pc);
        } else {
            c_i_[i] = 0.0;
        }
    }
}

double PRVTEOS::volumeTranslation(int i) const {
    if (i < 0 || i >= static_cast<int>(c_i_.size())) {
        throw std::out_of_range("Component index out of range");
    }
    return c_i_[i];
}

double PRVTEOS::cMix(const std::vector<double>& x) const {
    if (!mixture_.isValidComposition(x)) {
        throw std::invalid_argument("Invalid composition vector");
    }
    if (x.size() != c_i_.size()) {
        throw std::invalid_argument("Composition size mismatch");
    }
    double c = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        c += x[i] * c_i_[i];
    }
    return c;
}

double PRVTEOS::rhoPrimeFromRho(double rho, double c_mix) const {
    if (!(std::isfinite(rho) && rho > 0.0)) {
        throw std::invalid_argument("Density must be positive");
    }

    const double denom = 1.0 - c_mix * rho; // v' = (1 - c*rho)/rho
    if (!(std::isfinite(denom) && denom > 0.0)) {
        throw std::runtime_error("PR-VT: invalid translated volume (1 - c_mix*rho <= 0)");
    }
    return rho / denom;
}

double PRVTEOS::pressure(double T, double rho, const std::vector<double>& x) const {
    const double c_mix = cMix(x);
    const double rho_p = rhoPrimeFromRho(rho, c_mix);
    return pr_.pressure(T, rho_p, x);
}

std::vector<double> PRVTEOS::fugacityCoefficients(double T, double rho, const std::vector<double>& x) const {
    const double c_mix = cMix(x);
    const double rho_p = rhoPrimeFromRho(rho, c_mix);
    return pr_.fugacityCoefficients(T, rho_p, x);
}

double PRVTEOS::compressibility(double T, double rho, const std::vector<double>& x) const {
    const double P = pressure(T, rho, x);
    return P / (rho * Constants::GAS_CONSTANT * T);
}

double PRVTEOS::residualHelmholtz(double T, double rho, const std::vector<double>& x) const {
    const double c_mix = cMix(x);
    const double rho_p = rhoPrimeFromRho(rho, c_mix);
    return pr_.residualHelmholtz(T, rho_p, x);
}

double PRVTEOS::dadt(double T, double rho, const std::vector<double>& x) const {
    const double c_mix = cMix(x);
    const double rho_p = rhoPrimeFromRho(rho, c_mix);
    return pr_.dadt(T, rho_p, x);
}

double PRVTEOS::dadrho(double T, double rho, const std::vector<double>& x) const {
    const double c_mix = cMix(x);
    const double rho_p = rhoPrimeFromRho(rho, c_mix);

    const double denom = 1.0 - c_mix * rho;
    const double drhop_drho = 1.0 / (denom * denom);

    return pr_.dadrho(T, rho_p, x) * drhop_drho;
}

double PRVTEOS::dPdrho(double T, double rho, const std::vector<double>& x) const {
    const double c_mix = cMix(x);
    const double rho_p = rhoPrimeFromRho(rho, c_mix);

    const double denom = 1.0 - c_mix * rho;
    const double drhop_drho = 1.0 / (denom * denom);

    return pr_.dPdrho(T, rho_p, x) * drhop_drho;
}

EOSResult PRVTEOS::calculate(
    double T,
    double rho,
    const std::vector<double>& x,
    const DerivativeSpec& deriv_spec) const
{
    EOSResult out;
    out.ln_phi.assign(static_cast<size_t>(nc_), 0.0);
    out.phi.assign(static_cast<size_t>(nc_), 0.0);

    try {
        const double c_mix = cMix(x);
        const double rho_p = rhoPrimeFromRho(rho, c_mix);

        // Ask PR for derivatives at (T, rho') and then transform density derivatives by chain rule.
        const auto prr = pr_.calculate(T, rho_p, x, deriv_spec);
        if (!prr.success) {
            throw std::runtime_error("PR-VT: base Peng-Robinson calculate() failed: " + prr.error_message);
        }

        out = prr;
        out.compressibility = prr.pressure / (rho * Constants::GAS_CONSTANT * T);

        const double denom = 1.0 - c_mix * rho;
        const double drhop_drho = 1.0 / (denom * denom);

        if (deriv_spec.density) {
            const double da_drhop = prr.da_drho ? *prr.da_drho : pr_.dadrho(T, rho_p, x);
            out.da_drho = da_drhop * drhop_drho;
        }

        if (deriv_spec.second_order) {
            const double da_drhop = prr.da_drho ? *prr.da_drho : pr_.dadrho(T, rho_p, x);
            const double d2a_drhop2 = prr.d2a_drho2 ? *prr.d2a_drho2 : std::numeric_limits<double>::quiet_NaN();
            const double gpp = (2.0 * c_mix) / (denom * denom * denom); // d²rho'/drho²

            if (std::isfinite(d2a_drhop2)) {
                out.d2a_drho2 = d2a_drhop2 * (drhop_drho * drhop_drho) + da_drhop * gpp;
            } else {
                out.d2a_drho2 = std::numeric_limits<double>::quiet_NaN();
            }

            if (prr.d2a_dTdrho) {
                out.d2a_dTdrho = (*prr.d2a_dTdrho) * drhop_drho;
            }
        }

        out.success = true;
        return out;
    } catch (const std::exception& e) {
        out.success = false;
        out.error_message = e.what();
        return out;
    }
}

std::vector<double> PRVTEOS::densityRootsTP(double T, double P, const std::vector<double>& x) const {
    const double c_mix = cMix(x);

    // Roots from PR in terms of the shifted volume v' (i.e. rho').
    const auto rho_p = pr_.densityRootsTP(T, P, x);

    std::vector<double> rho;
    rho.reserve(rho_p.size());
    for (double rp : rho_p) {
        if (!(std::isfinite(rp) && rp > 0.0)) continue;
        const double denom = 1.0 + c_mix * rp; // rho = rho'/(1 + c*rho')
        if (!(std::isfinite(denom) && denom > 0.0)) continue;
        const double r = rp / denom;
        if (std::isfinite(r) && r > 0.0) rho.push_back(r);
    }

    std::sort(rho.begin(), rho.end());
    rho.erase(std::unique(rho.begin(), rho.end(), [](double a, double b) {
        return std::abs(a - b) <= 1e-12 * std::max(1.0, std::max(std::abs(a), std::abs(b)));
    }), rho.end());
    return rho;
}

double PRVTEOS::criticalTemperature(int i) const {
    return pr_.criticalTemperature(i);
}

double PRVTEOS::criticalPressure(int i) const {
    return pr_.criticalPressure(i);
}

double PRVTEOS::acentricFactor(int i) const {
    return pr_.acentricFactor(i);
}

} // namespace Cubic
} // namespace DMThermo

