/**
 * @file SRK.cpp
 * @brief Implementation of Soave-Redlich-Kwong EOS (cubic)
 */

#include "thermo/eos/SRK.h"
#include "thermo/numerics/polynomial_roots.h"
#include <cmath>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <limits>

namespace DMThermo {
namespace Cubic {

namespace {
constexpr double SRK_OMEGA_A = 0.42747;
constexpr double SRK_OMEGA_B = 0.08664;
} // namespace

SRKEOS::SRKEOS(const Core::Mixture& mixture)
    : mixture_(mixture)
    , nc_(mixture_.numComponents())
{
    initializeCriticalCache();
}

void SRKEOS::initializeCriticalCache() {
    Tc_.assign(nc_, 0.0);
    Pc_.assign(nc_, 0.0);
    omega_.assign(nc_, 0.0);

    for (int i = 0; i < nc_; ++i) {
        const auto& comp = mixture_.component(i);

        double Tc = comp.Tc();
        double Pc = comp.Pc();
        double omega = comp.omega();

        if (!(std::isfinite(Tc) && std::isfinite(Pc) && std::isfinite(omega)) || Tc <= 0.0 || Pc <= 0.0) {
            throw std::runtime_error("SRK EOS requires Tc/Pc/omega from databanks for component '" + comp.name() + "'");
        }

        Tc_[i] = Tc;
        Pc_[i] = Pc;
        omega_[i] = omega;
    }
}

SRKEOS::MixParams SRKEOS::mixtureParams(double T, const std::vector<double>& x, bool compute_temperature_derivatives) const {
    if (!mixture_.isValidComposition(x)) {
        throw std::invalid_argument("Invalid composition vector");
    }

    MixParams mp;
    mp.a_i.assign(nc_, 0.0);
    mp.b_i.assign(nc_, 0.0);
    mp.a_i_mix.assign(nc_, 0.0);
    if (compute_temperature_derivatives) {
        mp.da_i_dT.assign(nc_, 0.0);
        mp.d2a_i_dT2.assign(nc_, 0.0);
        mp.da_mix_dT = 0.0;
        mp.d2a_mix_dT2 = 0.0;
    }

    const double R = Constants::GAS_CONSTANT;

    for (int i = 0; i < nc_; ++i) {
        const double Tc = Tc_[i];
        const double Pc = Pc_[i];
        const double omega = omega_[i];

        const double Tr = T / Tc;
        const double kappa = 0.480 + 1.574 * omega - 0.176 * omega * omega;
        const double sqrt_Tr = std::sqrt(std::max(Tr, 1e-16));
        const double alpha = (1.0 + kappa * (1.0 - sqrt_Tr));
        const double alpha2 = alpha * alpha;

        const double a0 = SRK_OMEGA_A * R * R * Tc * Tc / Pc;
        mp.a_i[i] = a0 * alpha2;
        if (compute_temperature_derivatives) {
            const double dalpha_dT = -kappa / (2.0 * Tc * sqrt_Tr);
            const double d2alpha_dT2 = kappa / (4.0 * Tc * Tc * sqrt_Tr * sqrt_Tr * sqrt_Tr);
            mp.da_i_dT[i] = a0 * 2.0 * alpha * dalpha_dT;
            mp.d2a_i_dT2[i] = a0 * 2.0 * (dalpha_dT * dalpha_dT + alpha * d2alpha_dT2);
        }

        mp.b_i[i] = SRK_OMEGA_B * R * Tc / Pc;
        mp.b_mix += x[i] * mp.b_i[i];
    }

    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            const double kij = mixture_.kij(i, j);
            const double aij = std::sqrt(mp.a_i[i] * mp.a_i[j]) * (1.0 - kij);
            mp.a_mix += x[i] * x[j] * aij;
            mp.a_i_mix[i] += x[j] * aij;
            if (compute_temperature_derivatives) {
                const double ai = mp.a_i[i];
                const double aj = mp.a_i[j];
                if (ai > 0.0 && aj > 0.0) {
                    const double dai = mp.da_i_dT[i];
                    const double daj = mp.da_i_dT[j];
                    const double si = dai / ai;
                    const double sj = daj / aj;
                    const double term = 0.5 * (si + sj);
                    mp.da_mix_dT += x[i] * x[j] * aij * term;

                    const double d2ai = mp.d2a_i_dT2[i];
                    const double d2aj = mp.d2a_i_dT2[j];
                    const double ti = d2ai / ai;
                    const double tj = d2aj / aj;
                    const double factor = 0.5 * (ti + tj) + 0.5 * si * sj - 0.25 * (si * si + sj * sj);
                    mp.d2a_mix_dT2 += x[i] * x[j] * aij * factor;
                }
            }
        }
    }

    return mp;
}

double SRKEOS::pressureFromMix(double T, double v, double a, double b) const {
    const double R = Constants::GAS_CONSTANT;
    if (v <= b) {
        return std::numeric_limits<double>::infinity();
    }

    const double denom = v * (v + b);
    return R * T / (v - b) - a / denom;
}

std::vector<double> SRKEOS::lnPhi(double T, double rho, const std::vector<double>& x, double& Z_out) const {
    const auto mp = mixtureParams(T, x);

    const double v = 1.0 / rho;
    const double P = pressureFromMix(T, v, mp.a_mix, mp.b_mix);
    const double R = Constants::GAS_CONSTANT;
    const double Z = P * v / (R * T);

    if (!std::isfinite(P) || !std::isfinite(Z) || P <= 0.0) {
        throw std::runtime_error("SRK EOS: invalid state (non-finite pressure/Z)");
    }

    const double A = mp.a_mix * P / (R * R * T * T);
    const double B = mp.b_mix * P / (R * T);

    if (Z <= B) {
        throw std::runtime_error("SRK EOS: invalid state (Z <= B)");
    }

    const double log_term = std::log(1.0 + B / Z);

    std::vector<double> ln_phi(nc_, 0.0);
    for (int i = 0; i < nc_; ++i) {
        const double bi_over_b = mp.b_i[i] / mp.b_mix;
        const double da_over_a = (mp.a_mix > 0.0) ? (2.0 * mp.a_i_mix[i] / mp.a_mix) : 0.0;
        ln_phi[i] = bi_over_b * (Z - 1.0)
                  - std::log(Z - B)
                  - (A / B) * (da_over_a - bi_over_b) * log_term;
    }

    Z_out = Z;
    return ln_phi;
}

double SRKEOS::residualHelmholtzFromLnPhi(
    double Z,
    const std::vector<double>& x,
    const std::vector<double>& ln_phi) const
{
    double g_res_over_rt = 0.0;
    for (int i = 0; i < nc_; ++i) {
        g_res_over_rt += x[i] * ln_phi[i];
    }
    // a_res/(RT) at fixed (T, rho) relates to the (T, P) fugacity coefficients via:
    //   g_res/(RT) = sum_i x_i ln(phi_i)
    //   a_res/(RT) = g_res/(RT) - Z + 1 + ln(Z)
    return g_res_over_rt - (Z - 1.0) + std::log(std::max(Z, 1e-300));
}

EOSResult SRKEOS::calculate(
    double T,
    double rho,
    const std::vector<double>& x,
    const DerivativeSpec& deriv_spec) const
{
    EOSResult result;
    result.ln_phi.assign(nc_, 0.0);
    result.phi.assign(nc_, 0.0);

    try {
        if (rho <= 0.0) {
            throw std::invalid_argument("Density must be positive");
        }
        if (!mixture_.isValidComposition(x)) {
            throw std::invalid_argument("Invalid composition vector");
        }

        double Z = 0.0;
        result.ln_phi = lnPhi(T, rho, x, Z);
        for (int i = 0; i < nc_; ++i) {
            result.phi[i] = std::exp(result.ln_phi[i]);
        }

        result.pressure = pressure(T, rho, x);
        result.compressibility = Z;
        result.a_residual = residualHelmholtzFromLnPhi(Z, x, result.ln_phi);

        if (deriv_spec.composition) {
            result.chemical_potential = result.ln_phi;
        }

        if (deriv_spec.temperature) {
            result.da_dT = dadt(T, rho, x);
        }
        if (deriv_spec.density) {
            result.da_drho = dadrho(T, rho, x);
        }
        if (deriv_spec.second_order) {
            const auto mp = mixtureParams(T, x, true);
            const double b = mp.b_mix;
            const double a = mp.a_mix;
            const double da_dT = mp.da_mix_dT;
            const double d2a_dT2 = mp.d2a_mix_dT2;
            const double R = Constants::GAS_CONSTANT;

            if (!(std::isfinite(a) && std::isfinite(b) && std::isfinite(da_dT) && std::isfinite(d2a_dT2)) || b <= 0.0) {
                throw std::runtime_error("SRKEOS::calculate: invalid mixture parameters for second derivatives");
            }

            const double brho = b * rho;
            const double one_minus_brho = 1.0 - brho;
            const double one_plus_brho = 1.0 + brho;
            if (!(std::isfinite(one_minus_brho) && std::isfinite(one_plus_brho) && one_minus_brho > 0.0 && one_plus_brho > 0.0)) {
                throw std::runtime_error("SRKEOS::calculate: invalid log arguments for second derivatives");
            }

            const double log_term = std::log(one_plus_brho);

            const double invT = 1.0 / T;
            const double invT2 = invT * invT;
            const double invT3 = invT2 * invT;
            const double g1 = da_dT * invT - a * invT2;
            const double g2 = d2a_dT2 * invT - 2.0 * da_dT * invT2 + 2.0 * a * invT3;

            result.d2a_dT2 = -(log_term / (b * R)) * g2;
            result.d2a_drho2 = (b * b) / (one_minus_brho * one_minus_brho) + (a * b / (R * T)) / (one_plus_brho * one_plus_brho);
            result.d2a_dTdrho = -(1.0 / (R * one_plus_brho)) * g1;
        }

        result.success = true;
        return result;
    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = e.what();
        return result;
    }
}

double SRKEOS::pressure(double T, double rho, const std::vector<double>& x) const {
    const auto mp = mixtureParams(T, x);
    const double v = 1.0 / rho;
    return pressureFromMix(T, v, mp.a_mix, mp.b_mix);
}

std::vector<double> SRKEOS::fugacityCoefficients(double T, double rho, const std::vector<double>& x) const {
    double Z = 0.0;
    const auto ln_phi = lnPhi(T, rho, x, Z);
    std::vector<double> phi(nc_, 0.0);
    for (int i = 0; i < nc_; ++i) {
        phi[i] = std::exp(ln_phi[i]);
    }
    return phi;
}

double SRKEOS::compressibility(double T, double rho, const std::vector<double>& x) const {
    const double R = Constants::GAS_CONSTANT;
    const double v = 1.0 / rho;
    const double P = pressure(T, rho, x);
    return P * v / (R * T);
}

double SRKEOS::residualHelmholtz(double T, double rho, const std::vector<double>& x) const {
    double Z = 0.0;
    const auto ln_phi = lnPhi(T, rho, x, Z);
    return residualHelmholtzFromLnPhi(Z, x, ln_phi);
}

double SRKEOS::dadt(double T, double rho, const std::vector<double>& x) const {
    if (!(std::isfinite(T) && std::isfinite(rho)) || T <= 0.0 || rho <= 0.0) {
        throw std::invalid_argument("SRKEOS::dadt: invalid T/rho");
    }

    const auto mp = mixtureParams(T, x, true);
    const double b = mp.b_mix;
    const double a = mp.a_mix;
    const double da_dT = mp.da_mix_dT;
    const double R = Constants::GAS_CONSTANT;

    if (!(std::isfinite(a) && std::isfinite(b) && std::isfinite(da_dT)) || b <= 0.0) {
        throw std::runtime_error("SRKEOS::dadt: invalid mixture parameters");
    }

    const double brho = b * rho;
    const double arg = 1.0 + brho;
    if (!(std::isfinite(arg) && arg > 0.0)) {
        throw std::runtime_error("SRKEOS::dadt: invalid log arguments");
    }

    const double log_term = std::log(arg);
    const double invT = 1.0 / T;
    const double denom = b * R;
    const double value = -(log_term / denom) * (da_dT * invT - a * invT * invT);
    if (!std::isfinite(value)) {
        throw std::runtime_error("SRKEOS::dadt: non-finite result");
    }
    return value;
}

double SRKEOS::dadrho(double T, double rho, const std::vector<double>& x) const {
    if (!(std::isfinite(T) && std::isfinite(rho)) || T <= 0.0 || rho <= 0.0) {
        throw std::invalid_argument("SRKEOS::dadrho: invalid T/rho");
    }

    const auto mp = mixtureParams(T, x);
    const double b = mp.b_mix;
    const double a = mp.a_mix;
    const double R = Constants::GAS_CONSTANT;

    if (!(std::isfinite(a) && std::isfinite(b)) || b <= 0.0) {
        throw std::runtime_error("SRKEOS::dadrho: invalid mixture parameters");
    }

    const double brho = b * rho;
    if (!(std::isfinite(brho) && brho < 1.0)) {
        throw std::runtime_error("SRKEOS::dadrho: invalid packing (b*rho >= 1)");
    }

    const double one_minus_brho = 1.0 - brho;
    const double one_plus_brho = 1.0 + brho;
    if (!(std::isfinite(one_minus_brho) && std::isfinite(one_plus_brho) && one_minus_brho > 0.0 && one_plus_brho > 0.0)) {
        throw std::runtime_error("SRKEOS::dadrho: invalid log arguments");
    }

    const double term1 = b / one_minus_brho;
    const double term2 = -(a / (R * T)) / one_plus_brho;
    const double value = term1 + term2;
    if (!std::isfinite(value)) {
        throw std::runtime_error("SRKEOS::dadrho: non-finite result");
    }
    return value;
}

double SRKEOS::dPdrho(double T, double rho, const std::vector<double>& x) const {
    const auto mp = mixtureParams(T, x);
    const double v = 1.0 / rho;
    const double b = mp.b_mix;
    const double a = mp.a_mix;
    const double R = Constants::GAS_CONSTANT;

    if (v <= b) {
        return std::numeric_limits<double>::infinity();
    }

    const double denom = v * (v + b);
    const double d_denom_dv = 2.0 * v + b;
    const double dPdv = -R * T / ((v - b) * (v - b)) + a * d_denom_dv / (denom * denom);
    const double dvdrho = -1.0 / (rho * rho);
    return dPdv * dvdrho;
}

std::vector<double> SRKEOS::compressibilityRootsTP(double T, double P, const std::vector<double>& x) const {
    if (!(std::isfinite(T) && std::isfinite(P) && T > 0.0 && P > 0.0)) {
        throw std::invalid_argument("SRKEOS::compressibilityRootsTP: invalid T/P");
    }

    const auto mp = mixtureParams(T, x);
    const double R = Constants::GAS_CONSTANT;

    const double A = mp.a_mix * P / (R * R * T * T);
    const double B = mp.b_mix * P / (R * T);

    // SRK cubic in Z:
    //   Z^3 - Z^2 + (A - B - B^2) Z - A B = 0
    const double c3 = 1.0;
    const double c2 = -1.0;
    const double c1 = A - B - B * B;
    const double c0 = -A * B;

    std::vector<double> Z = Numerics::RootFinding::realRootsCubic(c3, c2, c1, c0);
    const double eps = 1e-12 * std::max(1.0, std::abs(B));
    Z.erase(std::remove_if(Z.begin(), Z.end(), [&](double z) {
        return !(std::isfinite(z) && (z > B + eps));
    }), Z.end());
    std::sort(Z.begin(), Z.end());
    Z.erase(std::unique(Z.begin(), Z.end(), [](double a, double b) {
        return std::abs(a - b) <= 1e-12 * std::max(1.0, std::max(std::abs(a), std::abs(b)));
    }), Z.end());
    return Z;
}

std::vector<double> SRKEOS::densityRootsTP(double T, double P, const std::vector<double>& x) const {
    const double R = Constants::GAS_CONSTANT;
    const auto Z = compressibilityRootsTP(T, P, x);

    std::vector<double> rho;
    rho.reserve(Z.size());
    for (double z : Z) {
        const double r = P / (std::max(z, 1e-300) * R * T);
        if (std::isfinite(r) && r > 0.0) rho.push_back(r);
    }
    std::sort(rho.begin(), rho.end());
    rho.erase(std::unique(rho.begin(), rho.end(), [](double a, double b) {
        return std::abs(a - b) <= 1e-12 * std::max(1.0, std::max(std::abs(a), std::abs(b)));
    }), rho.end());
    return rho;
}

double SRKEOS::criticalTemperature(int i) const {
    if (i < 0 || i >= nc_) throw std::out_of_range("Component index out of range");
    return Tc_[i];
}

double SRKEOS::criticalPressure(int i) const {
    if (i < 0 || i >= nc_) throw std::out_of_range("Component index out of range");
    return Pc_[i];
}

double SRKEOS::acentricFactor(int i) const {
    if (i < 0 || i >= nc_) throw std::out_of_range("Component index out of range");
    return omega_[i];
}

} // namespace Cubic
} // namespace DMThermo
