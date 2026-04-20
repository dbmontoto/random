/**
 * @file VTPR.cpp
 * @brief Implementation of VTPR EOS (PR + Wong-Sandler + UNIFAC + volume translation)
 */

#include "thermo/eos/VTPR.h"

#include "thermo/activity/unifac.h"
#include "thermo/numerics/polynomial_roots.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DMThermo {
namespace Cubic {

namespace {
constexpr double PR_OMEGA_A = 0.45724;
constexpr double PR_OMEGA_B = 0.07780;
constexpr double SQRT2 = 1.4142135623730950488;

// Wong-Sandler constant for PR:
//   C* = (1/(2*sqrt(2))) * ln( (1+sqrt(2)) / (sqrt(2)-1) ) = ln(3+2*sqrt(2)) / (2*sqrt(2))
constexpr double C_STAR_PR = 0.6232252401402304402;
} // namespace

VTPREOS::VTPREOS(const Core::Mixture& mixture)
    : mixture_(mixture)
    , nc_(mixture_.numComponents())
{
    initializeCriticalCache();

    // VTPR requires UNIFAC tables and per-component VTPR inputs.
    if (!mixture_.unifacTables()) {
        throw std::runtime_error("VTPR EOS requires UNIFAC tables on the mixture (unifac_*.csv)");
    }
    for (int i = 0; i < nc_; ++i) {
        const auto& c = mixture_.component(i);
        if (!c.hasVtprVolumeTranslation()) {
            throw std::runtime_error("VTPR EOS requires c(T) volume translation for component '" + c.name() + "' (vtpr_pure.csv)");
        }
        if (c.unifacSubgroups().empty()) {
            throw std::runtime_error("VTPR EOS requires UNIFAC subgroup mapping for component '" + c.name() + "' (vtpr_groups.csv)");
        }
    }
}

void VTPREOS::initializeCriticalCache() {
    Tc_.assign(nc_, 0.0);
    Pc_.assign(nc_, 0.0);
    omega_.assign(nc_, 0.0);

    for (int i = 0; i < nc_; ++i) {
        const auto& comp = mixture_.component(i);
        const double Tc = comp.Tc();
        const double Pc = comp.Pc();
        const double omega = comp.omega();

        if (!(std::isfinite(Tc) && std::isfinite(Pc) && std::isfinite(omega)) || Tc <= 0.0 || Pc <= 0.0) {
            throw std::runtime_error("VTPR EOS requires Tc/Pc/omega for component '" + comp.name() + "'");
        }

        Tc_[i] = Tc;
        Pc_[i] = Pc;
        omega_[i] = omega;
    }
}

double VTPREOS::pressureFromMix(double T, double v, double a, double b) const {
    const double R = Constants::GAS_CONSTANT;
    if (v <= b) return std::numeric_limits<double>::infinity();
    const double denom = v * v + 2.0 * b * v - b * b;
    return R * T / (v - b) - a / denom;
}

VTPREOS::WSParams VTPREOS::wsParams(double T, const std::vector<double>& x, bool need_gamma) const {
    if (!mixture_.isValidComposition(x)) {
        throw std::invalid_argument("Invalid composition vector");
    }
    if (!(std::isfinite(T) && T > 0.0)) {
        throw std::invalid_argument("VTPR: invalid temperature");
    }

    WSParams p;
    p.a_i.assign(nc_, 0.0);
    p.b_i.assign(nc_, 0.0);
    if (need_gamma) {
        p.ln_gamma = Activity::lnGammaUNIFAC(T, mixture_, x);
        if (static_cast<int>(p.ln_gamma.size()) != nc_) {
            throw std::runtime_error("VTPR: UNIFAC ln(gamma) size mismatch");
        }
    } else {
        p.ln_gamma.assign(nc_, 0.0);
    }

    const double R = Constants::GAS_CONSTANT;

    // Volume translation c_mix(T).
    p.c_mix = 0.0;
    for (int i = 0; i < nc_; ++i) {
        p.c_mix += x[i] * mixture_.component(i).vtprC(T);
    }

    // PR pure-component parameters.
    for (int i = 0; i < nc_; ++i) {
        const double Tc = Tc_[i];
        const double Pc = Pc_[i];
        const double omega = omega_[i];

        const double Tr = T / Tc;
        const double kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        const double sqrt_Tr = std::sqrt(std::max(Tr, 1e-16));
        const double alpha = (1.0 + kappa * (1.0 - sqrt_Tr));
        const double alpha2 = alpha * alpha;

        const double a0 = PR_OMEGA_A * R * R * Tc * Tc / Pc;
        p.a_i[i] = a0 * alpha2;
        p.b_i[i] = PR_OMEGA_B * R * Tc / Pc;
    }

    // Wong-Sandler D = Σ x_i a_i/(b_i RT) + (Σ x_i ln γ_i)/C*
    double sum_a_over_brt = 0.0;
    double sum_x_lngamma = 0.0;
    for (int i = 0; i < nc_; ++i) {
        if (!(std::isfinite(p.b_i[i]) && p.b_i[i] > 0.0)) {
            throw std::runtime_error("VTPR: invalid PR b_i for component index " + std::to_string(i));
        }
        sum_a_over_brt += x[i] * (p.a_i[i] / (p.b_i[i] * R * T));
        sum_x_lngamma += x[i] * p.ln_gamma[i];
    }
    p.D = sum_a_over_brt + (sum_x_lngamma / C_STAR_PR);

    const double one_minus_D = 1.0 - p.D;
    if (!(std::isfinite(one_minus_D) && std::abs(one_minus_D) > 1e-14)) {
        throw std::runtime_error("VTPR: invalid 1-D in Wong-Sandler mixing");
    }

    // Q = Σ_i Σ_j x_i x_j [ b_ij - a_ij/(RT) ], where:
    //   b_ij = (b_i + b_j)/2
    //   a_ij = sqrt(a_i a_j) (1 - k_ij)
    p.Q = 0.0;
    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            const double bij = 0.5 * (p.b_i[i] + p.b_i[j]);
            const double kij = mixture_.kij(i, j);
            const double aij = std::sqrt(std::max(p.a_i[i] * p.a_i[j], 0.0)) * (1.0 - kij);
            const double Bij = bij - aij / (R * T);
            p.Q += x[i] * x[j] * Bij;
        }
    }

    p.b = p.Q / one_minus_D;
    p.a = R * T * p.Q * p.D / one_minus_D;

    if (!(std::isfinite(p.a) && std::isfinite(p.b)) || p.b <= 0.0) {
        throw std::runtime_error("VTPR: invalid mixed a/b from Wong-Sandler");
    }

    return p;
}

double VTPREOS::rhoEosFromRho(double rho, double c_mix) const {
    if (!(std::isfinite(rho) && rho > 0.0)) {
        throw std::invalid_argument("VTPR: density must be positive");
    }
    const double denom = 1.0 + c_mix * rho; // v_eos = v + c_mix
    if (!(std::isfinite(denom) && denom > 0.0)) {
        throw std::runtime_error("VTPR: invalid translated volume (1 + c_mix*rho <= 0)");
    }
    return rho / denom;
}

double VTPREOS::rhoFromRhoEos(double rho_eos, double c_mix) const {
    if (!(std::isfinite(rho_eos) && rho_eos > 0.0)) {
        throw std::invalid_argument("VTPR: rho_eos must be positive");
    }
    const double denom = 1.0 - c_mix * rho_eos; // v = v_eos - c_mix
    if (!(std::isfinite(denom) && denom > 0.0)) {
        throw std::runtime_error("VTPR: invalid back-translated volume (1 - c_mix*rho_eos <= 0)");
    }
    return rho_eos / denom;
}

std::vector<double> VTPREOS::lnPhi(double T, double rho_eos, const std::vector<double>& x, double& Z_out) const {
    const auto mp = wsParams(T, x, true);

    const double v = 1.0 / rho_eos;
    const double P = pressureFromMix(T, v, mp.a, mp.b);
    const double R = Constants::GAS_CONSTANT;
    const double Z = P * v / (R * T);

    if (!std::isfinite(P) || !std::isfinite(Z) || P <= 0.0) {
        throw std::runtime_error("VTPR: invalid state (non-finite pressure/Z)");
    }

    const double A = mp.a * P / (R * R * T * T);
    const double B = mp.b * P / (R * T);
    if (Z <= B) {
        throw std::runtime_error("VTPR: invalid state (Z <= B)");
    }

    const double num = Z + (1.0 + SQRT2) * B;
    const double den = Z + (1.0 - SQRT2) * B;
    if (!(std::isfinite(num) && std::isfinite(den) && num > 0.0 && den > 0.0)) {
        throw std::runtime_error("VTPR: invalid log arguments for Cbar");
    }
    const double Cbar = (1.0 / (2.0 * SQRT2)) * std::log(num / den);

    // Precompute Bij = b_ij - a_ij/(RT) for S_Q_i.
    std::vector<std::vector<double>> Bij(nc_, std::vector<double>(nc_, 0.0));
    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            const double bij = 0.5 * (mp.b_i[i] + mp.b_i[j]);
            const double kij = mixture_.kij(i, j);
            const double aij = std::sqrt(std::max(mp.a_i[i] * mp.a_i[j], 0.0)) * (1.0 - kij);
            Bij[i][j] = bij - aij / (R * T);
        }
    }

    const double one_minus_D = 1.0 - mp.D;
    if (!(std::isfinite(one_minus_D) && std::abs(one_minus_D) > 1e-14)) {
        throw std::runtime_error("VTPR: invalid 1-D");
    }

    std::vector<double> ln_phi(nc_, 0.0);
    for (int i = 0; i < nc_; ++i) {
        // ∂(ND)/∂N_i
        const double S_D_i = mp.a_i[i] / (mp.b_i[i] * R * T) + mp.ln_gamma[i] / C_STAR_PR;

        // (1/N) ∂(N^2 Q)/∂N_i = 2 Σ_j x_j (b - a/RT)_ij
        double S_Q_i = 0.0;
        for (int j = 0; j < nc_; ++j) {
            S_Q_i += x[j] * Bij[i][j];
        }
        S_Q_i *= 2.0;

        // f1 = S_Q_i/(1-D) - Q/(1-D)^2 * (1 - S_D_i)
        const double f1 = (S_Q_i / one_minus_D) - (mp.Q / (one_minus_D * one_minus_D)) * (1.0 - S_D_i);

        // (1/N) ∂(N^2 a_m)/∂N_i = RT (D b_i + b_m ∂(ND)/∂N_i)
        const double S_a_i = R * T * (mp.D * mp.b_i[i] + mp.b * S_D_i);

        const double f2 = (f1 / mp.b) - (S_a_i / mp.a);

        ln_phi[i] = (f1 / mp.b) * (Z - 1.0)
                  - std::log(Z - B)
                  + (A / B) * f2 * Cbar;
    }

    Z_out = Z;
    return ln_phi;
}

double VTPREOS::pressure(double T, double rho, const std::vector<double>& x) const {
    const auto mp = wsParams(T, x, false);
    const double rho_eos = rhoEosFromRho(rho, mp.c_mix);
    const double v = 1.0 / rho_eos;
    return pressureFromMix(T, v, mp.a, mp.b);
}

std::vector<double> VTPREOS::fugacityCoefficients(double T, double rho, const std::vector<double>& x) const {
    const auto mp = wsParams(T, x, true);
    const double rho_eos = rhoEosFromRho(rho, mp.c_mix);

    double Z = 0.0;
    const auto ln_phi = lnPhi(T, rho_eos, x, Z);
    std::vector<double> phi(nc_, 0.0);
    for (int i = 0; i < nc_; ++i) phi[i] = std::exp(ln_phi[i]);
    return phi;
}

double VTPREOS::compressibility(double T, double rho, const std::vector<double>& x) const {
    (void)x;
    const double P = pressure(T, rho, x);
    return P / (rho * Constants::GAS_CONSTANT * T);
}

double VTPREOS::residualHelmholtz(double T, double rho, const std::vector<double>& x) const {
    const auto mp = wsParams(T, x, false);
    const double rho_eos = rhoEosFromRho(rho, mp.c_mix);
    const double b = mp.b;
    const double a = mp.a;
    const double R = Constants::GAS_CONSTANT;

    const double brho = b * rho_eos;
    if (!(std::isfinite(brho) && brho < 1.0)) {
        throw std::runtime_error("VTPR: invalid packing (b*rho >= 1)");
    }

    const double one_minus_brho = 1.0 - brho;
    const double cp = 1.0 + SQRT2;
    const double cm = 1.0 - SQRT2;
    const double num = 1.0 + cp * brho;
    const double den = 1.0 + cm * brho;
    if (!(std::isfinite(num) && std::isfinite(den) && num > 0.0 && den > 0.0)) {
        throw std::runtime_error("VTPR: invalid log arguments");
    }

    const double term1 = -std::log(one_minus_brho);
    const double log_term = std::log(num / den);
    const double term2 = -(a / (2.0 * SQRT2 * b * R * T)) * log_term;
    return term1 + term2;
}

double VTPREOS::dadt(double T, double rho, const std::vector<double>& x) const {
    // Numerical derivative: VTPR a(T, rho, x) includes UNIFAC temperature dependence and
    // temperature-dependent volume translation c(T), so keep this robust for now.
    const double h = std::max(1e-6 * T, 1e-3);
    const double Tp = T + h;
    const double Tm = T - h;
    if (Tm <= 0.0) {
        const double ap = residualHelmholtz(Tp, rho, x);
        const double a0 = residualHelmholtz(T, rho, x);
        return (ap - a0) / (Tp - T);
    }
    const double ap = residualHelmholtz(Tp, rho, x);
    const double am = residualHelmholtz(Tm, rho, x);
    return (ap - am) / (Tp - Tm);
}

double VTPREOS::dadrho(double T, double rho, const std::vector<double>& x) const {
    const auto mp = wsParams(T, x, false);
    const double rho_eos = rhoEosFromRho(rho, mp.c_mix);
    const double b = mp.b;
    const double a = mp.a;
    const double R = Constants::GAS_CONSTANT;

    const double brho = b * rho_eos;
    if (!(std::isfinite(brho) && brho < 1.0)) {
        throw std::runtime_error("VTPR: invalid packing (b*rho >= 1)");
    }

    const double one_minus_brho = 1.0 - brho;
    const double num = 1.0 + (1.0 + SQRT2) * brho;
    const double den = 1.0 + (1.0 - SQRT2) * brho;
    if (!(std::isfinite(num) && std::isfinite(den) && num > 0.0 && den > 0.0)) {
        throw std::runtime_error("VTPR: invalid log arguments");
    }

    const double term1 = b / one_minus_brho;
    const double dlog = (1.0 + SQRT2) / num - (1.0 - SQRT2) / den;
    const double term2 = -(a / (2.0 * SQRT2 * R * T)) * dlog;
    const double dadrho_eos = term1 + term2;

    const double denom = 1.0 + mp.c_mix * rho;
    const double drhoeos_drho = 1.0 / (denom * denom);
    return dadrho_eos * drhoeos_drho;
}

double VTPREOS::dPdrho(double T, double rho, const std::vector<double>& x) const {
    const auto mp = wsParams(T, x, false);
    const double rho_eos = rhoEosFromRho(rho, mp.c_mix);

    const double v = 1.0 / rho_eos;
    const double b = mp.b;
    const double a = mp.a;
    const double R = Constants::GAS_CONSTANT;

    if (v <= b) return std::numeric_limits<double>::infinity();

    const double denom = v * v + 2.0 * b * v - b * b;
    const double dPdv = -R * T / ((v - b) * (v - b)) + 2.0 * a * (v + b) / (denom * denom);
    const double dvdrho_eos = -1.0 / (rho_eos * rho_eos);
    const double dPdrho_eos = dPdv * dvdrho_eos;

    const double denom_r = 1.0 + mp.c_mix * rho;
    const double drhoeos_drho = 1.0 / (denom_r * denom_r);
    return dPdrho_eos * drhoeos_drho;
}

EOSResult VTPREOS::calculate(
    double T,
    double rho,
    const std::vector<double>& x,
    const DerivativeSpec& deriv_spec) const
{
    EOSResult result;
    result.ln_phi.assign(static_cast<std::size_t>(nc_), 0.0);
    result.phi.assign(static_cast<std::size_t>(nc_), 0.0);

    try {
        const auto mp = wsParams(T, x, true);
        const double rho_eos = rhoEosFromRho(rho, mp.c_mix);

        double Z_eos = 0.0;
        result.ln_phi = lnPhi(T, rho_eos, x, Z_eos);
        for (int i = 0; i < nc_; ++i) {
            result.phi[i] = std::exp(result.ln_phi[i]);
        }

        result.pressure = pressureFromMix(T, 1.0 / rho_eos, mp.a, mp.b);
        result.compressibility = result.pressure / (rho * Constants::GAS_CONSTANT * T);
        result.a_residual = residualHelmholtz(T, rho, x);

        if (deriv_spec.composition) {
            result.chemical_potential = result.ln_phi;
        }
        if (deriv_spec.temperature) {
            result.da_dT = dadt(T, rho, x);
        }
        if (deriv_spec.density) {
            result.da_drho = dadrho(T, rho, x);
        }

        result.success = true;
        return result;
    } catch (const std::exception& e) {
        result.success = false;
        result.error_message = e.what();
        return result;
    }
}

std::vector<double> VTPREOS::densityRootsTP(double T, double P, const std::vector<double>& x) const {
    if (!(std::isfinite(T) && std::isfinite(P) && T > 0.0 && P > 0.0)) {
        throw std::invalid_argument("VTPREOS::densityRootsTP: invalid T/P");
    }

    const auto mp = wsParams(T, x, true);
    const double R = Constants::GAS_CONSTANT;
    const double A = mp.a * P / (R * R * T * T);
    const double B = mp.b * P / (R * T);

    // PR cubic in Z:
    //   Z^3 - (1-B)Z^2 + (A - 2B - 3B^2)Z - (AB - B^2 - B^3) = 0
    const double c3 = 1.0;
    const double c2 = B - 1.0;
    const double c1 = A - 2.0 * B - 3.0 * B * B;
    const double c0 = -A * B + B * B + B * B * B;

    std::vector<double> Z = Numerics::RootFinding::realRootsCubic(c3, c2, c1, c0);
    const double eps = 1e-12 * std::max(1.0, std::abs(B));
    Z.erase(std::remove_if(Z.begin(), Z.end(), [&](double z) {
        return !(std::isfinite(z) && (z > B + eps));
    }), Z.end());
    std::sort(Z.begin(), Z.end());
    Z.erase(std::unique(Z.begin(), Z.end(), [](double a, double b) {
        return std::abs(a - b) <= 1e-12 * std::max(1.0, std::max(std::abs(a), std::abs(b)));
    }), Z.end());

    // Convert to rho_eos, then to physical rho via volume translation (v = v_eos - c_mix).
    std::vector<double> rho;
    rho.reserve(Z.size());
    for (double z : Z) {
        const double rho_eos = P / (std::max(z, 1e-300) * R * T);
        const double r = rhoFromRhoEos(rho_eos, mp.c_mix);
        if (std::isfinite(r) && r > 0.0) rho.push_back(r);
    }
    std::sort(rho.begin(), rho.end());
    rho.erase(std::unique(rho.begin(), rho.end(), [](double a, double b) {
        return std::abs(a - b) <= 1e-12 * std::max(1.0, std::max(std::abs(a), std::abs(b)));
    }), rho.end());
    return rho;
}

double VTPREOS::criticalTemperature(int i) const {
    if (i < 0 || i >= nc_) throw std::out_of_range("Component index out of range");
    return Tc_[i];
}

double VTPREOS::criticalPressure(int i) const {
    if (i < 0 || i >= nc_) throw std::out_of_range("Component index out of range");
    return Pc_[i];
}

double VTPREOS::acentricFactor(int i) const {
    if (i < 0 || i >= nc_) throw std::out_of_range("Component index out of range");
    return omega_[i];
}

} // namespace Cubic
} // namespace DMThermo

