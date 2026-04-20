/**
 * @file pure_saturation.cpp
 * @brief Implementation of generic pure-component saturation solver
 */

#include "thermo/equilibrium/pure_saturation.h"
#include "thermo/core/constants.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Saturation {

PureSaturationSolver::PureSaturationSolver(EOSPtr eos, Stability::StabilityAnalyzerPtr stability)
    : eos_(std::move(eos))
    , stability_(std::move(stability))
{
    if (!eos_) {
        throw std::invalid_argument("PureSaturationSolver: EOS pointer is null");
    }
    if (!stability_) {
        throw std::invalid_argument("PureSaturationSolver: Stability analyzer pointer is null");
    }
    if (eos_->numComponents() != 1) {
        throw std::invalid_argument("PureSaturationSolver currently supports only pure systems (nc == 1)");
    }
}

PureSaturationSolver::PhiAtTP PureSaturationSolver::phiAtTP(double T, double P) const {
    PhiAtTP out;

    const std::vector<double> x = {1.0};
    auto roots = stability_->findDensityRoots(T, P, x);

    if (roots.num_roots < 2) {
        out.valid = false;
        out.message = "No two-phase density roots";
        return out;
    }

    const double rho_v = *std::min_element(roots.densities.begin(), roots.densities.end());
    const double rho_l = *std::max_element(roots.densities.begin(), roots.densities.end());

    const auto phi_v = eos_->fugacityCoefficients(T, rho_v, x);
    const auto phi_l = eos_->fugacityCoefficients(T, rho_l, x);

    if (phi_v.empty() || phi_l.empty()) {
        out.valid = false;
        out.message = "Empty fugacity coefficients";
        return out;
    }

    const double lnphi_v = std::log(std::max(phi_v[0], 1e-300));
    const double lnphi_l = std::log(std::max(phi_l[0], 1e-300));

    if (!std::isfinite(lnphi_v) || !std::isfinite(lnphi_l)) {
        out.valid = false;
        out.message = "Non-finite ln(phi)";
        return out;
    }

    out.lnphi_v = lnphi_v;
    out.lnphi_l = lnphi_l;
    out.rho_v = rho_v;
    out.rho_l = rho_l;
    out.valid = true;
    return out;
}

PureSaturationPoint PureSaturationSolver::solveAtT(double T, const PureSaturationConfig& config) const {
    PureSaturationPoint pt;
    pt.T = T;

    if (T <= 0.0) {
        pt.converged = false;
        pt.message = "Invalid temperature";
        return pt;
    }

    const double Pc = eos_->criticalPressure(0);
    const double Tc = eos_->criticalTemperature(0);

    double Pmin = std::max(config.P_min, 1.0);
    double Pmax = config.P_max;
    if (Pmax <= 0.0) {
        Pmax = (Pc > 0.0) ? 0.999 * Pc : 1e8;
    }
    if (Pc > 0.0) {
        Pmax = std::min(Pmax, 0.999 * Pc);
    }

    if (config.P_guess > 0.0 && std::isfinite(config.P_guess) && config.guess_span_decades > 0.0) {
        const double factor = std::pow(10.0, 0.5 * config.guess_span_decades);
        const double guess_min = config.P_guess / factor;
        const double guess_max = config.P_guess * factor;
        if (std::isfinite(guess_min) && std::isfinite(guess_max) && guess_min > 0.0 && guess_max > guess_min) {
            Pmin = std::max(Pmin, guess_min);
            Pmax = std::min(Pmax, guess_max);
        }
    }

    if (Tc > 0.0 && T >= Tc) {
        pt.converged = false;
        pt.message = "T >= Tc (no saturation)";
        return pt;
    }

    if (Pmax <= Pmin) {
        pt.converged = false;
        pt.message = "Invalid pressure bounds";
        return pt;
    }

    // Scan log(P) to bracket sign change of g(P) = lnphi_l - lnphi_v.
    const int nscan = std::max(10, config.max_scan_points);
    const double logPmin = std::log(Pmin);
    const double logPmax = std::log(Pmax);

    double P_lo = 0.0, P_hi = 0.0;
    double g_lo = 0.0, g_hi = 0.0;
    bool bracket_found = false;

    double prev_P = 0.0;
    double prev_g = 0.0;
    bool prev_valid = false;

    for (int i = 0; i < nscan; ++i) {
        const double t = (nscan <= 1) ? 0.0 : static_cast<double>(i) / (nscan - 1);
        const double P = std::exp(logPmin + t * (logPmax - logPmin));

        const auto phi = phiAtTP(T, P);
        if (!phi.valid) {
            prev_valid = false;
            continue;
        }

        const double g = phi.lnphi_l - phi.lnphi_v;
        if (!std::isfinite(g)) {
            prev_valid = false;
            continue;
        }

        if (prev_valid) {
            if ((prev_g < 0.0 && g > 0.0) || (prev_g > 0.0 && g < 0.0)) {
                P_lo = prev_P;
                P_hi = P;
                g_lo = prev_g;
                g_hi = g;
                bracket_found = true;
                break;
            }
        }

        prev_P = P;
        prev_g = g;
        prev_valid = true;
    }

    if (!bracket_found) {
        pt.converged = false;
        pt.message = "Failed to bracket saturation pressure";
        return pt;
    }

    // Bisection solve.
    double a = P_lo, b = P_hi;
    double ga = g_lo, gb = g_hi;

    for (int it = 0; it < config.max_iterations; ++it) {
        const double m = 0.5 * (a + b);
        const auto phi = phiAtTP(T, m);
        if (!phi.valid) {
            // Shrink interval conservatively.
            b = m;
            gb = 0.5 * (ga + gb);
            continue;
        }

        const double gm = phi.lnphi_l - phi.lnphi_v;
        if (!std::isfinite(gm)) {
            b = m;
            gb = 0.5 * (ga + gb);
            continue;
        }

        if (std::abs(gm) < config.tolerance) {
            pt.P_sat = m;
            pt.rho_v = phi.rho_v;
            pt.rho_l = phi.rho_l;
            pt.converged = true;
            pt.message = "Converged";
            return pt;
        }

        if ((ga < 0.0 && gm > 0.0) || (ga > 0.0 && gm < 0.0)) {
            b = m;
            gb = gm;
        } else {
            a = m;
            ga = gm;
        }

        if (std::abs(b - a) / std::max(1.0, m) < config.relP_tolerance) {
            const auto phi2 = phiAtTP(T, m);
            if (phi2.valid) {
                pt.P_sat = m;
                pt.rho_v = phi2.rho_v;
                pt.rho_l = phi2.rho_l;
                pt.converged = true;
                pt.message = "Converged (pressure tolerance)";
                return pt;
            }
        }
    }

    pt.converged = false;
    pt.message = "Bisection did not converge";
    return pt;
}

std::vector<PureSaturationPoint> PureSaturationSolver::solveCurve(
    double T_min,
    double T_max,
    int n_points,
    const PureSaturationConfig& config) const
{
    std::vector<PureSaturationPoint> curve;
    if (n_points <= 0) return curve;

    curve.reserve(static_cast<size_t>(n_points));
    for (int i = 0; i < n_points; ++i) {
        const double t = (n_points <= 1) ? 0.0 : static_cast<double>(i) / (n_points - 1);
        const double T = T_min + t * (T_max - T_min);
        curve.push_back(solveAtT(T, config));
    }
    return curve;
}

} // namespace Saturation
} // namespace Equilibrium
} // namespace DMThermo
