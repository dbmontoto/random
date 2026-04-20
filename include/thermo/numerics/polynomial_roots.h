/**
 * @file polynomial_roots.h
 * @brief Analytic real roots for low-order polynomials.
 *
 * This header is intentionally small and dependency-free (beyond the STL) so it can
 * be reused by EOS implementations and equilibrium utilities.
 */

#ifndef THERMO_NUMERICS_POLYNOMIAL_ROOTS_H
#define THERMO_NUMERICS_POLYNOMIAL_ROOTS_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace DMThermo {
namespace Numerics {
namespace RootFinding {

namespace Detail {

inline bool isNearlyZero(double v, double scale) {
    const double s = std::max(1.0, scale);
    return std::abs(v) <= 1e-14 * s;
}

inline void sortAndDedup(std::vector<double>& roots) {
    roots.erase(std::remove_if(roots.begin(), roots.end(), [](double r) { return !std::isfinite(r); }), roots.end());
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end(), [](double a, double b) {
        return std::abs(a - b) <= 1e-12 * std::max(1.0, std::max(std::abs(a), std::abs(b)));
    }), roots.end());
}

inline std::vector<double> realRootsLinear(double a1, double a0) {
    const double scale = std::max(std::abs(a1), std::abs(a0));
    if (isNearlyZero(a1, scale)) {
        return {};
    }
    return {-a0 / a1};
}

inline std::vector<double> realRootsQuadratic(double a2, double a1, double a0) {
    const double scale = std::max({std::abs(a2), std::abs(a1), std::abs(a0)});
    if (isNearlyZero(a2, scale)) {
        return realRootsLinear(a1, a0);
    }

    const double disc = a1 * a1 - 4.0 * a2 * a0;
    const double tol = 1e-14 * std::max(1.0, std::abs(a1 * a1) + std::abs(4.0 * a2 * a0));
    if (disc < -tol) {
        return {};
    }

    const double sd = std::sqrt(std::max(0.0, disc));
    // Stable quadratic formula.
    const double q = -0.5 * (a1 + std::copysign(sd, a1));
    std::vector<double> roots;
    if (!isNearlyZero(q, scale)) {
        roots.push_back(q / a2);
        roots.push_back(a0 / q);
    } else {
        roots.push_back(-a1 / (2.0 * a2));
    }
    sortAndDedup(roots);
    return roots;
}

} // namespace Detail

/**
 * @brief Real roots of a general cubic a3*x^3 + a2*x^2 + a1*x + a0 = 0.
 *
 * Returns unique real roots (sorted ascending). Degenerate cases (a3 ~ 0) fall
 * back to quadratic/linear handling.
 */
inline std::vector<double> realRootsCubic(double a3, double a2, double a1, double a0) {
    const double scale = std::max({std::abs(a3), std::abs(a2), std::abs(a1), std::abs(a0)});
    if (Detail::isNearlyZero(a3, scale)) {
        return Detail::realRootsQuadratic(a2, a1, a0);
    }

    // Normalize: x^3 + b*x^2 + c*x + d = 0.
    const double b = a2 / a3;
    const double c = a1 / a3;
    const double d = a0 / a3;

    // Depressed cubic: t^3 + p*t + q = 0 with x = t - b/3.
    const double one_third = 1.0 / 3.0;
    const double b_over_3 = b * one_third;
    const double p = c - (b * b) * one_third;
    const double q = (2.0 * b * b * b) / 27.0 - (b * c) / 3.0 + d;

    const double disc = (q * q) / 4.0 + (p * p * p) / 27.0;
    const double disc_tol = 1e-14 * std::max(1.0, std::abs((q * q) / 4.0) + std::abs((p * p * p) / 27.0));

    std::vector<double> roots;
    roots.reserve(3);

    if (disc > disc_tol) {
        // One real root.
        const double sd = std::sqrt(disc);
        const double u = std::cbrt(-q / 2.0 + sd);
        const double v = std::cbrt(-q / 2.0 - sd);
        const double t = u + v;
        roots.push_back(t - b_over_3);
    } else if (disc < -disc_tol) {
        // Three distinct real roots.
        if (Detail::isNearlyZero(p, std::max({1.0, std::abs(p), std::abs(q)}))) {
            // p ~ 0 => t^3 + q ~ 0.
            roots.push_back(std::cbrt(-q) - b_over_3);
        } else {
            const double r = 2.0 * std::sqrt(std::max(0.0, -p / 3.0));
            const double denom = std::sqrt(std::max(0.0, -(p * p * p) / 27.0));
            double cosphi = 0.0;
            if (denom > 0.0) {
                cosphi = std::clamp((-q / 2.0) / denom, -1.0, 1.0);
            }
            const double phi = std::acos(cosphi);
            constexpr double two_pi = 6.2831853071795864769;

            for (int k = 0; k < 3; ++k) {
                const double t = r * std::cos((phi + two_pi * k) / 3.0);
                roots.push_back(t - b_over_3);
            }
        }
    } else {
        // Multiple real roots (disc ~ 0).
        if (Detail::isNearlyZero(q, std::max({1.0, std::abs(p), std::abs(q)})) &&
            Detail::isNearlyZero(p, std::max({1.0, std::abs(p), std::abs(q)}))) {
            roots.push_back(-b_over_3);
        } else {
            const double u = std::cbrt(-q / 2.0);
            roots.push_back(2.0 * u - b_over_3);
            roots.push_back(-u - b_over_3);
        }
    }

    Detail::sortAndDedup(roots);
    return roots;
}

} // namespace RootFinding
} // namespace Numerics
} // namespace DMThermo

#endif // THERMO_NUMERICS_POLYNOMIAL_ROOTS_H

