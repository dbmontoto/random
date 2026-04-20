/**
 * @file rachford_rice.h
 * @brief Shared Rachford-Rice utilities (internal).
 */

#ifndef THERMO_EQUILIBRIUM_INTERNAL_RACHFORD_RICE_H
#define THERMO_EQUILIBRIUM_INTERNAL_RACHFORD_RICE_H

#include <cmath>
#include <limits>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Internal {

// Solves the standard 2-phase Rachford–Rice equation for vapor fraction beta in [0,1].
// Returns true when a bracket is found and bisection completes; false otherwise.
inline bool solveRachfordRiceBeta(
    const std::vector<double>& z,
    const std::vector<double>& K,
    double& beta_out,
    int max_iters = 80,
    double tol = 1e-12)
{
    if (z.size() != K.size() || z.empty()) {
        return false;
    }
    const int nc = static_cast<int>(z.size());

    auto f = [&](double beta) -> double {
        double s = 0.0;
        for (int i = 0; i < nc; ++i) {
            const double denom = 1.0 + beta * (K[static_cast<size_t>(i)] - 1.0);
            s += z[static_cast<size_t>(i)] * (K[static_cast<size_t>(i)] - 1.0) / denom;
        }
        return s;
    };

    double a = 0.0, b = 1.0;
    double fa = f(a), fb = f(b);
    if (!(std::isfinite(fa) && std::isfinite(fb)) || fa * fb > 0.0) {
        return false;
    }

    for (int it = 0; it < max_iters; ++it) {
        const double m = 0.5 * (a + b);
        const double fm = f(m);
        if (!std::isfinite(fm)) return false;
        if (std::abs(fm) < tol) {
            beta_out = m;
            return true;
        }
        if (fa * fm < 0.0) { b = m; fb = fm; } else { a = m; fa = fm; }
    }

    beta_out = 0.5 * (a + b);
    return true;
}

} // namespace Internal
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_INTERNAL_RACHFORD_RICE_H

