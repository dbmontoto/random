/**
 * @file k_values.h
 * @brief Shared K-value estimation helpers (internal).
 */

#ifndef THERMO_EQUILIBRIUM_INTERNAL_K_VALUES_H
#define THERMO_EQUILIBRIUM_INTERNAL_K_VALUES_H

#include "thermo/eos.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Internal {

inline std::vector<double> wilsonK(const EOS& eos, double T, double P) {
    const int nc = eos.numComponents();
    std::vector<double> K(static_cast<size_t>(nc), 1.0);
    if (!(std::isfinite(T) && std::isfinite(P) && T > 0.0 && P > 0.0)) {
        return K;
    }

    for (int i = 0; i < nc; ++i) {
        const double Pc = eos.criticalPressure(i);
        const double Tc = eos.criticalTemperature(i);
        const double omega = eos.acentricFactor(i);
        if (!(std::isfinite(Pc) && std::isfinite(Tc) && std::isfinite(omega)) || Pc <= 0.0 || Tc <= 0.0) {
            K[static_cast<size_t>(i)] = 1.0;
            continue;
        }
        const double exponent = 5.37 * (1.0 + omega) * (1.0 - Tc / T);
        double v = (Pc / P) * std::exp(exponent);
        if (!std::isfinite(v) || v <= 0.0) v = 1.0;
        K[static_cast<size_t>(i)] = std::clamp(v, 1e-12, 1e12);
    }
    return K;
}

} // namespace Internal
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_INTERNAL_K_VALUES_H

