/**
 * @file tv_support.h
 * @brief Shared helpers for TV outer solves (density continuity + phase post-processing).
 */

#ifndef THERMO_EQUILIBRIUM_INTERNAL_TV_SUPPORT_H
#define THERMO_EQUILIBRIUM_INTERNAL_TV_SUPPORT_H

#include "thermo/eos.h"
#include "thermo/equilibrium/density_selection.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Internal {

template <typename RhoGuessFn, typename PhaseRankFn, typename VolumeFn>
inline bool postProcessTVEvalAtTP(
    std::vector<PhaseState>& phases,
    const EOS& eos,
    const Stability::IStabilityAnalyzer& stability,
    double T,
    double P,
    Config::PhaseDetection phase_detection_space,
    int max_density_newton_iters,
    const Core::Units::PropertyVariable& pressure_tolerance,
    double beta_min,
    RhoGuessFn rhoGuess,
    PhaseRankFn phaseRank,
    VolumeFn volumeFromPhases,
    double& Vcalc_out)
{
    for (auto& ph : phases) {
        const double rho_guess = rhoGuess(ph.type, P, ph.density);
        const auto ev = Detail::evalAtTP(
            eos,
            stability,
            T,
            P,
            ph.x,
            ph.type,
            phase_detection_space,
            rho_guess,
            max_density_newton_iters,
            pressure_tolerance);
        ph.type = ev.phase;
        ph.density = ev.rho;
        ph.compressibility = eos.compressibility(T, ph.density, ph.x);
        ph.phi = ev.phi;
    }

    double bsum = 0.0;
    for (const auto& ph : phases) bsum += std::max(ph.fraction, beta_min);
    if (bsum > 0.0) {
        for (auto& ph : phases) ph.fraction = std::max(ph.fraction, beta_min) / bsum;
    }

    std::stable_sort(phases.begin(), phases.end(), [&](const PhaseState& a, const PhaseState& b) {
        return phaseRank(a.type) < phaseRank(b.type);
    });
    for (auto& ph : phases) {
        if (ph.type == PhaseType::Liquid) ph.type = PhaseType::Liquid1;
    }

    const double Vcalc = volumeFromPhases(phases);
    if (!(std::isfinite(Vcalc) && Vcalc > 0.0)) return false;
    Vcalc_out = Vcalc;
    return true;
}

} // namespace Internal
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_INTERNAL_TV_SUPPORT_H
