/**
 * @file phase_stability_outer_loop.h
 * @brief Shared helpers for stability-driven phase addition in outer loops.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_PHASE_STABILITY_OUTER_LOOP_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_PHASE_STABILITY_OUTER_LOOP_H

#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/core/constants.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/equilibrium/optimization/phase_allocations.h"

#include <optional>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

inline Config::StabilityConfig stabilityConfigFromEquilibriumConfig(const Config::EquilibriumOptimizerConfig& cfg) {
    auto stab_cfg = Config::StabilityConfig::defaults();
    stab_cfg.num_trial_compositions = cfg.num_stability_trials;
    stab_cfg.phase_detection_space = cfg.phase_detection_space;
    stab_cfg.max_density_newton_iters = cfg.max_density_newton_iters;
    stab_cfg.density_newton_tol = cfg.density_newton_tol;
    return stab_cfg;
}

inline std::optional<std::vector<double>> findUnstableComposition(
    const Stability::IStabilityAnalyzer& stability,
    double T,
    double P,
    const std::vector<PhaseState>& phases,
    const Config::EquilibriumOptimizerConfig& cfg,
    double x_min = Constants::MIN_MOLE_FRACTION)
{
    const auto stab_cfg = stabilityConfigFromEquilibriumConfig(cfg);
    for (const auto& ph : phases) {
        const auto stres = stability.analyze(T, P, ph.x, stab_cfg);
        if (!stres.is_stable && stres.unstable_composition.has_value()) {
            return clampAndNormalize(stres.unstable_composition.value(), x_min);
        }
    }
    return std::nullopt;
}

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_PHASE_STABILITY_OUTER_LOOP_H

