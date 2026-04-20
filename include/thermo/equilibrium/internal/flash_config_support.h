/**
 * @file flash_config_support.h
 * @brief Shared helpers to map flash configs to optimizer configs (internal).
 */

#ifndef THERMO_EQUILIBRIUM_INTERNAL_FLASH_CONFIG_SUPPORT_H
#define THERMO_EQUILIBRIUM_INTERNAL_FLASH_CONFIG_SUPPORT_H

#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/config/flash_config.h"

#include <algorithm>

namespace DMThermo {
namespace Equilibrium {
namespace Internal {

inline Config::EquilibriumOptimizerConfig equilibriumOptimizerConfigFrom(const Config::FlashConfig& config) {
    Config::EquilibriumOptimizerConfig ocfg = Config::EquilibriumOptimizerConfig::defaults();
    ocfg.max_phases = std::clamp(config.max_phases, 1, 3);
    ocfg.phase_detection_space = config.phase_detection_space;
    ocfg.max_density_newton_iters = config.max_density_newton_iters;
    ocfg.density_newton_tol = config.density_newton_tol;
    ocfg.use_phase_stability_outer_loop = config.perform_stability_test;
    ocfg.num_stability_trials = config.num_stability_trials;
    // Flash tolerances are typically much tighter than what a finite-difference BFGS
    // can reliably satisfy in gradient-norm terms. Keep function tolerance aligned
    // with the flash tolerance, but clamp the gradient-norm tolerance to a practical floor.
    ocfg.optimizer.tolerance = std::max(config.tolerance, 1e-6);
    ocfg.optimizer.function_tolerance = config.tolerance;
    ocfg.optimizer.max_iterations = config.max_iterations;
    return ocfg;
}

inline Config::EquilibriumOptimizerConfig equilibriumOptimizerConfigFrom(const Config::TVFlashConfig& config) {
    auto ocfg = equilibriumOptimizerConfigFrom(static_cast<const Config::FlashConfig&>(config));
    ocfg.tv_solve_strategy = config.tv_solve_strategy;
    ocfg.tv_bracket_steps = config.tv_bracket_steps;
    ocfg.tv_bracket_logp_step = config.tv_bracket_logp_step;
    ocfg.tv_p_min = config.tv_p_min;
    ocfg.tv_p_max = config.tv_p_max;
    return ocfg;
}

} // namespace Internal
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_INTERNAL_FLASH_CONFIG_SUPPORT_H
