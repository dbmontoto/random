/**
 * @file equilibrium_optimizer_config.h
 * @brief Configuration for general-purpose equilibrium optimizers.
 */

#ifndef THERMO_CONFIG_EQUILIBRIUM_OPTIMIZER_CONFIG_H
#define THERMO_CONFIG_EQUILIBRIUM_OPTIMIZER_CONFIG_H

#include "phase_detection.h"
#include "solver_config.h"
#include "tv_solve_strategy.h"

namespace DMThermo {
namespace Config {

struct EquilibriumOptimizerConfig {
    int max_phases = 3;

    PhaseDetection phase_detection_space = PhaseDetection::TP;
    int max_density_newton_iters = 40;
    double density_newton_tol = 1e-10;

    // TV solve strategy:
    // - BrentOuter: outer root in log(P) + inner Gibbs(T,P) equilibrium.
    // - DirectOptimizer: direct TV minimization using penalties (volume + pressure).
    TVSolveStrategy tv_solve_strategy = TVSolveStrategy::BrentOuter;

    // If true, uses TPD to add/remove phases in an outer loop and calls the optimizer
    // to refine a feasible equilibrium.
    bool use_phase_stability_outer_loop = true;
    int max_outer_iterations = 6;
    int num_stability_trials = 7;

    // Internal optimizer settings (unconstrained variables after transforms).
    OptimizerConfig optimizer = OptimizerConfig::defaults();

    // Feasibility / regularization:
    double beta_min = 1e-12;
    double composition_min = 1e-14;
    double infeasibility_penalty = 1e6;
    double duplicate_phase_penalty = 0.0; // 0 disables; non-zero can discourage identical phases.

    // Helmholtz (TV/TRho) specific: penalties guiding constrained equilibrium.
    double volume_penalty = 1e8;   // Enforces sum_k beta_k/rho_k = V [m3/mol]
    double pressure_penalty = 1e4; // Encourages mechanical equilibrium P_k = P_common in multiphase

    // Helmholtz (TV): outer solve configuration (P bracketing in log-space).
    int tv_bracket_steps = 50;
    double tv_bracket_logp_step = 0.5; // in decades (log10(P) step); 0.5 = half-decade
    double tv_p_min = 1e-3;            // Pa
    double tv_p_max = 1e12;            // Pa

    static EquilibriumOptimizerConfig defaults() { return EquilibriumOptimizerConfig{}; }
};

} // namespace Config
} // namespace DMThermo

#endif // THERMO_CONFIG_EQUILIBRIUM_OPTIMIZER_CONFIG_H
