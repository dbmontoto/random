/**
 * @file tv_solve_strategy.h
 * @brief Strategy selection for TV (Helmholtz) flash solving.
 */

#ifndef THERMO_CONFIG_TV_SOLVE_STRATEGY_H
#define THERMO_CONFIG_TV_SOLVE_STRATEGY_H

namespace DMThermo {
namespace Config {

/**
 * @brief Strategy for solving TV equilibrium (Helmholtz minimization).
 *
 * - BrentOuter: Solve an outer pressure root in log(P) to match V, calling TP (Gibbs) equilibrium internally.
 * - DirectOptimizer: Directly minimize a Helmholtz-based objective with explicit constraint penalties.
 */
enum class TVSolveStrategy {
    BrentOuter,
    DirectOptimizer
};

} // namespace Config
} // namespace DMThermo

#endif // THERMO_CONFIG_TV_SOLVE_STRATEGY_H

