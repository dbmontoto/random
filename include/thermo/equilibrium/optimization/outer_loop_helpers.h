/**
 * @file outer_loop_helpers.h
 * @brief Glue utilities between equilibrium configs and generic numerics runners.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_OUTER_LOOP_HELPERS_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_OUTER_LOOP_HELPERS_H

#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/numerics/outer_loop_runner.h"

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

inline Numerics::Optimization::OuterLoopConfig outerLoopConfigFrom(const Config::EquilibriumOptimizerConfig& cfg) {
    Numerics::Optimization::OuterLoopConfig out;
    out.enabled = cfg.use_phase_stability_outer_loop;
    out.max_outer_iterations = cfg.max_outer_iterations;
    out.optimizer = cfg.optimizer;
    return out;
}

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_OUTER_LOOP_HELPERS_H

