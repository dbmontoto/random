/**
 * @file equilibrium_state.h
 * @brief Shared equilibrium state for optimizers (flash + reactions).
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_EQUILIBRIUM_STATE_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_EQUILIBRIUM_STATE_H

#include <optional>
#include <string>
#include <vector>
#include "thermo/core/types.h"

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

struct EquilibriumState {
    bool converged = false;
    int iterations = 0;
    double objective_value = 0.0;
    double infeasibility = 0.0;
    double score_value = 0.0;
    double temperature = 0.0;
    double pressure = 0.0;

    std::vector<double> z;
    std::vector<PhaseState> phases;

    // Reaction extension point: ξ_r (optional, not yet implemented).
    std::optional<std::vector<double>> reaction_extents;

    std::string message;
    std::string method_used;
};

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_EQUILIBRIUM_STATE_H
