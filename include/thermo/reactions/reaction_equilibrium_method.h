/**
 * @file reaction_equilibrium_method.h
 * @brief Solver selection methods for reaction equilibrium problems.
 */

#ifndef THERMO_REACTIONS_REACTION_EQUILIBRIUM_METHOD_H
#define THERMO_REACTIONS_REACTION_EQUILIBRIUM_METHOD_H

namespace DMThermo {
namespace Reactions {

enum class ReactionEquilibriumMethod {
    Auto,
    GibbsExtent,       // single-phase extent-space minimization (fast)
    GibbsOptimizerTP   // general TP Gibbs optimizer with reaction extents (robust; can be multiphase)
};

} // namespace Reactions
} // namespace DMThermo

#endif // THERMO_REACTIONS_REACTION_EQUILIBRIUM_METHOD_H

