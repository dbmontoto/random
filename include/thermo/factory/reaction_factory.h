/**
 * @file reaction_factory.h
 * @brief Factory for creating reaction equilibrium solver objects
 */

#ifndef THERMO_FACTORY_REACTION_FACTORY_H
#define THERMO_FACTORY_REACTION_FACTORY_H

#include "thermo/eos.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/reactions/ireaction_equilibrium.h"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace DMThermo {
namespace Factory {

/**
 * @brief Factory for creating reaction equilibrium solver implementations
 *
 * Built-in types:
 * - "Auto": policy-based switching (fast extent vs robust optimizer backend)
 * - "GibbsExtent": single-phase extent-space Gibbs minimization
 * - "GibbsOptimizerTP": TP Gibbs optimizer with reaction extents (phase-capable)
 */
class ReactionFactory {
public:
    using ReactionCreator = std::function<Reactions::ReactionEquilibriumSolverPtr(
        EOSPtr,
        Equilibrium::Stability::StabilityAnalyzerPtr
    )>;

    static Reactions::ReactionEquilibriumSolverPtr create(
        const std::string& method,
        EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability = nullptr
    );

    static Reactions::ReactionEquilibriumSolverPtr create(EOSPtr eos);

    static void registerReactionSolver(const std::string& name, ReactionCreator creator);
    static std::vector<std::string> registeredTypes();

private:
    static std::unordered_map<std::string, ReactionCreator>& registry();
    static void initializeBuiltinSolvers();
    static bool initialized_;
};

} // namespace Factory
} // namespace DMThermo

#endif // THERMO_FACTORY_REACTION_FACTORY_H
