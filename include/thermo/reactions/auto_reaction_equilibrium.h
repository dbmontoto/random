/**
 * @file auto_reaction_equilibrium.h
 * @brief Auto reaction equilibrium solver that switches between specialized backends.
 */

#ifndef THERMO_REACTIONS_AUTO_REACTION_EQUILIBRIUM_H
#define THERMO_REACTIONS_AUTO_REACTION_EQUILIBRIUM_H

#include "ireaction_equilibrium.h"
#include "thermo/equilibrium/optimization/reactive_gibbs_tp_optimizer.h"
#include "reaction_equilibrium_method.h"

namespace DMThermo {
namespace Reactions {

class AutoReactionEquilibriumSolver final : public IReactionEquilibriumSolver {
public:
    explicit AutoReactionEquilibriumSolver(
        EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability = nullptr);

    ReactionEquilibriumResult solveTP(
        double T,
        double P,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const std::vector<double>& mu0,
        PhaseType phase = PhaseType::Vapor,
        const ReactionEquilibriumConfig& cfg = ReactionEquilibriumConfig{}) const override;

    ReactionEquilibriumResult solveTPFromDatabanks(
        double T,
        double P,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const Data::Databanks& databanks,
        const Core::Mixture& mixture,
        PhaseType phase = PhaseType::Vapor,
        const ReactionEquilibriumConfig& cfg = ReactionEquilibriumConfig{},
        Data::Diagnostics* diag = nullptr) const override;

    std::string algorithmName() const override { return "ReactionAuto"; }

private:
    EOSPtr eos_;
    Equilibrium::Stability::StabilityAnalyzerPtr stability_;
    GibbsReactionEquilibriumSolver extent_;
    Equilibrium::Optimization::ReactiveGibbsTPOptimizer optimizer_;
};

} // namespace Reactions
} // namespace DMThermo

#endif // THERMO_REACTIONS_AUTO_REACTION_EQUILIBRIUM_H
