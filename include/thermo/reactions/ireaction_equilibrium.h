/**
 * @file ireaction_equilibrium.h
 * @brief Interface for reaction equilibrium solvers.
 */

#ifndef THERMO_REACTIONS_IREACTION_EQUILIBRIUM_H
#define THERMO_REACTIONS_IREACTION_EQUILIBRIUM_H

#include "thermo/data/databanks.h"
#include "thermo/data/diagnostics.h"
#include "thermo/core/mixture.h"
#include "gibbs_reaction_equilibrium.h"
#include "reaction_system.h"

#include <memory>
#include <string>
#include <vector>

namespace DMThermo {
namespace Reactions {

class IReactionEquilibriumSolver {
public:
    virtual ~IReactionEquilibriumSolver() = default;

    virtual ReactionEquilibriumResult solveTP(
        double T,
        double P,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const std::vector<double>& mu0, // J/mol at cfg.standard_pressure
        PhaseType phase = PhaseType::Vapor,
        const ReactionEquilibriumConfig& cfg = ReactionEquilibriumConfig{}) const = 0;

    virtual ReactionEquilibriumResult solveTPFromDatabanks(
        double T,
        double P,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const Data::Databanks& databanks,
        const Core::Mixture& mixture,
        PhaseType phase = PhaseType::Vapor,
        const ReactionEquilibriumConfig& cfg = ReactionEquilibriumConfig{},
        Data::Diagnostics* diag = nullptr) const = 0;

    virtual std::string algorithmName() const = 0;
};

using ReactionEquilibriumSolverPtr = std::shared_ptr<IReactionEquilibriumSolver>;

} // namespace Reactions
} // namespace DMThermo

#endif // THERMO_REACTIONS_IREACTION_EQUILIBRIUM_H
