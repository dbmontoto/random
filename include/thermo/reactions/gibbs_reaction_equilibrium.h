/**
 * @file gibbs_reaction_equilibrium.h
 * @brief Single-phase reaction equilibrium at constant T,P via Gibbs minimization in extent space.
 */

#ifndef THERMO_REACTIONS_GIBBS_REACTION_EQUILIBRIUM_H
#define THERMO_REACTIONS_GIBBS_REACTION_EQUILIBRIUM_H

#include "thermo/eos.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/config/solver_config.h"
#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/data/databanks.h"
#include "thermo/data/diagnostics.h"
#include "thermo/core/types.h"
#include "reaction_system.h"
#include "reaction_equilibrium_method.h"

#include <optional>
#include <string>
#include <vector>

namespace DMThermo {
namespace Reactions {

struct ReactionEquilibriumConfig {
    double standard_pressure = 101325.0; // Pa
    double min_moles = 1e-30;
    double infeasibility_penalty = 1e8;
    Config::OptimizerConfig optimizer = Config::OptimizerConfig::defaults();

    // Thermochemistry policy for mu0(T) from databanks.
    // - Default: `DIPPR.GFOR` (constant mu0 from GFOR)
    // - Option:  `DIPPR.MU0_FORMATION(T)` (requires HFOR+GFOR+ICP)
    std::string mu0_model = "DIPPR.GFOR";
    double mu0_reference_temperature = 298.15;

    // Solver selection policy:
    // - GibbsExtent: fast single-phase extent-space minimization
    // - GibbsOptimizerTP: general-purpose TP Gibbs optimizer with extents (can be multiphase)
    // - Auto: use GibbsExtent unless max_phases > 1
    ReactionEquilibriumMethod method = ReactionEquilibriumMethod::Auto;
    int max_phases = 1;
    Config::EquilibriumOptimizerConfig equilibrium_optimizer = Config::EquilibriumOptimizerConfig::defaults();
};

struct ReactionEquilibriumResult {
    bool converged = false;
    int iterations = 0;
    std::string message;

    double temperature = 0.0;
    double pressure = 0.0;

    std::vector<double> n0;
    std::vector<double> n;
    std::vector<double> x;
    std::vector<double> extents;
    double objective_g_over_rt = 0.0;

    // Optional phase information (filled when the solver computes EOS-based phase state).
    std::vector<PhaseState> phases;
};

/**
 * @brief Reaction equilibrium solver (single-phase, EOS-agnostic via fugacity coefficients).
 *
 * Notes:
 * - Currently single-phase: computes equilibrium in a specified phase at fixed T,P.
 * - Standard chemical potentials must be provided (mu0_i at cfg.standard_pressure).
 * - Phase-split + reaction will be handled by the multi-phase Gibbs/Helmholtz optimizers (Epic 13).
 */
class GibbsReactionEquilibriumSolver final {
public:
    explicit GibbsReactionEquilibriumSolver(
        EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability = nullptr)
        : eos_(std::move(eos)), stability_(std::move(stability))
    {}

    ReactionEquilibriumResult solveTPFromDatabanks(
        double T,
        double P,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const Data::Databanks& databanks,
        const Core::Mixture& mixture,
        PhaseType phase = PhaseType::Vapor,
        const ReactionEquilibriumConfig& cfg = ReactionEquilibriumConfig{},
        Data::Diagnostics* diag = nullptr) const;

    ReactionEquilibriumResult solveTP(
        double T,
        double P,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const std::vector<double>& mu0, // J/mol at cfg.standard_pressure
        PhaseType phase = PhaseType::Vapor,
        const ReactionEquilibriumConfig& cfg = ReactionEquilibriumConfig{}) const;

private:
    EOSPtr eos_;
    Equilibrium::Stability::StabilityAnalyzerPtr stability_;
};

} // namespace Reactions
} // namespace DMThermo

#endif // THERMO_REACTIONS_GIBBS_REACTION_EQUILIBRIUM_H
