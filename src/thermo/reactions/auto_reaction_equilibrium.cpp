/**
 * @file auto_reaction_equilibrium.cpp
 * @brief Implementation of AutoReactionEquilibriumSolver.
 */

#include "thermo/reactions/auto_reaction_equilibrium.h"
#include "thermo/reactions/thermochemistry.h"

#include <algorithm>
#include <stdexcept>

namespace DMThermo {
namespace Reactions {

AutoReactionEquilibriumSolver::AutoReactionEquilibriumSolver(
    EOSPtr eos,
    Equilibrium::Stability::StabilityAnalyzerPtr stability)
    : eos_(std::move(eos))
    , stability_(std::move(stability))
    , extent_(eos_, stability_)
    , optimizer_(eos_, stability_)
{
    if (!eos_) {
        throw std::invalid_argument("AutoReactionEquilibriumSolver: EOS is null");
    }
}

ReactionEquilibriumResult AutoReactionEquilibriumSolver::solveTP(
    double T,
    double P,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const std::vector<double>& mu0,
    PhaseType phase,
    const ReactionEquilibriumConfig& cfg) const
{
    if (cfg.standard_pressure <= 0.0) {
        throw std::invalid_argument("AutoReactionEquilibriumSolver::solveTP: invalid standard_pressure");
    }

    ReactionEquilibriumMethod method = cfg.method;
    if (method == ReactionEquilibriumMethod::Auto) {
        method = (cfg.max_phases > 1) ? ReactionEquilibriumMethod::GibbsOptimizerTP : ReactionEquilibriumMethod::GibbsExtent;
    }

    if (method == ReactionEquilibriumMethod::GibbsExtent) {
        return extent_.solveTP(T, P, n0, system, mu0, phase, cfg);
    }

    if (method != ReactionEquilibriumMethod::GibbsOptimizerTP) {
        throw std::invalid_argument("AutoReactionEquilibriumSolver::solveTP: unknown reaction equilibrium method");
    }

    Equilibrium::Optimization::ReactiveGibbsTPConfig rcfg;
    rcfg.standard_pressure = cfg.standard_pressure;
    rcfg.min_moles = cfg.min_moles;
    rcfg.equilibrium = cfg.equilibrium_optimizer;
    rcfg.equilibrium.max_phases = std::clamp(cfg.max_phases, 1, 3);

    const auto st = optimizer_.minimizeTPReactive(T, P, n0, system, mu0, rcfg);
    if (!st.reaction_extents.has_value()) {
        throw std::runtime_error("AutoReactionEquilibriumSolver::solveTP: optimizer backend did not return extents");
    }

    ReactionEquilibriumResult out;
    out.temperature = T;
    out.pressure = P;
    out.n0 = n0;
    out.extents = st.reaction_extents.value();
    out.converged = st.converged;
    out.iterations = st.iterations;
    out.message = st.message;
    out.objective_g_over_rt = st.objective_value;
    out.phases = st.phases;

    out.n = system.molesFromExtents(n0, out.extents);
    double N = 0.0;
    for (double ni : out.n) N += std::max(ni, cfg.min_moles);
    out.x.assign(out.n.size(), 0.0);
    if (N > 0.0 && std::isfinite(N)) {
        for (std::size_t i = 0; i < out.n.size(); ++i) out.x[i] = std::max(out.n[i], cfg.min_moles) / N;
    }

    return out;
}

ReactionEquilibriumResult AutoReactionEquilibriumSolver::solveTPFromDatabanks(
    double T,
    double P,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    PhaseType phase,
    const ReactionEquilibriumConfig& cfg,
    Data::Diagnostics* diag) const
{
    const std::vector<double> mu0 = Thermochemistry::mu0DB(
        cfg.mu0_model,
        databanks,
        mixture,
        T,
        cfg.mu0_reference_temperature,
        diag);

    return solveTP(T, P, n0, system, mu0, phase, cfg);
}

} // namespace Reactions
} // namespace DMThermo
