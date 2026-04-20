/**
 * @file reaction_factory.cpp
 * @brief Implementation of reaction equilibrium solver factory
 */

#include "thermo/factory/reaction_factory.h"

#include "thermo/factory/stability_factory.h"
#include "thermo/reactions/auto_reaction_equilibrium.h"
#include "thermo/reactions/gibbs_reaction_equilibrium.h"
#include "thermo/reactions/thermochemistry.h"
#include "thermo/equilibrium/optimization/reactive_gibbs_tp_optimizer.h"

#include <algorithm>
#include <stdexcept>

namespace DMThermo {
namespace Factory {

bool ReactionFactory::initialized_ = false;

std::unordered_map<std::string, ReactionFactory::ReactionCreator>& ReactionFactory::registry() {
    static std::unordered_map<std::string, ReactionCreator> registry;
    return registry;
}

void ReactionFactory::registerReactionSolver(const std::string& name, ReactionCreator creator) {
    registry()[name] = std::move(creator);
}

void ReactionFactory::initializeBuiltinSolvers() {
    if (initialized_) return;

    registerReactionSolver("Auto", [](EOSPtr eos, Equilibrium::Stability::StabilityAnalyzerPtr stability) {
        return std::make_shared<Reactions::AutoReactionEquilibriumSolver>(std::move(eos), std::move(stability));
    });

    struct ExtentOnly final : public Reactions::IReactionEquilibriumSolver {
        Reactions::GibbsReactionEquilibriumSolver solver;
        explicit ExtentOnly(EOSPtr eos, Equilibrium::Stability::StabilityAnalyzerPtr stability)
            : solver(std::move(eos), std::move(stability)) {}

        Reactions::ReactionEquilibriumResult solveTP(
            double T, double P, const std::vector<double>& n0, const Reactions::ReactionSystem& system,
            const std::vector<double>& mu0, PhaseType phase, const Reactions::ReactionEquilibriumConfig& cfg) const override
        {
            return solver.solveTP(T, P, n0, system, mu0, phase, cfg);
        }

        Reactions::ReactionEquilibriumResult solveTPFromDatabanks(
            double T, double P, const std::vector<double>& n0, const Reactions::ReactionSystem& system,
            const Data::Databanks& databanks, const Core::Mixture& mixture,
            PhaseType phase, const Reactions::ReactionEquilibriumConfig& cfg, Data::Diagnostics* diag) const override
        {
            return solver.solveTPFromDatabanks(T, P, n0, system, databanks, mixture, phase, cfg, diag);
        }

        std::string algorithmName() const override { return "GibbsExtent"; }
    };

    registerReactionSolver("GibbsExtent", [](EOSPtr eos, Equilibrium::Stability::StabilityAnalyzerPtr stability) {
        return std::make_shared<ExtentOnly>(std::move(eos), std::move(stability));
    });

    struct OptimizerTP final : public Reactions::IReactionEquilibriumSolver {
        EOSPtr eos;
        Equilibrium::Stability::StabilityAnalyzerPtr stability;
        Equilibrium::Optimization::ReactiveGibbsTPOptimizer optimizer;

        explicit OptimizerTP(EOSPtr eos_in, Equilibrium::Stability::StabilityAnalyzerPtr stability_in)
            : eos(std::move(eos_in))
            , stability(std::move(stability_in))
            , optimizer(eos, stability)
        {}

        Reactions::ReactionEquilibriumResult solveTP(
            double T, double P, const std::vector<double>& n0, const Reactions::ReactionSystem& system,
            const std::vector<double>& mu0, PhaseType, const Reactions::ReactionEquilibriumConfig& cfg) const override
        {
            Equilibrium::Optimization::ReactiveGibbsTPConfig rcfg;
            rcfg.standard_pressure = cfg.standard_pressure;
            rcfg.min_moles = cfg.min_moles;
            rcfg.equilibrium = cfg.equilibrium_optimizer;
            rcfg.equilibrium.max_phases = std::clamp(cfg.max_phases, 1, 3);

            const auto st = optimizer.minimizeTPReactive(T, P, n0, system, mu0, rcfg);
            if (!st.reaction_extents.has_value()) {
                throw std::runtime_error("GibbsOptimizerTP: missing extents in optimizer result");
            }

            Reactions::ReactionEquilibriumResult out;
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

        Reactions::ReactionEquilibriumResult solveTPFromDatabanks(
            double T, double P, const std::vector<double>& n0, const Reactions::ReactionSystem& system,
            const Data::Databanks& databanks, const Core::Mixture& mixture,
            PhaseType phase, const Reactions::ReactionEquilibriumConfig& cfg, Data::Diagnostics* diag) const override
        {
            const std::vector<double> mu0 = Reactions::Thermochemistry::mu0DB(
                cfg.mu0_model, databanks, mixture, T, cfg.mu0_reference_temperature, diag);
            return solveTP(T, P, n0, system, mu0, phase, cfg);
        }

        std::string algorithmName() const override { return "GibbsOptimizerTP"; }
    };

    registerReactionSolver("GibbsOptimizerTP", [](EOSPtr eos, Equilibrium::Stability::StabilityAnalyzerPtr stability) {
        if (!stability) stability = StabilityFactory::create(eos);
        return std::make_shared<OptimizerTP>(std::move(eos), std::move(stability));
    });

    initialized_ = true;
}

Reactions::ReactionEquilibriumSolverPtr ReactionFactory::create(
    const std::string& method,
    EOSPtr eos,
    Equilibrium::Stability::StabilityAnalyzerPtr stability)
{
    initializeBuiltinSolvers();
    if (!stability) {
        stability = StabilityFactory::create(eos);
    }

    auto it = registry().find(method);
    if (it == registry().end()) {
        throw std::invalid_argument("Unknown reaction method: " + method);
    }
    return it->second(std::move(eos), std::move(stability));
}

Reactions::ReactionEquilibriumSolverPtr ReactionFactory::create(EOSPtr eos) {
    return create("Auto", std::move(eos), nullptr);
}

std::vector<std::string> ReactionFactory::registeredTypes() {
    initializeBuiltinSolvers();
    std::vector<std::string> types;
    for (const auto& pair : registry()) types.push_back(pair.first);
    std::sort(types.begin(), types.end());
    return types;
}

} // namespace Factory
} // namespace DMThermo
