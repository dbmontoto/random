/**
 * @file flash_factory.cpp
 * @brief Implementation of flash solver factory
 */

#include "thermo/factory/flash_factory.h"
#include "thermo/factory/stability_factory.h"
#include "thermo/equilibrium/auto_flash.h"
#include "thermo/equilibrium/gibbs_multi_flash.h"
#include "thermo/equilibrium/gibbs_optimizer_flash.h"
#include "thermo/equilibrium/helmholtz_optimizer_flash.h"
#include <stdexcept>
#include <algorithm>

namespace DMThermo {
namespace Factory {

bool FlashFactory::initialized_ = false;

std::unordered_map<std::string, FlashFactory::FlashCreator>& FlashFactory::registry() {
    static std::unordered_map<std::string, FlashCreator> registry;
    return registry;
}

void FlashFactory::initializeBuiltinSolvers() {
    if (initialized_) return;

    // Register Gibbs solver (EOS-agnostic multi-phase PT flash).
    registerFlashSolver("Gibbs", [](EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability) -> Equilibrium::Flash::FlashSolverPtr {
        return std::make_shared<Equilibrium::Flash::GibbsMultiPhaseFlashSolver>(eos, stability);
    });

    // Register Gibbs optimizer (general-purpose minimization backend; replaces SS phi-phi when selected).
    registerFlashSolver("GibbsOptimizer", [](EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability) -> Equilibrium::Flash::FlashSolverPtr {
        return std::make_shared<Equilibrium::Flash::GibbsOptimizerFlashSolver>(eos, stability);
    });

    // Register Helmholtz optimizer (general-purpose minimization backend for TV problems).
    registerFlashSolver("HelmholtzOptimizer", [](EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability) -> Equilibrium::Flash::FlashSolverPtr {
        return std::make_shared<Equilibrium::Flash::HelmholtzOptimizerFlashSolver>(eos, stability);
    });

    // Register Helmholtz solver (legacy alias -> optimizer backend).
    registerFlashSolver("Helmholtz", [](EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability) -> Equilibrium::Flash::FlashSolverPtr {
        return std::make_shared<Equilibrium::Flash::HelmholtzOptimizerFlashSolver>(eos, stability);
    });

    // Register Combined solver (legacy alias -> auto policy).
    registerFlashSolver("Combined", [](EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability) -> Equilibrium::Flash::FlashSolverPtr {
        return std::make_shared<Equilibrium::Flash::AutoFlashSolver>(eos, stability);
    });

    // Register Auto (defaults to Gibbs)
    registerFlashSolver("Auto", [](EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability) -> Equilibrium::Flash::FlashSolverPtr {
        // Default policy: prefer general-purpose optimizers (GibbsOptimizer / HelmholtzOptimizer).
        return std::make_shared<Equilibrium::Flash::AutoFlashSolver>(eos, stability);
    });

    initialized_ = true;
}

Equilibrium::Flash::FlashSolverPtr FlashFactory::create(
    const std::string& method,
    EOSPtr eos,
    Equilibrium::Stability::StabilityAnalyzerPtr stability)
{
    initializeBuiltinSolvers();

    // Create stability analyzer if not provided
    if (!stability) {
        stability = StabilityFactory::create(eos);
    }

    auto it = registry().find(method);
    if (it == registry().end()) {
        throw std::invalid_argument("Unknown flash method: " + method);
    }

    return it->second(eos, stability);
}

Equilibrium::Flash::FlashSolverPtr FlashFactory::create(EOSPtr eos) {
    return create("Auto", eos, nullptr);
}

Equilibrium::Flash::FlashSolverPtr FlashFactory::createGibbs(
    EOSPtr eos,
    Equilibrium::Stability::StabilityAnalyzerPtr stability)
{
    return create("Gibbs", eos, stability);
}

Equilibrium::Flash::FlashSolverPtr FlashFactory::createHelmholtz(
    EOSPtr eos,
    Equilibrium::Stability::StabilityAnalyzerPtr stability)
{
    return create("Helmholtz", eos, stability);
}

Equilibrium::Flash::FlashSolverPtr FlashFactory::createCombined(
    EOSPtr eos,
    Equilibrium::Stability::StabilityAnalyzerPtr stability)
{
    return create("Combined", eos, stability);
}

void FlashFactory::registerFlashSolver(const std::string& name, FlashCreator creator) {
    registry()[name] = creator;
}

std::vector<std::string> FlashFactory::registeredTypes() {
    initializeBuiltinSolvers();
    std::vector<std::string> types;
    for (const auto& pair : registry()) {
        types.push_back(pair.first);
    }
    std::sort(types.begin(), types.end());
    return types;
}

} // namespace Factory
} // namespace DMThermo
