/**
 * @file helmholtz_optimizer_flash.cpp
 * @brief Implementation of HelmholtzOptimizerFlashSolver.
 */

#include "thermo/equilibrium/helmholtz_optimizer_flash.h"
#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/equilibrium/internal/flash_config_support.h"
#include "thermo/equilibrium/internal/flash_result_support.h"
#include "thermo/equilibrium/internal/k_values.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

HelmholtzOptimizerFlashSolver::HelmholtzOptimizerFlashSolver(
    EOSPtr eos,
    Stability::StabilityAnalyzerPtr stability)
    : optimizer_(std::move(eos), std::move(stability))
{}

FlashResult HelmholtzOptimizerFlashSolver::solvePT(
    double,
    double,
    const std::vector<double>&,
    const Config::FlashConfig&) const
{
    throw std::runtime_error("HelmholtzOptimizerFlashSolver: solvePT not implemented (use Gibbs optimizer)");
}

FlashResult HelmholtzOptimizerFlashSolver::solveTV(
    double T,
    double V,
    const std::vector<double>& z,
    const Config::TVFlashConfig& config) const
{
    auto ocfg = Internal::equilibriumOptimizerConfigFrom(config);

    auto st = optimizer_.minimizeTV(T, V, z, ocfg);
    auto out = Internal::flashResultFromState(st);
    out.stability_test_performed = config.perform_stability_test;
    out.feed_was_stable = (out.num_phases == 1);
    return out;
}

FlashResult HelmholtzOptimizerFlashSolver::solvePH(
    double,
    double,
    const std::vector<double>&,
    const Config::PHFlashConfig&) const
{
    throw std::runtime_error("HelmholtzOptimizerFlashSolver: solvePH not implemented");
}

FlashResult HelmholtzOptimizerFlashSolver::solvePS(
    double,
    double,
    const std::vector<double>&,
    const Config::PSFlashConfig&) const
{
    throw std::runtime_error("HelmholtzOptimizerFlashSolver: solvePS not implemented");
}

FlashResult HelmholtzOptimizerFlashSolver::solveVLLE(
    double,
    double,
    const std::vector<double>&,
    const Config::FlashConfig&) const
{
    throw std::runtime_error("HelmholtzOptimizerFlashSolver: solveVLLE is only supported for TV via Helmholtz minimization");
}

std::vector<double> HelmholtzOptimizerFlashSolver::estimateKValues(double T, double P) const {
    return Internal::wilsonK(*optimizer_.eos(), T, P);
}

bool HelmholtzOptimizerFlashSolver::isValidResult(const FlashResult& result) const {
    if (!result.converged) return false;
    if (result.temperature <= 0.0) return false;
    if (result.num_phases < 1 || result.num_phases > 3) return false;
    if (result.phases.empty()) return false;
    return true;
}

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo
