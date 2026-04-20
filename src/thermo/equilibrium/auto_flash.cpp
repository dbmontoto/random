/**
 * @file auto_flash.cpp
 * @brief Implementation of AutoFlashSolver.
 */

#include "thermo/equilibrium/auto_flash.h"
#include "thermo/factory/stability_factory.h"

#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

namespace {

EOSPtr requireEOS(EOSPtr eos) {
    if (!eos) {
        throw std::invalid_argument("AutoFlashSolver: EOS is null");
    }
    return eos;
}

} // namespace

AutoFlashSolver::AutoFlashSolver(EOSPtr eos, Stability::StabilityAnalyzerPtr stability)
    : eos_(requireEOS(std::move(eos)))
    , stability_(stability ? std::move(stability) : Factory::StabilityFactory::create(eos_))
    , legacy_multi_(eos_, stability_)
    , gibbs_(eos_, stability_)
    , helmholtz_(eos_, stability_)
{}

FlashResult AutoFlashSolver::solvePT(
    double T,
    double P,
    const std::vector<double>& z,
    const Config::FlashConfig& config) const
{
    return gibbs_.solvePT(T, P, z, config);
}

FlashResult AutoFlashSolver::solveTV(
    double T,
    double V,
    const std::vector<double>& z,
    const Config::TVFlashConfig& config) const
{
    return helmholtz_.solveTV(T, V, z, config);
}

FlashResult AutoFlashSolver::solvePH(
    double P,
    double H,
    const std::vector<double>& z,
    const Config::PHFlashConfig& config) const
{
    // Prefer optimizer backend when ideal-gas correlations are available via EOS::coreMixture().
    if (eos_->coreMixture()) {
        return gibbs_.solvePH(P, H, z, config);
    }
    return legacy_multi_.solvePH(P, H, z, config);
}

FlashResult AutoFlashSolver::solvePS(
    double P,
    double S,
    const std::vector<double>& z,
    const Config::PSFlashConfig& config) const
{
    if (eos_->coreMixture()) {
        return gibbs_.solvePS(P, S, z, config);
    }
    return legacy_multi_.solvePS(P, S, z, config);
}

FlashResult AutoFlashSolver::solveVLLE(
    double T,
    double P,
    const std::vector<double>& z,
    const Config::FlashConfig& config) const
{
    return gibbs_.solveVLLE(T, P, z, config);
}

std::vector<double> AutoFlashSolver::estimateKValues(double T, double P) const {
    // Use the modern backend's Wilson estimate by default.
    return gibbs_.estimateKValues(T, P);
}

bool AutoFlashSolver::isValidResult(const FlashResult& result) const {
    if (!result.converged) return false;
    if (result.temperature <= 0.0) return false;
    if (result.num_phases < 1 || result.num_phases > 3) return false;
    if (result.phases.empty()) return false;
    return true;
}

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo
