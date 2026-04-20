/**
 * @file auto_flash.h
 * @brief Auto flash solver that routes to explicit backends (optimizers vs legacy).
 */

#ifndef THERMO_EQUILIBRIUM_AUTO_FLASH_H
#define THERMO_EQUILIBRIUM_AUTO_FLASH_H

#include "iflash.h"
#include "gibbs_multi_flash.h"
#include "gibbs_optimizer_flash.h"
#include "helmholtz_optimizer_flash.h"

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

/**
 * @brief Auto flash solver with explicit policy (no hidden defaults).
 *
 * Policy:
 * - PT/VLLE: uses GibbsOptimizer (general-purpose minimization backend).
 * - TV: always uses HelmholtzOptimizer.
 */
class AutoFlashSolver final : public IFlashSolver {
public:
    explicit AutoFlashSolver(
        EOSPtr eos,
        Stability::StabilityAnalyzerPtr stability = nullptr
    );

    FlashResult solvePT(
        double T,
        double P,
        const std::vector<double>& z,
        const Config::FlashConfig& config = Config::FlashConfig::defaults()
    ) const override;

    FlashResult solveTV(
        double T,
        double V,
        const std::vector<double>& z,
        const Config::TVFlashConfig& config = Config::TVFlashConfig{}
    ) const override;

    FlashResult solvePH(
        double P,
        double H,
        const std::vector<double>& z,
        const Config::PHFlashConfig& config = Config::PHFlashConfig{}
    ) const override;

    FlashResult solvePS(
        double P,
        double S,
        const std::vector<double>& z,
        const Config::PSFlashConfig& config = Config::PSFlashConfig{}
    ) const override;

    FlashResult solveVLLE(
        double T,
        double P,
        const std::vector<double>& z,
        const Config::FlashConfig& config = Config::FlashConfig::threePhase()
    ) const override;

    std::vector<double> estimateKValues(double T, double P) const override;
    bool isValidResult(const FlashResult& result) const override;

    std::string algorithmName() const override { return "Auto"; }
    EOSPtr eos() const override { return eos_; }

private:
    EOSPtr eos_;
    Stability::StabilityAnalyzerPtr stability_;

    GibbsMultiPhaseFlashSolver legacy_multi_;
    GibbsOptimizerFlashSolver gibbs_;
    HelmholtzOptimizerFlashSolver helmholtz_;
};

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_AUTO_FLASH_H
