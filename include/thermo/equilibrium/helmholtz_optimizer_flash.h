/**
 * @file helmholtz_optimizer_flash.h
 * @brief Flash solver wrapper using the general-purpose Helmholtz TV optimizer.
 */

#ifndef THERMO_EQUILIBRIUM_HELMHOLTZ_OPTIMIZER_FLASH_H
#define THERMO_EQUILIBRIUM_HELMHOLTZ_OPTIMIZER_FLASH_H

#include "iflash.h"
#include "optimization/helmholtz_tv_optimizer.h"

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

class HelmholtzOptimizerFlashSolver final : public IFlashSolver {
public:
    explicit HelmholtzOptimizerFlashSolver(
        EOSPtr eos,
        Stability::StabilityAnalyzerPtr stability = nullptr);

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

    std::string algorithmName() const override { return "HelmholtzOptimizer"; }

    EOSPtr eos() const override { return optimizer_.eos(); }

private:
    Optimization::HelmholtzTVOptimizer optimizer_;
};

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_HELMHOLTZ_OPTIMIZER_FLASH_H

