/**
 * @file gibbs_optimizer_flash.h
 * @brief Flash solver wrapper using the general-purpose Gibbs TP optimizer.
 */

#ifndef THERMO_EQUILIBRIUM_GIBBS_OPTIMIZER_FLASH_H
#define THERMO_EQUILIBRIUM_GIBBS_OPTIMIZER_FLASH_H

#include "iflash.h"
#include "optimization/gibbs_tp_optimizer.h"

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

class GibbsOptimizerFlashSolver final : public IFlashSolver {
public:
    explicit GibbsOptimizerFlashSolver(
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

    std::string algorithmName() const override { return "GibbsOptimizer"; }

    EOSPtr eos() const override { return optimizer_.eos(); }

private:
    Optimization::GibbsTPOptimizer optimizer_;
};

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_GIBBS_OPTIMIZER_FLASH_H

