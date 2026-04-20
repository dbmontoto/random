/**
 * @file gibbs_multi_flash.h
 * @brief EOS-agnostic multi-phase Gibbs flash solver (TPD-driven phase detection)
 */

#ifndef THERMO_EQUILIBRIUM_GIBBS_MULTI_FLASH_H
#define THERMO_EQUILIBRIUM_GIBBS_MULTI_FLASH_H

#include "iflash.h"
#include "istability.h"
#include "thermo/eos.h"

#include <memory>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

/**
 * @brief EOS-agnostic Gibbs minimization flash solver (PT, up to 3 phases)
 *
 * This solver replaces the 2-phase phi-phi method with a more general workflow:
 * - Phase detection via TPD stability tests
 * - Multi-phase split via generalized Rachford-Rice on fugacity-ratio K-values
 *
 * Notes:
 * - Works with any EOS that can compute fugacity coefficients and pressure at (T, rho, x).
 * - Uses the provided StabilityAnalyzer for density root selection and stability checks.
 */
class GibbsMultiPhaseFlashSolver final : public IFlashSolver {
public:
    explicit GibbsMultiPhaseFlashSolver(
        EOSPtr eos,
        Stability::StabilityAnalyzerPtr stability
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

    std::string algorithmName() const override { return "Gibbs-MultiPhase"; }
    EOSPtr eos() const override { return eos_; }

private:
    EOSPtr eos_;
    Stability::StabilityAnalyzerPtr stability_;
};

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_GIBBS_MULTI_FLASH_H
