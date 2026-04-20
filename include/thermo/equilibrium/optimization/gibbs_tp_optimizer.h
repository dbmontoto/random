/**
 * @file gibbs_tp_optimizer.h
 * @brief EOS-agnostic multi-phase Gibbs minimizer at constant T,P.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_GIBBS_TP_OPTIMIZER_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_GIBBS_TP_OPTIMIZER_H

#include "i_equilibrium_optimizer.h"

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

class GibbsTPOptimizer final : public IGibbsOptimizerTP {
public:
    explicit GibbsTPOptimizer(
        EOSPtr eos,
        Stability::StabilityAnalyzerPtr stability = nullptr);

    EquilibriumState minimizeTP(
        double T,
        double P,
        const std::vector<double>& z,
        const Config::EquilibriumOptimizerConfig& cfg = Config::EquilibriumOptimizerConfig::defaults()
    ) const override;

    EOSPtr eos() const override { return eos_; }

private:
    EOSPtr eos_;
    Stability::StabilityAnalyzerPtr stability_;
};

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_GIBBS_TP_OPTIMIZER_H

