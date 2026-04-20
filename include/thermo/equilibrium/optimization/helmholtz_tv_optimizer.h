/**
 * @file helmholtz_tv_optimizer.h
 * @brief Skeleton for EOS-agnostic Helmholtz (T,V) equilibrium optimizer.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_HELMHOLTZ_TV_OPTIMIZER_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_HELMHOLTZ_TV_OPTIMIZER_H

#include "i_equilibrium_optimizer.h"

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

class HelmholtzTVOptimizer final : public IHelmholtzOptimizerTV {
public:
    explicit HelmholtzTVOptimizer(
        EOSPtr eos,
        Stability::StabilityAnalyzerPtr stability = nullptr);

    EquilibriumState minimizeTV(
        double T,
        double V,
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

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_HELMHOLTZ_TV_OPTIMIZER_H
