/**
 * @file i_equilibrium_optimizer.h
 * @brief Interfaces for general-purpose equilibrium optimizers.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_I_EQUILIBRIUM_OPTIMIZER_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_I_EQUILIBRIUM_OPTIMIZER_H

#include <memory>
#include <vector>
#include "thermo/eos.h"
#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/equilibrium/istability.h"
#include "equilibrium_state.h"

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

class IGibbsOptimizerTP {
public:
    virtual ~IGibbsOptimizerTP() = default;

    virtual EquilibriumState minimizeTP(
        double T,
        double P,
        const std::vector<double>& z,
        const Config::EquilibriumOptimizerConfig& cfg = Config::EquilibriumOptimizerConfig::defaults()
    ) const = 0;

    virtual EOSPtr eos() const = 0;
};

class IHelmholtzOptimizerTV {
public:
    virtual ~IHelmholtzOptimizerTV() = default;

    virtual EquilibriumState minimizeTV(
        double T,
        double V,
        const std::vector<double>& z,
        const Config::EquilibriumOptimizerConfig& cfg = Config::EquilibriumOptimizerConfig::defaults()
    ) const = 0;

    virtual EOSPtr eos() const = 0;
};

using GibbsOptimizerTPPtr = std::shared_ptr<IGibbsOptimizerTP>;
using HelmholtzOptimizerTVPtr = std::shared_ptr<IHelmholtzOptimizerTV>;

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_I_EQUILIBRIUM_OPTIMIZER_H
