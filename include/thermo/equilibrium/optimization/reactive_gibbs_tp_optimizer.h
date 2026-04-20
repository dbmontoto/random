/**
 * @file reactive_gibbs_tp_optimizer.h
 * @brief EOS-agnostic multi-phase Gibbs minimizer at constant T,P with reaction extents.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_REACTIVE_GIBBS_TP_OPTIMIZER_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_REACTIVE_GIBBS_TP_OPTIMIZER_H

#include "equilibrium_state.h"
#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/reactions/reaction_system.h"

#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

struct ReactiveGibbsTPConfig {
    Config::EquilibriumOptimizerConfig equilibrium = Config::EquilibriumOptimizerConfig::defaults();
    double standard_pressure = 101325.0; // Pa
    double min_moles = 1e-30;
};

class ReactiveGibbsTPOptimizer final {
public:
    explicit ReactiveGibbsTPOptimizer(
        EOSPtr eos,
        Stability::StabilityAnalyzerPtr stability = nullptr);

    /**
     * @brief Minimize total Gibbs energy at constant T,P with reactions and optional phase split.
     *
     * Decision variables:
     * - Reaction extents (bounded via conservative per-reaction extent bounds)
     * - Component allocations across phases (softmax logits)
     *
     * Objective (dimensionless):
     *   G/(RT) = sum_i n_i * (mu0_i/(RT) + ln(P/P0)) + sum_k n_k * sum_i x_{ik}(ln x_{ik} + ln phi_{ik})
     *
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param n0 Initial component moles
     * @param system Reaction system (stoichiometry on components)
     * @param mu0 Standard chemical potentials at P0 [J/mol]
     * @param cfg Reactive config (phase split + constraints)
     */
    EquilibriumState minimizeTPReactive(
        double T,
        double P,
        const std::vector<double>& n0,
        const Reactions::ReactionSystem& system,
        const std::vector<double>& mu0,
        const ReactiveGibbsTPConfig& cfg = ReactiveGibbsTPConfig{}) const;

    EOSPtr eos() const { return eos_; }

private:
    EOSPtr eos_;
    Stability::StabilityAnalyzerPtr stability_;
};

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_REACTIVE_GIBBS_TP_OPTIMIZER_H
