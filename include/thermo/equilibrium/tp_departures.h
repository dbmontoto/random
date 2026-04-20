/**
 * @file tp_departures.h
 * @brief TP departure (residual) properties per density root.
 *
 * This provides a process-simulator style "departure function" interface at (T,P,x):
 * - Z
 * - ln(phi_i)
 * - h_dep = H - H_ig(T)                  [J/mol]
 * - s_dep = S - S_ig(T,P,x)              [J/(mol*K)]
 * - g_dep = G - G_ig(T,P,x)              [J/mol]
 *
 * Roots are density roots of P(T,rho,x) = P, typically corresponding to vapor/liquid
 * branches for cubic EOS and other non-ideal models.
 */

#ifndef THERMO_EQUILIBRIUM_TP_DEPARTURES_H
#define THERMO_EQUILIBRIUM_TP_DEPARTURES_H

#include "thermo/config/stability_config.h"
#include "thermo/core/types.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/istability.h"

#include <optional>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Departures {

struct TPDepartureRoot {
    double density = 0.0;                 // mol/m^3
    double Z = 0.0;                       // [-]
    PhaseType phase = PhaseType::Unknown; // Vapor/Liquid hint from root classification
    bool mechanically_stable = false;     // dP/drho > 0
    std::vector<double> ln_phi;           // [-]

    double h_dep = 0.0; // J/mol
    double s_dep = 0.0; // J/(mol*K)
    double g_dep = 0.0; // J/mol
};

struct TPDepartureResult {
    double temperature = 0.0;             // K
    double pressure = 0.0;                // Pa
    std::vector<double> x;                // composition used for evaluation
    std::vector<TPDepartureRoot> roots;   // one entry per density root
    std::optional<int> stable_root_index; // stable single-phase root index (if provided by root finder)
};

enum class RootSelection {
    Vapor,   // minimum-density root (largest Z)
    Liquid,  // maximum-density root (smallest Z)
    Stable   // thermodynamically stable single-phase root (min g among mech-stable)
};

TPDepartureResult departureRootsTP(
    const Stability::IStabilityAnalyzer& stability,
    double T,
    double P,
    const std::vector<double>& x,
    const Config::DensityRootConfig& config = Config::DensityRootConfig{}
);

TPDepartureResult departureRootsTP(
    const Stability::StabilityAnalyzerPtr& stability,
    double T,
    double P,
    const std::vector<double>& x,
    const Config::DensityRootConfig& config = Config::DensityRootConfig{}
);

TPDepartureResult departureRootsTP(
    EOSPtr eos,
    double T,
    double P,
    const std::vector<double>& x,
    const Config::DensityRootConfig& config = Config::DensityRootConfig{}
);

std::optional<int> selectRootIndex(const TPDepartureResult& result, RootSelection selection);
const TPDepartureRoot* selectRoot(const TPDepartureResult& result, RootSelection selection);

} // namespace Departures
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_TP_DEPARTURES_H
