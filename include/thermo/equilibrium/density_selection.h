/**
 * @file density_selection.h
 * @brief Utilities for selecting/refining densities at (T,P,x).
 */

#ifndef THERMO_EQUILIBRIUM_DENSITY_SELECTION_H
#define THERMO_EQUILIBRIUM_DENSITY_SELECTION_H

#include <vector>
#include "thermo/core/types.h"
#include "thermo/core/units.h"
#include "thermo/config/phase_detection.h"

namespace DMThermo {

class EOS;

namespace Equilibrium {
namespace Stability {
class IStabilityAnalyzer;
} // namespace Stability

namespace Detail {

struct PhiEval {
    double rho = 0.0;
    std::vector<double> phi;
    std::vector<double> lnphi;
    PhaseType phase = PhaseType::Unknown;
};

PhiEval evalAtTP(
    const EOS& eos,
    const Stability::IStabilityAnalyzer& stability,
    double T,
    double P,
    const std::vector<double>& x,
    PhaseType phase_hint,
    Config::PhaseDetection space,
    double rho_guess,
    int max_newton_iters,
    const Core::Units::PropertyVariable& pressure_tolerance);

} // namespace Detail
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_DENSITY_SELECTION_H
