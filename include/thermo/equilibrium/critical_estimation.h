/**
 * @file critical_estimation.h
 * @brief EOS-based critical point estimation for pure-component limits.
 */

#ifndef THERMO_EQUILIBRIUM_CRITICAL_ESTIMATION_H
#define THERMO_EQUILIBRIUM_CRITICAL_ESTIMATION_H

#include "thermo/core/units.h"
#include "thermo/eos.h"

#include <limits>
#include <string>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Critical {

struct EstimateInputs {
    double Tc_hint = std::numeric_limits<double>::quiet_NaN();
    double Pc_hint = std::numeric_limits<double>::quiet_NaN();
    double Vc_hint_m3_per_kmol = std::numeric_limits<double>::quiet_NaN();

    Core::Units::PropertyVariable temperatureHint() const {
        return Core::Units::PropertyVariable(Tc_hint, Core::Units::Unit::K);
    }

    Core::Units::PropertyVariable pressureHint() const {
        return Core::Units::PropertyVariable(Pc_hint, Core::Units::Unit::PA_ABS);
    }

    Core::Units::PropertyVariable criticalVolumeHint() const {
        return Core::Units::PropertyVariable(Vc_hint_m3_per_kmol, Core::Units::Unit::M3_PER_KMOL);
    }
};

struct EstimateResult {
    bool converged = false;
    int iterations = 0;

    double Tc = std::numeric_limits<double>::quiet_NaN();
    double Pc = std::numeric_limits<double>::quiet_NaN();
    double rho = std::numeric_limits<double>::quiet_NaN();
    double Zc = std::numeric_limits<double>::quiet_NaN();
    double Vc_m3_per_kmol = std::numeric_limits<double>::quiet_NaN();

    // "newton" or "spinodal" when converged.
    std::string method;
    std::string message;

    Core::Units::PropertyVariable temperatureValue() const {
        return Core::Units::PropertyVariable(Tc, Core::Units::Unit::K);
    }

    Core::Units::PropertyVariable pressureValue() const {
        return Core::Units::PropertyVariable(Pc, Core::Units::Unit::PA_ABS);
    }

    Core::Units::PropertyVariable densityValue() const {
        return Core::Units::PropertyVariable(rho, Core::Units::Unit::MOL_PER_M3);
    }

    Core::Units::PropertyVariable criticalVolumeValue() const {
        return Core::Units::PropertyVariable(Vc_m3_per_kmol, Core::Units::Unit::M3_PER_KMOL);
    }
};

/**
 * @brief Estimate critical point from an EOS and composition.
 *
 * The caller should pass a pure-component composition vector (one-hot for
 * multi-component EOS instances). For single-component EOS objects, x can be
 * omitted and defaults to {1.0}.
 */
EstimateResult estimateFromEOS(
    const EOS& eos,
    const std::vector<double>& x,
    const EstimateInputs& inputs = EstimateInputs{}
);

} // namespace Critical
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_CRITICAL_ESTIMATION_H
