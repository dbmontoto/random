/**
 * @file ideal_gas_cp.h
 * @brief Ideal-gas Cp correlation evaluation helpers (outside Core).
 */

#ifndef THERMO_PROPERTIES_IDEAL_GAS_CP_H
#define THERMO_PROPERTIES_IDEAL_GAS_CP_H

#include "thermo/core/component.h"

namespace DMThermo {
namespace Properties {

double evaluateIdealGasCp(const Core::IdealGasCpCorrelation& corr, double T);
double integrateIdealGasCp(const Core::IdealGasCpCorrelation& corr, double T, double T_ref);
double integrateIdealGasCpOverT(const Core::IdealGasCpCorrelation& corr, double T, double T_ref);

} // namespace Properties
} // namespace DMThermo

#endif // THERMO_PROPERTIES_IDEAL_GAS_CP_H
