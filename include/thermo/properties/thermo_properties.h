/**
 * @file thermo_properties.h
 * @brief Thermodynamic properties computed from EOS + ideal-gas correlations
 */

#ifndef THERMO_PROPERTIES_THERMO_PROPERTIES_H
#define THERMO_PROPERTIES_THERMO_PROPERTIES_H

#include "thermo/core/mixture.h"
#include "thermo/eos.h"
#include "reference_state.h"

#include <optional>
#include <string>
#include <vector>

namespace DMThermo {
namespace Properties {

struct ThermoProperties {
    double T = 0.0;
    double P = 0.0;
    double rho = 0.0;
    double v = 0.0;
    double Z = 0.0;

    double H = 0.0;   // J/mol
    double S = 0.0;   // J/(mol*K)
    double U = 0.0;   // J/mol
    double A = 0.0;   // J/mol
    double G = 0.0;   // J/mol

    std::optional<double> Cp; // J/(mol*K)
    std::optional<double> Cv; // J/(mol*K)

    // Transport properties (optional; not all models provide these)
    std::optional<double> mu; // Pa*s
    std::optional<double> k;  // W/(m*K)

    std::string model;        // EOS name
    bool success = false;
    std::string message;
};

/**
 * @brief Compute properties at (T, rho) for composition x using EOS residuals + ideal contributions.
 *
 * Requirements:
 * - `eos.calculate()` must provide `a_residual` and `da_dT` to compute S_res and U_res.
 * - Ideal-gas Cp correlations are taken from `mixture.component(i).idealGasCpCorrelation()`;
 *   correlation parsing/evaluation is handled in `Properties` (not `Core`).
 */
ThermoProperties calculateTV(
    const EOS& eos,
    const Core::Mixture& mixture,
    double T,
    double rho,
    const std::vector<double>& x,
    const ReferenceState& ref = defaultReferenceState()
);

} // namespace Properties
} // namespace DMThermo

#endif // THERMO_PROPERTIES_THERMO_PROPERTIES_H
