/**
 * @file thermochemistry.h
 * @brief Helpers to derive standard-state chemical potentials from databanks.
 */

#ifndef THERMO_REACTIONS_THERMOCHEMISTRY_H
#define THERMO_REACTIONS_THERMOCHEMISTRY_H

#include "thermo/core/mixture.h"
#include "thermo/data/databanks.h"
#include "thermo/data/diagnostics.h"

#include <string>
#include <vector>

namespace DMThermo {
namespace Reactions {

namespace Thermochemistry {

/**
 * @brief String model IDs for thermochemistry sources, used for config/selection.
 *
 * Naming convention: `SOURCE.FIELD` (and `SOURCE.FIELD(T)` for temperature-dependent models).
 */
namespace Models {
inline constexpr const char* DIPPR_GFOR = "DIPPR.GFOR";
inline constexpr const char* DIPPR_MU0_FORMATION_T = "DIPPR.MU0_FORMATION(T)";
} // namespace Models

namespace DIPPR {

/**
 * @brief Build a per-component vector of standard-state chemical potentials from DIPPR `GFOR`.
 *
 * - Uses `Core::Component::cas()` when available, otherwise falls back to unique NAME matching.
 * - `GFOR` is expected in `J/kmol` and is converted to `J/mol`.
 * - Throws if any component cannot be resolved or is missing `GFOR`.
 */
std::vector<double> mu0GFOR(
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    Data::Diagnostics* diag = nullptr);

/**
 * @brief Build a per-component vector of formation chemical potentials at temperature T using DIPPR.
 *
 * Uses:
 * - `GFOR` and `HFOR` at `T_ref` to construct a reference entropy via `s_ref = (h_ref - g_ref)/T_ref`
 * - DIPPR ideal-gas Cp (`ICP_*`) to integrate enthalpy/entropy from `T_ref` to `T`
 *
 * Returned values are in J/mol and represent a consistent formation-based standard chemical potential.
 *
 * Throws if any component is missing `GFOR`, `HFOR`, or `ICP_*` (or if `T`/`T_ref` are out of the Cp validity range).
 */
std::vector<double> mu0FormT(
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    double T,
    double T_ref = 298.15,
    Data::Diagnostics* diag = nullptr);

} // namespace DIPPR

/**
 * @brief Compute standard-state chemical potentials `mu0_i(T)` for a mixture from databanks.
 *
 * Supported models:
 * - `DIPPR.GFOR`: constant mu0 from `GFOR`
 * - `DIPPR.MU0_FORMATION(T)`: formation-basis mu0(T) using `HFOR`,`GFOR`, and ideal-gas `ICP_*`
 */
std::vector<double> mu0DB(
    const std::string& model,
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    double T,
    double T_ref = 298.15,
    Data::Diagnostics* diag = nullptr);

} // namespace Thermochemistry

} // namespace Reactions
} // namespace DMThermo

#endif // THERMO_REACTIONS_THERMOCHEMISTRY_H
