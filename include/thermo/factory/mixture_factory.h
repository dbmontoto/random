/**
 * @file mixture_factory.h
 * @brief Construct DMThermo::Core types from CSV databanks
 */

#ifndef THERMO_FACTORY_MIXTURE_FACTORY_H
#define THERMO_FACTORY_MIXTURE_FACTORY_H

#include "thermo/core/component.h"
#include "thermo/core/mixture.h"
#include "thermo/data/databanks.h"
#include "thermo/data/diagnostics.h"

#include <string>
#include <vector>

namespace DMThermo {
namespace Factory {

struct MixtureBuildOptions {
    Data::BipModel bip_model = Data::BipModel::PengRobinson;
    double bip_reference_temperature = 298.15; // K (BIPs are stored as constants today)
    bool require_pcsaft_params = false;
    bool require_critical_properties = false; // Tc/Pc/omega (needed for cubic EOS + Wilson-based init)
    bool require_ideal_gas_cp = false;         // DIPPR ICP_* (needed for PH/PS and mu0(T) formation models)
    bool require_vtpr_params = false;          // VTPR: requires vtpr_pure.csv + vtpr_groups.csv + UNIFAC tables
};

/**
 * @brief Build a single component from databanks using CAS or name identifier.
 *
 * @param databanks CSV-backed databanks
 * @param identifier CAS (preferred) or component name
 * @param options Build options (only `require_pcsaft_params` is relevant here)
 * @param diag Optional diagnostics sink
 */
Core::Component buildComponentFromDatabanks(
    const Data::Databanks& databanks,
    const std::string& identifier,
    const MixtureBuildOptions& options,
    Data::Diagnostics* diag = nullptr
);

/**
 * @brief Build a mixture from databanks using a list of CAS/name identifiers.
 *
 * The databanks are allowed to be incomplete. Missing BIPs default to zero.
 * Components may be built without PC-SAFT parameters (for cubic EOS usage).
 * For PC-SAFT usage, set `require_pcsaft_params = true` (missing params will throw).
 */
Core::Mixture buildMixtureFromDatabanks(
    const Data::Databanks& databanks,
    const std::vector<std::string>& identifiers,
    const MixtureBuildOptions& options,
    Data::Diagnostics* diag = nullptr
);

} // namespace Factory
} // namespace DMThermo

#endif // THERMO_FACTORY_MIXTURE_FACTORY_H
