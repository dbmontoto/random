/**
 * @file dippr_adapter.h
 * @brief Adapt CSV DIPPR records into the existing DIPPR::Parameters model
 */

#ifndef THERMO_DATA_DIPPR_ADAPTER_H
#define THERMO_DATA_DIPPR_ADAPTER_H

#include "thermo/data/databanks.h"
#include "thermo/data/diagnostics.h"
#include "thermo/data/records.h"
#include "thermo/dippr/equations.h"

#include <optional>

namespace DMThermo {
namespace Data {

/**
 * @brief Convert a CSV DIPPR record to DIPPR::Parameters (best-effort).
 *
 * Only equation forms supported by DIPPR::EquationType are converted.
 * Unsupported equation forms are skipped with warnings.
 */
std::optional<DIPPR::Parameters> toDipprParameters(const DipprRecord& record, Diagnostics* diag = nullptr);

} // namespace Data
} // namespace DMThermo

#endif // THERMO_DATA_DIPPR_ADAPTER_H
