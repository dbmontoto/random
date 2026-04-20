/**
 * @file unifac.h
 * @brief UNIFAC activity coefficient model (original-form combinatorial + residual)
 */

#ifndef THERMO_ACTIVITY_UNIFAC_H
#define THERMO_ACTIVITY_UNIFAC_H

#include "thermo/core/mixture.h"
#include <vector>

namespace DMThermo {
namespace Activity {

/**
 * @brief Compute UNIFAC ln(gamma_i) for a mixture at (T,x).
 *
 * Uses subgroup R/Q parameters and main-group interaction parameters from
 * `mixture.unifacTables()`. Interactions missing from the table default to 0
 * (i.e. tau=1).
 */
std::vector<double> lnGammaUNIFAC(
    double T,
    const Core::Mixture& mixture,
    const std::vector<double>& x
);

/**
 * @brief Compute UNIFAC excess Gibbs energy G^E/(RT) = Σ x_i ln(gamma_i).
 */
double excessGibbsOverRTUNIFAC(
    double T,
    const Core::Mixture& mixture,
    const std::vector<double>& x,
    std::vector<double>* ln_gamma_out = nullptr
);

} // namespace Activity
} // namespace DMThermo

#endif // THERMO_ACTIVITY_UNIFAC_H

