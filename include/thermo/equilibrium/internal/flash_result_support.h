/**
 * @file flash_result_support.h
 * @brief Shared helpers to convert optimizer states to flash results (internal).
 */

#ifndef THERMO_EQUILIBRIUM_INTERNAL_FLASH_RESULT_SUPPORT_H
#define THERMO_EQUILIBRIUM_INTERNAL_FLASH_RESULT_SUPPORT_H

#include "thermo/equilibrium/flash_result.h"
#include "thermo/equilibrium/optimization/equilibrium_state.h"

namespace DMThermo {
namespace Equilibrium {
namespace Internal {

inline Flash::FlashResult flashResultFromState(const Optimization::EquilibriumState& st) {
    Flash::FlashResult out;
    out.converged = st.converged;
    out.iterations = st.iterations;
    out.residual = 0.0;
    out.message = st.message;
    out.temperature = st.temperature;
    out.pressure = st.pressure;
    out.z = st.z;
    out.method_used = st.method_used;
    out.phases = st.phases;
    out.num_phases = static_cast<int>(out.phases.size());
    out.is_two_phase = (out.num_phases == 2);
    out.is_three_phase = (out.num_phases == 3);

    double vf = 0.0;
    for (const auto& ph : out.phases) {
        if (ph.type == PhaseType::Vapor) vf += ph.fraction;
    }
    out.vapor_fraction = vf;
    return out;
}

} // namespace Internal
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_INTERNAL_FLASH_RESULT_SUPPORT_H

