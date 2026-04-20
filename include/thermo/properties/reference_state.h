/**
 * @file reference_state.h
 * @brief Reference state conventions for thermodynamic properties
 */

#ifndef THERMO_PROPERTIES_REFERENCE_STATE_H
#define THERMO_PROPERTIES_REFERENCE_STATE_H

namespace DMThermo {
namespace Properties {

struct ReferenceState {
    double T_ref = 298.15;   // K
    double P_ref = 1e5;      // Pa
};

inline ReferenceState defaultReferenceState() {
    return ReferenceState{};
}

} // namespace Properties
} // namespace DMThermo

#endif // THERMO_PROPERTIES_REFERENCE_STATE_H

