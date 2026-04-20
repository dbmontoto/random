#ifndef DMTHERM_DMTHERM_H
#define DMTHERM_DMTHERM_H

#include "dmtherm_export.h"
#include "dmtherm_types.h"

#ifdef __cplusplus
extern "C" {
#endif

// --------------------------------------------------------------------------------------
// Version / memory
// --------------------------------------------------------------------------------------

DMTHERM_API uint32_t dmtherm_abi_version(void);
DMTHERM_API const char* dmtherm_version_string(void);

// Free memory allocated by the library (strings, arrays, etc.).
DMTHERM_API void dmtherm_free(void* p);

// Helpers to initialize versioned config structs with safe defaults.
// Callers should then fill required fields (databanks_dir, eos_type, components, etc.).
DMTHERM_API void dmtherm_system_config_set_defaults(dmtherm_system_config_t* cfg);
DMTHERM_API void dmtherm_flash_config_set_defaults(dmtherm_flash_config_t* cfg);
DMTHERM_API void dmtherm_flash_hs_config_set_defaults(dmtherm_flash_hs_config_t* cfg);

// Convert pressure values between supported pressure units.
// Supported unit tags: PA, PA_ABS, PA_GAUGE, KPA, KPA_ABS, KPA_GAUGE,
// MPA, MPA_ABS, MPA_GAUGE, BAR, BARG, PSI, PSIA, PSIG.
// Returns DMTHERM_STATUS_INVALID_ARGUMENT for unknown/incompatible unit combinations.
DMTHERM_API dmtherm_status_t dmtherm_convert_pressure(
    double value,
    dmtherm_pressure_unit_t from_unit,
    dmtherm_pressure_unit_t to_unit,
    double* out_value
);

// --------------------------------------------------------------------------------------
// System lifetime + errors
// --------------------------------------------------------------------------------------

// Create a new thermo system (databanks + mixture + EOS + solvers).
// On failure, returns non-OK and (optionally) allocates an error message into out_error_message.
DMTHERM_API dmtherm_status_t dmtherm_system_create(
    const dmtherm_system_config_t* cfg,
    dmtherm_system_t** out_system,
    char** out_error_message
);

DMTHERM_API void dmtherm_system_destroy(dmtherm_system_t* sys);

// Returns a pointer to an internal, null-terminated error string owned by the system.
// Valid until the next failing call on the same system, or system destruction.
DMTHERM_API const char* dmtherm_system_last_error(dmtherm_system_t* sys);

// --------------------------------------------------------------------------------------
// Mixture / component queries
// --------------------------------------------------------------------------------------

// Write per-component molecular weights into a caller-provided array.
// Units: g/mol
// `out_mw_g_per_mol` length must equal the system's number of components.
DMTHERM_API dmtherm_status_t dmtherm_system_get_component_mws(
    dmtherm_system_t* sys,
    double* out_mw_g_per_mol,
    size_t out_len
);

// --------------------------------------------------------------------------------------
// TP evaluation
// --------------------------------------------------------------------------------------
// Pressure inputs/outputs are absolute pressure in Pa(abs).

DMTHERM_API dmtherm_status_t dmtherm_system_evaluate_tp_phase(
    dmtherm_system_t* sys,
    double T,
    double P,
    const double* x,
    size_t x_len,
    dmtherm_root_selection_t root,
    dmtherm_tp_phase_eval_t* out
);

DMTHERM_API void dmtherm_tp_phase_eval_destroy(dmtherm_tp_phase_eval_t* out);

DMTHERM_API dmtherm_status_t dmtherm_system_evaluate_tp_phase_batch(
    dmtherm_system_t* sys,
    const double* Ts,
    const double* Ps,
    size_t n_points,
    const double* x,
    size_t x_len,
    dmtherm_root_selection_t root,
    int32_t include_messages,
    dmtherm_tp_phase_eval_batch_t* out
);

DMTHERM_API void dmtherm_tp_phase_eval_batch_destroy(dmtherm_tp_phase_eval_batch_t* out);

DMTHERM_API dmtherm_status_t dmtherm_system_evaluate_tp_phase_batch_into(
    dmtherm_system_t* sys,
    const double* Ts,
    const double* Ps,
    size_t n_points,
    const double* x,
    size_t x_len,
    dmtherm_root_selection_t root,
    int32_t include_messages,
    dmtherm_tp_phase_eval_batch_buffers_t* io
);

DMTHERM_API dmtherm_status_t dmtherm_system_departure_roots_tp(
    dmtherm_system_t* sys,
    double T,
    double P,
    const double* x,
    size_t x_len,
    dmtherm_tp_departure_roots_t* out
);

DMTHERM_API void dmtherm_tp_departure_roots_destroy(dmtherm_tp_departure_roots_t* out);

// --------------------------------------------------------------------------------------
// Critical estimation
// --------------------------------------------------------------------------------------

DMTHERM_API dmtherm_status_t dmtherm_system_estimate_critical(
    dmtherm_system_t* sys,
    int32_t component_index,
    dmtherm_critical_estimate_t* out
);

DMTHERM_API void dmtherm_critical_estimate_destroy(dmtherm_critical_estimate_t* out);

// --------------------------------------------------------------------------------------
// Flash PT
// --------------------------------------------------------------------------------------
// Pressure inputs/outputs are absolute pressure in Pa(abs).

DMTHERM_API dmtherm_status_t dmtherm_system_flash_pt(
    dmtherm_system_t* sys,
    double T,
    double P,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    dmtherm_flash_pt_result_t* out
);

DMTHERM_API dmtherm_status_t dmtherm_system_flash_pt_batch(
    dmtherm_system_t* sys,
    const double* Ts,
    const double* Ps,
    size_t n_points,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    int32_t include_messages,
    dmtherm_flash_pt_batch_t* out
);

DMTHERM_API dmtherm_status_t dmtherm_system_flash_pt_batch_into(
    dmtherm_system_t* sys,
    const double* Ts,
    const double* Ps,
    size_t n_points,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    int32_t include_messages,
    dmtherm_flash_pt_batch_buffers_t* io
);

DMTHERM_API dmtherm_status_t dmtherm_system_flash_tv(
    dmtherm_system_t* sys,
    double T,
    double V,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    dmtherm_flash_pt_result_t* out
);

DMTHERM_API dmtherm_status_t dmtherm_system_flash_ph(
    dmtherm_system_t* sys,
    double P,
    double H,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    dmtherm_flash_pt_result_t* out
);

DMTHERM_API dmtherm_status_t dmtherm_system_flash_ps(
    dmtherm_system_t* sys,
    double P,
    double S,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    dmtherm_flash_pt_result_t* out
);

DMTHERM_API void dmtherm_flash_pt_result_destroy(dmtherm_flash_pt_result_t* out);
DMTHERM_API void dmtherm_flash_pt_batch_destroy(dmtherm_flash_pt_batch_t* out);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // DMTHERM_DMTHERM_H
