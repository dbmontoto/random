#ifndef DMTHERM_DMTHERM_TYPES_H
#define DMTHERM_DMTHERM_TYPES_H

#include <stddef.h>
#include <stdint.h>

// C ABI versioning:
// - Increment DMTHERM_ABI_VERSION when making any breaking change to structs or functions.
// - Versioned structs include `struct_size` to allow forward-compatible extension.
#define DMTHERM_ABI_VERSION 2u

// --------------------------------------------------------------------------------------
// Status / enums
// --------------------------------------------------------------------------------------

typedef enum dmtherm_status {
    DMTHERM_STATUS_OK = 0,
    DMTHERM_STATUS_INVALID_ARGUMENT = 1,
    DMTHERM_STATUS_OUT_OF_MEMORY = 2,
    DMTHERM_STATUS_RUNTIME_ERROR = 3,
    DMTHERM_STATUS_INTERNAL_ERROR = 4
} dmtherm_status_t;

typedef enum dmtherm_phase {
    DMTHERM_PHASE_UNKNOWN = 0,
    DMTHERM_PHASE_VAPOR = 1,
    DMTHERM_PHASE_LIQUID = 2,
    DMTHERM_PHASE_LIQUID1 = 3,
    DMTHERM_PHASE_LIQUID2 = 4,
    DMTHERM_PHASE_SUPERCRITICAL = 5,
    DMTHERM_PHASE_SOLID = 6
} dmtherm_phase_t;

typedef enum dmtherm_root_selection {
    DMTHERM_ROOT_VAPOR = 0,
    DMTHERM_ROOT_LIQUID = 1,
    DMTHERM_ROOT_STABLE = 2
} dmtherm_root_selection_t;

typedef enum dmtherm_property_method {
    DMTHERM_PROPERTY_METHOD_AUTO = 0,
    DMTHERM_PROPERTY_METHOD_EOS_ONLY = 1,
    DMTHERM_PROPERTY_METHOD_DIPPR_ONLY = 2,
    DMTHERM_PROPERTY_METHOD_HYBRID_VAPOR_EOS_LIQUID_DIPPR = 3
} dmtherm_property_method_t;

typedef enum dmtherm_failure_mode {
    DMTHERM_FAILURE_MODE_STRICT = 0,
    DMTHERM_FAILURE_MODE_LENIENT = 1
} dmtherm_failure_mode_t;

// Pressure unit tags used by conversion helpers.
// Conventions:
// - *_ABS are absolute-pressure units.
// - *_GAUGE are gauge-pressure units (offset by 1 atm).
// - PA/KPA/MPA/PSI are differential-pressure units.
// - BAR is treated as absolute by convention; BARG is gauge.
typedef enum dmtherm_pressure_unit {
    DMTHERM_PRESSURE_UNIT_PA = 0,
    DMTHERM_PRESSURE_UNIT_PA_ABS = 1,
    DMTHERM_PRESSURE_UNIT_PA_GAUGE = 2,
    DMTHERM_PRESSURE_UNIT_KPA = 3,
    DMTHERM_PRESSURE_UNIT_KPA_ABS = 4,
    DMTHERM_PRESSURE_UNIT_KPA_GAUGE = 5,
    DMTHERM_PRESSURE_UNIT_MPA = 6,
    DMTHERM_PRESSURE_UNIT_MPA_ABS = 7,
    DMTHERM_PRESSURE_UNIT_MPA_GAUGE = 8,
    DMTHERM_PRESSURE_UNIT_BAR = 9,
    DMTHERM_PRESSURE_UNIT_BARG = 10,
    DMTHERM_PRESSURE_UNIT_PSI = 11,
    DMTHERM_PRESSURE_UNIT_PSIA = 12,
    DMTHERM_PRESSURE_UNIT_PSIG = 13
} dmtherm_pressure_unit_t;

typedef struct dmtherm_optional_double {
    int32_t has_value; // 0 or 1
    double value;
} dmtherm_optional_double_t;

// --------------------------------------------------------------------------------------
// Opaque handles
// --------------------------------------------------------------------------------------

typedef struct dmtherm_system dmtherm_system_t;

// --------------------------------------------------------------------------------------
// System configuration
// --------------------------------------------------------------------------------------

typedef struct dmtherm_system_config {
    uint32_t struct_size; // must be set to sizeof(dmtherm_system_config_t)

    const char* databanks_dir;  // path to CSV databanks directory
    const char* eos_type;       // e.g. "Peng-Robinson", "SRK", "PC-SAFT"
    const char* flash_method;   // e.g. "Auto", "GibbsOptimizer", "HelmholtzOptimizer", "Combined"

    const char* const* components; // array of component identifiers (CAS/CHEMID/name)
    size_t num_components;

    int32_t require_critical_properties; // enforce Tc/Pc/omega presence at build time (recommended for cubics)
    int32_t require_ideal_gas_cp;         // enforce ICP_* presence at build time (required for H/S/G)

    dmtherm_property_method_t property_method; // property routing policy

    // Reserved for future extension.
    // reserved_u32[0]: dmtherm_failure_mode_t (defaults to DMTHERM_FAILURE_MODE_STRICT)
    // remaining reserved slots should be zero-initialized.
    uint32_t reserved_u32[8];
    void* reserved_ptr[8];
} dmtherm_system_config_t;

// --------------------------------------------------------------------------------------
// TP phase evaluation (single selected root)
// --------------------------------------------------------------------------------------

typedef struct dmtherm_tp_phase_eval {
    uint32_t struct_size; // must be set to sizeof(dmtherm_tp_phase_eval_t)

    double temperature_k; // K
    double pressure_pa_abs; // Pa(abs)

    dmtherm_phase_t phase;
    double molar_density_mol_per_m3; // mol/m^3
    double compressibility;

    double enthalpy_departure_j_per_mol; // J/mol
    double entropy_departure_j_per_mol_k; // J/(mol*K)
    double gibbs_departure_j_per_mol; // J/mol

    double enthalpy_j_per_mol; // J/mol
    double entropy_j_per_mol_k; // J/(mol*K)
    double gibbs_j_per_mol; // J/mol

    double* ln_fugacity_coeff; // length = num_components
    size_t ln_fugacity_coeff_len;

    int32_t success; // 0 or 1 (evaluation succeeded at the thermo level)
    char* message;   // optional message (allocated by library; free via dmtherm_free)

    // Heat capacities (total, including ideal + residual when available).
    // Units: J/(mol*K)
    double cp_j_per_mol_k;
    double cv_j_per_mol_k;
} dmtherm_tp_phase_eval_t;

// --------------------------------------------------------------------------------------
// TP phase evaluation (batch: many (T,P) points, one selected root policy)
// --------------------------------------------------------------------------------------

// Batch result buffers are allocated by the library and must be freed by calling
// `dmtherm_tp_phase_eval_batch_destroy()`.
//
// Notes:
// - `status[i]` is the per-point dmtherm_status_t (C ABI call-level) for point i.
// - `success[i]` is the per-point thermo-level success flag returned by the engine.
// - `ln_fugacity_coeff` is a flattened array in row-major order with stride `ln_fugacity_coeff_stride`
//   (typically `num_components`).
// - Messages (optional) are stored in a single null-separated buffer. The message for point i
//   starts at `messages + message_offsets[i]`.
typedef struct dmtherm_tp_phase_eval_batch {
    uint32_t struct_size; // must be set to sizeof(dmtherm_tp_phase_eval_batch_t)

    size_t n_points;
    size_t num_components;

    int32_t* status; // length n_points (values are dmtherm_status_t)
    int32_t* success; // length n_points (0/1)
    dmtherm_phase_t* phase; // length n_points

    double* molar_density_mol_per_m3; // length n_points
    double* compressibility; // length n_points

    double* enthalpy_departure_j_per_mol; // length n_points
    double* entropy_departure_j_per_mol_k; // length n_points
    double* gibbs_departure_j_per_mol; // length n_points

    double* enthalpy_j_per_mol; // length n_points
    double* entropy_j_per_mol_k; // length n_points
    double* gibbs_j_per_mol; // length n_points

    double* ln_fugacity_coeff; // flattened array, length ln_fugacity_coeff_len
    size_t ln_fugacity_coeff_len;
    size_t ln_fugacity_coeff_stride; // typically num_components

    size_t* message_offsets; // length n_points (optional; may be null)
    char* messages;          // null-separated messages buffer (optional; may be null)
    size_t messages_len;     // bytes in messages buffer

    // Reserved for future extension (must be zero-initialized).
    uint32_t reserved_u32[8];
    void* reserved_ptr[8];

    // Heat capacities (total, including ideal + residual when available).
    // Arrays are length n_points.
    double* cp_j_per_mol_k; // J/(mol*K)
    double* cv_j_per_mol_k; // J/(mol*K)
} dmtherm_tp_phase_eval_batch_t;

// Caller-provided buffers variant (no allocations in the DLL).
// The caller must allocate all arrays and pass pointers + capacities via this struct.
// No destroy function is required (and none should be called) for the caller-owned buffers.
//
// Messages are optional:
// - If include_messages=0, message pointers may be null and will not be written.
// - If include_messages=1, the caller must provide `message_offsets` (length n_points) and a writable `messages` buffer.
//   Each point's message is written as a null-terminated string at `messages + message_offsets[i]`.
//   If the messages buffer is too small, messages will be truncated and `messages_truncated` will be set to 1.
typedef struct dmtherm_tp_phase_eval_batch_buffers {
    uint32_t struct_size; // must be set to sizeof(dmtherm_tp_phase_eval_batch_buffers_t)

    // These must match the call arguments.
    size_t n_points;
    size_t num_components;

    // Caller-provided output buffers (length n_points unless otherwise noted).
    int32_t* status; // dmtherm_status_t as int32
    int32_t* success; // 0/1
    dmtherm_phase_t* phase;

    double* molar_density_mol_per_m3;
    double* compressibility;

    double* enthalpy_departure_j_per_mol;
    double* entropy_departure_j_per_mol_k;
    double* gibbs_departure_j_per_mol;

    double* enthalpy_j_per_mol;
    double* entropy_j_per_mol_k;
    double* gibbs_j_per_mol;

    // Flattened array in row-major order with stride `ln_fugacity_coeff_stride` (must be >= num_components).
    // Capacity is ln_fugacity_coeff_len doubles.
    double* ln_fugacity_coeff;
    size_t ln_fugacity_coeff_len;
    size_t ln_fugacity_coeff_stride;

    // Optional packed messages buffer (caller-provided).
    size_t* message_offsets; // length n_points
    char* messages;
    size_t messages_capacity; // bytes

    // Output (written by DLL).
    size_t messages_len;      // bytes used (including null terminators)
    int32_t messages_truncated; // 0/1

    uint32_t reserved_u32[8];
    void* reserved_ptr[8];

    // Optional heat capacity outputs (arrays length n_points).
    // If pointers are null, the DLL will skip writing these.
    double* cp_j_per_mol_k; // J/(mol*K)
    double* cv_j_per_mol_k; // J/(mol*K)
} dmtherm_tp_phase_eval_batch_buffers_t;

// --------------------------------------------------------------------------------------
// TP departure roots (all density roots at TP)
// --------------------------------------------------------------------------------------

typedef struct dmtherm_tp_departure_root {
    dmtherm_phase_t phase;
    int32_t mechanically_stable; // 0 or 1

    double molar_density_mol_per_m3; // mol/m^3
    double compressibility;

    double enthalpy_departure_j_per_mol; // J/mol
    double entropy_departure_j_per_mol_k; // J/(mol*K)
    double gibbs_departure_j_per_mol; // J/mol

    double enthalpy_j_per_mol; // J/mol
    double entropy_j_per_mol_k; // J/(mol*K)
    double gibbs_j_per_mol; // J/mol

    double* ln_fugacity_coeff; // length = num_components
    size_t ln_fugacity_coeff_len;
} dmtherm_tp_departure_root_t;

typedef struct dmtherm_tp_departure_roots {
    uint32_t struct_size; // must be set to sizeof(dmtherm_tp_departure_roots_t)

    double temperature_k; // K
    double pressure_pa_abs; // Pa(abs)

    int32_t stable_root_index; // -1 if unknown/unavailable

    dmtherm_tp_departure_root_t* roots; // array length roots_len
    size_t roots_len;

    char* message; // optional message (allocated by library; free via dmtherm_free)
} dmtherm_tp_departure_roots_t;

// --------------------------------------------------------------------------------------
// Critical point estimation (pure-component limit)
// --------------------------------------------------------------------------------------

typedef struct dmtherm_critical_estimate {
    uint32_t struct_size; // must be set to sizeof(dmtherm_critical_estimate_t)

    int32_t component_index; // input component index used for estimation
    int32_t converged;       // 0/1
    int32_t iterations;      // solver iterations (0 for non-iterative fallback branches)

    double critical_temperature_k; // K
    double critical_pressure_pa_abs; // Pa(abs)
    double critical_density_mol_per_m3; // mol/m^3
    double critical_compressibility; // [-]
    double critical_volume_m3_per_kmol; // m^3/kmol

    char* method;  // optional ("newton" or "spinodal")
    char* message; // optional detail message
} dmtherm_critical_estimate_t;

// --------------------------------------------------------------------------------------
// Flash PT (minimal C ABI surface; expands over time)
// --------------------------------------------------------------------------------------

typedef struct dmtherm_flash_config {
    uint32_t struct_size; // must be set to sizeof(dmtherm_flash_config_t)
    int32_t max_phases;   // 1..3
    int32_t compute_properties; // 0/1 (populate per-phase/overall H/S/G and departures when available)
    uint32_t reserved_u32[8];
    void* reserved_ptr[8];
} dmtherm_flash_config_t;

// Optional PH/PS temperature solve controls (tolerance + bracket).
// To use with the flash APIs:
//   - Set `dmtherm_flash_config_t.reserved_ptr[0] = &hs_cfg`
//   - Set `hs_cfg.struct_size = sizeof(dmtherm_flash_hs_config_t)`
//
// Interpretation:
//   - For PH flash, `tolerance` is in J/mol.
//   - For PS flash, `tolerance` is in J/(mol*K).
typedef struct dmtherm_flash_hs_config {
    uint32_t struct_size; // must be set to sizeof(dmtherm_flash_hs_config_t)
    double tolerance;
    double temp_bracket_low;   // K
    double temp_bracket_high;  // K
    int32_t max_iterations;    // <=0 means use engine default
    int32_t use_newton;        // 0/1 (currently best-effort; solver may ignore)
    uint32_t reserved_u32[8];
    void* reserved_ptr[8];
} dmtherm_flash_hs_config_t;

typedef struct dmtherm_flash_phase {
    dmtherm_phase_t phase;
    double fraction;

    double molar_density_mol_per_m3; // mol/m^3
    double compressibility;

    double* composition; // composition (length = num_components)
    size_t composition_len;

    double* fugacity_coeff; // fugacity coeffs (length = num_components), if available
    size_t fugacity_coeff_len;

    dmtherm_optional_double_t enthalpy_departure_j_per_mol;
    dmtherm_optional_double_t entropy_departure_j_per_mol_k;
    dmtherm_optional_double_t gibbs_departure_j_per_mol;

    dmtherm_optional_double_t enthalpy_j_per_mol;
    dmtherm_optional_double_t entropy_j_per_mol_k;
    dmtherm_optional_double_t gibbs_j_per_mol;
} dmtherm_flash_phase_t;

typedef struct dmtherm_flash_pt_result {
    uint32_t struct_size; // must be set to sizeof(dmtherm_flash_pt_result_t)

    int32_t converged;
    int32_t iterations;
    double residual;
    char* message; // optional

    double temperature_k; // K
    double pressure_pa_abs; // Pa(abs)

    double vapor_fraction;
    int32_t num_phases;

    dmtherm_flash_phase_t* phases; // array length phases_len
    size_t phases_len;

    dmtherm_optional_double_t mixture_enthalpy_j_per_mol;
    dmtherm_optional_double_t mixture_entropy_j_per_mol_k;
    dmtherm_optional_double_t mixture_gibbs_j_per_mol;

    char* method_used; // optional
} dmtherm_flash_pt_result_t;

// --------------------------------------------------------------------------------------
// Flash PT (batch; fixed strides for phases and compositions)
// --------------------------------------------------------------------------------------

// Notes:
// - Results are stored in fixed-stride arrays to avoid per-point/per-phase heap allocations.
// - Phases for point i occupy indices [i*max_phases, i*max_phases + max_phases).
//   Only the first `num_phases[i]` entries are valid; the rest are filled with zeros/Unknown.
// - Composition arrays are flattened with stride `num_components`:
//   `phase_composition[(i*max_phases + j)*num_components + k]`.
typedef struct dmtherm_flash_pt_batch {
    uint32_t struct_size; // must be set to sizeof(dmtherm_flash_pt_batch_t)

    size_t n_points;
    size_t num_components;
    int32_t max_phases; // fixed stride per point

    // Per-point outputs (length n_points).
    int32_t* status;     // dmtherm_status_t as int32
    int32_t* converged;  // 0/1
    int32_t* iterations;
    double* residual;
    double* temperature_k; // K
    double* pressure_pa_abs; // Pa(abs)
    double* vapor_fraction;
    int32_t* num_phases;

    dmtherm_optional_double_t* mixture_enthalpy_j_per_mol;
    dmtherm_optional_double_t* mixture_entropy_j_per_mol_k;
    dmtherm_optional_double_t* mixture_gibbs_j_per_mol;

    // Per-phase outputs (length phases_len = n_points * max_phases).
    size_t phases_len;
    dmtherm_phase_t* phase;
    double* phase_fraction;
    double* phase_molar_density_mol_per_m3;
    double* phase_compressibility;

    // Flattened per-phase compositions (length phases_len * num_components).
    size_t phase_composition_len;
    double* phase_composition;

    // Flattened per-phase fugacity coefficients (length phases_len * num_components).
    // If `phase_fugacity_coeff_available[idx] == 0`, the corresponding coeff slice is not valid.
    int32_t* phase_fugacity_coeff_available; // length phases_len
    size_t phase_fugacity_coeff_len;
    double* phase_fugacity_coeff;

    dmtherm_optional_double_t* phase_enthalpy_departure_j_per_mol;
    dmtherm_optional_double_t* phase_entropy_departure_j_per_mol_k;
    dmtherm_optional_double_t* phase_gibbs_departure_j_per_mol;
    dmtherm_optional_double_t* phase_enthalpy_j_per_mol;
    dmtherm_optional_double_t* phase_entropy_j_per_mol_k;
    dmtherm_optional_double_t* phase_gibbs_j_per_mol;

    // Optional packed messages per point.
    size_t* message_offsets; // length n_points (optional; may be null)
    char* messages;          // null-separated buffer (optional; may be null)
    size_t messages_len;

    uint32_t reserved_u32[8];
    void* reserved_ptr[8];
} dmtherm_flash_pt_batch_t;

// Caller-provided buffers variant (no allocations in the DLL).
// The caller must allocate all arrays and pass pointers + capacities via this struct.
// No destroy function is required (and none should be called) for the caller-owned buffers.
//
// Messages are optional:
// - If include_messages=0, message pointers may be null and will not be written.
// - If include_messages=1, the caller must provide `message_offsets` (length n_points) and a writable `messages` buffer.
//   Each point's message is written as a null-terminated string at `messages + message_offsets[i]`.
//   If the messages buffer is too small, messages will be truncated and `messages_truncated` will be set to 1.
typedef struct dmtherm_flash_pt_batch_buffers {
    uint32_t struct_size; // must be set to sizeof(dmtherm_flash_pt_batch_buffers_t)

    // These must match the call arguments.
    size_t n_points;
    size_t num_components;
    int32_t max_phases; // fixed stride per point

    // Caller-provided per-point output buffers (length n_points).
    int32_t* status;     // dmtherm_status_t as int32
    int32_t* converged;  // 0/1
    int32_t* iterations;
    double* residual;
    double* temperature_k; // K
    double* pressure_pa_abs; // Pa(abs)
    double* vapor_fraction;
    int32_t* num_phases;

    dmtherm_optional_double_t* mixture_enthalpy_j_per_mol;
    dmtherm_optional_double_t* mixture_entropy_j_per_mol_k;
    dmtherm_optional_double_t* mixture_gibbs_j_per_mol;

    // Caller-provided per-phase output buffers.
    // Capacity is phases_len entries, typically n_points * max_phases.
    size_t phases_len;
    dmtherm_phase_t* phase;
    double* phase_fraction;
    double* phase_molar_density_mol_per_m3;
    double* phase_compressibility;

    // Flattened per-phase compositions (length phases_len * num_components).
    size_t phase_composition_len;
    double* phase_composition;

    // Flattened per-phase fugacity coefficients (length phases_len * num_components).
    // If `phase_fugacity_coeff_available[idx] == 0`, the corresponding coeff slice is not valid.
    int32_t* phase_fugacity_coeff_available; // length phases_len
    size_t phase_fugacity_coeff_len;
    double* phase_fugacity_coeff;

    dmtherm_optional_double_t* phase_enthalpy_departure_j_per_mol;
    dmtherm_optional_double_t* phase_entropy_departure_j_per_mol_k;
    dmtherm_optional_double_t* phase_gibbs_departure_j_per_mol;
    dmtherm_optional_double_t* phase_enthalpy_j_per_mol;
    dmtherm_optional_double_t* phase_entropy_j_per_mol_k;
    dmtherm_optional_double_t* phase_gibbs_j_per_mol;

    // Optional packed messages buffer (caller-provided).
    size_t* message_offsets; // length n_points
    char* messages;
    size_t messages_capacity; // bytes

    // Output (written by DLL).
    size_t messages_len;      // bytes used (including null terminators)
    int32_t messages_truncated; // 0/1

    uint32_t reserved_u32[8];
    void* reserved_ptr[8];
} dmtherm_flash_pt_batch_buffers_t;

#endif // DMTHERM_DMTHERM_TYPES_H
