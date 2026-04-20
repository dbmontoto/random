#include "dmtherm/dmtherm.h"

#include "thermo/core/units.h"
#include "thermo/engine/thermo_system.h"
#include "thermo/equilibrium/thermo_contributions.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <mutex>
#include <optional>
#include <new>
#include <string>
#include <vector>

using namespace DMThermo;

namespace {

constexpr uint32_t requiredSizeTpPhaseEvalV1() {
    return static_cast<uint32_t>(offsetof(dmtherm_tp_phase_eval_t, cp_j_per_mol_k));
}

constexpr uint32_t requiredSizeTpPhaseEvalBatchV1() {
    return static_cast<uint32_t>(offsetof(dmtherm_tp_phase_eval_batch_t, cp_j_per_mol_k));
}

constexpr uint32_t requiredSizeTpPhaseEvalBatchBuffersV1() {
    return static_cast<uint32_t>(offsetof(dmtherm_tp_phase_eval_batch_buffers_t, cp_j_per_mol_k));
}

constexpr uint32_t requiredSizeCriticalEstimateV1() {
    return static_cast<uint32_t>(offsetof(dmtherm_critical_estimate_t, method));
}

bool hasTpPhaseEvalHeatCaps(const dmtherm_tp_phase_eval_t* out) {
    if (!out) return false;
    const uint32_t need = static_cast<uint32_t>(offsetof(dmtherm_tp_phase_eval_t, cv_j_per_mol_k) + sizeof(double));
    return out->struct_size >= need;
}

bool hasTpPhaseEvalBatchHeatCaps(const dmtherm_tp_phase_eval_batch_t* out) {
    if (!out) return false;
    const uint32_t need = static_cast<uint32_t>(offsetof(dmtherm_tp_phase_eval_batch_t, cv_j_per_mol_k) + sizeof(double*));
    return out->struct_size >= need;
}

bool hasTpPhaseEvalBatchBuffersHeatCaps(const dmtherm_tp_phase_eval_batch_buffers_t* out) {
    if (!out) return false;
    const uint32_t need = static_cast<uint32_t>(offsetof(dmtherm_tp_phase_eval_batch_buffers_t, cv_j_per_mol_k) + sizeof(double*));
    return out->struct_size >= need;
}

struct SystemHandle {
    Engine::ThermoSystem sys;
    std::mutex mu;
    std::string last_error;
};

dmtherm_phase_t toPhase(PhaseType ph) {
    switch (ph) {
        case PhaseType::Vapor: return DMTHERM_PHASE_VAPOR;
        case PhaseType::Liquid: return DMTHERM_PHASE_LIQUID;
        case PhaseType::Liquid1: return DMTHERM_PHASE_LIQUID1;
        case PhaseType::Liquid2: return DMTHERM_PHASE_LIQUID2;
        case PhaseType::Supercritical: return DMTHERM_PHASE_SUPERCRITICAL;
        case PhaseType::Solid: return DMTHERM_PHASE_SOLID;
        case PhaseType::Unknown:
        default: return DMTHERM_PHASE_UNKNOWN;
    }
}

Engine::PropertyMethod toPropertyMethod(dmtherm_property_method_t m) {
    switch (m) {
        case DMTHERM_PROPERTY_METHOD_AUTO: return Engine::PropertyMethod::AUTO;
        case DMTHERM_PROPERTY_METHOD_EOS_ONLY: return Engine::PropertyMethod::EOS_ONLY;
        case DMTHERM_PROPERTY_METHOD_DIPPR_ONLY: return Engine::PropertyMethod::DIPPR_ONLY;
        case DMTHERM_PROPERTY_METHOD_HYBRID_VAPOR_EOS_LIQUID_DIPPR: return Engine::PropertyMethod::HYBRID_VAPOR_EOS_LIQUID_DIPPR;
        default: return Engine::PropertyMethod::EOS_ONLY;
    }
}

std::optional<Config::FailureMode> toFailureMode(uint32_t raw) {
    switch (raw) {
        case DMTHERM_FAILURE_MODE_STRICT: return Config::FailureMode::Strict;
        case DMTHERM_FAILURE_MODE_LENIENT: return Config::FailureMode::Lenient;
        default: return std::nullopt;
    }
}

Equilibrium::Departures::RootSelection toRootSelection(dmtherm_root_selection_t sel) {
    switch (sel) {
        case DMTHERM_ROOT_VAPOR: return Equilibrium::Departures::RootSelection::Vapor;
        case DMTHERM_ROOT_LIQUID: return Equilibrium::Departures::RootSelection::Liquid;
        case DMTHERM_ROOT_STABLE:
        default: return Equilibrium::Departures::RootSelection::Stable;
    }
}

std::optional<Core::Units::Unit> toCorePressureUnit(dmtherm_pressure_unit_t unit) {
    using Core::Units::Unit;
    switch (unit) {
        case DMTHERM_PRESSURE_UNIT_PA: return Unit::PA;
        case DMTHERM_PRESSURE_UNIT_PA_ABS: return Unit::PA_ABS;
        case DMTHERM_PRESSURE_UNIT_PA_GAUGE: return Unit::PA_GAUGE;
        case DMTHERM_PRESSURE_UNIT_KPA: return Unit::KPA;
        case DMTHERM_PRESSURE_UNIT_KPA_ABS: return Unit::KPA_ABS;
        case DMTHERM_PRESSURE_UNIT_KPA_GAUGE: return Unit::KPA_GAUGE;
        case DMTHERM_PRESSURE_UNIT_MPA: return Unit::MPA;
        case DMTHERM_PRESSURE_UNIT_MPA_ABS: return Unit::MPA_ABS;
        case DMTHERM_PRESSURE_UNIT_MPA_GAUGE: return Unit::MPA_GAUGE;
        case DMTHERM_PRESSURE_UNIT_BAR: return Unit::BAR;
        case DMTHERM_PRESSURE_UNIT_BARG: return Unit::BARG;
        case DMTHERM_PRESSURE_UNIT_PSI: return Unit::PSI;
        case DMTHERM_PRESSURE_UNIT_PSIA: return Unit::PSIA;
        case DMTHERM_PRESSURE_UNIT_PSIG: return Unit::PSIG;
        default: return std::nullopt;
    }
}

void setLastError(SystemHandle* h, const std::string& msg) {
    if (!h) return;
    std::lock_guard<std::mutex> lock(h->mu);
    h->last_error = msg;
}

char* dupCString(const std::string& s) {
    const size_t n = s.size();
    char* p = static_cast<char*>(std::malloc(n + 1));
    if (!p) return nullptr;
    std::memcpy(p, s.c_str(), n);
    p[n] = '\0';
    return p;
}

template <class T>
T* allocArray(size_t n) {
    if (n == 0) return nullptr;
    if (n > (static_cast<size_t>(-1) / sizeof(T))) return nullptr;
    return static_cast<T*>(std::malloc(n * sizeof(T)));
}

template <class T>
bool copyArray(const T* src, size_t n, T*& dst) {
    dst = allocArray<T>(n);
    if (n != 0 && !dst) return false;
    if (n != 0 && src) {
        std::memcpy(dst, src, n * sizeof(T));
    }
    return true;
}

dmtherm_optional_double_t opt(const std::optional<double>& v) {
    dmtherm_optional_double_t out{};
    if (v.has_value()) {
        out.has_value = 1;
        out.value = v.value();
    }
    return out;
}

dmtherm_status_t fillFlashResult(
    SystemHandle* h,
    const Equilibrium::Flash::FlashResult& fr,
    dmtherm_flash_pt_result_t* out)
{
    out->converged = fr.converged ? 1 : 0;
    out->iterations = fr.iterations;
    out->residual = fr.residual;
    out->temperature_k = fr.temperature;
    out->pressure_pa_abs = fr.pressure;
    out->vapor_fraction = fr.vapor_fraction;
    out->num_phases = fr.num_phases;

    if (!fr.message.empty()) {
        out->message = dupCString(fr.message);
        if (!out->message) {
            setLastError(h, "flash: out of memory (message)");
            return DMTHERM_STATUS_OUT_OF_MEMORY;
        }
    }

    if (!fr.method_used.empty()) {
        out->method_used = dupCString(fr.method_used);
        if (!out->method_used) {
            setLastError(h, "flash: out of memory (method_used)");
            return DMTHERM_STATUS_OUT_OF_MEMORY;
        }
    }

    out->mixture_enthalpy_j_per_mol = opt(fr.enthalpy);
    out->mixture_entropy_j_per_mol_k = opt(fr.entropy);
    out->mixture_gibbs_j_per_mol = opt(fr.gibbs_energy);

    out->phases_len = fr.phases.size();
    out->phases = allocArray<dmtherm_flash_phase_t>(out->phases_len);
    if (out->phases_len != 0 && !out->phases) {
        setLastError(h, "flash: out of memory (phases)");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    }
    if (out->phases_len != 0) {
        std::memset(out->phases, 0, out->phases_len * sizeof(dmtherm_flash_phase_t));
    }

    for (size_t i = 0; i < fr.phases.size(); ++i) {
        const auto& ph = fr.phases[i];
        auto& dst = out->phases[i];

        dst.phase = toPhase(ph.type);
        dst.fraction = ph.fraction;
        dst.molar_density_mol_per_m3 = ph.density;
        dst.compressibility = ph.compressibility;

        if (!copyArray(ph.x.data(), ph.x.size(), dst.composition)) {
            setLastError(h, "flash: out of memory (phase x)");
            return DMTHERM_STATUS_OUT_OF_MEMORY;
        }
        dst.composition_len = ph.x.size();

        if (!copyArray(ph.phi.data(), ph.phi.size(), dst.fugacity_coeff)) {
            setLastError(h, "flash: out of memory (phase phi)");
            return DMTHERM_STATUS_OUT_OF_MEMORY;
        }
        dst.fugacity_coeff_len = ph.phi.size();

        dst.enthalpy_departure_j_per_mol = opt(ph.H_res);
        dst.entropy_departure_j_per_mol_k = opt(ph.S_res);
        dst.gibbs_departure_j_per_mol = opt(ph.G_res);
        dst.enthalpy_j_per_mol = opt(ph.H);
        dst.entropy_j_per_mol_k = opt(ph.S);
        dst.gibbs_j_per_mol = opt(ph.G);
    }

    return DMTHERM_STATUS_OK;
}

} // namespace

extern "C" {

uint32_t dmtherm_abi_version(void) {
    return DMTHERM_ABI_VERSION;
}

const char* dmtherm_version_string(void) {
    return "DMThermEngine (C ABI)";
}

void dmtherm_free(void* p) {
    std::free(p);
}

dmtherm_status_t dmtherm_convert_pressure(
    double value,
    dmtherm_pressure_unit_t from_unit,
    dmtherm_pressure_unit_t to_unit,
    double* out_value)
{
    if (!out_value) return DMTHERM_STATUS_INVALID_ARGUMENT;
    const auto from = toCorePressureUnit(from_unit);
    const auto to = toCorePressureUnit(to_unit);
    if (!from.has_value() || !to.has_value()) return DMTHERM_STATUS_INVALID_ARGUMENT;

    try {
        *out_value = Core::Units::UnitConverter::convert(value, from.value(), to.value());
        return DMTHERM_STATUS_OK;
    } catch (const std::invalid_argument&) {
        return DMTHERM_STATUS_INVALID_ARGUMENT;
    } catch (...) {
        return DMTHERM_STATUS_RUNTIME_ERROR;
    }
}

void dmtherm_system_config_set_defaults(dmtherm_system_config_t* cfg) {
    if (!cfg) return;
    std::memset(cfg, 0, sizeof(*cfg));
    cfg->struct_size = static_cast<uint32_t>(sizeof(*cfg));
    cfg->flash_method = "Auto";
    cfg->require_ideal_gas_cp = 1;
    cfg->require_critical_properties = 0;
    cfg->property_method = DMTHERM_PROPERTY_METHOD_EOS_ONLY;
    cfg->reserved_u32[0] = DMTHERM_FAILURE_MODE_STRICT;
}

void dmtherm_flash_config_set_defaults(dmtherm_flash_config_t* cfg) {
    if (!cfg) return;
    std::memset(cfg, 0, sizeof(*cfg));
    cfg->struct_size = static_cast<uint32_t>(sizeof(*cfg));
    cfg->max_phases = 2;
    cfg->compute_properties = 0;
}

void dmtherm_flash_hs_config_set_defaults(dmtherm_flash_hs_config_t* cfg) {
    if (!cfg) return;
    std::memset(cfg, 0, sizeof(*cfg));
    cfg->struct_size = static_cast<uint32_t>(sizeof(*cfg));
    cfg->tolerance = 0.0; // use engine defaults unless overridden by caller
    cfg->temp_bracket_low = 0.0;
    cfg->temp_bracket_high = 0.0;
    cfg->max_iterations = 0;
    cfg->use_newton = 1;
}

dmtherm_status_t dmtherm_system_create(
    const dmtherm_system_config_t* cfg,
    dmtherm_system_t** out_system,
    char** out_error_message)
{
    if (out_error_message) *out_error_message = nullptr;
    if (!out_system) return DMTHERM_STATUS_INVALID_ARGUMENT;
    *out_system = nullptr;
    if (!cfg) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (cfg->struct_size < sizeof(dmtherm_system_config_t)) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (!cfg->databanks_dir || !cfg->eos_type) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (!cfg->components || cfg->num_components == 0) return DMTHERM_STATUS_INVALID_ARGUMENT;

    try {
        Engine::ThermoSystemConfig tcfg;
        tcfg.databanks_dir = cfg->databanks_dir;
        tcfg.eos_type = cfg->eos_type;
        tcfg.flash_method = (cfg->flash_method && cfg->flash_method[0] != '\0') ? cfg->flash_method : "Auto";

        tcfg.components.reserve(cfg->num_components);
        for (size_t i = 0; i < cfg->num_components; ++i) {
            const char* c = cfg->components[i];
            if (!c || c[0] == '\0') {
                return DMTHERM_STATUS_INVALID_ARGUMENT;
            }
            tcfg.components.emplace_back(c);
        }

        if (tcfg.eos_type == "PC-SAFT") {
            tcfg.mixture_options.bip_model = Data::BipModel::PCSaft;
        } else if (tcfg.eos_type == "VTPR") {
            tcfg.mixture_options.bip_model = Data::BipModel::VTPR;
        } else {
            tcfg.mixture_options.bip_model = Data::BipModel::PengRobinson;
        }
        tcfg.mixture_options.require_pcsaft_params = (tcfg.eos_type == "PC-SAFT");
        tcfg.mixture_options.require_critical_properties = (cfg->require_critical_properties != 0);
        tcfg.mixture_options.require_ideal_gas_cp = (cfg->require_ideal_gas_cp != 0);
        tcfg.mixture_options.require_vtpr_params = (tcfg.eos_type == "VTPR");
        tcfg.property_method = toPropertyMethod(cfg->property_method);
        {
            const auto failure_mode = toFailureMode(cfg->reserved_u32[0]);
            if (!failure_mode.has_value()) {
                return DMTHERM_STATUS_INVALID_ARGUMENT;
            }
            tcfg.failure_mode = failure_mode.value();
        }

        Data::Diagnostics diag;
        auto sys = Engine::ThermoSystem::build(tcfg, &diag);

        auto* handle = new SystemHandle();
        handle->sys = std::move(sys);

        *out_system = reinterpret_cast<dmtherm_system_t*>(handle);
        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc& e) {
        if (out_error_message) {
            *out_error_message = dupCString(e.what());
        }
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        if (out_error_message) {
            *out_error_message = dupCString(e.what());
        }
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        if (out_error_message) {
            *out_error_message = dupCString("dmtherm_system_create: unknown error");
        }
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

void dmtherm_system_destroy(dmtherm_system_t* sys) {
    if (!sys) return;
    auto* h = reinterpret_cast<SystemHandle*>(sys);
    delete h;
}

const char* dmtherm_system_last_error(dmtherm_system_t* sys) {
    thread_local std::string tls;
    if (!sys) {
        tls = "dmtherm_system_last_error: system is null";
        return tls.c_str();
    }
    auto* h = reinterpret_cast<SystemHandle*>(sys);
    std::lock_guard<std::mutex> lock(h->mu);
    tls = h->last_error;
    return tls.c_str();
}

dmtherm_status_t dmtherm_system_get_component_mws(
    dmtherm_system_t* sys,
    double* out_mw_g_per_mol,
    size_t out_len)
{
    if (!sys || !out_mw_g_per_mol) return DMTHERM_STATUS_INVALID_ARGUMENT;
    auto* h = reinterpret_cast<SystemHandle*>(sys);
    try {
        const int nc = h->sys.mixture().numComponents();
        if (out_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;
        for (int i = 0; i < nc; ++i) {
            out_mw_g_per_mol[static_cast<size_t>(i)] = h->sys.mixture().component(i).MW();
        }
        return DMTHERM_STATUS_OK;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_get_component_mws: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

void dmtherm_tp_phase_eval_destroy(dmtherm_tp_phase_eval_t* out) {
    if (!out) return;
    if (out->ln_fugacity_coeff) dmtherm_free(out->ln_fugacity_coeff);
    if (out->message) dmtherm_free(out->message);
    out->ln_fugacity_coeff = nullptr;
    out->ln_fugacity_coeff_len = 0;
    out->message = nullptr;
}

void dmtherm_tp_phase_eval_batch_destroy(dmtherm_tp_phase_eval_batch_t* out) {
    if (!out) return;
    if (out->status) dmtherm_free(out->status);
    if (out->success) dmtherm_free(out->success);
    if (out->phase) dmtherm_free(out->phase);
    if (out->molar_density_mol_per_m3) dmtherm_free(out->molar_density_mol_per_m3);
    if (out->compressibility) dmtherm_free(out->compressibility);
    if (out->enthalpy_departure_j_per_mol) dmtherm_free(out->enthalpy_departure_j_per_mol);
    if (out->entropy_departure_j_per_mol_k) dmtherm_free(out->entropy_departure_j_per_mol_k);
    if (out->gibbs_departure_j_per_mol) dmtherm_free(out->gibbs_departure_j_per_mol);
    if (out->enthalpy_j_per_mol) dmtherm_free(out->enthalpy_j_per_mol);
    if (out->entropy_j_per_mol_k) dmtherm_free(out->entropy_j_per_mol_k);
    if (out->gibbs_j_per_mol) dmtherm_free(out->gibbs_j_per_mol);
    if (out->ln_fugacity_coeff) dmtherm_free(out->ln_fugacity_coeff);
    if (out->message_offsets) dmtherm_free(out->message_offsets);
    if (out->messages) dmtherm_free(out->messages);
    if (hasTpPhaseEvalBatchHeatCaps(out)) {
        if (out->cp_j_per_mol_k) dmtherm_free(out->cp_j_per_mol_k);
        if (out->cv_j_per_mol_k) dmtherm_free(out->cv_j_per_mol_k);
    }

    out->n_points = 0;
    out->num_components = 0;
    out->status = nullptr;
    out->success = nullptr;
    out->phase = nullptr;
    out->molar_density_mol_per_m3 = nullptr;
    out->compressibility = nullptr;
    out->enthalpy_departure_j_per_mol = nullptr;
    out->entropy_departure_j_per_mol_k = nullptr;
    out->gibbs_departure_j_per_mol = nullptr;
    out->enthalpy_j_per_mol = nullptr;
    out->entropy_j_per_mol_k = nullptr;
    out->gibbs_j_per_mol = nullptr;
    out->ln_fugacity_coeff = nullptr;
    out->ln_fugacity_coeff_len = 0;
    out->ln_fugacity_coeff_stride = 0;
    out->message_offsets = nullptr;
    out->messages = nullptr;
    out->messages_len = 0;
    if (hasTpPhaseEvalBatchHeatCaps(out)) {
        out->cp_j_per_mol_k = nullptr;
        out->cv_j_per_mol_k = nullptr;
    }
}

dmtherm_status_t dmtherm_system_evaluate_tp_phase(
    dmtherm_system_t* sys,
    double T,
    double P,
    const double* x,
    size_t x_len,
    dmtherm_root_selection_t root,
    dmtherm_tp_phase_eval_t* out)
{
    if (!sys || !x || !out) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->struct_size < requiredSizeTpPhaseEvalV1()) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->ln_fugacity_coeff || out->message) return DMTHERM_STATUS_INVALID_ARGUMENT;

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (x_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        std::vector<double> xv(x, x + x_len);
        const auto sel = toRootSelection(root);
        const auto ev = h->sys.evaluateTPPhase(T, P, xv, sel);

        out->temperature_k = ev.T;
        out->pressure_pa_abs = ev.P;
        out->phase = toPhase(ev.phase);
        out->molar_density_mol_per_m3 = ev.rho;
        out->compressibility = ev.Z;
        out->enthalpy_departure_j_per_mol = ev.h_dep;
        out->entropy_departure_j_per_mol_k = ev.s_dep;
        out->gibbs_departure_j_per_mol = ev.g_dep;
        out->enthalpy_j_per_mol = ev.H;
        out->entropy_j_per_mol_k = ev.S;
        out->gibbs_j_per_mol = ev.G;
        if (hasTpPhaseEvalHeatCaps(out)) {
            out->cp_j_per_mol_k = ev.Cp;
            out->cv_j_per_mol_k = ev.Cv;
        }
        out->success = ev.success ? 1 : 0;

        if (!copyArray(ev.ln_phi.data(), ev.ln_phi.size(), out->ln_fugacity_coeff)) {
            setLastError(h, "dmtherm_system_evaluate_tp_phase: out of memory (ln_phi)");
            return DMTHERM_STATUS_OUT_OF_MEMORY;
        }
        out->ln_fugacity_coeff_len = ev.ln_phi.size();

        if (!ev.message.empty()) {
            out->message = dupCString(ev.message);
            if (!out->message) {
                setLastError(h, "dmtherm_system_evaluate_tp_phase: out of memory (message)");
                return DMTHERM_STATUS_OUT_OF_MEMORY;
            }
        }

        if (!ev.success) {
            setLastError(h, ev.message);
        } else {
            setLastError(h, "");
        }

        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_evaluate_tp_phase: out of memory");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_evaluate_tp_phase: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

dmtherm_status_t dmtherm_system_evaluate_tp_phase_batch(
    dmtherm_system_t* sys,
    const double* Ts,
    const double* Ps,
    size_t n_points,
    const double* x,
    size_t x_len,
    dmtherm_root_selection_t root,
    int32_t include_messages,
    dmtherm_tp_phase_eval_batch_t* out)
{
    if (!sys || !Ts || !Ps || !x || !out) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->struct_size < requiredSizeTpPhaseEvalBatchV1()) return DMTHERM_STATUS_INVALID_ARGUMENT;
    const bool include_heat_caps = hasTpPhaseEvalBatchHeatCaps(out);
    if (out->status || out->success || out->phase || out->molar_density_mol_per_m3 || out->compressibility || out->enthalpy_departure_j_per_mol || out->entropy_departure_j_per_mol_k || out->gibbs_departure_j_per_mol || out->enthalpy_j_per_mol || out->entropy_j_per_mol_k || out->gibbs_j_per_mol || out->ln_fugacity_coeff ||
        out->message_offsets || out->messages ||
        (include_heat_caps && (out->cp_j_per_mol_k || out->cv_j_per_mol_k)))
    {
        return DMTHERM_STATUS_INVALID_ARGUMENT;
    }

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (x_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        const size_t ln_phi_stride = x_len;
        if (n_points > 0 && ln_phi_stride > 0) {
            const size_t max = static_cast<size_t>(-1);
            if (n_points > max / ln_phi_stride) return DMTHERM_STATUS_INVALID_ARGUMENT;
        }
        const size_t ln_phi_len = n_points * ln_phi_stride;

        out->n_points = n_points;
        out->num_components = x_len;
        out->ln_fugacity_coeff_len = ln_phi_len;
        out->ln_fugacity_coeff_stride = ln_phi_stride;

        out->status = allocArray<int32_t>(n_points);
        out->success = allocArray<int32_t>(n_points);
        out->phase = allocArray<dmtherm_phase_t>(n_points);
        out->molar_density_mol_per_m3 = allocArray<double>(n_points);
        out->compressibility = allocArray<double>(n_points);
        out->enthalpy_departure_j_per_mol = allocArray<double>(n_points);
        out->entropy_departure_j_per_mol_k = allocArray<double>(n_points);
        out->gibbs_departure_j_per_mol = allocArray<double>(n_points);
        out->enthalpy_j_per_mol = allocArray<double>(n_points);
        out->entropy_j_per_mol_k = allocArray<double>(n_points);
        out->gibbs_j_per_mol = allocArray<double>(n_points);
        out->ln_fugacity_coeff = allocArray<double>(ln_phi_len);
        if (include_heat_caps) {
            out->cp_j_per_mol_k = allocArray<double>(n_points);
            out->cv_j_per_mol_k = allocArray<double>(n_points);
        }

        if ((n_points && (!out->status || !out->success || !out->phase || !out->molar_density_mol_per_m3 || !out->compressibility || !out->enthalpy_departure_j_per_mol || !out->entropy_departure_j_per_mol_k || !out->gibbs_departure_j_per_mol || !out->enthalpy_j_per_mol || !out->entropy_j_per_mol_k ||
                          !out->gibbs_j_per_mol || (include_heat_caps && (!out->cp_j_per_mol_k || !out->cv_j_per_mol_k)))) ||
            (ln_phi_len && !out->ln_fugacity_coeff))
        {
            setLastError(h, "dmtherm_system_evaluate_tp_phase_batch: out of memory (arrays)");
            dmtherm_tp_phase_eval_batch_destroy(out);
            return DMTHERM_STATUS_OUT_OF_MEMORY;
        }

        std::vector<double> xv(x, x + x_len);
        const auto sel = toRootSelection(root);

        std::vector<std::string> msgs;
        if (include_messages) {
            msgs.resize(n_points);
            out->message_offsets = allocArray<size_t>(n_points);
            if (n_points && !out->message_offsets) {
                setLastError(h, "dmtherm_system_evaluate_tp_phase_batch: out of memory (message_offsets)");
                dmtherm_tp_phase_eval_batch_destroy(out);
                return DMTHERM_STATUS_OUT_OF_MEMORY;
            }
        }

        std::string last_err;
        bool any_err = false;

        for (size_t i = 0; i < n_points; ++i) {
            const double T = Ts[i];
            const double P = Ps[i];

            // Default failure fill (matches CLI: zeros/Unknown when not OK).
            out->status[i] = DMTHERM_STATUS_INTERNAL_ERROR;
            out->success[i] = 0;
            out->phase[i] = DMTHERM_PHASE_UNKNOWN;
            out->molar_density_mol_per_m3[i] = 0.0;
            out->compressibility[i] = 0.0;
            out->enthalpy_departure_j_per_mol[i] = 0.0;
            out->entropy_departure_j_per_mol_k[i] = 0.0;
            out->gibbs_departure_j_per_mol[i] = 0.0;
            out->enthalpy_j_per_mol[i] = 0.0;
            out->entropy_j_per_mol_k[i] = 0.0;
            out->gibbs_j_per_mol[i] = 0.0;
            if (include_heat_caps) {
                out->cp_j_per_mol_k[i] = 0.0;
                out->cv_j_per_mol_k[i] = 0.0;
            }
            if (ln_phi_stride) {
                std::fill(out->ln_fugacity_coeff + i * ln_phi_stride, out->ln_fugacity_coeff + (i + 1) * ln_phi_stride, 0.0);
            }

            try {
                const auto ev = h->sys.evaluateTPPhase(T, P, xv, sel);

                out->status[i] = DMTHERM_STATUS_OK;
                out->success[i] = ev.success ? 1 : 0;
                out->phase[i] = toPhase(ev.phase);
                out->molar_density_mol_per_m3[i] = ev.rho;
                out->compressibility[i] = ev.Z;
                out->enthalpy_departure_j_per_mol[i] = ev.h_dep;
                out->entropy_departure_j_per_mol_k[i] = ev.s_dep;
                out->gibbs_departure_j_per_mol[i] = ev.g_dep;
                out->enthalpy_j_per_mol[i] = ev.H;
                out->entropy_j_per_mol_k[i] = ev.S;
                out->gibbs_j_per_mol[i] = ev.G;
                if (include_heat_caps) {
                    out->cp_j_per_mol_k[i] = ev.Cp;
                    out->cv_j_per_mol_k[i] = ev.Cv;
                }

                if (ln_phi_stride && ev.ln_phi.size() == ln_phi_stride) {
                    std::memcpy(out->ln_fugacity_coeff + i * ln_phi_stride, ev.ln_phi.data(), ln_phi_stride * sizeof(double));
                }

                if (include_messages) {
                    msgs[i] = ev.message;
                }

                if (!ev.success) {
                    any_err = true;
                    last_err = ev.message;
                }
            } catch (const std::bad_alloc&) {
                any_err = true;
                out->status[i] = DMTHERM_STATUS_OUT_OF_MEMORY;
                last_err = "out of memory";
                if (include_messages) msgs[i] = last_err;
            } catch (const std::exception& e) {
                any_err = true;
                out->status[i] = DMTHERM_STATUS_RUNTIME_ERROR;
                last_err = e.what();
                if (include_messages) msgs[i] = last_err;
            } catch (...) {
                any_err = true;
                out->status[i] = DMTHERM_STATUS_INTERNAL_ERROR;
                last_err = "unknown error";
                if (include_messages) msgs[i] = last_err;
            }
        }

        if (include_messages) {
            size_t total = 0;
            for (size_t i = 0; i < msgs.size(); ++i) {
                total += msgs[i].size() + 1; // null terminator
            }
            out->messages = allocArray<char>(total);
            if (total && !out->messages) {
                setLastError(h, "dmtherm_system_evaluate_tp_phase_batch: out of memory (messages)");
                dmtherm_tp_phase_eval_batch_destroy(out);
                return DMTHERM_STATUS_OUT_OF_MEMORY;
            }
            out->messages_len = total;

            size_t cursor = 0;
            for (size_t i = 0; i < msgs.size(); ++i) {
                out->message_offsets[i] = cursor;
                const auto& s = msgs[i];
                if (!s.empty()) {
                    std::memcpy(out->messages + cursor, s.data(), s.size());
                }
                cursor += s.size();
                out->messages[cursor++] = '\0';
            }
        }

        if (any_err) {
            setLastError(h, last_err);
        } else {
            setLastError(h, "");
        }

        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_evaluate_tp_phase_batch: out of memory");
        dmtherm_tp_phase_eval_batch_destroy(out);
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        dmtherm_tp_phase_eval_batch_destroy(out);
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_evaluate_tp_phase_batch: unknown error");
        dmtherm_tp_phase_eval_batch_destroy(out);
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

dmtherm_status_t dmtherm_system_evaluate_tp_phase_batch_into(
    dmtherm_system_t* sys,
    const double* Ts,
    const double* Ps,
    size_t n_points,
    const double* x,
    size_t x_len,
    dmtherm_root_selection_t root,
    int32_t include_messages,
    dmtherm_tp_phase_eval_batch_buffers_t* io)
{
    if (!sys || !Ts || !Ps || !x || !io) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (io->struct_size < requiredSizeTpPhaseEvalBatchBuffersV1()) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (io->n_points != n_points) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (io->num_components != x_len) return DMTHERM_STATUS_INVALID_ARGUMENT;

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (x_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        const bool include_heat_caps = hasTpPhaseEvalBatchBuffersHeatCaps(io);

        if (n_points != 0) {
            if (!io->status || !io->success || !io->phase) return DMTHERM_STATUS_INVALID_ARGUMENT;
            if (!io->molar_density_mol_per_m3 || !io->compressibility) return DMTHERM_STATUS_INVALID_ARGUMENT;
            if (!io->enthalpy_departure_j_per_mol || !io->entropy_departure_j_per_mol_k || !io->gibbs_departure_j_per_mol) return DMTHERM_STATUS_INVALID_ARGUMENT;
            if (!io->enthalpy_j_per_mol || !io->entropy_j_per_mol_k || !io->gibbs_j_per_mol) return DMTHERM_STATUS_INVALID_ARGUMENT;
            if (!io->ln_fugacity_coeff) return DMTHERM_STATUS_INVALID_ARGUMENT;
        }

        if (io->ln_fugacity_coeff_stride < x_len) return DMTHERM_STATUS_INVALID_ARGUMENT;
        const size_t required_lnphi = n_points * io->ln_fugacity_coeff_stride;
        if (required_lnphi > io->ln_fugacity_coeff_len) return DMTHERM_STATUS_INVALID_ARGUMENT;

        io->messages_len = 0;
        io->messages_truncated = 0;

        if (include_messages) {
            if (n_points != 0 && (!io->message_offsets || !io->messages)) return DMTHERM_STATUS_INVALID_ARGUMENT;
            if (io->messages_capacity < n_points) return DMTHERM_STATUS_INVALID_ARGUMENT;
        }

        std::vector<double> xv(x, x + x_len);
        const auto sel = toRootSelection(root);

        std::string last_err;
        bool any_err = false;

        size_t msg_cursor = 0;
        for (size_t i = 0; i < n_points; ++i) {
            const double T = Ts[i];
            const double P = Ps[i];

            io->status[i] = DMTHERM_STATUS_INTERNAL_ERROR;
            io->success[i] = 0;
            io->phase[i] = DMTHERM_PHASE_UNKNOWN;
            io->molar_density_mol_per_m3[i] = 0.0;
            io->compressibility[i] = 0.0;
            io->enthalpy_departure_j_per_mol[i] = 0.0;
            io->entropy_departure_j_per_mol_k[i] = 0.0;
            io->gibbs_departure_j_per_mol[i] = 0.0;
            io->enthalpy_j_per_mol[i] = 0.0;
            io->entropy_j_per_mol_k[i] = 0.0;
            io->gibbs_j_per_mol[i] = 0.0;
            if (include_heat_caps) {
                if (io->cp_j_per_mol_k) io->cp_j_per_mol_k[i] = 0.0;
                if (io->cv_j_per_mol_k) io->cv_j_per_mol_k[i] = 0.0;
            }
            if (io->ln_fugacity_coeff_stride) {
                std::fill(io->ln_fugacity_coeff + i * io->ln_fugacity_coeff_stride, io->ln_fugacity_coeff + (i + 1) * io->ln_fugacity_coeff_stride, 0.0);
            }

            std::string msg;
            try {
                const auto ev = h->sys.evaluateTPPhase(T, P, xv, sel);

                io->status[i] = DMTHERM_STATUS_OK;
                io->success[i] = ev.success ? 1 : 0;
                io->phase[i] = toPhase(ev.phase);
                io->molar_density_mol_per_m3[i] = ev.rho;
                io->compressibility[i] = ev.Z;
                io->enthalpy_departure_j_per_mol[i] = ev.h_dep;
                io->entropy_departure_j_per_mol_k[i] = ev.s_dep;
                io->gibbs_departure_j_per_mol[i] = ev.g_dep;
                io->enthalpy_j_per_mol[i] = ev.H;
                io->entropy_j_per_mol_k[i] = ev.S;
                io->gibbs_j_per_mol[i] = ev.G;
                if (include_heat_caps) {
                    if (io->cp_j_per_mol_k) io->cp_j_per_mol_k[i] = ev.Cp;
                    if (io->cv_j_per_mol_k) io->cv_j_per_mol_k[i] = ev.Cv;
                }

                const size_t ln_stride = io->ln_fugacity_coeff_stride;
                if (ln_stride && ev.ln_phi.size() == x_len) {
                    std::memcpy(io->ln_fugacity_coeff + i * ln_stride, ev.ln_phi.data(), x_len * sizeof(double));
                }

                msg = ev.message;

                if (!ev.success) {
                    any_err = true;
                    last_err = ev.message;
                }
            } catch (const std::bad_alloc&) {
                any_err = true;
                io->status[i] = DMTHERM_STATUS_OUT_OF_MEMORY;
                last_err = "out of memory";
                msg = last_err;
            } catch (const std::exception& e) {
                any_err = true;
                io->status[i] = DMTHERM_STATUS_RUNTIME_ERROR;
                last_err = e.what();
                msg = last_err;
            } catch (...) {
                any_err = true;
                io->status[i] = DMTHERM_STATUS_INTERNAL_ERROR;
                last_err = "unknown error";
                msg = last_err;
            }

            if (include_messages) {
                if (msg_cursor >= io->messages_capacity) {
                    io->messages_truncated = 1;
                    const size_t off = (io->messages_capacity != 0) ? (io->messages_capacity - 1) : 0;
                    io->message_offsets[i] = off;
                    if (io->messages_capacity != 0) {
                        io->messages[off] = '\0';
                    }
                } else {
                    const size_t off = msg_cursor;
                    io->message_offsets[i] = off;
                    const size_t remaining = io->messages_capacity - off;
                    const size_t copy_n = (msg.size() < (remaining - 1)) ? msg.size() : (remaining - 1);
                    if (copy_n != 0) {
                        std::memcpy(io->messages + off, msg.data(), copy_n);
                    }
                    io->messages[off + copy_n] = '\0';
                    msg_cursor += (copy_n + 1);
                    if (copy_n < msg.size()) io->messages_truncated = 1;
                }
            }
        }

        if (include_messages) {
            io->messages_len = (msg_cursor <= io->messages_capacity) ? msg_cursor : io->messages_capacity;
        } else {
            io->messages_len = 0;
        }

        if (any_err) {
            setLastError(h, last_err);
        } else {
            setLastError(h, "");
        }

        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_evaluate_tp_phase_batch_into: out of memory");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_evaluate_tp_phase_batch_into: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

void dmtherm_tp_departure_roots_destroy(dmtherm_tp_departure_roots_t* out) {
    if (!out) return;
    if (out->roots) {
        for (size_t i = 0; i < out->roots_len; ++i) {
            if (out->roots[i].ln_fugacity_coeff) dmtherm_free(out->roots[i].ln_fugacity_coeff);
            out->roots[i].ln_fugacity_coeff = nullptr;
            out->roots[i].ln_fugacity_coeff_len = 0;
        }
        dmtherm_free(out->roots);
    }
    if (out->message) dmtherm_free(out->message);
    out->roots = nullptr;
    out->roots_len = 0;
    out->message = nullptr;
}

void dmtherm_critical_estimate_destroy(dmtherm_critical_estimate_t* out) {
    if (!out) return;
    if (out->method) dmtherm_free(out->method);
    if (out->message) dmtherm_free(out->message);
    out->method = nullptr;
    out->message = nullptr;
}

dmtherm_status_t dmtherm_system_departure_roots_tp(
    dmtherm_system_t* sys,
    double T,
    double P,
    const double* x,
    size_t x_len,
    dmtherm_tp_departure_roots_t* out)
{
    if (!sys || !x || !out) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->struct_size < sizeof(dmtherm_tp_departure_roots_t)) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->roots || out->message) return DMTHERM_STATUS_INVALID_ARGUMENT;

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (x_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        std::vector<double> xv(x, x + x_len);
        const auto dep = h->sys.departureRootsTP(T, P, xv);

        out->temperature_k = dep.temperature;
        out->pressure_pa_abs = dep.pressure;
        out->stable_root_index = dep.stable_root_index.has_value() ? dep.stable_root_index.value() : -1;

        const auto ig = Equilibrium::Contributions::idealGasProps(h->sys.mixture(), T, P, xv);

        out->roots_len = dep.roots.size();
        out->roots = allocArray<dmtherm_tp_departure_root_t>(out->roots_len);
        if (out->roots_len != 0 && !out->roots) {
            setLastError(h, "dmtherm_system_departure_roots_tp: out of memory (roots)");
            return DMTHERM_STATUS_OUT_OF_MEMORY;
        }
        if (out->roots_len != 0) {
            std::memset(out->roots, 0, out->roots_len * sizeof(dmtherm_tp_departure_root_t));
        }

        for (size_t i = 0; i < dep.roots.size(); ++i) {
            const auto& r = dep.roots[i];
            auto& dst = out->roots[i];
            dst.phase = toPhase(r.phase);
            dst.mechanically_stable = r.mechanically_stable ? 1 : 0;
            dst.molar_density_mol_per_m3 = r.density;
            dst.compressibility = r.Z;
            dst.enthalpy_departure_j_per_mol = r.h_dep;
            dst.entropy_departure_j_per_mol_k = r.s_dep;
            dst.gibbs_departure_j_per_mol = r.g_dep;
            dst.enthalpy_j_per_mol = ig.h_ig + r.h_dep;
            dst.entropy_j_per_mol_k = ig.s_ig + r.s_dep;
            dst.gibbs_j_per_mol = dst.enthalpy_j_per_mol - T * dst.entropy_j_per_mol_k;

            if (!copyArray(r.ln_phi.data(), r.ln_phi.size(), dst.ln_fugacity_coeff)) {
                setLastError(h, "dmtherm_system_departure_roots_tp: out of memory (ln_phi)");
                return DMTHERM_STATUS_OUT_OF_MEMORY;
            }
            dst.ln_fugacity_coeff_len = r.ln_phi.size();
        }

        setLastError(h, "");
        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_departure_roots_tp: out of memory");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_departure_roots_tp: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

dmtherm_status_t dmtherm_system_estimate_critical(
    dmtherm_system_t* sys,
    int32_t component_index,
    dmtherm_critical_estimate_t* out)
{
    if (!sys || !out) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->struct_size < requiredSizeCriticalEstimateV1()) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->method || out->message) return DMTHERM_STATUS_INVALID_ARGUMENT;

    auto* h = reinterpret_cast<SystemHandle*>(sys);
    try {
        const auto result = h->sys.estimateCritical(static_cast<int>(component_index));

        out->component_index = component_index;
        out->converged = result.converged ? 1 : 0;
        out->iterations = result.iterations;
        out->critical_temperature_k = result.Tc;
        out->critical_pressure_pa_abs = result.Pc;
        out->critical_density_mol_per_m3 = result.rho;
        out->critical_compressibility = result.Zc;
        const Core::Units::PropertyVariable critical_volume = result.criticalVolumeValue();
        out->critical_volume_m3_per_kmol = critical_volume.as(Core::Units::Unit::M3_PER_KMOL);

        if (!result.method.empty()) {
            out->method = dupCString(result.method);
            if (!out->method) {
                setLastError(h, "dmtherm_system_estimate_critical: out of memory (method)");
                return DMTHERM_STATUS_OUT_OF_MEMORY;
            }
        }
        if (!result.message.empty()) {
            out->message = dupCString(result.message);
            if (!out->message) {
                dmtherm_critical_estimate_destroy(out);
                setLastError(h, "dmtherm_system_estimate_critical: out of memory (message)");
                return DMTHERM_STATUS_OUT_OF_MEMORY;
            }
        }

        if (!result.converged) {
            setLastError(h, result.message);
        } else {
            setLastError(h, "");
        }
        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_estimate_critical: out of memory");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_estimate_critical: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

void dmtherm_flash_pt_result_destroy(dmtherm_flash_pt_result_t* out) {
    if (!out) return;
    if (out->phases) {
        for (size_t i = 0; i < out->phases_len; ++i) {
            auto& ph = out->phases[i];
            if (ph.composition) dmtherm_free(ph.composition);
            if (ph.fugacity_coeff) dmtherm_free(ph.fugacity_coeff);
            ph.composition = nullptr;
            ph.composition_len = 0;
            ph.fugacity_coeff = nullptr;
            ph.fugacity_coeff_len = 0;
        }
        dmtherm_free(out->phases);
    }
    if (out->message) dmtherm_free(out->message);
    if (out->method_used) dmtherm_free(out->method_used);
    out->phases = nullptr;
    out->phases_len = 0;
    out->message = nullptr;
    out->method_used = nullptr;
}

void dmtherm_flash_pt_batch_destroy(dmtherm_flash_pt_batch_t* out) {
    if (!out) return;

    if (out->status) dmtherm_free(out->status);
    if (out->converged) dmtherm_free(out->converged);
    if (out->iterations) dmtherm_free(out->iterations);
    if (out->residual) dmtherm_free(out->residual);
    if (out->temperature_k) dmtherm_free(out->temperature_k);
    if (out->pressure_pa_abs) dmtherm_free(out->pressure_pa_abs);
    if (out->vapor_fraction) dmtherm_free(out->vapor_fraction);
    if (out->num_phases) dmtherm_free(out->num_phases);
    if (out->mixture_enthalpy_j_per_mol) dmtherm_free(out->mixture_enthalpy_j_per_mol);
    if (out->mixture_entropy_j_per_mol_k) dmtherm_free(out->mixture_entropy_j_per_mol_k);
    if (out->mixture_gibbs_j_per_mol) dmtherm_free(out->mixture_gibbs_j_per_mol);

    if (out->phase) dmtherm_free(out->phase);
    if (out->phase_fraction) dmtherm_free(out->phase_fraction);
    if (out->phase_molar_density_mol_per_m3) dmtherm_free(out->phase_molar_density_mol_per_m3);
    if (out->phase_compressibility) dmtherm_free(out->phase_compressibility);
    if (out->phase_composition) dmtherm_free(out->phase_composition);
    if (out->phase_fugacity_coeff_available) dmtherm_free(out->phase_fugacity_coeff_available);
    if (out->phase_fugacity_coeff) dmtherm_free(out->phase_fugacity_coeff);

    if (out->phase_enthalpy_departure_j_per_mol) dmtherm_free(out->phase_enthalpy_departure_j_per_mol);
    if (out->phase_entropy_departure_j_per_mol_k) dmtherm_free(out->phase_entropy_departure_j_per_mol_k);
    if (out->phase_gibbs_departure_j_per_mol) dmtherm_free(out->phase_gibbs_departure_j_per_mol);
    if (out->phase_enthalpy_j_per_mol) dmtherm_free(out->phase_enthalpy_j_per_mol);
    if (out->phase_entropy_j_per_mol_k) dmtherm_free(out->phase_entropy_j_per_mol_k);
    if (out->phase_gibbs_j_per_mol) dmtherm_free(out->phase_gibbs_j_per_mol);

    if (out->message_offsets) dmtherm_free(out->message_offsets);
    if (out->messages) dmtherm_free(out->messages);

    out->n_points = 0;
    out->num_components = 0;
    out->max_phases = 0;
    out->status = nullptr;
    out->converged = nullptr;
    out->iterations = nullptr;
    out->residual = nullptr;
    out->temperature_k = nullptr;
    out->pressure_pa_abs = nullptr;
    out->vapor_fraction = nullptr;
    out->num_phases = nullptr;
    out->mixture_enthalpy_j_per_mol = nullptr;
    out->mixture_entropy_j_per_mol_k = nullptr;
    out->mixture_gibbs_j_per_mol = nullptr;

    out->phases_len = 0;
    out->phase = nullptr;
    out->phase_fraction = nullptr;
    out->phase_molar_density_mol_per_m3 = nullptr;
    out->phase_compressibility = nullptr;
    out->phase_composition_len = 0;
    out->phase_composition = nullptr;
    out->phase_fugacity_coeff_available = nullptr;
    out->phase_fugacity_coeff_len = 0;
    out->phase_fugacity_coeff = nullptr;

    out->phase_enthalpy_departure_j_per_mol = nullptr;
    out->phase_entropy_departure_j_per_mol_k = nullptr;
    out->phase_gibbs_departure_j_per_mol = nullptr;
    out->phase_enthalpy_j_per_mol = nullptr;
    out->phase_entropy_j_per_mol_k = nullptr;
    out->phase_gibbs_j_per_mol = nullptr;

    out->message_offsets = nullptr;
    out->messages = nullptr;
    out->messages_len = 0;
}

dmtherm_status_t dmtherm_system_flash_pt(
    dmtherm_system_t* sys,
    double T,
    double P,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    dmtherm_flash_pt_result_t* out)
{
    if (!sys || !z || !out) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->struct_size < sizeof(dmtherm_flash_pt_result_t)) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->phases || out->message || out->method_used) return DMTHERM_STATUS_INVALID_ARGUMENT;

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (z_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        Config::FlashConfig fcfg = Config::FlashConfig::defaults();
        if (cfg && cfg->struct_size >= sizeof(dmtherm_flash_config_t)) {
            if (cfg->max_phases > 0) fcfg.max_phases = cfg->max_phases;
            fcfg.compute_properties = (cfg->compute_properties != 0);
        }

        std::vector<double> zv(z, z + z_len);
        const auto fr = h->sys.flashPT(T, P, zv, fcfg);
        const auto st = fillFlashResult(h, fr, out);
        if (st != DMTHERM_STATUS_OK) return st;

        setLastError(h, "");
        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_flash_pt: out of memory");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_flash_pt: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

dmtherm_status_t dmtherm_system_flash_pt_batch(
    dmtherm_system_t* sys,
    const double* Ts,
    const double* Ps,
    size_t n_points,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    int32_t include_messages,
    dmtherm_flash_pt_batch_t* out)
{
    if (!sys || !Ts || !Ps || !z || !out) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->struct_size < sizeof(dmtherm_flash_pt_batch_t)) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->status || out->converged || out->iterations || out->residual || out->temperature_k || out->pressure_pa_abs || out->vapor_fraction || out->num_phases || out->mixture_enthalpy_j_per_mol || out->mixture_entropy_j_per_mol_k ||
        out->mixture_gibbs_j_per_mol || out->phase || out->phase_fraction || out->phase_molar_density_mol_per_m3 || out->phase_compressibility || out->phase_composition || out->phase_fugacity_coeff_available || out->phase_fugacity_coeff ||
        out->phase_enthalpy_departure_j_per_mol || out->phase_entropy_departure_j_per_mol_k || out->phase_gibbs_departure_j_per_mol || out->phase_enthalpy_j_per_mol || out->phase_entropy_j_per_mol_k || out->phase_gibbs_j_per_mol || out->message_offsets || out->messages)
    {
        return DMTHERM_STATUS_INVALID_ARGUMENT;
    }

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (z_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        Config::FlashConfig fcfg = Config::FlashConfig::defaults();
        if (cfg && cfg->struct_size >= sizeof(dmtherm_flash_config_t)) {
            if (cfg->max_phases > 0) fcfg.max_phases = cfg->max_phases;
            fcfg.compute_properties = (cfg->compute_properties != 0);
        }

        const int32_t max_phases = static_cast<int32_t>(fcfg.max_phases);
        if (max_phases <= 0) return DMTHERM_STATUS_INVALID_ARGUMENT;

        const size_t max = static_cast<size_t>(-1);
        if (n_points > 0 && static_cast<size_t>(max_phases) > 0 && n_points > max / static_cast<size_t>(max_phases)) {
            return DMTHERM_STATUS_INVALID_ARGUMENT;
        }
        const size_t phases_len = n_points * static_cast<size_t>(max_phases);
        if (phases_len > 0 && z_len > 0 && phases_len > max / z_len) return DMTHERM_STATUS_INVALID_ARGUMENT;
        const size_t phase_vec_len = phases_len * z_len;

        out->n_points = n_points;
        out->num_components = z_len;
        out->max_phases = max_phases;
        out->phases_len = phases_len;
        out->phase_composition_len = phase_vec_len;
        out->phase_fugacity_coeff_len = phase_vec_len;
        out->messages_len = 0;

        out->status = allocArray<int32_t>(n_points);
        out->converged = allocArray<int32_t>(n_points);
        out->iterations = allocArray<int32_t>(n_points);
        out->residual = allocArray<double>(n_points);
        out->temperature_k = allocArray<double>(n_points);
        out->pressure_pa_abs = allocArray<double>(n_points);
        out->vapor_fraction = allocArray<double>(n_points);
        out->num_phases = allocArray<int32_t>(n_points);
        out->mixture_enthalpy_j_per_mol = allocArray<dmtherm_optional_double_t>(n_points);
        out->mixture_entropy_j_per_mol_k = allocArray<dmtherm_optional_double_t>(n_points);
        out->mixture_gibbs_j_per_mol = allocArray<dmtherm_optional_double_t>(n_points);

        out->phase = allocArray<dmtherm_phase_t>(phases_len);
        out->phase_fraction = allocArray<double>(phases_len);
        out->phase_molar_density_mol_per_m3 = allocArray<double>(phases_len);
        out->phase_compressibility = allocArray<double>(phases_len);
        out->phase_composition = allocArray<double>(phase_vec_len);
        out->phase_fugacity_coeff_available = allocArray<int32_t>(phases_len);
        out->phase_fugacity_coeff = allocArray<double>(phase_vec_len);

        out->phase_enthalpy_departure_j_per_mol = allocArray<dmtherm_optional_double_t>(phases_len);
        out->phase_entropy_departure_j_per_mol_k = allocArray<dmtherm_optional_double_t>(phases_len);
        out->phase_gibbs_departure_j_per_mol = allocArray<dmtherm_optional_double_t>(phases_len);
        out->phase_enthalpy_j_per_mol = allocArray<dmtherm_optional_double_t>(phases_len);
        out->phase_entropy_j_per_mol_k = allocArray<dmtherm_optional_double_t>(phases_len);
        out->phase_gibbs_j_per_mol = allocArray<dmtherm_optional_double_t>(phases_len);

        if ((n_points && (!out->status || !out->converged || !out->iterations || !out->residual || !out->temperature_k || !out->pressure_pa_abs || !out->vapor_fraction || !out->num_phases ||
                          !out->mixture_enthalpy_j_per_mol || !out->mixture_entropy_j_per_mol_k || !out->mixture_gibbs_j_per_mol)) ||
            (phases_len && (!out->phase || !out->phase_fraction || !out->phase_molar_density_mol_per_m3 || !out->phase_compressibility || !out->phase_fugacity_coeff_available || !out->phase_enthalpy_departure_j_per_mol || !out->phase_entropy_departure_j_per_mol_k ||
                            !out->phase_gibbs_departure_j_per_mol || !out->phase_enthalpy_j_per_mol || !out->phase_entropy_j_per_mol_k || !out->phase_gibbs_j_per_mol)) ||
            (phase_vec_len && (!out->phase_composition || !out->phase_fugacity_coeff)))
        {
            setLastError(h, "dmtherm_system_flash_pt_batch: out of memory (arrays)");
            dmtherm_flash_pt_batch_destroy(out);
            return DMTHERM_STATUS_OUT_OF_MEMORY;
        }

        std::vector<double> zv(z, z + z_len);

        std::vector<std::string> msgs;
        if (include_messages) {
            msgs.resize(n_points);
            out->message_offsets = allocArray<size_t>(n_points);
            if (n_points && !out->message_offsets) {
                setLastError(h, "dmtherm_system_flash_pt_batch: out of memory (message_offsets)");
                dmtherm_flash_pt_batch_destroy(out);
                return DMTHERM_STATUS_OUT_OF_MEMORY;
            }
        }

        std::string last_err;
        bool any_err = false;

        for (size_t i = 0; i < n_points; ++i) {
            const double T_in = Ts[i];
            const double P_in = Ps[i];

            out->status[i] = DMTHERM_STATUS_INTERNAL_ERROR;
            out->converged[i] = 0;
            out->iterations[i] = 0;
            out->residual[i] = 0.0;
            out->temperature_k[i] = T_in;
            out->pressure_pa_abs[i] = P_in;
            out->vapor_fraction[i] = 0.0;
            out->num_phases[i] = 0;
            out->mixture_enthalpy_j_per_mol[i] = dmtherm_optional_double_t{};
            out->mixture_entropy_j_per_mol_k[i] = dmtherm_optional_double_t{};
            out->mixture_gibbs_j_per_mol[i] = dmtherm_optional_double_t{};

            // Default phase fill (fixed stride).
            for (int32_t j = 0; j < max_phases; ++j) {
                const size_t idx = i * static_cast<size_t>(max_phases) + static_cast<size_t>(j);
                out->phase[idx] = DMTHERM_PHASE_UNKNOWN;
                out->phase_fraction[idx] = 0.0;
                out->phase_molar_density_mol_per_m3[idx] = 0.0;
                out->phase_compressibility[idx] = 0.0;
                out->phase_fugacity_coeff_available[idx] = 0;
                out->phase_enthalpy_departure_j_per_mol[idx] = dmtherm_optional_double_t{};
                out->phase_entropy_departure_j_per_mol_k[idx] = dmtherm_optional_double_t{};
                out->phase_gibbs_departure_j_per_mol[idx] = dmtherm_optional_double_t{};
                out->phase_enthalpy_j_per_mol[idx] = dmtherm_optional_double_t{};
                out->phase_entropy_j_per_mol_k[idx] = dmtherm_optional_double_t{};
                out->phase_gibbs_j_per_mol[idx] = dmtherm_optional_double_t{};

                if (z_len != 0) {
                    std::fill(out->phase_composition + idx * z_len, out->phase_composition + (idx + 1) * z_len, 0.0);
                    std::fill(out->phase_fugacity_coeff + idx * z_len, out->phase_fugacity_coeff + (idx + 1) * z_len, 0.0);
                }
            }

            try {
                const auto fr = h->sys.flashPT(T_in, P_in, zv, fcfg);

                out->status[i] = DMTHERM_STATUS_OK;
                out->converged[i] = fr.converged ? 1 : 0;
                out->iterations[i] = fr.iterations;
                out->residual[i] = fr.residual;
                out->temperature_k[i] = fr.temperature;
                out->pressure_pa_abs[i] = fr.pressure;
                out->vapor_fraction[i] = fr.vapor_fraction;
                out->num_phases[i] = fr.num_phases;
                out->mixture_enthalpy_j_per_mol[i] = opt(fr.enthalpy);
                out->mixture_entropy_j_per_mol_k[i] = opt(fr.entropy);
                out->mixture_gibbs_j_per_mol[i] = opt(fr.gibbs_energy);

                const size_t nph = std::min(fr.phases.size(), static_cast<size_t>(max_phases));
                for (size_t j = 0; j < nph; ++j) {
                    const auto& ph = fr.phases[j];
                    const size_t idx = i * static_cast<size_t>(max_phases) + j;

                    out->phase[idx] = toPhase(ph.type);
                    out->phase_fraction[idx] = ph.fraction;
                    out->phase_molar_density_mol_per_m3[idx] = ph.density;
                    out->phase_compressibility[idx] = ph.compressibility;

                    // x (always expected to be size num_components).
                    if (z_len != 0) {
                        const size_t copy_n = std::min(z_len, ph.x.size());
                        if (copy_n) {
                            std::memcpy(out->phase_composition + idx * z_len, ph.x.data(), copy_n * sizeof(double));
                        }
                        if (copy_n < z_len) {
                            std::fill(out->phase_composition + idx * z_len + copy_n, out->phase_composition + (idx + 1) * z_len, 0.0);
                        }
                    }

                    // phi (may be unavailable depending on property method).
                    if (ph.phi.size() == z_len && z_len != 0) {
                        std::memcpy(out->phase_fugacity_coeff + idx * z_len, ph.phi.data(), z_len * sizeof(double));
                        out->phase_fugacity_coeff_available[idx] = 1;
                    } else {
                        out->phase_fugacity_coeff_available[idx] = 0;
                        if (z_len != 0) {
                            std::fill(out->phase_fugacity_coeff + idx * z_len, out->phase_fugacity_coeff + (idx + 1) * z_len, 0.0);
                        }
                    }

                    out->phase_enthalpy_departure_j_per_mol[idx] = opt(ph.H_res);
                    out->phase_entropy_departure_j_per_mol_k[idx] = opt(ph.S_res);
                    out->phase_gibbs_departure_j_per_mol[idx] = opt(ph.G_res);
                    out->phase_enthalpy_j_per_mol[idx] = opt(ph.H);
                    out->phase_entropy_j_per_mol_k[idx] = opt(ph.S);
                    out->phase_gibbs_j_per_mol[idx] = opt(ph.G);
                }

                if (include_messages) {
                    msgs[i] = fr.message;
                }
            } catch (const std::bad_alloc&) {
                any_err = true;
                out->status[i] = DMTHERM_STATUS_OUT_OF_MEMORY;
                last_err = "out of memory";
                if (include_messages) msgs[i] = last_err;
            } catch (const std::exception& e) {
                any_err = true;
                out->status[i] = DMTHERM_STATUS_RUNTIME_ERROR;
                last_err = e.what();
                if (include_messages) msgs[i] = last_err;
            } catch (...) {
                any_err = true;
                out->status[i] = DMTHERM_STATUS_INTERNAL_ERROR;
                last_err = "unknown error";
                if (include_messages) msgs[i] = last_err;
            }
        }

        if (include_messages) {
            size_t total = 0;
            for (size_t i = 0; i < msgs.size(); ++i) {
                total += msgs[i].size() + 1;
            }
            out->messages = allocArray<char>(total);
            if (total && !out->messages) {
                setLastError(h, "dmtherm_system_flash_pt_batch: out of memory (messages)");
                dmtherm_flash_pt_batch_destroy(out);
                return DMTHERM_STATUS_OUT_OF_MEMORY;
            }
            out->messages_len = total;

            size_t cursor = 0;
            for (size_t i = 0; i < msgs.size(); ++i) {
                out->message_offsets[i] = cursor;
                const auto& s = msgs[i];
                if (!s.empty()) {
                    std::memcpy(out->messages + cursor, s.data(), s.size());
                }
                cursor += s.size();
                out->messages[cursor++] = '\0';
            }
        }

        if (any_err) {
            setLastError(h, last_err);
        } else {
            setLastError(h, "");
        }
        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_flash_pt_batch: out of memory");
        dmtherm_flash_pt_batch_destroy(out);
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        dmtherm_flash_pt_batch_destroy(out);
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_flash_pt_batch: unknown error");
        dmtherm_flash_pt_batch_destroy(out);
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

dmtherm_status_t dmtherm_system_flash_pt_batch_into(
    dmtherm_system_t* sys,
    const double* Ts,
    const double* Ps,
    size_t n_points,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    int32_t include_messages,
    dmtherm_flash_pt_batch_buffers_t* io)
{
    if (!sys || !Ts || !Ps || !z || !io) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (io->struct_size < sizeof(dmtherm_flash_pt_batch_buffers_t)) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (io->n_points != n_points) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (io->num_components != z_len) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (io->max_phases <= 0) return DMTHERM_STATUS_INVALID_ARGUMENT;

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (z_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        Config::FlashConfig fcfg = Config::FlashConfig::defaults();
        if (cfg && cfg->struct_size >= sizeof(dmtherm_flash_config_t)) {
            if (cfg->max_phases > 0) fcfg.max_phases = cfg->max_phases;
            fcfg.compute_properties = (cfg->compute_properties != 0);
        }

        const int32_t max_phases = io->max_phases;
        if (max_phases <= 0) return DMTHERM_STATUS_INVALID_ARGUMENT;
        if (static_cast<size_t>(max_phases) < fcfg.max_phases) return DMTHERM_STATUS_INVALID_ARGUMENT;

        const size_t max = static_cast<size_t>(-1);
        const size_t stride = static_cast<size_t>(max_phases);
        if (n_points > 0 && stride > 0 && n_points > max / stride) return DMTHERM_STATUS_INVALID_ARGUMENT;
        const size_t required_phases_len = n_points * stride;
        if (required_phases_len > io->phases_len) return DMTHERM_STATUS_INVALID_ARGUMENT;
        if (required_phases_len > 0 && z_len > 0 && required_phases_len > max / z_len) return DMTHERM_STATUS_INVALID_ARGUMENT;
        const size_t required_phase_vec_len = required_phases_len * z_len;
        if (required_phase_vec_len > io->phase_composition_len) return DMTHERM_STATUS_INVALID_ARGUMENT;
        if (required_phase_vec_len > io->phase_fugacity_coeff_len) return DMTHERM_STATUS_INVALID_ARGUMENT;

        if (n_points != 0) {
            if (!io->status || !io->converged || !io->iterations || !io->residual || !io->temperature_k || !io->pressure_pa_abs || !io->vapor_fraction || !io->num_phases) {
                return DMTHERM_STATUS_INVALID_ARGUMENT;
            }
            if (!io->mixture_enthalpy_j_per_mol || !io->mixture_entropy_j_per_mol_k || !io->mixture_gibbs_j_per_mol) return DMTHERM_STATUS_INVALID_ARGUMENT;
        }

        if (required_phases_len != 0) {
            if (!io->phase || !io->phase_fraction || !io->phase_molar_density_mol_per_m3 || !io->phase_compressibility || !io->phase_fugacity_coeff_available) return DMTHERM_STATUS_INVALID_ARGUMENT;
            if (!io->phase_enthalpy_departure_j_per_mol || !io->phase_entropy_departure_j_per_mol_k || !io->phase_gibbs_departure_j_per_mol || !io->phase_enthalpy_j_per_mol || !io->phase_entropy_j_per_mol_k || !io->phase_gibbs_j_per_mol) {
                return DMTHERM_STATUS_INVALID_ARGUMENT;
            }
        }

        if (required_phase_vec_len != 0) {
            if (!io->phase_composition || !io->phase_fugacity_coeff) return DMTHERM_STATUS_INVALID_ARGUMENT;
        }

        io->messages_len = 0;
        io->messages_truncated = 0;

        if (include_messages) {
            if (n_points != 0 && (!io->message_offsets || !io->messages)) return DMTHERM_STATUS_INVALID_ARGUMENT;
            if (io->messages_capacity < n_points) return DMTHERM_STATUS_INVALID_ARGUMENT;
        }

        std::vector<double> zv(z, z + z_len);

        std::string last_err;
        bool any_err = false;
        size_t msg_cursor = 0;

        for (size_t i = 0; i < n_points; ++i) {
            const double T_in = Ts[i];
            const double P_in = Ps[i];

            io->status[i] = DMTHERM_STATUS_INTERNAL_ERROR;
            io->converged[i] = 0;
            io->iterations[i] = 0;
            io->residual[i] = 0.0;
            io->temperature_k[i] = T_in;
            io->pressure_pa_abs[i] = P_in;
            io->vapor_fraction[i] = 0.0;
            io->num_phases[i] = 0;
            io->mixture_enthalpy_j_per_mol[i] = dmtherm_optional_double_t{};
            io->mixture_entropy_j_per_mol_k[i] = dmtherm_optional_double_t{};
            io->mixture_gibbs_j_per_mol[i] = dmtherm_optional_double_t{};

            // Default phase fill (fixed stride).
            for (int32_t j = 0; j < max_phases; ++j) {
                const size_t idx = i * stride + static_cast<size_t>(j);
                io->phase[idx] = DMTHERM_PHASE_UNKNOWN;
                io->phase_fraction[idx] = 0.0;
                io->phase_molar_density_mol_per_m3[idx] = 0.0;
                io->phase_compressibility[idx] = 0.0;
                io->phase_fugacity_coeff_available[idx] = 0;
                io->phase_enthalpy_departure_j_per_mol[idx] = dmtherm_optional_double_t{};
                io->phase_entropy_departure_j_per_mol_k[idx] = dmtherm_optional_double_t{};
                io->phase_gibbs_departure_j_per_mol[idx] = dmtherm_optional_double_t{};
                io->phase_enthalpy_j_per_mol[idx] = dmtherm_optional_double_t{};
                io->phase_entropy_j_per_mol_k[idx] = dmtherm_optional_double_t{};
                io->phase_gibbs_j_per_mol[idx] = dmtherm_optional_double_t{};

                if (z_len != 0) {
                    std::fill(io->phase_composition + idx * z_len, io->phase_composition + (idx + 1) * z_len, 0.0);
                    std::fill(io->phase_fugacity_coeff + idx * z_len, io->phase_fugacity_coeff + (idx + 1) * z_len, 0.0);
                }
            }

            std::string msg;
            try {
                const auto fr = h->sys.flashPT(T_in, P_in, zv, fcfg);

                io->status[i] = DMTHERM_STATUS_OK;
                io->converged[i] = fr.converged ? 1 : 0;
                io->iterations[i] = fr.iterations;
                io->residual[i] = fr.residual;
                io->temperature_k[i] = fr.temperature;
                io->pressure_pa_abs[i] = fr.pressure;
                io->vapor_fraction[i] = fr.vapor_fraction;
                io->num_phases[i] = fr.num_phases;
                io->mixture_enthalpy_j_per_mol[i] = opt(fr.enthalpy);
                io->mixture_entropy_j_per_mol_k[i] = opt(fr.entropy);
                io->mixture_gibbs_j_per_mol[i] = opt(fr.gibbs_energy);

                const size_t nph = std::min(fr.phases.size(), stride);
                for (size_t j = 0; j < nph; ++j) {
                    const auto& ph = fr.phases[j];
                    const size_t idx = i * stride + j;

                    io->phase[idx] = toPhase(ph.type);
                    io->phase_fraction[idx] = ph.fraction;
                    io->phase_molar_density_mol_per_m3[idx] = ph.density;
                    io->phase_compressibility[idx] = ph.compressibility;

                    // x (always expected to be size num_components).
                    if (z_len != 0) {
                        const size_t copy_n = std::min(z_len, ph.x.size());
                        if (copy_n) {
                            std::memcpy(io->phase_composition + idx * z_len, ph.x.data(), copy_n * sizeof(double));
                        }
                        if (copy_n < z_len) {
                            std::fill(io->phase_composition + idx * z_len + copy_n, io->phase_composition + (idx + 1) * z_len, 0.0);
                        }
                    }

                    // phi (may be unavailable depending on property method).
                    if (ph.phi.size() == z_len && z_len != 0) {
                        std::memcpy(io->phase_fugacity_coeff + idx * z_len, ph.phi.data(), z_len * sizeof(double));
                        io->phase_fugacity_coeff_available[idx] = 1;
                    } else {
                        io->phase_fugacity_coeff_available[idx] = 0;
                        if (z_len != 0) {
                            std::fill(io->phase_fugacity_coeff + idx * z_len, io->phase_fugacity_coeff + (idx + 1) * z_len, 0.0);
                        }
                    }

                    io->phase_enthalpy_departure_j_per_mol[idx] = opt(ph.H_res);
                    io->phase_entropy_departure_j_per_mol_k[idx] = opt(ph.S_res);
                    io->phase_gibbs_departure_j_per_mol[idx] = opt(ph.G_res);
                    io->phase_enthalpy_j_per_mol[idx] = opt(ph.H);
                    io->phase_entropy_j_per_mol_k[idx] = opt(ph.S);
                    io->phase_gibbs_j_per_mol[idx] = opt(ph.G);
                }

                msg = fr.message;
            } catch (const std::bad_alloc&) {
                any_err = true;
                io->status[i] = DMTHERM_STATUS_OUT_OF_MEMORY;
                last_err = "out of memory";
                msg = last_err;
            } catch (const std::exception& e) {
                any_err = true;
                io->status[i] = DMTHERM_STATUS_RUNTIME_ERROR;
                last_err = e.what();
                msg = last_err;
            } catch (...) {
                any_err = true;
                io->status[i] = DMTHERM_STATUS_INTERNAL_ERROR;
                last_err = "unknown error";
                msg = last_err;
            }

            if (include_messages) {
                if (msg_cursor >= io->messages_capacity) {
                    io->messages_truncated = 1;
                    const size_t off = (io->messages_capacity != 0) ? (io->messages_capacity - 1) : 0;
                    io->message_offsets[i] = off;
                    if (io->messages_capacity != 0) {
                        io->messages[off] = '\0';
                    }
                } else {
                    const size_t off = msg_cursor;
                    io->message_offsets[i] = off;
                    const size_t remaining = io->messages_capacity - off;
                    const size_t copy_n = (msg.size() < (remaining - 1)) ? msg.size() : (remaining - 1);
                    if (copy_n != 0) {
                        std::memcpy(io->messages + off, msg.data(), copy_n);
                    }
                    io->messages[off + copy_n] = '\0';
                    msg_cursor += (copy_n + 1);
                    if (copy_n < msg.size()) io->messages_truncated = 1;
                }
            }
        }

        if (include_messages) {
            io->messages_len = (msg_cursor <= io->messages_capacity) ? msg_cursor : io->messages_capacity;
        } else {
            io->messages_len = 0;
        }

        if (any_err) {
            setLastError(h, last_err);
        } else {
            setLastError(h, "");
        }
        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_flash_pt_batch_into: out of memory");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_flash_pt_batch_into: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

dmtherm_status_t dmtherm_system_flash_tv(
    dmtherm_system_t* sys,
    double T,
    double V,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    dmtherm_flash_pt_result_t* out)
{
    if (!sys || !z || !out) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->struct_size < sizeof(dmtherm_flash_pt_result_t)) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->phases || out->message || out->method_used) return DMTHERM_STATUS_INVALID_ARGUMENT;

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (z_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        Config::TVFlashConfig fcfg = Config::TVFlashConfig{};
        if (cfg && cfg->struct_size >= sizeof(dmtherm_flash_config_t)) {
            if (cfg->max_phases > 0) fcfg.max_phases = cfg->max_phases;
            fcfg.compute_properties = (cfg->compute_properties != 0);
        }

        std::vector<double> zv(z, z + z_len);
        const auto fr = h->sys.flashTV(T, V, zv, fcfg);
        const auto st = fillFlashResult(h, fr, out);
        if (st != DMTHERM_STATUS_OK) return st;

        setLastError(h, "");
        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_flash_tv: out of memory");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_flash_tv: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

dmtherm_status_t dmtherm_system_flash_ph(
    dmtherm_system_t* sys,
    double P,
    double H,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    dmtherm_flash_pt_result_t* out)
{
    if (!sys || !z || !out) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->struct_size < sizeof(dmtherm_flash_pt_result_t)) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->phases || out->message || out->method_used) return DMTHERM_STATUS_INVALID_ARGUMENT;

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (z_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        Config::PHFlashConfig fcfg = Config::PHFlashConfig{};
        if (cfg && cfg->struct_size >= sizeof(dmtherm_flash_config_t)) {
            if (cfg->max_phases > 0) fcfg.max_phases = cfg->max_phases;
            fcfg.compute_properties = (cfg->compute_properties != 0);

            const auto* hs_cfg = static_cast<const dmtherm_flash_hs_config_t*>(cfg->reserved_ptr[0]);
            if (hs_cfg && hs_cfg->struct_size >= sizeof(dmtherm_flash_hs_config_t)) {
                if (hs_cfg->max_iterations > 0) fcfg.max_iterations = hs_cfg->max_iterations;
                if (hs_cfg->tolerance > 0.0) fcfg.enthalpy_tolerance = hs_cfg->tolerance;
                if (hs_cfg->temp_bracket_low > 0.0) fcfg.temp_bracket_low = hs_cfg->temp_bracket_low;
                if (hs_cfg->temp_bracket_high > 0.0) fcfg.temp_bracket_high = hs_cfg->temp_bracket_high;
                fcfg.use_newton = (hs_cfg->use_newton != 0);
            }
        }

        std::vector<double> zv(z, z + z_len);
        const auto fr = h->sys.flashPH(P, H, zv, fcfg);
        const auto st = fillFlashResult(h, fr, out);
        if (st != DMTHERM_STATUS_OK) return st;

        setLastError(h, "");
        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_flash_ph: out of memory");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_flash_ph: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

dmtherm_status_t dmtherm_system_flash_ps(
    dmtherm_system_t* sys,
    double P,
    double S,
    const double* z,
    size_t z_len,
    const dmtherm_flash_config_t* cfg,
    dmtherm_flash_pt_result_t* out)
{
    if (!sys || !z || !out) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->struct_size < sizeof(dmtherm_flash_pt_result_t)) return DMTHERM_STATUS_INVALID_ARGUMENT;
    if (out->phases || out->message || out->method_used) return DMTHERM_STATUS_INVALID_ARGUMENT;

    auto* h = reinterpret_cast<SystemHandle*>(sys);

    try {
        const int nc = h->sys.mixture().numComponents();
        if (z_len != static_cast<size_t>(nc)) return DMTHERM_STATUS_INVALID_ARGUMENT;

        Config::PSFlashConfig fcfg = Config::PSFlashConfig{};
        if (cfg && cfg->struct_size >= sizeof(dmtherm_flash_config_t)) {
            if (cfg->max_phases > 0) fcfg.max_phases = cfg->max_phases;
            fcfg.compute_properties = (cfg->compute_properties != 0);

            const auto* hs_cfg = static_cast<const dmtherm_flash_hs_config_t*>(cfg->reserved_ptr[0]);
            if (hs_cfg && hs_cfg->struct_size >= sizeof(dmtherm_flash_hs_config_t)) {
                if (hs_cfg->max_iterations > 0) fcfg.max_iterations = hs_cfg->max_iterations;
                if (hs_cfg->tolerance > 0.0) fcfg.entropy_tolerance = hs_cfg->tolerance;
                if (hs_cfg->temp_bracket_low > 0.0) fcfg.temp_bracket_low = hs_cfg->temp_bracket_low;
                if (hs_cfg->temp_bracket_high > 0.0) fcfg.temp_bracket_high = hs_cfg->temp_bracket_high;
                fcfg.use_newton = (hs_cfg->use_newton != 0);
            }
        }

        std::vector<double> zv(z, z + z_len);
        const auto fr = h->sys.flashPS(P, S, zv, fcfg);
        const auto st = fillFlashResult(h, fr, out);
        if (st != DMTHERM_STATUS_OK) return st;

        setLastError(h, "");
        return DMTHERM_STATUS_OK;
    } catch (const std::bad_alloc&) {
        setLastError(h, "dmtherm_system_flash_ps: out of memory");
        return DMTHERM_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& e) {
        setLastError(h, e.what());
        return DMTHERM_STATUS_RUNTIME_ERROR;
    } catch (...) {
        setLastError(h, "dmtherm_system_flash_ps: unknown error");
        return DMTHERM_STATUS_INTERNAL_ERROR;
    }
}

} // extern "C"
