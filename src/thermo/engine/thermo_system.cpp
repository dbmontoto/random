/**
 * @file thermo_system.cpp
 * @brief Implementation of ThermoSystem
 */

#include "thermo/engine/thermo_system.h"

#include "thermo/engine/property_engine.h"
#include "thermo/factory/eos_factory.h"
#include "thermo/factory/flash_factory.h"
#include "thermo/factory/stability_factory.h"
#include "thermo/core/constants.h"
#include "thermo/core/units.h"
#include "thermo/data/component_resolver.h"
#include "thermo/data/dippr_adapter.h"
#include "thermo/equilibrium/thermo_contributions.h"
#include "thermo/reactions/reactive_flash.h"
#include "thermo/dippr/dippr.h"

#include <string>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace DMThermo {
namespace Engine {

namespace {

const char* methodToString(PropertyMethod m) {
    switch (m) {
        case PropertyMethod::AUTO: return "AUTO";
        case PropertyMethod::EOS_ONLY: return "EOS_ONLY";
        case PropertyMethod::DIPPR_ONLY: return "DIPPR_ONLY";
        case PropertyMethod::HYBRID_VAPOR_EOS_LIQUID_DIPPR: return "HYBRID_VAPOR_EOS_LIQUID_DIPPR";
        default: return "UNKNOWN";
    }
}

const char* dipprLiquidKind(int nc) {
    return (nc == 1) ? "pure" : "ideal-mix";
}

void appendDiagnostic(
    Data::Diagnostics* diag,
    Data::DiagnosticLevel level,
    Config::ErrorCategory category,
    const std::string& code,
    const std::string& message)
{
    if (!diag) return;
    switch (level) {
        case Data::DiagnosticLevel::Info:
            diag->info(message, category, code);
            break;
        case Data::DiagnosticLevel::Warning:
            diag->warn(message, category, code);
            break;
        case Data::DiagnosticLevel::Error:
            diag->error(message, category, code);
            break;
    }
}

Properties::ThermoProperties tpFailureResult(
    double T,
    double P,
    const EOSPtr& eos,
    const std::string& message)
{
    Properties::ThermoProperties out;
    out.T = T;
    out.P = P;
    out.model = eos ? eos->name() : std::string{};
    out.success = false;
    out.message = message;
    return out;
}

std::optional<double> estimateVcHintFromCriticals(const Core::Component& component) {
    const double Tc = component.Tc();
    const double Pc = component.Pc();
    const double Zc = component.Zc();
    if (!(std::isfinite(Tc) && Tc > 0.0 && std::isfinite(Pc) && Pc > 0.0 && std::isfinite(Zc) && Zc > 0.0)) {
        return std::nullopt;
    }
    const Core::Units::PropertyVariable critical_volume(
        (Zc * Constants::GAS_CONSTANT * Tc) / Pc,
        Core::Units::Unit::M3_PER_MOL);
    const double Vc_hint = critical_volume.as(Core::Units::Unit::M3_PER_KMOL);
    if (!(std::isfinite(Vc_hint) && Vc_hint > 0.0)) {
        return std::nullopt;
    }
    return Vc_hint;
}

void populateFlashPropertiesTP(
    const ThermoSystem& sys,
    double T,
    double P,
    const Config::FlashConfig& cfg,
    Equilibrium::Flash::FlashResult& fr)
{
    if (!cfg.compute_properties) return;
    if (!fr.converged || fr.phases.empty()) return;

    double H_total = 0.0;
    double S_total = 0.0;
    double G_total = 0.0;
    bool all_ok = true;

    for (auto& ph : fr.phases) {
        Equilibrium::Departures::RootSelection sel = Equilibrium::Departures::RootSelection::Stable;
        if (ph.type == PhaseType::Vapor) {
            sel = Equilibrium::Departures::RootSelection::Vapor;
        } else if (ph.type == PhaseType::Liquid || ph.type == PhaseType::Liquid1 || ph.type == PhaseType::Liquid2) {
            sel = Equilibrium::Departures::RootSelection::Liquid;
        }

        std::vector<double> x_phase;
        try {
            x_phase = Core::Mixture::normalizeComposition(ph.x);
        } catch (...) {
            all_ok = false;
            continue;
        }

        const auto ev = sys.evaluateTPPhase(T, P, x_phase, sel);
        if (!ev.success) {
            all_ok = false;
            continue;
        }

        ph.type = (ph.type == PhaseType::Unknown) ? ev.phase : ph.type;
        ph.density = ev.rho;
        ph.compressibility = ev.Z;
        ph.x = x_phase;

        ph.phi.assign(ev.ln_phi.size(), 1.0);
        for (size_t i = 0; i < ev.ln_phi.size(); ++i) {
            const double lnphi = ev.ln_phi[i];
            const double clamped = std::clamp(lnphi, -700.0, 700.0);
            ph.phi[i] = std::exp(clamped);
        }

        ph.H_res = ev.h_dep;
        ph.S_res = ev.s_dep;
        ph.G_res = ev.g_dep;
        ph.H = ev.H;
        ph.S = ev.S;
        ph.G = ev.G;

        H_total += ph.fraction * ev.H;
        S_total += ph.fraction * ev.S;
        G_total += ph.fraction * ev.G;
    }

    if (all_ok) {
        fr.enthalpy = H_total;
        fr.entropy = S_total;
        fr.gibbs_energy = G_total;
    }
}

} // namespace

ThermoSystem ThermoSystem::build(const ThermoSystemConfig& cfg, Data::Diagnostics* diag) {
    if (cfg.databanks_dir.empty()) {
        throw std::invalid_argument("ThermoSystem::build: databanks_dir is empty");
    }
    if (cfg.components.empty()) {
        throw std::invalid_argument("ThermoSystem::build: components list is empty");
    }
    if (cfg.eos_type.empty()) {
        throw std::invalid_argument("ThermoSystem::build: eos_type is empty");
    }

    ThermoSystem sys;
    const std::size_t diag_start = diag ? diag->items().size() : 0;
    const bool ok = sys.databanks_.loadAllFromDirectory(cfg.databanks_dir, diag);
    bool new_errors = false;
    if (diag) {
        const auto& items = diag->items();
        for (std::size_t i = diag_start; i < items.size(); ++i) {
            if (items[i].level == Data::DiagnosticLevel::Error) {
                new_errors = true;
                break;
            }
        }
    }
    if (!ok || new_errors) {
        appendDiagnostic(
            diag,
            Data::DiagnosticLevel::Error,
            Config::ErrorCategory::MissingData,
            "DATABANK_LOAD_FAILED",
            "ThermoSystem::build: failed to load databanks from " + cfg.databanks_dir.string());
        throw std::runtime_error("ThermoSystem::build: failed to load databanks from " + cfg.databanks_dir.string());
    }

    sys.mixture_ = Factory::buildMixtureFromDatabanks(sys.databanks_, cfg.components, cfg.mixture_options, diag);
    sys.eos_ = Factory::EOSFactory::create(cfg.eos_type, sys.mixture_);

    sys.stability_ = Factory::StabilityFactory::create(sys.eos_);
    sys.flash_ = Factory::FlashFactory::create(cfg.flash_method, sys.eos_, sys.stability_);
    sys.property_method_ = cfg.property_method;
    sys.dippr_liquid_requirements_ = cfg.dippr_liquid_requirements;
    sys.liquid_k_mixing_ = cfg.liquid_k_mixing;
    sys.liquid_mu_mixing_ = cfg.liquid_mu_mixing;
    sys.enable_liquid_k_high_pressure_correction_9g1_ = cfg.enable_liquid_k_high_pressure_correction_9g1;
    sys.failure_mode_ = cfg.failure_mode;
    return sys;
}

Properties::ThermoProperties ThermoSystem::propertiesTV(double T, double rho, const std::vector<double>& x) const {
    if (!eos_) {
        throw std::runtime_error("ThermoSystem::propertiesTV: EOS is null");
    }
    Engine::PropertyEngine engine(eos_, mixture_);
    engine.setLiquidThermalConductivityMixing(liquid_k_mixing_);
    engine.setLiquidViscosityMixing(liquid_mu_mixing_);
    engine.setLiquidKHighPressureCorrection9G1Enabled(enable_liquid_k_high_pressure_correction_9g1_);
    return engine.calculateTV(T, rho, x);
}

Properties::ThermoProperties ThermoSystem::propertiesTP(
    double T,
    double P,
    const std::vector<double>& x,
    PhaseType phase_hint) const
{
    return propertiesTPWithMethod(T, P, x, phase_hint, property_method_, dippr_liquid_requirements_, nullptr);
}

Properties::ThermoProperties ThermoSystem::propertiesTP(
    double T,
    double P,
    const std::vector<double>& x,
    PhaseType phase_hint,
    Data::Diagnostics* diag) const
{
    return propertiesTPWithMethod(T, P, x, phase_hint, property_method_, dippr_liquid_requirements_, diag);
}

Properties::ThermoProperties ThermoSystem::propertiesTPWithMethod(
    double T,
    double P,
    const std::vector<double>& x,
    PhaseType phase_hint,
    PropertyMethod method,
    const DipprLiquidRequirements& dippr_req,
    Data::Diagnostics* diag) const
{
    auto failOrThrow = [&](Config::ErrorCategory category, const std::string& code, const std::string& message)
        -> Properties::ThermoProperties
    {
        appendDiagnostic(diag, Data::DiagnosticLevel::Error, category, code, message);
        if (failure_mode_ == Config::FailureMode::Strict) {
            throw std::runtime_error(message);
        }
        return tpFailureResult(T, P, eos_, message);
    };

    if (!eos_ || !stability_) {
        return failOrThrow(
            Config::ErrorCategory::NotInitialized,
            "EOS_OR_STABILITY_NOT_INITIALIZED",
            "ThermoSystem::propertiesTP: EOS/stability not initialized");
    }

    auto isLiquidHint = [](PhaseType ph) {
        return ph == PhaseType::Liquid || ph == PhaseType::Liquid1 || ph == PhaseType::Liquid2;
    };

    auto resolveAutoHint = [&](PhaseType hint) -> PhaseType {
        if (hint != PhaseType::Unknown) return hint;
        const auto roots = stability_->findDensityRoots(T, P, x);
        if (roots.stable_root_index.has_value()) {
            const int k = roots.stable_root_index.value();
            if (k >= 0 && k < static_cast<int>(roots.phase_types.size())) {
                return roots.phase_types[k];
            }
            if (k >= 0 && k < static_cast<int>(roots.densities.size())) {
                return stability_->classifyPhase(T, roots.densities[k], x);
            }
        }
        return PhaseType::Unknown;
    };

    auto dipprLiquidProps = [&]() -> Properties::ThermoProperties {
        Engine::PropertyEngine engine(eos_, mixture_);
        engine.setLiquidThermalConductivityMixing(liquid_k_mixing_);
        engine.setLiquidViscosityMixing(liquid_mu_mixing_);
        engine.setLiquidKHighPressureCorrection9G1Enabled(enable_liquid_k_high_pressure_correction_9g1_);

        if (mixture_.numComponents() == 1) {
            const auto& c = mixture_.component(0);
            std::string id = !c.cas().empty() ? c.cas() : c.name();
            if (id.empty()) {
                throw std::runtime_error("ThermoSystem::propertiesTP: cannot resolve identifier for pure component (missing CAS/name)");
            }
            return engine.calculatePureLiquidDIPPR(databanks_, id, T, P, dippr_req, diag);
        }

        return engine.calculateLiquidDIPPRIdealMix(databanks_, T, P, x, dippr_req, diag);
    };

    const int nc = mixture_.numComponents();
    const std::string hint_s = phaseTypeToString(phase_hint);

    // DIPPR-only: liquid-only properties from DIPPR (pure or ideal-mix).
    if (method == PropertyMethod::DIPPR_ONLY) {
        const PhaseType inferred = resolveAutoHint(phase_hint);
        appendDiagnostic(
            diag,
            Data::DiagnosticLevel::Info,
            Config::ErrorCategory::None,
            "ROUTING_DECISION",
            std::string("Routing: method=") + methodToString(method) +
            ", hint=" + hint_s +
            ", inferred=" + phaseTypeToString(inferred) +
            " -> DIPPR(" + dipprLiquidKind(nc) + ")"
        );
        if (!isLiquidHint(inferred)) {
            return failOrThrow(
                Config::ErrorCategory::InvalidInput,
                "ROUTING_REJECTED_NON_LIQUID",
                std::string("ThermoSystem::propertiesTP: DIPPR_ONLY requires a liquid phase hint (or a TP state classified as liquid); hint=") +
                    hint_s + ", inferred=" + phaseTypeToString(inferred));
        }

        auto props = dipprLiquidProps();
        if (!props.success) {
            return failOrThrow(
                Config::ErrorCategory::MissingData,
                "DIPPR_ROUTING_FAILED",
                std::string("ThermoSystem::propertiesTP: DIPPR_ONLY liquid routing failed: ") + props.message);
        }
        return props;
    }

    // HYBRID: vapor EOS, liquid DIPPR (pure or ideal-mix).
    if (method == PropertyMethod::HYBRID_VAPOR_EOS_LIQUID_DIPPR) {
        const PhaseType inferred = resolveAutoHint(phase_hint);
        if (isLiquidHint(inferred)) {
            appendDiagnostic(
                diag,
                Data::DiagnosticLevel::Info,
                Config::ErrorCategory::None,
                "ROUTING_DECISION",
                std::string("Routing: method=") + methodToString(method) +
                ", hint=" + hint_s +
                ", inferred=" + phaseTypeToString(inferred) +
                " -> DIPPR(" + dipprLiquidKind(nc) + ")"
            );
            auto props = dipprLiquidProps();
            if (!props.success) {
                return failOrThrow(
                    Config::ErrorCategory::MissingData,
                    "DIPPR_ROUTING_FAILED",
                    std::string("ThermoSystem::propertiesTP: HYBRID liquid DIPPR failed: ") + props.message);
            }
            return props;
        }
        appendDiagnostic(
            diag,
            Data::DiagnosticLevel::Info,
            Config::ErrorCategory::None,
            "ROUTING_DECISION",
            std::string("Routing: method=") + methodToString(method) +
            ", hint=" + hint_s +
            ", inferred=" + phaseTypeToString(inferred) +
            " -> EOS"
        );
        // Otherwise fall through to EOS (vapor or unknown).
    }

    // AUTO: attempt DIPPR liquid routing if the state is classified as liquid; otherwise EOS.
    if (method == PropertyMethod::AUTO) {
        const PhaseType inferred = resolveAutoHint(phase_hint);
        if (isLiquidHint(inferred)) {
            auto props = dipprLiquidProps();
            if (props.success) {
                appendDiagnostic(
                    diag,
                    Data::DiagnosticLevel::Info,
                    Config::ErrorCategory::None,
                    "ROUTING_DECISION",
                    std::string("Routing: method=") + methodToString(method) +
                    ", hint=" + hint_s +
                    ", inferred=" + phaseTypeToString(inferred) +
                    " -> DIPPR(" + dipprLiquidKind(nc) + ")"
                );
                return props;
            }
            appendDiagnostic(
                diag,
                Data::DiagnosticLevel::Warning,
                Config::ErrorCategory::MissingData,
                "ROUTING_FALLBACK_TO_EOS",
                std::string("Routing fallback: method=") + methodToString(method) +
                ", hint=" + hint_s +
                ", inferred=" + phaseTypeToString(inferred) +
                " DIPPR(" + dipprLiquidKind(nc) + ") failed: " + props.message +
                " -> EOS"
            );
        }
        appendDiagnostic(
            diag,
            Data::DiagnosticLevel::Info,
            Config::ErrorCategory::None,
            "ROUTING_DECISION",
            std::string("Routing: method=") + methodToString(method) +
            ", hint=" + hint_s +
            ", inferred=" + phaseTypeToString(inferred) +
            " -> EOS"
        );
        // Fall through to EOS.
    }

    if (method == PropertyMethod::EOS_ONLY) {
        appendDiagnostic(
            diag,
            Data::DiagnosticLevel::Info,
            Config::ErrorCategory::None,
            "ROUTING_DECISION",
            std::string("Routing: method=") + methodToString(method) + ", hint=" + hint_s + " -> EOS");
    }

    double rho = 0.0;
    try {
        if (phase_hint == PhaseType::Vapor) {
            rho = stability_->findVaporDensity(T, P, x);
        } else if (phase_hint == PhaseType::Liquid || phase_hint == PhaseType::Liquid1 || phase_hint == PhaseType::Liquid2) {
            rho = stability_->findLiquidDensity(T, P, x);
        } else {
            const auto roots = stability_->findDensityRoots(T, P, x);
            if (roots.stable_root_index.has_value()) {
                rho = roots.densities[roots.stable_root_index.value()];
            } else if (!roots.densities.empty()) {
                rho = roots.densities.front();
            }
        }
    } catch (const std::invalid_argument& e) {
        return failOrThrow(
            Config::ErrorCategory::InvalidInput,
            "DENSITY_ROOT_INPUT_ERROR",
            std::string("ThermoSystem::propertiesTP: density root evaluation failed due to invalid input: ") + e.what());
    } catch (const std::exception& e) {
        return failOrThrow(
            Config::ErrorCategory::NumericalFailure,
            "DENSITY_ROOT_EVALUATION_FAILED",
            std::string("ThermoSystem::propertiesTP: density root evaluation failed: ") + e.what());
    }

    if (!(std::isfinite(rho) && rho > 0.0)) {
        return failOrThrow(
            Config::ErrorCategory::NumericalFailure,
            "DENSITY_ROOT_NOT_FOUND",
            "ThermoSystem::propertiesTP: no density root found for requested TP state");
    }

    auto props = propertiesTV(T, rho, x);
    if (props.success) {
        // Overwrite reported pressure with the target pressure for clarity.
        props.P = P;
        props.Z = (rho > 0.0) ? (P / (rho * Constants::GAS_CONSTANT * T)) : props.Z;

        // Ensure thermodynamic identities remain consistent with the reported TP state.
        if (props.v > 0.0) {
            props.U = props.H - props.P * props.v;
            props.G = props.H - T * props.S;
            props.A = props.U - T * props.S;
        }
    }
    return props;
}

TPPhaseEvaluation ThermoSystem::evaluateTPPhase(
    double T,
    double P,
    const std::vector<double>& x,
    Equilibrium::Departures::RootSelection selection,
    const Config::DensityRootConfig& config) const
{
    TPPhaseEvaluation out;
    out.T = T;
    out.P = P;
    out.x = x;

    try {
        if (!eos_ || !stability_) {
            throw std::runtime_error("ThermoSystem::evaluateTPPhase: EOS/stability not initialized");
        }
        if (!mixture_.isValidComposition(x)) {
            throw std::invalid_argument("ThermoSystem::evaluateTPPhase: invalid composition vector");
        }

        const auto dep = Equilibrium::Departures::departureRootsTP(stability_, T, P, x, config);
        const auto idx = Equilibrium::Departures::selectRootIndex(dep, selection);
        if (!idx.has_value()) {
            out.success = false;
            out.message = "No valid TP departure roots available for requested state";
            return out;
        }

        const auto& r = dep.roots[static_cast<size_t>(idx.value())];
        out.phase = r.phase;
        out.rho = r.density;
        out.Z = r.Z;
        out.ln_phi = r.ln_phi;
        out.h_dep = r.h_dep;
        out.s_dep = r.s_dep;
        out.g_dep = r.g_dep;

        const auto ig = Equilibrium::Contributions::idealGasProps(mixture_, T, P, x);
        out.H = ig.h_ig + out.h_dep;
        out.S = ig.s_ig + out.s_dep;
        out.G = out.H - T * out.S;

        const auto props = propertiesTV(T, out.rho, x);
        if (props.success) {
            if (props.Cp.has_value()) out.Cp = props.Cp.value();
            if (props.Cv.has_value()) out.Cv = props.Cv.value();
        }

        out.success = true;
        return out;
    } catch (const std::exception& e) {
        out.success = false;
        out.message = e.what();
        return out;
    }
}

Equilibrium::Departures::TPDepartureResult ThermoSystem::departureRootsTP(
    double T,
    double P,
    const std::vector<double>& x,
    const Config::DensityRootConfig& config) const
{
    if (!stability_) {
        throw std::runtime_error("ThermoSystem::departureRootsTP: stability analyzer not initialized");
    }
    return Equilibrium::Departures::departureRootsTP(stability_, T, P, x, config);
}

Equilibrium::Critical::EstimateResult ThermoSystem::estimateCritical(
    int component_index,
    const Equilibrium::Critical::EstimateInputs& inputs) const
{
    if (!eos_) {
        throw std::runtime_error("ThermoSystem::estimateCritical: EOS not initialized");
    }

    const int nc = mixture_.numComponents();
    if (nc <= 0) {
        throw std::runtime_error("ThermoSystem::estimateCritical: mixture has no components");
    }
    if (component_index < 0 || component_index >= nc) {
        throw std::invalid_argument("ThermoSystem::estimateCritical: component index out of range");
    }

    std::vector<double> x(static_cast<size_t>(nc), 0.0);
    x[static_cast<size_t>(component_index)] = 1.0;

    Equilibrium::Critical::EstimateInputs effective_inputs = inputs;
    const auto& component = mixture_.component(component_index);

    if (!(std::isfinite(effective_inputs.Tc_hint) && effective_inputs.Tc_hint > 0.0)) {
        const double Tc = component.Tc();
        if (std::isfinite(Tc) && Tc > 0.0) {
            effective_inputs.Tc_hint = Tc;
        }
    }
    if (!(std::isfinite(effective_inputs.Pc_hint) && effective_inputs.Pc_hint > 0.0)) {
        const double Pc = component.Pc();
        if (std::isfinite(Pc) && Pc > 0.0) {
            effective_inputs.Pc_hint = Pc;
        }
    }
    const double critical_volume_hint = effective_inputs.criticalVolumeHint().as(Core::Units::Unit::M3_PER_KMOL);
    if (!(std::isfinite(critical_volume_hint) && critical_volume_hint > 0.0)) {
        const auto critical_volume_estimate = estimateVcHintFromCriticals(component);
        if (critical_volume_estimate.has_value()) {
            const Core::Units::PropertyVariable critical_volume(
                critical_volume_estimate.value(),
                Core::Units::Unit::M3_PER_KMOL);
            effective_inputs.Vc_hint_m3_per_kmol = critical_volume.as(Core::Units::Unit::M3_PER_KMOL);
        }
    }

    return Equilibrium::Critical::estimateFromEOS(*eos_, x, effective_inputs);
}

Equilibrium::Flash::FlashResult ThermoSystem::flashPT(
    double T,
    double P,
    const std::vector<double>& z,
    const Config::FlashConfig& cfg) const
{
    if (!flash_) {
        throw std::runtime_error("ThermoSystem::flashPT: flash solver not initialized");
    }

    Config::FlashConfig use_cfg = cfg;
    if (use_cfg.use_dippr_psat_k_init) {
        const int nc = mixture_.numComponents();
        if (static_cast<int>(z.size()) != nc) {
            throw std::invalid_argument("ThermoSystem::flashPT: composition size mismatch");
        }

        std::vector<double> K = flash_->estimateKValues(T, P); // fallback (Wilson)
        if (static_cast<int>(K.size()) != nc) {
            K.assign(static_cast<std::size_t>(nc), 1.0);
        }

        const Data::ComponentResolver resolver(databanks_);
        for (int i = 0; i < nc; ++i) {
            const auto& c = mixture_.component(i);
            const std::string id = !c.cas().empty() ? c.cas() : c.name();
            if (id.empty()) continue;

            const auto resolved = resolver.resolve(id, nullptr);
            if (!resolved.dippr) continue;

            auto params = Data::toDipprParameters(*resolved.dippr, nullptr);
            if (!params) continue;

            try {
                DIPPR::Calculator calc(*params);
                const auto& p = calc.getParameters();
                if (!p.has_vapor_pressure) continue;

                const double Psat = calc.vaporPressure(T);
                if (!(std::isfinite(Psat) && Psat > 0.0 && std::isfinite(P) && P > 0.0)) continue;

                const double Ki = Psat / P;
                if (std::isfinite(Ki) && Ki > 0.0) {
                    K[static_cast<std::size_t>(i)] = std::clamp(Ki, use_cfg.min_k_value, use_cfg.max_k_value);
                }
            } catch (...) {
                // Best-effort only; keep fallback K[i].
            }
        }

        use_cfg.k_init = Config::KValueInit::Custom;
        use_cfg.custom_k_values = std::move(K);
    }

    auto populateFlashProperties = [&](Equilibrium::Flash::FlashResult& fr) {
        populateFlashPropertiesTP(*this, T, P, use_cfg, fr);
    };

    Equilibrium::Flash::FlashResult out;
    if (use_cfg.max_phases >= 3) {
        out = flash_->solveVLLE(T, P, z, use_cfg);
    } else {
        out = flash_->solvePT(T, P, z, use_cfg);
    }

    populateFlashProperties(out);
    return out;
}

Equilibrium::Flash::FlashResult ThermoSystem::flashTV(
    double T,
    double V,
    const std::vector<double>& z,
    const Config::TVFlashConfig& cfg) const
{
    if (!flash_) {
        throw std::runtime_error("ThermoSystem::flashTV: flash solver not initialized");
    }

    auto out = flash_->solveTV(T, V, z, cfg);
    populateFlashPropertiesTP(*this, out.temperature > 0.0 ? out.temperature : T, out.pressure, cfg, out);
    return out;
}

Equilibrium::Flash::FlashResult ThermoSystem::flashPH(
    double P,
    double H,
    const std::vector<double>& z,
    const Config::PHFlashConfig& cfg) const
{
    if (!flash_) {
        throw std::runtime_error("ThermoSystem::flashPH: flash solver not initialized");
    }

    auto out = flash_->solvePH(P, H, z, cfg);
    populateFlashPropertiesTP(*this, out.temperature, out.pressure > 0.0 ? out.pressure : P, cfg, out);
    return out;
}

Equilibrium::Flash::FlashResult ThermoSystem::flashPS(
    double P,
    double S,
    const std::vector<double>& z,
    const Config::PSFlashConfig& cfg) const
{
    if (!flash_) {
        throw std::runtime_error("ThermoSystem::flashPS: flash solver not initialized");
    }

    auto out = flash_->solvePS(P, S, z, cfg);
    populateFlashPropertiesTP(*this, out.temperature, out.pressure > 0.0 ? out.pressure : P, cfg, out);
    return out;
}

Reactions::ReactiveFlashResult ThermoSystem::reactiveFlashTP(
    double T,
    double P,
    const std::vector<double>& n0,
    const Reactions::ReactionSystem& system,
    const Reactions::ReactiveFlashConfig& cfg,
    Data::Diagnostics* diag) const
{
    if (!eos_ || !stability_) {
        throw std::runtime_error("ThermoSystem::reactiveFlashTP: EOS/stability not initialized");
    }

    Reactions::ReactiveFlashSolver solver(eos_, stability_);
    return solver.solveTPFromDatabanks(T, P, n0, system, databanks_, mixture_, cfg, diag);
}

Reactions::ReactiveFlashResult ThermoSystem::reactiveFlashTV(
    double T,
    double V,
    const std::vector<double>& n0,
    const Reactions::ReactionSystem& system,
    const Reactions::ReactiveFlashConfig& cfg,
    Data::Diagnostics* diag) const
{
    if (!eos_ || !stability_) {
        throw std::runtime_error("ThermoSystem::reactiveFlashTV: EOS/stability not initialized");
    }

    Reactions::ReactiveFlashSolver solver(eos_, stability_);
    return solver.solveTVFromDatabanks(T, V, n0, system, databanks_, mixture_, cfg, diag);
}

Reactions::ReactiveFlashResult ThermoSystem::reactiveFlashPH(
    double P,
    double H,
    const std::vector<double>& n0,
    const Reactions::ReactionSystem& system,
    const Reactions::ReactivePHFlashConfig& cfg,
    Data::Diagnostics* diag) const
{
    if (!eos_ || !stability_) {
        throw std::runtime_error("ThermoSystem::reactiveFlashPH: EOS/stability not initialized");
    }

    Reactions::ReactiveFlashSolver solver(eos_, stability_);
    return solver.solvePHFromDatabanks(P, H, n0, system, databanks_, mixture_, cfg, diag);
}

Reactions::ReactiveFlashResult ThermoSystem::reactiveFlashPS(
    double P,
    double S,
    const std::vector<double>& n0,
    const Reactions::ReactionSystem& system,
    const Reactions::ReactivePSFlashConfig& cfg,
    Data::Diagnostics* diag) const
{
    if (!eos_ || !stability_) {
        throw std::runtime_error("ThermoSystem::reactiveFlashPS: EOS/stability not initialized");
    }

    Reactions::ReactiveFlashSolver solver(eos_, stability_);
    return solver.solvePSFromDatabanks(P, S, n0, system, databanks_, mixture_, cfg, diag);
}

} // namespace Engine
} // namespace DMThermo
