/**
 * @file thermo_system.h
 * @brief High-level system object: databanks + mixture + EOS + solvers
 */

#ifndef THERMO_ENGINE_THERMO_SYSTEM_H
#define THERMO_ENGINE_THERMO_SYSTEM_H

#include "thermo/data/databanks.h"
#include "thermo/data/diagnostics.h"
#include "thermo/eos.h"
#include "thermo/factory/mixture_factory.h"
#include "thermo/properties/thermo_properties.h"
#include "thermo/config/error_handling.h"

#include "thermo/equilibrium/iflash.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/equilibrium/tp_departures.h"
#include "thermo/equilibrium/critical_estimation.h"
#include "thermo/reactions/reaction_system.h"
#include "thermo/reactions/reactive_flash.h"
#include "property_engine.h"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace DMThermo {
namespace Engine {

/**
 * @brief Build-time configuration for ThermoSystem.
 *
 * This object exists to make it explicit what comes from CSV databanks vs user inputs:
 * - Pure-component constants/correlations come from databanks during mixture build.
 * - EOS implementations operate only on the built Core::Mixture (no CSV access).
 */
struct ThermoSystemConfig {
    std::filesystem::path databanks_dir;
    std::string eos_type;
    std::vector<std::string> components; // identifiers: CHEMID/CAS/name
    Factory::MixtureBuildOptions mixture_options = [] {
        Factory::MixtureBuildOptions o;
        // ThermoSystem defaults to a strict, "fully-featured" databank build:
        // PH/PS flashes and thermochemistry paths require ICP_* (ideal-gas Cp) for all components.
        o.require_ideal_gas_cp = true;
        return o;
    }();
    std::string flash_method = "Auto";
    // Property routing policy:
    // - EOS_ONLY: always EOS-based properties
    // - DIPPR_ONLY: liquid properties from DIPPR correlations (pure or ideal-mixed)
    // - HYBRID_VAPOR_EOS_LIQUID_DIPPR: vapor from EOS, liquid from DIPPR (pure or ideal-mixed)
    // - AUTO: attempt DIPPR liquid routing when classified as liquid; otherwise EOS
    PropertyMethod property_method = PropertyMethod::EOS_ONLY;
    DipprLiquidRequirements dippr_liquid_requirements{};

    // Transport (DIPPR liquid) policy: make mixing/corrections explicit (no hidden defaults).
    LiquidThermalConductivityMixing liquid_k_mixing = LiquidThermalConductivityMixing::Dippr9I;
    LiquidViscosityMixing liquid_mu_mixing = LiquidViscosityMixing::Log;
    bool enable_liquid_k_high_pressure_correction_9g1 = true;

    // Error handling policy:
    // - Strict: throw categorized exceptions for routing/solver failures.
    // - Lenient: return failure objects where possible.
    Config::FailureMode failure_mode = Config::FailureMode::Strict;
};

/**
 * @brief TP phase evaluation result for a selected density root.
 *
 * Provides a process-simulator style snapshot at (T,P,x) for a chosen root:
 * - rho, Z
 * - ln(phi_i)
 * - departure functions h/s/g (relative to ideal gas at the same T,P,x)
 * - total H/S/G (relative to the library reference state conventions)
 */
struct TPPhaseEvaluation {
    double T = 0.0;                         // K
    double P = 0.0;                         // Pa
    std::vector<double> x;                  // composition used for evaluation

    PhaseType phase = PhaseType::Unknown;   // phase hint/classification for this root
    double rho = 0.0;                       // mol/m^3
    double Z = 0.0;                         // [-]

    std::vector<double> ln_phi;             // [-]

    double h_dep = 0.0;                     // J/mol
    double s_dep = 0.0;                     // J/(mol*K)
    double g_dep = 0.0;                     // J/mol

    double H = 0.0;                         // J/mol
    double S = 0.0;                         // J/(mol*K)
    double G = 0.0;                         // J/mol

    double Cp = 0.0;                        // J/(mol*K)
    double Cv = 0.0;                        // J/(mol*K)

    bool success = false;
    std::string message;
};

/**
 * @brief A self-contained, databank-backed thermo system.
 *
 * Owns:
 * - Loaded databanks (DIPPR/PC-SAFT/BIPs)
 * - Built Core::Mixture
 * - Selected EOS instance
 * - Default stability analyzer and flash solver for the EOS
 */
class ThermoSystem {
public:
    static ThermoSystem build(const ThermoSystemConfig& cfg, Data::Diagnostics* diag = nullptr);

    const Data::Databanks& databanks() const { return databanks_; }
    const Core::Mixture& mixture() const { return mixture_; }
    EOSPtr eos() const { return eos_; }
    Equilibrium::Stability::StabilityAnalyzerPtr stability() const { return stability_; }
    Equilibrium::Flash::FlashSolverPtr flash() const { return flash_; }

    Properties::ThermoProperties propertiesTV(double T, double rho, const std::vector<double>& x) const;

    /**
     * @brief EOS-based properties at (T,P) for a given phase hint.
     *
     * Finds a density root using the stability analyzer, then computes TV properties at that density.
     * If no root is found, falls back to the ideal-gas density as a best-effort (success=false).
     */
    Properties::ThermoProperties propertiesTP(
        double T,
        double P,
        const std::vector<double>& x,
        PhaseType phase_hint = PhaseType::Unknown
    ) const;

    Properties::ThermoProperties propertiesTP(
        double T,
        double P,
        const std::vector<double>& x,
        PhaseType phase_hint,
        Data::Diagnostics* diag
    ) const;

    Properties::ThermoProperties propertiesTPWithMethod(
        double T,
        double P,
        const std::vector<double>& x,
        PhaseType phase_hint,
        PropertyMethod method,
        const DipprLiquidRequirements& dippr_req = DipprLiquidRequirements{},
        Data::Diagnostics* diag = nullptr
    ) const;

    /**
     * @brief Evaluate a single TP phase root (vapor/liquid/stable) with departures + H/S/G.
     *
     * This is EOS-based (requires ln(phi) + EOS derivatives for departures).
     */
    TPPhaseEvaluation evaluateTPPhase(
        double T,
        double P,
        const std::vector<double>& x,
        Equilibrium::Departures::RootSelection selection = Equilibrium::Departures::RootSelection::Stable,
        const Config::DensityRootConfig& config = Config::DensityRootConfig{}
    ) const;

    /**
     * @brief TP departure properties per density root (Z, ln(phi), h/s/g departures).
     *
     * This is a convenience wrapper around `Equilibrium::Departures::departureRootsTP`
     * using the system's default stability analyzer.
     */
    Equilibrium::Departures::TPDepartureResult departureRootsTP(
        double T,
        double P,
        const std::vector<double>& x,
        const Config::DensityRootConfig& config = Config::DensityRootConfig{}
    ) const;

    /**
     * @brief Estimate pure-component critical point for component index.
     *
     * Uses a pure-component limit composition (one-hot in system component space)
     * and an EOS-based critical solve.
     */
    Equilibrium::Critical::EstimateResult estimateCritical(
        int component_index = 0,
        const Equilibrium::Critical::EstimateInputs& inputs = Equilibrium::Critical::EstimateInputs{}
    ) const;

    Equilibrium::Flash::FlashResult flashPT(
        double T,
        double P,
        const std::vector<double>& z,
        const Config::FlashConfig& cfg = Config::FlashConfig::defaults()
    ) const;

    Equilibrium::Flash::FlashResult flashTV(
        double T,
        double V,
        const std::vector<double>& z,
        const Config::TVFlashConfig& cfg = Config::TVFlashConfig{}
    ) const;

    Equilibrium::Flash::FlashResult flashPH(
        double P,
        double H,
        const std::vector<double>& z,
        const Config::PHFlashConfig& cfg = Config::PHFlashConfig{}
    ) const;

    Equilibrium::Flash::FlashResult flashPS(
        double P,
        double S,
        const std::vector<double>& z,
        const Config::PSFlashConfig& cfg = Config::PSFlashConfig{}
    ) const;

    // -------------------------------------------------------------------------
    // Reactive flash (reaction + phase equilibrium)
    // -------------------------------------------------------------------------

    Reactions::ReactiveFlashResult reactiveFlashTP(
        double T,
        double P,
        const std::vector<double>& n0,
        const Reactions::ReactionSystem& system,
        const Reactions::ReactiveFlashConfig& cfg = Reactions::ReactiveFlashConfig{},
        Data::Diagnostics* diag = nullptr
    ) const;

    Reactions::ReactiveFlashResult reactiveFlashTV(
        double T,
        double V,
        const std::vector<double>& n0,
        const Reactions::ReactionSystem& system,
        const Reactions::ReactiveFlashConfig& cfg = Reactions::ReactiveFlashConfig{},
        Data::Diagnostics* diag = nullptr
    ) const;

    Reactions::ReactiveFlashResult reactiveFlashPH(
        double P,
        double H,
        const std::vector<double>& n0,
        const Reactions::ReactionSystem& system,
        const Reactions::ReactivePHFlashConfig& cfg = Reactions::ReactivePHFlashConfig{},
        Data::Diagnostics* diag = nullptr
    ) const;

    Reactions::ReactiveFlashResult reactiveFlashPS(
        double P,
        double S,
        const std::vector<double>& n0,
        const Reactions::ReactionSystem& system,
        const Reactions::ReactivePSFlashConfig& cfg = Reactions::ReactivePSFlashConfig{},
        Data::Diagnostics* diag = nullptr
    ) const;

private:
    Data::Databanks databanks_{};
    Core::Mixture mixture_{std::vector<Core::Component>{}, Core::BinaryParameters::zeros(0)};
    EOSPtr eos_{};
    Equilibrium::Stability::StabilityAnalyzerPtr stability_{};
    Equilibrium::Flash::FlashSolverPtr flash_{};
    PropertyMethod property_method_ = PropertyMethod::EOS_ONLY;
    DipprLiquidRequirements dippr_liquid_requirements_{};
    LiquidThermalConductivityMixing liquid_k_mixing_ = LiquidThermalConductivityMixing::Dippr9I;
    LiquidViscosityMixing liquid_mu_mixing_ = LiquidViscosityMixing::Log;
    bool enable_liquid_k_high_pressure_correction_9g1_ = true;
    Config::FailureMode failure_mode_ = Config::FailureMode::Strict;
};

} // namespace Engine
} // namespace DMThermo

#endif // THERMO_ENGINE_THERMO_SYSTEM_H
