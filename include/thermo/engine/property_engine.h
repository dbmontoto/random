/**
 * @file property_engine.h
 * @brief High-level thermophysical property routing (EOS/DIPPR/hybrid)
 */

#ifndef THERMO_ENGINE_PROPERTY_ENGINE_H
#define THERMO_ENGINE_PROPERTY_ENGINE_H

#include "thermo/core/mixture.h"
#include "thermo/data/databanks.h"
#include "thermo/data/diagnostics.h"
#include "thermo/eos.h"
#include "thermo/properties/thermo_properties.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace DMThermo {
namespace Engine {

enum class PropertyMethod {
    AUTO,
    EOS_ONLY,
    DIPPR_ONLY,
    HYBRID_VAPOR_EOS_LIQUID_DIPPR
};

enum class LiquidThermalConductivityMixing {
    SimpleLinear, // k = Σ x_i k_i
    Dippr9I       // DIPPR procedure 9I (volume-fraction + harmonic mean)
};

enum class LiquidViscosityMixing {
    Log,          // ln(mu) = Σ x_i ln(mu_i)
    SimpleLinear  // mu = Σ x_i mu_i
};

struct DipprLiquidRequirements {
    bool require_liquid_density = true;
    bool require_liquid_heat_capacity = false; // needed for H/S
    bool require_liquid_viscosity = false;
    bool require_liquid_thermal_conductivity = false;
};

class PropertyEngine {
public:
    PropertyEngine(EOSConstPtr eos, Core::Mixture mixture);

    void setMethod(PropertyMethod method) { method_ = method; }
    PropertyMethod method() const { return method_; }

    void setLiquidThermalConductivityMixing(LiquidThermalConductivityMixing mixing) { k_mixing_ = mixing; }
    LiquidThermalConductivityMixing liquidThermalConductivityMixing() const { return k_mixing_; }

    void setLiquidViscosityMixing(LiquidViscosityMixing mixing) { mu_mixing_ = mixing; }
    LiquidViscosityMixing liquidViscosityMixing() const { return mu_mixing_; }

    void setLiquidKHighPressureCorrection9G1Enabled(bool enabled) { enable_k_9g1_ = enabled; }
    bool liquidKHighPressureCorrection9G1Enabled() const { return enable_k_9g1_; }

    // EOS-based thermodynamics at (T, rho)
    Properties::ThermoProperties calculateTV(double T, double rho, const std::vector<double>& x) const;

    // DIPPR-backed pure-liquid properties
    Properties::ThermoProperties calculatePureLiquidDIPPR(
        const Data::Databanks& databanks,
        const std::string& identifier,
        double T,
        double P,
        const DipprLiquidRequirements& req = DipprLiquidRequirements{},
        Data::Diagnostics* diag = nullptr
    ) const;

    // DIPPR-backed pure-component scalar correlations (no phase routing).
    double calculatePureVaporPressureDIPPR(
        const Data::Databanks& databanks,
        const std::string& identifier,
        double T,
        Data::Diagnostics* diag = nullptr
    ) const;

    double calculatePureHeatOfVaporizationDIPPR(
        const Data::Databanks& databanks,
        const std::string& identifier,
        double T,
        Data::Diagnostics* diag = nullptr
    ) const;

    // DIPPR-backed mixture liquid properties using ideal mixing of pure-component DIPPR correlations.
    Properties::ThermoProperties calculateLiquidDIPPRIdealMix(
        const Data::Databanks& databanks,
        double T,
        double P,
        const std::vector<double>& x,
        const DipprLiquidRequirements& req = DipprLiquidRequirements{},
        Data::Diagnostics* diag = nullptr
    ) const;

private:
    EOSConstPtr eos_;
    Core::Mixture mixture_;
    PropertyMethod method_ = PropertyMethod::EOS_ONLY;
    LiquidThermalConductivityMixing k_mixing_ = LiquidThermalConductivityMixing::Dippr9I;
    LiquidViscosityMixing mu_mixing_ = LiquidViscosityMixing::Log;
    bool enable_k_9g1_ = true;
};

} // namespace Engine
} // namespace DMThermo

#endif // THERMO_ENGINE_PROPERTY_ENGINE_H
