/**
 * @file units.cpp
 * @brief Unit tags and conversion helpers for thermodynamic properties.
 */

#include "thermo/core/units.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace DMThermo {
namespace Core {
namespace Units {

namespace {

constexpr double K_ATM_PA = 101325.0;
constexpr double K_PSI_TO_PA = 6894.757293168;
constexpr double K_LB_TO_KG = 0.45359237;
constexpr double K_LB_MOL_TO_MOL = 453.59237;
constexpr double K_FT3_TO_M3 = 0.028316846592;
constexpr double K_US_GAL_TO_M3 = 0.003785411784;
constexpr double K_LB_PER_FT3_TO_KG_PER_M3 = 16.01846337396014;
constexpr double K_BTU_PER_HR_FT_F_TO_W_PER_M_K = 1.730735;
constexpr double K_BTU_TO_J = 1055.05585262;
constexpr double K_KGF_M_TO_J = 9.80665;

enum class UnitDimension {
    Unknown = 0,
    Dimensionless,
    Temperature,
    Pressure,
    MolarVolume,
    MolarDensity,
    MassDensity,
    MolarFlow,
    MassFlow,
    MolarEnergy,
    MolarEntropy,
    BulkPropertyPerTemperature,
    SpecificHeatMass,
    ThermalConductivity
};

enum class PressureReference {
    None = 0,
    Absolute,
    Gauge,
    Differential
};

struct UnitDef {
    Unit unit;
    const char* token;
    UnitDimension dimension;
    double scale_to_base;
    double offset_to_base;
    PressureReference pressure_reference;
};

constexpr std::array<UnitDef, 67> kUnitDefs = {
    UnitDef{Unit::UNKNOWN, "UNKNOWN", UnitDimension::Unknown, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::DIMENSIONLESS, "DIMENSIONLESS", UnitDimension::Dimensionless, 1.0, 0.0, PressureReference::None},

    // Temperature (base: K)
    UnitDef{Unit::K, "K", UnitDimension::Temperature, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::DEG_C, "DEG_C", UnitDimension::Temperature, 1.0, 273.15, PressureReference::None},
    UnitDef{Unit::DEG_F, "DEG_F", UnitDimension::Temperature, 5.0 / 9.0, 255.3722222222222, PressureReference::None},
    UnitDef{Unit::R, "R", UnitDimension::Temperature, 5.0 / 9.0, 0.0, PressureReference::None},

    // Pressure (base: Pa with explicit abs/gauge/differential semantics)
    UnitDef{Unit::PA, "PA", UnitDimension::Pressure, 1.0, 0.0, PressureReference::Differential},
    UnitDef{Unit::PA_ABS, "PA_ABS", UnitDimension::Pressure, 1.0, 0.0, PressureReference::Absolute},
    UnitDef{Unit::PA_GAUGE, "PA_GAUGE", UnitDimension::Pressure, 1.0, K_ATM_PA, PressureReference::Gauge},
    UnitDef{Unit::KPA, "KPA", UnitDimension::Pressure, 1.0e3, 0.0, PressureReference::Differential},
    UnitDef{Unit::KPA_ABS, "KPA_ABS", UnitDimension::Pressure, 1.0e3, 0.0, PressureReference::Absolute},
    UnitDef{Unit::KPA_GAUGE, "KPA_GAUGE", UnitDimension::Pressure, 1.0e3, K_ATM_PA, PressureReference::Gauge},
    UnitDef{Unit::MPA, "MPA", UnitDimension::Pressure, 1.0e6, 0.0, PressureReference::Differential},
    UnitDef{Unit::MPA_ABS, "MPA_ABS", UnitDimension::Pressure, 1.0e6, 0.0, PressureReference::Absolute},
    UnitDef{Unit::MPA_GAUGE, "MPA_GAUGE", UnitDimension::Pressure, 1.0e6, K_ATM_PA, PressureReference::Gauge},
    UnitDef{Unit::ATM_ABS, "ATM_ABS", UnitDimension::Pressure, K_ATM_PA, 0.0, PressureReference::Absolute},
    UnitDef{Unit::ATM_GAUGE, "ATM_GAUGE", UnitDimension::Pressure, K_ATM_PA, K_ATM_PA, PressureReference::Gauge},
    UnitDef{Unit::BAR, "BAR", UnitDimension::Pressure, 1.0e5, 0.0, PressureReference::Absolute},
    UnitDef{Unit::BARG, "BARG", UnitDimension::Pressure, 1.0e5, K_ATM_PA, PressureReference::Gauge},
    UnitDef{Unit::PSI, "PSI", UnitDimension::Pressure, K_PSI_TO_PA, 0.0, PressureReference::Differential},
    UnitDef{Unit::PSIA, "PSIA", UnitDimension::Pressure, K_PSI_TO_PA, 0.0, PressureReference::Absolute},
    UnitDef{Unit::PSIG, "PSIG", UnitDimension::Pressure, K_PSI_TO_PA, K_ATM_PA, PressureReference::Gauge},

    // Molar volume (base: m^3/mol)
    UnitDef{Unit::M3_PER_MOL, "M3_PER_MOL", UnitDimension::MolarVolume, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::M3_PER_KMOL, "M3_PER_KMOL", UnitDimension::MolarVolume, 1.0e-3, 0.0, PressureReference::None},
    UnitDef{Unit::L_PER_MOL, "L_PER_MOL", UnitDimension::MolarVolume, 1.0e-3, 0.0, PressureReference::None},
    UnitDef{Unit::CM3_PER_MOL, "CM3_PER_MOL", UnitDimension::MolarVolume, 1.0e-6, 0.0, PressureReference::None},
    UnitDef{Unit::FT3_PER_LB_MOL, "FT3_PER_LB_MOL", UnitDimension::MolarVolume, K_FT3_TO_M3 / K_LB_MOL_TO_MOL, 0.0, PressureReference::None},

    // Molar density (base: mol/m^3)
    UnitDef{Unit::MOL_PER_M3, "MOL_PER_M3", UnitDimension::MolarDensity, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::MOL_PER_L, "MOL_PER_L", UnitDimension::MolarDensity, 1.0e3, 0.0, PressureReference::None},
    UnitDef{Unit::KMOL_PER_M3, "KMOL_PER_M3", UnitDimension::MolarDensity, 1.0e3, 0.0, PressureReference::None},
    UnitDef{Unit::LB_MOL_PER_FT3, "LB_MOL_PER_FT3", UnitDimension::MolarDensity, K_LB_MOL_TO_MOL / K_FT3_TO_M3, 0.0, PressureReference::None},

    // Mass density (base: kg/m^3)
    UnitDef{Unit::KG_PER_M3, "KG_PER_M3", UnitDimension::MassDensity, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::KG_PER_L, "KG_PER_L", UnitDimension::MassDensity, 1.0e3, 0.0, PressureReference::None},
    UnitDef{Unit::G_PER_CM3, "G_PER_CM3", UnitDimension::MassDensity, 1.0e3, 0.0, PressureReference::None},
    UnitDef{Unit::LB_PER_FT3, "LB_PER_FT3", UnitDimension::MassDensity, K_LB_PER_FT3_TO_KG_PER_M3, 0.0, PressureReference::None},
    UnitDef{Unit::LB_PER_GAL, "LB_PER_GAL", UnitDimension::MassDensity, K_LB_TO_KG / K_US_GAL_TO_M3, 0.0, PressureReference::None},

    // Molar flow (base: mol/s)
    UnitDef{Unit::MOL_PER_S, "MOL_PER_S", UnitDimension::MolarFlow, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::MOL_PER_HR, "MOL_PER_HR", UnitDimension::MolarFlow, 1.0 / 3600.0, 0.0, PressureReference::None},
    UnitDef{Unit::KMOL_PER_S, "KMOL_PER_S", UnitDimension::MolarFlow, 1.0e3, 0.0, PressureReference::None},
    UnitDef{Unit::KMOL_PER_HR, "KMOL_PER_HR", UnitDimension::MolarFlow, 1.0e3 / 3600.0, 0.0, PressureReference::None},
    UnitDef{Unit::LB_MOL_PER_S, "LB_MOL_PER_S", UnitDimension::MolarFlow, K_LB_MOL_TO_MOL, 0.0, PressureReference::None},
    UnitDef{Unit::LB_MOL_PER_HR, "LB_MOL_PER_HR", UnitDimension::MolarFlow, K_LB_MOL_TO_MOL / 3600.0, 0.0, PressureReference::None},

    // Mass flow (base: kg/s)
    UnitDef{Unit::KG_PER_S, "KG_PER_S", UnitDimension::MassFlow, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::KG_PER_HR, "KG_PER_HR", UnitDimension::MassFlow, 1.0 / 3600.0, 0.0, PressureReference::None},
    UnitDef{Unit::LB_PER_S, "LB_PER_S", UnitDimension::MassFlow, K_LB_TO_KG, 0.0, PressureReference::None},
    UnitDef{Unit::LB_PER_HR, "LB_PER_HR", UnitDimension::MassFlow, K_LB_TO_KG / 3600.0, 0.0, PressureReference::None},

    // Molar energy (base: J/mol)
    UnitDef{Unit::J_PER_MOL, "J_PER_MOL", UnitDimension::MolarEnergy, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::KJ_PER_MOL, "KJ_PER_MOL", UnitDimension::MolarEnergy, 1.0e3, 0.0, PressureReference::None},
    UnitDef{Unit::J_PER_KMOL, "J_PER_KMOL", UnitDimension::MolarEnergy, 1.0e-3, 0.0, PressureReference::None},
    UnitDef{Unit::KJ_PER_KMOL, "KJ_PER_KMOL", UnitDimension::MolarEnergy, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::BTU_PER_MOL, "BTU_PER_MOL", UnitDimension::MolarEnergy, K_BTU_TO_J, 0.0, PressureReference::None},
    UnitDef{Unit::BTU_PER_LB_MOL, "BTU_PER_LB_MOL", UnitDimension::MolarEnergy, K_BTU_TO_J / K_LB_MOL_TO_MOL, 0.0, PressureReference::None},

    // Molar entropy / molar heat capacity (base: J/mol/K)
    UnitDef{Unit::J_PER_MOL_K, "J_PER_MOL_K", UnitDimension::MolarEntropy, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::KJ_PER_MOL_K, "KJ_PER_MOL_K", UnitDimension::MolarEntropy, 1.0e3, 0.0, PressureReference::None},
    UnitDef{Unit::J_PER_KMOL_K, "J_PER_KMOL_K", UnitDimension::MolarEntropy, 1.0e-3, 0.0, PressureReference::None},
    UnitDef{Unit::KJ_PER_KMOL_K, "KJ_PER_KMOL_K", UnitDimension::MolarEntropy, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::BTU_PER_MOL_R, "BTU_PER_MOL_R", UnitDimension::MolarEntropy, K_BTU_TO_J * (9.0 / 5.0), 0.0, PressureReference::None},
    UnitDef{Unit::BTU_PER_LB_MOL_R, "BTU_PER_LB_MOL_R", UnitDimension::MolarEntropy, (K_BTU_TO_J * (9.0 / 5.0)) / K_LB_MOL_TO_MOL, 0.0, PressureReference::None},

    // Volumetric property per temperature (base: kg/m^3/K)
    UnitDef{Unit::KG_PER_M3_K, "KG_PER_M3_K", UnitDimension::BulkPropertyPerTemperature, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::LB_PER_FT3_R, "LB_PER_FT3_R", UnitDimension::BulkPropertyPerTemperature,
            K_LB_PER_FT3_TO_KG_PER_M3 * (9.0 / 5.0), 0.0, PressureReference::None},

    // Specific heat capacity (base: J/kg/K)
    UnitDef{Unit::J_PER_KG_K, "J_PER_KG_K", UnitDimension::SpecificHeatMass, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::KJ_PER_KG_K, "KJ_PER_KG_K", UnitDimension::SpecificHeatMass, 1.0e3, 0.0, PressureReference::None},
    UnitDef{Unit::KGF_M_PER_KG_K, "KGF_M_PER_KG_K", UnitDimension::SpecificHeatMass, K_KGF_M_TO_J, 0.0, PressureReference::None},
    UnitDef{Unit::BTU_PER_LB_R, "BTU_PER_LB_R", UnitDimension::SpecificHeatMass, K_BTU_TO_J / (K_LB_TO_KG * (5.0 / 9.0)), 0.0, PressureReference::None},

    // Thermal conductivity (base: W/m/K)
    UnitDef{Unit::W_PER_M_K, "W_PER_M_K", UnitDimension::ThermalConductivity, 1.0, 0.0, PressureReference::None},
    UnitDef{Unit::BTU_PER_HR_FT_F, "BTU_PER_HR_FT_F", UnitDimension::ThermalConductivity, K_BTU_PER_HR_FT_F_TO_W_PER_M_K, 0.0, PressureReference::None},
    UnitDef{Unit::BTU_PER_HR_FT_R, "BTU_PER_HR_FT_R", UnitDimension::ThermalConductivity, K_BTU_PER_HR_FT_F_TO_W_PER_M_K, 0.0, PressureReference::None},
};

const UnitDef* findUnitDef(Unit unit) {
    for (const auto& def : kUnitDefs) {
        if (def.unit == unit) {
            return &def;
        }
    }
    return nullptr;
}

bool canConvertPressure(PressureReference from_ref, PressureReference to_ref) {
    if (from_ref == PressureReference::Differential || to_ref == PressureReference::Differential) {
        return from_ref == PressureReference::Differential && to_ref == PressureReference::Differential;
    }
    return true;
}

} // namespace

const char* toString(Unit unit) {
    const UnitDef* def = findUnitDef(unit);
    return def ? def->token : "UNKNOWN";
}

Unit fromString(const std::string& token) {
    static const std::unordered_map<std::string, Unit> kCanonicalTokens = [] {
        std::unordered_map<std::string, Unit> tokens;
        tokens.reserve(kUnitDefs.size());
        for (const auto& def : kUnitDefs) {
            tokens.emplace(def.token, def.unit);
        }
        return tokens;
    }();

    const auto it = kCanonicalTokens.find(token);
    if (it != kCanonicalTokens.end()) {
        return it->second;
    }
    return Unit::UNKNOWN;
}

bool UnitConverter::canConvert(Unit from, Unit to) {
    if (from == to) return true;

    const UnitDef* from_def = findUnitDef(from);
    const UnitDef* to_def = findUnitDef(to);
    if (!from_def || !to_def) return false;

    if (from_def->dimension == UnitDimension::Unknown || to_def->dimension == UnitDimension::Unknown) return false;
    if (from_def->dimension != to_def->dimension) return false;

    if (from_def->dimension == UnitDimension::Pressure) {
        return canConvertPressure(from_def->pressure_reference, to_def->pressure_reference);
    }

    return true;
}

bool UnitConverter::areCompatible(Unit from, Unit to) {
    return canConvert(from, to);
}

double UnitConverter::convert(double value, Unit from, Unit to) {
    if (from == to) return value;
    if (!std::isfinite(value)) return value;

    if (!canConvert(from, to)) {
        throw std::invalid_argument(
            std::string("Incompatible unit conversion: ") + toString(from) + " -> " + toString(to));
    }

    const UnitDef* from_def = findUnitDef(from);
    const UnitDef* to_def = findUnitDef(to);
    if (!from_def || !to_def) {
        throw std::invalid_argument(
            std::string("Unknown unit conversion: ") + toString(from) + " -> " + toString(to));
    }

    const double base_value = value * from_def->scale_to_base + from_def->offset_to_base;
    const double out = (base_value - to_def->offset_to_base) / to_def->scale_to_base;
    if (!std::isfinite(out)) {
        throw std::runtime_error(
            std::string("Non-finite conversion result: ") + toString(from) + " -> " + toString(to));
    }
    return out;
}

} // namespace Units
} // namespace Core
} // namespace DMThermo
