#include "thermo/data/dippr_adapter.h"

#include <cmath>
#include <limits>
#include <string>

namespace DMThermo {
namespace Data {

namespace {

bool applyCorrelation(
    const std::optional<Correlation>& src,
    DIPPR::Coefficients& dst,
    bool& dst_has_flag,
    const std::string& label,
    const std::string& record_label,
    Diagnostics* diag)
{
    if (!src) return false;

    const auto eq = DIPPR::equationTypeFromForm(src->eq_form);
    if (!eq) {
        if (diag) {
            diag->warn(
                "Unsupported DIPPR equation form " + std::to_string(src->eq_form) + " for " + label + " in record '" + record_label + "'; skipping"
            );
        }
        return false;
    }

    dst = DIPPR::Coefficients(
        *eq,
        src->A, src->B, src->C, src->D, src->E, src->F, src->G,
        src->Tmin.value_or(std::numeric_limits<double>::quiet_NaN()),
        src->Tmax.value_or(std::numeric_limits<double>::quiet_NaN())
    );
    dst_has_flag = dst.isValid();
    if (!dst_has_flag && diag) {
        diag->warn("Invalid DIPPR coefficients for " + label + " in record '" + record_label + "'; skipping");
    }
    return true;
}

void requireTcForEvaluation(
    const DIPPR::Coefficients& coeffs,
    bool& has_flag,
    double Tc,
    const std::string& label,
    const std::string& record_label,
    Diagnostics* diag)
{
    if (!has_flag) return;

    const bool needs_tc =
        coeffs.equation == DIPPR::EquationType::EQ_106 ||
        coeffs.equation == DIPPR::EquationType::EQ_114 ||
        coeffs.equation == DIPPR::EquationType::EQ_123;

    if (!needs_tc) return;

    if (!(std::isfinite(Tc) && Tc > 0.0)) {
        has_flag = false;
        if (diag) {
            diag->warn("Missing Tc in record '" + record_label + "'; cannot use " + label + " correlation " +
                       std::string(DIPPR::equationTypeName(coeffs.equation)));
        }
    }
}

} // namespace

std::optional<DIPPR::Parameters> toDipprParameters(const DipprRecord& record, Diagnostics* diag) {
    DIPPR::Parameters p;

    p.name = record.key.name;
    p.cas = record.key.cas;

    if (record.Tc) p.Tc = *record.Tc;
    if (record.Pc) p.Pc = *record.Pc;
    if (record.omega) p.omega = *record.omega;

    std::string record_label = record.key.name;
    if (!record.key.cas.empty()) {
        record_label += " (" + record.key.cas + ")";
    } else if (record.key.chemid.has_value()) {
        record_label += " (CHEMID " + std::to_string(record.key.chemid.value()) + ")";
    }

    // Map correlation groups supported by the current DIPPR module.
    applyCorrelation(record.liquid_density, p.liquid_density, p.has_liquid_density, "liquid density", record_label, diag);
    applyCorrelation(record.vapor_pressure, p.vapor_pressure, p.has_vapor_pressure, "vapor pressure", record_label, diag);
    applyCorrelation(record.heat_of_vaporization, p.heat_of_vaporization, p.has_heat_of_vaporization, "heat of vaporization", record_label, diag);
    applyCorrelation(record.liquid_heat_capacity, p.liquid_heat_capacity, p.has_liquid_heat_capacity, "liquid heat capacity", record_label, diag);
    applyCorrelation(record.liquid_viscosity, p.liquid_viscosity, p.has_liquid_viscosity, "liquid viscosity", record_label, diag);
    applyCorrelation(record.vapor_viscosity, p.vapor_viscosity, p.has_vapor_viscosity, "vapor viscosity", record_label, diag);
    applyCorrelation(record.liquid_thermal_cond, p.liquid_thermal_conductivity, p.has_liquid_thermal_conductivity, "liquid thermal conductivity", record_label, diag);
    applyCorrelation(record.surface_tension, p.surface_tension, p.has_surface_tension, "surface tension", record_label, diag);
    applyCorrelation(record.ideal_gas_cp, p.ideal_gas_cp, p.has_ideal_gas_cp, "ideal gas Cp", record_label, diag);

    // Some DIPPR equations require Tc at evaluation time; if Tc is missing, treat those correlations as unavailable.
    requireTcForEvaluation(p.heat_of_vaporization, p.has_heat_of_vaporization, p.Tc, "heat of vaporization", record_label, diag);
    requireTcForEvaluation(p.liquid_heat_capacity, p.has_liquid_heat_capacity, p.Tc, "liquid heat capacity", record_label, diag);
    requireTcForEvaluation(p.liquid_thermal_conductivity, p.has_liquid_thermal_conductivity, p.Tc, "liquid thermal conductivity", record_label, diag);
    requireTcForEvaluation(p.surface_tension, p.has_surface_tension, p.Tc, "surface tension", record_label, diag);

    // If nothing was mapped, return nullopt (not useful for DIPPR::Calculator).
    const bool any =
        p.has_liquid_density ||
        p.has_vapor_pressure ||
        p.has_heat_of_vaporization ||
        p.has_liquid_heat_capacity ||
        p.has_liquid_viscosity ||
        p.has_vapor_viscosity ||
        p.has_liquid_thermal_conductivity ||
        p.has_surface_tension ||
        p.has_ideal_gas_cp;

    if (!any) {
        if (diag) diag->warn("No supported DIPPR correlations available for record '" + record_label + "'");
        return std::nullopt;
    }

    return p;
}

} // namespace Data
} // namespace DMThermo
