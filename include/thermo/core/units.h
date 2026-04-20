/**
 * @file units.h
 * @brief Unit tags and conversion helpers for thermodynamic properties.
 */

#ifndef THERMO_CORE_UNITS_H
#define THERMO_CORE_UNITS_H

#include <string>

namespace DMThermo {
namespace Core {
namespace Units {

enum class Unit {
    UNKNOWN = 0,
    DIMENSIONLESS,

    // Temperature
    K,
    DEG_C,
    DEG_F,
    R,

    // Pressure
    PA,       // differential
    PA_ABS,
    PA_GAUGE,

    KPA,      // differential
    KPA_ABS,
    KPA_GAUGE,

    MPA,      // differential
    MPA_ABS,
    MPA_GAUGE,

    ATM_ABS,
    ATM_GAUGE,

    // Absolute / gauge pair retained for bar by industry convention
    BAR,      // absolute
    BARG,     // gauge

    // PSI conventions: PSI (differential), PSIA/PSIG (absolute/gauge)
    PSI,      // differential
    PSIA,
    PSIG,

    // Molar volume
    M3_PER_MOL,
    M3_PER_KMOL,
    L_PER_MOL,
    CM3_PER_MOL,
    FT3_PER_LB_MOL,

    // Molar density
    MOL_PER_M3,
    MOL_PER_L,
    KMOL_PER_M3,
    LB_MOL_PER_FT3,

    // Mass density
    KG_PER_M3,
    KG_PER_L,
    G_PER_CM3,
    LB_PER_FT3,
    LB_PER_GAL,

    // Molar flow
    MOL_PER_S,
    MOL_PER_HR,
    KMOL_PER_S,
    KMOL_PER_HR,
    LB_MOL_PER_S,
    LB_MOL_PER_HR,

    // Mass flow
    KG_PER_S,
    KG_PER_HR,
    LB_PER_S,
    LB_PER_HR,

    // Molar energy
    J_PER_MOL,
    KJ_PER_MOL,
    J_PER_KMOL,
    KJ_PER_KMOL,
    BTU_PER_MOL,
    BTU_PER_LB_MOL,

    // Molar entropy / molar heat capacity
    J_PER_MOL_K,
    KJ_PER_MOL_K,
    J_PER_KMOL_K,
    KJ_PER_KMOL_K,
    BTU_PER_MOL_R,
    BTU_PER_LB_MOL_R,

    // Volumetric property per temperature
    KG_PER_M3_K,
    LB_PER_FT3_R,

    // Specific heat capacity (mass basis)
    J_PER_KG_K,
    KJ_PER_KG_K,
    KGF_M_PER_KG_K,
    BTU_PER_LB_R,

    // Thermal conductivity
    W_PER_M_K,
    BTU_PER_HR_FT_F,
    BTU_PER_HR_FT_R
};

/**
 * @brief Convert unit enum to stable string token (e.g. "M3_PER_KMOL").
 */
const char* toString(Unit unit);

/**
 * @brief Parse stable string token into Unit (unknown => Unit::UNKNOWN).
 */
Unit fromString(const std::string& token);

/**
 * @brief Runtime unit conversion helper.
 */
class UnitConverter {
public:
    static bool canConvert(Unit from, Unit to);
    static bool areCompatible(Unit from, Unit to);
    static double convert(double value, Unit from, Unit to);
};

/**
 * @brief Value + unit wrapper with explicit conversion.
 */
class PropertyVariable {
public:
    PropertyVariable() = default;
    PropertyVariable(double value, Unit unit)
        : value_(value)
        , unit_(unit)
    {
    }

    double value() const { return value_; }
    Unit unit() const { return unit_; }

    double as(Unit target_unit) const {
        return UnitConverter::convert(value_, unit_, target_unit);
    }

private:
    double value_ = 0.0;
    Unit unit_ = Unit::UNKNOWN;
};

} // namespace Units
} // namespace Core
} // namespace DMThermo

#endif // THERMO_CORE_UNITS_H
