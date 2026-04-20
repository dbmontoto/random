/**
 * @file records.h
 * @brief Typed records parsed from databank CSV files
 */

#ifndef THERMO_DATA_RECORDS_H
#define THERMO_DATA_RECORDS_H

#include <optional>
#include <string>

namespace DMThermo {
namespace Data {

struct ChemicalKey {
    std::optional<int> chemid;
    std::string cas;
    std::string name;
};

struct Correlation {
    int eq_form = 0; // DIPPR equation form number (future-proof: not a closed enum)
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double D = 0.0;
    double E = 0.0;
    double F = 0.0;
    double G = 0.0;
    std::optional<double> Tmin;
    std::optional<double> Tmax;
    std::string data_type;
    std::string error;
};

struct DipprRecord {
    ChemicalKey key;

    // Common constants (optional because databank may be incomplete)
    std::optional<double> MW;
    std::optional<double> Tc;
    std::optional<double> Pc;
    std::optional<double> Vc;
    std::optional<double> Zc;
    std::optional<double> omega;

    // Standard-state / formation properties (optional; units per CSV header key)
    // Note: DIPPR databank uses J/kmol and J/kmol/K for these columns.
    std::optional<double> HFOR; // Formation enthalpy [J/kmol]
    std::optional<double> GFOR; // Gibbs energy of formation [J/kmol]
    std::optional<double> ENT;  // Entropy [J/kmol/K]
    std::optional<double> HSTD; // Standard enthalpy [J/kmol]
    std::optional<double> SSTD; // Standard entropy [J/kmol/K]

    // Correlation groups (optional)
    std::optional<Correlation> liquid_density;         // LDN_*
    std::optional<Correlation> vapor_pressure;         // VP_*
    std::optional<Correlation> liquid_viscosity;       // LVS_*
    std::optional<Correlation> heat_of_vaporization;   // HVP/HVAP_*
    std::optional<Correlation> liquid_thermal_cond;    // LTC_*
    std::optional<Correlation> liquid_heat_capacity;   // LCP_*
    std::optional<Correlation> vapor_viscosity;        // VVS_*
    std::optional<Correlation> ideal_gas_cp;           // ICP_*
    std::optional<Correlation> surface_tension;        // ST_*
};

struct PCSaftRecord {
    ChemicalKey key;

    std::optional<double> MW;
    std::optional<double> m;
    std::optional<double> sigma;
    std::optional<double> epsilon_k;
    std::optional<double> kappa_ab;
    std::optional<double> epsilon_k_ab;
    std::optional<std::string> assoc_scheme;
    std::optional<int> assoc_num_sites;
    std::optional<double> Tmin;
    std::optional<double> Tmax;
};

struct BipRecord {
    ChemicalKey i;
    ChemicalKey j;
    std::optional<double> kij_pr;
    std::optional<double> kij_vtpr;
    std::optional<double> kij_pcsaft;
    std::optional<double> Tmin;
    std::optional<double> Tmax;
};

} // namespace Data
} // namespace DMThermo

#endif // THERMO_DATA_RECORDS_H
