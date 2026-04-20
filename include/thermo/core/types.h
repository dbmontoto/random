/**
 * @file types.h
 * @brief Core type definitions for the DMThermo namespace
 *
 * Contains enums, type aliases, and simple data structures used
 * throughout the thermodynamic library.
 */

#ifndef THERMO_CORE_TYPES_H
#define THERMO_CORE_TYPES_H

#include <vector>
#include <string>
#include <optional>

namespace DMThermo {

/**
 * @brief Enumeration of thermodynamic phase types
 */
enum class PhaseType {
    Unknown,
    Vapor,
    Liquid,
    Liquid1,     ///< First liquid phase (in LLE/VLLE)
    Liquid2,     ///< Second liquid phase (in LLE/VLLE)
    Supercritical,
    Solid
};

/**
 * @brief Convert PhaseType to string representation
 */
inline std::string phaseTypeToString(PhaseType phase) {
    switch (phase) {
        case PhaseType::Unknown: return "Unknown";
        case PhaseType::Vapor: return "Vapor";
        case PhaseType::Liquid: return "Liquid";
        case PhaseType::Liquid1: return "Liquid1";
        case PhaseType::Liquid2: return "Liquid2";
        case PhaseType::Supercritical: return "Supercritical";
        case PhaseType::Solid: return "Solid";
        default: return "Unknown";
    }
}

/**
 * @brief Result structure from EOS calculations
 *
 * Contains all thermodynamic properties calculated by an equation of state.
 * All derivative fields are optional and only populated when requested.
 */
struct EOSResult {
    /// Reduced residual Helmholtz energy a_res/(RT) [-]
    double a_residual = 0.0;

    /// Pressure [Pa]
    double pressure = 0.0;

    /// Compressibility factor Z = PV/(nRT) [-]
    double compressibility = 0.0;

    /// Natural logarithm of fugacity coefficients [-]
    std::vector<double> ln_phi;

    /// Fugacity coefficients (exp of ln_phi) [-]
    std::vector<double> phi;

    // Optional derivatives (only computed when requested)

    /// d(a_res/RT)/dT at constant rho, x [1/K]
    std::optional<double> da_dT;

    /// d(a_res/RT)/drho at constant T, x [m³/mol]
    std::optional<double> da_drho;

    /// d²(a_res/RT)/dT² [1/K²]
    std::optional<double> d2a_dT2;

    /// d²(a_res/RT)/drho² [(m³/mol)²]
    std::optional<double> d2a_drho2;

    /// d²(a_res/RT)/(dT drho) [m³/(mol·K)]
    std::optional<double> d2a_dTdrho;

    /// d(a_res/RT)/dx_i at constant T, rho [-]
    std::optional<std::vector<double>> da_dx;

    /// Residual chemical potential (mu_i^res/RT) [-].
    ///
    /// For EOS implementations providing fugacity coefficients:
    ///   mu_i^res/(RT) = ln(phi_i)
    ///
    /// Total chemical potentials require an ideal/standard-state model and are computed
    /// outside the EOS layer (e.g., in reaction thermochemistry tooling).
    std::optional<std::vector<double>> chemical_potential;

    /// Whether calculation was successful
    bool success = true;

    /// Error message if calculation failed
    std::string error_message;
};

/**
 * @brief State of a single thermodynamic phase
 */
struct PhaseState {
    /// Phase type identifier
    PhaseType type = PhaseType::Unknown;

    /// Molar density [mol/m³]
    double density = 0.0;

    /// Molar volume [m³/mol]
    double molarVolume() const { return (density > 0) ? 1.0 / density : 0.0; }

    /// Compressibility factor [-]
    double compressibility = 0.0;

    /// Mole fractions in this phase
    std::vector<double> x;

    /// Fugacity coefficients in this phase
    std::vector<double> phi;

    /// Phase fraction (moles in this phase / total moles) [-]
    double fraction = 0.0;

    /// Residual enthalpy [J/mol]
    std::optional<double> H_res;

    /// Residual entropy [J/(mol·K)]
    std::optional<double> S_res;

    /// Residual Gibbs energy [J/mol]
    std::optional<double> G_res;

    /// Total enthalpy [J/mol] (relative to the library reference state)
    std::optional<double> H;

    /// Total entropy [J/(mol*K)] (relative to the library reference state)
    std::optional<double> S;

    /// Total Gibbs energy [J/mol] (relative to the library reference state)
    std::optional<double> G;
};

/**
 * @brief Derivative specification for EOS calculations
 */
struct DerivativeSpec {
    bool temperature = false;      ///< Compute T derivatives
    bool density = false;          ///< Compute rho derivatives
    bool composition = false;      ///< Compute x derivatives
    bool second_order = false;     ///< Include second derivatives

    /// Create spec for no derivatives
    static DerivativeSpec none() { return DerivativeSpec{}; }

    /// Create spec for all first derivatives
    static DerivativeSpec first() {
        return DerivativeSpec{true, true, true, false};
    }

    /// Create spec for all derivatives
    static DerivativeSpec all() {
        return DerivativeSpec{true, true, true, true};
    }
};

} // namespace DMThermo

#endif // THERMO_CORE_TYPES_H
