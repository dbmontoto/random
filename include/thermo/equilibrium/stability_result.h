/**
 * @file stability_result.h
 * @brief Result structures for phase stability analysis
 */

#ifndef THERMO_EQUILIBRIUM_STABILITY_RESULT_H
#define THERMO_EQUILIBRIUM_STABILITY_RESULT_H

#include <vector>
#include <string>
#include <optional>
#include "thermo/core/types.h"

namespace DMThermo {
namespace Equilibrium {
namespace Stability {

/**
 * @brief Result from a single TPD minimization
 */
struct TPDTrialResult {
    /// Trial composition (initial)
    std::vector<double> w_initial;

    /// Converged composition
    std::vector<double> w_converged;

    /// TPD value at converged composition
    double tpd = 0.0;

    /// Did this trial converge?
    bool converged = false;

    /// Number of iterations
    int iterations = 0;

    /// Is this trial composition indicating instability?
    bool indicates_instability = false;
};

/**
 * @brief Complete result from stability analysis
 */
struct StabilityResult {
    // =========================================================================
    // Main Result
    // =========================================================================

    /// Is the phase stable?
    bool is_stable = true;

    /// Minimum TPD value found (negative = unstable)
    double min_tpd = 0.0;

    /// Trial composition that gave minimum TPD
    std::optional<std::vector<double>> unstable_composition;

    /// Predicted phase type of incipient phase
    std::optional<PhaseType> incipient_phase_type;

    // =========================================================================
    // Analysis Details
    // =========================================================================

    /// Number of TPD trials performed
    int num_trials = 0;

    /// Number of converged trials
    int num_converged = 0;

    /// Results from all trials (if requested)
    std::vector<TPDTrialResult> trial_results;

    /// Method used for analysis
    std::string method_used;

    // =========================================================================
    // Spinodal Information
    // =========================================================================

    /// Was spinodal analysis performed?
    bool spinodal_analyzed = false;

    /// Vapor-side spinodal density [mol/m³]
    std::optional<double> rho_spinodal_vapor;

    /// Liquid-side spinodal density [mol/m³]
    std::optional<double> rho_spinodal_liquid;

    /// Is the feed between spinodal points (absolutely unstable)?
    bool in_spinodal_region = false;

    // =========================================================================
    // Debug Information
    // =========================================================================

    /// Total computation time [ms]
    double computation_time_ms = 0.0;

    /// Any warnings or notes
    std::string notes;
};

/**
 * @brief Result from density root search
 */
struct DensityRootResult {
    /// Found density values [mol/m³]
    std::vector<double> densities;

    /// Phase classification for each root
    std::vector<PhaseType> phase_types;

    /// Mechanical stability (dP/drho > 0) for each root
    std::vector<bool> is_mechanically_stable;

    /// Gibbs energy for each root (for selecting equilibrium)
    std::vector<double> gibbs_energies;

    /// Index of thermodynamically stable root
    std::optional<int> stable_root_index;

    /// Number of roots found
    int num_roots = 0;
};

/**
 * @brief Result from spinodal detection
 */
struct SpinodalResult {
    /// Was spinodal detection successful?
    bool found = false;

    /// Vapor-side spinodal density [mol/m³]
    double rho_vapor_spinodal = 0.0;

    /// Liquid-side spinodal density [mol/m³]
    double rho_liquid_spinodal = 0.0;

    /// Pressure at spinodal points [Pa]
    double pressure_spinodal = 0.0;

    /// Temperature [K]
    double temperature = 0.0;

    /// Is there a spinodal gap (two-phase possible)?
    bool has_gap = false;
};

} // namespace Stability
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_STABILITY_RESULT_H
