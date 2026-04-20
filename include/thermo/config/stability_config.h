/**
 * @file stability_config.h
 * @brief Configuration structures for phase stability analysis
 */

#ifndef THERMO_CONFIG_STABILITY_CONFIG_H
#define THERMO_CONFIG_STABILITY_CONFIG_H

#include <vector>
#include "phase_detection.h"

namespace DMThermo {
namespace Config {

/**
 * @brief Method for stability test
 */
enum class StabilityMethod {
    TPD,            ///< Tangent Plane Distance minimization
    Spinodal,       ///< Spinodal criterion (dP/drho < 0)
    Combined        ///< Unstable if either TPD or spinodal criterion indicates instability
};

/**
 * @brief Trial composition generation method
 */
enum class TrialMethod {
    Wilson,         ///< K-values from Wilson correlation
    IdealSolution,  ///< Pure component limits
    Random,         ///< Random compositions
    Custom          ///< User-provided trials
};

/**
 * @brief Configuration for phase stability analysis
 */
struct StabilityConfig {
    // Convergence criteria
    double tpd_tolerance = 1e-10;       ///< TPD convergence tolerance
    double spinodal_tolerance = 1e-8;   ///< Spinodal detection tolerance
    int max_iterations = 100;           ///< Maximum TPD iterations

    // Method selection
    StabilityMethod method = StabilityMethod::Combined;
    TrialMethod trial_method = TrialMethod::Wilson;

    // Trial composition settings
    int num_trial_compositions = 5;     ///< Number of trial compositions
    std::vector<std::vector<double>> custom_trials;  ///< Custom trial compositions

    // Stability thresholds
    double tpd_instability_threshold = -1e-8;  ///< TPD below this = unstable
    bool use_relative_tpd = true;       ///< Scale TPD by feed reduced Gibbs magnitude before thresholding

    // Spinodal options
    bool use_spinodal_filtering = true; ///< Use spinodal information to guide/filter TP density root selection
    bool detect_spinodal_gap = true;    ///< Detect spinodal points

    // Density root finding
    double density_search_precision = 1e-8;  ///< Density solver tolerance
    int max_density_roots = 5;          ///< Maximum TP density roots retained during stability analysis

    // Phase detection / density selection
    PhaseDetection phase_detection_space = PhaseDetection::TP;
    int max_density_newton_iters = 40;   ///< Used when phase_detection_space == TRho
    double density_newton_tol = 1e-10;   ///< |P(T,rho,x)-P| tolerance [Pa] when refining rho

    // Output control
    bool verbose = false;               ///< Print analysis details
    bool store_all_trials = false;      ///< Keep results for all trials

    // =========================================================================
    // Factory Methods
    // =========================================================================

    /**
     * @brief Default configuration
     */
    static StabilityConfig defaults() {
        return StabilityConfig{};
    }

    /**
     * @brief Explicit TP-root-search density selection (robust cold-start default).
     */
    static StabilityConfig forTP() {
        StabilityConfig config;
        config.phase_detection_space = PhaseDetection::TP;
        return config;
    }

    /**
     * @brief Explicit TRho density tracking (Newton refinement in rho; prefers a good rho guess).
     */
    static StabilityConfig forTRho() {
        StabilityConfig config;
        config.phase_detection_space = PhaseDetection::TRho;
        return config;
    }

    /**
     * @brief High accuracy configuration
     */
    static StabilityConfig highAccuracy() {
        StabilityConfig config;
        config.tpd_tolerance = 1e-14;
        config.num_trial_compositions = 10;
        config.max_iterations = 200;
        return config;
    }

    /**
     * @brief Configuration for LLE (liquid-liquid) systems
     */
    static StabilityConfig forLLE() {
        StabilityConfig config;
        config.num_trial_compositions = 15;
        config.trial_method = TrialMethod::IdealSolution;
        config.use_spinodal_filtering = true;
        return config;
    }

    /**
     * @brief Configuration for near-critical systems
     */
    static StabilityConfig nearCritical() {
        StabilityConfig config;
        config.tpd_tolerance = 1e-12;
        config.spinodal_tolerance = 1e-10;
        config.use_spinodal_filtering = true;
        config.detect_spinodal_gap = true;
        return config;
    }
};

/**
 * @brief Configuration for density root finding
 */
struct DensityRootConfig {
    double tolerance = 1e-10;           ///< Root-finding tolerance
    int max_iterations = 50;            ///< Maximum iterations
    double min_density = 1e-6;          ///< Minimum search density [mol/m³]
    double max_density = 50000.0;       ///< Maximum search density [mol/m³]
    bool use_spinodal_bounds = true;    ///< Bound search by spinodal
    bool filter_metastable = true;      ///< Exclude mechanically unstable roots (dP/drho <= 0)
    int max_roots = 0;                  ///< Maximum roots retained after filtering (<=0 means unlimited)
};

/**
 * @brief Configuration for spinodal detection
 */
struct SpinodalConfig {
    double tolerance = 1e-8;            ///< Detection tolerance
    int max_iterations = 50;            ///< Maximum iterations
    double initial_step = 100.0;        ///< Initial density step [mol/m³]
    bool find_both_sides = true;        ///< Find both vapor and liquid spinodals
};

} // namespace Config
} // namespace DMThermo

#endif // THERMO_CONFIG_STABILITY_CONFIG_H
