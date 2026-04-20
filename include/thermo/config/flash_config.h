/**
 * @file flash_config.h
 * @brief Configuration structures for flash calculations
 *
 * Provides flexible configuration for flash solvers with sensible defaults.
 */

#ifndef THERMO_CONFIG_FLASH_CONFIG_H
#define THERMO_CONFIG_FLASH_CONFIG_H

#include <string>
#include <vector>
#include "phase_detection.h"
#include "tv_solve_strategy.h"

namespace DMThermo {
namespace Config {

/**
 * @brief Flash calculation method selection
 */
enum class FlashMethod {
    Auto,       ///< Automatically select based on conditions
    Gibbs,      ///< Gibbs energy minimization (successive substitution)
    Helmholtz,  ///< Helmholtz energy minimization (Newton-based)
    Combined    ///< Start with Gibbs, refine with Helmholtz
};

/**
 * @brief K-value initialization method
 */
enum class KValueInit {
    Wilson,         ///< Wilson correlation
    FromFugacity,   ///< Calculate from fugacity coefficients
    Custom          ///< User-provided initial K-values
};

/**
 * @brief Configuration for flash calculations
 */
struct FlashConfig {
    // Convergence criteria
    double tolerance = 1e-8;            ///< Convergence tolerance
    double rel_tolerance = 1e-8;        ///< Relative tolerance
    int max_iterations = 100;           ///< Maximum iterations

    // Method selection
    FlashMethod method = FlashMethod::Auto;
    KValueInit k_init = KValueInit::Wilson;
    std::vector<double> custom_k_values{}; ///< Used when k_init == Custom (size must equal nc)
    bool use_dippr_psat_k_init = false;    ///< If true, ThermoSystem may seed custom_k_values from DIPPR Psat/P

    // Stability test options
    bool perform_stability_test = true;  ///< Run stability test first
    int num_stability_trials = 5;        ///< Number of TPD trial compositions

    // Phase handling
    int max_phases = 2;                  ///< Maximum phases to consider (1, 2, or 3)

    // Phase detection / density selection
    PhaseDetection phase_detection_space = PhaseDetection::TP;
    int max_density_newton_iters = 40;   ///< Used when phase_detection_space == TRho
    double density_newton_tol = 1e-10;   ///< |P(T,rho,x)-P| tolerance [Pa] when refining rho

    // Acceleration options
    bool use_acceleration = true;        ///< Use GDEM/Wegstein acceleration
    double acceleration_factor = 0.8;    ///< Acceleration damping factor
    double max_acceleration = 5.0;       ///< Maximum acceleration factor

    // K-value bounds
    double min_k_value = 1e-6;          ///< Minimum K-value
    double max_k_value = 1e6;           ///< Maximum K-value

    // Near-critical handling
    bool handle_near_critical = true;   ///< Special handling near critical
    double critical_distance = 0.15;    ///< |ln(K)| threshold for near-critical

    // Output control
    bool verbose = false;               ///< Print iteration details
    bool compute_properties = false;    ///< Compute H, S in result

    // =========================================================================
    // Factory Methods for Common Configurations
    // =========================================================================

    /**
     * @brief Default configuration for general use
     */
    static FlashConfig defaults() {
        return FlashConfig{};
    }

    /**
     * @brief Explicit TP-root-search phase detection (robust cold-start default).
     */
    static FlashConfig forTP() {
        FlashConfig config;
        config.phase_detection_space = PhaseDetection::TP;
        return config;
    }

    /**
     * @brief Explicit TRho density tracking (Newton refinement in rho; prefers a good rho guess).
     */
    static FlashConfig forTRho() {
        FlashConfig config;
        config.phase_detection_space = PhaseDetection::TRho;
        return config;
    }

    /**
     * @brief High accuracy configuration (tighter tolerances)
     */
    static FlashConfig highAccuracy() {
        FlashConfig config;
        config.tolerance = 1e-12;
        config.rel_tolerance = 1e-10;
        config.max_iterations = 200;
        config.num_stability_trials = 10;
        return config;
    }

    /**
     * @brief Fast convergence configuration (looser tolerances)
     */
    static FlashConfig fastConvergence() {
        FlashConfig config;
        config.tolerance = 1e-6;
        config.rel_tolerance = 1e-5;
        config.max_iterations = 50;
        config.num_stability_trials = 3;
        return config;
    }

    /**
     * @brief Configuration for near-critical calculations
     */
    static FlashConfig nearCritical() {
        FlashConfig config;
        config.method = FlashMethod::Helmholtz;
        config.handle_near_critical = true;
        config.critical_distance = 0.10;
        config.tolerance = 1e-10;
        config.use_acceleration = false;  // More stable without acceleration
        return config;
    }

    /**
     * @brief Configuration for three-phase (VLLE) calculations
     */
    static FlashConfig threePhase() {
        FlashConfig config;
        config.num_stability_trials = 15;
        config.max_iterations = 200;
        config.tolerance = 1e-9;
        config.max_phases = 3;
        return config;
    }
};

/**
 * @brief Configuration for specific flash types
 */
struct PHFlashConfig : public FlashConfig {
    double enthalpy_tolerance = 100.0;  ///< Enthalpy matching tolerance [J/mol]
    double temp_bracket_low = 100.0;    ///< Temperature lower bound [K]
    double temp_bracket_high = 1000.0;  ///< Temperature upper bound [K]
    bool use_newton = true;             ///< Use Newton for T iteration
};

struct PSFlashConfig : public FlashConfig {
    double entropy_tolerance = 0.1;     ///< Entropy matching tolerance [J/mol/K]
    double temp_bracket_low = 100.0;    ///< Temperature lower bound [K]
    double temp_bracket_high = 1000.0;  ///< Temperature upper bound [K]
    bool use_newton = true;             ///< Use Newton for T iteration
};

struct TVFlashConfig : public FlashConfig {
    bool helmholtz_only = true;         ///< Always use Helmholtz for TV flash

    // Strategy selection for TV solve.
    TVSolveStrategy tv_solve_strategy = TVSolveStrategy::BrentOuter;

    // Outer pressure bracketing for Helmholtz TV solve.
    // The TV solver uses an outer root solve in log(P) to match the target V.
    int tv_bracket_steps = 50;          ///< Max samples on each side of initial guess
    double tv_bracket_logp_step = 0.5;  ///< log10(P) step per sample (0.5 = half-decade)
    double tv_p_min = 1e-3;             ///< Minimum P considered by bracket scan [Pa]
    double tv_p_max = 1e12;             ///< Maximum P considered by bracket scan [Pa]
};

} // namespace Config
} // namespace DMThermo

#endif // THERMO_CONFIG_FLASH_CONFIG_H
