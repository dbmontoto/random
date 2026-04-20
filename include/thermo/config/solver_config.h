/**
 * @file solver_config.h
 * @brief Configuration structures for numerical solvers
 */

#ifndef THERMO_CONFIG_SOLVER_CONFIG_H
#define THERMO_CONFIG_SOLVER_CONFIG_H

namespace DMThermo {
namespace Config {

/**
 * @brief Configuration for scalar root finding
 */
struct RootFinderConfig {
    double tolerance = 1e-10;           ///< Convergence tolerance
    double rel_tolerance = 1e-10;       ///< Relative tolerance
    int max_iterations = 50;            ///< Maximum iterations
    double damping_factor = 1.0;        ///< Newton damping factor
    bool use_line_search = true;        ///< Use backtracking line search
    double line_search_alpha = 1e-4;    ///< Armijo condition parameter
    double line_search_beta = 0.5;      ///< Step reduction factor

    static RootFinderConfig defaults() {
        return RootFinderConfig{};
    }

    static RootFinderConfig highAccuracy() {
        RootFinderConfig config;
        config.tolerance = 1e-14;
        config.rel_tolerance = 1e-12;
        config.max_iterations = 100;
        return config;
    }

    static RootFinderConfig robust() {
        RootFinderConfig config;
        config.damping_factor = 0.5;
        config.use_line_search = true;
        config.max_iterations = 100;
        return config;
    }
};

/**
 * @brief Configuration for multidimensional optimization
 */
struct OptimizerConfig {
    double tolerance = 1e-8;            ///< Gradient norm tolerance
    double function_tolerance = 1e-10;  ///< Function value tolerance
    int max_iterations = 100;           ///< Maximum iterations
    double initial_trust_radius = 1.0;  ///< Initial trust region radius
    double max_trust_radius = 100.0;    ///< Maximum trust region radius
    bool use_numerical_jacobian = true; ///< Use numerical Jacobian
    double jacobian_step = 1e-7;        ///< Finite difference step

    static OptimizerConfig defaults() {
        return OptimizerConfig{};
    }

    static OptimizerConfig highAccuracy() {
        OptimizerConfig config;
        config.tolerance = 1e-12;
        config.function_tolerance = 1e-14;
        config.jacobian_step = 1e-9;
        return config;
    }
};

/**
 * @brief Configuration for fixed-point iteration
 */
struct FixedPointConfig {
    double tolerance = 1e-8;            ///< Convergence tolerance
    int max_iterations = 100;           ///< Maximum iterations
    bool use_acceleration = true;       ///< Use GDEM/Wegstein
    double acceleration_delay = 3;      ///< Iterations before accelerating
    double max_acceleration = 5.0;      ///< Maximum acceleration factor
    double damping_factor = 1.0;        ///< Direct substitution damping

    static FixedPointConfig defaults() {
        return FixedPointConfig{};
    }

    static FixedPointConfig stable() {
        FixedPointConfig config;
        config.use_acceleration = false;
        config.damping_factor = 0.5;
        config.max_iterations = 200;
        return config;
    }
};

/**
 * @brief Configuration for successive substitution in association
 */
struct AssociationSolverConfig {
    double tolerance = 1e-10;           ///< X_A convergence tolerance
    int max_iterations = 100;           ///< Maximum iterations
    double damping_factor = 1.0;        ///< Damping for difficult cases
    bool use_analytical_jacobian = false; ///< Use Newton if available

    static AssociationSolverConfig defaults() {
        return AssociationSolverConfig{};
    }
};

} // namespace Config
} // namespace DMThermo

#endif // THERMO_CONFIG_SOLVER_CONFIG_H
