#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

#include <vector>
#include <functional>
#include <stdexcept>
#include <string>

namespace pcsaft {

/**
 * @brief Convergence result for iterative methods
 */
struct ConvergenceResult {
    bool converged;
    int iterations;
    double final_error;
    std::string message;

    ConvergenceResult() : converged(false), iterations(0), final_error(0.0) {}
};

/**
 * @brief Numerical methods class for root finding and optimization
 */
class NumericalMethods {
public:
    // ========================================================================
    // 1D Root Finding
    // ========================================================================

    /**
     * @brief Newton-Raphson method for finding roots
     * @param f Function to find root of
     * @param df Derivative of f
     * @param x0 Initial guess
     * @param tol Convergence tolerance
     * @param max_iter Maximum iterations
     * @return Root value
     */
    static double newtonRaphson(
        std::function<double(double)> f,
        std::function<double(double)> df,
        double x0,
        double tol = 1e-8,
        int max_iter = 100
    );

    /**
     * @brief Brent's method for finding roots (bracketing method)
     * Robust method that doesn't require derivatives
     * @param f Function to find root of
     * @param a Lower bound (f(a) and f(b) must have opposite signs)
     * @param b Upper bound
     * @param tol Convergence tolerance
     * @param max_iter Maximum iterations
     * @return Root value
     */
    static double brentMethod(
        std::function<double(double)> f,
        double a,
        double b,
        double tol = 1e-8,
        int max_iter = 100
    );

    /**
     * @brief Secant method for finding roots
     * @param f Function to find root of
     * @param x0 First initial guess
     * @param x1 Second initial guess
     * @param tol Convergence tolerance
     * @param max_iter Maximum iterations
     * @return Root value
     */
    static double secantMethod(
        std::function<double(double)> f,
        double x0,
        double x1,
        double tol = 1e-8,
        int max_iter = 100
    );

    /**
     * @brief Find bracket [a, b] where f(a) * f(b) < 0
     * @param f Function
     * @param x0 Starting point
     * @param a Output: lower bound
     * @param b Output: upper bound
     * @param max_steps Maximum expansion steps
     * @return True if bracket found
     */
    static bool findBracket(
        std::function<double(double)> f,
        double x0,
        double& a,
        double& b,
        int max_steps = 50
    );

    // ========================================================================
    // Multi-dimensional Newton Solver
    // ========================================================================

    /**
     * @brief Multi-dimensional Newton's method
     * @param F Function returning residual vector
     * @param J Jacobian matrix function (returns flattened matrix)
     * @param x0 Initial guess vector
     * @param tol Convergence tolerance
     * @param max_iter Maximum iterations
     * @return Solution vector
     */
    static std::vector<double> newtonMultiD(
        std::function<std::vector<double>(const std::vector<double>&)> F,
        std::function<std::vector<std::vector<double>>(const std::vector<double>&)> J,
        const std::vector<double>& x0,
        double tol = 1e-8,
        int max_iter = 100
    );

    /**
     * @brief Multi-dimensional Newton with numerical Jacobian
     * @param F Function returning residual vector
     * @param x0 Initial guess vector
     * @param tol Convergence tolerance
     * @param max_iter Maximum iterations
     * @return Solution vector
     */
    static std::vector<double> newtonMultiD_numerical(
        std::function<std::vector<double>(const std::vector<double>&)> F,
        const std::vector<double>& x0,
        double tol = 1e-8,
        int max_iter = 100
    );

    // ========================================================================
    // Successive Substitution
    // ========================================================================

    /**
     * @brief Simple successive substitution
     * Solves x = g(x) by iteration
     * @param g Fixed-point function
     * @param x0 Initial guess vector
     * @param tol Convergence tolerance
     * @param max_iter Maximum iterations
     * @param result Output: convergence information
     * @return Solution vector
     */
    static std::vector<double> successiveSubstitution(
        std::function<std::vector<double>(const std::vector<double>&)> g,
        const std::vector<double>& x0,
        double tol = 1e-8,
        int max_iter = 1000,
        ConvergenceResult* result = nullptr
    );

    /**
     * @brief Successive substitution with GDEM acceleration
     * GDEM = Generalized Dominant Eigenvalue Method
     * @param g Fixed-point function
     * @param x0 Initial guess vector
     * @param tol Convergence tolerance
     * @param max_iter Maximum iterations
     * @param result Output: convergence information
     * @return Solution vector
     */
    static std::vector<double> successiveSubstitutionGDEM(
        std::function<std::vector<double>(const std::vector<double>&)> g,
        const std::vector<double>& x0,
        double tol = 1e-8,
        int max_iter = 1000,
        ConvergenceResult* result = nullptr
    );

    // ========================================================================
    // Utility Functions
    // ========================================================================

    /**
     * @brief Compute numerical Jacobian using finite differences
     * @param F Function
     * @param x Point at which to evaluate Jacobian
     * @param h Step size (default: adaptive based on x)
     * @return Jacobian matrix
     */
    static std::vector<std::vector<double>> numericalJacobian(
        std::function<std::vector<double>(const std::vector<double>&)> F,
        const std::vector<double>& x,
        double h = 1e-7
    );

    /**
     * @brief Compute numerical gradient
     * @param f Function
     * @param x Point at which to evaluate gradient
     * @param h Step size
     * @return Gradient vector
     */
    static std::vector<double> numericalGradient(
        std::function<double(const std::vector<double>&)> f,
        const std::vector<double>& x,
        double h = 1e-7
    );

    /**
     * @brief Solve linear system Ax = b using Gaussian elimination
     * @param A Coefficient matrix (n x n)
     * @param b Right-hand side vector (n)
     * @return Solution vector x
     */
    static std::vector<double> solveLinearSystem(
        const std::vector<std::vector<double>>& A,
        const std::vector<double>& b
    );

    /**
     * @brief Compute L2 norm of vector
     */
    static double norm(const std::vector<double>& v);

    /**
     * @brief Compute max norm (infinity norm) of vector
     */
    static double normInf(const std::vector<double>& v);

    /**
     * @brief Vector dot product
     */
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
};

} // namespace pcsaft

#endif // NUMERICAL_METHODS_H
