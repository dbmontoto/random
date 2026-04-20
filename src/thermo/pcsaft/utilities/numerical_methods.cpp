#include "thermo/pcsaft/utilities/numerical_methods.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

namespace pcsaft {

// ============================================================================
// 1D Root Finding Methods
// ============================================================================

double NumericalMethods::newtonRaphson(
    std::function<double(double)> f,
    std::function<double(double)> df,
    double x0,
    double tol,
    int max_iter
) {
    double x = x0;

    // Parameters for Armijo backtracking line search
    const double c1 = 1e-4;        // Armijo condition parameter
    const double rho = 0.5;        // Step reduction factor
    const double alpha_min = 1e-10; // Minimum step size
    const int max_backtrack = 15;  // Maximum backtracking iterations

    for (int iter = 0; iter < max_iter; ++iter) {
        double fx = f(x);
        double dfx = df(x);

        if (std::abs(fx) < tol) {
            return x;
        }

        if (std::abs(dfx) < 1e-14) {
            throw std::runtime_error("Newton-Raphson: derivative too small");
        }

        double dx = fx / dfx;  // Newton step (note: we subtract dx, so this is positive when fx/dfx > 0)

        // Armijo backtracking line search
        double alpha = 1.0;
        double f_abs = std::abs(fx);

        for (int bt = 0; bt < max_backtrack; ++bt) {
            double x_trial = x - alpha * dx;
            double f_trial_abs = std::abs(f(x_trial));

            // Armijo condition: sufficient decrease
            if (f_trial_abs < f_abs * (1.0 - c1 * alpha)) {
                x = x_trial;
                break;
            }

            alpha *= rho;

            if (alpha < alpha_min) {
                // Accept reduced step
                x = x - alpha * dx;
                break;
            }

            // Last resort: use full step
            if (bt == max_backtrack - 1) {
                x = x - dx;
            }
        }

        if (std::abs(dx) * alpha < tol) {
            return x;
        }
    }

    throw std::runtime_error("Newton-Raphson: failed to converge");
}

double NumericalMethods::brentMethod(
    std::function<double(double)> f,
    double a,
    double b,
    double tol,
    int max_iter
) {
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0) {
        throw std::invalid_argument("Brent's method: f(a) and f(b) must have opposite signs");
    }

    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    double c = a;
    double fc = fa;
    bool mflag = true;
    double d = 0.0;

    for (int iter = 0; iter < max_iter; ++iter) {
        if (std::abs(fb) < tol || std::abs(b - a) < tol) {
            return b;
        }

        double s;

        if (fa != fc && fb != fc) {
            // Inverse quadratic interpolation
            s = a * fb * fc / ((fa - fb) * (fa - fc))
              + b * fa * fc / ((fb - fa) * (fb - fc))
              + c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            // Secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        // Conditions to use bisection instead
        bool cond1 = (s < (3*a + b)/4 || s > b);
        bool cond2 = mflag && std::abs(s - b) >= std::abs(b - c) / 2;
        bool cond3 = !mflag && std::abs(s - b) >= std::abs(c - d) / 2;
        bool cond4 = mflag && std::abs(b - c) < tol;
        bool cond5 = !mflag && std::abs(c - d) < tol;

        if (cond1 || cond2 || cond3 || cond4 || cond5) {
            s = (a + b) / 2;
            mflag = true;
        } else {
            mflag = false;
        }

        double fs = f(s);
        d = c;
        c = b;
        fc = fb;

        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    throw std::runtime_error("Brent's method: failed to converge");
}

double NumericalMethods::secantMethod(
    std::function<double(double)> f,
    double x0,
    double x1,
    double tol,
    int max_iter
) {
    double f0 = f(x0);
    double f1 = f(x1);

    for (int iter = 0; iter < max_iter; ++iter) {
        if (std::abs(f1) < tol) {
            return x1;
        }

        if (std::abs(f1 - f0) < 1e-14) {
            throw std::runtime_error("Secant method: function values too close");
        }

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        if (std::abs(x2 - x1) < tol) {
            return x2;
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    throw std::runtime_error("Secant method: failed to converge");
}

bool NumericalMethods::findBracket(
    std::function<double(double)> f,
    double x0,
    double& a,
    double& b,
    int max_steps
) {
    const double factor = 1.6;
    double step = 0.1;

    a = x0 - step;
    b = x0 + step;

    double fa = f(a);
    double fb = f(b);

    for (int i = 0; i < max_steps; ++i) {
        if (fa * fb < 0.0) {
            return true;
        }

        // Expand bracket
        if (std::abs(fa) < std::abs(fb)) {
            a -= step * factor;
            step *= factor;
            fa = f(a);
        } else {
            b += step * factor;
            step *= factor;
            fb = f(b);
        }
    }

    return false;
}

// ============================================================================
// Multi-dimensional Newton Solver
// ============================================================================

std::vector<double> NumericalMethods::newtonMultiD(
    std::function<std::vector<double>(const std::vector<double>&)> F,
    std::function<std::vector<std::vector<double>>(const std::vector<double>&)> J,
    const std::vector<double>& x0,
    double tol,
    int max_iter
) {
    std::vector<double> x = x0;
    int n = x0.size();

    // Parameters for Armijo backtracking line search
    const double c1 = 1e-4;        // Armijo condition parameter
    const double rho = 0.5;        // Step reduction factor
    const double alpha_min = 1e-8; // Minimum step size
    const int max_backtrack = 20;  // Maximum backtracking iterations

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> f = F(x);
        std::vector<std::vector<double>> jac = J(x);

        // Check convergence
        double error = norm(f);
        if (error < tol) {
            return x;
        }

        // Solve J * dx = -f
        std::vector<double> neg_f(n);
        for (int i = 0; i < n; ++i) {
            neg_f[i] = -f[i];
        }

        std::vector<double> dx = solveLinearSystem(jac, neg_f);

        // Armijo backtracking line search
        // Find step size alpha such that ||F(x + alpha*dx)|| < ||F(x)|| * (1 - c1*alpha)
        double alpha = 1.0;
        double f_norm = error;
        bool found_step = false;

        for (int bt = 0; bt < max_backtrack; ++bt) {
            // Trial point
            std::vector<double> x_trial(n);
            for (int i = 0; i < n; ++i) {
                x_trial[i] = x[i] + alpha * dx[i];
            }

            // Evaluate function at trial point
            std::vector<double> f_trial = F(x_trial);
            double f_trial_norm = norm(f_trial);

            // Armijo condition: sufficient decrease
            if (f_trial_norm < f_norm * (1.0 - c1 * alpha)) {
                // Accept step
                x = x_trial;
                found_step = true;
                break;
            }

            // Reduce step size
            alpha *= rho;

            if (alpha < alpha_min) {
                // Step too small, accept current step anyway
                x = x_trial;
                found_step = true;
                break;
            }
        }

        // Fallback: if line search fails, use full step
        if (!found_step) {
            for (int i = 0; i < n; ++i) {
                x[i] += dx[i];
            }
        }

        if (norm(dx) * alpha < tol) {
            return x;
        }
    }

    throw std::runtime_error("Multi-dimensional Newton: failed to converge");
}

std::vector<double> NumericalMethods::newtonMultiD_numerical(
    std::function<std::vector<double>(const std::vector<double>&)> F,
    const std::vector<double>& x0,
    double tol,
    int max_iter
) {
    auto J = [&](const std::vector<double>& x) {
        return numericalJacobian(F, x);
    };

    return newtonMultiD(F, J, x0, tol, max_iter);
}

// ============================================================================
// Successive Substitution
// ============================================================================

std::vector<double> NumericalMethods::successiveSubstitution(
    std::function<std::vector<double>(const std::vector<double>&)> g,
    const std::vector<double>& x0,
    double tol,
    int max_iter,
    ConvergenceResult* result
) {
    std::vector<double> x = x0;
    int n = x0.size();

    ConvergenceResult conv_result;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> x_new = g(x);

        // Calculate error
        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            error = std::max(error, std::abs(x_new[i] - x[i]));
        }

        conv_result.iterations = iter + 1;
        conv_result.final_error = error;

        if (error < tol) {
            conv_result.converged = true;
            conv_result.message = "Converged";
            if (result) *result = conv_result;
            return x_new;
        }

        x = x_new;
    }

    conv_result.converged = false;
    conv_result.message = "Maximum iterations reached";
    if (result) *result = conv_result;

    return x;
}

std::vector<double> NumericalMethods::successiveSubstitutionGDEM(
    std::function<std::vector<double>(const std::vector<double>&)> g,
    const std::vector<double>& x0,
    double tol,
    int max_iter,
    ConvergenceResult* result
) {
    // GDEM acceleration: periodically use previous iterations to accelerate
    std::vector<double> x = x0;
    int n = x0.size();

    std::vector<double> x_prev = x;
    std::vector<double> x_prev2 = x;

    ConvergenceResult conv_result;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> x_new = g(x);

        // Calculate error
        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            error = std::max(error, std::abs(x_new[i] - x[i]));
        }

        conv_result.iterations = iter + 1;
        conv_result.final_error = error;

        if (error < tol) {
            conv_result.converged = true;
            conv_result.message = "Converged";
            if (result) *result = conv_result;
            return x_new;
        }

        // GDEM acceleration every 3 iterations
        if (iter > 2 && iter % 3 == 0) {
            // Compute acceleration factor
            std::vector<double> delta1(n), delta2(n);
            for (int i = 0; i < n; ++i) {
                delta1[i] = x[i] - x_prev[i];
                delta2[i] = x_prev[i] - x_prev2[i];
            }

            double num = 0.0, denom = 0.0;
            for (int i = 0; i < n; ++i) {
                num += delta1[i] * delta1[i];
                denom += (delta1[i] - delta2[i]) * delta1[i];
            }

            if (std::abs(denom) > 1e-12) {
                double lambda = num / denom;
                lambda = std::max(0.0, std::min(1.0, lambda)); // Clamp

                // Accelerated update
                for (int i = 0; i < n; ++i) {
                    x_new[i] = x[i] + lambda * (x_new[i] - x[i]);
                }
            }
        }

        x_prev2 = x_prev;
        x_prev = x;
        x = x_new;
    }

    conv_result.converged = false;
    conv_result.message = "Maximum iterations reached";
    if (result) *result = conv_result;

    return x;
}

// ============================================================================
// Utility Functions
// ============================================================================

std::vector<std::vector<double>> NumericalMethods::numericalJacobian(
    std::function<std::vector<double>(const std::vector<double>&)> F,
    const std::vector<double>& x,
    double h
) {
    int n = x.size();
    std::vector<double> f0 = F(x);
    int m = f0.size();

    std::vector<std::vector<double>> J(m, std::vector<double>(n));

    for (int j = 0; j < n; ++j) {
        std::vector<double> x_plus = x;
        double hj = h * std::max(1.0, std::abs(x[j]));
        x_plus[j] += hj;

        std::vector<double> f_plus = F(x_plus);

        for (int i = 0; i < m; ++i) {
            J[i][j] = (f_plus[i] - f0[i]) / hj;
        }
    }

    return J;
}

std::vector<double> NumericalMethods::numericalGradient(
    std::function<double(const std::vector<double>&)> f,
    const std::vector<double>& x,
    double h
) {
    int n = x.size();
    std::vector<double> grad(n);

    for (int i = 0; i < n; ++i) {
        std::vector<double> x_plus = x;
        std::vector<double> x_minus = x;

        double hi = h * std::max(1.0, std::abs(x[i]));
        x_plus[i] += hi;
        x_minus[i] -= hi;

        grad[i] = (f(x_plus) - f(x_minus)) / (2.0 * hi);
    }

    return grad;
}

std::vector<double> NumericalMethods::solveLinearSystem(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& b
) {
    int n = b.size();

    // Create augmented matrix
    std::vector<std::vector<double>> aug(n, std::vector<double>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            aug[i][j] = A[i][j];
        }
        aug[i][n] = b[i];
    }

    // Gaussian elimination with partial pivoting
    for (int k = 0; k < n; ++k) {
        // Find pivot
        int pivot_row = k;
        double max_val = std::abs(aug[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(aug[i][k]) > max_val) {
                max_val = std::abs(aug[i][k]);
                pivot_row = i;
            }
        }

        if (max_val < 1e-14) {
            throw std::runtime_error("Linear system is singular");
        }

        // Swap rows
        if (pivot_row != k) {
            std::swap(aug[k], aug[pivot_row]);
        }

        // Eliminate
        for (int i = k + 1; i < n; ++i) {
            double factor = aug[i][k] / aug[k][k];
            for (int j = k; j <= n; ++j) {
                aug[i][j] -= factor * aug[k][j];
            }
        }
    }

    // Back substitution
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = aug[i][n];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }

    return x;
}

double NumericalMethods::norm(const std::vector<double>& v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

double NumericalMethods::normInf(const std::vector<double>& v) {
    double max_val = 0.0;
    for (double val : v) {
        max_val = std::max(max_val, std::abs(val));
    }
    return max_val;
}

double NumericalMethods::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    int n = std::min(a.size(), b.size());
    for (int i = 0; i < n; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

} // namespace pcsaft
