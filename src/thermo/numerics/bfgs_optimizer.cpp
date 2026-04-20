/**
 * @file bfgs_optimizer.cpp
 * @brief BFGS optimizer implementation (finite-difference gradient by default).
 */

#include "thermo/numerics/ioptimizer.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DMThermo {
namespace Numerics {
namespace Optimization {
namespace {

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("dot: size mismatch");
    }
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        s += a[i] * b[i];
    }
    return s;
}

double norm2(const std::vector<double>& v) {
    return std::sqrt(std::max(0.0, dot(v, v)));
}

std::vector<double> addScaled(const std::vector<double>& x, const std::vector<double>& p, double a) {
    if (x.size() != p.size()) {
        throw std::invalid_argument("addScaled: size mismatch");
    }
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        y[i] = x[i] + a * p[i];
    }
    return y;
}

std::vector<double> numericalGradient(
    const std::function<double(const std::vector<double>&)>& f,
    const std::vector<double>& x,
    double h_base,
    int& evals)
{
    std::vector<double> g(x.size(), 0.0);
    for (size_t i = 0; i < x.size(); ++i) {
        const double hi = h_base * std::max(1.0, std::abs(x[i]));
        std::vector<double> xp = x;
        std::vector<double> xm = x;
        xp[i] += hi;
        xm[i] -= hi;
        const double fp = f(xp); ++evals;
        const double fm = f(xm); ++evals;
        g[i] = (fp - fm) / (2.0 * hi);
    }
    return g;
}

class BFGSOptimizer final : public IOptimizer {
public:
    OptimizationResult minimize(
        std::function<double(const std::vector<double>&)> f,
        const std::vector<double>& x0,
        const Config::OptimizerConfig& config) const override
    {
        return minimizeWithGradient(std::move(f), {}, x0, config);
    }

    OptimizationResult minimizeWithGradient(
        std::function<double(const std::vector<double>&)> f,
        std::function<std::vector<double>(const std::vector<double>&)> grad,
        const std::vector<double>& x0,
        const Config::OptimizerConfig& config) const override
    {
        OptimizationResult out;
        if (x0.empty()) {
            out.converged = true;
            out.x_optimal = x0;
            out.message = "Empty input";
            return out;
        }

        const int n = static_cast<int>(x0.size());
        std::vector<double> x = x0;
        int evals = 0;

        auto evalF = [&](const std::vector<double>& xx) -> double {
            const double v = f(xx);
            return v;
        };

        double fx = evalF(x);
        ++evals;

        std::vector<double> g =
            grad ? grad(x) : numericalGradient(evalF, x, config.jacobian_step, evals);

        // H = inverse Hessian approximation.
        std::vector<std::vector<double>> H(static_cast<size_t>(n), std::vector<double>(static_cast<size_t>(n), 0.0));
        for (int i = 0; i < n; ++i) H[i][i] = 1.0;

        double gnorm = norm2(g);
        if (!std::isfinite(fx) || !std::isfinite(gnorm)) {
            out.converged = false;
            out.x_optimal = x;
            out.f_optimal = fx;
            out.message = "Initial objective/gradient not finite";
            out.function_evaluations = evals;
            return out;
        }

        for (int it = 0; it < config.max_iterations; ++it) {
            if (gnorm < config.tolerance) {
                out.converged = true;
                out.iterations = it;
                break;
            }

            // p = -H g
            std::vector<double> p(static_cast<size_t>(n), 0.0);
            for (int i = 0; i < n; ++i) {
                double s = 0.0;
                for (int j = 0; j < n; ++j) {
                    s += H[i][j] * g[j];
                }
                p[i] = -s;
            }

            // Backtracking line search (Armijo).
            const double c1 = 1e-4;
            double alpha = 1.0;
            const double gtp = dot(g, p);
            if (!(std::isfinite(gtp) && gtp < 0.0)) {
                // Not a descent direction; reset.
                for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) H[i][j] = (i == j) ? 1.0 : 0.0;
                for (int i = 0; i < n; ++i) p[i] = -g[i];
            }

            std::vector<double> x_new = x;
            double f_new = fx;
            for (int ls = 0; ls < 25; ++ls) {
                x_new = addScaled(x, p, alpha);
                f_new = evalF(x_new);
                ++evals;
                if (std::isfinite(f_new) && f_new <= fx + c1 * alpha * dot(g, p)) {
                    break;
                }
                alpha *= 0.5;
                if (alpha < 1e-12) {
                    break;
                }
            }

            if (!std::isfinite(f_new)) {
                out.converged = false;
                out.iterations = it + 1;
                out.message = "Objective became non-finite";
                break;
            }

            std::vector<double> g_new =
                grad ? grad(x_new) : numericalGradient(evalF, x_new, config.jacobian_step, evals);

            const double gnorm_new = norm2(g_new);
            if (!std::isfinite(gnorm_new)) {
                out.converged = false;
                out.iterations = it + 1;
                out.message = "Gradient became non-finite";
                break;
            }

            // Convergence on function change.
            if (std::abs(f_new - fx) < config.function_tolerance && gnorm_new < config.tolerance * 10.0) {
                x = x_new;
                fx = f_new;
                g = g_new;
                gnorm = gnorm_new;
                out.converged = true;
                out.iterations = it + 1;
                break;
            }

            // BFGS update.
            std::vector<double> s(static_cast<size_t>(n), 0.0), y(static_cast<size_t>(n), 0.0);
            for (int i = 0; i < n; ++i) {
                s[i] = x_new[i] - x[i];
                y[i] = g_new[i] - g[i];
            }

            const double ys = dot(y, s);
            if (std::isfinite(ys) && ys > 1e-14) {
                const double rho = 1.0 / ys;

                // Compute (I - rho s y^T)
                std::vector<std::vector<double>> A(static_cast<size_t>(n), std::vector<double>(static_cast<size_t>(n), 0.0));
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        A[i][j] = -rho * s[i] * y[j];
                    }
                    A[i][i] += 1.0;
                }

                // H = A H A^T + rho s s^T
                std::vector<std::vector<double>> AH(static_cast<size_t>(n), std::vector<double>(static_cast<size_t>(n), 0.0));
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        double sum = 0.0;
                        for (int k = 0; k < n; ++k) sum += A[i][k] * H[k][j];
                        AH[i][j] = sum;
                    }
                }

                std::vector<std::vector<double>> Hnew(static_cast<size_t>(n), std::vector<double>(static_cast<size_t>(n), 0.0));
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        double sum = 0.0;
                        for (int k = 0; k < n; ++k) sum += AH[i][k] * A[j][k]; // A^T
                        Hnew[i][j] = sum;
                    }
                }

                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        Hnew[i][j] += rho * s[i] * s[j];
                    }
                }

                H = std::move(Hnew);
            } else {
                // Reset if curvature condition fails.
                for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) H[i][j] = (i == j) ? 1.0 : 0.0;
            }

            x = x_new;
            fx = f_new;
            g = g_new;
            gnorm = gnorm_new;

            out.iterations = it + 1;
        }

        out.x_optimal = x;
        out.f_optimal = fx;
        out.gradient_norm = gnorm;
        out.function_evaluations = evals;
        if (out.message.empty()) {
            out.message = out.converged ? "Converged" : "Did not converge";
        }
        return out;
    }

    std::string name() const override { return "BFGS"; }
};

} // namespace

OptimizerPtr makeBFGSOptimizer() {
    return std::make_shared<BFGSOptimizer>();
}

} // namespace Optimization
} // namespace Numerics
} // namespace DMThermo

