/**
 * @file fixed_point.cpp
 * @brief Implementations of fixed-point iteration solvers.
 */

#include "thermo/numerics/ifixed_point.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace DMThermo {
namespace Numerics {
namespace FixedPoint {
namespace {

double normInfDiff(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("normInfDiff: size mismatch");
    }
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        m = std::max(m, std::abs(a[i] - b[i]));
    }
    return m;
}

} // namespace

FixedPointResult DirectSubstitution::solve(
    std::function<std::vector<double>(const std::vector<double>&)> g,
    const std::vector<double>& x0,
    const Config::FixedPointConfig& config) const
{
    FixedPointResult out;
    if (x0.empty()) {
        out.converged = true;
        out.x = x0;
        out.message = "Empty input";
        return out;
    }

    std::vector<double> x = x0;
    std::vector<double> x_new = x0;
    out.residual_history.clear();

    for (int it = 0; it < config.max_iterations; ++it) {
        x_new = g(x);
        if (x_new.size() != x.size()) {
            throw std::invalid_argument("DirectSubstitution: g(x) size mismatch");
        }

        if (config.damping_factor != 1.0) {
            const double a = std::clamp(config.damping_factor, 0.0, 1.0);
            for (size_t i = 0; i < x.size(); ++i) {
                x_new[i] = (1.0 - a) * x[i] + a * x_new[i];
            }
        }

        const double resid = normInfDiff(x_new, x);
        out.iterations = it + 1;
        out.residual = resid;
        if (config.use_acceleration) {
            out.residual_history.push_back(resid);
        }

        x = x_new;
        if (resid < config.tolerance) {
            out.converged = true;
            out.x = x;
            out.message = "Converged";
            return out;
        }
    }

    out.converged = false;
    out.x = x;
    out.message = "Maximum iterations reached";
    return out;
}

FixedPointResult GDEMAccelerated::solve(
    std::function<std::vector<double>(const std::vector<double>&)> g,
    const std::vector<double>& x0,
    const Config::FixedPointConfig& config) const
{
    FixedPointResult out;
    if (x0.empty()) {
        out.converged = true;
        out.x = x0;
        out.message = "Empty input";
        return out;
    }

    std::vector<double> x = x0;
    std::vector<double> x_prev = x0;
    std::vector<double> x_prev2 = x0;
    std::vector<double> x_new = x0;

    for (int it = 0; it < config.max_iterations; ++it) {
        x_new = g(x);
        if (x_new.size() != x.size()) {
            throw std::invalid_argument("GDEM: g(x) size mismatch");
        }

        if (config.damping_factor != 1.0) {
            const double a = std::clamp(config.damping_factor, 0.0, 1.0);
            for (size_t i = 0; i < x.size(); ++i) {
                x_new[i] = (1.0 - a) * x[i] + a * x_new[i];
            }
        }

        const double resid = normInfDiff(x_new, x);
        out.iterations = it + 1;
        out.residual = resid;
        out.residual_history.push_back(resid);

        if (resid < config.tolerance) {
            out.converged = true;
            out.x = x_new;
            out.message = "Converged";
            return out;
        }

        const bool do_accel = config.use_acceleration &&
            it >= static_cast<int>(config.acceleration_delay) &&
            (it % 3 == 0);

        if (do_accel) {
            std::vector<double> delta1(x.size(), 0.0), delta2(x.size(), 0.0);
            for (size_t i = 0; i < x.size(); ++i) {
                delta1[i] = x[i] - x_prev[i];
                delta2[i] = x_prev[i] - x_prev2[i];
            }

            double num = 0.0;
            double denom = 0.0;
            for (size_t i = 0; i < x.size(); ++i) {
                num += delta1[i] * delta1[i];
                denom += (delta1[i] - delta2[i]) * delta1[i];
            }

            if (std::isfinite(num) && std::isfinite(denom) && std::abs(denom) > 1e-14) {
                double lambda = num / denom;
                lambda = std::clamp(lambda, 0.0, 1.0);
                for (size_t i = 0; i < x.size(); ++i) {
                    x_new[i] = x[i] + lambda * (x_new[i] - x[i]);
                }
            }
        }

        x_prev2 = x_prev;
        x_prev = x;
        x = x_new;
    }

    out.converged = false;
    out.x = x;
    out.message = "Maximum iterations reached";
    return out;
}

FixedPointResult WegsteinAccelerated::solve(
    std::function<std::vector<double>(const std::vector<double>&)> g,
    const std::vector<double>& x0,
    const Config::FixedPointConfig& config) const
{
    FixedPointResult out;
    if (x0.empty()) {
        out.converged = true;
        out.x = x0;
        out.message = "Empty input";
        return out;
    }

    std::vector<double> x = x0;
    std::vector<double> x_prev = x0;
    std::vector<double> gx = g(x);
    if (gx.size() != x.size()) {
        throw std::invalid_argument("Wegstein: g(x) size mismatch");
    }

    for (int it = 0; it < config.max_iterations; ++it) {
        std::vector<double> gx_new = g(x);
        if (gx_new.size() != x.size()) {
            throw std::invalid_argument("Wegstein: g(x) size mismatch");
        }

        std::vector<double> x_new = gx_new;
        const bool do_accel = config.use_acceleration &&
            it >= static_cast<int>(config.acceleration_delay);

        if (do_accel) {
            for (size_t i = 0; i < x.size(); ++i) {
                const double denom = (gx_new[i] - gx[i]) - (x[i] - x_prev[i]);
                if (std::abs(denom) > 1e-14) {
                    double q = (gx_new[i] - gx[i]) / denom;
                    q = std::clamp(q, -config.max_acceleration, config.max_acceleration);
                    x_new[i] = q * x[i] + (1.0 - q) * gx_new[i];
                }
            }
        }

        if (config.damping_factor != 1.0) {
            const double a = std::clamp(config.damping_factor, 0.0, 1.0);
            for (size_t i = 0; i < x.size(); ++i) {
                x_new[i] = (1.0 - a) * x[i] + a * x_new[i];
            }
        }

        const double resid = normInfDiff(x_new, x);
        out.iterations = it + 1;
        out.residual = resid;
        out.residual_history.push_back(resid);

        if (resid < config.tolerance) {
            out.converged = true;
            out.x = x_new;
            out.message = "Converged";
            return out;
        }

        x_prev = x;
        gx = gx_new;
        x = x_new;
    }

    out.converged = false;
    out.x = x;
    out.message = "Maximum iterations reached";
    return out;
}

} // namespace FixedPoint
} // namespace Numerics
} // namespace DMThermo

