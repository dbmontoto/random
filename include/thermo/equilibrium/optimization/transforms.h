/**
 * @file transforms.h
 * @brief Variable transforms for constrained equilibrium optimization.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_TRANSFORMS_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_TRANSFORMS_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

inline std::vector<double> softmaxWithFixedLast(
    const std::vector<double>& u,
    double min_value)
{
    // Output size = u.size()+1, last logit fixed at 0.0.
    const size_t n = u.size() + 1;
    std::vector<double> logits(n, 0.0);
    for (size_t i = 0; i < u.size(); ++i) logits[i] = u[i];

    const double max_logit = *std::max_element(logits.begin(), logits.end());
    std::vector<double> expv(n, 0.0);
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        expv[i] = std::exp(logits[i] - max_logit);
        sum += expv[i];
    }
    if (!(sum > 0.0 && std::isfinite(sum))) {
        return std::vector<double>(n, 1.0 / static_cast<double>(n));
    }

    std::vector<double> out(n, 0.0);
    for (size_t i = 0; i < n; ++i) out[i] = expv[i] / sum;

    if (min_value > 0.0) {
        for (auto& v : out) v = std::max(v, min_value);
        double s = 0.0;
        for (double v : out) s += v;
        for (auto& v : out) v /= s;
    }

    return out;
}

inline double safeLog(double x, double min_value) {
    return std::log(std::max(x, min_value));
}

inline double sigmoid(double u) {
    if (u >= 0.0) {
        const double e = std::exp(-u);
        return 1.0 / (1.0 + e);
    }
    const double e = std::exp(u);
    return e / (1.0 + e);
}

inline double boundedFromUnbounded(double u, double lo, double hi) {
    const double s = sigmoid(u);
    return lo + (hi - lo) * s;
}

inline double unboundedFromBounded(double x, double lo, double hi) {
    const double eps = 1e-15;
    const double t = (x - lo) / std::max(hi - lo, eps);
    const double tc = std::clamp(t, eps, 1.0 - eps);
    return std::log(tc / (1.0 - tc));
}

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_TRANSFORMS_H
