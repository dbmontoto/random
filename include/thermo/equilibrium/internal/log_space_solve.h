/**
 * @file log_space_solve.h
 * @brief Shared log-space bracketing + Brent root solve utilities (internal).
 */

#ifndef THERMO_EQUILIBRIUM_INTERNAL_LOG_SPACE_SOLVE_H
#define THERMO_EQUILIBRIUM_INTERNAL_LOG_SPACE_SOLVE_H

#include "thermo/numerics/iroot_finder.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Internal {

struct LogSpaceRootSolveResult {
    bool converged = false;
    double root = 0.0; // root in log-space
    int iterations = 0;
    std::string message;
};

template <typename F>
inline LogSpaceRootSolveResult solveLogSpaceRootBrent(
    const Numerics::RootFinding::IScalarRootFinder& root_finder,
    F f_logx,
    double logx0,
    double step,
    int max_steps,
    double x_min,
    double x_max,
    const Config::RootFinderConfig& config,
    const char* label)
{
    LogSpaceRootSolveResult out;

    std::vector<std::pair<double, double>> samples;
    samples.reserve(200);
    for (int k = 0; k <= max_steps; ++k) {
        const int signs[3] = {0, -1, +1};
        for (int sgn : signs) {
            if (k == 0 && sgn != 0) continue;
            if (k > 0 && sgn == 0) continue;

            const double logx = logx0 + static_cast<double>(sgn * k) * step;
            const double x = std::exp(logx);
            if (x < x_min || x > x_max) continue;

            const double fv = f_logx(logx);
            if (std::isfinite(fv)) samples.emplace_back(logx, fv);
        }
    }

    if (samples.size() < 2) {
        out.converged = false;
        out.message = std::string(label) + ": failed to bracket (no finite evaluations)";
        return out;
    }
    std::sort(samples.begin(), samples.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    bool have_bracket = false;
    double a = 0.0;
    double b = 0.0;
    for (size_t i = 1; i < samples.size(); ++i) {
        const auto [x0, f0] = samples[i - 1];
        const auto [x1, f1] = samples[i];
        if (f0 == 0.0) { have_bracket = true; a = b = x0; break; }
        if (f1 == 0.0) { have_bracket = true; a = b = x1; break; }
        if ((f0 < 0.0 && f1 > 0.0) || (f0 > 0.0 && f1 < 0.0)) { have_bracket = true; a = x0; b = x1; break; }
    }

    if (!have_bracket) {
        out.converged = false;
        out.message = std::string(label) + ": failed to bracket (no sign change found)";
        return out;
    }

    out.converged = true;
    out.root = a;
    out.iterations = 0;

    if (a != b) {
        const auto rr = root_finder.solveBracketed(f_logx, a, b, config);
        if (!rr.converged) {
            out.converged = false;
            out.message = std::string(label) + ": Brent solve failed: " + rr.message;
            return out;
        }
        out.root = rr.root;
        out.iterations = rr.iterations;
    }

    return out;
}

} // namespace Internal
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_INTERNAL_LOG_SPACE_SOLVE_H

