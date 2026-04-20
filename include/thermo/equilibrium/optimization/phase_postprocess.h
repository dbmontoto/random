/**
 * @file phase_postprocess.h
 * @brief Shared utilities for pruning, normalizing, merging, and ordering phases.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_PHASE_POSTPROCESS_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_PHASE_POSTPROCESS_H

#include "thermo/core/constants.h"
#include "thermo/core/types.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

inline int phaseTypeRank(PhaseType p) {
    if (p == PhaseType::Vapor) return 0;
    if (p == PhaseType::Liquid || p == PhaseType::Liquid1) return 1;
    if (p == PhaseType::Liquid2) return 2;
    return 3;
}

inline bool areDuplicatePhases(
    const PhaseState& a,
    const PhaseState& b,
    double density_rel_tol = 1e-8,
    double x_abs_tol = 1e-8)
{
    if (!(a.density > 0.0 && b.density > 0.0)) return false;
    if (!(std::isfinite(a.density) && std::isfinite(b.density))) return false;
    const double rel = std::abs(a.density - b.density) / std::max(a.density, b.density);
    if (rel > density_rel_tol) return false;

    if (a.x.size() != b.x.size()) return false;
    double maxdx = 0.0;
    for (size_t i = 0; i < a.x.size(); ++i) {
        maxdx = std::max(maxdx, std::abs(a.x[i] - b.x[i]));
    }
    return maxdx < x_abs_tol;
}

struct PhasePostProcessOptions {
    double prune_fraction = 0.0;
    bool allow_empty = false;
    bool merge_duplicates = false;
    double duplicate_density_rel_tol = 1e-8;
    double duplicate_x_abs_tol = 1e-8;
    bool sort_by_type = false;
    bool canonicalize_liquid = true;
};

inline bool postProcessPhases(std::vector<PhaseState>& phases, const PhasePostProcessOptions& opts) {
    if (opts.prune_fraction > 0.0) {
        phases.erase(std::remove_if(phases.begin(), phases.end(), [&](const PhaseState& ps) {
            return ps.fraction < opts.prune_fraction;
        }), phases.end());
    }

    if (opts.merge_duplicates && phases.size() >= 2) {
        for (size_t i = 0; i + 1 < phases.size(); ++i) {
            for (size_t j = i + 1; j < phases.size(); ++j) {
                if (areDuplicatePhases(phases[i], phases[j], opts.duplicate_density_rel_tol, opts.duplicate_x_abs_tol)) {
                    phases[i].fraction += phases[j].fraction;
                    phases.erase(phases.begin() + static_cast<std::ptrdiff_t>(j));
                    --j;
                }
            }
        }
    }

    if (phases.empty()) {
        return opts.allow_empty;
    }

    double sum = 0.0;
    for (const auto& ph : phases) sum += ph.fraction;
    const double denom = std::max(sum, Constants::MIN_MOLE_FRACTION);
    for (auto& ph : phases) ph.fraction /= denom;

    if (opts.sort_by_type) {
        std::stable_sort(phases.begin(), phases.end(), [](const PhaseState& a, const PhaseState& b) {
            return phaseTypeRank(a.type) < phaseTypeRank(b.type);
        });
    }

    if (opts.canonicalize_liquid) {
        for (auto& ph : phases) {
            if (ph.type == PhaseType::Liquid) ph.type = PhaseType::Liquid1;
        }
    }

    return true;
}

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_PHASE_POSTPROCESS_H

