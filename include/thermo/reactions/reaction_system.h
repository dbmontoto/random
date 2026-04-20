/**
 * @file reaction_system.h
 * @brief Basic reaction system types (stoichiometry over mixture components)
 */

#ifndef THERMO_REACTIONS_REACTION_SYSTEM_H
#define THERMO_REACTIONS_REACTION_SYSTEM_H

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <limits>

namespace DMThermo {
namespace Reactions {

struct Reaction {
    std::string name;
    // Stoichiometric coefficients for the mixture components.
    // Convention: nu_i < 0 reactants, nu_i > 0 products.
    std::vector<double> nu;
};

struct ExtentBounds {
    double min = 0.0;
    double max = 0.0;
};

class ReactionSystem {
public:
    ReactionSystem(int num_components, std::vector<Reaction> reactions)
        : num_components_(num_components), reactions_(std::move(reactions))
    {
        if (num_components_ <= 0) {
            throw std::invalid_argument("ReactionSystem: num_components must be > 0");
        }
        for (std::size_t r = 0; r < reactions_.size(); ++r) {
            if (static_cast<int>(reactions_[r].nu.size()) != num_components_) {
                throw std::invalid_argument("ReactionSystem: reaction stoichiometry size mismatch at r=" + std::to_string(r));
            }
        }
    }

    int numComponents() const { return num_components_; }
    int numReactions() const { return static_cast<int>(reactions_.size()); }

    const Reaction& reaction(int r) const {
        if (r < 0 || r >= numReactions()) {
            throw std::out_of_range("ReactionSystem: reaction index out of range");
        }
        return reactions_[static_cast<std::size_t>(r)];
    }

    /**
     * @brief Conservative box bounds for a single reaction extent given n0 and non-negativity constraints.
     *
     * Bounds enforce `n_i = n0_i + nu_i * xi >= min_moles` for all components i.
     *
     * Notes:
     * - For multiple simultaneous reactions, this is a conservative per-extent box that does not guarantee global feasibility.
     * - Throws if the reaction has no reactants or no products, or if the bounds are infeasible.
     */
    ExtentBounds extentBounds(int r, const std::vector<double>& n0, double min_moles) const {
        if (static_cast<int>(n0.size()) != num_components_) {
            throw std::invalid_argument("ReactionSystem::extentBounds: n0 size mismatch");
        }
        if (!(std::isfinite(min_moles) && min_moles >= 0.0)) {
            throw std::invalid_argument("ReactionSystem::extentBounds: min_moles must be finite and >= 0");
        }
        const auto& nu = reaction(r).nu;

        bool has_pos = false;
        bool has_neg = false;
        double lo = -std::numeric_limits<double>::infinity();
        double hi = std::numeric_limits<double>::infinity();

        for (int i = 0; i < num_components_; ++i) {
            const double nui = nu[static_cast<std::size_t>(i)];
            if (!std::isfinite(nui)) {
                throw std::invalid_argument("ReactionSystem::extentBounds: stoichiometry contains non-finite value");
            }
            if (nui > 0.0) {
                has_pos = true;
                const double bound = (min_moles - n0[static_cast<std::size_t>(i)]) / nui;
                lo = std::max(lo, bound);
            } else if (nui < 0.0) {
                has_neg = true;
                const double bound = (n0[static_cast<std::size_t>(i)] - min_moles) / (-nui);
                hi = std::min(hi, bound);
            }
        }

        if (!has_pos || !has_neg) {
            throw std::invalid_argument("ReactionSystem::extentBounds: each reaction must have both reactants and products");
        }
        if (!(std::isfinite(lo) && std::isfinite(hi)) || lo > hi) {
            throw std::invalid_argument("ReactionSystem::extentBounds: infeasible extent bounds for given n0/min_moles");
        }

        return ExtentBounds{lo, hi};
    }

    std::vector<double> molesFromExtents(const std::vector<double>& n0, const std::vector<double>& extents) const {
        if (static_cast<int>(n0.size()) != num_components_) {
            throw std::invalid_argument("ReactionSystem: n0 size mismatch");
        }
        if (static_cast<int>(extents.size()) != numReactions()) {
            throw std::invalid_argument("ReactionSystem: extents size mismatch");
        }

        std::vector<double> n = n0;
        for (int r = 0; r < numReactions(); ++r) {
            const double xi = extents[static_cast<std::size_t>(r)];
            const auto& nu = reactions_[static_cast<std::size_t>(r)].nu;
            for (int i = 0; i < num_components_; ++i) {
                n[static_cast<std::size_t>(i)] += nu[static_cast<std::size_t>(i)] * xi;
            }
        }
        return n;
    }

private:
    int num_components_ = 0;
    std::vector<Reaction> reactions_;
};

} // namespace Reactions
} // namespace DMThermo

#endif // THERMO_REACTIONS_REACTION_SYSTEM_H
