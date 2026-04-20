/**
 * @file variable_blocks.h
 * @brief Helpers for encoding/decoding constrained variable blocks.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_VARIABLE_BLOCKS_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_VARIABLE_BLOCKS_H

#include "thermo/equilibrium/optimization/transforms.h"
#include "thermo/equilibrium/optimization/variable_layout.h"
#include "thermo/reactions/reaction_system.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

inline std::vector<double> decodeBoundedVector(const VarSlice& u, double lo, double hi) {
    if (u.size < 0) {
        throw std::invalid_argument("decodeBoundedVector: negative size");
    }
    std::vector<double> x(static_cast<size_t>(u.size), 0.0);
    for (int i = 0; i < u.size; ++i) {
        x[static_cast<size_t>(i)] = boundedFromUnbounded(u[i], lo, hi);
    }
    return x;
}

inline std::vector<double> encodeBoundedVector(const std::vector<double>& x, double lo, double hi) {
    std::vector<double> u(x.size(), 0.0);
    for (size_t i = 0; i < x.size(); ++i) {
        u[i] = unboundedFromBounded(x[i], lo, hi);
    }
    return u;
}

inline std::vector<double> decodeExtents(
    const VarSlice& u_xi,
    const std::vector<Reactions::ExtentBounds>& bounds)
{
    if (u_xi.size != static_cast<int>(bounds.size())) {
        throw std::invalid_argument("decodeExtents: size mismatch");
    }
    std::vector<double> xi(static_cast<size_t>(u_xi.size), 0.0);
    for (int r = 0; r < u_xi.size; ++r) {
        const auto& b = bounds[static_cast<size_t>(r)];
        xi[static_cast<size_t>(r)] = boundedFromUnbounded(u_xi[r], b.min, b.max);
    }
    return xi;
}

inline std::vector<double> encodeExtents(
    const std::vector<double>& xi,
    const std::vector<Reactions::ExtentBounds>& bounds)
{
    if (xi.size() != bounds.size()) {
        throw std::invalid_argument("encodeExtents: size mismatch");
    }
    std::vector<double> u(xi.size(), 0.0);
    for (size_t r = 0; r < xi.size(); ++r) {
        const auto& b = bounds[r];
        u[r] = unboundedFromBounded(xi[r], b.min, b.max);
    }
    return u;
}

inline std::vector<double> clampToBounds(
    const std::vector<double>& x,
    double lo,
    double hi)
{
    std::vector<double> out = x;
    for (auto& v : out) {
        if (!std::isfinite(v)) v = lo;
        v = std::clamp(v, lo, hi);
    }
    return out;
}

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_VARIABLE_BLOCKS_H

