/**
 * @file phase_allocations.h
 * @brief Helpers for multi-phase allocation variables and composition reconstruction.
 *
 * These utilities are shared by equilibrium optimizers that represent phase splits
 * via per-component allocations (softmax logits) and reconstruct phase fractions
 * and phase compositions from them.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_PHASE_ALLOCATIONS_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_PHASE_ALLOCATIONS_H

#include "thermo/core/constants.h"
#include "thermo/core/types.h"
#include "thermo/equilibrium/optimization/transforms.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

inline std::vector<double> clampAndNormalize(const std::vector<double>& x, double x_min) {
    std::vector<double> tmp = x;
    for (auto& v : tmp) {
        if (!std::isfinite(v) || v < x_min) v = x_min;
    }
    double s = std::accumulate(tmp.begin(), tmp.end(), 0.0);
    if (!(std::isfinite(s) && s > 0.0)) {
        throw std::invalid_argument("clampAndNormalize: non-positive sum");
    }
    for (auto& v : tmp) v /= s;
    return tmp;
}

inline PhaseType hintForPhaseIndex(int k) {
    if (k == 0) return PhaseType::Vapor;
    if (k == 1) return PhaseType::Liquid1;
    return PhaseType::Liquid2;
}

struct AllocState {
    int nc = 0;
    int M = 0;
    // w_{i,k} laid out as i-major: w[i*M + k]
    std::vector<double> w;
};

inline AllocState unpackAllocations(const double* u, int u_size, int nc, int M, double min_value) {
    if (M < 1) {
        throw std::invalid_argument("unpackAllocations: M must be >= 1");
    }
    AllocState st;
    st.nc = nc;
    st.M = M;
    st.w.assign(static_cast<size_t>(nc * M), 0.0);
    if (M == 1) {
        for (int i = 0; i < nc; ++i) st.w[static_cast<size_t>(i * M)] = 1.0;
        return st;
    }

    const int per_i = M - 1;
    if (u_size != nc * per_i) {
        throw std::invalid_argument("unpackAllocations: variable size mismatch");
    }

    std::vector<double> logits(static_cast<size_t>(M), 0.0);
    std::vector<double> expv(static_cast<size_t>(M), 0.0);
    for (int i = 0; i < nc; ++i) {
        for (int k = 0; k < per_i; ++k) {
            logits[static_cast<size_t>(k)] = u[i * per_i + k];
        }
        logits[static_cast<size_t>(M - 1)] = 0.0;

        const double max_logit = *std::max_element(logits.begin(), logits.end());
        double sum = 0.0;
        for (int k = 0; k < M; ++k) {
            expv[static_cast<size_t>(k)] = std::exp(logits[static_cast<size_t>(k)] - max_logit);
            sum += expv[static_cast<size_t>(k)];
        }

        if (!(sum > 0.0 && std::isfinite(sum))) {
            const double w = 1.0 / static_cast<double>(M);
            for (int k = 0; k < M; ++k) st.w[static_cast<size_t>(i * M + k)] = w;
            continue;
        }

        double renorm = 0.0;
        for (int k = 0; k < M; ++k) {
            double wk = expv[static_cast<size_t>(k)] / sum;
            if (min_value > 0.0) wk = std::max(wk, min_value);
            st.w[static_cast<size_t>(i * M + k)] = wk;
            renorm += wk;
        }

        if (min_value > 0.0 && renorm > 0.0) {
            for (int k = 0; k < M; ++k) st.w[static_cast<size_t>(i * M + k)] /= renorm;
        }
    }
    return st;
}

inline AllocState unpackAllocations(const std::vector<double>& u, int nc, int M, double min_value) {
    const double* ptr = u.empty() ? nullptr : u.data();
    return unpackAllocations(ptr, static_cast<int>(u.size()), nc, M, min_value);
}

inline std::vector<double> packLogitsFromAllocations(const AllocState& st, double min_value) {
    if (st.M <= 1) return {};
    const int per_i = st.M - 1;
    std::vector<double> u(static_cast<size_t>(st.nc * per_i), 0.0);
    for (int i = 0; i < st.nc; ++i) {
        const double w_last = std::max(st.w[static_cast<size_t>(i * st.M + (st.M - 1))], min_value);
        for (int k = 0; k < st.M - 1; ++k) {
            const double wk = std::max(st.w[static_cast<size_t>(i * st.M + k)], min_value);
            u[static_cast<size_t>(i * per_i + k)] = std::log(wk / w_last);
        }
    }
    return u;
}

inline AllocState appendPhaseFromComposition(
    const AllocState& st_old,
    const std::vector<double>& z,
    const std::vector<double>& x_new_phase,
    double beta_new,
    double z_min = Constants::MIN_MOLE_FRACTION)
{
    if (st_old.M < 1) {
        throw std::invalid_argument("appendPhaseFromComposition: st_old.M must be >= 1");
    }
    if (st_old.nc <= 0) {
        throw std::invalid_argument("appendPhaseFromComposition: st_old.nc must be > 0");
    }
    if (static_cast<int>(z.size()) != st_old.nc) {
        throw std::invalid_argument("appendPhaseFromComposition: z size mismatch");
    }
    if (static_cast<int>(x_new_phase.size()) != st_old.nc) {
        throw std::invalid_argument("appendPhaseFromComposition: x_new_phase size mismatch");
    }
    if (!(std::isfinite(beta_new) && beta_new >= 0.0)) {
        throw std::invalid_argument("appendPhaseFromComposition: beta_new must be finite and >= 0");
    }

    const std::vector<double> x_new = clampAndNormalize(x_new_phase, Constants::MIN_MOLE_FRACTION);

    const int Mold = st_old.M;
    const int M = Mold + 1;
    AllocState st_new;
    st_new.nc = st_old.nc;
    st_new.M = M;
    st_new.w.assign(static_cast<size_t>(st_new.nc * st_new.M), 0.0);

    for (int i = 0; i < st_old.nc; ++i) {
        const double zi = std::max(z[static_cast<size_t>(i)], z_min);
        double w_new_i = beta_new * x_new[static_cast<size_t>(i)] / zi;
        w_new_i = std::clamp(w_new_i, 0.0, 1.0 - 1e-12);

        for (int k = 0; k < Mold; ++k) {
            st_new.w[static_cast<size_t>(i * M + k)] =
                st_old.w[static_cast<size_t>(i * Mold + k)] * (1.0 - w_new_i);
        }
        st_new.w[static_cast<size_t>(i * M + (M - 1))] = w_new_i;

        // Normalize per-component to sum=1 (for safety).
        double s = 0.0;
        for (int k = 0; k < M; ++k) s += st_new.w[static_cast<size_t>(i * M + k)];
        if (!(std::isfinite(s) && s > 0.0)) {
            throw std::runtime_error("appendPhaseFromComposition: per-component allocation sum is not positive");
        }
        for (int k = 0; k < M; ++k) st_new.w[static_cast<size_t>(i * M + k)] /= s;
    }

    return st_new;
}

inline void computeBetasAndCompositions(
    const std::vector<double>& z,
    const AllocState& st,
    double beta_min,
    std::vector<double>& betas_out,
    std::vector<std::vector<double>>& x_out,
    double& penalty_out)
{
    betas_out.assign(static_cast<size_t>(st.M), 0.0);
    x_out.assign(static_cast<size_t>(st.M), std::vector<double>(static_cast<size_t>(st.nc), 0.0));
    penalty_out = 0.0;

    for (int k = 0; k < st.M; ++k) {
        double beta = 0.0;
        for (int i = 0; i < st.nc; ++i) {
            beta += z[static_cast<size_t>(i)] * st.w[static_cast<size_t>(i * st.M + k)];
        }
        betas_out[static_cast<size_t>(k)] = beta;
    }

    for (int k = 0; k < st.M; ++k) {
        double beta = betas_out[static_cast<size_t>(k)];
        if (beta < beta_min) {
            penalty_out += (beta_min - beta) * (beta_min - beta);
            beta = beta_min;
        }
        for (int i = 0; i < st.nc; ++i) {
            x_out[static_cast<size_t>(k)][static_cast<size_t>(i)] =
                z[static_cast<size_t>(i)] * st.w[static_cast<size_t>(i * st.M + k)] / beta;
        }
        x_out[static_cast<size_t>(k)] = clampAndNormalize(x_out[static_cast<size_t>(k)], Constants::MIN_MOLE_FRACTION);
    }

    // Renormalize betas (only matters if beta_min was applied).
    double sum = 0.0;
    for (double b : betas_out) sum += std::max(b, beta_min);
    if (sum > 0.0) {
        for (auto& b : betas_out) b = std::max(b, beta_min) / sum;
    }
}

inline void computeBetasAndCompositions(
    const std::vector<double>& z,
    const AllocState& st,
    double beta_min,
    double x_min,
    std::vector<double>& betas_out,
    std::vector<std::vector<double>>& x_out,
    double& penalty_out)
{
    betas_out.assign(static_cast<size_t>(st.M), 0.0);
    x_out.assign(static_cast<size_t>(st.M), std::vector<double>(static_cast<size_t>(st.nc), 0.0));
    penalty_out = 0.0;

    for (int k = 0; k < st.M; ++k) {
        double beta = 0.0;
        for (int i = 0; i < st.nc; ++i) {
            beta += z[static_cast<size_t>(i)] * st.w[static_cast<size_t>(i * st.M + k)];
        }
        betas_out[static_cast<size_t>(k)] = beta;
    }

    for (int k = 0; k < st.M; ++k) {
        double beta = betas_out[static_cast<size_t>(k)];
        if (beta < beta_min) {
            penalty_out += (beta_min - beta) * (beta_min - beta);
            beta = beta_min;
        }
        for (int i = 0; i < st.nc; ++i) {
            x_out[static_cast<size_t>(k)][static_cast<size_t>(i)] =
                z[static_cast<size_t>(i)] * st.w[static_cast<size_t>(i * st.M + k)] / beta;
        }
        x_out[static_cast<size_t>(k)] = clampAndNormalize(x_out[static_cast<size_t>(k)], x_min);
    }

    // Renormalize betas (only matters if beta_min was applied).
    double sum = 0.0;
    for (double b : betas_out) sum += std::max(b, beta_min);
    if (sum > 0.0) {
        for (auto& b : betas_out) b = std::max(b, beta_min) / sum;
    }
}

inline double phaseGibbsRT(const std::vector<double>& x, const std::vector<double>& lnphi, double x_min) {
    double g = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        g += x[i] * (safeLog(x[i], x_min) + lnphi[i]);
    }
    return g;
}

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_PHASE_ALLOCATIONS_H
