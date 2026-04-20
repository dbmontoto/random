/**
 * @file temperature_solve.h
 * @brief Shared temperature bracketing + bisection utilities for PH/PS outer solves (internal).
 */

#ifndef THERMO_EQUILIBRIUM_INTERNAL_TEMPERATURE_SOLVE_H
#define THERMO_EQUILIBRIUM_INTERNAL_TEMPERATURE_SOLVE_H

#include <algorithm>
#include <cmath>

namespace DMThermo {
namespace Equilibrium {
namespace Internal {

template <typename Result, typename ComputeFn, typename FvalFn>
inline bool bracketTemperatureByExpansion(
    double& T_lo,
    double& T_hi,
    Result& fr_lo,
    Result& fr_hi,
    ComputeFn compute,
    FvalFn fval,
    int max_expand_steps,
    double min_width)
{
    bool ok_lo = compute(T_lo, fr_lo);
    bool ok_hi = compute(T_hi, fr_hi);

    auto finiteF = [&](const Result& fr, double& out_f) -> bool {
        out_f = fval(fr);
        return std::isfinite(out_f);
    };

    for (int k = 0; k < max_expand_steps; ++k) {
        if (T_hi - T_lo < min_width) break;

        if (ok_lo && ok_hi) {
            double f_lo = 0.0;
            double f_hi = 0.0;
            const bool fin_lo = finiteF(fr_lo, f_lo);
            const bool fin_hi = finiteF(fr_hi, f_hi);
            if (!fin_lo) ok_lo = false;
            if (!fin_hi) ok_hi = false;
            if (ok_lo && ok_hi && !(f_lo * f_hi > 0.0)) {
                return true;
            }

            if (ok_lo && ok_hi) {
                // Expand in the direction that can produce a sign change.
                if (f_lo > 0.0 && f_hi > 0.0) {
                    T_lo = std::max(1.0, 0.8 * T_lo);
                    ok_lo = compute(T_lo, fr_lo);
                } else if (f_lo < 0.0 && f_hi < 0.0) {
                    T_hi = 1.2 * T_hi;
                    ok_hi = compute(T_hi, fr_hi);
                } else {
                    T_lo = std::max(1.0, 0.8 * T_lo);
                    T_hi = 1.2 * T_hi;
                    ok_lo = compute(T_lo, fr_lo);
                    ok_hi = compute(T_hi, fr_hi);
                }
                continue;
            }
        }

        // If an endpoint is invalid (e.g., out of correlation range), move it *inward*
        // toward the other endpoint (not outward) until we can evaluate f(T).
        if (!ok_lo && ok_hi) {
            const double mid = std::sqrt(std::max(1.0, T_lo) * std::max(1.0, T_hi));
            T_lo = std::max(1.0, 0.5 * (T_lo + mid));
            ok_lo = compute(T_lo, fr_lo);
            continue;
        }
        if (ok_lo && !ok_hi) {
            const double mid = std::sqrt(std::max(1.0, T_lo) * std::max(1.0, T_hi));
            T_hi = std::max(T_lo + min_width, 0.5 * (T_hi + mid));
            ok_hi = compute(T_hi, fr_hi);
            continue;
        }

        // If both fail, shrink the interval around the log-midpoint and try again.
        const double mid = std::sqrt(std::max(1.0, T_lo) * std::max(1.0, T_hi));
        T_lo = std::max(1.0, 0.9 * mid);
        T_hi = std::max(T_lo + min_width, 1.1 * mid);
        ok_lo = compute(T_lo, fr_lo);
        ok_hi = compute(T_hi, fr_hi);
    }

    if (!(ok_lo && ok_hi)) return false;
    double f_lo = 0.0;
    double f_hi = 0.0;
    if (!finiteF(fr_lo, f_lo) || !finiteF(fr_hi, f_hi)) return false;
    return !(f_lo * f_hi > 0.0);
}

template <typename Result, typename ComputeFn, typename FvalFn>
inline bool bracketTemperatureByLogScan(
    double T_lo,
    double T_hi,
    Result& fr_lo,
    Result& fr_hi,
    ComputeFn compute,
    FvalFn fval,
    int num_samples)
{
    Result prev;
    double prev_f = 0.0;
    bool have_prev = false;

    const double log_lo = std::log(T_lo);
    const double log_hi = std::log(T_hi);
    for (int i = 0; i <= num_samples; ++i) {
        const double t = std::exp(log_lo + (log_hi - log_lo) * (static_cast<double>(i) / num_samples));
        Result cur;
        if (!compute(t, cur)) continue;
        const double fc = fval(cur);
        if (!have_prev) {
            prev = cur;
            prev_f = fc;
            have_prev = true;
            continue;
        }
        if (prev_f * fc < 0.0) {
            fr_lo = prev;
            fr_hi = cur;
            return true;
        }
        prev = cur;
        prev_f = fc;
    }
    return false;
}

template <typename Result>
struct TemperatureBisectionResult {
    Result result{};
    bool converged = false;
    int iterations = 0;
};

template <typename Result, typename ComputeFn, typename FvalFn>
inline TemperatureBisectionResult<Result> solveTemperatureBisection(
    double a,
    double b,
    Result fr_lo,
    ComputeFn compute,
    FvalFn fval,
    int max_iterations,
    double tolerance)
{
    Result fr_mid;
    for (int it = 0; it < max_iterations; ++it) {
        const double m = 0.5 * (a + b);
        if (!compute(m, fr_mid)) {
            b = m;
            continue;
        }

        const double fm = fval(fr_mid);
        if (std::abs(fm) < tolerance) {
            TemperatureBisectionResult<Result> out;
            out.result = fr_mid;
            out.converged = true;
            out.iterations = it + 1;
            return out;
        }

        const double fa = fval(fr_lo);
        if (fa * fm < 0.0) {
            b = m;
        } else {
            a = m;
            fr_lo = fr_mid;
        }
    }

    TemperatureBisectionResult<Result> out;
    out.result = fr_mid;
    out.converged = false;
    out.iterations = max_iterations;
    return out;
}

} // namespace Internal
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_INTERNAL_TEMPERATURE_SOLVE_H
