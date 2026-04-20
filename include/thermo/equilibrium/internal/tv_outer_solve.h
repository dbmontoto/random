/**
 * @file tv_outer_solve.h
 * @brief Shared TV outer-solve wrapper around log(P) Brent solves (internal).
 */

#ifndef THERMO_EQUILIBRIUM_INTERNAL_TV_OUTER_SOLVE_H
#define THERMO_EQUILIBRIUM_INTERNAL_TV_OUTER_SOLVE_H

#include "thermo/equilibrium/internal/log_space_solve.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>

namespace DMThermo {
namespace Equilibrium {
namespace Internal {

inline double tvInitialPressureGuess(double P_rho, double P_ig, double fallback) {
    if (std::isfinite(P_rho) && P_rho > 0.0) return P_rho;
    if (std::isfinite(P_ig) && P_ig > 0.0) return P_ig;
    if (std::isfinite(fallback) && fallback > 0.0) return fallback;
    return 1e5;
}

inline bool tvIsConverged(double Vcalc, double Vtarget, double volume_abs_tolerance, double slack_factor = 10.0) {
    return std::isfinite(Vcalc) && (std::abs(Vcalc - Vtarget) <= volume_abs_tolerance * slack_factor);
}

inline std::string tvConvergenceMessage(bool converged) {
    return converged ? "Converged" : "Did not converge";
}

template <typename State>
inline void tvFinalizeStateForTVSolve(
    State& st,
    double P_star,
    double Vcalc,
    double Vtarget,
    double volume_abs_tolerance,
    const char* method_used)
{
    st.pressure = P_star;
    st.method_used = method_used;
    st.converged = tvIsConverged(Vcalc, Vtarget, volume_abs_tolerance);
    st.message = tvConvergenceMessage(st.converged);
}

template <typename Eval>
struct TVOuterSolveResult {
    bool converged = false; // outer solve + final evaluation succeeded
    double P_star = 0.0;
    int iterations = 0;
    double volume_abs_tolerance = 0.0; // absolute |Vcalc - V| tolerance used by the root solver
    std::string message;
    Eval final_eval{};
};

template <typename EvalAtPFn>
inline auto solveTVPressureBrentOuter(
    const Numerics::RootFinding::IScalarRootFinder& root_finder,
    EvalAtPFn evalAtP,
    double V,
    double P0,
    double P0_min,
    double tv_bracket_logp_step,
    int tv_bracket_steps,
    double tv_p_min,
    double tv_p_max,
    double optimizer_tolerance,
    const char* label)
    -> TVOuterSolveResult<std::decay_t<decltype(evalAtP(P0))>>
{
    using Eval = std::decay_t<decltype(evalAtP(P0))>;
    TVOuterSolveResult<Eval> out;

    const double P0_safe = (std::isfinite(P0) && P0 > 0.0) ? P0 : P0_min;
    const double logP0 = std::log(std::max(P0_safe, P0_min));
    const double step = std::clamp(tv_bracket_logp_step, 0.05, 2.0) * std::log(10.0);
    const int max_steps = std::clamp(tv_bracket_steps, 5, 300);
    const double Pmin = std::max(tv_p_min, 1e-12);
    const double Pmax = std::max(Pmin * 10.0, tv_p_max);

    Config::RootFinderConfig rfc = Config::RootFinderConfig::robust();
    // The inner TP equilibrium solve can introduce numerical noise into V(P), especially
    // near phase boundaries. Using an extremely tight volume tolerance can cause the
    // outer Brent solve to stagnate even when the final TV closure is already within
    // practical engineering tolerances. Keep a reasonable floor here.
    const double vol_rel_tol = std::max(5e-4, optimizer_tolerance);
    rfc.tolerance = std::max(1e-12, vol_rel_tol * V);
    rfc.rel_tolerance = 0.0;
    rfc.max_iterations = 160;
    out.volume_abs_tolerance = rfc.tolerance;

    auto f_logP = [&](double logP) -> double {
        try {
            const double P = std::exp(logP);
            const auto e = evalAtP(P);
            if (!e.ok) return std::numeric_limits<double>::quiet_NaN();
            return e.Vcalc - V;
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    };

    const auto rr = solveLogSpaceRootBrent(root_finder, f_logP, logP0, step, max_steps, Pmin, Pmax, rfc, label);
    if (!rr.converged) {
        out.converged = false;
        out.message = rr.message;
        return out;
    }

    out.P_star = std::exp(rr.root);
    out.iterations = rr.iterations;

    try {
        out.final_eval = evalAtP(out.P_star);
    } catch (...) {
        out.converged = false;
        out.message = std::string(label) + ": final evaluation threw";
        return out;
    }
    if (!out.final_eval.ok) {
        out.converged = false;
        out.message = std::string(label) + ": final evaluation failed";
        return out;
    }

    out.converged = true;
    return out;
}

} // namespace Internal
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_INTERNAL_TV_OUTER_SOLVE_H
