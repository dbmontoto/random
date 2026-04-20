/**
 * @file helmholtz_tv_direct_outer_loop_problem.h
 * @brief Helmholtz TV direct-penalty outer-loop problem definition (internal).
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_INTERNAL_HELMHOLTZ_TV_DIRECT_OUTER_LOOP_PROBLEM_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_INTERNAL_HELMHOLTZ_TV_DIRECT_OUTER_LOOP_PROBLEM_H

#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/core/constants.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/density_selection.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/equilibrium/optimization/phase_allocations.h"
#include "thermo/equilibrium/optimization/phase_postprocess.h"
#include "thermo/equilibrium/optimization/phase_stability_outer_loop.h"
#include "thermo/equilibrium/optimization/variable_blocks.h"
#include "thermo/equilibrium/optimization/variable_layout.h"
#include "thermo/numerics/ioptimizer.h"
#include "thermo/numerics/optimization_problem.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {
namespace Internal {

inline double phaseHelmholtzRT(const EOS& eos, double T, double rho, const std::vector<double>& x, double x_min) {
    // A/RT = a_residual + ln(rho) + sum x_i ln x_i + const(T)
    // const(T) cancels between phases, so we omit it.
    double a = eos.residualHelmholtz(T, rho, x);
    double mix = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        mix += x[i] * safeLog(x[i], x_min);
    }
    return a + std::log(std::max(rho, Constants::MIN_DENSITY)) + mix;
}

inline double pressureCommonPenalty(const std::vector<double>& Pk) {
    if (Pk.empty()) return 0.0;
    double mean = 0.0;
    for (double v : Pk) mean += v;
    mean /= static_cast<double>(Pk.size());
    double var = 0.0;
    for (double v : Pk) {
        const double d = v - mean;
        var += d * d;
    }
    return var / static_cast<double>(Pk.size());
}

inline double volumeFromPhases(const std::vector<PhaseState>& phases) {
    double V = 0.0;
    for (const auto& ph : phases) {
        if (!(ph.density > 0.0)) continue;
        V += ph.fraction / ph.density;
    }
    return V;
}

inline double phasePressureSpread(const EOS& eos, double T, const std::vector<PhaseState>& phases) {
    if (phases.size() <= 1) return 0.0;
    std::vector<double> Pk;
    Pk.reserve(phases.size());
    for (const auto& ph : phases) {
        const double P = eos.pressure(T, ph.density, ph.x);
        if (std::isfinite(P)) Pk.push_back(P);
    }
    if (Pk.size() <= 1) return 0.0;
    double mean = 0.0;
    for (double p : Pk) mean += p;
    mean /= static_cast<double>(Pk.size());
    double maxabs = 0.0;
    for (double p : Pk) maxabs = std::max(maxabs, std::abs(p - mean));
    return maxabs;
}

struct DirectTVOuterLoopProblem {
    using State = EquilibriumState;

    const EOS& eos;
    const Stability::IStabilityAnalyzer& stability;
    double T = 0.0;
    double V = 0.0;
    double rho_total = 0.0;
    const std::vector<double>& z;
    const Config::EquilibriumOptimizerConfig& cfg;
    int nc = 0;

    std::vector<double> makeVars(const AllocState& a0, const std::vector<double>& rho_seed) const {
        const auto u0 = packLogitsFromAllocations(a0, cfg.composition_min);
        std::vector<double> v = u0;
        const auto urho = encodeBoundedVector(
            clampToBounds(rho_seed, Constants::MIN_DENSITY, Constants::MAX_DENSITY),
            Constants::MIN_DENSITY,
            Constants::MAX_DENSITY);
        v.insert(v.end(), urho.begin(), urho.end());
        return v;
    }

    Numerics::Optimization::ProblemEvaluation evaluate(int Mcur, const std::vector<double>& v) const {
        const auto layout = AllocRhoLayout::make(nc, Mcur);
        if (!layout.matches(v)) {
            return Numerics::Optimization::ProblemEvaluation{
                std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()
            };
        }

        const auto ualloc_u = layout.allocSlice(v);
        const auto urho_u = layout.rhoSlice(v);

        const AllocState st = unpackAllocations(ualloc_u.data, ualloc_u.size, nc, Mcur, cfg.composition_min);

        std::vector<double> betas;
        std::vector<std::vector<double>> x;
        double infeas = 0.0;
        computeBetasAndCompositions(z, st, cfg.beta_min, cfg.composition_min, betas, x, infeas);

        const auto rhos = decodeBoundedVector(urho_u, Constants::MIN_DENSITY, Constants::MAX_DENSITY);

        double A = 0.0;
        std::vector<double> Pk;
        Pk.reserve(static_cast<size_t>(Mcur));
        double Vcalc = 0.0;

        for (int k = 0; k < Mcur; ++k) {
            const double rho = std::max(rhos[static_cast<size_t>(k)], Constants::MIN_DENSITY);
            Vcalc += betas[static_cast<size_t>(k)] / rho;

            const double a = phaseHelmholtzRT(eos, T, rho, x[static_cast<size_t>(k)], cfg.composition_min);
            if (!std::isfinite(a)) {
                return Numerics::Optimization::ProblemEvaluation{
                    std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()
                };
            }
            A += betas[static_cast<size_t>(k)] * a;

            const double P = eos.pressure(T, rho, x[static_cast<size_t>(k)]);
            if (std::isfinite(P)) Pk.push_back(P);
        }

        if (!(std::isfinite(Vcalc) && Vcalc > 0.0)) {
            return Numerics::Optimization::ProblemEvaluation{
                std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()
            };
        }

        const double dv = Vcalc - V;
        const double vol_pen = cfg.volume_penalty * dv * dv;
        const double pres_pen = cfg.pressure_penalty * pressureCommonPenalty(Pk);

        Numerics::Optimization::ProblemEvaluation out;
        out.objective = A + vol_pen + pres_pen;
        out.infeasibility = infeas;
        return out;
    }

    double penaltyWeight() const { return cfg.infeasibility_penalty; }

    double score(const EquilibriumState& st) const { return st.score_value; }

    EquilibriumState decode(int Mcur, const Numerics::Optimization::OptimizationResult& opt) const {
        EquilibriumState st;
        st.temperature = T;
        st.z = z;
        st.method_used = "HelmholtzOptimizer(T,V) (direct penalty optimization)";
        st.message = opt.message;
        st.iterations = opt.iterations;
        st.converged = opt.converged;

        const auto layout = AllocRhoLayout::make(nc, Mcur);
        if (!layout.matches(opt.x_optimal)) {
            throw std::runtime_error("DirectTVOuterLoopProblem::decode: optimizer variable size mismatch");
        }

        const auto ualloc_u = layout.allocSlice(opt.x_optimal);
        const auto urho_u = layout.rhoSlice(opt.x_optimal);

        const AllocState a = unpackAllocations(ualloc_u.data, ualloc_u.size, nc, Mcur, cfg.composition_min);

        std::vector<double> betas;
        std::vector<std::vector<double>> xk;
        double penalty = 0.0;
        computeBetasAndCompositions(z, a, cfg.beta_min, cfg.composition_min, betas, xk, penalty);

        const auto rhos = decodeBoundedVector(urho_u, Constants::MIN_DENSITY, Constants::MAX_DENSITY);

        double Pmean = 0.0;
        int Pcount = 0;
        std::vector<PhaseState> phases;
        phases.reserve(static_cast<size_t>(Mcur));
        for (int k = 0; k < Mcur; ++k) {
            PhaseState ps;
            ps.x = xk[static_cast<size_t>(k)];
            ps.density = rhos[static_cast<size_t>(k)];
            ps.type = stability.classifyPhase(T, ps.density, ps.x);
            ps.compressibility = eos.compressibility(T, ps.density, ps.x);
            ps.phi = eos.fugacityCoefficients(T, ps.density, ps.x);
            ps.fraction = betas[static_cast<size_t>(k)];
            phases.push_back(std::move(ps));

            const double Pk = eos.pressure(T, phases.back().density, phases.back().x);
            if (std::isfinite(Pk)) {
                Pmean += Pk;
                ++Pcount;
            }
        }

        PhasePostProcessOptions pp;
        pp.prune_fraction = cfg.beta_min * 10.0;
        pp.allow_empty = false;
        pp.merge_duplicates = true;
        pp.sort_by_type = true;
        pp.canonicalize_liquid = true;
        if (!postProcessPhases(phases, pp)) {
            st.converged = false;
            st.message = "All phases collapsed";
            st.objective_value = std::numeric_limits<double>::infinity();
            st.infeasibility = std::numeric_limits<double>::infinity();
            st.score_value = std::numeric_limits<double>::infinity();
            return st;
        }

        st.phases = std::move(phases);
        st.pressure = (Pcount > 0) ? (Pmean / Pcount) : std::numeric_limits<double>::quiet_NaN();

        const auto eval = evaluate(Mcur, opt.x_optimal);
        st.objective_value = eval.objective;
        st.infeasibility = eval.infeasibility;
        st.score_value = eval.penalized(penaltyWeight());

        const double Vcalc = volumeFromPhases(st.phases);
        const double vol_rel =
            (std::isfinite(Vcalc) && Vcalc > 0.0) ? (std::abs(Vcalc - V) / V) : std::numeric_limits<double>::infinity();
        st.converged = st.converged && std::isfinite(Vcalc) && (vol_rel < 5e-3);
        if (st.converged) st.message = "Converged";

        return st;
    }

    std::optional<std::vector<double>> unstableComposition(int, const EquilibriumState& candidate) const {
        const double P_common = candidate.pressure;
        if (!(std::isfinite(P_common) && P_common > 0.0)) return std::nullopt;
        return findUnstableComposition(stability, T, P_common, candidate.phases, cfg);
    }

    std::vector<double> expandInitialGuess(
        int Mold,
        int Mnew,
        const Numerics::Optimization::OptimizationResult& opt,
        const EquilibriumState& candidate,
        const std::vector<double>& unstable_w) const
    {
        const auto layout_old = AllocRhoLayout::make(nc, Mold);
        if (!layout_old.matches(opt.x_optimal)) {
            throw std::runtime_error("DirectTVOuterLoopProblem::expandInitialGuess: optimizer variable size mismatch");
        }

        const auto ualloc_u = layout_old.allocSlice(opt.x_optimal);
        const auto urho_u = layout_old.rhoSlice(opt.x_optimal);

        const AllocState st_old = unpackAllocations(ualloc_u.data, ualloc_u.size, nc, Mold, cfg.composition_min);

        std::vector<double> rho_seed_new = decodeBoundedVector(urho_u, Constants::MIN_DENSITY, Constants::MAX_DENSITY);
        rho_seed_new.resize(static_cast<size_t>(Mnew), rho_total);

        const double P_common = candidate.pressure;
        double rho_new = rho_total;
        try {
            const Core::Units::PropertyVariable pressure_tolerance(cfg.density_newton_tol, Core::Units::Unit::PA);
            const auto ev = Detail::evalAtTP(
                eos, stability, T, P_common, unstable_w, hintForPhaseIndex(Mnew - 1),
                cfg.phase_detection_space, rho_total, cfg.max_density_newton_iters, pressure_tolerance
            );
            if (std::isfinite(ev.rho) && ev.rho > 0.0) rho_new = ev.rho;
        } catch (...) {
            const double rho_ig = std::max(P_common / (Constants::GAS_CONSTANT * T), Constants::MIN_DENSITY);
            rho_new = std::clamp(200.0 * rho_ig, Constants::MIN_DENSITY, Constants::MAX_DENSITY);
        }
        rho_seed_new[static_cast<size_t>(Mnew - 1)] = rho_new;

        const AllocState st_new = appendPhaseFromComposition(st_old, z, unstable_w, 1e-3);
        return makeVars(st_new, rho_seed_new);
    }
};

} // namespace Internal
} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_INTERNAL_HELMHOLTZ_TV_DIRECT_OUTER_LOOP_PROBLEM_H
