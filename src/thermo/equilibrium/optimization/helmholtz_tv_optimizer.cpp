/**
 * @file helmholtz_tv_optimizer.cpp
 * @brief Skeleton implementation for HelmholtzTVOptimizer.
 */

#include "thermo/equilibrium/optimization/helmholtz_tv_optimizer.h"
#include "thermo/core/constants.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/density_selection.h"
#include "thermo/equilibrium/internal/tv_outer_solve.h"
#include "thermo/equilibrium/internal/tv_support.h"
#include "thermo/equilibrium/optimization/phase_allocations.h"
#include "thermo/equilibrium/optimization/internal/helmholtz_tv_direct_outer_loop_problem.h"
#include "thermo/equilibrium/optimization/variable_blocks.h"
#include "thermo/equilibrium/optimization/outer_loop_helpers.h"
#include "thermo/equilibrium/optimization/phase_postprocess.h"
#include "thermo/equilibrium/optimization/phase_stability_outer_loop.h"
#include "thermo/equilibrium/optimization/variable_layout.h"
#include "thermo/factory/numerics_factory.h"
#include "thermo/factory/stability_factory.h"
#include "thermo/equilibrium/optimization/gibbs_tp_optimizer.h"
#include "thermo/equilibrium/optimization/transforms.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

HelmholtzTVOptimizer::HelmholtzTVOptimizer(EOSPtr eos, Stability::StabilityAnalyzerPtr stability)
    : eos_(std::move(eos)), stability_(std::move(stability))
{
    if (!eos_) {
        throw std::invalid_argument("HelmholtzTVOptimizer: EOS pointer is null");
    }
    if (!stability_) {
        stability_ = Factory::StabilityFactory::create(eos_);
    }
}

EquilibriumState HelmholtzTVOptimizer::minimizeTV(
    double T,
    double V,
    const std::vector<double>& z_in,
    const Config::EquilibriumOptimizerConfig& cfg) const
{
    EquilibriumState out;
    out.method_used = "HelmholtzOptimizer(T,V)";
    out.temperature = T;
    out.z = z_in;

    const int nc = eos_->numComponents();
    if (static_cast<int>(z_in.size()) != nc) {
        out.converged = false;
        out.message = "Composition size mismatch";
        return out;
    }
    if (!(T > 0.0 && V > 0.0 && std::isfinite(V))) {
        out.converged = false;
        out.message = "Invalid T/V";
        return out;
    }

    const std::vector<double> z = clampAndNormalize(z_in, Constants::MIN_MOLE_FRACTION);
    const double rho_total = 1.0 / V;
    const double P_ig = rho_total * Constants::GAS_CONSTANT * T;

    // Fast path: single-phase at the specified (T, rho_total, z).
    // For a TV problem, if the system is stable at this density, we can return directly.
    const double P_rho = eos_->pressure(T, rho_total, z);
    if (cfg.tv_solve_strategy == Config::TVSolveStrategy::BrentOuter &&
        std::isfinite(P_rho) && P_rho > 0.0 && cfg.use_phase_stability_outer_loop) {
        const auto stab_cfg = stabilityConfigFromEquilibriumConfig(cfg);
        const auto stab = stability_->analyze(T, P_rho, z, stab_cfg);
        if (stab.is_stable) {
            PhaseState ps;
            ps.type = stability_->classifyPhase(T, rho_total, z);
            ps.density = rho_total;
            ps.compressibility = eos_->compressibility(T, ps.density, z);
            ps.x = z;
            ps.phi = eos_->fugacityCoefficients(T, ps.density, z);
            ps.fraction = 1.0;

            out.phases = {ps};
            out.pressure = P_rho;
            out.objective_value = Internal::phaseHelmholtzRT(*eos_, T, ps.density, ps.x, cfg.composition_min);
            out.infeasibility = 0.0;
            out.score_value = out.objective_value;
            out.converged = true;
            out.message = "Single phase (stable)";
            return out;
        }
    }

    // Outer solve in P: find P such that V_calc(T,P,z) matches target V.
    if (cfg.tv_solve_strategy == Config::TVSolveStrategy::DirectOptimizer) {
        // Optional single-phase stability shortcut (still reports the selected strategy).
        if (std::isfinite(P_rho) && P_rho > 0.0 && cfg.use_phase_stability_outer_loop) {
            const auto stab_cfg = stabilityConfigFromEquilibriumConfig(cfg);
            const auto stab = stability_->analyze(T, P_rho, z, stab_cfg);
            if (stab.is_stable) {
                PhaseState ps;
                ps.type = stability_->classifyPhase(T, rho_total, z);
                ps.density = rho_total;
                ps.compressibility = eos_->compressibility(T, ps.density, z);
                ps.x = z;
                ps.phi = eos_->fugacityCoefficients(T, ps.density, z);
                ps.fraction = 1.0;

                out.phases = {ps};
                out.pressure = P_rho;
                out.objective_value = Internal::phaseHelmholtzRT(*eos_, T, ps.density, ps.x, cfg.composition_min);
                out.infeasibility = 0.0;
                out.score_value = out.objective_value;
                out.converged = true;
                out.message = "Single phase (stable)";
                out.method_used = "HelmholtzOptimizer(T,V) (direct penalty optimization)";
                return out;
            }
        }

        // Direct TV optimization: minimize a Helmholtz-based objective with explicit penalties.
        // This path is EOS-agnostic and uses penalty terms to satisfy volume/mechanical equilibrium.
        // Phase creation is optionally handled via an outer stability loop (using a provisional P).
        const int max_phases = std::clamp(cfg.max_phases, 1, 3);

        auto optimizer = Factory::NumericsFactory::createBFGS();

        // Seed allocations and densities from a single TP equilibrium at an estimated pressure.
        GibbsTPOptimizer gibbs_seed(eos_, stability_);
        const double P_seed = (std::isfinite(P_rho) && P_rho > 0.0) ? P_rho : std::max(P_ig, 1e5);

        int M = std::min(2, max_phases);
        AllocState alloc;
        alloc.nc = nc;
        alloc.M = M;
        alloc.w.assign(static_cast<size_t>(nc * M), 0.0);

        std::vector<double> rho0(static_cast<size_t>(M), rho_total);

        try {
            Config::EquilibriumOptimizerConfig tp_cfg = cfg;
            tp_cfg.phase_detection_space = Config::PhaseDetection::TP;
            tp_cfg.max_phases = M;
            const auto st = gibbs_seed.minimizeTP(T, P_seed, z, tp_cfg);
            if (st.converged && static_cast<int>(st.phases.size()) == M) {
                for (int k = 0; k < M; ++k) {
                    const auto& ph = st.phases[static_cast<size_t>(k)];
                    rho0[static_cast<size_t>(k)] = std::max(ph.density, Constants::MIN_DENSITY);
                    for (int i = 0; i < nc; ++i) {
                        const double zi = std::max(z[static_cast<size_t>(i)], Constants::MIN_MOLE_FRACTION);
                        const double wk = std::clamp((ph.fraction * ph.x[static_cast<size_t>(i)]) / zi, 0.0, 1.0);
                        alloc.w[static_cast<size_t>(i * M + k)] = wk;
                    }
                }
                // Renormalize allocations per component to sum to 1.
                for (int i = 0; i < nc; ++i) {
                    double s = 0.0;
                    for (int k = 0; k < M; ++k) s += alloc.w[static_cast<size_t>(i * M + k)];
                    if (s > 0.0) {
                        for (int k = 0; k < M; ++k) alloc.w[static_cast<size_t>(i * M + k)] /= s;
                    } else {
                        for (int k = 0; k < M; ++k) alloc.w[static_cast<size_t>(i * M + k)] = 1.0 / M;
                    }
                }

                // Scale densities to approximately satisfy the target volume.
                double V0 = 0.0;
                for (int k = 0; k < M; ++k) {
                    const double b = std::max(st.phases[static_cast<size_t>(k)].fraction, cfg.beta_min);
                    V0 += b / std::max(rho0[static_cast<size_t>(k)], Constants::MIN_DENSITY);
                }
                if (std::isfinite(V0) && V0 > 0.0) {
                    const double alpha = V0 / V;
                    if (std::isfinite(alpha) && alpha > 0.0) {
                        for (auto& r : rho0) r = std::clamp(r * alpha, Constants::MIN_DENSITY, Constants::MAX_DENSITY);
                    }
                }
            } else {
                // Symmetry-breaking seed: vapor-like and liquid-like densities around rho_total.
                if (M == 2) {
                    rho0[0] = std::clamp(0.2 * rho_total, Constants::MIN_DENSITY, Constants::MAX_DENSITY);
                    rho0[1] = std::clamp(5.0 * rho_total, Constants::MIN_DENSITY, Constants::MAX_DENSITY);
                }
                for (int i = 0; i < nc; ++i) {
                    alloc.w[static_cast<size_t>(i * M + 0)] = 0.5;
                    if (M == 2) alloc.w[static_cast<size_t>(i * M + 1)] = 0.5;
                }
            }
        } catch (...) {
            if (M == 2) {
                rho0[0] = std::clamp(0.2 * rho_total, Constants::MIN_DENSITY, Constants::MAX_DENSITY);
                rho0[1] = std::clamp(5.0 * rho_total, Constants::MIN_DENSITY, Constants::MAX_DENSITY);
            }
            for (int i = 0; i < nc; ++i) {
                alloc.w[static_cast<size_t>(i * M + 0)] = 0.5;
                if (M == 2) alloc.w[static_cast<size_t>(i * M + 1)] = 0.5;
            }
        }

        Internal::DirectTVOuterLoopProblem problem{*eos_, *stability_, T, V, rho_total, z, cfg, nc};
        std::vector<double> vars = problem.makeVars(alloc, rho0);
        const auto res = Numerics::Optimization::runPhaseStabilityOuterLoop(
            *optimizer, problem, M, max_phases, std::move(vars), outerLoopConfigFrom(cfg));
        out = res.best_state;
        if (out.phases.empty()) {
            out.converged = false;
            out.message = "HelmholtzTVOptimizer: all phases collapsed";
            out.method_used = "HelmholtzOptimizer(T,V) (direct penalty optimization)";
            return out;
        }

        out.method_used = "HelmholtzOptimizer(T,V) (direct penalty optimization)";
        out.message = out.converged ? "Converged" : (out.message.empty() ? "Did not converge" : out.message);
        return out;
    }

    GibbsTPOptimizer gibbs(eos_, stability_);
    auto root_finder = Factory::NumericsFactory::createBrent();

    struct Eval {
        bool ok = false;
        double Vcalc = 0.0;
        EquilibriumState state;
    };

    auto guessFor = [&](PhaseType t, double P, double fallback) -> double {
        (void)t;
        (void)P;
        if (cfg.phase_detection_space != Config::PhaseDetection::TRho) return std::numeric_limits<double>::quiet_NaN();

        if (std::isfinite(fallback) && fallback > 0.0) return fallback;
        return std::numeric_limits<double>::quiet_NaN();
    };

    auto evalAtP = [&](double P) -> Eval {
        Eval e;
        if (!(std::isfinite(P) && P > 0.0)) return e;
        const Core::Units::PropertyVariable pressure_tolerance(cfg.density_newton_tol, Core::Units::Unit::PA);

        Config::EquilibriumOptimizerConfig tp_cfg = cfg;
        tp_cfg.phase_detection_space = Config::PhaseDetection::TP; // density handled below for volume continuity
        tp_cfg.max_phases = std::clamp(cfg.max_phases, 1, 3);

        EquilibriumState st;
        try {
            st = gibbs.minimizeTP(T, P, z, tp_cfg);
        } catch (...) {
            return e;
        }
        if (st.phases.empty()) return e;

        double Vcalc = 0.0;
        try {
            const bool ok = Equilibrium::Internal::postProcessTVEvalAtTP(
                st.phases,
                *eos_,
                *stability_,
                T,
                P,
                cfg.phase_detection_space,
                cfg.max_density_newton_iters,
                pressure_tolerance,
                cfg.beta_min,
                [&](PhaseType type, double P_at, double rho_cur) { return guessFor(type, P_at, rho_cur); },
                [&](PhaseType t) { return phaseTypeRank(t); },
                [&](const std::vector<PhaseState>& phases) { return Internal::volumeFromPhases(phases); },
                Vcalc);
            if (!ok) return e;
        } catch (...) {
            return e;
        }

        e.ok = true;
        e.Vcalc = Vcalc;
        e.state = std::move(st);
        e.state.pressure = P;
        return e;
    };

    // Initial guess: use EOS at rho_total if finite; otherwise ideal-gas.
    const double P0 = Equilibrium::Internal::tvInitialPressureGuess(P_rho, P_ig, 1e5);

    const auto rr = Equilibrium::Internal::solveTVPressureBrentOuter(
        *root_finder,
        evalAtP,
        V,
        P0,
        1e-3,
        cfg.tv_bracket_logp_step,
        cfg.tv_bracket_steps,
        cfg.tv_p_min,
        cfg.tv_p_max,
        cfg.optimizer.tolerance,
        "HelmholtzTVOptimizer");
    if (!rr.converged) {
        out.converged = false;
        out.message = rr.message;
        return out;
    }
    const double P_star = rr.P_star;
    out.iterations = rr.iterations;
    const auto final_eval = rr.final_eval;

    out = final_eval.state;

    // Compute Helmholtz objective (dimensionless A/RT) from final phases.
    double A = 0.0;
    for (const auto& ph : out.phases) {
        A += ph.fraction * Internal::phaseHelmholtzRT(*eos_, T, ph.density, ph.x, cfg.composition_min);
    }
    out.objective_value = A;
    out.infeasibility = 0.0;
    out.score_value = out.objective_value;

    const double Vcalc = Internal::volumeFromPhases(out.phases);
    const double spread = Internal::phasePressureSpread(*eos_, T, out.phases);
    (void)spread;
    Equilibrium::Internal::tvFinalizeStateForTVSolve(
        out,
        P_star,
        Vcalc,
        V,
        rr.volume_abs_tolerance,
        "HelmholtzOptimizer(T,V) (Brent outer + Gibbs(T,P) inner)");
    return out;
}

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo
