/**
 * @file reactive_flash.cpp
 * @brief Implementation of ReactiveFlashSolver.
 */

#include "thermo/reactions/reactive_flash.h"

#include "thermo/core/constants.h"
#include "thermo/equilibrium/density_selection.h"
#include "thermo/equilibrium/internal/flash_result_support.h"
#include "thermo/equilibrium/internal/hs_support.h"
#include "thermo/equilibrium/internal/temperature_solve.h"
#include "thermo/equilibrium/internal/tv_outer_solve.h"
#include "thermo/equilibrium/internal/tv_support.h"
#include "thermo/equilibrium/thermo_contributions.h"
#include "thermo/equilibrium/optimization/reactive_gibbs_tp_optimizer.h"
#include "thermo/factory/numerics_factory.h"
#include "thermo/factory/stability_factory.h"
#include "thermo/reactions/thermochemistry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace DMThermo {
namespace Reactions {

namespace {

int rankPhase(PhaseType t) {
    switch (t) {
        case PhaseType::Vapor: return 0;
        case PhaseType::Liquid: return 1;
        case PhaseType::Liquid1: return 1;
        case PhaseType::Liquid2: return 2;
        default: return 10;
    }
}

double volumeFromPhases(const std::vector<PhaseState>& phases) {
    double V = 0.0;
    for (const auto& ph : phases) {
        if (!(std::isfinite(ph.fraction) && std::isfinite(ph.density) && ph.density > 0.0)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        V += ph.fraction / ph.density;
    }
    return V;
}

} // namespace

ReactiveFlashSolver::ReactiveFlashSolver(
    EOSPtr eos,
    Equilibrium::Stability::StabilityAnalyzerPtr stability)
    : eos_(std::move(eos))
    , stability_(std::move(stability))
{
    if (!eos_) {
        throw std::invalid_argument("ReactiveFlashSolver: EOS is null");
    }
    if (!stability_) {
        stability_ = Factory::StabilityFactory::create(eos_);
    }
}

ReactiveFlashResult ReactiveFlashSolver::solveTP(
    double T,
    double P,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const std::vector<double>& mu0,
    const ReactiveFlashConfig& cfg) const
{
    if (!(std::isfinite(T) && T > 0.0 && std::isfinite(P) && P > 0.0)) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTP: invalid T/P");
    }
    if (static_cast<int>(n0.size()) != eos_->numComponents()) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTP: n0 size mismatch");
    }
    if (static_cast<int>(mu0.size()) != eos_->numComponents()) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTP: mu0 size mismatch");
    }
    if (!(std::isfinite(cfg.standard_pressure) && cfg.standard_pressure > 0.0)) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTP: invalid standard_pressure");
    }
    if (!(std::isfinite(cfg.min_moles) && cfg.min_moles >= 0.0)) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTP: invalid min_moles");
    }

    Equilibrium::Optimization::ReactiveGibbsTPOptimizer opt(eos_, stability_);
    Equilibrium::Optimization::ReactiveGibbsTPConfig ocfg;
    ocfg.standard_pressure = cfg.standard_pressure;
    ocfg.min_moles = cfg.min_moles;
    ocfg.equilibrium = cfg.equilibrium;

    auto st = opt.minimizeTPReactive(T, P, n0, system, mu0, ocfg);

    ReactiveFlashResult out;
    out.temperature = T;
    out.pressure = P;
    out.converged = st.converged;
    out.iterations = st.iterations;
    out.message = st.message;
    out.n0 = n0;
    out.z = st.z;
    out.flash = Equilibrium::Internal::flashResultFromState(st);

    out.extents = st.reaction_extents.value_or(std::vector<double>{});
    if (!out.extents.empty()) {
        out.n = system.molesFromExtents(n0, out.extents);
    }

    return out;
}

ReactiveFlashResult ReactiveFlashSolver::solveTPFromDatabanks(
    double T,
    double P,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    const ReactiveFlashConfig& cfg,
    Data::Diagnostics* diag) const
{
    const std::vector<double> mu0 = Thermochemistry::mu0DB(
        cfg.mu0_model,
        databanks,
        mixture,
        T,
        cfg.mu0_reference_temperature,
        diag);
    return solveTP(T, P, n0, system, mu0, cfg);
}

ReactiveFlashResult ReactiveFlashSolver::solveTV(
    double T,
    double V,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const std::vector<double>& mu0,
    const ReactiveFlashConfig& cfg) const
{
    if (!(std::isfinite(T) && T > 0.0 && std::isfinite(V) && V > 0.0)) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTV: invalid T/V");
    }
    if (static_cast<int>(n0.size()) != eos_->numComponents()) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTV: n0 size mismatch");
    }
    if (static_cast<int>(mu0.size()) != eos_->numComponents()) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTV: mu0 size mismatch");
    }
    if (!(std::isfinite(cfg.standard_pressure) && cfg.standard_pressure > 0.0)) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTV: invalid standard_pressure");
    }

    const int nc = eos_->numComponents();
    double N0 = 0.0;
    for (double ni : n0) N0 += std::max(ni, cfg.min_moles);
    if (!(std::isfinite(N0) && N0 > 0.0)) {
        throw std::invalid_argument("ReactiveFlashSolver::solveTV: invalid n0");
    }
    std::vector<double> z0(static_cast<size_t>(nc), 0.0);
    for (int i = 0; i < nc; ++i) z0[static_cast<size_t>(i)] = std::max(n0[static_cast<size_t>(i)], cfg.min_moles) / N0;
    {
        double s = std::accumulate(z0.begin(), z0.end(), 0.0);
        if (!(std::isfinite(s) && s > 0.0)) s = 1.0;
        for (auto& zi : z0) zi = std::max(zi / s, Constants::MIN_MOLE_FRACTION);
    }

    // Initial pressure guess from a single-phase evaluation at rho_total = 1/V.
    const double rho_total = 1.0 / V;
    const double P_rho = eos_->pressure(T, rho_total, z0);
    const double P_ig = Constants::GAS_CONSTANT * T / V;
    const double P0 = Equilibrium::Internal::tvInitialPressureGuess(P_rho, P_ig, 1e5);

    auto root_finder = Factory::NumericsFactory::createBrent();
    Equilibrium::Optimization::ReactiveGibbsTPOptimizer inner(eos_, stability_);

    struct Eval {
        bool ok = false;
        double Vcalc = 0.0;
        Equilibrium::Optimization::EquilibriumState st;
    };

    auto evalAtP = [&](double P) -> Eval {
        Eval e;
        if (!(std::isfinite(P) && P > 0.0)) return e;
        const Core::Units::PropertyVariable pressure_tolerance(cfg.equilibrium.density_newton_tol, Core::Units::Unit::PA);

        Equilibrium::Optimization::ReactiveGibbsTPConfig ocfg;
        ocfg.standard_pressure = cfg.standard_pressure;
        ocfg.min_moles = cfg.min_moles;
        ocfg.equilibrium = cfg.equilibrium;
        ocfg.equilibrium.max_phases = std::clamp(cfg.equilibrium.max_phases, 1, 3);

        // Density continuity is handled after the inner solve (same approach as HelmholtzTVOptimizer).
        const auto requested_space = ocfg.equilibrium.phase_detection_space;
        ocfg.equilibrium.phase_detection_space = Config::PhaseDetection::TP;

        Equilibrium::Optimization::EquilibriumState st;
        try {
            st = inner.minimizeTPReactive(T, P, n0, system, mu0, ocfg);
        } catch (...) {
            return e;
        }
        if (st.phases.empty()) return e;

        double Vcalc = 0.0;
        const bool ok = Equilibrium::Internal::postProcessTVEvalAtTP(
            st.phases,
            *eos_,
            *stability_,
            T,
            P,
            requested_space,
            cfg.equilibrium.max_density_newton_iters,
            pressure_tolerance,
            cfg.equilibrium.beta_min,
            [&](PhaseType, double, double rho_cur) {
                return (requested_space == Config::PhaseDetection::TRho && std::isfinite(rho_cur) && rho_cur > 0.0)
                    ? rho_cur
                    : std::numeric_limits<double>::quiet_NaN();
            },
            [&](PhaseType t) { return rankPhase(t); },
            [&](const std::vector<PhaseState>& phases) { return volumeFromPhases(phases); },
            Vcalc);
        if (!ok) return e;

        e.ok = true;
        e.Vcalc = Vcalc;
        e.st = std::move(st);
        e.st.pressure = P;
        return e;
    };

    ReactiveFlashResult out;
    out.temperature = T;
    out.pressure = 0.0;
    out.n0 = n0;

    const auto rr = Equilibrium::Internal::solveTVPressureBrentOuter(
        *root_finder,
        evalAtP,
        V,
        P0,
        1e-12,
        cfg.equilibrium.tv_bracket_logp_step,
        cfg.equilibrium.tv_bracket_steps,
        cfg.equilibrium.tv_p_min,
        cfg.equilibrium.tv_p_max,
        cfg.equilibrium.optimizer.tolerance,
        "ReactiveFlashSolver::solveTV");
    if (!rr.converged) {
        out.converged = false;
        out.message = rr.message;
        return out;
    }
    const double P_star = rr.P_star;
    out.iterations = rr.iterations;
    const auto final_eval = rr.final_eval;

    out.pressure = P_star;
    out.z = final_eval.st.z;
    out.extents = final_eval.st.reaction_extents.value_or(std::vector<double>{});
    if (!out.extents.empty()) {
        out.n = system.molesFromExtents(n0, out.extents);
    }

    auto st = final_eval.st;
    Equilibrium::Internal::tvFinalizeStateForTVSolve(
        st,
        P_star,
        final_eval.Vcalc,
        V,
        rr.volume_abs_tolerance,
        "ReactiveHelmholtz(T,V) (Brent outer + ReactiveGibbs(T,P) inner)");

    out.converged = st.converged;
    out.message = st.message;
    out.flash = Equilibrium::Internal::flashResultFromState(st);
    return out;
}

ReactiveFlashResult ReactiveFlashSolver::solveTVFromDatabanks(
    double T,
    double V,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    const ReactiveFlashConfig& cfg,
    Data::Diagnostics* diag) const
{
    const std::vector<double> mu0 = Thermochemistry::mu0DB(
        cfg.mu0_model,
        databanks,
        mixture,
        T,
        cfg.mu0_reference_temperature,
        diag);
    return solveTV(T, V, n0, system, mu0, cfg);
}

ReactiveFlashResult ReactiveFlashSolver::solvePH(
    double P,
    double H,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const Mu0Provider& mu0AtT,
    const Core::Mixture& mixture,
    const ReactivePHFlashConfig& cfg) const
{
    if (!(std::isfinite(P) && P > 0.0)) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePH: invalid P");
    }
    if (!std::isfinite(H)) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePH: invalid H");
    }
    if (static_cast<int>(n0.size()) != eos_->numComponents()) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePH: n0 size mismatch");
    }
    Equilibrium::Internal::requireIdealGasCp(mixture, "ReactiveFlashSolver::solvePH");
    if (cfg.max_iterations <= 0) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePH: invalid max_iterations");
    }

    auto computeH = [&](double T, ReactiveFlashResult& out) -> bool {
        const std::vector<double> mu0 = mu0AtT(T);
        try {
            out = solveTP(T, P, n0, system, mu0, cfg.base);
            if (!out.converged || out.flash.phases.empty()) return false;

            double H_total = 0.0;
            for (const auto& phase : out.flash.phases) {
                const auto ig = Equilibrium::Contributions::idealGasProps(mixture, T, P, phase.x);
                const auto rr = Equilibrium::Internal::residualPropsForHS(*eos_, T, phase.density, phase.x, P, "ReactiveFlashSolver::solvePH");
                const double h_phase = ig.h_ig + rr.h_res;
                H_total += phase.fraction * h_phase;
            }
            if (!std::isfinite(H_total)) return false;
            out.flash.enthalpy = H_total;
            return true;
        } catch (...) {
            return false;
        }
    };

    double T_lo = cfg.temp_bracket_low;
    double T_hi = cfg.temp_bracket_high;
    if (!(std::isfinite(T_lo) && std::isfinite(T_hi) && T_lo > 0.0 && T_hi > T_lo)) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePH: invalid temperature bracket");
    }

    ReactiveFlashResult fr_lo, fr_hi;
    auto fval = [&](const ReactiveFlashResult& fr) -> double { return fr.flash.enthalpy.value() - H; };
    const bool bracketed = Equilibrium::Internal::bracketTemperatureByExpansion(
        T_lo, T_hi, fr_lo, fr_hi, computeH, fval, 8, 1e-9);
    if (!bracketed) {
        // Try a coarse log-spaced scan inside the bracket.
        const bool found = Equilibrium::Internal::bracketTemperatureByLogScan(
            T_lo, T_hi, fr_lo, fr_hi, computeH, fval, 50);
        if (!found) {
            ReactiveFlashResult fail;
            fail.converged = false;
            fail.temperature = 0.0;
            fail.pressure = P;
            fail.n0 = n0;
            fail.message = "ReactiveFlashSolver::solvePH: failed to bracket solution";
            fail.flash.pressure = P;
            fail.flash.method_used = "ReactiveFlash(PH)";
            return fail;
        }
    }

    const auto rr = Equilibrium::Internal::solveTemperatureBisection(
        fr_lo.temperature,
        fr_hi.temperature,
        fr_lo,
        computeH,
        fval,
        cfg.max_iterations,
        cfg.enthalpy_tolerance);
    auto fr_mid = rr.result;
    if (rr.converged) {
        fr_mid.converged = true;
        fr_mid.iterations = rr.iterations;
        fr_mid.message = "PH converged";
        fr_mid.flash.method_used = "ReactiveFlash(PH) outer + " + fr_mid.flash.method_used;
        return fr_mid;
    }

    fr_mid.converged = false;
    fr_mid.iterations = cfg.max_iterations;
    fr_mid.message = "PH did not converge";
    fr_mid.flash.method_used = "ReactiveFlash(PH) outer + " + fr_mid.flash.method_used;
    return fr_mid;
}

ReactiveFlashResult ReactiveFlashSolver::solvePHFromDatabanks(
    double P,
    double H,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    const ReactivePHFlashConfig& cfg,
    Data::Diagnostics* diag) const
{
    const auto mu0AtT = [&](double T) {
        return Thermochemistry::mu0DB(
            cfg.base.mu0_model,
            databanks,
            mixture,
            T,
            cfg.base.mu0_reference_temperature,
            diag);
    };
    return solvePH(P, H, n0, system, mu0AtT, mixture, cfg);
}

ReactiveFlashResult ReactiveFlashSolver::solvePS(
    double P,
    double S,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const Mu0Provider& mu0AtT,
    const Core::Mixture& mixture,
    const ReactivePSFlashConfig& cfg) const
{
    if (!(std::isfinite(P) && P > 0.0)) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePS: invalid P");
    }
    if (!std::isfinite(S)) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePS: invalid S");
    }
    if (static_cast<int>(n0.size()) != eos_->numComponents()) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePS: n0 size mismatch");
    }
    Equilibrium::Internal::requireIdealGasCp(mixture, "ReactiveFlashSolver::solvePS");
    if (cfg.max_iterations <= 0) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePS: invalid max_iterations");
    }

    auto computeS = [&](double T, ReactiveFlashResult& out) -> bool {
        const std::vector<double> mu0 = mu0AtT(T);
        try {
            out = solveTP(T, P, n0, system, mu0, cfg.base);
            if (!out.converged || out.flash.phases.empty()) return false;

            double S_total = 0.0;
            for (const auto& phase : out.flash.phases) {
                const auto ig = Equilibrium::Contributions::idealGasProps(mixture, T, P, phase.x);
                const auto rr = Equilibrium::Internal::residualPropsForHS(*eos_, T, phase.density, phase.x, P, "ReactiveFlashSolver::solvePS");
                const double s_phase = ig.s_ig + rr.s_res;
                S_total += phase.fraction * s_phase;
            }
            if (!std::isfinite(S_total)) return false;
            out.flash.entropy = S_total;
            return true;
        } catch (...) {
            return false;
        }
    };

    double T_lo = cfg.temp_bracket_low;
    double T_hi = cfg.temp_bracket_high;
    if (!(std::isfinite(T_lo) && std::isfinite(T_hi) && T_lo > 0.0 && T_hi > T_lo)) {
        throw std::invalid_argument("ReactiveFlashSolver::solvePS: invalid temperature bracket");
    }

    ReactiveFlashResult fr_lo, fr_hi;
    auto fval = [&](const ReactiveFlashResult& fr) -> double { return fr.flash.entropy.value() - S; };
    const bool bracketed = Equilibrium::Internal::bracketTemperatureByExpansion(
        T_lo, T_hi, fr_lo, fr_hi, computeS, fval, 8, 1e-9);
    if (!bracketed) {
        // Try a coarse log-spaced scan inside the bracket.
        const bool found = Equilibrium::Internal::bracketTemperatureByLogScan(
            T_lo, T_hi, fr_lo, fr_hi, computeS, fval, 50);
        if (!found) {
            ReactiveFlashResult fail;
            fail.converged = false;
            fail.temperature = 0.0;
            fail.pressure = P;
            fail.n0 = n0;
            fail.message = "ReactiveFlashSolver::solvePS: failed to bracket solution";
            fail.flash.pressure = P;
            fail.flash.method_used = "ReactiveFlash(PS)";
            return fail;
        }
    }

    const auto rr = Equilibrium::Internal::solveTemperatureBisection(
        fr_lo.temperature,
        fr_hi.temperature,
        fr_lo,
        computeS,
        fval,
        cfg.max_iterations,
        cfg.entropy_tolerance);
    auto fr_mid = rr.result;
    if (rr.converged) {
        fr_mid.converged = true;
        fr_mid.iterations = rr.iterations;
        fr_mid.message = "PS converged";
        fr_mid.flash.method_used = "ReactiveFlash(PS) outer + " + fr_mid.flash.method_used;
        return fr_mid;
    }

    fr_mid.converged = false;
    fr_mid.iterations = cfg.max_iterations;
    fr_mid.message = "PS did not converge";
    fr_mid.flash.method_used = "ReactiveFlash(PS) outer + " + fr_mid.flash.method_used;
    return fr_mid;
}

ReactiveFlashResult ReactiveFlashSolver::solvePSFromDatabanks(
    double P,
    double S,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    const ReactivePSFlashConfig& cfg,
    Data::Diagnostics* diag) const
{
    const auto mu0AtT = [&](double T) {
        return Thermochemistry::mu0DB(
            cfg.base.mu0_model,
            databanks,
            mixture,
            T,
            cfg.base.mu0_reference_temperature,
            diag);
    };
    return solvePS(P, S, n0, system, mu0AtT, mixture, cfg);
}

} // namespace Reactions
} // namespace DMThermo
