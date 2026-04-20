/**
 * @file gibbs_optimizer_flash.cpp
 * @brief Implementation of GibbsOptimizerFlashSolver.
 */

#include "thermo/equilibrium/gibbs_optimizer_flash.h"
#include "thermo/config/equilibrium_optimizer_config.h"
#include "thermo/equilibrium/internal/flash_config_support.h"
#include "thermo/equilibrium/internal/flash_result_support.h"
#include "thermo/equilibrium/internal/hs_support.h"
#include "thermo/equilibrium/internal/k_values.h"
#include "thermo/equilibrium/internal/temperature_solve.h"
#include "thermo/equilibrium/thermo_contributions.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

GibbsOptimizerFlashSolver::GibbsOptimizerFlashSolver(
    EOSPtr eos,
    Stability::StabilityAnalyzerPtr stability)
	    : optimizer_(std::move(eos), std::move(stability))
	{}

FlashResult GibbsOptimizerFlashSolver::solvePT(
    double T,
    double P,
    const std::vector<double>& z,
    const Config::FlashConfig& config) const
{
    auto ocfg = Internal::equilibriumOptimizerConfigFrom(config);

	    auto st = optimizer_.minimizeTP(T, P, z, ocfg);
	    auto out = Internal::flashResultFromState(st);
	    out.stability_test_performed = config.perform_stability_test;
	    out.feed_was_stable = (out.num_phases == 1);
	    return out;
	}

FlashResult GibbsOptimizerFlashSolver::solveTV(
    double,
    double,
    const std::vector<double>&,
    const Config::TVFlashConfig&) const
{
    throw std::runtime_error("GibbsOptimizerFlashSolver: solveTV not implemented (use Helmholtz optimizer)");
}

FlashResult GibbsOptimizerFlashSolver::solvePH(
    double P,
    double H,
    const std::vector<double>& z,
    const Config::PHFlashConfig& config) const
{
    const auto& mix = Internal::requireCoreMixture(*optimizer_.eos(), "GibbsOptimizerFlashSolver::solvePH");
    Internal::requireIdealGasCp(mix, "GibbsOptimizerFlashSolver::solvePH");

    auto computeH = [&](double T, FlashResult& fr_out) -> bool {
        try {
            fr_out = solvePT(T, P, z, config);
            if (!fr_out.converged || fr_out.phases.empty()) return false;

            double H_total = 0.0;
            for (const auto& phase : fr_out.phases) {
                const auto ig = Contributions::idealGasProps(mix, T, P, phase.x);
                const auto rr = Internal::residualPropsForHS(*optimizer_.eos(), T, phase.density, phase.x, P, "GibbsOptimizerFlashSolver::solvePH");
                const double h_phase = ig.h_ig + rr.h_res;
                H_total += phase.fraction * h_phase;
            }
            fr_out.enthalpy = H_total;
            return true;
        } catch (...) {
            return false;
        }
    };

    double T_lo = config.temp_bracket_low;
    double T_hi = config.temp_bracket_high;
    if (T_lo <= 0.0 || T_hi <= T_lo) {
        throw std::invalid_argument("GibbsOptimizerFlashSolver: invalid PH temperature bracket");
    }

    FlashResult fr_lo, fr_hi;
    auto fval = [&](const FlashResult& fr) -> double { return fr.enthalpy.value() - H; };
    const bool bracketed = Internal::bracketTemperatureByExpansion(
        T_lo, T_hi, fr_lo, fr_hi, computeH, fval, 8, 1e-6);
    if (!bracketed) {
        FlashResult fail;
        fail.converged = false;
        fail.pressure = P;
        fail.temperature = 0.0;
        fail.z = z;
        fail.method_used = algorithmName() + " (PH)";
        fail.message = "Failed to bracket PH solution";
        return fail;
    }

    const auto rr = Internal::solveTemperatureBisection(
        T_lo, T_hi, fr_lo, computeH, fval, config.max_iterations, config.enthalpy_tolerance);
    auto fr_mid = rr.result;
    if (rr.converged) {
        fr_mid.method_used = algorithmName() + " (PH)";
        fr_mid.message = "PH converged";
        return fr_mid;
    }

    fr_mid.converged = false;
    fr_mid.method_used = algorithmName() + " (PH)";
    fr_mid.message = "PH did not converge";
    return fr_mid;
}

FlashResult GibbsOptimizerFlashSolver::solvePS(
    double P,
    double S,
    const std::vector<double>& z,
    const Config::PSFlashConfig& config) const
{
    const auto& mix = Internal::requireCoreMixture(*optimizer_.eos(), "GibbsOptimizerFlashSolver::solvePS");
    Internal::requireIdealGasCp(mix, "GibbsOptimizerFlashSolver::solvePS");

    auto computeS = [&](double T, FlashResult& fr_out) -> bool {
        try {
            fr_out = solvePT(T, P, z, config);
            if (!fr_out.converged || fr_out.phases.empty()) return false;

            double S_total = 0.0;
            for (const auto& phase : fr_out.phases) {
                const auto ig = Contributions::idealGasProps(mix, T, P, phase.x);
                const auto rr = Internal::residualPropsForHS(*optimizer_.eos(), T, phase.density, phase.x, P, "GibbsOptimizerFlashSolver::solvePS");
                const double s_phase = ig.s_ig + rr.s_res;
                S_total += phase.fraction * s_phase;
            }
            fr_out.entropy = S_total;
            return true;
        } catch (...) {
            return false;
        }
    };

    double T_lo = config.temp_bracket_low;
    double T_hi = config.temp_bracket_high;
    if (T_lo <= 0.0 || T_hi <= T_lo) {
        throw std::invalid_argument("GibbsOptimizerFlashSolver: invalid PS temperature bracket");
    }

    FlashResult fr_lo, fr_hi;
    auto fval = [&](const FlashResult& fr) -> double { return fr.entropy.value() - S; };
    const bool bracketed = Internal::bracketTemperatureByExpansion(
        T_lo, T_hi, fr_lo, fr_hi, computeS, fval, 8, 1e-6);
    if (!bracketed) {
        FlashResult fail;
        fail.converged = false;
        fail.pressure = P;
        fail.temperature = 0.0;
        fail.z = z;
        fail.method_used = algorithmName() + " (PS)";
        fail.message = "Failed to bracket PS solution";
        return fail;
    }

    const auto rr = Internal::solveTemperatureBisection(
        T_lo, T_hi, fr_lo, computeS, fval, config.max_iterations, config.entropy_tolerance);
    auto fr_mid = rr.result;
    if (rr.converged) {
        fr_mid.method_used = algorithmName() + " (PS)";
        fr_mid.message = "PS converged";
        return fr_mid;
    }

    fr_mid.converged = false;
    fr_mid.method_used = algorithmName() + " (PS)";
    fr_mid.message = "PS did not converge";
    return fr_mid;
}

FlashResult GibbsOptimizerFlashSolver::solveVLLE(
    double T,
    double P,
    const std::vector<double>& z,
    const Config::FlashConfig& config) const
{
    Config::FlashConfig cfg = config;
    cfg.max_phases = std::max(cfg.max_phases, 3);
    return solvePT(T, P, z, cfg);
}

std::vector<double> GibbsOptimizerFlashSolver::estimateKValues(double T, double P) const {
    return Internal::wilsonK(*optimizer_.eos(), T, P);
}

bool GibbsOptimizerFlashSolver::isValidResult(const FlashResult& result) const {
    if (!result.converged) return false;
    if (result.temperature <= 0.0 || result.pressure <= 0.0) return false;
    if (result.num_phases < 1 || result.num_phases > 3) return false;
    if (result.phases.empty()) return false;
    return true;
}

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo
