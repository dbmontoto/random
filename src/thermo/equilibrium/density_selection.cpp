/**
 * @file density_selection.cpp
 * @brief Implementation of density selection/refinement at (T,P,x).
 */

#include "thermo/equilibrium/density_selection.h"
#include "thermo/core/constants.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/istability.h"
#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Detail {

PhiEval evalAtTP(
    const EOS& eos,
    const Stability::IStabilityAnalyzer& stability,
    double T,
    double P,
    const std::vector<double>& x,
    PhaseType phase_hint,
    Config::PhaseDetection space,
    double rho_guess,
    int max_newton_iters,
    const Core::Units::PropertyVariable& pressure_tolerance)
{
    PhiEval out;
    out.phi.assign(x.size(), 1.0);
    out.lnphi.assign(x.size(), 0.0);

    const double newton_tol = pressure_tolerance.as(Core::Units::Unit::PA);
    if (!(std::isfinite(newton_tol) && newton_tol > 0.0)) {
        throw std::invalid_argument("evalAtTP: pressure_tolerance must be finite and > 0");
    }

    auto refineRhoAtP = [&](double rho0) -> std::optional<double> {
        if (!(std::isfinite(rho0) && rho0 > 0.0)) return std::nullopt;
        double rho = std::clamp(rho0, Constants::MIN_DENSITY, Constants::MAX_DENSITY);
        for (int it = 0; it < max_newton_iters; ++it) {
            const double Pcalc = eos.pressure(T, rho, x);
            if (!std::isfinite(Pcalc)) return std::nullopt;
            const double f = Pcalc - P;
            if (std::abs(f) < newton_tol) return rho;

            const double dP = eos.dPdrho(T, rho, x);
            if (!(std::isfinite(dP) && std::abs(dP) > 0.0)) return std::nullopt;

            double step = f / dP;
            if (!std::isfinite(step)) return std::nullopt;
            step = std::clamp(step, -0.5 * rho, 0.5 * rho);
            const double rho_new = rho - step;
            if (!(std::isfinite(rho_new) && rho_new > 0.0)) return std::nullopt;
            rho = std::clamp(rho_new, Constants::MIN_DENSITY, Constants::MAX_DENSITY);
        }
        const double Pcalc = eos.pressure(T, rho, x);
        if (std::isfinite(Pcalc) && std::abs(Pcalc - P) < 100.0 * newton_tol) {
            return rho;
        }
        return std::nullopt;
    };

    const double pressure_tol = std::max(newton_tol, 1e-6 * std::max(1.0, std::abs(P)));

    double rho = 0.0;

    if (space == Config::PhaseDetection::TRho) {
        const double rho_ig = std::max(P / (Constants::GAS_CONSTANT * T), Constants::MIN_DENSITY);
        double seed = rho_guess;
        if (!(std::isfinite(seed) && seed > 0.0)) {
            seed = (phase_hint == PhaseType::Vapor) ? rho_ig : std::min(Constants::MAX_DENSITY, 500.0 * rho_ig);
        }
        if (auto r = refineRhoAtP(seed)) {
            rho = *r;
        }
    }

    if (!(std::isfinite(rho) && rho > 0.0)) {
        const auto roots = stability.findDensityRoots(T, P, x);
        if (!roots.densities.empty()) {
            if (phase_hint == PhaseType::Vapor) {
                rho = *std::min_element(roots.densities.begin(), roots.densities.end());
            } else if (phase_hint == PhaseType::Liquid || phase_hint == PhaseType::Liquid1 || phase_hint == PhaseType::Liquid2) {
                rho = *std::max_element(roots.densities.begin(), roots.densities.end());
            } else if (roots.stable_root_index.has_value()) {
                const int k = roots.stable_root_index.value();
                if (k >= 0 && k < static_cast<int>(roots.densities.size())) rho = roots.densities[k];
            } else {
                rho = roots.densities.front();
            }
        }
        if (auto r = refineRhoAtP(rho)) {
            rho = *r;
        }
    }

    if (!(std::isfinite(rho) && rho > 0.0)) {
        throw std::runtime_error("evalAtTP: unable to find TP-consistent density root");
    }

    const double Pcalc = eos.pressure(T, rho, x);
    if (!(std::isfinite(Pcalc) && std::abs(Pcalc - P) <= pressure_tol)) {
        std::ostringstream oss;
        oss << "evalAtTP: selected density is off-pressure (|Pcalc-P|=" << std::abs(Pcalc - P)
            << " Pa, tol=" << pressure_tol << " Pa)";
        throw std::runtime_error(oss.str());
    }

    out.rho = rho;
    out.phi = eos.fugacityCoefficients(T, rho, x);
    if (out.phi.size() != x.size()) {
        throw std::runtime_error("evalAtTP: EOS returned wrong phi size");
    }
    for (size_t i = 0; i < x.size(); ++i) {
        out.lnphi[i] = std::log(std::max(out.phi[i], 1e-300));
    }
    out.phase = stability.classifyPhase(T, rho, x);
    return out;
}

} // namespace Detail
} // namespace Equilibrium
} // namespace DMThermo
