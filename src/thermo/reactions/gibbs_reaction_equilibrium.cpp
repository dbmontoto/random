/**
 * @file gibbs_reaction_equilibrium.cpp
 * @brief Implementation of single-phase TP reaction equilibrium solver.
 */

#include "thermo/reactions/gibbs_reaction_equilibrium.h"
#include "thermo/reactions/thermochemistry.h"
#include "thermo/factory/numerics_factory.h"
#include "thermo/core/constants.h"
#include "thermo/equilibrium/optimization/transforms.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace DMThermo {
namespace Reactions {

namespace {

double safeLog(double x, double min_value) {
    return std::log(std::max(x, min_value));
}

double rhoAtTP(
    const EOS& eos,
    const Equilibrium::Stability::StabilityAnalyzerPtr& stability,
    double T,
    double P,
    const std::vector<double>& x,
    PhaseType phase)
{
    if (stability) {
        if (phase == PhaseType::Liquid || phase == PhaseType::Liquid1 || phase == PhaseType::Liquid2) {
            return stability->findLiquidDensity(T, P, x);
        }
        if (phase == PhaseType::Vapor) {
            return stability->findVaporDensity(T, P, x);
        }
        const auto roots = stability->findDensityRoots(T, P, x);
        if (roots.stable_root_index.has_value()) {
            const int k = roots.stable_root_index.value();
            if (k >= 0 && k < static_cast<int>(roots.densities.size())) {
                return roots.densities[static_cast<std::size_t>(k)];
            }
        }
        if (!roots.densities.empty()) return roots.densities.front();
    }

    // If no stability analyzer is available, let EOS implementations provide TP roots when possible.
    try {
        const auto roots = eos.densityRootsTP(T, P, x);
        if (!roots.empty()) {
            if (phase == PhaseType::Liquid || phase == PhaseType::Liquid1 || phase == PhaseType::Liquid2) {
                return *std::max_element(roots.begin(), roots.end());
            }
            if (phase == PhaseType::Vapor) {
                return *std::min_element(roots.begin(), roots.end());
            }
            return roots.front();
        }
    } catch (...) {
        // Fall through to ideal-gas density.
    }

    // Fallback: treat EOS as ideal-gas-like for density initialization.
    const double R = Constants::GAS_CONSTANT;
    return (R > 0.0 && T > 0.0) ? (P / (R * T)) : 0.0;
}

} // namespace

ReactionEquilibriumResult GibbsReactionEquilibriumSolver::solveTPFromDatabanks(
    double T,
    double P,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    PhaseType phase,
    const ReactionEquilibriumConfig& cfg,
    Data::Diagnostics* diag) const
{
    const std::vector<double> mu0 = Thermochemistry::mu0DB(
        cfg.mu0_model,
        databanks,
        mixture,
        T,
        cfg.mu0_reference_temperature,
        diag);
    return solveTP(T, P, n0, system, mu0, phase, cfg);
}

ReactionEquilibriumResult GibbsReactionEquilibriumSolver::solveTP(
    double T,
    double P,
    const std::vector<double>& n0,
    const ReactionSystem& system,
    const std::vector<double>& mu0,
    PhaseType phase,
    const ReactionEquilibriumConfig& cfg) const
{
    if (!eos_) {
        throw std::invalid_argument("GibbsReactionEquilibriumSolver::solveTP: EOS is null");
    }
    if (!(T > 0.0) || !(P > 0.0) || !std::isfinite(T) || !std::isfinite(P)) {
        throw std::invalid_argument("GibbsReactionEquilibriumSolver::solveTP: invalid T/P");
    }
    if (system.numComponents() != static_cast<int>(n0.size())) {
        throw std::invalid_argument("GibbsReactionEquilibriumSolver::solveTP: n0 size mismatch");
    }
    if (system.numComponents() != static_cast<int>(mu0.size())) {
        throw std::invalid_argument("GibbsReactionEquilibriumSolver::solveTP: mu0 size mismatch");
    }
    if (system.numComponents() != eos_->numComponents()) {
        throw std::invalid_argument("GibbsReactionEquilibriumSolver::solveTP: EOS component count mismatch");
    }
    for (double ni : n0) {
        if (!(std::isfinite(ni) && ni >= 0.0)) {
            throw std::invalid_argument("GibbsReactionEquilibriumSolver::solveTP: invalid n0 (must be finite, >=0)");
        }
    }

    const double R = Constants::GAS_CONSTANT;
    const double invRT = 1.0 / (R * T);
    if (!std::isfinite(invRT) || invRT <= 0.0) {
        throw std::runtime_error("GibbsReactionEquilibriumSolver::solveTP: invalid 1/(R*T)");
    }

    const int nr = system.numReactions();
    std::vector<double> u0(static_cast<std::size_t>(nr), 0.0);

    // Conservative per-reaction extent bounds, used to keep the optimizer in a feasible region.
    // For multi-reaction systems this is only a box approximation; keep infeasibility penalties as a backstop.
    std::vector<ExtentBounds> bounds;
    bounds.reserve(static_cast<std::size_t>(nr));
    for (int r = 0; r < nr; ++r) {
        bounds.push_back(system.extentBounds(r, n0, cfg.min_moles));
        double xi0 = 0.0;
        if (xi0 < bounds.back().min || xi0 > bounds.back().max) {
            xi0 = 0.5 * (bounds.back().min + bounds.back().max);
        }
        xi0 = std::clamp(xi0, bounds.back().min, bounds.back().max);
        u0[static_cast<std::size_t>(r)] = Equilibrium::Optimization::unboundedFromBounded(
            xi0,
            bounds.back().min,
            bounds.back().max);
    }

    const double P0 = cfg.standard_pressure;
    if (!(P0 > 0.0) || !std::isfinite(P0)) {
        throw std::invalid_argument("GibbsReactionEquilibriumSolver::solveTP: invalid standard_pressure");
    }

    auto objective = [&](const std::vector<double>& u) -> double {
        if (static_cast<int>(u.size()) != nr) return std::numeric_limits<double>::infinity();

        std::vector<double> extents(static_cast<std::size_t>(nr), 0.0);
        for (int r = 0; r < nr; ++r) {
            extents[static_cast<std::size_t>(r)] = Equilibrium::Optimization::boundedFromUnbounded(
                u[static_cast<std::size_t>(r)],
                bounds[static_cast<std::size_t>(r)].min,
                bounds[static_cast<std::size_t>(r)].max);
        }

        const std::vector<double> n = system.molesFromExtents(n0, extents);
        double N = 0.0;
        double infeas = 0.0;
        for (double ni : n) {
            if (!std::isfinite(ni)) return std::numeric_limits<double>::infinity();
            if (ni < cfg.min_moles) {
                const double d = cfg.min_moles - ni;
                infeas += d * d;
            }
            N += std::max(ni, cfg.min_moles);
        }
        if (!(N > 0.0) || !std::isfinite(N)) return std::numeric_limits<double>::infinity();

        std::vector<double> x(n.size(), 0.0);
        for (std::size_t i = 0; i < n.size(); ++i) {
            x[i] = std::max(n[i], cfg.min_moles) / N;
        }

        const double rho = rhoAtTP(*eos_, stability_, T, P, x, phase);
        if (!(std::isfinite(rho) && rho > 0.0)) {
            return cfg.infeasibility_penalty * (1.0 + infeas);
        }

        std::vector<double> phi;
        try {
            phi = eos_->fugacityCoefficients(T, rho, x);
        } catch (...) {
            return cfg.infeasibility_penalty * (1.0 + infeas);
        }
        if (phi.size() != x.size()) return cfg.infeasibility_penalty * (1.0 + infeas);

        const double lnP = safeLog(P / P0, cfg.min_moles);

        double g_over_rt = 0.0;
        for (std::size_t i = 0; i < n.size(); ++i) {
            const double ni = std::max(n[i], cfg.min_moles);
            const double ln_x = safeLog(x[i], cfg.min_moles);
            const double ln_phi = safeLog(phi[i], cfg.min_moles);
            const double mu0_over_rt = mu0[i] * invRT;
            g_over_rt += ni * (mu0_over_rt + ln_x + ln_phi + lnP);
        }

        if (!std::isfinite(g_over_rt)) return std::numeric_limits<double>::infinity();
        return g_over_rt + cfg.infeasibility_penalty * infeas;
    };

    auto optimizer = Factory::NumericsFactory::createBFGS();
    const auto opt = optimizer->minimize(objective, u0, cfg.optimizer);

    ReactionEquilibriumResult out;
    out.temperature = T;
    out.pressure = P;
    out.n0 = n0;
    out.extents.assign(static_cast<std::size_t>(nr), 0.0);
    for (int r = 0; r < nr; ++r) {
        out.extents[static_cast<std::size_t>(r)] = Equilibrium::Optimization::boundedFromUnbounded(
            opt.x_optimal[static_cast<std::size_t>(r)],
            bounds[static_cast<std::size_t>(r)].min,
            bounds[static_cast<std::size_t>(r)].max);
    }
    out.converged = opt.converged;
    out.iterations = opt.iterations;
    out.message = opt.message;
    out.objective_g_over_rt = opt.f_optimal;

    out.n = system.molesFromExtents(n0, out.extents);
    double N = 0.0;
    for (double ni : out.n) N += std::max(ni, cfg.min_moles);
    out.x.assign(out.n.size(), 0.0);
    if (N > 0.0 && std::isfinite(N)) {
        for (std::size_t i = 0; i < out.n.size(); ++i) {
            out.x[i] = std::max(out.n[i], cfg.min_moles) / N;
        }
    }

    // Populate a single-phase EOS state (best-effort; throws on EOS errors).
    if (!out.x.empty()) {
        const double rho = rhoAtTP(*eos_, stability_, T, P, out.x, phase);
        if (std::isfinite(rho) && rho > 0.0) {
            PhaseState ps;
            ps.type = phase;
            ps.density = rho;
            ps.compressibility = eos_->compressibility(T, rho, out.x);
            ps.x = out.x;
            ps.phi = eos_->fugacityCoefficients(T, rho, out.x);
            ps.fraction = 1.0;
            out.phases = {ps};
        }
    }

    return out;
}

} // namespace Reactions
} // namespace DMThermo
