/**
 * @file tp_departures.cpp
 * @brief TP departure properties per density root.
 */

#include "thermo/equilibrium/tp_departures.h"

#include "thermo/core/constants.h"
#include "thermo/equilibrium/tpd_stability.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Departures {

namespace {

void requireValidComposition(const EOS& eos, const std::vector<double>& x) {
    const int nc = eos.numComponents();
    if (static_cast<int>(x.size()) != nc) {
        throw std::invalid_argument("departureRootsTP: composition size mismatch");
    }
}

double gDepFromLnPhi(double T, const std::vector<double>& x, const std::vector<double>& ln_phi) {
    if (ln_phi.size() != x.size()) {
        throw std::runtime_error("departureRootsTP: ln_phi/composition size mismatch");
    }
    double g_over_rt = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        g_over_rt += x[i] * ln_phi[i];
    }
    return Constants::GAS_CONSTANT * T * g_over_rt;
}

bool refineDensityNewton(
    const EOS& eos,
    double T,
    double P_target,
    const std::vector<double>& x,
    double& rho,
    int max_iters,
    double tol_P)
{
    if (!(std::isfinite(rho) && rho > 0.0)) return false;
    rho = std::max(rho, Constants::MIN_DENSITY);

    for (int it = 0; it < max_iters; ++it) {
        const double Pcalc = eos.pressure(T, rho, x);
        if (!std::isfinite(Pcalc)) return false;

        const double f = Pcalc - P_target;
        if (!std::isfinite(f)) return false;
        if (std::abs(f) <= tol_P) return true;

        const double dP = eos.dPdrho(T, rho, x);
        if (!std::isfinite(dP) || std::abs(dP) < 1e-18) return false;

        double step = f / dP;
        if (!std::isfinite(step)) return false;
        step = std::clamp(step, -0.5 * rho, 0.5 * rho);

        const double rho_new = rho - step;
        if (!(std::isfinite(rho_new) && rho_new > 0.0)) return false;
        if (rho_new == rho) return false;
        rho = std::max(rho_new, Constants::MIN_DENSITY);
    }

    const double Pfinal = eos.pressure(T, rho, x);
    return std::isfinite(Pfinal) && std::abs(Pfinal - P_target) <= 10.0 * tol_P;
}

} // namespace

TPDepartureResult departureRootsTP(
    const Stability::IStabilityAnalyzer& stability,
    double T,
    double P,
    const std::vector<double>& x,
    const Config::DensityRootConfig& config)
{
    auto eos = stability.eos();
    if (!eos) {
        throw std::runtime_error("departureRootsTP: stability analyzer returned null EOS");
    }
    requireValidComposition(*eos, x);
    if (!(std::isfinite(T) && std::isfinite(P) && T > 0.0 && P > 0.0)) {
        throw std::invalid_argument("departureRootsTP: invalid T/P");
    }

    const auto roots = stability.findDensityRoots(T, P, x, config);

    TPDepartureResult out;
    out.temperature = T;
    out.pressure = P;
    out.x = x;

    const double R = Constants::GAS_CONSTANT;
    out.roots.reserve(roots.densities.size());
    std::vector<int> old_to_new(roots.densities.size(), -1);

    const double refine_tol_P = std::max(config.tolerance, 1e-6 * P);
    const int refine_max_iters = std::min(20, std::max(1, config.max_iterations));

    for (size_t i = 0; i < roots.densities.size(); ++i) {
        double rho = roots.densities[i];
        if (!(std::isfinite(rho) && rho > 0.0)) {
            continue;
        }

        TPDepartureRoot r;
        // Ensure rho is consistent with the requested TP state before computing departures.
        {
            double rho_refined = rho;
            if (refineDensityNewton(*eos, T, P, x, rho_refined, refine_max_iters, refine_tol_P)) {
                rho = rho_refined;
            }
        }
        r.density = rho;

        if (i < roots.phase_types.size()) {
            r.phase = roots.phase_types[i];
        }
        {
            // Prefer to recompute after refinement (if any); fall back to root-finder metadata when available.
            const double dPdrho = eos->dPdrho(T, rho, x);
            if (std::isfinite(dPdrho)) {
                r.mechanically_stable = (dPdrho > 0.0);
            } else if (i < roots.is_mechanically_stable.size()) {
                r.mechanically_stable = roots.is_mechanically_stable[i];
            }
        }

        const auto eos_res = eos->calculate(T, rho, x, DerivativeSpec{true, false, false, false});
        if (!eos_res.success || !eos_res.da_dT.has_value() || eos_res.ln_phi.size() != x.size()) {
            // Discard unphysical roots (e.g., Z<=B or v<=b for cubics) and any EOS that
            // cannot provide the derivatives required for h_dep/s_dep.
            continue;
        }

        // Explicit TP definition for departures: Z = P/(rho*R*T) using the target pressure.
        r.Z = P / (rho * R * T);
        r.ln_phi = eos_res.ln_phi;

        // Departure functions relative to ideal gas at the same (T,P).
        const double a = eos_res.a_residual;         // A_res/(RT), relative to ideal gas at same (T,rho)
        const double da_dT = eos_res.da_dT.value();  // d(A_res/RT)/dT at constant rho
        r.s_dep = -R * (a + T * da_dT) + R * std::log(std::max(r.Z, 1e-300));
        r.h_dep = -R * T * T * da_dT + R * T * (r.Z - 1.0);
        r.g_dep = gDepFromLnPhi(T, x, r.ln_phi);

        old_to_new[i] = static_cast<int>(out.roots.size());
        out.roots.push_back(std::move(r));
    }

    // Map the root-finder stable index onto the retained roots (if possible).
    if (roots.stable_root_index.has_value()) {
        const int k = roots.stable_root_index.value();
        if (k >= 0 && k < static_cast<int>(old_to_new.size()) && old_to_new[static_cast<size_t>(k)] >= 0) {
            out.stable_root_index = old_to_new[static_cast<size_t>(k)];
        }
    }

    return out;
}

TPDepartureResult departureRootsTP(
    const Stability::StabilityAnalyzerPtr& stability,
    double T,
    double P,
    const std::vector<double>& x,
    const Config::DensityRootConfig& config)
{
    if (!stability) {
        throw std::invalid_argument("departureRootsTP: stability analyzer is null");
    }
    return departureRootsTP(*stability, T, P, x, config);
}

TPDepartureResult departureRootsTP(
    EOSPtr eos,
    double T,
    double P,
    const std::vector<double>& x,
    const Config::DensityRootConfig& config)
{
    if (!eos) {
        throw std::invalid_argument("departureRootsTP: EOS is null");
    }
    Stability::TPDStabilityAnalyzer stability(std::move(eos));
    return departureRootsTP(stability, T, P, x, config);
}

std::optional<int> selectRootIndex(const TPDepartureResult& result, RootSelection selection) {
    const int n = static_cast<int>(result.roots.size());
    if (n <= 0) return std::nullopt;

    if (selection == RootSelection::Vapor) {
        bool any_mech = false;
        for (const auto& r : result.roots) {
            if (r.mechanically_stable) {
                any_mech = true;
                break;
            }
        }

        int best = -1;
        double best_rho = std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            const auto& r = result.roots[static_cast<size_t>(i)];
            if (any_mech && !r.mechanically_stable) continue;
            if (!std::isfinite(r.density)) continue;
            if (r.density < best_rho) {
                best_rho = r.density;
                best = i;
            }
        }
        if (best >= 0) return best;
        return std::nullopt;
    }

    if (selection == RootSelection::Liquid) {
        bool any_mech = false;
        for (const auto& r : result.roots) {
            if (r.mechanically_stable) {
                any_mech = true;
                break;
            }
        }

        int best = -1;
        double best_rho = 0.0;
        for (int i = 0; i < n; ++i) {
            const auto& r = result.roots[static_cast<size_t>(i)];
            if (any_mech && !r.mechanically_stable) continue;
            if (!std::isfinite(r.density)) continue;
            if (best < 0 || r.density > best_rho) {
                best_rho = r.density;
                best = i;
            }
        }
        if (best >= 0) return best;
        return std::nullopt;
    }

    // Stable single-phase root.
    if (result.stable_root_index.has_value()) {
        const int k = result.stable_root_index.value();
        if (k >= 0 && k < n) return k;
    }

    bool any_mech = false;
    for (const auto& r : result.roots) {
        if (r.mechanically_stable) {
            any_mech = true;
            break;
        }
    }

    int best = -1;
    double best_g = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
        const auto& r = result.roots[static_cast<size_t>(i)];
        if (any_mech && !r.mechanically_stable) continue;
        if (!std::isfinite(r.g_dep)) continue;
        if (r.g_dep < best_g) {
            best_g = r.g_dep;
            best = i;
        }
    }
    if (best >= 0) return best;
    return std::nullopt;
}

const TPDepartureRoot* selectRoot(const TPDepartureResult& result, RootSelection selection) {
    const auto idx = selectRootIndex(result, selection);
    if (!idx.has_value()) return nullptr;
    const int k = idx.value();
    if (k < 0 || k >= static_cast<int>(result.roots.size())) return nullptr;
    return &result.roots[static_cast<size_t>(k)];
}

} // namespace Departures
} // namespace Equilibrium
} // namespace DMThermo
