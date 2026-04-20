/**
 * @file tpd_stability.cpp
 * @brief Implementation of EOS-based TPD stability analyzer
 */

#include "thermo/equilibrium/tpd_stability.h"
#include "thermo/core/constants.h"
#include "thermo/core/units.h"
#include "thermo/equilibrium/internal/k_values.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Stability {
namespace {

bool refineDensityNewton(
    const EOS& eos,
    double T,
    double P_target,
    const std::vector<double>& x,
    double& rho,
    int max_iters,
    double pressure_tol,
    double rho_min,
    double rho_max)
{
    if (!(rho > 0.0)) return false;
    rho = std::clamp(rho, rho_min, rho_max);

    for (int it = 0; it < max_iters; ++it) {
        const double Pcalc = eos.pressure(T, rho, x);
        if (!std::isfinite(Pcalc)) return false;
        const double f = Pcalc - P_target;
        if (!std::isfinite(f)) return false;
        if (std::abs(f) <= pressure_tol) return true;

        const double dP = eos.dPdrho(T, rho, x);
        if (!std::isfinite(dP) || std::abs(dP) < 1e-18) return false;

        double step = -f / dP;
        if (!std::isfinite(step)) return false;

        // Basic damping / backtracking to keep rho in bounds and evaluations finite.
        double rho_new = rho + step;
        for (int bt = 0; bt < 12; ++bt) {
            if (std::isfinite(rho_new) && rho_new > rho_min && rho_new < rho_max) {
                const double Ptry = eos.pressure(T, rho_new, x);
                if (std::isfinite(Ptry)) {
                    break;
                }
            }
            step *= 0.5;
            rho_new = rho + step;
        }

        if (!(std::isfinite(rho_new) && rho_new > rho_min && rho_new < rho_max)) {
            return false;
        }
        if (rho_new == rho) return false;
        rho = rho_new;
    }

    const double Pfinal = eos.pressure(T, rho, x);
    if (std::isfinite(Pfinal) && std::abs(Pfinal - P_target) <= pressure_tol) return true;
    return false;
}

double pressureTolerance(double P_target, const Core::Units::PropertyVariable& configured_tolerance) {
    const double configured_tol = configured_tolerance.as(Core::Units::Unit::PA);
    return std::max(configured_tol, 1e-6 * std::max(1.0, std::abs(P_target)));
}

Config::DensityRootConfig densityRootConfigFromStability(const Config::StabilityConfig& config) {
    Config::DensityRootConfig roots_cfg;
    roots_cfg.tolerance = config.density_search_precision;
    roots_cfg.max_iterations = std::max(1, config.max_iterations);
    roots_cfg.use_spinodal_bounds = config.use_spinodal_filtering;
    roots_cfg.filter_metastable = config.use_spinodal_filtering;
    roots_cfg.max_roots = std::max(0, config.max_density_roots);
    return roots_cfg;
}

double tpdReferenceScale(const std::vector<double>& z, const std::vector<double>& lnphi_z) {
    if (z.size() != lnphi_z.size() || z.empty()) return 1.0;

    double g_ref = 0.0;
    for (size_t i = 0; i < z.size(); ++i) {
        const double zi = std::max(z[i], Constants::MIN_MOLE_FRACTION);
        g_ref += zi * (std::log(zi) + lnphi_z[i]);
    }
    return std::max(1.0, std::abs(g_ref));
}

bool tpdIndicatesInstability(double tpd, double tpd_scale, const Config::StabilityConfig& config) {
    if (!std::isfinite(tpd)) return false;
    const double metric = config.use_relative_tpd ? (tpd / std::max(1.0, tpd_scale)) : tpd;
    return metric < config.tpd_instability_threshold;
}

} // namespace

TPDStabilityAnalyzer::TPDStabilityAnalyzer(EOSPtr eos)
    : eos_(std::move(eos))
{
    if (!eos_) {
        throw std::invalid_argument("TPDStabilityAnalyzer: EOS pointer is null");
    }
}

std::vector<double> TPDStabilityAnalyzer::normalize(const std::vector<double>& x) {
    double sum = std::accumulate(x.begin(), x.end(), 0.0);
    if (sum <= 0.0) {
        throw std::invalid_argument("Cannot normalize composition with non-positive sum");
    }
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        y[i] = x[i] / sum;
    }
    return y;
}

std::vector<double> TPDStabilityAnalyzer::safeLog(const std::vector<double>& x) {
    std::vector<double> out(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        out[i] = std::log(std::max(x[i], Constants::MIN_MOLE_FRACTION));
    }
    return out;
}

double TPDStabilityAnalyzer::tpdFromLnPhi(
    const std::vector<double>& z,
    const std::vector<double>& lnphi_z,
    const std::vector<double>& w,
    const std::vector<double>& lnphi_w)
{
    const auto lz = safeLog(z);
    const auto lw = safeLog(w);
    double tpd = 0.0;
    for (size_t i = 0; i < w.size(); ++i) {
        tpd += w[i] * ((lw[i] + lnphi_w[i]) - (lz[i] + lnphi_z[i]));
    }
    return tpd;
}

std::vector<double> TPDStabilityAnalyzer::wilsonK(const EOS& eos, double T, double P) {
    return Internal::wilsonK(eos, T, P);
}

TPDStabilityAnalyzer::LnPhiAtTP TPDStabilityAnalyzer::lnPhiAtTP(
    double T,
    double P,
    const std::vector<double>& x,
    PhaseType phase_hint,
    const Config::StabilityConfig& config,
    double rho_guess) const
{
    LnPhiAtTP out;
    out.ln_phi.assign(x.size(), 0.0);

    const int nc = eos_->numComponents();
    if (static_cast<int>(x.size()) != nc) {
        throw std::invalid_argument("lnPhiAtTP: composition size mismatch");
    }

    const double rho_ig = std::max(P / (Constants::GAS_CONSTANT * T), Constants::MIN_DENSITY);
    const double rho_min = std::max(Constants::MIN_DENSITY, 1e-6 * rho_ig);
    const double rho_max = std::max(50000.0, rho_ig * 10.0);
    const Core::Units::PropertyVariable configured_tolerance(config.density_newton_tol, Core::Units::Unit::PA);
    const double pressure_tol = pressureTolerance(P, configured_tolerance);

    double rho = 0.0;
    if (config.phase_detection_space == Config::PhaseDetection::TRho && rho_guess > 0.0) {
        rho = rho_guess;
        if (!refineDensityNewton(*eos_, T, P, x, rho, config.max_density_newton_iters, pressure_tol, rho_min, rho_max)) {
            rho = 0.0;
        }
    }

    if (rho <= 0.0) {
        auto roots_cfg = densityRootConfigFromStability(config);
        DensityRootResult roots = findDensityRoots(T, P, x, roots_cfg);
        if (roots.num_roots > 0) {
            // Choose root based on hint.
            if (phase_hint == PhaseType::Vapor) {
                rho = *std::min_element(roots.densities.begin(), roots.densities.end());
            } else if (phase_hint == PhaseType::Liquid || phase_hint == PhaseType::Liquid1) {
                rho = *std::max_element(roots.densities.begin(), roots.densities.end());
            } else if (roots.stable_root_index.has_value()) {
                rho = roots.densities[roots.stable_root_index.value()];
            } else {
                rho = roots.densities.front();
            }
        }

        if (rho <= 0.0) {
            throw std::runtime_error("lnPhiAtTP: unable to locate TP-consistent density root");
        }

        double rho_refined = rho;
        if (refineDensityNewton(*eos_, T, P, x, rho_refined, config.max_density_newton_iters, pressure_tol, rho_min, rho_max)) {
            rho = rho_refined;
        }
    }

    const double Pcalc = eos_->pressure(T, rho, x);
    if (!(std::isfinite(Pcalc) && std::abs(Pcalc - P) <= pressure_tol)) {
        std::ostringstream oss;
        oss << "lnPhiAtTP: selected density is off-pressure (|Pcalc-P|=" << std::abs(Pcalc - P)
            << " Pa, tol=" << pressure_tol << " Pa)";
        throw std::runtime_error(oss.str());
    }

    const auto r = eos_->calculate(T, rho, x, DerivativeSpec::none());
    if (!r.success) {
        throw std::runtime_error("lnPhiAtTP: EOS calculate() failed: " + r.error_message);
    }
    if (r.ln_phi.size() != static_cast<size_t>(nc)) {
        throw std::runtime_error("lnPhiAtTP: EOS returned invalid ln_phi size");
    }

    out.rho = rho;
    out.Z = r.compressibility;
    out.ln_phi = r.ln_phi;
    return out;
}

DensityRootResult TPDStabilityAnalyzer::findDensityRoots(
    double T,
    double P,
    const std::vector<double>& x,
    const Config::DensityRootConfig& config) const
{
    DensityRootResult result;

    const int nc = eos_->numComponents();
    if (static_cast<int>(x.size()) != nc) {
        throw std::invalid_argument("findDensityRoots: composition size mismatch");
    }

    auto f = [&](double rho) -> double {
        if (rho <= 0.0) return -P;
        const double Pcalc = eos_->pressure(T, rho, x);
        if (!std::isfinite(Pcalc)) {
            return std::numeric_limits<double>::infinity();
        }
        return Pcalc - P;
    };

    const double rho_ig = std::max(P / (Constants::GAS_CONSTANT * T), config.min_density);
    const double rho_min = std::max(config.min_density, 1e-6 * rho_ig);
    const double rho_max = std::max(config.max_density, rho_ig * 10.0);

    std::vector<double> roots;
    roots.reserve(8);

    // Fast path: allow EOS implementations to provide TP roots directly (e.g., cubic EOS).
    try {
        roots = eos_->densityRootsTP(T, P, x);
    } catch (...) {
        roots.clear();
    }

    if (!roots.empty()) {
        roots.erase(std::remove_if(roots.begin(), roots.end(), [&](double rho) {
            return !(std::isfinite(rho) && rho >= rho_min && rho <= rho_max);
        }), roots.end());
    }

    auto scanInterval = [&](double a_rho, double b_rho, int n_scan, std::vector<double>& roots) {
        if (a_rho <= 0.0 || b_rho <= a_rho) return;

        std::vector<double> grid(n_scan);
        const double log_min = std::log(a_rho);
        const double log_max = std::log(b_rho);
        for (int i = 0; i < n_scan; ++i) {
            const double t = static_cast<double>(i) / (n_scan - 1);
            grid[i] = std::exp(log_min + t * (log_max - log_min));
        }

        double prev_rho = grid[0];
        double prev_f = f(prev_rho);
        for (int i = 1; i < n_scan; ++i) {
            const double cur_rho = grid[i];
            const double cur_f = f(cur_rho);

            if (!std::isfinite(prev_f) || !std::isfinite(cur_f)) {
                prev_rho = cur_rho;
                prev_f = cur_f;
                continue;
            }

            if (prev_f == 0.0) {
                roots.push_back(prev_rho);
            } else if (cur_f == 0.0) {
                roots.push_back(cur_rho);
            } else if ((prev_f < 0.0 && cur_f > 0.0) || (prev_f > 0.0 && cur_f < 0.0)) {
                double a = prev_rho, b = cur_rho;
                double fa = prev_f;

                for (int it = 0; it < config.max_iterations; ++it) {
                    const double m = 0.5 * (a + b);
                    const double fm = f(m);
                    if (!std::isfinite(fm)) break;

                    if (std::abs(fm) < config.tolerance) {
                        a = b = m;
                        break;
                    }
                    if ((fa < 0.0 && fm > 0.0) || (fa > 0.0 && fm < 0.0)) {
                        b = m;
                    } else {
                        a = m;
                        fa = fm;
                    }
                }
                roots.push_back(0.5 * (a + b));
            }

            prev_rho = cur_rho;
            prev_f = cur_f;
        }
    };

    if (roots.empty()) {
        // If requested, use spinodal bounds to focus searches.
        if (config.use_spinodal_bounds) {
            const auto spin = findSpinodalPoints(T, x, Config::SpinodalConfig{});
            if (spin.found && spin.has_gap &&
                spin.rho_vapor_spinodal > rho_min && spin.rho_liquid_spinodal < rho_max &&
                spin.rho_vapor_spinodal < spin.rho_liquid_spinodal) {

                const double rho_v_hi = std::max(rho_min, 0.999 * spin.rho_vapor_spinodal);
                const double rho_l_lo = std::min(rho_max, 1.001 * spin.rho_liquid_spinodal);

                scanInterval(rho_min, rho_v_hi, 350, roots);
                scanInterval(rho_l_lo, rho_max, 350, roots);
            } else {
                scanInterval(rho_min, rho_max, 500, roots);
            }
        } else {
            scanInterval(rho_min, rho_max, 500, roots);
        }

        // Always try a direct bracket around the ideal-gas density for the vapor root.
        {
            double a = std::max(rho_min, rho_ig * 1e-3);
            double b = std::min(rho_max, rho_ig * 1e3);
            double fa = f(a);
            double fb = f(b);
            if (std::isfinite(fa) && std::isfinite(fb) && fa * fb < 0.0) {
                // Bisection refine.
                for (int it = 0; it < config.max_iterations; ++it) {
                    const double m = 0.5 * (a + b);
                    const double fm = f(m);
                    if (!std::isfinite(fm)) break;
                    if (std::abs(fm) < config.tolerance) {
                        a = b = m;
                        break;
                    }
                    if (fa * fm < 0.0) {
                        b = m;
                        fb = fm;
                    } else {
                        a = m;
                        fa = fm;
                    }
                }
                roots.push_back(0.5 * (a + b));
            }
        }
    }

    // Sort and deduplicate.
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end(), [](double a, double b) {
        return std::abs(a - b) <= 1e-6 * std::max(1.0, std::max(a, b));
    }), roots.end());

    std::vector<double> kept_roots;
    std::vector<PhaseType> kept_phase_types;
    std::vector<bool> kept_mech;
    std::vector<double> kept_gibbs;
    kept_roots.reserve(roots.size());
    kept_phase_types.reserve(roots.size());
    kept_mech.reserve(roots.size());
    kept_gibbs.reserve(roots.size());

    for (size_t i = 0; i < roots.size(); ++i) {
        const double rho = roots[i];
        const double dPdrho = eos_->dPdrho(T, rho, x);
        const bool mech_stable = std::isfinite(dPdrho) && (dPdrho > 0.0);

        if (config.filter_metastable && !mech_stable) {
            continue;
        }

        PhaseType phase = (rho < 2.0 * rho_ig) ? PhaseType::Vapor : PhaseType::Liquid;
        kept_roots.push_back(rho);
        kept_phase_types.push_back(phase);
        kept_mech.push_back(mech_stable);

        const double g_res_over_rt = eos_->residualGibbs(T, rho, x);
        kept_gibbs.push_back(g_res_over_rt);
    }

    if (config.max_roots > 0 && static_cast<int>(kept_roots.size()) > config.max_roots) {
        const size_t keep = static_cast<size_t>(config.max_roots);
        kept_roots.resize(keep);
        kept_phase_types.resize(keep);
        kept_mech.resize(keep);
        kept_gibbs.resize(keep);
    }

    result.num_roots = static_cast<int>(kept_roots.size());
    result.densities = std::move(kept_roots);
    result.phase_types = std::move(kept_phase_types);
    result.is_mechanically_stable = std::move(kept_mech);
    result.gibbs_energies = std::move(kept_gibbs);

    int best_idx = -1;
    double best_g = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < result.densities.size(); ++i) {
        if (result.is_mechanically_stable[i] && result.gibbs_energies[i] < best_g) {
            best_g = result.gibbs_energies[i];
            best_idx = static_cast<int>(i);
        }
    }
    if (best_idx >= 0) {
        result.stable_root_index = best_idx;
    }

    return result;
}

double TPDStabilityAnalyzer::findVaporDensity(double T, double P, const std::vector<double>& x) const {
    auto roots = findDensityRoots(T, P, x);
    if (roots.densities.empty()) return 0.0;

    // Prefer the minimum mechanically-stable root; otherwise return the minimum root.
    double best = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < roots.densities.size(); ++i) {
        if (roots.is_mechanically_stable[i]) {
            best = std::min(best, roots.densities[i]);
        }
    }
    if (std::isfinite(best)) return best;
    return *std::min_element(roots.densities.begin(), roots.densities.end());
}

double TPDStabilityAnalyzer::findLiquidDensity(double T, double P, const std::vector<double>& x) const {
    auto roots = findDensityRoots(T, P, x);
    if (roots.densities.empty()) return 0.0;

    // Prefer the maximum mechanically-stable root; otherwise return the maximum root.
    double best = 0.0;
    bool found = false;
    for (size_t i = 0; i < roots.densities.size(); ++i) {
        if (roots.is_mechanically_stable[i]) {
            best = std::max(best, roots.densities[i]);
            found = true;
        }
    }
    if (found) return best;
    return *std::max_element(roots.densities.begin(), roots.densities.end());
}

SpinodalResult TPDStabilityAnalyzer::findSpinodalPoints(
    double T,
    const std::vector<double>& x,
    const Config::SpinodalConfig& config) const
{
    SpinodalResult result;
    result.temperature = T;

    auto g = [&](double rho) -> double {
        return eos_->dPdrho(T, rho, x);
    };

    const double rho_min = Constants::MIN_DENSITY;
    const double rho_max = Constants::MAX_DENSITY;

    // Scan derivative sign changes.
    constexpr int n_scan = 600;
    std::vector<double> grid(n_scan);
    const double log_min = std::log(rho_min);
    const double log_max = std::log(rho_max);
    for (int i = 0; i < n_scan; ++i) {
        const double t = static_cast<double>(i) / (n_scan - 1);
        grid[i] = std::exp(log_min + t * (log_max - log_min));
    }

    std::vector<double> spinodals;
    double prev_rho = grid[0];
    double prev_g = g(prev_rho);
    for (int i = 1; i < n_scan && static_cast<int>(spinodals.size()) < (config.find_both_sides ? 2 : 1); ++i) {
        const double cur_rho = grid[i];
        const double cur_g = g(cur_rho);
        if (!std::isfinite(prev_g) || !std::isfinite(cur_g)) {
            prev_rho = cur_rho;
            prev_g = cur_g;
            continue;
        }
        if ((prev_g < 0.0 && cur_g > 0.0) || (prev_g > 0.0 && cur_g < 0.0)) {
            double a = prev_rho, b = cur_rho;
            double ga = prev_g;
            for (int it = 0; it < config.max_iterations; ++it) {
                const double m = 0.5 * (a + b);
                const double gm = g(m);
                if (!std::isfinite(gm)) break;
                if (std::abs(gm) < config.tolerance) {
                    a = b = m;
                    break;
                }
                if ((ga < 0.0 && gm > 0.0) || (ga > 0.0 && gm < 0.0)) {
                    b = m;
                } else {
                    a = m;
                    ga = gm;
                }
            }
            spinodals.push_back(0.5 * (a + b));
        }
        prev_rho = cur_rho;
        prev_g = cur_g;
    }

    if (!spinodals.empty()) {
        std::sort(spinodals.begin(), spinodals.end());
        result.found = true;
        result.rho_vapor_spinodal = spinodals.front();
        result.rho_liquid_spinodal = spinodals.back();
        result.has_gap = spinodals.size() >= 2;
        result.pressure_spinodal = eos_->pressure(T, result.rho_vapor_spinodal, x);
    }

    return result;
}

bool TPDStabilityAnalyzer::isInSpinodalRegion(double T, double rho, const std::vector<double>& x) const {
    const double dPdrho = eos_->dPdrho(T, rho, x);
    return std::isfinite(dPdrho) && (dPdrho < 0.0);
}

PhaseType TPDStabilityAnalyzer::classifyPhase(double T, double rho, const std::vector<double>& x) const {
    const double P = eos_->pressure(T, rho, x);
    const double rho_ig = (P > 0.0) ? P / (Constants::GAS_CONSTANT * T) : 0.0;
    return (rho_ig > 0.0 && rho < 2.0 * rho_ig) ? PhaseType::Vapor : PhaseType::Liquid;
}

std::vector<std::vector<double>> TPDStabilityAnalyzer::generateTrialCompositions(
    double T,
    double P,
    const std::vector<double>& z,
    int num_trials) const
{
    const int nc = eos_->numComponents();
    if (static_cast<int>(z.size()) != nc) {
        throw std::invalid_argument("generateTrialCompositions: composition size mismatch");
    }

    std::vector<std::vector<double>> trials;
    trials.reserve(static_cast<size_t>(num_trials));

    // Always include feed composition.
    trials.push_back(z);

    // Wilson-based vapor-like and liquid-like trials.
    const auto K = wilsonK(*eos_, T, P);
    std::vector<double> w_vap(nc), w_liq(nc);
    for (int i = 0; i < nc; ++i) {
        w_vap[i] = z[i] / std::max(K[i], 1e-12);
        w_liq[i] = z[i] * std::max(K[i], 1e-12);
    }
    trials.push_back(normalize(w_vap));
    trials.push_back(normalize(w_liq));

    // Pure-component limits.
    for (int i = 0; i < nc && static_cast<int>(trials.size()) < num_trials; ++i) {
        std::vector<double> w(nc, (1.0 - 0.999) / std::max(1, nc - 1));
        w[i] = 0.999;
        trials.push_back(w);
    }

    // Random trials around z.
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    while (static_cast<int>(trials.size()) < num_trials) {
        std::vector<double> w(nc);
        for (int i = 0; i < nc; ++i) {
            w[i] = std::max(uni(rng), Constants::MIN_MOLE_FRACTION);
        }
        trials.push_back(normalize(w));
    }

    return trials;
}

double TPDStabilityAnalyzer::tangentPlaneDistance(
    double T,
    double P,
    const std::vector<double>& z,
    const std::vector<double>& w) const
{
    const auto cfg = Config::StabilityConfig::defaults();
    const auto z_state = lnPhiAtTP(T, P, z, PhaseType::Unknown, cfg, 0.0);
    const auto w_vap = lnPhiAtTP(T, P, w, PhaseType::Vapor, cfg, 0.0);
    const auto w_liq = lnPhiAtTP(T, P, w, PhaseType::Liquid, cfg, 0.0);

    const double tpd_v = tpdFromLnPhi(z, z_state.ln_phi, w, w_vap.ln_phi);
    const double tpd_l = tpdFromLnPhi(z, z_state.ln_phi, w, w_liq.ln_phi);
    return std::min(tpd_v, tpd_l);
}

TPDTrialResult TPDStabilityAnalyzer::minimizeTPDWithPhaseHint(
    double T,
    double P,
    const std::vector<double>& z,
    const std::vector<double>& lnphi_z,
    double tpd_scale,
    const std::vector<double>& w_initial,
    PhaseType phase_hint,
    const Config::StabilityConfig& config) const
{
    TPDTrialResult out;
    out.w_initial = w_initial;

    std::vector<double> w = normalize(w_initial);
    std::vector<double> w_new = w;
    double rho_guess = 0.0;

    for (int it = 0; it < config.max_iterations; ++it) {
        const auto w_state = lnPhiAtTP(T, P, w, phase_hint, config, rho_guess);
        rho_guess = w_state.rho;

        std::vector<double> unnormalized(w.size(), 0.0);
        for (size_t i = 0; i < w.size(); ++i) {
            const double zi = std::max(z[i], Constants::MIN_MOLE_FRACTION);
            unnormalized[i] = zi * std::exp(lnphi_z[i] - w_state.ln_phi[i]);
            unnormalized[i] = std::max(unnormalized[i], Constants::MIN_MOLE_FRACTION);
        }
        w_new = normalize(unnormalized);

        double max_delta = 0.0;
        for (size_t i = 0; i < w.size(); ++i) {
            max_delta = std::max(max_delta, std::abs(w_new[i] - w[i]));
        }

        w = w_new;
        out.iterations = it + 1;

        if (max_delta < config.tpd_tolerance) {
            out.converged = true;
            break;
        }
    }

    const auto w_final = lnPhiAtTP(T, P, w, phase_hint, config, rho_guess);
    out.w_converged = w;
    out.tpd = tpdFromLnPhi(z, lnphi_z, w, w_final.ln_phi);
    out.indicates_instability = tpdIndicatesInstability(out.tpd, tpd_scale, config);
    return out;
}

TPDTrialResult TPDStabilityAnalyzer::minimizeTPD(
    double T,
    double P,
    const std::vector<double>& z,
    const std::vector<double>& w_initial,
    const Config::StabilityConfig& config) const
{
    const auto z_state = lnPhiAtTP(T, P, z, PhaseType::Unknown, config, 0.0);
    const double tpd_scale = tpdReferenceScale(z, z_state.ln_phi);

    // Run minimization for both vapor-like and liquid-like branches; return the best (lowest TPD).
    auto vap = minimizeTPDWithPhaseHint(T, P, z, z_state.ln_phi, tpd_scale, w_initial, PhaseType::Vapor, config);
    auto liq = minimizeTPDWithPhaseHint(T, P, z, z_state.ln_phi, tpd_scale, w_initial, PhaseType::Liquid, config);

    if (vap.tpd <= liq.tpd) return vap;
    return liq;
}

StabilityResult TPDStabilityAnalyzer::analyze(
    double T,
    double P,
    const std::vector<double>& z,
    const Config::StabilityConfig& config) const
{
    StabilityResult result;
    if (config.method == Config::StabilityMethod::Spinodal) {
        result.method_used = "Spinodal (EOS)";
    } else if (config.method == Config::StabilityMethod::Combined) {
        result.method_used = "TPD + Spinodal (EOS)";
    } else {
        result.method_used = "TPD (EOS)";
    }

    const int nc = eos_->numComponents();
    if (static_cast<int>(z.size()) != nc) {
        result.is_stable = true;
        result.notes = "Invalid composition size";
        return result;
    }

    // Optional spinodal analysis.
    if (config.detect_spinodal_gap) {
        auto sp_cfg = Config::SpinodalConfig{};
        sp_cfg.tolerance = config.spinodal_tolerance;
        sp_cfg.max_iterations = config.max_iterations;
        auto spin = findSpinodalPoints(T, z, sp_cfg);
        result.spinodal_analyzed = spin.found;
        if (spin.found) {
            result.rho_spinodal_vapor = spin.rho_vapor_spinodal;
            if (spin.has_gap) {
                result.rho_spinodal_liquid = spin.rho_liquid_spinodal;
            }
        }
    }

    bool has_feed_mech_assessment = false;
    if (config.method == Config::StabilityMethod::Spinodal || config.method == Config::StabilityMethod::Combined) {
        auto roots_cfg = densityRootConfigFromStability(config);
        roots_cfg.filter_metastable = false;
        const auto roots = findDensityRoots(T, P, z, roots_cfg);
        if (!roots.densities.empty()) {
            const bool any_mech_stable = std::any_of(
                roots.is_mechanically_stable.begin(),
                roots.is_mechanically_stable.end(),
                [](bool stable) { return stable; });
            result.in_spinodal_region = !any_mech_stable;
            has_feed_mech_assessment = true;
        }
    }

    if (config.method == Config::StabilityMethod::Spinodal) {
        if (!has_feed_mech_assessment) {
            // Fallback for rare cases where no TP roots were detected.
            const double rho = std::max(P / (Constants::GAS_CONSTANT * T), Constants::MIN_DENSITY);
            result.in_spinodal_region = isInSpinodalRegion(T, rho, z);
        }
        result.is_stable = !result.in_spinodal_region;
        return result;
    }

    // Run TPD trials.
    std::vector<std::vector<double>> trials;
    const int target_trials = std::max(1, config.num_trial_compositions);
    trials.reserve(static_cast<size_t>(target_trials));

    auto addTrial = [&](const std::vector<double>& w_in) {
        if (static_cast<int>(trials.size()) >= target_trials) return;
        if (static_cast<int>(w_in.size()) != nc) return;
        trials.push_back(normalize(w_in));
    };

    addTrial(z);

    Config::TrialMethod trial_method = config.trial_method;
    if (trial_method == Config::TrialMethod::Custom) {
        if (!config.custom_trials.empty()) {
            for (const auto& w : config.custom_trials) {
                addTrial(w);
            }
        } else {
            trial_method = Config::TrialMethod::Wilson;
        }
    }

    if (trial_method == Config::TrialMethod::Wilson) {
        const auto K = wilsonK(*eos_, T, P);
        std::vector<double> w_vap(nc, 0.0), w_liq(nc, 0.0);
        for (int i = 0; i < nc; ++i) {
            w_vap[i] = z[i] / std::max(K[i], 1e-12);
            w_liq[i] = z[i] * std::max(K[i], 1e-12);
        }
        addTrial(w_vap);
        addTrial(w_liq);
    } else if (trial_method == Config::TrialMethod::IdealSolution) {
        for (int i = 0; i < nc && static_cast<int>(trials.size()) < target_trials; ++i) {
            std::vector<double> w(static_cast<size_t>(nc), (1.0 - 0.999) / std::max(1, nc - 1));
            w[static_cast<size_t>(i)] = 0.999;
            addTrial(w);
        }
    }

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    while (static_cast<int>(trials.size()) < target_trials) {
        std::vector<double> w(static_cast<size_t>(nc), 0.0);
        for (int i = 0; i < nc; ++i) {
            w[static_cast<size_t>(i)] = std::max(uni(rng), Constants::MIN_MOLE_FRACTION);
        }
        addTrial(w);
    }

    if (trials.size() > static_cast<size_t>(target_trials)) {
        trials.resize(static_cast<size_t>(target_trials));
    }

    const auto z_state = lnPhiAtTP(T, P, z, PhaseType::Unknown, config, 0.0);
    const double tpd_scale = tpdReferenceScale(z, z_state.ln_phi);

    double best_tpd = std::numeric_limits<double>::infinity();
    std::vector<double> best_w;
    bool best_converged = false;
    bool have_best = false;
    int converged = 0;

    for (const auto& w0 : trials) {
        TPDTrialResult trial;
        try {
            trial = minimizeTPD(T, P, z, w0, config);
        } catch (...) {
            continue;
        }
        if (trial.converged) {
            ++converged;
        }
        if (config.store_all_trials) {
            result.trial_results.push_back(trial);
        }

        if (!std::isfinite(trial.tpd)) {
            continue;
        }

        if (!have_best) {
            best_tpd = trial.tpd;
            best_w = trial.w_converged;
            best_converged = trial.converged;
            have_best = true;
        } else {
            const bool prefer = (trial.converged && !best_converged) ||
                ((trial.converged == best_converged) && (trial.tpd < best_tpd));
            if (prefer) {
                best_tpd = trial.tpd;
                best_w = trial.w_converged;
                best_converged = trial.converged;
            }
        }
    }

    if (!have_best) {
        throw std::runtime_error("TPDStabilityAnalyzer::analyze: failed to evaluate any valid TPD trial");
    }

    result.num_trials = static_cast<int>(trials.size());
    result.num_converged = converged;
    result.min_tpd = std::isfinite(best_tpd) ? best_tpd : 0.0;
    const bool tpd_unstable = tpdIndicatesInstability(best_tpd, tpd_scale, config);
    result.is_stable = !tpd_unstable;
    if (config.method == Config::StabilityMethod::Combined && result.in_spinodal_region) {
        result.is_stable = false;
    }

    if (tpd_unstable && !best_w.empty()) {
        result.unstable_composition = best_w;

        // Identify whether the incipient phase is vapor-like or liquid-like by re-evaluating
        // TPD using density roots consistent with each phase hint.
        const auto w_vap = lnPhiAtTP(T, P, best_w, PhaseType::Vapor, config, 0.0);
        const auto w_liq = lnPhiAtTP(T, P, best_w, PhaseType::Liquid, config, 0.0);

        const double tpd_v = tpdFromLnPhi(z, z_state.ln_phi, best_w, w_vap.ln_phi);
        const double tpd_l = tpdFromLnPhi(z, z_state.ln_phi, best_w, w_liq.ln_phi);
        if (std::isfinite(tpd_v) && std::isfinite(tpd_l)) {
            result.incipient_phase_type = (tpd_v <= tpd_l) ? PhaseType::Vapor : PhaseType::Liquid;
        }
    }

    return result;
}

bool TPDStabilityAnalyzer::isStable(double T, double P, const std::vector<double>& z) const {
    auto cfg = Config::StabilityConfig::defaults();
    cfg.num_trial_compositions = 3;
    cfg.max_iterations = 50;
    auto res = analyze(T, P, z, cfg);
    return res.is_stable;
}

} // namespace Stability
} // namespace Equilibrium
} // namespace DMThermo
