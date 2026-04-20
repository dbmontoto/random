/**
 * @file gibbs_multi_flash.cpp
 * @brief Implementation of EOS-agnostic multi-phase Gibbs flash solver
 */

#include "thermo/equilibrium/gibbs_multi_flash.h"

#include "thermo/core/constants.h"
#include "thermo/equilibrium/density_selection.h"
#include "thermo/equilibrium/internal/hs_support.h"
#include "thermo/equilibrium/internal/k_values.h"
#include "thermo/equilibrium/internal/rachford_rice.h"
#include "thermo/equilibrium/internal/temperature_solve.h"
#include "thermo/equilibrium/thermo_contributions.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace DMThermo {
namespace Equilibrium {
namespace Flash {

namespace {

constexpr double X_MIN = Constants::MIN_MOLE_FRACTION;
constexpr double BETA_MIN = 1e-12;

std::vector<double> normalize(const std::vector<double>& x) {
    double sum = std::accumulate(x.begin(), x.end(), 0.0);
    if (!(std::isfinite(sum) && sum > 0.0)) {
        throw std::invalid_argument("normalize: non-positive sum");
    }
    std::vector<double> out(x.size());
    for (size_t i = 0; i < x.size(); ++i) out[i] = x[i] / sum;
    return out;
}

std::vector<double> clampAndNormalize(const std::vector<double>& x) {
    std::vector<double> tmp = x;
    for (auto& v : tmp) {
        if (!std::isfinite(v) || v < X_MIN) v = X_MIN;
    }
    return normalize(tmp);
}

double safeLog(double v) {
    return std::log(std::max(v, X_MIN));
}

struct PhaseWork {
    PhaseType hint = PhaseType::Unknown;
    double beta = 0.0;
    std::vector<double> x;
    double rho = 0.0;
    std::vector<double> phi;
    std::vector<double> lnphi;
    PhaseType classified = PhaseType::Unknown;

    PhaseWork() = default;
    PhaseWork(PhaseType hint_in, double beta_in, std::vector<double> x_in)
        : hint(hint_in), beta(beta_in), x(std::move(x_in)) {}
};

bool gaussSolve(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& x) {
    const int n = static_cast<int>(b.size());
    x.assign(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        // Pivot
        int piv = i;
        double best = std::abs(A[i][i]);
        for (int r = i + 1; r < n; ++r) {
            const double v = std::abs(A[r][i]);
            if (v > best) {
                best = v;
                piv = r;
            }
        }
        if (!(best > 0.0) || !std::isfinite(best)) return false;
        if (piv != i) {
            std::swap(A[piv], A[i]);
            std::swap(b[piv], b[i]);
        }

        const double diag = A[i][i];
        for (int c = i; c < n; ++c) A[i][c] /= diag;
        b[i] /= diag;

        for (int r = 0; r < n; ++r) {
            if (r == i) continue;
            const double f = A[r][i];
            if (f == 0.0) continue;
            for (int c = i; c < n; ++c) A[r][c] -= f * A[i][c];
            b[r] -= f * b[i];
        }
    }
    x = b;
    return true;
}

bool solveBetasNewton(
    const std::vector<double>& z,
    const std::vector<std::vector<double>>& K, // [nc][M], K[:,ref]=1
    std::vector<double>& beta)
{
    const int nc = static_cast<int>(z.size());
    const int M = static_cast<int>(beta.size());
    if (M < 1) return false;

    auto normalizeBeta = [&]() {
        for (auto& b : beta) {
            if (!std::isfinite(b) || b < 0.0) b = 0.0;
        }
        double s = std::accumulate(beta.begin(), beta.end(), 0.0);
        if (!(s > 0.0)) {
            beta.assign(static_cast<size_t>(M), 1.0 / std::max(1, M));
            return;
        }
        for (auto& b : beta) b /= s;
    };

    normalizeBeta();

    auto compute = [&](std::vector<double>& D, std::vector<double>& E) -> bool {
        D.assign(static_cast<size_t>(nc), 0.0);
        for (int i = 0; i < nc; ++i) {
            double di = 0.0;
            for (int m = 0; m < M; ++m) di += beta[m] * K[i][m];
            if (!(std::isfinite(di) && di > 0.0)) return false;
            D[i] = di;
        }

        E.assign(static_cast<size_t>(M), 0.0);
        for (int k = 0; k < M - 1; ++k) {
            double s = 0.0;
            for (int i = 0; i < nc; ++i) s += z[i] * K[i][k] / D[i];
            E[k] = s - 1.0;
        }
        E[M - 1] = std::accumulate(beta.begin(), beta.end(), 0.0) - 1.0;
        return true;
    };

    for (int it = 0; it < 50; ++it) {
        std::vector<double> D, E;
        if (!compute(D, E)) return false;

        double maxE = 0.0;
        for (double e : E) maxE = std::max(maxE, std::abs(e));
        if (maxE < 1e-14) return true;

        std::vector<std::vector<double>> J(static_cast<size_t>(M), std::vector<double>(static_cast<size_t>(M), 0.0));
        for (int k = 0; k < M - 1; ++k) {
            for (int m = 0; m < M; ++m) {
                double s = 0.0;
                for (int i = 0; i < nc; ++i) {
                    const double di = D[i];
                    s += z[i] * K[i][k] * K[i][m] / (di * di);
                }
                J[k][m] = -s;
            }
        }
        for (int m = 0; m < M; ++m) {
            J[M - 1][m] = 1.0;
        }

        std::vector<double> rhs(static_cast<size_t>(M));
        for (int r = 0; r < M; ++r) rhs[r] = -E[r];

        std::vector<double> delta;
        if (!gaussSolve(J, rhs, delta)) return false;

        // Damped step to keep D_i positive.
        double alpha = 1.0;
        for (int ls = 0; ls < 20; ++ls) {
            std::vector<double> trial = beta;
            for (int m = 0; m < M; ++m) trial[m] += alpha * delta[m];

            // Project to simplex quickly.
            for (auto& b : trial) if (!std::isfinite(b) || b < 0.0) b = 0.0;
            double s = std::accumulate(trial.begin(), trial.end(), 0.0);
            if (!(s > 0.0)) { alpha *= 0.5; continue; }
            for (auto& b : trial) b /= s;

            bool ok = true;
            for (int i = 0; i < nc; ++i) {
                double di = 0.0;
                for (int m = 0; m < M; ++m) di += trial[m] * K[i][m];
                if (!(std::isfinite(di) && di > 0.0)) { ok = false; break; }
            }
            if (ok) {
                beta = std::move(trial);
                break;
            }
            alpha *= 0.5;
        }

        normalizeBeta();
    }

    return true;
}

FlashResult finalizeResult(
    const EOS& eos,
    const Stability::IStabilityAnalyzer& stability,
    double T,
    double P,
    const std::vector<double>& z,
    std::vector<PhaseWork> phases,
    Config::PhaseDetection space,
    int max_density_newton_iters,
    const Core::Units::PropertyVariable& pressure_tolerance,
    const std::string& method_used,
    int iterations,
    double residual)
{
    FlashResult out;
    out.temperature = T;
    out.pressure = P;
    out.z = z;
    out.method_used = method_used;
    out.iterations = iterations;
    out.residual = residual;

    // Re-evaluate final densities/phis and classify.
    for (auto& ph : phases) {
        ph.x = clampAndNormalize(ph.x);
        auto ev = Detail::evalAtTP(
            eos, stability, T, P, ph.x, ph.hint,
            space, ph.rho, max_density_newton_iters, pressure_tolerance
        );
        ph.rho = ev.rho;
        ph.phi = std::move(ev.phi);
        ph.lnphi = std::move(ev.lnphi);
        ph.classified = ev.phase;
    }

    // Sort phases: Vapor, Liquid/Liquid1, Liquid2.
    std::stable_sort(phases.begin(), phases.end(), [](const PhaseWork& a, const PhaseWork& b) {
        auto rank = [](PhaseType p) {
            if (p == PhaseType::Vapor) return 0;
            if (p == PhaseType::Liquid || p == PhaseType::Liquid1) return 1;
            if (p == PhaseType::Liquid2) return 2;
            return 3;
        };
        return rank(a.classified) < rank(b.classified);
    });

    out.num_phases = static_cast<int>(phases.size());
    out.is_two_phase = (out.num_phases == 2);
    out.is_three_phase = (out.num_phases == 3);

    double vf = 0.0;
    out.phases.clear();
    for (size_t k = 0; k < phases.size(); ++k) {
        PhaseState ps;
        ps.type = phases[k].classified;
        if (ps.type == PhaseType::Liquid) ps.type = PhaseType::Liquid1;
        ps.density = phases[k].rho;
        ps.compressibility = eos.compressibility(T, ps.density, phases[k].x);
        ps.x = phases[k].x;
        ps.phi = phases[k].phi;
        ps.fraction = phases[k].beta;
        if (ps.type == PhaseType::Vapor) vf = ps.fraction;
        out.phases.push_back(std::move(ps));
    }
    out.vapor_fraction = vf;

    out.converged = true;
    out.message = (out.num_phases == 1) ? "Single phase" : "Multi-phase Gibbs converged";
    return out;
}

} // namespace

GibbsMultiPhaseFlashSolver::GibbsMultiPhaseFlashSolver(EOSPtr eos, Stability::StabilityAnalyzerPtr stability)
    : eos_(std::move(eos))
    , stability_(std::move(stability))
{
    if (!eos_) throw std::invalid_argument("GibbsMultiPhaseFlashSolver: EOS is null");
    if (!stability_) throw std::invalid_argument("GibbsMultiPhaseFlashSolver: stability analyzer is null");
}

std::vector<double> GibbsMultiPhaseFlashSolver::estimateKValues(double T, double P) const {
    return Internal::wilsonK(*eos_, T, P);
}

bool GibbsMultiPhaseFlashSolver::isValidResult(const FlashResult& result) const {
    if (!result.converged) return false;
    if (result.temperature <= 0.0 || result.pressure <= 0.0) return false;
    if (result.num_phases < 1 || result.num_phases > 3) return false;
    if (static_cast<int>(result.phases.size()) != result.num_phases) return false;
    double sum = 0.0;
    for (const auto& ph : result.phases) {
        sum += ph.fraction;
        if (!(ph.density > 0.0)) return false;
        double s = 0.0;
        for (double xi : ph.x) {
            if (!(xi >= -1e-12 && xi <= 1.0 + 1e-12)) return false;
            s += xi;
        }
        if (std::abs(s - 1.0) > 1e-6) return false;
    }
    return std::abs(sum - 1.0) < 1e-6;
}

FlashResult GibbsMultiPhaseFlashSolver::solvePT(
    double T,
    double P,
    const std::vector<double>& z_in,
    const Config::FlashConfig& config) const
{
    FlashResult fail;
    fail.temperature = T;
    fail.pressure = P;
    fail.z = z_in;
    fail.method_used = algorithmName();

    const int nc = eos_->numComponents();
    if (static_cast<int>(z_in.size()) != nc) {
        fail.converged = false;
        fail.message = "Composition size mismatch";
        return fail;
    }
    if (!(T > 0.0 && P > 0.0)) {
        fail.converged = false;
        fail.message = "Invalid T/P";
        return fail;
    }

    const std::vector<double> z = clampAndNormalize(z_in);
    const int max_phases = std::clamp(config.max_phases, 1, 3);
    const Core::Units::PropertyVariable pressure_tolerance(config.density_newton_tol, Core::Units::Unit::PA);

    // Fast path: single-phase or explicitly limited to one phase.
    if (max_phases == 1) {
        PhaseWork ph;
        ph.beta = 1.0;
        ph.hint = PhaseType::Unknown;
        ph.x = z;
        return finalizeResult(
            *eos_, *stability_, T, P, z, {ph},
            config.phase_detection_space, config.max_density_newton_iters, pressure_tolerance,
            algorithmName(), 0, 0.0
        );
    }

    // Stability test on feed.
    if (config.perform_stability_test) {
        auto stab_cfg = Config::StabilityConfig::defaults();
        stab_cfg.num_trial_compositions = config.num_stability_trials;
        stab_cfg.phase_detection_space = config.phase_detection_space;
        stab_cfg.max_density_newton_iters = config.max_density_newton_iters;
        stab_cfg.density_newton_tol = config.density_newton_tol;
        auto stab = stability_->analyze(T, P, z, stab_cfg);
        if (stab.is_stable) {
            PhaseWork ph;
            ph.beta = 1.0;
            ph.hint = PhaseType::Unknown;
            ph.x = z;
            auto out = finalizeResult(
                *eos_, *stability_, T, P, z, {ph},
                config.phase_detection_space, config.max_density_newton_iters, pressure_tolerance,
                algorithmName(), 0, 0.0
            );
            out.stability_test_performed = true;
            out.feed_was_stable = true;
            return out;
        }
    }

    // Initial two-phase guess from K-values.
    std::vector<double> K0;
    if (config.k_init == Config::KValueInit::Custom && !config.custom_k_values.empty()) {
        if (static_cast<int>(config.custom_k_values.size()) != nc) {
            throw std::invalid_argument("GibbsMultiPhaseFlashSolver: custom_k_values size mismatch");
        }
        K0 = config.custom_k_values;
    } else {
        K0 = estimateKValues(T, P);
    }

    // Compute a 2-phase initial split using standard RR.
    double beta_v = 0.5;
    (void)Internal::solveRachfordRiceBeta(z, K0, beta_v);
    beta_v = std::clamp(beta_v, 1e-6, 1.0 - 1e-6);

    std::vector<double> x_liq(nc, 0.0), y_vap(nc, 0.0);
    for (int i = 0; i < nc; ++i) {
        const double denom = 1.0 + beta_v * (K0[i] - 1.0);
        x_liq[i] = z[i] / std::max(denom, X_MIN);
        y_vap[i] = K0[i] * x_liq[i];
    }
    x_liq = clampAndNormalize(x_liq);
    y_vap = clampAndNormalize(y_vap);

    std::vector<PhaseWork> phases;
    phases.emplace_back(PhaseType::Vapor, beta_v, std::move(y_vap));
    phases.emplace_back(PhaseType::Liquid, 1.0 - beta_v, std::move(x_liq));

    auto stab_cfg = Config::StabilityConfig::defaults();
    stab_cfg.num_trial_compositions = std::max(5, config.num_stability_trials);
    stab_cfg.phase_detection_space = config.phase_detection_space;
    stab_cfg.max_density_newton_iters = config.max_density_newton_iters;
    stab_cfg.density_newton_tol = config.density_newton_tol;

    double last_resid = std::numeric_limits<double>::infinity();
    int last_iters = 0;

    for (int outer = 0; outer < 6; ++outer) {
        // Inner equilibrium iteration.
        bool converged = false;
        for (int it = 0; it < config.max_iterations; ++it) {
            const int M = static_cast<int>(phases.size());

            // Evaluate lnphi for each phase (and classify).
            for (auto& ph : phases) {
                ph.x = clampAndNormalize(ph.x);
                auto ev = Detail::evalAtTP(
                    *eos_, *stability_, T, P, ph.x, ph.hint,
                    config.phase_detection_space, ph.rho,
                    config.max_density_newton_iters, pressure_tolerance
                );
                ph.rho = ev.rho;
                ph.phi = std::move(ev.phi);
                ph.lnphi = std::move(ev.lnphi);
                ph.classified = ev.phase;
            }

            // Choose reference phase: vapor if present; otherwise phase 0.
            int ref = 0;
            for (int k = 0; k < M; ++k) {
                if (phases[k].classified == PhaseType::Vapor) { ref = k; break; }
            }

            // Build K matrix [nc][M], with K[:,ref]=1.
            std::vector<std::vector<double>> K(static_cast<size_t>(nc), std::vector<double>(static_cast<size_t>(M), 1.0));
            for (int i = 0; i < nc; ++i) {
                const double lnphi_ref = phases[ref].lnphi[i];
                for (int k = 0; k < M; ++k) {
                    if (k == ref) { K[i][k] = 1.0; continue; }
                    const double lnK = lnphi_ref - phases[k].lnphi[i];
                    double v = std::exp(lnK);
                    if (!std::isfinite(v) || v <= 0.0) v = 1.0;
                    K[i][k] = std::clamp(v, config.min_k_value, config.max_k_value);
                }
            }

            // Solve for beta via generalized RR.
            std::vector<double> beta(static_cast<size_t>(M), 0.0);
            for (int k = 0; k < M; ++k) beta[k] = std::max(0.0, phases[k].beta);
            if (!solveBetasNewton(z, K, beta)) {
                break;
            }

            // Compute phase compositions from beta and K.
            std::vector<double> D(static_cast<size_t>(nc), 0.0);
            for (int i = 0; i < nc; ++i) {
                double di = 0.0;
                for (int k = 0; k < M; ++k) di += beta[k] * K[i][k];
                D[i] = di;
            }

            for (int k = 0; k < M; ++k) {
                phases[k].beta = beta[k];
                for (int i = 0; i < nc; ++i) {
                    phases[k].x[i] = z[i] * K[i][k] / std::max(D[i], X_MIN);
                }
                phases[k].x = clampAndNormalize(phases[k].x);
            }

            // Remove negligible phases.
            phases.erase(std::remove_if(phases.begin(), phases.end(), [](const PhaseWork& ph) {
                return ph.beta < BETA_MIN;
            }), phases.end());
            if (phases.empty()) {
                break;
            }
            // Renormalize betas.
            double bsum = 0.0;
            for (const auto& ph : phases) bsum += ph.beta;
            if (bsum > 0.0) {
                for (auto& ph : phases) ph.beta /= bsum;
            }

            // Residual: max deviation of ln(f_i) between phases.
            const int M2 = static_cast<int>(phases.size());
            // Re-evaluate lnphi after possible phase removal.
            for (auto& ph : phases) {
                auto ev = Detail::evalAtTP(
                    *eos_, *stability_, T, P, ph.x, ph.hint,
                    config.phase_detection_space, ph.rho,
                    config.max_density_newton_iters, pressure_tolerance
                );
                ph.rho = ev.rho;
                ph.phi = std::move(ev.phi);
                ph.lnphi = std::move(ev.lnphi);
                ph.classified = ev.phase;
            }
            int ref2 = 0;
            for (int k = 0; k < M2; ++k) {
                if (phases[k].classified == PhaseType::Vapor) { ref2 = k; break; }
            }

            double resid = 0.0;
            for (int i = 0; i < nc; ++i) {
                const double lnfr = safeLog(phases[ref2].x[i]) + phases[ref2].lnphi[i];
                for (int k = 0; k < M2; ++k) {
                    const double lnfk = safeLog(phases[k].x[i]) + phases[k].lnphi[i];
                    resid = std::max(resid, std::abs(lnfk - lnfr));
                }
            }

            last_resid = resid;
            last_iters = it + 1;

            if (resid < config.tolerance) {
                converged = true;
                break;
            }
        }

        if (phases.size() == 1) {
            return finalizeResult(
                *eos_, *stability_, T, P, z, phases,
                config.phase_detection_space, config.max_density_newton_iters, pressure_tolerance,
                algorithmName(), last_iters, last_resid
            );
        }

        // Phase stability checks and phase addition (TPD-driven).
        bool added = false;
        if (max_phases > static_cast<int>(phases.size())) {
            for (const auto& ph : phases) {
                auto st = stability_->analyze(T, P, ph.x, stab_cfg);
                if (!st.is_stable && st.unstable_composition.has_value()) {
                    auto w = clampAndNormalize(st.unstable_composition.value());
                    // Skip if close to existing phases.
                    bool distinct = true;
                    for (const auto& ex : phases) {
                        double l1 = 0.0;
                        for (int i = 0; i < nc; ++i) l1 += std::abs(w[i] - ex.x[i]);
                        if (l1 < 1e-4) { distinct = false; break; }
                    }
                    if (!distinct) continue;

                    PhaseType hint = st.incipient_phase_type.value_or(PhaseType::Liquid);
                    if (hint == PhaseType::Liquid) hint = PhaseType::Liquid2;
                    phases.emplace_back(hint, 1e-3, std::move(w));
                    // Renormalize betas.
                    double s = 0.0;
                    for (auto& p : phases) s += p.beta;
                    for (auto& p : phases) p.beta /= s;
                    added = true;
                    break;
                }
            }
        }

        if (!added) {
            // Finalize even if not perfectly converged: return best effort with converged=false.
            auto out = finalizeResult(
                *eos_, *stability_, T, P, z, phases,
                config.phase_detection_space, config.max_density_newton_iters, pressure_tolerance,
                algorithmName(), last_iters, last_resid
            );
            if (!converged) {
                out.converged = false;
                out.message = "Multi-phase Gibbs did not converge";
            }
            return out;
        }
    }

    auto out = finalizeResult(
        *eos_, *stability_, T, P, z, phases,
        config.phase_detection_space, config.max_density_newton_iters, pressure_tolerance,
        algorithmName(), last_iters, last_resid
    );
    out.converged = false;
    out.message = "Multi-phase phase-addition loop did not converge";
    return out;
}

FlashResult GibbsMultiPhaseFlashSolver::solveVLLE(
    double T,
    double P,
    const std::vector<double>& z,
    const Config::FlashConfig& config) const
{
    Config::FlashConfig cfg = config;
    cfg.max_phases = std::max(cfg.max_phases, 3);
    return solvePT(T, P, z, cfg);
}

FlashResult GibbsMultiPhaseFlashSolver::solveTV(
    double,
    double,
    const std::vector<double>&,
    const Config::TVFlashConfig&) const
{
    throw std::runtime_error("GibbsMultiPhaseFlashSolver: TV flash not implemented");
}

FlashResult GibbsMultiPhaseFlashSolver::solvePH(
    double P,
    double H,
    const std::vector<double>& z,
    const Config::PHFlashConfig& config) const
{
    const auto& mix = Internal::requireCoreMixture(*eos_, "GibbsMultiPhaseFlashSolver::solvePH", "(PR/SRK)");
    Internal::requireIdealGasCp(mix, "GibbsMultiPhaseFlashSolver::solvePH");

    auto computeH = [&](double T, FlashResult& fr_out) -> bool {
        try {
            fr_out = solvePT(T, P, z, config);
            if (!fr_out.converged || fr_out.phases.empty()) return false;

            double H_total = 0.0;
            for (const auto& phase : fr_out.phases) {
                const auto ig = Contributions::idealGasProps(mix, T, P, phase.x);
                const auto rr = Internal::residualPropsForHS(*eos_, T, phase.density, phase.x, P, "GibbsMultiPhaseFlashSolver::solvePH");
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
        throw std::invalid_argument("GibbsMultiPhaseFlashSolver: invalid PH temperature bracket");
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

FlashResult GibbsMultiPhaseFlashSolver::solvePS(
    double P,
    double S,
    const std::vector<double>& z,
    const Config::PSFlashConfig& config) const
{
    const auto& mix = Internal::requireCoreMixture(*eos_, "GibbsMultiPhaseFlashSolver::solvePS", "(PR/SRK)");
    Internal::requireIdealGasCp(mix, "GibbsMultiPhaseFlashSolver::solvePS");

    auto computeS = [&](double T, FlashResult& fr_out) -> bool {
        try {
            fr_out = solvePT(T, P, z, config);
            if (!fr_out.converged || fr_out.phases.empty()) return false;

            double S_total = 0.0;
            for (const auto& phase : fr_out.phases) {
                const auto ig = Contributions::idealGasProps(mix, T, P, phase.x);
                const auto rr = Internal::residualPropsForHS(*eos_, T, phase.density, phase.x, P, "GibbsMultiPhaseFlashSolver::solvePS");
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
        throw std::invalid_argument("GibbsMultiPhaseFlashSolver: invalid PS temperature bracket");
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

} // namespace Flash
} // namespace Equilibrium
} // namespace DMThermo
