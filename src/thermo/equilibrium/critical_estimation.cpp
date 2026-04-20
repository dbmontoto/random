/**
 * @file critical_estimation.cpp
 * @brief EOS-based critical point estimation implementation.
 */

#include "thermo/equilibrium/critical_estimation.h"

#include "thermo/core/constants.h"
#include "thermo/core/units.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace DMThermo {
namespace Equilibrium {
namespace Critical {

namespace {

double norm2(double a, double b) {
    return a * a + b * b;
}

} // namespace

EstimateResult estimateFromEOS(
    const EOS& eos,
    const std::vector<double>& x,
    const EstimateInputs& inputs)
{
    EstimateResult out;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double R = Constants::GAS_CONSTANT;
    const auto temperature_hint = inputs.temperatureHint();
    const auto pressure_hint = inputs.pressureHint();
    const auto critical_volume_hint = inputs.criticalVolumeHint();

    const int nc = eos.numComponents();
    if (nc <= 0) {
        out.message = "EOS has no components";
        return out;
    }

    std::vector<double> composition = x;
    if (composition.empty() && nc == 1) {
        composition = {1.0};
    }
    if (static_cast<int>(composition.size()) != nc) {
        out.message = "Composition size mismatch for critical estimation";
        return out;
    }
    double sum_x = 0.0;
    for (double value : composition) {
        if (!(std::isfinite(value) && value >= 0.0)) {
            out.message = "Composition contains invalid entries";
            return out;
        }
        sum_x += value;
    }
    if (!(std::isfinite(sum_x) && sum_x > 0.0)) {
        out.message = "Composition sum must be finite and > 0";
        return out;
    }
    for (double& value : composition) {
        value /= sum_x;
    }

    double T0 = temperature_hint.as(Core::Units::Unit::K);
    if (!(std::isfinite(T0) && T0 > 0.0)) {
        T0 = 300.0;
    }

    double rho0 = nan;
    const double Vc_hint = critical_volume_hint.as(Core::Units::Unit::M3_PER_KMOL);
    const double Pc_hint = pressure_hint.as(Core::Units::Unit::PA_ABS);
    if (std::isfinite(Vc_hint) && Vc_hint > 0.0) {
        const double critical_volume_molar = Core::Units::PropertyVariable(Vc_hint, Core::Units::Unit::M3_PER_KMOL)
            .as(Core::Units::Unit::M3_PER_MOL);
        if (std::isfinite(critical_volume_molar) && critical_volume_molar > 0.0) {
            rho0 = 1.0 / critical_volume_molar;
        }
    } else if (std::isfinite(Pc_hint) && Pc_hint > 0.0) {
        constexpr double z_guess = 0.30;
        rho0 = Pc_hint / (z_guess * R * T0);
    } else {
        rho0 = 1000.0;
    }
    if (!(std::isfinite(rho0) && rho0 > 0.0)) {
        rho0 = 1000.0;
    }

    auto evalTRho = [&](double T, double rho, double& f1, double& f2) -> bool {
        if (!(std::isfinite(T) && T > 0.0 && std::isfinite(rho) && rho > 0.0)) return false;
        try {
            const double dPdrho = eos.dPdrho(T, rho, composition);
            if (!std::isfinite(dPdrho)) return false;
            f1 = dPdrho / (R * T);

            double h = std::max(rho * 1e-6, 1e-12);
            if (rho <= h) h = 0.5 * rho;
            if (!(h > 0.0)) return false;

            const double dPdrho_p = eos.dPdrho(T, rho + h, composition);
            const double dPdrho_m = eos.dPdrho(T, rho - h, composition);
            if (!(std::isfinite(dPdrho_p) && std::isfinite(dPdrho_m))) return false;

            const double d2Pdrho2 = (dPdrho_p - dPdrho_m) / (2.0 * h);
            if (!std::isfinite(d2Pdrho2)) return false;

            f2 = (d2Pdrho2 * rho) / (R * T);
            return std::isfinite(f1) && std::isfinite(f2);
        } catch (...) {
            return false;
        }
    };

    // ---------------------------------------------------------------------
    // Attempt 1: damped Newton solve in log(T), log(rho)
    // ---------------------------------------------------------------------
    {
        double u = std::log(T0);
        double v = std::log(rho0);

        auto eval = [&](double u_in, double v_in, double& f1, double& f2) -> bool {
            return evalTRho(std::exp(u_in), std::exp(v_in), f1, f2);
        };

        constexpr int max_it = 60;
        constexpr double tol = 1e-10;

        double f1 = 0.0;
        double f2 = 0.0;
        if (eval(u, v, f1, f2)) {
            for (int it = 0; it < max_it; ++it) {
                if (!eval(u, v, f1, f2)) {
                    out.iterations = it;
                    out.message = "Newton: failed evaluation";
                    break;
                }

                const double n2 = norm2(f1, f2);
                if (std::isfinite(n2) && n2 < tol * tol) {
                    out.converged = true;
                    out.iterations = it;
                    break;
                }

                constexpr double du = 1e-4;
                constexpr double dv = 1e-4;

                double f1_up = 0.0;
                double f2_up = 0.0;
                double f1_um = 0.0;
                double f2_um = 0.0;
                double f1_vp = 0.0;
                double f2_vp = 0.0;
                double f1_vm = 0.0;
                double f2_vm = 0.0;
                if (!eval(u + du, v, f1_up, f2_up) || !eval(u - du, v, f1_um, f2_um) ||
                    !eval(u, v + dv, f1_vp, f2_vp) || !eval(u, v - dv, f1_vm, f2_vm))
                {
                    out.iterations = it;
                    out.message = "Newton: failed Jacobian evaluation";
                    break;
                }

                const double a = (f1_up - f1_um) / (2.0 * du);
                const double c = (f2_up - f2_um) / (2.0 * du);
                const double b = (f1_vp - f1_vm) / (2.0 * dv);
                const double d = (f2_vp - f2_vm) / (2.0 * dv);

                const double det = a * d - b * c;
                if (!std::isfinite(det) || std::abs(det) < 1e-18) {
                    out.iterations = it;
                    out.message = "Newton: singular Jacobian";
                    break;
                }

                double delta_u = (b * f2 - d * f1) / det;
                double delta_v = (c * f1 - a * f2) / det;

                constexpr double max_step = 0.5;
                if (std::abs(delta_u) > max_step) delta_u = (delta_u > 0.0 ? 1.0 : -1.0) * max_step;
                if (std::abs(delta_v) > max_step) delta_v = (delta_v > 0.0 ? 1.0 : -1.0) * max_step;

                const double cur_n2 = n2;
                double alpha = 1.0;
                bool accepted = false;
                for (int line_search = 0; line_search < 20; ++line_search) {
                    const double ut = u + alpha * delta_u;
                    const double vt = v + alpha * delta_v;
                    double tf1 = 0.0;
                    double tf2 = 0.0;
                    if (eval(ut, vt, tf1, tf2)) {
                        const double trial_n2 = norm2(tf1, tf2);
                        if (std::isfinite(trial_n2) && trial_n2 < cur_n2) {
                            u = ut;
                            v = vt;
                            accepted = true;
                            break;
                        }
                    }
                    alpha *= 0.5;
                }
                if (!accepted) {
                    out.iterations = it;
                    out.message = "Newton: line search failed";
                    break;
                }
            }

            if (out.converged) {
                const double Tc = std::exp(u);
                const double rho_c = std::exp(v);
                try {
                    const double Pc = eos.pressure(Tc, rho_c, composition);
                    if (!(std::isfinite(Tc) && std::isfinite(rho_c) && std::isfinite(Pc))) {
                        out.converged = false;
                        out.message = "Newton: non-finite critical state";
                    } else {
                        out.Tc = Tc;
                        out.rho = rho_c;
                        out.Pc = Pc;
                        out.Zc = Pc / (rho_c * R * Tc);
                        const Core::Units::PropertyVariable Vc(1.0 / rho_c, Core::Units::Unit::M3_PER_MOL);
                        out.Vc_m3_per_kmol = Vc.as(Core::Units::Unit::M3_PER_KMOL);
                        out.method = "newton";
                        out.message = "Converged (Newton)";
                        return out;
                    }
                } catch (...) {
                    out.converged = false;
                    out.message = "Newton: failed to evaluate Pc";
                }
            }
        } else {
            out.message = "Newton: failed initial evaluation";
        }
    }

    // ---------------------------------------------------------------------
    // Attempt 2: spinodal-bracketing on T (robust for non-cubic models)
    // ---------------------------------------------------------------------
    auto findSpinodalRoots = [&](double T, double rho_min, double rho_max, int n_scan) -> std::vector<double> {
        std::vector<double> roots;
        if (!(std::isfinite(T) && T > 0.0)) return roots;
        if (!(std::isfinite(rho_min) && rho_min > 0.0)) return roots;
        if (!(std::isfinite(rho_max) && rho_max > rho_min)) return roots;
        n_scan = std::max(10, n_scan);

        auto g = [&](double rho) -> double {
            double f1 = 0.0;
            double f2 = 0.0;
            if (!evalTRho(T, rho, f1, f2)) return nan;
            return f1;
        };

        const double log_min = std::log(rho_min);
        const double log_max = std::log(rho_max);

        double rho_prev = std::exp(log_min);
        double g_prev = g(rho_prev);
        for (int i = 1; i <= n_scan; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(n_scan);
            const double rho = std::exp(log_min + t * (log_max - log_min));
            const double g_cur = g(rho);

            if (std::isfinite(g_prev) && std::isfinite(g_cur)) {
                if (g_prev == 0.0) {
                    roots.push_back(rho_prev);
                } else if (g_cur == 0.0) {
                    roots.push_back(rho);
                } else if (g_prev * g_cur < 0.0) {
                    double a = rho_prev;
                    double b = rho;
                    double fa = g_prev;
                    for (int it = 0; it < 80; ++it) {
                        const double mid = 0.5 * (a + b);
                        const double fm = g(mid);
                        if (!std::isfinite(fm)) break;
                        if (std::abs(fm) < 1e-12) {
                            a = mid;
                            b = mid;
                            break;
                        }
                        if (fa * fm < 0.0) {
                            b = mid;
                        } else {
                            a = mid;
                            fa = fm;
                        }
                    }
                    roots.push_back(0.5 * (a + b));
                }
            }

            rho_prev = rho;
            g_prev = g_cur;
        }

        std::sort(roots.begin(), roots.end());
        roots.erase(std::unique(roots.begin(), roots.end(), [](double a, double b) {
            return std::abs(a - b) <= 1e-10 * std::max(1.0, std::max(std::abs(a), std::abs(b)));
        }), roots.end());
        return roots;
    };

    const double rho_ref = rho0;
    const double rho_min = std::max(rho_ref * 0.02, 1e-6);
    const double rho_max = std::max(rho_ref * 50.0, rho_min * 10.0);

    auto hasTwoSpinodals = [&](double T) -> bool {
        const auto roots = findSpinodalRoots(T, rho_min, rho_max, 240);
        return roots.size() >= 2;
    };

    const double T_ref = T0;
    double T_low = T_ref * 0.95;
    double T_high = T_ref * 1.05;

    std::vector<double> roots_low;
    bool ok_low = false;
    bool ok_high = false;

    for (int i = 0; i < 30; ++i) {
        roots_low = findSpinodalRoots(T_low, rho_min, rho_max, 240);
        if (roots_low.size() >= 2) {
            ok_low = true;
            break;
        }
        T_low *= 0.90;
        if (T_low < 1.0) break;
    }

    for (int i = 0; i < 30; ++i) {
        if (!hasTwoSpinodals(T_high)) {
            ok_high = true;
            break;
        }
        T_high *= 1.10;
    }

    if (!(ok_low && ok_high && T_high > T_low)) {
        if (out.message.empty()) out.message = "Failed to bracket critical temperature";
        return out;
    }

    for (int it = 0; it < 60; ++it) {
        const double T_mid = 0.5 * (T_low + T_high);
        const auto roots = findSpinodalRoots(T_mid, rho_min, rho_max, 240);
        if (roots.size() >= 2) {
            T_low = T_mid;
            roots_low = roots;
        } else {
            T_high = T_mid;
        }
        if ((T_high - T_low) < 1e-12 * std::max(1.0, T_low)) break;
    }

    if (roots_low.size() < 2) {
        if (out.message.empty()) out.message = "Failed to find spinodals near Tc";
        return out;
    }

    const double rho1 = roots_low.front();
    const double rho2 = roots_low.back();
    const double Tc = T_low;
    const double rho_c = 0.5 * (rho1 + rho2);

    try {
        const double Pc = eos.pressure(Tc, rho_c, composition);
        if (!(std::isfinite(Tc) && std::isfinite(rho_c) && std::isfinite(Pc))) {
            out.converged = false;
            out.message = "Spinodal: non-finite critical state";
            return out;
        }
        out.converged = true;
        out.iterations = 0;
        out.Tc = Tc;
        out.rho = rho_c;
        out.Pc = Pc;
        out.Zc = Pc / (rho_c * R * Tc);
        const Core::Units::PropertyVariable Vc(1.0 / rho_c, Core::Units::Unit::M3_PER_MOL);
        out.Vc_m3_per_kmol = Vc.as(Core::Units::Unit::M3_PER_KMOL);
        out.method = "spinodal";
        out.message = "Converged (spinodal)";
        return out;
    } catch (...) {
        out.message = "Spinodal: failed to evaluate Pc";
        out.converged = false;
        return out;
    }
}

} // namespace Critical
} // namespace Equilibrium
} // namespace DMThermo
