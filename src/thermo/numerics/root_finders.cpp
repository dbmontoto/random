/**
 * @file root_finders.cpp
 * @brief Root finder implementations for the numerics interfaces.
 */

#include "thermo/numerics/iroot_finder.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DMThermo {
namespace Numerics {
namespace RootFinding {
namespace {

class NewtonRaphson final : public IScalarRootFinder {
public:
    ScalarRootResult solve(
        std::function<double(double)> f,
        double x0,
        const Config::RootFinderConfig& config) const override
    {
        ScalarRootResult out;
        double x = x0;
        double fx = f(x);

        for (int it = 0; it < config.max_iterations; ++it) {
            fx = f(x);
            if (!std::isfinite(fx)) {
                out.converged = false;
                out.iterations = it + 1;
                out.root = x;
                out.residual = std::numeric_limits<double>::infinity();
                out.message = "Non-finite f(x)";
                return out;
            }
            if (std::abs(fx) < config.tolerance) {
                out.converged = true;
                out.iterations = it + 1;
                out.root = x;
                out.residual = std::abs(fx);
                out.message = "Converged";
                return out;
            }

            const double h = 1e-7 * std::max(1.0, std::abs(x));
            const double dfx = (f(x + h) - f(x - h)) / (2.0 * h);
            if (!std::isfinite(dfx) || std::abs(dfx) < 1e-18) {
                break;
            }

            double dx = fx / dfx;
            dx *= config.damping_factor;
            if (!std::isfinite(dx)) break;

            if (config.use_line_search) {
                double alpha = 1.0;
                const double c1 = config.line_search_alpha;
                const double beta = config.line_search_beta;
                const double f_abs = std::abs(fx);
                for (int ls = 0; ls < 20; ++ls) {
                    const double x_try = x - alpha * dx;
                    const double f_try = f(x_try);
                    if (std::isfinite(f_try) && std::abs(f_try) < f_abs * (1.0 - c1 * alpha)) {
                        x = x_try;
                        break;
                    }
                    alpha *= beta;
                    if (alpha < 1e-12) {
                        x = x - alpha * dx;
                        break;
                    }
                }
            } else {
                x -= dx;
            }

            if (std::abs(dx) < config.tolerance) {
                out.converged = true;
                out.iterations = it + 1;
                out.root = x;
                out.residual = std::abs(f(x));
                out.message = "Converged";
                return out;
            }
        }

        out.converged = false;
        out.iterations = config.max_iterations;
        out.root = x;
        out.residual = std::abs(fx);
        out.message = "Maximum iterations reached";
        return out;
    }

    ScalarRootResult solveBracketed(
        std::function<double(double)> f,
        double a,
        double b,
        const Config::RootFinderConfig& config) const override
    {
        // Bisection fallback (robust).
        ScalarRootResult out;
        double fa = f(a);
        double fb = f(b);
        if (!(std::isfinite(fa) && std::isfinite(fb)) || fa * fb > 0.0) {
            out.converged = false;
            out.root = 0.5 * (a + b);
            out.iterations = 0;
            out.message = "Invalid bracket";
            return out;
        }

        double lo = a, hi = b;
        double flo = fa, fhi = fb;
        for (int it = 0; it < config.max_iterations; ++it) {
            const double mid = 0.5 * (lo + hi);
            const double fmid = f(mid);
            if (!std::isfinite(fmid)) break;
            if (std::abs(fmid) < config.tolerance ||
                std::abs(hi - lo) < config.rel_tolerance * std::max(1.0, std::abs(mid))) {
                out.converged = true;
                out.iterations = it + 1;
                out.root = mid;
                out.residual = std::abs(fmid);
                out.message = "Converged";
                return out;
            }
            if (flo * fmid < 0.0) {
                hi = mid; fhi = fmid;
            } else {
                lo = mid; flo = fmid;
            }
        }

        out.converged = false;
        out.iterations = config.max_iterations;
        out.root = 0.5 * (lo + hi);
        out.residual = std::abs(f(out.root));
        out.message = "Maximum iterations reached";
        return out;
    }

    std::string name() const override { return "NewtonRaphson"; }
};

class Brent final : public IScalarRootFinder {
public:
    ScalarRootResult solve(
        std::function<double(double)> f,
        double x0,
        const Config::RootFinderConfig& config) const override
    {
        // No bracket given: try to bracket around x0, then call bracketed.
        double a = x0, b = x0;
        double fa = f(a), fb = fa;
        double step = 0.1 * std::max(1.0, std::abs(x0));
        bool bracketed = false;
        for (int i = 0; i < 40; ++i) {
            a = x0 - step;
            b = x0 + step;
            fa = f(a);
            fb = f(b);
            if (std::isfinite(fa) && std::isfinite(fb) && fa * fb < 0.0) {
                bracketed = true;
                break;
            }
            step *= 1.6;
        }
        if (!bracketed) {
            ScalarRootResult out;
            out.converged = false;
            out.root = x0;
            out.message = "Failed to bracket root";
            return out;
        }
        return solveBracketed(std::move(f), a, b, config);
    }

    ScalarRootResult solveBracketed(
        std::function<double(double)> f,
        double a,
        double b,
        const Config::RootFinderConfig& config) const override
    {
        // Brent implementation (simplified, robust).
        ScalarRootResult out;
        double fa = f(a);
        double fb = f(b);
        if (!(std::isfinite(fa) && std::isfinite(fb)) || fa * fb > 0.0) {
            out.converged = false;
            out.root = 0.5 * (a + b);
            out.message = "Invalid bracket";
            return out;
        }
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }

        double c = a;
        double fc = fa;
        bool mflag = true;
        double d = 0.0;

        for (int it = 0; it < config.max_iterations; ++it) {
            if (std::abs(fb) < config.tolerance ||
                std::abs(b - a) < config.rel_tolerance * std::max(1.0, std::abs(b))) {
                out.converged = true;
                out.iterations = it + 1;
                out.root = b;
                out.residual = std::abs(fb);
                out.message = "Converged";
                return out;
            }

            double s = b;
            if (fa != fc && fb != fc) {
                s = a * fb * fc / ((fa - fb) * (fa - fc))
                  + b * fa * fc / ((fb - fa) * (fb - fc))
                  + c * fa * fb / ((fc - fa) * (fc - fb));
            } else {
                s = b - fb * (b - a) / (fb - fa);
            }

            const bool cond1 = (s < (3.0 * a + b) / 4.0 || s > b);
            const bool cond2 = mflag && (std::abs(s - b) >= std::abs(b - c) / 2.0);
            const bool cond3 = !mflag && (std::abs(s - b) >= std::abs(c - d) / 2.0);
            const bool cond4 = mflag && (std::abs(b - c) < config.tolerance);
            const bool cond5 = !mflag && (std::abs(c - d) < config.tolerance);
            if (cond1 || cond2 || cond3 || cond4 || cond5) {
                s = 0.5 * (a + b);
                mflag = true;
            } else {
                mflag = false;
            }

            double fs = f(s);
            if (!std::isfinite(fs)) {
                // Try safer bisection-style fallbacks within the bracket.
                double s_try = 0.5 * (a + b);
                double f_try = f(s_try);
                if (std::isfinite(f_try)) {
                    s = s_try;
                    fs = f_try;
                } else {
                    s_try = 0.5 * (s_try + b);
                    f_try = f(s_try);
                    if (std::isfinite(f_try)) {
                        s = s_try;
                        fs = f_try;
                    } else {
                        s_try = 0.5 * (a + s_try);
                        f_try = f(s_try);
                        if (std::isfinite(f_try)) {
                            s = s_try;
                            fs = f_try;
                        } else {
                            out.converged = false;
                            out.iterations = it + 1;
                            out.root = b;
                            out.residual = std::numeric_limits<double>::infinity();
                            out.message = "Non-finite f(x) within bracket";
                            return out;
                        }
                    }
                }
            }
            d = c;
            c = b;
            fc = fb;

            if (fa * fs < 0.0) {
                b = s;
                fb = fs;
            } else {
                a = s;
                fa = fs;
            }

            if (std::abs(fa) < std::abs(fb)) {
                std::swap(a, b);
                std::swap(fa, fb);
            }
        }

        out.converged = false;
        out.iterations = config.max_iterations;
        out.root = b;
        out.residual = std::abs(fb);
        out.message = "Maximum iterations reached";
        return out;
    }

    std::string name() const override { return "Brent"; }
};

class Secant final : public IScalarRootFinder {
public:
    ScalarRootResult solve(
        std::function<double(double)> f,
        double x0,
        const Config::RootFinderConfig& config) const override
    {
        ScalarRootResult out;
        double x1 = x0 + 1e-2 * std::max(1.0, std::abs(x0));
        double f0 = f(x0);
        double f1 = f(x1);
        for (int it = 0; it < config.max_iterations; ++it) {
            if (std::abs(f1) < config.tolerance) {
                out.converged = true;
                out.iterations = it + 1;
                out.root = x1;
                out.residual = std::abs(f1);
                out.message = "Converged";
                return out;
            }
            const double denom = (f1 - f0);
            if (!std::isfinite(denom) || std::abs(denom) < 1e-18) break;
            const double x2 = x1 - f1 * (x1 - x0) / denom;
            x0 = x1; f0 = f1;
            x1 = x2; f1 = f(x1);
        }
        out.converged = false;
        out.iterations = config.max_iterations;
        out.root = x1;
        out.residual = std::abs(f1);
        out.message = "Maximum iterations reached";
        return out;
    }

    ScalarRootResult solveBracketed(
        std::function<double(double)> f,
        double a,
        double b,
        const Config::RootFinderConfig& config) const override
    {
        // Bracketed secant with bisection fallback.
        ScalarRootResult out;
        double fa = f(a);
        double fb = f(b);
        if (!(std::isfinite(fa) && std::isfinite(fb)) || fa * fb > 0.0) {
            out.converged = false;
            out.root = 0.5 * (a + b);
            out.message = "Invalid bracket";
            return out;
        }
        double x0 = a, x1 = b;
        for (int it = 0; it < config.max_iterations; ++it) {
            if (std::abs(fb) < config.tolerance) {
                out.converged = true;
                out.iterations = it + 1;
                out.root = x1;
                out.residual = std::abs(fb);
                out.message = "Converged";
                return out;
            }
            double x2 = x1;
            const double denom = (fb - fa);
            if (std::abs(denom) > 1e-18) {
                x2 = x1 - fb * (x1 - x0) / denom;
            } else {
                x2 = 0.5 * (x0 + x1);
            }
            const double f2 = f(x2);
            if (!std::isfinite(f2)) break;
            if (fa * f2 < 0.0) {
                x1 = x2; fb = f2;
            } else {
                x0 = x2; fa = f2;
            }
            if (std::abs(x1 - x0) < config.rel_tolerance * std::max(1.0, std::abs(x1))) break;
        }
        out.converged = false;
        out.iterations = config.max_iterations;
        out.root = x1;
        out.residual = std::abs(fb);
        out.message = "Maximum iterations reached";
        return out;
    }

    std::string name() const override { return "Secant"; }
};

} // namespace

std::shared_ptr<IScalarRootFinder> makeNewtonRaphson() { return std::make_shared<NewtonRaphson>(); }
std::shared_ptr<IScalarRootFinder> makeBrent() { return std::make_shared<Brent>(); }
std::shared_ptr<IScalarRootFinder> makeSecant() { return std::make_shared<Secant>(); }

} // namespace RootFinding
} // namespace Numerics
} // namespace DMThermo
