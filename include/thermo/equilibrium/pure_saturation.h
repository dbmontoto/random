/**
 * @file pure_saturation.h
 * @brief Generic pure-component saturation solver using EOS + density roots
 */

#ifndef THERMO_EQUILIBRIUM_PURE_SATURATION_H
#define THERMO_EQUILIBRIUM_PURE_SATURATION_H

#include "thermo/eos.h"
#include "istability.h"
#include <memory>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Saturation {

struct PureSaturationPoint {
    double T = 0.0;
    double P_sat = 0.0;
    double rho_v = 0.0;
    double rho_l = 0.0;
    bool converged = false;
    std::string message;
};

struct PureSaturationConfig {
    double P_min = 1.0;              ///< Minimum pressure search bound [Pa]
    double P_max = 0.0;              ///< Maximum pressure search bound [Pa] (0 => use 0.999*Pc)
    double P_guess = 0.0;            ///< Optional initial guess [Pa] (0 => none)
    double guess_span_decades = 6.0; ///< If P_guess provided, scan within P_guess*10^(±span/2)
    int max_scan_points = 80;        ///< Pressure scan points to bracket saturation
    int max_iterations = 80;         ///< Bisection iterations
    double tolerance = 1e-10;        ///< |ln(phi_l)-ln(phi_v)| target
    double relP_tolerance = 1e-10;   ///< Relative pressure tolerance
};

/**
 * @brief Pure-component saturation solver (Psat, rhoV, rhoL vs T)
 *
 * Works with any DMThermo::EOS, using the attached StabilityAnalyzer to locate
 * vapor/liquid density roots at a given (T, P).
 */
class PureSaturationSolver {
public:
    PureSaturationSolver(EOSPtr eos, Stability::StabilityAnalyzerPtr stability);

    PureSaturationPoint solveAtT(
        double T,
        const PureSaturationConfig& config = PureSaturationConfig{}
    ) const;

    std::vector<PureSaturationPoint> solveCurve(
        double T_min,
        double T_max,
        int n_points,
        const PureSaturationConfig& config = PureSaturationConfig{}
    ) const;

private:
    EOSPtr eos_;
    Stability::StabilityAnalyzerPtr stability_;

    struct PhiAtTP {
        double lnphi_v = 0.0;
        double lnphi_l = 0.0;
        double rho_v = 0.0;
        double rho_l = 0.0;
        bool valid = false;
        std::string message;
    };

    PhiAtTP phiAtTP(double T, double P) const;
};

} // namespace Saturation
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_PURE_SATURATION_H
