/**
 * @file phase_detection.h
 * @brief Shared phase-detection / density-selection configuration types.
 */

#ifndef THERMO_CONFIG_PHASE_DETECTION_H
#define THERMO_CONFIG_PHASE_DETECTION_H

namespace DMThermo {
namespace Config {

/**
 * @brief Primary state space used for phase detection / density selection.
 *
 * - TP: determine densities by solving `P(T,rho,x)=P` (scan/bracket + root refine). This is the most
 *   robust "cold start" and naturally supports multiple density roots.
 * - TRho: track/refine densities using Newton steps in `rho` to enforce pressure, optionally using
 *   a previous `rho_guess` to stay on the same phase branch during iterative solvers.
 *
 * Guidance:
 * - Cubic EOS (PR/SRK): typically use `TP` (cheap evaluations, robust bracketing).
 * - Density-based EOS (PC-SAFT, translated models): consider `TRho` inside iterative solvers (TV/Helmholtz,
 *   multiphase optimization) for performance and branch consistency; use `TP` for cold-start robustness.
 *
 * See `docs/implementation/PHASE_DETECTION_POLICY.md` for details.
 */
enum class PhaseDetection {
    TP,
    TRho
};

} // namespace Config
} // namespace DMThermo

#endif // THERMO_CONFIG_PHASE_DETECTION_H
