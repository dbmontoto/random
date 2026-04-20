/**
 * @file reactive_flash.h
 * @brief Reactive flash solver (TP and TV) using generalized Gibbs/Helmholtz optimization backends.
 */

#ifndef THERMO_REACTIONS_REACTIVE_FLASH_H
#define THERMO_REACTIONS_REACTIVE_FLASH_H

#include "thermo/data/databanks.h"
#include "thermo/data/diagnostics.h"
#include "thermo/eos.h"
#include "thermo/equilibrium/flash_result.h"
#include "thermo/equilibrium/istability.h"
#include "thermo/config/equilibrium_optimizer_config.h"
#include "reaction_system.h"

#include <functional>
#include <string>
#include <vector>

namespace DMThermo {
namespace Reactions {

struct ReactiveFlashConfig {
    Config::EquilibriumOptimizerConfig equilibrium = Config::EquilibriumOptimizerConfig::defaults();
    double standard_pressure = 101325.0; // Pa
    double min_moles = 1e-30;

    // Databank thermochemistry policy (used only by *FromDatabanks methods).
    std::string mu0_model = "DIPPR.GFOR";
    double mu0_reference_temperature = 298.15;
};

struct ReactivePHFlashConfig {
    ReactiveFlashConfig base{};
    double enthalpy_tolerance = 100.0;  // J/mol
    double temp_bracket_low = 100.0;    // K
    double temp_bracket_high = 1000.0;  // K
    int max_iterations = 80;            // outer temperature iterations
};

struct ReactivePSFlashConfig {
    ReactiveFlashConfig base{};
    double entropy_tolerance = 0.1;     // J/(mol*K)
    double temp_bracket_low = 100.0;    // K
    double temp_bracket_high = 1000.0;  // K
    int max_iterations = 80;            // outer temperature iterations
};

struct ReactiveFlashResult {
    bool converged = false;
    int iterations = 0;
    std::string message;

    double temperature = 0.0;
    double pressure = 0.0;

    std::vector<double> n0;
    std::vector<double> n;
    std::vector<double> z;
    std::vector<double> extents;

    Equilibrium::Flash::FlashResult flash;
};

/**
 * @brief Reactive flash solver that couples phase equilibrium with reaction equilibrium.
 *
 * - TP (reactive Gibbs): minimizes total Gibbs energy at (T,P) with reaction extents and phase allocations.
 * - TV (reactive Helmholtz strategy): solves outer log(P) for the target molar volume V using Brent, with
 *   an inner reactive TP minimization at each P.
 *
 * Notes:
 * - Uses EOS fugacity coefficients only (EOS-agnostic).
 * - If you want sequential "react then flash", use the reaction equilibrium solver first and then call a flash solver.
 */
class ReactiveFlashSolver final {
public:
    explicit ReactiveFlashSolver(
        EOSPtr eos,
        Equilibrium::Stability::StabilityAnalyzerPtr stability = nullptr);

    using Mu0Provider = std::function<std::vector<double>(double T)>;

    ReactiveFlashResult solveTP(
        double T,
        double P,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const std::vector<double>& mu0, // J/mol at cfg.standard_pressure
        const ReactiveFlashConfig& cfg = ReactiveFlashConfig{}) const;

    ReactiveFlashResult solveTPFromDatabanks(
        double T,
        double P,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const Data::Databanks& databanks,
        const Core::Mixture& mixture,
        const ReactiveFlashConfig& cfg = ReactiveFlashConfig{},
        Data::Diagnostics* diag = nullptr) const;

    ReactiveFlashResult solveTV(
        double T,
        double V,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const std::vector<double>& mu0, // J/mol at cfg.standard_pressure
        const ReactiveFlashConfig& cfg = ReactiveFlashConfig{}) const;

    ReactiveFlashResult solveTVFromDatabanks(
        double T,
        double V,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const Data::Databanks& databanks,
        const Core::Mixture& mixture,
        const ReactiveFlashConfig& cfg = ReactiveFlashConfig{},
        Data::Diagnostics* diag = nullptr) const;

    ReactiveFlashResult solvePH(
        double P,
        double H,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const Mu0Provider& mu0AtT,
        const Core::Mixture& mixture,
        const ReactivePHFlashConfig& cfg = ReactivePHFlashConfig{}) const;

    ReactiveFlashResult solvePHFromDatabanks(
        double P,
        double H,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const Data::Databanks& databanks,
        const Core::Mixture& mixture,
        const ReactivePHFlashConfig& cfg = ReactivePHFlashConfig{},
        Data::Diagnostics* diag = nullptr) const;

    ReactiveFlashResult solvePS(
        double P,
        double S,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const Mu0Provider& mu0AtT,
        const Core::Mixture& mixture,
        const ReactivePSFlashConfig& cfg = ReactivePSFlashConfig{}) const;

    ReactiveFlashResult solvePSFromDatabanks(
        double P,
        double S,
        const std::vector<double>& n0,
        const ReactionSystem& system,
        const Data::Databanks& databanks,
        const Core::Mixture& mixture,
        const ReactivePSFlashConfig& cfg = ReactivePSFlashConfig{},
        Data::Diagnostics* diag = nullptr) const;

private:
    EOSPtr eos_;
    Equilibrium::Stability::StabilityAnalyzerPtr stability_;
};

} // namespace Reactions
} // namespace DMThermo

#endif // THERMO_REACTIONS_REACTIVE_FLASH_H
