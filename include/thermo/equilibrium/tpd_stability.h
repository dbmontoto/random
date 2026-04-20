/**
 * @file tpd_stability.h
 * @brief EOS-based TPD stability analyzer (works with any DMThermo::EOS)
 */

#ifndef THERMO_EQUILIBRIUM_TPD_STABILITY_H
#define THERMO_EQUILIBRIUM_TPD_STABILITY_H

#include "istability.h"
#include "thermo/eos.h"
#include <memory>

namespace DMThermo {
namespace Equilibrium {
namespace Stability {

/**
 * @brief Tangent Plane Distance (TPD) stability analyzer using only the EOS interface
 *
 * This implementation supports cubic EOS baselines (PR/SRK) and any other EOS
 * that can provide pressure and fugacity coefficients at a given (T, rho, x).
 */
class TPDStabilityAnalyzer : public IStabilityAnalyzer {
public:
    explicit TPDStabilityAnalyzer(EOSPtr eos);

    StabilityResult analyze(
        double T,
        double P,
        const std::vector<double>& z,
        const Config::StabilityConfig& config = Config::StabilityConfig::defaults()
    ) const override;

    bool isStable(double T, double P, const std::vector<double>& z) const override;

    double tangentPlaneDistance(
        double T,
        double P,
        const std::vector<double>& z,
        const std::vector<double>& w
    ) const override;

    TPDTrialResult minimizeTPD(
        double T,
        double P,
        const std::vector<double>& z,
        const std::vector<double>& w_initial,
        const Config::StabilityConfig& config = Config::StabilityConfig::defaults()
    ) const override;

    DensityRootResult findDensityRoots(
        double T,
        double P,
        const std::vector<double>& x,
        const Config::DensityRootConfig& config = Config::DensityRootConfig{}
    ) const override;

    double findVaporDensity(double T, double P, const std::vector<double>& x) const override;
    double findLiquidDensity(double T, double P, const std::vector<double>& x) const override;

    SpinodalResult findSpinodalPoints(
        double T,
        const std::vector<double>& x,
        const Config::SpinodalConfig& config = Config::SpinodalConfig{}
    ) const override;

    bool isInSpinodalRegion(double T, double rho, const std::vector<double>& x) const override;

    PhaseType classifyPhase(double T, double rho, const std::vector<double>& x) const override;

    std::vector<std::vector<double>> generateTrialCompositions(
        double T,
        double P,
        const std::vector<double>& z,
        int num_trials
    ) const override;

    EOSPtr eos() const override { return eos_; }

private:
    EOSPtr eos_;

    struct LnPhiAtTP {
        double rho = 0.0;
        double Z = 0.0;
        std::vector<double> ln_phi;
    };

    LnPhiAtTP lnPhiAtTP(
        double T,
        double P,
        const std::vector<double>& x,
        PhaseType phase_hint,
        const Config::StabilityConfig& config,
        double rho_guess) const;

    static std::vector<double> normalize(const std::vector<double>& x);
    static std::vector<double> safeLog(const std::vector<double>& x);
    static double tpdFromLnPhi(
        const std::vector<double>& z,
        const std::vector<double>& lnphi_z,
        const std::vector<double>& w,
        const std::vector<double>& lnphi_w);

    static std::vector<double> wilsonK(const EOS& eos, double T, double P);

    TPDTrialResult minimizeTPDWithPhaseHint(
        double T,
        double P,
        const std::vector<double>& z,
        const std::vector<double>& lnphi_z,
        double tpd_scale,
        const std::vector<double>& w_initial,
        PhaseType phase_hint,
        const Config::StabilityConfig& config) const;
};

} // namespace Stability
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_TPD_STABILITY_H
