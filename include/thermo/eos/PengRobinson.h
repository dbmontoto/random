/**
 * @file PengRobinson.h
 * @brief Peng-Robinson EOS implementation (cubic)
 *
 * Provides a modular implementation of the Peng-Robinson (PR) EOS through the
 * DMThermo::EOS interface for baseline comparisons against PC-SAFT.
 */

#ifndef THERMO_EOS_PENG_ROBINSON_H
#define THERMO_EOS_PENG_ROBINSON_H

#include "thermo/eos.h"
#include "thermo/core/constants.h"
#include "thermo/core/mixture.h"
#include <vector>
#include <string>

namespace DMThermo {
namespace Cubic {

/**
 * @brief Peng-Robinson equation of state (mixtures supported via vdW mixing rules)
 *
 * Mixing rules:
 * - b = Σ x_i b_i
 * - a = Σ_i Σ_j x_i x_j a_ij, a_ij = sqrt(a_i a_j) (1 - k_ij)
 *
 * Component parameters use Tc, Pc, omega (provided by databanks).
 */
class PengRobinsonEOS : public ::DMThermo::EOS {
public:
    explicit PengRobinsonEOS(const Core::Mixture& mixture);

    EOSResult calculate(
        double T,
        double rho,
        const std::vector<double>& x,
        const DerivativeSpec& deriv_spec = DerivativeSpec::none()
    ) const override;

    double pressure(double T, double rho, const std::vector<double>& x) const override;
    std::vector<double> fugacityCoefficients(double T, double rho, const std::vector<double>& x) const override;
    double compressibility(double T, double rho, const std::vector<double>& x) const override;
    double residualHelmholtz(double T, double rho, const std::vector<double>& x) const override;

    double dadt(double T, double rho, const std::vector<double>& x) const override;
    double dadrho(double T, double rho, const std::vector<double>& x) const override;
    double dPdrho(double T, double rho, const std::vector<double>& x) const override;

    std::string name() const override { return "Peng-Robinson"; }
    int numComponents() const override { return nc_; }
    bool hasAssociation() const override { return false; }

    double criticalTemperature(int i) const override;
    double criticalPressure(int i) const override;
    double acentricFactor(int i) const override;
    const Core::Mixture* coreMixture() const override { return &mixture_; }

    /// Access underlying mixture (for property/evaluation tooling)
    const Core::Mixture& mixture() const { return mixture_; }

    /**
     * @brief Real compressibility-factor roots at (T, P, x) from the PR cubic in Z.
     *
     * Returns unique real roots (sorted ascending) that satisfy the physical constraint Z > B.
     */
    std::vector<double> compressibilityRootsTP(double T, double P, const std::vector<double>& x) const;

    /**
     * @brief Real density roots at (T, P, x) from the PR cubic in Z.
     *
     * Returns unique real roots (sorted ascending) expressed as molar densities [mol/m^3].
     */
    std::vector<double> densityRootsTP(double T, double P, const std::vector<double>& x) const override;

private:
    Core::Mixture mixture_;
    int nc_ = 0;
    std::vector<double> Tc_;
    std::vector<double> Pc_;
    std::vector<double> omega_;

    struct MixParams {
        std::vector<double> a_i;
        std::vector<double> b_i;
        std::vector<double> a_i_mix; // a_i_mix = Σ_j x_j a_ij
        double a_mix = 0.0;
        double b_mix = 0.0;
        std::vector<double> da_i_dT;
        double da_mix_dT = 0.0;
        std::vector<double> d2a_i_dT2;
        double d2a_mix_dT2 = 0.0;
    };

    void initializeCriticalCache();

    MixParams mixtureParams(double T, const std::vector<double>& x, bool compute_temperature_derivatives = false) const;
    double pressureFromMix(double T, double v, double a, double b) const;

    std::vector<double> lnPhi(double T, double rho, const std::vector<double>& x, double& Z_out) const;
    double residualHelmholtzFromLnPhi(double Z, const std::vector<double>& x, const std::vector<double>& ln_phi) const;

};

} // namespace Cubic
} // namespace DMThermo

#endif // THERMO_EOS_PENG_ROBINSON_H
