/**
 * @file VTPR.h
 * @brief VTPR EOS: Peng-Robinson + Wong-Sandler mixing rule + UNIFAC + volume translation
 */

#ifndef THERMO_EOS_VTPR_H
#define THERMO_EOS_VTPR_H

#include "thermo/eos.h"
#include "thermo/core/constants.h"
#include "thermo/core/mixture.h"

#include <string>
#include <vector>

namespace DMThermo {
namespace Cubic {

/**
 * @brief VTPR (Ahlers/Gmehling style): PR with Wong-Sandler mixing using UNIFAC g^E.
 *
 * Requirements (provided via databank-built Core::Mixture):
 * - Tc/Pc/omega (DIPPR) for PR pure-component a_i(T), b_i.
 * - UNIFAC subgroup decomposition per component (vtpr_groups.csv).
 * - UNIFAC subgroup parameters + main-group interactions (unifac_*.csv).
 * - Temperature-dependent volume translation c_i(T) in m^3/mol (vtpr_pure.csv).
 * - Optional binary interaction k_ij from bip.csv column `kij_vtpr` (defaults to 0).
 */
class VTPREOS : public ::DMThermo::EOS {
public:
    explicit VTPREOS(const Core::Mixture& mixture);

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

    std::string name() const override { return "VTPR"; }
    int numComponents() const override { return nc_; }
    bool hasAssociation() const override { return false; }

    double criticalTemperature(int i) const override;
    double criticalPressure(int i) const override;
    double acentricFactor(int i) const override;
    const Core::Mixture* coreMixture() const override { return &mixture_; }

    std::vector<double> densityRootsTP(double T, double P, const std::vector<double>& x) const override;

private:
    Core::Mixture mixture_;
    int nc_ = 0;
    std::vector<double> Tc_;
    std::vector<double> Pc_;
    std::vector<double> omega_;

    struct WSParams {
        std::vector<double> a_i;     // PR a_i(T)
        std::vector<double> b_i;     // PR b_i
        std::vector<double> ln_gamma;
        double Q = 0.0;
        double D = 0.0;
        double a = 0.0;
        double b = 0.0;
        double c_mix = 0.0;          // mixture volume translation at T
    };

    void initializeCriticalCache();

    WSParams wsParams(double T, const std::vector<double>& x, bool need_gamma) const;

    double rhoEosFromRho(double rho, double c_mix) const;
    double rhoFromRhoEos(double rho_eos, double c_mix) const;

    double pressureFromMix(double T, double v, double a, double b) const;

    std::vector<double> lnPhi(double T, double rho_eos, const std::vector<double>& x, double& Z_out) const;
};

} // namespace Cubic
} // namespace DMThermo

#endif // THERMO_EOS_VTPR_H

