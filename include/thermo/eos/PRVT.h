/**
 * @file PRVT.h
 * @brief Peng-Robinson EOS with Peneloux-style constant volume translation (PR-VT)
 *
 * This implementation applies a constant volume translation to Peng-Robinson:
 *   v = v_PR + c(x), where c(x) = Σ_i x_i c_i.
 *
 * The translation constants are computed per component from critical properties:
 *   c_i = (Zc_i - Zc_PR) * (R Tc_i / Pc_i)
 * with Zc_PR = 0.307 for Peng-Robinson.
 *
 * If Zc is missing (NaN) for a component, c_i defaults to 0 (no translation).
 */

#ifndef THERMO_EOS_PRVT_H
#define THERMO_EOS_PRVT_H

#include "thermo/eos.h"
#include "thermo/eos/PengRobinson.h"
#include "thermo/core/mixture.h"
#include <vector>

namespace DMThermo {
namespace Cubic {

/**
 * @brief PR-VT EOS (Peng-Robinson with constant volume translation)
 */
class PRVTEOS : public ::DMThermo::EOS {
public:
    explicit PRVTEOS(const Core::Mixture& mixture);

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

    std::string name() const override { return "PR-VT"; }
    int numComponents() const override { return nc_; }
    bool hasAssociation() const override { return false; }

    std::vector<double> densityRootsTP(double T, double P, const std::vector<double>& x) const override;

    /// Per-component translation constant c_i [m^3/mol] (0.0 when unavailable).
    double volumeTranslation(int i) const;

    double criticalTemperature(int i) const override;
    double criticalPressure(int i) const override;
    double acentricFactor(int i) const override;
    const Core::Mixture* coreMixture() const override { return &mixture_; }

    /// Access underlying mixture (for property/evaluation tooling)
    const Core::Mixture& mixture() const { return mixture_; }

private:
    Core::Mixture mixture_;
    int nc_ = 0;
    PengRobinsonEOS pr_;
    std::vector<double> c_i_;

    double cMix(const std::vector<double>& x) const;
    double rhoPrimeFromRho(double rho, double c_mix) const;
};

} // namespace Cubic
} // namespace DMThermo

#endif // THERMO_EOS_PRVT_H

