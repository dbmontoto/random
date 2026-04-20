/**
 * @file PCSAFT.h
 * @brief PC-SAFT equation of state implementation
 *
 * This class implements the PC-SAFT equation of state as a concrete
 * implementation of the DMThermo::EOS interface.
 */

#ifndef THERMO_EOS_PCSAFT_H
#define THERMO_EOS_PCSAFT_H

#include "thermo/eos.h"
#include "thermo/core/mixture.h"
#include "thermo/core/constants.h"
#include <memory>

namespace DMThermo {
namespace PCSaft {

/**
 * @brief PC-SAFT equation of state implementation
 *
 * PC-SAFT equation of state exposed through the DMThermo::EOS interface.
 *
 * This type deliberately hides legacy implementation details from the public
 * header surface.
 */
class PCSaftEOS : public ::DMThermo::EOS {
public:
    /**
     * @brief Construct from new-architecture Core::Mixture
     *
     * Stores a copy of the Core::Mixture for ideal-gas correlations, and builds the
	     * legacy pcsaft::Mixture internally for the underlying model.
	     */
	    explicit PCSaftEOS(const Core::Mixture& mixture);
	    ~PCSaftEOS() override;

	    PCSaftEOS(const PCSaftEOS&) = delete;
	    PCSaftEOS& operator=(const PCSaftEOS&) = delete;
	    PCSaftEOS(PCSaftEOS&&) noexcept;
	    PCSaftEOS& operator=(PCSaftEOS&&) noexcept;

    // =========================================================================
    // EOS Interface Implementation
    // =========================================================================

    EOSResult calculate(
        double T,
        double rho,
        const std::vector<double>& x,
        const DerivativeSpec& deriv_spec = DerivativeSpec::none()
    ) const override;

    double pressure(
        double T,
        double rho,
        const std::vector<double>& x
    ) const override;

    std::vector<double> fugacityCoefficients(
        double T,
        double rho,
        const std::vector<double>& x
    ) const override;

    double compressibility(
        double T,
        double rho,
        const std::vector<double>& x
    ) const override;

    double residualHelmholtz(
        double T,
        double rho,
        const std::vector<double>& x
    ) const override;

    double dadt(
        double T,
        double rho,
        const std::vector<double>& x
    ) const override;

    double dadrho(
        double T,
        double rho,
        const std::vector<double>& x
    ) const override;

    double dPdrho(
        double T,
        double rho,
        const std::vector<double>& x
    ) const override;

    std::string name() const override { return "PC-SAFT"; }

    int numComponents() const override { return nc_; }

    bool hasAssociation() const override { return has_association_; }

    double criticalTemperature(int i) const override;
    double criticalPressure(int i) const override;
    double acentricFactor(int i) const override;
    const Core::Mixture* coreMixture() const override;

    // =========================================================================
    // PC-SAFT Specific Methods
    // =========================================================================

    /**
     * @brief Get packing fraction at given state
     * @return eta = pi/6 * rho * sum(x_i m_i d_i^3)
     */
    double packingFraction(double T, double rho, const std::vector<double>& x) const;

    /**
     * @brief Get hard-chain contribution to Helmholtz energy
     */
    double hardChainContribution(double T, double rho, const std::vector<double>& x) const;

    /**
     * @brief Get dispersion contribution to Helmholtz energy
     */
    double dispersionContribution(double T, double rho, const std::vector<double>& x) const;

    /**
     * @brief Get association contribution to Helmholtz energy
     */
    double associationContribution(double T, double rho, const std::vector<double>& x) const;

	private:
	    struct Impl;
	    std::unique_ptr<Impl> impl_;

	    /// New-architecture mixture (used for ideal-gas correlations and metadata)
	    Core::Mixture core_mixture_{std::vector<Core::Component>{}, Core::BinaryParameters::zeros(0)};

    /// Cached composition for avoiding unnecessary rebuilds
    mutable std::vector<double> cached_x_;

    /// Number of components
    int nc_;

    /// Does mixture have associating components?
    bool has_association_;

    /// Critical properties (optional; may be NaN if unavailable)
    std::vector<double> Tc_;
    std::vector<double> Pc_;
    std::vector<double> omega_;

    /**
     * @brief Create EOS with updated composition (thread-safe)
     */
    void updateComposition(const std::vector<double>& x) const;

    /**
     * @brief Initialize cached properties from mixture
     */
    void initializeCachedProperties();
};

/**
 * @brief Create PC-SAFT EOS from component names
 *
 * @param names Vector of component names (looked up in database)
 * @return Shared pointer to PC-SAFT EOS
 */
} // namespace PCSaft
} // namespace DMThermo

#endif // THERMO_EOS_PCSAFT_H
