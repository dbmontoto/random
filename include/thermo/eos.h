/**
 * @file eos.h
 * @brief Abstract interface for equations of state
 *
 * This interface defines the contract that all equation of state
 * implementations must fulfill. The design prioritizes:
 * - Thread safety via const methods
 * - Flexibility via optional derivative computation
 * - Extensibility via virtual methods
 */

#ifndef THERMO_EOS_H
#define THERMO_EOS_H

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "core/types.h"

namespace DMThermo {
namespace Core {
class Mixture;
} // namespace Core

/**
 * @brief Abstract interface for equations of state
 *
 * All methods are const to ensure thread safety. Implementations
 * should store immutable parameters and use thread-local caching
 * for intermediate calculations.
 *
 * Example usage:
 * @code
 * auto eos = DMThermo::Factory::EOSFactory::create("PC-SAFT", mixture);
 * EOSResult result = eos->calculate(300.0, 1000.0, {0.5, 0.5});
 * double P = result.pressure;
 * @endcode
 */
class EOS {
public:
    virtual ~EOS() = default;

    // =========================================================================
    // Core Calculation Methods (Thread-Safe)
    // =========================================================================

    /**
     * @brief Complete EOS calculation at given state
     *
     * @param T Temperature [K]
     * @param rho Molar density [mol/m^3]
     * @param x Mole fractions (must sum to 1)
     * @param deriv_spec Which derivatives to compute
     * @return EOSResult containing all calculated properties
     *
     * @note This method is thread-safe when called with the same or different
     *       state parameters from multiple threads.
     */
    virtual EOSResult calculate(
        double T,
        double rho,
        const std::vector<double>& x,
        const DerivativeSpec& deriv_spec = DerivativeSpec::none()
    ) const = 0;

    // =========================================================================
    // Convenience Methods (All Thread-Safe)
    // =========================================================================

    /**
     * @brief Calculate pressure at given state
     *
     * @param T Temperature [K]
     * @param rho Molar density [mol/m^3]
     * @param x Mole fractions
     * @return Pressure [Pa]
     */
    virtual double pressure(
        double T,
        double rho,
        const std::vector<double>& x
    ) const = 0;

    /**
     * @brief Calculate fugacity coefficients at given state
     *
     * @param T Temperature [K]
     * @param rho Molar density [mol/m^3]
     * @param x Mole fractions
     * @return Vector of fugacity coefficients phi_i
     */
    virtual std::vector<double> fugacityCoefficients(
        double T,
        double rho,
        const std::vector<double>& x
    ) const = 0;

    /**
     * @brief Calculate compressibility factor at given state
     *
     * @param T Temperature [K]
     * @param rho Molar density [mol/m^3]
     * @param x Mole fractions
     * @return Compressibility factor Z = P/(rho*R*T)
     */
    virtual double compressibility(
        double T,
        double rho,
        const std::vector<double>& x
    ) const = 0;

    /**
     * @brief Calculate residual Helmholtz energy a_res/(RT)
     *
     * @param T Temperature [K]
     * @param rho Molar density [mol/m^3]
     * @param x Mole fractions
     * @return Reduced residual Helmholtz energy a_res/(RT)
     */
    virtual double residualHelmholtz(
        double T,
        double rho,
        const std::vector<double>& x
    ) const = 0;

    /**
     * @brief Optional: real density roots at (T, P, x).
     *
     * This hook exists to support EOS implementations that can provide multiple TP
     * density roots directly (e.g., cubic EOS via analytic Z roots), enabling
     * stability/flash utilities to avoid expensive rho-scans.
     *
     * Implementations should return unique, positive roots sorted ascending. The
     * default implementation returns an empty vector (not supported).
     */
    virtual std::vector<double> densityRootsTP(
        double T,
        double P,
        const std::vector<double>& x
    ) const
    {
        (void)T;
        (void)P;
        (void)x;
        return {};
    }

    /**
     * @brief Calculate residual Gibbs energy g_res/(RT) at the given state.
     *
     * For any EOS returning fugacity coefficients:
     *   g_res/(RT) = Σ_i x_i ln(φ_i)
     *
     * This is a TP-style residual (departure from ideal gas at the same T,P).
     */
    virtual double residualGibbs(
        double T,
        double rho,
        const std::vector<double>& x
    ) const
    {
        const auto r = calculate(T, rho, x, DerivativeSpec::none());
        if (!r.success) {
            throw std::runtime_error("EOS::residualGibbs: calculate() failed: " + r.error_message);
        }
        if (r.ln_phi.size() != x.size()) {
            throw std::runtime_error("EOS::residualGibbs: ln_phi/composition size mismatch");
        }
        double g_res_over_rt = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            g_res_over_rt += x[i] * r.ln_phi[i];
        }
        return g_res_over_rt;
    }

    // =========================================================================
    // Derivative Methods
    // =========================================================================

    /**
     * @brief Calculate d(a_res/RT)/dT at constant rho, x
     *
     * @param T Temperature [K]
     * @param rho Molar density [mol/m^3]
     * @param x Mole fractions
     * @return Temperature derivative [1/K]
     */
    virtual double dadt(
        double T,
        double rho,
        const std::vector<double>& x
    ) const = 0;

    /**
     * @brief Calculate d(a_res/RT)/drho at constant T, x
     *
     * @param T Temperature [K]
     * @param rho Molar density [mol/m^3]
     * @param x Mole fractions
     * @return Density derivative [m^3/mol]
     */
    virtual double dadrho(
        double T,
        double rho,
        const std::vector<double>& x
    ) const = 0;

    /**
     * @brief Calculate dP/drho at constant T, x (for stability analysis)
     *
     * @param T Temperature [K]
     * @param rho Molar density [mol/m^3]
     * @param x Mole fractions
     * @return Pressure-density derivative [Pa*m^3/mol]
     */
    virtual double dPdrho(
        double T,
        double rho,
        const std::vector<double>& x
    ) const = 0;

    // =========================================================================
    // Metadata Methods
    // =========================================================================

    /**
     * @brief Get the name of this EOS implementation
     * @return Human-readable name (e.g., "PC-SAFT", "Peng-Robinson")
     */
    virtual std::string name() const = 0;

    /**
     * @brief Get number of components in the mixture
     * @return Number of components
     */
    virtual int numComponents() const = 0;

    /**
     * @brief Check if this EOS handles association
     * @return true if association contribution is included
     */
    virtual bool hasAssociation() const = 0;

    /**
     * @brief Get critical temperature estimate for component i
     * @param i Component index
     * @return Critical temperature [K] (estimate)
     */
    virtual double criticalTemperature(int i) const = 0;

    /**
     * @brief Get critical pressure estimate for component i
     * @param i Component index
     * @return Critical pressure [Pa] (estimate)
     */
    virtual double criticalPressure(int i) const = 0;

    /**
     * @brief Get acentric factor for component i
     * @param i Component index
     * @return Acentric factor [-]
     */
    virtual double acentricFactor(int i) const = 0;

    /**
     * @brief Optional access to the underlying Core::Mixture (if available).
     *
     * This enables EOS-agnostic tooling (e.g., PH/PS flash outer solves) to compute
     * ideal-gas contributions using the mixture's Cp correlations.
     *
     * @return Pointer to an immutable Core::Mixture, or nullptr if unavailable.
     */
    virtual const Core::Mixture* coreMixture() const { return nullptr; }
};

/**
 * @brief Shared pointer type alias for EOS
 */
using EOSPtr = std::shared_ptr<EOS>;

/**
 * @brief Const shared pointer type alias for EOS
 */
using EOSConstPtr = std::shared_ptr<const EOS>;

} // namespace DMThermo

#endif // THERMO_EOS_H
