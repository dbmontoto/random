/**
 * @file state.h
 * @brief Thermodynamic state representation
 *
 * Provides an immutable representation of a thermodynamic state
 * for passing to property calculations.
 */

#ifndef THERMO_CORE_STATE_H
#define THERMO_CORE_STATE_H

#include <vector>
#include <optional>
#include <memory>
#include "types.h"
#include "thermo/core/mixture.h"

namespace DMThermo {
namespace Core {

/**
 * @brief Immutable thermodynamic state
 *
 * Represents a complete thermodynamic state for property calculations.
 * All state variables are immutable after construction.
 */
class State {
public:
    // =========================================================================
    // Factory Methods (Named constructors for clarity)
    // =========================================================================

    /**
     * @brief Create state from T, P, x specification
     *
     * @param mixture Mixture definition
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param x Mole fractions
     * @param phase Phase type (optional, for density selection)
     * @return State object
     */
    static State fromTP(
        MixturePtr mixture,
        double T,
        double P,
        const std::vector<double>& x,
        PhaseType phase = PhaseType::Unknown
    );

    /**
     * @brief Create state from T, rho, x specification
     *
     * @param mixture Mixture definition
     * @param T Temperature [K]
     * @param rho Molar density [mol/m³]
     * @param x Mole fractions
     * @return State object
     */
    static State fromTRho(
        MixturePtr mixture,
        double T,
        double rho,
        const std::vector<double>& x
    );

    /**
     * @brief Create state from T, V, x specification
     *
     * @param mixture Mixture definition
     * @param T Temperature [K]
     * @param V Molar volume [m³/mol]
     * @param x Mole fractions
     * @return State object
     */
    static State fromTV(
        MixturePtr mixture,
        double T,
        double V,
        const std::vector<double>& x
    );

    // Default copy/move operations
    State(const State&) = default;
    State(State&&) = default;
    State& operator=(const State&) = default;
    State& operator=(State&&) = default;

    // =========================================================================
    // Accessors (All const, thread-safe)
    // =========================================================================

    /// Temperature [K]
    double T() const { return T_; }

    /// Molar density [mol/m³]
    double rho() const { return rho_; }

    /// Molar volume [m³/mol]
    double V() const { return (rho_ > 0) ? 1.0 / rho_ : 0.0; }

    /// Pressure [Pa] (if known from specification)
    std::optional<double> P() const { return P_; }

    /// Mole fractions
    const std::vector<double>& x() const { return x_; }

    /// Mole fraction of component i
    double x(int i) const { return x_[i]; }

    /// Phase type
    PhaseType phase() const { return phase_; }

    /// Mixture reference
    MixturePtr mixture() const { return mixture_; }

    /// Number of components
    int numComponents() const { return static_cast<int>(x_.size()); }

    // =========================================================================
    // Derived Quantities
    // =========================================================================

    /**
     * @brief Calculate average molecular weight
     * @return Average MW [g/mol]
     */
    double averageMW() const;

    /**
     * @brief Calculate mass density
     * @return Mass density [kg/m³]
     */
    double massDensity() const;

    /**
     * @brief Calculate specific volume
     * @return Specific volume [m³/kg]
     */
    double specificVolume() const;

    // =========================================================================
    // Builder Methods (Return new immutable State)
    // =========================================================================

    /**
     * @brief Create state with updated temperature
     */
    State withT(double T) const;

    /**
     * @brief Create state with updated density
     */
    State withRho(double rho) const;

    /**
     * @brief Create state with updated composition
     */
    State withX(const std::vector<double>& x) const;

    /**
     * @brief Create state with phase type set
     */
    State withPhase(PhaseType phase) const;

    /**
     * @brief Create state with pressure set (for caching)
     */
    State withP(double P) const;

private:
    State() = default;  // Private, use factory methods

    MixturePtr mixture_;
    double T_ = 0.0;
    double rho_ = 0.0;
    std::optional<double> P_;
    std::vector<double> x_;
    PhaseType phase_ = PhaseType::Unknown;
};

} // namespace Core
} // namespace DMThermo

#endif // THERMO_CORE_STATE_H
