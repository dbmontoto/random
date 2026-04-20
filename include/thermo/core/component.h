/**
 * @file component.h
 * @brief Thread-safe component definitions for thermodynamic calculations
 *
 * Components are immutable after construction, ensuring thread safety.
 */

#ifndef THERMO_CORE_COMPONENT_H
#define THERMO_CORE_COMPONENT_H

#include <string>
#include <vector>
#include <optional>
#include <memory>
#include <limits>
#include <stdexcept>
#include <cmath>

namespace DMThermo {
namespace Core {

/**
 * @brief Association scheme types for SAFT-based equations
 */
enum class AssociationScheme {
    None,           ///< No association
    Scheme1,        ///< 1 association site
    Scheme2A,       ///< 2 sites type A (e.g., simple acids)
    Scheme2B,       ///< 2 sites type B (e.g., alcohols, amines)
    Scheme3B,       ///< 3 sites (e.g., ammonia)
    Scheme4C        ///< 4 sites (e.g., water, glycols)
};

/**
 * @brief Association parameters for associating components
 */
struct AssociationParams {
    AssociationScheme scheme = AssociationScheme::None;
    double epsilon_AB = 0.0;     ///< Association energy [K]
    double kappa_AB = 0.0;       ///< Association volume [-]
    int num_sites = 0;           ///< Total number of association sites

    bool isAssociating() const { return scheme != AssociationScheme::None; }
};

/**
 * @brief Polymer-specific parameters
 */
struct PolymerParams {
    bool is_polymer = false;
    double m_per_mw = 0.0;       ///< Segment number per molecular weight [1/(g/mol)]
    double mw_monomer = 0.0;    ///< Monomer molecular weight [g/mol]

    /// Calculate segment number for given molecular weight
    double segmentNumber(double MW) const {
        return m_per_mw * MW;
    }
};

/**
 * @brief UNIFAC subgroup decomposition entry for a component
 */
struct UnifacSubgroupCount {
    int subgroup_id = 0;
    int count = 0;
};

/**
 * @brief VTPR volume-translation polynomial c(T) [m^3/mol]
 *
 * c(T) = c0 + c1*T + c2*T^2 + c3*T^3, with optional validity range.
 */
struct VtprVolumeTranslationPoly {
    double c0 = 0.0;
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0;
    std::optional<double> Tmin;
    std::optional<double> Tmax;

    double c(double T) const {
        return ((c3 * T + c2) * T + c1) * T + c0;
    }
};

/**
 * @brief Core-level ideal-gas Cp correlation coefficients (DIPPR form + coefficients only).
 *
 * This is a storage POD only. Parsing/evaluation is handled outside Core.
 */
struct IdealGasCpCorrelation {
    int eq_form = 0;
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double D = 0.0;
    double E = 0.0;
    double F = 0.0;
    double G = 0.0;
    double Tmin = 0.0;
    double Tmax = 1000.0;

    bool isValid() const {
        const bool finite_range = std::isfinite(Tmin) && std::isfinite(Tmax) && (Tmax > Tmin);
        const bool finite_coeffs =
            std::isfinite(A) && std::isfinite(B) && std::isfinite(C) &&
            std::isfinite(D) && std::isfinite(E) && std::isfinite(F) && std::isfinite(G);
        const bool any_nonzero = (A != 0.0) || (B != 0.0) || (C != 0.0) ||
                                 (D != 0.0) || (E != 0.0) || (F != 0.0) || (G != 0.0);
        if (!(eq_form > 0 && finite_range && finite_coeffs && any_nonzero)) {
            return false;
        }

        switch (eq_form) {
            case 105:
                return (B > 0.0) && (C > 0.0);
            case 102:
                return true;
            case 107:
                return true;
            case 114:
                return true;
            case 123:
                return true;
            case 100:
            case 101:
            case 106:
                return true;
            default:
                return false;
        }
    }

    static IdealGasCpCorrelation polynomial(
        double A, double B, double C, double D, double E,
        double Tmin = 200.0, double Tmax = 1500.0)
    {
        IdealGasCpCorrelation c;
        c.eq_form = 100;
        c.A = A;
        c.B = B;
        c.C = C;
        c.D = D;
        c.E = E;
        c.Tmin = Tmin;
        c.Tmax = Tmax;
        return c;
    }
};

/**
 * @brief Immutable component definition
 *
 * This class is fully thread-safe because all data members are
 * initialized at construction and never modified afterward.
 */
class Component {
public:
    // =========================================================================
    // Constructors
    // =========================================================================

    /**
     * @brief Construct component without PC-SAFT parameters
     *
     * Useful for cubic EOS baselines and databank-driven construction when
     * PC-SAFT parameters are not available.
     */
    explicit Component(const std::string& name);

    /**
     * @brief Construct non-associating component with PC-SAFT parameters
     *
     * @param name Component name
     * @param m Segment number [-]
     * @param sigma Segment diameter [Å]
     * @param epsilon_k Dispersion energy parameter [K]
     */
    Component(
        const std::string& name,
        double m,
        double sigma,
        double epsilon_k
    );

    /**
     * @brief Construct associating component with PC-SAFT parameters
     *
     * @param name Component name
     * @param m Segment number [-]
     * @param sigma Segment diameter [Å]
     * @param epsilon_k Dispersion energy parameter [K]
     * @param assoc Association parameters
     */
    Component(
        const std::string& name,
        double m,
        double sigma,
        double epsilon_k,
        const AssociationParams& assoc
    );

    /**
     * @brief Construct polymer component
     *
     * @param name Component name
     * @param sigma Segment diameter [Å]
     * @param epsilon_k Dispersion energy parameter [K]
     * @param poly Polymer parameters
     * @param MW Molecular weight [g/mol]
     */
    Component(
        const std::string& name,
        double sigma,
        double epsilon_k,
        const PolymerParams& poly,
        double MW
    );

    // Default copy/move operations are safe for immutable objects
    Component(const Component&) = default;
    Component(Component&&) = default;
    Component& operator=(const Component&) = default;
    Component& operator=(Component&&) = default;

    // =========================================================================
    // Accessors (All const, thread-safe)
    // =========================================================================

    /// Component name
    const std::string& name() const { return name_; }

    /// CAS registry number
    const std::string& cas() const { return cas_; }

    /// Segment number [-]
    double m() const { return m_; }

    /// Segment diameter [Å]
    double sigma() const { return sigma_; }

    /// Dispersion energy parameter [K]
    double epsilonK() const { return epsilon_k_; }

    /// Does this component have PC-SAFT parameters?
    bool hasPCSaftParams() const { return m_ > 0.0 && sigma_ > 0.0 && epsilon_k_ > 0.0; }

    /// Molecular weight [g/mol]
    double MW() const { return mw_; }

    /// Critical temperature [K]
    double Tc() const { return tc_; }

    /// Critical pressure [Pa]
    double Pc() const { return pc_; }

    /// Critical compressibility factor Zc [-] (optional; NaN when unset)
    double Zc() const { return zc_; }

    /// Acentric factor [-]
    double omega() const { return omega_; }

    /// Association parameters
    const AssociationParams& associationParams() const { return assoc_; }

    /// Is this an associating component?
    bool isAssociating() const { return assoc_.isAssociating(); }

    /// Polymer parameters
    const PolymerParams& polymerParams() const { return poly_; }

    /// Is this a polymer?
    bool isPolymer() const { return poly_.is_polymer; }

    /// Core-level ideal-gas Cp correlation coefficients
    const IdealGasCpCorrelation& idealGasCpCorrelation() const {
        if (!ideal_gas_cp_) throw std::runtime_error("Ideal-gas Cp correlation not set for component '" + name_ + "'");
        return *ideal_gas_cp_;
    }

    /// Does this component have ideal-gas Cp correlation data?
    bool hasIdealGasCp() const { return ideal_gas_cp_.has_value() && ideal_gas_cp_->isValid(); }

    /// Does this component have VTPR volume translation parameters?
    bool hasVtprVolumeTranslation() const { return vtpr_c_.has_value(); }

    /// VTPR volume translation c(T) [m^3/mol] (throws if not set)
    double vtprC(double T) const {
        if (!vtpr_c_) throw std::runtime_error("VTPR volume translation not set for component '" + name_ + "'");
        return vtpr_c_->c(T);
    }

    /// UNIFAC subgroup decomposition (may be empty)
    const std::vector<UnifacSubgroupCount>& unifacSubgroups() const { return unifac_groups_; }

    // =========================================================================
    // Builder Methods (Return new immutable Component)
    // =========================================================================

    /**
     * @brief Create copy with different molecular weight (for polymers)
     *
     * @param MW New molecular weight [g/mol]
     * @return New component with adjusted segment number
     */
    Component withMW(double MW) const;

    /**
     * @brief Create copy with CAS number set
     */
    Component withCAS(const std::string& cas) const;

    /**
     * @brief Create copy with critical properties set
     */
    Component withCritical(double Tc, double Pc, double omega) const;

    /**
     * @brief Create copy with critical compressibility set
     */
    Component withCriticalCompressibility(double Zc) const;

    /**
     * @brief Create copy with VTPR volume translation polynomial set (m^3/mol)
     */
    Component withVtprVolumeTranslation(const VtprVolumeTranslationPoly& c) const;

    /**
     * @brief Create copy with UNIFAC subgroup decomposition set
     */
    Component withUnifacSubgroups(const std::vector<UnifacSubgroupCount>& groups) const;

    /**
     * @brief Create copy with ideal-gas Cp correlation coefficients set
     */
    Component withIdealGasCp(const IdealGasCpCorrelation& coeffs) const;

    // =========================================================================
    // Utility Methods
    // =========================================================================

    /**
     * @brief Print component information
     */
    void print() const;

    /**
     * @brief Validate component parameters
     * @return true if parameters are physically reasonable
     */
    bool isValid() const;

private:
    std::string name_;
    std::string cas_;
    double m_;
    double sigma_;
    double epsilon_k_;
    double mw_ = 0.0;
    double tc_ = 0.0;
    double pc_ = 0.0;
    double zc_ = std::numeric_limits<double>::quiet_NaN();
    double omega_ = 0.0;
    std::optional<VtprVolumeTranslationPoly> vtpr_c_;
    std::vector<UnifacSubgroupCount> unifac_groups_;
    AssociationParams assoc_;
    PolymerParams poly_;
    std::optional<IdealGasCpCorrelation> ideal_gas_cp_;
};

/**
 * @brief Shared pointer type alias
 */
using ComponentPtr = std::shared_ptr<const Component>;

} // namespace Core
} // namespace DMThermo

#endif // THERMO_CORE_COMPONENT_H
