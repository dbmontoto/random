/**
 * @file mixture.h
 * @brief Thread-safe mixture definitions for thermodynamic calculations
 *
 * Mixtures are immutable after construction, ensuring thread safety.
 */

#ifndef THERMO_CORE_MIXTURE_H
#define THERMO_CORE_MIXTURE_H

#include <vector>
#include <string>
#include <memory>
#include <optional>
#include "thermo/core/component.h"
#include "thermo/data/vtpr.h"

namespace DMThermo {
namespace Core {

/**
 * @brief Binary interaction parameters storage
 */
struct BinaryParameters {
    /// kij interaction parameters (symmetric matrix stored as upper triangular)
    std::vector<double> kij;

    /// Number of components
    int n = 0;

    /// Get kij for component pair
    double get(int i, int j) const;

    /// Set kij for component pair (modifies during construction only)
    void set(int i, int j, double value);

    /// Create with all zeros
    static BinaryParameters zeros(int n);
};

/**
 * @brief Immutable mixture definition
 *
 * This class is fully thread-safe because all data members are
 * initialized at construction and never modified afterward.
 */
class Mixture {
public:
    // =========================================================================
    // Constructors
    // =========================================================================

    /**
     * @brief Construct mixture from components
     *
     * @param components Vector of components
     * @param kij Binary interaction parameters (optional, defaults to zero)
     */
    explicit Mixture(
        const std::vector<Component>& components,
        const BinaryParameters& kij = BinaryParameters{},
        std::shared_ptr<const Data::UnifacTables> unifac_tables = nullptr
    );

    /**
     * @brief Construct mixture from component pointers
     *
     * @param components Vector of component shared pointers
     * @param kij Binary interaction parameters (optional)
     */
    explicit Mixture(
        const std::vector<ComponentPtr>& components,
        const BinaryParameters& kij = BinaryParameters{},
        std::shared_ptr<const Data::UnifacTables> unifac_tables = nullptr
    );

    // Default copy/move operations are safe for immutable objects
    Mixture(const Mixture&) = default;
    Mixture(Mixture&&) = default;
    Mixture& operator=(const Mixture&) = default;
    Mixture& operator=(Mixture&&) = default;

    // =========================================================================
    // Accessors (All const, thread-safe)
    // =========================================================================

    /// Number of components
    int numComponents() const { return static_cast<int>(components_.size()); }

    /// Get component by index
    const Component& component(int i) const { return *components_[i]; }

    /// Get all components
    const std::vector<ComponentPtr>& components() const { return components_; }

    /// Get binary interaction parameter kij
    double kij(int i, int j) const { return kij_.get(i, j); }

    /// Get all binary parameters
    const BinaryParameters& binaryParameters() const { return kij_; }

    /// UNIFAC parameter tables (optional; required for VTPR/UNIFAC-based models)
    std::shared_ptr<const Data::UnifacTables> unifacTables() const { return unifac_tables_; }

    /// Does mixture have any associating components?
    bool hasAssociation() const { return has_association_; }

    /// Does mixture have any polymer components?
    bool hasPolymers() const { return has_polymers_; }

    /// Total number of association sites in mixture
    int totalAssociationSites() const { return total_sites_; }

    // =========================================================================
    // Mixture Property Calculations (Pure functions, thread-safe)
    // =========================================================================

    /**
     * @brief Calculate average molecular weight
     *
     * @param x Mole fractions
     * @return Average MW [g/mol]
     */
    double averageMW(const std::vector<double>& x) const;

    /**
     * @brief Calculate average segment number
     *
     * @param x Mole fractions
     * @return m_avg [-]
     */
    double averageSegmentNumber(const std::vector<double>& x) const;

    /**
     * @brief Calculate mixing rule for sigma (segment diameter)
     *
     * @param i Component index
     * @param j Component index
     * @return sigma_ij [Å]
     */
    double sigmaMix(int i, int j) const;

    /**
     * @brief Calculate mixing rule for epsilon (with kij correction)
     *
     * @param i Component index
     * @param j Component index
     * @return epsilon_ij/k [K]
     */
    double epsilonMix(int i, int j) const;

    // =========================================================================
    // Builder Methods (Return new immutable Mixture)
    // =========================================================================

    /**
     * @brief Create mixture with modified kij parameter
     *
     * @param i Component index
     * @param j Component index
     * @param value New kij value
     * @return New mixture with updated parameter
     */
    Mixture withKij(int i, int j, double value) const;

    /**
     * @brief Create mixture with all new kij parameters
     *
     * @param kij New binary parameters
     * @return New mixture with updated parameters
     */
    Mixture withBinaryParameters(const BinaryParameters& kij) const;

    // =========================================================================
    // Utility Methods
    // =========================================================================

    /**
     * @brief Validate composition vector
     *
     * @param x Mole fractions to validate
     * @return true if valid (all positive, sums to ~1)
     */
    bool isValidComposition(const std::vector<double>& x) const;

    /**
     * @brief Normalize composition to sum to 1
     *
     * @param x Mole fractions
     * @return Normalized mole fractions
     */
    static std::vector<double> normalizeComposition(const std::vector<double>& x);

    /**
     * @brief Print mixture information
     */
    void print() const;

private:
    std::vector<ComponentPtr> components_;
    BinaryParameters kij_;
    std::shared_ptr<const Data::UnifacTables> unifac_tables_;
    bool has_association_ = false;
    bool has_polymers_ = false;
    int total_sites_ = 0;

    void initializeDerivedProperties();
};

/**
 * @brief Shared pointer type alias
 */
using MixturePtr = std::shared_ptr<const Mixture>;

} // namespace Core
} // namespace DMThermo

#endif // THERMO_CORE_MIXTURE_H
