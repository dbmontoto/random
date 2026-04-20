#ifndef COMPONENT_H
#define COMPONENT_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cmath>
#include <algorithm>
#include "thermo/dippr/equations.h"

namespace pcsaft {

/**
 * @brief Association scheme types for PC-SAFT with association
 *
 * Common schemes:
 * - NONE: No association (simple PC-SAFT)
 * - SCHEME_2B: Two association sites (e.g., ethanol, acetone)
 * - SCHEME_3B: Three association sites
 * - SCHEME_4C: Four association sites, two types (e.g., water, alcohols)
 */
enum class AssociationScheme {
    NONE,
    SCHEME_1A,  // One site type A
    SCHEME_2B,  // Two sites type B (e.g., acetone)
    SCHEME_3B,  // Three sites type B
    SCHEME_4C   // Four sites: 2A + 2B (e.g., water, alcohols)
};

/**
 * @brief Association site parameters for PC-SAFT
 */
struct AssociationParams {
    AssociationScheme scheme;
    double epsilon_AB;  // Association energy [K]
    double kappa_AB;    // Association volume [-]
    int num_sites;      // Total number of association sites

    AssociationParams()
        : scheme(AssociationScheme::NONE), epsilon_AB(0.0), kappa_AB(0.0), num_sites(0) {}

    AssociationParams(AssociationScheme s, double eps, double kap, int n)
        : scheme(s), epsilon_AB(eps), kappa_AB(kap), num_sites(n) {}

    bool hasAssociation() const { return scheme != AssociationScheme::NONE; }
};

/**
 * @brief Type alias for DIPPR ideal gas heat capacity coefficients
 *
 * Uses DIPPR::Coefficients from the standalone DIPPR module.
 * Supports polynomial (EQ_100) and Aly-Lee (EQ_107) equations.
 *
 * For polynomial form (default):
 *   Cp_ig = A + B*T + C*T^2 + D*T^3 + E*T^4  [J/(mol*K)]
 *
 * Methods available: calculate(), integrate(), integrateOverT(), isValid()
 */
using IdealGasCpParams = DIPPR::Coefficients;

/**
 * @brief Polymer parameters for PC-SAFT
 *
 * For polymers, segment number scales with molecular weight:
 *   m = m_per_mw * MW
 */
struct PolymerParams {
    double m_per_mw;    // Segment number per unit MW [mol/g]
    double Mn;          // Number-average molecular weight [g/mol]
    double Mw;          // Weight-average molecular weight [g/mol]

    PolymerParams() : m_per_mw(0.0), Mn(0.0), Mw(0.0) {}
    PolymerParams(double mpmw, double mn = 0.0, double mw = 0.0)
        : m_per_mw(mpmw), Mn(mn), Mw(mw) {}

    double getPDI() const { return (Mn > 0) ? Mw / Mn : 1.0; }
    bool isValid() const { return m_per_mw > 0.0; }
};

/**
 * @brief Component class containing PC-SAFT parameters
 */
class Component {
public:
    // Constructors
    Component();
    Component(const std::string& name, double m, double sigma, double epsilon_k);
    Component(const std::string& name, double m, double sigma, double epsilon_k,
              const AssociationParams& assoc);

    // Polymer constructor: uses m_per_mw * MW for segment number
    Component(const std::string& name, const PolymerParams& poly,
              double sigma, double epsilon_k, double mw);

    // Getters
    std::string getName() const { return name_; }
    std::string getCAS() const { return cas_; }
    double getM() const { return m_; }              // Segment number [-]
    double getSigma() const { return sigma_; }       // Segment diameter [Angstrom]
    double getEpsilonK() const { return epsilon_k_; } // Dispersion energy [K]
    double getMW() const { return mw_; }             // Molecular weight [g/mol]
    double getVolumeTranslation() const { return c_; } // Volume translation parameter [cm³/mol]

    const AssociationParams& getAssociationParams() const { return assoc_params_; }
    bool hasAssociation() const { return assoc_params_.hasAssociation(); }

    const IdealGasCpParams& getIdealGasCpParams() const { return cp_ig_params_; }
    bool hasIdealGasCp() const { return cp_ig_params_.isValid(); }

    // Polymer getters
    bool isPolymer() const { return is_polymer_; }
    const PolymerParams& getPolymerParams() const { return poly_params_; }
    double getMPerMW() const { return poly_params_.m_per_mw; }
    Component withMW(double new_mw) const;  // Create polymer with different MW

    // Setters
    void setName(const std::string& name) { name_ = name; }
    void setCAS(const std::string& cas) { cas_ = cas; }
    void setM(double m) { m_ = m; }
    void setSigma(double sigma) { sigma_ = sigma; }
    void setEpsilonK(double epsilon_k) { epsilon_k_ = epsilon_k; }
    void setMW(double mw) { mw_ = mw; }
    void setVolumeTranslation(double c) { c_ = c; }
    void setAssociationParams(const AssociationParams& params) { assoc_params_ = params; }
    void setIdealGasCpParams(const IdealGasCpParams& params) { cp_ig_params_ = params; }
    void setPolymerParams(const PolymerParams& params);

    // Display
    void print() const;

private:
    std::string name_;           // Component name
    std::string cas_;            // CAS number
    double m_;                   // Segment number [-]
    double sigma_;               // Segment diameter [Angstrom]
    double epsilon_k_;           // Dispersion energy parameter epsilon/k [K]
    double mw_;                  // Molecular weight [g/mol]
    double c_;                   // Volume translation parameter [cm³/mol] (for I-PC-SAFT)
    bool is_polymer_;            // True if this is a polymer component
    AssociationParams assoc_params_;   // Association parameters
    IdealGasCpParams cp_ig_params_;    // Ideal gas heat capacity parameters
    PolymerParams poly_params_;        // Polymer-specific parameters
};

/**
 * @brief Component database for loading common components
 */
class ComponentDatabase {
public:
    static ComponentDatabase& getInstance();

    // Load database from CSV databanks directory (pass any file within it, e.g. data/pcsaft.csv)
    bool loadFromFile(const std::string& filename);

    // Load default database (embedded or from default path)
    void loadDefault();

    // Get component by name or CAS
    Component getByName(const std::string& name) const;
    Component getByCAS(const std::string& cas) const;

    // Check if component exists
    bool hasComponent(const std::string& name) const;

    // Get all available component names
    std::vector<std::string> getAvailableComponents() const;

    // Add custom component to database
    void addComponent(const Component& comp);

private:
    ComponentDatabase() { loadDefault(); }
    ComponentDatabase(const ComponentDatabase&) = delete;
    ComponentDatabase& operator=(const ComponentDatabase&) = delete;

    void initializeDefaultDatabase();

    std::map<std::string, Component> components_by_name_;
    std::map<std::string, Component> components_by_cas_;
};

// Convenience function for getting components
inline Component getComponent(const std::string& name) {
    return ComponentDatabase::getInstance().getByName(name);
}

} // namespace pcsaft

#endif // COMPONENT_H
