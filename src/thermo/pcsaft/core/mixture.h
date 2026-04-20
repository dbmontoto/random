#ifndef MIXTURE_H
#define MIXTURE_H

#include "thermo/pcsaft/core/component.h"
#include <vector>
#include <map>
#include <string>

namespace pcsaft {

/**
 * @brief Mixture class for handling multi-component systems
 *
 * Manages components, mole fractions, and binary interaction parameters.
 */
class Mixture {
public:
    // Constructors
    Mixture();
    Mixture(const std::vector<Component>& components, const std::vector<double>& mole_fractions);

    // Component management
    void addComponent(const Component& comp, double mole_fraction);
    void setMoleFractions(const std::vector<double>& mole_fractions);
    void normalize(); // Normalize mole fractions to sum to 1.0

    // Getters
    int getNumComponents() const { return static_cast<int>(components_.size()); }
    const std::vector<Component>& getComponents() const { return components_; }
    const std::vector<double>& getMoleFractions() const { return mole_fractions_; }
    const Component& getComponent(int i) const { return components_[i]; }
    double getMoleFraction(int i) const { return mole_fractions_[i]; }

    // Binary interaction parameters (kij)
    void setBinaryParameter(int i, int j, double kij);
    void setBinaryParameter(const std::string& comp1, const std::string& comp2, double kij);
    double getBinaryParameter(int i, int j) const;

    // Mixture properties (average values)
    double getMixtureMW() const; // Molecular weight [g/mol]

    // Check if mixture contains associating components
    bool hasAssociation() const;

    // Validation
    bool isValid() const;
    void validate() const; // Throws if invalid

    // Display
    void print() const;

private:
    std::vector<Component> components_;
    std::vector<double> mole_fractions_;

    // Binary interaction parameters stored as symmetric matrix
    // kij_[i][j] = kij_[j][i]
    std::vector<std::vector<double>> kij_;

    void initializeKij();
    int getComponentIndex(const std::string& name) const;
};

/**
 * @brief Helper function to create a mixture from component names
 */
Mixture createMixture(const std::vector<std::string>& component_names,
                      const std::vector<double>& mole_fractions);

/**
 * @brief Helper function to create a binary mixture
 */
Mixture createBinaryMixture(const std::string& comp1, const std::string& comp2,
                           double x1, double kij = 0.0);

} // namespace pcsaft

#endif // MIXTURE_H
