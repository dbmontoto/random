#include "thermo/pcsaft/core/mixture.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <mutex>
#include <optional>

#include "thermo/data/databanks.h"

namespace {

std::optional<std::filesystem::path> findDefaultBipCsv() {
    std::filesystem::path p = std::filesystem::current_path();
    for (int i = 0; i < 6; ++i) {
        const auto candidate = p / "data" / "bip.csv";
        if (std::filesystem::exists(candidate)) return candidate;
        if (!p.has_parent_path()) break;
        p = p.parent_path();
    }
    return std::nullopt;
}

const DMThermo::Data::BipCsvDatabank* defaultBipDatabank() {
    static std::once_flag once;
    static DMThermo::Data::BipCsvDatabank bip;
    static bool loaded = false;

    std::call_once(once, [&] {
        const auto csv = findDefaultBipCsv();
        if (!csv) {
            loaded = false;
            return;
        }
        DMThermo::Data::Diagnostics diag;
        loaded = bip.load(*csv, &diag);
    });

    return loaded ? &bip : nullptr;
}

} // namespace

namespace pcsaft {

// ============================================================================
// Mixture Implementation
// ============================================================================

Mixture::Mixture() {}

Mixture::Mixture(const std::vector<Component>& components,
                 const std::vector<double>& mole_fractions)
    : components_(components), mole_fractions_(mole_fractions) {

    if (components.size() != mole_fractions.size()) {
        throw std::invalid_argument("Number of components must match number of mole fractions");
    }

    initializeKij();
    // Apply CSV databank BIPs if available. Missing BIPs default to 0.
    const auto* bip = defaultBipDatabank();
    if (bip) {
        constexpr double Tref = 298.15;
        const int n = getNumComponents();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const std::string& cas_i = components_[i].getCAS();
                const std::string& cas_j = components_[j].getCAS();
                if (cas_i.empty() || cas_j.empty()) continue;
                const auto k = bip->kij(DMThermo::Data::BipModel::PCSaft, cas_i, cas_j, Tref, nullptr);
                if (k) {
                    setBinaryParameter(i, j, *k);
                }
            }
        }
    }
    validate();
}

void Mixture::addComponent(const Component& comp, double mole_fraction) {
    components_.push_back(comp);
    mole_fractions_.push_back(mole_fraction);
    initializeKij();
    // Re-apply databank BIPs after reshaping kij matrix.
    const auto* bip = defaultBipDatabank();
    if (bip) {
        constexpr double Tref = 298.15;
        const int n = getNumComponents();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const std::string& cas_i = components_[i].getCAS();
                const std::string& cas_j = components_[j].getCAS();
                if (cas_i.empty() || cas_j.empty()) continue;
                const auto k = bip->kij(DMThermo::Data::BipModel::PCSaft, cas_i, cas_j, Tref, nullptr);
                if (k) {
                    setBinaryParameter(i, j, *k);
                }
            }
        }
    }
}

void Mixture::setMoleFractions(const std::vector<double>& mole_fractions) {
    if (mole_fractions.size() != components_.size()) {
        throw std::invalid_argument("Number of mole fractions must match number of components");
    }
    mole_fractions_ = mole_fractions;
}

void Mixture::normalize() {
    double sum = 0.0;
    for (double x : mole_fractions_) {
        sum += x;
    }

    if (sum <= 0.0) {
        throw std::runtime_error("Cannot normalize: sum of mole fractions is zero or negative");
    }

    for (double& x : mole_fractions_) {
        x /= sum;
    }
}

void Mixture::initializeKij() {
    int n = getNumComponents();
    kij_.resize(n);
    for (int i = 0; i < n; ++i) {
        kij_[i].resize(n, 0.0); // Default kij = 0
    }
}

void Mixture::setBinaryParameter(int i, int j, double kij) {
    if (i < 0 || i >= getNumComponents() || j < 0 || j >= getNumComponents()) {
        throw std::out_of_range("Component indices out of range");
    }

    kij_[i][j] = kij;
    kij_[j][i] = kij; // Symmetric
}

void Mixture::setBinaryParameter(const std::string& comp1, const std::string& comp2, double kij) {
    int i = getComponentIndex(comp1);
    int j = getComponentIndex(comp2);
    setBinaryParameter(i, j, kij);
}

double Mixture::getBinaryParameter(int i, int j) const {
    if (i < 0 || i >= getNumComponents() || j < 0 || j >= getNumComponents()) {
        throw std::out_of_range("Component indices out of range");
    }
    return kij_[i][j];
}

int Mixture::getComponentIndex(const std::string& name) const {
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    for (int i = 0; i < getNumComponents(); ++i) {
        std::string comp_name = components_[i].getName();
        std::transform(comp_name.begin(), comp_name.end(), comp_name.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (comp_name == name_lower) {
            return i;
        }
    }

    throw std::runtime_error("Component '" + name + "' not found in mixture");
}

double Mixture::getMixtureMW() const {
    double mw = 0.0;
    for (int i = 0; i < getNumComponents(); ++i) {
        mw += mole_fractions_[i] * components_[i].getMW();
    }
    return mw;
}

bool Mixture::hasAssociation() const {
    for (const auto& comp : components_) {
        if (comp.hasAssociation()) {
            return true;
        }
    }
    return false;
}

bool Mixture::isValid() const {
    if (components_.empty()) {
        return false;
    }

    if (components_.size() != mole_fractions_.size()) {
        return false;
    }

    double sum = 0.0;
    for (double x : mole_fractions_) {
        if (x < 0.0) {
            return false; // Negative mole fraction
        }
        sum += x;
    }

    // Check if sum is approximately 1.0 (within tolerance)
    const double tolerance = 1e-6;
    if (std::abs(sum - 1.0) > tolerance) {
        return false;
    }

    return true;
}

void Mixture::validate() const {
    if (!isValid()) {
        throw std::runtime_error("Invalid mixture: check components and mole fractions");
    }
}

void Mixture::print() const {
    std::cout << "Mixture with " << getNumComponents() << " components:\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nComponent           Mole Fraction\n";
    std::cout << "----------------------------------------\n";

    for (int i = 0; i < getNumComponents(); ++i) {
        std::cout << std::setw(20) << std::left << components_[i].getName()
                  << std::setw(15) << std::right << mole_fractions_[i] << "\n";
    }

    // Print non-zero kij values
    bool has_kij = false;
    for (int i = 0; i < getNumComponents(); ++i) {
        for (int j = i+1; j < getNumComponents(); ++j) {
            if (std::abs(kij_[i][j]) > 1e-10) {
                has_kij = true;
                break;
            }
        }
        if (has_kij) break;
    }

    if (has_kij) {
        std::cout << "\nBinary Interaction Parameters:\n";
        std::cout << "----------------------------------------\n";
        for (int i = 0; i < getNumComponents(); ++i) {
            for (int j = i+1; j < getNumComponents(); ++j) {
                if (std::abs(kij_[i][j]) > 1e-10) {
                    std::cout << components_[i].getName() << " - "
                              << components_[j].getName() << ": kij = "
                              << kij_[i][j] << "\n";
                }
            }
        }
    }

    if (hasAssociation()) {
        std::cout << "\nMixture contains associating components.\n";
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

Mixture createMixture(const std::vector<std::string>& component_names,
                      const std::vector<double>& mole_fractions) {
    if (component_names.size() != mole_fractions.size()) {
        throw std::invalid_argument("Number of component names must match number of mole fractions");
    }

    std::vector<Component> components;
    components.reserve(component_names.size());

    for (const auto& name : component_names) {
        components.push_back(getComponent(name));
    }

    return Mixture(components, mole_fractions);
}

Mixture createBinaryMixture(const std::string& comp1, const std::string& comp2,
                           double x1, double kij) {
    Component c1 = getComponent(comp1);
    Component c2 = getComponent(comp2);

    std::vector<Component> components = {c1, c2};
    std::vector<double> mole_fractions = {x1, 1.0 - x1};

    Mixture mix(components, mole_fractions);

    if (std::abs(kij) > 1e-10) {
        mix.setBinaryParameter(0, 1, kij);
    }

    return mix;
}

} // namespace pcsaft
