/**
 * @file mixture.cpp
 * @brief Implementation of DMThermo::Core::Mixture
 */

#include "thermo/core/mixture.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <numeric>

namespace DMThermo {
namespace Core {

// =============================================================================
// BinaryParameters Implementation
// =============================================================================

double BinaryParameters::get(int i, int j) const {
    if (n == 0) return 0.0;
    if (i > j) std::swap(i, j);  // Ensure upper triangular access
    int idx = i * n - i * (i + 1) / 2 + j;
    if (idx < 0 || idx >= static_cast<int>(kij.size())) return 0.0;
    return kij[idx];
}

void BinaryParameters::set(int i, int j, double value) {
    if (n == 0) return;
    if (i > j) std::swap(i, j);  // Ensure upper triangular storage
    int idx = i * n - i * (i + 1) / 2 + j;
    if (idx >= 0 && idx < static_cast<int>(kij.size())) {
        kij[idx] = value;
    }
}

BinaryParameters BinaryParameters::zeros(int n) {
    BinaryParameters bp;
    bp.n = n;
    // Upper triangular storage: n*(n+1)/2 elements
    bp.kij.resize(n * (n + 1) / 2, 0.0);
    return bp;
}

// =============================================================================
// Mixture Constructors
// =============================================================================

Mixture::Mixture(
    const std::vector<Component>& components,
    const BinaryParameters& kij,
    std::shared_ptr<const Data::UnifacTables> unifac_tables)
    : kij_(kij.n > 0 ? kij : BinaryParameters::zeros(static_cast<int>(components.size())))
    , unifac_tables_(std::move(unifac_tables))
{
    components_.reserve(components.size());
    for (const auto& comp : components) {
        components_.push_back(std::make_shared<const Component>(comp));
    }
    initializeDerivedProperties();
}

Mixture::Mixture(
    const std::vector<ComponentPtr>& components,
    const BinaryParameters& kij,
    std::shared_ptr<const Data::UnifacTables> unifac_tables)
    : components_(components)
    , kij_(kij.n > 0 ? kij : BinaryParameters::zeros(static_cast<int>(components.size())))
    , unifac_tables_(std::move(unifac_tables))
{
    initializeDerivedProperties();
}

void Mixture::initializeDerivedProperties() {
    has_association_ = false;
    has_polymers_ = false;
    total_sites_ = 0;

    for (const auto& comp : components_) {
        if (comp->isAssociating()) {
            has_association_ = true;
            total_sites_ += comp->associationParams().num_sites;
        }
        if (comp->isPolymer()) {
            has_polymers_ = true;
        }
    }
}

// =============================================================================
// Mixture Property Calculations
// =============================================================================

double Mixture::averageMW(const std::vector<double>& x) const {
    double mw_avg = 0.0;
    for (size_t i = 0; i < components_.size(); ++i) {
        mw_avg += x[i] * components_[i]->MW();
    }
    return mw_avg;
}

double Mixture::averageSegmentNumber(const std::vector<double>& x) const {
    double m_avg = 0.0;
    for (size_t i = 0; i < components_.size(); ++i) {
        m_avg += x[i] * components_[i]->m();
    }
    return m_avg;
}

double Mixture::sigmaMix(int i, int j) const {
    // Lorentz combining rule: sigma_ij = (sigma_i + sigma_j) / 2
    return 0.5 * (components_[i]->sigma() + components_[j]->sigma());
}

double Mixture::epsilonMix(int i, int j) const {
    // Berthelot combining rule with kij correction:
    // epsilon_ij = sqrt(epsilon_i * epsilon_j) * (1 - kij)
    double eps_i = components_[i]->epsilonK();
    double eps_j = components_[j]->epsilonK();
    return std::sqrt(eps_i * eps_j) * (1.0 - kij_.get(i, j));
}

// =============================================================================
// Builder Methods
// =============================================================================

Mixture Mixture::withKij(int i, int j, double value) const {
    BinaryParameters new_kij = kij_;
    new_kij.set(i, j, value);
    return Mixture(components_, new_kij, unifac_tables_);
}

Mixture Mixture::withBinaryParameters(const BinaryParameters& kij) const {
    return Mixture(components_, kij, unifac_tables_);
}

// =============================================================================
// Utility Methods
// =============================================================================

bool Mixture::isValidComposition(const std::vector<double>& x) const {
    if (x.size() != components_.size()) return false;

    double sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] < 0.0 || x[i] > 1.0) return false;
        sum += x[i];
    }

    // Allow small tolerance for sum to 1
    return std::abs(sum - 1.0) < 1e-10;
}

std::vector<double> Mixture::normalizeComposition(const std::vector<double>& x) {
    double sum = std::accumulate(x.begin(), x.end(), 0.0);
    if (sum <= 0.0) {
        throw std::invalid_argument("Cannot normalize composition with non-positive sum");
    }

    std::vector<double> normalized(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        normalized[i] = x[i] / sum;
    }
    return normalized;
}

void Mixture::print() const {
    std::cout << "Mixture with " << numComponents() << " components:\n";
    std::cout << std::string(60, '-') << "\n";

    for (size_t i = 0; i < components_.size(); ++i) {
        const auto& comp = components_[i];
        std::cout << "  [" << i << "] " << comp->name();
        std::cout << " (m=" << std::fixed << std::setprecision(4) << comp->m();
        std::cout << ", sigma=" << comp->sigma() << " Å";
        std::cout << ", eps/k=" << comp->epsilonK() << " K";
        if (comp->isAssociating()) std::cout << ", ASSOC";
        if (comp->isPolymer()) std::cout << ", POLYMER";
        std::cout << ")\n";
    }

    // Print non-zero binary parameters
    bool has_nonzero_kij = false;
    for (int i = 0; i < numComponents(); ++i) {
        for (int j = i + 1; j < numComponents(); ++j) {
            double k = kij_.get(i, j);
            if (std::abs(k) > 1e-15) {
                if (!has_nonzero_kij) {
                    std::cout << "  Binary interaction parameters:\n";
                    has_nonzero_kij = true;
                }
                std::cout << "    k[" << i << "," << j << "] = " << k << "\n";
            }
        }
    }

    std::cout << std::string(60, '-') << "\n";
}

} // namespace Core
} // namespace DMThermo
