/**
 * @file conversions.cpp
 * @brief Implementation of conversion utilities between legacy and new types
 */

#include "thermo/legacy/pcsaft/conversions.h"
#include "thermo/legacy/dippr/equations.h"
#include "thermo/legacy/pcsaft/core/component.h"   // legacy pcsaft::Component
#include "thermo/legacy/pcsaft/core/mixture.h"     // legacy pcsaft::Mixture
#include <string>
#include <stdexcept>

namespace DMThermo {
namespace Core {

namespace {

int equationFormFromLegacyType(DIPPR::EquationType eq) {
    switch (eq) {
        case DIPPR::EquationType::EQ_100: return 100;
        case DIPPR::EquationType::EQ_101: return 101;
        case DIPPR::EquationType::EQ_102: return 102;
        case DIPPR::EquationType::EQ_105: return 105;
        case DIPPR::EquationType::EQ_106: return 106;
        case DIPPR::EquationType::EQ_107: return 107;
        case DIPPR::EquationType::EQ_114: return 114;
        case DIPPR::EquationType::EQ_123: return 123;
        default: return 0;
    }
}

DIPPR::Coefficients toLegacyIdealGasCp(const IdealGasCpCorrelation& corr) {
    const auto eq = DIPPR::equationTypeFromForm(corr.eq_form);
    if (!eq) {
        throw std::runtime_error("Unsupported ideal-gas Cp equation form " + std::to_string(corr.eq_form));
    }
    DIPPR::Coefficients coeffs(
        *eq,
        corr.A, corr.B, corr.C, corr.D, corr.E, corr.F, corr.G,
        corr.Tmin, corr.Tmax);
    if (!coeffs.isValid()) {
        throw std::runtime_error("Invalid ideal-gas Cp coefficients for eq_form " + std::to_string(corr.eq_form));
    }
    return coeffs;
}

IdealGasCpCorrelation fromLegacyIdealGasCp(const DIPPR::Coefficients& coeffs) {
    IdealGasCpCorrelation out;
    out.eq_form = equationFormFromLegacyType(coeffs.equation);
    out.A = coeffs.A;
    out.B = coeffs.B;
    out.C = coeffs.C;
    out.D = coeffs.D;
    out.E = coeffs.E;
    out.F = coeffs.F;
    out.G = coeffs.G;
    out.Tmin = coeffs.T_min;
    out.Tmax = coeffs.T_max;
    return out;
}

} // namespace

// =============================================================================
// Association Scheme Conversion
// =============================================================================

AssociationScheme convertScheme(pcsaft::AssociationScheme legacy) {
    switch (legacy) {
        case pcsaft::AssociationScheme::NONE:      return AssociationScheme::None;
        case pcsaft::AssociationScheme::SCHEME_1A: return AssociationScheme::Scheme1;
        case pcsaft::AssociationScheme::SCHEME_2B: return AssociationScheme::Scheme2B;
        case pcsaft::AssociationScheme::SCHEME_3B: return AssociationScheme::Scheme3B;
        case pcsaft::AssociationScheme::SCHEME_4C: return AssociationScheme::Scheme4C;
        default: return AssociationScheme::None;
    }
}

pcsaft::AssociationScheme convertScheme(AssociationScheme scheme) {
    switch (scheme) {
        case AssociationScheme::None:        return pcsaft::AssociationScheme::NONE;
        case AssociationScheme::Scheme1:     return pcsaft::AssociationScheme::SCHEME_1A;
        case AssociationScheme::Scheme2A:    return pcsaft::AssociationScheme::SCHEME_2B;  // Closest match
        case AssociationScheme::Scheme2B:    return pcsaft::AssociationScheme::SCHEME_2B;
        case AssociationScheme::Scheme3B:    return pcsaft::AssociationScheme::SCHEME_3B;
        case AssociationScheme::Scheme4C:    return pcsaft::AssociationScheme::SCHEME_4C;
        default: return pcsaft::AssociationScheme::NONE;
    }
}

// =============================================================================
// Component Conversion
// =============================================================================

Component fromLegacy(const pcsaft::Component& legacy) {
    // Create base component
    Component comp(
        legacy.getName(),
        legacy.getM(),
        legacy.getSigma(),
        legacy.getEpsilonK()
    );

    // Set CAS if available
    if (!legacy.getCAS().empty()) {
        comp = comp.withCAS(legacy.getCAS());
    }

    // Set molecular weight if available
    if (legacy.getMW() > 0) {
        comp = comp.withMW(legacy.getMW());
    }

    // Handle association parameters
    if (legacy.hasAssociation()) {
        const auto& legacy_assoc = legacy.getAssociationParams();

        // Create new association params
        AssociationParams assoc;
        assoc.scheme = convertScheme(legacy_assoc.scheme);
        assoc.epsilon_AB = legacy_assoc.epsilon_AB;
        assoc.kappa_AB = legacy_assoc.kappa_AB;
        assoc.num_sites = legacy_assoc.num_sites;

        // Need to reconstruct with association
        comp = Component(
            legacy.getName(),
            legacy.getM(),
            legacy.getSigma(),
            legacy.getEpsilonK(),
            assoc
        );

        if (!legacy.getCAS().empty()) {
            comp = comp.withCAS(legacy.getCAS());
        }
        if (legacy.getMW() > 0) {
            comp = comp.withMW(legacy.getMW());
        }
    }

    // Handle DIPPR Cp parameters
    if (legacy.hasIdealGasCp()) {
        comp = comp.withIdealGasCp(fromLegacyIdealGasCp(legacy.getIdealGasCpParams()));
    }

    return comp;
}

pcsaft::Component toLegacy(const Component& component) {
    pcsaft::Component legacy;

    legacy.setName(component.name());
    legacy.setCAS(component.cas());
    if (!component.hasPCSaftParams()) {
        throw std::runtime_error("Cannot convert component '" + component.name() + "' to legacy PC-SAFT without PC-SAFT parameters");
    }
    legacy.setM(component.m());
    legacy.setSigma(component.sigma());
    legacy.setEpsilonK(component.epsilonK());
    legacy.setMW(component.MW());

    // Handle association parameters
    if (component.isAssociating()) {
        const auto& assoc = component.associationParams();
        pcsaft::AssociationParams legacy_assoc(
            convertScheme(assoc.scheme),
            assoc.epsilon_AB,
            assoc.kappa_AB,
            assoc.num_sites
        );
        legacy.setAssociationParams(legacy_assoc);
    }

    // Handle Cp parameters
    if (component.hasIdealGasCp()) {
        legacy.setIdealGasCpParams(toLegacyIdealGasCp(component.idealGasCpCorrelation()));
    }

    return legacy;
}

// =============================================================================
// Mixture Conversion
// =============================================================================

Mixture fromLegacy(const pcsaft::Mixture& legacy) {
    int nc = legacy.getNumComponents();

    // Convert components
    std::vector<Component> components;
    components.reserve(nc);
    for (int i = 0; i < nc; ++i) {
        components.push_back(fromLegacy(legacy.getComponent(i)));
    }

    // Convert binary parameters
    BinaryParameters kij = BinaryParameters::zeros(nc);
    for (int i = 0; i < nc; ++i) {
        for (int j = i + 1; j < nc; ++j) {
            double k = legacy.getBinaryParameter(i, j);
            if (std::abs(k) > 1e-15) {
                kij.set(i, j, k);
            }
        }
    }

    return Mixture(components, kij);
}

pcsaft::Mixture toLegacy(const Mixture& mixture, const std::vector<double>& composition) {
    int nc = mixture.numComponents();

    // Convert components
    std::vector<pcsaft::Component> components;
    components.reserve(nc);
    for (int i = 0; i < nc; ++i) {
        components.push_back(toLegacy(mixture.component(i)));
    }

    // Create legacy mixture
    pcsaft::Mixture legacy(components, composition);

    // Set binary parameters
    for (int i = 0; i < nc; ++i) {
        for (int j = i + 1; j < nc; ++j) {
            double k = mixture.kij(i, j);
            if (std::abs(k) > 1e-15) {
                legacy.setBinaryParameter(i, j, k);
            }
        }
    }

    return legacy;
}

// =============================================================================
// Convenience Functions
// =============================================================================

Component getComponentByName(const std::string& name) {
    pcsaft::ComponentDatabase& db = pcsaft::ComponentDatabase::getInstance();
    pcsaft::Component legacy = db.getByName(name);
    return fromLegacy(legacy);
}

Mixture createMixture(const std::vector<std::string>& names) {
    std::vector<Component> components;
    components.reserve(names.size());
    for (const auto& name : names) {
        components.push_back(getComponentByName(name));
    }
    return Mixture(components);
}

Mixture createMixture(const std::vector<std::string>& names, const BinaryParameters& kij) {
    std::vector<Component> components;
    components.reserve(names.size());
    for (const auto& name : names) {
        components.push_back(getComponentByName(name));
    }
    return Mixture(components, kij);
}

} // namespace Core
} // namespace DMThermo
