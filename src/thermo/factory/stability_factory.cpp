/**
 * @file stability_factory.cpp
 * @brief Implementation of stability analyzer factory
 */

#include "thermo/factory/stability_factory.h"
#include "thermo/equilibrium/tpd_stability.h"
#include <cmath>
#include <stdexcept>

namespace DMThermo {
namespace Factory {

bool StabilityFactory::initialized_ = false;

void StabilityFactory::initializeBuiltinAnalyzers() {
    if (initialized_) return;
    initialized_ = true;
}

Equilibrium::Stability::StabilityAnalyzerPtr StabilityFactory::create(EOSPtr eos) {
    initializeBuiltinAnalyzers();

    // Stability analysis requires per-component critical properties for Wilson K-value initialization.
    const int nc = eos ? eos->numComponents() : 0;
    if (nc <= 0) {
        throw std::invalid_argument("StabilityFactory::create: EOS has no components");
    }

    for (int i = 0; i < nc; ++i) {
        const double Tc = eos->criticalTemperature(i);
        const double Pc = eos->criticalPressure(i);
        const double omega = eos->acentricFactor(i);
        if (!(std::isfinite(Tc) && std::isfinite(Pc) && std::isfinite(omega)) || Tc <= 0.0 || Pc <= 0.0) {
            throw std::invalid_argument("StabilityFactory::create: EOS is missing required critical properties (Tc/Pc/omega)");
        }
    }

    return std::make_shared<Equilibrium::Stability::TPDStabilityAnalyzer>(eos);
}

Equilibrium::Stability::StabilityAnalyzerPtr StabilityFactory::createTPD(EOSPtr eos) {
    return std::make_shared<Equilibrium::Stability::TPDStabilityAnalyzer>(eos);
}

Equilibrium::Stability::StabilityAnalyzerPtr StabilityFactory::createWithSpinodal(EOSPtr eos) {
    // Same as default for now (spinodal is not a required dependency for TPD-based stability).
    return create(eos);
}

Equilibrium::Stability::StabilityAnalyzerPtr StabilityFactory::createForLLE(EOSPtr eos) {
    // Same as default for now - could customize trial compositions
    return create(eos);
}

} // namespace Factory
} // namespace DMThermo
