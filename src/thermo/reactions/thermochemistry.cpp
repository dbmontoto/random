/**
 * @file thermochemistry.cpp
 * @brief Implementation of databank-backed thermochemistry helpers.
 */

#include "thermo/reactions/thermochemistry.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "thermo/data/dippr_adapter.h"
#include "thermo/dippr/dippr.h"

namespace DMThermo {
namespace Reactions {

namespace Thermochemistry {

std::vector<double> DIPPR::mu0GFOR(
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    Data::Diagnostics* diag)
{
    const int nc = mixture.numComponents();
    std::vector<double> mu0(static_cast<std::size_t>(nc), 0.0);

    for (int i = 0; i < nc; ++i) {
        const auto& c = mixture.component(i);
        const std::string& cas = c.cas();
        const std::string& name = c.name();

        const Data::DipprRecord* rec = nullptr;
        if (!cas.empty()) {
            rec = databanks.dippr.findByCAS(cas);
        } else if (!name.empty()) {
            rec = databanks.dippr.findByNameUnique(name, diag);
        }

        if (!rec) {
            throw std::runtime_error(
                "Thermochemistry::DIPPR::mu0GFOR: DIPPR record not found for component '" +
                (!cas.empty() ? cas : name) + "'");
        }
        if (!rec->GFOR.has_value() || !std::isfinite(*rec->GFOR)) {
            throw std::runtime_error(
                "Thermochemistry::DIPPR::mu0GFOR: missing/invalid GFOR for component '" +
                (!rec->key.cas.empty() ? rec->key.cas : rec->key.name) + "'");
        }

        // DIPPR stores GFOR in J/kmol; convert to J/mol.
        mu0[static_cast<std::size_t>(i)] = (*rec->GFOR) / 1000.0;
    }

    return mu0;
}

std::vector<double> DIPPR::mu0FormT(
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    double T,
    double T_ref,
    Data::Diagnostics* diag)
{
    if (!(std::isfinite(T) && T > 0.0)) {
        throw std::invalid_argument("Thermochemistry::DIPPR::mu0FormT: T must be finite and > 0");
    }
    if (!(std::isfinite(T_ref) && T_ref > 0.0)) {
        throw std::invalid_argument("Thermochemistry::DIPPR::mu0FormT: T_ref must be finite and > 0");
    }

    const int nc = mixture.numComponents();
    std::vector<double> mu0(static_cast<std::size_t>(nc), 0.0);

    for (int i = 0; i < nc; ++i) {
        const auto& c = mixture.component(i);
        const std::string& cas = c.cas();
        const std::string& name = c.name();

        const Data::DipprRecord* rec = nullptr;
        if (!cas.empty()) {
            rec = databanks.dippr.findByCAS(cas);
        } else if (!name.empty()) {
            rec = databanks.dippr.findByNameUnique(name, diag);
        }

        if (!rec) {
            throw std::runtime_error(
                "Thermochemistry::DIPPR::mu0FormT: DIPPR record not found for component '" +
                (!cas.empty() ? cas : name) + "'");
        }

        if (!rec->GFOR.has_value() || !std::isfinite(*rec->GFOR)) {
            throw std::runtime_error(
                "Thermochemistry::DIPPR::mu0FormT: missing/invalid GFOR for component '" +
                (!rec->key.cas.empty() ? rec->key.cas : rec->key.name) + "'");
        }
        if (!rec->HFOR.has_value() || !std::isfinite(*rec->HFOR)) {
            throw std::runtime_error(
                "Thermochemistry::DIPPR::mu0FormT: missing/invalid HFOR for component '" +
                (!rec->key.cas.empty() ? rec->key.cas : rec->key.name) + "'");
        }

        // Convert from J/kmol to J/mol.
        const double g_ref = (*rec->GFOR) / 1000.0;
        const double h_ref = (*rec->HFOR) / 1000.0;
        const double s_ref = (h_ref - g_ref) / T_ref; // J/mol/K, consistent with formation basis

        // Require ICP (ideal-gas Cp) to compute Δh and Δs.
        const auto params = Data::toDipprParameters(*rec, diag);
        if (!params || !params->has_ideal_gas_cp) {
            throw std::runtime_error(
                "Thermochemistry::DIPPR::mu0FormT: missing ICP (ideal-gas Cp) for component '" +
                (!rec->key.cas.empty() ? rec->key.cas : rec->key.name) + "'");
        }

        ::DIPPR::Calculator calc(*params);
        const double dh = calc.idealGasEnthalpy(T, T_ref); // J/mol
        const double ds = calc.idealGasEntropy(T, T_ref);  // J/mol/K

        const double h = h_ref + dh;
        const double s = s_ref + ds;
        mu0[static_cast<std::size_t>(i)] = h - T * s;
    }

    return mu0;
}

std::vector<double> mu0DB(
    const std::string& model,
    const Data::Databanks& databanks,
    const Core::Mixture& mixture,
    double T,
    double T_ref,
    Data::Diagnostics* diag)
{
    if (model == Models::DIPPR_GFOR) {
        return DIPPR::mu0GFOR(databanks, mixture, diag);
    }
    if (model == Models::DIPPR_MU0_FORMATION_T) {
        return DIPPR::mu0FormT(databanks, mixture, T, T_ref, diag);
    }

    throw std::invalid_argument("Thermochemistry::mu0DB: unknown model '" + model + "'");
}

} // namespace Thermochemistry

} // namespace Reactions
} // namespace DMThermo
