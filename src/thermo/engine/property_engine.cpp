#include "thermo/engine/property_engine.h"

#include "thermo/dippr/dippr.h"
#include "thermo/data/component_resolver.h"
#include "thermo/data/dippr_adapter.h"
#include "thermo/core/constants.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <limits>

namespace DMThermo {
namespace Engine {

namespace {

std::string formatDouble(double v) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(6) << v;
    return oss.str();
}

std::string formatRangeK(const DIPPR::Coefficients& c) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << "[" << c.T_min << ", " << c.T_max << "] K";
    return oss.str();
}

std::string labelFor(const Data::DipprRecord& r) {
    std::string label = r.key.name;
    if (!r.key.cas.empty()) {
        label += " (" + r.key.cas + ")";
    } else if (r.key.chemid.has_value()) {
        label += " (CHEMID " + std::to_string(r.key.chemid.value()) + ")";
    }
    return label;
}

double dippr9G1LiquidThermalConductivityCorrection(double Tr, double Pr) {
    // DIPPR procedure 9G-1 (COCO/TEA documentation):
    // k_hp = F(Tr, Pr) * k_low
    // F = 0.63*Tr^1.2*Pr/(30+Pr) + 0.98 + 0.0079*Pr*Tr^1.4
    return 0.63 * std::pow(Tr, 1.2) * Pr / (30.0 + Pr) + 0.98 + 0.0079 * Pr * std::pow(Tr, 1.4);
}

constexpr double DIPPR_K_HP_THRESHOLD_PA = 35.0e5;

} // namespace

PropertyEngine::PropertyEngine(EOSConstPtr eos, Core::Mixture mixture)
    : eos_(std::move(eos))
    , mixture_(std::move(mixture))
{
    if (!eos_) {
        throw std::invalid_argument("PropertyEngine: eos is null");
    }
}

Properties::ThermoProperties PropertyEngine::calculateTV(double T, double rho, const std::vector<double>& x) const {
    return Properties::calculateTV(*eos_, mixture_, T, rho, x);
}

Properties::ThermoProperties PropertyEngine::calculatePureLiquidDIPPR(
    const Data::Databanks& databanks,
    const std::string& identifier,
    double T,
    double P,
    const DipprLiquidRequirements& req,
    Data::Diagnostics* diag) const
{
    Properties::ThermoProperties out;
    out.T = T;
    out.P = P;
    out.model = "DIPPR (CSV)";

    const Data::ComponentResolver resolver(databanks);
    const auto resolved = resolver.resolve(identifier, diag);
    if (!resolved.dippr) {
        out.success = false;
        out.message = "No DIPPR record available for '" + identifier + "'";
        return out;
    }

    const std::string record_label = labelFor(*resolved.dippr);

    auto params = Data::toDipprParameters(*resolved.dippr, diag);
    if (!params) {
        out.success = false;
        out.message = "No supported DIPPR correlations available for '" + record_label + "'";
        return out;
    }

    try {
        DIPPR::Calculator calc(*params);
        const auto& p = calc.getParameters();

        if (!p.has_liquid_density && req.require_liquid_density) {
            throw std::runtime_error("DIPPR liquid density correlation required but missing for '" + record_label + "'");
        }
        if (!p.has_liquid_heat_capacity && req.require_liquid_heat_capacity) {
            throw std::runtime_error("DIPPR liquid heat capacity correlation required but missing for '" + record_label + "'");
        }
        if (!p.has_liquid_viscosity && req.require_liquid_viscosity) {
            throw std::runtime_error("DIPPR liquid viscosity correlation required but missing for '" + record_label + "'");
        }
        if (!p.has_liquid_thermal_conductivity && req.require_liquid_thermal_conductivity) {
            throw std::runtime_error("DIPPR liquid thermal conductivity correlation required but missing for '" + record_label + "'");
        }

        if (p.has_liquid_density) {
            try {
                out.rho = calc.liquidDensity(T);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "DIPPR liquid density evaluation failed for '" + record_label + "' at T=" + formatDouble(T) +
                    " K (valid " + formatRangeK(p.liquid_density) + "): " + e.what()
                );
            }
            out.v = (out.rho > 0.0) ? 1.0 / out.rho : 0.0;
        } else if (diag) {
            diag->warn(
                "DIPPR liquid density correlation not available for '" + record_label + "'",
                Config::ErrorCategory::MissingData,
                "DIPPR_MISSING_LIQUID_DENSITY");
        }

        const double R = Constants::GAS_CONSTANT;
        out.Z = (out.v > 0.0) ? (P * out.v / (R * T)) : 0.0;

        if (p.has_liquid_heat_capacity) {
            try {
                out.Cp = calc.liquidCp(T);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "DIPPR liquid Cp evaluation failed for '" + record_label + "' at T=" + formatDouble(T) +
                    " K (valid " + formatRangeK(p.liquid_heat_capacity) + "): " + e.what()
                );
            }
            // For incompressible liquid correlations, Cp≈Cv is a reasonable approximation.
            out.Cv = out.Cp;
            try {
                out.H = calc.liquidEnthalpy(T);
                out.S = calc.liquidEntropy(T);
            } catch (const std::exception& e) {
                if (req.require_liquid_heat_capacity) {
                    throw std::runtime_error(
                        "DIPPR liquid H/S evaluation failed for '" + record_label + "' at T=" + formatDouble(T) +
                        " K (Cp valid " + formatRangeK(p.liquid_heat_capacity) + "): " + e.what()
                    );
                }
                if (diag) {
                    diag->warn(
                        "DIPPR liquid H/S not computed for '" + record_label + "' at T=" + formatDouble(T) +
                        " K (Cp valid " + formatRangeK(p.liquid_heat_capacity) + "): " + e.what(),
                        Config::ErrorCategory::NumericalFailure,
                        "DIPPR_LIQUID_HS_EVALUATION_FAILED"
                    );
                }
            }
        } else if (diag) {
            diag->warn(
                "DIPPR liquid heat capacity correlation not available for '" + record_label + "'; H/S not computed",
                Config::ErrorCategory::MissingData,
                "DIPPR_MISSING_LIQUID_HEAT_CAPACITY");
        }

        if (p.has_liquid_viscosity) {
            try {
                out.mu = calc.liquidViscosity(T);
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "DIPPR liquid viscosity evaluation failed for '" + record_label + "' at T=" + formatDouble(T) +
                    " K (valid " + formatRangeK(p.liquid_viscosity) + "): " + e.what()
                );
            }
        } else if (diag) {
            diag->warn(
                "DIPPR liquid viscosity correlation not available for '" + record_label + "'",
                Config::ErrorCategory::MissingData,
                "DIPPR_MISSING_LIQUID_VISCOSITY");
        }

        if (p.has_liquid_thermal_conductivity) {
            try {
                double k_low = calc.liquidThermalConductivity(T);
                if (enable_k_9g1_ && (P > DIPPR_K_HP_THRESHOLD_PA)) {
                    if (!(std::isfinite(p.Tc) && p.Tc > 0.0 && std::isfinite(p.Pc) && p.Pc > 0.0)) {
                        if (req.require_liquid_thermal_conductivity) {
                            throw std::runtime_error("DIPPR 9G-1 high-pressure correction requires Tc/Pc in databank for '" + record_label + "'");
                        }
                        if (diag) {
                            diag->warn(
                                "Skipping DIPPR 9G-1 high-pressure correction for '" + record_label + "' (missing Tc/Pc)",
                                Config::ErrorCategory::MissingData,
                                "DIPPR_9G1_MISSING_TC_PC");
                        }
                    } else {
                        const double Tr = T / p.Tc;
                        const double Pr = P / p.Pc;
                        const double factor = dippr9G1LiquidThermalConductivityCorrection(Tr, Pr);
                        if (!(std::isfinite(factor) && factor > 0.0)) {
                            if (req.require_liquid_thermal_conductivity) {
                                throw std::runtime_error("DIPPR 9G-1 high-pressure correction produced invalid correction factor for '" + record_label + "'");
                            }
                            if (diag) {
                                diag->warn(
                                    "Skipping DIPPR 9G-1 high-pressure correction for '" + record_label + "' (invalid correction factor)",
                                    Config::ErrorCategory::NumericalFailure,
                                    "DIPPR_9G1_INVALID_FACTOR");
                            }
                        } else {
                            k_low *= factor;
                        }
                    }
                }
                out.k = k_low;
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "DIPPR liquid thermal conductivity evaluation failed for '" + record_label + "' at T=" + formatDouble(T) +
                    " K (valid " + formatRangeK(p.liquid_thermal_conductivity) + "): " + e.what()
                );
            }
        } else if (diag) {
            diag->warn(
                "DIPPR liquid thermal conductivity correlation not available for '" + record_label + "'",
                Config::ErrorCategory::MissingData,
                "DIPPR_MISSING_LIQUID_THERMAL_CONDUCTIVITY");
        }

        out.U = out.H - P * out.v;
        out.G = out.H - T * out.S;
        out.A = out.U - T * out.S;

        out.success = true;
        return out;
    } catch (const std::exception& e) {
        out.success = false;
        out.message = e.what();
        return out;
    }
}

double PropertyEngine::calculatePureVaporPressureDIPPR(
    const Data::Databanks& databanks,
    const std::string& identifier,
    double T,
    Data::Diagnostics* diag) const
{
    const Data::ComponentResolver resolver(databanks);
    const auto resolved = resolver.resolve(identifier, diag);
    if (!resolved.dippr) {
        throw std::runtime_error("No DIPPR record available for '" + identifier + "'");
    }

    const std::string record_label = labelFor(*resolved.dippr);
    auto params = Data::toDipprParameters(*resolved.dippr, diag);
    if (!params) {
        throw std::runtime_error("No supported DIPPR correlations available for '" + record_label + "'");
    }

    DIPPR::Calculator calc(*params);
    const auto& p = calc.getParameters();
    if (!p.has_vapor_pressure) {
        throw std::runtime_error("DIPPR vapor pressure correlation required but missing for '" + record_label + "'");
    }
    return calc.vaporPressure(T);
}

double PropertyEngine::calculatePureHeatOfVaporizationDIPPR(
    const Data::Databanks& databanks,
    const std::string& identifier,
    double T,
    Data::Diagnostics* diag) const
{
    const Data::ComponentResolver resolver(databanks);
    const auto resolved = resolver.resolve(identifier, diag);
    if (!resolved.dippr) {
        throw std::runtime_error("No DIPPR record available for '" + identifier + "'");
    }

    const std::string record_label = labelFor(*resolved.dippr);
    auto params = Data::toDipprParameters(*resolved.dippr, diag);
    if (!params) {
        throw std::runtime_error("No supported DIPPR correlations available for '" + record_label + "'");
    }

    DIPPR::Calculator calc(*params);
    const auto& p = calc.getParameters();
    if (!p.has_heat_of_vaporization) {
        throw std::runtime_error("DIPPR heat of vaporization correlation required but missing for '" + record_label + "'");
    }
    return calc.heatOfVaporization(T);
}

Properties::ThermoProperties PropertyEngine::calculateLiquidDIPPRIdealMix(
    const Data::Databanks& databanks,
    double T,
    double P,
    const std::vector<double>& x,
    const DipprLiquidRequirements& req,
    Data::Diagnostics* diag) const
{
    Properties::ThermoProperties out;
    out.T = T;
    out.P = P;
    out.model = "DIPPR (CSV ideal mix)";

    try {
        if (!mixture_.isValidComposition(x)) {
            throw std::invalid_argument("Invalid composition vector");
        }

        const int nc = mixture_.numComponents();
        if (static_cast<int>(x.size()) != nc) {
            throw std::invalid_argument("Composition size mismatch");
        }

        const Data::ComponentResolver resolver(databanks);

        std::vector<double> rho_i(nc, 0.0);
        std::vector<double> mu_i(nc, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> k_i(nc, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> Tc_i(nc, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> Pc_i(nc, std::numeric_limits<double>::quiet_NaN());
        bool all_have_cp = true;
        bool all_have_density = true;
        bool all_have_mu = true;
        bool all_have_k = true;
        bool all_have_critical = true;
        std::vector<std::string> missing_density;
        std::vector<std::string> density_failures;
        std::vector<std::string> missing_cp;
        std::vector<std::string> cp_failures;
        std::vector<std::string> missing_mu;
        std::vector<std::string> mu_failures;
        std::vector<std::string> missing_k;
        std::vector<std::string> k_failures;

        double v_mix = 0.0;  // molar volume [m^3/mol]
        double Cp_mix = 0.0;
        double H_mix = 0.0;
        double S_mix_pure = 0.0;
        double ln_mu_mix = 0.0;
        double mu_lin = 0.0;

        for (int i = 0; i < nc; ++i) {
            if (x[i] == 0.0) continue;

            const auto& c = mixture_.component(i);
            const std::string id = !c.cas().empty() ? c.cas() : c.name();
            if (id.empty()) {
                throw std::runtime_error("Cannot resolve identifier for component at index " + std::to_string(i) + " (missing CAS/name)");
            }

            const auto resolved = resolver.resolve(id, diag);
            if (!resolved.dippr) {
                throw std::runtime_error("No DIPPR record available for component '" + id + "'");
            }

            const std::string record_label = labelFor(*resolved.dippr);

            auto params = Data::toDipprParameters(*resolved.dippr, diag);
            if (!params) {
                throw std::runtime_error("No supported DIPPR correlations available for component '" + record_label + "'");
            }

            DIPPR::Calculator calc(*params);
            const auto& p = calc.getParameters();

            // Critical properties for 9G-1 high-pressure correction (optional unless required by P threshold).
            if (p.Tc > 0.0 && std::isfinite(p.Tc) && p.Pc > 0.0 && std::isfinite(p.Pc)) {
                Tc_i[i] = p.Tc;
                Pc_i[i] = p.Pc;
            } else {
                all_have_critical = false;
            }

            if (p.has_liquid_density) {
                try {
                    rho_i[i] = calc.liquidDensity(T);
                    if (!(std::isfinite(rho_i[i]) && rho_i[i] > 0.0)) {
                        throw std::runtime_error("DIPPR returned invalid liquid density");
                    }
                    v_mix += x[i] / rho_i[i];
                } catch (const std::exception& e) {
                    all_have_density = false;
                    density_failures.push_back(
                        record_label + " @ T=" + formatDouble(T) + " K (valid " + formatRangeK(p.liquid_density) + "): " + e.what()
                    );
                }
            } else {
                all_have_density = false;
                missing_density.push_back(record_label);
            }

            if (p.has_liquid_heat_capacity) {
                try {
                    const double Cp_i = calc.liquidCp(T);
                    const double H_i = calc.liquidEnthalpy(T);
                    const double S_i = calc.liquidEntropy(T);
                    Cp_mix += x[i] * Cp_i;
                    H_mix += x[i] * H_i;
                    S_mix_pure += x[i] * S_i;
                } catch (const std::exception& e) {
                    all_have_cp = false;
                    cp_failures.push_back(
                        record_label + " @ T=" + formatDouble(T) + " K (valid " + formatRangeK(p.liquid_heat_capacity) + "): " + e.what()
                    );
                }
            } else {
                all_have_cp = false;
                missing_cp.push_back(record_label);
            }

            if (p.has_liquid_viscosity) {
                try {
                    const double mu_val = calc.liquidViscosity(T);
                    if (!(std::isfinite(mu_val) && mu_val > 0.0)) {
                        throw std::runtime_error("DIPPR returned invalid liquid viscosity");
                    }
                    mu_i[i] = mu_val;
                    if (mu_mixing_ == LiquidViscosityMixing::SimpleLinear) {
                        mu_lin += x[i] * mu_val;
                    } else {
                        ln_mu_mix += x[i] * std::log(mu_val);
                    }
                } catch (const std::exception& e) {
                    all_have_mu = false;
                    mu_failures.push_back(
                        record_label + " @ T=" + formatDouble(T) + " K (valid " + formatRangeK(p.liquid_viscosity) + "): " + e.what()
                    );
                }
            } else {
                all_have_mu = false;
                missing_mu.push_back(record_label);
            }

            if (p.has_liquid_thermal_conductivity) {
                try {
                    const double k_val = calc.liquidThermalConductivity(T);
                    if (!(std::isfinite(k_val) && k_val > 0.0)) {
                        throw std::runtime_error("DIPPR returned invalid liquid thermal conductivity");
                    }
                    k_i[i] = k_val;
                } catch (const std::exception& e) {
                    all_have_k = false;
                    k_failures.push_back(
                        record_label + " @ T=" + formatDouble(T) + " K (valid " + formatRangeK(p.liquid_thermal_conductivity) + "): " + e.what()
                    );
                }
            } else {
                all_have_k = false;
                missing_k.push_back(record_label);
            }
        }

    if (req.require_liquid_density && !all_have_density) {
        std::string msg = "DIPPR mixture liquid density could not be computed:";
        if (!missing_density.empty()) {
            msg += " missing liquid density for {";
            for (size_t k = 0; k < missing_density.size(); ++k) {
                if (k) msg += ", ";
                msg += missing_density[k];
            }
            msg += "}";
        }
        if (!density_failures.empty()) {
            msg += " density evaluation failed for {";
            for (size_t k = 0; k < density_failures.size(); ++k) {
                if (k) msg += "; ";
                msg += density_failures[k];
            }
            msg += "}";
        }
        throw std::runtime_error(msg);
    }

    if (req.require_liquid_viscosity && !all_have_mu) {
        std::string msg = "DIPPR mixture liquid viscosity could not be computed:";
        if (!missing_mu.empty()) {
            msg += " missing liquid viscosity for {";
            for (size_t k = 0; k < missing_mu.size(); ++k) {
                if (k) msg += ", ";
                msg += missing_mu[k];
            }
            msg += "}";
        }
        if (!mu_failures.empty()) {
            msg += " viscosity evaluation failed for {";
            for (size_t k = 0; k < mu_failures.size(); ++k) {
                if (k) msg += "; ";
                msg += mu_failures[k];
            }
            msg += "}";
        }
        throw std::runtime_error(msg);
    }

    if (req.require_liquid_thermal_conductivity && !all_have_k) {
        std::string msg = "DIPPR mixture liquid thermal conductivity could not be computed:";
        if (!missing_k.empty()) {
            msg += " missing liquid thermal conductivity for {";
            for (size_t k = 0; k < missing_k.size(); ++k) {
                if (k) msg += ", ";
                msg += missing_k[k];
            }
            msg += "}";
        }
        if (!k_failures.empty()) {
            msg += " thermal conductivity evaluation failed for {";
            for (size_t k = 0; k < k_failures.size(); ++k) {
                if (k) msg += "; ";
                msg += k_failures[k];
            }
            msg += "}";
        }
        throw std::runtime_error(msg);
    }

    if (all_have_density) {
        if (!(std::isfinite(v_mix) && v_mix > 0.0)) {
            throw std::runtime_error("DIPPR mixture liquid density computation produced invalid results");
        }

        out.v = v_mix;
        out.rho = 1.0 / v_mix;
        out.Z = (P * v_mix) / (Constants::GAS_CONSTANT * T);
    } else if (diag) {
        diag->warn(
            "DIPPR mixture liquid density not available (one or more components missing/invalid density)",
            Config::ErrorCategory::MissingData,
            "DIPPR_MIXTURE_DENSITY_UNAVAILABLE");
    }

        if (all_have_mu) {
            if (mu_mixing_ == LiquidViscosityMixing::SimpleLinear) {
                if (!(std::isfinite(mu_lin) && mu_lin > 0.0)) {
                    throw std::runtime_error("DIPPR mixture liquid viscosity computation produced invalid results");
                }
                out.mu = mu_lin;
            } else {
                out.mu = std::exp(ln_mu_mix);
            }
        } else if (diag) {
            diag->warn(
                "DIPPR mixture liquid viscosity not available (one or more components missing/invalid viscosity)",
                Config::ErrorCategory::MissingData,
                "DIPPR_MIXTURE_VISCOSITY_UNAVAILABLE");
        }

        if (all_have_k) {
            auto apply9g1Mix = [&](double k_low) -> double {
                if (!enable_k_9g1_ || !(P > DIPPR_K_HP_THRESHOLD_PA)) return k_low;

                double Tc_mix = 0.0;
                double Pc_mix = 0.0;
                for (int i = 0; i < nc; ++i) {
                    if (x[i] == 0.0) continue;
                    const double tci = Tc_i[i];
                    const double pci = Pc_i[i];
                    if (!(std::isfinite(tci) && tci > 0.0 && std::isfinite(pci) && pci > 0.0)) {
                        if (req.require_liquid_thermal_conductivity) {
                            throw std::runtime_error("DIPPR 9G-1 high-pressure correction requires Tc/Pc for all components with x_i>0");
                        }
                        if (diag) {
                            diag->warn(
                                "Skipping DIPPR 9G-1 high-pressure correction (missing Tc/Pc in databank)",
                                Config::ErrorCategory::MissingData,
                                "DIPPR_9G1_MISSING_TC_PC");
                        }
                        return k_low;
                    }
                    Tc_mix += x[i] * tci;
                    Pc_mix += x[i] * pci;
                }

                if (!(Tc_mix > 0.0 && Pc_mix > 0.0)) return k_low;
                const double Tr = T / Tc_mix;
                const double Pr = P / Pc_mix;
                if (!(std::isfinite(Tr) && std::isfinite(Pr) && Tr > 0.0 && Pr > 0.0)) return k_low;
                const double factor = dippr9G1LiquidThermalConductivityCorrection(Tr, Pr);
                if (!(std::isfinite(factor) && factor > 0.0)) return k_low;
                return k_low * factor;
            };

            if (k_mixing_ == LiquidThermalConductivityMixing::SimpleLinear) {
                double k_lin = 0.0;
                for (int i = 0; i < nc; ++i) {
                    if (x[i] == 0.0) continue;
                    const double ki = k_i[i];
                    if (!(std::isfinite(ki) && ki > 0.0)) {
                        throw std::runtime_error("DIPPR mixture liquid thermal conductivity mixing encountered invalid k for a component");
                    }
                    k_lin += x[i] * ki;
                }
                out.k = apply9g1Mix(k_lin);
            } else {
                // DIPPR 9I requires volume fractions, hence density for all components with x_i>0.
                if (!all_have_density) {
                    if (req.require_liquid_thermal_conductivity) {
                        throw std::runtime_error("DIPPR mixture liquid thermal conductivity mixing (9I) requires liquid density for all components");
                    }
                    if (diag) {
                        diag->warn(
                            "DIPPR mixture liquid thermal conductivity not available (9I mixing requires density for all components)",
                            Config::ErrorCategory::MissingData,
                            "DIPPR_MIXTURE_THERMAL_CONDUCTIVITY_UNAVAILABLE");
                    }
                } else {
                    std::vector<double> phi(nc, 0.0);
                    for (int i = 0; i < nc; ++i) {
                        if (x[i] == 0.0) continue;
                        if (!(rho_i[i] > 0.0)) {
                            throw std::runtime_error(
                                "DIPPR mixture liquid thermal conductivity mixing encountered invalid density for component index " +
                                std::to_string(i)
                            );
                        }
                        phi[i] = (x[i] / rho_i[i]) / v_mix;
                    }

                    double k_mix_9i = 0.0;
                    for (int i = 0; i < nc; ++i) {
                        if (x[i] == 0.0) continue;
                        for (int j = 0; j < nc; ++j) {
                            if (x[j] == 0.0) continue;
                            const double ki = k_i[i];
                            const double kj = k_i[j];
                            if (!(std::isfinite(ki) && ki > 0.0 && std::isfinite(kj) && kj > 0.0)) {
                                throw std::runtime_error("DIPPR mixture liquid thermal conductivity mixing encountered invalid k for a component");
                            }
                            const double kij = (i == j) ? ki : (2.0 / (1.0 / ki + 1.0 / kj));
                            k_mix_9i += phi[i] * phi[j] * kij;
                        }
                    }

                    out.k = apply9g1Mix(k_mix_9i);
                }
            }
        } else if (diag) {
            diag->warn(
                "DIPPR mixture liquid thermal conductivity not available (one or more components missing/invalid thermal conductivity)",
                Config::ErrorCategory::MissingData,
                "DIPPR_MIXTURE_THERMAL_CONDUCTIVITY_UNAVAILABLE");
        }

    if (req.require_liquid_heat_capacity && !all_have_cp) {
        std::string msg = "DIPPR mixture liquid heat capacity could not be computed:";
        if (!missing_cp.empty()) {
            msg += " missing liquid heat capacity for {";
            for (size_t k = 0; k < missing_cp.size(); ++k) {
                if (k) msg += ", ";
                msg += missing_cp[k];
            }
            msg += "}";
        }
        if (!cp_failures.empty()) {
            msg += " Cp/H/S evaluation failed for {";
            for (size_t k = 0; k < cp_failures.size(); ++k) {
                if (k) msg += "; ";
                msg += cp_failures[k];
            }
            msg += "}";
        }
        throw std::runtime_error(msg);
    }

    if (all_have_cp) {
        double mixing_entropy = 0.0;
        for (int i = 0; i < nc; ++i) {
            if (x[i] > 0.0) {
                mixing_entropy -= Constants::GAS_CONSTANT * x[i] * std::log(x[i]);
            }
        }

        out.Cp = Cp_mix;
        // For incompressible liquid correlations, Cp≈Cv is a reasonable approximation.
        out.Cv = out.Cp;
        out.H = H_mix;
        out.S = S_mix_pure + mixing_entropy;
        if (out.v > 0.0) {
            out.U = out.H - P * out.v;
            out.G = out.H - T * out.S;
            out.A = out.U - T * out.S;
        }
    } else if (diag) {
        diag->warn(
            "DIPPR liquid heat capacity correlation not available for one or more components; H/S not computed",
            Config::ErrorCategory::MissingData,
            "DIPPR_MIXTURE_HEAT_CAPACITY_UNAVAILABLE");
    }

        out.success = true;
        return out;
    } catch (const std::exception& e) {
        out.success = false;
        out.message = e.what();
        return out;
    }
}

} // namespace Engine
} // namespace DMThermo
