#include "thermo/factory/mixture_factory.h"

#include "thermo/data/component_resolver.h"
#include "thermo/data/csv.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace DMThermo {
namespace Factory {

namespace {

std::optional<Core::AssociationScheme> parseAssociationScheme(const std::string& raw) {
    std::string s = Data::toLower(Data::trim(raw));
    if (s.empty() || s == "none" || s == "0") return Core::AssociationScheme::None;

    // Accept both short ("2b") and enum-ish ("scheme_2b") spellings.
    if (s == "1a" || s == "scheme_1a") return Core::AssociationScheme::Scheme1;
    if (s == "2a" || s == "scheme_2a") return Core::AssociationScheme::Scheme2A;
    if (s == "2b" || s == "scheme_2b") return Core::AssociationScheme::Scheme2B;
    if (s == "3b" || s == "scheme_3b") return Core::AssociationScheme::Scheme3B;
    if (s == "4c" || s == "scheme_4c") return Core::AssociationScheme::Scheme4C;

    return std::nullopt;
}

bool isSupportedDipprEquationForm(int eq_form) {
    switch (eq_form) {
        case 100:
        case 101:
        case 102:
        case 105:
        case 106:
        case 107:
        case 114:
        case 123:
            return true;
        default:
            return false;
    }
}

Core::IdealGasCpCorrelation convertIdealGasCpToPerMol(const Data::Correlation& c) {
    if (!isSupportedDipprEquationForm(c.eq_form)) {
        throw std::runtime_error("Unsupported ideal-gas Cp equation form " + std::to_string(c.eq_form));
    }

    // DIPPR databanks commonly store Cp correlations in J/(kmol*K). Convert to J/(mol*K).
    // EQ_107 uses C and E as temperature parameters (K) and must not be scaled.
    double A = c.A;
    double B = c.B;
    double C = c.C;
    double D = c.D;
    double E = c.E;
    double F = c.F;
    double G = c.G;

    if (c.eq_form == 107) {
        A /= 1000.0;
        B /= 1000.0;
        D /= 1000.0;
        F /= 1000.0;
        G /= 1000.0;
    } else {
        A /= 1000.0;
        B /= 1000.0;
        C /= 1000.0;
        D /= 1000.0;
        E /= 1000.0;
        F /= 1000.0;
        G /= 1000.0;
    }

    Core::IdealGasCpCorrelation coeffs;
    coeffs.eq_form = c.eq_form;
    coeffs.A = A;
    coeffs.B = B;
    coeffs.C = C;
    coeffs.D = D;
    coeffs.E = E;
    coeffs.F = F;
    coeffs.G = G;
    coeffs.Tmin = c.Tmin.value_or(0.0);
    coeffs.Tmax = c.Tmax.value_or(1000.0);
    if (!coeffs.isValid()) {
        throw std::runtime_error("Invalid ideal-gas Cp correlation (ICP_*) for eq_form " + std::to_string(c.eq_form));
    }
    return coeffs;
}

} // namespace

Core::Component buildComponentFromDatabanks(
    const Data::Databanks& databanks,
    const std::string& identifier,
    const MixtureBuildOptions& options,
    Data::Diagnostics* diag)
{
    const Data::ComponentResolver resolver(databanks);
    const auto resolved = resolver.resolve(identifier, diag);

    const Data::DipprRecord* dippr = resolved.dippr;
    const Data::PCSaftRecord* pcsaft = resolved.pcsaft;

    if (!dippr && !pcsaft) {
        throw std::runtime_error("Component '" + identifier + "' not found in DIPPR or PC-SAFT databanks");
    }

    const std::string name = !resolved.key.name.empty() ? resolved.key.name : identifier;
    const std::string cas = resolved.key.cas;

    std::optional<double> m = pcsaft ? pcsaft->m : std::nullopt;
    std::optional<double> sigma = pcsaft ? pcsaft->sigma : std::nullopt;
    std::optional<double> eps_k = pcsaft ? pcsaft->epsilon_k : std::nullopt;

    Core::Component comp = [&]() -> Core::Component {
        if (m && sigma && eps_k && *m > 0.0 && *sigma > 0.0 && *eps_k > 0.0) {
            const bool has_assoc_params =
                pcsaft &&
                pcsaft->kappa_ab && pcsaft->epsilon_k_ab &&
                (*pcsaft->kappa_ab > 0.0) && (*pcsaft->epsilon_k_ab > 0.0);

            if (has_assoc_params) {
                if (!pcsaft->assoc_scheme || !pcsaft->assoc_num_sites) {
                    throw std::runtime_error("PC-SAFT association parameters present but ASSOC_SCHEME/ASSOC_NUM_SITES missing for component '" + identifier + "'");
                }
                const auto scheme = parseAssociationScheme(*pcsaft->assoc_scheme);
                if (!scheme) {
                    throw std::runtime_error("Unknown ASSOC_SCHEME '" + *pcsaft->assoc_scheme + "' for component '" + identifier + "'");
                }
                if (*pcsaft->assoc_num_sites <= 0) {
                    throw std::runtime_error("ASSOC_NUM_SITES must be positive for component '" + identifier + "'");
                }
                Core::AssociationParams assoc;
                assoc.scheme = *scheme;
                assoc.epsilon_AB = *pcsaft->epsilon_k_ab;
                assoc.kappa_AB = *pcsaft->kappa_ab;
                assoc.num_sites = *pcsaft->assoc_num_sites;
                return Core::Component(name, *m, *sigma, *eps_k, assoc);
            }

            if ((pcsaft && (pcsaft->assoc_scheme || pcsaft->assoc_num_sites)) &&
                !(pcsaft->kappa_ab && pcsaft->epsilon_k_ab && (*pcsaft->kappa_ab > 0.0) && (*pcsaft->epsilon_k_ab > 0.0)))
            {
                throw std::runtime_error("ASSOC_SCHEME/ASSOC_NUM_SITES provided but KAPPA_AB/EPSILON_K_AB missing for component '" + identifier + "'");
            }

            return Core::Component(name, *m, *sigma, *eps_k);
        }
        if (options.require_pcsaft_params) {
            throw std::runtime_error("PC-SAFT parameters missing for component '" + identifier + "'");
        }
        if (diag) {
            diag->warn("PC-SAFT parameters missing for '" + identifier + "'; building component without PC-SAFT params");
        }
        return Core::Component(name);
    }();

    if (!cas.empty()) {
        comp = comp.withCAS(cas);
    }

    // Molecular weight (optional)
    if (dippr && dippr->MW) {
        comp = comp.withMW(*dippr->MW);
    } else if (pcsaft && pcsaft->MW) {
        comp = comp.withMW(*pcsaft->MW);
    }

    // Critical properties (optional unless caller requires them).
    if (dippr) {
        const bool has_Tc = dippr->Tc.has_value() && std::isfinite(*dippr->Tc) && (*dippr->Tc > 0.0);
        const bool has_Pc = dippr->Pc.has_value() && std::isfinite(*dippr->Pc) && (*dippr->Pc > 0.0);
        const bool has_Zc = dippr->Zc.has_value() && std::isfinite(*dippr->Zc) && (*dippr->Zc > 0.0);
        const bool has_omega = dippr->omega.has_value() && std::isfinite(*dippr->omega);

        if (has_Tc && has_Pc) {
            const double omega = has_omega ? *dippr->omega : std::numeric_limits<double>::quiet_NaN();
            comp = comp.withCritical(*dippr->Tc, *dippr->Pc, omega);
            if (has_Zc) {
                comp = comp.withCriticalCompressibility(*dippr->Zc);
            }
        }

        if (options.require_critical_properties && !(has_Tc && has_Pc && has_omega)) {
            std::string missing;
            if (!has_Tc) missing += " Tc";
            if (!has_Pc) missing += " Pc";
            if (!has_omega) missing += " omega";
            throw std::runtime_error("Critical properties required but missing for component '" + identifier + "':" + missing);
        }
    } else if (options.require_critical_properties) {
        throw std::runtime_error("Critical properties required but no DIPPR record found for component '" + identifier + "'");
    }

    // Ideal-gas Cp correlation (optional; stored as Core POD coefficients).
    if (dippr && dippr->ideal_gas_cp) {
        const auto& c = *dippr->ideal_gas_cp;
        try {
            comp = comp.withIdealGasCp(convertIdealGasCpToPerMol(c));
        } catch (const std::exception& e) {
            if (options.require_ideal_gas_cp) {
                throw std::runtime_error(std::string("Ideal-gas Cp (ICP_*) required but could not be parsed for '") + identifier + "': " + e.what());
            }
            if (diag) {
                diag->warn(std::string("Ideal-gas Cp not set for '") + identifier + "': " + e.what());
            }
        }
    } else if (options.require_ideal_gas_cp) {
        throw std::runtime_error("Ideal-gas Cp (ICP_*) required but missing for component '" + identifier + "'");
    }

    // VTPR inputs (optional unless caller requires them): c(T) and UNIFAC subgroup decomposition.
    {
        const Data::VtprPureRecord* vp = nullptr;
        std::vector<Data::VtprGroupRecord> vg;

        if (resolved.key.chemid.has_value()) {
            vp = databanks.vtpr.findPureByCHEMID(*resolved.key.chemid);
            vg = databanks.vtpr.groupsByCHEMID(*resolved.key.chemid);
        }
        if (!vp && !resolved.key.cas.empty()) {
            vp = databanks.vtpr.findPureByCAS(resolved.key.cas);
            if (vg.empty()) vg = databanks.vtpr.groupsByCAS(resolved.key.cas);
        }
        if (!vp && !resolved.key.name.empty()) {
            vp = databanks.vtpr.findPureByNameUnique(resolved.key.name, diag);
            if (vg.empty()) vg = databanks.vtpr.groupsByNameUnique(resolved.key.name, diag);
        }

        if (vp) {
            Core::VtprVolumeTranslationPoly c;
            c.c0 = vp->c0.value_or(0.0);
            c.c1 = vp->c1.value_or(0.0);
            c.c2 = vp->c2.value_or(0.0);
            c.c3 = vp->c3.value_or(0.0);
            c.Tmin = vp->Tmin;
            c.Tmax = vp->Tmax;
            // If all coefficients are zero and no explicit Tmin/Tmax, treat as "absent".
            const bool nonzero = (c.c0 != 0.0) || (c.c1 != 0.0) || (c.c2 != 0.0) || (c.c3 != 0.0);
            if (nonzero || c.Tmin.has_value() || c.Tmax.has_value()) {
                comp = comp.withVtprVolumeTranslation(c);
            }
        }

        if (!vg.empty()) {
            std::vector<Core::UnifacSubgroupCount> groups;
            groups.reserve(vg.size());
            for (const auto& row : vg) {
                if (row.subgroup_id <= 0 || row.count <= 0) continue;
                groups.push_back(Core::UnifacSubgroupCount{row.subgroup_id, row.count});
            }
            if (!groups.empty()) {
                comp = comp.withUnifacSubgroups(groups);
            }
        }

        if (options.require_vtpr_params) {
            if (!comp.hasVtprVolumeTranslation()) {
                throw std::runtime_error("VTPR parameters required but missing c(T) for component '" + identifier + "' (vtpr_pure.csv)");
            }
            if (comp.unifacSubgroups().empty()) {
                throw std::runtime_error("VTPR parameters required but missing UNIFAC subgroup mapping for component '" + identifier + "' (vtpr_groups.csv)");
            }
            if (!databanks.vtpr.hasUnifacTables()) {
                throw std::runtime_error("VTPR parameters required but UNIFAC tables are missing (unifac_subgroups.csv / unifac_interactions.csv)");
            }
        }
    }

    return comp;
}

Core::Mixture buildMixtureFromDatabanks(
    const Data::Databanks& databanks,
    const std::vector<std::string>& identifiers,
    const MixtureBuildOptions& options,
    Data::Diagnostics* diag)
{
    if (identifiers.empty()) {
        throw std::invalid_argument("Mixture must have at least one component");
    }

    std::vector<Core::Component> components;
    components.reserve(identifiers.size());
    for (const auto& id : identifiers) {
        components.push_back(buildComponentFromDatabanks(databanks, id, options, diag));
    }

    const int n = static_cast<int>(components.size());
    Core::BinaryParameters kij = Core::BinaryParameters::zeros(n);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const std::string& cas_i = components[i].cas();
            const std::string& cas_j = components[j].cas();
            if (cas_i.empty() || cas_j.empty()) continue;

            const auto k = databanks.bip.kij(options.bip_model, cas_i, cas_j, options.bip_reference_temperature, diag);
            if (k) {
                kij.set(i, j, *k);
            }
        }
    }

    std::shared_ptr<const Data::UnifacTables> unifac_tables;
    if (databanks.vtpr.hasUnifacTables()) {
        unifac_tables = std::make_shared<Data::UnifacTables>(databanks.vtpr.unifac());
    }
    if (options.require_vtpr_params && !unifac_tables) {
        throw std::runtime_error("VTPR parameters required but UNIFAC tables were not loaded");
    }

    return Core::Mixture(components, kij, std::move(unifac_tables));
}

} // namespace Factory
} // namespace DMThermo
