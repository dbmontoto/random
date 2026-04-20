#include "thermo/pcsaft/core/component.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>

#include "thermo/data/databanks.h"
#include "thermo/data/csv.h"

namespace pcsaft {

namespace {

std::optional<std::filesystem::path> findDefaultDataDir() {
    std::filesystem::path p = std::filesystem::current_path();
    for (int i = 0; i < 6; ++i) {
        const auto candidate = p / "data" / "pcsaft.csv";
        if (std::filesystem::exists(candidate)) {
            return p / "data";
        }
        if (!p.has_parent_path()) break;
        p = p.parent_path();
    }
    return std::nullopt;
}

std::optional<AssociationScheme> parseAssociationScheme(const std::string& raw) {
    std::string s = DMThermo::Data::toLower(DMThermo::Data::trim(raw));
    if (s.empty() || s == "none" || s == "0") return AssociationScheme::NONE;
    if (s == "1a" || s == "scheme_1a") return AssociationScheme::SCHEME_1A;
    if (s == "2b" || s == "scheme_2b") return AssociationScheme::SCHEME_2B;
    if (s == "3b" || s == "scheme_3b") return AssociationScheme::SCHEME_3B;
    if (s == "4c" || s == "scheme_4c") return AssociationScheme::SCHEME_4C;
    return std::nullopt;
}

std::optional<DIPPR::EquationType> mapDipprEquation(int eq_form) {
    switch (eq_form) {
        case 100: return DIPPR::EquationType::EQ_100;
        case 101: return DIPPR::EquationType::EQ_101;
        case 105: return DIPPR::EquationType::EQ_105;
        case 106: return DIPPR::EquationType::EQ_106;
        case 107: return DIPPR::EquationType::EQ_107;
        default: return std::nullopt;
    }
}

IdealGasCpParams idealGasCpPerMol(const DMThermo::Data::Correlation& c) {
    const auto eq = mapDipprEquation(c.eq_form);
    if (!eq) {
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

    if (*eq == DIPPR::EquationType::EQ_107) {
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

    return IdealGasCpParams(
        *eq,
        A, B, C, D, E, F, G,
        c.Tmin.value_or(0.0),
        c.Tmax.value_or(1000.0)
    );
}

} // namespace

// ============================================================================
// Component Implementation
// ============================================================================

Component::Component()
    : name_(""), cas_(""), m_(0.0), sigma_(0.0), epsilon_k_(0.0), mw_(0.0), c_(0.0),
      is_polymer_(false), assoc_params_(), cp_ig_params_(), poly_params_() {}

Component::Component(const std::string& name, double m, double sigma, double epsilon_k)
    : name_(name), cas_(""), m_(m), sigma_(sigma), epsilon_k_(epsilon_k), mw_(0.0), c_(0.0),
      is_polymer_(false), assoc_params_(), cp_ig_params_(), poly_params_() {}

Component::Component(const std::string& name, double m, double sigma, double epsilon_k,
                     const AssociationParams& assoc)
    : name_(name), cas_(""), m_(m), sigma_(sigma), epsilon_k_(epsilon_k), mw_(0.0), c_(0.0),
      is_polymer_(false), assoc_params_(assoc), cp_ig_params_(), poly_params_() {}

// Polymer constructor: segment number = m_per_mw * MW
Component::Component(const std::string& name, const PolymerParams& poly,
                     double sigma, double epsilon_k, double mw)
    : name_(name), cas_(""), m_(poly.m_per_mw * mw), sigma_(sigma), epsilon_k_(epsilon_k),
      mw_(mw), c_(0.0), is_polymer_(true), assoc_params_(), cp_ig_params_(), poly_params_(poly) {}

void Component::setPolymerParams(const PolymerParams& params) {
    poly_params_ = params;
    is_polymer_ = params.isValid();
    if (is_polymer_ && mw_ > 0) {
        m_ = poly_params_.m_per_mw * mw_;
    }
}

Component Component::withMW(double new_mw) const {
    if (!is_polymer_) {
        // For non-polymers, just copy and update MW
        Component copy = *this;
        copy.mw_ = new_mw;
        return copy;
    }
    // For polymers, recalculate segment number from m_per_mw
    Component poly(name_, poly_params_, sigma_, epsilon_k_, new_mw);
    poly.cas_ = cas_;
    poly.c_ = c_;
    poly.assoc_params_ = assoc_params_;
    poly.cp_ig_params_ = cp_ig_params_;
    return poly;
}

void Component::print() const {
    std::cout << "Component: " << name_;
    if (is_polymer_) {
        std::cout << " (Polymer)";
    }
    std::cout << "\n";
    if (!cas_.empty()) {
        std::cout << "  CAS: " << cas_ << "\n";
    }
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  m      = " << m_ << " [-]\n";
    std::cout << "  sigma  = " << sigma_ << " [Angstrom]\n";
    std::cout << "  eps/k  = " << epsilon_k_ << " [K]\n";
    if (mw_ > 0.0) {
        std::cout << "  MW     = " << mw_ << " [g/mol]\n";
    }
    if (is_polymer_) {
        std::cout << "  Polymer parameters:\n";
        std::cout << "    m/MW   = " << poly_params_.m_per_mw << " [mol/g]\n";
        if (poly_params_.Mn > 0) {
            std::cout << "    Mn     = " << poly_params_.Mn << " [g/mol]\n";
        }
        if (poly_params_.Mw > 0) {
            std::cout << "    Mw     = " << poly_params_.Mw << " [g/mol]\n";
            std::cout << "    PDI    = " << poly_params_.getPDI() << " [-]\n";
        }
    }
    if (hasAssociation()) {
        std::cout << "  Association:\n";
        std::cout << "    eps_AB/k = " << assoc_params_.epsilon_AB << " [K]\n";
        std::cout << "    kappa_AB = " << assoc_params_.kappa_AB << " [-]\n";
        std::cout << "    sites    = " << assoc_params_.num_sites << "\n";
    }
}

// ============================================================================
// ComponentDatabase Implementation
// ============================================================================

ComponentDatabase& ComponentDatabase::getInstance() {
    static ComponentDatabase instance;
    return instance;
}

void ComponentDatabase::loadDefault() {
    initializeDefaultDatabase();
}

void ComponentDatabase::initializeDefaultDatabase() {
    components_by_name_.clear();
    components_by_cas_.clear();

    const auto data_dir = findDefaultDataDir();
    if (!data_dir) {
        throw std::runtime_error("pcsaft::ComponentDatabase: could not locate data/pcsaft.csv (searched up from current working directory)");
    }

    DMThermo::Data::Diagnostics diag;
    DMThermo::Data::Databanks databanks;
    const bool ok = databanks.loadAllFromDirectory(*data_dir, &diag);
    if (!ok || diag.hasErrors()) {
        std::string msg = "pcsaft::ComponentDatabase: failed to load databanks from " + data_dir->string();
        for (const auto& d : diag.items()) {
            if (d.level == DMThermo::Data::DiagnosticLevel::Error) {
                msg += "\n  - " + d.message;
            }
        }
        throw std::runtime_error(msg);
    }

    for (const auto& r : databanks.pcsaft.records()) {
        const std::string name = DMThermo::Data::normalizeName(r.key.name);
        if (name.empty()) {
            throw std::runtime_error("pcsaft::ComponentDatabase: pcsaft.csv contains a row with empty NAME");
        }
        if (!r.m || !r.sigma || !r.epsilon_k || *r.m <= 0.0 || *r.sigma <= 0.0 || *r.epsilon_k <= 0.0) {
            throw std::runtime_error("pcsaft::ComponentDatabase: missing required PC-SAFT parameters (M/SIGMA/EPSILON_K) for '" + name + "'");
        }

        const bool has_assoc_params = r.kappa_ab && r.epsilon_k_ab && (*r.kappa_ab > 0.0) && (*r.epsilon_k_ab > 0.0);
        const bool has_assoc_meta = r.assoc_scheme.has_value() || r.assoc_num_sites.has_value();

        Component comp;
        if (has_assoc_params) {
            if (!r.assoc_scheme || !r.assoc_num_sites) {
                throw std::runtime_error("pcsaft::ComponentDatabase: association parameters present but ASSOC_SCHEME/ASSOC_NUM_SITES missing for '" + name + "'");
            }
            const auto scheme = parseAssociationScheme(*r.assoc_scheme);
            if (!scheme) {
                throw std::runtime_error("pcsaft::ComponentDatabase: unknown ASSOC_SCHEME '" + *r.assoc_scheme + "' for '" + name + "'");
            }
            if (*r.assoc_num_sites <= 0) {
                throw std::runtime_error("pcsaft::ComponentDatabase: ASSOC_NUM_SITES must be positive for '" + name + "'");
            }
            AssociationParams assoc(*scheme, *r.epsilon_k_ab, *r.kappa_ab, *r.assoc_num_sites);
            comp = Component(name, *r.m, *r.sigma, *r.epsilon_k, assoc);
        } else {
            if (has_assoc_meta) {
                throw std::runtime_error("pcsaft::ComponentDatabase: ASSOC_SCHEME/ASSOC_NUM_SITES provided but KAPPA_AB/EPSILON_K_AB missing for '" + name + "'");
            }
            comp = Component(name, *r.m, *r.sigma, *r.epsilon_k);
        }

        if (!r.key.cas.empty()) {
            comp.setCAS(r.key.cas);
        }
        if (r.MW) {
            comp.setMW(*r.MW);
        }

        const DMThermo::Data::DipprRecord* dippr = nullptr;
        if (!r.key.cas.empty()) {
            dippr = databanks.dippr.findByCAS(r.key.cas);
        }
        if (!dippr) {
            dippr = databanks.dippr.findByNameUnique(r.key.name, &diag);
        }
        if (dippr && dippr->ideal_gas_cp) {
            try {
                comp.setIdealGasCpParams(idealGasCpPerMol(*dippr->ideal_gas_cp));
            } catch (...) {
                // Ideal-gas Cp is optional for legacy PC-SAFT calculations.
            }
        }

        addComponent(comp);
    }

    if (components_by_name_.empty()) {
        throw std::runtime_error("pcsaft::ComponentDatabase: loaded 0 components from CSV databanks");
    }
}

bool ComponentDatabase::loadFromFile(const std::string& filename) {
    const std::filesystem::path p(filename);
    std::filesystem::path dir = p;
    if (!std::filesystem::is_directory(dir)) {
        dir = p.has_parent_path() ? p.parent_path() : std::filesystem::current_path();
    }

    components_by_name_.clear();
    components_by_cas_.clear();

    DMThermo::Data::Diagnostics diag;
    DMThermo::Data::Databanks databanks;
    const bool ok = databanks.loadAllFromDirectory(dir, &diag);
    if (!ok || diag.hasErrors()) {
        return false;
    }

    for (const auto& r : databanks.pcsaft.records()) {
        const std::string name = DMThermo::Data::normalizeName(r.key.name);
        if (name.empty()) continue;
        if (!r.m || !r.sigma || !r.epsilon_k || *r.m <= 0.0 || *r.sigma <= 0.0 || *r.epsilon_k <= 0.0) continue;

        const bool has_assoc_params = r.kappa_ab && r.epsilon_k_ab && (*r.kappa_ab > 0.0) && (*r.epsilon_k_ab > 0.0);
        const bool has_assoc_meta = r.assoc_scheme.has_value() || r.assoc_num_sites.has_value();

        Component comp;
        if (has_assoc_params && r.assoc_scheme && r.assoc_num_sites) {
            const auto scheme = parseAssociationScheme(*r.assoc_scheme);
            if (scheme && *r.assoc_num_sites > 0) {
                AssociationParams assoc(*scheme, *r.epsilon_k_ab, *r.kappa_ab, *r.assoc_num_sites);
                comp = Component(name, *r.m, *r.sigma, *r.epsilon_k, assoc);
            } else {
                comp = Component(name, *r.m, *r.sigma, *r.epsilon_k);
            }
        } else {
            (void)has_assoc_meta;
            comp = Component(name, *r.m, *r.sigma, *r.epsilon_k);
        }

        if (!r.key.cas.empty()) comp.setCAS(r.key.cas);
        if (r.MW) comp.setMW(*r.MW);
        addComponent(comp);
    }

    return !components_by_name_.empty();
}

Component ComponentDatabase::getByName(const std::string& name) const {
    // Convert to lowercase for case-insensitive search
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    auto it = components_by_name_.find(name_lower);
    if (it != components_by_name_.end()) {
        return it->second;
    }

    throw std::runtime_error("Component '" + name + "' not found in database. "
                           "Use getAvailableComponents() to see available components.");
}

Component ComponentDatabase::getByCAS(const std::string& cas) const {
    auto it = components_by_cas_.find(cas);
    if (it != components_by_cas_.end()) {
        return it->second;
    }

    throw std::runtime_error("Component with CAS '" + cas + "' not found in database.");
}

bool ComponentDatabase::hasComponent(const std::string& name) const {
    std::string name_lower = name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return components_by_name_.find(name_lower) != components_by_name_.end();
}

std::vector<std::string> ComponentDatabase::getAvailableComponents() const {
    std::vector<std::string> names;
    names.reserve(components_by_name_.size());
    for (const auto& pair : components_by_name_) {
        names.push_back(pair.second.getName());
    }
    return names;
}

void ComponentDatabase::addComponent(const Component& comp) {
    std::string name_lower = comp.getName();
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    components_by_name_[name_lower] = comp;

    if (!comp.getCAS().empty()) {
        components_by_cas_[comp.getCAS()] = comp;
    }
}

} // namespace pcsaft
