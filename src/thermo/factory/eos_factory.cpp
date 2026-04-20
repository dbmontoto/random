/**
 * @file eos_factory.cpp
 * @brief Implementation of EOS factory
 */

#include "thermo/factory/eos_factory.h"
#include "thermo/eos/PCSAFT.h"
#include "thermo/eos/PengRobinson.h"
#include "thermo/eos/PRVT.h"
#include "thermo/eos/SRK.h"
#include "thermo/eos/VTPR.h"
#include "thermo/legacy/pcsaft/conversions.h"
#include "thermo/factory/mixture_factory.h"
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <limits>

namespace DMThermo {
namespace Factory {

bool EOSFactory::initialized_ = false;

std::unordered_map<std::string, EOSFactory::EOSCreator>& EOSFactory::registry() {
    static std::unordered_map<std::string, EOSCreator> registry;
    return registry;
}

void EOSFactory::initializeBuiltinEOS() {
    if (initialized_) return;

    // Register PC-SAFT - uses conversion to legacy mixture
    registerEOS("PC-SAFT", [](const Core::Mixture& mix) -> EOSPtr {
        return createPCSaft(mix);
    });

    // Placeholder for future EOS implementations
    registerEOS("Peng-Robinson", [](const Core::Mixture& mix) -> EOSPtr {
        return createPengRobinson(mix);
    });

    registerEOS("PR-VT", [](const Core::Mixture& mix) -> EOSPtr {
        return std::make_shared<Cubic::PRVTEOS>(mix);
    });

    registerEOS("VTPR", [](const Core::Mixture& mix) -> EOSPtr {
        return std::make_shared<Cubic::VTPREOS>(mix);
    });

    registerEOS("SRK", [](const Core::Mixture& mix) -> EOSPtr {
        return createSRK(mix);
    });

    initialized_ = true;
}

EOSPtr EOSFactory::create(const std::string& eos_type, const Core::Mixture& mixture) {
    initializeBuiltinEOS();

    auto it = registry().find(eos_type);
    if (it == registry().end()) {
        throw std::invalid_argument("Unknown EOS type: " + eos_type);
    }

    return it->second(mixture);
}

EOSPtr EOSFactory::create(const std::string& eos_type, Core::MixturePtr mixture) {
    if (!mixture) {
        throw std::invalid_argument("Mixture pointer is null");
    }
    return create(eos_type, *mixture);
}

EOSPtr EOSFactory::createPCSaft(const Core::Mixture& mixture) {
    return std::make_shared<PCSaft::PCSaftEOS>(mixture);
}

EOSPtr EOSFactory::createPCSaft(Core::MixturePtr mixture) {
    if (!mixture) {
        throw std::invalid_argument("Mixture pointer is null");
    }
    return createPCSaft(*mixture);
}

EOSPtr EOSFactory::createPengRobinson(const Core::Mixture& mixture) {
    return std::make_shared<Cubic::PengRobinsonEOS>(mixture);
}

EOSPtr EOSFactory::createPengRobinson(Core::MixturePtr mixture) {
    if (!mixture) {
        throw std::invalid_argument("Mixture pointer is null");
    }
    return createPengRobinson(*mixture);
}

EOSPtr EOSFactory::createSRK(const Core::Mixture& mixture) {
    return std::make_shared<Cubic::SRKEOS>(mixture);
}

EOSPtr EOSFactory::createSRK(Core::MixturePtr mixture) {
    if (!mixture) {
        throw std::invalid_argument("Mixture pointer is null");
    }
    return createSRK(*mixture);
}

void EOSFactory::registerEOS(const std::string& name, EOSCreator creator) {
    registry()[name] = creator;
}

bool EOSFactory::isRegistered(const std::string& name) {
    initializeBuiltinEOS();
    return registry().find(name) != registry().end();
}

std::vector<std::string> EOSFactory::registeredTypes() {
    initializeBuiltinEOS();
    std::vector<std::string> types;
    for (const auto& pair : registry()) {
        types.push_back(pair.first);
    }
    std::sort(types.begin(), types.end());
    return types;
}

EOSPtr EOSFactory::createFromDatabanks(
    const std::string& eos_type,
    const Data::Databanks& databanks,
    const std::vector<std::string>& identifiers,
    MixtureBuildOptions options,
    Data::Diagnostics* diag)
{
    std::string type_lc = eos_type;
    std::transform(type_lc.begin(), type_lc.end(), type_lc.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    // Default model-specific options if caller didn't override.
    if (type_lc.find("pc") != std::string::npos && type_lc.find("saft") != std::string::npos) {
        options.bip_model = Data::BipModel::PCSaft;
        options.require_pcsaft_params = true;
        options.require_critical_properties = false;
        options.require_vtpr_params = false;
    } else if (type_lc.find("srk") != std::string::npos) {
        options.bip_model = Data::BipModel::PengRobinson;
        options.require_pcsaft_params = false;
        options.require_critical_properties = true;
        options.require_vtpr_params = false;
    } else if (type_lc.find("pr-vt") != std::string::npos || type_lc.find("prvt") != std::string::npos) {
        options.bip_model = Data::BipModel::PengRobinson; // current databank uses KIJ_PR for cubics
        options.require_pcsaft_params = false;
        options.require_critical_properties = true;
        options.require_vtpr_params = false;
    } else if (type_lc.find("vtpr") != std::string::npos) {
        options.bip_model = Data::BipModel::VTPR;
        options.require_pcsaft_params = false;
        options.require_critical_properties = true;
        options.require_vtpr_params = true;
    } else if (type_lc.find("peng") != std::string::npos || type_lc.find("robinson") != std::string::npos) {
        options.bip_model = Data::BipModel::PengRobinson;
        options.require_pcsaft_params = false;
        options.require_critical_properties = true;
        options.require_vtpr_params = false;
    }

    Core::Mixture mix = buildMixtureFromDatabanks(databanks, identifiers, options, diag);
    return create(eos_type, mix);
}

// Convenience function accepting legacy mixture directly - declaration only
// Implementation is above

} // namespace Factory
} // namespace DMThermo
