/**
 * @file component_resolver.h
 * @brief Resolve components across multiple CSV databanks (CHEMID/CAS/name)
 */

#ifndef THERMO_DATA_COMPONENT_RESOLVER_H
#define THERMO_DATA_COMPONENT_RESOLVER_H

#include "databanks.h"
#include "diagnostics.h"
#include "csv.h"

#include <string>

namespace DMThermo {
namespace Data {

struct ResolvedComponent {
    ChemicalKey key;
    const DipprRecord* dippr = nullptr;
    const PCSaftRecord* pcsaft = nullptr;
};

class ComponentResolver {
public:
    explicit ComponentResolver(const Databanks& databanks)
        : databanks_(databanks) {}

    ResolvedComponent resolve(const std::string& identifier, Diagnostics* diag = nullptr) const {
        ResolvedComponent out;

        // Heuristics:
        // - digits only => CHEMID
        // - contains '-' and digits => CAS
        // - otherwise => name
        const bool looks_chemid = !identifier.empty() && identifier.find_first_not_of("0123456789") == std::string::npos;
        const bool looks_cas =
            !looks_chemid &&
            (identifier.find('-') != std::string::npos) &&
            (identifier.find_first_of("0123456789") != std::string::npos);

        if (looks_chemid) {
            const int chemid = std::stoi(identifier);
            out.key.chemid = chemid;
            out.dippr = databanks_.dippr.findByCHEMID(chemid);
            out.pcsaft = databanks_.pcsaft.findByCHEMID(chemid);
        } else if (looks_cas) {
            out.key.cas = normalizeCAS(identifier);
            out.dippr = databanks_.dippr.findByCAS(out.key.cas);
            out.pcsaft = databanks_.pcsaft.findByCAS(out.key.cas);
        } else {
            out.key.name = identifier;
            out.dippr = databanks_.dippr.findByNameUnique(identifier, diag);
            out.pcsaft = databanks_.pcsaft.findByNameUnique(identifier, diag);
        }

        // Cross-databank "join" pass:
        // If one record was found, try the other databank using the resolved key fields.
        if (out.dippr && !out.pcsaft) {
            if (out.dippr->key.chemid) {
                out.pcsaft = databanks_.pcsaft.findByCHEMID(*out.dippr->key.chemid);
            }
            if (!out.pcsaft && !out.dippr->key.cas.empty()) {
                out.pcsaft = databanks_.pcsaft.findByCAS(out.dippr->key.cas);
            }
            if (!out.pcsaft && !out.dippr->key.name.empty()) {
                out.pcsaft = databanks_.pcsaft.findByNameUnique(out.dippr->key.name, diag);
            }
        }
        if (out.pcsaft && !out.dippr) {
            if (out.pcsaft->key.chemid) {
                out.dippr = databanks_.dippr.findByCHEMID(*out.pcsaft->key.chemid);
            }
            if (!out.dippr && !out.pcsaft->key.cas.empty()) {
                out.dippr = databanks_.dippr.findByCAS(out.pcsaft->key.cas);
            }
            if (!out.dippr && !out.pcsaft->key.name.empty()) {
                out.dippr = databanks_.dippr.findByNameUnique(out.pcsaft->key.name, diag);
            }
        }

        if (out.dippr) {
            if (!out.dippr->key.name.empty()) out.key.name = out.dippr->key.name;
            if (!out.dippr->key.cas.empty()) out.key.cas = out.dippr->key.cas;
            if (out.dippr->key.chemid) out.key.chemid = out.dippr->key.chemid;
        }
        if (out.pcsaft) {
            if (out.key.name.empty() && !out.pcsaft->key.name.empty()) out.key.name = out.pcsaft->key.name;
            if (out.key.cas.empty() && !out.pcsaft->key.cas.empty()) out.key.cas = out.pcsaft->key.cas;
            if (!out.key.chemid && out.pcsaft->key.chemid) out.key.chemid = out.pcsaft->key.chemid;
        }

        if (!out.dippr && !out.pcsaft) {
            if (diag) {
                diag->error(
                    "Component '" + identifier + "' not found in DIPPR or PC-SAFT databanks",
                    Config::ErrorCategory::MissingData,
                    "COMPONENT_NOT_FOUND");
            }
        }

        // Cross-check CAS consistency when both exist.
        if (diag && out.dippr && out.pcsaft) {
            const auto a = normalizeCAS(out.dippr->key.cas);
            const auto b = normalizeCAS(out.pcsaft->key.cas);
            if (!a.empty() && !b.empty() && a != b) {
                diag->warn(
                    "Resolved component '" + identifier + "' has mismatched CAS between databanks: " + a + " vs " + b,
                    Config::ErrorCategory::MissingData,
                    "COMPONENT_CAS_MISMATCH");
            }
        }

        return out;
    }

private:
    const Databanks& databanks_;
};

} // namespace Data
} // namespace DMThermo

#endif // THERMO_DATA_COMPONENT_RESOLVER_H
