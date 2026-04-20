#include "thermo/data/vtpr.h"

#include <algorithm>
#include <sstream>

namespace DMThermo {
namespace Data {

namespace {

bool requireHeaders(const CsvTable& t, const std::vector<std::string>& required_lc, const std::string& label, Diagnostics* diag) {
    bool ok = true;
    for (const auto& h : required_lc) {
        if (t.header_index.find(h) == t.header_index.end()) {
            ok = false;
            if (diag) diag->error(label + " missing required column: " + h);
        }
    }
    return ok;
}

std::optional<std::string> getField(const CsvTable& t, const std::vector<std::string>& row, const std::string& header_lc) {
    auto it = t.header_index.find(header_lc);
    if (it == t.header_index.end()) return std::nullopt;
    const std::size_t idx = it->second;
    if (idx >= row.size()) return std::nullopt;
    return row[idx];
}

std::string getString(const CsvTable& t, const std::vector<std::string>& row, const std::string& header_lc) {
    auto s = getField(t, row, header_lc);
    return s ? *s : std::string{};
}

std::optional<int> getInt(const CsvTable& t, const std::vector<std::string>& row, const std::string& header_lc) {
    auto s = getField(t, row, header_lc);
    if (!s) return std::nullopt;
    return parseInt(*s);
}

std::optional<double> getDouble(const CsvTable& t, const std::vector<std::string>& row, const std::string& header_lc) {
    auto s = getField(t, row, header_lc);
    if (!s) return std::nullopt;
    return parseDouble(*s);
}

template <typename RecordT>
std::string summarizeAmbiguity(const std::vector<RecordT>& records, const std::vector<std::size_t>& idxs) {
    std::ostringstream oss;
    for (std::size_t k = 0; k < idxs.size(); ++k) {
        if (k) oss << ", ";
        const auto& r = records[idxs[k]];
        if (!r.key.cas.empty()) {
            oss << r.key.cas;
        } else if (r.key.chemid.has_value()) {
            oss << "CHEMID " << r.key.chemid.value();
        } else {
            oss << "<unkeyed>";
        }
    }
    return oss.str();
}

} // namespace

long long VtprCsvDatabank::keyPair(int i, int j) {
    return (static_cast<long long>(i) << 32) | static_cast<unsigned long long>(j);
}

bool VtprCsvDatabank::loadFromDirectory(const std::filesystem::path& dir, Diagnostics* diag) {
    bool ok = true;
    // All vtpr-related files are optional at the Databanks level; EOS/mixture build will enforce requirements.
    const auto pure = dir / FILE_VTPR_PURE;
    const auto groups = dir / FILE_VTPR_GROUPS;
    const auto subgroups = dir / FILE_UNIFAC_SUBGROUPS;
    const auto interactions = dir / FILE_UNIFAC_INTERACTIONS;

    if (std::filesystem::exists(pure)) ok = loadPure(pure, diag) && ok;
    if (std::filesystem::exists(groups)) ok = loadGroups(groups, diag) && ok;
    if (std::filesystem::exists(subgroups)) ok = loadUnifacSubgroups(subgroups, diag) && ok;
    if (std::filesystem::exists(interactions)) ok = loadUnifacInteractions(interactions, diag) && ok;

    return ok;
}

bool VtprCsvDatabank::loadPure(const std::filesystem::path& path, Diagnostics* diag) {
    pure_.clear();
    pure_by_chemid_.clear();
    pure_by_cas_.clear();
    pure_by_name_.clear();

    CsvTable t = readCsvFile(path);
    if (!requireHeaders(t, {"name"}, "vtpr_pure.csv", diag)) {
        return false;
    }

    for (const auto& row : t.rows) {
        VtprPureRecord r;
        r.key.chemid = getInt(t, row, "chemid");
        r.key.name = getString(t, row, "name");
        r.key.cas = normalizeCAS(getString(t, row, "cas"));

        r.c0 = getDouble(t, row, "c0");
        r.c1 = getDouble(t, row, "c1");
        r.c2 = getDouble(t, row, "c2");
        r.c3 = getDouble(t, row, "c3");
        r.Tmin = getDouble(t, row, "t_min");
        r.Tmax = getDouble(t, row, "t_max");

        const std::size_t idx = pure_.size();
        pure_.push_back(std::move(r));

        if (pure_[idx].key.chemid) {
            const int chemid = *pure_[idx].key.chemid;
            auto [it, inserted] = pure_by_chemid_.emplace(chemid, idx);
            if (!inserted && diag) {
                diag->warn("vtpr_pure.csv has duplicate CHEMID " + std::to_string(chemid) + " (keeping first)");
            }
        }
        if (!pure_[idx].key.cas.empty()) {
            auto [it, inserted] = pure_by_cas_.emplace(pure_[idx].key.cas, idx);
            if (!inserted && diag) {
                diag->warn("vtpr_pure.csv has duplicate CAS " + pure_[idx].key.cas + " (keeping first)");
            }
        }
        if (!pure_[idx].key.name.empty()) {
            pure_by_name_[normalizeName(pure_[idx].key.name)].push_back(idx);
        }
    }

    if (diag) {
        for (const auto& pair : pure_by_name_) {
            if (pair.second.size() <= 1) continue;
            diag->warn(
                "vtpr_pure.csv has ambiguous NAME '" + pair.first + "' mapping to " + std::to_string(pair.second.size()) +
                " records (" + summarizeAmbiguity(pure_, pair.second) + "); prefer CAS/CHEMID identifiers"
            );
        }
    }

    return true;
}

bool VtprCsvDatabank::loadGroups(const std::filesystem::path& path, Diagnostics* diag) {
    groups_.clear();
    groups_by_chemid_.clear();
    groups_by_cas_.clear();
    groups_by_name_.clear();

    CsvTable t = readCsvFile(path);
    if (!requireHeaders(t, {"name", "subgroup_id", "count"}, "vtpr_groups.csv", diag)) {
        return false;
    }

    for (const auto& row : t.rows) {
        VtprGroupRecord r;
        r.key.chemid = getInt(t, row, "chemid");
        r.key.name = getString(t, row, "name");
        r.key.cas = normalizeCAS(getString(t, row, "cas"));
        r.subgroup_id = getInt(t, row, "subgroup_id").value_or(0);
        r.count = getInt(t, row, "count").value_or(0);

        const std::size_t idx = groups_.size();
        groups_.push_back(std::move(r));

        if (groups_[idx].key.chemid) {
            groups_by_chemid_[*groups_[idx].key.chemid].push_back(idx);
        }
        if (!groups_[idx].key.cas.empty()) {
            groups_by_cas_[groups_[idx].key.cas].push_back(idx);
        }
        if (!groups_[idx].key.name.empty()) {
            groups_by_name_[normalizeName(groups_[idx].key.name)].push_back(idx);
        }
    }

    return true;
}

bool VtprCsvDatabank::loadUnifacSubgroups(const std::filesystem::path& path, Diagnostics* diag) {
    (void)diag;
    unifac_.subgroups.clear();

    CsvTable t = readCsvFile(path);
    if (!requireHeaders(t, {"subgroup_id", "main_group_id", "r", "q"}, "unifac_subgroups.csv", diag)) {
        return false;
    }

    for (const auto& row : t.rows) {
        UnifacSubgroupRecord r;
        r.subgroup_id = getInt(t, row, "subgroup_id").value_or(0);
        r.main_group_id = getInt(t, row, "main_group_id").value_or(0);
        r.R = getDouble(t, row, "r").value_or(0.0);
        r.Q = getDouble(t, row, "q").value_or(0.0);
        r.name = getString(t, row, "name");
        if (r.subgroup_id > 0) {
            unifac_.subgroups[r.subgroup_id] = std::move(r);
        }
    }

    return true;
}

bool VtprCsvDatabank::loadUnifacInteractions(const std::filesystem::path& path, Diagnostics* diag) {
    unifac_.interactions.clear();

    CsvTable t = readCsvFile(path);
    if (!requireHeaders(t, {"main_i", "main_j", "a"}, "unifac_interactions.csv", diag)) {
        return false;
    }

    for (const auto& row : t.rows) {
        UnifacInteractionRecord r;
        r.main_i = getInt(t, row, "main_i").value_or(0);
        r.main_j = getInt(t, row, "main_j").value_or(0);
        r.a = getDouble(t, row, "a").value_or(0.0);
        r.b = getDouble(t, row, "b").value_or(0.0);
        r.c = getDouble(t, row, "c").value_or(0.0);
        if (r.main_i > 0 && r.main_j > 0) {
            unifac_.interactions[keyPair(r.main_i, r.main_j)] = std::move(r);
        }
    }

    return true;
}

const VtprPureRecord* VtprCsvDatabank::findPureByCHEMID(int chemid) const {
    auto it = pure_by_chemid_.find(chemid);
    if (it == pure_by_chemid_.end()) return nullptr;
    return &pure_[it->second];
}

const VtprPureRecord* VtprCsvDatabank::findPureByCAS(const std::string& cas) const {
    auto it = pure_by_cas_.find(normalizeCAS(cas));
    if (it == pure_by_cas_.end()) return nullptr;
    return &pure_[it->second];
}

const VtprPureRecord* VtprCsvDatabank::findPureByNameUnique(const std::string& name, Diagnostics* diag) const {
    const std::string key = normalizeName(name);
    auto it = pure_by_name_.find(key);
    if (it == pure_by_name_.end() || it->second.empty()) return nullptr;
    if (it->second.size() == 1) return &pure_[it->second.front()];
    if (diag) {
        diag->error(
            "Ambiguous VTPR NAME '" + name + "' (" + key + ") matched " + std::to_string(it->second.size()) +
            " records (" + summarizeAmbiguity(pure_, it->second) + "); use CAS or CHEMID"
        );
    }
    return nullptr;
}

std::vector<VtprGroupRecord> VtprCsvDatabank::groupsByCHEMID(int chemid) const {
    auto it = groups_by_chemid_.find(chemid);
    if (it == groups_by_chemid_.end() || it->second.empty()) return {};
    std::vector<VtprGroupRecord> out;
    out.reserve(it->second.size());
    for (std::size_t idx : it->second) out.push_back(groups_[idx]);
    return out;
}

std::vector<VtprGroupRecord> VtprCsvDatabank::groupsByCAS(const std::string& cas) const {
    auto it = groups_by_cas_.find(normalizeCAS(cas));
    if (it == groups_by_cas_.end() || it->second.empty()) return {};
    std::vector<VtprGroupRecord> out;
    out.reserve(it->second.size());
    for (std::size_t idx : it->second) out.push_back(groups_[idx]);
    return out;
}

std::vector<VtprGroupRecord> VtprCsvDatabank::groupsByNameUnique(const std::string& name, Diagnostics* diag) const {
    const std::string key = normalizeName(name);
    auto it = groups_by_name_.find(key);
    if (it == groups_by_name_.end() || it->second.empty()) return {};
    if (it->second.size() > 1 && diag) {
        diag->warn("VTPR groups NAME '" + name + "' (" + key + ") matched multiple rows; prefer CAS/CHEMID");
    }
    std::vector<VtprGroupRecord> out;
    out.reserve(it->second.size());
    for (std::size_t idx : it->second) out.push_back(groups_[idx]);
    return out;
}

} // namespace Data
} // namespace DMThermo
