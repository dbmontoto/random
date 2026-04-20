#include "thermo/data/databanks.h"
#include "thermo/data/csv.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace DMThermo {
namespace Data {

namespace {

bool requireHeaders(const CsvTable& t, const std::vector<std::string>& required_lc, const std::string& label, Diagnostics* diag) {
    bool ok = true;
    for (const auto& h : required_lc) {
        if (t.header_index.find(h) == t.header_index.end()) {
            ok = false;
            if (diag) {
                diag->error(
                    label + " missing required column: " + h,
                    Config::ErrorCategory::MissingData,
                    "CSV_MISSING_REQUIRED_COLUMN");
            }
        }
    }
    return ok;
}

void warnHeaderConflicts(const CsvTable& t, const std::string& label, Diagnostics* diag) {
    if (!diag) return;
    if (t.header_conflicts.empty()) return;

    for (const auto& c : t.header_conflicts) {
        std::ostringstream oss;
        oss << label << " header key conflict for '" << c.key << "' at columns [";
        for (std::size_t k = 0; k < c.columns.size(); ++k) {
            if (k) oss << ", ";
            oss << c.columns[k];
        }
        oss << "] (using first). Raw headers: ";
        for (std::size_t k = 0; k < c.columns.size(); ++k) {
            const std::size_t col = c.columns[k];
            if (k) oss << " | ";
            if (col < t.headers.size()) {
                oss << "'" << t.headers[col] << "'";
            } else {
                oss << "<out-of-range>";
            }
        }
        diag->warn(
            oss.str(),
            Config::ErrorCategory::InvalidInput,
            "CSV_HEADER_CONFLICT");
    }
}

std::optional<std::string> getField(const CsvTable& t, const std::vector<std::string>& row, const std::string& header_lc) {
    auto it = t.header_index.find(header_lc);
    if (it == t.header_index.end()) return std::nullopt;
    const std::size_t idx = it->second;
    if (idx >= row.size()) return std::nullopt;
    return row[idx];
}

std::optional<double> getDouble(const CsvTable& t, const std::vector<std::string>& row, const std::string& header_lc) {
    auto s = getField(t, row, header_lc);
    if (!s) return std::nullopt;
    return parseDouble(*s);
}

std::optional<int> getInt(const CsvTable& t, const std::vector<std::string>& row, const std::string& header_lc) {
    auto s = getField(t, row, header_lc);
    if (!s) return std::nullopt;
    return parseInt(*s);
}

std::string getString(const CsvTable& t, const std::vector<std::string>& row, const std::string& header_lc) {
    auto s = getField(t, row, header_lc);
    return s ? *s : std::string{};
}

std::optional<Correlation> parseCorrelation(
    const CsvTable& t,
    const std::vector<std::string>& row,
    const std::string& prefix_lc)
{
    // Require an equation form number; if missing/invalid, treat as absent.
    auto eq = getInt(t, row, prefix_lc + "_dippr_eqn");
    if (!eq || *eq <= 0) return std::nullopt;

    Correlation c;
    c.eq_form = *eq;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // A-E are treated as required once an equation form is present. Do not default to 0.0:
    // missing coefficients should flow through as invalid and be caught at evaluation time.
    c.A = getDouble(t, row, prefix_lc + "_a").value_or(nan);
    c.B = getDouble(t, row, prefix_lc + "_b").value_or(nan);
    c.C = getDouble(t, row, prefix_lc + "_c").value_or(nan);
    c.D = getDouble(t, row, prefix_lc + "_d").value_or(nan);
    // Some equation forms don't use (or allow omitting) E; treat it as optional there.
    // Note: EQ_123 does use E, so do not default it.
    if (*eq == 102 || *eq == 105 || *eq == 114) {
        c.E = getDouble(t, row, prefix_lc + "_e").value_or(0.0);
    } else {
        c.E = getDouble(t, row, prefix_lc + "_e").value_or(nan);
    }

    // F/G exist in the schema for some groups (e.g., ICP). When absent, default to 0.0.
    c.F = getDouble(t, row, prefix_lc + "_f").value_or(0.0);
    c.G = getDouble(t, row, prefix_lc + "_g").value_or(0.0);

    c.Tmin = getDouble(t, row, prefix_lc + "_min_t");
    c.Tmax = getDouble(t, row, prefix_lc + "_max_t");

    c.data_type = getString(t, row, prefix_lc + "_data_type");
    c.error = getString(t, row, prefix_lc + "_error");
    return c;
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

bool DipprCsvDatabank::load(const std::filesystem::path& csv_path, Diagnostics* diag) {
    records_.clear();
    by_chemid_.clear();
    by_cas_.clear();
    by_name_.clear();

    CsvTable t = readCsvFile(csv_path);
    warnHeaderConflicts(t, "dippr.csv", diag);

    if (!requireHeaders(t, {"name"}, "dippr.csv", diag)) {
        return false;
    }

    for (const auto& row : t.rows) {
        DipprRecord r;

        r.key.chemid = getInt(t, row, "chemid");
        r.key.name = getString(t, row, "name");
        r.key.cas = normalizeCAS(getString(t, row, "cas"));

        r.MW = getDouble(t, row, "mw");
        r.Tc = getDouble(t, row, "tc");
        r.Pc = getDouble(t, row, "pc");
        r.Vc = getDouble(t, row, "vc");
        r.Zc = getDouble(t, row, "zc");
        r.omega = getDouble(t, row, "omega");

        // Standard/formation properties (optional)
        r.HFOR = getDouble(t, row, "hfor");
        r.GFOR = getDouble(t, row, "gfor");
        r.ENT = getDouble(t, row, "ent");
        r.HSTD = getDouble(t, row, "hstd");
        r.SSTD = getDouble(t, row, "sstd");

        r.liquid_density = parseCorrelation(t, row, "ldn");
        r.vapor_pressure = parseCorrelation(t, row, "vp");
        r.liquid_viscosity = parseCorrelation(t, row, "lvs");
        r.heat_of_vaporization = parseCorrelation(t, row, "hvp");
        r.liquid_thermal_cond = parseCorrelation(t, row, "ltc");
        r.liquid_heat_capacity = parseCorrelation(t, row, "lcp");
        r.vapor_viscosity = parseCorrelation(t, row, "vvs");
        r.ideal_gas_cp = parseCorrelation(t, row, "icp");
        r.surface_tension = parseCorrelation(t, row, "st");

        const std::size_t idx = records_.size();
        records_.push_back(std::move(r));

        if (records_[idx].key.chemid) {
            const int chemid = *records_[idx].key.chemid;
            auto [it, inserted] = by_chemid_.emplace(chemid, idx);
            if (!inserted && diag) {
                diag->warn(
                    "dippr.csv has duplicate CHEMID " + std::to_string(chemid) +
                    " (keeping first, ignoring later record)"
                );
            }
        }
        if (!records_[idx].key.cas.empty()) {
            auto [it, inserted] = by_cas_.emplace(records_[idx].key.cas, idx);
            if (!inserted && diag) {
                diag->warn(
                    "dippr.csv has duplicate CAS " + records_[idx].key.cas +
                    " (keeping first, ignoring later record)"
                );
            }
        }
        if (!records_[idx].key.name.empty()) {
            const std::string key = normalizeName(records_[idx].key.name);
            by_name_[key].push_back(idx);
        }
    }

    // Detect ambiguous NAME keys (normalized) and warn: prefer CAS/CHEMID for disambiguation.
    if (diag) {
        for (const auto& pair : by_name_) {
            if (pair.second.size() <= 1) continue;
            diag->warn(
                "dippr.csv has ambiguous NAME '" + pair.first + "' mapping to " + std::to_string(pair.second.size()) +
                " records (" + summarizeAmbiguity(records_, pair.second) + "); prefer CAS/CHEMID identifiers"
            );
        }
    }

    if (diag && records_.empty()) {
        diag->warn("DIPPR databank loaded 0 records from " + csv_path.string());
    }
    return true;
}

const DipprRecord* DipprCsvDatabank::findByCHEMID(int chemid) const {
    auto it = by_chemid_.find(chemid);
    if (it == by_chemid_.end()) return nullptr;
    return &records_[it->second];
}

const DipprRecord* DipprCsvDatabank::findByCAS(const std::string& cas) const {
    auto it = by_cas_.find(normalizeCAS(cas));
    if (it == by_cas_.end()) return nullptr;
    return &records_[it->second];
}

const DipprRecord* DipprCsvDatabank::findByName(const std::string& name) const {
    auto it = by_name_.find(normalizeName(name));
    if (it == by_name_.end() || it->second.empty()) return nullptr;
    return &records_[it->second.front()];
}

const DipprRecord* DipprCsvDatabank::findByNameUnique(const std::string& name, Diagnostics* diag) const {
    const std::string key = normalizeName(name);
    auto it = by_name_.find(key);
    if (it == by_name_.end() || it->second.empty()) return nullptr;
    if (it->second.size() == 1) return &records_[it->second.front()];

    if (diag) {
        diag->error(
            "Ambiguous DIPPR NAME '" + name + "' (" + key + ") matched " + std::to_string(it->second.size()) +
            " records (" + summarizeAmbiguity(records_, it->second) + "); use CAS or CHEMID"
        );
    }
    return nullptr;
}

bool PCSaftCsvDatabank::load(const std::filesystem::path& csv_path, Diagnostics* diag) {
    records_.clear();
    by_chemid_.clear();
    by_cas_.clear();
    by_name_.clear();

    CsvTable t = readCsvFile(csv_path);
    warnHeaderConflicts(t, "pcsaft.csv", diag);

    if (!requireHeaders(t, {"name", "m", "sigma", "epsilon_k"}, "pcsaft.csv", diag)) {
        return false;
    }

    for (const auto& row : t.rows) {
        PCSaftRecord r;
        r.key.chemid = getInt(t, row, "chemid");
        r.key.name = getString(t, row, "name");
        r.key.cas = normalizeCAS(getString(t, row, "cas"));

        r.MW = getDouble(t, row, "mw");
        r.m = getDouble(t, row, "m");
        r.sigma = getDouble(t, row, "sigma");
        r.epsilon_k = getDouble(t, row, "epsilon_k");
        r.kappa_ab = getDouble(t, row, "kappa_ab");
        r.epsilon_k_ab = getDouble(t, row, "epsilon_k_ab");
        {
            const std::string scheme = getString(t, row, "assoc_scheme");
            if (!scheme.empty()) {
                r.assoc_scheme = scheme;
            }
        }
        r.assoc_num_sites = getInt(t, row, "assoc_num_sites");
        r.Tmin = getDouble(t, row, "t_min");
        r.Tmax = getDouble(t, row, "t_max");

        const std::size_t idx = records_.size();
        records_.push_back(std::move(r));

        if (records_[idx].key.chemid) {
            const int chemid = *records_[idx].key.chemid;
            auto [it, inserted] = by_chemid_.emplace(chemid, idx);
            if (!inserted && diag) {
                diag->warn(
                    "pcsaft.csv has duplicate CHEMID " + std::to_string(chemid) +
                    " (keeping first, ignoring later record)"
                );
            }
        }
        if (!records_[idx].key.cas.empty()) {
            auto [it, inserted] = by_cas_.emplace(records_[idx].key.cas, idx);
            if (!inserted && diag) {
                diag->warn(
                    "pcsaft.csv has duplicate CAS " + records_[idx].key.cas +
                    " (keeping first, ignoring later record)"
                );
            }
        }
        if (!records_[idx].key.name.empty()) {
            const std::string key = normalizeName(records_[idx].key.name);
            by_name_[key].push_back(idx);
        }
    }

    // Detect ambiguous NAME keys (normalized) and warn: prefer CAS/CHEMID for disambiguation.
    if (diag) {
        for (const auto& pair : by_name_) {
            if (pair.second.size() <= 1) continue;
            diag->warn(
                "pcsaft.csv has ambiguous NAME '" + pair.first + "' mapping to " + std::to_string(pair.second.size()) +
                " records (" + summarizeAmbiguity(records_, pair.second) + "); prefer CAS/CHEMID identifiers"
            );
        }
    }

    if (diag && records_.empty()) {
        diag->warn("PC-SAFT databank loaded 0 records from " + csv_path.string());
    }
    return true;
}

const PCSaftRecord* PCSaftCsvDatabank::findByCHEMID(int chemid) const {
    auto it = by_chemid_.find(chemid);
    if (it == by_chemid_.end()) return nullptr;
    return &records_[it->second];
}

const PCSaftRecord* PCSaftCsvDatabank::findByCAS(const std::string& cas) const {
    auto it = by_cas_.find(normalizeCAS(cas));
    if (it == by_cas_.end()) return nullptr;
    return &records_[it->second];
}

const PCSaftRecord* PCSaftCsvDatabank::findByName(const std::string& name) const {
    auto it = by_name_.find(normalizeName(name));
    if (it == by_name_.end() || it->second.empty()) return nullptr;
    return &records_[it->second.front()];
}

const PCSaftRecord* PCSaftCsvDatabank::findByNameUnique(const std::string& name, Diagnostics* diag) const {
    const std::string key = normalizeName(name);
    auto it = by_name_.find(key);
    if (it == by_name_.end() || it->second.empty()) return nullptr;
    if (it->second.size() == 1) return &records_[it->second.front()];

    if (diag) {
        diag->error(
            "Ambiguous PC-SAFT NAME '" + name + "' (" + key + ") matched " + std::to_string(it->second.size()) +
            " records (" + summarizeAmbiguity(records_, it->second) + "); use CAS or CHEMID"
        );
    }
    return nullptr;
}

std::size_t BipCsvDatabank::PairKeyHash::operator()(const PairKey& k) const {
    return std::hash<std::string>{}(k.a) ^ (std::hash<std::string>{}(k.b) << 1);
}

bool BipCsvDatabank::PairKeyEq::operator()(const PairKey& x, const PairKey& y) const {
    return x.a == y.a && x.b == y.b;
}

BipCsvDatabank::PairKey BipCsvDatabank::canonicalPair(const std::string& cas_i, const std::string& cas_j) {
    PairKey k{normalizeCAS(cas_i), normalizeCAS(cas_j)};
    if (k.a <= k.b) return k;
    return PairKey{k.b, k.a};
}

bool BipCsvDatabank::load(const std::filesystem::path& csv_path, Diagnostics* diag) {
    records_.clear();
    by_pair_.clear();

    CsvTable t = readCsvFile(csv_path);
    warnHeaderConflicts(t, "bip.csv", diag);

    if (!requireHeaders(t, {"i_cas", "j_cas"}, "bip.csv", diag)) {
        return false;
    }
    const bool has_kij_pr = (t.header_index.find("kij_pr") != t.header_index.end());
    const bool has_kij_vtpr = (t.header_index.find("kij_vtpr") != t.header_index.end());
    const bool has_kij_pcsaft = (t.header_index.find("kij_pcsaft") != t.header_index.end());
    if (!has_kij_pr && !has_kij_vtpr && !has_kij_pcsaft) {
        if (diag) {
            diag->error(
                "bip.csv missing required column: kij_pr, kij_vtpr, or kij_pcsaft",
                Config::ErrorCategory::MissingData,
                "CSV_MISSING_REQUIRED_COLUMN");
        }
        return false;
    }

    for (const auto& row : t.rows) {
        BipRecord r;

        r.i.name = getString(t, row, "i_name");
        r.i.cas = normalizeCAS(getString(t, row, "i_cas"));
        r.i.chemid = std::nullopt;

        r.j.name = getString(t, row, "j_name");
        r.j.cas = normalizeCAS(getString(t, row, "j_cas"));
        r.j.chemid = std::nullopt;

        r.kij_pr = getDouble(t, row, "kij_pr");
        r.kij_vtpr = getDouble(t, row, "kij_vtpr");
        r.kij_pcsaft = getDouble(t, row, "kij_pcsaft");
        r.Tmin = getDouble(t, row, "t_min");
        r.Tmax = getDouble(t, row, "t_max");

        const std::size_t idx = records_.size();
        records_.push_back(std::move(r));

        if (!records_[idx].i.cas.empty() && !records_[idx].j.cas.empty()) {
            PairKey k = canonicalPair(records_[idx].i.cas, records_[idx].j.cas);
            by_pair_[k].push_back(idx);
        }
    }

    if (diag && records_.empty()) {
        diag->warn("BIP databank loaded 0 records from " + csv_path.string());
    }
    return true;
}

std::optional<double> BipCsvDatabank::kij(
    BipModel model,
    const std::string& cas_i,
    const std::string& cas_j,
    double T,
    Diagnostics* diag) const
{
    PairKey k = canonicalPair(cas_i, cas_j);
    auto it = by_pair_.find(k);
    if (it == by_pair_.end()) return std::nullopt;

    const auto& idxs = it->second;
    if (idxs.empty()) return std::nullopt;

    auto valueFor = [&](const BipRecord& r) -> std::optional<double> {
        switch (model) {
            case BipModel::PengRobinson: return r.kij_pr;
            case BipModel::VTPR: return r.kij_vtpr;
            case BipModel::PCSaft: return r.kij_pcsaft;
        }
        return std::nullopt;
    };

    // Select best record for temperature: in-range and narrowest range.
    std::optional<std::size_t> best;
    double best_span = 0.0;
    for (std::size_t idx : idxs) {
        const auto& r = records_[idx];
        const auto v = valueFor(r);
        if (!v) continue;

        if (r.Tmin && r.Tmax) {
            if (T < *r.Tmin || T > *r.Tmax) continue;
            const double span = *r.Tmax - *r.Tmin;
            if (!best || span < best_span) {
                best = idx;
                best_span = span;
            }
        } else {
            // No validity range: accept as fallback if nothing better exists.
            if (!best) {
                best = idx;
                best_span = std::numeric_limits<double>::infinity();
            }
        }
    }

    if (!best) {
        // No in-range record found; return the first available and warn.
        for (std::size_t idx : idxs) {
            const auto v = valueFor(records_[idx]);
            if (v) {
                if (diag) {
                    diag->warn("No in-range BIP record found for pair " + k.a + "/" + k.b + " at T=" + std::to_string(T) + " K; using first available");
                }
                return v;
            }
        }
        return std::nullopt;
    }

    return valueFor(records_[*best]);
}

bool Databanks::loadAllFromDirectory(const std::filesystem::path& dir, Diagnostics* diag) {
    bool ok = true;
    const auto dippr_path = dir / "dippr.csv";
    const auto pcsaft_path = dir / "pcsaft.csv";
    const auto bip_path = dir / "bip.csv";

    try {
        ok = dippr.load(dippr_path, diag) && ok;
    } catch (const std::exception& e) {
        ok = false;
        if (diag) {
            diag->error(
                std::string("Failed to load dippr.csv: ") + e.what(),
                Config::ErrorCategory::MissingData,
                "CSV_LOAD_FAILED");
        }
    }
    try {
        ok = pcsaft.load(pcsaft_path, diag) && ok;
    } catch (const std::exception& e) {
        ok = false;
        if (diag) {
            diag->error(
                std::string("Failed to load pcsaft.csv: ") + e.what(),
                Config::ErrorCategory::MissingData,
                "CSV_LOAD_FAILED");
        }
    }
    try {
        ok = bip.load(bip_path, diag) && ok;
    } catch (const std::exception& e) {
        ok = false;
        if (diag) {
            diag->error(
                std::string("Failed to load bip.csv: ") + e.what(),
                Config::ErrorCategory::MissingData,
                "CSV_LOAD_FAILED");
        }
    }
    try {
        ok = vtpr.loadFromDirectory(dir, diag) && ok;
    } catch (const std::exception& e) {
        ok = false;
        if (diag) {
            diag->error(
                std::string("Failed to load vtpr/unifac databanks: ") + e.what(),
                Config::ErrorCategory::MissingData,
                "CSV_LOAD_FAILED");
        }
    }

    return ok;
}

} // namespace Data
} // namespace DMThermo
