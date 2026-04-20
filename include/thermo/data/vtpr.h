/**
 * @file vtpr.h
 * @brief CSV-backed data structures for VTPR (UNIFAC + volume translation)
 */

#ifndef THERMO_DATA_VTPR_H
#define THERMO_DATA_VTPR_H

#include "thermo/data/csv.h"
#include "thermo/data/diagnostics.h"
#include "thermo/data/records.h"

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace DMThermo {
namespace Data {

struct UnifacSubgroupRecord {
    int subgroup_id = 0;
    int main_group_id = 0;
    double R = 0.0;
    double Q = 0.0;
    std::string name;
};

struct UnifacInteractionRecord {
    int main_i = 0;
    int main_j = 0;
    double a = 0.0;  // K
    double b = 0.0;  // K^2 (optional)
    double c = 0.0;  // K^3 (optional)
};

struct VtprGroupRecord {
    ChemicalKey key;
    int subgroup_id = 0;
    int count = 0;
};

struct VtprPureRecord {
    ChemicalKey key;

    // Volume-translation polynomial coefficients (m^3/mol):
    //   c(T) = c0 + c1*T + c2*T^2 + c3*T^3
    std::optional<double> c0;
    std::optional<double> c1;
    std::optional<double> c2;
    std::optional<double> c3;
    std::optional<double> Tmin;
    std::optional<double> Tmax;
};

struct UnifacTables {
    std::unordered_map<int, UnifacSubgroupRecord> subgroups;
    std::unordered_map<long long, UnifacInteractionRecord> interactions; // key(main_i, main_j) (directional)
};

class VtprCsvDatabank {
public:
    // File names (relative to databanks dir)
    static constexpr const char* FILE_VTPR_PURE = "vtpr_pure.csv";
    static constexpr const char* FILE_VTPR_GROUPS = "vtpr_groups.csv";
    static constexpr const char* FILE_UNIFAC_SUBGROUPS = "unifac_subgroups.csv";
    static constexpr const char* FILE_UNIFAC_INTERACTIONS = "unifac_interactions.csv";

    bool loadFromDirectory(const std::filesystem::path& dir, Diagnostics* diag = nullptr);

    // Per-component VTPR inputs
    const VtprPureRecord* findPureByCHEMID(int chemid) const;
    const VtprPureRecord* findPureByCAS(const std::string& cas) const;
    const VtprPureRecord* findPureByNameUnique(const std::string& name, Diagnostics* diag) const;

    std::vector<VtprGroupRecord> groupsByCHEMID(int chemid) const;
    std::vector<VtprGroupRecord> groupsByCAS(const std::string& cas) const;
    std::vector<VtprGroupRecord> groupsByNameUnique(const std::string& name, Diagnostics* diag) const;

    // Global UNIFAC tables (subgroups + interactions)
    const UnifacTables& unifac() const { return unifac_; }
    bool hasUnifacTables() const { return !unifac_.subgroups.empty(); }

private:
    std::vector<VtprPureRecord> pure_;
    std::unordered_map<int, std::size_t> pure_by_chemid_;
    std::unordered_map<std::string, std::size_t> pure_by_cas_;
    std::unordered_map<std::string, std::vector<std::size_t>> pure_by_name_;

    std::vector<VtprGroupRecord> groups_;
    std::unordered_map<int, std::vector<std::size_t>> groups_by_chemid_;
    std::unordered_map<std::string, std::vector<std::size_t>> groups_by_cas_;
    std::unordered_map<std::string, std::vector<std::size_t>> groups_by_name_;

    UnifacTables unifac_;

    static long long keyPair(int i, int j);

    bool loadPure(const std::filesystem::path& path, Diagnostics* diag);
    bool loadGroups(const std::filesystem::path& path, Diagnostics* diag);
    bool loadUnifacSubgroups(const std::filesystem::path& path, Diagnostics* diag);
    bool loadUnifacInteractions(const std::filesystem::path& path, Diagnostics* diag);
};

} // namespace Data
} // namespace DMThermo

#endif // THERMO_DATA_VTPR_H
