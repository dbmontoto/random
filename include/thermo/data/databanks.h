/**
 * @file databanks.h
 * @brief CSV-backed databanks for DIPPR, PC-SAFT, and BIPs
 */

#ifndef THERMO_DATA_DATABANKS_H
#define THERMO_DATA_DATABANKS_H

#include "diagnostics.h"
#include "records.h"
#include "vtpr.h"
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace DMThermo {
namespace Data {

class DipprCsvDatabank {
public:
    bool load(const std::filesystem::path& csv_path, Diagnostics* diag = nullptr);

    const DipprRecord* findByCHEMID(int chemid) const;
    const DipprRecord* findByCAS(const std::string& cas) const;
    const DipprRecord* findByName(const std::string& name) const;
    const DipprRecord* findByNameUnique(const std::string& name, Diagnostics* diag) const;

    std::size_t size() const { return records_.size(); }
    const std::vector<DipprRecord>& records() const { return records_; }

private:
    std::vector<DipprRecord> records_;
    std::unordered_map<int, std::size_t> by_chemid_;
    std::unordered_map<std::string, std::size_t> by_cas_;
    std::unordered_map<std::string, std::vector<std::size_t>> by_name_;
};

class PCSaftCsvDatabank {
public:
    bool load(const std::filesystem::path& csv_path, Diagnostics* diag = nullptr);

    const PCSaftRecord* findByCHEMID(int chemid) const;
    const PCSaftRecord* findByCAS(const std::string& cas) const;
    const PCSaftRecord* findByName(const std::string& name) const;
    const PCSaftRecord* findByNameUnique(const std::string& name, Diagnostics* diag) const;

    std::size_t size() const { return records_.size(); }
    const std::vector<PCSaftRecord>& records() const { return records_; }

private:
    std::vector<PCSaftRecord> records_;
    std::unordered_map<int, std::size_t> by_chemid_;
    std::unordered_map<std::string, std::size_t> by_cas_;
    std::unordered_map<std::string, std::vector<std::size_t>> by_name_;
};

enum class BipModel {
    PengRobinson,
    VTPR,
    PCSaft
};

class BipCsvDatabank {
public:
    bool load(const std::filesystem::path& csv_path, Diagnostics* diag = nullptr);

    std::optional<double> kij(
        BipModel model,
        const std::string& cas_i,
        const std::string& cas_j,
        double T,
        Diagnostics* diag = nullptr
    ) const;

    std::size_t size() const { return records_.size(); }
    const std::vector<BipRecord>& records() const { return records_; }

private:
    struct PairKey {
        std::string a;
        std::string b;
    };

    struct PairKeyHash {
        std::size_t operator()(const PairKey& k) const;
    };

    struct PairKeyEq {
        bool operator()(const PairKey& x, const PairKey& y) const;
    };

    static PairKey canonicalPair(const std::string& cas_i, const std::string& cas_j);

    std::vector<BipRecord> records_;
    std::unordered_map<PairKey, std::vector<std::size_t>, PairKeyHash, PairKeyEq> by_pair_;
};

struct Databanks {
    DipprCsvDatabank dippr;
    PCSaftCsvDatabank pcsaft;
    BipCsvDatabank bip;
    VtprCsvDatabank vtpr;

    bool loadAllFromDirectory(const std::filesystem::path& dir, Diagnostics* diag = nullptr);
};

} // namespace Data
} // namespace DMThermo

#endif // THERMO_DATA_DATABANKS_H
