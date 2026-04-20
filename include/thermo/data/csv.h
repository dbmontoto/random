/**
 * @file csv.h
 * @brief Minimal CSV reader with basic RFC4180-style quoting
 */

#ifndef THERMO_DATA_CSV_H
#define THERMO_DATA_CSV_H

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace DMThermo {
namespace Data {

struct CsvTable {
    struct HeaderConflict {
        std::string key;
        std::vector<std::size_t> columns; // 0-based column indices in `headers`
    };

    std::vector<std::string> headers;
    std::unordered_map<std::string, std::size_t> header_index;
    std::vector<HeaderConflict> header_conflicts;
    std::vector<std::vector<std::string>> rows;
};

CsvTable readCsvFile(const std::filesystem::path& path, char delimiter = ',');

std::string trim(std::string s);
std::string toLower(std::string s);

std::optional<int> parseInt(const std::string& s);
std::optional<double> parseDouble(const std::string& s);

std::string normalizeName(const std::string& name);
std::string normalizeCAS(const std::string& cas);

} // namespace Data
} // namespace DMThermo

#endif // THERMO_DATA_CSV_H
