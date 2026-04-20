#include "thermo/data/csv.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace DMThermo {
namespace Data {

namespace {

std::string rtrim(std::string s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) {
        s.pop_back();
    }
    return s;
}

std::vector<std::string> parseCsvLine(const std::string& line, char delimiter) {
    std::vector<std::string> out;
    std::string field;
    bool in_quotes = false;

    for (std::size_t i = 0; i < line.size(); ++i) {
        char c = line[i];

        if (in_quotes) {
            if (c == '"') {
                // Escaped quote
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    field.push_back('"');
                    ++i;
                } else {
                    in_quotes = false;
                }
            } else {
                field.push_back(c);
            }
            continue;
        }

        if (c == '"') {
            in_quotes = true;
            continue;
        }

        if (c == delimiter) {
            out.push_back(field);
            field.clear();
            continue;
        }

        // Ignore CR if file has CRLF and we read raw lines.
        if (c == '\r') {
            continue;
        }

        field.push_back(c);
    }

    out.push_back(field);
    return out;
}

bool isBracket(char c) {
    return c == '(' || c == '[' || c == '{';
}

std::string stripUtf8Bom(std::string s) {
    // UTF-8 BOM: EF BB BF
    if (s.size() >= 3 && static_cast<unsigned char>(s[0]) == 0xEF && static_cast<unsigned char>(s[1]) == 0xBB &&
        static_cast<unsigned char>(s[2]) == 0xBF) {
        return s.substr(3);
    }
    return s;
}

std::string stripHeaderSuffixDecoration(std::string s) {
    // If headers include units/comments like "MW (g/mol)" or "Tc [K]", strip the suffix starting at the first bracket.
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (isBracket(s[i])) {
            s.resize(i);
            break;
        }
    }
    return s;
}

std::string canonicalizeHeaderKey(const std::string& raw) {
    std::string s = stripUtf8Bom(trim(raw));
    s = stripHeaderSuffixDecoration(std::move(s));
    s = toLower(trim(std::move(s)));

    std::string out;
    out.reserve(s.size());
    bool prev_us = false;

    for (char ch : s) {
        const unsigned char c = static_cast<unsigned char>(ch);
        if ((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9')) {
            out.push_back(static_cast<char>(c));
            prev_us = false;
            continue;
        }
        if (ch == '_') {
            if (!out.empty() && !prev_us) {
                out.push_back('_');
                prev_us = true;
            }
            continue;
        }
        // Treat any other punctuation/whitespace as a separator.
        if (!out.empty() && !prev_us) {
            out.push_back('_');
            prev_us = true;
        }
    }

    // Trim trailing '_'.
    while (!out.empty() && out.back() == '_') out.pop_back();
    return out;
}

std::string removeUnderscores(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c != '_') out.push_back(c);
    }
    return out;
}

bool matchBipSideField(const std::string& key, std::string& side, std::string& field) {
    // Supports aliases like cas_i <-> i_cas and name_j <-> j_name.
    const auto pos = key.find('_');
    if (pos == std::string::npos) return false;
    const std::string a = key.substr(0, pos);
    const std::string b = key.substr(pos + 1);

    const auto isSide = [](const std::string& s) { return s == "i" || s == "j"; };
    const auto isField = [](const std::string& s) { return s == "cas" || s == "name" || s == "chemid"; };

    if (isSide(a) && isField(b)) {
        side = a;
        field = b;
        return true;
    }
    if (isField(a) && isSide(b)) {
        side = b;
        field = a;
        return true;
    }
    return false;
}

std::vector<std::string> headerKeyAliases(const std::string& key) {
    std::vector<std::string> out;
    out.reserve(6);

    if (!key.empty()) out.push_back(key);

    const std::string no_us = removeUnderscores(key);
    if (!no_us.empty() && no_us != key) out.push_back(no_us);

    // Common Tmin/Tmax spellings.
    if (key == "t_min") out.push_back("min_t");
    if (key == "min_t") out.push_back("t_min");
    if (key == "t_max") out.push_back("max_t");
    if (key == "max_t") out.push_back("t_max");
    if (key == "tmin") {
        out.push_back("t_min");
        out.push_back("min_t");
    }
    if (key == "tmax") {
        out.push_back("t_max");
        out.push_back("max_t");
    }

    // BIP side-field swapping (i_cas <-> cas_i, j_name <-> name_j, etc.).
    std::string side;
    std::string field;
    if (matchBipSideField(key, side, field)) {
        out.push_back(side + "_" + field);
        out.push_back(field + "_" + side);
    }

    // De-duplicate while preserving order.
    std::vector<std::string> uniq;
    uniq.reserve(out.size());
    for (const auto& k : out) {
        if (k.empty()) continue;
        if (std::find(uniq.begin(), uniq.end(), k) == uniq.end()) {
            uniq.push_back(k);
        }
    }
    return uniq;
}

} // namespace

std::string trim(std::string s) {
    s = rtrim(std::move(s));
    std::size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) {
        ++i;
    }
    return s.substr(i);
}

std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::optional<int> parseInt(const std::string& s) {
    std::string t = trim(s);
    if (t.empty()) return std::nullopt;
    try {
        std::size_t pos = 0;
        int v = std::stoi(t, &pos);
        if (pos != t.size()) return std::nullopt;
        return v;
    } catch (...) {
        return std::nullopt;
    }
}

std::optional<double> parseDouble(const std::string& s) {
    std::string t = trim(s);
    if (t.empty()) return std::nullopt;
    try {
        std::size_t pos = 0;
        double v = std::stod(t, &pos);
        if (pos != t.size()) return std::nullopt;
        return v;
    } catch (...) {
        return std::nullopt;
    }
}

std::string normalizeName(const std::string& name) {
    std::string s = toLower(trim(name));
    // Collapse whitespace
    std::string out;
    out.reserve(s.size());
    bool in_space = false;
    for (char c : s) {
        const bool is_space = std::isspace(static_cast<unsigned char>(c)) != 0;
        if (is_space) {
            if (!in_space) out.push_back(' ');
            in_space = true;
            continue;
        }
        in_space = false;
        out.push_back(c);
    }
    // Normalize common punctuation variants.
    std::replace(out.begin(), out.end(), '_', '-');
    return out;
}

std::string normalizeCAS(const std::string& cas) {
    // CAS is typically "NNNNN-NN-N" (digits and dashes). Keep digits/dashes only.
    std::string s = trim(cas);
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if ((c >= '0' && c <= '9') || c == '-') {
            out.push_back(c);
        }
    }
    return out;
}

CsvTable readCsvFile(const std::filesystem::path& path, char delimiter) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open CSV file: " + path.string());
    }

    CsvTable table;
    std::string line;

    // Header
    if (!std::getline(in, line)) {
        throw std::runtime_error("CSV file is empty: " + path.string());
    }

    table.headers = parseCsvLine(line, delimiter);
    for (auto& h : table.headers) {
        h = trim(h);
    }

    // Build a canonical header index robust to common header decoration (units/comments) and punctuation variants.
    // Also expose conflicts so callers can warn (duplicate headers after canonicalization).
    std::unordered_map<std::string, std::vector<std::size_t>> indices_by_key;
    indices_by_key.reserve(table.headers.size() * 2);

    for (std::size_t i = 0; i < table.headers.size(); ++i) {
        const std::string key = canonicalizeHeaderKey(table.headers[i]);
        if (key.empty()) continue;
        for (const auto& alias : headerKeyAliases(key)) {
            indices_by_key[alias].push_back(i);
        }
    }

    for (auto& pair : indices_by_key) {
        auto& idxs = pair.second;
        std::sort(idxs.begin(), idxs.end());
        idxs.erase(std::unique(idxs.begin(), idxs.end()), idxs.end());

        if (!idxs.empty()) {
            table.header_index[pair.first] = idxs.front();
        }
        if (idxs.size() > 1) {
            CsvTable::HeaderConflict c;
            c.key = pair.first;
            c.columns = idxs;
            table.header_conflicts.push_back(std::move(c));
        }
    }

    std::sort(table.header_conflicts.begin(), table.header_conflicts.end(),
        [](const CsvTable::HeaderConflict& a, const CsvTable::HeaderConflict& b) {
            return a.key < b.key;
        });

    // Rows
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto fields = parseCsvLine(line, delimiter);
        // Allow short rows; pad to header size.
        if (fields.size() < table.headers.size()) {
            fields.resize(table.headers.size());
        }
        // Ignore extra columns beyond header.
        if (fields.size() > table.headers.size()) {
            fields.resize(table.headers.size());
        }
        for (auto& f : fields) {
            f = trim(f);
        }
        table.rows.push_back(std::move(fields));
    }

    return table;
}

} // namespace Data
} // namespace DMThermo
