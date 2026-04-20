#include "dmtherm/dmtherm.h"
#include "thermo/core/constants.h"
#include "thermo/core/units.h"
#include "thermo/data/databanks.h"
#include "thermo/data/component_resolver.h"
#include "thermo/equilibrium/critical_estimation.h"
#include "thermo/factory/eos_factory.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <deque>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <limits>

#ifdef _WIN32
#include <windows.h>
#endif

namespace {

struct Options {
    std::unordered_map<std::string, std::string> kv;
    std::unordered_set<std::string> flags;
};

std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string toUpper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

Options parseOptions(int argc, char** argv, int start_index) {
    Options out;
    for (int i = start_index; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") {
            out.flags.insert("help");
            continue;
        }

        if (a.rfind("--", 0) != 0) {
            throw std::invalid_argument("Unexpected positional argument: '" + a + "'");
        }

        auto eq = a.find('=');
        if (eq != std::string::npos) {
            const std::string key = a.substr(2, eq - 2);
            const std::string val = a.substr(eq + 1);
            out.kv[key] = val;
            continue;
        }

        const std::string key = a.substr(2);
        if (i + 1 < argc) {
            const std::string next = argv[i + 1];
            if (next.rfind("--", 0) != 0) {
                out.kv[key] = next;
                ++i;
                continue;
            }
        }

        out.flags.insert(key);
    }
    return out;
}

bool hasFlag(const Options& o, const std::string& k) {
    return o.flags.find(k) != o.flags.end();
}

std::string getStr(const Options& o, const std::string& k, const std::string& def = "") {
    auto it = o.kv.find(k);
    if (it != o.kv.end()) return it->second;
    return def;
}

double getDouble(const Options& o, const std::string& k) {
    auto it = o.kv.find(k);
    if (it == o.kv.end()) throw std::invalid_argument("Missing required option: --" + k);
    return std::stod(it->second);
}

int getInt(const Options& o, const std::string& k, int def) {
    auto it = o.kv.find(k);
    if (it == o.kv.end()) return def;
    return std::stoi(it->second);
}

std::vector<std::string> splitCsv(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        while (!item.empty() && std::isspace(static_cast<unsigned char>(item.front()))) item.erase(item.begin());
        while (!item.empty() && std::isspace(static_cast<unsigned char>(item.back()))) item.pop_back();
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

std::vector<double> parseCsvDoubles(const std::string& s) {
    std::vector<double> out;
    for (const auto& tok : splitCsv(s)) {
        out.push_back(std::stod(tok));
    }
    return out;
}

std::vector<double> normalizeComposition(std::vector<double> x) {
    double sum = 0.0;
    for (double v : x) sum += v;
    if (!(sum > 0.0) || !std::isfinite(sum)) {
        throw std::invalid_argument("Composition sum must be finite and > 0");
    }
    for (double& v : x) v /= sum;
    return x;
}

std::string csvEscape(const std::string& s) {
    bool needs_quotes = false;
    for (char c : s) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            needs_quotes = true;
            break;
        }
    }
    if (!needs_quotes) return s;
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out.push_back(c);
    }
    out.push_back('"');
    return out;
}

std::string jsonEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '"': out += "\\\""; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    std::ostringstream ss;
                    ss << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c);
                    out += ss.str();
                } else {
                    out.push_back(static_cast<char>(c));
                }
        }
    }
    return out;
}

const char* phaseToString(dmtherm_phase_t ph) {
    switch (ph) {
        case DMTHERM_PHASE_VAPOR: return "Vapor";
        case DMTHERM_PHASE_LIQUID: return "Liquid";
        case DMTHERM_PHASE_LIQUID1: return "Liquid1";
        case DMTHERM_PHASE_LIQUID2: return "Liquid2";
        case DMTHERM_PHASE_SUPERCRITICAL: return "Supercritical";
        case DMTHERM_PHASE_SOLID: return "Solid";
        case DMTHERM_PHASE_UNKNOWN:
        default: return "Unknown";
    }
}

dmtherm_root_selection_t parseRootSelection(std::string s) {
    s = toLower(std::move(s));
    if (s == "vapor") return DMTHERM_ROOT_VAPOR;
    if (s == "liquid") return DMTHERM_ROOT_LIQUID;
    if (s == "stable") return DMTHERM_ROOT_STABLE;
    throw std::invalid_argument("Invalid --root value (expected vapor|liquid|stable)");
}

dmtherm_property_method_t parsePropertyMethod(std::string s) {
    s = toLower(std::move(s));
    if (s == "auto") return DMTHERM_PROPERTY_METHOD_AUTO;
    if (s == "eos_only" || s == "eos") return DMTHERM_PROPERTY_METHOD_EOS_ONLY;
    if (s == "dippr_only" || s == "dippr") return DMTHERM_PROPERTY_METHOD_DIPPR_ONLY;
    if (s == "hybrid_vapor_eos_liquid_dippr" || s == "hybrid") return DMTHERM_PROPERTY_METHOD_HYBRID_VAPOR_EOS_LIQUID_DIPPR;
    throw std::invalid_argument("Invalid --property-method (expected auto|eos_only|dippr_only|hybrid_vapor_eos_liquid_dippr)");
}

const char* pressureUnitToken(dmtherm_pressure_unit_t unit) {
    switch (unit) {
        case DMTHERM_PRESSURE_UNIT_PA: return "PA";
        case DMTHERM_PRESSURE_UNIT_PA_ABS: return "PA_ABS";
        case DMTHERM_PRESSURE_UNIT_PA_GAUGE: return "PA_GAUGE";
        case DMTHERM_PRESSURE_UNIT_KPA: return "KPA";
        case DMTHERM_PRESSURE_UNIT_KPA_ABS: return "KPA_ABS";
        case DMTHERM_PRESSURE_UNIT_KPA_GAUGE: return "KPA_GAUGE";
        case DMTHERM_PRESSURE_UNIT_MPA: return "MPA";
        case DMTHERM_PRESSURE_UNIT_MPA_ABS: return "MPA_ABS";
        case DMTHERM_PRESSURE_UNIT_MPA_GAUGE: return "MPA_GAUGE";
        case DMTHERM_PRESSURE_UNIT_BAR: return "BAR";
        case DMTHERM_PRESSURE_UNIT_BARG: return "BARG";
        case DMTHERM_PRESSURE_UNIT_PSI: return "PSI";
        case DMTHERM_PRESSURE_UNIT_PSIA: return "PSIA";
        case DMTHERM_PRESSURE_UNIT_PSIG: return "PSIG";
        default: return "UNKNOWN";
    }
}

dmtherm_pressure_unit_t parsePressureUnit(std::string s) {
    s = toUpper(std::move(s));
    if (s == "PA") return DMTHERM_PRESSURE_UNIT_PA;
    if (s == "PA_ABS") return DMTHERM_PRESSURE_UNIT_PA_ABS;
    if (s == "PA_GAUGE") return DMTHERM_PRESSURE_UNIT_PA_GAUGE;
    if (s == "KPA") return DMTHERM_PRESSURE_UNIT_KPA;
    if (s == "KPA_ABS") return DMTHERM_PRESSURE_UNIT_KPA_ABS;
    if (s == "KPA_GAUGE") return DMTHERM_PRESSURE_UNIT_KPA_GAUGE;
    if (s == "MPA") return DMTHERM_PRESSURE_UNIT_MPA;
    if (s == "MPA_ABS") return DMTHERM_PRESSURE_UNIT_MPA_ABS;
    if (s == "MPA_GAUGE") return DMTHERM_PRESSURE_UNIT_MPA_GAUGE;
    if (s == "BAR") return DMTHERM_PRESSURE_UNIT_BAR;
    if (s == "BARG") return DMTHERM_PRESSURE_UNIT_BARG;
    if (s == "PSI") return DMTHERM_PRESSURE_UNIT_PSI;
    if (s == "PSIA") return DMTHERM_PRESSURE_UNIT_PSIA;
    if (s == "PSIG") return DMTHERM_PRESSURE_UNIT_PSIG;
    throw std::invalid_argument(
        "Invalid --P-unit (expected PA|PA_ABS|PA_GAUGE|KPA|KPA_ABS|KPA_GAUGE|"
        "MPA|MPA_ABS|MPA_GAUGE|BAR|BARG|PSI|PSIA|PSIG)");
}

double convertPressure(double value, dmtherm_pressure_unit_t from_unit, dmtherm_pressure_unit_t to_unit) {
    if (from_unit == to_unit) return value;
    double out = 0.0;
    const auto st = dmtherm_convert_pressure(value, from_unit, to_unit, &out);
    if (st == DMTHERM_STATUS_OK) return out;
    if (st == DMTHERM_STATUS_INVALID_ARGUMENT) {
        throw std::invalid_argument(
            std::string("Incompatible pressure conversion: ") + pressureUnitToken(from_unit) + " -> " + pressureUnitToken(to_unit));
    }
    throw std::runtime_error("dmtherm_convert_pressure failed with status " + std::to_string(static_cast<int>(st)));
}

double convertPressureToPaAbs(double value, dmtherm_pressure_unit_t from_unit) {
    return convertPressure(value, from_unit, DMTHERM_PRESSURE_UNIT_PA_ABS);
}

std::optional<double> optDouble(const dmtherm_optional_double_t& v) {
    if (v.has_value) return v.value;
    return std::nullopt;
}

std::vector<double> getComponentMWs(dmtherm_system_t* sys, size_t num_components) {
    std::vector<double> mws(num_components, 0.0);
    if (!sys || num_components == 0) return mws;
    const auto st = dmtherm_system_get_component_mws(sys, mws.data(), mws.size());
    if (st != DMTHERM_STATUS_OK) {
        throw std::runtime_error("dmtherm_system_get_component_mws failed (status=" + std::to_string(st) + "): " +
                                 std::string(dmtherm_system_last_error(sys)));
    }
    return mws;
}

double averageMW(const std::vector<double>& x, const std::vector<double>& component_mws) {
    if (x.size() != component_mws.size()) {
        throw std::invalid_argument("averageMW: composition length must match component_mws length");
    }
    double mw = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        mw += x[i] * component_mws[i];
    }
    return mw;
}

std::optional<double> molarToMass(double v_molar, double mw_g_per_mol) {
    if (!(std::isfinite(mw_g_per_mol) && mw_g_per_mol > 0.0)) return std::nullopt;
    return v_molar * 1000.0 / mw_g_per_mol;
}

const char* propertyMethodToString(dmtherm_property_method_t m) {
    switch (m) {
        case DMTHERM_PROPERTY_METHOD_AUTO: return "auto";
        case DMTHERM_PROPERTY_METHOD_EOS_ONLY: return "eos_only";
        case DMTHERM_PROPERTY_METHOD_DIPPR_ONLY: return "dippr_only";
        case DMTHERM_PROPERTY_METHOD_HYBRID_VAPOR_EOS_LIQUID_DIPPR: return "hybrid_vapor_eos_liquid_dippr";
        default: return "unknown";
    }
}

void jsonWriteOptionalDouble(std::ostream& os, const dmtherm_optional_double_t& v) {
    if (v.has_value) os << v.value;
    else os << "null";
}

void jsonWriteOptionalDouble(std::ostream& os, const std::optional<double>& v) {
    if (v.has_value()) os << v.value();
    else os << "null";
}

void jsonWriteStringOrNull(std::ostream& os, const char* s) {
    if (s && s[0] != '\0') os << "\"" << jsonEscape(s) << "\"";
    else os << "null";
}

void jsonWriteDoubleArray(std::ostream& os, const double* a, size_t n) {
    os << "[";
    for (size_t i = 0; i < n; ++i) {
        if (i) os << ",";
        os << a[i];
    }
    os << "]";
}

void jsonWriteStringArray(std::ostream& os, const std::vector<std::string>& a) {
    os << "[";
    for (size_t i = 0; i < a.size(); ++i) {
        if (i) os << ",";
        os << "\"" << jsonEscape(a[i]) << "\"";
    }
    os << "]";
}

const char* canonicalFlashInputKey(const char* key) {
    if (key == nullptr) return nullptr;
    const std::string token = toUpper(std::string(key));
    if (token == "T") return "temperature_k";
    if (token == "P") return "pressure_pa_abs";
    if (token == "V") return "molar_volume_m3_per_mol";
    if (token == "H") return "enthalpy_j_per_mol";
    if (token == "S") return "entropy_j_per_mol_k";
    return nullptr;
}

void jsonWriteTPPointCanonical(
    std::ostream& os,
    double T,
    double P,
    const std::string& root_name,
    dmtherm_status_t st,
    bool success,
    dmtherm_phase_t ph,
    double rho,
    double Z,
    double h_dep,
    double s_dep,
    double g_dep,
    double H,
    double S,
    double G,
    const double* ln_phi,
    size_t ln_phi_len,
    const std::string& message)
{
    os << "{";
    os << "\"temperature_k\":" << T << ",";
    os << "\"pressure_pa_abs\":" << P << ",";
    os << "\"root\":\"" << jsonEscape(root_name) << "\",";
    os << "\"status\":" << static_cast<int>(st) << ",";
    os << "\"success\":" << (success ? "true" : "false") << ",";
    os << "\"phase\":\"" << jsonEscape(phaseToString(ph)) << "\",";
    os << "\"molar_density_mol_per_m3\":" << rho << ",";
    os << "\"compressibility\":" << Z << ",";
    os << "\"enthalpy_departure_j_per_mol\":" << h_dep << ",";
    os << "\"entropy_departure_j_per_mol_k\":" << s_dep << ",";
    os << "\"gibbs_departure_j_per_mol\":" << g_dep << ",";
    os << "\"enthalpy_j_per_mol\":" << H << ",";
    os << "\"entropy_j_per_mol_k\":" << S << ",";
    os << "\"gibbs_j_per_mol\":" << G << ",";
    os << "\"ln_fugacity_coeff\":";
    if (ln_phi && ln_phi_len) jsonWriteDoubleArray(os, ln_phi, ln_phi_len);
    else os << "[]";
    os << ",";
    os << "\"message\":\"" << jsonEscape(message) << "\"";
    os << "}";
}

void printGlobalHelp(std::ostream& os) {
    os << "dmtherm - DMTherm console tool (C ABI)\n\n";
    os << "Usage:\n";
    os << "  dmtherm                         (interactive shell)\n";
    os << "  dmtherm --version\n";
    os << "  dmtherm shell                   (interactive mode)\n";
    os << "  dmtherm list components         (browse available components from CSV databanks)\n";
    os << "  dmtherm list models             (list available EOS models)\n";
    os << "  dmtherm eval critical --components <csv> [options]\n";
    os << "  dmtherm eval tp-phase --T <K> --P <value> --components <csv> --z <csv> [options]\n";
    os << "  dmtherm eval tp-roots --T <K> --P <value> --components <csv> --z <csv> [options]\n";
    os << "  dmtherm flash pt --T <K> --P <value> --components <csv> --z <csv> [options]\n";
    os << "  dmtherm flash tv --T <K> --V <m^3/mol> --components <csv> --z <csv> [options]\n";
    os << "  dmtherm flash ph --P <value> --H <J/mol> --components <csv> --z <csv> [options]\n";
    os << "  dmtherm flash ps --P <value> --S <J/mol/K> --components <csv> --z <csv> [options]\n\n";
    os << "  dmtherm sweep tp-grid --components <csv> --z <csv> --T-min <K> --T-max <K> --nT <N> --P-min <value> --P-max <value> --nP <N> [options]\n\n";
    os << "  dmtherm sweep isotherm --T <K> --components <csv> --z <csv> --P-min <value> --P-max <value> --nP <N> [options]\n";
    os << "  dmtherm sweep isobar --P <value> --components <csv> --z <csv> --T-min <K> --T-max <K> --nT <N> [options]\n\n";
    os << "Common options:\n";
    os << "  --databanks-dir <dir>        (default: data)\n";
    os << "  --eos <name>                 (default: Peng-Robinson)\n";
    os << "  --flash-method <name>        (default: Auto)\n";
    os << "  --property-method <name>     (default: eos_only)\n";
    os << "  --P-unit <unit>              (pressure input unit; default: PA_ABS)\n";
    os << "  --require-critical           (enforce Tc/Pc/omega)\n";
    os << "  --no-require-ideal-gas-cp    (allow missing ideal Cp; H/S/G may be unavailable)\n";
    os << "\n";
}

void printListComponentsHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm list components [options]\n\n";
    os << "Options:\n";
    os << "  --databanks-dir <dir>        (default: data)\n";
    os << "  --filter <text>              (case-insensitive; matches name/CAS/CHEMID)\n";
    os << "  --limit <N>                  (default: 50)\n";
    os << "  --csv                        (print CSV instead of text)\n";
}

void printListModelsHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm list models [options]\n\n";
    os << "Options:\n";
    os << "  --csv                        (print CSV instead of text)\n";
}

// --------------------------------------------------------------------------------------
// Interactive shell
// --------------------------------------------------------------------------------------

struct ShellSession {
    std::string databanks_dir = "data";
    std::string eos_type = "Peng-Robinson";
    std::string flash_method = "Auto";
    std::string property_method = "eos_only";
    std::string pressure_unit = "PA_ABS";
    bool require_critical = false;
    bool require_ideal_gas_cp = true;
    bool always_show_menu = true;
    bool sticky_menu = true;
    int ui_log_limit = 500;
    int ui_tail_lines = 40;
    std::vector<std::string> components;
    std::vector<double> z;

    std::deque<std::string> ui_log;

    bool databanks_loaded = false;
    std::string databanks_loaded_dir;
    DMThermo::Data::Databanks databanks;
    DMThermo::Data::Diagnostics databanks_diag;
};

std::string joinCsv(const std::vector<std::string>& items) {
    std::ostringstream oss;
    for (size_t i = 0; i < items.size(); ++i) {
        if (i) oss << ",";
        oss << items[i];
    }
    return oss.str();
}

std::string joinCsvDoubles(const std::vector<double>& items) {
    std::ostringstream oss;
    oss << std::setprecision(17);
    for (size_t i = 0; i < items.size(); ++i) {
        if (i) oss << ",";
        oss << items[i];
    }
    return oss.str();
}

std::string trimCopy(std::string s) {
    auto isSpace = [](unsigned char c) { return std::isspace(c) != 0; };
    while (!s.empty() && isSpace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
    while (!s.empty() && isSpace(static_cast<unsigned char>(s.back()))) s.pop_back();
    return s;
}

std::vector<std::string> tokenizeCommandLine(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    enum class Quote { None, Single, Double } q = Quote::None;

    auto flush = [&]() {
        if (!cur.empty()) {
            out.push_back(cur);
            cur.clear();
        }
    };

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (q == Quote::None) {
            if (std::isspace(static_cast<unsigned char>(c))) {
                flush();
                continue;
            }
            if (c == '\'') {
                q = Quote::Single;
                continue;
            }
            if (c == '"') {
                q = Quote::Double;
                continue;
            }
            cur.push_back(c);
            continue;
        }

        const char qch = (q == Quote::Single) ? '\'' : '"';
        if (c == '\\' && i + 1 < line.size()) {
            const char next = line[i + 1];
            if (next == qch || next == '\\') {
                cur.push_back(next);
                ++i;
                continue;
            }
        }
        if (c == qch) {
            q = Quote::None;
            continue;
        }
        cur.push_back(c);
    }
    flush();
    return out;
}

bool hasOpt(const std::vector<std::string>& args, const std::string& opt_with_dashes) {
    for (const auto& a : args) {
        if (a == opt_with_dashes) return true;
        if (a.rfind(opt_with_dashes + "=", 0) == 0) return true;
    }
    return false;
}

std::vector<std::string> injectSessionDefaults(std::vector<std::string> args, const ShellSession& s) {
    if (args.empty()) return args;

    const std::string cmd1 = toLower(args[0]);
    const bool needs_system = (cmd1 == "eval" || cmd1 == "flash" || cmd1 == "sweep");
    if (!needs_system) return args;

    if (!hasOpt(args, "--databanks-dir") && !s.databanks_dir.empty()) {
        args.push_back("--databanks-dir");
        args.push_back(s.databanks_dir);
    }
    if (!hasOpt(args, "--eos") && !s.eos_type.empty()) {
        args.push_back("--eos");
        args.push_back(s.eos_type);
    }
    if (!hasOpt(args, "--flash-method") && !s.flash_method.empty()) {
        args.push_back("--flash-method");
        args.push_back(s.flash_method);
    }
    if (!hasOpt(args, "--property-method") && !s.property_method.empty()) {
        args.push_back("--property-method");
        args.push_back(s.property_method);
    }
    if (!hasOpt(args, "--P-unit") && !s.pressure_unit.empty()) {
        args.push_back("--P-unit");
        args.push_back(s.pressure_unit);
    }
    if (s.require_critical && !hasOpt(args, "--require-critical")) {
        args.push_back("--require-critical");
    }
    if (!s.require_ideal_gas_cp && !hasOpt(args, "--no-require-ideal-gas-cp")) {
        args.push_back("--no-require-ideal-gas-cp");
    }

    if (!hasOpt(args, "--components") && !s.components.empty()) {
        args.push_back("--components");
        args.push_back(joinCsv(s.components));
    }
    if (!hasOpt(args, "--z") && !s.z.empty()) {
        args.push_back("--z");
        args.push_back(joinCsvDoubles(s.z));
    }
    return args;
}

void printShellMenu(std::ostream& os, const ShellSession& s) {
    os << "Menu:\n";
    os << "  Session: "
       << "db=" << (s.databanks_dir.empty() ? "<unset>" : s.databanks_dir)
       << " eos=" << (s.eos_type.empty() ? "<unset>" : s.eos_type)
       << " comps=" << (s.components.empty() ? "0" : std::to_string(s.components.size()))
       << " z=" << (s.z.empty() ? "0" : std::to_string(s.z.size()))
       << " flash=" << (s.flash_method.empty() ? "<unset>" : s.flash_method)
       << " props=" << (s.property_method.empty() ? "<unset>" : s.property_method)
       << " p-unit=" << (s.pressure_unit.empty() ? "<unset>" : s.pressure_unit)
       << "\n";
    os << "  Set:   [2]DB  [3]EOS  [4]Comps  [5]z   Eval:  [6]TPphase  [7]TProots\n";
    os << "  Flash: [8]PT  [9]TV  [10]PH  [11]PS   Sweep: [12]TPgrid  [13]Isotherm  [14]Isobar\n";
    os << "  Util:  [1]Show  [15]Run  [16]Models  [17]Components  [18]Critical  [0]Exit\n";
}

void printShellHelp(std::ostream& os) {
    os << "Shell commands:\n";
    os << "  menu                 Show interactive menu\n";
    os << "  show                 Show current session defaults\n";
    os << "  set <key> <value>    Set a session default\n";
    os << "  clear                Clear the on-screen output log (sticky UI mode)\n";
    os << "  list components      Print available components from CSV databanks\n";
    os << "  list models          Print available EOS models (values for --eos)\n";
    os << "  run <cmd...>         Run a dmtherm command (session defaults injected)\n";
    os << "  <cmd...>             You can also type dmtherm commands directly\n";
    os << "  exit | quit          Exit the shell\n\n";
    os << "Keys for 'set':\n";
    os << "  databanks-dir | eos | flash-method | property-method | p-unit\n";
    os << "  require-critical (true/false)\n";
    os << "  require-ideal-gas-cp (true/false)\n";
    os << "  always-show-menu (true/false)\n";
    os << "  sticky-menu (true/false)\n";
    os << "  ui-log-limit (N lines)\n";
    os << "  ui-tail-lines (N lines)\n";
    os << "  components (csv) | z (csv)\n";
}

void printShellSession(std::ostream& os, const ShellSession& s) {
    os << "session:\n";
    os << "  databanks-dir: " << (s.databanks_dir.empty() ? "<empty>" : s.databanks_dir) << "\n";
    os << "  eos: " << (s.eos_type.empty() ? "<empty>" : s.eos_type) << "\n";
    os << "  flash-method: " << (s.flash_method.empty() ? "<empty>" : s.flash_method) << "\n";
    os << "  property-method: " << (s.property_method.empty() ? "<empty>" : s.property_method) << "\n";
    os << "  p-unit: " << (s.pressure_unit.empty() ? "<empty>" : s.pressure_unit) << "\n";
    os << "  require-critical: " << (s.require_critical ? "true" : "false") << "\n";
    os << "  require-ideal-gas-cp: " << (s.require_ideal_gas_cp ? "true" : "false") << "\n";
    os << "  always-show-menu: " << (s.always_show_menu ? "true" : "false") << "\n";
    os << "  sticky-menu: " << (s.sticky_menu ? "true" : "false") << "\n";
    os << "  ui-log-limit: " << s.ui_log_limit << "\n";
    os << "  ui-tail-lines: " << s.ui_tail_lines << "\n";
    os << "  components: " << (s.components.empty() ? "<unset>" : joinCsv(s.components)) << "\n";
    os << "  z: " << (s.z.empty() ? "<unset>" : joinCsvDoubles(s.z)) << "\n";
}

bool parseBool(const std::string& s, bool& out) {
    const std::string v = toLower(trimCopy(s));
    if (v == "1" || v == "true" || v == "yes" || v == "y" || v == "on") {
        out = true;
        return true;
    }
    if (v == "0" || v == "false" || v == "no" || v == "n" || v == "off") {
        out = false;
        return true;
    }
    return false;
}

std::string shellPressureUnitToken(const ShellSession& s) {
    if (s.pressure_unit.empty()) return "PA_ABS";
    return s.pressure_unit;
}

double pressureDefaultInSessionUnit(double pa_abs_default, const ShellSession& s) {
    try {
        return convertPressure(
            pa_abs_default,
            DMTHERM_PRESSURE_UNIT_PA_ABS,
            parsePressureUnit(shellPressureUnitToken(s)));
    } catch (...) {
        return pa_abs_default;
    }
}

std::string pressurePromptLabel(const std::string& quantity, const ShellSession& s) {
    return quantity + " [" + shellPressureUnitToken(s) + "]";
}

std::string promptLine(const std::string& label, const std::string& def = "") {
    std::cout << label;
    if (!def.empty()) std::cout << " [" << def << "]";
    std::cout << ": ";
    std::string line;
    if (!std::getline(std::cin, line)) return def;
    line = trimCopy(line);
    return line.empty() ? def : line;
}

double promptDouble(const std::string& label, double def) {
    for (;;) {
        const std::string s = promptLine(label, std::to_string(def));
        try {
            return std::stod(s);
        } catch (...) {
            std::cout << "Invalid number. Try again.\n";
        }
    }
}

int promptInt(const std::string& label, int def) {
    for (;;) {
        const std::string s = promptLine(label, std::to_string(def));
        try {
            return std::stoi(s);
        } catch (...) {
            std::cout << "Invalid integer. Try again.\n";
        }
    }
}

bool promptYesNo(const std::string& label, bool def) {
    for (;;) {
        const std::string s = promptLine(label + " (y/n)", def ? "y" : "n");
        bool v = def;
        if (parseBool(s, v)) return v;
        std::cout << "Invalid value. Use y/n/true/false.\n";
    }
}

std::string promptEOSFromList(const std::string& current) {
    const auto types = DMThermo::Factory::EOSFactory::registeredTypes();
    std::cout << "Available EOS models:\n";
    for (size_t i = 0; i < types.size(); ++i) {
        std::cout << "  [" << (i + 1) << "] " << types[i] << "\n";
    }
    std::cout << "  [0] Enter custom\n";

    for (;;) {
        std::cout << "Current: " << (current.empty() ? "<unset>" : current) << "\n";
        const std::string pick = promptLine("Select EOS model # (blank=keep current)", "");
        if (pick.empty()) return current;

        try {
            size_t pos = 0;
            const int idx = std::stoi(pick, &pos);
            if (pos != pick.size()) {
                // Allow typing a model name directly.
                return pick;
            }
            if (idx == 0) {
                return promptLine("eos", current);
            }
            if (idx > 0 && static_cast<size_t>(idx) <= types.size()) {
                return types[static_cast<size_t>(idx) - 1];
            }
        } catch (...) {
            // Allow typing a model name directly.
            return pick;
        }

        std::cout << "Invalid selection.\n";
    }
}

std::string promptModelsArgFromList() {
    const auto types = DMThermo::Factory::EOSFactory::registeredTypes();
    std::cout << "Available EOS models:\n";
    for (size_t i = 0; i < types.size(); ++i) {
        std::cout << "  [" << (i + 1) << "] " << types[i] << "\n";
    }
    std::cout << "  [a] all\n";

    for (;;) {
        const std::string pick = promptLine("Select models (csv indices or 'a')", "a");
        const std::string pl = toLower(trimCopy(pick));
        if (pl.empty() || pl == "a" || pl == "all") return "all";

        std::vector<std::string> selected;
        for (const auto& tok : splitCsv(pick)) {
            try {
                const int idx = std::stoi(tok);
                if (idx <= 0) continue;
                if (static_cast<size_t>(idx) > types.size()) continue;
                const auto& name = types[static_cast<size_t>(idx) - 1];
                if (std::find(selected.begin(), selected.end(), name) == selected.end()) {
                    selected.push_back(name);
                }
            } catch (...) {
                // ignore
            }
        }

        if (selected.empty()) {
            std::cout << "No valid indices selected.\n";
            continue;
        }
        return joinCsv(selected);
    }
}

bool shellUsesStickyMenu(const ShellSession& s) {
    return s.always_show_menu && s.sticky_menu;
}

void uiLogClear(ShellSession& s) {
    s.ui_log.clear();
}

void uiLogAppendLine(ShellSession& s, std::string line) {
    if (!line.empty() && line.back() == '\r') line.pop_back();

    // Drop trailing empty lines to reduce spam.
    if (line.empty() && !s.ui_log.empty() && s.ui_log.back().empty()) return;

    s.ui_log.push_back(std::move(line));
    while (static_cast<int>(s.ui_log.size()) > std::max(1, s.ui_log_limit)) {
        s.ui_log.pop_front();
    }
}

void uiLogAppendText(ShellSession& s, const std::string& text) {
    std::string line;
    line.reserve(256);
    for (char c : text) {
        if (c == '\n') {
            uiLogAppendLine(s, std::move(line));
            line.clear();
            continue;
        }
        line.push_back(c);
    }
    if (!line.empty()) uiLogAppendLine(s, std::move(line));
}

void shellPrintln(ShellSession& s, const std::string& line) {
    if (shellUsesStickyMenu(s)) {
        uiLogAppendLine(s, line);
    } else {
        std::cout << line << "\n";
    }
}

bool enableAnsiEscapeCodes() {
#ifdef _WIN32
    static bool attempted = false;
    static bool enabled = false;
    if (attempted) return enabled;
    attempted = true;

    HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    if (h == INVALID_HANDLE_VALUE) return false;

    DWORD mode = 0;
    if (!GetConsoleMode(h, &mode)) return false;

    // Try to enable ANSI processing; if it fails, we fall back to printing newlines.
    mode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    if (!SetConsoleMode(h, mode)) return false;
    enabled = true;
    return true;
#else
    return true;
#endif
}

void clearScreen() {
    if (enableAnsiEscapeCodes()) {
        // Clear screen + move cursor home.
        std::cout << "\x1b[2J\x1b[H";
        return;
    }

#ifdef _WIN32
    // Try Win32 console clearing (works even without ANSI).
    HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    if (h != INVALID_HANDLE_VALUE) {
        CONSOLE_SCREEN_BUFFER_INFO csbi;
        if (GetConsoleScreenBufferInfo(h, &csbi)) {
            const DWORD cell_count = static_cast<DWORD>(csbi.dwSize.X) * static_cast<DWORD>(csbi.dwSize.Y);
            DWORD written = 0;
            const COORD home{0, 0};
            const bool ok_chars = FillConsoleOutputCharacterA(h, ' ', cell_count, home, &written) != 0;
            const bool ok_attr = FillConsoleOutputAttribute(h, csbi.wAttributes, cell_count, home, &written) != 0;
            if (ok_chars && ok_attr) {
                (void)SetConsoleCursorPosition(h, home);
                return;
            }
        }
    }
#endif

    // Last resort fallback (keeps things usable in non-console terminals / redirected output).
    for (int i = 0; i < 80; ++i) std::cout << "\n";
}

void renderStickyShellUI(const ShellSession& s) {
    clearScreen();

    const int tail = std::max(0, s.ui_tail_lines);
    const std::size_t total = s.ui_log.size();
    const std::size_t start = (tail > 0 && total > static_cast<std::size_t>(tail)) ? (total - static_cast<std::size_t>(tail)) : 0;

    for (std::size_t i = start; i < total; ++i) {
        std::cout << s.ui_log[i] << "\n";
    }
    if (total > 0) std::cout << "\n";

    printShellMenu(std::cout, s);
}

bool loadDatabanksForShell(ShellSession& s) {
    if (s.databanks_loaded && s.databanks_loaded_dir == s.databanks_dir) return true;

    s.databanks_loaded = false;
    s.databanks_loaded_dir.clear();
    s.databanks = DMThermo::Data::Databanks{};
    s.databanks_diag = DMThermo::Data::Diagnostics{};

    if (s.databanks_dir.empty() || !std::filesystem::exists(s.databanks_dir)) {
        shellPrintln(s, "databanks dir does not exist: " + s.databanks_dir);
        return false;
    }

    const bool ok = s.databanks.loadAllFromDirectory(s.databanks_dir, &s.databanks_diag);
    if (!ok || s.databanks_diag.hasErrors()) {
        shellPrintln(s, "Failed to load databanks from: " + s.databanks_dir);
        for (const auto& d : s.databanks_diag.items()) {
            if (d.level == DMThermo::Data::DiagnosticLevel::Error) {
                shellPrintln(s, "  error: " + d.message);
            }
        }
        return false;
    }

    s.databanks_loaded = true;
    s.databanks_loaded_dir = s.databanks_dir;
    return true;
}

bool isPCSaftEos(const std::string& eos_type) {
    std::string s = toLower(eos_type);
    s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c) {
        return c == '-' || c == '_' || std::isspace(c);
    }), s.end());
    return (s == "pcsaft");
}

bool recordMatchesSearch(const DMThermo::Data::DipprRecord& r, const std::string& q_lower) {
    if (q_lower.empty()) return true;
    auto contains = [&](const std::string& s) {
        const std::string sl = toLower(s);
        return sl.find(q_lower) != std::string::npos;
    };

    if (contains(r.key.name)) return true;
    if (contains(r.key.cas)) return true;
    if (r.key.chemid.has_value()) {
        const std::string cs = std::to_string(r.key.chemid.value());
        if (contains(cs)) return true;
    }
    return false;
}

bool recordSatisfiesSessionRequirements(const DMThermo::Data::Databanks& db, const DMThermo::Data::DipprRecord& r, const ShellSession& s) {
    if (s.require_ideal_gas_cp && !r.ideal_gas_cp.has_value()) return false;
    if (s.require_critical) {
        if (!r.Tc.has_value() || !r.Pc.has_value() || !r.omega.has_value()) return false;
    }
    if (isPCSaftEos(s.eos_type)) {
        const DMThermo::Data::PCSaftRecord* pr = nullptr;
        if (r.key.chemid.has_value()) pr = db.pcsaft.findByCHEMID(r.key.chemid.value());
        if (!pr && !r.key.cas.empty()) pr = db.pcsaft.findByCAS(r.key.cas);
        if (!pr && !r.key.name.empty()) pr = db.pcsaft.findByName(r.key.name);
        if (!pr) return false;
    }
    return true;
}

std::string recordIdentifier(const DMThermo::Data::DipprRecord& r) {
    if (!r.key.cas.empty()) return r.key.cas;
    if (r.key.chemid.has_value()) return std::to_string(r.key.chemid.value());
    return r.key.name;
}

void setEqualComposition(ShellSession& s) {
    s.z.clear();
    if (s.components.empty()) return;
    const double v = 1.0 / static_cast<double>(s.components.size());
    s.z.assign(s.components.size(), v);
}

void browseAndSelectComponents(ShellSession& s) {
    if (!loadDatabanksForShell(s)) return;

    if (!s.components.empty()) {
        const bool replace = promptYesNo("Replace existing components", true);
        if (replace) {
            s.components.clear();
            s.z.clear();
        }
    }

    std::cout << "Search matches are filtered by current session requirements:\n";
    std::cout << "  require-critical=" << (s.require_critical ? "true" : "false") << "\n";
    std::cout << "  require-ideal-gas-cp=" << (s.require_ideal_gas_cp ? "true" : "false") << "\n";
    std::cout << "  eos=" << s.eos_type << (isPCSaftEos(s.eos_type) ? " (PC-SAFT parameters required)" : "") << "\n";

    for (;;) {
        const std::string query = promptLine("Search (name/CAS/CHEMID), blank=all, 'q' to quit", "");
        if (toLower(query) == "q") break;

        const std::string ql = toLower(query);
        std::vector<const DMThermo::Data::DipprRecord*> matches;
        matches.reserve(64);

        for (const auto& r : s.databanks.dippr.records()) {
            if (!recordMatchesSearch(r, ql)) continue;
            if (!recordSatisfiesSessionRequirements(s.databanks, r, s)) continue;
            matches.push_back(&r);
        }

        if (matches.empty()) {
            std::cout << "No matches.\n";
            continue;
        }

        const size_t limit = 20;
        const size_t shown = std::min(limit, matches.size());
        std::cout << "\nMatches (showing " << shown << " of " << matches.size() << "):\n";
        for (size_t i = 0; i < shown; ++i) {
            const auto& r = *matches[i];
            std::cout << "  [" << (i + 1) << "] " << r.key.name;
            if (!r.key.cas.empty()) std::cout << "  CAS=" << r.key.cas;
            if (r.key.chemid.has_value()) std::cout << "  CHEMID=" << r.key.chemid.value();
            if (r.MW.has_value()) std::cout << "  molecular_weight=" << r.MW.value();
            std::cout << "\n";
        }
        if (matches.size() > shown) {
            std::cout << "Refine your search to see more.\n";
        }

        const std::string pick = promptLine("Select indices to add (csv), 'a' for all shown, blank to search again", "");
        if (pick.empty()) continue;

        std::vector<size_t> sel;
        const std::string pick_l = toLower(pick);
        if (pick_l == "a" || pick_l == "all") {
            for (size_t i = 0; i < shown; ++i) sel.push_back(i + 1);
        } else {
            for (const auto& tok : splitCsv(pick)) {
                try {
                    const int v = std::stoi(tok);
                    if (v <= 0) continue;
                    sel.push_back(static_cast<size_t>(v));
                } catch (...) {
                    // ignore token
                }
            }
        }

        if (sel.empty()) {
            std::cout << "No valid indices selected.\n";
            continue;
        }

        size_t added = 0;
        for (size_t idx1 : sel) {
            if (idx1 == 0 || idx1 > shown) continue;
            const auto& r = *matches[idx1 - 1];
            const std::string id = recordIdentifier(r);
            if (id.empty()) continue;
            if (std::find(s.components.begin(), s.components.end(), id) != s.components.end()) continue;
            s.components.push_back(id);
            ++added;
        }

        if (added == 0) {
            std::cout << "No new components added.\n";
        } else {
            setEqualComposition(s);
            std::cout << "Added " << added << " component(s). components=" << joinCsv(s.components) << "\n";
            std::cout << "Updated z to equal fractions. Use menu [5] or 'set z ...' to change.\n";
        }

        const bool more = promptYesNo("Add more components", false);
        if (!more) break;
    }
}

bool ensureSessionHasComponents(ShellSession& s) {
    if (!s.components.empty()) return true;

    const bool browse = promptYesNo("Browse components from databanks", true);
    if (browse) {
        browseAndSelectComponents(s);
    }
    if (s.components.empty()) {
        const std::string comps = promptLine("components (csv)", "");
        if (!comps.empty()) {
            s.components = splitCsv(comps);
        }
    }
    if (s.components.empty()) {
        shellPrintln(s, "components are required.");
        return false;
    }
    if (!s.z.empty() && s.z.size() != s.components.size()) {
        s.z.clear();
    }
    return true;
}

bool ensureSessionHasMixture(ShellSession& s) {
    if (!s.components.empty() && !s.z.empty() && s.z.size() == s.components.size()) return true;

    if (s.components.empty()) {
        const bool browse = promptYesNo("Browse components from databanks", true);
        if (browse) {
            browseAndSelectComponents(s);
        }
        if (s.components.empty()) {
            const std::string comps = promptLine("components (csv)", "");
            if (!comps.empty()) {
                s.components = splitCsv(comps);
            }
        }
    }
    if (s.components.empty()) {
        shellPrintln(s, "components are required.");
        return false;
    }

    if (s.z.empty() || s.z.size() != s.components.size()) {
        const std::string defz = (s.z.empty()) ? "" : joinCsvDoubles(s.z);
        const std::string zcsv = promptLine("z (csv, same length as components)", defz);
        if (!zcsv.empty()) {
            s.z = normalizeComposition(parseCsvDoubles(zcsv));
        }
    }
    if (s.z.empty() || s.z.size() != s.components.size()) {
        shellPrintln(s, "z must be provided and match components length.");
        return false;
    }
    return true;
}

int runCommand(int argc, char** argv); // forward

int dispatchShellCommand(const std::vector<std::string>& args) {
    std::vector<std::string> av;
    av.reserve(args.size() + 1);
    av.push_back("dmtherm");
    for (const auto& a : args) av.push_back(a);

    std::vector<char*> argv;
    argv.reserve(av.size());
    for (auto& s : av) argv.push_back(const_cast<char*>(s.c_str()));
    return runCommand(static_cast<int>(argv.size()), argv.data());
}

std::string quoteForDisplay(const std::string& s) {
    if (s.empty()) return "\"\"";

    const bool needs_quotes = (s.find_first_of(" \t\r\n\"") != std::string::npos);
    if (!needs_quotes) return s;

    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (const char c : s) {
        if (c == '"' || c == '\\') out.push_back('\\');
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}

std::string joinArgsForDisplay(const std::vector<std::string>& args) {
    std::ostringstream oss;
    for (size_t i = 0; i < args.size(); ++i) {
        if (i) oss << " ";
        oss << quoteForDisplay(args[i]);
    }
    return oss.str();
}

struct ScopedStdCapture {
    std::ostringstream capture;
    std::streambuf* old_cout = nullptr;
    std::streambuf* old_cerr = nullptr;

    ScopedStdCapture() {
        old_cout = std::cout.rdbuf(capture.rdbuf());
        old_cerr = std::cerr.rdbuf(capture.rdbuf());
    }

    ~ScopedStdCapture() {
        std::cout.flush();
        std::cerr.flush();
        std::cout.rdbuf(old_cout);
        std::cerr.rdbuf(old_cerr);
    }
};

int dispatchShellCommandToUser(ShellSession& s, const std::vector<std::string>& args) {
    if (!shellUsesStickyMenu(s)) {
        return dispatchShellCommand(args);
    }

    uiLogAppendLine(s, "$ dmtherm " + joinArgsForDisplay(args));
    ScopedStdCapture cap;
    const int rc = dispatchShellCommand(args);
    uiLogAppendText(s, cap.capture.str());
    if (rc != 0) uiLogAppendLine(s, "(exit code " + std::to_string(rc) + ")");
    return rc;
}

void handleShellSet(ShellSession& s, const std::vector<std::string>& tokens) {
    if (tokens.size() < 3) {
        shellPrintln(s, "Usage: set <key> <value>");
        return;
    }
    const std::string key = toLower(tokens[1]);
    std::string value;
    for (size_t i = 2; i < tokens.size(); ++i) {
        if (i != 2) value += " ";
        value += tokens[i];
    }
    value = trimCopy(value);

    if (key == "databanks-dir") {
        s.databanks_dir = value;
        s.databanks_loaded = false;
        s.databanks_loaded_dir.clear();
        return;
    }
    if (key == "eos") {
        s.eos_type = value;
        return;
    }
    if (key == "flash-method") {
        s.flash_method = value;
        return;
    }
    if (key == "property-method") {
        // validate early
        (void)parsePropertyMethod(value);
        s.property_method = value;
        return;
    }
    if (key == "p-unit" || key == "pressure-unit") {
        s.pressure_unit = pressureUnitToken(parsePressureUnit(value));
        return;
    }
    if (key == "require-critical") {
        bool v = s.require_critical;
        if (!parseBool(value, v)) {
            shellPrintln(s, "Invalid boolean. Use true/false.");
            return;
        }
        s.require_critical = v;
        return;
    }
    if (key == "require-ideal-gas-cp") {
        bool v = s.require_ideal_gas_cp;
        if (!parseBool(value, v)) {
            shellPrintln(s, "Invalid boolean. Use true/false.");
            return;
        }
        s.require_ideal_gas_cp = v;
        return;
    }
    if (key == "always-show-menu") {
        bool v = s.always_show_menu;
        if (!parseBool(value, v)) {
            shellPrintln(s, "Invalid boolean. Use true/false.");
            return;
        }
        s.always_show_menu = v;
        return;
    }
    if (key == "sticky-menu") {
        bool v = s.sticky_menu;
        if (!parseBool(value, v)) {
            shellPrintln(s, "Invalid boolean. Use true/false.");
            return;
        }
        s.sticky_menu = v;
        return;
    }
    if (key == "ui-log-limit") {
        int v = s.ui_log_limit;
        try {
            v = std::stoi(value);
        } catch (...) {
            shellPrintln(s, "Invalid integer.");
            return;
        }
        if (v <= 0) {
            shellPrintln(s, "ui-log-limit must be > 0.");
            return;
        }
        s.ui_log_limit = v;
        while (static_cast<int>(s.ui_log.size()) > s.ui_log_limit) {
            s.ui_log.pop_front();
        }
        return;
    }
    if (key == "ui-tail-lines") {
        int v = s.ui_tail_lines;
        try {
            v = std::stoi(value);
        } catch (...) {
            shellPrintln(s, "Invalid integer.");
            return;
        }
        if (v < 0) {
            shellPrintln(s, "ui-tail-lines must be >= 0.");
            return;
        }
        s.ui_tail_lines = v;
        return;
    }
    if (key == "components") {
        s.components = splitCsv(value);
        if (!s.z.empty() && s.z.size() != s.components.size()) {
            s.z.clear();
        }
        return;
    }
    if (key == "z") {
        s.z = normalizeComposition(parseCsvDoubles(value));
        return;
    }

    shellPrintln(s, "Unknown key: " + key);
}

int runShell() {
    ShellSession s{};

    bool was_sticky = shellUsesStickyMenu(s);
    if (was_sticky) {
        uiLogAppendLine(s, "DMTherm interactive shell");
        uiLogAppendLine(s, "Type 'menu' to see options, 'help' for shell commands, 'exit' to quit.");
        uiLogAppendLine(s, "");
        std::ostringstream ss;
        printShellSession(ss, s);
        uiLogAppendText(s, ss.str());
    } else {
        std::cout << "DMTherm interactive shell\n";
        std::cout << "Type 'menu' to see options, 'help' for shell commands, 'exit' to quit.\n\n";
        printShellSession(std::cout, s);
        std::cout << "\n";
        if (!s.always_show_menu) {
            printShellMenu(std::cout, s);
        }
    }

    for (;;) {
        const bool sticky = shellUsesStickyMenu(s);
        if (sticky != was_sticky) {
            if (sticky) {
                uiLogAppendLine(s, "Sticky UI enabled.");
            } else {
                std::cout << "\nSticky UI disabled.\n";
            }
            was_sticky = sticky;
        }
        if (sticky) {
            renderStickyShellUI(s);
        } else if (s.always_show_menu) {
            std::cout << "\n";
            printShellMenu(std::cout, s);
        }
        if (sticky) {
            std::cout << "dmtherm> ";
        } else {
            std::cout << "\ndmtherm> ";
        }
        std::string line;
        if (!std::getline(std::cin, line)) {
            if (!sticky) std::cout << "\n";
            return 0;
        }
        line = trimCopy(line);
        if (line.empty()) continue;

        auto tokens = tokenizeCommandLine(line);
        if (tokens.empty()) continue;

        const std::string cmd = toLower(tokens[0]);

        if (cmd == "exit" || cmd == "quit") return 0;
        if (cmd == "clear" || cmd == "cls") {
            if (sticky) {
                uiLogClear(s);
            } else {
                clearScreen();
            }
            continue;
        }
        if (cmd == "help" || cmd == "?") {
            if (sticky) {
                std::ostringstream oss;
                printShellHelp(oss);
                uiLogAppendText(s, oss.str());
            } else {
                printShellHelp(std::cout);
            }
            continue;
        }
        if (cmd == "menu") {
            if (!sticky) {
                printShellMenu(std::cout, s);
            }
            continue;
        }
        if (cmd == "show") {
            if (sticky) {
                std::ostringstream oss;
                printShellSession(oss, s);
                uiLogAppendText(s, oss.str());
            } else {
                printShellSession(std::cout, s);
            }
            continue;
        }
        if (cmd == "set") {
            try {
                handleShellSet(s, tokens);
            } catch (const std::exception& e) {
                shellPrintln(s, std::string("Error: ") + e.what());
            }
            continue;
        }

        // Menu selection by number.
        bool parsed_num = false;
        int choice = -1;
        try {
            size_t pos = 0;
            choice = std::stoi(tokens[0], &pos);
            parsed_num = (pos == tokens[0].size());
        } catch (...) {
            parsed_num = false;
        }
        if (parsed_num) {
            try {
                if (choice == 0) return 0;
                if (choice == 1) {
                    if (sticky) {
                        std::ostringstream oss;
                        printShellSession(oss, s);
                        uiLogAppendText(s, oss.str());
                    } else {
                        printShellSession(std::cout, s);
                    }
                } else if (choice == 2) {
                    s.databanks_dir = promptLine("databanks-dir", s.databanks_dir);
                    s.databanks_loaded = false;
                    s.databanks_loaded_dir.clear();
                    if (sticky) uiLogAppendLine(s, "databanks-dir: " + (s.databanks_dir.empty() ? "<empty>" : s.databanks_dir));
                } else if (choice == 3) {
                    s.eos_type = promptEOSFromList(s.eos_type);
                    if (sticky) uiLogAppendLine(s, "eos: " + (s.eos_type.empty() ? "<empty>" : s.eos_type));
                } else if (choice == 4) {
                    browseAndSelectComponents(s);
                    if (sticky && !s.components.empty()) {
                        uiLogAppendLine(s, "components (" + std::to_string(s.components.size()) + "): " + joinCsv(s.components));
                        if (!s.z.empty()) uiLogAppendLine(s, "z: " + joinCsvDoubles(s.z));
                    }
                } else if (choice == 5) {
                    if (s.components.empty()) {
                        shellPrintln(s, "Set components first.");
                    } else {
                        const std::string zcsv = promptLine("z (csv, same length as components)", joinCsvDoubles(s.z));
                        s.z = normalizeComposition(parseCsvDoubles(zcsv));
                        if (sticky && !s.z.empty()) uiLogAppendLine(s, "z: " + joinCsvDoubles(s.z));
                    }
                } else if (choice == 6) {
                    if (!ensureSessionHasMixture(s)) continue;
                    const double T = promptDouble("T [K]", 298.15);
                    const double P = promptDouble(pressurePromptLabel("P", s), pressureDefaultInSessionUnit(101325.0, s));
                    const std::string root = promptLine("root (stable|vapor|liquid)", "stable");
                    auto args = std::vector<std::string>{"eval", "tp-phase", "--T", std::to_string(T), "--P", std::to_string(P), "--root", root};
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 7) {
                    if (!ensureSessionHasMixture(s)) continue;
                    const double T = promptDouble("T [K]", 298.15);
                    const double P = promptDouble(pressurePromptLabel("P", s), pressureDefaultInSessionUnit(101325.0, s));
                    auto args = std::vector<std::string>{"eval", "tp-roots", "--T", std::to_string(T), "--P", std::to_string(P)};
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 8) {
                    if (!ensureSessionHasMixture(s)) continue;
                    const double T = promptDouble("T [K]", 298.15);
                    const double P = promptDouble(pressurePromptLabel("P", s), pressureDefaultInSessionUnit(101325.0, s));
                    const int max_phases = promptInt("max-phases", 2);
                    const bool props = promptYesNo("compute properties", true);
                    const bool json = promptYesNo("json output", false);
                    std::vector<std::string> args{"flash", "pt", "--T", std::to_string(T), "--P", std::to_string(P), "--max-phases", std::to_string(max_phases)};
                    if (props) args.push_back("--properties");
                    if (json) args.push_back("--json");
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 9) {
                    if (!ensureSessionHasMixture(s)) continue;
                    const double T = promptDouble("T [K]", 298.15);
                    const double V = promptDouble("V [m^3/mol]", 1e-3);
                    const int max_phases = promptInt("max-phases", 2);
                    const bool props = promptYesNo("compute properties", true);
                    const bool json = promptYesNo("json output", false);
                    std::vector<std::string> args{"flash", "tv", "--T", std::to_string(T), "--V", std::to_string(V), "--max-phases", std::to_string(max_phases)};
                    if (props) args.push_back("--properties");
                    if (json) args.push_back("--json");
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 10) {
                    if (!ensureSessionHasMixture(s)) continue;
                    const double P = promptDouble(pressurePromptLabel("P", s), pressureDefaultInSessionUnit(101325.0, s));
                    const double H = promptDouble("H [J/mol]", 0.0);
                    const int max_phases = promptInt("max-phases", 2);
                    const bool props = promptYesNo("compute properties", true);
                    const bool json = promptYesNo("json output", false);
                    std::vector<std::string> args{"flash", "ph", "--P", std::to_string(P), "--H", std::to_string(H), "--max-phases", std::to_string(max_phases)};
                    if (props) args.push_back("--properties");
                    if (json) args.push_back("--json");
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 11) {
                    if (!ensureSessionHasMixture(s)) continue;
                    const double P = promptDouble(pressurePromptLabel("P", s), pressureDefaultInSessionUnit(101325.0, s));
                    const double Sval = promptDouble("S [J/mol/K]", 0.0);
                    const int max_phases = promptInt("max-phases", 2);
                    const bool props = promptYesNo("compute properties", true);
                    const bool json = promptYesNo("json output", false);
                    std::vector<std::string> args{"flash", "ps", "--P", std::to_string(P), "--S", std::to_string(Sval), "--max-phases", std::to_string(max_phases)};
                    if (props) args.push_back("--properties");
                    if (json) args.push_back("--json");
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 12) {
                    if (!ensureSessionHasMixture(s)) continue;
                    const double T_min = promptDouble("T-min [K]", 300.0);
                    const double T_max = promptDouble("T-max [K]", 310.0);
                    const int nT = promptInt("nT", 2);
                    const double P_min = promptDouble(pressurePromptLabel("P-min", s), pressureDefaultInSessionUnit(1e5, s));
                    const double P_max = promptDouble(pressurePromptLabel("P-max", s), pressureDefaultInSessionUnit(2e5, s));
                    const int nP = promptInt("nP", 2);
                    const bool p_log = promptYesNo("P-log", false);
                    const std::string root = promptLine("root (stable|vapor|liquid)", "stable");
                    std::vector<std::string> args{
                        "sweep", "tp-grid",
                        "--T-min", std::to_string(T_min),
                        "--T-max", std::to_string(T_max),
                        "--nT", std::to_string(nT),
                        "--P-min", std::to_string(P_min),
                        "--P-max", std::to_string(P_max),
                        "--nP", std::to_string(nP),
                        "--root", root
                    };
                    if (p_log) args.push_back("--P-log");
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 13) {
                    if (!ensureSessionHasMixture(s)) continue;
                    const double T = promptDouble("T [K]", 300.0);
                    const double P_min = promptDouble(pressurePromptLabel("P-min", s), pressureDefaultInSessionUnit(1e5, s));
                    const double P_max = promptDouble(pressurePromptLabel("P-max", s), pressureDefaultInSessionUnit(2e5, s));
                    const int nP = promptInt("nP", 10);
                    const bool p_log = promptYesNo("P-log", false);
                    const std::string root = promptLine("root (stable|vapor|liquid)", "stable");
                    std::vector<std::string> args{
                        "sweep", "isotherm",
                        "--T", std::to_string(T),
                        "--P-min", std::to_string(P_min),
                        "--P-max", std::to_string(P_max),
                        "--nP", std::to_string(nP),
                        "--root", root
                    };
                    if (p_log) args.push_back("--P-log");
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 14) {
                    if (!ensureSessionHasMixture(s)) continue;
                    const double P = promptDouble(pressurePromptLabel("P", s), pressureDefaultInSessionUnit(1e5, s));
                    const double T_min = promptDouble("T-min [K]", 300.0);
                    const double T_max = promptDouble("T-max [K]", 310.0);
                    const int nT = promptInt("nT", 10);
                    const std::string root = promptLine("root (stable|vapor|liquid)", "stable");
                    std::vector<std::string> args{
                        "sweep", "isobar",
                        "--P", std::to_string(P),
                        "--T-min", std::to_string(T_min),
                        "--T-max", std::to_string(T_max),
                        "--nT", std::to_string(nT),
                        "--root", root
                    };
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 15) {
                    const std::string raw = promptLine("command (e.g. eval tp-phase --T 298.15 --P 101325)", "");
                    auto toks = tokenizeCommandLine(raw);
                    toks = injectSessionDefaults(std::move(toks), s);
                    (void)dispatchShellCommandToUser(s, toks);
                } else if (choice == 16) {
                    std::vector<std::string> args{"list", "models"};
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 17) {
                    const std::string filter = promptLine("filter (blank=all)", "");
                    const int limit = promptInt("limit", 50);
                    const bool csv = promptYesNo("csv output", false);
                    std::vector<std::string> args{"list", "components", "--databanks-dir", s.databanks_dir};
                    if (!filter.empty()) {
                        args.push_back("--filter");
                        args.push_back(filter);
                    }
                    args.push_back("--limit");
                    args.push_back(std::to_string(limit));
                    if (csv) args.push_back("--csv");
                    (void)dispatchShellCommandToUser(s, args);
                } else if (choice == 18) {
                    if (!ensureSessionHasComponents(s)) continue;
                    const std::string models = promptModelsArgFromList();
                    const bool csv = promptYesNo("csv output", false);
                    std::vector<std::string> args{"eval", "critical"};
                    if (!models.empty()) {
                        args.push_back("--models");
                        args.push_back(models);
                    }
                    if (csv) args.push_back("--csv");
                    args = injectSessionDefaults(std::move(args), s);
                    (void)dispatchShellCommandToUser(s, args);
                } else {
                    shellPrintln(s, "Unknown menu selection.");
                }
            } catch (const std::exception& e) {
                shellPrintln(s, std::string("Error: ") + e.what());
            }
            continue;
        }

        // "run <cmd...>" is a convenience wrapper.
        if (cmd == "run") {
            tokens.erase(tokens.begin());
        }

        // Treat as a dmtherm command line (inject session defaults).
        try {
            auto args = injectSessionDefaults(std::move(tokens), s);
            const int rc = dispatchShellCommandToUser(s, args);
            if (!sticky && rc != 0) std::cout << "(exit code " << rc << ")\n";
        } catch (const std::exception& e) {
            shellPrintln(s, std::string("Error: ") + e.what());
        }
    }
}

void printTPPhaseHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm eval tp-phase --T <K> --P <value> --components <csv> --z <csv> [options]\n\n";
    os << "Options:\n";
    os << "  --root vapor|liquid|stable   (default: stable)\n";
    os << "  --P-unit <unit>              (default: PA_ABS)\n";
    os << "  --csv                        (print a single CSV row)\n";
}

void printTPRootsHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm eval tp-roots --T <K> --P <value> --components <csv> --z <csv> [options]\n\n";
    os << "Options:\n";
    os << "  --P-unit <unit>              (default: PA_ABS)\n";
    os << "  --csv                        (print CSV rows; one per root)\n";
}

void printEvalCriticalHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm eval critical --components <csv> [options]\n\n";
    os << "Options:\n";
    os << "  --databanks-dir <dir>        (default: data)\n";
    os << "  --models <csv|all>           (default: all)\n";
    os << "  --csv                        (print CSV rows; one per (component, model))\n";
}

void printFlashPTHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm flash pt --T <K> --P <value> --components <csv> --z <csv> [options]\n\n";
    os << "Options:\n";
    os << "  --P-unit <unit>              (default: PA_ABS)\n";
    os << "  --max-phases <N>             (default: 2)\n";
    os << "  --properties                 (populate per-phase/overall H/S/G when available)\n";
    os << "  --json                       (write JSON instead of text)\n";
}

void printFlashTVHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm flash tv --T <K> --V <m^3/mol> --components <csv> --z <csv> [options]\n\n";
    os << "Options:\n";
    os << "  --max-phases <N>             (default: 2)\n";
    os << "  --properties                 (populate per-phase/overall H/S/G when available)\n";
    os << "  --json                       (write JSON instead of text)\n";
}

void printFlashPHHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm flash ph --P <value> --H <J/mol> --components <csv> --z <csv> [options]\n\n";
    os << "Options:\n";
    os << "  --P-unit <unit>              (default: PA_ABS)\n";
    os << "  --max-phases <N>             (default: 2)\n";
    os << "  --properties                 (populate per-phase/overall H/S/G when available)\n";
    os << "  --enthalpy-tol <J/mol>       (default: engine default)\n";
    os << "  --temp-low <K>               (default: engine default)\n";
    os << "  --temp-high <K>              (default: engine default)\n";
    os << "  --max-iterations <N>         (default: engine default)\n";
    os << "  --json                       (write JSON instead of text)\n";
}

void printFlashPSHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm flash ps --P <value> --S <J/mol/K> --components <csv> --z <csv> [options]\n\n";
    os << "Options:\n";
    os << "  --P-unit <unit>              (default: PA_ABS)\n";
    os << "  --max-phases <N>             (default: 2)\n";
    os << "  --properties                 (populate per-phase/overall H/S/G when available)\n";
    os << "  --entropy-tol <J/mol/K>      (default: engine default)\n";
    os << "  --temp-low <K>               (default: engine default)\n";
    os << "  --temp-high <K>              (default: engine default)\n";
    os << "  --max-iterations <N>         (default: engine default)\n";
    os << "  --json                       (write JSON instead of text)\n";
}

void printSweepTPGridHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm sweep tp-grid --components <csv> --z <csv> --T-min <K> --T-max <K> --nT <N> --P-min <value> --P-max <value> --nP <N> [options]\n\n";
    os << "Options:\n";
    os << "  --root vapor|liquid|stable   (default: stable)\n";
    os << "  --P-unit <unit>              (default: PA_ABS)\n";
    os << "  --P-log                      (log-space pressures between P-min and P-max)\n";
    os << "  --json                       (write JSON instead of CSV)\n";
    os << "  --out <path>                 (write to file instead of stdout)\n";
    os << "  --no-header                  (omit CSV header)\n";
}

void printSweepIsothermHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm sweep isotherm --T <K> --components <csv> --z <csv> --P-min <value> --P-max <value> --nP <N> [options]\n\n";
    os << "Options:\n";
    os << "  --root vapor|liquid|stable   (default: stable)\n";
    os << "  --P-unit <unit>              (default: PA_ABS)\n";
    os << "  --P-log                      (log-space pressures between P-min and P-max)\n";
    os << "  --json                       (write JSON instead of CSV)\n";
    os << "  --out <path>                 (write to file instead of stdout)\n";
    os << "  --no-header                  (omit CSV header)\n";
}

void printSweepIsobarHelp(std::ostream& os) {
    os << "Usage:\n";
    os << "  dmtherm sweep isobar --P <value> --components <csv> --z <csv> --T-min <K> --T-max <K> --nT <N> [options]\n\n";
    os << "Options:\n";
    os << "  --root vapor|liquid|stable   (default: stable)\n";
    os << "  --P-unit <unit>              (default: PA_ABS)\n";
    os << "  --json                       (write JSON instead of CSV)\n";
    os << "  --out <path>                 (write to file instead of stdout)\n";
    os << "  --no-header                  (omit CSV header)\n";
}

int cmdListComponents(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printListComponentsHelp(std::cout);
        return 0;
    }

    const std::string databanks_dir = getStr(o, "databanks-dir", "data");
    const std::string ql = toLower(getStr(o, "filter", ""));
    const int limit = std::max(0, getInt(o, "limit", 50));
    const bool csv = hasFlag(o, "csv");

    if (!std::filesystem::exists(databanks_dir)) {
        std::cerr << "databanks dir does not exist: " << databanks_dir << "\n";
        return 2;
    }

    DMThermo::Data::Databanks db;
    DMThermo::Data::Diagnostics diag;
    const bool ok = db.loadAllFromDirectory(databanks_dir, &diag);
    if (!ok || diag.hasErrors()) {
        std::cerr << "Failed to load databanks from: " << databanks_dir << "\n";
        for (const auto& d : diag.items()) {
            if (d.level == DMThermo::Data::DiagnosticLevel::Error) {
                std::cerr << "  error: " << d.message << "\n";
            }
        }
        return 2;
    }

    std::vector<const DMThermo::Data::DipprRecord*> matches;
    matches.reserve(256);
    for (const auto& r : db.dippr.records()) {
        if (!recordMatchesSearch(r, ql)) continue;
        matches.push_back(&r);
    }

    std::sort(matches.begin(), matches.end(), [](const auto* a, const auto* b) {
        const std::string an = toLower(a->key.name);
        const std::string bn = toLower(b->key.name);
        if (an != bn) return an < bn;
        return recordIdentifier(*a) < recordIdentifier(*b);
    });

    const size_t shown = std::min(matches.size(), static_cast<size_t>(limit));
    if (csv) {
        std::cout << "identifier,name,cas,chemid,molecular_weight_g_per_mol,critical_temperature_k,critical_pressure_pa_abs,acentric_factor,critical_volume_m3_per_kmol,critical_compressibility,has_ideal_gas_cp,has_critical\n";
        for (size_t i = 0; i < shown; ++i) {
            const auto& r = *matches[i];
            const std::string id = recordIdentifier(r);
            const bool has_critical = r.Tc.has_value() && r.Pc.has_value() && r.omega.has_value();
            std::cout << csvEscape(id) << ","
                      << csvEscape(r.key.name) << ","
                      << csvEscape(r.key.cas) << ",";
            if (r.key.chemid.has_value()) std::cout << r.key.chemid.value();
            std::cout << ",";
            if (r.MW.has_value()) std::cout << r.MW.value();
            std::cout << ",";
            if (r.Tc.has_value()) std::cout << r.Tc.value();
            std::cout << ",";
            if (r.Pc.has_value()) std::cout << r.Pc.value();
            std::cout << ",";
            if (r.omega.has_value()) std::cout << r.omega.value();
            std::cout << ",";
            if (r.Vc.has_value()) std::cout << r.Vc.value();
            std::cout << ",";
            if (r.Zc.has_value()) std::cout << r.Zc.value();
            std::cout << ","
                      << (r.ideal_gas_cp.has_value() ? "1" : "0") << ","
                      << (has_critical ? "1" : "0") << "\n";
        }
    } else {
        std::cout << "Found " << matches.size() << " components";
        if (!ql.empty()) std::cout << " matching filter '" << ql << "'";
        std::cout << " (showing " << shown << ")\n";
        for (size_t i = 0; i < shown; ++i) {
            const auto& r = *matches[i];
            std::cout << "  [" << (i + 1) << "] " << r.key.name;
            const std::string id = recordIdentifier(r);
            if (!id.empty()) std::cout << "  id=" << id;
            if (!r.key.cas.empty()) std::cout << "  CAS=" << r.key.cas;
            if (r.key.chemid.has_value()) std::cout << "  CHEMID=" << r.key.chemid.value();
            std::cout << "\n";
        }
    }

    return 0;
}

int cmdListModels(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printListModelsHelp(std::cout);
        return 0;
    }

    const bool csv = hasFlag(o, "csv");

    const auto types = DMThermo::Factory::EOSFactory::registeredTypes();
    if (csv) {
        std::cout << "eos_type,kind,notes\n";
        for (const auto& t : types) {
            std::string kind;
            std::string notes;
            if (t == "PC-SAFT") {
                kind = "pcsaft";
                notes = "requires PC-SAFT parameters (pcsaft.csv)";
            } else if (t == "SRK" || t == "Peng-Robinson" || t == "PR-VT" || t == "VTPR") {
                kind = "cubic";
                if (t == "PR-VT") {
                    notes = "Peneloux-style volume-translated Peng-Robinson (uses DIPPR Zc when available)";
                } else if (t == "VTPR") {
                    notes = "PR + Wong-Sandler mixing (UNIFAC) + volume translation (vtpr_*.csv + unifac_*.csv)";
                }
            } else {
                kind = "unknown";
            }

            std::cout << csvEscape(t) << "," << csvEscape(kind) << "," << csvEscape(notes) << "\n";
        }
    } else {
        std::cout << "Available EOS models (use with --eos or in shell via 'set eos ...'):\n";
        for (const auto& t : types) {
            std::cout << "  - " << t;
            if (t == "PC-SAFT") std::cout << " (requires PC-SAFT parameters in pcsaft.csv)";
            std::cout << "\n";
        }
    }

    return 0;
}

struct SystemBuildInputs {
    std::string databanks_dir;
    std::string eos_type;
    std::string flash_method;
    dmtherm_property_method_t property_method;
    bool require_critical_properties;
    bool require_ideal_gas_cp;
    std::vector<std::string> components;
};

SystemBuildInputs parseSystemBuildInputs(const Options& o) {
    SystemBuildInputs in{};
    in.databanks_dir = getStr(o, "databanks-dir", "data");
    in.eos_type = getStr(o, "eos", "Peng-Robinson");
    in.flash_method = getStr(o, "flash-method", "Auto");
    in.property_method = parsePropertyMethod(getStr(o, "property-method", "eos_only"));
    in.require_critical_properties = hasFlag(o, "require-critical");
    in.require_ideal_gas_cp = !hasFlag(o, "no-require-ideal-gas-cp");

    const std::string components_csv = getStr(o, "components", "");
    if (components_csv.empty()) {
        throw std::invalid_argument("Missing required option: --components");
    }
    in.components = splitCsv(components_csv);
    if (in.components.empty()) {
        throw std::invalid_argument("--components must not be empty");
    }
    return in;
}

std::vector<double> parseComposition(const Options& o, const std::vector<std::string>& components) {
    const std::string z_csv = getStr(o, "z", "");
    if (z_csv.empty()) throw std::invalid_argument("Missing required option: --z");
    auto z = parseCsvDoubles(z_csv);
    if (z.size() != components.size()) {
        std::ostringstream msg;
        msg << "--z length (" << z.size() << ") must match --components length (" << components.size() << ")";
        throw std::invalid_argument(msg.str());
    }
    return normalizeComposition(std::move(z));
}

dmtherm_system_t* buildSystem(const SystemBuildInputs& in) {
    if (!std::filesystem::exists(in.databanks_dir)) {
        throw std::invalid_argument("databanks dir does not exist: " + in.databanks_dir);
    }

    std::vector<const char*> comp_ptrs;
    comp_ptrs.reserve(in.components.size());
    for (const auto& c : in.components) comp_ptrs.push_back(c.c_str());

    dmtherm_system_config_t cfg{};
    dmtherm_system_config_set_defaults(&cfg);
    cfg.databanks_dir = in.databanks_dir.c_str();
    cfg.eos_type = in.eos_type.c_str();
    cfg.flash_method = in.flash_method.c_str();
    cfg.components = comp_ptrs.data();
    cfg.num_components = comp_ptrs.size();
    cfg.require_critical_properties = in.require_critical_properties ? 1 : 0;
    cfg.require_ideal_gas_cp = in.require_ideal_gas_cp ? 1 : 0;
    cfg.property_method = in.property_method;

    dmtherm_system_t* sys = nullptr;
    char* err = nullptr;
    const auto st = dmtherm_system_create(&cfg, &sys, &err);
    if (st != DMTHERM_STATUS_OK) {
        std::string msg = err ? err : "";
        if (err) dmtherm_free(err);
        throw std::runtime_error("dmtherm_system_create failed (status=" + std::to_string(st) + "): " + msg);
    }
    if (err) dmtherm_free(err);
    return sys;
}

void printFlashResultJson(
    std::ostream& os,
    const char* mode_name,
    const SystemBuildInputs& in,
    const std::vector<double>& z,
    int max_phases,
    bool compute_properties,
    const char* a_key,
    double a,
    const char* b_key,
    double b,
    const dmtherm_flash_pt_result_t& out)
{
    os << std::setprecision(17);

    std::optional<double> v_mix;
    if (out.phases && out.phases_len != 0) {
        double v_acc = 0.0;
        double frac_acc = 0.0;
        for (size_t k = 0; k < out.phases_len; ++k) {
            const auto& ph = out.phases[k];
            if (std::isfinite(ph.fraction) && std::isfinite(ph.molar_density_mol_per_m3) && ph.molar_density_mol_per_m3 > 0.0) {
                v_acc += ph.fraction * (1.0 / ph.molar_density_mol_per_m3);
                frac_acc += ph.fraction;
            }
        }
        if (frac_acc > 0.0 && std::isfinite(v_acc)) {
            v_mix = v_acc / frac_acc;
        }
    }

    std::optional<double> U_mix;
    if (out.mixture_enthalpy_j_per_mol.has_value && v_mix.has_value() && std::isfinite(out.pressure_pa_abs)) {
        U_mix = out.mixture_enthalpy_j_per_mol.value - out.pressure_pa_abs * v_mix.value();
    }
    std::optional<double> A_mix;
    if (U_mix.has_value() && out.mixture_entropy_j_per_mol_k.has_value && std::isfinite(out.temperature_k)) {
        A_mix = U_mix.value() - out.temperature_k * out.mixture_entropy_j_per_mol_k.value;
    }

    os << "{";
    os << "\"schema\":\"dmtherm.flash." << jsonEscape(mode_name) << ".v1\",";
    os << "\"meta\":{";
    os << "\"databanks_dir\":\"" << jsonEscape(in.databanks_dir) << "\",";
    os << "\"eos\":\"" << jsonEscape(in.eos_type) << "\",";
    os << "\"flash_method\":\"" << jsonEscape(in.flash_method) << "\",";
    os << "\"property_method\":\"" << jsonEscape(propertyMethodToString(in.property_method)) << "\",";
    os << "\"require_critical\":" << (in.require_critical_properties ? "true" : "false") << ",";
    os << "\"require_ideal_gas_cp\":" << (in.require_ideal_gas_cp ? "true" : "false") << ",";
    os << "\"max_phases\":" << max_phases << ",";
    os << "\"compute_properties\":" << (compute_properties ? "true" : "false") << ",";
    os << "\"components\":";
    jsonWriteStringArray(os, in.components);
    os << ",";
    os << "\"z\":";
    jsonWriteDoubleArray(os, z.data(), z.size());
    os << "},";
    const char* a_canonical_key = canonicalFlashInputKey(a_key);
    const char* b_canonical_key = canonicalFlashInputKey(b_key);
    const std::string a_input_key = (a_canonical_key != nullptr) ? a_canonical_key : (a_key ? a_key : "value_a");
    const std::string b_input_key = (b_canonical_key != nullptr) ? b_canonical_key : (b_key ? b_key : "value_b");
    os << "\"input\":{";
    os << "\"" << jsonEscape(a_input_key) << "\":" << a << ",";
    os << "\"" << jsonEscape(b_input_key) << "\":" << b;
    os << "},";
    os << "\"result\":{";
    os << "\"converged\":" << (out.converged ? "true" : "false") << ",";
    os << "\"iterations\":" << out.iterations << ",";
    os << "\"residual\":" << out.residual << ",";
    os << "\"temperature_k\":" << out.temperature_k << ",";
    os << "\"pressure_pa_abs\":" << out.pressure_pa_abs << ",";
    os << "\"vapor_fraction\":" << out.vapor_fraction << ",";
    os << "\"phase_count\":" << out.num_phases << ",";
    os << "\"method_used\":";
    jsonWriteStringOrNull(os, out.method_used);
    os << ",";
    os << "\"message\":";
    jsonWriteStringOrNull(os, out.message);
    os << ",";
    os << "\"mixture_enthalpy_j_per_mol\":";
    jsonWriteOptionalDouble(os, out.mixture_enthalpy_j_per_mol);
    os << ",";
    os << "\"mixture_entropy_j_per_mol_k\":";
    jsonWriteOptionalDouble(os, out.mixture_entropy_j_per_mol_k);
    os << ",";
    os << "\"mixture_gibbs_j_per_mol\":";
    jsonWriteOptionalDouble(os, out.mixture_gibbs_j_per_mol);
    os << ",";
    os << "\"mixture_molar_volume_m3_per_mol\":";
    jsonWriteOptionalDouble(os, v_mix);
    os << ",";
    os << "\"mixture_internal_energy_j_per_mol\":";
    jsonWriteOptionalDouble(os, U_mix);
    os << ",";
    os << "\"mixture_helmholtz_j_per_mol\":";
    jsonWriteOptionalDouble(os, A_mix);
    os << ",";
    os << "\"phases\":[";
    for (size_t k = 0; k < out.phases_len; ++k) {
        if (k) os << ",";
        const auto& ph = out.phases[k];

        const double v = (std::isfinite(ph.molar_density_mol_per_m3) && ph.molar_density_mol_per_m3 > 0.0) ? (1.0 / ph.molar_density_mol_per_m3) : 0.0;
        std::optional<double> U;
        if (ph.enthalpy_j_per_mol.has_value && std::isfinite(ph.molar_density_mol_per_m3) && ph.molar_density_mol_per_m3 > 0.0) {
            U = ph.enthalpy_j_per_mol.value - out.pressure_pa_abs * v;
        }
        std::optional<double> A;
        if (U.has_value() && ph.entropy_j_per_mol_k.has_value) {
            A = U.value() - out.temperature_k * ph.entropy_j_per_mol_k.value;
        }

        os << "{";
        os << "\"phase\":\"" << jsonEscape(phaseToString(ph.phase)) << "\",";
        os << "\"fraction\":" << ph.fraction << ",";
        os << "\"molar_density_mol_per_m3\":" << ph.molar_density_mol_per_m3 << ",";
        os << "\"molar_volume_m3_per_mol\":" << v << ",";
        os << "\"compressibility\":" << ph.compressibility << ",";
        os << "\"composition\":";
        if (ph.composition && ph.composition_len) jsonWriteDoubleArray(os, ph.composition, ph.composition_len);
        else os << "[]";
        os << ",";
        os << "\"fugacity_coeff\":";
        if (ph.fugacity_coeff && ph.fugacity_coeff_len) jsonWriteDoubleArray(os, ph.fugacity_coeff, ph.fugacity_coeff_len);
        else os << "[]";
        os << ",";
        os << "\"enthalpy_departure_j_per_mol\":";
        jsonWriteOptionalDouble(os, ph.enthalpy_departure_j_per_mol);
        os << ",";
        os << "\"entropy_departure_j_per_mol_k\":";
        jsonWriteOptionalDouble(os, ph.entropy_departure_j_per_mol_k);
        os << ",";
        os << "\"gibbs_departure_j_per_mol\":";
        jsonWriteOptionalDouble(os, ph.gibbs_departure_j_per_mol);
        os << ",";
        os << "\"enthalpy_j_per_mol\":";
        jsonWriteOptionalDouble(os, ph.enthalpy_j_per_mol);
        os << ",";
        os << "\"entropy_j_per_mol_k\":";
        jsonWriteOptionalDouble(os, ph.entropy_j_per_mol_k);
        os << ",";
        os << "\"gibbs_j_per_mol\":";
        jsonWriteOptionalDouble(os, ph.gibbs_j_per_mol);
        os << ",";
        os << "\"internal_energy_j_per_mol\":";
        jsonWriteOptionalDouble(os, U);
        os << ",";
        os << "\"helmholtz_j_per_mol\":";
        jsonWriteOptionalDouble(os, A);
        os << "}";
    }
    os << "]";
    os << "}";
    os << "}";
    os << "\n";
}

void printFlashResult(
    const dmtherm_flash_pt_result_t& out,
    const std::vector<std::string>& components,
    const std::vector<double>& z,
    const std::vector<double>& component_mws)
{
    std::cout << std::setprecision(17);
    std::cout << "temperature_k=" << out.temperature_k << ", pressure_pa_abs=" << out.pressure_pa_abs << "\n";
    std::cout << "converged=" << (out.converged ? "true" : "false")
              << ", iterations=" << out.iterations
              << ", residual=" << out.residual
              << ", vapor_fraction=" << out.vapor_fraction
              << ", phase_count=" << out.num_phases << "\n";
    if (out.method_used && out.method_used[0] != '\0') {
        std::cout << "method: " << out.method_used << "\n";
    }
    if (out.message && out.message[0] != '\0') {
        std::cout << "message: " << out.message << "\n";
    }

    std::optional<double> mw_mix;
    if (!z.empty() && z.size() == component_mws.size()) {
        mw_mix = averageMW(z, component_mws);
        if (mw_mix.has_value()) {
            std::cout << "mixture_molecular_weight_g_per_mol=" << mw_mix.value() << "\n";
        }
    }

    if (auto v = optDouble(out.mixture_enthalpy_j_per_mol)) std::cout << "mixture_enthalpy_j_per_mol=" << v.value() << "\n";
    if (auto v = optDouble(out.mixture_entropy_j_per_mol_k)) std::cout << "mixture_entropy_j_per_mol_k=" << v.value() << "\n";
    if (auto v = optDouble(out.mixture_gibbs_j_per_mol)) std::cout << "mixture_gibbs_j_per_mol=" << v.value() << "\n";

    std::optional<double> v_mix;
    {
        double v_acc = 0.0;
        double frac_acc = 0.0;
        for (size_t k = 0; k < out.phases_len; ++k) {
            const auto& ph = out.phases[k];
            if (std::isfinite(ph.fraction) && std::isfinite(ph.molar_density_mol_per_m3) && ph.molar_density_mol_per_m3 > 0.0) {
                v_acc += ph.fraction * (1.0 / ph.molar_density_mol_per_m3);
                frac_acc += ph.fraction;
            }
        }
        if (frac_acc > 0.0 && std::isfinite(v_acc)) {
            v_mix = v_acc / frac_acc;
        }
    }
    if (v_mix.has_value()) std::cout << "mixture_molar_volume_m3_per_mol=" << v_mix.value() << "\n";

    std::optional<double> U_mix;
    if (auto H = optDouble(out.mixture_enthalpy_j_per_mol); H && v_mix.has_value()) {
        U_mix = H.value() - out.pressure_pa_abs * v_mix.value();
        std::cout << "mixture_internal_energy_j_per_mol=" << U_mix.value() << "\n";
    }
    if (auto S = optDouble(out.mixture_entropy_j_per_mol_k); S && U_mix.has_value()) {
        const double A_mix = U_mix.value() - out.temperature_k * S.value();
        std::cout << "mixture_helmholtz_j_per_mol=" << A_mix << "\n";
    }

    if (mw_mix.has_value()) {
        if (auto v = optDouble(out.mixture_enthalpy_j_per_mol)) {
            if (auto vm = molarToMass(v.value(), mw_mix.value())) std::cout << "mixture_enthalpy_j_per_kg=" << vm.value() << "\n";
        }
        if (auto v = optDouble(out.mixture_entropy_j_per_mol_k)) {
            if (auto vm = molarToMass(v.value(), mw_mix.value())) std::cout << "mixture_entropy_j_per_kg_k=" << vm.value() << "\n";
        }
    }

    for (size_t k = 0; k < out.phases_len; ++k) {
        const auto& ph = out.phases[k];
        std::cout << "\nphase[" << k << "] type=" << phaseToString(ph.phase)
                  << " fraction=" << ph.fraction
                  << " molar_density_mol_per_m3=" << ph.molar_density_mol_per_m3
                  << " compressibility=" << ph.compressibility << "\n";

        if (std::isfinite(ph.molar_density_mol_per_m3) && ph.molar_density_mol_per_m3 > 0.0) {
            std::cout << "  molar_volume_m3_per_mol=" << (1.0 / ph.molar_density_mol_per_m3) << "\n";
        }

        if (auto v = optDouble(ph.enthalpy_j_per_mol)) std::cout << "  enthalpy_j_per_mol=" << v.value() << "\n";
        if (auto v = optDouble(ph.entropy_j_per_mol_k)) std::cout << "  entropy_j_per_mol_k=" << v.value() << "\n";
        if (auto v = optDouble(ph.gibbs_j_per_mol)) std::cout << "  gibbs_j_per_mol=" << v.value() << "\n";

        if (std::isfinite(ph.molar_density_mol_per_m3) && ph.molar_density_mol_per_m3 > 0.0) {
            const double vv = 1.0 / ph.molar_density_mol_per_m3;
            if (auto H = optDouble(ph.enthalpy_j_per_mol)) {
                const double U = H.value() - out.pressure_pa_abs * vv;
                std::cout << "  internal_energy_j_per_mol=" << U << "\n";
                if (auto S = optDouble(ph.entropy_j_per_mol_k)) {
                    const double A = U - out.temperature_k * S.value();
                    std::cout << "  helmholtz_j_per_mol=" << A << "\n";
                }
            }
        }

        if (ph.composition && ph.composition_len == component_mws.size() && ph.composition_len == components.size()) {
            std::vector<double> x(ph.composition, ph.composition + ph.composition_len);
            const double mw = averageMW(x, component_mws);
            std::cout << "  molecular_weight_g_per_mol=" << mw << "\n";
            if (auto v = optDouble(ph.enthalpy_j_per_mol)) {
                if (auto vm = molarToMass(v.value(), mw)) std::cout << "  enthalpy_j_per_kg=" << vm.value() << "\n";
            }
            if (auto v = optDouble(ph.entropy_j_per_mol_k)) {
                if (auto vm = molarToMass(v.value(), mw)) std::cout << "  entropy_j_per_kg_k=" << vm.value() << "\n";
            }
        }

        if (auto v = optDouble(ph.enthalpy_departure_j_per_mol)) std::cout << "  enthalpy_departure_j_per_mol=" << v.value() << "\n";
        if (auto v = optDouble(ph.entropy_departure_j_per_mol_k)) std::cout << "  entropy_departure_j_per_mol_k=" << v.value() << "\n";
        if (auto v = optDouble(ph.gibbs_departure_j_per_mol)) std::cout << "  gibbs_departure_j_per_mol=" << v.value() << "\n";

        std::cout << "  composition =";
        for (size_t i = 0; i < ph.composition_len; ++i) std::cout << " " << ph.composition[i];
        std::cout << "\n";

        if (ph.fugacity_coeff && ph.fugacity_coeff_len) {
            std::cout << "  fugacity_coeff =";
            for (size_t i = 0; i < ph.fugacity_coeff_len; ++i) std::cout << " " << ph.fugacity_coeff[i];
            std::cout << "\n";
        }

        if (ph.composition && ph.composition_len == components.size()) {
            // no-op: composition already printed; placeholder for future per-component formatted output
        }
    }
}

using FlashFn = dmtherm_status_t (*)(
    dmtherm_system_t*,
    double,
    double,
    const double*,
    size_t,
    const dmtherm_flash_config_t*,
    dmtherm_flash_pt_result_t*);

int cmdFlashGeneric(
    int argc,
    char** argv,
    int start_index,
    const char* mode_name,
    FlashFn fn,
    void (*help_printer)(std::ostream&),
    const char* a_key,
    const char* b_key)
{
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        help_printer(std::cout);
        return 0;
    }

    const auto in = parseSystemBuildInputs(o);
    const auto a = getDouble(o, a_key);
    double b = getDouble(o, b_key);
    const auto z = parseComposition(o, in.components);
    const int max_phases = getInt(o, "max-phases", 2);
    const bool compute_properties = hasFlag(o, "properties");
    const bool json = hasFlag(o, "json");
    const bool b_is_pressure = (std::string(b_key) == "P");
    dmtherm_pressure_unit_t p_unit = DMTHERM_PRESSURE_UNIT_PA_ABS;
    if (b_is_pressure) {
        p_unit = parsePressureUnit(getStr(o, "P-unit", "PA_ABS"));
        b = convertPressureToPaAbs(b, p_unit);
    }

    dmtherm_system_t* sys = nullptr;
    try {
        sys = buildSystem(in);
        std::vector<double> component_mws;
        if (!json) {
            component_mws = getComponentMWs(sys, in.components.size());
        }

        dmtherm_flash_config_t fcfg{};
        dmtherm_flash_config_set_defaults(&fcfg);
        fcfg.max_phases = max_phases;
        fcfg.compute_properties = compute_properties ? 1 : 0;

        dmtherm_flash_pt_result_t out{};
        out.struct_size = static_cast<uint32_t>(sizeof(out));

        const auto st = fn(sys, a, b, z.data(), z.size(), &fcfg, &out);
        const std::string last_error = dmtherm_system_last_error(sys);
        if (st != DMTHERM_STATUS_OK) {
            dmtherm_flash_pt_result_destroy(&out);
            throw std::runtime_error(std::string("flash ") + mode_name + " failed (status=" + std::to_string(st) + "): " + last_error);
        }

        if (json) {
            printFlashResultJson(std::cout, mode_name, in, z, max_phases, compute_properties, a_key, a, b_key, b, out);
        } else {
            std::cout << "mode=" << mode_name << "  " << a_key << "=" << a << "  " << b_key << "=" << b << "\n";
            if (b_is_pressure) {
                std::cout << "pressure input unit=" << pressureUnitToken(p_unit) << ", solver unit=PA_ABS\n";
            }
            printFlashResult(out, in.components, z, component_mws);
        }

        const bool ok = (out.converged != 0);
        dmtherm_flash_pt_result_destroy(&out);
        dmtherm_system_destroy(sys);
        return ok ? 0 : 2;
    } catch (...) {
        if (sys) dmtherm_system_destroy(sys);
        throw;
    }
}

int cmdEvalTPPhase(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printTPPhaseHelp(std::cout);
        return 0;
    }

    const auto in = parseSystemBuildInputs(o);
    const auto T = getDouble(o, "T");
    const auto P_input = getDouble(o, "P");
    const auto p_unit = parsePressureUnit(getStr(o, "P-unit", "PA_ABS"));
    const auto P = convertPressureToPaAbs(P_input, p_unit);
    const auto z = parseComposition(o, in.components);
    const auto root = parseRootSelection(getStr(o, "root", "stable"));
    const bool csv = hasFlag(o, "csv");

    dmtherm_system_t* sys = nullptr;
    try {
        sys = buildSystem(in);

        dmtherm_tp_phase_eval_t out{};
        out.struct_size = static_cast<uint32_t>(sizeof(out));

        const auto st = dmtherm_system_evaluate_tp_phase(sys, T, P, z.data(), z.size(), root, &out);
        const std::string last_error = dmtherm_system_last_error(sys);
        if (st != DMTHERM_STATUS_OK) {
            dmtherm_tp_phase_eval_destroy(&out);
            throw std::runtime_error("evaluate_tp_phase failed (status=" + std::to_string(st) + "): " + last_error);
        }

        if (csv) {
            std::cout << std::setprecision(17);
            std::cout << "temperature_k,pressure_pa_abs,root,success,phase,molar_density_mol_per_m3,molar_volume_m3_per_mol,compressibility,enthalpy_departure_j_per_mol,entropy_departure_j_per_mol_k,gibbs_departure_j_per_mol,enthalpy_j_per_mol,entropy_j_per_mol_k,gibbs_j_per_mol,internal_energy_j_per_mol,helmholtz_j_per_mol,molecular_weight_g_per_mol,enthalpy_j_per_kg,entropy_j_per_kg_k,internal_energy_j_per_kg,helmholtz_j_per_kg,ln_fugacity_coeff,message\n";
            std::ostringstream lnphi;
            lnphi << std::setprecision(17);
            for (size_t i = 0; i < out.ln_fugacity_coeff_len; ++i) {
                if (i) lnphi << ",";
                lnphi << out.ln_fugacity_coeff[i];
            }

            std::optional<double> mw_mix;
            std::optional<double> H_mass;
            std::optional<double> S_mass;
            std::optional<double> U_mass;
            std::optional<double> A_mass;
            try {
                const auto component_mws = getComponentMWs(sys, in.components.size());
                mw_mix = averageMW(z, component_mws);
                H_mass = molarToMass(out.enthalpy_j_per_mol, mw_mix.value());
                S_mass = molarToMass(out.entropy_j_per_mol_k, mw_mix.value());
                if (std::isfinite(out.molar_density_mol_per_m3) && out.molar_density_mol_per_m3 > 0.0) {
                    const double v = 1.0 / out.molar_density_mol_per_m3;
                    const double U = out.enthalpy_j_per_mol - out.pressure_pa_abs * v;
                    const double A = U - out.temperature_k * out.entropy_j_per_mol_k;
                    U_mass = molarToMass(U, mw_mix.value());
                    A_mass = molarToMass(A, mw_mix.value());
                }
            } catch (...) {
                // Best-effort only for CLI convenience output.
            }

            std::optional<double> v;
            std::optional<double> U;
            std::optional<double> A;
            if (std::isfinite(out.molar_density_mol_per_m3) && out.molar_density_mol_per_m3 > 0.0) {
                v = 1.0 / out.molar_density_mol_per_m3;
                U = out.enthalpy_j_per_mol - out.pressure_pa_abs * v.value();
                A = U.value() - out.temperature_k * out.entropy_j_per_mol_k;
            }

            auto printOpt = [](const std::optional<double>& v) {
                if (v.has_value()) std::cout << v.value();
            };

            std::cout << out.temperature_k << "," << out.pressure_pa_abs << "," << getStr(o, "root", "stable") << ","
                      << (out.success ? "1" : "0") << ","
                      << phaseToString(out.phase) << ","
                      << out.molar_density_mol_per_m3 << ",";
            printOpt(v);
            std::cout << "," << out.compressibility << ","
                      << out.enthalpy_departure_j_per_mol << "," << out.entropy_departure_j_per_mol_k << "," << out.gibbs_departure_j_per_mol << ","
                      << out.enthalpy_j_per_mol << "," << out.entropy_j_per_mol_k << "," << out.gibbs_j_per_mol << ",";
            printOpt(U);
            std::cout << ",";
            printOpt(A);
            std::cout << ",";
            printOpt(mw_mix);
            std::cout << ",";
            printOpt(H_mass);
            std::cout << ",";
            printOpt(S_mass);
            std::cout << ",";
            printOpt(U_mass);
            std::cout << ",";
            printOpt(A_mass);
            std::cout << ","
                      << "\"" << lnphi.str() << "\","
                      << "\"" << (out.message ? out.message : "") << "\"\n";
        } else {
            std::cout << std::setprecision(17);
            std::cout << "temperature_k=" << out.temperature_k << ", pressure_pa_abs=" << out.pressure_pa_abs << "\n";
            std::cout << "pressure input unit=" << pressureUnitToken(p_unit) << ", solver unit=PA_ABS\n";
            std::cout << "root=" << getStr(o, "root", "stable") << ", success=" << (out.success ? "true" : "false")
                      << ", phase=" << phaseToString(out.phase) << "\n";
            std::cout << "molar_density_mol_per_m3=" << out.molar_density_mol_per_m3 << ", compressibility=" << out.compressibility << "\n";
            if (std::isfinite(out.molar_density_mol_per_m3) && out.molar_density_mol_per_m3 > 0.0) {
                const double v = 1.0 / out.molar_density_mol_per_m3;
                std::cout << "molar_volume_m3_per_mol=" << v << "\n";
                const double U = out.enthalpy_j_per_mol - out.pressure_pa_abs * v;
                const double A = U - out.temperature_k * out.entropy_j_per_mol_k;
                std::cout << "internal_energy_j_per_mol=" << U << ", helmholtz_j_per_mol=" << A << "\n";
            }
            std::cout << "enthalpy_departure_j_per_mol=" << out.enthalpy_departure_j_per_mol << ", entropy_departure_j_per_mol_k=" << out.entropy_departure_j_per_mol_k
                      << ", gibbs_departure_j_per_mol=" << out.gibbs_departure_j_per_mol << "\n";
            std::cout << "enthalpy_j_per_mol=" << out.enthalpy_j_per_mol << ", entropy_j_per_mol_k=" << out.entropy_j_per_mol_k << ", gibbs_j_per_mol=" << out.gibbs_j_per_mol
                      << "\n";

            try {
                const auto component_mws = getComponentMWs(sys, in.components.size());
                const double mw = averageMW(z, component_mws);
                std::cout << "molecular_weight_g_per_mol=" << mw << "\n";
                if (auto v = molarToMass(out.enthalpy_j_per_mol, mw)) std::cout << "enthalpy_j_per_kg=" << v.value() << "\n";
                if (auto v = molarToMass(out.entropy_j_per_mol_k, mw)) std::cout << "entropy_j_per_kg_k=" << v.value() << "\n";
                if (std::isfinite(out.molar_density_mol_per_m3) && out.molar_density_mol_per_m3 > 0.0) {
                    const double vv = 1.0 / out.molar_density_mol_per_m3;
                    const double U = out.enthalpy_j_per_mol - out.pressure_pa_abs * vv;
                    const double A = U - out.temperature_k * out.entropy_j_per_mol_k;
                    if (auto vm = molarToMass(U, mw)) std::cout << "internal_energy_j_per_kg=" << vm.value() << "\n";
                    if (auto vm = molarToMass(A, mw)) std::cout << "helmholtz_j_per_kg=" << vm.value() << "\n";
                }
            } catch (...) {
                // Best-effort only for CLI convenience output.
            }

            if (out.ln_fugacity_coeff && out.ln_fugacity_coeff_len == in.components.size()) {
                std::cout << "ln_fugacity_coeff:\n";
                for (size_t i = 0; i < in.components.size(); ++i) {
                    std::cout << "  [" << i << "] " << in.components[i] << " = " << out.ln_fugacity_coeff[i] << "\n";
                }
            }
            if (out.message && out.message[0] != '\0') {
                std::cout << "message: " << out.message << "\n";
            }
        }

        const bool ok = (out.success != 0);
        dmtherm_tp_phase_eval_destroy(&out);
        dmtherm_system_destroy(sys);
        return ok ? 0 : 2;
    } catch (...) {
        if (sys) dmtherm_system_destroy(sys);
        throw;
    }
}

int cmdEvalTPRoots(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printTPRootsHelp(std::cout);
        return 0;
    }

    const auto in = parseSystemBuildInputs(o);
    const auto T = getDouble(o, "T");
    const auto P_input = getDouble(o, "P");
    const auto p_unit = parsePressureUnit(getStr(o, "P-unit", "PA_ABS"));
    const auto P = convertPressureToPaAbs(P_input, p_unit);
    const auto z = parseComposition(o, in.components);
    const bool csv = hasFlag(o, "csv");

    dmtherm_system_t* sys = nullptr;
    try {
        sys = buildSystem(in);

        dmtherm_tp_departure_roots_t out{};
        out.struct_size = static_cast<uint32_t>(sizeof(out));

        const auto st = dmtherm_system_departure_roots_tp(sys, T, P, z.data(), z.size(), &out);
        const std::string last_error = dmtherm_system_last_error(sys);
        if (st != DMTHERM_STATUS_OK) {
            dmtherm_tp_departure_roots_destroy(&out);
            throw std::runtime_error("departure_roots_tp failed (status=" + std::to_string(st) + "): " + last_error);
        }

        if (csv) {
            std::cout << std::setprecision(17);
            std::cout << "temperature_k,pressure_pa_abs,root_index,mechanically_stable,phase,molar_density_mol_per_m3,compressibility,enthalpy_departure_j_per_mol,entropy_departure_j_per_mol_k,gibbs_departure_j_per_mol,enthalpy_j_per_mol,entropy_j_per_mol_k,gibbs_j_per_mol,ln_fugacity_coeff\n";
            for (size_t i = 0; i < out.roots_len; ++i) {
                const auto& r = out.roots[i];
                std::ostringstream lnphi;
                lnphi << std::setprecision(17);
                for (size_t j = 0; j < r.ln_fugacity_coeff_len; ++j) {
                    if (j) lnphi << ",";
                    lnphi << r.ln_fugacity_coeff[j];
                }
                std::cout << out.temperature_k << "," << out.pressure_pa_abs << "," << i << ","
                          << (r.mechanically_stable ? "1" : "0") << ","
                          << phaseToString(r.phase) << ","
                          << r.molar_density_mol_per_m3 << "," << r.compressibility << ","
                          << r.enthalpy_departure_j_per_mol << "," << r.entropy_departure_j_per_mol_k << "," << r.gibbs_departure_j_per_mol << ","
                          << r.enthalpy_j_per_mol << "," << r.entropy_j_per_mol_k << "," << r.gibbs_j_per_mol << ","
                          << "\"" << lnphi.str() << "\"\n";
            }
            if (out.message && out.message[0] != '\0') {
                std::cerr << "message: " << out.message << "\n";
            }
            std::cerr << "stable_root_index=" << out.stable_root_index << "\n";
        } else {
            std::cout << std::setprecision(17);
            std::cout << "temperature_k=" << out.temperature_k << ", pressure_pa_abs=" << out.pressure_pa_abs << "\n";
            std::cout << "pressure input unit=" << pressureUnitToken(p_unit) << ", solver unit=PA_ABS\n";
            std::cout << "stable_root_index=" << out.stable_root_index << "\n";
            if (out.message && out.message[0] != '\0') {
                std::cout << "message: " << out.message << "\n";
            }

            for (size_t i = 0; i < out.roots_len; ++i) {
                const auto& r = out.roots[i];
                std::cout << "\nroot[" << i << "] phase=" << phaseToString(r.phase)
                          << ", mechanically_stable=" << (r.mechanically_stable ? "true" : "false") << "\n";
                std::cout << "molar_density_mol_per_m3=" << r.molar_density_mol_per_m3 << ", compressibility=" << r.compressibility << "\n";
                std::cout << "enthalpy_departure_j_per_mol=" << r.enthalpy_departure_j_per_mol << ", entropy_departure_j_per_mol_k=" << r.entropy_departure_j_per_mol_k
                          << ", gibbs_departure_j_per_mol=" << r.gibbs_departure_j_per_mol << "\n";
                std::cout << "enthalpy_j_per_mol=" << r.enthalpy_j_per_mol << ", entropy_j_per_mol_k=" << r.entropy_j_per_mol_k << ", gibbs_j_per_mol=" << r.gibbs_j_per_mol
                          << "\n";
                if (r.ln_fugacity_coeff && r.ln_fugacity_coeff_len == in.components.size()) {
                    std::cout << "ln_fugacity_coeff:\n";
                    for (size_t j = 0; j < in.components.size(); ++j) {
                        std::cout << "  [" << j << "] " << in.components[j] << " = " << r.ln_fugacity_coeff[j] << "\n";
                    }
                }
            }
        }

        dmtherm_tp_departure_roots_destroy(&out);
        dmtherm_system_destroy(sys);
        return 0;
    } catch (...) {
        if (sys) dmtherm_system_destroy(sys);
        throw;
    }
}

int cmdEvalCritical(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printEvalCriticalHelp(std::cout);
        return 0;
    }

    const std::string databanks_dir = getStr(o, "databanks-dir", "data");
    const std::string components_csv = getStr(o, "components", "");
    if (components_csv.empty()) {
        throw std::invalid_argument("Missing required option: --components");
    }

    std::vector<std::string> components = splitCsv(components_csv);
    if (components.empty()) {
        throw std::invalid_argument("--components must not be empty");
    }

    std::vector<std::string> models;
    const std::string models_csv = trimCopy(getStr(o, "models", "all"));
    if (toLower(models_csv) == "all" || models_csv.empty()) {
        models = DMThermo::Factory::EOSFactory::registeredTypes();
    } else {
        models = splitCsv(models_csv);
    }

    const bool csv = hasFlag(o, "csv");

    if (!std::filesystem::exists(databanks_dir)) {
        throw std::invalid_argument("databanks dir does not exist: " + databanks_dir);
    }

    DMThermo::Data::Databanks db;
    DMThermo::Data::Diagnostics diag;
    const bool ok = db.loadAllFromDirectory(databanks_dir, &diag);
    if (!ok || diag.hasErrors()) {
        std::string msg = "Failed to load databanks from: " + databanks_dir + "\n";
        for (const auto& d : diag.items()) {
            if (d.level == DMThermo::Data::DiagnosticLevel::Error) {
                msg += "  error: " + d.message + "\n";
            }
        }
        throw std::runtime_error(msg);
    }

    // Pre-validate models.
    {
        const auto available = DMThermo::Factory::EOSFactory::registeredTypes();
        for (const auto& m : models) {
            if (std::find(available.begin(), available.end(), m) == available.end()) {
                throw std::invalid_argument("Unknown EOS model '" + m + "' (use 'dmtherm list models')");
            }
        }
    }

    const DMThermo::Data::ComponentResolver resolver(db);
    const double R = DMThermo::Constants::GAS_CONSTANT;

    auto optToDouble = [](const std::optional<double>& v) {
        return v.has_value() ? v.value() : std::numeric_limits<double>::quiet_NaN();
    };

    auto computeZcFromTPV = [&](double critical_temperature, double critical_pressure, double critical_volume) -> std::optional<double> {
        if (!(std::isfinite(critical_temperature) && critical_temperature > 0.0)) return std::nullopt;
        if (!(std::isfinite(critical_pressure) && critical_pressure > 0.0)) return std::nullopt;
        if (!(std::isfinite(critical_volume) && critical_volume > 0.0)) return std::nullopt;
        const DMThermo::Core::Units::PropertyVariable critical_volume_var(
            critical_volume,
            DMThermo::Core::Units::Unit::M3_PER_KMOL);
        const double critical_volume_molar = critical_volume_var.as(DMThermo::Core::Units::Unit::M3_PER_MOL);
        return (critical_pressure * critical_volume_molar) / (R * critical_temperature);
    };

    auto solveCriticalFromEOS = [&](const DMThermo::EOS& eos, double Tc_hint, double Pc_hint, double critical_volume_hint) {
        DMThermo::Equilibrium::Critical::EstimateInputs inputs;
        inputs.Tc_hint = Tc_hint;
        inputs.Pc_hint = Pc_hint;
        const DMThermo::Core::Units::PropertyVariable critical_volume(
            critical_volume_hint,
            DMThermo::Core::Units::Unit::M3_PER_KMOL);
        inputs.Vc_hint_m3_per_kmol = critical_volume.as(DMThermo::Core::Units::Unit::M3_PER_KMOL);
        return DMThermo::Equilibrium::Critical::estimateFromEOS(eos, {1.0}, inputs);
    };

    int failures = 0;

    if (csv) {
        std::cout << std::setprecision(17);
        std::cout << "component,model,critical_temperature_actual_k,Tc_actual_K,critical_pressure_actual_pa_abs,Pc_actual_Pa,critical_volume_actual_m3_per_kmol,Vc_actual_m3_per_kmol,critical_compressibility_actual,Zc_actual,critical_temperature_est_k,Tc_est_K,critical_pressure_est_pa_abs,Pc_est_Pa,critical_volume_est_m3_per_kmol,Vc_est_m3_per_kmol,critical_compressibility_est,Zc_est,note\n";
    }

    for (const auto& id : components) {
        auto resolved = resolver.resolve(id, &diag);
        if (!resolved.dippr) {
            std::cerr << "Component '" << id << "' not found in dippr.csv (required for actual critical properties)\n";
            ++failures;
            continue;
        }

        const auto& r = *resolved.dippr;
        const std::string name = !r.key.name.empty() ? r.key.name : id;
        const std::string ident = !r.key.cas.empty() ? r.key.cas : (!r.key.name.empty() ? r.key.name : id);

        const double Tc_actual = optToDouble(r.Tc);
        const double Pc_actual = optToDouble(r.Pc);
        const double critical_volume_actual = optToDouble(r.Vc);

        double Zc_actual = optToDouble(r.Zc);
        if (!std::isfinite(Zc_actual)) {
            if (auto z = computeZcFromTPV(Tc_actual, Pc_actual, critical_volume_actual)) {
                Zc_actual = z.value();
            }
        }

        if (!csv) {
            std::cout << "\ncomponent=" << name;
            if (!r.key.cas.empty()) std::cout << "  CAS=" << r.key.cas;
            if (r.key.chemid.has_value()) std::cout << "  CHEMID=" << r.key.chemid.value();
            std::cout << "\n";
            std::cout << "  actual: ";
            std::cout << "Tc=" << Tc_actual << " K, Pc=" << Pc_actual << " Pa";
            std::cout << ", critical_volume=" << critical_volume_actual << " m^3/kmol";
            std::cout << ", Zc=" << Zc_actual << "\n";
        }

        for (const auto& model : models) {
            double Tc_est = std::numeric_limits<double>::quiet_NaN();
            double Pc_est = std::numeric_limits<double>::quiet_NaN();
            double Zc_est = std::numeric_limits<double>::quiet_NaN();
            double critical_volume_est = std::numeric_limits<double>::quiet_NaN();
            std::string note;

            DMThermo::EOSPtr eos;
            DMThermo::Data::Diagnostics eos_diag;
            try {
                eos = DMThermo::Factory::EOSFactory::createFromDatabanks(model, db, {id}, {}, &eos_diag);
            } catch (const std::exception& e) {
                note = std::string("eos_create_failed: ") + e.what();
            }

            if (eos) {
                const auto est = solveCriticalFromEOS(*eos, Tc_actual, Pc_actual, critical_volume_actual);
                if (est.converged) {
                    Tc_est = est.Tc;
                    Pc_est = est.Pc;
                    critical_volume_est = est.Vc_m3_per_kmol;
                    Zc_est = est.Zc;
                    note = "critical_from_eos";
                    if (!est.message.empty()) {
                        note += ";";
                        note += est.message;
                    }
                } else {
                    note = std::string("critical_solve_failed: ") + (est.message.empty() ? "unknown" : est.message);
                }
            }

            // Fallback for cubic EOS: Tc/Pc are parameters, Zc is model-constant.
            if (!std::isfinite(Tc_est) || !std::isfinite(Pc_est) || !std::isfinite(critical_volume_est) || !std::isfinite(Zc_est)) {
                if ((model == "Peng-Robinson" || model == "SRK") &&
                    (std::isfinite(Tc_actual) && Tc_actual > 0.0) &&
                    (std::isfinite(Pc_actual) && Pc_actual > 0.0))
                {
                    Tc_est = Tc_actual;
                    Pc_est = Pc_actual;
                    Zc_est = (model == "Peng-Robinson") ? 0.307 : (1.0 / 3.0);
                    const double critical_volume_molar_est = (Zc_est * R * Tc_est) / Pc_est;
                    const DMThermo::Core::Units::PropertyVariable critical_volume_var(
                        critical_volume_molar_est,
                        DMThermo::Core::Units::Unit::M3_PER_MOL);
                    critical_volume_est = critical_volume_var.as(DMThermo::Core::Units::Unit::M3_PER_KMOL);
                    if (!note.empty()) note += ";";
                    note += "fallback_Zc_const";
                } else if (model == "PR-VT" &&
                    (std::isfinite(Tc_actual) && Tc_actual > 0.0) &&
                    (std::isfinite(Pc_actual) && Pc_actual > 0.0))
                {
                    Tc_est = Tc_actual;
                    Pc_est = Pc_actual;
                    Zc_est = (std::isfinite(Zc_actual) && Zc_actual > 0.0) ? Zc_actual : 0.307;
                    const double critical_volume_molar_est = (Zc_est * R * Tc_est) / Pc_est;
                    const DMThermo::Core::Units::PropertyVariable critical_volume_var(
                        critical_volume_molar_est,
                        DMThermo::Core::Units::Unit::M3_PER_MOL);
                    critical_volume_est = critical_volume_var.as(DMThermo::Core::Units::Unit::M3_PER_KMOL);
                    if (!note.empty()) note += ";";
                    note += "fallback_Zc_from_dippr_or_pr_const";
                }
            }

            if (csv) {
                std::cout << csvEscape(ident) << ","
                          << csvEscape(model) << ","
                          << Tc_actual << ","
                          << Tc_actual << ","
                          << Pc_actual << ","
                          << Pc_actual << ","
                          << critical_volume_actual << ","
                          << critical_volume_actual << ","
                          << Zc_actual << ","
                          << Zc_actual << ","
                          << Tc_est << ","
                          << Tc_est << ","
                          << Pc_est << ","
                          << Pc_est << ","
                          << critical_volume_est << ","
                          << critical_volume_est << ","
                          << Zc_est << ","
                          << Zc_est << ","
                          << csvEscape(note)
                          << "\n";
            } else {
                std::cout << "  est(" << model << "): ";
                std::cout << "Tc=" << Tc_est << " K, Pc=" << Pc_est << " Pa";
                std::cout << ", critical_volume=" << critical_volume_est << " m^3/kmol";
                std::cout << ", Zc=" << Zc_est;
                if (!note.empty()) std::cout << " (" << note << ")";
                std::cout << "\n";
            }
        }
    }

    return (failures == 0) ? 0 : 2;
}

int cmdFlashPT(int argc, char** argv, int start_index) {
    return cmdFlashGeneric(
        argc,
        argv,
        start_index,
        "pt",
        &dmtherm_system_flash_pt,
        &printFlashPTHelp,
        "T",
        "P");
}

int cmdFlashTV(int argc, char** argv, int start_index) {
    return cmdFlashGeneric(
        argc,
        argv,
        start_index,
        "tv",
        &dmtherm_system_flash_tv,
        &printFlashTVHelp,
        "T",
        "V");
}

int cmdFlashPH(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printFlashPHHelp(std::cout);
        return 0;
    }

    const auto in = parseSystemBuildInputs(o);
    const auto P_input = getDouble(o, "P");
    const auto p_unit = parsePressureUnit(getStr(o, "P-unit", "PA_ABS"));
    const auto P = convertPressureToPaAbs(P_input, p_unit);
    const auto H = getDouble(o, "H");
    const auto z = parseComposition(o, in.components);
    const int max_phases = getInt(o, "max-phases", 2);
    const bool compute_properties = hasFlag(o, "properties");
    const bool json = hasFlag(o, "json");

    dmtherm_flash_hs_config_t hs_cfg{};
    dmtherm_flash_hs_config_set_defaults(&hs_cfg);
    hs_cfg.tolerance = std::stod(getStr(o, "enthalpy-tol", "0"));
    hs_cfg.temp_bracket_low = std::stod(getStr(o, "temp-low", "0"));
    hs_cfg.temp_bracket_high = std::stod(getStr(o, "temp-high", "0"));
    hs_cfg.max_iterations = getInt(o, "max-iterations", 0);

    dmtherm_system_t* sys = nullptr;
    try {
        sys = buildSystem(in);
        std::vector<double> component_mws;
        if (!json) {
            component_mws = getComponentMWs(sys, in.components.size());
        }

        dmtherm_flash_config_t fcfg{};
        dmtherm_flash_config_set_defaults(&fcfg);
        fcfg.max_phases = max_phases;
        fcfg.compute_properties = compute_properties ? 1 : 0;
        fcfg.reserved_ptr[0] = &hs_cfg;

        dmtherm_flash_pt_result_t out{};
        out.struct_size = static_cast<uint32_t>(sizeof(out));

        const auto st = dmtherm_system_flash_ph(sys, P, H, z.data(), z.size(), &fcfg, &out);
        const std::string last_error = dmtherm_system_last_error(sys);
        if (st != DMTHERM_STATUS_OK) {
            dmtherm_flash_pt_result_destroy(&out);
            throw std::runtime_error(std::string("flash ph failed (status=") + std::to_string(st) + "): " + last_error);
        }

        if (json) {
            printFlashResultJson(std::cout, "ph", in, z, max_phases, compute_properties, "P", P, "H", H, out);
        } else {
            std::cout << "mode=ph  P=" << P << "  H=" << H << "\n";
            std::cout << "pressure input unit=" << pressureUnitToken(p_unit) << ", solver unit=PA_ABS\n";
            printFlashResult(out, in.components, z, component_mws);
        }

        const bool ok = (out.converged != 0);
        dmtherm_flash_pt_result_destroy(&out);
        dmtherm_system_destroy(sys);
        return ok ? 0 : 2;
    } catch (...) {
        if (sys) dmtherm_system_destroy(sys);
        throw;
    }
}

int cmdFlashPS(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printFlashPSHelp(std::cout);
        return 0;
    }

    const auto in = parseSystemBuildInputs(o);
    const auto P_input = getDouble(o, "P");
    const auto p_unit = parsePressureUnit(getStr(o, "P-unit", "PA_ABS"));
    const auto P = convertPressureToPaAbs(P_input, p_unit);
    const auto S = getDouble(o, "S");
    const auto z = parseComposition(o, in.components);
    const int max_phases = getInt(o, "max-phases", 2);
    const bool compute_properties = hasFlag(o, "properties");
    const bool json = hasFlag(o, "json");

    dmtherm_flash_hs_config_t hs_cfg{};
    dmtherm_flash_hs_config_set_defaults(&hs_cfg);
    hs_cfg.tolerance = std::stod(getStr(o, "entropy-tol", "0"));
    hs_cfg.temp_bracket_low = std::stod(getStr(o, "temp-low", "0"));
    hs_cfg.temp_bracket_high = std::stod(getStr(o, "temp-high", "0"));
    hs_cfg.max_iterations = getInt(o, "max-iterations", 0);

    dmtherm_system_t* sys = nullptr;
    try {
        sys = buildSystem(in);
        std::vector<double> component_mws;
        if (!json) {
            component_mws = getComponentMWs(sys, in.components.size());
        }

        dmtherm_flash_config_t fcfg{};
        dmtherm_flash_config_set_defaults(&fcfg);
        fcfg.max_phases = max_phases;
        fcfg.compute_properties = compute_properties ? 1 : 0;
        fcfg.reserved_ptr[0] = &hs_cfg;

        dmtherm_flash_pt_result_t out{};
        out.struct_size = static_cast<uint32_t>(sizeof(out));

        const auto st = dmtherm_system_flash_ps(sys, P, S, z.data(), z.size(), &fcfg, &out);
        const std::string last_error = dmtherm_system_last_error(sys);
        if (st != DMTHERM_STATUS_OK) {
            dmtherm_flash_pt_result_destroy(&out);
            throw std::runtime_error(std::string("flash ps failed (status=") + std::to_string(st) + "): " + last_error);
        }

        if (json) {
            printFlashResultJson(std::cout, "ps", in, z, max_phases, compute_properties, "P", P, "S", S, out);
        } else {
            std::cout << "mode=ps  P=" << P << "  S=" << S << "\n";
            std::cout << "pressure input unit=" << pressureUnitToken(p_unit) << ", solver unit=PA_ABS\n";
            printFlashResult(out, in.components, z, component_mws);
        }

        const bool ok = (out.converged != 0);
        dmtherm_flash_pt_result_destroy(&out);
        dmtherm_system_destroy(sys);
        return ok ? 0 : 2;
    } catch (...) {
        if (sys) dmtherm_system_destroy(sys);
        throw;
    }
}

int cmdSweepTPGrid(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printSweepTPGridHelp(std::cout);
        return 0;
    }

    const auto in = parseSystemBuildInputs(o);
    const auto z = parseComposition(o, in.components);
    const std::string root_name = getStr(o, "root", "stable");
    const auto root = parseRootSelection(root_name);
    const auto p_unit = parsePressureUnit(getStr(o, "P-unit", "PA_ABS"));

    const double T_min = getDouble(o, "T-min");
    const double T_max = getDouble(o, "T-max");
    const int nT = getInt(o, "nT", -1);
    const double P_min = convertPressureToPaAbs(getDouble(o, "P-min"), p_unit);
    const double P_max = convertPressureToPaAbs(getDouble(o, "P-max"), p_unit);
    const int nP = getInt(o, "nP", -1);
    const bool p_log = hasFlag(o, "P-log");
    const bool no_header = hasFlag(o, "no-header");
    const bool json = hasFlag(o, "json");
    const std::string out_path = getStr(o, "out", "");

    if (!(std::isfinite(T_min) && std::isfinite(T_max) && T_min > 0.0 && T_max >= T_min)) {
        throw std::invalid_argument("Invalid T range: require T-min > 0 and T-max >= T-min");
    }
    if (!(std::isfinite(P_min) && std::isfinite(P_max) && P_min > 0.0 && P_max >= P_min)) {
        throw std::invalid_argument("Invalid P range: require P-min > 0 and P-max >= P-min");
    }
    if (nT <= 0 || nP <= 0) {
        throw std::invalid_argument("nT and nP must be > 0");
    }
    if (p_log && (P_min <= 0.0 || P_max <= 0.0)) {
        throw std::invalid_argument("P-log requires P-min and P-max > 0");
    }

    std::vector<double> Ts(static_cast<size_t>(nT));
    if (nT == 1) {
        Ts[0] = T_min;
    } else {
        for (int i = 0; i < nT; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(nT - 1);
            Ts[static_cast<size_t>(i)] = T_min + t * (T_max - T_min);
        }
    }

    std::vector<double> Ps(static_cast<size_t>(nP));
    if (nP == 1) {
        Ps[0] = P_min;
    } else if (p_log) {
        const double a = std::log10(P_min);
        const double b = std::log10(P_max);
        for (int i = 0; i < nP; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(nP - 1);
            Ps[static_cast<size_t>(i)] = std::pow(10.0, a + t * (b - a));
        }
    } else {
        for (int i = 0; i < nP; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(nP - 1);
            Ps[static_cast<size_t>(i)] = P_min + t * (P_max - P_min);
        }
    }

    std::ofstream out_file;
    std::ostream* os = &std::cout;
    if (!out_path.empty()) {
        out_file.open(out_path, std::ios::out | std::ios::binary);
        if (!out_file) {
            throw std::runtime_error("Failed to open --out file: " + out_path);
        }
        os = &out_file;
    }

    dmtherm_system_t* sys = nullptr;
    try {
        sys = buildSystem(in);

        (*os) << std::setprecision(17);
        if (!json && !no_header) {
            (*os) << "temperature_k,pressure_pa_abs,root,success,phase,molar_density_mol_per_m3,compressibility,enthalpy_departure_j_per_mol,entropy_departure_j_per_mol_k,gibbs_departure_j_per_mol,enthalpy_j_per_mol,entropy_j_per_mol_k,gibbs_j_per_mol,ln_fugacity_coeff,message\n";
        }
        if (json) {
            (*os) << "{";
            (*os) << "\"schema\":\"dmtherm.sweep.tp_grid.v1\",";
            (*os) << "\"meta\":{";
            (*os) << "\"databanks_dir\":\"" << jsonEscape(in.databanks_dir) << "\",";
            (*os) << "\"eos\":\"" << jsonEscape(in.eos_type) << "\",";
            (*os) << "\"flash_method\":\"" << jsonEscape(in.flash_method) << "\",";
            (*os) << "\"property_method\":\"" << jsonEscape(getStr(o, "property-method", "eos_only")) << "\",";
            (*os) << "\"require_critical\":" << (in.require_critical_properties ? "true" : "false") << ",";
            (*os) << "\"require_ideal_gas_cp\":" << (in.require_ideal_gas_cp ? "true" : "false") << ",";
            (*os) << "\"root\":\"" << jsonEscape(root_name) << "\",";
            (*os) << "\"temperature_min_k\":" << T_min << ",";
            (*os) << "\"temperature_max_k\":" << T_max << ",";
            (*os) << "\"nT\":" << nT << ",";
            (*os) << "\"pressure_min_pa_abs\":" << P_min << ",";
            (*os) << "\"pressure_max_pa_abs\":" << P_max << ",";
            (*os) << "\"nP\":" << nP << ",";
            (*os) << "\"pressure_log\":" << (p_log ? "true" : "false") << ",";
            (*os) << "\"P_unit\":\"PA_ABS\",";
            (*os) << "\"P_input_unit\":\"" << pressureUnitToken(p_unit) << "\",";
            (*os) << "\"components\":[";
            for (size_t i = 0; i < in.components.size(); ++i) {
                if (i) (*os) << ",";
                (*os) << "\"" << jsonEscape(in.components[i]) << "\"";
            }
            (*os) << "],";
            (*os) << "\"z\":[";
            for (size_t i = 0; i < z.size(); ++i) {
                if (i) (*os) << ",";
                (*os) << z[i];
            }
            (*os) << "]";
            (*os) << "},";
            (*os) << "\"points\":[";
        }

        int failures = 0;
        dmtherm_tp_phase_eval_t out{};
        out.struct_size = static_cast<uint32_t>(sizeof(out));
        bool first_json = true;

        for (double T : Ts) {
            for (double P : Ps) {
                const auto st = dmtherm_system_evaluate_tp_phase(sys, T, P, z.data(), z.size(), root, &out);
                const std::string last_error = dmtherm_system_last_error(sys);

                std::string msg;
                if (out.message && out.message[0] != '\0') msg = out.message;
                else if (st != DMTHERM_STATUS_OK && !last_error.empty()) msg = last_error;

                std::ostringstream lnphi;
                lnphi << std::setprecision(17);
                if (out.ln_fugacity_coeff && out.ln_fugacity_coeff_len) {
                    for (size_t i = 0; i < out.ln_fugacity_coeff_len; ++i) {
                        if (i) lnphi << ",";
                        lnphi << out.ln_fugacity_coeff[i];
                    }
                }

                const int ok_eval = (st == DMTHERM_STATUS_OK && out.success == 1) ? 1 : 0;
                if (!ok_eval) ++failures;

                const dmtherm_phase_t ph = (st == DMTHERM_STATUS_OK) ? out.phase : DMTHERM_PHASE_UNKNOWN;
                const double rho = (st == DMTHERM_STATUS_OK) ? out.molar_density_mol_per_m3 : 0.0;
                const double Z = (st == DMTHERM_STATUS_OK) ? out.compressibility : 0.0;
                const double h_dep = (st == DMTHERM_STATUS_OK) ? out.enthalpy_departure_j_per_mol : 0.0;
                const double s_dep = (st == DMTHERM_STATUS_OK) ? out.entropy_departure_j_per_mol_k : 0.0;
                const double g_dep = (st == DMTHERM_STATUS_OK) ? out.gibbs_departure_j_per_mol : 0.0;
                const double H = (st == DMTHERM_STATUS_OK) ? out.enthalpy_j_per_mol : 0.0;
                const double S = (st == DMTHERM_STATUS_OK) ? out.entropy_j_per_mol_k : 0.0;
                const double G = (st == DMTHERM_STATUS_OK) ? out.gibbs_j_per_mol : 0.0;

                if (json) {
                    if (!first_json) (*os) << ",";
                    first_json = false;
                    jsonWriteTPPointCanonical(
                        *os,
                        T,
                        P,
                        root_name,
                        st,
                        ok_eval != 0,
                        ph,
                        rho,
                        Z,
                        h_dep,
                        s_dep,
                        g_dep,
                        H,
                        S,
                        G,
                        (st == DMTHERM_STATUS_OK) ? out.ln_fugacity_coeff : nullptr,
                        (st == DMTHERM_STATUS_OK) ? out.ln_fugacity_coeff_len : 0,
                        msg);
                } else {
                    (*os) << T << "," << P << "," << root_name << ","
                          << ok_eval << ","
                          << phaseToString(ph) << ","
                          << rho << "," << Z << ","
                          << h_dep << "," << s_dep << "," << g_dep << ","
                          << H << "," << S << "," << G << ","
                          << csvEscape(lnphi.str()) << ","
                          << csvEscape(msg) << "\n";
                }

                dmtherm_tp_phase_eval_destroy(&out);
            }
        }

        if (json) {
            (*os) << "]}";
            (*os) << "\n";
        }

        dmtherm_system_destroy(sys);
        return (failures == 0) ? 0 : 2;
    } catch (...) {
        if (sys) dmtherm_system_destroy(sys);
        throw;
    }
}

int cmdSweepIsotherm(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printSweepIsothermHelp(std::cout);
        return 0;
    }

    const auto in = parseSystemBuildInputs(o);
    const auto z = parseComposition(o, in.components);
    const std::string root_name = getStr(o, "root", "stable");
    const auto root = parseRootSelection(root_name);
    const auto p_unit = parsePressureUnit(getStr(o, "P-unit", "PA_ABS"));
    const double T = getDouble(o, "T");
    const double P_min = convertPressureToPaAbs(getDouble(o, "P-min"), p_unit);
    const double P_max = convertPressureToPaAbs(getDouble(o, "P-max"), p_unit);
    const int nP = getInt(o, "nP", -1);
    const bool p_log = hasFlag(o, "P-log");
    const bool no_header = hasFlag(o, "no-header");
    const bool json = hasFlag(o, "json");
    const std::string out_path = getStr(o, "out", "");

    if (!(std::isfinite(T) && T > 0.0)) throw std::invalid_argument("Invalid T (must be > 0)");
    if (!(std::isfinite(P_min) && std::isfinite(P_max) && P_min > 0.0 && P_max >= P_min)) {
        throw std::invalid_argument("Invalid P range: require P-min > 0 and P-max >= P-min");
    }
    if (nP <= 0) throw std::invalid_argument("nP must be > 0");
    if (p_log && (P_min <= 0.0 || P_max <= 0.0)) throw std::invalid_argument("P-log requires P-min and P-max > 0");

    std::vector<double> Ps(static_cast<size_t>(nP));
    if (nP == 1) {
        Ps[0] = P_min;
    } else if (p_log) {
        const double a = std::log10(P_min);
        const double b = std::log10(P_max);
        for (int i = 0; i < nP; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(nP - 1);
            Ps[static_cast<size_t>(i)] = std::pow(10.0, a + t * (b - a));
        }
    } else {
        for (int i = 0; i < nP; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(nP - 1);
            Ps[static_cast<size_t>(i)] = P_min + t * (P_max - P_min);
        }
    }

    std::ofstream out_file;
    std::ostream* os = &std::cout;
    if (!out_path.empty()) {
        out_file.open(out_path, std::ios::out | std::ios::binary);
        if (!out_file) throw std::runtime_error("Failed to open --out file: " + out_path);
        os = &out_file;
    }

    dmtherm_system_t* sys = nullptr;
    try {
        sys = buildSystem(in);

        (*os) << std::setprecision(17);
        if (!json && !no_header) {
            (*os) << "temperature_k,pressure_pa_abs,root,success,phase,molar_density_mol_per_m3,compressibility,enthalpy_departure_j_per_mol,entropy_departure_j_per_mol_k,gibbs_departure_j_per_mol,enthalpy_j_per_mol,entropy_j_per_mol_k,gibbs_j_per_mol,ln_fugacity_coeff,message\n";
        }
        if (json) {
            (*os) << "{";
            (*os) << "\"schema\":\"dmtherm.sweep.isotherm.v1\",";
            (*os) << "\"meta\":{";
            (*os) << "\"databanks_dir\":\"" << jsonEscape(in.databanks_dir) << "\",";
            (*os) << "\"eos\":\"" << jsonEscape(in.eos_type) << "\",";
            (*os) << "\"flash_method\":\"" << jsonEscape(in.flash_method) << "\",";
            (*os) << "\"property_method\":\"" << jsonEscape(getStr(o, "property-method", "eos_only")) << "\",";
            (*os) << "\"require_critical\":" << (in.require_critical_properties ? "true" : "false") << ",";
            (*os) << "\"require_ideal_gas_cp\":" << (in.require_ideal_gas_cp ? "true" : "false") << ",";
            (*os) << "\"root\":\"" << jsonEscape(root_name) << "\",";
            (*os) << "\"temperature_k\":" << T << ",";
            (*os) << "\"pressure_min_pa_abs\":" << P_min << ",";
            (*os) << "\"pressure_max_pa_abs\":" << P_max << ",";
            (*os) << "\"nP\":" << nP << ",";
            (*os) << "\"pressure_log\":" << (p_log ? "true" : "false") << ",";
            (*os) << "\"P_unit\":\"PA_ABS\",";
            (*os) << "\"P_input_unit\":\"" << pressureUnitToken(p_unit) << "\",";
            (*os) << "\"components\":[";
            for (size_t i = 0; i < in.components.size(); ++i) {
                if (i) (*os) << ",";
                (*os) << "\"" << jsonEscape(in.components[i]) << "\"";
            }
            (*os) << "],";
            (*os) << "\"z\":[";
            for (size_t i = 0; i < z.size(); ++i) {
                if (i) (*os) << ",";
                (*os) << z[i];
            }
            (*os) << "]";
            (*os) << "},";
            (*os) << "\"points\":[";
        }

        int failures = 0;
        dmtherm_tp_phase_eval_t out{};
        out.struct_size = static_cast<uint32_t>(sizeof(out));
        bool first_json = true;

        for (double P : Ps) {
            const auto st = dmtherm_system_evaluate_tp_phase(sys, T, P, z.data(), z.size(), root, &out);
            const std::string last_error = dmtherm_system_last_error(sys);

            std::string msg;
            if (out.message && out.message[0] != '\0') msg = out.message;
            else if (st != DMTHERM_STATUS_OK && !last_error.empty()) msg = last_error;

            std::ostringstream lnphi;
            lnphi << std::setprecision(17);
            if (out.ln_fugacity_coeff && out.ln_fugacity_coeff_len) {
                for (size_t i = 0; i < out.ln_fugacity_coeff_len; ++i) {
                    if (i) lnphi << ",";
                    lnphi << out.ln_fugacity_coeff[i];
                }
            }

            const int ok_eval = (st == DMTHERM_STATUS_OK && out.success == 1) ? 1 : 0;
            if (!ok_eval) ++failures;

            const dmtherm_phase_t ph = (st == DMTHERM_STATUS_OK) ? out.phase : DMTHERM_PHASE_UNKNOWN;
            const double rho = (st == DMTHERM_STATUS_OK) ? out.molar_density_mol_per_m3 : 0.0;
            const double Z = (st == DMTHERM_STATUS_OK) ? out.compressibility : 0.0;
            const double h_dep = (st == DMTHERM_STATUS_OK) ? out.enthalpy_departure_j_per_mol : 0.0;
            const double s_dep = (st == DMTHERM_STATUS_OK) ? out.entropy_departure_j_per_mol_k : 0.0;
            const double g_dep = (st == DMTHERM_STATUS_OK) ? out.gibbs_departure_j_per_mol : 0.0;
            const double H = (st == DMTHERM_STATUS_OK) ? out.enthalpy_j_per_mol : 0.0;
            const double S = (st == DMTHERM_STATUS_OK) ? out.entropy_j_per_mol_k : 0.0;
            const double G = (st == DMTHERM_STATUS_OK) ? out.gibbs_j_per_mol : 0.0;

            if (json) {
                if (!first_json) (*os) << ",";
                first_json = false;
                jsonWriteTPPointCanonical(
                    *os,
                    T,
                    P,
                    root_name,
                    st,
                    ok_eval != 0,
                    ph,
                    rho,
                    Z,
                    h_dep,
                    s_dep,
                    g_dep,
                    H,
                    S,
                    G,
                    (st == DMTHERM_STATUS_OK) ? out.ln_fugacity_coeff : nullptr,
                    (st == DMTHERM_STATUS_OK) ? out.ln_fugacity_coeff_len : 0,
                    msg);
            } else {
                (*os) << T << "," << P << "," << root_name << ","
                      << ok_eval << ","
                      << phaseToString(ph) << ","
                      << rho << "," << Z << ","
                      << h_dep << "," << s_dep << "," << g_dep << ","
                      << H << "," << S << "," << G << ","
                      << csvEscape(lnphi.str()) << ","
                      << csvEscape(msg) << "\n";
            }

            dmtherm_tp_phase_eval_destroy(&out);
        }

        if (json) {
            (*os) << "]}";
            (*os) << "\n";
        }

        dmtherm_system_destroy(sys);
        return (failures == 0) ? 0 : 2;
    } catch (...) {
        if (sys) dmtherm_system_destroy(sys);
        throw;
    }
}

int cmdSweepIsobar(int argc, char** argv, int start_index) {
    const auto o = parseOptions(argc, argv, start_index);
    if (hasFlag(o, "help")) {
        printSweepIsobarHelp(std::cout);
        return 0;
    }

    const auto in = parseSystemBuildInputs(o);
    const auto z = parseComposition(o, in.components);
    const std::string root_name = getStr(o, "root", "stable");
    const auto root = parseRootSelection(root_name);
    const auto p_unit = parsePressureUnit(getStr(o, "P-unit", "PA_ABS"));
    const double P = convertPressureToPaAbs(getDouble(o, "P"), p_unit);
    const double T_min = getDouble(o, "T-min");
    const double T_max = getDouble(o, "T-max");
    const int nT = getInt(o, "nT", -1);
    const bool no_header = hasFlag(o, "no-header");
    const bool json = hasFlag(o, "json");
    const std::string out_path = getStr(o, "out", "");

    if (!(std::isfinite(P) && P > 0.0)) throw std::invalid_argument("Invalid P (must be > 0)");
    if (!(std::isfinite(T_min) && std::isfinite(T_max) && T_min > 0.0 && T_max >= T_min)) {
        throw std::invalid_argument("Invalid T range: require T-min > 0 and T-max >= T-min");
    }
    if (nT <= 0) throw std::invalid_argument("nT must be > 0");

    std::vector<double> Ts(static_cast<size_t>(nT));
    if (nT == 1) {
        Ts[0] = T_min;
    } else {
        for (int i = 0; i < nT; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(nT - 1);
            Ts[static_cast<size_t>(i)] = T_min + t * (T_max - T_min);
        }
    }

    std::ofstream out_file;
    std::ostream* os = &std::cout;
    if (!out_path.empty()) {
        out_file.open(out_path, std::ios::out | std::ios::binary);
        if (!out_file) throw std::runtime_error("Failed to open --out file: " + out_path);
        os = &out_file;
    }

    dmtherm_system_t* sys = nullptr;
    try {
        sys = buildSystem(in);

        (*os) << std::setprecision(17);
        if (!json && !no_header) {
            (*os) << "temperature_k,pressure_pa_abs,root,success,phase,molar_density_mol_per_m3,compressibility,enthalpy_departure_j_per_mol,entropy_departure_j_per_mol_k,gibbs_departure_j_per_mol,enthalpy_j_per_mol,entropy_j_per_mol_k,gibbs_j_per_mol,ln_fugacity_coeff,message\n";
        }
        if (json) {
            (*os) << "{";
            (*os) << "\"schema\":\"dmtherm.sweep.isobar.v1\",";
            (*os) << "\"meta\":{";
            (*os) << "\"databanks_dir\":\"" << jsonEscape(in.databanks_dir) << "\",";
            (*os) << "\"eos\":\"" << jsonEscape(in.eos_type) << "\",";
            (*os) << "\"flash_method\":\"" << jsonEscape(in.flash_method) << "\",";
            (*os) << "\"property_method\":\"" << jsonEscape(getStr(o, "property-method", "eos_only")) << "\",";
            (*os) << "\"require_critical\":" << (in.require_critical_properties ? "true" : "false") << ",";
            (*os) << "\"require_ideal_gas_cp\":" << (in.require_ideal_gas_cp ? "true" : "false") << ",";
            (*os) << "\"root\":\"" << jsonEscape(root_name) << "\",";
            (*os) << "\"pressure_pa_abs\":" << P << ",";
            (*os) << "\"P_unit\":\"PA_ABS\",";
            (*os) << "\"P_input_unit\":\"" << pressureUnitToken(p_unit) << "\",";
            (*os) << "\"temperature_min_k\":" << T_min << ",";
            (*os) << "\"temperature_max_k\":" << T_max << ",";
            (*os) << "\"nT\":" << nT << ",";
            (*os) << "\"components\":[";
            for (size_t i = 0; i < in.components.size(); ++i) {
                if (i) (*os) << ",";
                (*os) << "\"" << jsonEscape(in.components[i]) << "\"";
            }
            (*os) << "],";
            (*os) << "\"z\":[";
            for (size_t i = 0; i < z.size(); ++i) {
                if (i) (*os) << ",";
                (*os) << z[i];
            }
            (*os) << "]";
            (*os) << "},";
            (*os) << "\"points\":[";
        }

        int failures = 0;
        dmtherm_tp_phase_eval_t out{};
        out.struct_size = static_cast<uint32_t>(sizeof(out));
        bool first_json = true;

        for (double T : Ts) {
            const auto st = dmtherm_system_evaluate_tp_phase(sys, T, P, z.data(), z.size(), root, &out);
            const std::string last_error = dmtherm_system_last_error(sys);

            std::string msg;
            if (out.message && out.message[0] != '\0') msg = out.message;
            else if (st != DMTHERM_STATUS_OK && !last_error.empty()) msg = last_error;

            std::ostringstream lnphi;
            lnphi << std::setprecision(17);
            if (out.ln_fugacity_coeff && out.ln_fugacity_coeff_len) {
                for (size_t i = 0; i < out.ln_fugacity_coeff_len; ++i) {
                    if (i) lnphi << ",";
                    lnphi << out.ln_fugacity_coeff[i];
                }
            }

            const int ok_eval = (st == DMTHERM_STATUS_OK && out.success == 1) ? 1 : 0;
            if (!ok_eval) ++failures;

            const dmtherm_phase_t ph = (st == DMTHERM_STATUS_OK) ? out.phase : DMTHERM_PHASE_UNKNOWN;
            const double rho = (st == DMTHERM_STATUS_OK) ? out.molar_density_mol_per_m3 : 0.0;
            const double Z = (st == DMTHERM_STATUS_OK) ? out.compressibility : 0.0;
            const double h_dep = (st == DMTHERM_STATUS_OK) ? out.enthalpy_departure_j_per_mol : 0.0;
            const double s_dep = (st == DMTHERM_STATUS_OK) ? out.entropy_departure_j_per_mol_k : 0.0;
            const double g_dep = (st == DMTHERM_STATUS_OK) ? out.gibbs_departure_j_per_mol : 0.0;
            const double H = (st == DMTHERM_STATUS_OK) ? out.enthalpy_j_per_mol : 0.0;
            const double S = (st == DMTHERM_STATUS_OK) ? out.entropy_j_per_mol_k : 0.0;
            const double G = (st == DMTHERM_STATUS_OK) ? out.gibbs_j_per_mol : 0.0;

            if (json) {
                if (!first_json) (*os) << ",";
                first_json = false;
                jsonWriteTPPointCanonical(
                    *os,
                    T,
                    P,
                    root_name,
                    st,
                    ok_eval != 0,
                    ph,
                    rho,
                    Z,
                    h_dep,
                    s_dep,
                    g_dep,
                    H,
                    S,
                    G,
                    (st == DMTHERM_STATUS_OK) ? out.ln_fugacity_coeff : nullptr,
                    (st == DMTHERM_STATUS_OK) ? out.ln_fugacity_coeff_len : 0,
                    msg);
            } else {
                (*os) << T << "," << P << "," << root_name << ","
                      << ok_eval << ","
                      << phaseToString(ph) << ","
                      << rho << "," << Z << ","
                      << h_dep << "," << s_dep << "," << g_dep << ","
                      << H << "," << S << "," << G << ","
                      << csvEscape(lnphi.str()) << ","
                      << csvEscape(msg) << "\n";
            }

            dmtherm_tp_phase_eval_destroy(&out);
        }

        if (json) {
            (*os) << "]}";
            (*os) << "\n";
        }

        dmtherm_system_destroy(sys);
        return (failures == 0) ? 0 : 2;
    } catch (...) {
        if (sys) dmtherm_system_destroy(sys);
        throw;
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            return runShell();
        }
        return runCommand(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

namespace {
int runCommand(int argc, char** argv) {
    if (argc < 2) {
        printGlobalHelp(std::cerr);
        return 2;
    }

    const std::string cmd1 = argv[1];
    if (cmd1 == "--help" || cmd1 == "help" || cmd1 == "-h") {
        printGlobalHelp(std::cout);
        return 0;
    }
    if (cmd1 == "--version") {
        std::cout << dmtherm_version_string() << "\n";
        std::cout << "ABI " << dmtherm_abi_version() << "\n";
        return 0;
    }
    if (cmd1 == "shell") {
        return runShell();
    }

    if (cmd1 == "list") {
        if (argc < 3) {
            printGlobalHelp(std::cerr);
            return 2;
        }
        const std::string cmd2 = argv[2];
        if (cmd2 == "components") {
            return cmdListComponents(argc, argv, 3);
        }
        if (cmd2 == "models") {
            return cmdListModels(argc, argv, 3);
        }
        std::cerr << "Unknown list subcommand: " << cmd2 << "\n";
        return 2;
    }

    if (cmd1 == "eval") {
        if (argc < 3) {
            printGlobalHelp(std::cerr);
            return 2;
        }
        const std::string cmd2 = argv[2];
        if (cmd2 == "critical") {
            return cmdEvalCritical(argc, argv, 3);
        }
        if (cmd2 == "tp-phase") {
            return cmdEvalTPPhase(argc, argv, 3);
        }
        if (cmd2 == "tp-roots") {
            return cmdEvalTPRoots(argc, argv, 3);
        }
        std::cerr << "Unknown eval subcommand: " << cmd2 << "\n";
        return 2;
    }

    if (cmd1 == "flash") {
        if (argc < 3) {
            printGlobalHelp(std::cerr);
            return 2;
        }
        const std::string cmd2 = argv[2];
        if (cmd2 == "pt") {
            return cmdFlashPT(argc, argv, 3);
        }
        if (cmd2 == "tv") {
            return cmdFlashTV(argc, argv, 3);
        }
        if (cmd2 == "ph") {
            return cmdFlashPH(argc, argv, 3);
        }
        if (cmd2 == "ps") {
            return cmdFlashPS(argc, argv, 3);
        }
        std::cerr << "Unknown flash subcommand: " << cmd2 << "\n";
        return 2;
    }

    if (cmd1 == "sweep") {
        if (argc < 3) {
            printGlobalHelp(std::cerr);
            return 2;
        }
        const std::string cmd2 = argv[2];
        if (cmd2 == "tp-grid") {
            return cmdSweepTPGrid(argc, argv, 3);
        }
        if (cmd2 == "isotherm") {
            return cmdSweepIsotherm(argc, argv, 3);
        }
        if (cmd2 == "isobar") {
            return cmdSweepIsobar(argc, argv, 3);
        }
        std::cerr << "Unknown sweep subcommand: " << cmd2 << "\n";
        return 2;
    }

    std::cerr << "Unknown command: " << cmd1 << "\n";
    printGlobalHelp(std::cerr);
    return 2;
}
} // namespace
