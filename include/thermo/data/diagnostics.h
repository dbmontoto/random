/**
 * @file diagnostics.h
 * @brief Simple diagnostics collector for databank loading and lookups
 */

#ifndef THERMO_DATA_DIAGNOSTICS_H
#define THERMO_DATA_DIAGNOSTICS_H

#include "thermo/config/error_handling.h"

#include <string>
#include <utility>
#include <vector>

namespace DMThermo {
namespace Data {

enum class DiagnosticLevel {
    Info,
    Warning,
    Error
};

struct Diagnostic {
    DiagnosticLevel level = DiagnosticLevel::Info;
    Config::ErrorCategory category = Config::ErrorCategory::None;
    std::string code;
    std::string message;
};

class Diagnostics {
public:
    void info(
        std::string message,
        Config::ErrorCategory category = Config::ErrorCategory::None,
        std::string code = {})
    {
        items_.push_back({DiagnosticLevel::Info, category, std::move(code), std::move(message)});
    }

    void warn(
        std::string message,
        Config::ErrorCategory category = Config::ErrorCategory::None,
        std::string code = {})
    {
        items_.push_back({DiagnosticLevel::Warning, category, std::move(code), std::move(message)});
    }

    void error(
        std::string message,
        Config::ErrorCategory category = Config::ErrorCategory::None,
        std::string code = {})
    {
        items_.push_back({DiagnosticLevel::Error, category, std::move(code), std::move(message)});
    }

    bool hasErrors() const;
    const std::vector<Diagnostic>& items() const { return items_; }

private:
    std::vector<Diagnostic> items_;
};

} // namespace Data
} // namespace DMThermo

#endif // THERMO_DATA_DIAGNOSTICS_H
