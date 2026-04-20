/**
 * @file error_handling.h
 * @brief Shared error taxonomy and failure policy controls.
 */

#ifndef THERMO_CONFIG_ERROR_HANDLING_H
#define THERMO_CONFIG_ERROR_HANDLING_H

namespace DMThermo {
namespace Config {

enum class ErrorCategory {
    None,
    InvalidInput,
    MissingData,
    NumericalFailure,
    NotInitialized,
    UnsupportedOperation
};

enum class FailureMode {
    Strict,  // throw on categorized failures
    Lenient  // return failure objects where possible
};

inline const char* errorCategoryToString(ErrorCategory category) {
    switch (category) {
        case ErrorCategory::None: return "None";
        case ErrorCategory::InvalidInput: return "InvalidInput";
        case ErrorCategory::MissingData: return "MissingData";
        case ErrorCategory::NumericalFailure: return "NumericalFailure";
        case ErrorCategory::NotInitialized: return "NotInitialized";
        case ErrorCategory::UnsupportedOperation: return "UnsupportedOperation";
        default: return "Unknown";
    }
}

inline const char* failureModeToString(FailureMode mode) {
    switch (mode) {
        case FailureMode::Strict: return "Strict";
        case FailureMode::Lenient: return "Lenient";
        default: return "Unknown";
    }
}

} // namespace Config
} // namespace DMThermo

#endif // THERMO_CONFIG_ERROR_HANDLING_H
