#include "thermo/data/diagnostics.h"

namespace DMThermo {
namespace Data {

bool Diagnostics::hasErrors() const {
    for (const auto& d : items_) {
        if (d.level == DiagnosticLevel::Error) {
            return true;
        }
    }
    return false;
}

} // namespace Data
} // namespace DMThermo

