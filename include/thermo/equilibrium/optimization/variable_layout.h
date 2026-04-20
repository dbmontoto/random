/**
 * @file variable_layout.h
 * @brief Shared helpers for slicing optimizer variable vectors into named blocks.
 */

#ifndef THERMO_EQUILIBRIUM_OPTIMIZATION_VARIABLE_LAYOUT_H
#define THERMO_EQUILIBRIUM_OPTIMIZATION_VARIABLE_LAYOUT_H

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Optimization {

struct VarSlice {
    const double* data = nullptr;
    int size = 0;

    bool empty() const { return size == 0; }
    const double& operator[](int i) const { return data[static_cast<size_t>(i)]; }
};

inline VarSlice sliceVector(const std::vector<double>& v, int offset, int size, const char* label) {
    if (offset < 0 || size < 0) {
        throw std::invalid_argument(std::string(label) + ": negative offset/size");
    }
    const size_t off = static_cast<size_t>(offset);
    const size_t len = static_cast<size_t>(size);
    if (off + len > v.size()) {
        throw std::invalid_argument(std::string(label) + ": slice out of range");
    }
    return VarSlice{v.data() + off, size};
}

inline int allocationVariableCount(int nc, int M) {
    if (M <= 1) return 0;
    return nc * (M - 1);
}

struct AllocRhoLayout {
    int nc = 0;
    int M = 0;
    int alloc_offset = 0;
    int alloc_size = 0;
    int rho_offset = 0;
    int rho_size = 0;
    int total_size = 0;

    static AllocRhoLayout make(int nc_, int M_) {
        AllocRhoLayout out;
        out.nc = nc_;
        out.M = M_;
        out.alloc_offset = 0;
        out.alloc_size = allocationVariableCount(nc_, M_);
        out.rho_offset = out.alloc_size;
        out.rho_size = M_;
        out.total_size = out.alloc_size + out.rho_size;
        return out;
    }

    bool matches(const std::vector<double>& v) const { return static_cast<int>(v.size()) == total_size; }

    VarSlice allocSlice(const std::vector<double>& v) const { return sliceVector(v, alloc_offset, alloc_size, "AllocRhoLayout::allocSlice"); }
    VarSlice rhoSlice(const std::vector<double>& v) const { return sliceVector(v, rho_offset, rho_size, "AllocRhoLayout::rhoSlice"); }
};

struct ExtentAllocLayout {
    int nr = 0;
    int nc = 0;
    int M = 0;
    int extent_offset = 0;
    int extent_size = 0;
    int alloc_offset = 0;
    int alloc_size = 0;
    int total_size = 0;

    static ExtentAllocLayout make(int nr_, int nc_, int M_) {
        ExtentAllocLayout out;
        out.nr = nr_;
        out.nc = nc_;
        out.M = M_;
        out.extent_offset = 0;
        out.extent_size = nr_;
        out.alloc_offset = out.extent_size;
        out.alloc_size = allocationVariableCount(nc_, M_);
        out.total_size = out.extent_size + out.alloc_size;
        return out;
    }

    bool matches(const std::vector<double>& v) const { return static_cast<int>(v.size()) == total_size; }

    VarSlice extentSlice(const std::vector<double>& v) const { return sliceVector(v, extent_offset, extent_size, "ExtentAllocLayout::extentSlice"); }
    VarSlice allocSlice(const std::vector<double>& v) const { return sliceVector(v, alloc_offset, alloc_size, "ExtentAllocLayout::allocSlice"); }
};

} // namespace Optimization
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_OPTIMIZATION_VARIABLE_LAYOUT_H

