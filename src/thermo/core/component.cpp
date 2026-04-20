/**
 * @file component.cpp
 * @brief Implementation of DMThermo::Core::Component
 */

#include "thermo/core/component.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

namespace DMThermo {
namespace Core {

// =============================================================================
// Component Constructors
// =============================================================================

Component::Component(
    const std::string& name)
    : name_(name)
    , cas_("")
    , m_(0.0)
    , sigma_(0.0)
    , epsilon_k_(0.0)
    , mw_(0.0)
    , tc_(0.0)
    , pc_(0.0)
    , zc_(std::numeric_limits<double>::quiet_NaN())
    , omega_(0.0)
    , vtpr_c_(std::nullopt)
    , unifac_groups_()
    , assoc_()
    , poly_()
    , ideal_gas_cp_(std::nullopt)
{
}

Component::Component(
    const std::string& name,
    double m,
    double sigma,
    double epsilon_k)
    : name_(name)
    , cas_("")
    , m_(m)
    , sigma_(sigma)
    , epsilon_k_(epsilon_k)
    , mw_(0.0)
    , tc_(0.0)
    , pc_(0.0)
    , zc_(std::numeric_limits<double>::quiet_NaN())
    , omega_(0.0)
    , vtpr_c_(std::nullopt)
    , unifac_groups_()
    , assoc_()
    , poly_()
    , ideal_gas_cp_(std::nullopt)
{
}

Component::Component(
    const std::string& name,
    double m,
    double sigma,
    double epsilon_k,
    const AssociationParams& assoc)
    : name_(name)
    , cas_("")
    , m_(m)
    , sigma_(sigma)
    , epsilon_k_(epsilon_k)
    , mw_(0.0)
    , tc_(0.0)
    , pc_(0.0)
    , zc_(std::numeric_limits<double>::quiet_NaN())
    , omega_(0.0)
    , vtpr_c_(std::nullopt)
    , unifac_groups_()
    , assoc_(assoc)
    , poly_()
    , ideal_gas_cp_(std::nullopt)
{
}

Component::Component(
    const std::string& name,
    double sigma,
    double epsilon_k,
    const PolymerParams& poly,
    double MW)
    : name_(name)
    , cas_("")
    , m_(poly.segmentNumber(MW))
    , sigma_(sigma)
    , epsilon_k_(epsilon_k)
    , mw_(MW)
    , tc_(0.0)
    , pc_(0.0)
    , zc_(std::numeric_limits<double>::quiet_NaN())
    , omega_(0.0)
    , vtpr_c_(std::nullopt)
    , unifac_groups_()
    , assoc_()
    , poly_(poly)
    , ideal_gas_cp_(std::nullopt)
{
    poly_.is_polymer = true;
}

// =============================================================================
// Builder Methods
// =============================================================================

Component Component::withMW(double MW) const {
    if (!poly_.is_polymer) {
        // For non-polymers, just set MW without changing m
        Component result = *this;
        result.mw_ = MW;
        return result;
    }

    // For polymers, recalculate segment number
    Component result = *this;
    result.mw_ = MW;
    result.m_ = poly_.segmentNumber(MW);
    return result;
}

Component Component::withCAS(const std::string& cas) const {
    Component result = *this;
    result.cas_ = cas;
    return result;
}

Component Component::withCritical(double Tc, double Pc, double omega) const {
    Component result = *this;
    result.tc_ = Tc;
    result.pc_ = Pc;
    result.omega_ = omega;
    return result;
}

Component Component::withCriticalCompressibility(double Zc) const {
    Component result = *this;
    result.zc_ = Zc;
    return result;
}

Component Component::withVtprVolumeTranslation(const VtprVolumeTranslationPoly& c) const {
    Component result = *this;
    result.vtpr_c_ = c;
    return result;
}

Component Component::withUnifacSubgroups(const std::vector<UnifacSubgroupCount>& groups) const {
    Component result = *this;
    result.unifac_groups_ = groups;
    return result;
}

Component Component::withIdealGasCp(const IdealGasCpCorrelation& coeffs) const {
    Component result = *this;
    result.ideal_gas_cp_ = coeffs;
    return result;
}

// =============================================================================
// Utility Methods
// =============================================================================

void Component::print() const {
    std::cout << "Component: " << name_ << "\n";
    if (!cas_.empty()) {
        std::cout << "  CAS: " << cas_ << "\n";
    }
    std::cout << "  PC-SAFT Parameters:\n";
    std::cout << "    m = " << std::fixed << std::setprecision(4) << m_ << "\n";
    std::cout << "    sigma = " << sigma_ << " Å\n";
    std::cout << "    epsilon/k = " << epsilon_k_ << " K\n";

    if (mw_ > 0) {
        std::cout << "  MW = " << mw_ << " g/mol\n";
    }

    if (tc_ > 0) {
        std::cout << "  Critical Properties:\n";
        std::cout << "    Tc = " << tc_ << " K\n";
        std::cout << "    Pc = " << pc_ / 1e5 << " bar\n";
        if (std::isfinite(zc_)) {
            std::cout << "    Zc = " << zc_ << "\n";
        }
        std::cout << "    omega = " << omega_ << "\n";
    }

    if (assoc_.isAssociating()) {
        std::cout << "  Association:\n";
        std::cout << "    epsilon_AB/k = " << assoc_.epsilon_AB << " K\n";
        std::cout << "    kappa_AB = " << assoc_.kappa_AB << "\n";
        std::cout << "    sites = " << assoc_.num_sites << "\n";
    }

    if (poly_.is_polymer) {
        std::cout << "  Polymer:\n";
        std::cout << "    m/MW = " << poly_.m_per_mw << " (g/mol)^-1\n";
        if (poly_.mw_monomer > 0) {
            std::cout << "    MW_monomer = " << poly_.mw_monomer << " g/mol\n";
        }
    }

    if (ideal_gas_cp_ && ideal_gas_cp_->isValid()) {
        std::cout << "  Ideal Gas Cp Correlation:\n";
        std::cout << "    eq_form: " << ideal_gas_cp_->eq_form << "\n";
        std::cout << "    T range: " << ideal_gas_cp_->Tmin << " - " << ideal_gas_cp_->Tmax << " K\n";
    }
}

bool Component::isValid() const {
    // Basic parameter checks
    if (name_.empty()) return false;

    // PC-SAFT parameters are optional for non-SAFT models; treat "all zero" as unset.
    const bool has_pcsaft = (m_ != 0.0) || (sigma_ != 0.0) || (epsilon_k_ != 0.0);
    if (has_pcsaft) {
        if (m_ <= 0) return false;
        if (sigma_ <= 0) return false;
        if (epsilon_k_ <= 0) return false;
    }

    // Association parameter checks
    if (assoc_.isAssociating()) {
        if (assoc_.epsilon_AB < 0) return false;
        if (assoc_.kappa_AB < 0) return false;
        if (assoc_.num_sites <= 0) return false;
    }

    // Polymer parameter checks
    if (poly_.is_polymer) {
        if (poly_.m_per_mw <= 0) return false;
    }

    // Critical property checks (if set)
    if (tc_ > 0) {
        if (!std::isfinite(tc_) || !std::isfinite(pc_) || pc_ <= 0) return false;
        if (std::isfinite(zc_)) {
            if (zc_ <= 0.0 || zc_ > 1.0) return false;
        }
        if (!std::isfinite(omega_)) return false;
        if (omega_ < -1 || omega_ > 2) return false;  // Reasonable acentric factor range
    }

    if (vtpr_c_.has_value()) {
        if (!std::isfinite(vtpr_c_->c0) || !std::isfinite(vtpr_c_->c1) || !std::isfinite(vtpr_c_->c2) || !std::isfinite(vtpr_c_->c3)) {
            return false;
        }
        if (vtpr_c_->Tmin.has_value() && !std::isfinite(*vtpr_c_->Tmin)) return false;
        if (vtpr_c_->Tmax.has_value() && !std::isfinite(*vtpr_c_->Tmax)) return false;
        if (vtpr_c_->Tmin && vtpr_c_->Tmax && (*vtpr_c_->Tmax < *vtpr_c_->Tmin)) return false;
    }

    for (const auto& g : unifac_groups_) {
        if (g.subgroup_id <= 0) return false;
        if (g.count <= 0) return false;
    }

    return true;
}

} // namespace Core
} // namespace DMThermo
