/**
 * @file PCSAFT.cpp
 * @brief Implementation of PC-SAFT equation of state
 */

#include "thermo/eos/PCSAFT.h"
#include "thermo/legacy/pcsaft/core/component.h"
#include "thermo/legacy/pcsaft/conversions.h"
#include "thermo/core/mixture.h"
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace DMThermo {
namespace PCSaft {

namespace {

pcsaft::AssociationScheme toPcsaftScheme(Core::AssociationScheme scheme) {
    switch (scheme) {
        case Core::AssociationScheme::None: return pcsaft::AssociationScheme::NONE;
        case Core::AssociationScheme::Scheme1: return pcsaft::AssociationScheme::SCHEME_1A;
        case Core::AssociationScheme::Scheme2A: return pcsaft::AssociationScheme::SCHEME_2B; // closest match
        case Core::AssociationScheme::Scheme2B: return pcsaft::AssociationScheme::SCHEME_2B;
        case Core::AssociationScheme::Scheme3B: return pcsaft::AssociationScheme::SCHEME_3B;
        case Core::AssociationScheme::Scheme4C: return pcsaft::AssociationScheme::SCHEME_4C;
        default: return pcsaft::AssociationScheme::NONE;
    }
}

pcsaft::Component toPcsaftComponent(const Core::Component& c) {
    pcsaft::Component out;
    out.setName(c.name());
    out.setCAS(c.cas());

    if (c.hasPCSaftParams()) {
        out.setM(c.m());
        out.setSigma(c.sigma());
        out.setEpsilonK(c.epsilonK());
    }

    if (c.MW() > 0.0 && std::isfinite(c.MW())) {
        out.setMW(c.MW());
    }

    if (c.isAssociating()) {
        const auto& assoc = c.associationParams();
        pcsaft::AssociationParams legacy_assoc(
            toPcsaftScheme(assoc.scheme),
            assoc.epsilon_AB,
            assoc.kappa_AB,
            assoc.num_sites
        );
        out.setAssociationParams(legacy_assoc);
    }

    if (c.hasIdealGasCp()) {
        const auto& cp = c.idealGasCpCoeffs();
        const auto eq = DIPPR::equationTypeFromForm(static_cast<int>(cp.form));
        if (eq) {
            out.setIdealGasCpParams(DIPPR::Coefficients(
                *eq,
                cp.A, cp.B, cp.C, cp.D, cp.E, cp.F, cp.G,
                cp.T_min, cp.T_max
            ));
        }
    }

    return out;
}

pcsaft::Mixture toPcsaftMixture(const Core::Mixture& mix, const std::vector<double>& x) {
    const int nc = mix.numComponents();
    std::vector<pcsaft::Component> components;
    components.reserve(static_cast<std::size_t>(nc));
    for (int i = 0; i < nc; ++i) {
        components.push_back(toPcsaftComponent(mix.component(i)));
    }

    pcsaft::Mixture out(components, x);
    for (int i = 0; i < nc; ++i) {
        for (int j = i + 1; j < nc; ++j) {
            const double kij = mix.kij(i, j);
            if (std::abs(kij) > 1e-15) {
                out.setBinaryParameter(i, j, kij);
            }
        }
    }
    return out;
}

} // namespace

struct PCSaftEOS::Impl {
    pcsaft::Mixture mixture_template;
    mutable std::mutex eos_mutex;
    mutable std::unique_ptr<pcsaft::PCSaftEOS> eos;
};

PCSaftEOS::PCSaftEOS(const Core::Mixture& mixture)
    : impl_(std::make_unique<Impl>())
    , core_mixture_(mixture)
    , nc_(mixture.numComponents())
{
    const std::vector<double> z(
        static_cast<size_t>(mixture.numComponents()),
        1.0 / std::max(1, mixture.numComponents())
    );
    impl_->mixture_template = toPcsaftMixture(mixture, z);

    // If critical properties exist on the Core::Mixture, expose them via the EOS interface.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Tc_.assign(nc_, nan);
    Pc_.assign(nc_, nan);
    omega_.assign(nc_, nan);
    for (int i = 0; i < nc_; ++i) {
        const auto& c = mixture.component(i);
        if (std::isfinite(c.Tc()) && c.Tc() > 0.0) Tc_[i] = c.Tc();
        if (std::isfinite(c.Pc()) && c.Pc() > 0.0) Pc_[i] = c.Pc();
        if (std::isfinite(c.omega())) omega_[i] = c.omega();
    }

    impl_->eos = std::make_unique<pcsaft::PCSaftEOS>(impl_->mixture_template);
    initializeCachedProperties();
}

PCSaftEOS::~PCSaftEOS() = default;

PCSaftEOS::PCSaftEOS(PCSaftEOS&&) noexcept = default;
PCSaftEOS& PCSaftEOS::operator=(PCSaftEOS&&) noexcept = default;

const Core::Mixture* PCSaftEOS::coreMixture() const {
    return (core_mixture_.numComponents() > 0) ? &core_mixture_ : nullptr;
}

void PCSaftEOS::initializeCachedProperties() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto ensureSized = [&](std::vector<double>& v, const char* label) {
        if (v.empty()) {
            v.assign(nc_, nan);
            return;
        }
        if (static_cast<int>(v.size()) != nc_) {
            throw std::invalid_argument(std::string("PCSaftEOS: ") + label + " size mismatch");
        }
    };

    ensureSized(Tc_, "Tc");
    ensureSized(Pc_, "Pc");
    ensureSized(omega_, "omega");
    has_association_ = false;
    cached_x_.resize(nc_, 0.0);

    for (int i = 0; i < nc_; ++i) {
        const auto& comp = impl_->mixture_template.getComponent(i);
        if (comp.hasAssociation()) {
            has_association_ = true;
        }
    }
}

void PCSaftEOS::updateComposition(const std::vector<double>& x) const {
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);

    // Check if composition has changed
    bool same = true;
    for (int i = 0; i < nc_ && same; ++i) {
        if (std::abs(x[i] - cached_x_[i]) > 1e-15) {
            same = false;
        }
    }
    if (same) return;  // No need to rebuild

    // Rebuild the mixture with new composition
    std::vector<pcsaft::Component> components;
    for (int i = 0; i < nc_; ++i) {
        components.push_back(impl_->mixture_template.getComponent(i));
    }

    // Create new mixture with updated composition
    pcsaft::Mixture new_mix(components, x);

    // Copy binary parameters
    for (int i = 0; i < nc_; ++i) {
        for (int j = i + 1; j < nc_; ++j) {
            double kij = impl_->mixture_template.getBinaryParameter(i, j);
            if (std::abs(kij) > 1e-15) {
                new_mix.setBinaryParameter(i, j, kij);
            }
        }
    }

    // Rebuild EOS
    impl_->eos = std::make_unique<pcsaft::PCSaftEOS>(new_mix);

    // Cache new composition
    cached_x_ = x;
}

EOSResult PCSaftEOS::calculate(
    double T,
    double rho,
    const std::vector<double>& x,
    const DerivativeSpec& deriv_spec) const
{
    EOSResult result;

    // Update composition if needed
    updateComposition(x);

    // Call legacy calculate (thread-safe via mutex in updateComposition)
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    bool calc_derivs = deriv_spec.temperature || deriv_spec.density || deriv_spec.second_order;
    pcsaft::PCSaftResults legacy = impl_->eos->calculate(T, rho, calc_derivs);

    // Map results
    result.a_residual = legacy.a_res;
    result.pressure = legacy.pressure;
    result.compressibility = legacy.Z;

    // Fugacity coefficients
    result.phi = legacy.fugacity_coeff;
    result.ln_phi.resize(nc_);
    for (int i = 0; i < nc_; ++i) {
        result.ln_phi[i] = std::log(legacy.fugacity_coeff[i]);
    }

    // Derivatives if requested
    if (deriv_spec.temperature) {
        result.da_dT = legacy.dadt;
    }
    if (deriv_spec.density) {
        result.da_drho = legacy.dadrho;
    }
    if (deriv_spec.second_order) {
        result.d2a_dT2 = legacy.d2adt2;
        result.d2a_drho2 = legacy.d2adrho2;
        result.d2a_dTdrho = legacy.d2adtdrho;
    }

    // Chemical potentials
    if (deriv_spec.composition) {
        result.chemical_potential = legacy.mu_res;
    }

    result.success = true;
    return result;
}

double PCSaftEOS::pressure(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    return impl_->eos->calculatePressure(T, rho);
}

std::vector<double> PCSaftEOS::fugacityCoefficients(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    return impl_->eos->calculateFugacityCoefficients(T, rho);
}

double PCSaftEOS::compressibility(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    return impl_->eos->calculateZ(T, rho);
}

double PCSaftEOS::residualHelmholtz(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    pcsaft::PCSaftResults res = impl_->eos->calculate(T, rho, false);
    return res.a_res;
}

double PCSaftEOS::dadt(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    pcsaft::PCSaftResults res = impl_->eos->calculate(T, rho, true);
    return res.dadt;
}

double PCSaftEOS::dadrho(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    pcsaft::PCSaftResults res = impl_->eos->calculate(T, rho, true);
    return res.dadrho;
}

double PCSaftEOS::dPdrho(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    // dP/drho = RT * (1 + 2*rho*da/drho + rho^2*d^2a/drho^2)
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    pcsaft::PCSaftResults res = impl_->eos->calculate(T, rho, true);

    double R = Constants::GAS_CONSTANT;
    return R * T * (1.0 + 2.0 * rho * res.dadrho + rho * rho * res.d2adrho2);
}

double PCSaftEOS::criticalTemperature(int i) const {
    if (i < 0 || i >= nc_) {
        throw std::out_of_range("Component index out of range");
    }
    return Tc_[i];
}

double PCSaftEOS::criticalPressure(int i) const {
    if (i < 0 || i >= nc_) {
        throw std::out_of_range("Component index out of range");
    }
    return Pc_[i];
}

double PCSaftEOS::acentricFactor(int i) const {
    if (i < 0 || i >= nc_) {
        throw std::out_of_range("Component index out of range");
    }
    return omega_[i];
}

double PCSaftEOS::packingFraction(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    return impl_->eos->getPackingFraction(T, rho);
}

double PCSaftEOS::hardChainContribution(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    return impl_->eos->calc_a_hc(T, rho);
}

double PCSaftEOS::dispersionContribution(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    return impl_->eos->calc_a_disp(T, rho);
}

double PCSaftEOS::associationContribution(
    double T,
    double rho,
    const std::vector<double>& x) const
{
    updateComposition(x);
    std::lock_guard<std::mutex> lock(impl_->eos_mutex);
    pcsaft::PCSaftResults res = impl_->eos->calculate(T, rho, false);
    return res.a_assoc;
}

// ============================================================================
// Factory Functions
// ============================================================================

} // namespace PCSaft
} // namespace DMThermo
