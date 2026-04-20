#include "thermo/pcsaft/eos/PCSAFT.h"
#include "thermo/pcsaft/core/constants.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/dual.hpp>

// Use dual2nd for true second-order automatic differentiation
using dual2nd = autodiff::HigherOrderDual<2, double>;

namespace pcsaft {

// Dispersion integrals coefficients (from Gross & Sadowski 2001, Table 3)
constexpr double a0[7] = { 0.910563145, 0.636128145, 2.686134789, -26.54736249,
                            97.75920878, -159.5915409, 91.29777408 };
constexpr double a1[7] = { -0.308401692, 0.186053116, -2.797738490, 21.41953503,
                            -65.25588533, 83.31868048, -33.74692293 };
constexpr double a2[7] = { -0.090614835, 0.452784281, 0.596270073, -1.72481385,
                            -4.13021158, 13.77663187, -8.67284703 };

constexpr double b0[7] = { 0.724094694, 2.238279186, -4.002584949, -21.00357682,
                            26.85564116, 206.5513384, -355.6023561 };
constexpr double b1[7] = { -0.575549808, 0.699509552, 3.892567339, -17.21547165,
                            192.6722645, -161.8264616, -165.2076935 };
constexpr double b2[7] = { 0.097688312, -0.255757498, -9.155856153, 20.64207597,
                            -38.80443005, 93.62677877, -29.66690559 };

// Gross & Vrabec (2006) universal constants for improved numerical stability
// These provide better convergence for short-chain molecules
constexpr double a_universal[7][3] = {
    { 0.9105631445, -0.3084016918, -0.0906148351 },
    { 0.6361281449,  0.1860531159,  0.4527842806 },
    { 2.6861347891, -2.5030047259,  0.5962700728 },
    {-26.547362491,  21.419793629, -1.7241829131 },
    { 97.759208784, -65.255885330, -4.1302112531 },
    {-159.59154087,  83.318680481,  13.776631870 },
    { 91.297774084, -33.746922930, -8.6728470368 }
};

constexpr double b_universal[7][3] = {
    { 0.7240946941, -0.5755498075,  0.0976883116 },
    { 2.2382791861,  0.6995095521, -0.2557574982 },
    {-4.0025849485,  3.8925673390, -9.1558561530 },
    {-21.003576815, -17.215471648,  20.642075974 },
    { 26.855641363,  192.67226447, -38.804430052 },
    { 206.55133841, -161.82646165,  93.626774077 },
    {-355.60235612, -165.20769346, -29.666905585 }
};

// ============================================================================
// Constructor
// ============================================================================

PCSaftEOS::PCSaftEOS(const Mixture& mixture)
    : mixture_(mixture), nc_(mixture.getNumComponents()) {

    // Extract pure component parameters
    m_.resize(nc_);
    sigma_.resize(nc_);
    epsilon_.resize(nc_);

    for (int i = 0; i < nc_; ++i) {
        const Component& comp = mixture_.getComponent(i);
        m_[i] = comp.getM();
        sigma_[i] = comp.getSigma();
        epsilon_[i] = comp.getEpsilonK();
    }
}

// ============================================================================
// Main calculation method
// ============================================================================

PCSaftResults PCSaftEOS::calculate(double T, double rho, bool calc_derivatives) {
    PCSaftResults res;

    res.temperature = T;
    res.density = rho;

    // Calculate each contribution
    res.a_hc = calc_a_hc(T, rho);
    res.a_disp = calc_a_disp(T, rho);
    res.a_assoc = mixture_.hasAssociation() ? calc_a_assoc(T, rho) : 0.0;
    res.a_res = res.a_hc + res.a_disp + res.a_assoc;

    // Calculate first derivatives
    res.dadrho = calc_dadrho(T, rho);
    res.dadt = calc_dadt(T, rho);

    // Calculate pressure and compressibility factor
    // Z = 1 + ρ(∂a_res/∂ρ)_T
    res.Z = 1.0 + rho * res.dadrho;
    res.pressure = res.Z * rho * R * T;

    // Calculate chemical potentials and fugacity coefficients
    res.mu_res = calc_mu_res(T, rho);
    res.fugacity_coeff.resize(nc_);
    for (int i = 0; i < nc_; ++i) {
        res.fugacity_coeff[i] = std::exp(res.mu_res[i]);
    }

    // Calculate second derivatives if requested
    if (calc_derivatives) {
        res.d2adrho2 = calc_d2adrho2(T, rho);
        res.d2adt2 = calc_d2adt2(T, rho);
        res.d2adtdrho = calc_d2adtdrho(T, rho);
    }

    return res;
}

double PCSaftEOS::calculatePressure(double T, double rho) const {
    double dadrho = calc_dadrho(T, rho);
    double Z = 1.0 + rho * dadrho;
    return Z * rho * R * T;
}

double PCSaftEOS::calculateZ(double T, double rho) const {
    double dadrho = calc_dadrho(T, rho);
    return 1.0 + rho * dadrho;
}

std::vector<double> PCSaftEOS::calculateFugacityCoefficients(double T, double rho) const {
    std::vector<double> mu_res = calc_mu_res(T, rho);
    std::vector<double> fugacity_coeff(nc_);

    for (int i = 0; i < nc_; ++i) {
        fugacity_coeff[i] = std::exp(mu_res[i]);
    }

    return fugacity_coeff;
}

// ============================================================================
// Hard-chain contribution
// ============================================================================

double PCSaftEOS::calc_a_hs(double eta) const {
    // Carnahan-Starling equation for hard spheres (simplified, single-component)
    // NOTE: This is a simplified version. The full mixture formula should be:
    // ã^{hs} = (1/ζ_0)[3ζ_1ζ_2/(1-ζ_3) + ζ_2³/(ζ_3(1-ζ_3)²) + (ζ_2³/ζ_3² - ζ_0) ln(1-ζ_3)]
    // This simplified version is only accurate for pure components!

    // Smooth clamping to avoid discontinuities in derivatives (autodiff-safe)
    // Allow small negative eta (will give ~0) and clamp upper bound smoothly
    constexpr double eta_min = 1e-15;  // Avoid exactly zero
    constexpr double eta_max = 0.999;  // Stay away from singularity at eta=1

    double eta_safe = smoothClamp(eta, eta_min, eta_max);

    double term1 = 4.0 * eta_safe - 3.0 * eta_safe * eta_safe;
    double term2 = (1.0 - eta_safe) * (1.0 - eta_safe);

    return term1 / term2;
}

double PCSaftEOS::calc_a_hs_mixture(const std::vector<double>& zeta) const {
    // Full hard-sphere Helmholtz energy for mixtures
    // From Gross & Sadowski (2001), Equation 8 (Mansoori et al. formula):    // ã^{hs} = (1/ζ_0)[3ζ_1ζ_2/(1-ζ_3) + ζ_2³/(ζ_3(1-ζ_3)²) + (ζ_2³/ζ_3² - ζ_0) ln(1-ζ_3)]

    // Smooth clamping to avoid discontinuities (autodiff-safe)
    constexpr double zeta_min = 1e-15;
    constexpr double zeta3_max = 0.999;  // Stay away from singularity at zeta3=1

    double zeta0_safe = smoothClampLower(zeta[0], zeta_min);
    double zeta3_safe = smoothClamp(zeta[3], zeta_min, zeta3_max);

    double one_minus_zeta3 = 1.0 - zeta3_safe;
    double one_minus_zeta3_sq = one_minus_zeta3 * one_minus_zeta3;

    // Use safe zeta3 in divisions to avoid singularities
    double term1 = 3.0 * zeta[1] * zeta[2] / one_minus_zeta3;
    double term2 = zeta[2] * zeta[2] * zeta[2] / (zeta3_safe * one_minus_zeta3_sq);
    double term3 = (zeta[2] * zeta[2] * zeta[2] / (zeta3_safe * zeta3_safe) - zeta[0]) * std::log(one_minus_zeta3);

    return (term1 + term2 + term3) / zeta0_safe;
}

double PCSaftEOS::calc_a_hc(double T, double rho) const {
    const std::vector<double>& x = mixture_.getMoleFractions();
    std::vector<double> d = calc_d(T);

    // Calculate packing fraction η
    double eta = calc_eta(T, rho);

    // Smooth clamping (autodiff-safe)
    constexpr double eta_min = 1e-15;
    constexpr double eta_max = 0.999;
    eta = smoothClamp(eta, eta_min, eta_max);

    // Calculate mean segment diameter
    double d_bar = 0.0;
    for (int i = 0; i < nc_; ++i) {
        d_bar += x[i] * m_[i] * d[i];
    }

    // Calculate ζ0, ζ1, ζ2, ζ3
    std::vector<double> zeta(4, 0.0);
    for (int i = 0; i < nc_; ++i) {
        double d_i = d[i];
        zeta[0] += x[i] * m_[i];
        zeta[1] += x[i] * m_[i] * d_i;
        zeta[2] += x[i] * m_[i] * d_i * d_i;
        zeta[3] += x[i] * m_[i] * d_i * d_i * d_i;
    }

    // Convert density from mol/m³ to molecules/Å³
    double rho_reduced = rho * N_A * 1e-30;

    // Multiply by π/6 * ρ_reduced
    double factor = PI / 6.0 * rho_reduced;
    for (int k = 0; k < 4; ++k) {
        zeta[k] *= factor;
    }

    // Hard-sphere contribution (use full mixture formula)
    double a_hs = calc_a_hs_mixture(zeta);

    // Radial distribution function at contact for reference
    double g_hs_ref = calc_g_hs(zeta[3]);

    // Hard-chain contribution
    // a_hc = m_bar * a_hs - Σ_i x_i (m_i - 1) ln(g_ii)
    double m_bar = calc_m_bar();

    double sum_chain = 0.0;
    for (int i = 0; i < nc_; ++i) {
        if (m_[i] > 1.0) {
            // Calculate component-specific g_ii using proper formula
            double g_ii = calc_g_ii(d[i], zeta[2], zeta[3]);
            sum_chain += x[i] * (m_[i] - 1.0) * std::log(g_ii);
        }
    }

    double a_hc = m_bar * a_hs - sum_chain;

    return a_hc;
}

double PCSaftEOS::calc_g_hs(double eta) const {
    // Simplified radial distribution function at contact
    // Smooth clamping (autodiff-safe)
    constexpr double eta_min = 1e-15;
    constexpr double eta_max = 0.999;
    double eta_safe = smoothClamp(eta, eta_min, eta_max);

    double one_minus_eta = 1.0 - eta_safe;
    return (1.0 - 0.5 * eta_safe) / (one_minus_eta * one_minus_eta * one_minus_eta);
}

double PCSaftEOS::calc_g_ii(double d_i, double zeta_2, double zeta_3) const {
    // Component-specific radial distribution function at contact
    // From Gross & Sadowski (2001), Equation 11 for i=i case
    // g_ii^hs = 1/(1-ζ_3) + (d_i/2) · 3ζ_2/(1-ζ_3)² + (d_i²/4) · 2ζ_2²/(1-ζ_3)³

    if (zeta_3 <= 0.0) return 1.0;
    if (zeta_3 >= 1.0) return 1e10;

    double one_minus_zeta3 = 1.0 - zeta_3;
    double one_minus_zeta3_sq = one_minus_zeta3 * one_minus_zeta3;
    double one_minus_zeta3_cube = one_minus_zeta3_sq * one_minus_zeta3;

    double term1 = 1.0 / one_minus_zeta3;
    double term2 = (d_i / 2.0) * 3.0 * zeta_2 / one_minus_zeta3_sq;
    double term3 = (d_i * d_i / 4.0) * 2.0 * zeta_2 * zeta_2 / one_minus_zeta3_cube;

    return term1 + term2 + term3;
}

// ============================================================================
// Dispersion contribution
// ============================================================================

double PCSaftEOS::calc_a_disp(double T, double rho) const {
    const std::vector<double>& x = mixture_.getMoleFractions();
    std::vector<double> d = calc_d(T);

    double eta = calc_eta(T, rho);

    // Smooth clamping (autodiff-safe) - only for numerical boundary protection
    constexpr double eta_min = 1e-15;
    constexpr double eta_max = 0.999;
    eta = smoothClamp(eta, eta_min, eta_max);

    double m_bar = calc_m_bar();

    // Calculate m²εσ³ and m²ε²σ³
    // Σ_i Σ_j x_i x_j m_i m_j ε_ij σ_ij³
    double m2es3 = 0.0;
    double m2e2s3 = 0.0;

    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            // Combining rules (Lorentz-Berthelot)
            double sigma_ij = 0.5 * (sigma_[i] + sigma_[j]);
            double epsilon_ij = std::sqrt(epsilon_[i] * epsilon_[j]);

            // Apply binary interaction parameter kij
            double kij = mixture_.getBinaryParameter(i, j);
            epsilon_ij *= (1.0 - kij);

            double sigma_ij3 = sigma_ij * sigma_ij * sigma_ij;

            m2es3 += x[i] * x[j] * m_[i] * m_[j] * epsilon_ij * sigma_ij3;
            m2e2s3 += x[i] * x[j] * m_[i] * m_[j] * epsilon_ij * epsilon_ij * sigma_ij3;
        }
    }

    // Convert density from mol/m³ to molecules/Å³
    double rho_reduced = rho * N_A * 1e-30;

    // Calculate dispersion integrals I1 and I2 from power series
    double I1, I2;
    calc_dispersion_integrals(eta, I1, I2);

    // Compressibility factor contribution from hard-chain term
    // C1 = (1 + Z_hc + rho * dZ_hc/drho)^-1
    // Standard PC-SAFT expression from Gross & Sadowski (2001)
    double eta2 = eta * eta;
    double eta3 = eta * eta * eta;
    double eta4 = eta2 * eta2;
    double one_minus_eta = 1.0 - eta;
    double one_minus_eta4 = one_minus_eta * one_minus_eta * one_minus_eta * one_minus_eta;
    double two_minus_eta = 2.0 - eta;

    double C1 = 1.0 / (1.0 + m_bar * (8.0 * eta - 2.0 * eta2) / one_minus_eta4
                       + (1.0 - m_bar) * (20.0 * eta - 27.0 * eta2 + 12.0 * eta3 - 2.0 * eta4)
                         / (one_minus_eta * one_minus_eta * two_minus_eta * two_minus_eta));

    // First-order perturbation term (mean-attractive energy)
    double a_1 = -2.0 * PI * rho_reduced * I1 * m2es3 / T;

    // Second-order perturbation term (fluctuation term)
    double a_2 = -PI * rho_reduced * m_bar * C1 * I2 * m2e2s3 / (T * T);

    return a_1 + a_2;
}

void PCSaftEOS::calc_dispersion_integrals(double eta, double& I1, double& I2) const {
    double m_bar = calc_m_bar();

    I1 = 0.0;
    I2 = 0.0;

    // Standard PC-SAFT dispersion integrals
    // Gross & Sadowski (2001) - Equation 12
    // Power series: I = Σ(i=0 to 6) ai * η^i
    // where ai = ai0 + (m-1)/m * ai1 + (m-1)(m-2)/m² * ai2

    for (int i = 0; i < 7; ++i) {
        double eta_pow = std::pow(eta, i);

        // I1 coefficient
        double I1_coeff = a0[i] + (m_bar - 1.0) / m_bar * a1[i] +
                         (m_bar - 1.0) * (m_bar - 2.0) / (m_bar * m_bar) * a2[i];
        I1 += I1_coeff * eta_pow;

        // I2 coefficient
        double I2_coeff = b0[i] + (m_bar - 1.0) / m_bar * b1[i] +
                         (m_bar - 1.0) * (m_bar - 2.0) / (m_bar * m_bar) * b2[i];
        I2 += I2_coeff * eta_pow;
    }
}

// ============================================================================
// Association contribution
// ============================================================================

double PCSaftEOS::calc_a_assoc(double T, double rho) const {
    if (!mixture_.hasAssociation()) {
        return 0.0;
    }

    const std::vector<double>& x = mixture_.getMoleFractions();

    // Solve for fraction of non-bonded sites X_A
    std::vector<double> X_A = solve_association(T, rho);

    // Calculate association contribution
    // a_assoc = Σ_i x_i Σ_A (ln(X_Ai) - X_Ai/2 + 1/2)
    double a_assoc = 0.0;

    int site_idx = 0;
    for (int i = 0; i < nc_; ++i) {
        const AssociationParams& assoc = mixture_.getComponent(i).getAssociationParams();

        if (assoc.hasAssociation()) {
            int num_sites = assoc.num_sites;
            for (int A = 0; A < num_sites; ++A) {
                double X_Ai = X_A[site_idx++];
                a_assoc += x[i] * (std::log(X_Ai) - X_Ai / 2.0 + 0.5);
            }
        }
    }

    return a_assoc;
}

std::vector<double> PCSaftEOS::solve_association(double T, double rho) const {
    // Count total number of association sites
    int total_sites = 0;
    for (int i = 0; i < nc_; ++i) {
        total_sites += mixture_.getComponent(i).getAssociationParams().num_sites;
    }

    if (total_sites == 0) {
        return std::vector<double>();
    }

    // Convert density from mol/m³ to molecules/Å³
    double rho_reduced = rho * N_A * 1e-30;

    // Initialize X_A = 1 (all sites non-bonded initially)
    std::vector<double> X_A(total_sites, 1.0);
    std::vector<double> X_A_new(total_sites);

    const std::vector<double>& x = mixture_.getMoleFractions();
    std::vector<double> d = calc_d(T);

    // Calculate association strengths Δ^AB for all site pairs
    std::vector<std::vector<double>> Delta(nc_, std::vector<double>(nc_, 0.0));
    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            if (mixture_.getComponent(i).hasAssociation() &&
                mixture_.getComponent(j).hasAssociation()) {
                Delta[i][j] = calc_delta_assoc(i, j, T, d, rho);
            }
        }
    }

    // Successive substitution to solve X_A equations
    const int max_iter = 100;
    const double tol = 1e-10;

    for (int iter = 0; iter < max_iter; ++iter) {
        int site_idx = 0;
        for (int i = 0; i < nc_; ++i) {
            int num_sites_i = mixture_.getComponent(i).getAssociationParams().num_sites;

            for (int A = 0; A < num_sites_i; ++A) {
                // X_Ai = 1 / (1 + ρ Σ_j x_j Σ_B X_Bj Δ^AB_ij)
                double sum = 0.0;

                int site_idx_j = 0;
                for (int j = 0; j < nc_; ++j) {
                    int num_sites_j = mixture_.getComponent(j).getAssociationParams().num_sites;

                    for (int B = 0; B < num_sites_j; ++B) {
                        sum += x[j] * X_A[site_idx_j++] * Delta[i][j];
                    }
                }

                X_A_new[site_idx++] = 1.0 / (1.0 + rho_reduced * sum);
            }
        }

        // Check convergence
        double max_diff = 0.0;
        for (int k = 0; k < total_sites; ++k) {
            max_diff = std::max(max_diff, std::abs(X_A_new[k] - X_A[k]));
        }

        X_A = X_A_new;

        if (max_diff < tol) {
            break;
        }
    }

    return X_A;
}

double PCSaftEOS::calc_delta_assoc(int i, int j, double T, const std::vector<double>& d, double rho) const {
    // Association strength Δ^AB = g(d_ij) * [exp(ε^AB/kT) - 1] * κ^AB * σ_ij³

    const AssociationParams& assoc_i = mixture_.getComponent(i).getAssociationParams();
    const AssociationParams& assoc_j = mixture_.getComponent(j).getAssociationParams();

    if (!assoc_i.hasAssociation() || !assoc_j.hasAssociation()) {
        return 0.0;
    }

    // Combining rules for cross-association
    double epsilon_AB = 0.5 * (assoc_i.epsilon_AB + assoc_j.epsilon_AB);
    double kappa_AB = std::sqrt(assoc_i.kappa_AB * assoc_j.kappa_AB);
    double sigma_ij = 0.5 * (sigma_[i] + sigma_[j]);

    double d_ij = 0.5 * (d[i] + d[j]);

    // Calculate packing fraction using full mixture composition and actual density
    double eta = calc_eta(T, rho);
    double g_hs = calc_g_hs(eta);

    double Delta = g_hs * (std::exp(epsilon_AB / T) - 1.0) * kappa_AB * sigma_ij * sigma_ij * sigma_ij;

    return Delta;
}

// ============================================================================
// Derivative calculations (numerical for now, can be analytical later)
// ============================================================================

double PCSaftEOS::calc_dadrho(double T, double rho) const {
    // AUTODIFF: ∂a/∂ρ at constant T using automatic differentiation
    // This provides machine-precision derivatives and is faster than numerical
    using namespace autodiff;

    const std::vector<double>& x = mixture_.getMoleFractions();
    std::vector<dual> x_dual(nc_);
    for (int i = 0; i < nc_; ++i) {
        x_dual[i] = x[i];
    }

    auto f = [&](dual rho_d) -> dual {
        dual result = calc_a_hc_templated(x_dual, dual(T), rho_d) +
                     calc_a_disp_templated(x_dual, dual(T), rho_d);
        if (mixture_.hasAssociation()) {
            result += calc_a_assoc_templated(x_dual, dual(T), rho_d);
        }
        return result;
    };

    dual rho_d = rho;
    return derivative(f, wrt(rho_d), at(rho_d));
}

double PCSaftEOS::calc_dadt(double T, double rho) const {
    // AUTODIFF: ∂a/∂T at constant ρ using automatic differentiation
    // This provides machine-precision derivatives and is faster than numerical
    using namespace autodiff;

    const std::vector<double>& x = mixture_.getMoleFractions();
    std::vector<dual> x_dual(nc_);
    for (int i = 0; i < nc_; ++i) {
        x_dual[i] = x[i];
    }

    auto f = [&](dual T_d) -> dual {
        dual result = calc_a_hc_templated(x_dual, T_d, dual(rho)) +
                     calc_a_disp_templated(x_dual, T_d, dual(rho));
        if (mixture_.hasAssociation()) {
            result += calc_a_assoc_templated(x_dual, T_d, dual(rho));
        }
        return result;
    };

    dual T_d = T;
    return derivative(f, wrt(T_d), at(T_d));
}

double PCSaftEOS::calc_d2adrho2(double T, double rho) const {
    // HYBRID AUTODIFF: Numerical derivative of exact (autodiff) first derivative
    // This is more reliable than nested autodiff and much better than double numerical
    // Error is O(h²) where the first derivative is exact to machine precision
    const double h = std::max(std::abs(rho) * 1e-5, 1e-8);
    double da_plus = calc_dadrho(T, rho + h);
    double da_minus = calc_dadrho(T, rho - h);
    return (da_plus - da_minus) / (2.0 * h);
}

double PCSaftEOS::calc_d2adt2(double T, double rho) const {
    // HYBRID AUTODIFF: Numerical derivative of exact (autodiff) first derivative
    // This is more reliable than nested autodiff and much better than double numerical
    const double h = std::max(std::abs(T) * 1e-5, 1e-6);
    double da_plus = calc_dadt(T + h, rho);
    double da_minus = calc_dadt(T - h, rho);
    return (da_plus - da_minus) / (2.0 * h);
}

double PCSaftEOS::calc_d2adtdrho(double T, double rho) const {
    // HYBRID AUTODIFF: Compute d(da/dT)/drho using numerical diff of exact first derivative
    // This gives the mixed partial ∂²a/(∂T∂ρ)
    const double h = std::max(std::abs(rho) * 1e-5, 1e-8);
    double da_plus = calc_dadt(T, rho + h);
    double da_minus = calc_dadt(T, rho - h);
    return (da_plus - da_minus) / (2.0 * h);
}

std::vector<double> PCSaftEOS::calc_mu_res(double T, double rho) const {
    // Residual chemical potential for fugacity coefficient calculation
    // ln(φ_i) = μ_i^res/(RT) - ln(Z)
    //
    // For PC-SAFT mixtures:
    // μ_i^res/(RT) = a^res + (∂a^res/∂x_i) - Σ_j x_j(∂a^res/∂x_j)
    //
    // This implementation uses AUTODIFF for hard-chain and dispersion (machine precision)
    // and numerical for association (if present)

    std::vector<double> mu_res(nc_);
    const std::vector<double>& x = mixture_.getMoleFractions();

    // Calculate baseline a_res
    double a_hc = calc_a_hc(T, rho);
    double a_disp = calc_a_disp(T, rho);
    double a_assoc = 0.0;
    if (mixture_.hasAssociation()) {
        a_assoc = calc_a_assoc(T, rho);
    }
    double a_res = a_hc + a_disp + a_assoc;


    // For pure components: simplified
    if (nc_ == 1) {
        double Z = calculateZ(T, rho);
        mu_res[0] = a_res - std::log(Z);
        return mu_res;
    }

    // For mixtures: use autodiff for both hard-chain and dispersion
    // (machine precision accurate, replaces broken analytical formulas)
    // Calculate ∂a^res/∂x_i for each component
    std::vector<double> da_hc_dx = calc_da_hc_dx_autodiff(T, rho);    // AUTODIFF: machine precision
    std::vector<double> da_disp_dx = calc_da_disp_dx_autodiff(T, rho); // AUTODIFF: machine precision

    std::vector<double> da_dx(nc_);
    for (int i = 0; i < nc_; ++i) {
        da_dx[i] = da_hc_dx[i] + da_disp_dx[i];
    }

    // Add association derivatives if present
    if (mixture_.hasAssociation()) {
        // Use autodiff for association (machine precision accurate)
        std::vector<double> da_assoc_dx = calc_da_assoc_dx_autodiff(T, rho);
        for (int i = 0; i < nc_; ++i) {
            da_dx[i] += da_assoc_dx[i];
        }
    }

    // Calculate Σ_j x_j (∂a^res/∂x_j)
    double sum_x_dadx = 0.0;
    for (int j = 0; j < nc_; ++j) {
        sum_x_dadx += x[j] * da_dx[j];
    }

    // Apply formula: μ_i^res/(RT) = a^res + (∂a^res/∂x_i) - Σ_j x_j(∂a^res/∂x_j)
    // Then convert to ln(φ_i) = μ_i^res/(RT) - ln(Z)
    double Z = calculateZ(T, rho);
    double ln_Z = std::log(Z);

    // DEBUG: Check derivative values
    static int deriv_debug = 0;
    if (deriv_debug < 3 && rho > 10000.0) {
        std::cout << "\n[DEBUG Composition Derivatives]\n";
        std::cout << "  T=" << T << " K, rho=" << rho << " mol/m³\n";
        std::cout << "  a_res=" << a_res << "\n";
        std::cout << "  da_hc_dx: [" << da_hc_dx[0];
        for (int i = 1; i < nc_; ++i) std::cout << ", " << da_hc_dx[i];
        std::cout << "]\n";
        std::cout << "  da_disp_dx: [" << da_disp_dx[0];
        for (int i = 1; i < nc_; ++i) std::cout << ", " << da_disp_dx[i];
        std::cout << "]\n";
        std::cout << "  da_dx (total): [" << da_dx[0];
        for (int i = 1; i < nc_; ++i) std::cout << ", " << da_dx[i];
        std::cout << "]\n";
        std::cout << "  sum_x_dadx=" << sum_x_dadx << "\n";
        std::cout << "  Z=" << Z << ", ln(Z)=" << ln_Z << "\n";
        for (int i = 0; i < nc_; ++i) {
            double mu_res_over_RT = a_res + da_dx[i] - sum_x_dadx;
            std::cout << "  Component " << i << ": mu_res/RT=" << mu_res_over_RT
                      << ", ln(φ)=" << (mu_res_over_RT - ln_Z)
                      << ", φ=" << std::exp(mu_res_over_RT - ln_Z) << "\n";
        }
        deriv_debug++;
    }

    for (int i = 0; i < nc_; ++i) {
        double mu_res_over_RT = a_res + da_dx[i] - sum_x_dadx;
        mu_res[i] = mu_res_over_RT - ln_Z;  // This is ln(φ_i)
    }

    return mu_res;
}

// ============================================================================
// Analytical composition derivatives - VALIDATION ONLY
// ============================================================================
// WARNING: These analytical derivative implementations contain fundamental errors
// (sign errors, factor-of-2 errors, magnitude errors) discovered during testing.
//
// PRODUCTION CODE USES: calc_da_hc_dx_autodiff(), calc_da_disp_dx_autodiff()
//
// These functions are KEPT FOR VALIDATION PURPOSES ONLY to allow comparison
// between analytical formulas from literature and machine-precision autodiff.
// See docs/debugging/ANALYTICAL_DERIVATIVES_STATUS.md for full error analysis.
// ============================================================================

// [VALIDATION ONLY] Hard-chain derivative - contains fundamental errors
std::vector<double> PCSaftEOS::calc_da_hc_dx(double T, double rho) const {
    // Calculate ∂a^{hc}/∂x_k at constant T, ρ, x_{j≠k}
    // From Gross & Sadowski (2001), Appendix A, Equation A.35:
    //
    // (∂ã^{hc}/∂x_k) = m_k ã^{hs} + m̄ (∂ã^{hs}/∂x_k) -
    //                  Σ_i x_i(m_i - 1)(g_{ii}^{hs})^{-1}(∂g_{ii}^{hs}/∂x_k)

    std::vector<double> da_dx(nc_);
    const std::vector<double>& x = mixture_.getMoleFractions();
    std::vector<double> d = calc_d(T);

    double eta = calc_eta(T, rho);
    if (eta <= 0.0 || eta >= 1.0) {
        return da_dx;  // Return zeros
    }

    double m_bar = calc_m_bar();
    double rho_reduced = rho * N_A * 1e-30;

    // Calculate ζ_n values (packing fraction moments)
    std::vector<double> zeta(4, 0.0);
    for (int i = 0; i < nc_; ++i) {
        double d_i3 = d[i] * d[i] * d[i];
        zeta[0] += x[i] * m_[i];
        zeta[1] += x[i] * m_[i] * d[i];
        zeta[2] += x[i] * m_[i] * d[i] * d[i];
        zeta[3] += x[i] * m_[i] * d_i3;
    }
    for (int n = 0; n < 4; ++n) {
        zeta[n] *= (PI / 6.0) * rho_reduced;
    }

    // Calculate hard-sphere contribution using full mixture formula
    double a_hs = calc_a_hs_mixture(zeta);
    double one_minus_eta = 1.0 - zeta[3];

    for (int k = 0; k < nc_; ++k) {
        // Calculate ζ_{n,xk} = ∂ζ_n/∂x_k (Equation A.34)
        std::vector<double> zeta_xk(4);
        double d_k3 = d[k] * d[k] * d[k];
        zeta_xk[0] = (PI / 6.0) * rho_reduced * m_[k];
        zeta_xk[1] = (PI / 6.0) * rho_reduced * m_[k] * d[k];
        zeta_xk[2] = (PI / 6.0) * rho_reduced * m_[k] * d[k] * d[k];
        zeta_xk[3] = (PI / 6.0) * rho_reduced * m_[k] * d_k3;

        // Calculate (∂ã^{hs}/∂x_k) (Equation A.36)
        double da_hs_dxk = -(zeta_xk[0] / zeta[0]) * a_hs + (1.0 / zeta[0]) * (
            3.0 * (zeta_xk[1] * zeta[2] + zeta[1] * zeta_xk[2]) / one_minus_eta +
            3.0 * zeta[1] * zeta[2] * zeta_xk[3] / (one_minus_eta * one_minus_eta) +
            3.0 * zeta[2] * zeta[2] * zeta_xk[2] / (zeta[3] * one_minus_eta * one_minus_eta) +
            zeta[2] * zeta[2] * zeta[2] * zeta_xk[3] * (3.0 * zeta[3] - 1.0) /
                (zeta[3] * zeta[3] * one_minus_eta * one_minus_eta * one_minus_eta) +
            (3.0 * zeta[2] * zeta[2] * zeta_xk[2] * zeta[3] - 2.0 * zeta[2] * zeta[2] * zeta[2] * zeta_xk[3]) /
                (zeta[3] * zeta[3] * zeta[3]) * std::log(one_minus_eta) +
            (zeta[0] - zeta[2] * zeta[2] * zeta[2] / (zeta[3] * zeta[3])) * zeta_xk[3] / one_minus_eta
        );

        // First term: m_k ã^{hs}
        double term1 = m_[k] * a_hs;

        // Second term: m̄ (∂ã^{hs}/∂x_k)
        double term2 = m_bar * da_hs_dxk;

        // Third term: -(m_k - 1) ln g_kk^{hs}
        // This comes from the product rule: ∂/∂x_k [Σ_i x_i(m_i-1)ln g_ii] includes a direct term
        double term3 = 0.0;
        if (m_[k] > 1.0) {
            double g_kk = calc_g_ii(d[k], zeta[2], zeta[3]);
            term3 = -(m_[k] - 1.0) * std::log(g_kk);
        }

        // Fourth term: - Σ_i x_i(m_i - 1)(g_{ii}^{hs})^{-1}(∂g_{ii}^{hs}/∂x_k)
        // This is the indirect term from how g_ii depends on x_k through ζ_n
        double term4 = 0.0;
        for (int i = 0; i < nc_; ++i) {
            if (m_[i] <= 1.0) continue;

            // Calculate component-specific g_{ii}^{hs} for i-i pair
            double g_ii = calc_g_ii(d[i], zeta[2], zeta[3]);

            // Calculate (∂g_{ii}^{hs}/∂x_k) (Equation A.37, for i=i case)
            // From G&S: (d_i/2) and (d_i²/4) factors for i=i
            double dg_ii_dxk = zeta_xk[3] / (one_minus_eta * one_minus_eta) +
                (d[i] / 2.0) * (3.0 * zeta_xk[2] / (one_minus_eta * one_minus_eta) +
                                6.0 * zeta[2] * zeta_xk[3] / (one_minus_eta * one_minus_eta * one_minus_eta)) +
                (d[i] * d[i] / 4.0) * (4.0 * zeta[2] * zeta_xk[2] / (one_minus_eta * one_minus_eta * one_minus_eta) +
                                       6.0 * zeta[2] * zeta[2] * zeta_xk[3] /
                                           (one_minus_eta * one_minus_eta * one_minus_eta * one_minus_eta));

            term4 -= x[i] * (m_[i] - 1.0) * dg_ii_dxk / g_ii;
        }

        da_dx[k] = term1 + term2 + term3 + term4;
    }

    return da_dx;
}

// [VALIDATION ONLY] Dispersion derivative - contains fundamental errors
std::vector<double> PCSaftEOS::calc_da_disp_dx(double T, double rho) const {
    // Calculate ∂a^{disp}/∂x_k at constant T, ρ, x_{j≠k}
    // From Gross & Sadowski (2001), Appendix A, Equation A.38:
    //
    // (∂ã^{disp}/∂x_k) = -2πρ[I_{1,xk} m̄²εσ³ + I_1(m̄²εσ³)_{xk}] -
    //                     πρ[m_k C_1 I_2 + m̄ C_{1,xk} I_2 + m̄ C_1 I_{2,xk}] m̄²ε²σ³ +
    //                     m̄ C_1 I_2(m̄²ε²σ³)_{xk}

    std::vector<double> da_dx(nc_);
    const std::vector<double>& x = mixture_.getMoleFractions();
    std::vector<double> d = calc_d(T);

    double eta = calc_eta(T, rho);
    if (eta <= 0.0) {
        return da_dx;
    }

    double m_bar = calc_m_bar();
    double rho_reduced = rho * N_A * 1e-30;

    // Calculate m̄²εσ³ and m̄²ε²σ³ (Equations A.12, A.13)
    double m2es3 = 0.0;
    double m2e2s3 = 0.0;

    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            double sigma_ij = 0.5 * (sigma_[i] + sigma_[j]);
            double epsilon_ij = std::sqrt(epsilon_[i] * epsilon_[j]);
            double kij = mixture_.getBinaryParameter(i, j);
            epsilon_ij *= (1.0 - kij);
            double sigma_ij3 = sigma_ij * sigma_ij * sigma_ij;
            double eps_over_T = epsilon_ij / T;

            m2es3 += x[i] * x[j] * m_[i] * m_[j] * eps_over_T * sigma_ij3;
            m2e2s3 += x[i] * x[j] * m_[i] * m_[j] * eps_over_T * eps_over_T * sigma_ij3;
        }
    }

    // Calculate I_1, I_2
    double I1, I2;
    calc_dispersion_integrals(eta, I1, I2);

    // Calculate C_1 (Equation A.11)
    double C1 = calc_C1(eta);

    for (int k = 0; k < nc_; ++k) {
        // Calculate ζ_{3,xk} = ∂ζ_3/∂x_k (Equation A.34, n=3)
        double d_k3 = d[k] * d[k] * d[k];
        double zeta_3_xk = (PI / 6.0) * rho_reduced * m_[k] * d_k3;

        // Calculate (m̄²εσ³)_{xk} (Equation A.39)
        double m2es3_xk = 0.0;
        for (int j = 0; j < nc_; ++j) {
            double sigma_kj = 0.5 * (sigma_[k] + sigma_[j]);
            double epsilon_kj = std::sqrt(epsilon_[k] * epsilon_[j]);
            double kkj = mixture_.getBinaryParameter(k, j);
            epsilon_kj *= (1.0 - kkj);
            double sigma_kj3 = sigma_kj * sigma_kj * sigma_kj;
            double eps_kj_over_T = epsilon_kj / T;

            m2es3_xk += 2.0 * x[j] * m_[k] * m_[j] * eps_kj_over_T * sigma_kj3;
        }

        // Calculate (m̄²ε²σ³)_{xk} (Equation A.40)
        double m2e2s3_xk = 0.0;
        for (int j = 0; j < nc_; ++j) {
            double sigma_kj = 0.5 * (sigma_[k] + sigma_[j]);
            double epsilon_kj = std::sqrt(epsilon_[k] * epsilon_[j]);
            double kkj = mixture_.getBinaryParameter(k, j);
            epsilon_kj *= (1.0 - kkj);
            double sigma_kj3 = sigma_kj * sigma_kj * sigma_kj;
            double eps_kj_over_T = epsilon_kj / T;

            m2e2s3_xk += 2.0 * x[j] * m_[k] * m_[j] * eps_kj_over_T * eps_kj_over_T * sigma_kj3;
        }

        // Calculate C_{1,xk} (Equation A.41)
        double C2 = calc_C2(eta);
        double one_minus_eta = 1.0 - eta;
        double two_minus_eta = 2.0 - eta;
        double C1_xk = C2 * zeta_3_xk - C1 * C1 * (
            m_[k] * (8.0 * eta - 2.0 * eta * eta) / std::pow(one_minus_eta, 4) -
            m_[k] * (20.0 * eta - 27.0 * eta * eta + 12.0 * eta * eta * eta - 2.0 * eta * eta * eta * eta) /
                (one_minus_eta * one_minus_eta * two_minus_eta * two_minus_eta)
        );

        // Calculate I_{1,xk} and I_{2,xk} (Equations A.42-A.45)
        double I1_xk = calc_I_derivative(eta, m_bar, m_[k], zeta_3_xk, 1);
        double I2_xk = calc_I_derivative(eta, m_bar, m_[k], zeta_3_xk, 2);

        // Assemble the full derivative (Equation A.38)
        double term1 = -2.0 * PI * rho_reduced * (I1_xk * m2es3 + I1 * m2es3_xk);

        double term2 = -PI * rho_reduced * (
            m_[k] * C1 * I2 * m2e2s3 +
            m_bar * C1_xk * I2 * m2e2s3 +
            m_bar * C1 * I2_xk * m2e2s3 +
            m_bar * C1 * I2 * m2e2s3_xk
        );

        da_dx[k] = term1 + term2;
    }

    return da_dx;
}

// Numerical association composition derivative (fixes buggy analytical version)
std::vector<double> PCSaftEOS::calc_da_assoc_dx_numerical(double T, double rho) const {
    // Use central finite differences to calculate ∂a^{assoc}/∂x_k
    const std::vector<double>& x = mixture_.getMoleFractions();
    int nc = nc_;
    std::vector<double> da_dx_unc(nc);  // Unconstrained derivatives

    const double h = 1e-7;  // Perturbation size

    // Step 1: Calculate constrained derivatives ∂a/∂x_k maintaining Σx_i = 1
    // Use perturbation: increase x_k, decrease x_(nc-1) by same amount
    std::vector<double> da_dx(nc);

    for (int k = 0; k < nc - 1; ++k) {
        // Create perturbed compositions maintaining Σx_i = 1
        std::vector<double> x_plus = x;
        std::vector<double> x_minus = x;

        x_plus[k] += h;
        x_plus[nc-1] -= h;  // Maintain sum = 1

        x_minus[k] -= h;
        x_minus[nc-1] += h;  // Maintain sum = 1

        // Create temporary mixtures with perturbed compositions
        Mixture mix_plus(mixture_.getComponents(), x_plus);
        Mixture mix_minus(mixture_.getComponents(), x_minus);

        // Copy binary interaction parameters
        for (int i = 0; i < nc; ++i) {
            for (int j = i+1; j < nc; ++j) {
                double kij = mixture_.getBinaryParameter(i, j);
                mix_plus.setBinaryParameter(i, j, kij);
                mix_minus.setBinaryParameter(i, j, kij);
            }
        }

        // Calculate association contributions
        PCSaftEOS eos_plus(mix_plus);
        PCSaftEOS eos_minus(mix_minus);

        double a_plus = eos_plus.calc_a_assoc(T, rho);
        double a_minus = eos_minus.calc_a_assoc(T, rho);

        // Central difference gives constrained derivative directly
        da_dx[k] = (a_plus - a_minus) / (2.0 * h);
    }

    // Last component: by constraint Σ(∂a/∂x_i)_constr = 0
    // Therefore: ∂a/∂x_{nc-1} = -Σ_{i=0}^{nc-2} (∂a/∂x_i)
    double sum = 0.0;
    for (int k = 0; k < nc - 1; ++k) {
        sum += da_dx[k];
    }
    da_dx[nc-1] = -sum;

    return da_dx;
}

// [VALIDATION ONLY] Association derivative - contains fundamental errors
std::vector<double> PCSaftEOS::calc_da_assoc_dx(double T, double rho) const {
    // Calculate ∂a^{assoc}/∂x_k at constant T, ρ
    //
    // a^{assoc} = Σ_i x_i [Σ_A (ln(X_A_i) - X_A_i/2) + M_i/2]
    //
    // This is complex and requires solving for ∂X_A/∂x_k
    // For now: numerical derivative as fallback

    std::vector<double> da_dx(nc_);

    // Numerical derivative (placeholder for analytical)
    const std::vector<double>& x = mixture_.getMoleFractions();
    double a_base = calc_a_assoc(T, rho);

    const double dx_pert = 1e-6;
    for (int k = 0; k < nc_; ++k) {
        std::vector<double> x_pert = x;
        int j_reduce = (k == 0) ? 1 : 0;
        for (int j = 0; j < nc_; ++j) {
            if (j != k && x[j] > x[j_reduce]) j_reduce = j;
        }

        double delta = std::min(dx_pert, x[j_reduce] * 0.5);
        x_pert[k] += delta;
        x_pert[j_reduce] -= delta;

        std::vector<Component> comps = mixture_.getComponents();
        Mixture mix_pert(comps, x_pert);
        PCSaftEOS eos_pert(mix_pert);

        double a_pert = eos_pert.calc_a_assoc(T, rho);
        da_dx[k] = (a_pert - a_base) / delta;
    }

    return da_dx;
}

// ============================================================================
// Helper functions
// ============================================================================

double PCSaftEOS::calc_eta(double T, double rho) const {
    const std::vector<double>& x = mixture_.getMoleFractions();
    std::vector<double> d = calc_d(T);

    double zeta3 = 0.0;
    for (int i = 0; i < nc_; ++i) {
        double d_i = d[i];
        zeta3 += x[i] * m_[i] * d_i * d_i * d_i;
    }

    // Convert density from mol/m³ to molecules/Å³
    // Factor: N_A [molecules/mol] × 10^-30 [Å³/m³]
    double rho_reduced = rho * N_A * 1e-30;

    double eta = PI / 6.0 * rho_reduced * zeta3;

    return eta;
}

std::vector<double> PCSaftEOS::calc_d(double T) const {
    // Temperature-dependent hard-sphere diameter
    // d_i = σ_i [1 - 0.12 exp(-3 ε_i/kT)]
    std::vector<double> d(nc_);

    for (int i = 0; i < nc_; ++i) {
        d[i] = sigma_[i] * (1.0 - 0.12 * std::exp(-3.0 * epsilon_[i] / T));
    }

    return d;
}

double PCSaftEOS::calc_m_bar() const {
    const std::vector<double>& x = mixture_.getMoleFractions();

    double m_bar = 0.0;
    for (int i = 0; i < nc_; ++i) {
        m_bar += x[i] * m_[i];
    }

    return m_bar;
}

double PCSaftEOS::calc_C1(double eta) const {
    // Calculate compressibility factor C_1
    // From Gross & Sadowski (2001), Appendix A, Equation A.11
    //
    // C_1 = (1 + Z^{hc} + ρ ∂Z^{hc}/∂ρ)^{-1}
    //     = 1 / (1 + m̄ (8η-2η²)/(1-η)⁴ + (1-m̄) (20η-27η²+12η³-2η⁴)/[(1-η)(2-η)]²)

    double m_bar = calc_m_bar();
    double one_minus_eta = 1.0 - eta;
    double two_minus_eta = 2.0 - eta;

    if (one_minus_eta <= 0.0) return 0.0;

    double term1 = m_bar * (8.0 * eta - 2.0 * eta * eta) / std::pow(one_minus_eta, 4);
    double term2 = (1.0 - m_bar) * (20.0 * eta - 27.0 * eta * eta + 12.0 * eta * eta * eta - 2.0 * eta * eta * eta * eta) /
                   (one_minus_eta * one_minus_eta * two_minus_eta * two_minus_eta);

    double C1 = 1.0 / (1.0 + term1 + term2);
    return C1;
}

double PCSaftEOS::calc_C2(double eta) const {
    // Calculate derivative of C_1 with respect to eta
    // From Gross & Sadowski (2001), Appendix A, Equation A.31
    //
    // C_2 = ∂C_1/∂η = -C_1² [m̄ (-4η²+20η+8)/(1-η)⁵ +
    //                         (1-m̄) (2η³+12η²-48η+40)/[(1-η)(2-η)]³]

    double m_bar = calc_m_bar();
    double one_minus_eta = 1.0 - eta;
    double two_minus_eta = 2.0 - eta;

    if (one_minus_eta <= 0.0) return 0.0;

    double C1 = calc_C1(eta);

    double term1 = m_bar * (-4.0 * eta * eta + 20.0 * eta + 8.0) / std::pow(one_minus_eta, 5);

    double eta2 = eta * eta;
    double eta3 = eta2 * eta;
    double term2 = (1.0 - m_bar) * (2.0 * eta3 + 12.0 * eta2 - 48.0 * eta + 40.0) /
                   (one_minus_eta * one_minus_eta * one_minus_eta *
                    two_minus_eta * two_minus_eta * two_minus_eta);

    double C2 = -C1 * C1 * (term1 + term2);
    return C2;
}

double PCSaftEOS::calc_I_derivative(double eta, double m_bar, double m_k, double zeta_3_xk, int order) const {
    // Calculate derivative of dispersion integrals I_1 or I_2 with respect to composition
    // From Gross & Sadowski (2001), Appendix A, Equations A.42-A.45
    //
    // I_{n,xk} = Σ_{i=0}^6 [a_i(m̄) i ζ_{3,xk} η^{i-1} + a_{i,xk} η^i]  for n=1
    // I_{n,xk} = Σ_{i=0}^6 [b_i(m̄) i ζ_{3,xk} η^{i-1} + b_{i,xk} η^i]  for n=2
    //
    // where a_{i,xk} = (m_k/m̄²) a_{1i} + (m_k/m̄²)(3 - 4/m̄) a_{2i}
    //       b_{i,xk} = (m_k/m̄²) b_{1i} + (m_k/m̄²)(3 - 4/m̄) b_{2i}

    double I_xk = 0.0;
    double m_bar2 = m_bar * m_bar;

    // Coefficient derivatives (Equations A.44, A.45)
    double coeff_factor1 = m_k / m_bar2;
    double coeff_factor2 = (m_k / m_bar2) * (3.0 - 4.0 / m_bar);

    for (int i = 0; i < 7; ++i) {
        double eta_power_i = std::pow(eta, i);
        double eta_power_i_minus_1 = (i > 0) ? std::pow(eta, i - 1) : 0.0;

        // Calculate a_i(m_bar) or b_i(m_bar)
        double ai_or_bi;
        if (order == 1) {
            // For I_1 derivative: use a coefficients
            ai_or_bi = a0[i] + ((m_bar - 1.0) / m_bar) * a1[i] +
                       ((m_bar - 1.0) / m_bar) * ((m_bar - 2.0) / m_bar) * a2[i];
        } else {
            // For I_2 derivative: use b coefficients
            ai_or_bi = b0[i] + ((m_bar - 1.0) / m_bar) * b1[i] +
                       ((m_bar - 1.0) / m_bar) * ((m_bar - 2.0) / m_bar) * b2[i];
        }

        // Calculate coefficient derivative (a_{i,xk} or b_{i,xk})
        double coeff_xk;
        if (order == 1) {
            coeff_xk = coeff_factor1 * a1[i] + coeff_factor2 * a2[i];
        } else {
            coeff_xk = coeff_factor1 * b1[i] + coeff_factor2 * b2[i];
        }

        // Equation A.42 or A.43
        I_xk += ai_or_bi * i * zeta_3_xk * eta_power_i_minus_1 + coeff_xk * eta_power_i;
    }

    return I_xk;
}

// ============================================================================
// AUTODIFF Implementation - Phase 1: Parallel AD-based derivatives
// ============================================================================

// Template helper functions to work with both double and autodiff::dual types

/**
 * @brief Templated version of calc_d - temperature-dependent diameter
 */
template<typename T>
std::vector<T> PCSaftEOS::calc_d_templated(T T_val) const {
    std::vector<T> d(nc_);
    for (int i = 0; i < nc_; ++i) {
        d[i] = sigma_[i] * (1.0 - 0.12 * exp(-3.0 * epsilon_[i] / T_val));
    }
    return d;
}

/**
 * @brief Templated version of calc_m_bar - mean segment number
 */
template<typename T>
T PCSaftEOS::calc_m_bar_templated(const std::vector<T>& x) const {
    T m_bar = 0.0;
    for (int i = 0; i < nc_; ++i) {
        m_bar += x[i] * m_[i];
    }
    return m_bar;
}

/**
 * @brief Templated version of calc_g_ii - component-specific RDF
 * Uses smooth clamping consistent with non-templated version
 */
template<typename T>
T PCSaftEOS::calc_g_ii_templated(T d_i, T zeta_2, T zeta_3) const {
    // Smooth clamping consistent with non-templated version
    constexpr double zeta_min = 1e-15;
    constexpr double zeta3_max = 0.999;

    T zeta3_safe = smoothClampT(zeta_3, zeta_min, zeta3_max);

    T one_minus_zeta3 = 1.0 - zeta3_safe;
    T one_minus_zeta3_sq = one_minus_zeta3 * one_minus_zeta3;
    T one_minus_zeta3_cube = one_minus_zeta3_sq * one_minus_zeta3;

    T term1 = 1.0 / one_minus_zeta3;
    T term2 = (d_i / 2.0) * 3.0 * zeta_2 / one_minus_zeta3_sq;
    T term3 = (d_i * d_i / 4.0) * 2.0 * zeta_2 * zeta_2 / one_minus_zeta3_cube;

    return term1 + term2 + term3;
}

/**
 * @brief Templated version of calc_a_hs_mixture - Mansoori formula
 * Uses smooth clamping consistent with non-templated version for autodiff compatibility
 */
template<typename T>
T PCSaftEOS::calc_a_hs_mixture_templated(const std::vector<T>& zeta) const {
    // Smooth clamping consistent with non-templated version
    constexpr double zeta_min = 1e-15;
    constexpr double zeta3_max = 0.999;

    T zeta0_safe = smoothClampLowerT(zeta[0], zeta_min);
    T zeta3_safe = smoothClampT(zeta[3], zeta_min, zeta3_max);

    T one_minus_zeta3 = 1.0 - zeta3_safe;
    T one_minus_zeta3_sq = one_minus_zeta3 * one_minus_zeta3;

    T term1 = 3.0 * zeta[1] * zeta[2] / one_minus_zeta3;
    T term2 = zeta[2] * zeta[2] * zeta[2] / (zeta3_safe * one_minus_zeta3_sq);
    T term3 = (zeta[2] * zeta[2] * zeta[2] / (zeta3_safe * zeta3_safe) - zeta[0]) * log(one_minus_zeta3);

    return (term1 + term2 + term3) / zeta0_safe;
}

/**
 * @brief Templated version of calc_a_hc - hard-chain contribution
 *
 * This is the core function that autodiff will differentiate.
 * Works with both double and autodiff::dual types.
 */
template<typename T>
T PCSaftEOS::calc_a_hc_templated(const std::vector<T>& x, T T_val, T rho) const {
    // Calculate temperature-dependent diameters
    std::vector<T> d = calc_d_templated(T_val);

    // Calculate packing fractions ζ0, ζ1, ζ2, ζ3
    std::vector<T> zeta(4, T(0.0));
    for (int i = 0; i < nc_; ++i) {
        T d_i = d[i];
        T d_i_pow = 1.0;
        for (int n = 0; n < 4; ++n) {
            zeta[n] += x[i] * m_[i] * d_i_pow;
            d_i_pow *= d_i;
        }
    }

    // Convert density from mol/m³ to molecules/Å³
    T rho_reduced = rho * N_A * 1e-30;

    // Multiply by π/6 * ρ
    T factor = PI / 6.0 * rho_reduced;
    for (int k = 0; k < 4; ++k) {
        zeta[k] *= factor;
    }

    // Hard-sphere contribution (Mansoori formula)
    T a_hs = calc_a_hs_mixture_templated(zeta);

    // Mean segment number
    T m_bar = calc_m_bar_templated(x);

    // Hard-chain contribution: a_hc = m_bar * a_hs - Σ_i x_i(m_i - 1) ln(g_ii)
    T sum_chain = 0.0;
    for (int i = 0; i < nc_; ++i) {
        if (m_[i] > 1.0) {
            T g_ii = calc_g_ii_templated(d[i], zeta[2], zeta[3]);
            sum_chain += x[i] * (m_[i] - 1.0) * log(g_ii);
        }
    }

    T a_hc = m_bar * a_hs - sum_chain;

    return a_hc;
}

/**
 * @brief Templated version of calc_eta - packing fraction
 */
template<typename T>
T PCSaftEOS::calc_eta_templated(const std::vector<T>& x, T T_val, T rho) const {
    std::vector<T> d = calc_d_templated(T_val);

    T zeta3 = 0.0;
    for (int i = 0; i < nc_; ++i) {
        T d_i = d[i];
        zeta3 += x[i] * m_[i] * d_i * d_i * d_i;
    }

    // Convert density from mol/m³ to molecules/Å³
    T rho_reduced = rho * N_A * 1e-30;

    return PI / 6.0 * rho_reduced * zeta3;
}

template<typename T>
T PCSaftEOS::calc_g_hs_templated(T eta) const {
    // Simplified radial distribution function at contact
    // Uses smooth clamping consistent with non-templated version
    constexpr double eta_min = 1e-15;
    constexpr double eta_max = 0.999;

    T eta_safe = smoothClampT(eta, eta_min, eta_max);

    T one_minus_eta = T(1.0) - eta_safe;
    return (T(1.0) - T(0.5) * eta_safe) / (one_minus_eta * one_minus_eta * one_minus_eta);
}

/**
 * @brief Templated version of calc_dispersion_integrals
 */
template<typename T>
void PCSaftEOS::calc_dispersion_integrals_templated(T eta, T m_bar, T& I1, T& I2) const {
    I1 = 0.0;
    I2 = 0.0;

    for (int i = 0; i < 7; ++i) {
        T eta_pow = pow(eta, i);
        T factor1 = (m_bar - 1.0) / m_bar;
        T factor2 = (m_bar - 1.0) * (m_bar - 2.0) / (m_bar * m_bar);

        I1 += (a0[i] + factor1 * a1[i] + factor2 * a2[i]) * eta_pow;
        I2 += (b0[i] + factor1 * b1[i] + factor2 * b2[i]) * eta_pow;
    }
}

/**
 * @brief Templated version of calc_a_disp - dispersion contribution
 * Updated to include full C1 expression consistent with non-templated version
 */
template<typename T>
T PCSaftEOS::calc_a_disp_templated(const std::vector<T>& x, T T_val, T rho) const {
    std::vector<T> d = calc_d_templated(T_val);

    T eta = calc_eta_templated(x, T_val, rho);

    // Apply smooth clamping consistent with non-templated version
    constexpr double eta_min = 1e-15;
    constexpr double eta_max = 0.999;
    eta = smoothClampT(eta, eta_min, eta_max);

    T m_bar = calc_m_bar_templated(x);

    // Calculate m²εσ³ and m²ε²σ³
    T m2es3 = 0.0;   // Σ_i Σ_j x_i x_j m_i m_j ε_ij σ_ij³
    T m2e2s3 = 0.0;  // Σ_i Σ_j x_i x_j m_i m_j ε_ij² σ_ij³

    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            // Combining rules
            double sigma_ij = 0.5 * (sigma_[i] + sigma_[j]);
            double epsilon_ij = std::sqrt(epsilon_[i] * epsilon_[j]);

            // Apply binary interaction parameter
            double kij = mixture_.getBinaryParameter(i, j);
            epsilon_ij *= (1.0 - kij);

            double sigma_ij3 = sigma_ij * sigma_ij * sigma_ij;

            m2es3 += x[i] * x[j] * m_[i] * m_[j] * epsilon_ij * sigma_ij3;
            m2e2s3 += x[i] * x[j] * m_[i] * m_[j] * epsilon_ij * epsilon_ij * sigma_ij3;
        }
    }

    // Calculate dispersion integrals I1 and I2
    T I1, I2;
    calc_dispersion_integrals_templated(eta, m_bar, I1, I2);

    // Convert density from mol/m³ to molecules/Å³
    T rho_reduced = rho * N_A * 1e-30;

    // C1 compressibility correction - FULL expression matching non-templated version
    // From Gross & Sadowski (2001)
    T eta2 = eta * eta;
    T eta3 = eta * eta * eta;
    T eta4 = eta2 * eta2;
    T one_minus_eta = T(1.0) - eta;
    T one_minus_eta2 = one_minus_eta * one_minus_eta;
    T one_minus_eta4 = one_minus_eta2 * one_minus_eta2;
    T two_minus_eta = T(2.0) - eta;
    T two_minus_eta2 = two_minus_eta * two_minus_eta;

    T C1_denom = T(1.0) + m_bar * (T(8.0) * eta - T(2.0) * eta2) / one_minus_eta4
                 + (T(1.0) - m_bar) * (T(20.0) * eta - T(27.0) * eta2 + T(12.0) * eta3 - T(2.0) * eta4)
                   / (one_minus_eta2 * two_minus_eta2);
    T C1 = T(1.0) / C1_denom;

    // First-order perturbation term (mean-attractive energy)
    T a_1 = T(-2.0) * PI * rho_reduced * I1 * m2es3 / T_val;

    // Second-order perturbation term (fluctuation term)
    T a_2 = T(-1.0) * PI * rho_reduced * m_bar * C1 * I2 * m2e2s3 / (T_val * T_val);

    return a_1 + a_2;
}

/**
 * @brief Calculate hard-chain composition derivatives using autodiff
 *
 * This computes CONSTRAINED derivatives respecting Σx_i = 1.
 * The constrained derivative is: ∂a/∂x_k|_{Σx=1} = ∂a/∂x_k|_{unc} - ∂a/∂x_{nc-1}|_{unc}
 */
std::vector<double> PCSaftEOS::calc_da_hc_dx_autodiff(double T, double rho) const {
    using namespace autodiff;

    const std::vector<double>& x_double = mixture_.getMoleFractions();
    std::vector<double> da_dx_unc(nc_);  // Unconstrained derivatives

    // Step 1: Compute unconstrained derivatives for all components
    for (int k = 0; k < nc_; ++k) {
        // Create dual number composition vector
        std::vector<dual> x_dual(nc_);
        for (int i = 0; i < nc_; ++i) {
            x_dual[i] = x_double[i];
        }

        // Define function: a_hc as a function of x[k]
        auto f = [&](dual xk) -> dual {
            std::vector<dual> x_temp = x_dual;
            x_temp[k] = xk;
            return calc_a_hc_templated(x_temp, dual(T), dual(rho));
        };

        // Compute unconstrained derivative using autodiff
        dual xk = x_double[k];
        da_dx_unc[k] = derivative(f, wrt(xk), at(xk));
    }

    // Step 2: Apply constraint transformation
    // For constrained derivatives respecting Σx_i = 1:
    // ∂a/∂x_k|_{Σx=1} = ∂a/∂x_k|_{unc} - ∂a/∂x_ref|_{unc}
    // where ref is chosen as:
    //   - For k < nc-1: ref = nc-1 (last component is dependent)
    //   - For k = nc-1: ref = 0 (first component is dependent)
    std::vector<double> da_dx(nc_);

    // DEBUG: Check unconstrained derivatives
    static int autodiff_debug = 0;
    if (autodiff_debug < 2 && rho > 10000.0) {
        std::cout << "\n[DEBUG calc_da_hc_dx_autodiff]\n";
        std::cout << "  Unconstrained derivatives: [";
        for (int i = 0; i < nc_; ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << da_dx_unc[i];
        }
        std::cout << "]\n";
    }

    // Apply correct constraint transformation for composition derivatives
    // For constrained derivatives respecting Σx_i = 1:
    // ∂a/∂x_k|_const = ∂a/∂x_k|_unc - Σ_j x_j (∂a/∂x_j|_unc)
    // This ensures Gibbs-Duhem relation: Σ_k x_k (∂a/∂x_k|_const) = 0
    double sum_x_dadx_unc = 0.0;
    for (int j = 0; j < nc_; ++j) {
        sum_x_dadx_unc += x_double[j] * da_dx_unc[j];
    }

    for (int k = 0; k < nc_; ++k) {
        da_dx[k] = da_dx_unc[k] - sum_x_dadx_unc;
    }

    if (autodiff_debug < 2 && rho > 10000.0) {
        std::cout << "  Constrained derivatives: [";
        for (int i = 0; i < nc_; ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << da_dx[i];
        }
        std::cout << "]\n";
        autodiff_debug++;
    }

    return da_dx;
}

/**
 * @brief Calculate dispersion composition derivatives using autodiff
 *
 * This computes CONSTRAINED derivatives respecting Σx_i = 1.
 * Same approach as hard-chain: compute unconstrained, then apply constraint transformation.
 */
std::vector<double> PCSaftEOS::calc_da_disp_dx_autodiff(double T, double rho) const {
    using namespace autodiff;

    const std::vector<double>& x_double = mixture_.getMoleFractions();
    std::vector<double> da_dx_unc(nc_);  // Unconstrained derivatives

    // Step 1: Compute unconstrained derivatives for all components
    for (int k = 0; k < nc_; ++k) {
        // Create dual number composition vector
        std::vector<dual> x_dual(nc_);
        for (int i = 0; i < nc_; ++i) {
            x_dual[i] = x_double[i];
        }

        // Define function: a_disp as a function of x[k]
        auto f = [&](dual xk) -> dual {
            std::vector<dual> x_temp = x_dual;
            x_temp[k] = xk;
            return calc_a_disp_templated(x_temp, dual(T), dual(rho));
        };

        // Compute unconstrained derivative using autodiff
        dual xk = x_double[k];
        da_dx_unc[k] = derivative(f, wrt(xk), at(xk));
    }

    // Step 2: Apply correct constraint transformation for composition derivatives
    // For constrained derivatives respecting Σx_i = 1:
    // ∂a/∂x_k|_const = ∂a/∂x_k|_unc - Σ_j x_j (∂a/∂x_j|_unc)
    // This ensures Gibbs-Duhem relation: Σ_k x_k (∂a/∂x_k|_const) = 0
    std::vector<double> da_dx(nc_);

    double sum_x_dadx_unc = 0.0;
    for (int j = 0; j < nc_; ++j) {
        sum_x_dadx_unc += x_double[j] * da_dx_unc[j];
    }

    for (int k = 0; k < nc_; ++k) {
        da_dx[k] = da_dx_unc[k] - sum_x_dadx_unc;
    }

    return da_dx;
}

// ============================================================================
// AUTODIFF: Templated Association Functions
// ============================================================================

template<typename T>
T PCSaftEOS::calc_delta_assoc_templated(int i, int j, T T_val, const std::vector<T>& d, T eta) const {
    // Association strength Δ^AB = g(d_ij) * [exp(ε^AB/kT) - 1] * κ^AB * σ_ij³

    const AssociationParams& assoc_i = mixture_.getComponent(i).getAssociationParams();
    const AssociationParams& assoc_j = mixture_.getComponent(j).getAssociationParams();

    if (!assoc_i.hasAssociation() || !assoc_j.hasAssociation()) {
        return T(0.0);
    }

    // Combining rules for cross-association
    double epsilon_AB = 0.5 * (assoc_i.epsilon_AB + assoc_j.epsilon_AB);
    double kappa_AB = std::sqrt(assoc_i.kappa_AB * assoc_j.kappa_AB);
    double sigma_ij = 0.5 * (sigma_[i] + sigma_[j]);

    T d_ij = 0.5 * (d[i] + d[j]);

    // Hard-sphere radial distribution function at contact
    T g_hs = calc_g_hs_templated(eta);

    T Delta = g_hs * (exp(epsilon_AB / T_val) - 1.0) * kappa_AB * sigma_ij * sigma_ij * sigma_ij;

    return Delta;
}

template<typename T>
std::vector<T> PCSaftEOS::solve_association_templated(const std::vector<T>& x, T T_val, T rho) const {
    // Count total number of association sites
    int total_sites = 0;
    for (int i = 0; i < nc_; ++i) {
        total_sites += mixture_.getComponent(i).getAssociationParams().num_sites;
    }

    if (total_sites == 0) {
        return std::vector<T>();
    }

    // Convert density from mol/m³ to molecules/Å³
    T rho_reduced = rho * T(N_A * 1e-30);

    // Initialize X_A = 1 (all sites non-bonded initially)
    std::vector<T> X_A(total_sites, T(1.0));
    std::vector<T> X_A_new(total_sites);

    std::vector<T> d = calc_d_templated(T_val);
    T eta = calc_eta_templated(x, T_val, rho);

    // Calculate association strengths Δ^AB for all site pairs
    std::vector<std::vector<T>> Delta(nc_, std::vector<T>(nc_, T(0.0)));
    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            if (mixture_.getComponent(i).hasAssociation() &&
                mixture_.getComponent(j).hasAssociation()) {
                Delta[i][j] = calc_delta_assoc_templated(i, j, T_val, d, eta);
            }
        }
    }

    // Successive substitution to solve X_A equations
    const int max_iter = 100;
    const double tol = 1e-10;

    for (int iter = 0; iter < max_iter; ++iter) {
        int site_idx = 0;
        for (int i = 0; i < nc_; ++i) {
            int num_sites_i = mixture_.getComponent(i).getAssociationParams().num_sites;

            for (int A = 0; A < num_sites_i; ++A) {
                // X_Ai = 1 / (1 + ρ Σ_j x_j Σ_B X_Bj Δ^AB_ij)
                T sum = T(0.0);

                int site_idx_j = 0;
                for (int j = 0; j < nc_; ++j) {
                    int num_sites_j = mixture_.getComponent(j).getAssociationParams().num_sites;

                    for (int B = 0; B < num_sites_j; ++B) {
                        sum += x[j] * X_A[site_idx_j++] * Delta[i][j];
                    }
                }

                X_A_new[site_idx++] = T(1.0) / (T(1.0) + rho_reduced * sum);
            }
        }

        // Check convergence
        double max_diff = 0.0;
        for (int k = 0; k < total_sites; ++k) {
            max_diff = std::max(max_diff, std::abs(val(X_A_new[k] - X_A[k])));
        }

        X_A = X_A_new;

        if (max_diff < tol) {
            break;
        }
    }

    return X_A;
}

template<typename T>
T PCSaftEOS::calc_a_assoc_templated(const std::vector<T>& x, T T_val, T rho) const {
    if (!mixture_.hasAssociation()) {
        return T(0.0);
    }

    // Solve for fraction of non-bonded sites X_A
    std::vector<T> X_A = solve_association_templated(x, T_val, rho);

    // Calculate association contribution
    // a_assoc = Σ_i x_i Σ_A (ln(X_Ai) - X_Ai/2 + 1/2)
    T a_assoc = T(0.0);

    int site_idx = 0;
    for (int i = 0; i < nc_; ++i) {
        const AssociationParams& assoc = mixture_.getComponent(i).getAssociationParams();

        if (assoc.hasAssociation()) {
            int num_sites = assoc.num_sites;

            for (int A = 0; A < num_sites; ++A) {
                a_assoc += x[i] * (log(X_A[site_idx]) - X_A[site_idx] / T(2.0) + T(0.5));
                site_idx++;
            }
        }
    }

    return a_assoc;
}

std::vector<double> PCSaftEOS::calc_da_assoc_dx_autodiff(double T, double rho) const {
    using namespace autodiff;

    const std::vector<double>& x_double = mixture_.getMoleFractions();
    std::vector<double> da_dx_unc(nc_);  // Unconstrained derivatives

    // Step 1: Compute unconstrained derivatives for all components
    for (int k = 0; k < nc_; ++k) {
        // Create dual number composition vector
        std::vector<dual> x_dual(nc_);
        for (int i = 0; i < nc_; ++i) {
            x_dual[i] = x_double[i];
        }

        // Define function: a_assoc as a function of x[k]
        auto f = [&](dual xk) -> dual {
            std::vector<dual> x_temp = x_dual;
            x_temp[k] = xk;
            return calc_a_assoc_templated(x_temp, dual(T), dual(rho));
        };

        // Compute unconstrained derivative using autodiff
        dual xk = x_double[k];
        da_dx_unc[k] = derivative(f, wrt(xk), at(xk));
    }

    // Step 2: Apply correct constraint transformation for composition derivatives
    // For constrained derivatives respecting Σx_i = 1:
    // ∂a/∂x_k|_const = ∂a/∂x_k|_unc - Σ_j x_j (∂a/∂x_j|_unc)
    // This ensures Gibbs-Duhem relation: Σ_k x_k (∂a/∂x_k|_const) = 0
    std::vector<double> da_dx(nc_);

    double sum_x_dadx_unc = 0.0;
    for (int j = 0; j < nc_; ++j) {
        sum_x_dadx_unc += x_double[j] * da_dx_unc[j];
    }

    for (int k = 0; k < nc_; ++k) {
        da_dx[k] = da_dx_unc[k] - sum_x_dadx_unc;
    }

    return da_dx;
}

} // namespace pcsaft
