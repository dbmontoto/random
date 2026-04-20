#ifndef PCSAFT_H
#define PCSAFT_H

#include "thermo/pcsaft/core/mixture.h"
#include <vector>
#include <memory>

namespace pcsaft {

/**
 * @brief Structure to hold PC-SAFT calculation results
 */
struct PCSaftResults {
    // State variables
    double temperature;          // [K]
    double density;              // Total molar density [mol/m³]
    double pressure;             // [Pa]

    // Helmholtz energy contributions [dimensionless, A/(nRT)]
    double a_res;                // Total residual Helmholtz energy
    double a_hc;                 // Hard-chain contribution
    double a_disp;               // Dispersion contribution
    double a_assoc;              // Association contribution

    // Compressibility factor
    double Z;                    // Z = PV/(nRT)

    // Chemical potentials and fugacity coefficients
    std::vector<double> mu_res;  // Residual chemical potentials [dimensionless]
    std::vector<double> fugacity_coeff; // Fugacity coefficients [-]

    // Derivatives (for property calculations)
    double dadt;                 // ∂(a_res)/∂T at constant ρ, x
    double dadrho;               // ∂(a_res)/∂ρ at constant T, x
    double d2adt2;               // ∂²(a_res)/∂T²
    double d2adrho2;             // ∂²(a_res)/∂ρ²
    double d2adtdrho;            // ∂²(a_res)/∂T∂ρ

    PCSaftResults() : temperature(0), density(0), pressure(0),
                      a_res(0), a_hc(0), a_disp(0), a_assoc(0), Z(0),
                      dadt(0), dadrho(0), d2adt2(0), d2adrho2(0), d2adtdrho(0) {}
};

/**
 * @brief PC-SAFT Equation of State with association
 *
 * Implements the full PC-SAFT equation including:
 * - Hard-chain contribution
 * - Dispersion contribution
 * - Association contribution
 * - All necessary derivatives for thermodynamic properties
 */
class PCSaftEOS {
public:
    // Constructor
    explicit PCSaftEOS(const Mixture& mixture);

    // Main calculation methods
    /**
     * @brief Calculate all PC-SAFT properties at given T and ρ
     * @param T Temperature [K]
     * @param rho Total molar density [mol/m³]
     * @param calc_derivatives Calculate second derivatives (slower)
     * @return PCSaftResults structure with all calculated properties
     */
    PCSaftResults calculate(double T, double rho, bool calc_derivatives = true);

    /**
     * @brief Calculate pressure at given T and ρ
     * @param T Temperature [K]
     * @param rho Total molar density [mol/m³]
     * @return Pressure [Pa]
     */
    double calculatePressure(double T, double rho) const;

    /**
     * @brief Calculate compressibility factor at given T and ρ
     * @param T Temperature [K]
     * @param rho Total molar density [mol/m³]
     * @return Z = PV/(nRT)
     */
    double calculateZ(double T, double rho) const;

    /**
     * @brief Calculate fugacity coefficients at given T and ρ
     * @param T Temperature [K]
     * @param rho Total molar density [mol/m³]
     * @return Vector of fugacity coefficients for each component
     */
    std::vector<double> calculateFugacityCoefficients(double T, double rho) const;

    // Getters
    const Mixture& getMixture() const { return mixture_; }

    // ========================================================================
    // PUBLIC for testing/validation - Core contributions and derivatives
    // ========================================================================

    /**
     * @brief Calculate hard-chain contribution to Helmholtz energy
     * @param T Temperature [K]
     * @param rho Molar density [mol/m³]
     * @return a_hc [dimensionless]
     */
    double calc_a_hc(double T, double rho) const;

    /**
     * @brief Calculate dispersion contribution to Helmholtz energy
     * @param T Temperature [K]
     * @param rho Molar density [mol/m³]
     * @return a_disp [dimensionless]
     */
    double calc_a_disp(double T, double rho) const;

    // ========================================================================
    // ANALYTICAL DERIVATIVES - VALIDATION ONLY
    // ========================================================================
    // NOTE: These analytical derivative implementations contain fundamental errors
    // (sign errors, factor-of-2 errors, magnitude errors) as documented in
    // docs/debugging/ANALYTICAL_DERIVATIVES_STATUS.md
    //
    // PRODUCTION CODE USES: calc_da_hc_dx_autodiff(), calc_da_disp_dx_autodiff()
    //
    // These functions are KEPT FOR VALIDATION PURPOSES ONLY to allow comparison
    // between analytical formulas from literature and machine-precision autodiff.
    // DO NOT USE IN PRODUCTION CALCULATIONS.
    // ========================================================================

    /**
     * @brief [VALIDATION ONLY] Calculate analytical composition derivative of hard-chain term
     * ∂a_hc/∂x_i at constant T, ρ
     *
     * @warning Contains fundamental errors - use calc_da_hc_dx_autodiff() instead
     */
    std::vector<double> calc_da_hc_dx(double T, double rho) const;

    /**
     * @brief [VALIDATION ONLY] Calculate analytical composition derivative of dispersion term
     * ∂a_disp/∂x_i at constant T, ρ
     *
     * @warning Contains fundamental errors - use calc_da_disp_dx_autodiff() instead
     */
    std::vector<double> calc_da_disp_dx(double T, double rho) const;

    /**
     * @brief [VALIDATION ONLY] Calculate analytical composition derivative of association term
     * ∂a_assoc/∂x_i at constant T, ρ
     *
     * @warning Contains fundamental errors - use autodiff version when available
     */
    std::vector<double> calc_da_assoc_dx(double T, double rho) const;

    /**
     * @brief Numerical association composition derivative (TEMPORARY FIX)
     * Uses finite differences instead of buggy analytical version
     */
    std::vector<double> calc_da_assoc_dx_numerical(double T, double rho) const;

    // ========================================================================
    // AUTODIFF-based derivatives (Phase 1: parallel implementation)
    // ========================================================================

    /**
     * @brief Calculate hard-chain composition derivative using automatic differentiation
     * ∂a_hc/∂x_i at constant T, ρ
     *
     * This provides a reference implementation using autodiff library to validate
     * the analytical derivatives. Should give machine-precision accurate results.
     */
    std::vector<double> calc_da_hc_dx_autodiff(double T, double rho) const;

    /**
     * @brief Calculate dispersion composition derivative using automatic differentiation
     * ∂a_disp/∂x_i at constant T, ρ
     */
    std::vector<double> calc_da_disp_dx_autodiff(double T, double rho) const;

    /**
     * @brief Calculate association composition derivative using automatic differentiation
     * ∂a_assoc/∂x_i at constant T, ρ (PRODUCTION VERSION)
     */
    std::vector<double> calc_da_assoc_dx_autodiff(double T, double rho) const;

    /**
     * @brief Public method to get packing fraction (eta) for analysis
     * @param T Temperature [K]
     * @param rho Density [mol/m³]
     * @return Packing fraction η
     */
    double getPackingFraction(double T, double rho) const {
        return calc_eta(T, rho);
    }

private:
    Mixture mixture_;
    int nc_;  // Number of components

    // Mixture-averaged parameters (calculated once in constructor)
    std::vector<double> m_;       // Segment numbers
    std::vector<double> sigma_;   // Segment diameters [Angstrom]
    std::vector<double> epsilon_; // Dispersion energies [K]

    // ========================================================================
    // Core PC-SAFT contribution calculations
    // ========================================================================

    /**
     * @brief Calculate hard-sphere contribution
     * Used as building block for hard-chain term
     */
    double calc_a_hs(double eta) const;

    /**
     * @brief Calculate hard-sphere Helmholtz energy for mixtures (Mansoori et al.)
     * @param zeta Packing fraction vector [ζ_0, ζ_1, ζ_2, ζ_3]
     * @return a_hs [dimensionless]
     * Full mixture formula from Equation 8 of Gross & Sadowski (2001)
     */
    double calc_a_hs_mixture(const std::vector<double>& zeta) const;

    /**
     * @brief Calculate association contribution to Helmholtz energy
     * @param T Temperature [K]
     * @param rho Molar density [mol/m³]
     * @return a_assoc [dimensionless]
     */
    double calc_a_assoc(double T, double rho) const;

    // ========================================================================
    // Derivative calculations
    // ========================================================================

    /**
     * @brief Calculate ∂(a_res)/∂ρ at constant T
     */
    double calc_dadrho(double T, double rho) const;

    /**
     * @brief Calculate ∂(a_res)/∂T at constant ρ
     */
    double calc_dadt(double T, double rho) const;

    /**
     * @brief Calculate ∂²(a_res)/∂ρ²
     */
    double calc_d2adrho2(double T, double rho) const;

    /**
     * @brief Calculate ∂²(a_res)/∂T²
     */
    double calc_d2adt2(double T, double rho) const;

    /**
     * @brief Calculate ∂²(a_res)/∂T∂ρ
     */
    double calc_d2adtdrho(double T, double rho) const;

    /**
     * @brief Calculate residual chemical potentials (for fugacity)
     * Returns ln(φ_i) for each component
     */
    std::vector<double> calc_mu_res(double T, double rho) const;

    // ========================================================================
    // Helper functions
    // ========================================================================

    /**
     * @brief Calculate reduced density η = π/6 * ρ * d³
     * where d is temperature-dependent hard-sphere diameter
     */
    double calc_eta(double T, double rho) const;

    /**
     * @brief Calculate temperature-dependent hard-sphere diameter
     * d = σ [1 - 0.12 exp(-3ε/kT)]
     */
    std::vector<double> calc_d(double T) const;

    /**
     * @brief Calculate mean segment number
     */
    double calc_m_bar() const;

    /**
     * @brief Calculate power series integrals I1 and I2 for dispersion
     */
    void calc_dispersion_integrals(double eta, double& I1, double& I2) const;

    /**
     * @brief Solve association equations for fraction of non-bonded sites
     * Returns X_A for each association site
     */
    std::vector<double> solve_association(double T, double rho) const;

    /**
     * @brief Calculate association strength Δ^AB
     */
    double calc_delta_assoc(int i, int j, double T, const std::vector<double>& d, double rho) const;

    /**
     * @brief Calculate radial distribution function g_hs at contact (simplified)
     */
    double calc_g_hs(double eta) const;

    /**
     * @brief Calculate component-specific radial distribution function g_ii^hs
     * @param d_i Diameter of component i
     * @param zeta_2 Packing fraction ζ_2
     * @param zeta_3 Packing fraction ζ_3
     * @return g_ii^hs value
     */
    double calc_g_ii(double d_i, double zeta_2, double zeta_3) const;

    /**
     * @brief Calculate compressibility factor C_1 (Equation A.11 from Gross & Sadowski)
     */
    double calc_C1(double eta) const;

    /**
     * @brief Calculate derivative of C_1 with respect to eta (Equation A.31)
     */
    double calc_C2(double eta) const;

    /**
     * @brief Calculate derivative of dispersion integrals I_1 or I_2 with respect to composition
     * @param eta Packing fraction
     * @param m_bar Mean segment number
     * @param m_k Segment number of component k
     * @param zeta_3_xk Derivative of ζ_3 with respect to x_k
     * @param order 1 for I_1 derivative, 2 for I_2 derivative
     * @return I_{order,xk} derivative value
     */
    double calc_I_derivative(double eta, double m_bar, double m_k, double zeta_3_xk, int order) const;

    // ========================================================================
    // Templated helper functions for autodiff
    // ========================================================================

    /**
     * @brief Templated version of calc_d - works with dual numbers
     */
    template<typename T>
    std::vector<T> calc_d_templated(T T_val) const;

    /**
     * @brief Templated version of calc_m_bar - works with dual numbers
     */
    template<typename T>
    T calc_m_bar_templated(const std::vector<T>& x) const;

    /**
     * @brief Templated version of calc_g_ii - works with dual numbers
     */
    template<typename T>
    T calc_g_ii_templated(T d_i, T zeta_2, T zeta_3) const;

    /**
     * @brief Templated version of calc_a_hs_mixture - works with dual numbers
     */
    template<typename T>
    T calc_a_hs_mixture_templated(const std::vector<T>& zeta) const;

    /**
     * @brief Templated version of calc_a_hc - works with dual numbers
     */
    template<typename T>
    T calc_a_hc_templated(const std::vector<T>& x, T T_val, T rho) const;

    /**
     * @brief Templated version of calc_eta - works with dual numbers
     */
    template<typename T>
    T calc_eta_templated(const std::vector<T>& x, T T_val, T rho) const;

    /**
     * @brief Templated version of calc_g_hs - works with dual numbers
     */
    template<typename T>
    T calc_g_hs_templated(T eta) const;

    /**
     * @brief Templated version of calc_dispersion_integrals - works with dual numbers
     */
    template<typename T>
    void calc_dispersion_integrals_templated(T eta, T m_bar, T& I1, T& I2) const;

    /**
     * @brief Templated version of calc_a_disp - works with dual numbers
     */
    template<typename T>
    T calc_a_disp_templated(const std::vector<T>& x, T T_val, T rho) const;

    /**
     * @brief Templated version of association functions - works with dual numbers
     */
    template<typename T>
    T calc_delta_assoc_templated(int i, int j, T T_val, const std::vector<T>& d, T eta) const;

    template<typename T>
    std::vector<T> solve_association_templated(const std::vector<T>& x, T T_val, T rho) const;

    template<typename T>
    T calc_a_assoc_templated(const std::vector<T>& x, T T_val, T rho) const;
};

} // namespace pcsaft

#endif // PCSAFT_H
