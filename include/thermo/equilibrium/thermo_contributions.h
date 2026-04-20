/**
 * @file thermo_contributions.h
 * @brief Shared ideal-gas and residual thermodynamic contributions used by flash solvers.
 */

#ifndef THERMO_EQUILIBRIUM_THERMO_CONTRIBUTIONS_H
#define THERMO_EQUILIBRIUM_THERMO_CONTRIBUTIONS_H

#include "thermo/core/mixture.h"
#include "thermo/eos.h"

#include <vector>

namespace DMThermo {
namespace Equilibrium {
namespace Contributions {

struct IdealProps {
    double h_ig = 0.0; // J/mol (relative to T_ref)
    double s_ig = 0.0; // J/mol/K (relative to T_ref, includes mixing and -R ln(P/P_ref))
};

struct ResidualProps {
    double h_res = 0.0; // J/mol
    double s_res = 0.0; // J/mol/K (departure from ideal gas at same T,P)
    double g_res = 0.0; // J/mol (departure from ideal gas at same T,P)
};

IdealProps idealGasProps(const Core::Mixture& mix, double T, double P, const std::vector<double>& x);
ResidualProps residualProps(const EOS& eos, double T, double rho, const std::vector<double>& x, double P_target = 0.0);

} // namespace Contributions
} // namespace Equilibrium
} // namespace DMThermo

#endif // THERMO_EQUILIBRIUM_THERMO_CONTRIBUTIONS_H
