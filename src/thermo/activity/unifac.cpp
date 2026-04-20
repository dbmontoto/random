#include "thermo/activity/unifac.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace DMThermo {
namespace Activity {

namespace {

struct SubgroupTotals {
    double total_count = 0.0;
    double total_QX = 0.0;
};

double interactionA(const Data::UnifacTables& tables, int main_i, int main_j, double T) {
    (void)T;
    if (main_i <= 0 || main_j <= 0) return 0.0;
    const long long key = (static_cast<long long>(main_i) << 32) | static_cast<unsigned long long>(main_j);
    auto it = tables.interactions.find(key);
    if (it == tables.interactions.end()) return 0.0;
    const auto& p = it->second;
    return p.a + p.b * T + p.c * T * T;
}

double tau(const Data::UnifacTables& tables, int main_i, int main_j, double T) {
    const double A = interactionA(tables, main_i, main_j, T);
    return std::exp(-A / T);
}

struct GroupInfo {
    int subgroup_id = 0;
    int main_group_id = 0;
    double Q = 0.0;
};

} // namespace

std::vector<double> lnGammaUNIFAC(double T, const Core::Mixture& mixture, const std::vector<double>& x) {
    return [&]() {
        std::vector<double> ln_gamma;
        (void)excessGibbsOverRTUNIFAC(T, mixture, x, &ln_gamma);
        return ln_gamma;
    }();
}

double excessGibbsOverRTUNIFAC(
    double T,
    const Core::Mixture& mixture,
    const std::vector<double>& x,
    std::vector<double>* ln_gamma_out)
{
    if (!(std::isfinite(T) && T > 0.0)) {
        throw std::invalid_argument("UNIFAC: invalid temperature");
    }
    if (!mixture.isValidComposition(x)) {
        throw std::invalid_argument("UNIFAC: invalid composition vector");
    }
    const auto tables_ptr = mixture.unifacTables();
    if (!tables_ptr) {
        throw std::runtime_error("UNIFAC: mixture has no UNIFAC tables");
    }
    const auto& tables = *tables_ptr;

    const int nc = mixture.numComponents();
    std::vector<double> ln_gamma(nc, 0.0);

    // Component r_i, q_i from subgroup sums.
    std::vector<double> r(nc, 0.0), q(nc, 0.0);
    std::vector<int> total_groups_count(nc, 0);
    for (int i = 0; i < nc; ++i) {
        const auto& comp = mixture.component(i);
        for (const auto& g : comp.unifacSubgroups()) {
            auto it = tables.subgroups.find(g.subgroup_id);
            if (it == tables.subgroups.end()) {
                throw std::runtime_error("UNIFAC: missing subgroup_id=" + std::to_string(g.subgroup_id));
            }
            r[i] += static_cast<double>(g.count) * it->second.R;
            q[i] += static_cast<double>(g.count) * it->second.Q;
            total_groups_count[i] += g.count;
        }
    }

    // Combinatorial part (original UNIFAC form).
    const double z = 10.0;
    double r_mix = 0.0, q_mix = 0.0;
    for (int i = 0; i < nc; ++i) {
        r_mix += x[i] * r[i];
        q_mix += x[i] * q[i];
    }
    if (!(std::isfinite(r_mix) && std::isfinite(q_mix)) || r_mix <= 0.0 || q_mix <= 0.0) {
        throw std::runtime_error("UNIFAC: invalid r_mix/q_mix");
    }

    std::vector<double> V(nc, 0.0), F(nc, 0.0), ln_gamma_C(nc, 0.0);
    for (int i = 0; i < nc; ++i) {
        V[i] = r[i] / r_mix;
        F[i] = q[i] / q_mix;
        const double term = V[i] / std::max(F[i], 1e-300);
        ln_gamma_C[i] = 1.0 - V[i] + std::log(std::max(V[i], 1e-300))
                      - (z / 2.0) * q[i] * (1.0 - term + std::log(std::max(term, 1e-300)));
    }

    // Collect distinct subgroups present in the mixture.
    std::vector<GroupInfo> groups;
    groups.reserve(64);
    std::unordered_map<int, std::size_t> subgroup_index;
    subgroup_index.reserve(64);

    auto addGroup = [&](int subgroup_id) {
        auto it = subgroup_index.find(subgroup_id);
        if (it != subgroup_index.end()) return it->second;
        auto jt = tables.subgroups.find(subgroup_id);
        if (jt == tables.subgroups.end()) {
            throw std::runtime_error("UNIFAC: missing subgroup_id=" + std::to_string(subgroup_id));
        }
        GroupInfo gi;
        gi.subgroup_id = subgroup_id;
        gi.main_group_id = jt->second.main_group_id;
        gi.Q = jt->second.Q;
        const std::size_t idx = groups.size();
        groups.push_back(gi);
        subgroup_index[subgroup_id] = idx;
        return idx;
    };

    // Mixture subgroup fractions X_k and theta_k.
    std::vector<double> Xg;
    Xg.assign(0, 0.0);
    for (int i = 0; i < nc; ++i) {
        for (const auto& g : mixture.component(i).unifacSubgroups()) {
            (void)addGroup(g.subgroup_id);
        }
    }
    Xg.assign(groups.size(), 0.0);

    double total_groups_moles = 0.0;
    for (int i = 0; i < nc; ++i) {
        total_groups_moles += x[i] * static_cast<double>(total_groups_count[i]);
        for (const auto& g : mixture.component(i).unifacSubgroups()) {
            const std::size_t k = subgroup_index[g.subgroup_id];
            Xg[k] += x[i] * static_cast<double>(g.count);
        }
    }
    if (!(std::isfinite(total_groups_moles) && total_groups_moles > 0.0)) {
        throw std::runtime_error("UNIFAC: total subgroup count is zero");
    }
    for (double& v : Xg) v /= total_groups_moles;

    std::vector<double> theta(groups.size(), 0.0);
    double sum_QX = 0.0;
    for (std::size_t k = 0; k < groups.size(); ++k) {
        sum_QX += groups[k].Q * Xg[k];
    }
    if (!(std::isfinite(sum_QX) && sum_QX > 0.0)) {
        throw std::runtime_error("UNIFAC: invalid sum_QX");
    }
    for (std::size_t k = 0; k < groups.size(); ++k) {
        theta[k] = (groups[k].Q * Xg[k]) / sum_QX;
    }

    auto lnGammaGroup = [&](const std::vector<double>& theta_in) -> std::vector<double> {
        const std::size_t ng = groups.size();
        std::vector<double> S(ng, 0.0);
        for (std::size_t k = 0; k < ng; ++k) {
            double s = 0.0;
            for (std::size_t m = 0; m < ng; ++m) {
                s += theta_in[m] * tau(tables, groups[m].main_group_id, groups[k].main_group_id, T);
            }
            S[k] = s;
        }

        std::vector<double> lnG(ng, 0.0);
        for (std::size_t k = 0; k < ng; ++k) {
            const double Sk = std::max(S[k], 1e-300);
            double sum_term = 0.0;
            for (std::size_t m = 0; m < ng; ++m) {
                const double Sm = std::max(S[m], 1e-300);
                sum_term += theta_in[m] * tau(tables, groups[k].main_group_id, groups[m].main_group_id, T) / Sm;
            }
            lnG[k] = groups[k].Q * (1.0 - std::log(Sk) - sum_term);
        }
        return lnG;
    };

    const auto lnG_mix = lnGammaGroup(theta);

    // Residual part per component.
    std::vector<double> ln_gamma_R(nc, 0.0);
    for (int i = 0; i < nc; ++i) {
        // theta_k for pure component i.
        std::vector<double> theta_i(groups.size(), 0.0);
        double sum_QX_i = 0.0;
        const double total_cnt = std::max(1.0, static_cast<double>(total_groups_count[i]));
        for (const auto& g : mixture.component(i).unifacSubgroups()) {
            const std::size_t k = subgroup_index[g.subgroup_id];
            const double Xki = static_cast<double>(g.count) / total_cnt;
            sum_QX_i += groups[k].Q * Xki;
        }
        if (!(std::isfinite(sum_QX_i) && sum_QX_i > 0.0)) {
            throw std::runtime_error("UNIFAC: invalid sum_QX_i for component " + std::to_string(i));
        }
        for (const auto& g : mixture.component(i).unifacSubgroups()) {
            const std::size_t k = subgroup_index[g.subgroup_id];
            const double Xki = static_cast<double>(g.count) / total_cnt;
            theta_i[k] = (groups[k].Q * Xki) / sum_QX_i;
        }

        const auto lnG_i = lnGammaGroup(theta_i);

        double sum = 0.0;
        for (const auto& g : mixture.component(i).unifacSubgroups()) {
            const std::size_t k = subgroup_index[g.subgroup_id];
            sum += static_cast<double>(g.count) * (lnG_mix[k] - lnG_i[k]);
        }
        ln_gamma_R[i] = sum;
    }

    for (int i = 0; i < nc; ++i) {
        ln_gamma[i] = ln_gamma_C[i] + ln_gamma_R[i];
    }

    if (ln_gamma_out) *ln_gamma_out = ln_gamma;

    double GE_over_RT = 0.0;
    for (int i = 0; i < nc; ++i) {
        GE_over_RT += x[i] * ln_gamma[i];
    }
    return GE_over_RT;
}

} // namespace Activity
} // namespace DMThermo
