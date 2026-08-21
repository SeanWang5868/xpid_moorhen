#pragma once

#include <cmath>
#include <optional>
#include <tuple>
#include <vector>

#include <gemmi/model.hpp>

#include "geometry.h"
#include "cone_acceptors.h"
#include "rotatable_groups.h"

namespace xhpi {

struct HydrogenConformer {
    double phi = 0.0;
    std::vector<gemmi::Position> hydrogen_positions;
};

struct ConePositiveEvidence {
    double phi = 0.0;
    size_t hydrogen_index = 0;
    gemmi::Position hydrogen_position;
    double theta = -1.0;
    double xh_pi_angle = -1.0;
    int is_hudson = 0;
    int is_plevin = 0;
};

inline constexpr double CONE_H_VDW_RADIUS = 1.20;
inline constexpr double CONE_DEFAULT_VDW_RADIUS = 1.70;
inline constexpr double CONE_CLASH_SCALE = 0.75;
inline constexpr double CONE_ABSOLUTE_MIN_HA = 1.30;
inline constexpr double CONE_HBOND_CONTACT_HA_MAX = 2.80;
inline constexpr double CONE_HBOND_CONTACT_DHA_MIN = 120.0;
inline constexpr double CONE_ACCEPTOR_MIN_OCCUPANCY = 0.20;

struct ConeEnvironmentAtom {
    gemmi::Position position;
    const gemmi::Atom* atom = nullptr;
    const gemmi::Residue* residue = nullptr;
    double vdw_radius = CONE_DEFAULT_VDW_RADIUS;
    bool is_hbond_acceptor = false;
    double occupancy = 1.0;
};

inline double cone_vdw_radius(const std::string& element) {
    if (element == "H" || element == "D") return 1.20;
    if (element == "C") return 1.70;
    if (element == "N") return 1.55;
    if (element == "O") return 1.52;
    if (element == "S" || element == "P") return 1.80;
    if (element == "SE") return 1.90;
    return CONE_DEFAULT_VDW_RADIUS;
}

inline ConeEnvironmentAtom make_cone_environment_atom(
        const gemmi::Position& position,
        const gemmi::Atom& atom,
        const gemmi::Residue& residue) {
    const std::string element = gemmi::to_upper(std::string(atom.element.name()));
    return {position, &atom, &residue, cone_vdw_radius(element),
            is_cone_hbond_acceptor(residue, atom), atom.occ};
}

inline double cone_dha_angle(const gemmi::Position& x_pos,
                             const gemmi::Position& h_pos,
                             const gemmi::Position& acceptor_pos) {
    return calculate_angle_vectors(x_pos - h_pos, acceptor_pos - h_pos);
}

inline bool is_valid_cone_hbond_contact(const gemmi::Position& x_pos,
                                        const gemmi::Position& h_pos,
                                        const ConeEnvironmentAtom& environment) {
    if (!environment.is_hbond_acceptor ||
        environment.occupancy < CONE_ACCEPTOR_MIN_OCCUPANCY) {
        return false;
    }
    const double distance = h_pos.dist(environment.position);
    if (distance < CONE_ABSOLUTE_MIN_HA || distance > CONE_HBOND_CONTACT_HA_MAX) {
        return false;
    }
    const double angle = cone_dha_angle(x_pos, h_pos, environment.position);
    return angle >= CONE_HBOND_CONTACT_DHA_MIN;
}

inline bool cone_hydrogen_is_sterically_valid(
        const gemmi::Position& x_pos,
        const gemmi::Position& h_pos,
        const std::vector<ConeEnvironmentAtom>& environment) {
    for (const ConeEnvironmentAtom& env : environment) {
        const double distance = h_pos.dist(env.position);
        if (is_valid_cone_hbond_contact(x_pos, h_pos, env)) continue;
        if (distance < CONE_ABSOLUTE_MIN_HA) return false;
        const double clash_distance = CONE_CLASH_SCALE *
            (CONE_H_VDW_RADIUS + env.vdw_radius);
        if (distance < clash_distance) return false;
    }
    return true;
}

inline bool cone_conformer_is_sterically_valid(
        const HydrogenConformer& conformer,
        const gemmi::Position& x_pos,
        const std::vector<ConeEnvironmentAtom>& environment) {
    for (const gemmi::Position& h_pos : conformer.hydrogen_positions) {
        if (!cone_hydrogen_is_sterically_valid(x_pos, h_pos, environment)) return false;
    }
    return true;
}

inline std::vector<HydrogenConformer> generate_group_conformers(
        const gemmi::Position& parent_pos,
        const gemmi::Position& x_pos,
        const RotatableGroupDefinition& definition,
        int step_degrees = 1) {
    std::vector<HydrogenConformer> conformers;
    if (step_degrees <= 0) return conformers;

    gemmi::Position axis = x_pos - parent_pos;
    gemmi::Position u = normalize(axis);
    if (magnitude(u) < 1e-8) return conformers;

    gemmi::Position reference(1.0, 0.0, 0.0);
    if (std::abs(dot_product(u, reference)) > 0.99) {
        reference = gemmi::Position(0.0, 1.0, 0.0);
    }
    gemmi::Position v = normalize(cross_product(u, reference));
    if (magnitude(v) < 1e-8) return conformers;
    gemmi::Position w = cross_product(u, v);

    const double beta = (180.0 - definition.parent_x_h_angle) * M_PI / 180.0;
    const double axial = definition.xh_bond_length * std::cos(beta);
    const double radial = definition.xh_bond_length * std::sin(beta);

    for (int phi = 0; phi < definition.rotation_period; phi += step_degrees) {
        HydrogenConformer conformer;
        conformer.phi = static_cast<double>(phi);
        conformer.hydrogen_positions.reserve(
            static_cast<size_t>(definition.hydrogen_count));
        for (int h_index = 0; h_index < definition.hydrogen_count; ++h_index) {
            const double angle_deg = static_cast<double>(phi + 120 * h_index);
            const double angle = angle_deg * M_PI / 180.0;
            conformer.hydrogen_positions.push_back(
                x_pos + u * axial + v * (radial * std::cos(angle)) +
                w * (radial * std::sin(angle)));
        }
        conformers.push_back(std::move(conformer));
    }
    return conformers;
}

inline auto cone_evidence_rank(const ConePositiveEvidence& evidence) {
    const int positive_systems = evidence.is_hudson + evidence.is_plevin;
    const double theta_for_rank = evidence.theta >= 0.0 ? evidence.theta : 999.0;
    return std::make_tuple(
        evidence.is_hudson && evidence.is_plevin ? 1 : 0,
        positive_systems,
        -theta_for_rank,
        evidence.xh_pi_angle,
        -evidence.phi,
        -static_cast<long long>(evidence.hydrogen_index));
}

inline std::optional<ConePositiveEvidence> evaluate_binary_cone(
        const gemmi::Position& parent_pos,
        const gemmi::Position& x_pos,
        const RotatableGroupDefinition& definition,
        const std::vector<ConeEnvironmentAtom>& environment,
        const gemmi::Position& pi_center,
        const gemmi::Position& pi_normal,
        double dist_x_pi,
        double max_dist,
        double xpcn_angle,
        double proj_dist,
        size_t ring_size,
        int step_degrees = 1) {
    std::optional<ConePositiveEvidence> best;
    const double ring_radius = ring_size == 6 ? 2.0 : 1.6;
    for (const HydrogenConformer& conformer :
         generate_group_conformers(parent_pos, x_pos, definition, step_degrees)) {
        if (!cone_conformer_is_sterically_valid(conformer, x_pos, environment)) continue;
        for (size_t h_index = 0; h_index < conformer.hydrogen_positions.size(); ++h_index) {
            const gemmi::Position& h_pos = conformer.hydrogen_positions[h_index];
            const double xh_pi_angle = calculate_xh_picenter_angle(x_pos, h_pos, pi_center);
            const double theta = calculate_hudson_theta(pi_center, x_pos, h_pos, pi_normal);
            const int is_plevin =
                dist_x_pi < max_dist && xpcn_angle >= 0.0 && xpcn_angle < 25.0 &&
                xh_pi_angle >= 120.0 ? 1 : 0;
            const int is_hudson =
                dist_x_pi <= max_dist && std::isfinite(proj_dist) &&
                proj_dist <= ring_radius && theta >= 0.0 && theta <= 40.0 ? 1 : 0;
            if (!is_plevin && !is_hudson) continue;

            ConePositiveEvidence evidence{
                conformer.phi, h_index, h_pos, theta, xh_pi_angle,
                is_hudson, is_plevin};
            if (!best || cone_evidence_rank(evidence) > cone_evidence_rank(*best)) {
                best = evidence;
            }
        }
    }
    return best;
}

}  // namespace xhpi
