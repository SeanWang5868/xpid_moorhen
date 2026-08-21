#pragma once

#include <cmath>
#include <vector>

#include <gemmi/model.hpp>

#include "geometry.h"
#include "rotatable_groups.h"

namespace xhpi {

struct HydrogenConformer {
    double phi = 0.0;
    std::vector<gemmi::Position> hydrogen_positions;
};

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

}  // namespace xhpi
