#pragma once

#include <string>
#include <unordered_map>

namespace xhpi {

enum class RotatableGroupKind {
    SingleHydrogen,
    Methyl
};

struct RotatableGroupDefinition {
    std::string parent_atom_name;
    RotatableGroupKind kind = RotatableGroupKind::SingleHydrogen;
    std::string element;
    int hydrogen_count = 1;
    double xh_bond_length = 1.0;
    double parent_x_h_angle = 109.5;
    int rotation_period = 360;
};

inline RotatableGroupDefinition single_hydrogen_group(
        const std::string& parent,
        const std::string& element,
        double length,
        double angle) {
    return {parent, RotatableGroupKind::SingleHydrogen, element, 1,
            length, angle, 360};
}

inline RotatableGroupDefinition methyl_group(
        const std::string& parent,
        double length,
        double angle) {
    return {parent, RotatableGroupKind::Methyl, "C", 3,
            length, angle, 120};
}

// Nuclear X-H distances and parent-X-H angles are the frozen production
// geometry used by XPID. Cation-pi groups are deliberately absent.
inline const std::unordered_map<std::string, RotatableGroupDefinition> ROTATABLE_GROUPS = {
    {"SER:OG", single_hydrogen_group("CB", "O", 0.972, 108.539)},
    {"THR:OG1", single_hydrogen_group("CB", "O", 0.972, 109.544)},
    {"TYR:OH", single_hydrogen_group("CZ", "O", 0.966, 109.970)},
    {"CYS:SG", single_hydrogen_group("CB", "S", 1.338, 97.543)},

    {"ALA:CB", methyl_group("CA", 1.092, 109.742)},
    {"VAL:CG1", methyl_group("CB", 1.092, 109.5)},
    {"VAL:CG2", methyl_group("CB", 1.092, 109.5)},
    {"LEU:CD1", methyl_group("CG", 1.092, 109.5)},
    {"LEU:CD2", methyl_group("CG", 1.092, 109.5)},
    {"ILE:CD1", methyl_group("CG1", 1.092, 109.5)},
    {"ILE:CG2", methyl_group("CB", 1.092, 109.5)},
    {"MET:CE", methyl_group("SD", 1.092, 109.5)},
    {"MSE:CE", methyl_group("SE", 1.092, 109.5)},
    {"THR:CG2", methyl_group("CB", 1.092, 109.532)}
};

inline const RotatableGroupDefinition* get_rotatable_group_definition(
        const std::string& residue_name,
        const std::string& atom_name) {
    auto it = ROTATABLE_GROUPS.find(residue_name + ":" + atom_name);
    return it == ROTATABLE_GROUPS.end() ? nullptr : &it->second;
}

}  // namespace xhpi
