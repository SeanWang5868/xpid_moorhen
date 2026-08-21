#pragma once

#include <algorithm>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include <gemmi/model.hpp>

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

struct DonorConformer {
    const gemmi::Atom* x_atom = nullptr;
    char x_altloc = '\0';
    const gemmi::Atom* parent_atom = nullptr;
    char parent_altloc = '\0';

    char active_altloc() const {
        return x_altloc != '\0' ? x_altloc : parent_altloc;
    }

    double occupancy() const {
        if (!x_atom || !parent_atom) return 0.0;
        return std::min(static_cast<double>(x_atom->occ),
                        static_cast<double>(parent_atom->occ));
    }
};

struct DonorConformerResolution {
    std::vector<DonorConformer> conformers;
    std::string issue;
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

inline DonorConformerResolution resolve_donor_conformers(
        const gemmi::Residue& residue,
        const gemmi::Atom& x_atom,
        const RotatableGroupDefinition& definition) {
    std::map<char, std::vector<const gemmi::Atom*>> parents_by_altloc;
    for (const gemmi::Atom& atom : residue.atoms) {
        if (atom.name == definition.parent_atom_name) {
            parents_by_altloc[atom.altloc].push_back(&atom);
        }
    }
    if (parents_by_altloc.empty()) return {{}, "missing_parent"};

    for (const auto& item : parents_by_altloc) {
        if (item.second.size() != 1) {
            const std::string label = item.first == '\0'
                ? "<blank>" : std::string(1, item.first);
            return {{}, "duplicate_parent_altloc:" + label};
        }
    }

    DonorConformerResolution resolution;
    const char x_altloc = x_atom.altloc;
    if (x_altloc != '\0') {
        auto exact = parents_by_altloc.find(x_altloc);
        if (exact != parents_by_altloc.end()) {
            resolution.conformers.push_back(
                {&x_atom, x_altloc, exact->second.front(), x_altloc});
            return resolution;
        }
        auto blank = parents_by_altloc.find('\0');
        if (blank != parents_by_altloc.end()) {
            resolution.conformers.push_back(
                {&x_atom, x_altloc, blank->second.front(), '\0'});
            return resolution;
        }
        std::string available;
        for (const auto& item : parents_by_altloc) {
            if (!available.empty()) available += ',';
            available += item.first == '\0' ? "<blank>" : std::string(1, item.first);
        }
        resolution.issue = "incompatible_parent_altloc:" + available;
        return resolution;
    }

    auto blank = parents_by_altloc.find('\0');
    if (blank != parents_by_altloc.end()) {
        resolution.conformers.push_back(
            {&x_atom, '\0', blank->second.front(), '\0'});
        return resolution;
    }
    for (const auto& item : parents_by_altloc) {
        resolution.conformers.push_back(
            {&x_atom, '\0', item.second.front(), item.first});
    }
    return resolution;
}

}  // namespace xhpi
