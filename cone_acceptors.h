#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>

#include <gemmi/model.hpp>
#include <gemmi/util.hpp>

namespace xhpi {

inline bool cone_acceptor_has_local_hydrogen(const gemmi::Residue& residue,
                                             const gemmi::Atom& atom) {
    for (const gemmi::Atom& candidate : residue.atoms) {
        const std::string element = gemmi::to_upper(std::string(candidate.element.name()));
        if (element != "H" && element != "D") continue;
        if (!gemmi::is_same_conformer(atom.altloc, candidate.altloc)) continue;
        if (atom.pos.dist(candidate.pos) <= 1.45) return true;
    }
    return false;
}

inline bool is_cone_hbond_acceptor(const gemmi::Residue& residue,
                                   const gemmi::Atom& atom) {
    static const std::unordered_set<std::string> backbone_acceptors = {"O", "OXT"};
    static const std::unordered_map<std::string, std::unordered_set<std::string>> sidechain_acceptors = {
        {"ASP", {"OD1", "OD2"}},
        {"GLU", {"OE1", "OE2"}},
        {"ASN", {"OD1"}},
        {"GLN", {"OE1"}},
        {"SER", {"OG"}},
        {"THR", {"OG1"}},
        {"TYR", {"OH"}},
        {"MET", {"SD"}},
        {"MSE", {"SE"}}
    };
    static const std::unordered_set<std::string> hydroxyl_acceptors = {
        "SER:OG", "THR:OG1", "TYR:OH"
    };
    static const std::unordered_set<std::string> non_acceptor_nitrogens = {
        "LYS:NZ", "ARG:NE", "ARG:NH1", "ARG:NH2", "ASN:ND2", "GLN:NE2"
    };

    const std::string name = gemmi::to_upper(gemmi::trim_str(atom.name));
    const std::string res_name = gemmi::to_upper(gemmi::trim_str(residue.name));
    const std::string element = gemmi::to_upper(std::string(atom.element.name()));
    const std::string key = res_name + ":" + name;

    if (element == "O" && backbone_acceptors.count(name) != 0) return true;
    if (hydroxyl_acceptors.count(key) != 0) return true;

    auto sidechain_it = sidechain_acceptors.find(res_name);
    if (sidechain_it != sidechain_acceptors.end() &&
        sidechain_it->second.count(name) != 0) {
        return !cone_acceptor_has_local_hydrogen(residue, atom);
    }
    if (non_acceptor_nitrogens.count(key) != 0) return false;
    if (res_name == "HIS" && (name == "ND1" || name == "NE2")) {
        return !cone_acceptor_has_local_hydrogen(residue, atom);
    }
    if (residue.is_water() && element == "O") return true;
    if (element == "O" || element == "S") {
        return !cone_acceptor_has_local_hydrogen(residue, atom);
    }
    return false;
}

}  // namespace xhpi
