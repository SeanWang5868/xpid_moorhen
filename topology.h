#pragma once

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <gemmi/chemcomp.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/util.hpp>

#include "config.h"

namespace xhpi {

struct AromaticRing {
    std::vector<std::string> atoms;
    std::string source;
};

struct DonorTopology {
    std::vector<std::string> heavy_neighbors;
    std::vector<std::string> bonded_hydrogens;
};

struct MonomerTopology {
    std::vector<AromaticRing> rings;
    std::unordered_map<std::string, std::string> atom_elements;
    std::unordered_map<std::string, DonorTopology> donors;
};

inline std::string upper_copy(std::string s) {
    return gemmi::to_upper(gemmi::trim_str(s));
}

inline std::string default_monomer_library_path() {
    if (const char* env = std::getenv("GEMMI_MON_LIB_PATH")) {
        if (*env != '\0') return std::string(env);
    }
    if (const char* ccp4_lib = std::getenv("CCP4_LIB")) {
        if (*ccp4_lib != '\0') {
            std::filesystem::path monomers = std::filesystem::path(ccp4_lib) / "data" / "monomers";
            return monomers.string();
        }
    }
    return "";
}

inline std::filesystem::path find_monomer_cif_path(const std::string& res_upper,
                                                    const std::string& mon_lib_path) {
    if (res_upper.empty() || mon_lib_path.empty()) return {};

    std::filesystem::path root(mon_lib_path);
    std::string first_letter(1, static_cast<char>(std::tolower(static_cast<unsigned char>(res_upper[0]))));
    std::vector<std::filesystem::path> candidates = {
        root / first_letter / (res_upper + ".cif"),
        root / first_letter / (res_upper + ".cif.gz"),
        root / (res_upper + ".cif"),
        root / (res_upper + ".cif.gz")
    };

    for (const auto& path : candidates) {
        if (std::filesystem::exists(path)) return path;
    }
    return {};
}

inline const gemmi::cif::Block* find_component_block(const gemmi::cif::Document& doc,
                                                      const std::string& res_upper) {
    const std::string target = "COMP_" + res_upper;
    for (const auto& block : doc.blocks) {
        if (upper_copy(block.name) == target) return &block;
    }
    if (doc.blocks.size() == 1) return &doc.blocks.front();
    return nullptr;
}

inline std::string ring_key(std::vector<std::string> atoms) {
    std::sort(atoms.begin(), atoms.end());
    std::string key;
    for (const std::string& atom : atoms) {
        key += atom;
        key += '\x1f';
    }
    return key;
}

inline bool is_ring_element(const std::unordered_map<std::string, std::string>& atom_elements,
                            const std::string& atom_name) {
    auto it = atom_elements.find(atom_name);
    if (it == atom_elements.end()) return false;
    return it->second == "C" || it->second == "N" || it->second == "O" || it->second == "S";
}

inline void add_ring_if_new(std::vector<AromaticRing>& rings,
                            std::unordered_set<std::string>& seen,
                            std::vector<std::string> atoms,
                            const std::string& source,
                            const std::unordered_map<std::string, std::string>* atom_elements = nullptr) {
    atoms.erase(std::remove_if(atoms.begin(), atoms.end(),
                               [](const std::string& atom) { return atom.empty(); }),
                atoms.end());
    if (atoms.size() != 5 && atoms.size() != 6) return;
    if (atom_elements) {
        for (const std::string& atom : atoms) {
            if (!is_ring_element(*atom_elements, atom)) return;
        }
    }

    std::string key = ring_key(atoms);
    if (seen.insert(key).second) {
        rings.push_back({std::move(atoms), source});
    }
}

inline std::vector<AromaticRing> rings_from_fallback(const std::string& res_upper) {
    std::vector<AromaticRing> rings;
    auto it = FALLBACK_RINGS.find(res_upper);
    if (it == FALLBACK_RINGS.end()) return rings;
    for (const auto& atom_names : it->second) {
        rings.push_back({atom_names, "fallback"});
    }
    return rings;
}

inline void dfs_aromatic_cycles(const std::map<std::string, std::set<std::string>>& graph,
                                const std::string& start,
                                std::vector<std::string>& path,
                                std::set<std::vector<std::string>>& found,
                                int& visits) {
    if (++visits > 500) return;
    if (path.size() > 8) return;

    auto graph_it = graph.find(path.back());
    if (graph_it == graph.end()) return;

    for (const std::string& nb : graph_it->second) {
        if (nb == start && (path.size() == 5 || path.size() == 6)) {
            std::vector<std::string> cycle = path;
            std::sort(cycle.begin(), cycle.end());
            found.insert(std::move(cycle));
        } else if (std::find(path.begin(), path.end(), nb) == path.end()) {
            path.push_back(nb);
            dfs_aromatic_cycles(graph, start, path, found, visits);
            path.pop_back();
        }
    }
}

inline void find_cycles_in_graph(const std::map<std::string, std::set<std::string>>& graph,
                                 std::vector<std::vector<std::string>>& cycles) {
    std::set<std::vector<std::string>> found_cycles;
    int visits = 0;
    for (const auto& item : graph) {
        std::vector<std::string> path = {item.first};
        dfs_aromatic_cycles(graph, item.first, path, found_cycles, visits);
        if (visits > 500) break;
    }
    cycles.assign(found_cycles.begin(), found_cycles.end());
}

inline bool is_standard_polymer_residue(const std::string& res_upper) {
    static const std::unordered_set<std::string> residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "A", "C", "G", "U", "DA", "DC", "DG", "DT"
    };
    return residues.count(res_upper) != 0;
}

inline MonomerTopology build_monomer_topology(const std::string& res_upper,
                                              const gemmi::ChemComp& cc) {
    MonomerTopology topology;
    std::unordered_set<std::string> seen;

    for (const gemmi::ChemComp::Atom& atom : cc.atoms) {
        topology.atom_elements[gemmi::trim_str(atom.id)] = gemmi::to_upper(std::string(atom.el.name()));
    }

    std::map<std::string, std::set<std::string>> all_bond_graph;
    std::map<std::string, std::set<std::string>> aromatic_bond_graph;
    std::unordered_map<std::string, std::vector<std::string>> bonded_h;
    std::unordered_map<std::string, std::vector<std::string>> heavy_neighbors;

    for (const gemmi::Restraints::Bond& bond : cc.rt.bonds) {
        std::string a1 = gemmi::trim_str(bond.id1.atom);
        std::string a2 = gemmi::trim_str(bond.id2.atom);
        if (a1.empty() || a2.empty()) continue;

        std::string e1 = topology.atom_elements.count(a1) ? topology.atom_elements[a1] : "";
        std::string e2 = topology.atom_elements.count(a2) ? topology.atom_elements[a2] : "";
        bool h1 = e1 == "H" || e1 == "D";
        bool h2 = e2 == "H" || e2 == "D";

        if (!h1 && !h2) {
            heavy_neighbors[a1].push_back(a2);
            heavy_neighbors[a2].push_back(a1);
            if (is_ring_element(topology.atom_elements, a1) && is_ring_element(topology.atom_elements, a2)) {
                all_bond_graph[a1].insert(a2);
                all_bond_graph[a2].insert(a1);
                if (bond.aromatic || bond.type == gemmi::BondType::Aromatic) {
                    aromatic_bond_graph[a1].insert(a2);
                    aromatic_bond_graph[a2].insert(a1);
                }
            }
        } else if (h1 && !h2) {
            bonded_h[a2].push_back(a1);
        } else if (!h1 && h2) {
            bonded_h[a1].push_back(a2);
        }
    }

    for (const auto& item : bonded_h) {
        DonorTopology donor;
        donor.bonded_hydrogens = item.second;
        auto heavy_it = heavy_neighbors.find(item.first);
        if (heavy_it != heavy_neighbors.end()) donor.heavy_neighbors = heavy_it->second;
        topology.donors[item.first] = std::move(donor);
    }

    for (const gemmi::Restraints::Plane& plane : cc.rt.planes) {
        std::vector<std::string> atoms;
        atoms.reserve(plane.ids.size());
        for (const gemmi::Restraints::AtomId& atom_id : plane.ids) {
            atoms.push_back(gemmi::trim_str(atom_id.atom));
        }
        add_ring_if_new(topology.rings, seen, atoms, "monomer-plane", &topology.atom_elements);

        if (atoms.size() > 6) {
            std::unordered_set<std::string> plane_atoms(atoms.begin(), atoms.end());
            std::map<std::string, std::set<std::string>> plane_graph;
            for (const auto& item : all_bond_graph) {
                if (plane_atoms.count(item.first) == 0) continue;
                for (const std::string& nb : item.second) {
                    if (plane_atoms.count(nb) != 0) {
                        plane_graph[item.first].insert(nb);
                    }
                }
            }
            std::vector<std::vector<std::string>> plane_cycles;
            find_cycles_in_graph(plane_graph, plane_cycles);
            for (const auto& cycle : plane_cycles) {
                add_ring_if_new(topology.rings, seen, cycle, "monomer-plane-cycle", &topology.atom_elements);
            }
        }
    }

    if (!topology.rings.empty()) return topology;

    std::vector<std::vector<std::string>> aromatic_cycles;
    find_cycles_in_graph(aromatic_bond_graph, aromatic_cycles);
    for (const auto& cycle : aromatic_cycles) {
        add_ring_if_new(topology.rings, seen, cycle, "monomer-aromatic-bond", &topology.atom_elements);
    }

    if (!topology.rings.empty()) return topology;

    if (!is_standard_polymer_residue(res_upper)) {
        std::vector<std::vector<std::string>> broad_cycles;
        find_cycles_in_graph(all_bond_graph, broad_cycles);
        for (const auto& cycle : broad_cycles) {
            add_ring_if_new(topology.rings, seen, cycle, "monomer-bond-cycle", &topology.atom_elements);
        }
    }

    return topology;
}

inline MonomerTopology read_monomer_topology(const std::string& res_upper,
                                             const std::filesystem::path& cif_path) {
    gemmi::cif::Document doc = gemmi::cif::read_file(cif_path.string());
    const gemmi::cif::Block* block = find_component_block(doc, res_upper);
    if (!block) return {};
    gemmi::ChemComp cc = gemmi::make_chemcomp_from_block(*block);
    return build_monomer_topology(res_upper, cc);
}

inline const MonomerTopology& get_monomer_topology(const std::string& res_name,
                                                   const std::string& mon_lib_path = "") {
    static std::unordered_map<std::string, MonomerTopology> cache;

    std::string res_upper = upper_copy(res_name);
    std::string lib_path = mon_lib_path.empty() ? default_monomer_library_path() : mon_lib_path;
    std::string cache_key = lib_path + "::" + res_upper;

    MonomerTopology topology;
    std::filesystem::path cif_path = find_monomer_cif_path(res_upper, lib_path);
    auto cached = cache.find(cache_key);
    if (cached != cache.end()) {
        const bool has_topology = !cached->second.rings.empty() ||
                                  !cached->second.donors.empty() ||
                                  !cached->second.atom_elements.empty();
        if (has_topology || cif_path.empty()) return cached->second;
        cache.erase(cached);
    }

    if (!cif_path.empty()) {
        try {
            topology = read_monomer_topology(res_upper, cif_path);
        } catch (...) {
            topology = {};
        }
    }

    if (topology.rings.empty()) {
        topology.rings = rings_from_fallback(res_upper);
    }

    auto inserted = cache.emplace(cache_key, std::move(topology));
    return inserted.first->second;
}

inline const std::vector<AromaticRing>& get_aromatic_rings(const std::string& res_name,
                                                           const std::string& mon_lib_path = "") {
    return get_monomer_topology(res_name, mon_lib_path).rings;
}

inline std::string get_monomer_parent_atom(const std::string& res_name,
                                           const std::string& atom_name,
                                           const std::string& mon_lib_path = "") {
    const MonomerTopology& topology = get_monomer_topology(res_name, mon_lib_path);
    auto donor_it = topology.donors.find(atom_name);
    if (donor_it == topology.donors.end()) return "";
    const DonorTopology& donor = donor_it->second;
    if (donor.bonded_hydrogens.empty()) return "";
    if (donor.heavy_neighbors.size() != 1) return "";
    return donor.heavy_neighbors.front();
}

inline bool monomer_atom_has_dictionary_hydrogens(const std::string& res_name,
                                                  const std::string& atom_name,
                                                  const std::string& mon_lib_path = "") {
    const MonomerTopology& topology = get_monomer_topology(res_name, mon_lib_path);
    auto donor_it = topology.donors.find(atom_name);
    return donor_it != topology.donors.end() && !donor_it->second.bonded_hydrogens.empty();
}

} // namespace xhpi
