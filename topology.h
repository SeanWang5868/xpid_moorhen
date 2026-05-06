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

inline void add_ring_if_new(std::vector<AromaticRing>& rings,
                            std::unordered_set<std::string>& seen,
                            std::vector<std::string> atoms,
                            const std::string& source) {
    atoms.erase(std::remove_if(atoms.begin(), atoms.end(),
                               [](const std::string& atom) { return atom.empty(); }),
                atoms.end());
    if (atoms.size() != 5 && atoms.size() != 6) return;

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

inline std::vector<AromaticRing> read_monomer_aromatic_rings(const std::string& res_upper,
                                                             const std::filesystem::path& cif_path) {
    std::vector<AromaticRing> rings;
    std::unordered_set<std::string> seen;

    gemmi::cif::Document doc = gemmi::cif::read_file(cif_path.string());
    const gemmi::cif::Block* block = find_component_block(doc, res_upper);
    if (!block) return rings;

    gemmi::ChemComp cc = gemmi::make_chemcomp_from_block(*block);

    for (const gemmi::Restraints::Plane& plane : cc.rt.planes) {
        std::vector<std::string> atoms;
        atoms.reserve(plane.ids.size());
        for (const gemmi::Restraints::AtomId& atom_id : plane.ids) {
            atoms.push_back(gemmi::trim_str(atom_id.atom));
        }
        add_ring_if_new(rings, seen, std::move(atoms), "monomer-plane");
    }

    if (!rings.empty()) return rings;

    std::map<std::string, std::set<std::string>> graph;
    for (const gemmi::Restraints::Bond& bond : cc.rt.bonds) {
        if (!bond.aromatic && bond.type != gemmi::BondType::Aromatic) continue;
        std::string a1 = gemmi::trim_str(bond.id1.atom);
        std::string a2 = gemmi::trim_str(bond.id2.atom);
        if (a1.empty() || a2.empty()) continue;
        graph[a1].insert(a2);
        graph[a2].insert(a1);
    }

    if (graph.empty()) return rings;

    std::set<std::vector<std::string>> found_cycles;
    int visits = 0;
    for (const auto& item : graph) {
        std::vector<std::string> path = {item.first};
        dfs_aromatic_cycles(graph, item.first, path, found_cycles, visits);
        if (visits > 500) break;
    }

    for (const std::vector<std::string>& cycle : found_cycles) {
        add_ring_if_new(rings, seen, cycle, "monomer-aromatic-bond");
    }

    return rings;
}

inline const std::vector<AromaticRing>& get_aromatic_rings(const std::string& res_name,
                                                           const std::string& mon_lib_path = "") {
    static std::unordered_map<std::string, std::vector<AromaticRing>> cache;

    std::string res_upper = upper_copy(res_name);
    std::string lib_path = mon_lib_path.empty() ? default_monomer_library_path() : mon_lib_path;
    std::string cache_key = lib_path + "::" + res_upper;

    auto cached = cache.find(cache_key);
    if (cached != cache.end()) return cached->second;

    std::vector<AromaticRing> rings;
    std::filesystem::path cif_path = find_monomer_cif_path(res_upper, lib_path);
    if (!cif_path.empty()) {
        try {
            rings = read_monomer_aromatic_rings(res_upper, cif_path);
        } catch (...) {
            rings.clear();
        }
    }

    if (rings.empty()) {
        rings = rings_from_fallback(res_upper);
    }

    auto inserted = cache.emplace(cache_key, std::move(rings));
    return inserted.first->second;
}

} // namespace xhpi
