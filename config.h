#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace xhpi {

inline constexpr double DIST_SEARCH_LIMIT = 6.0;
inline constexpr double DIST_CUTOFF_H = 1.3;
inline constexpr double MIN_COVALENT_XH = 0.5;
inline constexpr double METAL_BLOCKING_RADIUS = 2.6;
inline constexpr double PLANARITY_CUTOFF = 0.5;
inline constexpr double MIN_ATOM_OCCUPANCY = 0.10;
inline constexpr double TETRAHEDRAL_ANGLE = 109.5;

inline const std::unordered_set<std::string> BLOCKING_METALS = {
    "ZN", "FE", "CU", "MN", "MG", "CO", "NI", "CA", "CD", "HG",
    "NA", "K", "PT", "AU", "AG", "FE2", "FE3"
};

inline const std::unordered_set<std::string> TARGET_ELEMENTS_X = {"C", "N", "O", "S"};
inline const std::unordered_set<std::string> TARGET_ELEMENTS_H = {"H", "D"};

inline const std::unordered_map<std::string, std::vector<std::vector<std::string>>> FALLBACK_RINGS = {
    {"TRP", {
        {"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
        {"CD1", "CD2", "NE1", "CG", "CE2"}
    }},
    {"TYR", {{"CD1", "CD2", "CE1", "CE2", "CZ", "CG"}}},
    {"PTR", {{"CD1", "CD2", "CE1", "CE2", "CZ", "CG"}}},
    {"PHE", {{"CD1", "CD2", "CE1", "CE2", "CZ", "CG"}}},
    {"HIS", {{"CE1", "ND1", "NE2", "CG", "CD2"}}},
    {"BER", {
        {"C1", "N1", "C3", "C6", "C8", "C12"},
        {"C8", "C12", "C13", "C15", "C16", "C18"},
        {"C2", "C4", "C5", "C9", "C11", "C14"}
    }},
    {"4PO", {{"N2", "C6", "C7", "C8", "C9", "C10"}}}
};

inline const std::unordered_map<std::string, std::unordered_map<std::string, std::string>> ROTATABLE_MAPPING = {
    {"ALA", {{"CB", "CA"}}},
    {"VAL", {{"CG1", "CB"}, {"CG2", "CB"}}},
    {"LEU", {{"CD1", "CG"}, {"CD2", "CG"}}},
    {"ILE", {{"CD1", "CG1"}, {"CG2", "CB"}}},
    {"MET", {{"CE", "SD"}}},
    {"MSE", {{"CE", "SE"}}},
    {"THR", {{"CG2", "CB"}, {"OG1", "CB"}}},
    {"SER", {{"OG", "CB"}}},
    {"TYR", {{"OH", "CZ"}}},
    {"CYS", {{"SG", "CB"}}},
    {"LYS", {{"NZ", "CE"}}}
};

inline const std::unordered_set<std::string> FLEXIBLE_DONORS = {"OG", "OG1", "OH", "SG"};
inline const std::unordered_set<std::string> RIGID_DONORS = {"CB", "CG1", "CG2", "CD1", "CD2", "CE", "NZ"};

inline const std::unordered_map<std::string, double> BOND_LENGTHS = {
    {"C", 1.09},
    {"N", 1.01},
    {"O", 0.96},
    {"S", 1.33}
};

inline const std::unordered_set<std::string> CATION_DONORS = {
    "LYS:NZ",
    "ARG:NH1",
    "ARG:NH2",
    "ARG:NE"
};

inline constexpr double PI_PI_DIST_MAX = 5.5;
inline constexpr double PI_PI_ANGLE_PARALLEL_MAX = 35.0;
inline constexpr double PI_PI_ANGLE_TSHAPED_MIN = 60.0;

inline double get_dynamic_threshold(const std::string& elem) {
    if (elem == "N" || elem == "O") return 4.3;
    if (elem == "C") return 4.5;
    if (elem == "S") return 4.8;
    return 4.5;
}

inline double get_bond_length(const std::string& elem) {
    auto it = BOND_LENGTHS.find(elem);
    return it == BOND_LENGTHS.end() ? 1.09 : it->second;
}

inline std::string get_cone_parent_atom(const std::string& res_name, const std::string& atom_name) {
    auto res_it = ROTATABLE_MAPPING.find(res_name);
    if (res_it == ROTATABLE_MAPPING.end()) return "";
    auto atom_it = res_it->second.find(atom_name);
    return atom_it == res_it->second.end() ? "" : atom_it->second;
}

inline bool is_flexible_donor(const std::string& res_name, const std::string& atom_name) {
    return (res_name == "SER" && atom_name == "OG") ||
           (res_name == "THR" && atom_name == "OG1") ||
           (res_name == "TYR" && atom_name == "OH") ||
           (res_name == "CYS" && atom_name == "SG");
}

inline bool is_cation_donor(const std::string& res_name, const std::string& atom_name) {
    return CATION_DONORS.count(res_name + ":" + atom_name) != 0;
}

inline bool is_cone_scan_suppressed(const std::string& res_name, const std::string& atom_name) {
    return (res_name == "ARG" && (atom_name == "NH1" || atom_name == "NH2" || atom_name == "NE")) ||
           ((res_name == "ASN" || res_name == "GLN") && (atom_name == "ND2" || atom_name == "NE2")) ||
           (res_name == "HIS" && (atom_name == "ND1" || atom_name == "NE2")) ||
           (res_name == "TRP" && atom_name == "NE1");
}

} // namespace xhpi
