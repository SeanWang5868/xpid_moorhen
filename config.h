#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace xhpi {

inline constexpr double DIST_SEARCH_LIMIT = 6.0;
inline constexpr double MIN_COVALENT_XH = 0.5;
inline constexpr double METAL_BLOCKING_RADIUS = 2.6;
inline constexpr double PLANARITY_CUTOFF = 0.5;
inline constexpr double MIN_ATOM_OCCUPANCY = 0.10;

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
    {"PHE", {{"CD1", "CD2", "CE1", "CE2", "CZ", "CG"}}},
    {"HIS", {{"CE1", "ND1", "NE2", "CG", "CD2"}}}
};

inline const std::unordered_map<std::string, double> COVALENT_XH_MAX = {
    {"C", 1.25},
    {"N", 1.25},
    {"O", 1.20},
    {"S", 1.55}
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

inline double get_covalent_xh_max(const std::string& elem) {
    auto it = COVALENT_XH_MAX.find(elem);
    return it == COVALENT_XH_MAX.end() ? 0.0 : it->second;
}

inline bool is_cation_donor(const std::string& res_name, const std::string& atom_name) {
    return CATION_DONORS.count(res_name + ":" + atom_name) != 0;
}

} // namespace xhpi
