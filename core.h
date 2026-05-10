#pragma once

#include <algorithm>
#include <cmath>
#include <exception>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gemmi/model.hpp>
#include <gemmi/monlib.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/topo.hpp>
#include <gemmi/util.hpp>

#include "config.h"
#include "geometry.h"
#include "topology.h"

namespace xhpi {

inline std::string element_name_upper(const gemmi::Element& element) {
    return gemmi::to_upper(std::string(element.name()));
}

inline bool is_target_x_element(const std::string& elem) {
    return TARGET_ELEMENTS_X.count(elem) != 0;
}

inline bool is_hydrogen_element(const std::string& elem) {
    return TARGET_ELEMENTS_H.count(elem) != 0;
}

inline int seqid_int_or_zero(const gemmi::Residue& residue) {
    return residue.seqid.num.has_value() ? residue.seqid.num.value : 0;
}

inline std::string seqid_string(const gemmi::Residue& residue) {
    return residue.seqid.str();
}

inline std::string model_id_for(const gemmi::Model& model, size_t model_idx) {
    return model.num != 0 ? std::to_string(model.num) : std::to_string(model_idx + 1);
}

inline bool same_residue(const gemmi::const_CRA& a, const gemmi::const_CRA& b) {
    return a.chain == b.chain && a.residue == b.residue;
}

inline void normalize_atom_and_residue_names(gemmi::Structure& structure) {
    for (gemmi::Model& model : structure.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                residue.name = upper_copy(residue.name);
                for (gemmi::Atom& atom : residue.atoms) {
                    atom.name = gemmi::trim_str(atom.name);
                }
            }
        }
    }
}

inline void select_best_altconf(gemmi::Structure& structure) {
    for (gemmi::Model& model : structure.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                std::set<char> altlocs;
                for (const gemmi::Atom& atom : residue.atoms) {
                    if (atom.altloc != '\0') altlocs.insert(atom.altloc);
                }
                if (altlocs.empty()) continue;

                if (altlocs.size() == 1) {
                    for (gemmi::Atom& atom : residue.atoms) {
                        if (atom.altloc != '\0') atom.altloc = '\0';
                    }
                    continue;
                }

                std::map<char, double> occ_sum;
                std::map<char, int> occ_count;
                for (char alt : altlocs) {
                    occ_sum[alt] = 0.0;
                    occ_count[alt] = 0;
                }
                for (const gemmi::Atom& atom : residue.atoms) {
                    if (atom.altloc != '\0') {
                        occ_sum[atom.altloc] += atom.occ;
                        occ_count[atom.altloc] += 1;
                    }
                }

                char best_alt = *altlocs.begin();
                double best_avg = -1.0;
                for (char alt : altlocs) {
                    double avg = occ_count[alt] > 0 ? occ_sum[alt] / static_cast<double>(occ_count[alt]) : 0.0;
                    if (avg > best_avg || (avg == best_avg && alt < best_alt)) {
                        best_avg = avg;
                        best_alt = alt;
                    }
                }

                std::vector<gemmi::Atom> kept;
                kept.reserve(residue.atoms.size());
                for (gemmi::Atom atom : residue.atoms) {
                    if (atom.altloc != '\0' && atom.altloc != best_alt) continue;
                    if (atom.altloc == best_alt) atom.altloc = '\0';
                    kept.push_back(std::move(atom));
                }
                residue.atoms.swap(kept);
            }
        }
    }
}

inline void prepare_structure_for_detection(gemmi::Structure& structure) {
    normalize_atom_and_residue_names(structure);
    select_best_altconf(structure);
    gemmi::remove_waters(structure);
    structure.remove_empty_chains();
}

inline std::vector<std::string> collect_residue_names(const gemmi::Structure& structure) {
    std::set<std::string> seen;
    std::vector<std::string> resnames;
    for (const gemmi::Model& model : structure.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& residue : chain.residues) {
                std::string name = upper_copy(residue.name);
                if (!name.empty() && seen.insert(name).second) {
                    resnames.push_back(std::move(name));
                }
            }
        }
    }
    return resnames;
}

inline bool remove_residue_from_topology_error(gemmi::Structure& structure,
                                               const std::string& message) {
    static const std::regex bonded_pattern(R"(bonded to ([^/]+)/([^ ]+) ([^/]+)/([^ ]+) failed)");
    std::smatch match;
    if (!std::regex_search(message, match, bonded_pattern)) return false;

    std::string bad_chain = match[1].str();
    std::string bad_seqid = gemmi::trim_str(match[3].str());

    for (gemmi::Model& model : structure.models) {
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != bad_chain) continue;
            for (auto it = chain.residues.begin(); it != chain.residues.end(); ++it) {
                if (gemmi::trim_str(it->seqid.str()) == bad_seqid) {
                    chain.residues.erase(it);
                    return true;
                }
            }
        }
    }
    return false;
}

inline void add_hydrogens_from_monomer_library(gemmi::Structure& structure) {
    std::string mon_lib_path = default_monomer_library_path();
    if (mon_lib_path.empty() || structure.models.empty()) return;

    std::vector<std::string> resnames = collect_residue_names(structure);
    if (resnames.empty()) return;

    gemmi::MonLib monlib;
    gemmi::Logger logger;
    try {
        monlib.read_monomer_lib(mon_lib_path, resnames, logger);
    } catch (...) {
        return;
    }

    constexpr int max_attempts = 10;
    for (size_t model_idx = 0; model_idx < structure.models.size(); ++model_idx) {
        for (int attempt = 0; attempt < max_attempts; ++attempt) {
            try {
                gemmi::prepare_topology(structure, monlib, model_idx,
                                        gemmi::HydrogenChange::ReAddButWater,
                                        false, logger, true);
                break;
            } catch (const std::exception& err) {
                std::string message = err.what();
                if (message.find("link") != std::string::npos) {
                    structure.connections.clear();
                    continue;
                }
                if (remove_residue_from_topology_error(structure, message)) {
                    continue;
                }
                // Keep the original heavy-atom model usable if hydrogen placement
                // cannot be completed for a problematic component or link.
                break;
            } catch (...) {
                break;
            }
        }
    }
}

inline void prepare_structure_for_detection(gemmi::Structure& structure, bool add_hydrogens) {
    normalize_atom_and_residue_names(structure);
    select_best_altconf(structure);
    if (add_hydrogens) {
        add_hydrogens_from_monomer_library(structure);
    }
    gemmi::remove_waters(structure);
    structure.remove_empty_chains();
}

inline bool is_excluded_donor_atom(const std::string& res_name, const std::string& atom_name) {
    return ((res_name == "ASP" || res_name == "GLU") &&
            (atom_name == "OD1" || atom_name == "OD2" || atom_name == "OE1" || atom_name == "OE2")) ||
           atom_name == "OXT";
}

inline bool is_donor_blocked(const gemmi::Atom& x_atom,
                             gemmi::Model& model,
                             gemmi::NeighborSearch& ns,
                             const gemmi::Position& search_pos) {
    std::string x_elem = element_name_upper(x_atom.element);
    auto close_atoms = ns.find_atoms(search_pos, '\0', 0.0, METAL_BLOCKING_RADIUS);

    for (const gemmi::NeighborSearch::Mark* mark : close_atoms) {
        double dist = mark->pos.dist(search_pos);
        if (dist < 0.01) continue;

        gemmi::const_CRA cra = mark->to_cra(model);
        if (!cra.atom) continue;

        std::string neighbor_el = element_name_upper(cra.atom->element);
        if (x_elem == "S" && neighbor_el == "S" && cra.atom->name == "SG" &&
            dist >= 1.8 && dist <= 2.2) {
            return true;
        }

        if (BLOCKING_METALS.count(neighbor_el) != 0 && dist <= METAL_BLOCKING_RADIUS) {
            return true;
        }
    }
    return false;
}

inline std::vector<gemmi::Position> collect_positions(const std::vector<const gemmi::Atom*>& atoms) {
    std::vector<gemmi::Position> positions;
    positions.reserve(atoms.size());
    for (const gemmi::Atom* atom : atoms) positions.push_back(atom->pos);
    return positions;
}

inline std::string join_remark_parts(const std::vector<std::string>& parts) {
    std::string remark;
    for (size_t i = 0; i < parts.size(); ++i) {
        if (i != 0) remark += ", ";
        remark += parts[i];
    }
    return remark;
}

inline void maybe_add_pi_pi_remark(std::vector<std::string>& remarks,
                                   const gemmi::Residue& donor_residue,
                                   const gemmi::Position& pi_center,
                                   const gemmi::Position& pi_normal) {
    const auto& donor_rings = get_aromatic_rings(donor_residue.name);
    for (const AromaticRing& donor_ring : donor_rings) {
        std::vector<const gemmi::Atom*> donor_atoms;
        donor_atoms.reserve(donor_ring.atoms.size());
        for (const std::string& atom_name : donor_ring.atoms) {
            if (const gemmi::Atom* atom = donor_residue.find_atom(atom_name, '*')) {
                donor_atoms.push_back(atom);
            }
        }
        if (donor_atoms.size() != donor_ring.atoms.size()) continue;

        std::vector<gemmi::Position> positions = collect_positions(donor_atoms);
        gemmi::Position donor_center = calculate_center(positions);
        gemmi::Position donor_normal = calculate_normal(positions);
        auto [pp_dist, pp_angle, pp_offset] = calculate_pi_pi_geometry(
            pi_center, pi_normal, donor_center, donor_normal);
        (void) pp_offset;

        if (pp_dist < 3.0 || pp_dist > PI_PI_DIST_MAX) continue;
        std::ostringstream ss;
        ss.setf(std::ios::fixed);
        ss.precision(1);
        if (pp_angle <= PI_PI_ANGLE_PARALLEL_MAX) {
            ss << "Pi-Pi Parallel d=" << pp_dist;
            remarks.push_back(ss.str());
            return;
        }
        if (pp_angle >= PI_PI_ANGLE_TSHAPED_MIN) {
            ss << "Pi-Pi T-shaped d=" << pp_dist;
            remarks.push_back(ss.str());
            return;
        }
    }
}

inline void record_hit(std::vector<InteractionResult>& results,
                       const gemmi::Structure& structure,
                       const std::string& model_id,
                       const gemmi::Chain& pi_chain,
                       const gemmi::Residue& pi_res,
                       const AromaticRing& ring,
                       int ring_index,
                       const gemmi::Position& pi_center,
                       const gemmi::Position& pi_normal,
                       double pi_b_mean,
                       double avg_pi_occ,
                       const gemmi::const_CRA& x_cra,
                       const gemmi::Atom& x_atom,
                       const std::string& h_name,
                       double dist_x_pi,
                       int is_plevin,
                       int is_hudson,
                       const gemmi::Position& x_pos,
                       double theta,
                       double xh_pi_angle,
                       double xpcn_angle,
                       double proj_dist,
                       bool is_cone,
                       double combined_occ,
                       int sym_op,
                       double min_occ) {
    if (combined_occ < min_occ) return;

    std::vector<std::string> remarks;
    if (is_cone) remarks.push_back("Cone");
    if (pi_res.name == "TRP" && ring.atoms.size() == 5) remarks.push_back("TRP 5-ring acceptor");
    if (sym_op != 0) remarks.push_back("SymContact op " + std::to_string(sym_op));
    if (is_cation_donor(x_cra.residue->name, x_atom.name)) remarks.push_back("Cation-pi");
    maybe_add_pi_pi_remark(remarks, *x_cra.residue, pi_center, pi_normal);

    InteractionResult result;
    result.pdb = structure.name.empty() ? "UNKNOWN" : structure.name;
    result.model = model_id;
    result.resolution = structure.resolution;

    result.pi_chain = pi_chain.name;
    result.pi_res = pi_res.name;
    result.pi_id = seqid_int_or_zero(pi_res);
    result.pi_seqid = seqid_string(pi_res);
    result.pi_ring = "ring" + std::to_string(ring_index + 1) + ":" + ring.source;

    result.X_chain = x_cra.chain->name;
    result.X_res = x_cra.residue->name;
    result.X_id = seqid_int_or_zero(*x_cra.residue);
    result.X_seqid = seqid_string(*x_cra.residue);
    result.X_atom = x_atom.name;
    result.H_atom = h_name;
    result.method = is_cone ? "Implicit/Cone" : "Explicit";

    result.dist_X_Pi = dist_x_pi;
    result.is_plevin = is_plevin;
    result.is_hudson = is_hudson;
    result.remark = join_remark_parts(remarks);

    result.pi_avg_b = pi_b_mean;
    result.pi_center_x = pi_center.x;
    result.pi_center_y = pi_center.y;
    result.pi_center_z = pi_center.z;
    result.X_b = x_atom.b_iso;
    result.X_xyz_x = x_pos.x;
    result.X_xyz_y = x_pos.y;
    result.X_xyz_z = x_pos.z;

    if (pi_chain.name == x_cra.chain->name &&
        pi_res.seqid.num.has_value() &&
        x_cra.residue->seqid.num.has_value()) {
        result.seq_sep = pi_res.seqid.num.value - x_cra.residue->seqid.num.value;
    }
    result.combined_occ = combined_occ;
    result.theta = theta;
    result.angle_xh_pi = xh_pi_angle;
    result.angle_xpcn = xpcn_angle;
    result.proj_dist = proj_dist;
    result.sym_op = sym_op;

    (void) avg_pi_occ;
    results.push_back(std::move(result));
}

inline void collect_cone_environment(gemmi::Model& model,
                                     gemmi::NeighborSearch& ns,
                                     const gemmi::const_CRA& x_cra,
                                     const gemmi::Position& x_pos,
                                     std::vector<gemmi::Position>& env_coords,
                                     std::vector<gemmi::Position>& acceptor_coords) {
    auto neighbors = ns.find_atoms(x_pos, '\0', 0.0, 4.0);
    for (const gemmi::NeighborSearch::Mark* mark : neighbors) {
        if (mark->pos.dist(x_pos) < 0.01) continue;
        gemmi::const_CRA n_cra = mark->to_cra(model);
        if (!n_cra.atom || same_residue(n_cra, x_cra)) continue;

        std::string elem = element_name_upper(n_cra.atom->element);
        if (elem == "H" || elem == "D" || elem.empty()) continue;

        double dist = mark->pos.dist(x_pos);
        if (dist <= 4.0) env_coords.push_back(mark->pos);
        if (dist <= 3.5 && (elem == "O" || elem == "N" || elem == "S")) {
            acceptor_coords.push_back(mark->pos);
        }
    }
}

inline std::vector<gemmi::Position> generate_wobbled_hydrogens(const gemmi::Position& parent_pos,
                                                               const gemmi::Position& x_pos,
                                                               const std::vector<gemmi::Position>& orig_h_positions,
                                                               const std::vector<gemmi::Position>& env_coords) {
    std::vector<gemmi::Position> candidates;
    gemmi::Position axis = normalize(x_pos - parent_pos);
    if (magnitude(axis) < 1e-8) return candidates;

    const int angles[] = {-20, -15, -10, -5, 5, 10, 15, 20};
    for (const gemmi::Position& h_orig : orig_h_positions) {
        gemmi::Position vec_xh = h_orig - x_pos;
        for (int angle_deg : angles) {
            double angle_rad = static_cast<double>(angle_deg) * M_PI / 180.0;
            gemmi::Position h_pos = x_pos + rotate_vector_around_axis(vec_xh, axis, angle_rad);
            if (!has_clash(h_pos, env_coords, 2.0)) candidates.push_back(h_pos);
        }
    }
    return candidates;
}

inline bool run_explicit_track(std::vector<InteractionResult>& results,
                               std::vector<gemmi::Position>& orig_h_positions,
                               const gemmi::Structure& structure,
                               gemmi::Model& model,
                               gemmi::NeighborSearch& ns,
                               const std::string& model_id,
                               const gemmi::Chain& pi_chain,
                               const gemmi::Residue& pi_res,
                               const AromaticRing& ring,
                               int ring_index,
                               const gemmi::Position& pi_center,
                               const gemmi::Position& pi_normal,
                               double pi_b_mean,
                               double avg_pi_occ,
                               const gemmi::const_CRA& x_cra,
                               const gemmi::Atom& x_atom,
                               const gemmi::Position& x_pos,
                               double dist_x_pi,
                               double max_dist,
                               double xpcn_angle,
                               double proj_dist,
                               double combined_occ,
                               int sym_op,
                               double min_occ) {
    bool found = false;
    auto h_marks = ns.find_atoms(x_pos, x_atom.altloc, 0.0, DIST_CUTOFF_H);

    for (const gemmi::NeighborSearch::Mark* h_mark : h_marks) {
        if (h_mark->image_idx != 0) continue;

        gemmi::const_CRA h_cra = h_mark->to_cra(model);
        if (!h_cra.atom || !same_residue(h_cra, x_cra)) continue;

        const gemmi::Atom& h_atom = *h_cra.atom;
        std::string h_elem = element_name_upper(h_atom.element);
        if (!is_hydrogen_element(h_elem)) continue;

        double xh_dist = x_pos.dist(h_atom.pos);
        if (xh_dist <= MIN_COVALENT_XH || xh_dist > DIST_CUTOFF_H) continue;

        orig_h_positions.push_back(h_atom.pos);

        if (h_atom.altloc != '\0' && x_atom.altloc != '\0' && h_atom.altloc != x_atom.altloc) continue;

        double xh_pi_angle = calculate_xh_picenter_angle(x_pos, h_atom.pos, pi_center);
        double theta = calculate_hudson_theta(pi_center, x_pos, h_atom.pos, pi_normal);
        if (xh_pi_angle < 0.0 || theta < 0.0 || xpcn_angle < 0.0) continue;

        int is_plevin = (dist_x_pi < max_dist && xh_pi_angle >= 120.0 && xpcn_angle < 25.0) ? 1 : 0;
        int is_hudson = (dist_x_pi <= max_dist && std::isfinite(proj_dist) &&
                         proj_dist <= (ring.atoms.size() == 6 ? 2.0 : 1.6) && theta <= 40.0) ? 1 : 0;

        if (is_plevin || is_hudson) {
            found = true;
            double h_combined_occ = std::min(combined_occ, static_cast<double>(h_atom.occ));
            record_hit(results, structure, model_id, pi_chain, pi_res, ring, ring_index,
                       pi_center, pi_normal, pi_b_mean, avg_pi_occ, x_cra, x_atom,
                       h_atom.name, dist_x_pi, is_plevin, is_hudson, x_pos, theta,
                       xh_pi_angle, xpcn_angle, proj_dist, false, h_combined_occ, sym_op, min_occ);
        }
    }

    return found;
}

inline void run_cone_track(std::vector<InteractionResult>& results,
                           const std::vector<gemmi::Position>& orig_h_positions,
                           const gemmi::Structure& structure,
                           gemmi::Model& model,
                           gemmi::NeighborSearch& ns,
                           const std::string& model_id,
                           const gemmi::Chain& pi_chain,
                           const gemmi::Residue& pi_res,
                           const AromaticRing& ring,
                           int ring_index,
                           const gemmi::Position& pi_center,
                           const gemmi::Position& pi_normal,
                           double pi_b_mean,
                           double avg_pi_occ,
                           const gemmi::const_CRA& x_cra,
                           const gemmi::Atom& x_atom,
                           const gemmi::Position& x_pos,
                           double dist_x_pi,
                           double max_dist,
                           double xpcn_angle,
                           double proj_dist,
                           double combined_occ,
                           int sym_op,
                           double min_occ,
                           bool generate_missing_h_cone) {
    (void) generate_missing_h_cone;
    std::string parent_name = get_cone_parent_atom(x_cra.residue->name, x_atom.name);
    if (parent_name.empty()) return;

    const gemmi::Atom* parent_atom = x_cra.residue->find_atom(parent_name, '*');
    if (!parent_atom) return;

    std::vector<gemmi::Position> env_coords;
    std::vector<gemmi::Position> acceptor_coords;
    collect_cone_environment(model, ns, x_cra, x_pos, env_coords, acceptor_coords);

    bool locked = check_hbond_locked(x_pos, orig_h_positions, acceptor_coords);
    if (locked) return;

    std::string x_elem = element_name_upper(x_atom.element);
    std::vector<gemmi::Position> h_candidates;
    if (is_flexible_donor(x_cra.residue->name, x_atom.name)) {
        h_candidates = generate_rotated_hydrogens(parent_atom->pos, x_pos, x_elem, env_coords, 2.0, 72);
    } else if (!orig_h_positions.empty()) {
        h_candidates = generate_wobbled_hydrogens(parent_atom->pos, x_pos, orig_h_positions, env_coords);
    }

    double best_xh_angle = -1.0;
    double best_theta = 0.0;
    double best_angle = 180.0;
    int best_plevin = 0;
    int best_hudson = 0;

    for (const gemmi::Position& h_pos : h_candidates) {
        double theta = calculate_hudson_theta(pi_center, x_pos, h_pos, pi_normal);
        double xh_pi_angle = calculate_xh_picenter_angle(x_pos, h_pos, pi_center);
        if (theta < 0.0 || xh_pi_angle < 0.0 || xpcn_angle < 0.0) continue;

        int is_plevin = (dist_x_pi < max_dist && xpcn_angle < 25.0 && xh_pi_angle >= 120.0) ? 1 : 0;
        int is_hudson = (dist_x_pi <= max_dist && std::isfinite(proj_dist) &&
                         proj_dist <= (ring.atoms.size() == 6 ? 2.0 : 1.6) && theta <= 40.0) ? 1 : 0;

        if ((is_plevin || is_hudson) && xh_pi_angle > best_xh_angle) {
            best_xh_angle = xh_pi_angle;
            best_theta = theta;
            best_angle = xh_pi_angle;
            best_plevin = is_plevin;
            best_hudson = is_hudson;
        }
    }

    if (best_xh_angle >= 0.0) {
        record_hit(results, structure, model_id, pi_chain, pi_res, ring, ring_index,
                   pi_center, pi_normal, pi_b_mean, avg_pi_occ, x_cra, x_atom, "(virt)",
                   dist_x_pi, best_plevin, best_hudson, x_pos, best_theta, best_angle,
                   xpcn_angle, proj_dist, true, combined_occ, sym_op, min_occ);
    }
}

inline std::vector<InteractionResult> detect_interactions(const gemmi::Structure& input_structure,
                                                          bool use_cone = true,
                                                          double min_occ = 0.0,
                                                          bool generate_missing_h_cone = true,
                                                          bool add_hydrogens = true) {
    std::vector<InteractionResult> results;
    if (input_structure.models.empty()) return results;

    gemmi::Structure structure = input_structure;
    prepare_structure_for_detection(structure, add_hydrogens);

    for (size_t model_idx = 0; model_idx < structure.models.size(); ++model_idx) {
        gemmi::Model& model = structure.models[model_idx];
        std::string model_id = model_id_for(model, model_idx);

        gemmi::NeighborSearch ns(model, structure.cell, DIST_SEARCH_LIMIT);
        ns.populate(true);

        for (const gemmi::Chain& pi_chain : model.chains) {
            for (const gemmi::Residue& pi_res : pi_chain.residues) {
                const auto& rings = get_aromatic_rings(pi_res.name);
                if (rings.empty()) continue;

                for (size_t ring_idx = 0; ring_idx < rings.size(); ++ring_idx) {
                    const AromaticRing& ring = rings[ring_idx];

                    std::vector<const gemmi::Atom*> pi_atoms;
                    pi_atoms.reserve(ring.atoms.size());
                    for (const std::string& atom_name : ring.atoms) {
                        if (const gemmi::Atom* atom = pi_res.find_atom(atom_name, '*')) {
                            pi_atoms.push_back(atom);
                        }
                    }
                    if (pi_atoms.size() != ring.atoms.size()) continue;

                    std::vector<gemmi::Position> ring_positions = collect_positions(pi_atoms);
                    gemmi::Position pi_center = calculate_center(ring_positions);
                    gemmi::Position pi_normal = calculate_normal(ring_positions);
                    double planarity = calculate_planarity_deviation(ring_positions, pi_center, pi_normal);
                    if (planarity > PLANARITY_CUTOFF) continue;

                    double avg_pi_occ = 0.0;
                    for (const gemmi::Atom* atom : pi_atoms) avg_pi_occ += atom->occ;
                    avg_pi_occ /= static_cast<double>(pi_atoms.size());
                    if (avg_pi_occ < MIN_ATOM_OCCUPANCY) continue;

                    double pi_b_mean = calculate_mean_b(pi_atoms);
                    char pi_alt = pi_atoms.empty() ? '\0' : pi_atoms.front()->altloc;
                    auto candidates = ns.find_atoms(pi_center, pi_alt, 0.0, DIST_SEARCH_LIMIT);

                    for (const gemmi::NeighborSearch::Mark* x_mark : candidates) {
                        if (x_mark->image_idx != 0) continue;

                        gemmi::const_CRA x_cra = x_mark->to_cra(model);
                        if (!x_cra.atom || !x_cra.residue || !x_cra.chain) continue;
                        if (x_cra.chain == &pi_chain && x_cra.residue == &pi_res) continue;

                        const gemmi::Atom& x_atom = *x_cra.atom;
                        std::string x_elem = element_name_upper(x_atom.element);
                        if (!is_target_x_element(x_elem)) continue;
                        if (is_excluded_donor_atom(x_cra.residue->name, x_atom.name)) continue;
                        if (x_atom.occ < MIN_ATOM_OCCUPANCY) continue;
                        if (pi_alt != '\0' && x_atom.altloc != '\0' && pi_alt != x_atom.altloc) continue;
                        if (is_donor_blocked(x_atom, model, ns, x_atom.pos)) continue;

                        gemmi::Position x_pos = x_atom.pos;
                        double dist_x_pi = calculate_distance(x_pos, pi_center);
                        double max_dist = get_dynamic_threshold(x_elem);
                        if (dist_x_pi > max_dist) continue;

                        double xpcn_angle = calculate_xpcn_angle(pi_normal, x_pos, pi_center);
                        double proj_dist = calculate_projection_dist(pi_normal, pi_center, x_pos);
                        if (xpcn_angle < 0.0 || !std::isfinite(proj_dist)) continue;

                        double combined_occ = std::min(avg_pi_occ, static_cast<double>(x_atom.occ));
                        int sym_op = 0;

                        std::vector<gemmi::Position> orig_h_positions;
                        bool found_explicit = run_explicit_track(
                            results, orig_h_positions, structure, model, ns, model_id,
                            pi_chain, pi_res, ring, static_cast<int>(ring_idx), pi_center, pi_normal,
                            pi_b_mean, avg_pi_occ, x_cra, x_atom, x_pos, dist_x_pi, max_dist,
                            xpcn_angle, proj_dist, combined_occ, sym_op, min_occ);

                        if (!found_explicit && use_cone &&
                            !is_cone_scan_suppressed(x_cra.residue->name, x_atom.name)) {
                            run_cone_track(
                                results, orig_h_positions, structure, model, ns, model_id,
                                pi_chain, pi_res, ring, static_cast<int>(ring_idx), pi_center, pi_normal,
                                pi_b_mean, avg_pi_occ, x_cra, x_atom, x_pos, dist_x_pi, max_dist,
                                xpcn_angle, proj_dist, combined_occ, sym_op, min_occ,
                                generate_missing_h_cone);
                        }
                    }
                }
            }
        }
    }

    return results;
}

} // namespace xhpi
