#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include <gemmi/model.hpp>

#include "config.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace xhpi {

struct InteractionResult {
    std::string pdb = "UNKNOWN";
    std::string model = "1";
    double resolution = 0.0;

    std::string pi_chain;
    std::string pi_res;
    int pi_id = 0;
    std::string pi_seqid;
    std::string pi_ring;

    std::string X_chain;
    std::string X_res;
    int X_id = 0;
    std::string X_seqid;
    std::string X_atom;
    std::string H_atom;
    std::string method;

    double dist_X_Pi = 0.0;
    int is_plevin = 0;
    int is_hudson = 0;
    std::string remark;

    std::string pi_ss_type;
    std::string pi_ss_id;
    std::string X_ss_type;
    std::string X_ss_id;

    double pi_avg_b = 0.0;
    double pi_center_x = 0.0;
    double pi_center_y = 0.0;
    double pi_center_z = 0.0;
    double X_b = 0.0;
    double X_xyz_x = 0.0;
    double X_xyz_y = 0.0;
    double X_xyz_z = 0.0;

    int seq_sep = 0;
    double combined_occ = 1.0;
    double theta = 0.0;
    double angle_xh_pi = 180.0;
    double angle_xpcn = 0.0;
    double proj_dist = 0.0;
    int sym_op = 0;
};

inline double dot_product(const gemmi::Position& v1, const gemmi::Position& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline gemmi::Position cross_product(const gemmi::Position& v1, const gemmi::Position& v2) {
    return gemmi::Position(
        v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x
    );
}

inline double magnitude(const gemmi::Position& v) {
    return std::sqrt(dot_product(v, v));
}

inline gemmi::Position normalize(const gemmi::Position& v) {
    double mag = magnitude(v);
    if (mag < 1e-8) return gemmi::Position(0, 0, 0);
    return v / mag;
}

inline double calculate_angle_vectors(const gemmi::Position& v1, const gemmi::Position& v2) {
    double mag1 = magnitude(v1);
    double mag2 = magnitude(v2);
    if (mag1 < 1e-8 || mag2 < 1e-8) return -1.0;
    double cos_val = dot_product(v1, v2) / (mag1 * mag2);
    cos_val = std::max(-1.0, std::min(1.0, cos_val));
    return std::acos(cos_val) * (180.0 / M_PI);
}

inline double calculate_distance(const gemmi::Position& p1, const gemmi::Position& p2) {
    return p1.dist(p2);
}

inline gemmi::Position calculate_center(const std::vector<gemmi::Position>& positions) {
    if (positions.empty()) return gemmi::Position(0, 0, 0);
    gemmi::Position sum(0, 0, 0);
    for (const auto& pos : positions) sum += pos;
    return sum / static_cast<double>(positions.size());
}

inline double calculate_mean_b(const std::vector<const gemmi::Atom*>& atoms) {
    if (atoms.empty()) return 0.0;
    double total = 0.0;
    for (const gemmi::Atom* atom : atoms) total += atom->b_iso;
    return total / static_cast<double>(atoms.size());
}

// Jacobi covariance eigenvector; equivalent to the SVD normal used in xpid2.
inline gemmi::Position calculate_normal(const std::vector<gemmi::Position>& ring_atoms) {
    size_t n = ring_atoms.size();
    if (n < 3) return gemmi::Position(0, 0, 0);
    gemmi::Position center = calculate_center(ring_atoms);

    double cov[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
    for (const auto& pos : ring_atoms) {
        double dx = pos.x - center.x;
        double dy = pos.y - center.y;
        double dz = pos.z - center.z;
        cov[0][0] += dx * dx; cov[0][1] += dx * dy; cov[0][2] += dx * dz;
        cov[1][1] += dy * dy; cov[1][2] += dy * dz;
        cov[2][2] += dz * dz;
    }
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    double v[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
    double d[3] = {cov[0][0], cov[1][1], cov[2][2]};
    double a[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            a[i][j] = cov[i][j];

    for (int step = 0; step < 50; ++step) {
        double sm = 0.0;
        for (int i = 0; i < 2; ++i)
            for (int j = i + 1; j < 3; ++j)
                sm += std::abs(a[i][j]);
        if (sm < 1e-12) break;

        for (int i = 0; i < 2; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                double g = 100.0 * std::abs(a[i][j]);
                if (step > 3 &&
                    std::abs(d[i]) + g == std::abs(d[i]) &&
                    std::abs(d[j]) + g == std::abs(d[j])) {
                    a[i][j] = 0.0;
                } else if (std::abs(a[i][j]) > 0.0) {
                    double h = d[j] - d[i];
                    double t;
                    if (std::abs(h) + g == std::abs(h)) {
                        t = a[i][j] / h;
                    } else {
                        double theta = 0.5 * h / a[i][j];
                        t = 1.0 / (std::abs(theta) + std::sqrt(1.0 + theta * theta));
                        if (theta < 0.0) t = -t;
                    }
                    double c = 1.0 / std::sqrt(1.0 + t * t);
                    double s = t * c;
                    double tau = s / (1.0 + c);
                    h = t * a[i][j];
                    d[i] -= h;
                    d[j] += h;
                    a[i][j] = 0.0;

                    for (int k = 0; k <= i - 1; ++k) {
                        double g2 = a[k][i], h2 = a[k][j];
                        a[k][i] = g2 - s * (h2 + g2 * tau);
                        a[k][j] = h2 + s * (g2 - h2 * tau);
                    }
                    for (int k = i + 1; k <= j - 1; ++k) {
                        double g2 = a[i][k], h2 = a[k][j];
                        a[i][k] = g2 - s * (h2 + g2 * tau);
                        a[k][j] = h2 + s * (g2 - h2 * tau);
                    }
                    for (int k = j + 1; k < 3; ++k) {
                        double g2 = a[i][k], h2 = a[j][k];
                        a[i][k] = g2 - s * (h2 + g2 * tau);
                        a[j][k] = h2 + s * (g2 - h2 * tau);
                    }
                    for (int k = 0; k < 3; ++k) {
                        double g2 = v[k][i], h2 = v[k][j];
                        v[k][i] = g2 - s * (h2 + g2 * tau);
                        v[k][j] = h2 + s * (g2 - h2 * tau);
                    }
                }
            }
        }
    }

    int min_idx = 0;
    if (d[1] < d[min_idx]) min_idx = 1;
    if (d[2] < d[min_idx]) min_idx = 2;
    return normalize(gemmi::Position(v[0][min_idx], v[1][min_idx], v[2][min_idx]));
}

inline double calculate_planarity_deviation(const std::vector<gemmi::Position>& ring_atoms,
                                            const gemmi::Position& center,
                                            const gemmi::Position& normal) {
    if (ring_atoms.size() < 3 || magnitude(normal) < 1e-8) return 999.0;
    gemmi::Position n = normalize(normal);
    double max_dev = 0.0;
    for (const auto& pos : ring_atoms) {
        double dev = std::abs(dot_product(pos - center, n));
        max_dev = std::max(max_dev, dev);
    }
    return max_dev;
}

inline double calculate_hudson_theta(const gemmi::Position& pi_center,
                                     const gemmi::Position& x_pos,
                                     const gemmi::Position& h_pos,
                                     const gemmi::Position& normal) {
    gemmi::Position v_x_pi = pi_center - x_pos;
    gemmi::Position v_xh = h_pos - x_pos;
    double norm_xpi = magnitude(v_x_pi);
    if (norm_xpi < 1e-8) return -1.0;

    double projection_length = dot_product(v_xh, v_x_pi) / norm_xpi;
    if (projection_length <= 0.0) return -1.0;

    double angle = calculate_angle_vectors(normal, v_xh);
    if (angle < 0.0) return -1.0;
    if (angle >= 90.0) angle = 180.0 - angle;
    return angle;
}

inline double calculate_xpcn_angle(const gemmi::Position& normal,
                                   const gemmi::Position& x_pos,
                                   const gemmi::Position& pi_center) {
    gemmi::Position v_x_pi = pi_center - x_pos;
    double angle = calculate_angle_vectors(normal, v_x_pi);
    if (angle < 0.0) return -1.0;
    if (angle > 90.0) angle = 180.0 - angle;
    return angle;
}

inline double calculate_xh_picenter_angle(const gemmi::Position& x_pos,
                                          const gemmi::Position& h_pos,
                                          const gemmi::Position& pi_center) {
    return calculate_angle_vectors(x_pos - h_pos, pi_center - h_pos);
}

inline double calculate_projection_dist(const gemmi::Position& normal,
                                        const gemmi::Position& pi_center,
                                        const gemmi::Position& x_pos) {
    double denominator = dot_product(normal, normal);
    if (denominator < 1e-12) return std::numeric_limits<double>::quiet_NaN();
    double t = dot_product(normal, pi_center - x_pos) / denominator;
    gemmi::Position projection_point = x_pos + normal * t;
    return projection_point.dist(pi_center);
}

inline bool check_h_direction(const gemmi::Position& x_pos,
                              const gemmi::Position& h_pos,
                              const gemmi::Position& pi_center) {
    return dot_product(h_pos - x_pos, pi_center - x_pos) > 0.0;
}

inline bool has_clash(const gemmi::Position& h_pos,
                      const std::vector<gemmi::Position>& env_coords,
                      double clash_cutoff) {
    for (const gemmi::Position& env : env_coords) {
        if (h_pos.dist(env) < clash_cutoff) return true;
    }
    return false;
}

inline bool check_hbond_locked(const gemmi::Position& x_pos,
                               const std::vector<gemmi::Position>& orig_h_positions,
                               const std::vector<gemmi::Position>& acceptor_coords,
                               double dist_cutoff = 3.5,
                               double angle_cutoff_deg = 120.0) {
    if (orig_h_positions.empty() || acceptor_coords.empty()) return false;

    for (const gemmi::Position& h_pos : orig_h_positions) {
        gemmi::Position v_xh = h_pos - x_pos;
        gemmi::Position xh_norm = normalize(v_xh);
        if (magnitude(xh_norm) < 1e-8) continue;

        for (const gemmi::Position& acc : acceptor_coords) {
            if (h_pos.dist(acc) > dist_cutoff) continue;
            gemmi::Position ha_norm = normalize(acc - h_pos);
            if (magnitude(ha_norm) < 1e-8) continue;
            double angle = calculate_angle_vectors(-xh_norm, ha_norm);
            if (angle >= angle_cutoff_deg) return true;
        }
    }
    return false;
}

inline gemmi::Position rotate_vector_around_axis(const gemmi::Position& vec,
                                                 const gemmi::Position& axis,
                                                 double angle_rad) {
    gemmi::Position u = normalize(axis);
    double cos_t = std::cos(angle_rad);
    double sin_t = std::sin(angle_rad);
    return vec * cos_t + cross_product(u, vec) * sin_t + u * dot_product(u, vec) * (1.0 - cos_t);
}

inline std::vector<gemmi::Position> generate_rotated_hydrogens(const gemmi::Position& parent_pos,
                                                               const gemmi::Position& x_pos,
                                                               const std::string& element,
                                                               const std::vector<gemmi::Position>& env_coords = {},
                                                               double clash_cutoff = 2.0,
                                                               int num_samples = 72) {
    std::vector<gemmi::Position> out;
    gemmi::Position axis = x_pos - parent_pos;
    gemmi::Position u = normalize(axis);
    if (magnitude(u) < 1e-8) return out;

    gemmi::Position arbitrary(1.0, 0.0, 0.0);
    if (std::abs(dot_product(u, arbitrary)) > 0.99) arbitrary = gemmi::Position(0.0, 1.0, 0.0);

    gemmi::Position v = normalize(cross_product(u, arbitrary));
    if (magnitude(v) < 1e-8) return out;
    gemmi::Position w = cross_product(u, v);

    double bond_length = get_bond_length(element);
    double theta_rad = TETRAHEDRAL_ANGLE * M_PI / 180.0;
    double h_proj_u = bond_length * std::cos(M_PI - theta_rad);
    double h_radius = bond_length * std::sin(M_PI - theta_rad);

    out.reserve(static_cast<size_t>(num_samples));
    for (int i = 0; i < num_samples; ++i) {
        double phi = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(num_samples);
        gemmi::Position h_pos = x_pos + u * h_proj_u + v * (h_radius * std::cos(phi)) +
                                w * (h_radius * std::sin(phi));
        if (!has_clash(h_pos, env_coords, clash_cutoff)) out.push_back(h_pos);
    }
    return out;
}

inline std::tuple<double, double, double> calculate_pi_pi_geometry(const gemmi::Position& center1,
                                                                   const gemmi::Position& normal1,
                                                                   const gemmi::Position& center2,
                                                                   const gemmi::Position& normal2) {
    gemmi::Position vec = center2 - center1;
    double dist = magnitude(vec);
    gemmi::Position n1 = normalize(normal1);
    gemmi::Position n2 = normalize(normal2);
    if (magnitude(n1) < 1e-8 || magnitude(n2) < 1e-8) return {dist, 90.0, dist};

    double cos_angle = std::abs(dot_product(n1, n2));
    cos_angle = std::max(0.0, std::min(1.0, cos_angle));
    double angle = std::acos(cos_angle) * 180.0 / M_PI;
    double proj_along = dist > 0.0 ? std::abs(dot_product(vec, n1)) : 0.0;
    double offset = std::sqrt(std::max(dist * dist - proj_along * proj_along, 0.0));
    return {dist, angle, offset};
}

} // namespace xhpi
