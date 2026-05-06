#pragma once

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "core.h"

namespace xhpi {

inline std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (unsigned char c : value) {
        switch (c) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20) {
                    out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c)
                        << std::dec << std::setfill(' ');
                } else {
                    out << static_cast<char>(c);
                }
        }
    }
    return out.str();
}

inline std::string json_string(const std::string& value) {
    return "\"" + json_escape(value) + "\"";
}

inline std::string json_number(double value, int precision) {
    if (!std::isfinite(value)) return "null";
    std::ostringstream out;
    out.setf(std::ios::fixed);
    out << std::setprecision(precision) << value;
    return out.str();
}

inline void write_json_field(std::ostringstream& json,
                             const std::string& key,
                             const std::string& value,
                             bool comma = true) {
    json << "    " << json_string(key) << ": " << json_string(value);
    if (comma) json << ",";
    json << "\n";
}

inline void write_json_field(std::ostringstream& json,
                             const std::string& key,
                             int value,
                             bool comma = true) {
    json << "    " << json_string(key) << ": " << value;
    if (comma) json << ",";
    json << "\n";
}

inline void write_json_field(std::ostringstream& json,
                             const std::string& key,
                             double value,
                             int precision,
                             bool comma = true) {
    json << "    " << json_string(key) << ": " << json_number(value, precision);
    if (comma) json << ",";
    json << "\n";
}

inline std::string format_results_as_json(const std::vector<InteractionResult>& results) {
    std::ostringstream json;
    json << "[\n";

    for (size_t i = 0; i < results.size(); ++i) {
        const InteractionResult& r = results[i];
        json << "  {\n";
        write_json_field(json, "pdb", r.pdb);
        write_json_field(json, "model", r.model);
        write_json_field(json, "resolution", r.resolution, 3);

        write_json_field(json, "pi_chain", r.pi_chain);
        write_json_field(json, "pi_res", r.pi_res);
        write_json_field(json, "pi_id", r.pi_id);
        write_json_field(json, "pi_seqid", r.pi_seqid);
        write_json_field(json, "pi_ring", r.pi_ring);

        write_json_field(json, "X_chain", r.X_chain);
        write_json_field(json, "X_res", r.X_res);
        write_json_field(json, "X_id", r.X_id);
        write_json_field(json, "X_seqid", r.X_seqid);
        write_json_field(json, "X_atom", r.X_atom);
        write_json_field(json, "H_atom", r.H_atom);

        write_json_field(json, "dist_X_Pi", r.dist_X_Pi, 3);
        write_json_field(json, "is_plevin", r.is_plevin);
        write_json_field(json, "is_hudson", r.is_hudson);
        write_json_field(json, "remark", r.remark);

        write_json_field(json, "pi_ss_type", r.pi_ss_type);
        write_json_field(json, "pi_ss_id", r.pi_ss_id);
        write_json_field(json, "X_ss_type", r.X_ss_type);
        write_json_field(json, "X_ss_id", r.X_ss_id);

        write_json_field(json, "pi_avg_b", r.pi_avg_b, 2);
        write_json_field(json, "pi_center_x", r.pi_center_x, 3);
        write_json_field(json, "pi_center_y", r.pi_center_y, 3);
        write_json_field(json, "pi_center_z", r.pi_center_z, 3);
        write_json_field(json, "X_b", r.X_b, 2);
        write_json_field(json, "X_xyz_x", r.X_xyz_x, 3);
        write_json_field(json, "X_xyz_y", r.X_xyz_y, 3);
        write_json_field(json, "X_xyz_z", r.X_xyz_z, 3);

        write_json_field(json, "seq_sep", r.seq_sep);
        write_json_field(json, "combined_occ", r.combined_occ, 3);
        write_json_field(json, "method", r.method);
        write_json_field(json, "theta", r.theta, 2);
        write_json_field(json, "angle_xh_pi", r.angle_xh_pi, 2);
        write_json_field(json, "angle_xpcn", r.angle_xpcn, 2);
        write_json_field(json, "angle_XH_Pi", r.angle_xh_pi, 2);
        write_json_field(json, "angle_XPCN", r.angle_xpcn, 2);
        write_json_field(json, "proj_dist", r.proj_dist, 3);
        write_json_field(json, "sym_op", r.sym_op, false);
        json << "  }";
        if (i + 1 < results.size()) json << ",";
        json << "\n";
    }

    json << "]\n";
    return json.str();
}

inline std::string detect_xhpi_interactions_json(const gemmi::Structure& st) {
    std::vector<InteractionResult> raw_results = detect_interactions(st);
    return format_results_as_json(raw_results);
}

inline std::string xpid(const gemmi::Structure& st) {
    return detect_xhpi_interactions_json(st);
}

} // namespace xhpi
