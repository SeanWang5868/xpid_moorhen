#include <cmath>
#include <functional>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "core.h"

namespace {

int failures = 0;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void run_test(const std::string& name, const std::function<void()>& test) {
    try {
        test();
        std::cout << "[PASS] " << name << "\n";
    } catch (const std::exception& error) {
        ++failures;
        std::cerr << "[FAIL] " << name << ": " << error.what() << "\n";
    }
}

gemmi::Atom atom(const std::string& name,
                 const std::string& element,
                 double x,
                 double y,
                 double z,
                 double occupancy = 1.0,
                 char altloc = '\0') {
    gemmi::Atom value;
    value.name = name;
    value.element = gemmi::Element(element);
    value.pos = gemmi::Position(x, y, z);
    value.occ = static_cast<float>(occupancy);
    value.altloc = altloc;
    return value;
}

gemmi::Residue residue(const std::string& name, int seqid) {
    gemmi::Residue value;
    value.name = name;
    value.seqid = gemmi::SeqId(seqid, ' ');
    return value;
}

gemmi::Residue phe_ring(int seqid = 1) {
    gemmi::Residue phe = residue("PHE", seqid);
    phe.atoms = {
        atom("CG",  "C", 11.4, 10.0, 10.0),
        atom("CD1", "C", 10.7, 11.212, 10.0),
        atom("CE1", "C",  9.3, 11.212, 10.0),
        atom("CZ",  "C",  8.6, 10.0, 10.0),
        atom("CE2", "C",  9.3,  8.788, 10.0),
        atom("CD2", "C", 10.7,  8.788, 10.0),
    };
    return phe;
}

gemmi::Structure structure_with_donor(gemmi::Residue donor) {
    gemmi::Structure structure;
    structure.name = "contract";
    structure.cell = gemmi::UnitCell(30, 30, 30, 90, 90, 90);
    structure.models.emplace_back(1);
    structure.models.back().chains.emplace_back("A");
    gemmi::Chain& chain = structure.models.back().chains.back();
    chain.residues.push_back(phe_ring());
    chain.residues.push_back(std::move(donor));
    return structure;
}

std::vector<xhpi::InteractionResult> explicit_hits(
        gemmi::Structure structure,
        const std::string& monomer_library = "") {
    return xhpi::detect_interactions(
        structure,
        false,  // no legacy cone path in explicit-H contract tests
        0.0,
        false,
        false,
        monomer_library);
}

void test_fused_graph_returns_only_minimal_chordless_cycles() {
    std::vector<std::pair<std::string, std::string>> edges = {
        {"A1", "A2"}, {"A2", "A3"}, {"A3", "A4"},
        {"A4", "A5"}, {"A5", "A6"}, {"A6", "A1"},
        {"A3", "B1"}, {"B1", "B2"}, {"B2", "B3"},
        {"B3", "B4"}, {"B4", "A4"},
    };
    std::map<std::string, std::set<std::string>> graph;
    for (const auto& edge : edges) {
        graph[edge.first].insert(edge.second);
        graph[edge.second].insert(edge.first);
    }
    std::vector<std::vector<std::string>> cycles;
    xhpi::find_cycles_in_graph(graph, cycles);

    require(cycles.size() == 2, "fused graph must contain exactly two minimal rings");
    std::set<std::string> keys;
    for (const auto& cycle : cycles) keys.insert(xhpi::ring_key(cycle));
    require(keys.count(xhpi::ring_key({"A1", "A2", "A3", "A4", "A5", "A6"})) == 1,
            "first six-membered ring missing");
    require(keys.count(xhpi::ring_key({"A3", "A4", "B1", "B2", "B3", "B4"})) == 1,
            "second six-membered ring missing");

    std::map<std::string, std::set<std::string>> chorded;
    const std::vector<std::pair<std::string, std::string>> chorded_edges = {
        {"C1", "C2"}, {"C2", "C3"}, {"C3", "C4"},
        {"C4", "C5"}, {"C5", "C6"}, {"C6", "C1"},
        {"C1", "C4"},
    };
    for (const auto& edge : chorded_edges) {
        chorded[edge.first].insert(edge.second);
        chorded[edge.second].insert(edge.first);
    }
    cycles.clear();
    xhpi::find_cycles_in_graph(chorded, cycles);
    require(cycles.empty(), "a six-membered perimeter with an internal chord is not one aromatic ring");
}

void test_nonstandard_components_have_no_hardcoded_ring_fallback() {
    require(xhpi::rings_from_fallback("PTR").empty(), "PTR requires a usable dictionary");
    require(xhpi::rings_from_fallback("BER").empty(), "BER requires a usable dictionary");
    require(xhpi::rings_from_fallback("4PO").empty(), "4PO requires a usable dictionary");

    std::vector<xhpi::AromaticRing> incomplete_trp = {
        {{"CD1", "CD2", "NE1", "CG", "CE2"}, "dictionary"},
    };
    require(!xhpi::standard_rings_are_complete("TRP", incomplete_trp),
            "a partial TRP dictionary must trigger the complete two-ring fallback");
}

void test_ring_altloc_variants_never_mix_labels() {
    gemmi::Residue phe = residue("PHE", 1);
    const std::vector<std::string> names = {"CG", "CD1", "CE1", "CZ", "CE2", "CD2"};
    for (size_t i = 0; i < names.size(); ++i) {
        phe.atoms.push_back(atom(names[i], "C", static_cast<double>(i), 0, 0, 0.7, 'A'));
        phe.atoms.push_back(atom(names[i], "C", static_cast<double>(i), 1, 0, 0.3, 'B'));
    }
    auto variants = xhpi::atom_variants_for_names(phe, names);
    require(variants.size() == 2, "complete A and B rings must both be retained");
    require(variants[0].altloc == 'A' && variants[1].altloc == 'B',
            "ring variants must be deterministic by altloc");
    for (const auto& variant : variants) {
        for (const gemmi::Atom* ring_atom : variant.atoms) {
            require(ring_atom->altloc == variant.altloc, "ring variant mixed incompatible altlocs");
        }
    }

    gemmi::Residue incomplete = residue("PHE", 2);
    incomplete.atoms.push_back(atom("CG", "C", 0, 0, 0, 1.0, 'A'));
    incomplete.atoms.push_back(atom("CD1", "C", 1, 0, 0, 1.0, 'B'));
    for (const std::string& name : {"CE1", "CZ", "CE2", "CD2"}) {
        incomplete.atoms.push_back(atom(name, "C", 0, 0, 0));
    }
    require(xhpi::atom_variants_for_names(incomplete, names).empty(),
            "incomplete A/B atoms must not be combined into a ring");
}

void test_lower_occupancy_donor_altloc_can_contribute() {
    gemmi::Residue ser = residue("SER", 2);
    ser.atoms = {
        atom("OG", "O", 10, 10, 15.2, 0.8, 'A'),
        atom("HG", "H", 10, 10, 14.2, 0.8, 'A'),
        atom("OG", "O", 10, 10, 13.0, 0.2, 'B'),
        atom("HG", "H", 10, 10, 12.0, 0.2, 'B'),
    };
    auto hits = explicit_hits(structure_with_donor(std::move(ser)));
    require(hits.size() == 1, "lower-occupancy B conformer should produce one hit");
    require(hits.front().X_atom == "OG", "unexpected donor atom");
    require(std::abs(hits.front().combined_occ - 0.2) < 1e-6,
            "hit occupancy must retain the contributing altloc occupancy");
}

void test_incompatible_hydrogen_altloc_is_rejected() {
    gemmi::Residue ser = residue("SER", 2);
    ser.atoms = {
        atom("OG", "O", 10, 10, 13.0, 1.0, 'B'),
        atom("HG", "H", 10, 10, 12.0, 1.0, 'A'),
    };
    require(explicit_hits(structure_with_donor(std::move(ser))).empty(),
            "A hydrogen cannot be assigned to a B donor");
}

void test_cation_pi_atoms_are_not_xh_pi_donors() {
    gemmi::Residue lys = residue("LYS", 2);
    lys.atoms = {
        atom("NZ", "N", 10, 10, 13.0),
        atom("HZ1", "H", 10, 10, 12.0),
    };
    require(explicit_hits(structure_with_donor(std::move(lys))).empty(),
            "LYS NZ belongs to cation-pi, not XH-pi");
}

void test_unknown_component_rejects_ambiguous_hydrogen_owner() {
    gemmi::Residue ligand = residue("UNK", 2);
    ligand.atoms = {
        atom("O1", "O", 10.0, 10, 13.0),
        atom("N1", "N", 10.2, 10, 13.0),
        atom("H1", "H", 10.0, 10, 12.0),
    };
    require(explicit_hits(structure_with_donor(std::move(ligand))).empty(),
            "a nearby H with ambiguous heavy-atom ownership must be rejected");
}

void test_dictionary_bonds_are_authoritative_for_explicit_hydrogens() {
    const std::string monomer_library =
        (std::filesystem::path(XPID_TEST_SOURCE_DIR) / "test" / "fixtures" / "monomers").string();
    require(xhpi::dictionary_hydrogen_match(
                "BMA", "O6", "DO6", monomer_library) ==
                xhpi::DictionaryHydrogenMatch::Match,
            "dictionary H name must match its neutron D equivalent");
    require(xhpi::dictionary_hydrogen_match(
                "BMA", "O6", "H1", monomer_library) ==
                xhpi::DictionaryHydrogenMatch::Mismatch,
            "an unbonded nearby H must not be assigned to O6");

    gemmi::Residue bma = residue("BMA", 2);
    bma.atoms = {
        atom("O6", "O", 10, 10, 13.0),
        atom("H1", "H", 10, 10, 12.0),
    };
    require(explicit_hits(
                structure_with_donor(std::move(bma)), monomer_library).empty(),
            "dictionary mismatch must reject an otherwise plausible X-H geometry");
}

void test_sulfur_deuterium_uses_element_specific_cutoff() {
    gemmi::Residue cys = residue("CYS", 2);
    cys.atoms = {
        atom("SG", "S", 10, 10, 13.5),
        atom("DG", "D", 10, 10, 12.162),
    };
    auto hits = explicit_hits(structure_with_donor(std::move(cys)));
    require(hits.size() == 1, "1.338 A S-D bond must survive explicit-H detection");
    require(hits.front().H_atom == "DG", "neutron deuterium identity must be retained");
    require(xhpi::canonical_hydrogen_name("DG") == "HG",
            "deuterium name must canonicalize to dictionary hydrogen name");
}

}  // namespace

int main() {
    run_test("minimal chordless aromatic cycles", test_fused_graph_returns_only_minimal_chordless_cycles);
    run_test("nonstandard fallback exclusion", test_nonstandard_components_have_no_hardcoded_ring_fallback);
    run_test("ring altloc compatibility", test_ring_altloc_variants_never_mix_labels);
    run_test("lower occupancy donor altloc", test_lower_occupancy_donor_altloc_can_contribute);
    run_test("incompatible hydrogen altloc", test_incompatible_hydrogen_altloc_is_rejected);
    run_test("cation-pi exclusion", test_cation_pi_atoms_are_not_xh_pi_donors);
    run_test("ambiguous hydrogen ownership", test_unknown_component_rejects_ambiguous_hydrogen_owner);
    run_test("dictionary hydrogen ownership", test_dictionary_bonds_are_authoritative_for_explicit_hydrogens);
    run_test("sulfur-deuterium cutoff", test_sulfur_deuterium_uses_element_specific_cutoff);

    if (failures != 0) {
        std::cerr << failures << " contract test(s) failed\n";
        return 1;
    }
    std::cout << "All XPID Moorhen contract tests passed\n";
    return 0;
}
