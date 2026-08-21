#include <cmath>
#include <functional>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "core.h"
#include "cone.h"

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

void test_group_specific_single_hydrogen_geometry() {
    const gemmi::Position parent(0.0, 0.0, 0.0);
    const gemmi::Position x_pos(1.42, 0.0, 0.0);
    for (const auto& expected : std::vector<std::tuple<std::string, std::string, double, double>>{
             {"SER", "OG", 0.972, 108.539},
             {"CYS", "SG", 1.338, 97.543}}) {
        const auto* definition = xhpi::get_rotatable_group_definition(
            std::get<0>(expected), std::get<1>(expected));
        require(definition != nullptr, "rotatable single-H definition missing");
        auto conformers = xhpi::generate_group_conformers(parent, x_pos, *definition);
        require(conformers.size() == 360, "single-H group requires a full 360 degree scan");
        require(conformers.front().hydrogen_positions.size() == 1,
                "single-H group generated the wrong hydrogen count");
        const gemmi::Position& h_pos = conformers.front().hydrogen_positions.front();
        require(std::abs(h_pos.dist(x_pos) - std::get<2>(expected)) < 1e-10,
                "group-specific X-H nuclear distance changed");
        const double angle = xhpi::calculate_angle_vectors(parent - x_pos, h_pos - x_pos);
        require(std::abs(angle - std::get<3>(expected)) < 1e-10,
                "group-specific parent-X-H angle changed");
    }
}

void test_methyl_conformer_contains_three_coupled_hydrogens() {
    const auto* definition = xhpi::get_rotatable_group_definition("ALA", "CB");
    require(definition != nullptr, "ALA methyl definition missing");
    const gemmi::Position parent(0.0, 0.0, 0.0);
    const gemmi::Position x_pos(1.53, 0.0, 0.0);
    auto conformers = xhpi::generate_group_conformers(parent, x_pos, *definition);
    require(conformers.size() == 120, "methyl symmetry requires one 120 degree period");
    require(conformers.front().hydrogen_positions.size() == 3,
            "a methyl conformer must contain all three hydrogens");

    const auto& hydrogens = conformers.front().hydrogen_positions;
    for (const gemmi::Position& h_pos : hydrogens) {
        require(std::abs(h_pos.dist(x_pos) - 1.092) < 1e-10,
                "methyl C-H nuclear distance changed");
        const double angle = xhpi::calculate_angle_vectors(parent - x_pos, h_pos - x_pos);
        require(std::abs(angle - 109.742) < 1e-10,
                "methyl parent-C-H angle changed");
    }
    require(std::abs(hydrogens[0].dist(hydrogens[1]) -
                     hydrogens[1].dist(hydrogens[2])) < 1e-10,
            "methyl hydrogens are not rotationally coupled");
    require(xhpi::get_rotatable_group_definition("LYS", "NZ") == nullptr,
            "cationic LYS NZ must not enter the XH-pi cone model");
}

void test_any_hydrogen_clash_invalidates_complete_methyl_conformer() {
    const auto* definition = xhpi::get_rotatable_group_definition("ALA", "CB");
    const gemmi::Position parent(0.0, 0.0, 0.0);
    const gemmi::Position x_pos(1.53, 0.0, 0.0);
    auto conformer = xhpi::generate_group_conformers(parent, x_pos, *definition).front();

    gemmi::Residue ligand = residue("LIG", 9);
    ligand.atoms.push_back(atom("C1", "C", 0.0, 0.0, 0.0));
    const gemmi::Position blocker_pos = conformer.hydrogen_positions.front();
    std::vector<xhpi::ConeEnvironmentAtom> environment = {
        xhpi::make_cone_environment_atom(blocker_pos, ligand.atoms.front(), ligand)
    };
    require(!xhpi::cone_conformer_is_sterically_valid(conformer, x_pos, environment),
            "one clashing methyl H must invalidate the complete CH3 conformer");
}

void test_valid_hbond_contact_is_not_a_steric_clash() {
    const auto* definition = xhpi::get_rotatable_group_definition("SER", "OG");
    const gemmi::Position parent(0.0, 0.0, 0.0);
    const gemmi::Position x_pos(1.42, 0.0, 0.0);
    auto conformer = xhpi::generate_group_conformers(parent, x_pos, *definition).front();
    const gemmi::Position h_pos = conformer.hydrogen_positions.front();
    const gemmi::Position direction = xhpi::normalize(h_pos - x_pos);
    const gemmi::Position contact_pos = h_pos + direction * 1.8;

    gemmi::Residue asp = residue("ASP", 9);
    asp.atoms.push_back(atom("OD1", "O", contact_pos.x, contact_pos.y, contact_pos.z));
    std::vector<xhpi::ConeEnvironmentAtom> acceptor_environment = {
        xhpi::make_cone_environment_atom(contact_pos, asp.atoms.front(), asp)
    };
    require(xhpi::cone_conformer_is_sterically_valid(
                conformer, x_pos, acceptor_environment),
            "a valid 1.8 A, 180 degree H-bond contact must not count as a clash");

    gemmi::Residue lys = residue("LYS", 10);
    lys.atoms.push_back(atom("NZ", "N", contact_pos.x, contact_pos.y, contact_pos.z));
    std::vector<xhpi::ConeEnvironmentAtom> non_acceptor_environment = {
        xhpi::make_cone_environment_atom(contact_pos, lys.atoms.front(), lys)
    };
    require(!xhpi::cone_conformer_is_sterically_valid(
                conformer, x_pos, non_acceptor_environment),
            "the same close contact to non-acceptor LYS NZ must remain a clash");
}

void test_minimal_cone_acceptor_typing() {
    gemmi::Residue ala = residue("ALA", 1);
    ala.atoms.push_back(atom("O", "O", 0, 0, 0));
    ala.atoms.push_back(atom("N", "N", 1, 0, 0));
    require(xhpi::is_cone_hbond_acceptor(ala, ala.atoms[0]),
            "backbone carbonyl oxygen must be a Cone acceptor");
    require(!xhpi::is_cone_hbond_acceptor(ala, ala.atoms[1]),
            "backbone nitrogen must not be a Cone acceptor");

    gemmi::Residue ser = residue("SER", 2);
    ser.atoms.push_back(atom("OG", "O", 0, 0, 0));
    ser.atoms.push_back(atom("HG", "H", 0.96, 0, 0));
    require(xhpi::is_cone_hbond_acceptor(ser, ser.atoms[0]),
            "neutral hydroxyl oxygen remains an acceptor when protonated");
}

void test_cone_parent_altloc_resolution() {
    const auto* definition = xhpi::get_rotatable_group_definition("LEU", "CD1");
    require(definition != nullptr, "LEU CD1 rotatable definition missing");

    gemmi::Residue labelled = residue("LEU", 2);
    labelled.atoms = {
        atom("CG", "C", 0, 0, 0, 0.65, 'A'),
        atom("CG", "C", 0, 1, 0, 0.35, 'B'),
        atom("CD1", "C", 1, 1, 0, 0.35, 'B'),
    };
    auto exact = xhpi::resolve_donor_conformers(labelled, labelled.atoms[2], *definition);
    require(exact.issue.empty() && exact.conformers.size() == 1,
            "labelled donor must resolve one compatible parent");
    require(exact.conformers.front().parent_altloc == 'B',
            "labelled donor must prefer the exact parent altloc");
    require(std::abs(exact.conformers.front().occupancy() - 0.35) < 1e-6,
            "donor conformer occupancy must include its parent");

    gemmi::Residue shared_parent = residue("LEU", 3);
    shared_parent.atoms = {
        atom("CG", "C", 0, 0, 0),
        atom("CD1", "C", 1, 0, 0, 0.4, 'B'),
    };
    auto shared = xhpi::resolve_donor_conformers(
        shared_parent, shared_parent.atoms[1], *definition);
    require(shared.issue.empty() && shared.conformers.size() == 1 &&
            shared.conformers.front().parent_altloc == '\0',
            "labelled donor may use one shared blank parent");

    gemmi::Residue blank_donor = residue("LEU", 4);
    blank_donor.atoms = {
        atom("CG", "C", 0, 0, 0, 0.7, 'A'),
        atom("CG", "C", 0, 1, 0, 0.3, 'B'),
        atom("CD1", "C", 1, 0, 0),
    };
    auto split = xhpi::resolve_donor_conformers(
        blank_donor, blank_donor.atoms[2], *definition);
    require(split.issue.empty() && split.conformers.size() == 2,
            "blank donor with labelled parents must produce both conformers");
    require(split.conformers[0].parent_altloc == 'A' &&
            split.conformers[1].parent_altloc == 'B',
            "split parent conformers must be deterministic by altloc");

    gemmi::Residue incompatible = residue("LEU", 5);
    incompatible.atoms = {
        atom("CG", "C", 0, 0, 0, 1.0, 'A'),
        atom("CD1", "C", 1, 0, 0, 1.0, 'B'),
    };
    auto rejected = xhpi::resolve_donor_conformers(
        incompatible, incompatible.atoms[1], *definition);
    require(rejected.conformers.empty() && !rejected.issue.empty(),
            "incompatible labelled parent must be rejected explicitly");
}

void test_binary_cone_uses_deterministic_positive_evidence() {
    const auto* definition = xhpi::get_rotatable_group_definition("SER", "OG");
    require(definition != nullptr, "SER OG rotatable definition missing");
    const gemmi::Position parent(11.42, 10.0, 13.0);
    const gemmi::Position x_pos(10.0, 10.0, 13.0);
    const gemmi::Position pi_center(10.0, 10.0, 10.0);
    const gemmi::Position pi_normal(0.0, 0.0, 1.0);
    auto evidence = xhpi::evaluate_binary_cone(
        parent, x_pos, *definition, {}, pi_center, pi_normal,
        3.0, 4.3, 0.0, 0.0, 6);
    require(evidence.has_value(), "unblocked SER cone must find positive evidence");
    require(evidence->is_hudson == 1 && evidence->is_plevin == 1,
            "best SER evidence must satisfy both published criteria");
    require(evidence->phi >= 0.0 && evidence->phi < 360.0,
            "representative virtual hydrogen must retain deterministic phi");
}

void test_auto_cone_detects_methyl_without_explicit_hydrogen() {
    gemmi::Residue ala = residue("ALA", 2);
    ala.atoms = {
        atom("CA", "C", 11.53, 10.0, 13.0),
        atom("CB", "C", 10.0, 10.0, 13.0),
    };
    gemmi::Structure structure = structure_with_donor(std::move(ala));
    auto auto_hits = xhpi::detect_interactions(
        structure, true, 0.0, true, false);
    require(auto_hits.size() == 1, "auto Cone must detect one implicit ALA methyl hit");
    require(auto_hits.front().X_atom == "CB" && auto_hits.front().H_atom == "virt",
            "implicit methyl hit must retain donor identity and virtual-H marker");
    require(auto_hits.front().method == "Implicit/Cone",
            "methyl hit must be reported through the Cone route");

    auto explicit_only_hits = xhpi::detect_interactions(
        structure, false, 0.0, false, false);
    require(explicit_only_hits.empty(),
            "a hydrogen-free methyl group cannot pass explicit-only detection");
}

void test_auto_cone_is_independent_of_riding_hydrogen_direction() {
    gemmi::Residue ser_without_h = residue("SER", 2);
    ser_without_h.atoms = {
        atom("CB", "C", 11.42, 10.0, 13.0),
        atom("OG", "O", 10.0, 10.0, 13.0),
    };
    gemmi::Residue ser_with_away_h = ser_without_h;
    ser_with_away_h.atoms.push_back(atom("HG", "H", 10.0, 10.0, 13.972));

    auto no_h_hits = xhpi::detect_interactions(
        structure_with_donor(std::move(ser_without_h)), true, 0.0, true, false);
    auto away_h_hits = xhpi::detect_interactions(
        structure_with_donor(std::move(ser_with_away_h)), true, 0.0, true, false);
    require(no_h_hits.size() == 1 && away_h_hits.size() == 1,
            "auto Cone must evaluate SER with or without a riding hydrogen");
    require(no_h_hits.front().is_hudson == away_h_hits.front().is_hudson &&
            no_h_hits.front().is_plevin == away_h_hits.front().is_plevin &&
            std::abs(no_h_hits.front().theta - away_h_hits.front().theta) < 1e-10 &&
            std::abs(no_h_hits.front().angle_xh_pi - away_h_hits.front().angle_xh_pi) < 1e-10,
            "riding-H coordinates must not bias the binary Cone result");
}

void test_periodic_metal_contact_blocks_donor() {
    gemmi::Structure structure;
    structure.cell = gemmi::UnitCell(30, 30, 30, 90, 90, 90);
    structure.models.emplace_back(1);
    structure.models.back().chains.emplace_back("A");
    gemmi::Residue cys = residue("CYS", 1);
    cys.atoms.push_back(atom("SG", "S", 0.5, 10.0, 10.0));
    gemmi::Residue zinc = residue("ZN", 2);
    zinc.atoms.push_back(atom("ZN", "Zn", 29.0, 10.0, 10.0));
    structure.models.back().chains.back().residues.push_back(std::move(cys));
    structure.models.back().chains.back().residues.push_back(std::move(zinc));

    gemmi::Model& model = structure.models.back();
    gemmi::Atom& sulfur = model.chains.front().residues.front().atoms.front();
    gemmi::NeighborSearch ns(model, structure.cell, xhpi::DIST_SEARCH_LIMIT);
    ns.populate(true);
    require(xhpi::is_donor_blocked(sulfur, model, ns, sulfur.pos),
            "a periodic 1.5 A metal contact must block XH-pi donor assignment");
}

void test_equivalent_explicit_altloc_hits_are_deduplicated() {
    gemmi::Residue pro = residue("PRO", 2);
    pro.atoms = {
        atom("CD", "C", 10.0, 10.0, 13.0, 0.5, 'A'),
        atom("HD2", "H", 10.0, 10.0, 12.0, 0.5, 'A'),
        atom("CD", "C", 10.0, 10.0, 13.0, 0.5, 'B'),
        atom("HD2", "H", 10.0, 10.0, 12.0, 0.5, 'B'),
    };
    auto hits = explicit_hits(structure_with_donor(std::move(pro)));
    require(hits.size() == 1,
            "equivalent explicit altloc hits must produce one Moorhen interaction");
    require(std::abs(hits.front().combined_occ - 0.5) < 1e-6,
            "deduplication must retain the selected conformer occupancy");
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
    run_test("group-specific single-H geometry", test_group_specific_single_hydrogen_geometry);
    run_test("coupled methyl geometry", test_methyl_conformer_contains_three_coupled_hydrogens);
    run_test("complete methyl sterics", test_any_hydrogen_clash_invalidates_complete_methyl_conformer);
    run_test("valid H-bond contact sterics", test_valid_hbond_contact_is_not_a_steric_clash);
    run_test("minimal Cone acceptor typing", test_minimal_cone_acceptor_typing);
    run_test("Cone parent altloc resolution", test_cone_parent_altloc_resolution);
    run_test("deterministic binary Cone evidence", test_binary_cone_uses_deterministic_positive_evidence);
    run_test("auto Cone methyl routing", test_auto_cone_detects_methyl_without_explicit_hydrogen);
    run_test("auto Cone riding-H independence", test_auto_cone_is_independent_of_riding_hydrogen_direction);
    run_test("periodic metal donor blocking", test_periodic_metal_contact_blocks_donor);
    run_test("explicit altloc result deduplication", test_equivalent_explicit_altloc_hits_are_deduplicated);

    if (failures != 0) {
        std::cerr << failures << " contract test(s) failed\n";
        return 1;
    }
    std::cout << "All XPID Moorhen contract tests passed\n";
    return 0;
}
