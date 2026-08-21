# XPID Moorhen detection contract

This C++ implementation is a focused detector for transferring XH-pi hits to
Moorhen. It intentionally follows only the chemistry that changes the identity
or geometry label of a hit.

## Required behaviour

- Aromatic acceptors are complete chordless five- or six-membered cycles from
  CCP4 monomer aromatic bonds. PHE, TYR, TRP and HIS have complete built-in
  fallbacks; non-standard components do not.
- Alternate conformations are evaluated as chemically compatible complete
  states. Blank atoms may be shared, but labelled A/B states are never mixed.
  A lower-occupancy state is not removed before detection.
- Explicit H/D atoms require an element-appropriate covalent distance and an
  unambiguous owner. When a usable monomer dictionary is available, its bond
  table is authoritative. Neutron D names are matched to their H equivalents.
- Lys/Arg cation-pi atoms are excluded from XH-pi detection.
- Default output remains Hudson/Plevin XH-pi data for consumption by Moorhen.

## Deliberately out of scope

SASA, cooperativity, P-slab development metrics, batch processing, reporting,
PDB download/resolution, provenance files and extended H-bond annotations are
not part of this C++ detector.
