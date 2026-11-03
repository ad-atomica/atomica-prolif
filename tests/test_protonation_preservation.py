"""A prepared source graph, including deliberately absent H, is authoritative."""

import pickle
from pathlib import Path
from typing import TypedDict

import numpy as np
import pytest
from rdkit import Chem

from prolif import Molecule
from prolif.interactions import Anionic, Cationic, ImplicitHBAcceptor, ImplicitHBDonor
from prolif.io import MoleculeStandardizer
from prolif.pickling import PROLIF_PICKLE_OPTIONS, RDKitPickleHandler

from .protonation_utils import (
    HYDROGEN_REPRESENTATIONS,
    atoms_by_name,
    graph_signature,
    local_state,
    peptide,
    represent_hydrogens,
    state_peptide,
)


class HistidineEvidence(TypedDict, total=False):
    histidine_states: dict[str, str]
    complete_histidine_hydrogens: list[str]


# Literal (formal charge, total H) expectations; no template-derived oracle.
STATES = [
    ("ASP:acid1", ("OD1", "OD2"), ((0, 1), (0, 0))),
    ("ASP:acid2", ("OD1", "OD2"), ((0, 0), (0, 1))),
    ("ASP:anion2", ("OD1", "OD2"), ((0, 0), (-1, 0))),
    ("GLU:acid2", ("OE1", "OE2"), ((0, 0), (0, 1))),
    ("GLU:anion1", ("OE1", "OE2"), ((-1, 0), (0, 0))),
    ("LYS:0", ("NZ",), ((0, 2),)),
    ("LYS:+", ("NZ",), ((1, 3),)),
    ("ARG:0", ("NE", "NH1", "NH2"), ((0, 1), (0, 1), (0, 2))),
    ("ARG:+", ("NE", "NH1", "NH2"), ((0, 1), (1, 2), (0, 2))),
    ("CYS:0", ("SG",), ((0, 1),)),
    ("CYS:-", ("SG",), ((-1, 0),)),
    ("TYR:-", ("OH",), ((-1, 0),)),
    ("HIS:HID", ("ND1", "NE2"), ((0, 1), (0, 0))),
    ("HIS:HIE", ("ND1", "NE2"), ((0, 0), (0, 1))),
    ("HIS:HIP", ("ND1", "NE2"), ((1, 1), (0, 1))),
]


def preserve(source: Chem.Mol) -> Molecule:
    return MoleculeStandardizer(preserve_protonation=True)(
        source, source_is_prepared=True
    )


@pytest.mark.parametrize("hydrogens", [2, 4])
@pytest.mark.parametrize("representation", HYDROGEN_REPRESENTATIONS)
def test_hypervalent_methionine_is_not_a_standard_residue(
    hydrogens: int, representation: str
) -> None:
    source = peptide("GMG")
    sulfur = atoms_by_name(source)["SD"]
    sulfur.SetNumExplicitHs(hydrogens)
    sulfur.SetNoImplicit(True)
    Chem.SanitizeMol(source)
    source = represent_hydrogens(source, representation)
    before = source.ToBinary()
    with pytest.raises(ValueError, match=r"MET2.A.*SD"):
        preserve(source)
    assert source.ToBinary() == before


@pytest.mark.parametrize("representation", HYDROGEN_REPRESENTATIONS)
def test_arg_imine_ammonium_is_not_a_guanidinium_resonance_form(
    representation: str,
) -> None:
    source = state_peptide("ARG:+")
    atoms = atoms_by_name(source)
    atoms["NH1"].SetFormalCharge(0)
    atoms["NH1"].SetNumExplicitHs(1)
    atoms["NH2"].SetFormalCharge(1)
    atoms["NH2"].SetNumExplicitHs(3)
    Chem.SanitizeMol(source)
    source = represent_hydrogens(source, representation)
    before = source.ToBinary()
    with pytest.raises(ValueError, match=r"ARG2.A.*arginine"):
        preserve(source)
    assert source.ToBinary() == before


@pytest.mark.parametrize("name", ["H3", "3H", "HN3", "HT3"])
def test_terminal_h_name_on_a_sidechain_is_rejected(name: str) -> None:
    source = state_peptide("HIS:HID", "graph")
    carbon = atoms_by_name(source)["CB"]
    hydrogen = next(a for a in carbon.GetNeighbors() if a.GetAtomicNum() == 1)
    hydrogen.GetPDBResidueInfo().SetName(name)
    before = source.ToBinary()
    with pytest.raises(ValueError, match=r"HIS2.A.*hydrogen name H3.*expected N"):
        preserve(source)
    assert source.ToBinary() == before


@pytest.mark.parametrize(
    ("state", "alias", "parent", "name"),
    [
        ("ASP:acid1", "ASP", "OD1", "HD2"),
        ("ASP:acid1", "ASH", "OD1", "HD2"),
        ("ASP:acid2", "ASH", "OD2", "HD1"),
        ("GLU:acid1", "GLU", "OE1", "HE2"),
        ("GLU:acid1", "GLH", "OE1", "HE2"),
        ("GLU:acid2", "GLH", "OE2", "HE1"),
    ],
)
def test_conflicting_acid_h_name_is_rejected(
    state: str, alias: str, parent: str, name: str
) -> None:
    source = state_peptide(state, "graph", alias)
    oxygen = atoms_by_name(source)[parent]
    hydrogen = next(a for a in oxygen.GetNeighbors() if a.GetAtomicNum() == 1)
    hydrogen.GetPDBResidueInfo().SetName(name)
    before = source.ToBinary()
    with pytest.raises(ValueError, match=r"conflicting hydrogen name"):
        preserve(source)
    assert source.ToBinary() == before
    # The chemically matching name must preserve this same proton location.
    hydrogen.GetPDBResidueInfo().SetName("H" + parent[1:])
    assert graph_signature(preserve(source)) == graph_signature(source)


@pytest.mark.parametrize("state", ["HIS:HID", "HIS:HIE", "HIS:HIP"])
def test_equivalent_kekule_input_retains_chemical_assignments(state: str) -> None:
    expected = state_peptide(state)
    source = Chem.Mol(expected)
    Chem.Kekulize(source, clearAromaticFlags=True)
    before = source.ToBinary()
    assert graph_signature(preserve(source)) == graph_signature(expected)
    assert source.ToBinary() == before


@pytest.mark.parametrize("representation", HYDROGEN_REPRESENTATIONS)
@pytest.mark.parametrize(
    ("state", "charges"),
    [
        ("HID", {"ND1": 1, "NE2": -1}),
        ("HIE", {"ND1": -1, "NE2": 1}),
        ("HIE", {"CE1": 1, "CD2": -1}),
        ("HIP", {"ND1": 0, "CE1": 1}),
    ],
)
def test_nonstandard_histidine_charge_distribution_is_rejected(
    state: str, charges: dict[str, int], representation: str
) -> None:
    source = state_peptide(f"HIS:{state}", alias=state)
    atoms = atoms_by_name(source)
    for name, charge in charges.items():
        atoms[name].SetFormalCharge(charge)
    Chem.SanitizeMol(source)
    source = represent_hydrogens(source, representation)
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in source.GetAtoms())
    ring = ("CG", "ND1", "CE1", "NE2", "CD2")
    assert sum(q for q, _ in local_state(source, ring)) == int(state == "HIP")
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match=rf"{state}2.A.*histidine"):
        preserve(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


@pytest.mark.parametrize("representation", HYDROGEN_REPRESENTATIONS)
@pytest.mark.parametrize("charged_nitrogen", ["ND1", "NE2"])
def test_histidinium_charge_on_either_nitrogen_is_preserved(
    charged_nitrogen: str, representation: str
) -> None:
    source = state_peptide("HIS:HIP", alias="HIP")
    atoms = atoms_by_name(source)
    for name in ("ND1", "NE2"):
        atoms[name].SetFormalCharge(int(name == charged_nitrogen))
    Chem.SanitizeMol(source)
    source = represent_hydrogens(source, representation)
    assert local_state(source, (charged_nitrogen,)) == ((1, 1),)
    assert graph_signature(preserve(source)) == graph_signature(source)


def test_preparation_declaration_is_required() -> None:
    source = Chem.MolFromSequence("GKG")
    with pytest.raises(ValueError, match=r"source_is_prepared.*True"):
        MoleculeStandardizer(preserve_protonation=True)(source)


@pytest.mark.parametrize("source", ["prepared.pdb", Path("prepared.pdb")])
def test_file_input_cannot_declare_a_prepared_graph(source: str | Path) -> None:
    with pytest.raises(ValueError, match=r"(?i)(prepared|RDKit|ProLIF)"):
        MoleculeStandardizer(preserve_protonation=True)(source, source_is_prepared=True)


def test_preparation_declaration_requires_preservation_mode() -> None:
    with pytest.raises(ValueError, match="preserve_protonation"):
        MoleculeStandardizer()(peptide("GKG"), source_is_prepared=True)


def test_broad_and_histidine_only_modes_conflict() -> None:
    with pytest.raises(ValueError, match="preserve_histidine"):
        MoleculeStandardizer(preserve_protonation=True, preserve_histidine=True)


@pytest.mark.parametrize(
    "evidence", ["histidine_states", "complete_histidine_hydrogens"]
)
def test_histidine_only_evidence_conflicts_with_prepared_graph(evidence: str) -> None:
    kwargs: HistidineEvidence = (
        {"histidine_states": {"HIS2.A": "HID"}}
        if evidence == "histidine_states"
        else {"complete_histidine_hydrogens": ["HIS2.A"]}
    )
    with pytest.raises(ValueError, match=r"(?i)(histidine|conflict)"):
        MoleculeStandardizer(preserve_protonation=True)(
            state_peptide("HIS:HID"), source_is_prepared=True, **kwargs
        )


@pytest.mark.parametrize("representation", HYDROGEN_REPRESENTATIONS)
@pytest.mark.parametrize(
    ("state", "names", "expected"), STATES, ids=[s[0] for s in STATES]
)
def test_prepared_sidechain_states_are_preserved(
    state: str, names: tuple[str, ...], expected: tuple, representation: str
) -> None:
    source = state_peptide(state, representation)
    assert local_state(source, names) == expected
    before = graph_signature(source)
    result = preserve(source)
    assert local_state(result, names) == expected
    assert local_state(result[1], names) == expected
    assert graph_signature(result) == before
    assert graph_signature(source) == before


@pytest.mark.parametrize(
    ("state", "names", "expected"),
    [
        ("GLU:acid1", ("OE1", "OE2"), ((0, 1), (0, 0))),
        ("ASP:anion1", ("OD1", "OD2"), ((-1, 0), (0, 0))),
        ("GLU:anion2", ("OE1", "OE2"), ((0, 0), (-1, 0))),
        ("ARG:0NE", ("NE", "NH1", "NH2"), ((0, 0), (0, 2), (0, 2))),
        ("ARG:+NE", ("NE", "NH1", "NH2"), ((1, 1), (0, 2), (0, 2))),
        ("SER:-", ("OG",), ((-1, 0),)),
        ("THR:-", ("OG1",), ((-1, 0),)),
    ],
)
def test_alternative_resonance_and_less_common_states(
    state: str, names: tuple[str, ...], expected: tuple
) -> None:
    source = state_peptide(state)
    assert local_state(source, names) == expected
    result = preserve(source)
    assert local_state(result, names) == expected
    assert graph_signature(result) == graph_signature(source)


@pytest.mark.parametrize("code", list("ACDEFGHIKLMNPQRSTVWY"))
def test_every_standard_amino_acid_is_supported(code: str) -> None:
    source = represent_hydrogens(peptide(f"G{code}G"), "counted")
    result = preserve(source)
    assert graph_signature(result) == graph_signature(source)


@pytest.mark.parametrize(
    ("state", "alias"),
    [
        ("LYS:0", "LYN"),
        ("LYS:0", "LSN"),
        ("ARG:0", "ARGN"),
        ("ASP:acid1", "ASH"),
        ("GLU:acid2", "GLH"),
        ("CYS:-", "CYM"),
        ("CYS:-", "CYN"),
        ("CYS:-", "CYX"),
        ("TYR:-", "TYM"),
        ("HIS:HID", "HSD"),
        ("HIS:HIE", "HISB"),
        ("HIS:HIP", "HISH"),
    ],
)
def test_explicit_alias_agrees_with_prepared_state(state: str, alias: str) -> None:
    source = state_peptide(state, alias=alias)
    result = preserve(source)
    assert str(result[1].resid) == f"{alias}2.A"
    assert graph_signature(result) == graph_signature(source)


@pytest.mark.parametrize(
    ("state", "alias"),
    [
        ("LYS:+", "LYN"),
        ("LYS:+", "LSN"),
        ("ARG:+", "ARGN"),
        ("ASP:anion2", "ASH"),
        ("GLU:anion1", "GLH"),
        ("CYS:0", "CYM"),
        ("CYS:0", "CYN"),
        ("CYS:0", "CYX"),
        ("TYR:0", "TYM"),
        ("HIS:HIE", "HID"),
    ],
)
def test_alias_conflict_rejects_without_mutation(state: str, alias: str) -> None:
    source = Molecule(state_peptide(state, alias=alias))
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match=rf"{alias}2.A"):
        preserve(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


@pytest.mark.parametrize("input_kind", ["rdkit", "prolif"])
def test_prepared_graph_is_returned_as_an_independent_copy(input_kind: str) -> None:
    source = state_peptide("LYS:0", "graph", "LYN")
    atoms_by_name(source)["CA"].SetIsotope(13)
    # Non-template atom order and multiple conformers must survive.
    source = Chem.RenumberAtoms(source, list(reversed(range(source.GetNumAtoms()))))
    conf = Chem.Conformer(source.GetConformer())
    conf.SetAtomPosition(0, (5.25, 2.5, -1.75))
    source.AddConformer(conf, assignId=True)
    if input_kind == "prolif":
        source = Molecule(source)
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    signature = graph_signature(source)
    result = preserve(source)
    assert result is not source
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before
    assert graph_signature(result) == signature
    for residue in result.residues.values():
        for atom in residue.GetAtoms():
            parent = result.GetAtomWithIdx(atom.GetUnsignedProp("mapindex"))
            assert atom.GetFormalCharge() == parent.GetFormalCharge()
            assert atom.GetTotalNumHs(includeNeighbors=True) == parent.GetTotalNumHs(
                includeNeighbors=True
            )


@pytest.mark.parametrize("representation", HYDROGEN_REPRESENTATIONS)
def test_mixed_peptide_with_disulfide_and_charged_termini(representation: str) -> None:
    source = peptide(
        "GDEKRCCYHSTPG",
        {
            2: ("ASP:acid1", "ASH"),
            3: ("GLU:anion2", None),
            4: ("LYS:0", "LSN"),
            5: ("ARG:+NE", None),
            6: ("CYS:0", "CYX"),
            7: ("CYS:0", None),
            8: ("TYR:-", "TYM"),
            9: ("HIS:HID", "HSD"),
            10: ("SER:-", None),
            11: ("THR:0", None),
        },
    )
    source = Chem.RWMol(source)
    sulfurs = [atoms_by_name(source, number)["SG"] for number in (6, 7)]
    for atom in sulfurs:
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
    source.AddBond(sulfurs[0].GetIdx(), sulfurs[1].GetIdx(), Chem.BondType.SINGLE)
    nterm = atoms_by_name(source, 1)["N"]
    nterm.SetFormalCharge(1)
    nterm.SetNoImplicit(True)
    nterm.SetNumExplicitHs(3)
    cterm = atoms_by_name(source, 13)["OXT"]
    cterm.SetFormalCharge(-1)
    cterm.SetNoImplicit(True)
    cterm.SetNumExplicitHs(0)
    Chem.SanitizeMol(source)
    source = represent_hydrogens(source, representation)
    result = preserve(source)
    assert graph_signature(result) == graph_signature(source)
    for number in (6, 7):
        assert local_state(result, ("SG",), number) == ((0, 0),)
    assert local_state(result, ("N",), 1) == ((1, 3),)
    assert local_state(result, ("OXT",), 13) == ((-1, 0),)
    assert local_state(result[3], ("NZ",), 4) == ((0, 2),)


@pytest.mark.parametrize("charge", [0, 1])
def test_proline_n_terminus_preserves_neutral_and_cationic_states(charge: int) -> None:
    source = peptide("PG")
    nitrogen = atoms_by_name(source, 1)["N"]
    nitrogen.SetFormalCharge(charge)
    nitrogen.SetNoImplicit(True)
    nitrogen.SetNumExplicitHs(1 + charge)
    Chem.SanitizeMol(source)
    result = preserve(source)
    assert local_state(result, ("N",), 1) == ((charge, 1 + charge),)
    assert graph_signature(result) == graph_signature(source)


@pytest.mark.parametrize("single_oxygen", ["O", "OXT"])
@pytest.mark.parametrize("charge", [0, -1])
def test_c_terminal_acid_and_anion_preserve_oxygen_position(
    single_oxygen: str, charge: int
) -> None:
    source = peptide("AG")
    atoms = atoms_by_name(source)
    for name in ("O", "OXT"):
        single = name == single_oxygen
        atom = atoms[name]
        atom.SetFormalCharge(charge if single else 0)
        atom.SetNumExplicitHs(int(single and charge == 0))
        atom.SetNoImplicit(True)
        source.GetBondBetweenAtoms(atoms["C"].GetIdx(), atom.GetIdx()).SetBondType(
            Chem.BondType.SINGLE if single else Chem.BondType.DOUBLE
        )
    Chem.SanitizeMol(source)
    result = preserve(source)
    assert local_state(result, (single_oxygen,)) == ((charge, int(charge == 0)),)
    assert graph_signature(result) == graph_signature(source)


def test_missing_terminal_oxygen_is_not_treated_as_prepared_aldehyde() -> None:
    source = Chem.RWMol(peptide("AG"))
    source.RemoveAtom(atoms_by_name(source)["OXT"].GetIdx())
    Chem.SanitizeMol(source)
    with pytest.raises(ValueError, match=r"GLY2.A"):
        preserve(source)


@pytest.mark.parametrize("unsupported_name", ["HOH", "ACE", "NME"])
def test_nonprotein_residues_require_a_separate_preparation_path(
    unsupported_name: str,
) -> None:
    source = peptide("AG")
    for atom in atoms_by_name(source, 1).values():
        atom.GetPDBResidueInfo().SetResidueName(unsupported_name)
    with pytest.raises(ValueError, match=rf"{unsupported_name}1.A"):
        preserve(source)


def test_sanitization_cannot_create_radicals_in_the_prepared_source() -> None:
    source = state_peptide("LYS:0", "fixed-counted")
    nitrogen = atoms_by_name(source)["NZ"]
    nitrogen.SetNumExplicitHs(1)
    source.UpdatePropertyCache(strict=True)
    assert nitrogen.GetNumRadicalElectrons() == 0
    assert local_state(source, ("NZ",)) == ((0, 1),)
    # Establish the defect on a private copy without normalizing the declared source.
    sanitized = Chem.Mol(source)
    Chem.SanitizeMol(sanitized)
    assert atoms_by_name(sanitized)["NZ"].GetNumRadicalElectrons() == 1
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match=r"(?i)(sanitiz|radical)"):
        preserve(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


@pytest.mark.parametrize(
    ("code", "resname", "carbon", "oxygen", "nitrogen"),
    [("N", "ASN", "CG", "OD1", "ND2"), ("Q", "GLN", "CD", "OE1", "NE2")],
)
def test_sidechain_amide_cannot_be_replaced_by_an_imidic_acid(
    code: str, resname: str, carbon: str, oxygen: str, nitrogen: str
) -> None:
    source = represent_hydrogens(peptide(f"G{code}G"), "fixed-counted")
    atoms = atoms_by_name(source)
    for name, order in (
        (oxygen, Chem.BondType.SINGLE),
        (nitrogen, Chem.BondType.DOUBLE),
    ):
        atoms[name].SetNumExplicitHs(1)
        source.GetBondBetweenAtoms(
            atoms[carbon].GetIdx(), atoms[name].GetIdx()
        ).SetBondType(order)
    Chem.SanitizeMol(source)
    assert local_state(source, (oxygen, nitrogen)) == ((0, 1), (0, 1))
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in source.GetAtoms())
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match=rf"{resname}2.A"):
        preserve(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


def test_peptide_bond_cannot_be_replaced_by_an_imidic_acid() -> None:
    source = represent_hydrogens(peptide("GAG"), "fixed-counted")
    first, second = atoms_by_name(source, 1), atoms_by_name(source, 2)
    first["O"].SetNumExplicitHs(1)
    second["N"].SetNumExplicitHs(0)
    source.GetBondBetweenAtoms(first["C"].GetIdx(), first["O"].GetIdx()).SetBondType(
        Chem.BondType.SINGLE
    )
    source.GetBondBetweenAtoms(first["C"].GetIdx(), second["N"].GetIdx()).SetBondType(
        Chem.BondType.DOUBLE
    )
    Chem.SanitizeMol(source)
    assert local_state(source, ("O",), 1) == ((0, 1),)
    assert local_state(source, ("N",), 2) == ((0, 0),)
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in source.GetAtoms())
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match=r"GLY1.A.*ALA2.A"):
        preserve(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


def test_generic_cysteine_state_is_independent_of_an_unrelated_cym_alias() -> None:
    source = represent_hydrogens(
        peptide("GCCG", {2: ("CYS:0", None), 3: ("CYS:-", "CYM")}), "graph"
    )
    result = preserve(source)
    assert local_state(result, ("SG",), 2) == ((0, 1),)
    assert local_state(result, ("SG",), 3) == ((-1, 0),)
    assert graph_signature(result) == graph_signature(source)


def test_named_hydrogens_on_protonated_n_terminal_arginine() -> None:
    source = peptide("RG", {1: ("ARG:+", None)})
    nitrogen = atoms_by_name(source, 1)["N"]
    nitrogen.SetFormalCharge(1)
    nitrogen.SetNumExplicitHs(3)
    nitrogen.SetNoImplicit(True)
    Chem.SanitizeMol(source)
    source = represent_hydrogens(source, "graph")
    nitrogen = atoms_by_name(source, 1)["N"]
    hydrogens = [atom for atom in nitrogen.GetNeighbors() if atom.GetAtomicNum() == 1]
    for atom, name in zip(hydrogens, ("H", "H2", "H3"), strict=True):
        atom.GetPDBResidueInfo().SetName(name)
    assert local_state(source, ("N",), 1) == ((1, 3),)
    assert local_state(source, ("NE", "NH1", "NH2"), 1) == ((0, 1), (1, 2), (0, 2))
    result = preserve(source)
    assert graph_signature(result) == graph_signature(source)


@pytest.mark.parametrize(
    "damage",
    [
        "missing-heavy",
        "extra-heavy",
        "wrong-bond-order",
        "wrong-topology",
        "extra-bond",
        "crosslink",
        "radical",
        "unknown-residue",
        "metal-link",
    ],
)
def test_unsupported_chemistry_rejects_without_mutating_source(damage: str) -> None:
    source = Chem.RWMol(state_peptide("LYS:0"))
    atoms = atoms_by_name(source)
    if damage == "missing-heavy":
        source.RemoveAtom(atoms["NZ"].GetIdx())
    elif damage == "extra-heavy":
        fluorine = Chem.Atom("F")
        fluorine.SetMonomerInfo(
            Chem.AtomPDBResidueInfo(
                " FX ", residueName="LYS", residueNumber=2, chainId="A"
            )
        )
        index = source.AddAtom(fluorine)
        source.AddBond(atoms["CE"].GetIdx(), index, Chem.BondType.SINGLE)
        atoms["CE"].SetNumExplicitHs(1)
    elif damage == "wrong-bond-order":
        source.GetBondBetweenAtoms(
            atoms["CE"].GetIdx(), atoms["NZ"].GetIdx()
        ).SetBondType(Chem.BondType.DOUBLE)
        atoms["CE"].SetNumExplicitHs(1)
        atoms["NZ"].SetNumExplicitHs(1)
    elif damage == "wrong-topology":
        source.RemoveBond(atoms["CE"].GetIdx(), atoms["NZ"].GetIdx())
        source.AddBond(atoms["CD"].GetIdx(), atoms["NZ"].GetIdx(), Chem.BondType.SINGLE)
        atoms["CE"].SetNumExplicitHs(3)
        atoms["CD"].SetNumExplicitHs(1)
    elif damage in {"extra-bond", "crosslink"}:
        first = atoms["CB"]
        second = (
            atoms["CD"] if damage == "extra-bond" else atoms_by_name(source, 1)["CA"]
        )
        first.SetNumExplicitHs(1)
        second.SetNumExplicitHs(1)
        source.AddBond(first.GetIdx(), second.GetIdx(), Chem.BondType.SINGLE)
    elif damage == "radical":
        atoms["NZ"].SetNumExplicitHs(1)
        atoms["NZ"].SetNumRadicalElectrons(1)
    elif damage == "unknown-residue":
        for atom in atoms.values():
            atom.GetPDBResidueInfo().SetResidueName("UNK")
    elif damage == "metal-link":
        zinc = Chem.Atom("Zn")
        zinc.SetFormalCharge(2)
        zinc.SetMonomerInfo(
            Chem.AtomPDBResidueInfo(
                " ZN ", residueName="ZN", residueNumber=9, chainId="A"
            )
        )
        source.AddBond(atoms["NZ"].GetIdx(), source.AddAtom(zinc), Chem.BondType.DATIVE)
    source.UpdatePropertyCache(strict=False)
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match=r"(LYS2.A|UNK2.A|ZN9.A)"):
        preserve(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


@pytest.mark.parametrize("name", ["HG", "HB2"])
def test_recognized_hydrogen_name_must_agree_with_attachment(name: str) -> None:
    source = state_peptide("CYS:0", "graph")
    heavy = atoms_by_name(source)["SG" if name == "HB2" else "CB"]
    hydrogen = next(atom for atom in heavy.GetNeighbors() if atom.GetAtomicNum() == 1)
    hydrogen.GetPDBResidueInfo().SetName(name)
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match=r"CYS2.A.*(?i:hydrogen|attachment|conflict)"):
        preserve(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


def test_unbonded_hydrogen_is_rejected() -> None:
    source = Chem.RWMol(state_peptide("CYS:0", "graph"))
    sulfur = atoms_by_name(source)["SG"]
    hydrogen = next(atom for atom in sulfur.GetNeighbors() if atom.GetAtomicNum() == 1)
    source.RemoveBond(sulfur.GetIdx(), hydrogen.GetIdx())
    with pytest.raises(ValueError, match=r"CYS2.A"):
        preserve(source)


def test_named_hydrogens_in_peptide_with_termini_and_disulfide() -> None:
    source = Chem.RWMol(peptide("CKDC", {2: ("LYS:0", "LYN"), 3: ("ASP:acid1", "ASH")}))
    sulfurs = [atoms_by_name(source, number)["SG"] for number in (1, 4)]
    for atom in sulfurs:
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
    source.AddBond(sulfurs[0].GetIdx(), sulfurs[1].GetIdx(), Chem.BondType.SINGLE)
    nitrogen = atoms_by_name(source, 1)["N"]
    nitrogen.SetFormalCharge(1)
    nitrogen.SetNumExplicitHs(3)
    nitrogen.SetNoImplicit(True)
    Chem.SanitizeMol(source)
    source = represent_hydrogens(source, "graph")
    targets = [
        (1, "N", ("H1", "H2", "H3")),
        (2, "NZ", ("HZ1", "HZ2")),
        (3, "OD1", ("HD1",)),
        (4, "N", ("H",)),
    ]
    for number, heavy_name, hydrogen_names in targets:
        heavy = atoms_by_name(source, number)[heavy_name]
        hydrogens = [atom for atom in heavy.GetNeighbors() if atom.GetAtomicNum() == 1]
        for atom, name in zip(hydrogens, hydrogen_names, strict=True):
            atom.GetPDBResidueInfo().SetName(name)
    result = preserve(source)
    assert graph_signature(result) == graph_signature(source)
    assert local_state(result, ("NZ",), 2) == ((0, 2),)
    assert local_state(result, ("OD1", "OD2"), 3) == ((0, 1), (0, 0))


def contact_ligand(source: Chem.Mol, target: str, smiles: str = "CO") -> Molecule:
    """Two-atom ligand placed at a hand-selected ideal contact angle, 2.8 A away."""
    atom = atoms_by_name(source)[target]
    neighbor = next(a for a in atom.GetNeighbors() if a.GetAtomicNum() != 1)
    positions = source.GetConformer().GetPositions()
    center = positions[atom.GetIdx()]
    toward_neighbor = positions[neighbor.GetIdx()] - center
    toward_neighbor /= np.linalg.norm(toward_neighbor)
    perpendicular = np.array([-toward_neighbor[1], toward_neighbor[0], 0.0])
    theta = np.radians(
        120 if atom.GetHybridization() == Chem.HybridizationType.SP2 else 109.5
    )
    outward = np.cos(theta) * toward_neighbor + np.sin(theta) * perpendicular
    probe = center + 2.8 * outward
    toward_target = -outward
    probe_perpendicular = np.array([-toward_target[1], toward_target[0], 0.0])
    carbon = probe + 1.4 * (
        np.cos(np.radians(109.5)) * toward_target
        + np.sin(np.radians(109.5)) * probe_perpendicular
    )
    ligand = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(ligand.GetNumAtoms())
    conf.SetAtomPosition(0, carbon)
    conf.SetAtomPosition(1, probe)
    ligand.AddConformer(conf)
    return Molecule.from_rdkit(ligand)


@pytest.mark.parametrize(
    ("state", "alias", "target", "accepts", "donates", "ionic"),
    [
        ("LYS:0", "LYN", "NZ", True, True, False),
        ("LYS:+", "LYS", "NZ", False, True, True),
        ("CYS:0", "CYS", "SG", False, True, False),
        ("CYS:-", "CYS", "SG", False, False, True),
        ("ASP:acid1", "ASH", "OD1", True, True, False),
        ("ASP:anion1", "ASP", "OD1", True, False, True),
    ],
)
def test_detector_geometry_and_ionic_metadata_survive_h_removal(
    state: str, alias: str, target: str, accepts: bool, donates: bool, ionic: bool
) -> None:
    source = state_peptide(state, "graph", alias)
    methanol = contact_ligand(source, target)[0]
    charged_ligand = contact_ligand(
        source, target, "C[O-]" if state.startswith("LYS") else "C[NH3+]"
    )[0]
    detectors = (
        (ImplicitHBDonor(), methanol, accepts),
        (ImplicitHBAcceptor(), methanol, donates),
        (Anionic() if state.startswith("LYS") else Cationic(), charged_ligand, ionic),
    )
    prepared = preserve(source)
    proteins = (
        Molecule(source),
        prepared,
        Molecule(Chem.RemoveHs(prepared)),
        Molecule(Chem.RemoveHs(prepared, sanitize=False)),
    )
    for detector, ligand, expected in detectors:
        outcomes = []
        for protein in proteins:
            residue = protein[1]
            target_index = atoms_by_name(residue)[target].GetIdx()
            contacts = [
                metadata
                for metadata in detector.detect(ligand, residue)
                if metadata["indices"]["protein"] == (target_index,)
                and metadata["indices"]["ligand"] == (1,)
            ]
            assert bool(contacts) is expected
            for metadata in contacts:
                assert metadata["distance"] == pytest.approx(2.8)
                if isinstance(detector, ImplicitHBAcceptor | ImplicitHBDonor):
                    assert "ideal_donor_angle" in metadata
                    assert "ideal_acceptor_angle" in metadata
                    assert "vina_hbond_potential" in metadata
            outcomes.append(contacts)
        assert outcomes[1:] == [outcomes[0]] * 3


def test_serialization_preserves_prepared_state_and_detector_geometry() -> None:
    result = preserve(state_peptide("LYS:0", "graph", "LYN"))
    ligand = contact_ligand(result, "NZ")[0]
    before = list(ImplicitHBAcceptor().detect(ligand, result[1]))
    assert before
    handler = RDKitPickleHandler(PROLIF_PICKLE_OPTIONS)
    handler.set()
    try:
        restored = pickle.loads(pickle.dumps(result))
    finally:
        handler.reset()
    assert local_state(restored[1], ("NZ",)) == ((0, 2),)
    after = list(ImplicitHBAcceptor().detect(ligand, restored[1]))
    assert len(after) == len(before)
    for expected, actual in zip(before, after, strict=True):
        assert actual.keys() == expected.keys()
        for key, value in expected.items():
            if key in {"indices", "parent_indices"}:
                assert actual[key] == value
            else:
                # RDKit's float32 pickle coordinates introduce angle differences
                # of order 1e-5 degrees; identities and contact counts stay exact.
                assert actual[key] == pytest.approx(value, abs=1e-4)
