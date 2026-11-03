"""Source evidence precedes template selection, including in implicit analyses."""

import pickle
from pathlib import Path
from typing import TypedDict

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

from prolif import Molecule
from prolif.interactions import Anionic, ImplicitHBAcceptor, ImplicitHBDonor
from prolif.interactions.base import Interaction
from prolif.io import MoleculeStandardizer
from prolif.pickling import PROLIF_PICKLE_OPTIONS, RDKitPickleHandler
from prolif.residue import Residue, ResidueId


class EvidenceInputs(TypedDict, total=False):
    histidine_states: dict[str | ResidueId, str]
    complete_histidine_hydrogens: list[str | ResidueId]


def histidine(
    state: str = "HIE", name: str = "HIS", *, hydrogens: bool = True, number: int = 1
) -> Chem.Mol:
    """Free histidine with a planar imidazole and only ring-N graph hydrogens."""
    mol = Chem.MolFromSequence("H")
    counts = {"HID": (1, 0), "HIE": (0, 1), "HIP": (1, 1)}[state]
    for index, count in zip((6, 9), counts, strict=True):
        atom = mol.GetAtomWithIdx(index)
        atom.SetNumExplicitHs(count)
        atom.SetFormalCharge(int(state == "HIP" and index == 6))
    Chem.SanitizeMol(mol)
    if hydrogens:
        mol = Chem.AddHs(mol, onlyOnAtoms=[6, 9], addResidueInfo=True)
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        info.SetResidueName(name)
        info.SetResidueNumber(number)
        info.SetChainId("A")
        info.SetSerialNumber(atom.GetIdx() + 1)
        if atom.GetAtomicNum() == 1:
            nitrogen = atom.GetNeighbors()[0].GetPDBResidueInfo().GetName().strip()
            info.SetName(" HD1" if nitrogen == "ND1" else " HE2")
    Compute2DCoords(mol)
    return mol


def ring_state(mol: Chem.Mol) -> tuple[int, int, int]:
    atoms = {a.GetPDBResidueInfo().GetName().strip(): a for a in mol.GetAtoms()}
    return (
        atoms["ND1"].GetTotalNumHs(includeNeighbors=True),
        atoms["NE2"].GetTotalNumHs(includeNeighbors=True),
        sum(atoms[n].GetFormalCharge() for n in ("CG", "ND1", "CE1", "NE2", "CD2")),
    )


EXPECTED = {"HID": (1, 0, 0), "HIE": (0, 1, 0), "HIP": (1, 1, 1)}


def test_complete_mixed_his_retains_source_states() -> None:
    mols = [histidine(state, number=i) for i, state in enumerate(EXPECTED, 1)]
    source = Molecule(Chem.CombineMols(Chem.CombineMols(*mols[:2]), mols[2]))
    result = MoleculeStandardizer(preserve_histidine=True)(
        source, complete_histidine_hydrogens=["HIS1.A", "HIS2.A", "HIS3.A"]
    )
    assert [ring_state(r) for r in result.residues.values()] == list(EXPECTED.values())
    for i, state in enumerate(EXPECTED, 1):
        parent_atoms = [
            a.GetIdx()
            for a in result.GetAtoms()
            if a.GetPDBResidueInfo().GetResidueNumber() == i
        ]
        parent_fragment = Chem.PathToSubmol(
            result,
            [
                b.GetIdx()
                for b in result.GetBonds()
                if b.GetBeginAtomIdx() in parent_atoms
                and b.GetEndAtomIdx() in parent_atoms
            ],
        )
        assert ring_state(parent_fragment) == EXPECTED[state]


@pytest.mark.parametrize(
    ("name", "state"),
    [
        ("HID", "HID"),
        ("HSD", "HID"),
        ("HIE", "HIE"),
        ("HSE", "HIE"),
        ("HIP", "HIP"),
        ("HSP", "HIP"),
        ("HISA", "HID"),
        ("HISB", "HIE"),
        ("HISH", "HIP"),
    ],
)
@pytest.mark.parametrize("hydrogens", [True, False])
def test_explicit_alias_is_evidence(name: str, state: str, hydrogens: bool) -> None:
    result = MoleculeStandardizer(preserve_histidine=True)(
        histidine(state, name, hydrogens=hydrogens)
    )
    assert ring_state(result[0]) == EXPECTED[state]
    assert str(result[0].resid) == f"{name}1.A"


@pytest.mark.parametrize("state", EXPECTED)
def test_explicit_caller_assignment(state: str) -> None:
    result = MoleculeStandardizer(preserve_histidine=True)(
        histidine(state, hydrogens=False),
        histidine_states={ResidueId("HIS", 1, "A"): state},
    )
    assert ring_state(result["HIS1.A"]) == EXPECTED[state]


@pytest.mark.parametrize("hydrogens", [True, False])
def test_generic_his_requires_evidence(hydrogens: bool) -> None:
    with pytest.raises(
        ValueError,
        match=r"HIS1.A.*unresolved.*histidine_states.*complete_histidine_hydrogens",
    ):
        MoleculeStandardizer(preserve_histidine=True)(histidine(hydrogens=hydrogens))


@pytest.mark.parametrize("evidence", ["name", "caller", "complete"])
def test_conflicting_evidence_is_rejected(evidence: str) -> None:
    mol = histidine("HIE", "HID" if evidence == "name" else "HIS")
    kwargs: EvidenceInputs = (
        {"histidine_states": {"HIS1.A": "HID"}} if evidence != "name" else {}
    )
    if evidence == "complete":
        kwargs["complete_histidine_hydrogens"] = ["HIS1.A"]
    with pytest.raises(ValueError, match=r"HI[SD]1.A.*conflict"):
        MoleculeStandardizer(preserve_histidine=True)(mol, **kwargs)


def test_alias_and_caller_conflict_without_graph_h() -> None:
    with pytest.raises(ValueError, match="conflict"):
        MoleculeStandardizer(preserve_histidine=True)(
            histidine("HID", "HID", hydrogens=False), histidine_states={"HID1.A": "HIE"}
        )


def test_partial_inventory_can_corroborate_hip_but_not_establish_hid() -> None:
    source = histidine("HID")
    result = MoleculeStandardizer(preserve_histidine=True)(
        source, histidine_states={"HIS1.A": "HIP"}
    )
    assert ring_state(result[0]) == EXPECTED["HIP"]
    with pytest.raises(ValueError, match="conflict"):
        MoleculeStandardizer(preserve_histidine=True)(
            source,
            histidine_states={"HIS1.A": "HIP"},
            complete_histidine_hydrogens=["HIS1.A"],
        )


@pytest.mark.parametrize(
    "damage", ["wrong-attachment", "unbonded", "duplicate", "missing-ring"]
)
def test_invalid_inventory_is_rejected(damage: str) -> None:
    mol = Chem.RWMol(histidine("HID"))
    h = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
    if damage in {"wrong-attachment", "unbonded"}:
        mol.RemoveBond(6, h)
        if damage == "wrong-attachment":
            mol.AddBond(9, h, Chem.BondType.SINGLE)
    elif damage == "duplicate":
        mol.GetAtomWithIdx(8).GetPDBResidueInfo().SetName(" ND1")
    else:
        mol.RemoveAtom(8)
    with pytest.raises(ValueError, match=r"HIS1.A.*(attachment|duplicate|ring)"):
        MoleculeStandardizer(preserve_histidine=True)(
            mol, complete_histidine_hydrogens=["HIS1.A"]
        )


@pytest.mark.parametrize("hname", ["HND", "HND1", "HD1"])
def test_supported_hydrogen_name_aliases(hname: str) -> None:
    mol = histidine("HID")
    mol.GetAtomWithIdx(11).GetPDBResidueInfo().SetName(hname)
    assert (
        ring_state(
            MoleculeStandardizer(preserve_histidine=True)(
                mol, complete_histidine_hydrogens=["HIS1.A"]
            )[0]
        )
        == EXPECTED["HID"]
    )


@pytest.mark.parametrize("input_kind", ["rdkit", "prolif"])
def test_source_copy_identity_mapping_and_coordinates(input_kind: str) -> None:
    source = histidine("HID", "HSD")
    if input_kind == "prolif":
        source = Molecule(source)
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    result = MoleculeStandardizer(preserve_histidine=True)(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before
    assert result is not source
    assert list(map(str, result.residues)) == ["HSD1.A"]
    np.testing.assert_array_equal(
        result.GetConformer().GetPositions(), source.GetConformer().GetPositions()
    )
    assert [
        (a.GetAtomicNum(), a.GetPDBResidueInfo().GetName()) for a in result.GetAtoms()
    ] == [
        (a.GetAtomicNum(), a.GetPDBResidueInfo().GetName()) for a in source.GetAtoms()
    ]
    assert [a.GetUnsignedProp("mapindex") for a in result[0].GetAtoms()] == list(
        range(source.GetNumAtoms())
    )


def test_failure_does_not_mutate_source() -> None:
    source = Molecule(histidine("HIE", "HID"))
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match="conflict"):
        MoleculeStandardizer(preserve_histidine=True)(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before
    assert str(source[0].resid) == "HID1.A"


@pytest.mark.parametrize("name", ["HISA", "HISB", "HISH"])
def test_four_letter_pdb_alias_survives_owned_parser(tmp_path: Path, name: str) -> None:
    state = {"HISA": "HID", "HISB": "HIE", "HISH": "HIP"}[name]
    mol = histidine(state, hydrogens=False)
    # Extended PDB uses the otherwise blank column after the residue name.
    lines = Chem.MolToPDBBlock(mol).splitlines(keepends=True)
    block = "".join(
        line[:17] + name + line[21:] if line.startswith(("ATOM", "HETATM")) else line
        for line in lines
    )
    path = tmp_path / "histidine.pdb"
    path.write_text(block)
    result = MoleculeStandardizer(preserve_histidine=True)(path)
    assert str(result[0].resid) == f"{name}1.A"
    assert ring_state(result[0]) == EXPECTED[state]


def test_evidence_arguments_cannot_be_silently_ignored() -> None:
    with pytest.raises(ValueError, match="preserve_histidine=True"):
        MoleculeStandardizer()(histidine(), histidine_states={"HIS1.A": "HIE"})


def test_legacy_default_keeps_name_normalization_and_in_place_contract() -> None:
    source = Molecule(histidine("HID", "HSD"))
    assert MoleculeStandardizer()(source) is source
    assert str(source[0].resid) == "HID1.A"
    assert ring_state(MoleculeStandardizer()(histidine("HIE"))[0]) == EXPECTED["HIP"]


def contact_ligand(protein: Molecule, smiles: str) -> Molecule:
    """Put the terminal O 2.9 A along ND1's outward in-plane bisector."""
    mol = Chem.MolFromSmiles(smiles)
    mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()))
    nitrogen = next(
        a
        for a in protein.GetAtoms()
        if a.GetPDBResidueInfo().GetName().strip() == "ND1"
    )
    positions = protein.GetConformer().GetPositions()
    origin = positions[nitrogen.GetIdx()]
    bonds = np.array(
        [
            positions[a.GetIdx()] - origin
            for a in nitrogen.GetNeighbors()
            if a.GetAtomicNum() != 1
        ]
    )
    bonds /= np.linalg.norm(bonds, axis=1)[:, None]
    direction = -bonds.sum(axis=0)
    direction /= np.linalg.norm(direction)
    perpendicular = np.cross(direction, [0, 0, 1])
    oxygen = origin + 2.9 * direction
    for i in range(mol.GetNumAtoms()):
        position = oxygen
        if i != mol.GetNumAtoms() - 1:
            position = oxygen + 1.43 * (
                -np.cos(np.radians(109.5)) * direction
                + np.sin(np.radians(109.5)) * perpendicular
            )
        mol.GetConformer().SetAtomPosition(i, position)
    return Molecule.from_rdkit(mol)


@pytest.mark.parametrize("remove_h", [False, True])
@pytest.mark.parametrize("state", EXPECTED)
def test_real_contacts_preserve_assignment_after_h_removal(
    state: str, remove_h: bool
) -> None:
    protein = MoleculeStandardizer(preserve_histidine=True)(
        histidine(state), complete_histidine_hydrogens=["HIS1.A"]
    )
    if remove_h:
        protein = Molecule(Chem.RemoveHs(protein))
    assert ring_state(protein[0]) == EXPECTED[state]
    methanol = contact_ligand(protein, "CO")[0]
    formate = contact_ligand(protein, "O=C[O-]")[0]
    nd1 = next(
        a.GetIdx()
        for a in protein[0].GetAtoms()
        if a.GetPDBResidueInfo().GetName().strip() == "ND1"
    )

    def at_nd1(detector: Interaction, ligand: Residue) -> bool:
        return any(
            c["indices"]["protein"] == (nd1,)
            for c in detector.detect(ligand, protein[0])
        )

    assert at_nd1(ImplicitHBDonor(), methanol) is (state == "HIE")
    assert at_nd1(ImplicitHBAcceptor(), methanol) is (state != "HIE")
    assert at_nd1(Anionic(), formate) is (state == "HIP")


def test_serialization_keeps_preserved_parent_chemistry() -> None:
    source = MoleculeStandardizer(preserve_histidine=True)(
        histidine(), complete_histidine_hydrogens=["HIS1.A"]
    )
    handler = RDKitPickleHandler(PROLIF_PICKLE_OPTIONS)
    handler.set()
    try:
        result = pickle.loads(pickle.dumps(source))
    finally:
        handler.reset()
    assert ring_state(result) == EXPECTED["HIE"]
    assert list(ImplicitHBDonor().detect(contact_ligand(result, "CO")[0], result[0]))


@pytest.mark.parametrize("name", ["HE1", "HD2", "HE2"])
def test_named_hydrogen_must_agree_with_connectivity(name: str) -> None:
    source = histidine("HID")
    source.GetAtomWithIdx(11).GetPDBResidueInfo().SetName(name)
    with pytest.raises(ValueError, match="conflicting hydrogen attachment"):
        MoleculeStandardizer(preserve_histidine=True)(
            source, complete_histidine_hydrogens=["HIS1.A"]
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"histidine_states": {"HIS1.A": "HIS"}}, "HID, HIE or HIP"),
        ({"histidine_states": {"HIS9.A": "HID"}}, "Unknown source"),
        ({"complete_histidine_hydrogens": ["HIS9.A"]}, "Unknown source"),
        ({"complete_histidine_hydrogens": ["HIS1.A"]}, "invalid complete"),
    ],
)
def test_invalid_evidence(kwargs: EvidenceInputs, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        MoleculeStandardizer(preserve_histidine=True)(
            histidine(hydrogens=False), **kwargs
        )


def test_no_inferred_inventory_even_if_both_ring_h_are_observed() -> None:
    with pytest.raises(ValueError, match="unresolved"):
        MoleculeStandardizer(preserve_histidine=True)(histidine("HIP"))


def test_wrong_custom_template_cannot_override_resolved_state() -> None:
    standardizer = MoleculeStandardizer(preserve_histidine=True)
    standardizer.engines["HIE"] = standardizer.engines["HIS"]
    with pytest.raises(ValueError, match="did not preserve HIE"):
        standardizer(histidine("HIE", "HIE"))


def test_equivalent_source_resonance_is_allowed() -> None:
    source = histidine("HIP", "HIP", hydrogens=False)
    resonances = list(Chem.ResonanceMolSupplier(source))
    assert len(resonances) >= 2
    for resonance in resonances:
        result = MoleculeStandardizer(preserve_histidine=True)(resonance)
        assert ring_state(result[0]) == EXPECTED["HIP"]


@pytest.mark.parametrize("remove_h", [False, True])
def test_connected_peptide_retains_histidine_and_backbone_chemistry(
    remove_h: bool,
) -> None:
    source = Chem.MolFromSequence("AHG")
    ne2 = next(
        a.GetIdx()
        for a in source.GetAtoms()
        if a.GetPDBResidueInfo().GetName().strip() == "NE2"
    )
    source = Chem.AddHs(source, onlyOnAtoms=[ne2], addResidueInfo=True)
    source.GetAtomWithIdx(source.GetNumAtoms() - 1).GetPDBResidueInfo().SetName(" HE2")
    Compute2DCoords(source)
    before = source.GetConformer().GetPositions().copy()
    result = MoleculeStandardizer(preserve_histidine=True)(
        source, complete_histidine_hydrogens=["HIS2.A"]
    )
    assert len(Chem.GetMolFrags(result)) == 1
    if remove_h:
        result = Molecule(Chem.RemoveHs(result))
    assert ring_state(result["HIS2.A"]) == EXPECTED["HIE"]
    np.testing.assert_array_equal(
        result.GetConformer().GetPositions()[: len(before) - 1], before[:-1]
    )
    for atom in result.GetAtoms():
        assert atom.GetNumRadicalElectrons() == 0
    acceptor_pattern = ImplicitHBAcceptor().lig_pattern
    acceptors = {m[0] for m in result.GetSubstructMatches(acceptor_pattern)}
    for residue in result.residues.values():
        for atom in residue.GetAtoms():
            if (
                atom.GetPDBResidueInfo().GetName().strip() == "N"
                and residue.resid.number > 1
            ):
                assert atom.GetUnsignedProp("mapindex") not in acceptors


@pytest.mark.parametrize(
    "damage", ["ring-chord", "moved-sidechain", "external-ring-bond"]
)
def test_malformed_ring_connectivity_is_not_silently_repaired(damage: str) -> None:
    source = Chem.RWMol(histidine("HIE"))
    if damage == "ring-chord":
        source.AddBond(5, 8, Chem.BondType.SINGLE)
    elif damage == "moved-sidechain":
        source.RemoveBond(4, 5)
        source.AddBond(4, 7, Chem.BondType.SINGLE)
    else:
        metal = Chem.Atom(30)
        metal.SetFormalCharge(2)
        metal.SetMonomerInfo(
            Chem.AtomPDBResidueInfo(
                " ZN ", residueName="ZN", residueNumber=2, chainId="A"
            )
        )
        index = source.AddAtom(metal)
        source.AddBond(6, index, Chem.BondType.DATIVE)
    with pytest.raises(ValueError, match=r"HIS1.A.*ring connectivity"):
        MoleculeStandardizer(preserve_histidine=True)(
            source, complete_histidine_hydrogens=["HIS1.A"]
        )


def test_radical_histidine_is_not_an_ordinary_tautomer() -> None:
    source = histidine("HIE", "HIE")
    carbon = source.GetAtomWithIdx(7)
    carbon.SetNumRadicalElectrons(1)
    carbon.SetNoImplicit(True)
    Chem.SanitizeMol(source)
    with pytest.raises(ValueError, match=r"HIE1.A.*ring.*radical"):
        MoleculeStandardizer(preserve_histidine=True)(source)


def test_template_must_preserve_source_heavy_connectivity(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    standardizer = MoleculeStandardizer(preserve_histidine=True)
    engine = standardizer.engines["HIE"]
    apply = engine.apply

    def rewired(residue: Residue) -> Residue:
        fixed = Chem.RWMol(apply(residue))
        fixed.RemoveBond(0, 1)
        fixed.AddBond(0, 4, Chem.BondType.SINGLE)
        return Residue(fixed)

    monkeypatch.setattr(engine, "apply", rewired)
    with pytest.raises(
        ValueError, match="template changed source heavy-atom connectivity"
    ):
        standardizer(histidine("HIE", "HIE"))
