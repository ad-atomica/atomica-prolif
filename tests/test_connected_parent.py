"""Fragmented residues retain connected chemistry without sharing mutable owners."""

import gc
import pickle
import weakref

import numpy as np
import pytest
from prolif import Molecule
from prolif.interactions import ImplicitHBAcceptor, ImplicitHBDonor
from prolif.pickling import PROLIF_PICKLE_OPTIONS, RDKitPickleHandler
from prolif.residue import Residue
from rdkit import Chem


def molecule_with_coordinates(smiles: str) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(smiles)
    molecule.AddConformer(Chem.Conformer(molecule.GetNumAtoms()))
    return molecule


@pytest.fixture
def peptide() -> Molecule:
    molecule = molecule_with_coordinates("CC(=O)NCC(=O)O")
    # The amide N's peptide plane differs from its fragmented N-CA-C plane.
    positions = [
        (-2.4, -0.6, 0),
        (-1.2, 0, 0),
        (-1.2, 1.2, 0),
        (0, 0, 0),
        (0.6, 1.0392304845, 0),
        (1.1, 1.6, 1.2),
        (1.4, 2.8, 1.2),
        (1.8, 0.8, 1.7),
    ]
    for atom, position in zip(molecule.GetAtoms(), positions, strict=True):
        atom.SetMonomerInfo(
            Chem.AtomPDBResidueInfo(
                f" {atom.GetSymbol():<3}",
                residueName="ACE" if atom.GetIdx() < 3 else "GLY",
                residueNumber=1 if atom.GetIdx() < 3 else 2,
                chainId="A",
            )
        )
        molecule.GetConformer().SetAtomPosition(atom.GetIdx(), position)
    return Molecule(molecule)


def oxygen_ligand(peptide: Molecule, elevation: float = 0) -> Residue:
    coordinates = peptide.GetConformer().GetPositions()
    origin = coordinates[3]
    bonds = coordinates[[1, 4]] - origin
    bonds /= np.linalg.norm(bonds, axis=1)[:, None]
    direction = -bonds.sum(axis=0)
    direction /= np.linalg.norm(direction)
    normal = np.cross(*bonds)
    normal /= np.linalg.norm(normal)
    direction = (
        np.cos(np.radians(elevation)) * direction
        + np.sin(np.radians(elevation)) * normal
    )
    ligand = molecule_with_coordinates("CC=O")
    perpendicular = np.cross(direction, bonds[0])
    perpendicular /= np.linalg.norm(perpendicular)
    for index, position in (
        (2, origin + 2.9 * direction),
        (1, origin + 4.13 * direction),
        (0, origin + 4.88 * direction + 1.299 * perpendicular),
    ):
        ligand.GetConformer().SetAtomPosition(index, position)
    return next(iter(Molecule.from_rdkit(ligand).residues.values()))


@pytest.mark.parametrize(("elevation", "expected"), [(0, True), (60, False)])
def test_backbone_donor_uses_connected_peptide_plane(
    peptide: Molecule, elevation: float, expected: bool
) -> None:
    before = peptide.ToBinary()
    contacts = [
        contact
        for contact in ImplicitHBAcceptor().detect(
            oxygen_ligand(peptide, elevation), peptide.residues["GLY2.A"]
        )
        if contact["parent_indices"]["protein"] == (3,)
    ]
    assert bool(contacts) is expected
    if expected:
        assert contacts[0]["donor_plane_angle"] == pytest.approx(0, abs=1e-6)
        assert contacts[0]["ideal_donor_angle"] == 120
        assert contacts[0]["indices"]["protein"] == (0,)
    assert peptide.ToBinary() == before


@pytest.mark.parametrize("mode", ["connected", "standalone", "custom-pattern"])
def test_backbone_amide_is_not_an_acceptor(peptide: Molecule, mode: str) -> None:
    ligand = molecule_with_coordinates("CO")
    ligand.GetConformer().SetAtomPosition(1, (0, 0, 2.9))
    ligand.GetConformer().SetAtomPosition(0, (0, 0, 4.33))
    ligand_residue = next(iter(Molecule.from_rdkit(ligand).residues.values()))
    protein_residue = peptide.residues["GLY2.A"]
    if mode == "standalone":
        protein_residue = Residue(protein_residue)
    detector = (
        ImplicitHBDonor(acceptor="[N]", ignore_geometry_checks=True)
        if mode == "custom-pattern"
        else ImplicitHBDonor(ignore_geometry_checks=True)
    )
    contacts = [
        contact
        for contact in detector.detect(ligand_residue, protein_residue)
        if contact["parent_indices"]["protein"] == (3,)
    ]
    assert bool(contacts) is (mode != "connected")


@pytest.mark.parametrize("copy_mode", ["pickle", "supplied-residues"])
def test_parent_geometry_survives_copy_and_serialization(
    peptide: Molecule, copy_mode: str
) -> None:
    handler = RDKitPickleHandler(PROLIF_PICKLE_OPTIONS)
    handler.set()
    try:
        receptor = pickle.loads(pickle.dumps(peptide))
    finally:
        handler.reset()
    if copy_mode == "supplied-residues":
        shifted = Chem.Mol(receptor)
        for index, position in enumerate(shifted.GetConformer().GetPositions() + 100):
            shifted.GetConformer().SetAtomPosition(index, position)
        Molecule(shifted, residues=list(receptor.residues.values()))
    contacts = [
        contact
        for contact in ImplicitHBAcceptor().detect(
            oxygen_ligand(receptor), receptor.residues["GLY2.A"]
        )
        if contact["parent_indices"]["protein"] == (3,)
    ]
    assert len(contacts) == 1
    assert contacts[0]["donor_plane_angle"] == pytest.approx(0, abs=1e-3)


def test_parent_does_not_retain_discarded_molecule_wrappers(peptide: Molecule) -> None:
    receptor = Molecule(Chem.Mol(peptide))
    receptor_ref = weakref.ref(receptor)
    residue_ref = weakref.ref(receptor.residues["GLY2.A"])
    del receptor
    gc.collect()
    assert receptor_ref() is None
    assert residue_ref() is None


@pytest.mark.parametrize(
    ("smiles", "expected_count"),
    [
        ("OCC(O)" * 501 + "CN(C)C", 1),
        ("NCC" * 1001 + "N(C)C", 1),
        ("NCC(=O)" * 1002 + "O", 0),
    ],
    ids=["total-matches", "recursive-matches", "amide-exclusions"],
)
def test_connected_matching_is_not_truncated(smiles: str, expected_count: int) -> None:
    receptor = molecule_with_coordinates(smiles)
    nitrogens = [a.GetIdx() for a in receptor.GetAtoms() if a.GetAtomicNum() == 7]
    parent_index = nitrogens[-1]
    for atom in receptor.GetAtoms():
        atom.SetMonomerInfo(
            Chem.AtomPDBResidueInfo(
                atom.GetSymbol(),
                residueName="LIG",
                residueNumber=1 + (atom.GetIdx() >= parent_index),
            )
        )
    ligand = molecule_with_coordinates("CO")
    ligand.GetConformer().SetAtomPosition(1, (2.9, 0, 0))
    ligand.GetConformer().SetAtomPosition(0, (4.33, 0, 0))
    ligand_residue = next(iter(Molecule.from_rdkit(ligand).residues.values()))
    contacts = [
        contact
        for contact in ImplicitHBDonor(ignore_geometry_checks=True).detect(
            ligand_residue, Molecule(receptor).residues["LIG2"]
        )
        if contact["parent_indices"]["protein"] == (parent_index,)
    ]
    assert len(contacts) == expected_count
    if contacts:
        assert contacts[0]["indices"]["protein"] == (0,)
