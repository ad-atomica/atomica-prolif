"""Preserving-mode fragments must reflect the sanitized connected graph."""

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

from prolif import Molecule
from prolif.io import MoleculeStandardizer


def resonant_arginine() -> Chem.Mol:
    mol = Chem.MolFromSequence("GRG")
    atoms = {
        a.GetPDBResidueInfo().GetName().strip(): a
        for a in mol.GetAtoms()
        if a.GetPDBResidueInfo().GetResidueNumber() == 2
    }
    for name in ("NH1", "NH2"):
        atoms[name].SetFormalCharge(int(name == "NH1"))
        mol.GetBondBetweenAtoms(atoms["CZ"].GetIdx(), atoms[name].GetIdx()).SetBondType(
            Chem.BondType.DOUBLE if name == "NH1" else Chem.BondType.SINGLE
        )
    Chem.SanitizeMol(mol)
    Compute2DCoords(mol)
    return mol


@pytest.mark.parametrize("graph_h", [False, True])
def test_preserving_residues_use_final_parent_chemistry(graph_h: bool) -> None:
    source = resonant_arginine()
    if graph_h:
        source = Chem.AddHs(source, addCoords=True, addResidueInfo=True)
        for atom in source.GetAtoms():
            if atom.GetAtomicNum() == 1:
                atom.GetPDBResidueInfo().SetName(f"Q{atom.GetIdx():03d}")
    result = MoleculeStandardizer(preserve_histidine=True)(source)
    rebuilt = Molecule(Chem.RemoveHs(result))
    for residue in result.residues.values():
        for atom in residue.GetAtoms():
            parent = result.GetAtomWithIdx(atom.GetUnsignedProp("mapindex"))
            assert atom.GetFormalCharge() == parent.GetFormalCharge()
            assert atom.GetHybridization() == parent.GetHybridization()
            assert atom.GetTotalNumHs(includeNeighbors=True) == parent.GetTotalNumHs(
                includeNeighbors=True
            )
    before = next(
        a
        for a in result[1].GetAtoms()
        if a.GetPDBResidueInfo().GetName().strip() == "NH1"
    )
    after = next(
        a
        for a in rebuilt[1].GetAtoms()
        if a.GetPDBResidueInfo().GetName().strip() == "NH1"
    )
    assert before.GetHybridization() == after.GetHybridization()
    np.testing.assert_array_equal(
        result.GetConformer().GetPositions(), source.GetConformer().GetPositions()
    )
