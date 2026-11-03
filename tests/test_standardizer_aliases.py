"""Atom and extended PDB aliases must reach the same chemistry templates."""

from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

from prolif.io import MoleculeStandardizer
from prolif.io._histidine import read_preserving_pdb


def amino_acid(sequence: str = "A", name: str | None = None) -> Chem.Mol:
    mol = Chem.MolFromSequence(sequence)
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        info.SetSerialNumber(atom.GetIdx() + 1)
        if name is not None:
            info.SetResidueName(name)
    Compute2DCoords(mol)
    return mol


def atom_names(mol: Chem.Mol) -> list[str]:
    return [atom.GetPDBResidueInfo().GetName().strip() for atom in mol.GetAtoms()]


def bond_orders(mol: Chem.Mol) -> dict[tuple[int, int], float]:
    return {
        tuple(
            sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        ): bond.GetBondTypeAsDouble()
        for bond in mol.GetBonds()
    }


def extended_pdb(mol: Chem.Mol, name: str) -> str:
    return "".join(
        line[:17] + f"{name:<4}" + line[21:]
        if line.startswith(("ATOM  ", "HETATM"))
        else line
        for line in Chem.MolToPDBBlock(mol).splitlines(keepends=True)
    )


@pytest.mark.parametrize("preserve_histidine", [False, True])
@pytest.mark.parametrize(
    ("canonical", "alias", "order"),
    [
        ("O", "O", 2.0),
        ("O", "O1", 2.0),
        ("O", "OT1", 2.0),
        ("O", "OCT1", 2.0),
        ("O", "OC1", 2.0),
        ("OXT", "OXT", 1.0),
        ("OXT", "O2", 1.0),
        ("OXT", "OT2", 1.0),
        ("OXT", "OCT2", 1.0),
        ("OXT", "OC2", 1.0),
        ("OXT", "OT", 1.0),
    ],
)
def test_protein_oxygen_alias_retains_carbonyl_connectivity(
    canonical: str, alias: str, order: float, preserve_histidine: bool
) -> None:
    source = amino_acid()
    oxygen = source.GetAtomWithIdx(atom_names(source).index(canonical))
    oxygen.GetPDBResidueInfo().SetName(alias)
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    result = MoleculeStandardizer(preserve_histidine=preserve_histidine)(source)

    for mol in (result, result[0]):
        names = atom_names(mol)
        atom = mol.GetAtomWithIdx(names.index(alias))
        neighbors = [
            a.GetPDBResidueInfo().GetName().strip() for a in atom.GetNeighbors()
        ]
        assert neighbors == ["C"]
        assert atom.GetBonds()[0].GetBondTypeAsDouble() == order
        assert mol.GetNumBonds() == 5
        assert names == atom_names(source)
        np.testing.assert_array_equal(
            mol.GetConformer().GetPositions(), source.GetConformer().GetPositions()
        )
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


def test_common_aliases_merge_with_residue_specific_aliases() -> None:
    source = amino_acid("S")
    expected_bonds = bond_orders(source)
    for old, new in (("OG", "OG1"), ("O", "O1")):
        source.GetAtomWithIdx(
            atom_names(source).index(old)
        ).GetPDBResidueInfo().SetName(new)
    result = MoleculeStandardizer(preserve_histidine=True)(source)
    assert bond_orders(result) == expected_bonds
    assert bond_orders(result[0]) == expected_bonds
    assert atom_names(result) == atom_names(source)


@pytest.mark.parametrize(("first", "second"), [("O", "O1"), ("O", "O"), ("O2", "OT2")])
@pytest.mark.filterwarnings("ignore:Residue .* has a different number of heavy atoms")
def test_duplicate_canonical_names_fail_without_mutating_source(
    first: str, second: str
) -> None:
    source = amino_acid()
    source.GetAtomWithIdx(3).GetPDBResidueInfo().SetName(first)
    source.GetAtomWithIdx(5).GetPDBResidueInfo().SetName(second)
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match="duplicate canonical atom name"):
        MoleculeStandardizer(preserve_histidine=True)(source)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


@pytest.mark.parametrize(
    ("sequence", "name", "normalized", "site", "charge"),
    [
        ("D", "ASPH", "ASH", "OD2", 0),
        ("D", "ASPP", "ASH", "OD2", 0),
        ("E", "GLUH", "GLH", "OE2", 0),
        ("E", "GLUP", "GLH", "OE2", 0),
        ("C", "CYS1", "CYX", "SG", -1),
        ("A", "NALA", "ALA", "N", 0),
        ("A", "CALA", "ALA", "OXT", 0),
        ("K", "LYSH", "LYS", "NZ", 1),
        ("R", "ARGN", "ARG", "NH2", 1),
    ],
)
def test_extended_residue_alias_matches_in_memory_standardization(
    tmp_path: Path,
    sequence: str,
    name: str,
    normalized: str,
    site: str,
    charge: int,
) -> None:
    source = amino_acid(sequence, name)
    before = source.ToBinary(Chem.PropertyPickleOptions.AllProps)
    path = tmp_path / f"{name}.pdb"
    block = extended_pdb(source, name)
    path.write_text(block)

    parsed = read_preserving_pdb(path)
    assert {a.GetPDBResidueInfo().GetResidueName() for a in parsed.GetAtoms()} == {name}
    assert atom_names(parsed) == atom_names(source)
    np.testing.assert_allclose(
        parsed.GetConformer().GetPositions(),
        source.GetConformer().GetPositions(),
        atol=0.0005,
    )

    standardizer = MoleculeStandardizer(preserve_histidine=True)
    direct, from_pdb = standardizer(source), standardizer(path)
    for result in (direct, from_pdb):
        assert result[0].resid.name == normalized
        for mol in (result, result[0]):
            atom = mol.GetAtomWithIdx(atom_names(mol).index(site))
            assert atom.GetFormalCharge() == charge
            assert atom_names(mol) == atom_names(source)
            assert bond_orders(mol) == bond_orders(direct)
    assert source.ToBinary(Chem.PropertyPickleOptions.AllProps) == before
    assert path.read_text() == block


def test_unknown_extended_residue_name_is_not_guessed(tmp_path: Path) -> None:
    path = tmp_path / "unknown.pdb"
    path.write_text(extended_pdb(amino_acid(), "ALAX"))
    with pytest.raises(
        ValueError, match=r"unsupported extended PDB residue name.*ALAX"
    ):
        read_preserving_pdb(path)


def test_standard_pdb_residue_name_is_retained(tmp_path: Path) -> None:
    path = tmp_path / "tyrosinate.pdb"
    path.write_text(extended_pdb(amino_acid("Y", "TYM"), "TYM"))
    parsed = read_preserving_pdb(path)
    names = {a.GetPDBResidueInfo().GetResidueName() for a in parsed.GetAtoms()}
    assert names == {"TYM"}


@pytest.mark.parametrize(
    ("code", "name", "expected"), [("D", "ASPH", "ASH"), ("A", "NALA", "ALA")]
)
def test_legacy_reader_retains_extended_template_selection(
    tmp_path: Path, code: str, name: str, expected: str
) -> None:
    source = amino_acid(code, name)
    path = tmp_path / "source.pdb"
    path.write_text(extended_pdb(source, name))
    result = MoleculeStandardizer()(path)
    assert result[0].resid.name == expected
    assert bond_orders(result[0]) == bond_orders(MoleculeStandardizer()(source)[0])
