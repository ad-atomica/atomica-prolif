# ruff: noqa: T201
# Audit command-line tools intentionally print their results.
"""Audit-only, synthetic probes for atom aliases and extended PDB residue names.

Run from repository root:
Use run_audit.sh in this directory from the repository root.
No library files are changed. Outputs and deterministic PDB fixtures go to
.context/titratable-audit/.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import TypedDict

from rdkit import Chem, rdBase
from rdkit.Chem.rdDistGeom import EmbedMolecule

from prolif.io import MoleculeStandardizer
from prolif.io._histidine import read_preserving_pdb

OUT = Path(".context/titratable-audit")
OUT.mkdir(parents=True, exist_ok=True)


class AliasSnapshot(TypedDict):
    name: str
    site: str
    charge: int
    graph_h: int
    total_h: int
    radicals: int
    neighbors: list[str]
    atom_count: int
    bond_count: int


class AliasResult(TypedDict, total=False):
    parent: AliasSnapshot
    residue: AliasSnapshot
    error: str
    warnings: list[str]


class AtomAliasRow(TypedDict):
    canonical: str
    alias: str
    preserve_histidine: bool
    source: AliasSnapshot
    result: AliasResult


class PDBRow(TypedDict):
    alias: str
    preserve_histidine: bool
    source: AliasSnapshot
    parsed_pdb: AliasSnapshot | None
    direct_result: AliasResult
    pdb_result: AliasResult


class ParentContextRow(TypedDict):
    preserve_histidine: bool
    parent_NZ: AliasSnapshot
    residue_NZ: AliasSnapshot
    residue_context: dict[str, bool]


def atom_name(atom: Chem.Atom) -> str:
    return atom.GetPDBResidueInfo().GetName().strip()


def fixture(
    sequence: str,
    residue_name: str | None = None,
    charges: dict[str, int] | None = None,
) -> Chem.Mol:
    mol = Chem.MolFromSequence(sequence)
    for atom in mol.GetAtoms():
        if charges and atom_name(atom) in charges:
            atom.SetFormalCharge(charges[atom_name(atom)])
        info = atom.GetPDBResidueInfo()
        info.SetSerialNumber(atom.GetIdx() + 1)
        if residue_name:
            info.SetResidueName(residue_name)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addResidueInfo=True)
    assert EmbedMolecule(mol, randomSeed=42) == 0
    return Chem.RemoveHs(mol)


def snapshot(mol: Chem.Mol, site: str) -> AliasSnapshot:
    mol.UpdatePropertyCache(strict=False)
    atom = next(a for a in mol.GetAtoms() if atom_name(a) == site)
    return {
        "name": atom.GetPDBResidueInfo().GetResidueName().strip(),
        "site": site,
        "charge": atom.GetFormalCharge(),
        "graph_h": sum(a.GetAtomicNum() == 1 for a in atom.GetNeighbors()),
        "total_h": atom.GetTotalNumHs(includeNeighbors=True),
        "radicals": atom.GetNumRadicalElectrons(),
        "neighbors": sorted(atom_name(a) for a in atom.GetNeighbors()),
        "atom_count": mol.GetNumAtoms(),
        "bond_count": mol.GetNumBonds(),
    }


def run(source: Chem.Mol | Path, site: str, preserve: bool) -> AliasResult:
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        try:
            result = MoleculeStandardizer(preserve_histidine=preserve)(source)
            answer: AliasResult = {
                "parent": snapshot(result, site),
                "residue": snapshot(result[0], site),
            }
        except Exception as exc:
            answer = {"error": f"{type(exc).__name__}: {exc}"}
    answer["warnings"] = [str(w.message) for w in recorded]
    return answer


def atom_alias_probes() -> list[AtomAliasRow]:
    rows: list[AtomAliasRow] = []
    for canonical, aliases in [
        ("OXT", ("OXT", "OT", "OT2", "O2", "OCT2", "OC2")),
        ("O", ("O", "O1", "OT1", "OCT1", "OC1")),
    ]:
        for alias in aliases:
            mol = fixture("A")
            atom = next(a for a in mol.GetAtoms() if atom_name(a) == canonical)
            atom.GetPDBResidueInfo().SetName(alias)
            for preserve in (False, True):
                rows.append(
                    {
                        "canonical": canonical,
                        "alias": alias,
                        "preserve_histidine": preserve,
                        "source": snapshot(mol, alias),
                        "result": run(Chem.Mol(mol), alias, preserve),
                    }
                )
    return rows


def pdb_boundary_probes() -> list[PDBRow]:
    rows: list[PDBRow] = []
    cases = [
        ("D", "ASPH", "OD2", {"OD2": 0}),
        ("D", "ASPP", "OD2", {"OD2": 0}),
        ("E", "GLUH", "OE2", {"OE2": 0}),
        ("E", "GLUP", "OE2", {"OE2": 0}),
        ("K", "LYSH", "NZ", {"NZ": 1}),
        ("R", "ARGN", "NH2", {"NH2": 0}),
        ("C", "CYS1", "SG", {"SG": -1}),
        ("A", "NALA", "N", {"N": 1}),
        ("A", "CALA", "OXT", {"OXT": -1}),
    ]
    for sequence, name, site, charges in cases:
        mol = fixture(sequence, name, charges)
        # Extended four-letter residue name occupies columns 18-21, as used by
        # the supported histidine extended-name convention in the production reader.
        lines = Chem.MolToPDBBlock(mol).splitlines()
        block = (
            "\n".join(
                line[:17] + name + line[21:]
                if line.startswith(("ATOM  ", "HETATM"))
                else line
                for line in lines
            )
            + "\n"
        )
        path = OUT / f"alias-{name}.pdb"
        path.write_text(block)
        for preserve in (False, True):
            parser = (
                read_preserving_pdb(path)
                if preserve
                else Chem.MolFromPDBFile(str(path), removeHs=False)
            )
            # Record names and atom/bond counts immediately after reading to
            # distinguish parsing loss from standardizer-template behavior.
            rows.append(
                {
                    "alias": name,
                    "preserve_histidine": preserve,
                    "source": snapshot(mol, site),
                    "parsed_pdb": snapshot(parser, site) if parser else None,
                    "direct_result": run(Chem.Mol(mol), site, preserve),
                    "pdb_result": run(path, site, preserve),
                }
            )
    return rows


def parent_context_probe() -> list[ParentContextRow]:
    rows: list[ParentContextRow] = []
    mol = fixture("AKA")
    for preserve in (False, True):
        result = MoleculeStandardizer(preserve_histidine=preserve)(Chem.Mol(mol))
        rows.append(
            {
                "preserve_histidine": preserve,
                "parent_NZ": snapshot(result, "NZ"),
                "residue_NZ": snapshot(result[1], "NZ"),
                "residue_context": {
                    str(r.resid): hasattr(r, "_parent_molecule")
                    for r in result.residues.values()
                },
            }
        )
    return rows


if __name__ == "__main__":
    atom_rows = atom_alias_probes()
    pdb_rows = pdb_boundary_probes()
    result = {
        "rdkit_version": rdBase.rdkitVersion,
        "atom_aliases": atom_rows,
        "extended_residue_names": pdb_rows,
        "parent_context": parent_context_probe(),
    }
    (OUT / "alias-probe-results.json").write_text(json.dumps(result, indent=2) + "\n")
    for row in atom_rows:
        out = row["result"]
        print(
            "atom",
            row["alias"],
            "preserve",
            row["preserve_histidine"],
            "source_neighbors",
            row["source"]["neighbors"],
            "parent",
            out["parent"]["neighbors"] if "parent" in out else None,
            "residue",
            out["residue"]["neighbors"] if "residue" in out else None,
            "error",
            out.get("error"),
        )
    for pdb_row in pdb_rows:
        direct, pdb = pdb_row["direct_result"], pdb_row["pdb_result"]
        print(
            "pdb",
            pdb_row["alias"],
            "preserve",
            pdb_row["preserve_histidine"],
            "parsed",
            pdb_row["parsed_pdb"]["name"] if pdb_row["parsed_pdb"] else None,
            "direct_site_charge",
            direct["residue"]["charge"] if "residue" in direct else None,
            "pdb_site_charge",
            pdb["residue"]["charge"] if "residue" in pdb else None,
            "pdb_error",
            pdb.get("error"),
        )
