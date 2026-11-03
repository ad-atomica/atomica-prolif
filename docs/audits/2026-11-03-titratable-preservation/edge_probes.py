# ruff: noqa: T201
# Audit command-line tools intentionally print their results.
"""Additional controls for source naming, source mutation and ambiguous PDBs."""

import json
import warnings
from pathlib import Path
from typing import TypedDict

from audit import Snapshot, atoms_at, chemistry, fixture, heavy_info, snapshot
from rdkit import Chem, rdBase
from rdkit.Chem.rdDepictor import Compute2DCoords

from prolif import Molecule
from prolif.io import MoleculeStandardizer
from prolif.io._histidine import read_preserving_pdb

OUT = Path(".context/titratable-audit")


class Evaluation(TypedDict, total=False):
    parent: Snapshot
    residue: Snapshot
    error: str
    warnings: list[str]


EdgeValue = str | bool | tuple[str, str, str] | Snapshot | Evaluation


def evaluate(
    source: Chem.Mol | Path, number: int, sites: list[str], preserve: bool
) -> Evaluation:
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        try:
            result = MoleculeStandardizer(preserve_histidine=preserve)(source)
            res = next(r for r in result.residues.values() if r.resid.number == number)
            answer: Evaluation = {
                "parent": snapshot(result, number, sites),
                "residue": snapshot(res, number, sites),
            }
        except Exception as exc:
            answer = {"error": f"{type(exc).__name__}: {exc}"}
        answer["warnings"] = [str(w.message) for w in ws]
        return answer


def run() -> None:  # noqa: PLR0912
    records: list[dict[str, EdgeValue]] = []
    # Explicit conventional H names must reproduce the same state failures as
    # the unique arbitrary Q names in the primary matrix.
    specs = [
        ("LYS", "neutral", "LYN"),
        ("ARG", "neutral2", "ARGN"),
        ("ASP", "acid2", "ASP"),
        ("GLU", "acid2", "GLU"),
        ("CYS", "anion", "CYS"),
        ("CYS", "neutral", "CYM"),
    ]
    hydrogen_names = {
        "NZ": ["HZ1", "HZ2", "HZ3"],
        "NE": ["HE"],
        "NH1": ["HH11", "HH12"],
        "NH2": ["HH21", "HH22"],
        "OD2": ["HD2"],
        "OE2": ["HE2"],
        "SG": ["HG"],
    }
    for spec in specs:
        mol, number, sites = fixture(*spec)
        named = Chem.Mol(mol)
        for n in sites:
            a = atoms_at(named, number)[n]
            hs = [h for h in a.GetNeighbors() if h.GetAtomicNum() == 1]
            for h, name in zip(hs, hydrogen_names.get(n, []), strict=False):
                h.GetPDBResidueInfo().SetName(name)
        for preserve in (False, True):
            records.append(
                {
                    "probe": "canonical_h_names",
                    "spec": spec,
                    "preserve_histidine": preserve,
                    "arbitrary_names": evaluate(mol, number, sites, preserve),
                    "canonical_names": evaluate(named, number, sites, preserve),
                }
            )
    # True intact source ProLIF input; preserve=True must return without source
    # mutation, while legacy is documented to replace residues in place.
    for preserve in (False, True):
        mol, number, sites = fixture("LYS", "neutral", "LYN")
        source = Molecule(Chem.RemoveHs(mol))
        before_parent = chemistry(source)
        before_identity = heavy_info(source)
        before_res = snapshot(source[1], number, sites)
        result = MoleculeStandardizer(preserve_histidine=preserve)(source)
        records.append(
            {
                "probe": "mutation",
                "preserve_histidine": preserve,
                "same_object": source is result,
                "source_parent_chemistry_unchanged": before_parent == chemistry(source),
                "source_identity_unchanged": before_identity == heavy_info(source),
                "source_residue_before": before_res,
                "source_residue_after": snapshot(source[1], number, sites),
                "result_parent": snapshot(result, number, sites),
                "result_residue": snapshot(result[1], number, sites),
            }
        )
    # Discard charge fields from H-free PDBs, retaining CONECT records. This is
    # deliberately ambiguous evidence, not a source-assigned molecule anymore.
    for spec in [
        ("TYR", "anion", "TYR"),
        ("NTERM", "cation", "ALA"),
        ("CTERM", "anion", "ALA"),
        ("LYS", "neutral", "LYS"),
        ("ASP", "acid2", "ASP"),
    ]:
        mol, number, sites = fixture(*spec)
        mol = Chem.RemoveHs(mol)
        block = (
            "\n".join(
                line[:78] + "  " + line[80:]
                if line.startswith(("ATOM  ", "HETATM"))
                else line
                for line in Chem.MolToPDBBlock(mol).splitlines()
            )
            + "\n"
        )
        path = OUT / f"ambiguous-{spec[0]}.pdb"
        path.write_text(block)
        for preserve in (False, True):
            parsed = (
                read_preserving_pdb(path)
                if preserve
                else Chem.MolFromPDBFile(str(path), removeHs=False)
            )
            parsed.UpdatePropertyCache(strict=False)
            records.append(
                {
                    "probe": "ambiguous_pdb",
                    "spec": spec,
                    "preserve_histidine": preserve,
                    "original_source": snapshot(mol, number, sites),
                    "parsed": snapshot(parsed, number, sites),
                    "result": evaluate(path, number, sites, preserve),
                }
            )
    # A different residue name changes the guessed force field for the entire
    # receptor, hence the CYS template, despite identical CYS source evidence.
    for preserve in (False, True):
        mol, number, sites = fixture("CYS", "neutral", "CYS")
        mol = Chem.RemoveHs(mol)
        other = Chem.MolFromSequence("A")
        Compute2DCoords(other)
        for a in other.GetAtoms():
            a.GetPDBResidueInfo().SetResidueNumber(10)
            a.GetPDBResidueInfo().SetResidueName("DALA")
        source = Chem.CombineMols(mol, other)
        records.append(
            {
                "probe": "forcefield_context",
                "preserve_histidine": preserve,
                "CYS_alone": evaluate(mol, number, sites, preserve),
                "with_DALA": evaluate(source, number, sites, preserve),
            }
        )
    # A named proton on the other acid oxygen tests whether the template trusts
    # name over actual source attachment. Coordinates stay unchanged.
    for preserve in (False, True):
        mol, number, sites = fixture("ASP", "acid1", "ASH")
        proton = next(
            h
            for h in atoms_at(mol, number)["OD1"].GetNeighbors()
            if h.GetAtomicNum() == 1
        )
        proton.GetPDBResidueInfo().SetName("HD2")
        records.append(
            {
                "probe": "contradictory_h_name_attachment",
                "preserve_histidine": preserve,
                "source": snapshot(mol, number, sites),
                "result": evaluate(mol, number, sites, preserve),
            }
        )
    OUT.joinpath("edge-results.json").write_text(json.dumps(records, indent=2) + "\n")
    print("Wrote", len(records), "edge probes")
    for r in records:
        if r["probe"] == "canonical_h_names":
            print(
                r["probe"],
                r["spec"],
                r["preserve_histidine"],
                r["arbitrary_names"] == r["canonical_names"],
            )
        elif r["probe"] == "mutation":
            print(
                r["probe"],
                r["preserve_histidine"],
                r["same_object"],
                r["source_identity_unchanged"],
            )
        elif r["probe"] == "ambiguous_pdb":
            print(
                r["probe"],
                r["spec"],
                r["preserve_histidine"],
                isinstance(r["result"], dict) and "error" in r["result"],
            )
        else:
            print(r)


if __name__ == "__main__":
    rdBase.DisableLog("rdApp.warning")
    rdBase.DisableLog("rdApp.error")
    run()
