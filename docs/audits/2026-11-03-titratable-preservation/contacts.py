# ruff: noqa: T201
# Audit command-line tools intentionally print their results.
"""Geometry-enabled detector probes using the same ligand poses at each boundary."""

import json
import warnings
from collections.abc import Iterator
from itertools import product
from pathlib import Path
from typing import TypedDict

import numpy as np
from audit import REPRESENTATIONS, Snapshot, atom_name, fixture, represent, snapshot
from evidence import JSONValue, validate_json
from numpy.typing import NDArray
from rdkit import Chem, rdBase

from prolif import Molecule
from prolif.interactions import Anionic, Cationic, ImplicitHBAcceptor, ImplicitHBDonor
from prolif.interactions.base import Interaction
from prolif.io import MoleculeStandardizer
from prolif.residue import Residue


class ChangeExample(TypedDict):
    pose: int
    direction: list[float]
    ligand_positions: list[list[float]]
    contacts: dict[str, JSONValue]


class DetectorResult(TypedDict):
    hits: dict[str, list[int]]
    change_example: ChangeExample | dict[str, JSONValue]


class ContactRow(TypedDict):
    family: str
    state: str
    alias: str
    site: str
    preserve_histidine: bool
    representation: str


class ContactResult(ContactRow, total=False):
    error: str
    site_states: dict[str, Snapshot]
    detectors: dict[str, DetectorResult]


def ligand_at(
    point: NDArray[np.float64], direction: NDArray[np.float64], smiles: str
) -> Residue:
    mol = Chem.MolFromSmiles(smiles)
    mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()))
    anchor = mol.GetNumAtoms() - 1
    remote = point + 2.9 * direction
    perp = np.cross(direction, [0.0, 0.0, 1.0])
    if np.linalg.norm(perp) < 0.01:
        perp = np.cross(direction, [0.0, 1.0, 0.0])
    perp /= np.linalg.norm(perp)
    outward = -np.cos(np.radians(109.5)) * direction + np.sin(np.radians(109.5)) * perp
    for i in range(mol.GetNumAtoms()):
        position = remote if i == anchor else remote + (anchor - i) * 1.43 * outward
        mol.GetConformer().SetAtomPosition(i, position)
    return Molecule.from_rdkit(mol)[0]


def unit_directions(count: int = 96) -> Iterator[NDArray[np.float64]]:
    for i in range(count):
        z = 1 - 2 * (i + 0.5) / count
        theta = i * np.pi * (3 - np.sqrt(5))
        yield np.array(
            [np.sqrt(1 - z * z) * np.cos(theta), np.sqrt(1 - z * z) * np.sin(theta), z]
        )


def run() -> None:
    specs = [
        ("LYS", "neutral", "LYN", "NZ"),
        ("ARG", "neutral2", "ARGN", "NH2"),
        ("ARG", "cation1", "ARG", "NH1"),
        ("ASP", "acid2", "ASP", "OD2"),
        ("ASP", "anion2", "ASH", "OD2"),
        ("ASP", "acid1", "ASH", "OD1"),
        ("GLU", "acid2", "GLU", "OE2"),
        ("CYS", "anion", "CYS", "SG"),
        ("CYS", "neutral", "CYM", "SG"),
        ("SS", "disulfide", "CYX", "SG"),
        ("TYR", "anion", "TYR", "OH"),
        ("NTERM", "cation", "ALA", "N"),
        ("CTERM", "anion", "ALA", "OXT"),
    ]
    detectors: list[tuple[str, Interaction, str]] = [
        ("protein_donor", ImplicitHBAcceptor(), "CO"),
        ("protein_acceptor", ImplicitHBDonor(), "CO"),
        ("protein_cation", Anionic(), "O=C[O-]"),
        ("protein_anion", Cationic(), "C[NH3+]"),
    ]
    results: list[ContactResult] = []
    for family, state, alias, site in specs:
        for preserve, representation in product((False, True), REPRESENTATIONS):
            mol, number, _ = fixture(family, state, alias)
            mol = represent(mol, representation)
            source = Molecule(mol)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    standardized = MoleculeStandardizer(preserve_histidine=preserve)(
                        mol
                    )
                except Exception as exc:
                    results.append(
                        {
                            "family": family,
                            "state": state,
                            "alias": alias,
                            "site": site,
                            "preserve_histidine": preserve,
                            "representation": representation,
                            "error": f"{type(exc).__name__}: {exc}",
                        }
                    )
                    continue
            removed = Molecule(Chem.RemoveHs(standardized))
            stages = {
                "source": source,
                "standardized": standardized,
                "removed": removed,
            }
            residues = {
                stage: next(r for r in m.residues.values() if r.resid.number == number)
                for stage, m in stages.items()
            }
            indices = {
                stage: next(a.GetIdx() for a in r.GetAtoms() if atom_name(a) == site)
                for stage, r in residues.items()
            }
            point = residues["source"].GetConformer().GetPositions()[indices["source"]]
            row: ContactResult = {
                "family": family,
                "state": state,
                "alias": alias,
                "site": site,
                "preserve_histidine": preserve,
                "representation": representation,
                "site_states": {
                    s: snapshot(r, number, [site]) for s, r in residues.items()
                },
                "detectors": {},
            }
            for name, detector, smiles in detectors:
                hits: dict[str, list[int]] = {stage: [] for stage in stages}
                examples: ChangeExample | dict[str, JSONValue] = {}
                for pose, direction in enumerate(unit_directions()):
                    ligand = ligand_at(point, direction, smiles)
                    pose_hits: dict[str, JSONValue] = {}
                    for stage, res in residues.items():
                        contacts = [
                            c
                            for c in detector.detect(ligand, res)
                            if c["indices"]["protein"] == (indices[stage],)
                            and c["indices"]["ligand"] == (ligand.GetNumAtoms() - 1,)
                        ]
                        if contacts:
                            hits[stage].append(pose)
                            pose_hits[stage] = validate_json(
                                json.loads(json.dumps(contacts[0]))
                            )
                    if len(pose_hits) not in {0, len(stages)} and not examples:
                        examples = {
                            "pose": pose,
                            "direction": direction.tolist(),
                            "ligand_positions": ligand.GetConformer()
                            .GetPositions()
                            .tolist(),
                            "contacts": pose_hits,
                        }
                row["detectors"][name] = {"hits": hits, "change_example": examples}
            results.append(row)
    Path(".context/titratable-audit/contact-results.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    for r in results:
        if r["preserve_histidine"] and r["representation"] == "no_graph_h":
            print(
                r["family"],
                r["state"],
                r["alias"],
                r["site"],
                {
                    k: {s: len(v) for s, v in d["hits"].items()}
                    for k, d in r["detectors"].items()
                },
            )


if __name__ == "__main__":
    rdBase.DisableLog("rdApp.warning")
    run()
