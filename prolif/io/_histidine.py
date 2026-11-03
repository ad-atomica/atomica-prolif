"""Validate source histidine evidence before selecting a chemistry template."""

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TypeVar

import numpy as np
from rdkit import Chem

from prolif.io.constants import ATOMNAME_ALIASES, STANDARD_RESNAME_MAP
from prolif.molecule import Molecule
from prolif.residue import Residue, ResidueId

HistidineKey = TypeVar("HistidineKey", bound=str | ResidueId)

_COUNTS = {"HID": (1, 0), "HIE": (0, 1), "HIP": (1, 1)}
_TEMPLATES = {"HID": "HID", "HIE": "HIE", "HIP": "HIS"}
_RING = ("CG", "ND1", "CE1", "NE2", "CD2")
_ALIASES = {
    name: "HIP" if template == "HIS" else template
    for name, template in STANDARD_RESNAME_MAP.items()
    if template in {"HIS", "HID", "HIE"} and name != "HIS"
}


def _error(residue: Residue, message: str) -> ValueError:
    return ValueError(f"Histidine {residue.resid}: {message}")


def check_histidine_hydrogens(molecule: Chem.Mol) -> None:
    """Reject disconnected H before residue fragmentation can obscure the cause."""
    for atom in molecule.GetAtoms():
        resid = ResidueId.from_atom(atom)
        if (
            atom.GetAtomicNum() == 1
            and resid.name.upper() in {"HIS", *_ALIASES}
            and (
                atom.GetDegree() != 1
                or atom.GetNeighbors()[0].GetAtomicNum() == 1
                or ResidueId.from_atom(atom.GetNeighbors()[0]) != resid
            )
        ):
            raise ValueError(
                f"Histidine {resid}: invalid hydrogen attachment; "
                "supply bonded source hydrogens within the residue."
            )


def _ring_atoms(
    residue: Residue, parent: Chem.Mol | None = None
) -> dict[str, Chem.Atom]:
    atoms: dict[str, Chem.Atom] = {}
    for atom in residue.GetAtoms():
        name = atom.GetPDBResidueInfo().GetName().strip()
        if name in atoms:
            raise _error(
                residue, f"duplicate atom name {name!r}; provide unique names."
            )
        atoms[name] = atom
    for i, name in enumerate(_RING):
        if name not in atoms or atoms[name].GetAtomicNum() != (
            7 if name in {"ND1", "NE2"} else 6
        ):
            raise _error(residue, f"missing or invalid ring atom {name}.")
        if atoms[name].GetNumRadicalElectrons():
            raise _error(residue, f"unsupported ring radical at {name}.")
        next_name = _RING[(i + 1) % len(_RING)]
        if (
            next_name not in atoms
            or residue.GetBondBetweenAtoms(
                atoms[name].GetIdx(), atoms[next_name].GetIdx()
            )
            is None
        ):
            raise _error(residue, f"missing ring bond {name}-{next_name}.")
    if "CB" not in atoms or atoms["CB"].GetAtomicNum() != 6:
        raise _error(residue, "missing or invalid ring sidechain atom CB.")
    for i, name in enumerate(_RING):
        expected = {_RING[i - 1], _RING[(i + 1) % len(_RING)]}
        if name == "CG":
            expected.add("CB")
        atom = atoms[name]
        if parent is None:
            expected_indices = {atoms[n].GetIdx() for n in expected}
        else:
            expected_indices = {atoms[n].GetUnsignedProp("mapindex") for n in expected}
            atom = parent.GetAtomWithIdx(atom.GetUnsignedProp("mapindex"))
        heavy_indices = {
            a.GetIdx() for a in atom.GetNeighbors() if a.GetAtomicNum() != 1
        }
        if heavy_indices != expected_indices:
            raise _error(residue, f"unsupported ring connectivity at {name}.")
    return atoms


def _inventory(residue: Residue, atoms: dict[str, Chem.Atom]) -> tuple[int, int]:
    """Count graph H, never H inferred from the input's unstandardized valences."""
    counts = {"ND1": 0, "NE2": 0}
    aliases = ATOMNAME_ALIASES["HIS"]
    for atom in residue.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        name = atom.GetPDBResidueInfo().GetName().strip()
        canonical = aliases.get(name, name)
        if (
            atom.GetDegree() != 1
            or atom.GetNeighbors()[0].GetAtomicNum() == 1
            or atom.GetBonds()[0].GetBondType() != Chem.BondType.SINGLE
        ):
            raise _error(
                residue,
                f"invalid hydrogen attachment for {name}; "
                "supply bonded source hydrogens.",
            )
        parent = atom.GetNeighbors()[0]
        attached = parent.GetPDBResidueInfo().GetName().strip()
        expected = {"HD1": "ND1", "HE2": "NE2", "HD2": "CD2", "HE1": "CE1"}.get(
            canonical
        )
        if expected is not None and parent.GetIdx() != atoms[expected].GetIdx():
            raise _error(
                residue,
                f"conflicting hydrogen attachment: {name} is bonded to {attached}, "
                f"expected {expected}.",
            )
        if attached in counts:
            counts[attached] += 1
    if max(counts.values()) > 1:
        raise _error(
            residue, "invalid ring hydrogen inventory: more than one H on a ring N."
        )
    return counts["ND1"], counts["NE2"]


@dataclass
class HistidineEvidence:
    state: str
    source: Residue

    @property
    def template(self) -> str:
        return _TEMPLATES[self.state]

    def verify(self, fixed: Residue) -> None:
        """Compare chemistry and atom correspondence with the pre-template source."""
        source = self.source
        if fixed.GetNumAtoms() != source.GetNumAtoms():
            raise _error(source, "template changed source atom correspondence.")
        heavy_bonds = []
        for mol in (source, fixed):
            heavy_bonds.append(
                {
                    tuple(
                        sorted(
                            (
                                b.GetBeginAtom().GetUnsignedProp("mapindex"),
                                b.GetEndAtom().GetUnsignedProp("mapindex"),
                            )
                        )
                    )
                    for b in mol.GetBonds()
                    if b.GetBeginAtom().GetAtomicNum() != 1
                    and b.GetEndAtom().GetAtomicNum() != 1
                }
            )
        if heavy_bonds[0] != heavy_bonds[1]:
            raise _error(source, "template changed source heavy-atom connectivity.")
        for before, after in zip(source.GetAtoms(), fixed.GetAtoms(), strict=True):
            if (
                before.GetAtomicNum(),
                before.GetIsotope(),
                before.GetUnsignedProp("mapindex"),
                before.GetPDBResidueInfo().GetName(),
                ResidueId.from_atom(before, use_segid=source._use_segid),
            ) != (
                after.GetAtomicNum(),
                after.GetIsotope(),
                after.GetUnsignedProp("mapindex"),
                after.GetPDBResidueInfo().GetName(),
                ResidueId.from_atom(after, use_segid=source._use_segid),
            ):
                raise _error(
                    source,
                    "template changed source atom correspondence or residue identity.",
                )
            if before.GetAtomicNum() == 1 and (
                [a.GetUnsignedProp("mapindex") for a in before.GetNeighbors()]
                != [a.GetUnsignedProp("mapindex") for a in after.GetNeighbors()]
            ):
                raise _error(source, "template changed source hydrogen attachment.")
        if fixed.GetNumConformers() != source.GetNumConformers() or any(
            not np.array_equal(a.GetPositions(), b.GetPositions())
            for a, b in zip(source.GetConformers(), fixed.GetConformers(), strict=True)
        ):
            raise _error(source, "template changed source coordinates.")
        atoms = _ring_atoms(fixed)
        _inventory(fixed, atoms)
        counts = tuple(
            atoms[name].GetTotalNumHs(includeNeighbors=True) for name in ("ND1", "NE2")
        )
        charge = sum(atoms[name].GetFormalCharge() for name in _RING)
        if (
            counts != _COUNTS[self.state]
            or charge != int(self.state == "HIP")
            or not all(atoms[name].GetIsAromatic() for name in _RING)
        ):
            raise _error(
                source,
                f"template did not preserve {self.state}: "
                f"ring N H counts={counts}, ring charge={charge}.",
            )


def resolve_histidines(
    molecule: Molecule,
    states: Mapping[HistidineKey, str] | None,
    complete_hydrogens: Iterable[str | ResidueId],
) -> dict[str, HistidineEvidence]:
    """Resolve all histidines while names, graph H and bonds are still source data."""
    residues = {
        str(r.resid): r
        for r in molecule.residues.values()
        if r.resid.name.upper() in {"HIS", *_ALIASES}
    }
    assignments: dict[str, str] = {}
    for source_key, state in (states or {}).items():
        key = str(source_key)
        if key in assignments:
            raise ValueError(f"Duplicate histidine assignment for {key}.")
        if state not in _COUNTS:
            raise ValueError(
                f"Histidine {key}: state must be HID, HIE or HIP, got {state!r}."
            )
        assignments[key] = state
    complete = {str(key) for key in complete_hydrogens}
    unknown = (assignments.keys() | complete) - residues.keys()
    if unknown:
        raise ValueError(
            f"Unknown source histidine residue keys: {sorted(unknown)}. "
            "Use the original residue identities."
        )
    resolved = {}
    for key, residue in residues.items():
        atoms = _ring_atoms(residue, parent=molecule)
        observed = _inventory(residue, atoms)
        evidence = []
        if alias := _ALIASES.get(residue.resid.name.upper()):
            evidence.append(alias)
        if key in assignments:
            evidence.append(assignments[key])
        if key in complete:
            inventory_state = next(
                (s for s, counts in _COUNTS.items() if counts == observed), None
            )
            if inventory_state is None:
                raise _error(
                    residue,
                    f"invalid complete ring hydrogen inventory {observed}; "
                    "expected HID, HIE or HIP.",
                )
            evidence.append(inventory_state)
        if not evidence:
            raise _error(
                residue,
                "unresolved source assignment; supply histidine_states or declare "
                "complete_histidine_hydrogens after validating source preparation.",
            )
        state = evidence[0]
        if any(s != state for s in evidence) or any(
            n > expected for n, expected in zip(observed, _COUNTS[state], strict=True)
        ):
            raise _error(
                residue,
                f"conflicting source evidence: assignments={evidence}, "
                f"observed ring H={observed}. "
                "Correct the source or its preparation evidence.",
            )
        charge = sum(atoms[name].GetFormalCharge() for name in _RING)
        if charge != 0 and charge != int(state == "HIP"):
            raise _error(
                residue, f"conflicting source ring charge {charge} for {state}."
            )
        resolved[key] = HistidineEvidence(
            state, Residue(residue, use_segid=residue._use_segid)
        )
    return resolved


def read_preserving_pdb(path: str | Path) -> Chem.Mol:
    """Restore extended PDB residue names before RDKit's three-character boundary."""
    molecule = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False)
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError(f"Could not read PDB topology {path}.")
    source_names: dict[int, str] = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
        except ValueError as exc:
            raise ValueError(
                "Preserving PDB residue names requires decimal atom serial numbers; "
                "pass a molecule with intact residue metadata instead."
            ) from exc
        if serial in source_names:
            raise ValueError(
                "Preserving PDB residue names requires unique atom serial numbers "
                "in the first model."
            )
        # MDAnalysis/GROMOS extended names occupy the otherwise blank column 21.
        # RDKit and gemmi both expose only the standard three-character field.
        extended_name = line[17:21].strip()
        source_names[serial] = (
            extended_name if extended_name in _ALIASES else line[17:20].strip()
        )
    for atom in molecule.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info.GetSerialNumber() not in source_names:
            raise ValueError(
                "Could not map PDB source residue names by atom serial number."
            )
        info.SetResidueName(source_names[info.GetSerialNumber()])
    return molecule


def update_parent(parent: Chem.RWMol, fixed: Residue) -> None:
    """Transfer template chemistry without reordering atoms or cutting peptide bonds."""
    indices = {a.GetUnsignedProp("mapindex") for a in fixed.GetAtoms()}
    intra_bonds = {
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
        for index in indices
        for b in parent.GetAtomWithIdx(index).GetBonds()
        if b.GetBeginAtomIdx() in indices and b.GetEndAtomIdx() in indices
    }
    for i, j in intra_bonds:
        parent.RemoveBond(i, j)
    for atom in fixed.GetAtoms():
        # Copy before fragment-specific peptide radical/H adjustments. Implicit
        # valence will be recomputed against the full connected graph.
        parent.ReplaceAtom(atom.GetUnsignedProp("mapindex"), atom)
    for bond in fixed.GetBonds():
        i, j = (
            a.GetUnsignedProp("mapindex")
            for a in (bond.GetBeginAtom(), bond.GetEndAtom())
        )
        parent.AddBond(i, j, bond.GetBondType())
        parent.GetBondBetweenAtoms(i, j).SetIsAromatic(bond.GetIsAromatic())
