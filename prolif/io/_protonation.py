"""Validate and retain deliberately prepared, standard protein chemistry.

This path is a validator, not a protonation or topology inference engine. Source
bonding, charge and H inventories remain authoritative after explicit preparation
provenance is supplied through MoleculeStandardizer.
"""

from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from rdkit import Chem

from prolif.io._histidine import _ALIASES as HISTIDINE_ALIASES
from prolif.io.constants import ATOMNAME_ALIASES, STANDARD_AA, STANDARD_RESNAME_MAP
from prolif.molecule import Molecule
from prolif.residue import ResidueId

_FAMILIES = frozenset(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
)
_BASE = {
    "ASH": "ASP",
    "GLH": "GLU",
    "CYX": "CYS",
    "HID": "HIS",
    "HIE": "HIS",
    "TYM": "TYR",
}
_ORDER = {
    "SING": Chem.BondType.SINGLE,
    "DOUB": Chem.BondType.DOUBLE,
    "TRIP": Chem.BondType.TRIPLE,
}


@dataclass
class _Template:
    elements: dict[str, int]
    bonds: dict[frozenset[str], Chem.BondType]
    aromatic: set[str]
    hydrogen_parents: dict[str, str]
    aliases: dict[str, str]


@lru_cache(maxsize=20)
def _template(family: str) -> _Template:
    block = STANDARD_AA.find_block(family)
    assert block is not None
    elements = {
        row[0]: Chem.GetPeriodicTable().GetAtomicNumber(row[1])
        for row in block.find("_chem_comp_atom.", ["atom_id", "type_symbol"])
    }
    bonds = {}
    hydrogens = {}
    for row in block.find(
        "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]
    ):
        first, second, order = row
        # The ARG CIF includes the terminal N-H3 bond but omits its leaving-H
        # atom row. It must never become an expected heavy-atom bond.
        if family == "ARG" and "H3" in {first, second}:
            elements["H3"] = 1
        if elements.get(first) == 1:
            hydrogens[first] = second
        elif elements.get(second) == 1:
            hydrogens[second] = first
        else:
            bonds[frozenset((first, second))] = _ORDER[order]
    # A single CIF state does not enumerate every supported proton site. These
    # names still assert attachment even when the reference state lacks that H.
    hydrogens.update({"H": "N", "H2": "N", "H3": "N"})
    if family == "ASP":
        hydrogens.update({"HD1": "OD1", "HD2": "OD2"})
    elif family == "GLU":
        hydrogens.update({"HE1": "OE1", "HE2": "OE2"})
    aromatic = {
        row[0]
        for row in block.find("_chem_comp_atom.", ["atom_id", "pdbx_aromatic_flag"])
        if row[1] == "Y"
    }
    aliases = {**ATOMNAME_ALIASES["Protein"], **ATOMNAME_ALIASES.get(family, {})}
    return _Template(elements, bonds, aromatic, hydrogens, aliases)


@dataclass
class _Residue:
    resid: ResidueId
    family: str
    atoms: dict[str, Chem.Atom]
    template: _Template

    def fail(self, detail: str) -> ValueError:
        return ValueError(
            f"Residue {self.resid}: {detail} Supply a complete prepared protein graph."
        )


def _qh(atom: Chem.Atom) -> tuple[int, int]:
    return atom.GetFormalCharge(), atom.GetTotalNumHs(includeNeighbors=True)


def _canonical(name: str, template: _Template) -> str:
    return name if name in template.elements else template.aliases.get(name, name)


def _residues(
    mol: Chem.Mol, use_segid: bool
) -> tuple[list[_Residue], dict[int, tuple[_Residue, str]]]:
    grouped: dict[ResidueId, list[Chem.Atom]] = defaultdict(list)
    for atom in mol.GetAtoms():
        if atom.GetPDBResidueInfo() is None:
            raise ValueError(
                "Prepared protein atoms require source residue and atom names."
            )
        grouped[ResidueId.from_atom(atom, use_segid=use_segid)].append(atom)
    residues = []
    index = {}
    for resid, atoms in grouped.items():
        normalized = STANDARD_RESNAME_MAP.get(resid.name.upper(), resid.name.upper())
        family = _BASE.get(normalized, normalized)
        if family not in _FAMILIES:
            raise ValueError(
                f"Residue {resid}: unsupported residue for strict protonation "
                "preservation; select complete standard protein residues."
            )
        template = _template(family)
        residue = _Residue(resid, family, {}, template)
        for atom in atoms:
            source_name = atom.GetPDBResidueInfo().GetName().strip()
            name = _canonical(source_name, template)
            if not source_name or name in residue.atoms:
                raise residue.fail(f"empty or duplicate canonical atom name {name!r}.")
            residue.atoms[name] = atom
            index[atom.GetIdx()] = residue, name
        expected = {n: z for n, z in template.elements.items() if z != 1 and n != "OXT"}
        observed = {
            n: a.GetAtomicNum()
            for n, a in residue.atoms.items()
            if a.GetAtomicNum() != 1 and n != "OXT"
        }
        if expected != observed:
            raise residue.fail("missing, extra or incorrectly named heavy atoms.")
        if "OXT" in residue.atoms and residue.atoms["OXT"].GetAtomicNum() != 8:
            raise residue.fail("invalid terminal oxygen OXT.")
        residues.append(residue)
    return residues, index


def _hydrogens(residue: _Residue, index: dict[int, tuple[_Residue, str]]) -> None:
    for name, atom in residue.atoms.items():
        if atom.GetAtomicNum() != 1:
            continue
        if (
            atom.GetFormalCharge()
            or atom.GetDegree() != 1
            or atom.GetBonds()[0].GetBondType() != Chem.BondType.SINGLE
        ):
            raise residue.fail(f"invalid hydrogen attachment at {name}.")
        parent_res, parent_name = index[atom.GetNeighbors()[0].GetIdx()]
        if parent_res is not residue or atom.GetNeighbors()[0].GetAtomicNum() == 1:
            raise residue.fail(f"invalid hydrogen attachment at {name}.")
        expected_parent = residue.template.hydrogen_parents.get(name)
        if expected_parent is not None and parent_name != expected_parent:
            raise residue.fail(
                f"conflicting hydrogen name {name}: attached to {parent_name}, "
                f"expected {expected_parent}."
            )


def _topology(
    residue: _Residue, index: dict[int, tuple[_Residue, str]]
) -> dict[str, list[tuple[_Residue, str]]]:
    expected = {
        pair
        for pair in residue.template.bonds
        if all(name in residue.atoms for name in pair)
    }
    observed = set()
    external: dict[str, list[tuple[_Residue, str]]] = defaultdict(list)
    for name, atom in residue.atoms.items():
        if atom.GetNumRadicalElectrons():
            raise residue.fail(f"unsupported radical at {name}.")
        if atom.GetAtomicNum() == 1:
            continue
        if atom.GetIsAromatic() != (name in residue.template.aromatic):
            raise residue.fail(f"unsupported aromaticity at {name}.")
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if neighbor.GetAtomicNum() == 1:
                continue
            other, other_name = index[neighbor.GetIdx()]
            if other is residue:
                observed.add(frozenset((name, other_name)))
            else:
                allowed = (name, other_name) in {("C", "N"), ("N", "C")} or (
                    name == other_name == "SG"
                    and residue.family == other.family == "CYS"
                )
                if not allowed or bond.GetBondType() != Chem.BondType.SINGLE:
                    raise residue.fail(
                        f"unsupported inter-residue connectivity at {name} "
                        f"to {other.resid}:{other_name}."
                    )
                external[name].append((other, other_name))
                if len(external[name]) != 1:
                    raise residue.fail(f"multiple inter-residue bonds at {name}.")
    if expected != observed:
        raise residue.fail(
            "unsupported heavy-atom connectivity; "
            "source bonds do not match the residue."
        )
    return external


def _acid(
    mol: Chem.Mol, res: _Residue, carbon: str, oxygens: tuple[str, str]
) -> set[str]:
    states = []
    orders = []
    for name in oxygens:
        atom = res.atoms[name]
        order = mol.GetBondBetweenAtoms(
            res.atoms[carbon].GetIdx(), atom.GetIdx()
        ).GetBondType()
        state = _qh(atom)
        if order == Chem.BondType.DOUBLE:
            valid = state == (0, 0)
        elif order == Chem.BondType.SINGLE:
            valid = state in {(-1, 0), (0, 1)}
        else:
            valid = False
        if not valid:
            raise res.fail(
                f"unsupported acid/carboxylate state at {name}: charge/H={state}."
            )
        states.append(state)
        orders.append(order)
    if set(orders) != {Chem.BondType.SINGLE, Chem.BondType.DOUBLE}:
        raise res.fail(f"unsupported acid/carboxylate bond orders at {carbon}.")
    return set(oxygens)


def _sidechain(
    mol: Chem.Mol, res: _Residue, external: dict[str, list[tuple[_Residue, str]]]
) -> tuple[set[str], set[frozenset[str]]]:
    family, atoms = res.family, res.atoms
    variable = set()
    bonds: set[frozenset[str]] = set()
    state = ""
    if family in {"ASP", "GLU"}:
        carbon, oxygens = (
            ("CG", ("OD1", "OD2")) if family == "ASP" else ("CD", ("OE1", "OE2"))
        )
        variable = _acid(mol, res, carbon, oxygens)
        bonds = {frozenset((carbon, o)) for o in oxygens}
        state = "acid" if sum(_qh(atoms[o])[1] for o in oxygens) else "anion"
    elif family == "ARG":
        variable = {"NE", "NH1", "NH2"}
        nitrogen_orders = {
            n: mol.GetBondBetweenAtoms(
                atoms["CZ"].GetIdx(), atoms[n].GetIdx()
            ).GetBondType()
            for n in variable
        }
        orders = list(nitrogen_orders.values())
        q = sum(_qh(atoms[n])[0] for n in variable)
        h = sum(_qh(atoms[n])[1] for n in variable)
        if (
            q not in {0, 1}
            or h != 4 + q
            or orders.count(Chem.BondType.SINGLE) != 2
            or orders.count(Chem.BondType.DOUBLE) != 1
            or any(atoms[n].GetFormalCharge() not in {0, 1} for n in variable)
            or any(
                atoms[n].GetFormalCharge() != 0 and order != Chem.BondType.DOUBLE
                for n, order in nitrogen_orders.items()
            )
        ):
            raise res.fail(
                "unsupported arginine charge, hydrogen inventory or resonance bonding."
            )
        state = "neutral" if q == 0 else "cation"
        bonds = {frozenset(("CZ", n)) for n in variable}
    elif family == "HIS":
        variable = {"CG", "ND1", "CE1", "NE2", "CD2"}
        counts = (
            atoms["ND1"].GetTotalNumHs(includeNeighbors=True),
            atoms["NE2"].GetTotalNumHs(includeNeighbors=True),
        )
        state = {(1, 0): "HID", (0, 1): "HIE", (1, 1): "HIP"}.get(counts, "")
        nitrogen_charges = sorted(atoms[n].GetFormalCharge() for n in ("ND1", "NE2"))
        if (
            not state
            or nitrogen_charges != [0, int(state == "HIP")]
            or any(atoms[n].GetFormalCharge() != 0 for n in ("CG", "CE1", "CD2"))
            or any(
                atoms[n].GetTotalNumHs(includeNeighbors=True) != int(n != "CG")
                for n in ("CG", "CE1", "CD2")
            )
        ):
            raise res.fail("unsupported histidine charge or hydrogen inventory.")
    elif family in {"LYS", "CYS", "TYR", "SER", "THR"}:
        name = {"LYS": "NZ", "CYS": "SG", "TYR": "OH", "SER": "OG", "THR": "OG1"}[
            family
        ]
        value = _qh(atoms[name])
        if family == "LYS":
            state = {(0, 2): "neutral", (1, 3): "cation"}.get(value, "")
        elif family == "CYS" and external.get("SG"):
            state = "disulfide" if value == (0, 0) else ""
        else:
            state = {(0, 1): "neutral", (-1, 0): "anion"}.get(value, "")
        if not state:
            raise res.fail(f"unsupported {family} state at {name}: charge/H={value}.")
        variable = {name}
    alias = res.resid.name.upper()
    required = {
        "LYN": {"neutral"},
        "LSN": {"neutral"},
        "LYSH": {"cation"},
        "ARGN": {"neutral"},
        "ASH": {"acid"},
        "ASPH": {"acid"},
        "ASPP": {"acid"},
        "GLH": {"acid"},
        "GLUH": {"acid"},
        "GLUP": {"acid"},
        "CYSH": {"neutral"},
        "CYM": {"anion"},
        "CYN": {"anion"},
        "CYX": {"anion", "disulfide"},
        "CYS1": {"anion", "disulfide"},
        "CYS2": {"anion", "disulfide"},
        "NCYX": {"anion", "disulfide"},
        "CCYX": {"anion", "disulfide"},
        "TYM": {"anion"},
        **{n: {s} for n, s in HISTIDINE_ALIASES.items()},
    }.get(alias)
    if required is not None and state not in required:
        raise res.fail(f"conflicting residue alias {alias} and source state {state}.")
    return variable, bonds


def _backbone(
    mol: Chem.Mol, res: _Residue, external: dict[str, list[tuple[_Residue, str]]]
) -> tuple[set[str], set[frozenset[str]]]:
    variable = {"N"}
    variable_bonds: set[frozenset[str]] = set()
    atoms = res.atoms
    n = atoms["N"]
    n_heavy = sum(a.GetAtomicNum() != 1 for a in n.GetNeighbors())
    n_q, n_h = _qh(n)
    if external.get("N"):
        n_valid = n_q == 0 and n_h == 3 - n_heavy
    else:
        n_valid = n_q in {0, 1} and n_h == 3 + n_q - n_heavy
    if not n_valid:
        raise res.fail(f"unsupported backbone/N-terminal charge/H at N: {(n_q, n_h)}.")
    if external.get("C"):
        if "OXT" in atoms:
            raise res.fail("internal peptide carbon has an unexpected terminal OXT.")
    else:
        if "OXT" not in atoms:
            raise res.fail("missing terminal oxygen OXT at free C terminus.")
        variable |= _acid(mol, res, "C", ("O", "OXT"))
        variable_bonds |= {frozenset(("C", "O")), frozenset(("C", "OXT"))}
    alias = res.resid.name.upper()
    if len(alias) == 4 and alias[1:] in STANDARD_RESNAME_MAP:
        if alias.startswith("N") and external.get("N"):
            raise res.fail(
                "N-terminal residue alias conflicts with peptide connectivity."
            )
        if alias.startswith("C") and external.get("C"):
            raise res.fail(
                "C-terminal residue alias conflicts with peptide connectivity."
            )
    return variable, variable_bonds


def _validate_residue(
    mol: Chem.Mol, res: _Residue, index: dict[int, tuple[_Residue, str]]
) -> None:
    external = _topology(res, index)
    variable, variable_bonds = _sidechain(mol, res, external)
    backbone, backbone_bonds = _backbone(mol, res, external)
    variable |= backbone
    variable_bonds |= backbone_bonds
    atoms = res.atoms
    for name, atom in atoms.items():
        if atom.GetAtomicNum() == 1:
            continue
        if name not in variable and atom.GetFormalCharge() != 0:
            raise res.fail(f"unsupported formal charge at {name}.")
        if (
            name not in variable
            and atom.GetTotalValence()
            != {
                6: 4,
                7: 3,
                8: 2,
                16: 2,
            }[atom.GetAtomicNum()]
        ):
            # RDKit also sanitizes hypervalent sulfur, but MET is a thioether.
            # Fixed bond orders plus ordinary valence constrain the full H count.
            raise res.fail(f"unsupported hydrogen inventory or valence at {name}.")
        for neighbor in atom.GetNeighbors():
            other, other_name = index[neighbor.GetIdx()]
            pair = frozenset((name, other_name))
            if (
                other is not res
                or neighbor.GetAtomicNum() == 1
                or pair in variable_bonds
            ):
                continue
            if name in res.template.aromatic and other_name in res.template.aromatic:
                if not mol.GetBondBetweenAtoms(
                    atom.GetIdx(), neighbor.GetIdx()
                ).GetIsAromatic():
                    raise res.fail(f"unsupported aromatic bond {name}-{other_name}.")
            elif (
                mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType()
                != res.template.bonds[pair]
            ):
                raise res.fail(f"unsupported bond order {name}-{other_name}.")


def preserve_prepared(
    source: Chem.Mol | str | Path,
    *,
    source_is_prepared: bool,
    histidine_evidence: bool,
) -> Molecule:
    """Validate a private copy of a prepared protein, then bind its residue views."""
    if source_is_prepared is not True:
        raise ValueError(
            "Set source_is_prepared=True only after preparation has assigned source "
            "bond orders, charges and hydrogen counts, including deliberately absent H."
        )
    if not isinstance(source, Chem.Mol):
        raise ValueError(
            "Strict protonation preservation requires a prepared RDKit or ProLIF "
            "molecule; prepare PDB input before declaring its chemistry."
        )
    if histidine_evidence:
        raise ValueError(
            "Histidine-only evidence conflicts with the authoritative prepared graph."
        )
    use_segid = (
        source[0]._use_segid
        if isinstance(source, Molecule) and source.n_residues
        else False
    )
    mol = Chem.Mol(source)
    if not mol.GetNumAtoms() or not mol.GetNumConformers():
        raise ValueError("Prepared protein requires atoms and source coordinates.")
    if any(not np.isfinite(conf.GetPositions()).all() for conf in mol.GetConformers()):
        raise ValueError("Prepared protein contains non-finite coordinates.")
    residues, index = _residues(mol, use_segid)
    for residue in residues:
        _hydrogens(residue, index)
    try:
        mol.UpdatePropertyCache(strict=True)
        before = [
            (_qh(a), a.GetNumRadicalElectrons(), a.GetChiralTag())
            for a in mol.GetAtoms()
        ]
        Chem.SanitizeMol(mol)
    except Exception as exc:
        raise ValueError(f"Invalid prepared protein chemistry: {exc}") from exc
    after = [
        (_qh(a), a.GetNumRadicalElectrons(), a.GetChiralTag()) for a in mol.GetAtoms()
    ]
    if before != after:
        raise ValueError(
            "Sanitization changed source charge, hydrogen counts, radicals or "
            "stereochemistry; supply a validated prepared graph."
        )
    for residue in residues:
        _validate_residue(mol, residue, index)
    return Molecule(mol, use_segid=use_segid)
