# ruff: noqa: T201
# Audit command-line tools intentionally print their results.
"""Reproduce source-to-standardizer chemistry changes; no library modifications."""

import json
import sys
import warnings
from collections import Counter
from collections.abc import Callable, Iterator, Sequence
from pathlib import Path
from typing import Protocol, TypedDict

from rdkit import Chem, rdBase
from rdkit.Chem.rdDistGeom import EmbedMolecule

from prolif import Molecule
from prolif.interactions import Cationic, ImplicitHBAcceptor
from prolif.io import MoleculeStandardizer


class SiteSnapshot(TypedDict):
    q: int
    h: int
    rad: int
    hyb: str
    bonds: list[tuple[str, str]]
    donor: bool
    acceptor: bool
    cation: bool
    anion: bool


Snapshot = dict[str, SiteSnapshot]
HeavyInfo = list[tuple[str, str, int, list[float]]]
Chemistry = list[tuple[str, int, int, int]]


class CaseResult(TypedDict):
    family: str
    state: str
    alias: str
    representation: str
    preserve_histidine: bool
    source: Snapshot
    remove_h_source: Snapshot
    warnings: list[str]


class CaseOutcome(CaseResult, total=False):
    parent: Snapshot
    residue: Snapshot
    output_resname: str
    parent_preserved: bool
    residue_preserved: bool
    heavy_order_coords_preserved: bool
    identity_preserved: bool
    source_chemistry_unchanged: bool
    removed_parent: Snapshot
    removed_residue: Snapshot
    h_removal_preserves_parent: bool
    rebuild_preserves_residue: bool
    remove_h_error: str
    error: str


class StandardizerFactory(Protocol):
    def __call__(
        self, *, preserve_histidine: bool
    ) -> Callable[[Chem.Mol], Molecule]: ...


standardizer_factory: StandardizerFactory = MoleculeStandardizer


def atom_name(a: Chem.Atom) -> str:
    return a.GetPDBResidueInfo().GetName().strip()


def key(a: Chem.Atom) -> str:
    i = a.GetPDBResidueInfo()
    return f"{i.GetChainId()}:{i.GetResidueNumber()}:{atom_name(a)}"


def atoms_at(mol: Chem.Mol, number: int) -> dict[str, Chem.Atom]:
    return {
        atom_name(a): a
        for a in mol.GetAtoms()
        if a.GetPDBResidueInfo().GetResidueNumber() == number
    }


def charge(a: Chem.Atom, value: int) -> None:
    a.SetFormalCharge(value)
    a.SetNumExplicitHs(0)
    a.SetNoImplicit(False)
    a.SetNumRadicalElectrons(0)


def fixture(  # noqa: PLR0912
    family: str, state: str, alias: str | None = None
) -> tuple[Chem.Mol, int, list[str]]:
    one = {
        "ASP": "D",
        "GLU": "E",
        "LYS": "K",
        "ARG": "R",
        "CYS": "C",
        "TYR": "Y",
        "SER": "S",
        "THR": "T",
        "NTERM": "A",
        "NPRO": "P",
        "CTERM": "A",
        "SS": "C",
    }[family]
    sequence = (
        "CGC"
        if family == "SS"
        else (
            one + "G"
            if family in {"NTERM", "NPRO"}
            else "G" + one
            if family == "CTERM"
            else "G" + one + "G"
        )
    )
    mol = Chem.MolFromSequence(sequence)
    number = 1 if family in {"NTERM", "NPRO", "SS"} else 2
    aa = atoms_at(mol, number)
    for a in mol.GetAtoms():
        a.GetPDBResidueInfo().SetChainId("A")
    if family in {"ASP", "GLU"}:
        n1, n2, center = (
            ("OD1", "OD2", "CG") if family == "ASP" else ("OE1", "OE2", "CD")
        )
        swapped = state.endswith("1")
        single, double = (n1, n2) if swapped else (n2, n1)
        for n in (n1, n2):
            charge(aa[n], -int(state.startswith("anion") and n == single))
            mol.GetBondBetweenAtoms(aa[center].GetIdx(), aa[n].GetIdx()).SetBondType(
                Chem.BondType.SINGLE if n == single else Chem.BondType.DOUBLE
            )
        sites = [n1, n2]
    elif family == "ARG":
        double = "NH1" if state.endswith("1") else "NH2"
        for n in ("NE", "NH1", "NH2"):
            charge(aa[n], int(state.startswith("cation") and n == double))
            mol.GetBondBetweenAtoms(aa["CZ"].GetIdx(), aa[n].GetIdx()).SetBondType(
                Chem.BondType.DOUBLE if n == double else Chem.BondType.SINGLE
            )
        sites = ["NE", "NH1", "NH2"]
    elif family == "SS":
        b = atoms_at(mol, 3)["SG"]
        for a in (aa["SG"], b):
            charge(a, 0)
        rw = Chem.RWMol(mol)
        rw.AddBond(aa["SG"].GetIdx(), b.GetIdx(), Chem.BondType.SINGLE)
        mol = rw.GetMol()
        aa = atoms_at(mol, number)
        sites = ["SG"]
    else:
        site = {
            "LYS": "NZ",
            "CYS": "SG",
            "TYR": "OH",
            "SER": "OG",
            "THR": "OG1",
            "NTERM": "N",
            "NPRO": "N",
            "CTERM": "OXT",
        }[family]
        charge(aa[site], 1 if state == "cation" else -1 if state == "anion" else 0)
        sites = [site]
    if alias:
        numbers = {1, 3} if family == "SS" else {number}
        for a in mol.GetAtoms():
            if a.GetPDBResidueInfo().GetResidueNumber() in numbers:
                a.GetPDBResidueInfo().SetResidueName(alias)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addResidueInfo=True)
    # Distinct H names deliberately exercise connectivity recovery. They cannot
    # accidentally match template H names on a different atom. Coordinates are
    # embedded with the correct source bonding before any template is applied.
    for a in mol.GetAtoms():
        info = a.GetPDBResidueInfo()
        info.SetSerialNumber(a.GetIdx() + 1)
        if a.GetAtomicNum() == 1:
            info.SetName(f"Q{a.GetIdx():03d}")
    assert EmbedMolecule(mol, randomSeed=0xA70) == 0
    return mol, number, sites


HB = ImplicitHBAcceptor()
PATTERNS = {
    "donor": Chem.MolFromSmarts(HB.donor),
    "acceptor": Chem.MolFromSmarts(HB.acceptor),
    "cation": Cationic().lig_pattern,
    "anion": Cationic().prot_pattern,
}


def snapshot(mol: Chem.Mol, number: int, sites: Sequence[str]) -> Snapshot:
    matches = {
        kind: {m[0] for m in mol.GetSubstructMatches(p)} for kind, p in PATTERNS.items()
    }
    aa = atoms_at(mol, number)
    return {
        n: {
            "q": aa[n].GetFormalCharge(),
            "h": aa[n].GetTotalNumHs(includeNeighbors=True),
            "rad": aa[n].GetNumRadicalElectrons(),
            "hyb": str(aa[n].GetHybridization()),
            "bonds": sorted(
                (
                    key(b),
                    str(
                        mol.GetBondBetweenAtoms(
                            aa[n].GetIdx(), b.GetIdx()
                        ).GetBondType()
                    ),
                )
                for b in aa[n].GetNeighbors()
                if b.GetAtomicNum() != 1
            ),
            "donor": aa[n].GetIdx() in matches["donor"],
            "acceptor": aa[n].GetIdx() in matches["acceptor"],
            "cation": aa[n].GetIdx() in matches["cation"],
            "anion": aa[n].GetIdx() in matches["anion"],
        }
        for n in sites
    }


def state_sig(snap: Snapshot) -> dict[str, tuple[int, int, int]]:
    return {n: (s["q"], s["h"], s["rad"]) for n, s in snap.items()}


REPRESENTATIONS = ("all_graph_h", "no_graph_h", "counted_h", "fixed_h_counts")


def represent(source: Chem.Mol, representation: str) -> Chem.Mol:
    if representation != "all_graph_h":
        source = Chem.RemoveHs(
            source,
            updateExplicitCount=representation in {"counted_h", "fixed_h_counts"},
        )
    if representation == "fixed_h_counts":
        for atom in source.GetAtoms():
            atom.SetNoImplicit(True)
        Chem.SanitizeMol(source)
    assert all(a.GetNumRadicalElectrons() == 0 for a in source.GetAtoms())
    return source


def heavy_info(mol: Chem.Mol) -> HeavyInfo:
    return [
        (
            key(a),
            a.GetPDBResidueInfo().GetResidueName(),
            a.GetAtomicNum(),
            list(mol.GetConformer().GetAtomPosition(a.GetIdx())),
        )
        for a in mol.GetAtoms()
        if a.GetAtomicNum() != 1
    ]


def chemistry(mol: Chem.Mol) -> Chemistry:
    return [
        (
            key(a),
            a.GetFormalCharge(),
            a.GetTotalNumHs(includeNeighbors=True),
            a.GetNumRadicalElectrons(),
        )
        for a in mol.GetAtoms()
        if a.GetAtomicNum() != 1
    ]


def case_specs() -> Iterator[tuple[str, str, str]]:  # noqa: PLR0912
    for family, neutral_alias in [("ASP", "ASH"), ("GLU", "GLH")]:
        for state in ("acid2", "acid1", "anion2", "anion1"):
            for name in (family, neutral_alias):
                yield family, state, name
    for name in ("LYS", "LYN", "LSN", "LYSH"):
        for state in ("neutral", "cation"):
            yield "LYS", state, name
    for name in ("ARG", "ARGN"):
        for state in ("neutral2", "neutral1", "cation2", "cation1"):
            yield "ARG", state, name
    for name in ("CYS", "CYM", "CYX", "CYN", "CYSH", "CYS1", "CYS2"):
        for state in ("neutral", "anion"):
            yield "CYS", state, name
    for name in ("CYS", "CYX", "CYM", "CYS2"):
        yield "SS", "disulfide", name
    for name in ("TYR", "TYM"):
        for state in ("neutral", "anion"):
            yield "TYR", state, name
    for family in ("SER", "THR"):
        for state in ("neutral", "anion"):
            yield family, state, family
    for family, names in [
        ("NTERM", ("ALA", "NALA")),
        ("NPRO", ("PRO", "NPRO")),
        ("CTERM", ("ALA", "CALA")),
    ]:
        for name in names:
            for state in (
                ("neutral", "anion") if family == "CTERM" else ("neutral", "cation")
            ):
                yield family, state, name


def run_case(
    family: str, state: str, alias: str, representation: str, preserve: bool
) -> CaseOutcome:
    source, number, sites = fixture(family, state, alias)
    source = represent(source, representation)
    before = chemistry(source)
    row: CaseOutcome = {
        "family": family,
        "state": state,
        "alias": alias,
        "representation": representation,
        "preserve_histidine": preserve,
        "source": snapshot(source, number, sites),
        "remove_h_source": snapshot(Chem.RemoveHs(source), number, sites),
        "warnings": [],
    }
    with warnings.catch_warnings(record=True) as caught:
        try:
            result = standardizer_factory(preserve_histidine=preserve)(source)
            res = next(r for r in result.residues.values() if r.resid.number == number)
            row["parent"] = snapshot(result, number, sites)
            row["residue"] = snapshot(res, number, sites)
            row["output_resname"] = res.resid.name
            row["parent_preserved"] = state_sig(row["source"]) == state_sig(
                row["parent"]
            )
            row["residue_preserved"] = state_sig(row["source"]) == state_sig(
                row["residue"]
            )
            old, new = heavy_info(source), heavy_info(result)
            row["heavy_order_coords_preserved"] = [x[:1] + x[2:] for x in old] == [
                x[:1] + x[2:] for x in new
            ]
            row["identity_preserved"] = old == new
            row["source_chemistry_unchanged"] = before == chemistry(source)
            try:
                removed = Molecule(Chem.RemoveHs(result))
                removed_res = next(
                    r for r in removed.residues.values() if r.resid.number == number
                )
                row["removed_parent"] = snapshot(removed, number, sites)
                row["removed_residue"] = snapshot(removed_res, number, sites)
                row["h_removal_preserves_parent"] = state_sig(
                    row["parent"]
                ) == state_sig(row["removed_parent"])
                row["rebuild_preserves_residue"] = state_sig(
                    row["residue"]
                ) == state_sig(row["removed_residue"])
            except Exception as exc:
                row["remove_h_error"] = f"{type(exc).__name__}: {exc}"
        except Exception as exc:
            row["error"] = f"{type(exc).__name__}: {exc}"
        row["warnings"] = sorted({str(w.message) for w in caught})
    return row


def main() -> None:
    rdBase.DisableLog("rdApp.warning")
    rdBase.DisableLog("rdApp.error")
    results = [
        run_case(*spec, representation, preserve)
        for spec in case_specs()
        for representation in REPRESENTATIONS
        for preserve in (False, True)
    ]
    output = Path(
        sys.argv[1] if len(sys.argv) > 1 else ".context/titratable-audit/results.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps({"rdkit": rdBase.rdkitVersion, "cases": results}, indent=2) + "\n"
    )
    print(f"Wrote {len(results)} cases to {output}")
    for mode in (False, True):
        selected = [r for r in results if r["preserve_histidine"] == mode]
        print(
            "preserve_histidine",
            mode,
            Counter(
                "error"
                if "error" in r
                else "preserved"
                if r["residue_preserved"]
                else "changed"
                for r in selected
            ),
        )
    print(
        "Source RemoveHs changes:",
        sum(state_sig(r["source"]) != state_sig(r["remove_h_source"]) for r in results),
    )


if __name__ == "__main__":
    main()
