"""Prepared peptide fixtures independent of ProLIF residue templates."""

from collections.abc import Mapping

from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

HYDROGEN_REPRESENTATIONS = ("graph", "removed", "counted", "fixed-counted")


def atoms_by_name(mol: Chem.Mol, number: int = 2) -> dict[str, Chem.Atom]:
    return {
        atom.GetPDBResidueInfo().GetName().strip(): atom
        for atom in mol.GetAtoms()
        if atom.GetPDBResidueInfo().GetResidueNumber() == number
    }


def set_state(
    mol: Chem.Mol, state: str, *, number: int = 2, alias: str | None = None
) -> None:
    """Assign literal H/charge and bond-order evidence to one peptide residue."""
    atoms = atoms_by_name(mol, number)

    def charge_h(name: str, charge: int, hydrogens: int) -> None:
        atom = atoms[name]
        atom.SetFormalCharge(charge)
        atom.SetNumExplicitHs(hydrogens)
        atom.SetNoImplicit(True)

    def bond(first: str, second: str, order: Chem.BondType) -> None:
        mol.GetBondBetweenAtoms(
            atoms[first].GetIdx(), atoms[second].GetIdx()
        ).SetBondType(order)

    family, _, variant = state.partition(":")
    if family in {"ASP", "GLU"}:
        carbon, oxygens = (
            ("CG", ("OD1", "OD2")) if family == "ASP" else ("CD", ("OE1", "OE2"))
        )
        single = oxygens[0] if variant.endswith("1") else oxygens[1]
        for oxygen in oxygens:
            is_single = oxygen == single
            charge_h(
                oxygen,
                -int(is_single and variant.startswith("anion")),
                int(is_single and variant.startswith("acid")),
            )
            bond(
                carbon,
                oxygen,
                Chem.BondType.SINGLE if is_single else Chem.BondType.DOUBLE,
            )
    elif family == "LYS":
        cation = variant == "+"
        charge_h("NZ", int(cation), 3 if cation else 2)
    elif family == "ARG":
        cation = variant.startswith("+")
        double = "NE" if variant.endswith("NE") else "NH1"
        for nitrogen in ("NE", "NH1", "NH2"):
            is_double = nitrogen == double
            charge = int(cation and is_double)
            h = (1 if nitrogen == "NE" else 2) - int(is_double) + charge
            charge_h(nitrogen, charge, h)
            bond(
                "CZ",
                nitrogen,
                Chem.BondType.DOUBLE if is_double else Chem.BondType.SINGLE,
            )
    elif family in {"CYS", "TYR", "SER", "THR"}:
        target = {"CYS": "SG", "TYR": "OH", "SER": "OG", "THR": "OG1"}[family]
        charge_h(target, -int(variant == "-"), int(variant != "-"))
    elif family == "HIS":
        nd_h, ne_h, nd_charge = {
            "HID": (1, 0, 0),
            "HIE": (0, 1, 0),
            "HIP": (1, 1, 1),
        }[variant]
        charge_h("ND1", nd_charge, nd_h)
        charge_h("NE2", 0, ne_h)
    else:
        raise AssertionError(f"Unknown fixture state: {state}")
    if alias is not None:
        for atom in atoms.values():
            atom.GetPDBResidueInfo().SetResidueName(alias)


def peptide(
    sequence: str, states: Mapping[int, tuple[str, str | None]] | None = None
) -> Chem.Mol:
    mol = Chem.MolFromSequence(sequence)
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        info.SetChainId("A")
        info.SetSerialNumber(atom.GetIdx() + 1)
        # RDKit's sequence builder puts ILE CD1 on its CG2. Standard protein names
        # put the ethyl branch on CG1; fix names without changing the molecule.
        if info.GetResidueName() == "ILE":
            name = info.GetName().strip()
            if name in {"CG1", "CG2"}:
                info.SetName(" CG2" if name == "CG1" else " CG1")
    for number, (state, alias) in (states or {}).items():
        set_state(mol, state, number=number, alias=alias)
    Chem.SanitizeMol(mol)
    Compute2DCoords(mol)
    return mol


def represent_hydrogens(mol: Chem.Mol, representation: str) -> Chem.Mol:
    if representation in {"counted", "fixed-counted"}:
        result = Chem.Mol(mol)
        for atom in result.GetAtoms():
            hydrogens = atom.GetTotalNumHs()
            atom.SetNoImplicit(representation == "fixed-counted")
            atom.SetNumExplicitHs(hydrogens)
        Chem.SanitizeMol(result)
        return result
    result = Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    for atom in result.GetAtoms():
        if atom.GetAtomicNum() == 1:
            # Unknown names assert source connectivity without claiming template names.
            atom.GetPDBResidueInfo().SetName(f"Q{atom.GetIdx():03d}")
    if representation == "graph":
        return result
    if representation in {"removed", "removed-unsanitized"}:
        return Chem.RemoveHs(result, sanitize=representation == "removed")
    raise AssertionError(f"Unknown hydrogen representation: {representation}")


def state_peptide(
    state: str, representation: str = "counted", alias: str | None = None
) -> Chem.Mol:
    code = {
        "ASP": "D",
        "GLU": "E",
        "LYS": "K",
        "ARG": "R",
        "CYS": "C",
        "TYR": "Y",
        "HIS": "H",
        "SER": "S",
        "THR": "T",
    }[state.split(":")[0]]
    return represent_hydrogens(
        peptide(f"G{code}G", {2: (state, alias)}), representation
    )


def graph_signature(mol: Chem.Mol) -> tuple:
    """Observable chemistry, atom identity/order, and every conformer's coordinates."""
    atoms = []
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        atoms.append(
            (
                atom.GetAtomicNum(),
                atom.GetIsotope(),
                atom.GetFormalCharge(),
                atom.GetTotalNumHs(includeNeighbors=True),
                atom.GetChiralTag(),
                atom.GetIsAromatic(),
                atom.GetNumRadicalElectrons(),
                info.GetName(),
                info.GetResidueName(),
                info.GetResidueNumber(),
                info.GetChainId(),
                info.GetSerialNumber(),
                info.GetInsertionCode(),
                info.GetAltLoc(),
                info.GetOccupancy(),
                info.GetTempFactor(),
            )
        )
    bonds = tuple(
        (
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond.GetBondType(),
            bond.GetIsAromatic(),
            bond.GetStereo(),
            tuple(bond.GetStereoAtoms()),
        )
        for bond in mol.GetBonds()
    )
    coords = tuple(
        tuple(map(tuple, conf.GetPositions())) for conf in mol.GetConformers()
    )
    return tuple(atoms), bonds, coords


def local_state(mol: Chem.Mol, names: tuple[str, ...], number: int = 2) -> tuple:
    atoms = atoms_by_name(mol, number)
    return tuple(
        (
            atoms[name].GetFormalCharge(),
            atoms[name].GetTotalNumHs(includeNeighbors=True),
        )
        for name in names
    )
