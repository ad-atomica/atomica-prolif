Helper functions
================

.. automodule:: prolif.utils
    :members: get_residues_near_ligand, select_over_trajectory, pandas_series_to_bv, to_dataframe, to_bitvectors

.. automodule:: prolif.io.molecule_standardizer
    :members: MoleculeStandardizer

Preserving source histidine assignments
---------------------------------------

Template selection and hydrogen representation are separate decisions. The
standardizer's legacy default normalizes residue names and selects a cationic
template for generic HIS, regardless of supplied ring hydrogens. To preserve
an already prepared receptor's histidine states, opt in explicitly::

    from rdkit import Chem
    from prolif import Molecule
    from prolif.io import MoleculeStandardizer

    standardizer = MoleculeStandardizer(preserve_histidine=True)
    receptor = standardizer(
        source_receptor,
        histidine_states={"HIS100.A": "HIE", "HIS125.A": "HIP"},
        complete_histidine_hydrogens=["HIS42.A"],
    )
    implicit_receptor = Molecule(Chem.RemoveHs(receptor))

Evidence must cover **every histidine**, using any combination of:

* Explicit state names: HID/HSD (neutral, ND1 hydrogen), HIE/HSE (neutral,
  NE2 hydrogen), HIP/HSP (cationic, both), or the corresponding built-in
  forcefield aliases (including HISA/HISB/HISH).
* ``histidine_states``: original residue keys (strings or ``ResidueId``) mapped
  to canonical HID/HIE/HIP assignments from source preparation.
* ``complete_histidine_hydrogens``: original residue keys for which preparation
  guarantees a complete inventory of graph H attached to the ring nitrogens.
  This can encode mixed states even if every residue is called HIS. It does
  not require all protein hydrogens to be present.

Generic HIS alone is ambiguous. A single observed hydrogen in an input of
unknown completeness cannot prove the other ring nitrogen was intentionally
left unprotonated. Completeness is a caller assertion, never inferred from a
PDB file or unstandardized valence. All evidence must agree: observed graph H
may corroborate a partial inventory, but names, attachments and explicit
assignments cannot contradict each other. Missing evidence, malformed ring/H
connectivity and conflicting evidence raise residue-specific ``ValueError``
messages before templates run. No protonation state is predicted.

Preservation mode returns a copy and retains source histidine residue names
and keys, atom correspondence, order and coordinates. HID/HIE/HIS are internal
template choices; they do not replace a public HIS or HSD key. Output validation
checks hydrogen counts/attachments, source heavy-atom connectivity, aromatic
ring chemistry, absence of ring radicals and net ring charge,
allowing equivalent resonance representations. Custom template overrides must
also satisfy these checks. The connected parent and residues both receive the
standardized chemistry, allowing the illustrated ``RemoveHs`` analysis copy to
retain implicit donor/acceptor/ionic eligibility. Removing atoms manually without
retaining hydrogen counts is not the supported operation.

For PDB paths in either mode, recognized four-letter residue aliases in extended
columns 18-21 are restored before template selection. Unknown extensions raise
an error instead of being truncated. This uses the first model and requires
unique decimal atom serials; other serial encodings require a caller-provided
molecule with intact metadata. If an external parser already truncated HISB to
HIS, use explicit source evidence; the lost assignment cannot be recovered.

Common protein atom aliases (including O1/OT1 and O2/OT2/OXT variants) are applied
alongside residue-specific aliases. A custom template's explicit atom names take
precedence. Duplicate names after alias mapping raise an error instead of silently
dropping source bonds.

The default ``preserve_histidine=False`` remains compatible with existing calls,
including in-place mutation of ProLIF inputs and residue-name normalization.
Evidence arguments require the option so they cannot be silently ignored.
Other residues retain existing template/name normalization behavior in both
modes. This capability is scoped to histidine; it does not establish
protonation preservation for every titratable residue, modified histidine or
metal-coordinated histidine with unsupported ring-N connectivity.

.. automodule:: prolif.io.template_engine
    :members: TemplateEngine, RDKitMolTemplateEngine, CIFTemplateEngine, assign_bond_orders_from_template, strip_bonds, assign_intra_props

.. automodule:: prolif.io.cif
    :members: cif_template_reader

.. automodule:: prolif.exceptions
    :members: OptionalException, FragmentedResidueError, trigger
