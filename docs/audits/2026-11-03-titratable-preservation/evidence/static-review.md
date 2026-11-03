# Non-histidine protonation audit: independent static review

Audit scope: source inspection plus narrow synthetic probes of atom-name aliases,
extended PDB residue-name parsing, and legacy/opt-in connected-parent behavior.
No production code, tests, commits, or external systems were changed. The broader
chemical-state, disulfide, and interaction matrices are owned by the root audit.

Evidence labels below distinguish code guarantees from experimentally observed
behavior. Synthetic observations use RDKit 2025.09.6 in the existing `.venv`.

## Prioritized findings

### 1. Non-histidine templates overwrite selected source charges without resolving source evidence

**Proven from implementation.** Name normalization chooses the non-histidine
template at `prolif/io/molecule_standardizer.py:265-289`. Only histidines receive
evidence resolution and post-template verification (`:235-243`, `:314-315`).
Built-in CIF application removes all residue bonds and recreates them by atom
name (`prolif/io/template_engine.py:114-116`, `:182-190`, `:235-275`). Charge
overrides are applied without checking the source value (`:293-304`).

| Source names with intact metadata | Selected template | Forced site charge / fixed bond pattern | Preservation implication |
| --- | --- | --- | --- |
| ASP/NASP/CASP | ASP | OD2 = -1; CG=OD1, CG-OD2 | A neutral source acid still named ASP does not select ASH from H or charge. |
| ASH/ASPH/ASPP | ASH | OD2 = 0; CG=OD1, CG-OD2; OD2-HD2 if present | Supports a name-selected neutral state, but does not preserve alternative oxygen labeling/attachment automatically. |
| GLU/NGLU/CGLU | GLU | OE2 = -1; CD=OE1, CD-OE2 | Same limitation as ASP. |
| GLH/GLUH/GLUP | GLH | OE2 = 0; CD=OE1, CD-OE2; OE2-HE2 if present | Same limitation as ASH. |
| LYS/LYN/LSN/LYSH/NLYS/CLYS | LYS | NZ = +1 | Neutral-state names LYN and LSN collapse into the cationic template. |
| ARG/ARGN/NARG/CARG | ARG | NH2 = +1; CZ=NH2, CZ-NH1, CZ-NE | ARGN does not select a neutral template. |
| CYS/CYSH/NCYS/CCYS | CYS, except global GROMOS inference | SG = 0; CB-SG, SG-HG if present | Source thiolate charge is overwritten; context can change selection. |
| CYX/CYM/CYN/CYS1/CYS2/NCYX/CCYX | CYX | SG = -1; CB-SG | Thiolate and disulfide names share one charge assignment. |
| TYR/NTYR/CTYR | TYR | CZ-OH, OH-HH if present; no OH charge override | Existing OH formal charge survives this stage; it is incorrect to claim all TYR is forcibly neutralized. |

Alias evidence is in `prolif/io/constants.py:64-101`, `:111-119`, `:154-183`;
charge overrides are `:204-216`. Acid template bonds are in
`prolif/data/molecule_standardizer/templates/standard_aa.cif:218-219`,
`:267-269`, `:508-509`, `:563-565`; ARG `:112-119`; LYS `:931-936`;
CYS/CYX `:312-315`, `:357-359`; TYR `:1224-1225`.

**Important limit:** The table states assigned site charges, not guaranteed final
valid protonation states. Existing graph H, explicit-H counters, `NoImplicit`,
other atom charges, radicals and parent connectivity can produce extra inferred H,
radicals, or sanitization failures. The implementation does not reset all atom
charges, remove all H, or infer a consistent source state. A conflict can fail
sanitization rather than silently produce a conventional target state. Formal
charges on atoms absent from the override table remain carried over from source.

### 2. Terminal/backbone oxygen aliases can silently lose their bonds

**Experimentally confirmed.** All nine alternative names declared under the XML
`Protein` namespace failed the synthetic heavy-only alanine probe:

- OXT alternatives: OT, OT2, O2, OCT2, OC2.
- O alternatives: O1, OT1, OCT1, OC1.

The source oxygen had neighbor C. After standardization every aliased oxygen had
no neighbors in the residue, with inferred total H = 2. No warning or exception
was produced. Canonical OXT and O controls retained their C bond. Legacy mode
left the connected parent bond intact but returned the damaged residue. With
`preserve_histidine=True`, both returned connected parent and residue lost the
bond. These inputs contain no histidine.

**Mechanism:** `parse_altnames` stores each XML residue/type bucket separately and
does not inherit the `type="Protein"` namespace (`prolif/io/xml.py:20-35`).
`assign_intra_props` looks up only the residue-specific bucket
(`prolif/io/template_engine.py:235-240`). It never merges the common Protein
aliases declared at `standard_aa_name.xml:4-10`. Missing names then cause template
bonds to be skipped at `template_engine.py:261-264`, after all source bonds were
stripped. Terminal O aliases are recognized for the *heavy atom count exemption*
(`constants.py:293-296`, `molecule_standardizer.py:473-479`), which does not repair
their topology. Missing bonds do not change the heavy atom count, so the soft
count check does not catch this.

This is a broader connectivity problem than choosing a terminal protonation
state, and propagating residue chemistry to the parent makes it relevant to
the opt-in path even in histidine-free molecules.

### 3. Four-letter non-histidine aliases lose state evidence at the PDB boundary

**Experimentally confirmed, both standardizer modes.** Nine deterministic extended
PDB fixtures use columns 18-21 for the four-letter name, following the existing
histidine convention. Initial parsing kept atom/bond counts and the source site
charge in every fixture; the name was already truncated before standardization.

| Source name | Parsed name | Intact in-memory result at site | PDB-path result at site |
| --- | --- | --- | --- |
| ASPH, ASPP | ASP | OD2 charge 0 | OD2 charge -1 |
| GLUH, GLUP | GLU | OE2 charge 0 | OE2 charge -1 |
| CYS1 | CYS | SG charge -1 | SG charge 0, with unknown-forcefield CYS warning |
| LYSH | LYS | NZ charge +1 | NZ charge +1 |
| ARGN | ARG | NH2 charge +1 | NH2 charge +1 |
| NALA | NAL | N charge +1 | Missing-template ValueError for NAL |
| CALA | CAL | OXT charge -1 | Missing-template ValueError for CAL |

The source ARGN fixture was neutral; its intact-name result is already changed
by template choice. Therefore ARGN truncation is a name-loss observation, not an
additional charge discrepancy between direct and PDB input. CYS1 is a single
thiolate test here, not a disulfide test.

**Mechanism:** legacy mode delegates directly to RDKit's three-character PDB
reader (`molecule_standardizer.py:372-382`). The preserving reader restores an
extended name only when it belongs to histidine `_ALIASES`; all other names use
`line[17:20]` (`prolif/io/_histidine.py:313-325`). Thus the histidine preservation
fix intentionally does not extend non-histidine alias preservation. No claim is
made that four-letter names are standard PDB; these probes exercise the extended
format already supported by this project for histidines.

### 4. Legacy and opt-in mode expose different whole-molecule chemistry for every residue

**Proven from code and confirmed in AKA peptide probe.** Both modes call the same
non-histidine template engines. Legacy updates residue objects and parent atom
names, but does not transfer standardized bonds or charges to the connected
parent (`molecule_standardizer.py:275-284`, `:336-349`).

Opt-in creates a parent copy, applies `update_parent` for *every* residue, then
sanitizes and returns it (`:244`, `:316-317`, `:339-347`). `update_parent` copies
atom properties and reconstructs intra-residue bonds while retaining
inter-residue bonds (`_histidine.py:329-350`). Thus the option changes more than
histidine at the connected-graph representation level, even though the
non-histidine template choices remain the same.

In a heavy-only AKA peptide with neutral source Lys, legacy parent NZ charge was
0 while residue NZ was +1; opt-in parent and residue NZ were both +1. Rebuilding
`Molecule(Chem.RemoveHs(result))` from a legacy result therefore starts from
different charge/bond data than its standardized residue collection. This is not
evidence that H removal itself changed protonation.

Also, legacy returned standardized residues lack `_parent_molecule` in this
probe, whereas opt-in residues have it. Fresh engine results are `Residue` objects
(`template_engine.py:124`); legacy merely replaces the residue collection, while
opt-in constructs `Molecule(..., residues=new_residues)`, whose constructor binds
each copied residue to a shared geometry parent (`prolif/molecule.py:108-118`).
This can change implicit geometry handling; contact-level consequences are left
to the root audit's actual detector probes.

### 5. Disulfide and termini decisions are not source-state resolvers

**Static findings; defer exact chemistry outcomes to dynamic matrix.**

- There is no S-S bond/neighbor-aware charge choice in the CIF engine. CYX gets SG
  -1 whether source SG is a free thiolate or bonded to another residue. The
  `MAX_DISULPHIDE_LENGTH` constant exists but has no call sites in production
  code (`constants.py:30`; repository-wide symbol search).
- `is_peptide_bond` means *any bond between different residues*, explicitly
  including disulfides (`prolif/utils.py:217-229`). The fragment adjustment in
  `molecule_standardizer.py:319-332` therefore applies to S-S endpoints and other
  crosslinks when their rebuilt residue atoms have inferred H. It is not a
  chemically specific terminal-state assignment.
- `forcefield_guesser` uses the names of all residues, with fixed priority AMBER,
  CHARMM, GROMOS, OPLS (`molecule_standardizer.py:405-412`). GROMOS inference makes
  generic CYS choose CYX (`:453-454`); unknown forcefield chooses neutral CYS and
  warns (`:446-451`). Adding an unrelated forcefield-marker residue can therefore
  alter generic CYS template selection, irrespective of its local charge/H.
- N/C prefix aliases normalize to ordinary residue templates. No standard charge
  override exists for backbone N or terminal OXT. Their source formal charges
  remain, while conventional C=O/C-OXT and N-H bonds are reconstructed if named
  atoms exist. Charge-less inputs do not acquire an explicit terminal-state
  assignment merely from NALA/CALA names.
- Missing heavy atoms only produce a warning; unknown atom names or duplicate
  normalized names are not comprehensively validated. `nm_2_idx[nm] = index`
  overwrites earlier name collisions (`template_engine.py:248`) and missing
  template bonds are skipped. This bounds any preservation claim.
- All H bonds are stripped with other bonds. Recognized H names are rebound from
  templates; remaining unbound H are attached to the nearest heavy atom within
  2 A (`template_engine.py:277-278`, `:314-367`). This is topology repair, not
  verification of original H attachment. Non-histidines receive no inventory
  completeness or assignment-agreement checks.

## Existing test coverage and gaps

`tests/test_constants.py:50-75` checks mapping counts and a few charge dictionary
entries; `:83-86` checks terminal alias membership. These are not behavioral
preservation checks. `tests/test_io.py:190-216` tests CYS name conversion and HSD;
`:331-346` tests single-forcefield marker choices. Mixed-forcefield precedence
and local CYS source evidence are not tested there.

`tests/test_io.py:266-318` checks custom TPO template application/equivalence,
`:361-377` exercises missing-heavy-atom warnings, `:396-448` exercises routine
standardization and BEN bond topology, `:450-473` exercises nearest-heavy H
attachment, and `:475-493` checks peptide N exclusion. None asserts preservation
across non-histidine charge/H/name conflicts, neutral Lys/Arg aliases, protonated
Asp/Glu with generic names, disulfide vs thiolate, Tyr phenolate, or terminal
state variants.

`tests/test_histidine_preservation.py:63-108` covers mixed and aliased histidines;
`:383-416` tests AHG connected peptide and backbone behavior. Its non-histidine
residues in that peptide are alanine/glycine, so it does not cover copying
titratable non-histidine template changes to the connected parent. Repository
test-name/content searches did not find explicit ASH/GLH/LYN/ARGN/CYM/disulfide
standardizer regression matrices or terminal oxygen alias connectivity tests.

High-value coverage additions, without proposing an implementation here:
source charge/H/name agreement and conflict cases per residue; canonical and
alternative oxygen resonance labeling; graph-H/implicit-H/NoImplicit variants;
direct molecule vs PDB aliases; all Protein O/OXT aliases; free thiolate vs S-S;
mixed-forcefield marker contexts; termini and caps; both modes' residue/parent/
H-removed reconstructions; actual detector outputs after standardization.

## Reproduction and limits

Run from repository root:

```sh
MISE_UV_VERSION=0.11.24 uv run --no-sync python .context/titratable-audit/alias_probes.py
```

Artifacts: `alias_probes.py`, `alias-probe-results.json`,
`alias-probe-summary.txt`, and `alias-*.pdb` in this directory. JSON retains source,
parsed, parent, residue, site H/charge/radical counts, atom/bond counts, and
captured Python warnings/errors. The printed summary is abbreviated; JSON is the
full observation record. The conformer generator emitted an RDKit UFF warning
for the charged CYS fixture but returned success and valid coordinates; it did
not stop standardization.

No universal chemical preservation conclusion follows from matching formal
charges or donor eligibility alone. These probes establish the listed concrete
boundaries; broader biological contexts, modified residues, metal coordination,
missing atom inventories and arbitrary custom templates remain outside them.
