---
orphan: true
---

# Broader source protonation preservation

Alex invoked implement-suggestions after the extensive titratable-residue audit.
The approved work is to preserve the evidence in Git, repair parsing/connectivity,
extend strict source preservation and check graph/detector consistency. Work in
this Conductor checkout; retain existing commits and legacy API compatibility.
Do not release, merge, or modify the Atomica consumer as part of library work.

## Design and decisions

The caller's assigned chemical graph is authoritative. Introduce
`MoleculeStandardizer(preserve_protonation=True)` with per-call
`source_is_prepared=True`. This flag is a substantive declaration: the input
RDKit/ProLIF molecule has intentionally assigned bond orders, charges and total
hydrogen counts, including deliberately absent hydrogens. A parser's inferred
valences alone do not establish preparation. Reject missing declaration and
path/string input with preparation instructions. Do not infer protonation or
offer a best-effort fallback. Reject conflicting histidine-only evidence options
in this path. The existing histidine-only option continues to work as documented.

The strict path validates and returns a copy of the prepared connected graph;
it does not apply chemistry-replacing templates to it. Using a known graph
directly avoids reconstructing a valid state and subsequently checking whether
that reconstruction destroyed it. Retain residue names, atom names/order,
isotopes, coordinates, connectivity, charge, total H and stereochemistry. Build
residues from the validated connected graph so the exact same chemistry is used
before and after supported H removal. Preserve scientifically equivalent source
carboxylate/guanidinium resonance encodings instead of canonicalizing their charge.

Resolve/validate all supported titratable groups: Asp/Glu acid versus carboxylate
and proton location; Lys neutral/cation; Arg neutral/cation and resonance;
Cys thiol/thiolate/disulfide by connectivity; Tyr phenol/phenolate; histidine
HID/HIE/HIP; Ser/Thr neutral/alkoxide; neutral/cationic N termini (including proline)
and acid/anionic C termini. Explicit variant aliases must agree with the source
state; generic family names may contain any supported source state. Generic CYS
is interpreted locally, never by unrelated force-field marker residues.

Validate standard residue heavy-atom inventories and topology, supported peptide
and disulfide links, hydrogen attachment and recognized names, valence/radicals,
aromaticity and allowed local states. Unsupported modified/metal-bound residues,
missing atoms, unknown residues and custom chemistry must fail with a residue/
atom-specific error. Do not accept an arbitrary sanitized but chemically different
residue merely because its total charge is plausible. Unmodified standard protein
residues are the qualification boundary, not every possible biomolecule.

For existing template-based paths, merge common Protein atom aliases with the
residue-specific aliases; reject alias collisions instead of silently dropping
bonds. Preserve recognized extended residue names at library-owned PDB parsing
boundaries. Keep name normalization in legacy modes. Make preserving-mode residue
views derive from the sanitized connected parent, avoiding stale hybridization.

## Alternatives considered

- Per-group caller assignment/completeness dictionaries and new state-specific
  chemistry templates: flexible for raw PDB preparation, but duplicate preparation
  logic and expand the possibilities for inventing or overwriting an assignment.
- Accept sanitized graphs without declared provenance: simpler but confuses
  parser-imputed H with intentional source assignments; rejected.
- Require declared prepared chemistry, validate it, and preserve it directly:
  selected by the implementer under the user's implement-suggestions instruction.
  It makes ambiguity explicit and keeps this library from becoming a preparation
  engine. Users needing raw PDBs must prepare a typed molecular graph first.

## Acceptance

Known valid audited states must preserve chemistry/coordinates/identities in all
four H representations. Contradictory aliases, radicals, unsupported connectivity,
missing atoms and undeclared evidence must fail without mutating inputs. Tests
must include a mixed receptor, a representative named-H peptide with termini and
disulfide, actual implicit/ionic interactions, alternative resonance forms,
PDB/atom aliases, serialization, and equality of detector eligibility/geometry
before and after H removal. Existing tests, lint, mypy and build must pass, with
any pre-existing limitation reported. Changes are unreleased; version unchanged.
