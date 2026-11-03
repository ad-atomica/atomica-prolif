# Independent audit harness review

Reviewed `audit.py`, `contacts.py`, `edge_probes.py`, `baseline.py`, their JSON
results, and the relevant production methods. No root-owned harness or production
files were changed. Review includes the final four-representation matrix:
560 cases = 70 source specifications × 4 H representations × 2 modes.

## Assessment

The central findings are supported. I found no fixture or measurement error that
invalidates the demonstrated charge/H changes, PDB/atom alias defects, or detector
differences. The major reporting risk is equating selected-site q/H/radical
equality with full chemical or interaction preservation. A concrete remaining
parent/residue hybridization mismatch demonstrates why that distinction matters.

## Independently verified

- Final `results.json` contains 560 rows. Each mode has 140 unchanged selected-site
  q/H/radical signatures, 98 changed signatures, and 42 errors. The three-way
  counts are descriptive outcomes, not 560 passing regression tests.
- Parsed `baseline-results.json` has 280 rows and equals all 280 current legacy
  rows exactly after JSON normalization. The baseline class is extracted from
  `origin/main`; its supporting production modules (`constants`, template engine,
  molecule, residue, utils) have no diff from that ref in this workspace. This
  establishes equivalence for these fixtures rather than universal equivalence.
- All 70 source fixtures were independently rebuilt. They contained no radicals
  anywhere. Acid sources had the requested side-chain charge/H totals (neutral:
  0/1; anion: -1/0). Arginine sources had the requested guanidine charge/H totals
  (neutral: 0/4; cation: +1/5). Each source disulfide SG had charge 0, H 0 and two
  heavy neighbors, one the other SG. The fourth representation additionally
  asserts radical-free source chemistry after fixing explicit H counts.
- Every full selected-site `source` snapshot equals its `remove_h_source`
  snapshot in all 560 rows, including heavy bond orders, hybridization and SMARTS
  eligibility. Graph H atoms are deliberately removed; this is preservation of
  the recorded chemical state, not literal graph identity.
- Source-neutral LYN NZ is 0/H2. In the three representations permitting implicit
  H it becomes +1/H3 in standardized residues; opt-in parents carry that same
  charge/H change. With `fixed_h_counts`, NZ becomes +1/H2/radical 1 instead.
- NH1-charged ARG cation1 starts with side-chain total charge +1. The template
  retains NH1 +1 while assigning NH2 +1, giving side-chain +2. With implicit H
  allowed NH1 becomes H3; with fixed H counts NH1 remains H2 with a radical.
  This fixture is a valid initial charge/bond-order resonance assignment, not an
  accidental +2 source.
- Both SG atoms of the opt-in CYX disulfide were checked directly, beyond the
  main matrix's first-SG snapshot. With implicit H allowed each has charge -1,
  total H 1 and the original S-S bond. With fixed H counts the recorded first SG
  instead has charge -1, H 0 and radical 1 in the connected parent. The fragment
  has charge -1/H0/radical 0. Do not generalize -1/H1 across all representations.
- Generic CYS disulfide parent remains neutral/H0; its standardized residue gains
  a radical from cutting the cross-residue bond. That belongs to fragment
  bookkeeping, not a source protonation-state change. Matrix signature-change
  totals include this case, so those totals must not be called protonation-change
  counts without qualification.
- Canonical *target-site* H-name controls in all 12 edge rows produce the same
  recorded outcomes as the arbitrary Q-name controls. Other H names remain Q
  names, so this is not an exhaustive conventional PDB-H-name audit.

## Important additional representation difference

For opt-in ARG cation1/ARG with `no_graph_h`, NH1 states are:

| Stage | Charge | Total H | Radical electrons | Hybridization |
| --- | --- | --- | --- | --- |
| Source | +1 | 2 | 0 | SP2 |
| Standardized connected parent | +1 | 3 | 0 | SP3 |
| Standardized returned residue | +1 | 3 | 0 | SP2 |
| H-removed/rebuilt residue | +1 | 3 | 0 | SP3 |

The chemistry signatures compare equal between standardized and rebuilt
residues, but hybridization does not. The contact harness records 15/96 donor
poses at standardized NH1 and 39/96 after rebuilding. This is correctly measured
and should not be misattributed to an H-removal protonation change.

Production mechanism: `CIFTemplateEngine.apply` sanitizes via
`assign_intra_props`, then assigns source-geometry hybridization
(`prolif/io/template_engine.py:114-124`). Opt-in later sanitizes the connected
parent (`prolif/io/molecule_standardizer.py:339-346`) but returns the previously
built residue objects. Several CYS/disulfide fixtures also show parent/residue
hybridization differences. Therefore opt-in parent/residue agreement in charge
and H does not establish complete agreement in detector inputs.

## Detector and metadata review

`audit.py:93` now correctly gets cation and anion SMARTS from
`Cationic().lig_pattern` and `.prot_pattern`. Inverted interaction classes swap
arguments during detection; their stored pattern attributes are not swapped.

`contacts.py:55` uses the correct detector roles: ImplicitHBAcceptor for a protein
donor, ImplicitHBDonor for a protein acceptor, Anionic for a protein cation, and
Cationic for a protein anion. It filters returned metadata to the intended
protein site and ligand anchor rather than counting neighboring matches.

The contact matrix has 104 scenarios, with 18 standardization errors represented
as error rows. Each successful row uses identical ligand poses across source,
standardized and reconstructed stages. Implicit H-bond geometry checks remain
enabled. Charged-pair detection is distance based. The synthetic methanol anchor
geometry supports the intended comparisons; the charged ligand geometry is
sufficient for a distance detector.

Counts such as 92/96 acceptor poses or 96/96 ionic poses describe this designed
pose sample at 2.9 A; they are not probabilities, empirical interaction rates,
or a benchmark of biological specificity. No energy relaxation is claimed.
Static SMARTS eligibility and real detector hits should remain distinct.

## Scope and optional follow-ups

- `state_sig` is explicitly q/H/radical at selected sites (`audit.py:102-103`).
  It excludes bond orders, hybridization and neighboring states, even though
  selected-site bonds and hybridization are separately present in JSON. In the
  current rows classified unchanged, I found no selected-site heavy-bond,
  hybridization or SMARTS-eligibility differences from source; this observation
  still is not a whole-molecule preservation proof.
- Heavy-coordinate/order checks ignore graph-H coordinates and hydrogen identity.
  Source-copy checks in the main matrix compare heavy-atom q/H/radical values and
  run after successful standardization. They do not prove every source property
  remains unchanged on failure. The edge ProLIF-input control separately shows
  legacy in-place mutation and opt-in copying for one neutral-Lys fixture.
- Total H includes graph H, explicit atom counts and implicit H. An increase
  from H2 to H3 need not add a hydrogen graph atom or a new coordinate. Use
  “inferred total H” rather than asserting a new positioned proton was added.
- Intentionally contradictory residue-name/state combinations are useful
  adversarial cases, but distinguish them from accepted same-state aliases and
  source-generic naming. They test whether trusted chemistry gets overridden;
  they do not imply every such naming convention is chemically self-consistent.
- Charge-stripped H-free PDB edge fixtures are deliberately ambiguous. A changed
  state already present in their `parsed` snapshot is parser/evidence loss,
  whereas a later change belongs to standardization. The harness records both.
- The arbitrary H names test nearest-heavy reconstruction. The targeted canonical
  controls substantially reduce the risk that central findings are artifacts of
  those names. Full named-H inventories, atom-name collisions, arbitrary custom
  templates and modified/metal-bound residues remain outside the finite matrix.
- Useful optional future controls would assert whole-molecule heavy bond/H
  attachment fingerprints before and after source copying (including failures),
  and add an explicit measurement of parent/residue hybridization agreement.
  These are not required to substantiate the present counterexamples. The newly
  added fixed-H/NoImplicit representation closes the most important earlier gap.

The independent terminal/PDB alias script and JSON support the alias claims
without relying on root's matrix construction. No further controls are required
before reporting the concrete observed failures with the qualifications above.
