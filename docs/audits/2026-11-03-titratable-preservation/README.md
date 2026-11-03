# Audit of source protonation preservation beyond histidine

**Audit date:** 2026-11-03. **Conclusion:** general source protonation preservation
is not supported. The audit reproduces charge changes, proton relocation,
introduced radicals, lost connectivity, and changed interaction detection outside
histidine. `preserve_histidine=True` protects supported histidines; it does not
protect these other residues.

This is an investigation with executable evidence, not a fix. The observations
describe the recorded baseline; subsequent fixes are separate.
The archived evidence is immutable and can be regenerated at the audited commit.

## Scope and evidence

Audited checkout: `b1094973e2d194171a460857de3587116676adad`, branch
`conductor/prolif-preserve-histidine`. Baseline:
`origin/main` at `e876dd4ddcd8df9f3c1a0d9f782107d132294a5f`.
Experiments use the locked environment, RDKit **2025.09.6**, and this checkout's
ProLIF **2.2.3+atomica.1**. The version label alone does not identify this unreleased
branch; use the commit above. Exact dependency versions and raw-result hashes are
in the generated [manifest](evidence/manifest.json).

- **560 chemistry cases:** 70 source-state/name combinations, four H
  representations, and both standardizer modes. Coverage includes Asp, Glu, Lys,
  Arg, Cys, Tyr, neutral disulfides, N termini including proline, C termini, and
  Ser/Thr controls. It includes matching names, generic names, contradictory
  aliases, alternative acid proton positions, and carboxylate/guanidinium
  resonance encodings.
- **104 interaction scenarios:** 13 selected fixtures, four representations,
  both modes. 86 reach interaction detection; 18 stop at the reproduced
  standardization errors. Each completed scenario uses 96 fixed ligand poses,
  four real detectors, and three stages: source, standardized residues, and
  H-removed/rebuilt residues. Geometry checks remain enabled.
- **42 alias/context probes:** 22 oxygen-name probes, 18 direct-molecule versus
  extended-PDB comparisons, and two parent-context checks.
- **28 edge probes:** conventional H names, ProLIF source mutation, ambiguous
  PDBs, force-field context, and contradictory H naming/attachment.
- **280 baseline comparisons:** every current legacy-mode matrix result matches
  the result from `origin/main`'s standardizer exactly after JSON normalization.
  Its shared constants, template engine, molecule and utility code are unchanged
  between these commits. The baseline class is loaded separately; no checkout or
  production-code mutation is needed.
- **499 existing tests passed**, with 40 dependency/fixture warnings, covering
  histidine preservation, IO, connected-parent behavior, molecules, residues,
  interactions and serialization. These passing tests do not cover the newly
  demonstrated non-histidine failures.

The main matrix uses chemically valid, initially radical-free peptides made with
`Chem.MolFromSequence`: an internal residue in Gly-X-Gly, a dipeptide for termini,
and a Cys-Gly-Cys peptide with an explicit S-S bond. Source states are constructed
before embedding with a fixed seed. The four representations are all graph H,
ordinary `RemoveHs`, `RemoveHs(updateExplicitCount=True)`, and the latter with
`NoImplicit=True` on every atom. An H-free graph with retained chemistry is
different from an unassigned, charge-less PDB; the latter is tested separately.

Measurements distinguish the connected molecule from its residue fragments.
They record site formal charges, total H **including graph neighbors**, radicals,
hybridization, named heavy-atom bonds, donor/acceptor/ionic SMARTS eligibility,
identities, heavy-atom order and coordinates. Fragment radicals at cut bonds are
not automatically classified as changed source protonation.

## Findings by chemical group

| Group | Passing controls | Reproduced failures |
| --- | --- | --- |
| Asp | Conventionally encoded ASP anion and ASH acid on OD2 | Generic ASP overwrites a neutral acid; contradictory ASH can neutralize an anion; acid proton on OD1 moves or errors; opposite carboxylate resonance errors. |
| Glu | Conventionally encoded GLU anion and GLH acid on OE2 | Same behavior as Asp, on OE1/OE2. |
| Lys | Ordinary cationic LYS | Neutral LYS, LYN and LSN are assigned NZ +1. With ordinary H representations, H2 becomes H3; with fixed counts a radical is introduced. |
| Arg | Cationic source with the template's NH2-positive encoding | Neutral ARG/ARGN becomes charged; an NH1-positive guanidinium resonance form becomes a dication; fixed counts can create radicals. |
| Cys | Neutral CYS in the unknown-force-field control; isolated thiolate under CYM/CYX | CYS thiolate is neutralized; thiol under a CYX-family alias is assigned -1; global force-field context changes selection; H/radical behavior depends on representation. |
| Disulfide | Connected S-S graph named CYS remains neutral in the tested unknown-force-field context | CYX/CYM/CYS2 assign -1 to both sulfur atoms while retaining S-S connectivity. |
| Tyr | Both neutral and anionic source OH survive with intact TYR metadata and charge | TYM has no built-in template; lost charge evidence in a PDB is not recovered or rejected. |
| N terminus | Neutral and +1 N, including proline, with canonical atom names and explicit source chemistry | Ambiguous PDB input is accepted as inferred neutral chemistry; four-letter prefix names can be truncated into unsupported names. |
| C terminus | Neutral and -1 OXT with canonical atom names and explicit source chemistry | Alternative oxygen names can lose their bonds; ambiguous PDB input is accepted as an acid; four-letter prefix names can be truncated. |
| Ser/Thr controls | Neutral OH and explicitly charged alkoxide states survive in these fixtures | This is a tested control, not a prediction that these states occur in a particular receptor. |

These passing controls are conditional observations, not guarantees for all
preparations, aliases, stereochemistry, crosslinks, or metal complexes.

### 1. Non-histidine source evidence does not control template selection

Only histidine gets evidence resolution before template selection and verification
after it. Other residues use name normalization and a fixed template. The charge
override table assigns ARG NH2 +1, LYS NZ +1, ASP OD2 -1, GLU OE2 -1, CYS SG 0,
and CYX SG -1; ASH/GLH assign zero at their corresponding oxygen.
See [standardization](../../../prolif/io/molecule_standardizer.py),
[name/charge tables](../../../prolif/io/constants.py), and
[template application](../../../prolif/io/template_engine.py).

For neutral lysine with an intact `LYN` name, all-graph-H and ordinary H-removed
inputs produce **NZ charge 0/H2 → charge +1/H3**. The extra proton is inferred;
the library need not add a new hydrogen atom to change the chemical state.
Conventional HZ naming reproduces the same result as the arbitrary unique H
names used in the matrix. Thus this is not a fixture-name artifact. The meaning
of LYN as neutral lysine with two NZ hydrogens agrees with the
[OpenMM preparation documentation](https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html#openmm.app.modeller.Modeller.addHydrogens).

Neutral ARG/ARGN similarly becomes +1. More seriously, a valid guanidinium
resonance encoding with **NH1 +1** is rebuilt with a double bond to NH2 and an
additional **NH2 +1**, while the original NH1 charge is retained. The group charge
becomes **+2**. This is not an equivalent resonance representation: total charge
changes. With implicit H allowed, NH1 acquires an additional H; fixed H counts
instead expose a radical.

With `NoImplicit=True`, neutral LYN becomes NZ +1/H2 with one radical electron;
neutral ARG can likewise become a charged radical. CYS thiolate can become a
neutral sulfur radical. All these sources started radical-free. Sanitization
success therefore does not establish preservation, or even an ordinary intended
amino-acid protonation state.

### 2. Acids expose both silent overwrites and errors

A source acid named ASP or GLU with ordinary H-free chemistry silently becomes
an anion, loses its hydroxyl donor, and gains an ionic contact. Supplying the
actual acid H as a graph atom, or retaining a fixed explicit H count, instead
causes an oxygen-valence error in these cases. That incidental valence failure
is not an actionable source-evidence conflict check.

ASH and GLH support the conventional neutral templates, but do not verify the
input. An anion carrying one of these names can be silently neutralized. For a
neutral ASH acid whose proton belongs to OD1 rather than OD2:

- Ordinary H-free input moves the proton to OD2 and reverses the local C-O
  bond orders while keeping atom names and coordinates.
- Graph H attached to OD1 under a non-template name causes a valence error.
- If that same graph H is named HD2 despite being attached to OD1, the template
  silently reattaches it to OD2. Neither mode rejects the contradictory
  name/connectivity evidence.

The neutral acid's proton placement is chemically distinct at fixed coordinates;
it must not be dismissed as an interchangeable carboxylate resonance drawing.
Separately, genuine anionic carboxylate resonance alternatives with the negative
charge on OD1/OE1 fail: the template puts a double bond on the still-negative
oxygen. A preservation implementation needs resonance-aware matching rather
than comparing one arbitrary bond/charge localization.

### 3. Cysteine and disulfides require distinct treatment

Generic CYS chooses neutral SG and can turn a source thiolate into a thiol.
CYM, CYN, CYX, CYS1 and CYS2 share the CYX template and its SG -1 override. They
cannot represent both free thiolate and neutral disulfide correctly through
that single rule.

For the neutral Cys-Gly-Cys disulfide named CYX, both SG atoms receive **-1**.
With ordinary H representations, the opt-in connected graph retains the S-S bond
and reports an inferred H on each SG. With fixed H counts, it instead has
SG -1/H0 with a radical. The exact defective product depends on representation;
the charge overwrite is common. Real ionic detection gains a contact to a
cationic ligand at a sulfur that was neutral in the source.

For the same disulfide named CYS in the unknown-force-field context, the
connected source charge/H state is preserved. A temporary residue-fragment
radical appears at the severed S-S endpoint; this is recorded separately from
the true CYX charge defect.

Selection is also **nonlocal**: adding a separate residue named DALA causes
GROMOS inference, changing an otherwise identical generic CYS from the neutral
template to CYX. Its local source evidence did not change. The relevant logic
is `forcefield_guesser`/`convert_to_standard_resname`, not a bond-aware cysteine
state resolver.

Five opt-in contradictory thiol/CYX-family graph-H cases also change total sulfur H
during subsequent `RemoveHs`: the already-corrupted standardized graph reports
H2, while the removed graph reports H0. These are distinct from removal of H
from the valid sources, which preserved every selected-site snapshot.

### 4. Parsing and atom names can destroy evidence or connectivity

The preserving PDB reader restores four-letter names only for histidine aliases.
The extended format is already supported by this project for those histidines;
it is not the standard three-character PDB residue-name field. Under both modes:

| Extended source name | Parsed name | Direct in-memory input versus PDB-path result |
| --- | --- | --- |
| ASPH / ASPP | ASP | Neutral acid becomes anion only on the PDB path. |
| GLUH / GLUP | GLU | Neutral acid becomes anion only on the PDB path. |
| CYS1 | CYS | Isolated thiolate becomes neutral on the PDB path. |
| LYSH | LYS | Name is lost; this cation control selects the same final state. |
| ARGN | ARG | Name is lost; both paths already apply the charged ARG template. |
| NALA / CALA | NAL / CAL | PDB path raises a missing-template error. |

Initial parsing retained the tested atom/bond counts and source site charges.
The first three rows therefore isolate name loss followed by the wrong template,
rather than attributing all charge changes directly to the PDB parser.

Nine oxygen aliases declared in the repository's Protein atom-name namespace
silently lose their original C-O bond:

- OXT aliases: **OT, OT2, O2, OCT2, OC2**.
- Carbonyl O aliases: **O1, OT1, OCT1, OC1**.

All nine become disconnected oxygen atoms with inferred H2 in the neutral
alanine probes. Canonical O and OXT controls keep their bonds. Both modes damage
the residue; the opt-in mode also damages the returned connected molecule.
There is no warning in these probes. Source coordinates and heavy-atom order
still match exactly, showing why those checks alone are insufficient.

The engine strips source bonds, applies only residue-specific atom aliases, and
skips template bonds whose atom names cannot be found. It does not merge the
shared Protein aliases. Some of those names are recognized separately for the
heavy-atom-count exemption, which does not restore connectivity. The independent
[static review](evidence/static-review.md) records the
precise source locations.

### 5. Parent/residue consistency changes the practical effect of the option

Legacy standardization replaces residue chemistry but leaves the connected
molecule's original charges/bonds. For neutral LYN, the returned residue is +1
while its parent is still neutral. `Molecule(Chem.RemoveHs(result))` rebuilds
residues from that unchanged parent and exposes the neutral state again. This
is inconsistent representation, not protonation restoration by `RemoveHs`.

`preserve_histidine=True` copies standardized chemistry for **every residue** into
the connected graph. Thus LYN's unwanted +1 persists into the H-removed molecule.
The same propagation applies to the disulfide and oxygen-alias defects. This
is an integration-relevant consequence even for a receptor with no histidine.
It does not mean the histidine option originally introduced those templates.

The opt-in path also does not ensure every parent/residue property is identical.
The ARG NH1-positive resonance probe ends with parent NH1 hybridization SP3 but
returned residue NH1 SP2. After rewrapping, the residue becomes SP3. Its donor
pose count changes from **15 to 39** despite unchanged standardized q/H/rad.
The CIF engine assigns hybridization from source geometry after sanitizing the
fragment; the connected parent is later sanitized separately. Geometry metadata
must be included in interaction-equivalence checks.

The source-mutation probe agrees with the documented contract: legacy modifies
the supplied ProLIF object's names/residues in place; opt-in returns a separate
object and leaves that source unchanged. This is a targeted success-path check,
not exhaustive proof of every metadata field or failure-path atomicity.

## Actual interaction consequences

The following counts use the ordinary H-removed input and opt-in mode. Each row
uses identical receptor coordinates and the same **96** probe directions at each
stage. The ligand contact atom is 2.9 Å away; methanol probes implicit H bonds,
formate probes protein cations, and methylammonium probes protein anions. These
are synthetic contact tests, not experimental protonation assignments or docking
pose statistics. No geometry checks were disabled.

| Source fixture/site | Detector outcome | Source → standardized → H-removed |
| --- | --- | --- |
| Neutral LYN NZ | Protein accepts H bond | **92 → 0 → 0** |
| Neutral LYN NZ | Protein cation ionic contact | **0 → 96 → 96** |
| Neutral ARGN NH2 | Protein cation ionic contact | **0 → 96 → 96** |
| Neutral acid named ASP, OD2 | Protein donates H bond | **13 → 0 → 0** |
| Neutral acid named ASP, OD2 | Protein anion ionic contact | **0 → 96 → 96** |
| Anion named ASH, OD2 | Protein anion ionic contact | **96 → 0 → 0** |
| CYS thiolate, SG | Protein anion ionic contact | **96 → 0 → 0** |
| CYS thiolate, SG | Protein donates H bond | **0 → 39 → 39** |
| Neutral disulfide named CYX, SG | Protein anion ionic contact | **0 → 96 → 96** |
| TYR phenolate, OH | Protein anion ionic contact | **96 → 96 → 96** |
| Cationic N terminus, N | Protein cation ionic contact | **96 → 96 → 96** |
| Anionic C terminus, OXT | Protein anion ionic contact | **96 → 96 → 96** |

Fixed H counts can produce different detector errors. For example, neutral LYN
then loses donor contacts (38 → 0) while retaining acceptor matches and gaining
ionic contacts, because the output is a cationic radical rather than ordinary
protonated lysine. A single representation is not enough to qualify the path.
The [contact JSON](evidence/contact-results.json.gz) records
exact pose indices, coordinates, metadata examples and all tested variants.

## What the passing results establish

The entire selected-site snapshot of each valid source—charge, total H, radicals,
heavy bonds, hybridization and SMARTS eligibility—survived direct `RemoveHs` in
all **560** cases. This supports retaining implicit detection: the dominant
failures are earlier template/evidence changes, not an inherent need to choose
protonation when removing hydrogen graph atoms. It does not establish that every
possible molecule or H-removal configuration is safe; RDKit exposes several
[removal options](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RemoveHs).

For each mode, the matrix has **140 unchanged selected-site q/H/rad signatures,
98 changed signatures, and 42 standardization errors**. These are deliberately
constructed stress/control cases, not failure rates in real receptors. Four of
the changed signatures per mode are the neutral CYS disulfide's fragment radical
bookkeeping. Matching signatures are not full chemical preservation: they do not
include all other atoms, bonds, geometry metadata, or source identity.

Every successful main-matrix case preserved heavy-atom order and coordinates.
Source RDKit heavy-atom charge/H/rad signatures were unchanged in successful
calls. Names can still be normalized, chemistry can change at the same coordinates,
and the alias probes demonstrate bond loss despite unchanged positions.

Tyr and the termini have no forced charge override at their relevant sites, so
intact source charges survive the tested canonical inputs. When graph H and PDB
charge fields are deliberately removed, the source evidence is lost: phenolate
becomes inferred phenol, a cationic N terminus becomes inferred neutral N, and an
anionic C terminus becomes inferred acid. Both modes accept those ambiguous
non-histidine inputs. This is a parser/evidence limitation, not proof that a
known original state can be reconstructed from missing information.

## Recommended follow-up and limits

The evidence supports keeping implicit detection while extending source-state
resolution and verification before using a broad preservation claim. Priority:

1. Resolve reliable per-residue evidence for Asp/Glu/Lys/Arg/Cys/Tyr and termini;
   require preparation or caller evidence when ambiguous, and reject conflicts.
   Preserve observed H attachment and charge rather than treating names as an
   unconditional override. Give neutral Lys/Arg and disulfides distinct chemistry.
2. Verify source → standardized connected graph and residue views. Include group
   charge, proton position, radicals, heavy connectivity, atom mapping/identity,
   coordinates and resonance equivalence; avoid leaving source charges behind
   when canonicalizing a resonance form.
3. Repair shared atom-alias application and retain supported non-histidine
   extended names at library-owned parsing boundaries. Missing source evidence
   from external parsers needs explicit caller input, not guessing.
4. Check standardized → H-removed graph and actual detector behavior, including
   hybridization/geometry context. Add mixed-receptor and prepared-receptor
   integration fixtures before changing the consumer's preservation claims.

The canonical Tyr/terminus controls do not qualify every terminal atom name,
caps, custom template, modified residue, covalent adduct, metal coordination,
missing atom, alternate location, or force field. The audit does not infer
experimental states, perform pH prediction, or measure prevalence in Atomica's
receptor collection. It does not exercise the separate Atomica application or
its private preparation/removal logic. No claim is made about upstream versions
other than the recorded baseline or RDKit versions other than the installed one.

## Reproduction and artifacts

From the repository root with the existing locked environment installed:

```bash
bash docs/audits/2026-11-03-titratable-preservation/run_audit.sh
MISE_UV_VERSION=0.11.24 uv run --no-sync pytest -q \
  tests/test_histidine_preservation.py tests/test_io.py \
  tests/test_connected_parent.py tests/test_molecule.py tests/test_residues.py \
  tests/test_interactions.py tests/test_pickling.py
```

`run_audit.sh` rebuilds the matrices, compares the baseline, verifies the recorded
findings and regenerates [matrix.csv](matrix.csv). `validate.py` asserts the
observed defects as historical reproductions; these are not normal product
regressions asserting desired behavior, and will need updating after fixes.
The scripts do not change library files. They write raw evidence under the
gitignored `.context/titratable-audit/` directory and regenerate the local CSV.

- [Chemistry fixtures and matrix](audit.py)
- [Real detector probes](contacts.py)
- [Parsing/alias probes](alias_probes.py)
- [Evidence and mutation edge probes](edge_probes.py)
- [Baseline comparator](baseline.py) and [observation validator](validate.py)
- [JSON evidence validation helpers](evidence.py)
- [Full chemistry JSON](evidence/results.json.gz),
  [alias JSON](evidence/alias-probe-results.json.gz),
  [edge JSON](evidence/edge-results.json.gz), and
  [existing-suite output](evidence/pytest.log.gz)

The audit scripts pass repository Ruff lint, formatting and mypy checks. An independent
review checked the source-state construction, alias probes, state signatures,
resonance distinctions and boundary interpretation; its
[review notes](evidence/harness-review.md) distinguish
measured results from untested guarantees.

Historical evidence is checked in under `evidence/`; JSON files use deterministic
gzip compression. Their uncompressed SHA-256 hashes match the archived manifest.
The reproduction assertions describe the audited commit, not expected behavior
after later fixes. To reproduce historical results, use the recorded checkout.
