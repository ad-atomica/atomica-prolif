---
orphan: true
---

# Source-preserving histidine standardization

Preserve a source receptor's assigned histidine protonation and tautomer state
before template selection. This implements the task-only handoff from
[Atomica PR #1625](https://github.com/ad-atomica/atomica/pull/1625#discussion_r4447698611).
It does not predict protonation or audit other titratable residues.

`MoleculeStandardizer(preserve_histidine=True)` returns a new molecule and keeps
source histidine residue identities. The default remains the legacy in-place
ProLIF behavior, including name normalization and generic HIS's cationic template.
Other residues still use existing templates/name normalization in either mode.

Per-call `histidine_states` maps original residue keys (strings or `ResidueId`)
to HID/HIE/HIP. `complete_histidine_hydrogens` lists original residue keys for
which the caller guarantees the graph contains the complete ring-N hydrogen
inventory. This does not require all protein hydrogens. Completeness is never
inferred from PDB format, inferred valence, or one observed ring hydrogen.

Resolve before any source name changes. Explicit state aliases are evidence;
generic HIS is not. Check every observed ring hydrogen's actual attachment,
known hydrogen-name aliases, ring topology, and agreement between all stated
authorities. Missing, conflicting, duplicate, or unsupported evidence raises an
actionable residue-specific ValueError. Partial hydrogens can corroborate an
explicit assignment but cannot establish one without completeness. Reuse the
existing HID/HIE/HIS chemistry templates on a per-residue basis without changing
the public histidine keys. Verify hydrogen counts, aromatic ring chemistry and
net ring charge (allow resonance), atom correspondence, and coordinates afterward.

The preservation result must have consistent connected-parent and residue
chemistry so removing graph H with RDKit `RemoveHs` retains assigned chemistry
and implicit detector eligibility. Retain the source heavy-atom order and all
source coordinates; never mutate a source ProLIF or RDKit molecule, on success
or failure. Restore four-letter histidine aliases when this library reads PDB;
externally discarded names require explicit caller evidence.

No Atomica checkout/dependency changes, release tags, publishing, or source
review replies. Commit, push the existing feature branch, and create a draft
PR against this fork's main.
