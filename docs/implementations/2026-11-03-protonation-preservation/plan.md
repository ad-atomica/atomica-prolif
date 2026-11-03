---
orphan: true
---

# Broader source protonation preservation implementation plan

Execute under the user's implement-suggestions authorization. Existing workspace
is the requested isolated Conductor checkout; no extra implementation worktree.
Use TDD for fixes and an independent whole-change review. Spec: [spec.md](spec.md).

- [x] Archive the existing audit JSON/reviews and push the existing audit commits.
- [x] Add failing oxygen/residue alias and collision regressions. Fix shared
  Protein alias application and recognized extended PDB names; run IO/histidine
  regressions; commit independently.
- [x] Add failing parent/residue/H-removal consistency regressions; derive
  preserving-mode residues from the sanitized connected graph. Run connected,
  histidine, residue and interaction tests; commit independently.
- [x] Add strict prepared-graph API/evidence and source-copy failure tests. Write
  and verify regressions for valid neutral/cation/anion states and conflict cases.
  Implement topology and state validation without destructive templating.
- [x] Add mixed prepared receptor, H inventory/alias, resonance and real detector
  regressions. Verify both graph-H and H-removed forms and source immutability.
- [x] Document API, provenance contract and supported scope in API docs/ATOMICA.
- [x] Run full pytest, Ruff, mypy, package build, installed-wheel regressions and
  independent review. Address substantive findings with reproductions.
- [x] Commit/push coherent changes and open draft PRs in a reviewable stack.
  Leave the top branch checked out for Conductor; do not merge or release.

Completed implementation: `d9b4ee3`. Published draft chain:
[alias repair #2](https://github.com/ad-atomica/atomica-prolif/pull/2),
[parent/residue consistency #3](https://github.com/ad-atomica/atomica-prolif/pull/3),
[prepared preservation #4](https://github.com/ad-atomica/atomica-prolif/pull/4),
following the existing histidine/audit PR #1. See [verification.md](verification.md).

Review focus: source H versus inferred H; inter-residue bonds and termini;
explicit alias conflicts; equivalent resonance; no mutation on failure; residue
fragment chemistry versus connected-parent chemistry; unsupported species must
fail closed. Archive original evidence without rerunning historical expectations
against changed implementations.
