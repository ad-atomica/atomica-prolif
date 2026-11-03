---
orphan: true
---

# Histidine Preservation Implementation Plan

> Execute inline using superpowers:executing-plans and test-driven-development;
> request one independent whole-change review before delivery.

**Goal:** Preserve reliable source histidine assignments through standardization
and supported removal of hydrogen graph atoms.

**Architecture:** Add an opt-in copy path to MoleculeStandardizer. Resolve and
validate source evidence in a focused internal histidine module, select existing
templates per residue, and validate output against the original evidence. Keep
connected parent chemistry consistent with standardized residues.

**Tech stack:** Python >=3.10, RDKit >=2025.9.6, gemmi, uv, pytest.

**Spec:** [spec.md](spec.md).

## Global constraints

- Work in this Conductor checkout on `conductor/prolif-preserve-histidine`.
- Preserve legacy defaults; evidence arguments require the preservation option.
- No protonation predictor, all-residue preservation claim, or Atomica changes.
- Keep package version unchanged: this is an unreleased feature PR.

## Review focus

- Missing or partial graph H must never silently select HIP for generic HIS.
- Name/assignment/completeness evidence must agree with actual connectivity.
- Source identities, atom order, mapping and coordinates must survive copying.
- Connected parents and H-removed graphs must expose the resolved chemistry.
- Four-letter aliases must survive library-owned PDB parsing.

## Task 1: Source preservation and interaction regressions

**Files:** create `prolif/io/_histidine.py`, `tests/test_histidine_preservation.py`;
modify `prolif/io/molecule_standardizer.py`.

**Interfaces:** `MoleculeStandardizer(..., preserve_histidine=False)` and
`__call__(..., histidine_states=None, complete_histidine_hydrogens=())`.

- [x] Construct realistic small histidines and a mixed HIS receptor; demonstrate
  the original neutral-to-cation overwrite with a failing chemistry regression.
- [x] Add failing evidence, identity/copy, PDB parsing, and detector regressions.
- [x] Resolve source assignments before templates and reject unreliable evidence.
- [x] Apply existing templates without leaking template names into histidine keys.
- [x] Validate output state/correspondence and synchronize connected chemistry.
- [x] Verify actual implicit contacts and charge eligibility before/after RemoveHs.
- [ ] Run new tests and existing IO, connected-parent, molecule, residue,
  interaction, and serialization tests; commit the focused implementation.

## Task 2: Documentation and delivery

**Files:** `ATOMICA.md`, `docs/source/modules/utils.rst`, this plan.

- [x] Document API defaults, accepted/rejected evidence, mutation/name contracts,
  supported H removal and four-letter parsing, and histidine-only limits.
- [x] Run full pytest, relevant Ruff and mypy checks, and `uv build`.
- [x] Request independent review; address important findings with regressions.
- [ ] Record actual verification, commit, push, and open a draft PR against main.

## Execution notes

- The user supplied a detailed implementation handoff and authorized routine API
  decisions; proceed without another design approval or another workspace.
- The uv shim has no default version. Use installed uv 0.11.24 via
  `MISE_UV_VERSION=0.11.24`. Initial automatic Python 3.14 selection could not
  build locked Pillow (missing JPEG headers); use the supported Python 3.11.

## Verification and review

- Baseline relevant suites: 431 passed. Direct pre-fix reproduction showed HIE
  encoded as HIS change from `(ND1 H=0, NE2 H=1, ring charge=0)` to `(1, 1, +1)`.
  The new API regressions failed before implementation.
- Final focused source suite: 499 passed (68 new tests plus the existing IO,
  connected-parent, molecule, residue, interaction and serialization suites).
- Full source suite: 661 passed, 51 dependency/fixture warnings in 101.72 s,
  including the final rerun after the annotation-only API correction.
- Built wheel installed under an isolated target and tested from `/tmp`: all
  68 new regressions pass. Wheel and source archive build with `uv build`.
- Repository-wide Ruff lint and formatting pass; `git diff --check` passes.
- Scoped mypy (`--follow-imports=silent`) passes the three changed Python files
  and a consumer example using ordinary string-keyed and ResidueId-keyed dicts.
  The regular mypy invocation reports an unchanged `no-any-return` error at
  `prolif/interactions/utils.py:35`; reproduced using origin/main's identical
  helper. No unrelated type-check workaround or interaction change was added.
- Sphinx HTML builds with notebook execution disabled; warnings remain for
  existing docs, duplicate autodoc targets and unlisted legacy notebooks.
- Independent review found two important acceptance gaps: malformed ring-carbon
  connectivity and aromatic ring radicals. Both have failing-then-passing
  regressions, including an external ring-N bond check against the connected
  parent and a custom-template heavy-bond preservation check. Ring validation
  now rejects these unsupported species before template application, and the
  same validation runs afterward. Existing resonance coverage remains green.
- Also replaced the per-residue full-parent bond scan with incident-bond
  traversal to avoid quadratic work. No review findings remain deferred.
- Ruling: keep the package version unchanged because this is an unreleased
  feature PR; release/version/pin updates belong to subsequent authorized work.
- Ruling: use the connected source graph as the preservation boundary and return
  a separate standardized parent/residue graph. Source histidine names remain
  public identities; other residue names retain existing normalization.
- Ruling: preserve four-letter histidine PDB aliases from raw fixed columns by
  source serial number. Unsupported serial encodings require intact caller
  metadata rather than guessing a lost state.
