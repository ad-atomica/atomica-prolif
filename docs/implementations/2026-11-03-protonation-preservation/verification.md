---
orphan: true
---

# Verification of prepared protein preservation

Verified 2026-11-03 with Python 3.11.15, RDKit 2025.09.6 and the locked project
environment. This feature is unreleased; the package version is unchanged.
Strict implementation commit: `d9b4ee3`.

Published as draft PRs [#2](https://github.com/ad-atomica/atomica-prolif/pull/2)
(aliases), [#3](https://github.com/ad-atomica/atomica-prolif/pull/3) (consistency)
and [#4](https://github.com/ad-atomica/atomica-prolif/pull/4) (prepared preservation).

## Delivered changes

- Historical audit and essential raw evidence archived in commits `bdeeadc` and
  `8c0fd5c`. All five compressed result files match the SHA-256 values in the
  historical manifest after decompression. The historical matrix is unchanged.
- Alias parsing/connectivity repair: `2aab7e3`, including all nine oxygen aliases,
  extended residue names, canonical-name collision errors and custom-template
  precedence. Recorded test-first result: 32 failures became 37 passing tests;
  subsequent legacy-reader coverage is included in the current full suite.
- Final connected-parent chemistry supplies preserving residue views: `a2f0cd9`.
  This fixes stale hybridization independently of the strict API.
- The new `preserve_protonation=True` path validates and copies intentionally
  prepared graphs. It requires `source_is_prepared=True` and rejects conflicting
  aliases, unsupported topology, missing atoms, radicals and inconsistent named H.

## Automated checks

Commands use `MISE_UV_VERSION=0.11.24 uv run --no-sync` unless stated otherwise.

| Check | Result |
| --- | --- |
| `pytest -q` | 909 passed, 51 dependency/fixture warnings |
| `pytest -q tests/test_protonation_preservation.py tests/test_standardizer_aliases.py tests/test_standardizer_consistency.py` | 248 passed |
| `ruff check --preview .` | Passed |
| `ruff format --check .` | 67 files already formatted |
| `mypy .` | No issues in 66 source files |
| `uv build --out-dir .context/protonation-dist` | Wheel and sdist built |
| Installed wheel: strict, alias, consistency and histidine regression modules | 316 passed; import path checked to use the installed wheel |
| Isolated staged snapshot: same four regression modules | 316 passed; import path checked to use the detached verification worktree |
| Sphinx HTML build | Completed; existing documentation warnings remain |

The complete documentation build emitted 103 warnings, including pre-existing
notebook MIME output, duplicate changelog labels, docstring formatting and orphan
pages. Subsequent incremental builds reported 14. These are not a clean docs-lint
result; the changed API documentation rendered successfully.

The original API regressions were recorded failing before implementation. On
resuming this work, 10 additional named-H conflict tests failed: a single CIF
state omitted supported acid H names and the third terminal H. Supplementing
attachment validation made all 10 pass without moving source hydrogens.

## Independent review

A fresh reviewer inspected the strict implementation and the alias/consistency
commits, including independent chemical probes. One important finding was fixed:
the initial histidine validator checked ring totals but allowed charged ring
carbons or opposite charges on neutral-tautomer nitrogens. These passed RDKit
sanitization but changed ionic and donor/acceptor eligibility.

Sixteen regressions reproduced that gap across all four H representations and
then passed after per-atom charge validation. Eight positive controls retain
ordinary histidinium with +1 on either ring nitrogen. The final full-suite result
above includes these fixes.

The review's minor insertion-code scope gap was addressed in the API docs:
residue names/numbers/chains (or segments) must give unique ProLIF identifiers.
Insertion codes remain metadata; general insertion-code identity support is
outside this change. Existing collisions fail instead of merging valid residues.

Review scope decisions: caller preparation provenance remains an explicit caller
responsibility; modified residues/metals/caps and additional states are rejected;
unknown unique H names use their explicit connectivity; broader historical
histidine behavior and Atomica consumer integration remain separate work.

## Additional audit replay

The historical audit's 70 source-state/name fixtures were run in all four H
representations through the strict API: **212 preserved, 68 rejected conflicting
aliases, zero unexpected outcomes**. Accepted inputs retained whole-graph
chemistry, atom identity/order, coordinates and source immutability. Site charge,
H, radicals, hybridization and detector SMARTS eligibility agreed in connected
and residue views, before and after default hydrogen removal. Local replay code
and detailed results are under `.context/protonation-qualify-audit.py` and
`.context/protonation-qualified-audit.json`; permanent regression coverage is in
the test modules listed above.

An additional probe confirmed that unsanitized RDKit hydrogen removal without
updating H counts can discard graph-only H inventories. The documented supported
recipe keeps `Chem.RemoveHs`'s default sanitization; the replay also verified
`sanitize=False, updateExplicitCount=True` after initializing query ring metadata.
The library does not change RDKit's removal behavior.

These are controlled qualification fixtures, including mixed peptides, named H,
termini and disulfides. Atomica's own prepared receptor inputs still need consumer
qualification before adoption; no Atomica integration or release was performed.
