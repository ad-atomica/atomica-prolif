# Atomica ProLIF implementation plan

> Execute inline using the executing-plans, test-driven-development,
> atomica-test-design, and verification-before-completion skills.

**Goal:** Package Atomica's ProLIF fixes in an independently installable fork.

**Architecture:** Start at upstream v2.2.3, apply the exact existing two-file
patch, retain `prolif` imports and distribution metadata, and publish versioned
GitHub Release assets. Atomica integration is explicitly deferred.

**Tech stack:** Python >=3.10, RDKit, ProLIF, setuptools, uv, GitHub Actions.

**Spec:** [Fork specification](spec.md).

## Global constraints

- Base: `11a90d36ad0167aa6322bca250e61d65157e82cf` (upstream `v2.2.3`).
- Distribution: `prolif`; initial version: `2.2.3+atomica.1`.
- RDKit runtime floor: `2025.9.6`; keep the baseline lock on this qualified version.
- Preserve the exact patched chemistry implementation already used by Atomica.
- Do not modify the Atomica consumer checkout.

## Review focus

Exercise fragmentation across an amide bond, donor-plane rejection, configured
acceptor patterns, large-parent matching limits, shared-residue rebinding,
pickling, garbage collection, and imports from an installed wheel. Existing
upstream molecule, residue, interaction, and pickling tests cover compatibility.

## Task 1: Move and verify the chemistry fixes

**Files:** `prolif/molecule.py`, `prolif/interactions/interactions.py`,
`tests/test_connected_parent.py`.

- [x] Write standalone module-contract tests through `Molecule`
  and the implicit detectors. Use small synthetic molecules with assigned
  residue boundaries and coordinates; do not depend on Atomica imports.
- [x] Run against pristine upstream and observe the expected chemistry and
  geometry failures.
- [x] Apply Atomica's existing patch with `patch -p1` in `prolif/`.
- [x] Verify both source digests match Atomica's qualified patched digests.
- [x] Run the new tests and upstream molecule/residue/interaction/pickling tests.

## Task 2: Build independently installable artifacts

**Files:** `pyproject.toml`, `uv.lock`, `ATOMICA.md`, `.github/workflows/build.yml`,
`.github/workflows/_build-package.yml`, `.github/workflows/ci.yml`.

- [x] Set version `2.2.3+atomica.1`; update fork source and issue URLs while
  retaining upstream authorship, license, and documentation links.
- [x] Replace inherited PyPI publishing with Atomica-tag GitHub Release assets,
  validate tag/metadata agreement, and include SHA256 checksums.
- [x] Keep existing build and CI entry points working with the committed version.
- [x] Refresh `uv.lock`; build wheel and source archive using `uv build`.
- [x] Install the wheel into an isolated target, prove imports resolve there,
  and run the focused chemistry regressions against that artifact.
- [x] Build a wheel from the source archive and verify its package identity.
- [x] Document exact future consumer pins, updates, and the publication steps.
- [x] Commit the complete fork changes and export a git bundle for handoff.

## Task 3: Publish the fork when permitted

- [x] Create `ad-atomica/atomica-prolif` as a GitHub fork of upstream.
- [ ] Push the prepared `main` branch and set it as the fork's default branch.
- [ ] Push release tag `v2.2.3+atomica.1` and verify GitHub Release artifacts.
- [ ] Report exact outcomes and any remaining access blocker.

## Verification and deviations

- Pristine upstream: five expected regression failures and six passes.
- Patched source: 396 tests passed in 13.65 seconds, including all 11 new
  regressions and the upstream molecule/residue/interaction/pickling suites.
- Installed wheel outside the source checkout: all 11 regressions passed in
  2.23 seconds, with version and import origin checked explicitly. Final package
  metadata verification after the RDKit requirement change also passed; the
  final installed-wheel regression run passed all 11 tests in 2.39 seconds.
- Wheel, source archive, and wheel rebuilt from the source archive all built.
- Pinned Ruff 0.11.2 checks and formatting, actionlint, and `git diff --check`
  passed. Existing upstream/MDAnalysis deprecation warnings remain.
- Ruling: apply upstream Ruff formatting after importing and digest-checking
  the existing patch. An AST comparison with Atomica's qualified installation
  confirms that only formatting changed the patched sources.
- Review corrected the inherited RDKit compatibility claim: require the tested
  `2025.9.6` floor because older versions lack recursive-match control. Keep this
  floor in runtime requirements, optional/dev groups, the baseline lock, and CI.
- Review also moved the inherited type-check report from privileged
  `pull_request_target` comments to a read-only `pull_request` job summary,
  with checkout credentials disabled and SHA-based checkouts.
- The cloud integration could not create forks. Alex directed the Mac GitHub
  CLI, which successfully created the organization fork.
- Commit `68b21ba96990e729a93b4241f0ffb1659bb35bca` contains the tested fork.
  A standalone Git bundle and the wheel/source archive are saved beside this
  checkout in the parent Atomica workspace's ignored `.context/` directory.
- Publication remains blocked: the cloud credential also received HTTP 403 on
  push. Conductor's Mac publication command then returned `Timed out waiting
  for client`. A remote ref check confirmed that `main` was not created. A
  subsequent read-only Mac connectivity check also failed to return and was
  cancelled. No release tag was created or pushed, and no release is published.
- The global hk Git hook found the parent Atomica configuration and ran its
  checks against the wrong repository paths. Fork checks were run explicitly;
  `HK=0` disabled that unrelated global hook for commits in this separate repo.
- Atomica's tracked files and current branch remain unchanged.
