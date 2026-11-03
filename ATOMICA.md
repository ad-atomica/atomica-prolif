# Atomica ProLIF

This fork starts at upstream [ProLIF v2.2.3](https://github.com/chemosim-lab/ProLIF/tree/v2.2.3),
commit `11a90d36ad0167aa6322bca250e61d65157e82cf`, and preserves upstream's
Apache-2.0 license and authorship.

The distribution and Python imports remain `prolif`. The initial Atomica version
is `2.2.3+atomica.1`. Install this **instead of** the upstream distribution:
installing another distribution that also owns `prolif/` can leave overlapping
files. Packages are published to GitHub Releases; this fork does not publish to
the upstream PyPI project. PEP 440 local versions are not intended for PyPI.

## Changes

The existing Atomica patch changes `prolif/molecule.py` and
`prolif/interactions/interactions.py`:

- Residues share a plain copy of the connected molecule for geometry, preserving
  the bonds cut during fragmentation without adding interaction candidates.
- Implicit donor geometry uses the connected peptide plane and neighbor angles.
- Implicit acceptors at cut bonds recheck the configured SMARTS against the
  connected parent with unrestricted total and recursive match counts. This
  excludes backbone amide nitrogens that fragment-only matching can accept.
- Supplied residues are copied before rebinding their geometry context. The
  plain connected graph avoids retaining discarded molecule/residue wrappers.

Standalone residues and explicit hydrogen-bond detectors retain upstream
behavior. Focused regressions live in `tests/test_connected_parent.py`; upstream
molecule, residue, interaction, and pickling tests remain in place.

## Install a release

Once the release is published:

```sh
uv pip install 'prolif @ https://github.com/ad-atomica/atomica-prolif/releases/download/v2.2.3%2Batomica.1/prolif-2.2.3%2Batomica.1-py3-none-any.whl'
```

For Pixi's `[tool.pixi.pypi-dependencies]`:

```toml
prolif = { url = "https://github.com/ad-atomica/atomica-prolif/releases/download/v2.2.3%2Batomica.1/prolif-2.2.3%2Batomica.1-py3-none-any.whl" }
```

The consuming package should require `prolif==2.2.3+atomica.1`. Lock the URL and
hash using the consumer's resolver. An exact public pin `==2.2.3` also accepts
local versions under PEP 440, so use the full version to require this fork.
Check `SHA256SUMS` from the same release before installation. This fork requires
RDKit >=2025.9.6, the version qualified with Atomica's patch. Older upstream
RDKit versions lack the recursive-match control used by connected-parent
acceptor matching. Atomica's dependency and runtime patch setup have deliberately
not been changed by this fork-creation work.

## Develop and verify

```sh
uv sync --locked --only-group tests
uv run --no-sync pytest -q tests/test_connected_parent.py tests/test_molecule.py tests/test_residues.py tests/test_interactions.py tests/test_pickling.py
uv build
```

The release workflow installs the wheel into a fresh environment and runs the
connected-parent tests outside the source tree. The reusable upstream build
workflow also installs both the source archive and wheel and checks package data.

## Maintain and release

Keep `upstream` pointing to `https://github.com/chemosim-lab/ProLIF.git` and
`origin` pointing to `https://github.com/ad-atomica/atomica-prolif.git`. The fork's
`main` branch carries the Atomica patches; upstream tags remain intact.

1. Review and port changes from an upstream tag explicitly. Keep the
   connected-parent regressions green before removing or adjusting a patch.
2. Set a new version in `pyproject.toml`, such as `2.2.3+atomica.2`, and refresh
   the lock with `uv lock`. Build and run the focused suite.
3. Commit and push `main`, then create and push an annotated tag that exactly
   matches the package version: `git tag -a v2.2.3+atomica.1 -m 'ProLIF 2.2.3+atomica.1'`
   and `git push origin v2.2.3+atomica.1` for this initial release.
4. Enable GitHub Actions for the fork if GitHub disabled inherited workflows,
   and grant the release workflow permission to write repository contents.
   The Atomica release workflow publishes wheel, source archive, and checksums
   after installation and chemistry verification. No PyPI credentials are used.
5. Verify the remote assets before changing Atomica's pin or deleting its local
   patch installer. Consumer integration is a separate change.

Do not move a published release tag or replace its package assets with different
bytes. Increment the Atomica version for a corrected package.
