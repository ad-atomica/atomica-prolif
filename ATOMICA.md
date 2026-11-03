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

### Source-preserving histidine standardization (unreleased)

The template paths now apply shared protein atom aliases (including terminal
oxygen names) and reject duplicate canonical names. Library-owned PDB reading
retains recognized four-letter residue aliases in both modes; unsupported extended
names, ambiguous serial mappings and conflicting atom-name mappings raise errors.
This fixes lost C–O connectivity and non-histidine template selection at parsing.

Histidine preservation now derives residue views from the finalized connected
graph. This keeps hybridization and H inventories consistent when the analysis
molecule is rebuilt after H removal. It does not extend the histidine-only
option's source-state guarantee to other residues.

`MoleculeStandardizer(preserve_histidine=True)` resolves histidine states from
source evidence **before** selecting per-residue HID/HIE/HIS templates. It
preserves a mixture of HID, HIE and HIP even when a writer calls them all HIS.
This fixes an existing template-selection limitation: implicit H removal does
not itself require choosing a new protonation state.

Explicit state aliases (HID/HIE/HIP, HSD/HSE/HSP and built-in forcefield aliases)
are evidence. Generic HIS requires per-call `histidine_states` assignments or a
`complete_histidine_hydrogens` declaration for each such residue. The latter
asserts that ring-N graph hydrogens are complete, not that all receptor H are
present. PDB format and a single observed H cannot establish completeness.
Names, caller assignments and validated H attachment/inventory must agree;
ambiguous or conflicting inputs, unsupported ring connectivity and ring radicals
raise actionable errors. No pH-based prediction
or global HIS override is performed.

The opt-in path returns a copy, retains original histidine residue identities,
atom correspondence and coordinates, and verifies state after template
application. Its connected graph and residues share the resolved chemistry;
`Molecule(Chem.RemoveHs(result))` can be used for implicit detection while
retaining the assigned state. Supported four-letter aliases are retained when
the standardizer reads an extended PDB (first model, unique decimal serials).
Assignments lost by an external parser require explicit caller evidence.

Default calls retain legacy behavior, including in-place mutation of ProLIF
inputs, name normalization, and generic HIS's cationic template. Other residues
still use their existing templates and name normalization in either mode. This
change does not guarantee preservation for all titratable residues or modified/
metal-coordinated histidines with unsupported connectivity. See the
[API documentation](docs/source/modules/utils.rst) and
`tests/test_histidine_preservation.py` for examples and regressions.

This is an unreleased library change; package version and existing release
assets remain unchanged. Atomica consumer integration must supply trustworthy
source evidence and enable this option before making its H-removed analysis
copy; updating its dependency is a separate PR.

## Install a release

Install the released wheel:

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
uv run --no-sync pytest -q tests/test_histidine_preservation.py tests/test_io.py tests/test_connected_parent.py tests/test_molecule.py tests/test_residues.py tests/test_interactions.py tests/test_pickling.py
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
