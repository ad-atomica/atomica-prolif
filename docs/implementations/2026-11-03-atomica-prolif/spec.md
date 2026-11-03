# Atomica ProLIF fork

Prepare `ad-atomica/atomica-prolif` as a fork of `chemosim-lab/ProLIF` based on
upstream `v2.2.3` (`11a90d36ad0167aa6322bca250e61d65157e82cf`). The fork owns
Atomica's connected-parent implicit hydrogen-bond chemistry and geometry fixes.
Alex requested the fork and packaging, and explicitly deferred all changes to
the Atomica consumer repository.

## Approach

Keep the `prolif` distribution and import names. Use version
`2.2.3+atomica.1` to distinguish the fork in dependency resolution and analysis
metadata. This avoids installing two distributions that both own `prolif/`.
Publish wheels and source archives as GitHub Release assets; PEP 440 local
versions are unsuitable for PyPI publication. A Git dependency is a usable
alternative but requires a build during consumer installation. Renaming only
the distribution would risk overlapping installations and change metadata
lookups; renaming the import would expand this into an API migration.

These packaging choices are the implementer's proposal. The repository name,
organization, patch migration, and deferred Atomica integration are Alex's
requirements.

## Required behavior

- Preserve upstream history, authorship, Apache-2.0 license, and documentation.
- Apply the existing Atomica `prolif.patch`, verify the imported source digests,
  and preserve its behavior. Formatting may follow upstream's pinned Ruff style.
- Preserve connected-parent donor geometry, amide acceptor chemistry, residue
  indices, standalone residues, custom SMARTS, and unrestricted parent matches.
- Preserve cloning, pickling, and collection of discarded molecule wrappers.
- Add focused standalone regressions and verify the built wheel separately from
  the source checkout.
- Build wheel and source archive using standard PEP 517 tooling.
- Require RDKit >=2025.9.6, the qualified Atomica version with recursive-match
  controls needed by the patch, and reflect that floor in CI.
- Release only Atomica tags matching the package metadata; do not publish the
  fork to upstream's PyPI project.
- Leave every tracked file and branch in `ad-atomica/atomica` unchanged.

## Status

The cloud GitHub integration rejected fork creation with HTTP 403. Alex directed
use of the GitHub CLI on his Mac; this successfully created the organization fork.
The patch and package are implemented and verified locally. Pushing the fork
and publishing release assets remain pending: the cloud credential cannot push
this repository, and Conductor timed out waiting for the Mac client during
publication. The remote fork still has its upstream default branch, without the
prepared `main` branch. No production deployment or consumer migration is part
of this work; no production release actions are required.
