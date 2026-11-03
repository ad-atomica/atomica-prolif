# ruff: noqa: T201
# Audit command-line tools intentionally print their results.
"""Compare current legacy behavior to the same fixtures on origin/main's class."""

import importlib.util
import json
import subprocess
from collections.abc import Callable
from pathlib import Path

import audit
from rdkit import Chem

from prolif import Molecule


def run() -> None:
    out = Path(".context/titratable-audit")
    path = out / "baseline_standardizer.py"
    path.write_text(
        subprocess.check_output(
            ["git", "show", "origin/main:prolif/io/molecule_standardizer.py"], text=True
        )
    )
    spec = importlib.util.spec_from_file_location("prolif.io.audit_baseline", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load baseline standardizer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def baseline_factory(
        *, preserve_histidine: bool = False
    ) -> Callable[[Chem.Mol], Molecule]:
        assert not preserve_histidine

        def standardize(source: Chem.Mol) -> Molecule:
            result = module.MoleculeStandardizer()(source)
            if not isinstance(result, Molecule):
                raise TypeError("Baseline standardizer did not return a Molecule")
            return result

        return standardize

    audit.standardizer_factory = baseline_factory
    audit.rdBase.DisableLog("rdApp.warning")
    audit.rdBase.DisableLog("rdApp.error")
    result = [
        audit.run_case(*s, representation, False)
        for s in audit.case_specs()
        for representation in audit.REPRESENTATIONS
    ]
    path = out / "baseline-results.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    # Normalize tuples to JSON arrays before equality comparison.
    baseline = json.loads(path.read_text())
    current = [
        r
        for r in json.loads((out / "results.json").read_text())["cases"]
        if not r["preserve_histidine"]
    ]
    identical = sum(a == b for a, b in zip(baseline, current, strict=False))
    print("Baseline cases:", len(baseline), "identical to current legacy:", identical)
    assert len(baseline) == len(current) == identical


if __name__ == "__main__":
    run()
