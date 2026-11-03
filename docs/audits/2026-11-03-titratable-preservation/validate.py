# ruff: noqa: T201
# Audit command-line tools intentionally print their results.
"""Verify the audit's observations and emit a compact, reviewable CSV.

These are historical reproduction assertions, not desired product behavior.
After a fix, failed assertions can mean a reported defect no longer reproduces.
"""

import csv
import hashlib
import json
import platform
import subprocess
from importlib.metadata import version
from pathlib import Path

import numpy as np
from evidence import JSONRecord, JSONValue, array, integer, load_json, record, records
from rdkit import rdBase

import prolif

OUT = Path(".context/titratable-audit")


def load(name: str) -> JSONValue:
    return load_json(OUT / name)


def signature(snap: JSONValue) -> str:
    return "; ".join(
        f"{n}:q={record(a)['q']},H={record(a)['h']},rad={record(a)['rad']}"
        for n, a in record(snap).items()
    )


def site(row: JSONRecord, stage: str, name: str) -> JSONRecord:
    return record(record(row[stage])[name])


def total_charge(row: JSONRecord, stage: str) -> int:
    return sum(integer(record(a)["q"]) for a in record(row[stage]).values())


def run() -> None:
    rows = records(record(load("results.json"))["cases"])
    baseline = load("baseline-results.json")
    aliases = record(load("alias-probe-results.json"))
    edges = records(load("edge-results.json"))
    contacts = records(load("contact-results.json"))
    assert len(rows) == 560
    assert baseline == [r for r in rows if not r["preserve_histidine"]]
    assert all(r["source"] == r["remove_h_source"] for r in rows)
    succeeded = [r for r in rows if "error" not in r]
    assert all(r["heavy_order_coords_preserved"] for r in succeeded)
    assert all(r["source_chemistry_unchanged"] for r in succeeded)
    for mode in (False, True):
        selected = [r for r in rows if r["preserve_histidine"] == mode]
        assert sum("error" in r for r in selected) == 42
        assert sum(r.get("residue_preserved") is True for r in selected) == 140
        assert sum(r.get("residue_preserved") is False for r in selected) == 98

    def case(
        family: str,
        state: str,
        alias: str,
        representation: str = "no_graph_h",
        preserve: bool = True,
    ) -> JSONRecord:
        return next(
            r
            for r in rows
            if (
                r["family"],
                r["state"],
                r["alias"],
                r["representation"],
                r["preserve_histidine"],
            )
            == (family, state, alias, representation, preserve)
        )

    lys = case("LYS", "neutral", "LYN")
    assert (site(lys, "source", "NZ")["q"], site(lys, "source", "NZ")["h"]) == (0, 2)
    assert (site(lys, "parent", "NZ")["q"], site(lys, "parent", "NZ")["h"]) == (1, 3)
    arg = case("ARG", "cation1", "ARG")
    assert total_charge(arg, "source") == 1
    assert total_charge(arg, "parent") == 2
    assert site(case("SS", "disulfide", "CYX"), "parent", "SG")["q"] == -1
    assert (
        site(case("LYS", "neutral", "LYN", "fixed_h_counts"), "parent", "NZ")["rad"]
        == 1
    )
    assert case("TYR", "anion", "TYR")["parent_preserved"]
    assert case("NTERM", "cation", "ALA")["parent_preserved"]
    assert case("CTERM", "anion", "ALA")["parent_preserved"]
    for row in records(aliases["atom_aliases"]):
        neighbors = record(record(row["result"])["residue"])["neighbors"]
        assert neighbors == (["C"] if row["alias"] == row["canonical"] else [])
    assert all(
        r["canonical_names"] == r["arbitrary_names"]
        for r in edges
        if r["probe"] == "canonical_h_names"
    )
    assert len(contacts) == 104
    assert sum("error" in r for r in contacts) == 18
    contact = next(
        r
        for r in contacts
        if r["family"] == "LYS"
        and r["representation"] == "no_graph_h"
        and r["preserve_histidine"]
    )
    hits = record(record(record(contact["detectors"])["protein_acceptor"])["hits"])
    assert (
        len(array(hits["source"])) == 92
        and hits["standardized"] == hits["removed"] == []
    )
    hits = record(record(record(contact["detectors"])["protein_cation"])["hits"])
    assert (
        hits["source"] == []
        and len(array(hits["standardized"])) == len(array(hits["removed"])) == 96
    )

    fields = [
        "family",
        "state",
        "alias",
        "representation",
        "preserve_histidine",
        "source",
        "parent",
        "residue",
        "removed_parent",
        "removed_residue",
        "parent_preserved",
        "residue_preserved",
        "identity_preserved",
        "heavy_order_coords_preserved",
        "error",
    ]
    csv_path = Path(__file__).with_name("matrix.csv")
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: signature(row[k])
                    if isinstance(row.get(k), dict)
                    else row.get(k, "")
                    for k in fields
                }
            )

    meta = {
        "python": platform.python_version(),
        "rdkit": rdBase.rdkitVersion,
        "numpy": np.__version__,
        "gemmi": version("gemmi"),
        "prolif": prolif.__version__,
        "head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "base": subprocess.check_output(
            ["git", "rev-parse", "origin/main"], text=True
        ).strip(),
        "matrix_cases": len(rows),
        "baseline_cases": len(array(baseline)),
        "contact_scenarios": len(contacts),
        "contact_standardization_errors": 18,
        "sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in OUT.glob("*results.json")
        },
    }
    OUT.joinpath("manifest.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(
        "Verified 560 matrix cases, 280 baseline comparisons, "
        "104 contact scenarios, alias and edge controls."
    )
    print("All 560 selected-site source snapshots survive RemoveHs unchanged.")
    print(
        "Per mode: 140 unchanged site q/H/rad signatures, "
        "98 changed, 42 standardization errors."
    )
    print("These counts describe designed probes, not rates in real receptors.")
    print(f"Wrote {csv_path} and {OUT / 'manifest.json'}.")


if __name__ == "__main__":
    run()
