#!/usr/bin/env bash
# Run from the repository root with the locked uv environment already installed.
set -euo pipefail
export MISE_UV_VERSION="${MISE_UV_VERSION:-0.11.24}"
audit_dir=docs/audits/2026-11-03-titratable-preservation
audit_output=.context/titratable-audit
mkdir -p "$audit_output"
for audit_script in audit alias_probes edge_probes contacts baseline; do
  uv run --no-sync python "$audit_dir/$audit_script.py" > "$audit_output/$audit_script.log" 2>&1
done
uv run --no-sync python "$audit_dir/validate.py"
