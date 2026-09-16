#!/usr/bin/env bash

set -euo pipefail

if [[ "${CONDA_DEFAULT_ENV:-}" != "bear-load-practical" ]]; then
  echo "ERROR: activate bear-load-practical before running this script." >&2
  exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BIN_DIR="$SCRIPT_DIR/bin"
mkdir -p "$BIN_DIR"

for TOOL in minimap2 transanno liftOver bigWigToBedGraph; do
  TOOL_PATH=$(command -v "$TOOL" || true)
  if [[ -z "$TOOL_PATH" ]]; then
    echo "ERROR: $TOOL was not found in the active Conda environment." >&2
    exit 1
  fi

  ln -sfn "$TOOL_PATH" "$BIN_DIR/$TOOL"
  printf '%-20s -> %s\n' "$BIN_DIR/$TOOL" "$TOOL_PATH"
done

echo "All Day 2 GERP command links are ready."
