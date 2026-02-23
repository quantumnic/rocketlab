#!/usr/bin/env bash
# Build rocketlab WASM module and copy artefacts into web/pkg/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "🦀 Building rocketlab WASM module..."
cd "$PROJECT_DIR"

# Check for wasm-pack
if ! command -v wasm-pack &>/dev/null; then
  echo "⚠️  wasm-pack not found. Installing..."
  cargo install wasm-pack
fi

wasm-pack build --target web --features wasm --out-dir web/pkg --no-typescript

echo "✅ WASM build complete → web/pkg/"
ls -lh web/pkg/*.wasm web/pkg/*.js 2>/dev/null || true
