#!/bin/sh
set -eu
ATLAS_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ATLAS_BIN="$ATLAS_DIR/server/target/release/semantic-explorer-server"
if [ ! -x "$ATLAS_BIN" ] || [ "$ATLAS_DIR/server/Cargo.toml" -nt "$ATLAS_BIN" ] || [ "$ATLAS_DIR/server/Cargo.lock" -nt "$ATLAS_BIN" ] || [ "$ATLAS_DIR/../website/lib/atlas-status.json" -nt "$ATLAS_BIN" ] || [ -n "$(find "$ATLAS_DIR/server/src" -type f -name '*.rs' -newer "$ATLAS_BIN" -print -quit 2>/dev/null)" ]; then
  cargo build --release --quiet --manifest-path "$ATLAS_DIR/server/Cargo.toml"
fi
if [ -z "${OPENAI_API_KEY_FILE:-}" ] && [ -f "$ATLAS_DIR/local-config.json" ]; then
  OPENAI_API_KEY_FILE=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("key_file", ""))' "$ATLAS_DIR/local-config.json")
  export OPENAI_API_KEY_FILE
fi
if [ "$#" -eq 0 ]; then set -- serve; fi
exec "$ATLAS_BIN" --data "$ATLAS_DIR/data/knowledge.json" --repo "$ATLAS_DIR/.." --web "$ATLAS_DIR/web" "$@"
