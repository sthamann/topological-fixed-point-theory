#!/usr/bin/env python3
"""Fetch the official public weights once, with a recorded immutable revision."""
import json
import pathlib
import time
from huggingface_hub import HfApi, snapshot_download

ROOT = pathlib.Path(__file__).resolve().parent
MODEL = "Qwen/Qwen3-32B"
manifest = ROOT / "model_manifest.json"
if manifest.exists():
    revision = json.loads(manifest.read_text())["revision"]
else:
    revision = HfApi().model_info(MODEL).sha
    manifest.write_text(json.dumps({"model": MODEL, "revision": revision,
        "source": "https://huggingface.co/Qwen/Qwen3-32B",
        "created_unix": time.time()}, indent=2) + "\n")
print(f"Downloading {MODEL}@{revision} to {ROOT / 'model'}", flush=True)
path = snapshot_download(MODEL, revision=revision, local_dir=ROOT / "model",
    cache_dir=ROOT / "cache", max_workers=8)
print(f"READY {path}", flush=True)
