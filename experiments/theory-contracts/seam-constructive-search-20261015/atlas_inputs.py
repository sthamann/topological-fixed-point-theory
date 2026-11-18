"""Read source-pinned search inputs through the canonical Atlas MCP, without AI calls."""
from __future__ import annotations
import argparse
import hashlib
import json
import select
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PROBLEMS = [
    "problem:common-source-state-event-selection",
    "problem:native-charged-seam-selection",
    "problem:future-quotient-and-gns-reconstruction",
]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def save(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    tmp.replace(path)


class Atlas:
    def __init__(self):
        self.process = subprocess.Popen(
            [str(REPO / "semantic_explorer/run.sh"), "mcp"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, bufsize=1,
        )
        self.serial = 0

    def call(self, name, arguments):
        self.serial += 1
        request = {"jsonrpc": "2.0", "id": self.serial, "method": "tools/call",
                   "params": {"name": name, "arguments": arguments}}
        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()
        if not select.select([self.process.stdout], [], [], 120)[0]:
            raise TimeoutError(f"Atlas {name} exceeded 120s")
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError(f"Atlas exited: {self.process.poll()}")
        response = json.loads(line)
        if response.get("id") != self.serial or "error" in response:
            raise RuntimeError(f"Atlas protocol error: {response}")
        result = response["result"]
        if result.get("isError"):
            raise RuntimeError(f"Atlas rejected {name}: {result}")
        return json.loads(result["content"][0]["text"])

    def close(self):
        self.process.terminate()
        try:
            self.process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self.process.kill()
            self.process.wait()


def prepare(output: Path):
    atlas = Atlas()
    try:
        status = atlas.call("graph_status", {})
        overview = atlas.call("theory_overview", {})
        states = [atlas.call("research_state", {"query": p, "lang": "en", "max_chars": 60000})
                  for p in PROBLEMS]
        refs = {}
        for state in states:
            if state.get("budget", {}).get("truncated"):
                raise ValueError("Incomplete research packet: increase bounded packet budget")
            for problem in state.get("problems", []):
                for finding in problem.get("findings", []):
                    if finding.get("source_state") != "current":
                        raise ValueError(f"Stale finding: {finding['id']}")
                    for ref in finding.get("sources", []):
                        refs[(ref["source_id"], ref["line_start"], ref["line_end"])] = ref
        excerpts = []
        for (source_id, start, end), ref in refs.items():
            result = atlas.call("read_source", {"source_id": source_id,
                                               "line_start": start, "line_end": end})
            path = (REPO / result["path"]).resolve()
            if not path.is_relative_to(REPO.resolve()) or not path.is_file():
                raise ValueError("Invalid repository source")
            actual = digest(path.read_bytes())
            if actual != result["sha256"] or actual != ref["sha256"]:
                raise ValueError(f"Changed original source: {result['path']}")
            excerpts.append(result)
        packet = {"schema": 1, "collected_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                  "purpose": "constructive source/mark/state/time search; no automatic promotion",
                  "atlas_revision": states[0]["revision"], "graph_status": status,
                  "theory_overview": overview, "research_states": states, "sources": excerpts,
                  "proof_boundary": "Atlas statements retain their premises. Computational candidates do not close P1/P2 or T1-T8."}
        save(output, packet)
        print(json.dumps({"path": str(output), "sha256": digest(output.read_bytes()),
                          "atlas_revision": packet["atlas_revision"], "source_excerpts": len(excerpts)}))
        return packet
    finally:
        atlas.close()


def verify(packet_path: Path):
    packet = json.loads(packet_path.read_text())
    changed = []
    for source in packet["sources"]:
        path = (REPO / source["path"]).resolve()
        if not path.is_relative_to(REPO.resolve()) or not path.is_file() or digest(path.read_bytes()) != source["sha256"]:
            changed.append(source["path"])
    if changed:
        raise ValueError(f"Source packet changed: {changed}")
    return {"source_state": "current", "sha256": digest(packet_path.read_bytes())}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "atlas_packet.json")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(verify(args.output)))
    else:
        prepare(args.output)
