"""Discover native R113 source packets and replay the checked consumer."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

from source import evaluate
from verify import digest


SOURCE_DIRS = ("docs/research/fundamental-plan", "experiments/theory-contracts", "research/originlab")
PACKET_SCHEMA = "tfpt.original-blocked-source.v1"
PACKET_MARKERS = (PACKET_SCHEMA.encode(), b'"native_inputs_sha256"')


def packet_in(value):
    if isinstance(value,dict):
        candidate = value if value.get("schema")==PACKET_SCHEMA else value.get("whole_source_bundle")
        if isinstance(candidate,dict) and candidate.get("schema")==PACKET_SCHEMA:
            return candidate
    return None


def packets(root):
    """Read candidate data only. Discovery never executes repository files."""
    for group in SOURCE_DIRS:
        for path in sorted((root/group).rglob("*.json")):
            if not path.resolve().is_relative_to(root.resolve()):
                continue
            try:
                data = path.read_bytes()
                if not all(token in data for token in PACKET_MARKERS):
                    continue
                value = packet_in(json.loads(data))
            except (OSError, ValueError):
                continue
            if value is not None:
                yield path,value


def inventory(root, domain):
    if not domain.is_file():
        return {"matching":[],"rejected":[]}
    expected = digest(domain)
    name = json.loads(domain.read_bytes())["domain_id"]
    matching,rejected = [],[]
    for path,value in packets(root):
        entry = {"path":str(path.relative_to(root)),"sha256":digest(path)}
        if value.get("input_role")!="source_reviewed_original":
            rejected.append({**entry,"reason":"CONTROL_OR_UNREVIEWED_ROLE"})
        elif value.get("domain_id")!=name or value.get("native_inputs_sha256")!=expected:
            rejected.append({**entry,"reason":"DIFFERENT_ORIGINAL_DOMAIN"})
        else:
            dependencies = {}
            references = value.get("source_refs")
            for reference in references if isinstance(references,list) else ():
                if isinstance(reference,dict) and isinstance(reference.get("path"),str):
                    ref = (root/reference["path"]).resolve()
                    if ref.is_relative_to(root.resolve()):
                        dependencies[reference["path"]] = digest(ref) if ref.is_file() else None
            matching.append({**entry,"dependencies":dependencies})
    return {"matching":matching,"rejected":rejected}


def replay(root, domain, entry):
    path = (root/entry["path"]).resolve()
    if not path.is_relative_to(root.resolve()) or digest(path)!=entry["sha256"]:
        raise ValueError("discovered packet changed before replay")
    original = json.loads(domain.read_bytes())
    if original.get("whole_source_bundle") is not None:
        raise ValueError("automatic replay cannot replace an attached original source")
    packet = packet_in(json.loads(path.read_bytes()))
    if packet is None or packet.get("input_role")!="source_reviewed_original":
        raise ValueError("automatic replay requires a source-reviewed input role")
    source_sha = digest(domain)
    if packet.get("domain_id")!=original.get("domain_id") or packet.get("native_inputs_sha256")!=source_sha:
        raise ValueError("discovered packet is not bound to the selected native domain")
    references = packet.get("source_refs")
    if not isinstance(references,list) or any(not isinstance(ref,dict) or not isinstance(ref.get("path"),str)
           or not (root/ref["path"]).resolve().is_relative_to(root.resolve())
           for ref in references):
        raise ValueError("discovered packet source references must stay in the repository")
    overlay = {**original,"whole_source_bundle":packet,
               "original_domain_binding":{"path":str(domain.resolve()),"sha256":source_sha}}
    with TemporaryDirectory(prefix="originlab-native-source-") as work:
        candidate = Path(work)/"native.json"
        candidate.write_text(json.dumps(overlay),encoding="utf-8")
        result = evaluate(root,candidate)
    return {**result,"discovered_packet":entry,"original_domain_sha256":source_sha,
            "physical_solution":False}
