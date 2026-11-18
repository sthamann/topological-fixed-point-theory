"""Call both unchanged entry points of the original Q source."""
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

repo=Path(sys.argv[1]).resolve()
root=repo/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
domain=root/"domain/native-diagnostic-v1.json"
os.environ["TFPT_FIELD_DOMAIN"]=str(domain)
sys.path.insert(0,str(root/"runtime"))
import field_problem as source

candidate=source.canonicalize({"schema":source.SCHEMA,"domain_sha256":source.DOMAIN["_sha256"],"genes":[2,2,2]})
assembled=source.A.assemble(candidate)
direct=source.ADAPTER.evaluate_source(candidate,assembled,source.DOMAIN)
with TemporaryDirectory() as temp:
    process=subprocess.run([sys.executable,"-B",str(root/"runtime/launch.py"),"--domain",str(domain),
                            "--out",temp,"--mode","integrate"],capture_output=True,text=True)
    pre=json.loads((Path(temp)/"preflight.json").read_text())
    out=json.loads((Path(temp)/"result.json").read_text())
print(json.dumps({
    "domain_sha256":source.DOMAIN["_sha256"], "candidate_id":assembled["candidate_id"],
    "direct":direct,"stage":{
       "status":out["status"],"exit_code":process.returncode,
       "native_atoms":pre["original_atom_count"],
       "diagnostic_configurations":pre["finite_diagnostic_candidate_count"],
       "integration_plan":pre["integration_plan"],
    }},ensure_ascii=False,indent=2))
assert direct["status"]=="original_full_source_operation_not_supplied"
assert all(direct.get(k) is None for k in ("Z0","ZU","ZY","ZUY"))
assert process.returncode==2 and pre["original_atom_count"]==0
