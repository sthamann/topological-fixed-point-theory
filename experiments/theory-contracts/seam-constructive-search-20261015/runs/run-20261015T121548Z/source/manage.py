"""Deploy two GPU islands and collect fixed-template Lean receipts on this Mac."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from atlas_inputs import HERE, REPO, save, verify

REMOTE = "/home/aime/tfpt-constructive-search-20261015"
REMOTE_PYTHON = "/home/aime/venvs/torch/bin/python"
HOST = "fieldai"


def verify_inputs(run=None):
    verify(HERE / "atlas_packet.json")
    spec = json.loads((HERE / "domain_spec.json").read_text())
    for source in spec["source_manifest"]:
        path = (REPO / source["path"]).resolve()
        if not path.is_relative_to(REPO.resolve()) or hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError(f"Changed domain source: {source['path']}")
    if run:
        for filename, expected in {**run["deployment"], **run["proof_code"]}.items():
            if hashlib.sha256((HERE / filename).read_bytes()).hexdigest() != expected:
                raise ValueError(f"Changed running search input: {filename}")


def ssh_script(source, timeout=60):
    proc = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", HOST,
                           "python3", "-"], input=source, text=True, capture_output=True, timeout=timeout)
    if proc.returncode:
        raise RuntimeError(proc.stderr[-4000:])
    return json.loads(proc.stdout)


def snapshot(run):
    return ssh_script("""import json, pathlib, os
p=pathlib.Path(%r)
out={}
for i in (0,1):
    for kind in ('status','archive'):
        f=p/f'island-{i}-{kind}.json'
        out[f'{i}-{kind}']=json.loads(f.read_text()) if f.exists() else None
f=p/'launch.json'
out['launch']=json.loads(f.read_text()) if f.exists() else None
f=p/'classification-union.json'
out['classification']=json.loads(f.read_text()) if f.exists() else None
for item in (out.get('launch') or {}).get('workers',[]):
    try: os.kill(item['pid'],0); item['process_exists']=True
    except ProcessLookupError: item['process_exists']=False
print(json.dumps(out))
""" % run["remote_run"])


def deploy():
    verify_inputs()
    files = ["gpu_search.py", "problem.py", "domain_spec.json", "atlas_packet.json"]
    import importlib.util
    spec = importlib.util.spec_from_file_location("deploy_problem", HERE / "problem.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    files += list(getattr(module, "SOURCE_FILES", []))
    files = sorted(set(files))
    for name in files:
        if Path(name).name != name or not (HERE / name).is_file():
            raise ValueError(f"Deployment expects an explicit local basename: {name}")
    manifest = {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest() for f in files}
    revision = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()[:20]
    remote_code = REMOTE + "/code/" + revision
    ssh_script("import pathlib,json; pathlib.Path(%r).mkdir(parents=True,exist_ok=True); print('{}')" % remote_code)
    subprocess.run(["scp", "-q", *[str(HERE / f) for f in files], f"{HOST}:{remote_code}/"], check=True, timeout=120)
    observed = ssh_script("import pathlib,hashlib,json\np=pathlib.Path(%r)\nprint(json.dumps({f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in %r}))" % (remote_code, files))
    if manifest != observed:
        raise ValueError("Deployment hash mismatch")
    return manifest, remote_code


def start(args):
    if args.hours <= 0 or args.generations < 1 or args.population < 2 or args.batch < 1:
        raise ValueError("Invalid positive search bounds")
    manifest, remote_code = deploy()
    run_id = time.strftime("run-%Y%m%dT%H%M%SZ", time.gmtime())
    local_run = HERE / "runs" / run_id
    local_run.mkdir(parents=True, exist_ok=False)
    remote_run = REMOTE + "/runs/" + run_id
    run = {"run_id": run_id, "remote_run": remote_run, "local_run": str(local_run),
           "hours": args.hours, "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "deployment": manifest, "remote_code": remote_code,
           "proof_code": {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest()
                          for f in ("proof_check.py", "lean/FixedCertificate.lean", "manage.py", "atlas_inputs.py")},
           "status": "launching", "physical_solution": False}
    save(local_run / "run.json", run)
    source = """import pathlib,json,subprocess,os,time
p=pathlib.Path(%r); p.mkdir(parents=True,exist_ok=False)
root=pathlib.Path(%r)
workers=[]
for island,gpu in ((0,1),(1,2)):
    env=os.environ.copy(); env['CUDA_VISIBLE_DEVICES']=str(gpu); env['OMP_NUM_THREADS']='4'
    args=[%r,str(root/'gpu_search.py'),'--problem',str(root/'problem.py'),
          '--device','cuda:0','--island',str(island),'--seed',str(731+island),
          '--output',str(p),'--population',str(%r),'--batch',str(%r),
          '--generations',str(%r),'--hours',str(%r),'--source-manifest',str(root/'atlas_packet.json'),
          '--exhaustive-after','--exhaustive-shards','2']
    log=(p/f'island-{island}.log').open('ab')
    proc=subprocess.Popen(args,env=env,cwd=root,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    workers.append({'island':island,'gpu_index':gpu,'pid':proc.pid,'args':args})
data={'workers':workers,'started':time.time()}; (p/'launch.json').write_text(json.dumps(data,indent=2)); print(json.dumps(data))
""" % (remote_run, remote_code, REMOTE_PYTHON, args.population, args.batch, args.generations, args.hours)
    run["launch"] = ssh_script(source)
    run["status"] = "started_unverified"
    save(local_run / "run.json", run)
    save(HERE / "current_run.json", run)
    log = (local_run / "collector.log").open("ab")
    proc = subprocess.Popen([sys.executable, str(HERE / "manage.py"), "collect", "--run", str(local_run / "run.json")],
                            cwd=HERE, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                            start_new_session=True)
    run["collector_pid"] = proc.pid
    if sys.platform == "darwin":
        guard = subprocess.Popen(["caffeinate", "-i", "-w", str(proc.pid)],
                                 stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        run["sleep_guard_pid"] = guard.pid
    save(local_run / "run.json", run)
    save(HERE / "current_run.json", run)
    print(json.dumps(run, indent=2))


def collect(run):
    # The collector checks only fixed generated expressions; it never executes arbitrary candidate code.
    from proof_check import check_candidate
    import problem
    output = Path(run["local_run"])
    deadline = time.monotonic() + (run["hours"] + 0.5) * 3600
    state_file = output / "proof_state.json"
    state = json.loads(state_file.read_text()) if state_file.exists() else {"processed": {}, "status": "running"}
    while time.monotonic() < deadline:
        if (output / "STOP_COLLECTOR").exists():
            state["status"] = "stopped"
            break
        try:
            verify_inputs(run)
            snap = snapshot(run)
            save(output / "snapshot.json", snap)
            queue = {}
            for island in (0, 1):
                archive = snap.get(f"{island}-archive") or {}
                for row in archive.get("candidates", []):
                    if row["result"]["feasible"] and row["key"] not in state["processed"]:
                        queue[row["key"]] = row
            for row in (snap.get("classification") or {}).get("survivors", []):
                if row["result"]["feasible"] and row["key"] not in state["processed"]:
                    queue[row["key"]] = row
            for key, row in list(queue.items())[:4]:
                candidate = problem.proof_candidate(row["genotype"])
                dest = output / "proofs" / key
                dest.mkdir(parents=True, exist_ok=True)
                save(dest / "search_record.json", row)
                receipt = check_candidate(candidate, dest)
                state["processed"][key] = receipt
                save(state_file, state)
            statuses = [snap.get(f"{i}-status") or {} for i in (0, 1)]
            state["last_collected_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            state["numerical_evaluations"] = sum(x.get("total_evaluated_here", 0) for x in statuses)
            state["classified_unique"] = (snap.get("classification") or {}).get("classified_unique", 0)
            state["classification_complete"] = (snap.get("classification") or {}).get("complete", False)
            state["physical_solution"] = False
            state["certified_candidate_count"] = sum(bool(r.get("ok")) for r in state["processed"].values())
            pending = set(queue) - set(state["processed"])
            finished = all(x.get("reason") and x["reason"] != "running" for x in statuses)
            if finished and not pending and state["classification_complete"]:
                state["status"] = "finite_class_finished"
                break
            workers = (snap.get("launch") or {}).get("workers", [])
            if finished and any(x.get("reason", "").endswith(("error", "stopped", "sentinel")) or
                                x.get("reason") == "duration_limit" for x in statuses):
                state["status"] = "incomplete_gpu_search"
                break
            if workers and not any(w.get("process_exists") for w in workers) and not state["classification_complete"]:
                state["status"] = "workers_exited_without_complete_classification"
                break
            state.pop("last_error", None)
        except ValueError as exc:
            # An altered premise cannot silently continue into a certificate.
            state["last_error"] = str(exc)
            state["status"] = "source_changed_stopped"
            ssh_script("import pathlib,json; (pathlib.Path(%r)/'STOP').touch(); print('{}')" % run["remote_run"])
            break
        except Exception as exc:
            state["last_error"] = str(exc)
            state["last_error_at"] = time.time()
        save(state_file, state)
        time.sleep(20)
    else:
        state["status"] = "collection_window_finished"
    save(state_file, state)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["start", "status", "stop", "collect"])
    p.add_argument("--run", type=Path, default=HERE / "current_run.json")
    p.add_argument("--hours", type=float, default=6)
    p.add_argument("--population", type=int, default=128)
    p.add_argument("--batch", type=int, default=128)
    p.add_argument("--generations", type=int, default=20)
    args = p.parse_args()
    if args.command == "start":
        return start(args)
    run = json.loads(args.run.read_text())
    if args.command == "collect":
        return collect(run)
    if args.command == "stop":
        ssh_script("import pathlib,json; (pathlib.Path(%r)/'STOP').touch(); print('{}')" % run["remote_run"])
        (Path(run["local_run"]) / "STOP_COLLECTOR").touch()
    snap = snapshot(run)
    # The status command intentionally omits the large candidate population.
    print(json.dumps({k: v for k, v in snap.items() if not k.endswith("-archive")}, indent=2))


if __name__ == "__main__":
    main()
