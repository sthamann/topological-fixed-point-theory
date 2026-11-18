#!/usr/bin/env python3
"""Execute untrusted proposed checks only inside the restricted CPU container."""
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import selectors
import signal
import subprocess
import sys
import tempfile
import time
import uuid

IMAGE = "tfpt-seam-check:20261015"
OUTPUT_LIMIT = 32768
REQUEST_LIMIT = 20 * 1024 * 1024
ISOLATION = {
    "network": "none", "read_only_root": True, "capabilities": "none",
    "no_new_privileges": True, "pids_limit": 128, "memory_bytes": 2147483648,
    "cpus": 2, "tmpfs_tmp_bytes": 134217728, "user": "65534:65534",
    "gpus": False, "host_secrets": False, "docker_socket": False,
    "sources": "explicit_sha256_verified_read_only_files", "image": IMAGE,
}


def validate(request):
    if not isinstance(request, dict):
        raise ValueError("request must be an object")
    job_id = request.get("job_id")
    if not isinstance(job_id, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", job_id):
        raise ValueError("invalid job_id")
    code = request.get("code")
    if not isinstance(code, str) or len(code.encode("utf-8")) > 1024 * 1024:
        raise ValueError("code must be a UTF-8 string of at most 1 MiB")
    timeout = request.get("timeout_seconds", 90)
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or not 0 < timeout <= 90:
        raise ValueError("timeout_seconds must be > 0 and <= 90")
    sources = request.get("sources", [])
    if not isinstance(sources, list) or len(sources) > 256:
        raise ValueError("sources must be a list of at most 256 files")
    seen, checked, total = set(), [], 0
    for item in sources:
        if not isinstance(item, dict) or set(item) != {"path", "content", "sha256"}:
            raise ValueError("each source must contain only path, content, sha256; symlinks are forbidden")
        name, content, digest = item["path"], item["content"], item["sha256"]
        if not isinstance(name, str) or not name or len(name) > 1024 or "\\" in name:
            raise ValueError("invalid source path")
        parts = name.split("/")
        if PurePosixPath(name).is_absolute() or any(p in ("", ".", "..") or not re.fullmatch(r"[A-Za-z0-9_.-]+", p) for p in parts):
            raise ValueError("absolute, traversal, and noncanonical source paths are forbidden")
        if name in seen or any(name.startswith(old + "/") or old.startswith(name + "/") for old in seen):
            raise ValueError("duplicate or conflicting source path")
        if not isinstance(content, str) or not isinstance(digest, str) or not re.fullmatch(r"[a-fA-F0-9]{64}", digest):
            raise ValueError("invalid source content or sha256")
        data = content.encode("utf-8")
        actual = hashlib.sha256(data).hexdigest()
        if actual != digest.lower():
            raise ValueError("source SHA-256 mismatch: " + name)
        total += len(data)
        if total > 16 * 1024 * 1024:
            raise ValueError("sources exceed 16 MiB")
        seen.add(name)
        checked.append((name, data, actual))
    return job_id, code, checked, float(timeout)


def docker_command(name, work, image=IMAGE):
    return ["docker", "create", "--rm", "--name", name, "--network", "none",
            "--read-only", "--cap-drop", "ALL", "--security-opt", "no-new-privileges:true",
            "--pids-limit", "128", "--memory", "2g", "--memory-swap", "2g",
            "--cpus", "2", "--user", "65534:65534", "--log-driver", "none",
            "--tmpfs", "/tmp:rw,noexec,nosuid,nodev,size=128m", "--workdir", "/tmp",
            "--env", "PYTHONDONTWRITEBYTECODE=1", "--env", "PYTHONUNBUFFERED=1",
            "--mount", f"type=bind,source={work / 'check.py'},target=/job/check.py,readonly",
            "--mount", f"type=bind,source={work / 'sources'},target=/sources,readonly",
            image, "python3", "-I", "-u", "/job/check.py"]


def cleanup(name):
    # Never use broad container/prune operations: delete only this unpredictable name.
    try:
        done = subprocess.run(["docker", "rm", "-f", name], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, timeout=15, check=False)
        if done.returncode == 0 or b"No such container" in done.stderr:
            return True
    except (OSError, subprocess.TimeoutExpired):
        pass
    return False


def run(request):
    started = time.monotonic()
    job_id, code, sources, timeout = validate(request)
    root = Path(__file__).resolve().parent / "jobs"
    if root.is_symlink():
        raise ValueError("jobs directory may not be a symlink")
    root.mkdir(mode=0o700, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=job_id + "-", dir=root))
    source_root = work / "sources"
    source_root.mkdir()
    (work / "request.json").write_text(json.dumps(request, ensure_ascii=False), encoding="utf-8")
    (work / "check.py").write_text(code, encoding="utf-8")
    (work / "check.py").chmod(0o444)
    for name, data, _ in sources:
        path = source_root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        # The tree is fresh and never exposed writable to the untrusted process.
        if path.is_symlink() or any(p.is_symlink() for p in path.parents if p != work):
            raise ValueError("source symlink forbidden")
        with path.open("xb") as stream:
            stream.write(data)
        path.chmod(0o444)
    for path in sorted(source_root.rglob("*"), reverse=True):
        if path.is_dir():
            path.chmod(0o555)
    source_root.chmod(0o555)
    name = "tfpt-seam-check-" + uuid.uuid4().hex
    result = {"job_id": job_id, "exit_code": None, "stdout": "", "stderr": "",
              "timed_out": False, "container_isolation": ISOLATION.copy(),
              "source_hashes": {n: h for n, _, h in sources},
              "status": "computational_check_only", "output_truncated": False,
              "artifact_directory": str(work), "container_name": name}
    process = None
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    remaining = OUTPUT_LIMIT
    selector = selectors.DefaultSelector()
    old_handlers = {}

    def interrupted(signum, frame):
        raise InterruptedError("worker interrupted by signal " + str(signum))

    try:
        for sig in (signal.SIGTERM, signal.SIGINT):
            old_handlers[sig] = signal.signal(sig, interrupted)
        deadline = time.monotonic() + timeout
        inspected = subprocess.run(["docker", "image", "inspect", "--format", "{{.Id}}", IMAGE],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   timeout=min(timeout, 15), check=False)
        image_id = inspected.stdout.decode("ascii", errors="replace").strip()
        if inspected.returncode or not re.fullmatch(r"sha256:[a-f0-9]{64}", image_id):
            raise OSError("required checker image is unavailable")
        result["container_isolation"]["image_id"] = image_id
        created = subprocess.run(docker_command(name, work, image_id), stdin=subprocess.DEVNULL,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 timeout=max(0.001, min(deadline - time.monotonic(), 30)), check=False)
        if created.returncode:
            raise OSError("docker create failed: " + created.stderr[:4096].decode("utf-8", errors="replace"))
        process = subprocess.Popen(["docker", "start", "--attach", name], stdin=subprocess.DEVNULL,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for label in buffers:
            pipe = getattr(process, label)
            os.set_blocking(pipe.fileno(), False)
            selector.register(pipe, selectors.EVENT_READ, label)
        while selector.get_map():
            left = deadline - time.monotonic()
            if left <= 0:
                result["timed_out"] = True
                break
            for key, _ in selector.select(min(0.1, left)):
                chunk = os.read(key.fileobj.fileno(), 65536)
                if not chunk:
                    selector.unregister(key.fileobj)
                    continue
                take = min(remaining, len(chunk))
                buffers[key.data].extend(chunk[:take])
                remaining -= take
                if take < len(chunk):
                    result["output_truncated"] = True
        if not result["timed_out"]:
            try:
                result["exit_code"] = process.wait(timeout=max(0.001, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                result["timed_out"] = True
    except subprocess.TimeoutExpired:
        result["timed_out"] = True
    except (OSError, InterruptedError) as exc:
        result["worker_error"] = str(exc)
    finally:
        # docker-client termination alone does not stop its container.
        try:
            if process is not None:
                if process.poll() is None:
                    process.kill()
                process.wait(timeout=5)
        except (OSError, subprocess.TimeoutExpired) as exc:
            result["worker_error"] = "Docker client cleanup: " + str(exc)
        finally:
            result["container_cleanup_confirmed"] = cleanup(name)
            if process is not None:
                for pipe in (process.stdout, process.stderr):
                    pipe.close()
        selector.close()
        for sig, previous in old_handlers.items():
            signal.signal(sig, previous)
    remaining = OUTPUT_LIMIT
    for label, content in buffers.items():
        decoded = bytes(content).decode("utf-8", errors="replace").encode("utf-8")
        result[label] = decoded[:remaining].decode("utf-8", errors="ignore")
        remaining -= len(result[label].encode("utf-8"))
    result["wall_seconds"] = round(time.monotonic() - started, 6)
    (work / "result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def main():
    started = time.monotonic()
    try:
        raw = sys.stdin.buffer.read(REQUEST_LIMIT + 1)
        if len(raw) > REQUEST_LIMIT:
            raise ValueError("request exceeds 20 MiB")
        result = run(json.loads(raw))
    except Exception as exc:
        result = {"exit_code": None, "stdout": "", "stderr": "", "timed_out": False,
                  "container_isolation": ISOLATION.copy(), "source_hashes": {},
                  "status": "computational_check_only", "worker_error": str(exc),
                  "wall_seconds": round(time.monotonic() - started, 6)}
    print(json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
