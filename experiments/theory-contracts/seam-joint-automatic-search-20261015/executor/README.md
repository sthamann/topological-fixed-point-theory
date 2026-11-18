# Isolated computational checker

This worker executes proposed Python calculations inside a CPU-only Docker container. Every result has `status: computational_check_only`: passing calculations are not a proof of a derivation, physical source selection, or the full TFPT claim. The controller must review premises, assertions and original sources separately. Generated output is untrusted data.

The host worker uses only Python's standard library. It never executes generated code on the host, imports submitted modules on the host, installs packages in a host environment, or mounts the complete repository. The image contains Python 3.11, SymPy 1.14.0 and NumPy 2.3.5. The base `python:3.11-slim` tag is mutable; each run resolves the checker tag to its immutable local image ID, uses that ID for container creation, and records it under `container_isolation.image_id`.

Build and integration tests on the remote Docker host:

```sh
docker build -t tfpt-seam-check:20261015 /home/aime/tfpt-seam-search-20261015/executor
python3 /home/aime/tfpt-seam-search-20261015/executor/test_worker.py
```

The controller invokes `ssh fieldai python3 /home/aime/tfpt-seam-search-20261015/executor/worker.py` and writes one JSON document to stdin, then closes stdin. The worker writes exactly one JSON document to stdout. Example input:

```json
{"job_id":"candidate-0001","code":"from sympy import Rational\nassert Rational(1, 3) + Rational(1, 6) == Rational(1, 2)\nprint('identity holds')","sources":[],"timeout_seconds":90}
```

Optional explicit original sources are objects containing exactly `path`, UTF-8 `content` and `sha256`. The worker verifies hashes before creating a job, reconstructs regular files in a fresh private directory, and mounts them read-only at `/sources/<path>`. Absolute paths, traversal, duplicate/conflicting paths, noncanonical path components and symlink entries are rejected. Paths name the submitted source tree, never arbitrary host files. Source hashes establish content identity, not source authenticity; selecting genuine original sources is the controller's responsibility. Inputs are limited to 20 MiB total, code to 1 MiB, source contents to 16 MiB, and source count to 256.

Each container has network `none`, a read-only root, all Linux capabilities dropped, `no-new-privileges`, UID/GID 65534, 128-process limit, 2 GiB memory including swap, two CPUs, and a 128 MiB `/tmp` tmpfs with `noexec,nosuid,nodev`. Docker's default seccomp policy remains enabled. Only the submitted checker and explicit source tree are mounted. No GPU, secret environment values, Docker socket or other host directories are supplied. The Docker daemon and its standard container/kernel boundary remain trusted.

The deadline is at most 90 seconds (a smaller positive `timeout_seconds` is permitted). Container creation reserves a unique unpredictable name before attachment. Cleanup always targets only that name with `docker rm -f`, including after a timeout or catchable interrupt; killing only the Docker client would leave a workload running. `container_cleanup_confirmed` reports whether removal or absence was confirmed. An uncatchable host kill, host crash or failed Docker daemon can still prevent cleanup and must be handled operationally; this worker does not claim a VM-level security boundary.

Return fields include `exit_code`, `stdout`, `stderr`, `timed_out`, `container_isolation`, `source_hashes`, `wall_seconds`, `status`, `container_cleanup_confirmed`, `container_name`, `artifact_directory` and `output_truncated`. `exit_code` is null on timeout or worker failure; `worker_error` distinguishes malformed requests and infrastructure failures from a nonzero checker exit. Combined captured stdout/stderr is limited to 32 KiB and excess output is drained and discarded. Container logs are disabled to avoid unbounded daemon log storage.

Requests, source snapshots, checker code and result JSON persist in `executor/jobs/<safe-job-id>-<random>/` on the remote host; a retry never overwrites an earlier job. These are audit artifacts, not executable host scripts. The worker removes containers but retains those job directories for controller inspection. No broad Docker cleanup or service changes are performed.

The integration suite checks an exact symbolic identity and rational arithmetic, a false assertion, source/root write protection, nonroot/capability restrictions, outbound network denial, timeout removal, output truncation, hash/path/symlink validation, and the actual Docker configuration including memory/CPU/process limits and the mount allowlist.
