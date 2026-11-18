#!/usr/bin/env python3
"""Integration tests. Run on the Docker host beside worker.py."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import uuid

import worker


HERE = Path(__file__).resolve().parent


class WorkerTests(unittest.TestCase):
    def request(self, code="print('ok')", **kwargs):
        value = {"job_id": "test-" + uuid.uuid4().hex[:12], "code": code,
                 "sources": [], "timeout_seconds": 15}
        value.update(kwargs)
        proc = subprocess.run([sys.executable, str(HERE / "worker.py")],
                              input=json.dumps(value), text=True, capture_output=True, timeout=55)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stderr, "")
        result = json.loads(proc.stdout)
        self.assertEqual(result["status"], "computational_check_only")
        self.assertLessEqual(len(result["stdout"].encode()) + len(result["stderr"].encode()), 32768)
        return result

    def successful(self, result):
        self.assertEqual(result["exit_code"], 0, result)
        self.assertFalse(result["timed_out"], result)
        self.assertNotIn("worker_error", result)
        self.assertTrue(result["container_cleanup_confirmed"], result)

    def test_exact_sympy(self):
        result = self.request("""from sympy import Matrix, Rational, symbols, expand
import sympy, numpy
x = symbols('x')
assert expand((x+1)**3 - (x**3 + 3*x*x + 3*x + 1)) == 0
assert Matrix([[1, 2], [3, 4]]).det() == -2
assert Rational(1,3) + Rational(1,6) == Rational(1,2)
assert sympy.__version__ == '1.14.0'
assert numpy.__version__ == '2.3.5'
print('exact identities checked')
""")
        self.successful(result)
        self.assertIn("exact identities checked", result["stdout"])
        work = Path(result["artifact_directory"])
        self.assertTrue((work / "request.json").is_file())
        self.assertTrue((work / "check.py").is_file())
        self.assertEqual(json.loads((work / "result.json").read_text()), result)

    def test_false_assertion(self):
        result = self.request("assert 2 + 2 == 5, 'deliberately false'")
        self.assertEqual(result["exit_code"], 1, result)
        self.assertIn("AssertionError", result["stderr"])
        self.assertTrue(result["container_cleanup_confirmed"])

    def test_readonly_sources_and_root_and_privileges(self):
        content = "from sympy import Rational\nVALUE = Rational(2, 7)\n"
        digest = hashlib.sha256(content.encode()).hexdigest()
        result = self.request("""import os
from pathlib import Path
assert os.getuid() == 65534 and os.getgid() == 65534
status = Path('/proc/self/status').read_text()
assert 'CapEff:\\t0000000000000000' in status
assert 'NoNewPrivs:\\t1' in status
source = Path('/sources/experiments/original.py')
assert 'Rational(2, 7)' in source.read_text()
for path in (source, Path('/sources/new.py'), Path('/etc/untrusted-write')):
    try:
        path.write_text('modified')
    except OSError:
        pass
    else:
        raise AssertionError('write permitted: ' + str(path))
Path('/tmp/allowed').write_text('yes')
assert not Path('/var/run/docker.sock').exists()
assert not Path('/dev/nvidia0').exists()
print('source and root write protection enforced')
""", sources=[{"path": "experiments/original.py", "content": content, "sha256": digest}])
        self.successful(result)
        self.assertEqual(result["source_hashes"], {"experiments/original.py": digest})

    def test_no_network(self):
        result = self.request("""import socket
s = socket.socket()
s.settimeout(1)
try:
    s.connect(('1.1.1.1', 443))
except OSError:
    print('outbound network blocked')
else:
    raise AssertionError('unexpected outbound network')
finally:
    s.close()
""")
        self.successful(result)
        self.assertIn("network blocked", result["stdout"])

    def test_timeout_kills_container(self):
        result = self.request("import time\nprint('started', flush=True)\ntime.sleep(60)", timeout_seconds=1.5)
        self.assertTrue(result["timed_out"], result)
        self.assertIn("started", result["stdout"], result)
        self.assertTrue(result["container_cleanup_confirmed"], result)
        proc = subprocess.run(["docker", "inspect", result["container_name"]], capture_output=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertLess(result["wall_seconds"], 25)

    def test_output_budget(self):
        result = self.request("import sys\nsys.stdout.buffer.write(b'x' * 80000)\nsys.stderr.buffer.write(b'y' * 80000)")
        self.successful(result)
        self.assertTrue(result["output_truncated"], result)
        self.assertEqual(len(result["stdout"].encode()) + len(result["stderr"].encode()), 32768)

    def test_bad_hash_and_paths_rejected(self):
        for name, digest in [("safe.py", "0" * 64), ("../escape.py", hashlib.sha256(b"x").hexdigest()),
                             ("/etc/passwd", hashlib.sha256(b"x").hexdigest()),
                             ("a/../../escape", hashlib.sha256(b"x").hexdigest()),
                             ("a\\escape", hashlib.sha256(b"x").hexdigest())]:
            with self.subTest(path=name):
                result = self.request(sources=[{"path": name, "content": "x", "sha256": digest}])
                self.assertIn("worker_error", result)
                self.assertIsNone(result["exit_code"])
                self.assertNotIn("artifact_directory", result)
        result = self.request(sources=[{"path": "symlink", "content": "x", "sha256": hashlib.sha256(b"x").hexdigest(), "type": "symlink"}])
        self.assertIn("symlinks are forbidden", result["worker_error"])

    def test_actual_docker_configuration(self):
        with tempfile.TemporaryDirectory(dir=HERE) as name:
            work = Path(name)
            (work / "sources").mkdir()
            (work / "check.py").write_text("pass")
            container = "tfpt-seam-check-test-" + uuid.uuid4().hex
            try:
                subprocess.run(worker.docker_command(container, work), check=True, capture_output=True)
                info = json.loads(subprocess.check_output(["docker", "inspect", container]))[0]
                config = info["HostConfig"]
                self.assertEqual(config["NetworkMode"], "none")
                self.assertTrue(config["ReadonlyRootfs"])
                self.assertIn("ALL", config["CapDrop"])
                self.assertIn("no-new-privileges:true", config["SecurityOpt"])
                self.assertEqual(config["PidsLimit"], 128)
                self.assertEqual(config["Memory"], 2 * 1024**3)
                self.assertEqual(config["MemorySwap"], 2 * 1024**3)
                self.assertEqual(config["NanoCpus"], 2000000000)
                self.assertEqual(info["Config"]["User"], "65534:65534")
                self.assertFalse(config.get("DeviceRequests"))
                self.assertEqual(set(config["Tmpfs"]), {"/tmp"})
                self.assertIn("size=128m", config["Tmpfs"]["/tmp"])
                self.assertEqual({m["Destination"] for m in info["Mounts"]}, {"/sources", "/job/check.py"})
                self.assertTrue(all(not m["RW"] for m in info["Mounts"]))
            finally:
                self.assertTrue(worker.cleanup(container))


if __name__ == "__main__":
    unittest.main(verbosity=2)
