"""Optional independent theorem and inverse synthesis backends."""
from pathlib import Path
from tempfile import TemporaryDirectory
import re
import subprocess

from verify import digest


def inverse_rank(low, high):
    try:
        import z3
    except ImportError:
        return dict(status="UNAVAILABLE", backend="z3-solver")
    n = z3.Int("rank_parameter")
    quotient = z3.Int("positive_divisibility_witness")
    solver = z3.Solver()
    solver.set(timeout=15000)
    # The independently derived divisibility constraint is the interface;
    # the rank answer is absent from the solver input.
    solver.add(n >= low, n <= high, quotient > 0, n*quotient == 4)
    out = []
    while (result := solver.check()) == z3.sat:
        model = solver.model()
        found = model[n].as_long()
        witness = model[quotient].as_long()
        if found*witness != 4 or not low <= found <= high or witness <= 0:
            raise ArithmeticError("solver witness failed exact replay")
        out.append(found)
        solver.add(n != found)
    return dict(status="EXHAUSTED" if result == z3.unsat else "UNKNOWN",
                backend=z3.get_version_string(), solutions=sorted(out), bounds=[low,high],
                constraint="n*q=4, q>0 after primitive integral projection")


def lean_check(path, toolchain):
    if toolchain is None:
        return dict(status="NOT_RUN")
    command = ["elan","run",toolchain,"lean"]
    try:
        with TemporaryDirectory() as tmp:
            result = subprocess.run([*command,str(path.resolve())],cwd=tmp,text=True,
                                    capture_output=True,timeout=45,check=False)
            version = subprocess.check_output([*command,"--version"],text=True,timeout=10).strip()
    except (OSError,subprocess.SubprocessError) as error:
        return dict(status="UNAVAILABLE",error=str(error))
    output = result.stdout+result.stderr
    audit = re.search(r"originlab_rank_from_divisibility' depends on axioms: \[([^\]]*)\]",output)
    axioms = set(audit.group(1).split(", ")) if audit and audit.group(1) else set()
    valid = result.returncode == 0 and audit is not None and axioms <= {
        "propext","Classical.choice","Quot.sound",
    }
    return dict(status="KERNEL_CHECKED" if valid else "FAILED",sha256=digest(path),
                toolchain=toolchain,version=version,output=output,
                scope="rank from n dividing 4; the D/A discriminant argument is a premise")
