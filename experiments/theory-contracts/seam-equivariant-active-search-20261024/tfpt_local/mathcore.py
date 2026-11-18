"""Actual E8 sign-lift relations, using the existing source basis and cocycle.

T = D_a S_T0, j = D_b S_J. D_a acts on the OUTPUT lattice coordinate.
We test T^2 = 1 and T j T j = 1 on the integral lattice generators.
This is a finite reflection-lift subproblem, not the four-puncture classification.
"""
from __future__ import annotations

import ast
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as S

DEFAULT_REPO = Path("/Users/stefanhamann/Projekte/tfpt-theoryv4")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def parity(x):
    return x.bit_count() & 1


def bits(x, n=8):
    return [(x >> i) & 1 for i in range(n)]


def mask(values):
    if not isinstance(values, list) or len(values) != 8:
        raise ValueError("Expected exactly eight binary integers")
    if any(type(x) is not int or x not in (0, 1) for x in values):
        raise ValueError("Only the integers 0 and 1 are allowed")
    return sum(x << i for i, x in enumerate(values))


def exact_gf2(rows, rhs, n=8):
    """Independent Gaussian elimination; returns affine solution or contradiction."""
    work = [[int(r), int(b), 1 << i] for i, (r, b) in enumerate(zip(rows, rhs))]
    pivot = 0
    pivots = []
    for col in range(n):
        found = next((i for i in range(pivot, len(work)) if (work[i][0] >> col) & 1), None)
        if found is None:
            continue
        work[pivot], work[found] = work[found], work[pivot]
        for i in range(len(work)):
            if i != pivot and (work[i][0] >> col) & 1:
                work[i] = [a ^ b for a, b in zip(work[i], work[pivot])]
        pivots.append(col)
        pivot += 1
    for row, target, provenance in work:
        if row == 0 and target:
            return {"status": "inconsistent", "row_combination": provenance, "rank": pivot}
    answer = sum(work[i][1] << col for i, col in enumerate(pivots))
    null_basis = []
    for free in sorted(set(range(n)) - set(pivots)):
        v = 1 << free
        for i, col in enumerate(pivots):
            v |= ((work[i][0] >> free) & 1) << col
        null_basis.append(v)
    return {"status": "solved", "solution": answer, "rank": pivot,
            "null_basis": null_basis, "solution_count": 1 << len(null_basis)}


class Source:
    def __init__(self, repo=DEFAULT_REPO):
        self.repo = Path(repo)
        paths = {
            "basis_and_original_cocycle": self.repo / "experiments/theory-contracts/seam-joint-ward-state-20261014/operator_transport_check.py",
            "literal_J": self.repo / "verification/v634_st31_structure.py",
            "existing_exact_checker": self.repo / "experiments/theory-contracts/marked-seam-monodromy-20261023/checker.py",
            "scope": self.repo / "experiments/theory-contracts/marked-seam-monodromy-20261023/THEORY.txt",
        }
        self.sources = {name: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                        for name, path in paths.items()}
        tree = ast.parse(paths["basis_and_original_cocycle"].read_text())
        selected = [x for x in tree.body if isinstance(x, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id in {"e", "basis"} for t in x.targets)]
        fn = [x for x in tree.body if isinstance(x, ast.FunctionDef)
              and x.name in {"eps_car", "split_even_lattice", "eps_e8"}]
        if len(fn) != 3:
            raise ValueError("Original cocycle definitions unavailable")
        env = {"S": S, "q": S.ones(8, 1), "s": S.ones(8, 1)/2}
        exec(compile(ast.Module(body=selected + fn, type_ignores=[]), str(paths["basis_and_original_cocycle"]), "exec"), env)
        self.B = env["basis"]
        self.original_eps = env["eps_e8"]
        jtree = ast.parse(paths["literal_J"].read_text())
        jnodes = [x for x in jtree.body if isinstance(x, ast.FunctionDef) and x.name == "J_vec"]
        jenv = {"np": np}
        exec(compile(ast.Module(body=jnodes, type_ignores=[]), str(paths["literal_J"]), "exec"), jenv)
        self.J = S.Matrix.hstack(*[S.Matrix(jenv["J_vec"](np.array(list(S.eye(8)[:, i]), dtype=object))) for i in range(8)])
        self.T0 = S.diag(1, -1, 1, -1, 1, -1, 1, -1)
        self.matrices = {"J": self.B.inv()*self.J*self.B, "T": self.B.inv()*self.T0*self.B}
        self.G = self.B.T*self.B
        self.units = [tuple(int(x) for x in S.eye(8)[:, i]) for i in range(8)]
        self.images = {k: [tuple(int(x) for x in 2*self.B*m[:, i]) for i in range(8)]
                       for k, m in {"I": S.eye(8), **self.matrices}.items()}
        for name, m in self.matrices.items():
            if any(x.q != 1 for x in m) or m.T*self.G*m != self.G:
                raise ArithmeticError(f"{name} is not an integral lattice isometry")
        if self.matrices["T"]**2 != S.eye(8) or (self.matrices["T"]*self.matrices["J"])**2 != S.eye(8):
            raise ArithmeticError("Original reflected base relation changed")
        self.identity = digest(self.sources)

    @staticmethod
    @lru_cache(None)
    def eps(u, v):
        """Exact doubled-coordinate formula from the existing checker."""
        c, h = u[0] % 2, v[0] % 2
        if any(x % 2 != c for x in u) or any(x % 2 != h for x in v):
            raise ValueError("Mixed coordinate parity is not in E8")
        d = [(x-c)//2 for x in u]
        f = [(x-h)//2 for x in v]
        if sum(d) % 2 or sum(f) % 2:
            raise ValueError("Vector outside E8")
        exponent = c*(sum(f)//2) + c*h
        exponent += sum(d[i]*f[j] for i in range(8) for j in range(i))
        exponent += c*h*sum(i*(d[i]+f[i]) for i in range(8))
        return -1 if exponent % 2 else 1

    @lru_cache(None)
    def ordered(self, n, kind):
        point, phase = (0,)*8, 1
        for b, k in zip(self.images[kind], n):
            unit = self.eps(b, tuple(-x for x in b)) if k < 0 else 1
            if k < 0:
                b = tuple(-x for x in b)
            for _ in range(abs(k)):
                phase *= unit*self.eps(point, b)
                point = tuple(x+y for x, y in zip(point, b))
        return phase

    @lru_cache(None)
    def eta(self, kind, n):
        return self.ordered(n, kind)*self.ordered(n, "I")

    @lru_cache(None)
    def image(self, kind, n):
        return tuple(int(x) for x in self.matrices[kind]*S.Matrix(n))

    def apply(self, kind, char, n):
        out = self.image(kind, n)
        out_mask = sum((x % 2) << i for i, x in enumerate(out))
        sign = self.eta(kind, n)*(-1 if parity(char & out_mask) else 1)
        return sign, out

    def word(self, word, a, b, n):
        sign = 1
        for kind in reversed(word):
            phase, n = self.apply(kind, a if kind == "T" else b, n)
            sign *= phase
        return sign, n

    @lru_cache(None)
    def residual_bits(self, a, b):
        result = []
        for word in [("T", "T"), ("T", "J", "T", "J")]:
            for n in self.units:
                sign, out = self.word(word, a, b, n)
                if out != n:
                    raise ArithmeticError("A sign cannot repair a base-space error")
                result.append(int(sign == -1))
        return tuple(result)

    def task(self, b, fixed=None):
        if type(b) is not int or not 0 <= b < 256:
            raise ValueError("The J character must be an integer in 0..255")
        fixed = {str(k): v for k, v in (fixed or {}).items()}
        if any(not k.isdigit() or not 0 <= int(k) < 8 or type(v) is not int or v not in (0, 1)
               for k, v in fixed.items()):
            raise ValueError("Invalid explicitly declared character marking")
        zero = self.residual_bits(0, b)
        columns = [[x ^ y for x, y in zip(self.residual_bits(1 << j, b), zero)] for j in range(8)]
        rows = [sum(columns[j][i] << j for j in range(8)) for i in range(16)]
        labels = [f"{word}:basis[{i}]" for word in ["T*T=I", "T*J*T*J=I"] for i in range(8)]
        rhs = list(zero)
        for i, value in sorted(fixed.items()):
            rows.append(1 << int(i)); rhs.append(int(value)); labels.append(f"declared_character_mark[{i}]")
        task = {"kind": "E8-reflection-sign-lift", "source_id": self.identity, "j_character": b,
                "rows": rows, "rhs": rhs, "labels": labels,
                "scope": "fixed-J spatial reflection lift relations only; no puncture orbit classification",
                "declared_character_marks": fixed or {}}
        task["id"] = digest(task)
        task["equation_key"] = digest({"rows": rows, "rhs": rhs})
        return task

    def check(self, task, candidate):
        expected_task = self.task(task["j_character"], task["declared_character_marks"])
        if task != expected_task:
            raise ValueError("Task/source identity or exact equations have been modified")
        if not isinstance(candidate, dict):
            return {"status": "invalid", "error": "Expected a JSON object"}
        if candidate.get("status") == "inconsistent":
            rows = candidate.get("rows")
            if not isinstance(rows, list) or not rows or any(type(i) is not int or not 0 <= i < len(task['rows']) for i in rows) or len(rows) != len(set(rows)):
                return {"status": "invalid", "error": "Expected distinct valid equation indices"}
            left = right = 0
            for i in rows:
                left ^= task["rows"][i]; right ^= task["rhs"][i]
            if left == 0 and right == 1:
                return {"status": "verified_inconsistent", "rows": rows, "physical_completion": False}
            return {"status": "counterexample", "relation": "contradiction-certificate", "left": bits(left), "right": right}
        try:
            a = mask(candidate["bits"])
        except (KeyError, ValueError, TypeError) as e:
            return {"status": "invalid", "error": str(e)}
        for i, (row, target) in enumerate(zip(task["rows"], task["rhs"])):
            if parity(row & a) != target:
                return {"status": "counterexample", "relation": task["labels"][i],
                        "row": bits(row), "expected": target, "observed": parity(row & a)}
        # Recompute directly in the original central extension, not through the equations.
        residual = self.residual_bits(a, task["j_character"])
        if any(residual):
            raise ArithmeticError("Linear interface and direct lattice action disagree")
        return {"status": "verified", "bits": bits(a),
                "scope": task["scope"], "physical_completion": False}

    def audit(self):
        rng = np.random.default_rng(20261024)
        cocycles, lift_checks = 0, 0
        for _ in range(30):
            n, m = [tuple(map(int, rng.integers(-2, 3, 8))) for _ in range(2)]
            u, v = self.B*S.Matrix(n), self.B*S.Matrix(m)
            u2, v2 = tuple(int(2*x) for x in u), tuple(int(2*x) for x in v)
            if self.eps(u2, v2) != self.original_eps(u, v):
                raise ArithmeticError("Literal source cocycle mismatch")
            cocycles += 1
            for kind in ["J", "T"]:
                gn, gm = self.image(kind, n), self.image(kind, m)
                gu, gv = self.B*S.Matrix(gn), self.B*S.Matrix(gm)
                left = self.eta(kind, tuple(x+y for x, y in zip(n, m)))*self.eps(u2, v2)
                right = self.eta(kind, n)*self.eta(kind, m)*self.eps(tuple(int(2*x) for x in gu), tuple(int(2*x) for x in gv))
                if left != right:
                    raise ArithmeticError("Lift fails the original product law")
                lift_checks += 1
        counts, unique_equations = {}, set()
        for b in range(256):
            for n in self.units:
                phase, out = self.word(('J',)*4, 0, b, n)
                if phase != 1 or out != n:
                    raise ArithmeticError("Corrected J is not an order-four lattice-field lift")
            task = self.task(b)
            direct = [a for a in range(256) if not any(self.residual_bits(a, b))]
            linear = [a for a in range(256) if all(parity(r & a) == v for r, v in zip(task['rows'], task['rhs']))]
            if direct != linear:
                raise ArithmeticError("Linearized equations changed the actual valid character set")
            solution = exact_gf2(task["rows"], task["rhs"])
            if solution["status"] == "solved":
                if solution["solution_count"] != len(direct):
                    raise ArithmeticError("Elimination disagrees with exhaustive original action")
                self.check(task, {"bits": bits(solution["solution"])})
            else:
                if direct:
                    raise ArithmeticError("False inconsistent certificate")
            key = str(len(direct))
            counts[key] = counts.get(key, 0) + 1
            unique_equations.add(task["equation_key"])
        return {"status": "finite_check_passed", "sources": self.sources, "source_id": self.identity,
                "cocycle_cross_checks": cocycles, "product_law_sample_checks": lift_checks,
                "all_character_pairs_checked": 65536, "j_characters_by_reflection_solution_count": counts,
                "all_256_J_lifts_have_fourth_power_identity": True,
                "linearization_matches_full_direct_solution_sets": True,
                "distinct_unmarked_linear_systems": len(unique_equations),
                "scope": "fixed base T0 and J; all 256x256 character corrections, not full monodromy classification",
                "physical_completion": False}
