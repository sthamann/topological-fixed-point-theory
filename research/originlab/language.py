"""OriginDSL for finite typed rational processes.

A process is a matrix from its input boundary to its output boundary.
There is no physical interpretation of a dimension or iteration counter.
"""
from dataclasses import dataclass
from functools import lru_cache
import sympy as sp

ARITY = dict(compose=2, glue=2, mirror=1, reverse=1, branch=1, identify=2,
             trace=1, iterate=1, fix=1, lift=1, atom=0)
PARAMETER = {"branch", "iterate", "lift"}


class TypeMismatch(ValueError):
    pass


@dataclass(frozen=True)
class Term:
    op: str
    args: tuple = ()
    parameter: int | str | None = None

    @classmethod
    def parse(cls, data, depth=0):
        if depth > 20 or not isinstance(data, list) or not data:
            raise ValueError("invalid or overdeep OriginDSL term")
        op, *items = data
        if not isinstance(op, str) or op not in ARITY:
            raise ValueError("unknown OriginDSL operation")
        param = None
        if op in PARAMETER:
            if not items or type(items[0]) is not int or not 1 <= items[0] <= 8:
                raise ValueError("bounded positive integer parameter required")
            param, *items = items
        elif op == "atom":
            if len(items) != 1 or not isinstance(items[0], str):
                raise ValueError("atom must name a declared primitive")
            param, *items = items
        if len(items) != ARITY[op]:
            raise ValueError("wrong OriginDSL arity")
        return cls(op, tuple(cls.parse(x, depth+1) for x in items), param)

    def data(self):
        return [self.op, *([] if self.parameter is None else [self.parameter]),
                *(x.data() for x in self.args)]

    @property
    def size(self):
        return 1 + sum(x.size for x in self.args)


class Evaluator:
    def __init__(self, atoms, max_dimension=24, boundaries=None, seams=None):
        self.atoms = atoms
        self.max_dimension = max_dimension
        self.boundaries = boundaries
        self.seams = seams or {}

    @lru_cache(maxsize=16000)
    def boundary(self, term):
        if self.boundaries is None:
            return None, None
        op = term.op
        if op == "atom":
            if term.parameter not in self.boundaries:
                raise TypeMismatch("unknown primitive")
            return self.boundaries[term.parameter]
        before, after = self.boundary(term.args[0])
        if op in {"compose", "glue"}:
            return self.boundary(term.args[1])[0], after
        if op == "mirror":
            return after, before
        if op == "trace":
            return ("trace",), ("trace",)
        if op in {"branch", "lift"}:
            return (op,term.parameter,before), (op,term.parameter,after)
        return before, after

    @lru_cache(maxsize=16000)
    def evaluate(self, term):
        op = term.op
        if op == "atom":
            if term.parameter not in self.atoms:
                raise TypeMismatch("unknown primitive")
            atom = self.atoms[term.parameter]
            if max(atom.shape)>self.max_dimension:
                raise TypeMismatch("primitive exceeds boundary dimension")
            return atom
        matrices = [self.evaluate(arg) for arg in term.args]
        a = matrices[0]
        b = matrices[-1]
        p = term.parameter
        input_a,output_a = self.boundary(term.args[0])
        input_b,output_b = self.boundary(term.args[-1])
        if op in {"compose", "glue"}:
            seam = self.seams.get((output_b,input_a)) if op == "glue" else None
            if seam is not None:
                if seam.shape != (a.cols,b.rows):
                    raise TypeMismatch("declared seam has incompatible dimensions")
                return a*seam*b
            if a.cols != b.rows or input_a != output_b:
                raise TypeMismatch("incompatible boundaries")
            return a*b
        if op == "mirror":
            return a.H
        if op == "reverse":
            # Reflection of both ordered boundary frames.
            return a[::-1, ::-1]
        if op == "identify":
            if (input_a,output_a)!=(input_b,output_b) or a != b:
                raise TypeMismatch("exact equality in the declared boundary frames failed")
            return a
        if op == "trace":
            if a.rows != a.cols or input_a != output_a:
                raise TypeMismatch("trace requires a closed boundary")
            return sp.ImmutableMatrix([[a.trace()]])
        if op == "iterate":
            if a.rows != a.cols or input_a != output_a:
                raise TypeMismatch("iteration requires an endomorphism")
            return a**p
        if op == "fix":
            if a.rows != a.cols or input_a != output_a:
                raise TypeMismatch("fix requires an endomorphism")
            basis = (a-sp.eye(a.rows)).nullspace()
            if not basis:
                return sp.ImmutableMatrix.zeros(a.rows)
            v = sp.Matrix.hstack(*basis)
            return sp.ImmutableMatrix(v*(v.H*v).inv()*v.H)
        if op in {"branch", "lift"}:
            if max(a.shape)*p > self.max_dimension:
                raise TypeMismatch("boundary dimension exceeds search class")
            unit = sp.eye(p) if op == "branch" else sp.ones(p)/p
            return sp.ImmutableMatrix(sp.kronecker_product(unit, a))
        raise ValueError("unknown operation")


def positive(matrix):
    """Exact Hermitian PSD test by Schur elimination, including zero pivots."""
    if matrix.rows != matrix.cols or matrix != matrix.H:
        return False
    remaining = sp.Matrix(matrix)
    while remaining.rows:
        pivot = remaining[0, 0]
        column = remaining[1:, 0]
        if pivot < 0 or (pivot == 0 and any(column)):
            return False
        remainder = remaining[1:, 1:]
        remaining = remainder if pivot == 0 else remainder-column*column.H/pivot
    return True
