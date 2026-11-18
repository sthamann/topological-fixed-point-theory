"""Joint peripheral and compiler-clock transport of actual Phi12 intertwiners."""
import contextlib
import io
import json
from hashlib import sha256
from pathlib import Path
import runpy

import sympy as sp

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "e-series-compiler-closure-20261023/cover_dictionary.py"
digest = sha256(SOURCE.read_bytes()).hexdigest()
with contextlib.redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(SOURCE))


def need(condition, label):
    if not bool(condition):
        raise RuntimeError(label)


A = old["COMP12"]
BE, BC = old["B_NATIVE"], old["B_CURVE"]
J, deck = old["J8"], old["DECK6"]
need(J * BE == BE * A**9, "actual native peripheral J")
need(deck * BC == BC * A**4, "actual oriented Burau peripheral deck")
need(A**12 == sp.eye(4), "original exact clock")
X = sp.symbols("X")
phi12 = X**4 - X**2 + 1
need(sp.factor(A.charpoly(X).as_expr()) == phi12, "common actual Phi12 generator")

# In the two certified integral cyclic bases, the clock-intertwiner space
# is the commutant of A. Its integral basis is I,A,A^2,A^3.
basis = [A**k for k in range(4)]
vecs = sp.Matrix.hstack(*(sp.Matrix(16, 1, list(B)) for B in basis))
unknown = sp.Matrix(4, 4, sp.symbols("h0:16"))
equations = list(unknown * A - A * unknown)
linear, _ = sp.linear_eq_to_matrix(equations, list(unknown))
need(16 - linear.rank() == 4 and vecs.rank() == 4, "complete rational clock Hom space")

# A is cyclic with e0,Ae0,A^2e0,A^3e0 the standard basis. For any integral
# commuting matrix H, H e0 is integral and determines H as that polynomial.
cyclic = sp.Matrix.hstack(*(B[:, 0] for B in basis))
need(cyclic == sp.eye(4), "complete integral commutant is Z[A]")

JL, DB = A**9, A**4
relative = A**(-5)
need(relative == A**7, "relative meridian exponent")
for col, H in enumerate(basis):
    image = DB * H * JL.inv()
    coefficients = relative[:, col]
    reconstructed = sp.zeros(4)
    for coefficient, generator in zip(coefficients, basis):
        reconstructed += coefficient * generator
    need(image == reconstructed, "actual peripheral action on every clock Hom basis vector")

need(sp.factor(relative.charpoly(X).as_expr()) == phi12, "relative faithful Phi12 character")
need((relative - sp.eye(4)).det() != 0, "no nonzero peripheral-flat clock intertwiner")
need(relative**12 == sp.eye(4), "relative return12")
for n in (1, 2, 3, 4, 6):
    need(relative**n != sp.eye(4), "relative return has exact order12")
need(all((relative**n - sp.eye(4)).det() != 0 for n in range(1, 12)),
     "every nonzero rational intertwiner has twelve-element peripheral orbit")
need(relative**4 * A**8 == sp.eye(4), "finite plus infinity sphere relation")
need(A**8 == DB**(-4), "actual infinity transport")
need(sha256(SOURCE.read_bytes()).hexdigest() == digest, "source unchanged")

print(json.dumps({
    "research_id": "MATH.SEAM.JOINT.RELATIVE.05",
    "verdict": "PASS",
    "rational_clock_hom_dimension": 4,
    "integral_clock_hom": "Z[A]=Z[zeta12]",
    "integral_clock_isomorphisms": "Z[zeta12]^units torsor",
    "native_peripheral_exponent": 9,
    "burau_peripheral_exponent": 4,
    "relative_hom_peripheral_exponent": 7,
    "relative_hom_infinity_exponent": 8,
    "flat_clock_hom_dimension": 0,
    "relative_return_order": 12,
    "minimal_trivializing_cover": "the existing C12",
    "scope": "The retained Burau integral model and native Phi12 modules. Clock equality alone is weaker than local-system equality. The relative Hom system is determined by existing monodromies; no line twist is imposed to change either source.",
    "source_sha256": {str(SOURCE.relative_to(HERE.parents[2])): digest},
}, indent=2))

