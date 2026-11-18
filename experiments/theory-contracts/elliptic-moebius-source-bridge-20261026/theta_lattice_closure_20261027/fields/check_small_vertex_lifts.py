"""Exact standard-cocycle lift of the marked order-24 lattice subgroup.

The lattice basis, Gram matrix and marked J,V,sigma matrices are imported from
the already checked theta-lattice construction.  The only freedom enumerated
here is multiplication of each standard lattice-isometry lift by one of the
256 sign characters of Lambda/2 Lambda.
"""
from __future__ import annotations

import contextlib
import io
import itertools
import json
import runpy
from pathlib import Path

import numpy as np


def require(test, message):
    if not test:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
LATTICE_CHECK = HERE.parent / "lattice" / "check_theta_lattice.py"
with contextlib.redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(LATTICE_CHECK))

s = source["s"]
B = source["B"]
Gram_exact = source["Gram"]
real_matrix = source["real_matrix"]


def lattice_matrix(g):
    A = s.simplify(B.inv() * real_matrix(g) * B)
    require(all(x.is_Integer for x in A), "marked map is not integral")
    require(A.T * Gram_exact * A == Gram_exact, "marked map is not an isometry")
    return np.array(A.tolist(), dtype=np.int64)


Amat_int = {
    "J": lattice_matrix(s.I * s.eye(4)),
    "V": lattice_matrix(source["V"]),
    "sigma": lattice_matrix(source["sigma"]),
}
E = {name: s.Matrix(A.tolist()) for name, A in Amat_int.items()}
I8_exact = s.eye(8)
require(E["J"]**4 == I8_exact, "exact J^4")
require(E["V"]**2 == I8_exact, "exact V^2")
require(E["V"] * E["J"] == E["J"] * E["V"], "exact [V,J]")
require(E["sigma"]**3 == I8_exact, "exact sigma^3")
require(E["sigma"] * E["J"] == E["J"] * E["sigma"], "exact [sigma,J]")
require(E["V"] * E["sigma"] * E["V"]**-1 == E["sigma"]**-1,
        "exact V sigma V^-1")
Amat = {name: A % 2 for name, A in Amat_int.items()}
I8 = np.eye(8, dtype=np.uint8)

# Standard basis cocycle epsilon(m,n)=(-1)^(m^T F n).  The diagonal choice
# makes epsilon(m,m)=(-1)^(<m,m>/2), while F+F^T=Gram modulo two.
F = np.zeros((8, 8), dtype=np.uint8)
for i in range(8):
    F[i, i] = (int(Gram_exact[i, i]) // 2) % 2
    for j in range(i):
        F[i, j] = int(Gram_exact[i, j]) % 2
Gram2 = np.array(Gram_exact.tolist(), dtype=np.int64) % 2
require(np.array_equal((F + F.T) % 2, Gram2), "cocycle commutator")

vectors = np.array(
    [[(mask >> i) & 1 for i in range(8)] for mask in range(256)],
    dtype=np.uint8,
)


def mask_of(x):
    return sum(int(x[i]) << i for i in range(8))


def standard_q(A):
    """Zero-on-basis quadratic correction for the fixed cocycle."""
    D = (A.T @ F @ A + F) % 2
    require(np.array_equal(D, D.T), "cocycle correction symmetric")
    require(not np.any(np.diag(D)), "cocycle correction alternating")
    q = np.zeros(256, dtype=np.uint8)
    for mask, x in enumerate(vectors):
        value = 0
        for i in range(8):
            for j in range(i + 1, 8):
                value ^= int(D[i, j] & x[i] & x[j])
        q[mask] = value
    # Exhaustive polarization identity on Lambda/2Lambda.
    for xmask, x in enumerate(vectors):
        for ymask, y in enumerate(vectors):
            zmask = mask_of(x ^ y)
            lhs = int(q[zmask] ^ q[xmask] ^ q[ymask])
            rhs = int((x @ D @ y) % 2)
            require(lhs == rhs, "quadratic cocycle correction")
    return q


qzero = {name: standard_q(A) for name, A in Amat.items()}


def q_with_character(name, bits):
    t = np.array([(bits >> i) & 1 for i in range(8)], dtype=np.uint8)
    return qzero[name] ^ ((vectors @ t) % 2)


# A lift is (A,q), acting by e_m -> (-1)^q(m) e_(A m).
def compose(first, second):
    """Return first after second."""
    A1, q1 = first
    A2, q2 = second
    image_indices = np.array([mask_of((A2 @ x) % 2) for x in vectors])
    return ((A1 @ A2) % 2, q2 ^ q1[image_indices])


IDENTITY = (I8, np.zeros(256, dtype=np.uint8))


def power(op, exponent):
    result = IDENTITY
    for _ in range(exponent):
        result = compose(op, result)
    return result


def inverse(op):
    A, q = op
    # All matrices here are invertible over F2.
    Ainv_exact = s.Matrix(A.tolist()).inv()
    require(all(x.is_Integer for x in Ainv_exact), "inverse matrix")
    Ainv = np.array(Ainv_exact.tolist(), dtype=np.int64) % 2
    image_indices = np.array([mask_of((Ainv @ x) % 2) for x in vectors])
    return (Ainv, q[image_indices])


def equal(left, right):
    return np.array_equal(left[0], right[0]) and np.array_equal(left[1], right[1])


def operator_key(op):
    A, q = op
    return (tuple(int(x) for x in A.ravel()), tuple(int(x) for x in q))


def candidate_lifts(name, order):
    candidates = []
    for bits in range(256):
        op = (Amat[name], q_with_character(name, bits))
        if equal(power(op, order), IDENTITY):
            candidates.append((bits, op))
    return candidates


candidates = {
    "J": candidate_lifts("J", 4),
    "V": candidate_lifts("V", 2),
    "sigma": candidate_lifts("sigma", 3),
}
require([len(candidates[x]) for x in ("J", "V", "sigma")] == [256, 16, 16],
        "unexpected power-normalized lift counts")

solutions = []
for jbits, Jhat in candidates["J"]:
    for vbits, Vhat in candidates["V"]:
        if not equal(compose(Jhat, Vhat), compose(Vhat, Jhat)):
            continue
        for sbits, Shat in candidates["sigma"]:
            if not equal(compose(Shat, Jhat), compose(Jhat, Shat)):
                continue
            if not equal(compose(compose(Vhat, Shat), Vhat), power(Shat, 2)):
                continue
            solutions.append((jbits, vbits, sbits, Jhat, Vhat, Shat))

require(len(solutions) == 256, "strict marked lift census")
require(not any(jbits == 0 for jbits, *_ in solutions),
        "zero-on-basis J lift unexpectedly extends strictly")

# Lexicographically first exact section of C4 x S3.
jbits, vbits, sbits, Jhat, Vhat, Shat = solutions[0]
require((jbits, vbits, sbits) == (0x03, 0x35, 0x2A), "selected solution changed")

# Check the complete presentation on all 256 residue fields.
require(equal(power(Jhat, 4), IDENTITY), "J^4")
require(equal(power(Vhat, 2), IDENTITY), "V^2")
require(equal(compose(Jhat, Vhat), compose(Vhat, Jhat)), "[J,V]")
require(equal(power(Shat, 3), IDENTITY), "sigma^3")
require(equal(compose(Shat, Jhat), compose(Jhat, Shat)), "[sigma,J]")
require(equal(compose(compose(Vhat, Shat), Vhat), power(Shat, 2)),
        "V sigma V = sigma^-1")

# Closure is exactly the geometric group C4 x S3 of order 24.
generators = [Jhat, Vhat, Shat]
group = {operator_key(IDENTITY): IDENTITY}
frontier = [IDENTITY]
while frontier:
    old = frontier.pop()
    for generator in generators:
        new = compose(generator, old)
        key = operator_key(new)
        if key not in group:
            group[key] = new
            frontier.append(new)
require(len(group) == 24, "lifted geometry group order")

# The native clock c=J sigma is an honest order-12 vertex automorphism.
chat = compose(Jhat, Shat)
require(equal(power(chat, 12), IDENTITY), "c^12")
require(not any(equal(power(chat, n), IDENTITY) for n in (1, 2, 3, 4, 6)),
        "c exact order")
require(equal(power(chat, 4), Shat), "c^4=sigma")
require(equal(power(chat, 9), Jhat), "c^9=J")

# The weaker two-object construction needs no split of V.  Take the
# zero-on-basis cocycle lift T of the lattice isomorphism V between the two
# branch coordinate frames.  It deliberately fails the same-frame involution
# and commutation requirements.  Transport the second object's marks through
# T and define the return branch by beta_0 T^-1.
T = (Amat["V"], qzero["V"])
Tinv = inverse(T)
require(equal(compose(T, Tinv), IDENTITY) and equal(compose(Tinv, T), IDENTITY),
        "T inverse")
require(not equal(power(T, 2), IDENTITY), "T unexpectedly an involution")
require(not equal(compose(T, Jhat), compose(Jhat, T)),
        "T unexpectedly commutes with beta zero")

beta0, sigma0 = Jhat, Shat
beta1 = compose(compose(T, beta0), Tinv)
sigma1 = compose(compose(T, power(sigma0, 2)), Tinv)
require(equal(power(beta1, 4), IDENTITY), "transported beta order")
require(equal(power(sigma1, 3), IDENTITY), "transported sigma order")
require(equal(compose(beta1, sigma1), compose(sigma1, beta1)),
        "transported marks commute")

alpha_01 = T
alpha_10 = compose(beta0, Tinv)
require(equal(compose(alpha_10, alpha_01), Jhat), "alpha^2 on object zero")
require(equal(compose(alpha_01, alpha_10), beta1), "alpha^2 on object one")
require(equal(power(beta0, 4), IDENTITY) and equal(power(beta1, 4), IDENTITY),
        "alpha^8")
require(equal(compose(alpha_01, sigma0), compose(power(sigma1, 2), alpha_01)),
        "alpha cubic intertwining 0 to 1")
require(equal(compose(alpha_10, sigma1), compose(power(sigma0, 2), alpha_10)),
        "alpha cubic intertwining 1 to 0")

c1hat = compose(beta1, sigma1)
require(equal(power(c1hat, 12), IDENTITY), "branch one c^12")
require(equal(power(c1hat, 4), sigma1), "branch one c^4")
require(equal(power(c1hat, 9), beta1), "branch one c^9")

out = {
    "verdict": "TWO_OBJECT_VERTEX_GROUPOID_EXISTS_AND_STRICT_ORDER24_SECTION_IS_OPTIONAL",
    "source_checker": str(LATTICE_CHECK),
    "standard_cocycle": "epsilon(m,n)=(-1)^(m^T F n), lower-triangular Gram section",
    "characters_enumerated_per_generator": 256,
    "power_normalized_candidate_counts": {
        "J_order_4": len(candidates["J"]),
        "V_order_2": len(candidates["V"]),
        "sigma_order_3": len(candidates["sigma"]),
    },
    "strict_solution_count": len(solutions),
    "selected_character_masks": {
        "J": f"0x{jbits:02x}",
        "V": f"0x{vbits:02x}",
        "sigma": f"0x{sbits:02x}",
    },
    "selected_character_supports_zero_based": {
        "J": [i for i in range(8) if (jbits >> i) & 1],
        "V": [i for i in range(8) if (vbits >> i) & 1],
        "sigma": [i for i in range(8) if (sbits >> i) & 1],
    },
    "relations": [
        "J^4=1", "V^2=1", "[V,J]=1", "sigma^3=1",
        "[sigma,J]=1", "V sigma V^-1=sigma^-1",
    ],
    "lifted_group_order": len(group),
    "native_clock": ["order(c)=12", "c^4=sigma", "c^9=J"],
    "two_object_alpha": [
        "alpha_01=T", "alpha_10=beta0 T^-1",
        "alpha^2=diag(beta0,beta1)", "alpha^8=1",
        "sigma1=T sigma0^-1 T^-1",
    ],
    "groupoid_does_not_use": ["T^2=1", "[T,beta0]=1", "full G31 splitting"],
    "explicit_T_square_is_identity": False,
    "explicit_T_commutes_with_beta0": False,
    "canonical_zero_on_basis_J_extends_strictly": False,
    "scope": (
        "standard cocycle on the actual theta E8 basis; proves the marked order-24 subgroup "
        "and two-object alpha groupoid only, not a splitting of the full G31 extension"
    ),
}
print(json.dumps(out, ensure_ascii=False, sort_keys=True, indent=2))
