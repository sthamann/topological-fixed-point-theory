"""Test fixed internal charge clocks against the original whole v252 Higgs family.

The integer Yukawa and Majorana matrices are constructor controls, not a fit or
claimed values of the physical Yukawas. A fixed-H commutator is not a test of
the gauge action that simultaneously transports the Higgs field.
"""
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "verification"))
import v252_full_finite_triple as orig

Z3 = np.zeros((3, 3), complex)
Z2 = np.zeros((2, 2), complex)
# Nonzero values in every Yukawa matrix, with no implied physical selection.
Y = np.array([[1, 2, 3], [2, 5, 7], [3, 7, 11]], complex)
MR = np.array([[1, 2, 0], [2, 3, 1], [0, 1, 4]], complex)
D_y = orig.full_D(Y, Y, Y, Y, Z3)
D_m = orig.full_D(Z3, Z3, Z3, Z3, MR)
K, gamma = orig.Kswap(), orig.grading()
bp = orig.rep(1, Z2, Z3)
Q = [np.array([[0, 1], [-1, 0]], complex),
     np.array([[0, 1j], [1j, 0]], complex),
     np.eye(2, dtype=complex), np.diag([-1j, 1j])]

def fluctuation(q):
    a = orig.rep(0, q, Z3)
    om = a @ (D_y @ bp - bp @ D_y)
    A = om + om.conj().T
    return A + K @ A.conj() @ K

Hs = [fluctuation(q) for q in Q]
def full_charge(fn):
    part = np.array([fn(name, chir, color, gen) for name, chir, color, gen in orig.P])
    # Same particle/antiparticle convention as orig.Kswap and complex charge conjugation.
    return np.diag(np.r_[part, -part]).astype(complex)

def B(n, c, col, g): return 1 if n in ("u", "d") else 0  # 3B
def L(n, c, col, g): return 1 if n in ("nu", "e") else 0
def Qem(n, c, col, g): return {"nu": 0, "e": -3, "u": 2, "d": -1}[n]  # 3Q
def hyper(n, c, col, g):  # 6Y
    if n in ("nu", "e"):
        return -3 if c == "L" else (0 if n == "nu" else -6)
    return 1 if c == "L" else (4 if n == "u" else -2)
def color(n, c, col, g): return [1, -1, 0][col] if n in ("u", "d") else 0
def color2(n, c, col, g): return [0, 1, -1][col] if n in ("u", "d") else 0

charges = {"3B": full_charge(B), "L": full_charge(L),
           "3(B-L)": full_charge(lambda n,c,k,g: B(n,c,k,g)-3*L(n,c,k,g)),
           "3Qem": full_charge(Qem), "6Y": full_charge(hyper),
           "color_cartan": full_charge(color), "color_cartan_2": full_charge(color2)}
def nnz(A): return int(np.count_nonzero(A))
def comm(A, B): return A @ B - B @ A

assert all(np.array_equal(H, H.conj().T) for H in Hs)
assert all(nnz(H @ gamma + gamma @ H) == 0 for H in Hs)
assert nnz(D_m @ gamma + gamma @ D_m) == 0
for name, lam in charges.items():
    assert nnz(comm(lam, gamma)) == 0
    assert np.all(np.equal(np.diag(lam).real, np.rint(np.diag(lam).real)))
    print(name, "gamma_trace", np.trace(lam @ gamma).real,
          "Yukawa", nnz(comm(lam, D_y)),
          "Higgs4", [nnz(comm(lam, H)) for H in Hs],
          "Majorana", nnz(comm(lam, D_m)))

for name in ("3B", "color_cartan", "color_cartan_2"):
    lam=charges[name]
    assert nnz(comm(lam, D_y)) == 0
    assert all(nnz(comm(lam, H)) == 0 for H in Hs)
    assert nnz(comm(lam, D_m)) == 0
    assert np.trace(lam @ gamma) == 0
for name in ("L", "3(B-L)"):
    assert all(nnz(comm(charges[name], H)) == 0 for H in Hs)
    assert nnz(comm(charges[name], D_m)) > 0
assert nnz(comm(charges["3Qem"], D_y)) == 0
assert nnz(comm(charges["3Qem"], Hs[0])) > 0
assert nnz(comm(charges["6Y"], D_y)) > 0
assert nnz(comm(charges["6Y"], D_m)) == 0
# The reference Higgs direction is the original supplied Yukawa block.  In the
# five-dimensional named span the three rows below are triangular: the neutral
# Yukawa link fixes d=0 in a*3B+b*L+c*3Q+d*6Y+e*Cartan; a charged Higgs link
# then fixes c=0; the neutrino Majorana link fixes b=0.  B and Cartan remain.
assert np.array_equal(Hs[2], D_y)
for name in ("3B", "L", "3Qem", "color_cartan", "color_cartan_2"):
    assert nnz(comm(charges[name], D_y)) == 0
for name in ("3B", "L", "6Y", "color_cartan", "color_cartan_2"):
    if name != "6Y":
        assert nnz(comm(charges[name], Hs[0])) == 0
for name in ("3B", "3Qem", "6Y", "color_cartan", "color_cartan_2"):
    assert nnz(comm(charges[name], D_m)) == 0
print("PASS: exact integer matrix comparisons on the 96-dimensional constructor control")
