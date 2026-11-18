"""Replay the existing six induced insertions, not a new source search."""
from pathlib import Path
import json
import mpmath as mp
import sympy as sp

mp.mp.dps = 65
om = mp.exp(2j * mp.pi / 3)
ii = mp.j
hs = [-om, -om, om]
C = 2 * mp.power(3, -mp.mpf(3) / 4) * mp.exp(-5j * mp.pi / 12)
x = mp.root(3, 4) * mp.cos(5 * mp.pi / 12)
aa = mp.atan2(2 * x, 1 - mp.sqrt(3))
bb = mp.atan2(mp.sqrt(2) * x, 1 - mp.sqrt(3) / 2)
ell = aa - mp.sqrt(2) * bb
I02 = C * ell
I21 = -C * (ell + mp.pi * (mp.sqrt(2) - 1))
P = mp.matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
F = mp.matrix([[0, 0, I02], [0, 0, 0], [0, I21, 0]])
N = mp.matrix([[0, 0, 0], [0, 0, -ii * I02], [-ii * I21, 0, 0]])

def fmt(z):
    return {"re": mp.nstr(mp.re(z), 60), "im": mp.nstr(mp.im(z), 60)}

def matrix_fmt(m):
    return [[fmt(m[r, s]) for s in range(3)] for r in range(3)]

def maxnorm(m):
    return max(abs(v) for v in m)

records = []
all_error = mp.mpf(0)
for a, b in [(0, 1), (1, 2), (2, 0), (1, 0), (2, 1), (0, 2)]:
    m = mp.matrix(3, 3)
    entries = []
    for r in range(3):
        for s in range(3):
            t = (r + 1) % 3
            entry = {"r": r, "s": s, "q": s, "p": t}
            if r != s and r != (s + 1) % 3:
                entry["reason"] = "disjoint open supports"
            else:
                eq = -1 if s in (a, (a + 1) % 3) else 1
                ep = -1 if t in (b, (b + 1) % 3) else 1
                entry["kummer_eta_sign"] = eq * ep
                if eq * ep == 1:
                    entry["reason"] = "two eta arcs cancel"
                else:
                    th = ((5 if r == s else 7) * mp.pi / 6
                          + 2 * mp.pi * s / 3)
                    ray = mp.exp(ii * th)
                    z0 = -hs[(a - s) % 3] * hs[(b - t) % 3] / om**(a + b)
                    # Each individual square root continues from 1 on this ray.
                    # The whole squared differential must not be principal-rooted.
                    def integrand(v):
                        w = ray * v
                        zeta = z0 * mp.sqrt(1 + w / om**((a + 2) % 3))
                        zeta *= mp.sqrt(1 - w / om**((b + 2) % 3))
                        return zeta * ray / (1 - w**3)
                    period = mp.quad(integrand, [0, mp.mpf('0.25'), 1, 4, mp.inf])
                    m[r, s] = (-1 if r == s else 1) * period
                    entry.update({"reason": "nonzero signed continuous-root period",
                                  "ray_angle_pi": mp.nstr(th / mp.pi, 30),
                                  "zeta_at_0": fmt(z0), "period": fmt(period)})
            entry["scalar_reader_entry_without_4X_over_pi"] = fmt(m[r, s])
            entries.append(entry)
    common = P**(1 - b) * m * P**a
    k = a if b == (a + 1) % 3 else (a - 1) % 3
    expected = om**(-k) * (F if b == (a + 1) % 3 else N)
    error = maxnorm(common - expected)
    all_error = max(all_error, error)
    records.append({"a": a, "b": b, "entries": entries,
                    "scalar_cable_matrix": matrix_fmt(m),
                    "common_G_matrix": matrix_fmt(common),
                    "expected_phase": f"omega^(-{k})",
                    "expected_base": "F" if b == (a + 1) % 3 else "N",
                    "max_abs_residual": mp.nstr(error, 10)})

# Independent exact algebra, with coefficients kept symbolic.
u, v, w = sp.symbols("alpha beta w", nonzero=True)
os = -sp.Rational(1, 2) + sp.sqrt(3) * sp.I / 2
ps = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
def e(r, s):
    m = sp.zeros(3)
    m[r, s] = 1
    return m
fs = u * e(0, 2) + v * e(2, 1)
ns = -sp.I * (u * e(1, 2) + v * e(2, 0))
def zero(m):
    return all(sp.simplify(z) == 0 for z in m)
exact = {
    "F_squared": zero(fs**2 - u*v*e(0, 1)),
    "N_squared": zero(ns**2 + u*v*e(1, 0)),
    "E00_product": zero((fs**2/(u*v))*(-ns**2/(u*v)) - e(0, 0)),
    "E11_product": zero((-ns**2/(u*v))*(fs**2/(u*v)) - e(1, 1)),
    "E02_recovery": zero(e(0, 0)*fs*e(2, 2)/u-e(0, 2)),
    "E21_recovery": zero(e(2, 2)*fs*e(1, 1)/v-e(2, 1)),
    "E12_recovery": zero(sp.I*e(1, 1)*ns*e(2, 2)/u-e(1, 2)),
    "E20_recovery": zero(sp.I*e(2, 2)*ns*e(0, 0)/v-e(2, 0)),
}
for a, b in [(0, 1), (1, 2), (2, 0), (1, 0), (2, 1), (0, 2)]:
    z2 = (1+w**3)/(1-w**3)
    alpha2 = 1/((w+os**a)*(w+os**((a+1)%3))
                *(w-os**b)*(w-os**((b+1)%3)))
    proposed = -(w+os**((a+2)%3))*(w-os**((b+2)%3))/(1-w**3)**2
    exact[f"squared_differential_{a}{b}"] = sp.simplify(z2*alpha2-proposed) == 0
assert all(exact.values())
assert all_error < mp.mpf('1e-55')
out = {"scope": "Existing NS Phi=0 Hardy/G/J reader; all six existing offdiagonal induced mu inputs. Not whole physical Yukawa/source selection.",
       "normalization": "Multiply every reported scalar family matrix by 4 Xstar/pi and tensor P_gamma,+; Xstar=1/(2 cosh(ell_NS/2)).",
       "frame": "G_a=G Pc^(-a), JE^(b-1) output. Common=P^(1-b) Scalar P^a.",
       "working_decimal_digits": mp.mp.dps,
       "I02": fmt(I02), "I21": fmt(I21),
       "I21_over_I02": fmt(I21/I02),
       "all_54_entries": records,
       "total_entries": 54, "nonzero_entries": 12, "exact_zeros": 42,
       "max_common_frame_residual": mp.nstr(all_error, 15),
       "exact_symbolic_identities": exact,
       "algebra_scope": "Unital coefficient algebra generated by F,N is M3(C). This is not a CAR-product or physical-insertion identification."}
target = Path(__file__).with_name("HARDY-BMU-ALL-SIX-INPUT-READER.json")
target.write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps({"output": str(target), "total_entries": 54, "nonzero": 12,
                  "max_abs_residual": mp.nstr(all_error, 10), "exact_identities": exact}, indent=2))
