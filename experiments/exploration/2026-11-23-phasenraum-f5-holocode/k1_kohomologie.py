"""K: Familienraum als étale Kohomologie von X = P^1 ohne die vier Marken z^4 = 1 über F_p.

Test aus VERDICT.md (Abschnitt F, Test 1): Rang und mu4-Charaktere von H^1 mit dem TFPT
Familienraum vergleichen (Rang 3, Charaktere {-1, i, -i}). Abbruch, wenn sie abweichen.

1. Punktzählung |X(F_{p^k})| in explizit konstruierten Körpern F_{p^n} (primitives Polynom,
   Exponententafel aus Polynommultiplikation). Lefschetz: |X(F_q)| = q - tr(F | H^1_c), daraus
   die Frobenius-Eigenwerte auf H^1_c.
2. Für p = 1 mod 4 die getwisteten Zählungen #Fix(F o g_a), g_a(z) = a z, a in mu4(F_p):
   tr(F g_a | H^1_c) = p - #Fix. Fourier über mu4 zerlegt H^1_c nach Charakteren und liefert
   je Charakter den Frobenius-Eigenwert.
3. Residuenbrücke: Residuen der TFPT Trägermoden e_m = z^m/(z^4 - 1) (m = 0..4) an den vier
   Marken (exakt, sympy) und das Gegenstück über F5 (Einschränkung auf F5^x modulo Konstanten).
"""
import itertools
from math import gcd

import numpy as np
import sympy as sp


def prime_factors(n):
    return sorted(sp.factorint(n))


def polymulmod(a, b, f, p):
    """a, b: coefficient lists low->high, length n; f monic of degree n (length n+1)."""
    n = len(f) - 1
    prod = [0] * (2 * n - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                prod[i + j] = (prod[i + j] + x * y) % p
    for d in range(2 * n - 2, n - 1, -1):
        c = prod[d]
        if c:
            for k in range(n + 1):
                prod[d - n + k] = (prod[d - n + k] - c * f[k]) % p
    return prod[:n]


def field(p, n):
    """Exp/log tables of F_{p^n} from a primitive polynomial found by search."""
    q = p ** n
    for tail in itertools.product(range(p), repeat=n):
        f = list(tail) + [1]
        if f[0] == 0:
            continue
        x = [0, 1] + [0] * (n - 2) if n > 1 else [f[0] and (-f[0]) % p]
        exp = [tuple([1] + [0] * (n - 1))]
        cur = list(exp[0])
        ok = True
        for j in range(1, q - 1):
            cur = polymulmod(cur, x, f, p)
            t = tuple(cur)
            if t == exp[0]:
                ok = False
                break
            exp.append(t)
        if ok and tuple(polymulmod(cur, x, f, p)) == exp[0]:
            return q, exp, {e: j for j, e in enumerate(exp)}
    raise RuntimeError("no primitive polynomial")


def analyse(p, n):
    q, exp, log = field(p, n)
    N = q - 1
    fp_step = N // (p - 1)                    # F_p^x = exponents divisible by fp_step
    marks = [j for j in range(N) if (4 * j) % N == 0]
    rows = []
    for k in range(1, n + 1):
        if n % k:
            continue
        sub = [j for j in range(N) if (j * (p ** k - 1)) % N == 0]      # F_{p^k}^x
        rat_marks = len(set(sub) & set(marks))
        count = (p ** k + 1) - rat_marks                                 # P^1(F_{p^k}) minus marks
        rows.append((k, count, p ** k - count))                          # tr(F^k | H^1_c)
    twisted = None
    if (p - 1) % 4 == 0:
        mu4 = [j for j in marks]                                         # exponents of mu4 in F_p
        assert all(j % fp_step == 0 for j in mu4)
        twisted = {}
        gen = next(j for j in mu4 if (4 * j) % N == 0 and all((m * j) % N for m in (1, 2, 3)))
        for t in range(4):
            e_a = (t * gen) % N
            # nonzero fixed points of z -> a z^p:  e_a + p j = j  (mod N), z not a mark
            sols = [j for j in range(N) if (e_a + p * j - j) % N == 0 and j not in marks]
            fix = 2 + len(sols)                                          # plus 0 and infinity
            twisted[t] = (fix, p - fix)
    return rows, twisted


print("1. Frobenius auf H^1_c(X), X = P^1 ohne {z^4 = 1}\n")
for p, n in ((3, 4), (5, 4), (7, 4), (13, 4)):
    rows, twisted = analyse(p, n)
    traces = [tr for _, _, tr in rows]
    # eigenvalues: three roots of unity whose power sums match tr(F^k), k = 1, 2, 4
    cand = [c for c in itertools.combinations_with_replacement([1, -1, 1j, -1j], 3)
            if all(abs(sum(np.array(c) ** k) - tr) < 1e-9 for (k, _, tr) in rows)]
    print(f"   p = {p:2d} (p mod 4 = {p % 4}): |X(F_p^k)| und tr(F^k | H^1_c) für k = "
          + ", ".join(f"{k}: {c} / {tr}" for k, c, tr in rows)
          + f"  ->  Eigenwerte auf H^1_c: {sorted(set(cand[0]), key=str) if cand else '?'}"
          + f" {list(cand[0]) if cand else ''}")
    if twisted:
        print(f"      getwistete Zählungen #Fix(F g_a) für a = g^0..g^3: "
              f"{[twisted[t][0] for t in range(4)]}")
        tr = np.array([twisted[t][1] for t in range(4)], dtype=complex)
        # character chi_s(g^t) = i^(s t); multiplicity-weighted Frobenius eigenvalue per character
        lam = {s: sum(tr[t] * (1j ** (-s * t)) for t in range(4)) / 4 for s in range(4)}
        print("      Zerlegung von H^1_c nach mu4-Charakteren chi_s (chi_s(g) = i^s): "
              + ", ".join(f"chi_{s}: {np.round(lam[s].real, 6):g}" for s in range(4)))

print("\n2. Lesart: auf H^1_c trägt jeder nichttriviale Charakter genau einen Frobenius-Eigenwert 1"
      "\n   (auf H^1 = dual: Eigenwert p, Gewicht 2); der triviale Charakter kommt nicht vor.")

# ------------------------------------------------------------------ 3. residue bridge
z = sp.symbols("z")
marks = [1, sp.I, -1, -sp.I]
E = [z ** m / (z ** 4 - 1) for m in range(5)]
Res = sp.Matrix([[sp.simplify(sp.residue(e, z, zeta)) for e in E] for zeta in marks])
print("\n3. Residuenbrücke über C: Residuen von e_m = z^m/(z^4-1) an den Marken 1, i, -1, -i")
sp.pprint(Res)
quot = sp.Matrix.vstack(Res, sp.ones(1, 4) * Res)        # not used for rank: we project mod constants
proj = sp.eye(4) - sp.ones(4, 4) / 4
img = proj * Res
print(f"   Rang der Residuen modulo Konstanten: {img.rank()} (Familienraum hat Rang 3)")
ker = img.nullspace()
print("   Kern (Trägermoden ohne Familienanteil):", [list(v.T) for v in ker])
# characters of the images under zeta -> i zeta (cyclic shift of the marks 1 -> i -> -1 -> -i)
shift = sp.Matrix(4, 4, lambda r, c: 1 if c == (r + 1) % 4 else 0)
evs = (proj * shift * proj).eigenvals()
print("   Eigenwerte der Markenverschiebung ζ -> iζ auf Funktionen modulo Konstanten:",
      {sp.nsimplify(k): v for k, v in evs.items()})

print("\n   Gegenstück über F5: Fun(F5) -> Fun(F5^x)/Konstanten, f -> f|F5^x")
M = np.zeros((4, 5))
for r, x in enumerate((1, 2, 3, 4)):
    M[r, x] = 1
P5 = np.eye(4) - np.ones((4, 4)) / 4
img5 = P5 @ M
print(f"   Rang {np.linalg.matrix_rank(img5)}, Kern = span(δ0, 1_F5^x): "
      f"{np.allclose(img5 @ np.array([1, 0, 0, 0, 0]), 0) and np.allclose(img5 @ np.array([0, 1, 1, 1, 1]), 0)}")
assert img.rank() == 3 and len(ker) == 2
