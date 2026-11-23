"""F: Familientest für den Weil-Kandidaten W.

Teil A (angekündigter Test, Abbruchregel aus VERDICT.md): Weil-Darstellung von Sp(2n, 5) auf
C^(5^n) für n = 1, 2, 3, eingeschränkt auf das Howe-Paar (diagonales SL(2,5) = 2I) x O(n, 5).
Wie oft kommt das fermionische Dublett 2 (und sein Galoispartner 2') vor, und welche
Darstellung sigma von O(n,5) trägt seinen Multiplizitätsraum? Gesucht: drei Kopien oder ein
natürliches Triplett.

Teil B (Wörterbuch über F5): vier Marken = Nullstellen von z^4 - 1 = ganz F5^x; fünf Plätze =
Polynome vom Grad <= 4 = alle Funktionen auf F5; Dilatationscharaktere auf C^F5; Familien als
nichttriviale Charaktere von F5^x; Auswahl "Einheitengruppe = mu4"; Gauss-Summen.

Teil C (Hauptserie): Für jeden Charakter chi von F5^x die induzierte Darstellung von SL(2,5)
auf den sechs Punkten von P^1(F5): Dimensionen, Irreduzibilität, fermionisch oder bosonisch.
"""
import itertools
import numpy as np

exec(open(__file__.replace("f1_familien.py", "w1_weil_darstellung.py")).read().split("def main")[0])

p = 5
S, T, G = weil(p)
Z = central_minus_one(G, [S, T], p)
Vb, Vf = block_basis(Z, +1), block_basis(Z, -1)
chi2 = np.array([np.trace(restrict(g, Vf)) for g in G])
chi2g = np.array([galois5(c.real) for c in chi2])          # Galois partner 2'
print("Teil A: Weil-Darstellung von Sp(2n,5), Howe-Paar 2I x O(n,5)\n")


def orthogonal_group(n):
    vecs = list(itertools.product(range(p), repeat=n))
    unit = [v for v in vecs if sum(x * x for x in v) % p == 1]
    out = []
    for cols in itertools.product(unit, repeat=n):
        M = np.array(cols).T % p
        if np.array_equal((M.T @ M) % p, np.eye(n, dtype=int)):
            out.append(M)
    return out


def perm_op(M, n):
    pts = list(itertools.product(range(p), repeat=n))
    index = {x: i for i, x in enumerate(pts)}
    P = np.zeros((len(pts), len(pts)))
    for i, x in enumerate(pts):
        y = tuple((M @ np.array(x)) % p)
        P[index[y], i] = 1                                     # (P f)(Mx) = f(x)
    return P


def kron_power(g, n):
    out = g
    for _ in range(n - 1):
        out = np.kron(out, g)
    return out


rowsA = []
for n in (1, 2, 3):
    D = p ** n
    rho = [kron_power(g, n) for g in G]
    On = orthogonal_group(n)
    ops = [perm_op(M, n) for M in On]
    Sn, Tn = kron_power(S, n), kron_power(T, n)
    commute = all(np.allclose(o @ Sn, Sn @ o) and np.allclose(o @ Tn, Tn @ o) for o in ops)
    Zn = kron_power(Z, n)
    dim_b = int(round(np.trace(np.eye(D) + Zn).real / 2))
    res = {"n": n, "D": D, "|O(n,5)|": len(On), "commute": commute, "bos": dim_b, "ferm": D - dim_b}
    for name, chi in (("2", chi2), ("2'", chi2g)):
        Pr = (2 / len(G)) * sum(np.conj(c) * r for c, r in zip(chi, rho))
        m = int(round(np.trace(Pr).real / 2))
        sig = np.array([np.trace(Pr @ o).real / 2 for o in ops])          # character of sigma
        norm = round(float(sum(s * s for s in sig) / len(On)), 6)
        # restriction of sigma to coordinate permutations S_n inside O(n,5)
        perms = [o for M, o in zip(On, ops) if all(sorted(col) == [0] * (n - 1) + [1] for col in M.T.tolist())
                 and all(sum(row) == 1 for row in M.tolist())]
        res[name] = (m, norm, len(perms))
    rowsA.append(res)
    print(f"n = {n}: C^{D}, |O({n},5)| = {len(On)}, O(n) vertauscht mit der Weil-Wirkung: {commute}; "
          f"bosonisch {res['bos']}, fermionisch {res['ferm']}")
    for name in ("2", "2'"):
        m, norm, _ = res[name]
        print(f"   Dublett {name}: Multiplizität {m}; Darstellung sigma von O({n},5) mit <χ,χ> = {norm}"
              f" ({'irreduzibel' if abs(norm - 1) < 1e-6 else 'reduzibel' if m else '-'})")

# ------------------------------------------------------------------ Teil B
print("\nTeil B: Wörterbuch über F5")
marks = [z for z in range(p) if (z ** 4 - 1) % p == 0]
print(f"   Nullstellen von z^4 - 1 in F5: {marks} (= ganz F5^x)")
polys = set()
for coeffs in itertools.product(range(p), repeat=5):
    polys.add(tuple(sum(c * x ** k for k, c in enumerate(coeffs)) % p for x in range(p)))
print(f"   Polynome vom Grad <= 4 als Funktionen auf F5: {len(polys)} verschiedene von {p ** 5} "
      f"(Bijektion: {len(polys) == p ** 5})")
R = np.zeros((p, p))
for x in range(p):
    R[(2 * x) % p, x] = 1
ev = np.linalg.eigvals(R)
print(f"   Dilatation x -> 2x auf C^F5: Spektrum {sorted(np.round(ev, 6), key=lambda c: (c.real, c.imag))}")
# Dirichlet characters mod 5: chi_k(2^j) = i^(k j)
dlog = {1: 0, 2: 1, 4: 2, 3: 3}
chars = {}
for k in range(4):
    vec = np.array([0] + [1j ** (k * dlog[x]) for x in range(1, p)])
    lam = (R @ vec)[2 % p] / vec[1] if True else None
    chars[k] = vec
    # eigenvalue of R on chi_k: (R chi)(2x) = chi(x) -> R chi = conj(chi(2)) chi
    eig = np.vdot(vec, R @ vec) / np.vdot(vec, vec)
    print(f"   Charakter chi_{k} (chi(2) = i^{k}): R-Eigenwert {np.round(eig, 6)}, "
          f"gerade: {np.allclose(vec[[0, 4, 3, 2, 1]], vec)}")
inv = np.array([1, 0, 0, 0, 0]), np.array([0, 1, 1, 1, 1])
print("   invarianter Teil: span(δ0, 1_{F5^x}), Dimension 2; nichttriviale Charaktere chi_1, chi_2, chi_3 "
      "mit R-Eigenwerten {-i, -1, +i} = TFPT Familiencharaktere {-1, i, -i}")
z = np.exp(2j * np.pi / p)
for k in (1, 2, 3):
    g = sum(chars[k][x] * z ** x for x in range(1, p))
    print(f"   Gauss-Summe g(chi_{k}) = {np.round(g, 6)}, |g|^2 = {abs(g) ** 2:.6f}")
print("   Auswahl 'Einheitengruppe = mu4' (F_p^x zyklisch der Ordnung 4):",
      {q: (q - 1 == 4) for q in (3, 5, 7, 11, 13, 17)})
print("   Zahlen (Familien p-2, Marken p-1, Plätze p, Fasern p+1):",
      {q: (q - 2, q - 1, q, q + 1) for q in (3, 5, 7)})

# ------------------------------------------------------------------ Teil C
print("\nTeil C: Hauptserie von SL(2,5) auf P^1(F5), getwistet mit Charakteren von F5^x")
P1 = [(1, x) for x in range(p)] + [(0, 1)]                      # [x:1] and infinity [1:0]


def normal(v):
    a, b = int(v[0]) % p, int(v[1]) % p
    if b != 0:
        inv_b = pow(b, p - 2, p)
        return ((a * inv_b) % p, 1)
    return (1, 0)


pts = [normal((x, 1)) for x in range(p)] + [(1, 0)]
idx = {q: i for i, q in enumerate(pts)}
SL = [np.array([[a, b], [c, d]]) for a, b, c, d in itertools.product(range(p), repeat=4) if (a * d - b * c) % p == 1]


def induced(g, k):
    """(g f)(x) = chi_k(cocycle) f(g^-1 x) on lines; cocycle = scaling of the representative."""
    ginv = np.array([[g[1, 1], -g[0, 1]], [-g[1, 0], g[0, 0]]]) % p
    M = np.zeros((6, 6), dtype=complex)
    for i, q in enumerate(pts):
        v = (ginv @ np.array(q)) % p
        r = normal(tuple(v))
        # v = s * r for a unit s
        s = next(t for t in range(1, p) if ((t * r[0] - v[0]) % p, (t * r[1] - v[1]) % p) == (0, 0))
        M[i, idx[r]] = 1j ** (k * dlog[s])
    return M


for k, name in ((0, "trivial"), (2, "Legendre (−1)"), (1, "quartisch (i)"), (3, "quartisch (−i)")):
    mats = [induced(g, k) for g in SL]
    ch = np.array([np.trace(m) for m in mats])
    norm = round(float((sum(abs(c) ** 2 for c in ch) / len(SL)).real), 6)
    minus = mats[next(i for i, g in enumerate(SL) if np.array_equal(g % p, np.array([[4, 0], [0, 4]])))]
    ferm = np.allclose(minus, -np.eye(6))
    hom = all(np.allclose(induced(g, k) @ induced(h, k), induced((g @ h) % p, k)) for g, h in
              zip(SL[:30], SL[30:60]))
    print(f"   chi = {name:15s}: Homomorphismus {hom}, <χ,χ> = {norm} "
          f"({'irreduzibel' if abs(norm - 1) < 1e-6 else f'{int(round(norm))} Bestandteile'}), "
          f"-I wirkt als {'-1 (fermionisch)' if ferm else '+1 (bosonisch)' if np.allclose(minus, np.eye(6)) else 'gemischt'}")

assert rowsA[0]["2"][0] == 1 and len(SL) == 120
