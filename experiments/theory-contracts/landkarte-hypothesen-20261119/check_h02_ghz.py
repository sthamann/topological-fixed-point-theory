"""H2: symmetrische Stabilisatorzustände, Majorana Marken, und der tatsächliche Fünfer Fix_P(Sym^4 C^4)."""
import itertools, math, numpy as np
from collections import Counter
from common import Checks
c = Checks("h02_ghz")
H = np.array([[1, 1], [1, -1]]) / math.sqrt(2); S = np.diag([1, 1j])
def stabilizer_states(n):
    N = 2 ** n
    def one(g, q):
        M = np.eye(1)
        for i in range(n): M = np.kron(M, g if i == q else np.eye(2))
        return M
    def cz(a, b): return np.diag([-1.0 if (x >> (n - 1 - a)) & 1 and (x >> (n - 1 - b)) & 1 else 1.0 for x in range(N)])
    gens = [one(H, q) for q in range(n)] + [one(S, q) for q in range(n)] + [cz(a, b) for a, b in itertools.combinations(range(n), 2)]
    def key(v):
        i = np.flatnonzero(np.abs(v) > 1e-9)[0]; w = v / v[i]; w = w / np.linalg.norm(w)
        return tuple(np.round(w.real, 6)) + tuple(np.round(w.imag, 6))
    v0 = np.zeros(N, complex); v0[0] = 1; seen = {key(v0): v0}; front = [v0]
    while front:
        nxt = []
        for v in front:
            for g in gens:
                w = g @ v; k = key(w)
                if k not in seen: seen[k] = w; nxt.append(w)
        front = nxt
    return list(seen.values())
def symmetric(v, n):
    N = 2 ** n
    for p in [tuple([1, 0] + list(range(2, n))), tuple(list(range(1, n)) + [0])]:
        w = np.zeros(N, complex)
        for x in range(N):
            bits = [(x >> (n - 1 - i)) & 1 for i in range(n)]; y = 0
            for i in range(n): y = (y << 1) | bits[p[i]]
            w[y] = v[x]
        if not np.allclose(w, v, atol=1e-7): return False
    return True
def stars(v, n):                       # Majorana Punkte auf der Kugel
    N = 2 ** n; v = v / np.linalg.norm(v)
    a = [sum(v[x] for x in range(N) if bin(x).count("1") == k) / math.comb(n, k) for k in range(n + 1)]
    co = [(-1) ** k * math.comb(n, k) * a[k] for k in range(n + 1)]       # sum_k co_k z^(n-k)
    lead = next(i for i, x in enumerate(co) if abs(x) > 1e-9)
    rts = list(np.roots(co[lead:])) if lead < n else []
    pts = [np.array([0, 0, 1.0])] * lead
    for z in rts: pts.append(np.array([2 * z.real, 2 * z.imag, abs(z) ** 2 - 1]) / (abs(z) ** 2 + 1))
    return pts
def polygon(pts):                       # regelmäßiges n Eck auf einem Großkreis?
    P = np.array(pts); n = len(P)
    if np.linalg.matrix_rank(P, tol=1e-6) != 2: return False
    d = sorted(round(float(np.linalg.norm(P[i] - P[j])), 5) for i in range(n) for j in range(i + 1, n))
    Q = np.array([[math.cos(2 * math.pi * k / n), math.sin(2 * math.pi * k / n), 0] for k in range(n)])
    return d == sorted(round(float(np.linalg.norm(Q[i] - Q[j])), 5) for i in range(n) for j in range(i + 1, n))
res = {}
for n in (3, 4):
    st = stabilizer_states(n); sym = [v for v in st if symmetric(v, n)]
    kinds = Counter()
    for v in sym:
        P = stars(v, n); distinct = min(np.linalg.norm(P[i] - P[j]) for i in range(n) for j in range(i + 1, n)) > 0.05   # mehrfache Nullstellen sind numerisch nur auf 1e-3 genau
        kinds["zusammenfallend" if not distinct else ("regelmäßiges Vieleck auf Großkreis" if polygon(P) else "andere Lage")] += 1
    res[n] = (len(st), len(sym), dict(kinds))
c.ok("3 Qubits: 1080 Stabilisatorzustände; alle symmetrischen mit drei verschiedenen Marken sind gleichseitige Dreiecke auf einem Großkreis",
     res[3][0] == 1080 and "andere Lage" not in res[3][2], res[3])
c.ok("4 Qubits: 36720 Stabilisatorzustände, 18 symmetrisch, 12 mit vier verschiedenen Marken, alle zwölf sind Quadrate auf einem Großkreis",
     res[4][0] == 36720 and res[4][1] == 18 and res[4][2].get("regelmäßiges Vieleck auf Großkreis") == 12 and "andere Lage" not in res[4][2], res[4])
c.ok("GRENZE: das Prinzip wählt für n Qubits das regelmäßige n Eck, also die Lage und nicht die Vier", res[3][2].get("regelmäßiges Vieleck auf Großkreis", 0) > 0)
# der tatsächliche Fünfer des Originals: Fixraum der Zwei Qubit Pauligruppe in Sym^4(C^4)
X = np.array([[0, 1], [1, 0]]); Z = np.diag([1, -1]); I2 = np.eye(2)
paulis = [np.kron(a, b) for a in (I2, X, Z, X @ Z) for b in (I2, X, Z, X @ Z)]
def h4(M):
    e = np.linalg.eigvals(M); p = [sum(e ** k) for k in range(1, 5)]
    return (p[0] ** 4 + 6 * p[0] ** 2 * p[1] + 3 * p[1] ** 2 + 8 * p[0] * p[2] + 6 * p[3]) / 24
full = sum(np.trace(ph * P) ** 4 for P in paulis for ph in (1, 1j, -1, -1j)) / 64
sym4 = sum(h4(ph * P) for P in paulis for ph in (1, 1j, -1, -1j)) / 64
c.ok("Fixraum der Pauligruppe (Ordnung 64): 16 im vollen vierten Tensorgrad, 5 in Sym^4(C^4) (wie im Atlas)", abs(full - 16) < 1e-9 and abs(sym4 - 5) < 1e-9, (full.real, sym4.real))
c.ok("TYPGRENZE: der Fünfer des Originals liegt in Sym^4(C^4) (35 dim), die Majorana Lesart braucht Sym^4(C^2) (5 dim); das sind verschiedene Räume",
     math.comb(4 + 3, 4) == 35 and math.comb(4 + 1, 4) == 5)
# gemeinsame Darstellung nur über die Ikosaedergruppe: Charakter von Sym^4(C^2) gegen Standarddarstellung von S6 auf PSL(2,5)
chi5 = lambda th: 5 if th == 0 else math.sin(5 * th / 2) / math.sin(th / 2)
spin2 = [round(chi5(t), 9) for t in (0, math.pi, 2 * math.pi / 3, 2 * math.pi / 5, 4 * math.pi / 5)]
perm6 = [6 - 1, 2 - 1, 0 - 1, 1 - 1, 1 - 1]          # Fixpunkte von PSL(2,5) auf den 6 Punkten der projektiven Geraden, minus 1
c.ok("als Darstellung der Ikosaedergruppe A5 stimmen beide Fünfer überein (Charaktere 5, 1, -1, 0, 0)", spin2 == [5, 1, -1, 0, 0] and perm6 == [5, 1, -1, 0, 0], (spin2, perm6))
# Zerlegung von Spin 2 unter der Symmetrie des Quadrats (D4, Achse durch die Pole): m = 0 | +-1 | +-2
d4 = [(1, 0.0), (2, math.pi / 2), (1, math.pi), (2, math.pi), (2, math.pi)]          # Klassen von D4: e, 2 C4, C2, 2 C2', 2 C2''
chiE = [2, 0, -2, 0, 0]
multE = sum(k * chi5(t) * e for (k, t), e in zip(d4, chiE)) / 8
mult1 = sum(k * chi5(t) ** 2 for k, t in d4) / 8 - multE ** 2                      # Summe m_i^2 minus Dublettanteil
c.ok("unter der Quadratsymmetrie D4 zerfällt Spin 2 in genau ein Dublett und drei Einzelne: 2 + (1+1+1), kein Triplett",
     abs(multE - 1) < 1e-9 and abs(mult1 - 3) < 1e-9, (round(multE, 9), round(mult1, 9)))
c.done({"zaehlung": res, "urteil": "In der eigenen Klasse exakt: symmetrische Stabilisatorzustände mit n verschiedenen Marken sind regelmäßige n Ecke. Bindung an das Original nur über eine nicht gezeigte Wahl einer Ikosaedergruppe in S6."})
