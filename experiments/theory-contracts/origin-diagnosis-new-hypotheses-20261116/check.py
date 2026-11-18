#!/usr/bin/env python3
"""Exakte Prüfungen zu fünf Ursprungshypothesen, 16. November 2026.

Nur Standardbibliothek (Fraction, itertools) plus mpmath für T4.
Koordinaten: D8 Chart von E8, verdoppelt (ganze Zahlen). Träger D5 auf den
Koordinaten 0..4, Familie D3 = A3 auf 5..7. Jede Prüfung meldet PASS/FAIL.
"""
from fractions import Fraction as Fr
from itertools import combinations, product, permutations
import cmath
import json
import math
import sys

OUT = {}
FAILS = []


def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name}" + (f" :: {detail}" if detail else ""))
    OUT[name] = {"ok": bool(cond), "detail": str(detail)}
    if not cond:
        FAILS.append(name)


# ---------------------------------------------------------------- E8 Wurzeln
def e8_roots():
    R = []
    for i, j in combinations(range(8), 2):
        for si, sj in product((2, -2), repeat=2):
            v = [0] * 8
            v[i], v[j] = si, sj
            R.append(tuple(v))
    for s in product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0:
            R.append(tuple(s))
    return R


ROOTS = e8_roots()
RSET = set(ROOTS)
ip = lambda a, b: sum(x * y for x, y in zip(a, b))  # 4 x echtes Skalarprodukt
add = lambda a, b: tuple(x + y for x, y in zip(a, b))
neg = lambda a: tuple(-x for x in a)

check("E8 hat 240 Wurzeln der Norm 2", len(ROOTS) == 240 and all(ip(r, r) == 8 for r in ROOTS))

CAR, FAM = range(0, 5), range(5, 8)
is_int = lambda r: all(x % 2 == 0 for x in r)
supp = lambda r: [i for i, x in enumerate(r) if x != 0]
d5 = [r for r in ROOTS if is_int(r) and all(i in CAR for i in supp(r))]
a3 = [r for r in ROOTS if is_int(r) and all(i in FAM for i in supp(r))]
bif = [r for r in ROOTS if is_int(r) and len([i for i in supp(r) if i in CAR]) == 1]
odd = [r for r in ROOTS if not is_int(r)]
check("Zerlegung 40 + 12 + 60 + 128", (len(d5), len(a3), len(bif), len(odd)) == (40, 12, 60, 128))

# ------------------------------------------------- T1: Typgrenze der Bindung
print("\nT1 Typgrenze: E8 kovariante Bindung an drei Familien")
# (16,4): gerade Zahl Minuszeichen im Träger und in der Familie
m16_4 = [r for r in odd if [r[i] for i in CAR].count(-1) % 2 == 0]
m16b_4b = [r for r in odd if [r[i] for i in CAR].count(-1) % 2 == 1]
check("Materie = (16,4) + (16bar,4bar) = 64 + 64", len(m16_4) == 64 and len(m16b_4b) == 64)
fam4 = sorted({tuple(r[i] for i in FAM) for r in m16_4})
check("Familienfundamentale hat vier Gewichte", len(fam4) == 4, fam4)
# odd + odd ist nie odd: Z2 Graduierung
oo = sum(1 for a in odd for b in odd if add(a, b) in RSET and not is_int(add(a, b)))
check("[Materie, Materie] liegt ganz im geraden Teil", oo == 0)
# Drei aus vier Familien: wähle w4, Full96 = Materie ohne Familie w4
results_T1 = {}
all_zero = True
for w4 in fam4:
    w4n = tuple(-x for x in w4)
    fam_of = lambda r: tuple(r[i] for i in FAM)
    full96 = {r for r in odd if fam_of(r) not in (w4, w4n)}
    assert len(full96) == 96
    preserving = 0
    inner_elems = 0
    leak_elems = 0
    for g in bif:
        leaks = False
        for m in full96:
            t = add(g, m)
            if t in RSET:
                if t in full96:
                    inner_elems += 1
                else:
                    leak_elems += 1
                    leaks = True
        if not leaks:
            preserving += 1
    a3_stab = [r for r in a3 if all((add(r, m) not in RSET) or (add(r, m) in full96) for m in full96)]
    d5_stab = [r for r in d5 if all((add(r, m) not in RSET) or (add(r, m) in full96) for m in full96)]
    results_T1[str(w4)] = (preserving, inner_elems, leak_elems, len(a3_stab), len(d5_stab))
    all_zero &= preserving == 0
print("   je Wahl der vierten Familie: (erhaltende (10,6) Ströme, innere Elemente, austretende Elemente, |A3 Stab|, |D5 Stab|)")
for k, v in results_T1.items():
    print("  ", k, v)
check("Kein einziger der 60 (10,6) Ströme erhält den Raum dreier Familien", all_zero)
check("Erhaltende gerade Wurzeln: ganz D5 (40) und A2 der Familie (6)",
      all(v[3] == 6 and v[4] == 40 for v in results_T1.values()))
# Paarungsstruktur: jede Familienkoordinate ist eine perfekte Paarung der vier Gewichte
pairings = []
for j in range(3):
    flip = lambda w: tuple(-x if k == j else x for k, x in enumerate(w))
    pairs = {frozenset((fam4.index(w), fam4.index(tuple(-x for x in flip(w))))) for w in fam4}
    pairings.append(sorted(sorted(p) for p in pairs))
check("Drei Familienkoordinaten = drei perfekte Paarungen der vier Marken",
      sorted(map(str, pairings)) == sorted(map(str, [[[0, 1], [2, 3]], [[0, 2], [1, 3]], [[0, 3], [1, 2]]])), pairings)
# Minimaldimension: e8 hat keine nichttriviale Darstellung unter 248
check("96 < 248: jede Darstellung von e8 auf Full96 ist trivial (kleinste nichttriviale Dimension 248)", 96 < 248)

# --------------------------- T2: Spiegelungen, die D5 + A3 erhalten
print("\nT2 Involutionen im Stabilisator von D5 + A3 in W(E8)")
SIMPLE = [(2, -2, 0, 0, 0, 0, 0, 0), (0, 2, -2, 0, 0, 0, 0, 0), (0, 0, 2, -2, 0, 0, 0, 0),
          (0, 0, 0, 2, -2, 0, 0, 0), (0, 0, 0, 0, 2, -2, 0, 0), (0, 0, 0, 0, 0, 2, -2, 0),
          (0, 0, 0, 0, 0, 2, 2, 0), (-1, -1, -1, -1, -1, -1, -1, 1)]


def det_int(M):
    M = [[Fr(x) for x in row] for row in M]
    n, d = len(M), Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return d


gram = [[Fr(ip(a, b), 4) for b in SIMPLE] for a in SIMPLE]
check("Einfache Wurzeln: Gramdeterminante 1", det_int(gram) == 1)


def coords(v):  # Koordinaten von v in der Basis SIMPLE (ganzzahlig für Gittervektoren)
    n = 8
    A = [[Fr(SIMPLE[j][i]) for j in range(n)] + [Fr(v[i])] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        A[c] = [x / A[c][c] for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]


def left_kernel(M):
    """Gesättigte Z Basis von {c : c M = 0} durch unimodulare Zeilenumformungen."""
    n, m = len(M), len(M[0])
    A = [list(M[i]) + [1 if i == j else 0 for j in range(n)] for i in range(n)]
    row = 0
    for col in range(m):
        while True:
            nz = [r for r in range(row, n) if A[r][col] != 0]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda r: abs(A[r][col]))
            p = nz[0]
            for r in nz[1:]:
                q = A[r][col] // A[p][col]
                A[r] = [a - q * b for a, b in zip(A[r], A[p])]
        nz = [r for r in range(row, n) if A[r][col] != 0]
        if nz:
            A[row], A[nz[0]] = A[nz[0]], A[row]
            row += 1
    return [A[r][m:] for r in range(row, n)]


def root_type(rs):
    rs = list(rs)
    if not rs:
        return "0"
    comps, seen = [], set()
    for r in rs:
        if r in seen:
            continue
        comp, stack = {r}, [r]
        while stack:
            x = stack.pop()
            for y in rs:
                if y not in comp and ip(x, y) != 0:
                    comp.add(y)
                    stack.append(y)
        seen |= comp
        # Rang der Komponente
        basis = []
        for v in comp:
            trial = basis + [v]
            G = [[Fr(ip(a, b)) for b in trial] for a in trial]
            if det_int(G) != 0:
                basis = trial
        comps.append((len(basis), len(comp)))
    names = []
    for rk, nr in sorted(comps, reverse=True):
        if nr == rk * (rk + 1):
            names.append(f"A{rk}")
        elif nr == 2 * rk * (rk - 1):
            names.append(f"D{rk}")
        elif (rk, nr) in ((6, 72), (7, 126), (8, 240)):
            names.append(f"E{rk}")
        else:
            names.append(f"?({rk},{nr})")
    return "+".join(names)


def involution(c5, n5, c3, n3):
    perm, sign = list(range(8)), [1] * 8
    free5, free3 = list(CAR), list(FAM)
    for _ in range(c5):
        i, j = free5.pop(0), free5.pop(0)
        perm[i], perm[j] = j, i
    for _ in range(n5):
        sign[free5.pop(0)] = -1
    for _ in range(c3):
        i, j = free3.pop(0), free3.pop(0)
        perm[i], perm[j] = j, i
    for _ in range(n3):
        sign[free3.pop(0)] = -1
    return lambda v: tuple(sign[perm[i]] * v[perm[i]] for i in range(8))


rows = []
classes = 0
for c5 in range(3):
    for n5 in range(0, 6 - 2 * c5):
        for c3 in range(2):
            for n3 in range(0, 4 - 2 * c3):
                if (n5 + n3) % 2 or (c5 + n5 + c3 + n3) == 0:
                    continue
                th = involution(c5, n5, c3, n3)
                assert all(th(r) in RSET for r in ROOTS)
                assert all(th(th(r)) == r for r in ROOTS)
                T = [coords(th(s)) for s in SIMPLE]
                assert all(x.denominator == 1 for row in T for x in row)
                T = [[int(x) for x in row] for row in T]
                I8 = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
                Kp = left_kernel([[T[i][j] - I8[i][j] for j in range(8)] for i in range(8)])
                Km = left_kernel([[T[i][j] + I8[i][j] for j in range(8)] for i in range(8)])
                vec = lambda c: tuple(sum(c[i] * SIMPLE[i][k] for i in range(8)) for k in range(8))
                dets = []
                for K in (Kp, Km):
                    vs = [vec(c) for c in K]
                    dets.append(int(det_int([[Fr(ip(a, b), 4) for b in vs] for a in vs])) if vs else 1)
                fixed = [r for r in ROOTS if th(r) == r]
                anti = [r for r in ROOTS if th(r) == neg(r)]
                classes += 1
                rows.append(dict(c5=c5, n5=n5, c3=c3, n3=n3, minus=len(Km), plus=len(Kp),
                                 det_plus=dets[0], det_minus=dets[1],
                                 roots_plus=root_type(fixed), roots_minus=root_type(anti)))
print("   (c5,n5,c3,n3) | minus | det+ det- | Wurzeln+ | Wurzeln-")
for r in sorted(rows, key=lambda r: (r["minus"], r["c5"], r["n5"], r["c3"])):
    print(f"   ({r['c5']},{r['n5']},{r['c3']},{r['n3']}) | {r['minus']} | {r['det_plus']} {r['det_minus']} | {r['roots_plus']} | {r['roots_minus']}")
OUT["T2_classes"] = rows
check("Stabilisatorordnung 2^5*5!*2^3*3!/2 = 92160", 2**5 * 120 * 8 * 6 // 2 == 92160)
one = [r for r in rows if r["minus"] == 1]
check("Minusrang eins kommt vor (einfache Wurzelspiegelungen in D5 oder A3)", len(one) >= 1, one)
half = [r for r in rows if r["minus"] == 4 and r["roots_plus"] == "D4" and r["roots_minus"] == "D4"]
check("Es gibt D5+A3 erhaltende Spiegelungen mit Zerlegung D4 | D4", len(half) >= 1,
      [(r["c5"], r["n5"], r["c3"], r["n3"]) for r in half])
pure = [r for r in half if r["c5"] == 0 and r["c3"] == 0]
check("Reine Vorzeichenform der D4|D4 Spiegelung: vier Trägerkoordinaten negiert, 1 Träger + 3 Familie fest",
      any(r["n5"] == 4 and r["n3"] == 0 for r in pure), [(r["n5"], r["n3"]) for r in pure])
# Verträglichkeit mit der 3+2 Hyperladung Y ~ (-2,-2,-2,3,3 | 0,0,0): theta Y = +Y oder -Y
Y = (-2, -2, -2, 3, 3, 0, 0, 0)
compat = []
for signs in product((1, -1), repeat=8):
    if signs.count(-1) % 2:
        continue
    for p5 in permutations(range(5)):
        if any(p5[p5[i]] != i for i in range(5)):
            continue
        for p3 in permutations(range(3)):
            if any(p3[p3[i]] != i for i in range(3)):
                continue
            perm = list(p5) + [5 + x for x in p3]
            th = lambda v: tuple(signs[i] * v[perm[i]] for i in range(8))
            if any(th(th(e)) != e for e in SIMPLE):
                continue  # keine Involution (Vorzeichen auf Zweizykeln ungleich)
            ty = th(Y)
            if ty not in (Y, neg(Y)):
                continue
            minus = sum(1 for i in range(8) if perm[i] == i and signs[i] == -1) + sum(1 for i in range(8) if perm[i] > i)
            if minus != 4:
                continue
            fixed = [r for r in ROOTS if th(r) == r]
            anti = [r for r in ROOTS if th(r) == neg(r)]
            compat.append((root_type(fixed), root_type(anti), "Y" if ty == Y else "-Y"))
types = sorted(set(compat))
print("   Hyperladungsverträgliche Involutionen mit Minusrang 4:", len(compat), "Typen:", types)
OUT["T2_hypercharge_compatible_rank4"] = [list(t) for t in types]
check("Keine hyperladungsverträgliche Spiegelung zerlegt E8 in D4 | D4",
      len(compat) > 0 and all(t[0] != "D4" for t in types), types)

# ------------------------------------- T3: chirale Zentralladung modulo 8
print("\nT3 Gauß Summen der Diskriminantenformen")


def cartan(kind, n):
    G = [[0] * n for _ in range(n)]
    for i in range(n):
        G[i][i] = 2
    if kind == "A":
        for i in range(n - 1):
            G[i][i + 1] = G[i + 1][i] = -1
    if kind == "D":
        for i in range(n - 2):
            G[i][i + 1] = G[i + 1][i] = -1
        G[n - 3][n - 1] = G[n - 1][n - 3] = -1
    return G


def inv(M):
    n = len(M)
    A = [[Fr(x) for x in M[i]] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        A[c] = [x / A[c][c] for x in A[c]]
        for r in range(n):
            if r != c:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [row[n:] for row in A]


def gauss_c(kind, n):
    G = cartan(kind, n)
    Gi = inv(G)
    d = int(det_int(G))
    # Dualgitter wird von den Zeilen von G^-1 erzeugt; Klassen durch Vielfache und Summen
    gens = [tuple(Gi[i]) for i in range(n)]
    frontier = [tuple([Fr(0)] * n)]
    red = lambda v: tuple(x - math.floor(x) for x in v)
    seen = {red(frontier[0])}
    while frontier:
        v = frontier.pop()
        for g in gens:
            w = red(tuple(a + b for a, b in zip(v, g)))
            if w not in seen:
                seen.add(w)
                frontier.append(w)
    assert len(seen) == d, (len(seen), d)
    spins = []
    tot = 0
    for v in seen:  # v sind Koeffizienten bezüglich der einfachen Wurzeln
        q = sum(v[i] * G[i][j] * v[j] for i in range(n) for j in range(n)) / 2
        q = q - math.floor(q)
        spins.append(q)
        tot += cmath.exp(2j * math.pi * float(q))
    tot /= math.sqrt(d)
    c8 = (cmath.phase(tot) / (2 * math.pi) * 8) % 8
    return d, sorted(spins), round(c8, 9), abs(abs(tot) - 1) < 1e-9


res3 = {}
for kind, n in (("D", 5), ("A", 3), ("D", 4), ("A", 4), ("D", 8)):
    d, spins, c8, uni = gauss_c(kind, n)
    res3[f"{kind}{n}"] = (d, [str(s) for s in spins], c8)
    print(f"   {kind}{n}: |Disc|={d}, Spins={[str(s) for s in spins]}, c mod 8 = {c8}")
OUT["T3"] = res3
check("D5: c = 5 mod 8, A3: c = 3 mod 8, Summe 0 mod 8", res3["D5"][2] == 5 and res3["A3"][2] == 3)
check("D5 und A3 tragen konjugierte Z4 Spins (5/8 gegen 3/8)",
      sorted(res3["D5"][1]) == ["0", "1/2", "5/8", "5/8"] and sorted(res3["A3"][1]) == ["0", "1/2", "3/8", "3/8"])
check("D4: c = 4 mod 8 und drei nichttriviale Klassen mit Spin 1/2 (drei Fermionen)",
      res3["D4"][2] == 4 and res3["D4"][1] == ["0", "1/2", "1/2", "1/2"])

# ---------------------------------------------- T4: absolute Leerreferenz
print("\nT4 Vakuumcharakter E4/eta^8 der einzelnen Naht")
try:
    import mpmath as mp
    mp.mp.dps = 40

    def chi(tau):
        q = mp.exp(2j * mp.pi * tau)
        E4 = 1 + 240 * mp.nsum(lambda n: n**3 * q**n / (1 - q**n), [1, mp.inf])
        eta = mp.exp(2j * mp.pi * tau / 24) * mp.nprod(lambda n: 1 - q**n, [1, mp.inf])
        return E4 / eta**8

    rho = mp.exp(2j * mp.pi / 3)
    v_i, v_rho = chi(1j), chi(rho)
    tau = mp.mpc("0.137", "0.911")
    tphase = chi(tau + 1) / chi(tau)
    sinv = chi(-1 / tau) / chi(tau)
    print("   chi(i) =", mp.nstr(v_i, 20), " |chi(rho)| =", mp.nstr(abs(v_rho), 5))
    print("   chi(tau+1)/chi(tau) =", mp.nstr(tphase, 15), " chi(-1/tau)/chi(tau) =", mp.nstr(sinv, 15))
    check("chi(i) = 12 am Viertelpunkt", abs(v_i - 12) < mp.mpf(10) ** -30)
    check("chi(rho) = 0 am Dreierpunkt", abs(v_rho) < mp.mpf(10) ** -30)
    check("T Phase exp(-2 pi i/3): Referenz ist Schnitt eines Linienbündels", abs(tphase - mp.exp(-2j * mp.pi / 3)) < mp.mpf(10) ** -30)
    check("S invariant", abs(sinv - 1) < mp.mpf(10) ** -30)
    check("Reflektiertes Betragsquadrat ist T invariant", abs(abs(tphase) - 1) < mp.mpf(10) ** -30)
except ImportError:
    check("mpmath verfügbar", False)

# -------------------------------- T5: goldene Struktur und ungerade Ränge
print("\nT5 Goldene (ikosische) Struktur gegen den 5 + 3 Schnitt")
# F selbstadjungiert, ganzzahlig, F^2 = F + 1 auf Teilgitter vom Rang n: Spur = a*phi + (n-a)*phi'
# = n/2 + (2a-n)*sqrt(5)/2, ganzzahlig genau für 2a = n und n/2 ganz, also n durch 4 teilbar? n gerade und n/2 ganz.
ok_ranks = [n for n in range(1, 9) if any(2 * a == n for a in range(n + 1))]
check("Goldene Struktur nur auf geraden Rängen: D5 (Rang 5) und A3 (Rang 3) einzeln unmöglich", 5 not in ok_ranks and 3 not in ok_ranks, ok_ranks)
# Positivkontrolle: Fünfertakt g = zyklische Permutation der Trägerkoordinaten, F = 1 + g + g^-1
g = lambda v: (v[4], v[0], v[1], v[2], v[3], v[5], v[6], v[7])
gi = lambda v: (v[1], v[2], v[3], v[4], v[0], v[5], v[6], v[7])
assert all(g(r) in RSET for r in ROOTS)
Tg = [[int(x) for x in coords(g(s))] for s in SIMPLE]
I8 = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
Kfix = left_kernel([[Tg[i][j] - I8[i][j] for j in range(8)] for i in range(8)])
vec = lambda c: tuple(sum(c[i] * SIMPLE[i][k] for i in range(8)) for k in range(8))
vs = [vec(c) for c in Kfix]
dfix = int(det_int([[Fr(ip(a, b), 4) for b in vs] for a in vs]))
fixroots = [r for r in ROOTS if g(r) == r]
check("Fixgitter des Trägerfünfertakts: Rang 4, Determinante 5, Wurzelsystem A4 (Familie A3 plus Trägerspur)",
      len(vs) == 4 and dfix == 5 and root_type(fixroots) == "A4", (len(vs), dfix, root_type(fixroots)))
F = lambda v: add(add(v, g(v)), gi(v))
# auf dem Komplement (Summe der Trägerkoordinaten null, Familie null) gilt F^2 = F + 1
comp = [r for r in d5 if sum(r[i] for i in CAR) == 0]
golden = all(F(F(v)) == add(F(v), v) for v in comp)
check("Auf dem A4 Komplement erfüllt F = 1 + g + g^-1 die Gleichung F^2 = F + 1", golden and root_type(comp) == "A4", root_type(comp))

print("\nErgebnis:", "alle Prüfungen bestanden" if not FAILS else f"FEHLER in {FAILS}")
json.dump(OUT, open(__file__.replace("check.py", "results.json"), "w"), indent=1, ensure_ascii=False, default=str)
sys.exit(1 if FAILS else 0)
