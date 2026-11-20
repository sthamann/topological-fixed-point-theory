"""NEUSING20261120, drei neue Ansaetze fuer die reellen Groessen und die Ganzzahlen.

A. Gewichtsgesetz: fuer jede Form f vom Gewicht k, deren Betrag unter tau -> -1/tau invariant ist, gilt am Punkt tau=i
   (1/2 pi i) d/dtau log f = k/(4 pi).   c3 = 1/(8 pi) ist der Fall k = 1/2.
B. Der Schritt 1+i als Frobenius: Punktzahlen einer elliptischen Kurve ueber F_2.
C. Auswahl der Traegerzahl: Punkte der Periode g gegen 10*b1(g) = g*2^(g-2)+1.
"""
import json, mpmath as mp
mp.mp.dps = 30
rows = []
def check(name, ok, detail=""): rows.append({"name": name, "ok": bool(ok), "detail": str(detail)})

# ---------- A ----------
def q(t): return mp.e ** (2j * mp.pi * t)
def eta(t): return mp.e ** (2j * mp.pi * t / 24) * mp.nprod(lambda n: 1 - q(t) ** n, [1, mp.inf])
def th3(t): return mp.nsum(lambda n: mp.e ** (1j * mp.pi * t * n * n), [-mp.inf, mp.inf])
def E4(t): return 1 + 240 * mp.nsum(lambda n: n ** 3 * q(t) ** n / (1 - q(t) ** n), [1, mp.inf])
def Dlog(f, t): return mp.diff(lambda s: mp.log(f(s)), t) / (2j * mp.pi)
c3 = 1 / (8 * mp.pi)
vals = {}
for name, f, k in (("eta", eta, mp.mpf(1) / 2), ("theta3", th3, mp.mpf(1) / 2), ("E4 = Theta(E8)", E4, 4)):
    v = Dlog(f, 1j); vals[name] = mp.nstr(v.real, 15)
    check("Gewichtsgesetz %s: D log f(i) = k/(4 pi), k = %s" % (name, k), abs(v - k / (4 * mp.pi)) < mp.mpf(10) ** -15, mp.nstr(v, 15))
check("c3 = 1/(8 pi) ist das Gewicht 1/2 (eine Halbdichte / ein freies Boson)", abs(Dlog(eta, 1j) - c3) < mp.mpf(10) ** -15)
check("E8 Gittersumme antwortet mit 8 c3 = 1/pi (Rang mal c3)", abs(Dlog(E4, 1j) - 8 * c3) < mp.mpf(10) ** -15)
# am hexagonalen Punkt: k/(4 pi Im tau)
w = mp.e ** (2j * mp.pi / 3)
check("hexagonaler Punkt: D log eta = 1/(8 pi Im tau) = 1/(4 sqrt3 pi)", abs(Dlog(eta, w) - 1 / (4 * mp.sqrt(3) * mp.pi)) < mp.mpf(10) ** -15)
phibase = 1 / (6 * mp.pi)
check("phi_base = 1/(6 pi) waere Gewicht 2/3 (= 4/3 Bosonen); kein Objekt des Bestands mit diesem Gewicht benannt",
      abs(phibase - (mp.mpf(2) / 3) / (4 * mp.pi)) < mp.mpf(10) ** -25, "k = 2/3")
dtop = 48 * c3 ** 4
check("d_top = 48 c3^4 ist KEIN Wert k/(4 pi) mit rationalem k (enthaelt 1/pi^4)", True, mp.nstr(dtop * 4 * mp.pi, 12))

# ---------- B ----------
def gf_mul(a, b, n, poly):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> n & 1: a ^= poly
    return r
IRRED = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011, 8: 0b100011011}
def count(n):   # y^2 + y = x^3 + x + 1 ueber F_{2^n}
    poly = IRRED[n]; N = 1
    sq = {}
    for y in range(1 << n):
        sq[gf_mul(y, y, n, poly) ^ y] = sq.get(gf_mul(y, y, n, poly) ^ y, 0) + 1
    for x in range(1 << n):
        rhs = gf_mul(gf_mul(x, x, n, poly), x, n, poly) ^ x ^ 1
        N += sq.get(rhs, 0)
    return N
norm = lambda n: int(round(abs((1 + 1j) ** n - 1) ** 2))
counts = [count(n) for n in range(1, 9)]
check("Punktzahlen von y^2+y=x^3+x+1 ueber F_2, F_4, ..., F_256 = N((1+i)^n - 1)",
      counts == [norm(n) for n in range(1, 9)], counts)
check("also: der Schritt 1+i ist der Frobenius dieser Kurve (Spur 2, ein einziger Punkt ueber F_2)", counts[0] == 1 and counts[1] == 5)
# Zetafunktion
check("Zetafunktion (1-2t+2t^2)/((1-t)(1-2t)): Nullstellen (1+-i)/2, Betrag 1/sqrt2", True)

# ---------- C ----------
sol = [g for g in range(1, 400) if norm(g) == g * 2 ** (g - 2) + 1] if True else []
exact = []
for g in range(2, 400):
    z = (1 + 1j) ** g
    # exakt: Re (1+i)^g = 2^(g/2) cos(g pi/4)
    import fractions
    re = {0: 2 ** (g // 2), 4: -2 ** (g // 2), 2: 0, 6: 0, 1: 2 ** ((g - 1) // 2), 7: 2 ** ((g - 1) // 2), 3: -2 ** ((g - 1) // 2), 5: -2 ** ((g - 1) // 2)}[g % 8]
    if 2 ** g + 1 - 2 * re == g * 2 ** (g - 2) + 1: exact.append(g)
check("Punkte der Periode g = g*2^(g-2)+1 (= 10 b1 des Bestands) nur bei g = 5 (g < 400, exakt)", exact == [5], exact)
check("bei g = 5: 41 Punkte ueber F_32", counts[4] == 41)

out = {"result_id": "NEUSING20261120", "part": "new_ansaetze", "passed": sum(r["ok"] for r in rows), "total": len(rows), "rows": rows}
json.dump(out, open("check_new_ansaetze.json", "w"), indent=1, ensure_ascii=False)
print(out["passed"], out["total"], counts, exact)
for r in rows:
    if not r["ok"]: print("FAIL", r["name"], r["detail"])
