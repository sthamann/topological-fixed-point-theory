"""H10: Frobenius Spuren der Kurve mit Multiplikation durch i als Ereignisprotokoll. Festlegung vor der Auszählung, dann Vergleich mit Zufallskurven."""
import numpy as np, sympy as sp
from common import Checks
c = Checks("h10_frobenius")
# --- Festlegung (steht vor jeder Auszählung fest)
T = {2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 30, 41, 48}          # ganze Zahlen, die im Kern der Theorie getragen werden
P_TH = [5, 13, 41, 137]                                          # Primzahlen der Theorie
P_ALL = [p for p in sp.primerange(5, 1000)]                      # zweiter, breiter Satz
NORM = "Betrag der Spur |a_p|, a_p = p + 1 - Punktezahl"
def chi_table(p):
    t = -np.ones(p, dtype=np.int64); t[(np.arange(1, p, dtype=np.int64) ** 2) % p] = 1; t[0] = 0; return t
def traces(A, B, p):
    """a_p für viele Kurven y^2 = x^3 + A x + B zugleich (A, B Vektoren)."""
    x = np.arange(p, dtype=np.int64)[None, :]; f = (x ** 3 % p + (A[:, None] % p) * x + B[:, None]) % p
    return -chi_table(p)[f].sum(axis=1)
one = lambda A, B, p: int(traces(np.array([A]), np.array([B]), p)[0])
# --- 1. die Spur ist durch p = a^2 + b^2 erzwungen
def two_squares(p):
    for a in range(1, int(p ** .5) + 1, 2):
        b2 = p - a * a; b = int(round(b2 ** .5))
        if b * b == b2: return a, b
ok1 = all((one(-1, 0, p) == 0) if p % 4 == 3 else (abs(one(-1, 0, p)) == 2 * two_squares(p)[0]) for p in P_ALL)
c.ok("für y^2 = x^3 - x gilt bei allen 166 Primzahlen von 5 bis 1000: a_p = 0 für p = 3 mod 4, sonst |a_p| = 2a mit p = a^2 + b^2, a ungerade", ok1 and len(P_ALL) == 166, len(P_ALL))
th = {p: abs(one(-1, 0, p)) for p in P_TH}
c.ok("an den Primzahlen der Theorie: |a_5| = 2, |a_13| = 6, |a_41| = 10, |a_137| = 22", th == {5: 2, 13: 6, 41: 10, 137: 22}, th)
hits_th = sum(v in T for v in th.values())
c.ok("drei der vier Werte liegen in der festgelegten Menge (2, 6, 10), 22 nicht", hits_th == 3)
tw = {D: [abs(one(-D, 0, p)) for p in P_TH] for D in range(1, 41) if all(D % p for p in P_TH)}
vals = [sorted({v[i] for v in tw.values()}) for i in range(4)]
c.ok("alle Drehungen y^2 = x^3 - Dx derselben Kurve geben an diesen Primzahlen nur 2a oder 2b: {2,4}, {4,6}, {8,10}, {8,22}", vals == [[2, 4], [4, 6], [8, 10], [8, 22]], vals)
c.ok("FOLGE: ein Treffer sagt nur, dass 2a oder 2b in der Menge liegt. Das ist eine Aussage über die Zerlegung der Primzahl in zwei Quadrate, nicht über die Theorie", True)
# --- 2. Vergleich mit 1000 Zufallskurven, gleiche Primzahlen, gleiche Normierung
rng = np.random.default_rng(20261119); N = 1000
A = rng.integers(-500, 501, N); B = rng.integers(-500, 501, N)
good = np.array([(4 * int(a) ** 3 + 27 * int(b) ** 2) != 0 for a, b in zip(A, B)]); A, B = A[good], B[good]
Tarr = np.array(sorted(T))
def hits(A, B, primes):
    h = np.zeros(len(A), dtype=int)
    for p in primes: h += np.isin(np.abs(traces(A, B, p)), Tarr)
    return h
h_th = hits(A, B, P_TH); frac_th = float((h_th >= hits_th).mean())
c.ok("enger Satz: mindestens drei Treffer aus vier erreicht ein merklicher Teil der Zufallskurven (kein seltenes Ereignis)", frac_th > 0.05, f"Anteil {frac_th:.3f}, Mittel {h_th.mean():.2f}")
h_all = hits(A, B, P_ALL); cm_all = int(hits(np.array([-1]), np.array([0]), P_ALL)[0]); frac_all = float((h_all >= cm_all).mean())
c.ok("breiter Satz (166 Primzahlen): die Kurve mit i hat WENIGER Treffer als fast alle Zufallskurven, weil die Hälfte ihrer Spuren null ist", frac_all > 0.95, f"Kurve {cm_all}, Zufall Mittel {h_all.mean():.1f}, Anteil mit mindestens so vielen {frac_all:.3f}")
P1 = [p for p in P_ALL if p % 4 == 1]; cm1 = int(hits(np.array([-1]), np.array([0]), P1)[0]); h1 = hits(A, B, P1); frac1 = float((h1 >= cm1).mean())
c.ok("auch nur auf den Primzahlen p = 1 mod 4 liegt die Kurve im Bereich der Zufallskurven", 0.05 < frac1, f"Kurve {cm1} von {len(P1)}, Zufall Mittel {h1.mean():.1f}, Anteil {frac1:.3f}")
c.ok("ABBRUCHREGEL DER HYPOTHESE GREIFT: die Trefferzahl liegt im Bereich der Zufallskurven. Die Idee ist in dieser Form erledigt", frac_th > 0.05 and frac1 > 0.05)
# --- 3. was exakt übrig bleibt: die Vier als Galoisgruppe ordnet jeder Primzahl ein Element zu
cls = {p: p % 5 for p in P_TH}
c.ok("Frobenius in der Uhrengruppe (Z/5)^* = Z/4: 41 wirkt als 1, 13 als 3, 137 als 2, 5 ist verzweigt", cls == {5: 0, 13: 3, 41: 1, 137: 2}, cls)
cnt = {r: sum(1 for p in sp.primerange(7, 100000) if p % 5 == r) for r in (1, 2, 3, 4)}; tot = sum(cnt.values())
c.ok("über alle Primzahlen bis 100000 sind die vier Klassen gleich häufig (je ein Viertel, Dirichlet)", all(abs(v / tot - .25) < .005 for v in cnt.values()), cnt)
c.ok("GRENZE: diese Zuordnung ist ein Etikett je Primzahl. Sie hat keine Reihenfolge außer der Größe der Primzahl, keine Rate und keine Zusammensetzung von Ereignissen. Dynamik folgt daraus nicht", True)
c.done({"festlegung": {"T": sorted(T), "primzahlen": P_TH, "normierung": NORM, "zufallskurven": int(len(A)), "seed": 20261119},
        "treffer_theorie_primzahlen": hits_th, "anteil_zufall_mindestens_so_gut": frac_th, "breit": {"kurve": cm_all, "zufall_mittel": float(h_all.mean()), "anteil": frac_all},
        "p_gleich_1_mod_4": {"kurve": cm1, "von": len(P1), "zufall_mittel": float(h1.mean()), "anteil": frac1},
        "urteil": "Fällt in dieser Form. Die Spuren sind durch p = a^2 + b^2 erzwungen und tragen keine Auskunft über die Theorie; die Trefferzahl liegt im Bereich der Zufallskurven. Exakt bleibt nur das Etikett p mod 5 in der Vier."})
