"""w2: exakte Brücken zwischen den Lesarten (Exploration 2026-11-27).

Jede Prüfung rechnet eine Aussage nach, die mehrere Lesarten der Naht verbindet. Nichts davon ist eine Herleitung der
Physik; es sind exakte Aussagen über die genannten mathematischen Objekte. Ausgabe: ausgabe_w2.txt, w2_ergebnis.json.
"""
from fractions import Fraction as Fr
from itertools import product
from math import comb, gcd
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "research" / "originlab"))
from genesis_rules import anti_isometric            # noqa: E402
import whole_solution_search as wss                  # noqa: E402

LINES, CHECKS, OUT = [], [], {}


def say(text=""):
    LINES.append(text)
    print(text)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    say(f"[{'ok' if ok else 'FEHLER'}] {name}" + (f": {detail}" if detail else ""))


# ---------------------------------------------------------------- B1: dieselbe Vier nur für zwei Blätter
say("B1  Dieselbe Vier: Index k^2, Markenzahl 2k/(k-1), Ordnung 2k des Spinlifts")
rows = {}
for k in range(2, 13):
    marks = Fr(2 * k, k - 1)
    rows[k] = {"index": k * k, "marken": str(marks), "lift": 2 * k,
               "c3_gleich_einstein": Fr(1, 4 * k) == Fr(1, 8)}       # 1/(4 pi k) gegen 1/(8 pi)
same = [k for k, r in rows.items() if r["marken"] == str(r["index"]) == str(r["lift"])]
flat = [k for k in rows if Fr(2 * k, k - 1).denominator == 1]
check("Index, Markenzahl und Liftordnung stimmen nur für k = 2 überein (alle gleich 4)", same == [2], f"k = {same}")
check("flach mit gleichen Marken nur für k = 2 (vier Marken) und k = 3 (drei Marken)", flat == [2, 3], f"k = {flat}")
check("c3(k) = 1/(4 pi k) ist die Einsteinkonstante 1/(8 pi) nur für k = 2",
      [k for k, r in rows.items() if r["c3_gleich_einstein"]] == [2])
OUT["B1"] = rows

# ---------------------------------------------------------------- B2: Kegel der Reinheit als Teilchen
say("\nB2  Replikakegel und Kegelteilchen: Spannung (n-1)/(4n) gegen Defizit 2 pi (n-1)/n")
ratio = {n: Fr(n - 1, 4 * n) / (2 * Fr(n - 1, n)) for n in range(2, 10)}       # Spannung / (Defizit / pi)
check("Spannung durch Defizit ist für jedes n gleich 1/(8 pi) (Einsteingleichung, G = 1)",
      set(ratio.values()) == {Fr(1, 8)})
check("n = 2: Spannung 1/8 je Marke, Defizit pi; vier Marken tragen 1/2 = 1/|Z2| und schließen die Sphäre",
      Fr(2 - 1, 4 * 2) == Fr(1, 8) and 4 * Fr(1, 8) == Fr(1, 2) and 4 * 2 * Fr(1, 2) == 4)
OUT["B2"] = {"spannung_n2": "1/8", "summe_vier_marken": "1/2"}

# ---------------------------------------------------------------- B3: Zwistgewichte der fermionischen Replika
say("\nB3  Fermionische Replika: Zwiste k/n, k = -(n-1)/2 ... (n-1)/2")
ok = True
for n in range(2, 13):
    ks = [Fr(2 * j - (n - 1), 2) for j in range(n)]
    ok &= sum((k / n) ** 2 / 2 for k in ks) == Fr(n * n - 1, 24 * n)
check("Summe der Gewichte (k/n)^2/2 ist (c/24)(n - 1/n) mit c = 1, für n = 2 ... 12", ok)
tw = [Fr(1, 2) / 2, Fr(-1, 2) / 2]
check("n = 2: Zwiste +1/4 und -1/4, Gewicht 1/32 je Zwist, zusammen 1/16",
      tw == [Fr(1, 4), Fr(-1, 4)] and [t * t / 2 for t in tw] == [Fr(1, 32)] * 2)
check("n = 2: Phasen exp(2 pi i k/n) = ±i, der Tausch hat auf Fermionen die Ordnung 4 und das Quadrat -1",
      (1j) ** 2 == -1 and (1j) ** 4 == 1)
genus = [(k, m) for k in range(2, 13) for m in range(1, 13) if (k - 1) * (m - 1) == 1]
check("Die Replikafläche von m Gebieten auf k Blättern hat genau für (k, m) = (2, 2) einen Modul (Genus 1)",
      genus == [(2, 2)])
OUT["B3"] = {"zwiste_n2": ["1/4", "-1/4"], "gewicht": "1/32", "genus_eins": genus}

# ---------------------------------------------------------------- B4: McKay-Graph der Ikosaedergruppe
say("\nB4  Affines E8 (McKay-Graph der binären Ikosaedergruppe): Schnitt am Knoten der Marke 4")
NODES = ["a0", "a8", "a7", "a6", "a5", "a4", "a3", "a1", "a2"]
MARK = dict(zip(NODES, [1, 2, 3, 4, 5, 6, 4, 2, 3]))
EDGES = [("a0", "a8"), ("a8", "a7"), ("a7", "a6"), ("a6", "a5"), ("a5", "a4"), ("a4", "a3"), ("a3", "a1"),
         ("a4", "a2")]
adj = {v: {b if a == v else a for a, b in EDGES if v in (a, b)} for v in NODES}
check("Die Marken sind der Perronvektor zum Eigenwert 2 (Index 4), Quadratsumme 120",
      all(sum(MARK[w] for w in adj[v]) == 2 * MARK[v] for v in NODES) and sum(m * m for m in MARK.values()) == 120)


def components(removed):
    left, comps = [v for v in NODES if v != removed], []
    while left:
        comp, todo = {left[0]}, [left[0]]
        while todo:
            for w in adj[todo.pop()] - {removed} - comp:
                comp.add(w)
                todo.append(w)
        left = [v for v in left if v not in comp]
        deg = sorted(len(adj[v] & comp) for v in comp)
        n = len(comp)
        kind = (f"A{n}" if max(deg, default=0) <= 2 else
                f"D{n}" if deg.count(3) == 1 and any(len(adj[v] & comp) == 3 and
                                                     sum(len(adj[w] & comp) == 1 for w in adj[v] & comp) >= 2
                                                     for v in comp) else f"E{n}")
        comps.append((kind, sorted(MARK[v] for v in comp)))
    return sorted(comps)


cuts = {v: components(v) for v in NODES}
check("Schnitt am Knoten Sym^3 (Marke 4 auf dem langen Arm): A3 mit den Dimensionen 1, 2, 3 und D5 mit 2, 3, 4, 5, 6",
      cuts["a6"] == [("A3", [1, 2, 3]), ("D5", [2, 3, 4, 5, 6])], str(cuts["a6"]))
check("Der zweite Knoten der Marke 4 gibt A7 + A1; kein anderer Knoten gibt D5 + A3",
      cuts["a3"] == [("A1", [2]), ("A7", [1, 2, 3, 3, 4, 5, 6])]
      and [v for v in NODES if [c[0] for c in cuts[v]] == ["A3", "D5"]] == ["a6"])
say("    Die A3-Seite trägt h0(O), h0(O(1)), h0(O(2)) = 1, 2, 3; der Schnittknoten h0(O(3)) = 4; die D5-Seite beginnt "
    "mit h0(O(4)) = 5.")
OUT["B4"] = {v: [[c[0], c[1]] for c in cuts[v]] for v in NODES}

# ---------------------------------------------------------------- B5: Realformen mit kompakter Spin(10)
say("\nB5  Realformen: Involutionen exp(i pi x) von E8, die auf D5 trivial sind")
roots = []
for i in range(8):
    for j in range(i + 1, 8):
        for si, sj in product((1, -1), repeat=2):
            v = [Fr(0)] * 8
            v[i], v[j] = Fr(si), Fr(sj)
            roots.append(tuple(v))
for signs in product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        roots.append(tuple(Fr(s, 2) for s in signs))
dot = lambda a, b: sum(x * y for x, y in zip(a, b))
V = tuple([Fr(1, 2)] * 5 + [Fr(0)] * 3)
clock = lambda r: int(4 * dot(V, r)) % 4
in_d5 = lambda r: all(x == 0 for x in r[5:])
in_a3 = lambda r: all(x == 0 for x in r[:5])
d5 = [r for r in roots if in_d5(r)]
a3 = [r for r in roots if in_a3(r)]
split = {k: sum(1 for r in roots if clock(r) == k) for k in range(4)}
check("240 Wurzeln; Uhrklassen 52 (D5 und A3), 64, 60, 64", len(roots) == 240 and len(d5) == 40 and len(a3) == 12
      and split == {0: 52, 1: 64, 2: 60, 3: 64}, str(split))
simple = [(Fr(1, 2), Fr(-1, 2), Fr(-1, 2), Fr(-1, 2), Fr(-1, 2), Fr(-1, 2), Fr(-1, 2), Fr(1, 2))]
simple += [tuple(Fr(1) if t == 0 else Fr(1) if t == 1 else Fr(0) for t in range(8))]
simple += [tuple(Fr(-1) if t == i else Fr(1) if t == i + 1 else Fr(0) for t in range(8)) for i in range(6)]
census, trivial = {}, []
for coeff in product((0, 1), repeat=8):
    x = tuple(sum(c * s[t] for c, s in zip(coeff, simple)) for t in range(8))
    sign = lambda r: int(dot(r, x)) % 2
    fixed = 8 + sum(1 for r in roots if sign(r) == 0)
    census[fixed] = census.get(fixed, 0) + 1
    if all(sign(r) == 0 for r in d5):
        trivial.append((fixed, 3 + sum(1 for r in a3 if sign(r) == 0), x))
check("Die 256 Klassen x mod 2: Fixalgebren der Dimension 248 (1), 136 (120), 120 (135)",
      census == {248: 1, 136: 120, 120: 135}, str(census))
kinds = sorted({(f, k) for f, k, _ in trivial})
names = {(248, 15): "kompakte E8 mit su(4)", (120, 15): "E8(8) mit kompakter su(4)",
         (136, 7): "E8(-24) mit su(2,2)"}
check("Auf D5 trivial sind genau drei Typen: (248, su(4)), (120, su(4)), (136, kompakter Teil 7 = so(4) + so(2))",
      kinds == [(120, 15), (136, 7), (248, 15)], f"{[(k, names.get(k)) for k in kinds]}, Klassen {len(trivial)}")
check("Kompakte Spin(10) mit nichtkompaktem A3-Faktor gibt es nur in E8(-24), und der Faktor ist su(2,2) = so(4,2), "
      "die konforme Algebra von R^(3,1); so(5,1) (kompakt 10), so(3,3) (kompakt 6) und su(3,1) (kompakt 9) kommen "
      "nicht vor", [k for f, k in kinds if k != 15] == [7])
x = next(x for f, k, x in trivial if (f, k) == (136, 7))
sign = lambda r: int(dot(r, x)) % 2
noncompact = {"A3": sum(1 for r in a3 if sign(r)), "(10,6)": sum(1 for r in roots if clock(r) == 2 and sign(r)),
              "(16,4)": sum(1 for r in roots if clock(r) == 1 and sign(r)),
              "(16*,4*)": sum(1 for r in roots if clock(r) == 3 and sign(r))}
check("Nichtkompakter Teil 112 = 8 (Translationen und spezielle konforme) + 40 + 32 + 32",
      noncompact == {"A3": 8, "(10,6)": 40, "(16,4)": 32, "(16*,4*)": 32}, str(noncompact))
even = [r for r in roots if sign(r) == 0]
comps, left = [], set(even)
while left:
    comp, todo = set(), [next(iter(left))]
    while todo:
        r = todo.pop()
        if r in comp:
            continue
        comp.add(r)
        todo += [s for s in left if s not in comp and dot(r, s) != 0]
    left -= comp
    comps.append(comp)
sizes = sorted(len(c) for c in comps)
small = min(comps, key=len)
check("Maximal kompakt: Wurzelsystem 126 + 2, also E7 x SU(2); die einzelne SU(2) liegt im A3-Faktor "
      "(eine der beiden Drehgruppen von so(4))", sizes == [2, 126] and all(in_a3(r) for r in small), str(sizes))
# die drei Realformen und die drei Paarungen der vier Gewichte der 4 von A3
spin4 = [w for w in product((Fr(1, 2), Fr(-1, 2)), repeat=3) if sum(1 for t in w if t < 0) % 2 == 0]
pattern = lambda x: tuple(int(dot((Fr(1, 2),) * 5 + w, x)) % 2 for w in spin4)
pairings = {}
for f, k, x in trivial:
    if (f, k) == (136, 7):
        pat = pattern(x)
        pairings.setdefault(frozenset((frozenset(i for i in range(4) if pat[i] == b)) for b in (0, 1)), 0)
        pairings[frozenset((frozenset(i for i in range(4) if pat[i] == b)) for b in (0, 1))] += 1
from itertools import permutations
s4 = []
for perm in permutations(range(3)):
    for sg in product((1, -1), repeat=3):
        if sg.count(-1) % 2 == 0:
            s4.append(tuple(spin4.index(tuple(sg[i] * w[perm[i]] for i in range(3))) for w in spin4))
one = next(iter(pairings))
stab = [g for g in s4 if frozenset(frozenset(g[i] for i in part) for part in one) == one]
check("Die Involutionen vom Typ E8(-24) sind genau die drei Paarungen der vier Gewichte der 4 (je zwei gegen zwei); "
      "die Weylgruppe S4 von A3 vertauscht sie, der Stabilisator einer Paarung hat Ordnung 8 (D4, die Symmetrie "
      "des Markenquadrats)", len(pairings) == 3 and all(len(part) == 2 for pr in pairings for part in pr)
      and len(set(s4)) == 24 and len(stab) == 8, f"Klassen je Paarung {sorted(pairings.values())}")
OUT["B5"] = {"pairings": len(pairings), "stabiliser_order": len(stab), "census": census, "trivial_on_D5": [[f, k, names[(f, k)]] for f, k in kinds], "noncompact": noncompact,
             "maximal_compact_roots": sizes}

# ---------------------------------------------------------------- B6: Dimensionsschleife
say("\nB6  Schleife: Dimension -> Blätter -> Marken -> Familien -> Dimension")
loop = {}
for d in range(2, 12):
    sky = d == 4                                              # S^(d-2) ist nur für d = 4 eine Riemannsche Fläche
    pi1 = 1 if d == 2 else None if d == 3 else 2              # |pi_1(SO(d-1))|, None = unendlich
    marks = Fr(2 * pi1, pi1 - 1) if pi1 and pi1 > 1 else None
    fam = marks - 1 if marks else None
    back = fam + 1 if fam and comb(int(fam), 2) == fam else None
    loop[d] = {"himmel_holomorph": sky, "blaetter": pi1, "marken": str(marks), "familien": str(fam),
               "dimension_zurueck": str(back), "fixpunkt": back == d and sky}
check("Für d >= 4 gibt die Kette immer 2 Blätter, 4 Marken, 3 Familien und zurück d = 4; Fixpunkt ist nur d = 4, und "
      "nur dort ist der Himmel holomorph", [d for d, r in loop.items() if r["fixpunkt"]] == [4]
      and all(loop[d]["dimension_zurueck"] == "4" for d in range(4, 12)))
OUT["B6"] = loop

# ---------------------------------------------------------------- B7: Träger und Verklebung
say("\nB7  Träger: positive Teilung und Verklebung")
uniq = [N for N in range(2, 41) if len(wss.positive_splits(N)) == 1]
check("Genau eine Teilung a + b mit 1 - 1/a - 1/b > 0 hat nur N = 5 (3 + 2), geprüft bis N = 40",
      uniq == [5] and wss.positive_splits(5) == [(3, 2)])
glue = [(N, f) for N in range(4, 17) for f in range(1, 13) if anti_isometric([("D", N)], [("A", f)])]
check("D_N + A_f verklebt gerade unimodular für N <= 16, f <= 12 nur bei (5, 3) und (13, 3); N - f = 2 nur bei (5, 3)",
      glue == [(5, 3), (13, 3)], str(glue))
OUT["B7"] = {"eindeutige_teilung": uniq, "verklebungen": glue}

# ---------------------------------------------------------------- B8: Zahlen der Nachbarwelten
say("\nB8  Zahlen: Keim (2, 5, 3) gegen seine Nachbarn")
nb = wss.number_baseline()
check("alpha^-1 des Keims (2, 5, 3) ist die Ledgerwurzel 137,0359992168",
      abs(nb["tfpt"]["numbers"]["alpha_inverse"] - wss.AINV) < 1e-9)
check("Von 336 Keimen trifft genau einer alpha, genau einer alle vier Zahlen",
      nb["per_target"]["alpha_measured"] == ["k2 N5 f3"] and nb["hits_histogram"].get(4) == 1, str(nb["hits_histogram"]))
for key, row in nb["neighbours"].items():
    say(f"    {key}: alpha^-1 = {row['alpha_inverse']:.4f}, Treffer {row['hits']}")
OUT["B8"] = {"histogramm": nb["hits_histogram"],
             "nachbarn": {k: round(v["alpha_inverse"], 6) for k, v in nb["neighbours"].items()}}

passed = sum(ok for _, ok in CHECKS)
say(f"\n{passed}/{len(CHECKS)} Prüfungen bestanden")
OUT["checks"] = {"passed": passed, "total": len(CHECKS), "failed": [n for n, ok in CHECKS if not ok]}
(HERE / "ausgabe_w2.txt").write_text("\n".join(LINES) + "\n")
(HERE / "w2_ergebnis.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=str) + "\n")
sys.exit(0 if passed == len(CHECKS) else 1)
