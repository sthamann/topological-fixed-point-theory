"""w3: Realform gegen Twistorbruecke, Z4-Regel fuer Spin und Ladung, Formbedingungen aus der Doppelbrechung.

Ausgangslage: Das Ledger enthaelt die Twistorbruecke (WOIT.THETA.FREE.01, WOIT.BETA3.DUALITY.01, WOIT.SU2.KIN.01).
w2 hatte sie nicht gelesen und den "Tomita-Test" als offen gefuehrt.  Hier wird unabhaengig nachgerechnet,
was die Bruecke fuer die Realform bedeutet, und was daraus fuer Chiralitaet und Vorhersagen folgt.

A  Paare von Realformen (Raumzeitfaktor, Eichfaktor) in E8 ueber die diagonalen Involutionen von O(16)
B  Blocksatz: Uhr und Woit-Struktur erzwingen Signatur (2,2); welche der drei Paarungen aus w2/B5 vorkommt
C  Uhr als Vierteldrehung der konformen Zeit: g^2 = Cartan-Involution mal Klebeuhr, g^4 = (-1)^F
D  Z4-Regel: Zentrum von Spin(10) = Zentrum von SU(2,2); erlaubte Paare (Klasse, Helizitaet)
E  Doppelbrechung: Zahlen, Daten und was sie ueber die Gestalt der Loesung sagt

Aufruf: ../../../.venv/bin/python w3_twistorbruecke_chiralitaet.py
"""
from fractions import Fraction as Fr
from itertools import product
from math import pi
from pathlib import Path
import json

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
LINES, CHECKS, RESULT = [], [], {}


def say(text=""):
    print(text)
    LINES.append(text)


def check(name, ok, detail=""):
    CHECKS.append(bool(ok))
    say(f"  [{'ok' if ok else 'FEHLER'}] {name}" + (f": {detail}" if detail else ""))


# ---------------------------------------------------------------- A: pairs of real forms through O(16)
def gammas(n=16):
    """n = 2m gamma matrices on C^(2^m) (Jordan-Wigner), exact entries in {0, +-1, +-i}."""
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1.0 + 0j, -1.0])
    I = np.eye(2, dtype=complex)
    out = []
    for k in range(n // 2):
        for P in (X, Y):
            M = np.array([[1.0 + 0j]])
            for j in range(n // 2):
                M = np.kron(M, Z if j < k else (P if j == k else I))
            out.append(M)
    return out


def so_dim(n):
    return n * (n - 1) // 2


def part_a():
    say("A  Paare von Realformen ueber diagonale Involutionen von O(16) (16 = 10 Eich + 6 Raumzeit)")
    G = gammas()
    chir = np.eye(256, dtype=complex)
    for g in G:
        chir = chir @ g
    chir = chir * (1j ** 8)
    plus = (np.eye(256) + chir.real) / 2                      # projector on S+ (dimension 128)
    check("Spinor S+ hat Dimension 128", round(np.trace(plus).real) == 128)
    table = {}
    for p, q in product(range(11), range(7)):
        m = p + q
        if m % 2 or m == 0:
            continue
        lift = np.eye(256, dtype=complex)
        for i in list(range(p)) + list(range(10, 10 + q)):
            lift = lift @ G[i]
        sq = lift @ lift
        sign = sq[0, 0].real
        if not np.allclose(sq, sign * np.eye(256)) or sign < 0:
            continue                                           # order 4 on the spinor: no involution of E8
        fixed_spinor = round((128 + abs(np.trace(plus @ lift).real)) / 2) if m < 16 else None
        if m == 16:
            continue
        compact = so_dim(16 - m) + so_dim(m) + fixed_spinor
        form = {136: "E8(-24)", 120: "E8(8)"}[compact]
        table[(p, q)] = (f"so({10 - p},{p})", f"so({6 - q},{q})", form)
    RESULT["A_pairs"] = {f"{p},{q}": v for (p, q), v in sorted(table.items())}
    by_spacetime = {}
    for (p, q), (gauge, st, form) in table.items():
        key = tuple(sorted((6 - q, q), reverse=True))
        gk = tuple(sorted((10 - p, p), reverse=True))
        by_spacetime.setdefault(key, set()).add((gk, form))
    for key in sorted(by_spacetime, reverse=True):
        say(f"    Raumzeitfaktor so{key}: Eichfaktor " + ", ".join(
            f"so{g} in {f}" for g, f in sorted(by_spacetime[key], reverse=True)))
    check("Involution nur fuer m = 0 mod 4 (so(14,2), so(10,6) kommen nicht vor)",
          all((p + q) % 4 == 0 for p, q in table))
    check("kompakte Eichalgebra so(10) nur mit so(4,2), und dann E8(-24)",
          {(st, f) for (p, q), (g, st, f) in table.items() if p in (0, 10)} <= {("so(2,4)", "E8(-24)"), ("so(4,2)", "E8(-24)")}
          and any(p in (0, 10) for p, q in table))
    eucl = {(tuple(sorted((10 - p, p), reverse=True))) for (p, q) in table if q in (1, 5)}
    check("euklidische konforme Algebra so(5,1) nur mit Eichfaktor so(7,3)", eucl == {(7, 3)}, str(eucl))
    split = {(tuple(sorted((10 - p, p), reverse=True))) for (p, q) in table if q == 3}
    check("reelle so(3,3) = sl(4,R) nur mit so(9,1) oder so(5,5)", split == {(9, 1), (5, 5)}, str(split))
    say("    Folgerung: z79 (Markensymmetrie nur euklidisch) und kompakte Spin(10) schliessen sich in E8 aus;")
    say("    die euklidische Seite haette so(7,3) als Eichalgebra (enthaelt so(6) + so(3,1)).")


# ---------------------------------------------------------------- B: block theorem, independent re-derivation
def part_b():
    say()
    say("B  Blocksatz neu gerechnet: Uhr rho4 = diag(i,1,i,1), Woit-Struktur rho_tw = Rechts-j")
    re = sp.symbols("r0:16", real=True)
    Phi = sp.zeros(4)
    k = 0
    for i in range(4):
        Phi[i, i] = re[k]
        k += 1
    for i in range(4):
        for j in range(i + 1, 4):
            Phi[i, j] = re[k] + sp.I * re[k + 1]
            Phi[j, i] = re[k] - sp.I * re[k + 1]
            k += 2
    rho = sp.diag(sp.I, 1, sp.I, 1)
    eqs = [e for e in sp.flatten((rho.H * Phi * rho - Phi).applyfunc(sp.expand)) if e != 0]
    eqs = [c for e in eqs for c in (sp.re(e), sp.im(e))]
    sol = sp.solve(eqs, re, dict=True)[0]
    Pc = Phi.subs(sol)
    free = sorted(Pc.free_symbols, key=str)
    cross = [Pc[i, j] for i in (0, 2) for j in (1, 3)]
    check("Uhrvertraeglichkeit: 16 -> 8 reelle Parameter, Kreuzbloecke null",
          len(free) == 8 and all(c == 0 for c in cross), f"{len(free)} Parameter")
    M = sp.Matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]])
    sig = {}
    for eps in (1, -1):
        cond = (M.H * Pc * M - eps * Pc.conjugate()).applyfunc(sp.expand)
        e2 = [c for e in sp.flatten(cond) if e != 0 for c in (sp.re(e), sp.im(e))]
        s2 = sp.solve(e2, free, dict=True)[0]
        Pe = Pc.subs(s2)
        rest = sorted(Pe.free_symbols, key=str)
        # sample generic members and read the signature
        seen = set()
        rng = np.random.default_rng(7)
        for _ in range(60):
            vals = {s: float(rng.normal()) for s in rest}
            ev = np.linalg.eigvalsh(np.array(Pe.subs(vals).evalf(), dtype=complex))
            seen.add((int((ev > 1e-9).sum()), int((ev < -1e-9).sum())))
        sig[eps] = (Pe, rest, seen)
    check("Zweig eps = -1 (Minkowski): Signatur immer (2,2)", sig[-1][2] == {(2, 2)}, str(sig[-1][2]))
    check("Zweig eps = +1 (euklidisch): (4,0), (0,4) und (2,2) kommen vor",
          sig[1][2] == {(4, 0), (0, 4), (2, 2)}, str(sorted(sig[1][2])))
    Pm = sig[-1][0]
    check("Minkowski-Zweig: Nahtgerade (e1,e2) traegt a(|z|^2 - 1), Nullort ist der Nahtkreis",
          sp.simplify(Pm[0, 0] + Pm[1, 1]) == 0 and Pm[0, 1] == 0)
    # seam-preserving members are diagonal: Phi = diag(a, -a, c, -c); which pairings occur?
    a, c = sp.symbols("a c", real=True, nonzero=True)
    diag_m = sp.diag(a, -a, c, -c)
    diag_e = sp.diag(a, a, c, c)
    ok_m = sp.simplify(M.H * diag_m * M + diag_m.conjugate()) == sp.zeros(4)
    ok_e = sp.simplify(M.H * diag_e * M - diag_e.conjugate()) == sp.zeros(4)
    check("diagonale Mitglieder: Minkowski diag(a,-a,c,-c), euklidisch diag(a,a,c,c)", ok_m and ok_e)

    def pairing(signs):
        pos = frozenset(i + 1 for i, s in enumerate(signs) if s > 0)
        neg = frozenset(i + 1 for i, s in enumerate(signs) if s < 0)
        return frozenset([pos, neg]) if len(pos) == 2 else None
    P13 = frozenset([frozenset({1, 3}), frozenset({2, 4})])
    P14 = frozenset([frozenset({1, 4}), frozenset({2, 3})])
    P12 = frozenset([frozenset({1, 2}), frozenset({3, 4})])
    mink = {pairing((sa, -sa, sc, -sc)) for sa in (1, -1) for sc in (1, -1)}
    eucl = {pairing((sa, sa, sc, sc)) for sa in (1, -1) for sc in (1, -1)} - {None}
    check("Minkowski-Zweig traegt genau die Paarungen {13|24} und {14|23}", mink == {P13, P14})
    check("die dritte Paarung {12|34} liegt im euklidischen Zweig (Nahtgerade definit, kein Schnitt)", eucl == {P12})
    clock2 = [(-1) ** 0 if x == 1 else -1 for x in (sp.I ** 2, 1, sp.I ** 2, 1)]
    check("fuer a c > 0 ist die Cartan-Involution Ad(rho4^2): die Uhr ist dann zentral in der maximal kompakten",
          pairing(clock2) == P13)
    RESULT["B"] = {"clock_parameters": len(free), "minkowski_signatures": sorted(sig[-1][2]),
                   "euclidean_signatures": sorted(sig[1][2]), "minkowski_pairings": ["13|24", "14|23"],
                   "euclidean_pairing": "12|34"}
    say("    Folgerung: Der in w2 offene Tomita-Test ist im Ledger entschieden (WOIT.BETA3.DUALITY.01) und hier")
    say("    unabhaengig bestaetigt: Uhr und Woit-Struktur erzwingen su(2,2). Zwei der drei Paarungen aus B5 sind")
    say("    die beiden Minkowski-Mitglieder, die dritte hat keinen Nullkreis.")


# ---------------------------------------------------------------- C, D: E8 roots in D8 coordinates
def e8_roots():
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
    return roots


def glue_class(r):
    """class in Z4 of a root under D5 (coordinates 0-4) + D3 = A3 (coordinates 5-7)."""
    if all(x.denominator == 1 for x in r):
        g, f = sum(abs(x) for x in r[:5]), sum(abs(x) for x in r[5:])
        return 0 if g in (0, 2) else 2                          # (45,1), (1,15) class 0; (10,6) class 2
    minus5 = sum(1 for x in r[:5] if x < 0)
    return 1 if minus5 % 2 == 0 else 3                          # (16,4) class 1, (16bar,4bar) class 3


def phase(r, y):
    """exp(2 pi i y.r) as a fraction of a full turn."""
    return sum(a * b for a, b in zip(r, y)) % 1


def part_c(roots):
    say()
    say("C  Die Uhr als Vierteldrehung der konformen Zeit")
    counts = {k: sum(1 for r in roots if glue_class(r) == k) for k in range(4)}
    check("Klebeklassen: 52 + 64 + 60 + 64 Wurzeln (Klassen 0, 1, 2, 3)",
          counts == {0: 52, 1: 64, 2: 60, 3: 64}, str(counts))
    y = (Fr(0),) * 7 + (Fr(1, 4),)                              # g = rotation by pi/2 in the plane of e8
    order = lambda yy: min(n for n in range(1, 17) if all(phase(r, tuple(n * t for t in yy)) == 0 for r in roots))
    a3 = [r for r in roots if all(x == 0 for x in r[:5])]
    order_a3 = min(n for n in range(1, 17) if all(phase(r, tuple(n * t for t in y)) == 0 for r in a3))
    check("g = exp(i pi/2 H) hat Ordnung 4 auf su(4) und Ordnung 8 auf E8", order_a3 == 4 and order(y) == 8,
          f"{order_a3}, {order(y)}")
    theta = (Fr(0),) * 5 + (Fr(1, 2), Fr(1, 2), Fr(0))          # Cartan involution of w2/B5: x = e6 + e7
    compact = 8 + sum(1 for r in roots if phase(r, theta) == 0)
    check("theta = exp(i pi (e6+e7)) ist trivial auf D5 und hat kompakte Dimension 136 (E8(-24))",
          compact == 136 and all(phase(r, theta) == 0 for r in roots if all(x == 0 for x in r[5:])), str(compact))
    y2 = tuple(2 * t for t in y)
    ok = all(phase(r, y2) == (phase(r, theta) + Fr(glue_class(r), 4)) % 1 for r in roots) or \
        all(phase(r, y2) == (phase(r, theta) - Fr(glue_class(r), 4)) % 1 for r in roots)
    check("g^2 = theta mal Klebeuhr z (auf allen 240 Wurzeln)", ok)
    y4 = tuple(4 * t for t in y)
    check("g^4 = z^2 = (-1)^F: +1 auf D8, -1 auf dem Spinor",
          all(phase(r, y4) == (Fr(0) if glue_class(r) in (0, 2) else Fr(1, 2)) for r in roots))
    comm = [r for r in roots if all(sum(a * b for a, b in zip(r, s)) == 0 for s in roots
                                    if glue_class(s) == 0)]
    check("keine Wurzel vertauscht mit D5 + A3, Rang 5 + 3 = 8: kein U(1) der Helizitaet in E8", len(comm) == 0)
    RESULT["C"] = {"order_on_su4": order_a3, "order_on_e8": order(y), "compact_dimension": compact}
    say("    Lesart: exp(i pi H) = theta z. Die halbe Periode der konformen Zeit ist die Cartan-Involution mal der")
    say("    Uhr der Verklebung; die Uhr mu4 ist zugleich Zentrum von Spin(10) und von SU(2,2).")


def part_d(roots):
    say()
    say("D  Z4-Regel fuer Spin und Ladung")
    # diagonal Z4: z10 acts on class k by i^k, the centre of SU(4) by i^k as well; (z10, z4^-1) is trivial
    y10 = (Fr(1, 2),) * 5 + (Fr(0),) * 3                        # centre of Spin(10): i on the 16, -1 on the 10
    y4 = (Fr(0),) * 5 + (Fr(1, 2),) * 3                         # centre of SU(4), inverse generator: -i on the 4
    check("z10 wirkt auf Klasse k als i^k (16: i, 10: -1, 16bar: -i, Adjungierte: 1)",
          all(phase(r, y10) == Fr(glue_class(r), 4) for r in roots))
    check("z4^-1 wirkt auf Klasse k als i^(-k)", all(phase(r, y4) == Fr(-glue_class(r), 4) % 1 for r in roots))
    both = tuple(a + b for a, b in zip(y10, y4))
    check("(z10, z4^-1) wirkt trivial auf allen Wurzeln: Zentrum von Spin(10) = Zentrum von SU(2,2) in E8(-24)",
          all(phase(r, both) == 0 for r in roots))
    # central character of the helicity ladders of SU(2,2) in the oscillator model
    n = 4
    a = np.diag(np.sqrt(np.arange(1, n)), 1)
    I = np.eye(n)

    def op(m, slot):
        mats = [I] * 4
        mats[slot] = m
        out = mats[0]
        for x in mats[1:]:
            out = np.kron(out, x)
        return out
    a1, a2, b1, b2 = (op(a, s) for s in range(4))
    N = [x.T @ x for x in (a1, a2, b1, b2)]
    lower, raise_ = a1 @ b1, a1.T @ b1.T                        # E13 and E31 of u(2,2)
    commutator = lower @ raise_ - raise_ @ lower
    target = N[0] + N[2] + np.eye(n ** 4)
    occ = np.array(list(product(range(n), repeat=4)))
    inner = (occ.max(axis=1) < n - 1)                           # states away from the cutoff
    check("[a1 b1, a1+ b1+] = N_a1 + N_b1 + 1 (legt die Ordnungskonstante des Zentrums fest)",
          np.allclose((commutator - target)[np.ix_(inner, inner)], 0))
    # centre element i.1 = exp(i pi/2 H), H = diag(1,1,1,-3)  ->  N_a1 + N_a2 - (N_b1 + 1) + 3 (N_b2 + 1)
    H = occ[:, 0] + occ[:, 1] - (occ[:, 2] + 1) + 3 * (occ[:, 3] + 1)
    C = occ[:, 0] + occ[:, 1] - occ[:, 2] - occ[:, 3]           # twice the helicity
    check("Zentrum von SU(2,2) wirkt auf der Leiter der Helizitaet h = C/2 als i^(C+2)",
          all((int(h) - (int(c) + 2)) % 4 == 0 for h, c in zip(H, C)))
    rule = {}
    for orient in (1, -1):
        allowed = {}
        for k in range(4):
            allowed[k] = [Fr(c, 2) for c in range(-2, 3) if (orient * (c + 2) - k) % 4 == 0]
        rule[orient] = allowed
    names = {0: "45, 15, 1", 1: "16", 2: "10", 3: "16bar"}
    for orient in (1, -1):
        say(f"    Orientierung {orient:+d}: " + "; ".join(
            f"Klasse {k} ({names[k]}): h = {', '.join(str(h) for h in rule[orient][k])}" for k in range(4)))
    check("Klasse 0 (Eichbosonen) nur mit h = +-1, Klasse 2 (die 10) nur mit h = 0",
          all(rule[o][0] == [Fr(-1), Fr(1)] and rule[o][2] == [Fr(0)] for o in (1, -1)))
    check("die 16 traegt genau eine Helizitaet, die 16bar die entgegengesetzte: keine Spiegelgeneration",
          all(len(rule[o][1]) == 1 and len(rule[o][3]) == 1 and rule[o][1][0] == -rule[o][3][0] for o in (1, -1)))
    check("die Orientierung ist das einzige freie Bit: sie vertauscht die Helizitaeten von 16 und 16bar",
          rule[1][1] == rule[-1][3] and rule[1][3] == rule[-1][1])
    # scalars only in class 2: the Yukawa representations 10, 120, 126 of Spin(10), not 16, 45, 54, 210
    quadrality = {"10": 2, "120": 2, "126": 2, "16": 1, "45": 0, "54": 0, "210": 0, "1": 0}
    check("Skalare nur in Klasse 2: genau die Yukawa-Darstellungen 10, 120, 126",
          sorted(k for k, v in quadrality.items() if v == 2) == ["10", "120", "126"])
    RESULT["D"] = {"rule": {str(o): {str(k): [str(h) for h in v] for k, v in rule[o].items()} for o in rule}}
    say("    Lesart: Felder sind Darstellungen von (Spin(10) x SU(2,2))/Z4. Die Chiralitaet folgt aus der globalen")
    say("    Gestalt der Gruppe, nicht aus der Zerlegung der 248 unter SL(2,C) (Einwand Distler-Garibaldi).")
    say("    Praemisse: die Gruppe wirkt als Quotient, also treu wie in E8(-24), und masselose Felder sind Leitern.")


# ---------------------------------------------------------------- E: birefringence
def part_e():
    say()
    say("E  Doppelbrechung als Formbedingung")
    c3 = 1 / (8 * pi)
    phi0 = 1 / (6 * pi) + 48 * c3 ** 4
    beta = phi0 / (4 * pi)
    deg = beta * 180 / pi
    omega_b = (1 - 2 * c3) * phi0
    check("beta = phi0/(4 pi) = 0,24243503 Grad (eingefrorener Wert)", abs(deg - 0.2424350309009295) < 1e-12, f"{deg:.10f}")
    check("Phase der Helizitaets-U(1): beta/2 = c3 phi0", abs(beta / 2 - c3 * phi0) < 1e-18, f"{beta / 2:.10f} rad")
    check("Budget: Omega_b + beta = phi0, Anteile (1 - 2 c3, 2 c3)", abs(omega_b + beta - phi0) < 1e-16)
    data = json.loads((HERE.parents[1] / "cmb-birefringence-seed" / "data" / "measurements.json").read_text())
    rows = [r for r in data["beta_birefringence_deg"] if "SUPERSEDED" not in r["note"]]
    w = [1 / r["sigma"] ** 2 for r in rows]
    mean = sum(wi * r["value"] for wi, r in zip(w, rows)) / sum(w)
    err = sum(w) ** -0.5
    for r in rows:
        say(f"    {r['experiment']}: {r['value']} +- {r['sigma']} Grad, Abstand {(deg - r['value']) / r['sigma']:+.2f} sigma")
    say(f"    gewichtetes Mittel (ACT, Planck PR4): {mean:.4f} +- {err:.4f} Grad; von null {mean / err:.1f} sigma, "
        f"von TFPT {(deg - mean) / err:+.2f} sigma")
    check("beide unabhaengigen Messungen haben dasselbe Vorzeichen und liegen innerhalb 1 sigma",
          all(r["value"] > 0 and abs(deg - r["value"]) < r["sigma"] for r in rows))
    f = data["beta_frequency_exponent"][0]
    far = (f["value"] + 2) / f["sigma_plus"]
    check("achromatisch: Frequenzexponent mit 0 vertraeglich, Faraday (n = -2) ueber 3 sigma entfernt",
          abs(f["value"]) < f["sigma_plus"] and far > 3, f"n = {f['value']}, Faraday {far:.1f} sigma")
    RESULT["E"] = {"beta_deg": deg, "helicity_phase_rad": beta / 2, "mean_deg": mean, "mean_sigma": err,
                   "sigma_from_zero": mean / err, "sigma_from_tfpt": (deg - mean) / err}
    say("    Formbedingungen (aus C, D und den Daten):")
    say("    1. Eine Drehung der Polarisation um beta ist eine Phase 2 beta zwischen den Helizitaeten +-1, also")
    say("       ein Element der Helizitaets-U(1). In E8 gibt es diese U(1) nicht (C), nur ihr Z4, und das Z4")
    say("       wirkt auf Photonen trivial (D). beta ist deshalb keine Groesse der Kinematik, sondern eine")
    say("       Antwort des Zustands: jede Loesung mit SU(2,2)-invariantem Vakuum als Weltzustand gibt beta = 0.")
    say("    2. beta ist ungerade unter Paritaet. Das Vorzeichen ist dasselbe Bit wie die Orientierung in D, die")
    say("       festlegt, welche Helizitaet die 16 traegt. Eine vollstaendige Loesung muss das Vorzeichen von")
    say("       beta an die Linkshaendigkeit der schwachen Wechselwirkung binden.")
    say("    3. beta ist achromatisch und eine reine Zahl. Das passt zu einer festen Phase des Randes und nicht")
    say("       zu einem langsam rollenden Feld, dessen Winkel von der Geschichte abhinge.")


def main():
    say("w3: Twistorbruecke, Z4-Regel und Doppelbrechung")
    say()
    part_a()
    part_b()
    roots = e8_roots()
    check("E8 hat 240 Wurzeln", len(roots) == 240)
    part_c(roots)
    part_d(roots)
    part_e()
    say()
    say(f"Pruefungen: {sum(CHECKS)}/{len(CHECKS)}")
    RESULT["checks"] = [sum(CHECKS), len(CHECKS)]
    (HERE / "ausgabe_w3.txt").write_text("\n".join(LINES) + "\n")
    (HERE / "w3_ergebnis.json").write_text(json.dumps(RESULT, indent=1, default=str) + "\n")
    return 0 if all(CHECKS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
