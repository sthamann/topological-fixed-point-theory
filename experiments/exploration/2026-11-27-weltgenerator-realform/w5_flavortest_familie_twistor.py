"""w5: Flavortest der Realform-Lesart.  Traegt die Familienuhr des Ledgers (mu4 auf H_1(P^1 ohne mu4) = Z^3, mit
den Transportoperatoren Q, K, R, L und der Blattparitaet Sigma) die Aussage 'family_is_twistor' aus w4?

Geprueft wird an den Matrizen des Ledgers, exakt:
  A  die Objekte des Ledgers und die Familienuhr als Coxeter-Element von W(A3)
  B  Familienuhr gegen die Twistoruhr rho4 = diag(i,1,i,1) in sl(4)
  C  welche Realform von A3 das Quadrat der Familienuhr als Cartan-Involution hat; D4 als eingeschraenkte Weyl-Gruppe
  D  Sigma, die Exponentenparitaet, die ganzzahlige D4 des Ledgers (v98, v146) und die Operatoren Q, K, R, L
  E  die Familien im Twistor-Woerterbuch des Bestands: Moebius-Algebra der Nahtgeraden in su(2,2)

Aufruf: ../../../.venv/bin/python w5_flavortest_familie_twistor.py
"""
from itertools import permutations, product
from pathlib import Path
import json

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
LINES, CHECKS, RESULT = [], [], {}
t = sp.symbols("t")
I3 = sp.eye(3)


def say(text=""):
    print(text)
    LINES.append(text)


def check(name, ok, detail=""):
    CHECKS.append(bool(ok))
    say(f"  [{'ok' if ok else 'FEHLER'}] {name}" + (f": {detail}" if detail else ""))
    return bool(ok)


# ---------------------------------------------------------------- ledger objects
R_ROT = sp.Matrix([[0, 0, -1], [1, 0, -1], [0, 1, -1]])      # GATE.QGEOM.02 (v70): mu4 on the puncture homology
Q = sp.Matrix([[3, 1, 0], [3, 2, 0], [3, 2, 1]])             # FLAV.KQ.01 (v10), transport basis
R = sp.Matrix([[1, 3, 0], [1, 5, 2], [2, 5, 3]])
SIG = sp.diag(1, -1, -1)
K = R + Q * SIG
L = R + Q * (I3 + SIG)
QP = (Q + SIG * Q * SIG) / 2
QM = (Q - SIG * Q * SIG) / 2
RHO4 = [1, 0, 1, 0]                                          # WOIT.BETA3.DUALITY.01: diag(i,1,i,1) as powers of i


def perm_matrix(p):
    """4x4 matrix of the permutation i -> p[i] of the four marks."""
    m = sp.zeros(4, 4)
    for i, j in enumerate(p):
        m[j, i] = 1
    return m


def quotient(m4):
    """Action on Z^4/Z(1,1,1,1) in the basis g1, g2, g3 with g4 = -(g1+g2+g3): the puncture homology."""
    cols = []
    for j in range(3):
        v = m4[:, j]
        cols.append(sp.Matrix([v[i] - v[3] for i in range(3)]))
    return sp.Matrix.hstack(*cols)


def signature(m):
    """Eigenvalues of an involution or of a finite-order matrix, sorted, as strings."""
    out = []
    for value, mult in m.eigenvals().items():
        out += [sp.nsimplify(value)] * mult
    return sorted(str(x) for x in out)


def cycle_type(p):
    seen, out = set(), []
    for i in range(4):
        n, j = 0, i
        while j not in seen:
            seen.add(j)
            j = p[j]
            n += 1
        if n:
            out.append(n)
    return tuple(sorted(out, reverse=True))


def ad_spectrum(exponents):
    """Multiplicities of the eigenvalues 1, i, -1, -i of Ad(diag(i^a)) on sl(4)."""
    counts = [0, 0, 0, 0]
    for a, b in product(exponents, repeat=2):
        counts[(a - b) % 4] += 1
    counts[0] -= 1
    return tuple(counts)


def generated(gens):
    group, frontier = [sp.eye(gens[0].rows)], [sp.eye(gens[0].rows)]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                x = g * h
                if x not in group:
                    group.append(x)
                    new.append(x)
        frontier = new
    return group


def split_index(mats, dim=3):
    """Index of the sum of the common eigen-lines of commuting involutions in Z^dim (1 = diagonal over Z)."""
    lines = []
    for signs in product((1, -1), repeat=len(mats)):
        space = sp.Matrix.vstack(*[m - s * sp.eye(dim) for m, s in zip(mats, signs)]).nullspace()
        for v in space:
            v = v * sp.ilcm(*[x.q for x in v])
            lines.append(v / sp.igcd(*[int(x) for x in v]))
    return abs(sp.Matrix.hstack(*lines).det()) if len(lines) == dim else None


def part_a():
    say("A  Objekte des Ledgers und die Familienuhr")
    check("K = R + Q Sigma, L = R + Q (1 + Sigma) wie im Ledger; det (Q, K, R, L) = (3, 4, 8, 20)",
          K == sp.Matrix([[4, 2, 0], [4, 3, 2], [5, 3, 2]]) and [m.det() for m in (Q, K, R, L)] == [3, 4, 8, 20])
    cyc = perm_matrix((1, 2, 3, 0))                              # g1 -> g2 -> g3 -> g4 -> g1
    check("die Drehung der vier Marken auf Z^4/Z(1,1,1,1) ist die Matrix R_rot von v70", quotient(cyc) == R_ROT)
    simple = [perm_matrix(p) for p in ((1, 0, 2, 3), (0, 2, 1, 3), (0, 1, 3, 2))]
    cox = simple[0] * simple[1] * simple[2]
    check("sie ist ein Coxeter-Element von W(A3) = S4 (Produkt der drei einfachen Spiegelungen)",
          cox in (cyc, cyc.T), "bis auf die Richtung")
    check("char = (t+1)(t^2+1): Eigenwerte i^m mit m = 1, 2, 3 = Exponenten von A3, Coxeter-Zahl 4",
          sp.factor(R_ROT.charpoly(t).as_expr()) == (t + 1) * (t**2 + 1) and R_ROT**4 == I3)
    spec = sorted(QP.eigenvals().keys())
    clock_poly = sp.expand(sp.prod([t - sp.I**m for m in spec]))
    check("Spec(Q+) = {1,2,3}; exp(i pi Q+/2) hat das charakteristische Polynom der Uhr (auf Akte: FLAV.QGEO.04)",
          spec == [1, 2, 3] and clock_poly == sp.expand(R_ROT.charpoly(t).as_expr()))
    RESULT["A"] = {"clock": "Coxeter element of W(A3) on Z^4/Z(1,1,1,1)", "spec_Q_plus": [int(x) for x in spec]}


def part_b():
    say()
    say("B  Familienuhr gegen Twistoruhr in sl(4)")
    cox_exp = [1, 0, 3, 2]                                       # Z = diag(i, 1, -i, -1)
    s_rho, s_cox = ad_spectrum(RHO4), ad_spectrum(cox_exp)
    check("Ad(rho4) auf sl(4): Eigenwerte (1, i, -1, -i) mit Vielfachheit (7, 4, 0, 4)", s_rho == (7, 4, 0, 4))
    cyc = perm_matrix((1, 2, 3, 0))
    check("die Vierer-Drehung der Marken hat in C^4 die Eigenwerte {1, i, -1, -i}, wie Z = diag(i, 1, -i, -1)",
          sp.expand(cyc.charpoly(t).as_expr()) == t**4 - 1 and sorted(cox_exp) == [0, 1, 2, 3])
    check("Ad(Coxeter) auf sl(4): Vielfachheit (3, 4, 4, 4)", s_cox == (3, 4, 4, 4))
    outer = lambda s: (s[0], s[3], s[2], s[1])
    check("die Uhren sind in Aut(sl(4)) nicht konjugiert (auch nicht mit dem aeusseren Automorphismus)",
          s_rho not in (s_cox, outer(s_cox)))
    sq = lambda e: [(2 * a) % 4 for a in e]
    check("ihre Quadrate sind dasselbe Element theta = diag(-1, 1, -1, 1), Typ (2,2)",
          sq(RHO4) == sq(cox_exp) == [2, 0, 2, 0])
    j = [0, 0, 2, 2]
    check("Z = rho4 J mit J = diag(1, 1, -1, -1)", [(a + b) % 4 for a, b in zip(RHO4, j)] == cox_exp)
    pairing = lambda e: sorted(sorted(i + 1 for i in range(4) if e[i] == v) for v in sorted(set(e)))
    theta = sq(RHO4)
    theta_j = [(a + b) % 4 for a, b in zip(theta, j)]
    pairs = [pairing(theta), pairing(j), pairing(theta_j)]
    check("theta, J, theta J sind die drei Paarungen der vier Marken",
          pairs == [[[1, 3], [2, 4]], [[1, 2], [3, 4]], [[1, 4], [2, 3]]], str(pairs))
    fixed = lambda e: ad_spectrum(e)[0]
    check("rho4 ist ein Skalar je theta-Eigenraum: Zentralisator Dimension 7 = s(gl2 + gl2), zentral in K",
          fixed(RHO4) == 7 == fixed(theta) and all(RHO4[i] == RHO4[k] for i in range(4) for k in range(4)
                                                  if theta[i] == theta[k]))
    check("die Coxeter-Klasse ist regulaer: Zentralisator Dimension 3 = Cartan", fixed(cox_exp) == 3)
    sym3 = sorted((3 - a) * 0 + a * 1 for a in range(4))         # Sym^3 of diag(i, 1): i^a
    doublet = sorted([1, 1, 0, 0])                               # diag(i, 1) (x) 1_2
    check("eine Vierteldrehung diag(i, 1) auf C^2: Sym^3 liegt in der Coxeter-Klasse, C^2 (x) C^2 in der von rho4",
          sym3 == sorted(cox_exp) and doublet == sorted(RHO4))
    # the vector 6 = Lambda^2 C^4 of SO(4,2); an SU(4) lift carries an eighth root of unity, Lambda^2 the phase i^k
    pairs6 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    six = lambda e, k: {p: (e[p[0]] + e[p[1]] + k) % 4 for p in pairs6}
    time_plane = [p for p, v in six(theta, 0).items() if v == 0]
    space = [p for p in pairs6 if p not in time_plane]
    check("theta hat auf dem Vektor 6 die Signatur (4,2); die zeitartige Ebene sind die Paare 13 und 24",
          time_plane == [(0, 2), (1, 3)] and all(six(theta, 0)[p] == 2 for p in space))
    lifts = lambda e: [k for k in range(4) if (sum(e) + 2 * k) % 4 == 0]            # det = 1 in SU(4)
    turn = {}
    for name, e in (("rho4", RHO4), ("family", cox_exp)):
        for k in lifts(e):
            v = six(e, k)
            if sorted(v[p] for p in time_plane) == [1, 3]:
                turn.setdefault(name, []).append(sorted(v[p] for p in space))
    check("Twistoruhr: Vierteldrehung der zeitartigen Ebene, Raum fest (bis auf das Zentrum)",
          [0, 0, 0, 0] in turn["rho4"], str(turn["rho4"]))
    check("Familienuhr: dieselbe Vierteldrehung der Zeitebene und eine halbe Drehung in einer raumartigen Ebene",
          [0, 0, 2, 2] in turn["family"], str(turn["family"]))
    check("beide Uhren brauchen in SU(4) eine achte Einheitswurzel: Uhr^2 = theta mal Zentrum, Uhr^4 = -1 auf der 4",
          lifts(RHO4) == [1, 3] == lifts(cox_exp))
    RESULT["B"] = {"ad_spectrum_rho4": s_rho, "ad_spectrum_family_clock": s_cox, "conjugate": False,
                   "same_square": True, "pairings": pairs,
                   "on_vector_6": {"rho4": "quarter turn of the time-like plane", "family": "the same quarter turn and a half turn in one space-like plane"}}
    say("    Lesart: beide Uhren sind vierte Wurzeln derselben Involution theta. rho4 ist die zentrale Wurzel in K")
    say("    (Vierteldrehung der konformen Zeit), die Familienuhr die regulaere (Coxeter). Sie unterscheiden sich")
    say("    um die Involution J der dritten Paarung: eine raeumliche halbe Drehung.")


def part_c():
    say()
    say("C  Realform von A3 aus dem Quadrat der Familienuhr")
    s4 = {p: perm_matrix(p) for p in permutations(range(4))}
    weyl = {p: quotient(m) for p, m in s4.items()}
    squares = {}
    for p, m in weyl.items():
        sq = m * m
        squares.setdefault(cycle_type(p), set()).add(tuple(signature(sq)) if sq * sq == I3 else ("keine Involution",))
    table = {str(k): sorted(v)[0] for k, v in squares.items()}
    only = [k for k, v in squares.items() if v == {("-1", "-1", "1")}]
    check("nur die Vierer-Zykel haben als Quadrat eine nichttriviale Involution; ihre Signatur ist (1, -1, -1)",
          only == [(4,)] and all(v in ({("1", "1", "1")}, {("keine Involution",)}) for k, v in squares.items()
                                 if k != (4,)), json.dumps(table))
    forms = {}
    for q, name in ((0, "su(4)"), (1, "su(3,1)"), (2, "su(2,2)")):
        p = [0, 1, 2, 3]
        for a, b in ((0, 2), (1, 3))[:q]:
            p[a], p[b] = b, a
        m4 = perm_matrix(tuple(p))
        plus = sum(mult for val, mult in m4.eigenvals().items() if val == 1)
        forms[name] = {"type": (plus, 4 - plus), "theta_on_cartan": signature(quotient(m4)),
                       "compact_dimension": plus**2 + (4 - plus)**2 - 1, "real_rank": q, "inner": True}
    minus = -I3
    forms["sl(4,R)"] = {"theta_on_cartan": signature(minus), "real_rank": 3, "inner": False,
                        "in_weyl_group": minus in weyl.values()}
    dt = quotient(perm_matrix((2, 3, 0, 1)))
    forms["su*(4)"] = {"theta_on_cartan": signature(-dt), "real_rank": 1, "inner": False,
                       "in_weyl_group": (-dt) in weyl.values()}
    check("Cartan-Involution auf der maximal zerfallenden Cartan: su(4) (1,1,1), su(3,1) (1,1,-1), su(2,2) (1,-1,-1); "
          "sl(4,R) und su*(4) liegen nicht in W(A3)",
          forms["su(4)"]["theta_on_cartan"] == ["1", "1", "1"] and forms["su(3,1)"]["theta_on_cartan"] == ["-1", "1", "1"]
          and forms["su(2,2)"]["theta_on_cartan"] == ["-1", "-1", "1"] and forms["su(2,2)"]["type"] == (2, 2)
          and forms["su(2,2)"]["compact_dimension"] == 7
          and not forms["sl(4,R)"]["in_weyl_group"] and not forms["su*(4)"]["in_weyl_group"])
    clock = quotient(perm_matrix((1, 2, 3, 0)))
    theta = clock * clock
    match = [name for name, row in forms.items() if row["theta_on_cartan"] == signature(theta)
             and row.get("inner", False)]
    check("Uhr^2 = Cartan-Involution waehlt unter den fuenf Realformen von A3 genau su(2,2)", match == ["su(2,2)"])
    n = [0, 0, 2, 2]                                             # SU(4) lift of (13)(24): eigenvalues 1, 1, -1, -1
    fixed6 = sum(1 for a in range(4) for b in range(a + 1, 4) if (n[a] + n[b]) % 4 == 0)
    compact = 45 + ad_spectrum(n)[0] + 10 * fixed6 + 16 * n.count(0) + 16 * n.count(0)
    check("mit trivialer Wirkung auf D5 ist Uhr^2 die Cartan-Involution von E8(-24): 45 + 7 + 20 + 32 + 32 = 136",
          (ad_spectrum(n)[0], fixed6, compact) == (7, 2, 136))
    cent = [p for p, m in weyl.items() if m * theta == theta * m]
    types = sorted(cycle_type(p) for p in cent)
    check("Zentralisator von Uhr^2 in W(A3) ist D4 (Ordnung 8), erzeugt von der Uhr und einer Kantenspiegelung",
          len(cent) == 8 and types.count((4,)) == 2 and types.count((2, 2)) == 3 and types.count((2, 1, 1)) == 2,
          str(types))
    # restricted roots: project the 12 roots e_i - e_j onto the -1 eigenspace of theta = (13)(24)
    th4 = perm_matrix((2, 3, 0, 1))
    restricted = {}
    for i, k in permutations(range(4), 2):
        r = sp.zeros(4, 1)
        r[i], r[k] = 1, -1
        proj = (r - th4 * r) / 2
        restricted[tuple(proj)] = restricted.get(tuple(proj), 0) + 1
    by_norm = {}
    for vec, mult in restricted.items():
        by_norm.setdefault((sum(x * x for x in vec), mult), 0)
        by_norm[(sum(x * x for x in vec), mult)] += 1
    zero = sum(m for v, m in restricted.items() if all(x == 0 for x in v))
    check("eingeschraenkte Wurzeln: 4 lange (Norm 2, einfach) und 4 kurze (Norm 1, doppelt), keine Null: C2 mit "
          "Vielfachheiten (1, 2), m = u(1), reeller Rang 2",
          by_norm == {(2, 1): 4, (1, 2): 4} and zero == 0, str({str(k): v for k, v in by_norm.items()}))
    # D4 on the plane a = <e1 - e3, e2 - e4> against the group generated by the reflections in restricted roots
    u1, u2 = sp.Matrix([1, 0, -1, 0]), sp.Matrix([0, 1, 0, -1])

    def on_plane(m4):
        cols = []
        for u in (u1, u2):
            v = m4 * u
            cols.append(sp.Matrix([v[0], v[1]]))                # coordinates along u1 and u2
        return sp.Matrix.hstack(*cols)
    d4_plane = []
    for p in cent:
        x = on_plane(s4[p])
        if x not in d4_plane:
            d4_plane.append(x)
    refl = generated([sp.diag(-1, 1), sp.Matrix([[0, 1], [1, 0]])])
    check("D4 wirkt treu auf der Ebene a und ist dort die Weyl-Gruppe W(C2) der eingeschraenkten Wurzeln",
          len(d4_plane) == 8 and len(refl) == 8 and all(x in refl for x in d4_plane))
    c4 = perm_matrix((1, 2, 3, 0))
    fixed = sp.Matrix([1, -1, 1, -1])                            # theta = +1 inside the sum-zero space
    check("die Uhr ist -1 auf der kompakten Geraden t (B1) und eine Vierteldrehung auf a (E): Coxeter-Element von "
          "W(C2), Coxeter-Zahl 4", c4 * fixed == -fixed and on_plane(c4)**2 == -sp.eye(2))
    RESULT["C"] = {"squares_by_cycle_type": table, "real_forms": {k: {a: (list(b) if isinstance(b, tuple) else b)
                                                                    for a, b in v.items()} for k, v in forms.items()},
                   "selected": match, "centraliser_order": len(cent), "e8_compact_dimension": compact,
                   "restricted_roots": {"long": [2, 1, 4], "short": [1, 2, 4]}}
    say("    Lesart: die Familienaufspaltung B1 + E des Ledgers ist die Aufspaltung t + a der maximal zerfallenden")
    say("    Cartan von su(2,2), und die Flavorsymmetrie D4 ist ihre eingeschraenkte Weyl-Gruppe.")


def part_d():
    say()
    say("D  Sigma, Exponentenparitaet, die ganzzahlige D4 des Ledgers und die Transportoperatoren")
    T = I3 - 2 * (QP - 2 * I3)**2                               # (-1)^Q+ on the spectrum {1, 2, 3}
    check("T = (-1)^(Q+) ist ganzzahlig: T = [[-1,0,0],[0,1,0],[0,4,-1]], T^2 = 1",
          T == sp.Matrix([[-1, 0, 0], [0, 1, 0], [0, 4, -1]]) and T * T == I3)
    v4 = [SIG, T, SIG * T]
    check("Sigma und T vertauschen, Sigma != T, und Sigma, T, Sigma T haben alle die Signatur (1, -1, -1)",
          SIG * T == T * SIG and SIG != T and all(signature(m) == ["-1", "-1", "1"] for m in v4))
    basis = sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 2, 1]])        # columns e1, b = (0,1,2), e3: the cusp basis of v146
    b1, e1, e3 = sp.Matrix([0, 1, 2]), I3[:, 0], I3[:, 2]
    check("gemeinsame Eigenvektoren e1, b = (0,1,2), e3 mit Q+ = 3, 2, 1 bilden eine Basis von Z^3 (det 1)",
          QP * b1 == 2 * b1 and QP * e1 == 3 * e1 and QP * e3 == e3 and abs(basis.det()) == 1)
    # the integer deck of the ledger (v98, v141, v146): G = T_A Sigma
    TA = sp.Matrix([[0, 1, 0], [1, 0, 0], [2, -2, 1]])
    G = TA * SIG
    check("Deck des Ledgers G = T_A Sigma: Ordnung 4, char (t+1)(t^2+1), G^2 = Sigma T, <G, Sigma> = D4",
          G**4 == I3 and sp.expand(G.charpoly(t).as_expr()) == sp.expand((t + 1) * (t**2 + 1)) and G * G == SIG * T
          and len(generated([G, SIG])) == 8 and SIG * G * SIG == G.inv())
    check("die selbstkonjugierte Gerade von G (Eigenwert -1) ist Z e3 mit Q+ = 1; E = <e1, b> traegt Q+ = 3 und 2",
          G * e3 == -e3 and (G * G + I3) * e1 == sp.zeros(3, 1) and (G * G + I3) * b1 == sp.zeros(3, 1))
    # the other non-central choice: a clock r with r^2 = T, Sigma r Sigma = r^-1
    found = []
    for sign in (1, -1):
        for a, b, c, d in product(range(-3, 4), repeat=4):
            block = sp.Matrix([[a, 0, b], [0, sign, 0], [c, 0, d]])      # basis (e1, b, e3); r commutes with T
            r = basis * block * basis.inv()
            if r * r == T and SIG * r * SIG * r == I3 and all(x.is_integer for x in r):
                found.append(r)
    clocks = [r for r in found if sp.expand(r.charpoly(t).as_expr()) == sp.expand((t + 1) * (t**2 + 1))]
    check("zweite Wahl r^2 = T (Exponentenparitaet): genau r und r^-1 ganzzahlig; dann ist B1 = Z b mit Q+ = 2",
          len(clocks) == 2 and clocks[0] * clocks[1] == I3 and clocks[0] * b1 == -b1
          and len(generated([clocks[0], SIG])) == 8, str(clocks[0].tolist()) if clocks else "")
    r = clocks[0]
    # character grading against Q+: on omega_k the degree N = diag(1,2,3) is flipped by every reflection of D4
    N = sp.diag(1, 2, 3)
    iota = sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    delta = sp.diag(sp.I, -1, -sp.I)
    refl = [iota, delta * iota, delta**2 * iota, delta**3 * iota]
    check("im Charakterraum kippt jede Spiegelung von D4 den Grad: s N s^-1 = 4 - N; Sigma vertauscht mit Q+",
          all(sp.simplify(x * N * x.inv() - (4 * I3 - N)) == sp.zeros(3, 3) for x in refl) and SIG * QP == QP * SIG)
    commuting = lambda g: [x for x in generated([g, SIG]) if x * QP == QP * x]
    check("in beiden ganzzahligen D4 vertauscht genau die Vierergruppe {1, Sigma, T, Sigma T} mit Q+, die Uhr "
          "nicht: Q+ ist dort nicht der Charaktergrad (der mit der Uhr vertauscht und von Sigma gekippt wuerde)",
          all(sorted(map(str, commuting(g))) == sorted(map(str, [I3, SIG, T, SIG * T])) for g in (G, r))
          and G * QP != QP * G and r * QP != QP * r)
    # lattice types
    dts = [quotient(perm_matrix(q)) for q in ((2, 3, 0, 1), (1, 0, 3, 2))]
    root_basis = sp.Matrix([[1, 0, 0], [-1, 1, 0], [0, -1, 1], [0, 0, -1]])
    to_roots = lambda m4: (root_basis.T * root_basis).inv() * root_basis.T * m4 * root_basis
    on_roots = [to_roots(perm_matrix(q)) for q in ((2, 3, 0, 1), (1, 0, 3, 2))]
    idx = (split_index([SIG, T]), split_index(dts), split_index(on_roots))
    check("Vierergruppe: im Transportgitter ueber Z diagonal (Index 1), im Punktierungsgitter A3* Index 2, im "
          "Wurzelgitter A3 Index 4", idx == (1, 2, 4), str(idx))

    def clock_index(m):
        """Index of ker(m + 1) + ker(m^2 + 1) in Z^3, with saturated integer bases."""
        vecs = []
        for space in ((m + I3).nullspace(), (m * m + I3).nullspace()):
            mat = sp.Matrix.hstack(*[v * sp.ilcm(*[x.q for x in v]) for v in space])
            gram = smith_normal_form(mat, domain=ZZ)
            vecs.append((mat, sp.prod([abs(gram[i, i]) for i in range(mat.cols)])))
        full = abs(sp.Matrix.hstack(vecs[0][0], vecs[1][0]).det())
        return full / (vecs[0][1] * vecs[1][1])
    cyc = perm_matrix((1, 2, 3, 0))
    lam2 = sp.Matrix.hstack(*[sp.Matrix([1 if i in pair else 0 for i in range(3)])
                              for pair in ((0, 1), (0, 2), (1, 2))])     # g1+g2, g1+g3, g2+g3 in A3*
    middle = lam2.inv() * R_ROT * lam2
    ci = {"G": clock_index(G), "r": clock_index(r), "A3*": clock_index(R_ROT), "A3": clock_index(to_roots(cyc)),
          "Lambda2": clock_index(middle)}
    check("Uhrtyp B1 + E: Deck G Index 1; Uhr auf A3* (v70) Index 2, auf A3 Index 2: G ist ueber Q zur Uhr der "
          "Punktierungshomologie konjugiert, ueber Z nicht", (ci["G"], ci["r"], ci["A3*"], ci["A3"]) == (1, 1, 2, 2)
          and G.charpoly(t) == R_ROT.charpoly(t), str({k: str(v) for k, v in ci.items()}))
    check("das Gitter zwischen A3 und A3* (Gewichte von Lambda^2 C^4, Vektor 6 von SO(4,2)) ist ganzzahlig "
          "uhrstabil und hat den Typ des Transportgitters (Index 1)",
          abs(lam2.det()) == 2 and all(x.is_integer for x in middle) and ci["Lambda2"] == 1)
    ops = {"Q": Q, "K": K, "R": R, "L": L}
    disc = {n: sp.discriminant(m.charpoly(t).as_expr(), t) for n, m in ops.items()}
    fields, has_i = {}, {}
    for n, m in ops.items():
        fields[n], has_i[n] = [], False
        for f, _ in sp.factor_list(m.charpoly(t).as_expr(), t)[1]:
            deg = sp.degree(f, t)
            if deg == 1:
                fields[n].append("Q")
            elif deg == 2:
                d = sp.discriminant(f, t)
                fields[n].append(f"Q({sp.sqrt(sp.Integer(d) // max(k * k for k in range(1, 40) if d % (k * k) == 0))})")
                has_i[n] = has_i[n] or sp.sqrt(-d).is_rational
            else:
                fields[n].append(f"kubischer Koerper, Diskriminante {disc[n]}")      # holds no quadratic subfield
    check("Q, K, R, L haben einfache Eigenwerte, ihr rationaler Kommutant ist Q[M]; keiner der Faktorkoerper "
          "enthaelt i: keine rationale Matrix mit Eigenwerten +-i vertauscht mit einer von ihnen",
          all(d != 0 for d in disc.values()) and not any(has_i.values()), str(fields))
    check("R hat ein komplexes Eigenwertpaar (Diskriminante < 0), Q, K, L sind reell diagonalisierbar",
          disc["R"] < 0 and all(disc[n] > 0 for n in "QKL"), str({n: int(d) for n, d in disc.items()}))
    comm = {n: [m * g == g * m for g in (G, r, SIG, T)] for n, m in ops.items()}
    check("keiner der vier Operatoren vertauscht mit G, r, Sigma oder T", not any(any(v) for v in comm.values()))
    snf = {n: [int(smith_normal_form(m, domain=ZZ)[i, i]) for i in range(3)] for n, m in ops.items()}
    check("Kokerne zyklisch: Q Z3, K Z4, R Z8, L Z20",
          [abs(v[2]) for v in snf.values()] == [3, 4, 8, 20] and all(abs(v[0]) == abs(v[1]) == 1 for v in snf.values()),
          str(snf))
    RESULT["D"] = {"T": T.tolist(), "ledger_deck_G": G.tolist(), "G_squared": "Sigma T", "second_choice_r": r.tolist(),
                   "v4_split_index": list(idx), "clock_index": {k: str(v) for k, v in ci.items()},
                   "discriminants": {n: int(d) for n, d in disc.items()}, "factor_fields": fields, "snf": snf}
    say("    Lesart: die Vierergruppe {1, Sigma, T, Sigma T} enthaelt die Blattparitaet, die Exponentenparitaet und")
    say("    das Quadrat des Decks. Q+ ist in keiner ganzzahligen D4 der Charaktergrad, und Q, K, R, L sind keine")
    say("    D4-invarianten Operatoren: die Lesart liefert die Transportzahlen nicht.")


def part_e():
    say()
    say("E  Die Familien im Twistor-Woerterbuch des Bestands (h6h33: 4 = H^0(O(1)) (x) C^2, Familien = H^0(K(D)))")
    z = sp.symbols("z")
    omega = lambda k: z**(k - 1) / (z**4 - 1)                   # omega_k = z^(k-1) dz / (z^4 - 1)
    vec = lambda k: z**(k - 1)                                  # v_k = z^(k-1) d/dz: L_-1, L_0, L_1 of the seam line
    q = 1 / (z**4 - 1)                                          # q = dz^2 / (z^4 - 1)
    form = lambda f, m, w: sp.simplify(f.subs(z, m) * sp.diff(m, z)**w)  # pullback of a differential of weight w
    rot, inv = sp.I * z, 1 / z
    check("omega_k hat unter z -> iz den Charakter i^k und geht unter z -> 1/z in omega_(4-k) ueber (v137, v146)",
          all(sp.simplify(form(omega(k), rot, 1) - sp.I**k * omega(k)) == 0
              and sp.simplify(form(omega(k), inv, 1) - omega(4 - k)) == 0 for k in (1, 2, 3)))
    check("das quadratische Differential q = dz^2/(z^4-1) der vier Marken ist unter beiden Erzeugern ungerade",
          sp.simplify(form(q, rot, 2) + q) == 0 and sp.simplify(form(q, inv, 2) + q) == 0)
    check("omega_k = v_k eingesetzt in q, und v_k hat den Charakter i^(k-2): (-i, 1, i)",
          all(sp.simplify(omega(k) - vec(k) * q) == 0
              and sp.simplify(form(vec(k), rot, -1) - sp.I**(k - 2) * vec(k)) == 0 for k in (1, 2, 3)))
    say("    Also: H^0(K(D)) = Moebius-Algebra sl(2) der Nahtgeraden mal dem Vorzeichencharakter von q.")
    rho4 = sp.diag(sp.I, 1, sp.I, 1)
    c2 = sp.diag(sp.I, 1)
    kron = sp.kronecker_product
    sl2 = {"L0": sp.diag(1, -1) / 2, "L-": sp.Matrix([[0, 1], [0, 0]]), "L+": sp.Matrix([[0, 0], [1, 0]])}
    gens = {n: kron(sp.eye(2), x) for n, x in sl2.items()}
    eig = {n: sp.simplify((rho4 * g * rho4.inv())[list(g).index(next(x for x in g if x != 0))]
                          / next(x for x in g if x != 0)) for n, g in gens.items()}
    check("rho4 = 1 (x) diag(i,1); auf der Unteralgebra 1 (x) sl(2) wirkt Ad(rho4) mit (1, i, -i)",
          rho4 == kron(sp.eye(2), c2) and sorted(map(str, eig.values())) == ["-I", "1", "I"]
          and all(rho4 * g * rho4.inv() == eig[n] * g for n, g in gens.items()), str(eig))
    check("die Familienuhr ist -Ad(rho4) auf dieser sl(2): Eigenwerte (-1, i, -i), wie die Uhr des Ledgers",
          sp.expand(sp.prod([t + e for e in eig.values()])) == sp.expand(R_ROT.charpoly(t).as_expr()))
    h = kron(sp.eye(2), sp.diag(1, -1))                         # hermitian form, + on {e1,e3}, - on {e2,e4}
    su11 = {"compact": sp.I * sp.diag(1, -1), "boost1": sp.Matrix([[0, 1], [1, 0]]),
            "boost2": sp.Matrix([[0, -sp.I], [sp.I, 0]])}
    preserved = all((kron(sp.eye(2), x).H * h + h * kron(sp.eye(2), x)) == sp.zeros(4, 4) for x in su11.values())
    killing = [sp.sign((x * x).trace()) for x in su11.values()]
    check("mit der Form der Signatur (2,2) ist (1 (x) sl(2)) geschnitten su(2,2) die su(1,1) = so(2,1): eine "
          "kompakte Richtung, zwei nichtkompakte", preserved and h.eigenvals() == {1: 2, -1: 2} and killing == [-1, 1, 1])
    quarter = sp.exp(sp.I * sp.pi / 4) * sp.diag(sp.exp(sp.I * sp.pi / 4), sp.exp(-sp.I * sp.pi / 4))
    check("die kompakte Richtung L0 erzeugt die Uhr selbst: diag(i,1) = Phase mal exp(i pi/2 L0); im Ledger ist das "
          "die Familie omega_2 (Charakter -1, selbstkonjugiert)", sp.simplify(quarter - c2) == sp.zeros(2, 2)
          and eig["L0"] == 1)
    theta_sig = sorted(str(sp.simplify(e * e)) for e in eig.values())
    check("theta = rho4^2 hat auf den Familien die Signatur (1, -1, -1), wie in Teil C", theta_sig == ["-1", "-1", "1"])
    RESULT["E"] = {"families": "Moebius sl(2) of the seam line times the sign character of dz^2/(z^4-1)",
                   "embedding": "1 (x) sl(2) in sl(C^2 (x) C^2); real form su(1,1) = so(2,1) inside su(2,2)",
                   "family_clock": "-Ad(rho4)", "compact_direction": "L0 = omega_2"}
    say("    Lesart: im Woerterbuch des Bestands sind die drei Familien die konforme Algebra so(2,1) der Nahtgeraden")
    say("    in so(4,2), die mittlere Familie ist die Richtung der konformen Zeit, und Familienuhr und Twistoruhr")
    say("    sind dieselbe Vierteldrehung; das Vorzeichen ist der Charakter von q.")


def main():
    say("w5: Flavortest der Realform-Lesart (family_is_twistor)")
    say()
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    say()
    say("Urteil")
    say("  1. Als Cartan von sl(4) gelesen (Uhr = Weyl-Element) hat die Familienuhr die Coxeter-Klasse und ist nicht")
    say("     rho4. Im Woerterbuch des Bestands (Familien = Moebius-Algebra der Nahtgeraden) ist sie -Ad(rho4).")
    say("  2. In beiden Lesarten ist das Quadrat der Uhr eine Involution der Signatur (1, -1, -1): su(2,2), eine")
    say("     kompakte und zwei nichtkompakte Familienrichtungen.")
    say("  3. Q, K, R, L folgen nicht aus der Lesart: sie vertauschen mit keinem Element der ganzzahligen D4 ausser")
    say("     der Eins, und Q+ ist dort nicht der Charaktergrad.")
    ok = all(CHECKS)
    say()
    say(f"{sum(CHECKS)} von {len(CHECKS)} Pruefungen bestanden")
    RESULT["verdict"] = {"family_space_as_cartan_of_sl4": "clock is the Coxeter class, not rho4",
                         "family_space_as_moebius_sl2": "clock is -Ad(rho4), compact direction L0",
                         "clock_square_signature": [1, -1, -1], "real_form": "su(2,2)",
                         "transport_operators_explained": False, "checks": [sum(CHECKS), len(CHECKS)]}
    (HERE / "ausgabe_w5.txt").write_text("\n".join(LINES) + "\n")
    (HERE / "w5_ergebnis.json").write_text(json.dumps(RESULT, indent=1, default=str) + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
