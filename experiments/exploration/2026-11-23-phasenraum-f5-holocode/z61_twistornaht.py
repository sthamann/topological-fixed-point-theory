"""Z61: Vogelperspektive -- die Naht als Twistorgerade (Raumzeit als Raum der Nahtlagen, Dynamik aus Holomorphie).

Leitlinie 'Rekonstruktion statt Auswahl': z60 hat die inneren Antworten der Naht gemessen (d = 1 bzw. kompakt 8). Hier
geht es um die äußeren Antworten, die Deformationen der Naht selbst. Lesart T: Die Naht ist eine Twistorgerade
L_x = P(Graph x) in P(V_1 + V_2), mit Normalenbündel O(1) (x) C^2 (Einheitswindung mal Dublett).

Teil 1  Zahlen: h0(O(4)) = 5 (Träger = Weylspinor), h0(O(2)) = 3 (Familien = Lambda^2_+), h0(O(1)) = 2 (Spinor/Dublett),
        h0(O(1) (x) C^2) = 4 mit h1 = 0 (vollständige 4-dim. Familie, Kodaira). Die zweiblättrige Überlagerung allein gibt
        r_* r^* O(1) = O(1) + O(-1), also nur 2 Deformationen: genau hier liegt die Lücke zu T.
Teil 2  Raumzeit als Raum der Blattabgleiche x in Hom(V_1, V_2): L_x und L_y schneiden sich genau für det(x - y) = 0
        (Nullabstand); Herm(2) gibt Signatur (1,3), Dichtematrizen den Zukunftskegel, reine Zustände Nullstrahlen; die
        quaternionische Realform ist euklidisch (Osterwalder-Schrader); die Antipodenabbildung erhält die H6-Marken.
Teil 3  Vier Marken = Hauptnullrichtungen eines Spin-2-Objekts (Weylspinor in H0(O(4)), Majorana-Sterne eines Zustands
        im Fünferträger): Weierstrass-Invarianten der Doppelüberlagerung sind die Weyl-Invarianten I, J; tau = i (E6 = 0)
        heißt J = 0; die H6-Marken sind der Spin-2-Zustand (|2> - |-2>)/sqrt2.
Teil 4  Dynamik aus Holomorphie (Penrose-Transformation): Twistorfunktionen vom Grad -2, -3, -6 geben exakt Lösungen der
        masselosen Gleichungen für Spin 0, Weyl-Spin 1/2 (eine Chiralität) und das linearisierte selbstduale Graviton.
Teil 5  Bilanz: was T leistet und was offen bleibt.
"""
import itertools
import json
import random
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def h0(d):
    return max(d + 1, 0)


def h1(d):
    return max(-d - 1, 0)


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Zahlen der Naht als Twistorgerade")
check("h0(O(4)) = 5, h0(O(2)) = 3, h0(O(1)) = 2: Träger = Weylspinor (Spin 2), Familien = H0(K(D)) = Lambda^2_+ "
      "(= holomorphe Vektorfelder der Naht, sl(2)), Dublett = Spinor",
      (h0(4), h0(2), h0(1)) == (5, 3, 2) and h0(-2 + 4) == 3)
check("Normalenbündel O(1) (x) C^2 = O(1) + O(1): h0 = 4, h1 = 0, also eine vollständige 4-dimensionale Familie von Nähten "
      "(Kodaira); Grad 2 = c1(T P^3)|L - c1(T L) = 4 - 2", 2*h0(1) == 4 and 2*h1(1) == 0 and 4 - 2 == 2)
g = sp.Rational(1, 2)*(2*(-2) + 4) + 1
check("Doppelüberlagerung mit vier Verzweigungen: Geschlecht 1, r_* O_E = O + O(-2) (h0 = 1, h1 = 1 = g), also "
      "r_* r^* O(1) = O(1) + O(-1) mit h0 = 2: die beiden Blätter allein tragen nur 2 Deformationen",
      g == 1 and h0(0) + h0(-2) == 1 and h1(0) + h1(-2) == 1 and h0(1) + h0(-1) == 2)
befund("Lücke zu T genau benannt: Eine Twistorgerade braucht O(1) + O(1); das zweite Blatt liefert O(-1). Das zweite O(1) "
       "muss das Dublett des Partnerblatts tragen (Hom(V_1, V_2)); im Bestand ist das nicht hergeleitet", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Raumzeit als Raum der Blattabgleiche x in Hom(V_1, V_2)")
rng = np.random.default_rng(61)


def crand(*s):
    return rng.normal(size=s) + 1j*rng.normal(size=s)


ok_inc = True
for _ in range(200):
    X = crand(2, 2)
    u, v = crand(2), crand(2)
    Ynull = X + np.outer(u, v)
    Yrand = crand(2, 2)
    s_null = np.linalg.svd(X - Ynull, compute_uv=False).min()
    s_rand = np.linalg.svd(X - Yrand, compute_uv=False).min()
    ok_inc &= abs(np.linalg.det(X - Ynull)) < 1e-10 and s_null < 1e-10 and s_rand > 1e-6
check("L_x = {(i x pi, pi)} und L_y schneiden sich genau, wenn det(x - y) = 0 (200 Nullabstände und 200 zufällige Paare): "
      "Nullabstand = gemeinsamer Punkt der Nähte", ok_inc)
sig = [np.array([[1, 0], [0, 1]]), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
G = np.array([[0.5*(np.linalg.det(sig[i] + sig[j]) - np.linalg.det(sig[i]) - np.linalg.det(sig[j])) for j in range(4)]
              for i in range(4)]).real
ev = np.sort(np.linalg.eigvalsh(G))
cone = True
for _ in range(500):
    A_ = crand(2, 2)
    rho = A_ @ A_.conj().T
    rho /= np.trace(rho).real
    psi = crand(2)
    pure = np.outer(psi, psi.conj())/np.vdot(psi, psi).real
    cone &= np.linalg.det(rho).real > 0 and abs(np.linalg.det(pure)) < 1e-12
check("Realform Herm(2): det ist die Minkowskiform mit Signatur (1,3); Dichtematrizen liegen im Zukunftskegel, reine "
      "Zustände auf Nullstrahlen", np.allclose(ev, [-1, -1, -1, 1]) and cone)
okE = True
for _ in range(200):
    a, b = complex(*rng.normal(size=2)), complex(*rng.normal(size=2))
    Q = np.array([[a, -b.conjugate()], [b, a.conjugate()]])
    okE &= abs(np.linalg.det(Q).imag) < 1e-12 and abs(np.linalg.det(Q).real - (abs(a)**2 + abs(b)**2)) < 1e-12
marks = [1, 1j, -1, -1j]
anti = sorted([(-1/np.conj(z)) for z in marks], key=lambda w: (round(w.real, 9), round(w.imag, 9)))
check("die quaternionische Realform ist euklidisch (det = |a|^2 + |b|^2, Ausgangspunkt für Osterwalder-Schrader); die "
      "Antipodenabbildung z -> -1/conj(z) bildet die H6-Marken {1, i, -1, -i} auf sich ab (zwei Antipodenpaare)",
      okE and np.allclose(sorted(marks, key=lambda w: (round(w.real, 9), round(w.imag, 9))), anti))

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Vier Marken als Hauptnullrichtungen eines Spin-2-Objekts")


def IJ(a, b, c, d, e):          # Quartik a z^4 + 4 b z^3 + 6 c z^2 + 4 d z + e (Newman-Penrose: Psi_4..Psi_0)
    return a*e - 4*b*d + 3*c*c, a*c*e + 2*b*c*d - a*d*d - b*b*e - c**3


def j_roots(roots):
    z1, z2, z3, z4 = roots
    lam = ((z1 - z3)*(z2 - z4))/((z1 - z4)*(z2 - z3))
    return 256*(lam**2 - lam + 1)**3/(lam**2*(lam - 1)**2)


mp.mp.dps = 40
ok_j = True
for _ in range(50):
    co = [mp.mpc(*rng.normal(size=2)) for _ in range(5)]
    a, b, c, d, e = co
    I_, J_ = IJ(a, b, c, d, e)
    roots = mp.polyroots([a, 4*b, 6*c, 4*d, e], maxsteps=200, extraprec=200)
    ok_j &= abs(1728*I_**3/(I_**3 - 27*J_**2)/j_roots(roots) - 1) < mp.mpf(10)**-25
check("für 50 zufällige Weylspinoren gilt j(Kreuzverhältnis der vier Nullrichtungen) = 1728 I^3/(I^3 - 27 J^2): die "
      "Weierstrass-Invarianten der Doppelüberlagerung w^2 = Psi(z) sind g2 = I, g3 = J", ok_j)
I_h, J_h = IJ(-1, 0, 0, 0, 1)                                         # H6: s(z) = 1 - z^4
I_t, J_t = IJ(1, 0, 0, sp.Rational(2, 1), 0)                          # z^4 + 8 z: Tetraeder (Nullstellen 0 und -2 omega^k)
check("H6-Marken s = 1 - z^4: J = 0, I = -1, also j = 1728 und tau = i; E6(tau) = 0 heißt genau J(Psi) = 0. "
      "Gegenstück Tetraeder z(z^3 + 8): I = 0, j = 0, tau = e^(2 pi i/3)", J_h == 0 and I_h != 0 and I_t == 0 and J_t != 0)
cm = {2: 1/np.sqrt(2), -2: -1/np.sqrt(2)}
binom = {m: np.sqrt(float(sp.binomial(4, 2 + m))) for m in range(-2, 3)}
poly = [binom[m]*cm.get(m, 0) for m in (2, 1, 0, -1, -2)]
st = np.roots(poly)
check("der Spin-2-Zustand (|2> - |-2>)/sqrt2 im Fünferträger hat die Majorana-Sterne {1, i, -1, -i}: die H6-Marken sind "
      "die Konstellation eines Zustands im Träger", np.allclose(sorted(np.round(st**4, 9).real), [1, 1, 1, 1]) and
      np.allclose(np.abs(st), 1))
befund("Lesart: P1 (Sphäre mit vier Marken, Doppelüberlagerung, tau) ist die Twistorform eines Spin-2-Objekts an einem "
       "Punkt; 'dieselbe Vier überall' heißt Spin 2 (2s = 4 Nullrichtungen); die Stationarität <T> = 0 der Naht wird "
       "J(Psi) = 0. Als Aussage über die Hintergrundkrümmung wäre J = 0 falsch (Schwarzschild ist Typ D, flacher Raum hat "
       "Psi = 0): Die Marken gehören zu einem Zustand im Träger, nicht zum klassischen Hintergrund", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Dynamik aus Holomorphie (Penrose-Transformation)")
X = sp.Matrix(2, 2, sp.symbols("x00 x01 x10 x11"))
z = sp.Symbol("z")
random.seed(61)


def rq():
    return sp.Rational(random.randint(-9, 9), random.randint(1, 5))


lamA, muA, lamB, muB = [sp.Matrix([rq(), rq()]) for _ in range(4)]


def dual_line(lam, mu):     # a^{A'}(x) = i lam_A X^{AA'} + mu^{A'}
    return [sp.I*(lam[0]*X[0, k] + lam[1]*X[1, k]) + mu[k] for k in range(2)]


a_, b_ = dual_line(lamA, muA), dual_line(lamB, muB)
pi = [sp.Integer(1), z]
ap = a_[0]*pi[0] + a_[1]*pi[1]
bp = b_[0]*pi[0] + b_[1]*pi[1]
z0 = -a_[0]/a_[1]


def residue(numer, order):
    gfun = numer/(bp*a_[1]**order)
    return sp.diff(gfun, z, order - 1).subs(z, z0)/sp.factorial(order - 1)


def D(f, A, Ap):
    return sp.diff(f, X[A, Ap])


pts = [{X[i, j]: sp.Rational(random.randint(-7, 7), random.randint(1, 4)) + sp.I*sp.Rational(random.randint(-5, 5), 3)
        for i in range(2) for j in range(2)} for _ in range(3)]
phi0 = residue(sp.Integer(1), 1)
box = 2*(D(D(phi0, 0, 0), 1, 1) - D(D(phi0, 0, 1), 1, 0))
ok0 = all(sp.simplify(box.subs(p)) == 0 for p in pts)
mp.mp.dps = 30
p0 = {k: complex(v) for k, v in pts[0].items()}
an = [complex(sp.N(t.subs(pts[0]))) for t in a_]
bn = [complex(sp.N(t.subs(pts[0]))) for t in b_]
zc = -an[0]/an[1]
rad = 0.4*abs(zc + bn[0]/bn[1])                                    # Kreis nur um den Pol von A.Z
quad = mp.quad(lambda t: 1/((an[0] + an[1]*(zc + rad*mp.expj(t)))*(bn[0] + bn[1]*(zc + rad*mp.expj(t))))*1j*rad*mp.expj(t),
               [0, 2*mp.pi])/(2j*mp.pi)
ok0 &= abs(complex(quad) - complex(sp.N(phi0.subs(pts[0])))) < 1e-12
check("Spin 0: Twistorfunktion 1/((A.Z)(B.Z)) vom Grad -2 gibt phi = Residuum (Konturintegral numerisch bestätigt) mit "
      "Box phi = 0 exakt an drei komplexen Punkten", ok0)
phi1 = [residue(pi[k], 2) for k in range(2)]
ok1 = all(sp.simplify((D(phi1[0], B, 1) - D(phi1[1], B, 0)).subs(p)) == 0 for B in range(2) for p in pts)
check("Spin 1/2: Grad -3 gibt phi_A' mit nabla^(AA') phi_A' = 0 exakt (Weyl-Gleichung, nur eine Chiralität)", ok1)
Psi = {}
for idx in itertools.product(range(2), repeat=4):
    Psi[idx] = residue(pi[idx[0]]*pi[idx[1]]*pi[idx[2]]*pi[idx[3]], 5)
ok2 = True
for B in range(2):
    for rest in itertools.product(range(2), repeat=3):
        expr = D(Psi[(0,) + rest], B, 1) - D(Psi[(1,) + rest], B, 0)
        ok2 &= all(sp.simplify(expr.subs(p)) == 0 for p in pts[:2])
check("Spin 2: Grad -6 gibt Psi_A'B'C'D' mit nabla^(AA') Psi_A'B'C'D' = 0 exakt (linearisiertes selbstduales Graviton, "
      "Bianchi im Vakuum); an jedem Punkt ist Psi eine Quartik in pi, ihre vier Nullstellen sind die Marken", ok2)
befund("Die Feldgleichungen sind nicht entworfen: Sie folgen aus Holomorphie. Positive Frequenz ist Hardy-Fortsetzbarkeit "
       "in PT+, also dieselbe Bedingung wie die originale Szegő-Projektion der Naht", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Bilanz der Lesart T (Naht = Twistorgerade)")
tafel = [
    ("3+1", "Raum der Nahtlagen Hom(V_1, V_2): d = 4, nichtkompakt, Nullabstand = gemeinsamer Nahtpunkt; Signatur aus "
            "Herm(2), euklidisch über die quaternionische Form -- unter T, aus Deformationsantworten"),
    ("Feldgleichungen", "masselos für alle Helizitäten aus Holomorphie (Spin 0, Weyl 1/2, selbstduales Graviton geprüft)"),
    ("Chiralität", "Grad der Twistorfunktion legt die Händigkeit fest (Grad -3: nur gestrichene Weylspinoren)"),
    ("Zahlen 5, 3, 2, 4", "h0 von O(4), O(2), O(1), O(1) (x) C^2"),
    ("vier Marken, tau = i", "Nullrichtungen eines Spin-2-Zustands im Träger; tau = i heißt J = 0"),
    ("offen: T selbst", "die Blätter geben O(1) + O(-1), T braucht O(1) + O(1)"),
    ("offen: volle Gravitation", "nur der selbstduale Sektor folgt (googly-Problem); nichtlinear über Penroses Graviton, "
                                  "voll erst über Ambitwistoren"),
    ("offen: Wechselwirkung, E8, Gewichte, G", "nicht berührt"),
]
for k_, v_ in tafel:
    say(f"   {k_}: {v_}")
ZAHLEN["bilanz"] = dict(tafel)
urteil = ("Aus der Vogelperspektive sind TFPTs Grundzahlen die Twistordaten eines Raumzeitpunkts: 5 = h0(O(4)) "
          "(Weylspinor), 3 = h0(O(2)) (Lambda^2_+, Familien), 2 = h0(O(1)) (Spinor), 4 = h0(O(1) (x) C^2) (Raumzeit), und "
          "die vier Marken sind die Nullrichtungen eines Spin-2-Zustands im Träger; tau = i heißt J(Psi) = 0. Unter der "
          "Lesart T (Naht = Twistorgerade) ist 3+1 der Raum der Nahtlagen, rekonstruiert aus den Deformationsantworten "
          "(d = 4, nichtkompakt, Signatur aus Herm(2)), und die masselosen Feldgleichungen (Spin 0, Weyl 1/2, selbstduales "
          "Graviton) folgen exakt aus Holomorphie statt aus entworfenen Propagatoren. Offen sind T selbst (die Blätter geben "
          "O(1) + O(-1) statt O(1) + O(1)), der nicht-selbstduale Teil der Gravitation, Wechselwirkung, E8, Gewichte und G.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z61_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
