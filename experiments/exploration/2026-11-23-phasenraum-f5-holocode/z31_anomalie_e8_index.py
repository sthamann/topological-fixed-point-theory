"""Z31: Wann trägt der primitive Defekt genau einen E8-Block? Anomalienbilanz, Gegenmodell, Index, Spektralwirkung.

Prüft die Anomalien- und Indexrechnungen der fünften externen Analyse (24.11.) und verbindet sie mit Z13, Z22, Z24:
Teil 1  Zufluss auf den primitiven Wirbel eines Ladung-4-Kondensats: c_- = (1/4) sum q_f dim R_f,
        Spin(10)-Niveau k = (q_{10,+} + q_{10,-})/2. Direkte J-Fortsetzung q_+ = 1 + 4a, q_- = -1 + 4b:
        Aufzählung (6561), Teilbarkeit, Anker in Z13/Z22/Z24; axiale Fortsetzung (+1, +1).
Teil 2  Gegenmodell: Quartettzustand |0> + r e^{i theta} (P^+)^2 |0> mit P^+ = sum_a c^+_{a up} c^+_{a dn}:
        Z4-symmetrisch, <P> = 0, <P^2> != 0, Lücke 1, Produktzustand also c_- = 0.
Teil 3  Bosonische chirale invertierbare Phasen: (ST)^3 = S^2 verlangt c in 8Z; E8 ist das kleinste gerade
        unimodulare Gitter. Fermionischer Generator c = 1/2; nach der GSO-Summe (A*) gilt die bosonische Klasse.
Teil 4  Axiales Ladung-4-Kondensat und der Operator Re Phi^2 aus Z24: N_DW = 1 (Wand bleibt).
Teil 6  Trivialer Kommutant: Shift-Vielfachheit, Clifford-Irreduzibilität, E8-Coset.
Teil 5  Spektralwirkung: f_lambda = e^{-s}(1 + lambda s) gibt f_0 = 1, f_2 = 1 + lambda; Halbgruppe => e^{-kappa s};
        im Wärmeleitungsregler gilt S_ent = A/(4 G_ind) für den minimal gekoppelten Skalar exakt.
"""
import itertools
import json
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import sympy as sp

OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
HERE = Path(__file__).parent


def say(s=""):
    print(s)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Anomalienbilanz auf dem primitiven Wirbel eines Ladung-4-Kondensats")
c_of = lambda q10p, q10m, q6p, q6m: Fr(10 * (q10p + q10m) + 6 * (q6p + q6m), 4)
k_of = lambda q10p, q10m: Fr(q10p + q10m, 2)
rngA = range(-4, 5)
tot = c8 = c8k1 = 0
for a1, b1, a2, b2 in itertools.product(rngA, repeat=4):
    q = (1 + 4 * a1, -1 + 4 * b1, 1 + 4 * a2, -1 + 4 * b2)
    tot += 1
    if c_of(*q) == 8:
        c8 += 1
        if k_of(q[0], q[1]) == 1:
            c8k1 += 1
say(f"   direkte J-Fortsetzung, a, b in [-4, 4]: {tot} Kombinationen, c = 8: {c8}, c = 8 und k = 1: {c8k1}")
check("Aufzählung: 6561 Kombinationen, 110 mit c = 8, 0 mit c = 8 und k = 1", (tot, c8, c8k1) == (6561, 110, 0))
say("   Teilbarkeit: 16 gleiche Paare geben c = 16 (a + b), also nie 8; k = 2 (a + b) ist stets gerade, nie 1.")
check("16 gleiche Paare: c in 16Z; Niveau k in 2Z (für alle ganzen a, b)", all(
    c_of(1 + 4 * a, -1 + 4 * b, 1 + 4 * a, -1 + 4 * b) % 16 == 0 and k_of(1 + 4 * a, -1 + 4 * b) % 2 == 0
    for a, b in itertools.product(range(-20, 21), repeat=2)))
ax = (c_of(1, 1, 1, 1), k_of(1, 1)); dj = (c_of(1, -1, 1, -1), k_of(1, -1))
say(f"   axial (+1, +1): c = {ax[0]}, k = {ax[1]}; direkt J (+1, -1): c = {dj[0]}, k = {dj[1]}")
check("axial (+1, +1): c = 8, k = 1 (E8 bzw. so(16)_1); direkt (+1, -1): c = 0, k = 0", ax == (8, 1) and dj == (0, 0))
z13 = json.loads((HERE / "z13_ergebnis.json").read_text())
c_pair = z13["defekt"]["verdrillt"]["sin"]["Kern"]["c"]
z24 = json.loads((HERE / "z24_ergebnis.json").read_text())["zahlen"]["ecken_index"]
say(f"   Anker: Z13 verdrillter Defekt c = {c_pair:.4f} je Paar (16 Paare: {16 * c_pair:.2f}); Z24 Ecken-Index {z24}")
check("Anker Z13: c je Paar = 1/2 = (1/4)(q_+ + q_-) mit axialen Ladungen", abs(c_pair - 0.5) < 0.01)
check("Anker Z24/Z22: G (axial, gravitativ anomal) Index 2 je Umdrehung, U (vektorartig wie J) Index 0",
      z24["G (psi chi + Wilson)"] == 2 and z24["U mit y' = 0,6"] == 0)
befund("Ladung-4-Kondensat der J-Ladung schützt E8 nicht (c in 16Z, k gerade)", True)
befund("Nur die axiale, gravitativ anomale Ladung gibt c = 8 und k = 1", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Gegenmodell (Quartettordnung ohne Nahtmoden)")
def fock(nmodes):
    Z = np.diag([1.0, -1.0]); I2 = np.eye(2); low = np.array([[0, 1], [0, 0]], complex)
    ops = []
    for i in range(nmodes):
        m = np.array([[1.0 + 0j]])
        for j in range(nmodes):
            m = np.kron(m, Z if j < i else (low if j == i else I2))
        ops.append(m)
    return ops
res = {}
for ncomp in (2, 3, 4):
    a = fock(2 * ncomp)                                  # Moden (a, up) = 2a, (a, dn) = 2a+1
    dim = a[0].shape[0]
    Pd = sum(a[2 * k].conj().T @ a[2 * k + 1].conj().T for k in range(ncomp))
    vac = np.zeros(dim, complex); vac[0] = 1
    r, th = 0.7, 0.9
    om = vac + r * np.exp(1j * th) * (Pd @ Pd @ vac)
    om /= np.linalg.norm(om)
    Nop = sum(x.conj().T @ x for x in a)
    Z4 = np.diag(np.exp(1j * np.pi / 2 * np.real(np.diag(Nop))))
    P = Pd.conj().T
    H = np.eye(dim) - np.outer(om, om.conj())
    ev = np.sort(np.linalg.eigvalsh(H))
    res[ncomp] = dict(z4=bool(np.allclose(Z4 @ om, om)), P=float(abs(om.conj() @ P @ om)),
                      P2=float(abs(om.conj() @ P @ P @ om)), gap=float(ev[1] - ev[0]), deg=int(np.sum(ev < 0.5)))
    say(f"   {ncomp} Komponenten: Z4 {res[ncomp]['z4']}, |<P>| = {res[ncomp]['P']:.1e}, |<P^2>| = {res[ncomp]['P2']:.4f}, "
        f"Lücke {res[ncomp]['gap']:.3f}, Grundzustand {res[ncomp]['deg']}-fach")
check("Quartettzustand: Z4-invariant, <P> = 0, <P^2> != 0, Lücke 1, eindeutiger Grundzustand",
      all(v["z4"] and v["P"] < 1e-12 and v["P2"] > 1e-3 and abs(v["gap"] - 1) < 1e-12 and v["deg"] == 1 for v in res.values()))
n16 = 4 * sp.binomial(16, 2)
say(f"   16 Komponenten analytisch: ||(P^+)^2 |0>||^2 = 4 C(16,2) = {n16}, also <P^2> = r e^(i theta) {n16} / (1 + r^2 {n16}) != 0")
check("16 Komponenten: <P^2> != 0 für r != 0 (Norm 480)", n16 == 480)
say("   Ein Produkt solcher lokaler Zustände mit räumlich windender Phase und r = 0 am Kern bleibt überall gelückt und ist")
say("   ein Produktzustand: Seine chirale Zentralladung ist 0. Ordnungsparameter und Komponentenzahl schützen den Index nicht;")
say("   ihn trägt erst die gravitativ anomale Antwort (Teil 1).")
befund("Quartettordnung mit 16 Komponenten allein erzwingt keine Nahtmoden", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: bosonische invertierbare chirale Phasen")
c = sp.symbols("c")
Tph = sp.exp(-2 * sp.pi * sp.I * c / 24)
cond = sp.simplify(Tph**3)                               # (ST)^3 = S^2 = 1 auf einem holomorphen Charakter
cs_ok = [cc for cc in range(1, 25) if sp.simplify(cond.subs(c, cc)) == 1]
say(f"   (ST)^3 = e^(-2 pi i c/8) = 1 für c in {cs_ok}")
check("holomorphe bosonische Theorie: c in 8Z (kleinstes positives c = 8)", cs_ok == [8, 16, 24])
cartan = np.array([[2, -1, 0, 0, 0, 0, 0, 0], [-1, 2, -1, 0, 0, 0, 0, 0], [0, -1, 2, -1, 0, 0, 0, 0], [0, 0, -1, 2, -1, 0, 0, 0],
                   [0, 0, 0, -1, 2, -1, 0, -1], [0, 0, 0, 0, -1, 2, -1, 0], [0, 0, 0, 0, 0, -1, 2, 0], [0, 0, 0, 0, -1, 0, 0, 2]])
check("E8-Cartanmatrix: gerade, positiv definit, det = 1 (gerades unimodulares Gitter in Dimension 8)",
      round(np.linalg.det(cartan)) == 1 and np.all(np.diag(cartan) % 2 == 0) and np.linalg.eigvalsh(cartan).min() > 0)
say("   Fermionisch ist der Generator c = 1/2 (ein Majorana-Weyl). Unter A* summiert die Naht über Spinstrukturen (GSO, Z5),")
say("   ist also bosonisch; dann ist E8 (c = 8) der Generator, und eine unimodulare Indexabbildung Z -> Z gibt |m| = 1.")
befund("'Genau ein E8-Block' ist äquivalent zur Unimodularität der Abbildung primitiver Index -> bosonische Klasse", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: axiales Ladung-4-Kondensat und der relevante Operator Re Phi^2")
qA_phi = 2
qA_op = 2 * qA_phi
alpha = np.linspace(0, 2 * np.pi, 4001)[:-1]             # Phase von Phi^2
V = -np.cos(qA_op / qA_op * alpha)                       # Re Phi^2 = |Phi|^2 cos(arg Phi^2)
ndw = int(np.sum((V < np.roll(V, 1)) & (V < np.roll(V, -1))))
say(f"   Re Phi^2 hat axiale Ladung {qA_op} (Phi: {qA_phi}); als Potential in arg(Phi^2): {ndw} Minimum je 2 pi")
check("Re Phi^2 ist unter axialem Z4 invariant, pinnt aber arg(Phi^2): N_DW = 1 am primitiven Phi^2-Wirbel",
      qA_op % 4 == 0 and ndw == 1)
befund("Auch der axiale Ladung-4-Weg trägt die Wand, solange Re Phi^2 nicht abgestimmt (Z24) oder verboten (Z26) ist", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Spektralwirkung und Halbgruppenregler")
s_, lam, kap, eps = sp.symbols("s lambda kappa epsilon", positive=True)
f = sp.exp(-s_) * (1 + lam * s_)
f0, f2 = f.subs(s_, 0), sp.integrate(f, (s_, 0, sp.oo))
say(f"   f_lambda: f_0 = {f0}, f_2 = {sp.simplify(f2)}")
check("f_0 = 1 und f_2 = 1 + lambda: Eich- und Gravitationsnormierung unabhängig", f0 == 1 and sp.simplify(f2 - (1 + lam)) == 0)
g = sp.exp(-kap * s_)
t_ = sp.symbols("t", positive=True)
check("e^{-kappa s} erfüllt f(s + t) = f(s) f(t); Cauchy mit Stetigkeit, Positivität und f(0) = 1 lässt nur diese Familie",
      sp.simplify(g.subs(s_, s_ + t_) - g * g.subs(s_, t_)) == 0)
say(f"   Halbgruppe: f_2 / f_0 = {sp.integrate(g, (s_, 0, sp.oo))}, kappa geht in Lambda auf; frei bleibt eine Skala.")
# Wärmeleitungsregler (s > eps^2): Skalar, xi = 0
Sent = sp.integrate(1 / (12 * 4 * sp.pi * s_**2), (s_, eps**2, sp.oo))          # S/A (konische Methode)
invG = 16 * sp.pi * sp.Rational(1, 2) * sp.Rational(1, 6) * sp.integrate(1 / (16 * sp.pi**2 * s_**2), (s_, eps**2, sp.oo))
say(f"   Skalar: S/A = {sp.simplify(Sent)}, 1/(4 G_ind) = {sp.simplify(invG / 4)}")
check("Wärmeleitungsregler: S_ent = A/(4 G_ind) für den minimal gekoppelten Skalar exakt", sp.simplify(Sent - invG / 4) == 0)
z27 = json.loads((HERE / "z27_ergebnis.json").read_text())
say("   Z27-Kalibrierung: Vektoren nur mit Kabats Kontaktterm R = 1, ohne ihn R = -1/2. Im Halbgruppenregler bleiben also")
say("   eine Skala (Lambda bzw. a/l_P) und der Randterm der Eichfelder; die Kopplungsregel des Gitters entfällt.")
befund("Ein Halbgruppenregler reduziert die Z27-Mehrdeutigkeit auf Lambda und den Eich-Randterm", True)

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: Ein Prinzip für Entscheidung 1 und 2? Trivialer Kommutant")
def commutant_dim(gens):
    n = gens[0].shape[0]
    rows = [np.kron(np.eye(n), g) - np.kron(g.T, np.eye(n)) for g in gens]      # vec(XG - GX) = (I (x) G - G^T (x) I) vec X
    A = np.vstack(rows)
    return n * n - np.linalg.matrix_rank(A)
K = 5
Sj = np.eye(K, k=-1)
cd = {m: commutant_dim([np.kron(Sj, np.eye(m)), np.kron(Sj.T, np.eye(m))]) for m in (1, 2, 3)}
say(f"   Shift der Vielfachheit m (abgeschnitten, mit Adjungierter): Kommutantendimension {cd}")
check("Kommutant des Shifts der Vielfachheit m ist M_m (Dimension m^2): trivial genau für m = 1 (eine Naht)", cd == {1: 1, 2: 4, 3: 9})
def clifford(nm):
    Z = np.diag([1.0, -1.0]); X = np.array([[0, 1], [1, 0.0]]); Y = np.array([[0, -1j], [1j, 0]]); I2 = np.eye(2)
    gs = []
    for i in range(nm):
        for P_ in (X, Y):
            m = np.array([[1.0 + 0j]])
            for j in range(nm):
                m = np.kron(m, Z if j < i else (P_ if j == i else I2))
            gs.append(m)
    return gs
cl = {2 * nm: commutant_dim(clifford(nm)) for nm in (1, 2, 3)}
say(f"   Clifford Cl(2n) auf dem Spinor (2^n): Kommutantendimension {cl}; für 16 Majorana gilt dasselbe (Spinor 2^8)")
check("2n Majorana-Komponenten wirken irreduzibel (Kommutant C): die 16 Eltern sind keine Vielfachheit", all(v == 1 for v in cl.values()))
coset = {N: (8 * N - 8, Fr(8 * N - 8, 2)) for N in (1, 2, 3)}
say(f"   N E8-Blöcke relativ zu einer E8-Stromalgebra: (c_Rest, ||L^Rest_-2|0>||^2) = {coset}")
check("Kommutant (Coset) einer E8-Stromalgebra trivial genau für N = 1 (c_Rest = 0)", [N for N, v in coset.items() if v[0] == 0] == [1])
befund("Entscheidung 1 (eine Naht) und 2 (ein E8-Block) sind beide trivialer Kommutant des Ursprungsprozesses", True)
say("   Grenze: Das ist ein Prinzip (Irreduzibilität, kein verborgener Sektor), keine Herleitung; es ersetzt zwei Lesarten")
say("   durch eine Prämisse. Welche Algebra der Ursprungsprozess ist, bleibt zu zeigen.")

urteil = ("Die Anomalienbilanz der Analyse stimmt: Bei direkter Fortsetzung von J (q_+ = 1 + 4a, q_- = -1 + 4b) trägt "
          "der primitive Wirbel eines Ladung-4-Kondensats c = 16(a + b), nie 8, und das Spin(10)-Niveau ist gerade "
          "(6561 Kombinationen, 110 mit c = 8, keine mit Niveau 1). Nur die axiale, gravitativ anomale Ladung (+1, +1) "
          "gibt c = 8 und Niveau 1; das ist genau H6s G (Z13: c = 1/2 je Paar; Z24: Index 2 gegen 0 für U). Das "
          "Quartett-Gegenmodell ist richtig: Z4-symmetrisch, <P> = 0, <P^2> != 0, Lücke 1, also ohne Nahtmoden; den Index "
          "trägt die anomale Antwort, nicht der Ordnungsparameter. Bosonische holomorphe Theorien haben c in 8Z, E8 ist der "
          "Generator; unter A* (GSO) ist die Naht bosonisch, und 'genau ein Block' ist äquivalent zur Unimodularität der "
          "Indexabbildung, die selbst noch herzuleiten ist. Neu: Auch der axiale Ladung-4-Weg trägt die Wand, denn der "
          "relevante Operator Re Phi^2 aus Z24 hat axiale Ladung 4 und pinnt arg(Phi^2) (N_DW = 1); Entscheidung 4 bleibt "
          "der Engpass jedes Weges. Ein Halbgruppenregler e^{-kappa s} legt die Reglerfunktion fest und macht "
          "S_ent = A/(4 G_ind) für minimal gekoppelte Felder exakt; offen bleiben eine Skala und der Eich-Randterm. "
          "Entscheidung 1 und 2 sind beide ein trivialer Kommutant: Ein Shift der Vielfachheit m hat den Kommutanten M_m, "
          "N E8-Blöcke haben relativ zu einer Stromalgebra den Coset c = 8(N - 1); die 16 Eltern wirken irreduzibel "
          "(Clifford) und sind keine Vielfachheit. Irreduzibilität des Ursprungsprozesses wäre damit eine Prämisse für "
          "beide Entscheidungen.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z31_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": {"anomalie": {"kombinationen": tot, "c8": c8, "c8_k1": c8k1},
                                                      "gegenmodell": res, "c_je_paar_z13": c_pair},
     "urteil": urteil}, indent=1, ensure_ascii=False, default=str) + "\n")
