"""Z1: Prüfung der externen Kritik an H6h2 (23. November 2026) mit Rechnungen.

Die Kritik nennt vier Punkte, die hier einzeln gerechnet werden.
 1. Modulare Kovarianz wählt kein Quadrat: Z_AA(-1/tau) = Z_AA(tau) gilt auch bei tau = 2i.
    Geprüft wird, ob Zeitneutralität (gleiches thermisches System in beiden Rekonstruktionen)
    mehr verlangt als Kovarianz, und die vorgeschlagene Verschärfung (echte Drehung der Ordnung
    größer zwei mit erhaltener Spinstruktur ohne Nullmoden).
 2. Spinorwirkung der Vierteldrehung: vier Vierteldrehungen geben auf Spinoren -1. Hält die
    Viertelzählung c3 = 2/(4 pi 4) = 1/(8 pi) mit der gelifteten Drehung?
 3. Familienvorzeichen: S4 = Aut(Q8) auf Im H gibt {1, i, -i}; {-1, i, -i} verlangt das Signum.
    Woher kommt es, ohne es einzusetzen?
 4. E8 und die T Phase exp(-2 pi i c/24): was verlangt Basisunabhängigkeit genau?
"""
import itertools
import json
from pathlib import Path
import numpy as np
from mpmath import mp, mpc, jtheta, exp, pi, qp

mp.dps = 30
np.set_printoptions(precision=4, suppress=True)
OK = []


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


# ---------------------------------------------------------------------------------------------
print("1. Kovarianz gegen Zeitneutralität\n")


def eta(tau):
    return exp(1j * pi * tau / 12) * qp(exp(2j * pi * tau))


def z_aa(tau):
    return abs(jtheta(3, 0, exp(1j * pi * tau)) / eta(tau)) ** 2


for tau in (mpc(0, 2), mpc(0, 1)):
    print(f"   tau = {mp.nstr(tau, 3)}: Z_AA(tau) = {mp.nstr(z_aa(tau), 12)}, Z_AA(-1/tau) = {mp.nstr(z_aa(-1 / tau), 12)}")
print(f"   Kovarianz gilt bei tau = 2i: {check('Z_AA kovariant bei 2i', abs(z_aa(mpc(0, 2)) - z_aa(mpc(0, 0.5))) < 1e-20)}"
      " (die Kritik hat hier recht)")
# the two Osterwalder-Schrader reconstructions of the same torus
for name, (space, beta) in (("A: Raum 1, Zeit 2", (1.0, 2.0)), ("B: Raum 2, Zeit 1", (2.0, 1.0))):
    lam1 = np.pi / space                                    # lowest antiperiodic level on the circle
    occ = 1 / (1 + np.exp(beta * lam1))
    print(f"   Rekonstruktion {name}: beta/L = {beta / space:.2f}, Besetzung der tiefsten Mode {occ:.5f}")
print("   Gleiche Zustandssumme, aber verschiedene thermische Systeme (kalt gegen heiß). Zeitneutral im"
      "\n   Sinn von H6 ist erst die Gleichheit der Systeme; das ist die Fixpunktbedingung, nicht die Kovarianz.")


def lattice_has(tau, unit):
    """Does multiplication by unit map the lattice Z + Z tau to itself?"""
    for w in (1, tau):
        v = unit * w
        # solve v = a + b tau with integers a, b
        b = v.imag / tau.imag
        a = v.real - b * tau.real
        if abs(a - round(a)) > 1e-9 or abs(b - round(b)) > 1e-9:
            return False
    return True


extra = []
for x in np.linspace(-0.5, 0.5, 101):
    for y in sorted(set(np.linspace(0.8, 2.0, 121)) | {np.sqrt(3) / 2}):
        tau = complex(x, y)
        if abs(tau) < 1 - 1e-12:
            continue
        for order, unit in ((4, 1j), (6, np.exp(1j * np.pi / 3))):
            if lattice_has(tau, unit):
                extra.append((round(x, 4), round(y, 4), order))
found = sorted({(float(x), float(y), o) for x, y, o in extra})
print(f"\n   Tori mit Drehung der Ordnung > 2 im Fundamentalbereich: {found}"
      "\n   (die beiden Randpunkte -1/2 und +1/2 sind derselbe Sechseckstorus)")
check("Drehung >2 nur bei i und rho", found == [(-0.5, 0.866, 6), (0.0, 1.0, 4), (0.5, 0.866, 6)])


def spin_action(basis_image):
    """Spin structures as holonomy signs h(e1), h(e2); action h -> h o g."""
    out = {}
    for h1, h2 in itertools.product((-1, 1), repeat=2):
        imgs = []
        for (c1, c2) in basis_image:                       # g(e_i) = c1 e1 + c2 e2
            imgs.append((h1 ** (c1 % 2)) * (h2 ** (c2 % 2)))
        out[(h1, h2)] = tuple(imgs)
    return out


name = {(-1, -1): "AA", (-1, 1): "AP", (1, -1): "PA", (1, 1): "PP"}
sq = spin_action([(0, 1), (-1, 0)])                         # i*1 = i, i*i = -1
hx = spin_action([(0, 1), (-1, 1)])                         # rho*1 = rho, rho*rho = rho - 1
print("   Quadrat, Vierteldrehung:  " + ", ".join(f"{name[k]}->{name[v]}" for k, v in sq.items()))
print("   Sechseck, Sechsteldrehung: " + ", ".join(f"{name[k]}->{name[v]}" for k, v in hx.items()))
fix_sq = [name[k] for k, v in sq.items() if k == v]
fix_hx = [name[k] for k, v in hx.items() if k == v]
print(f"   invariant: Quadrat {fix_sq}, Sechseck {fix_hx}; PP hat Nullmoden")
check("Quadrat hält AA", "AA" in fix_sq)
check("Sechseck hält nur PP", fix_hx == ["PP"])
print("   Die verschärfte Fassung der Kritik (echte Drehung, Ordnung > 2, Spinstruktur ohne Nullmoden"
      "\n   erhalten) wählt genau tau = i mit (A,A). Das ist dieselbe Bedingung wie in N und deckt sich mit"
      "\n   dem Atlasbefund vom 16.11. (nur j = 1728 erlaubt Ordnung vier auf den Marken).")

# ---------------------------------------------------------------------------------------------
print("\n2. Gelifte Vierteldrehung g(z) = i z + 1/2 auf Spinoren, Spinstruktur (A,A)\n")
marks = [0, 0.5, 0.5 + 0.5j, 0.5j]


def g(z):
    return 1j * z + 0.5


def mod1(z):
    return complex(round(z.real % 1.0, 9) % 1.0, round(z.imag % 1.0, 9) % 1.0)


img = [mod1(g(m)) for m in marks]
cyc = [marks.index(next(n for n in marks if abs(mod1(n) - w) < 1e-9)) for w in img]
print(f"   Marken 0, 1/2, (1+i)/2, i/2 gehen auf Indizes {cyc}: Viererzyklus "
      f"{check('g zyklisch auf Marken', cyc == [1, 2, 3, 0])}")
z0 = 0.123 + 0.456j
print(f"   g^4 = Identität auf dem Torus: {check('g^4 = id', abs(mod1(g(g(g(g(z0))))) - mod1(z0)) < 1e-12)}")

sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
S = np.diag([np.exp(-1j * np.pi / 4), np.exp(1j * np.pi / 4)])   # spinor rotation by pi/2
shift = np.array([0.5, 0.0])


def rot(k):
    return np.array([-k[1], k[0]])


def block(k, a):
    ks = [np.array(k, float)]
    for _ in range(3):
        ks.append(rot(ks[-1]))
    G = np.zeros((8, 8), complex)
    D = np.zeros((8, 8), complex)
    for j, kk in enumerate(ks):
        jj = (j + 1) % 4
        phase = np.exp(-2j * np.pi * np.dot(ks[jj], a))
        G[2 * jj:2 * jj + 2, 2 * j:2 * j + 2] = phase * S
        D[2 * j:2 * j + 2, 2 * j:2 * j + 2] = 2 * np.pi * (kk[0] * sx + kk[1] * sy)
    return G, D


G, D = block((0.5, 1.5), shift)
comm = np.max(np.abs(G @ D - D @ G))
g4 = np.linalg.matrix_power(G, 4)
print(f"   [g~, D] = {comm:.1e}: {check('g~ vertauscht mit D', comm < 1e-12)};  g~^4 = -1: "
      f"{check('g~^4 = -1', np.allclose(g4, -np.eye(8)))}")
ev = np.linalg.eigvals(G)
phases = np.round(np.angle(ev) / (np.pi / 4)) % 8
counts = {int(p): int(np.sum(phases == p)) for p in sorted(set(phases))}
print(f"   Eigenwerte von g~ auf einer Bahn (in Einheiten pi/4): {counts}")
check("je Sektor zwei Zustände je Bahn", counts == {1: 2, 3: 2, 5: 2, 7: 2})
naive = sum(np.linalg.matrix_power(G, k) for k in range(4)) / 4
print(f"   naiver Projektor (1/4) Summe g~^k, k=0..3, ist idempotent: "
      f"{'ja' if np.allclose(naive @ naive, naive) else 'nein'}; invarianter Sektor hat Dimension "
      f"{int(np.sum(np.abs(ev - 1) < 1e-9))}")
print("   Mit Spinoren gibt es keinen invarianten Sektor (g~^4 = -1). Die Sektoren sind die vier Eigenwerte"
      "\n   exp(i pi (2m+1)/4), und jeder trägt auf jeder Bahn genau 2 von 8 Zuständen.")

t = 0.002
n = np.arange(-90, 90) + 0.5
K1, K2 = np.meshgrid(n, n)
total = 2 * np.sum(np.exp(-t * 4 * np.pi ** 2 * (K1 ** 2 + K2 ** 2)))
per_sector = total / 4                                       # exact: every orbit has length 4
fixed = int(np.sum((np.abs(K1 + K2) < 1e-12) & (np.abs(K1 - K2) < 1e-12)))
print(f"   Impulse mit Fixpunkt unter der Drehung: {fixed}; also Wärmespur je Sektor = Spur/4 exakt.")
print(f"   t * Spur je Sektor bei t = {t}: {t * per_sector:.12f}, 1/(8 pi) = {1 / (8 * np.pi):.12f}")
check("c3 je Spinorsektor", abs(t * per_sector - 1 / (8 * np.pi)) < 1e-9)
for a, label in ((np.zeros(2), "Drehung um 0"), (shift, "Drehung um (1+i)/4")):
    Gb, _ = block((0.5, 0.5), a)
    e = np.round(np.angle(np.linalg.eigvals(Gb)) / (np.pi / 4)) % 8
    print(f"   {label}: Sektorzählung auf der Bahn {{{', '.join(f'{int(p)}: {int(np.sum(e == p))}' for p in sorted(set(e)))}}}")
print("   Ergebnis: Die fermionische Drehphase ändert die Viertelzählung nicht. Jeder der vier spinoriellen"
      "\n   Sektoren hat exakt ein Viertel jedes Eigenraums, ohne Fixpunktkorrektur, also c3 = 1/(8 pi).")

# ---------------------------------------------------------------------------------------------
print("\n3. Familienvorzeichen\n")
R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
C = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
group = [np.eye(3, dtype=int)]
frontier = [np.eye(3, dtype=int)]
while frontier:
    new = []
    for x in frontier:
        for y in (R, C):
            z = y @ x
            if not any(np.array_equal(z, w) for w in group):
                group.append(z)
                new.append(z)
    frontier = new
print(f"   <R, C> auf Im H: Ordnung {len(group)} {check('Ordnung 24', len(group) == 24)}; R hat "
      f"Eigenwerte {np.round(np.linalg.eigvals(R), 6)}")
diag = [np.array(d) for d in ((1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1))]


def perm_of(m):
    P = np.zeros((4, 4))
    for a, d in enumerate(diag):
        v = m @ d
        b = next(b for b, e in enumerate(diag) if np.array_equal(v, e) or np.array_equal(v, -e))
        P[b, a] = 1
    return P


ok_twist, ok_wedge = True, True
for m in group:
    P = perm_of(m)
    sgn = round(np.linalg.det(P))
    chi_std = np.trace(P) - 1
    chi_rot = np.trace(m)
    P2 = perm_of(m @ m)
    chi_std_sq = np.trace(P2) - 1
    ok_twist &= abs(chi_rot - sgn * chi_std) < 1e-12
    ok_wedge &= abs(chi_rot - (chi_std ** 2 - chi_std_sq) / 2) < 1e-12
print(f"   Charakter auf Im H = Signum x Standard für alle 24 Elemente: {check('rot = sgn x std', ok_twist)}")
print(f"   Charakter auf Im H = Lambda^2(Standard) für alle 24 Elemente: {check('rot = Lambda2 std', ok_wedge)}")
print("   Also: Die chiralen 2-Formen (Im H) sind die zweite äußere Potenz der Standarddarstellung,"
      "\n   weil Lambda^2 V = V* x det V und det(Standard) = Signum.")

# the seam itself: characters of the mark group E[2] (= the three even spin structures)
weil = {}


def e2(h, m):
    """Weil pairing on E[2] for the lattice Z[i]: e2(h, m) = (-1)^(2 det(h, m))."""
    d = 2 * (h.real * m.imag - h.imag * m.real)
    return -1 if round(4 * d) % 4 == 2 else 1


halfp = [0.5, 0.5j, 0.5 + 0.5j]
labels = {0.5: "Halbperiode 1/2", 0.5j: "Halbperiode i/2", 0.5 + 0.5j: "Halbperiode (1+i)/2"}
M = np.zeros((3, 3), complex)
for col, h in enumerate(halfp):
    ih = mod1(1j * h)
    row = next(r for r, hh in enumerate(halfp) if abs(mod1(hh) - ih) < 1e-9)
    sign = e2(ih, 0.5)                                       # translation part t = 1/2 of g
    M[row, col] = sign
spec = np.sort_complex(np.round(np.linalg.eigvals(M), 9))
print(f"\n   g = i z + 1/2 auf den drei nichttrivialen Charakteren der Markengruppe E[2] (Weilpaarung):")
print(f"   Matrix {M.real.astype(int).tolist()}, Eigenwerte {spec}")
check("Familiencharaktere aus g", np.allclose(np.sort_complex(spec), np.sort_complex(np.array([-1, -1j, 1j]))))
M0 = np.zeros((3, 3))
for col, h in enumerate(halfp):
    ih = mod1(1j * h)
    row = next(r for r, hh in enumerate(halfp) if abs(mod1(hh) - ih) < 1e-9)
    M0[row, col] = 1
print(f"   Kontrolle Drehung um 0 (ohne Verschiebung): Eigenwerte {np.round(np.linalg.eigvals(M0), 9)}")
print("   Die nichttrivialen Charaktere von E[2] sind die drei geraden Spinstrukturen; (1+i)/2 gehört zu (A,A)."
      "\n   Die zyklische Markendrehung ist affin. Ihr Verschiebungsteil wirkt über die Weilpaarung als"
      "\n   Vorzeichen, und das ergibt genau {-1, i, -i}: (A,A) mit -1, AP und PA mit +-i. Das Vorzeichen"
      "\n   ist damit nicht eingesetzt; es ist der Verschiebungsteil der Drehung um den Viertelpunkt.")
marks_perm = np.zeros((4, 4))
for a, b in enumerate(cyc):
    marks_perm[b, a] = 1
Q = np.eye(4) - np.ones((4, 4)) / 4
h1 = np.linalg.eigvals(Q @ marks_perm @ Q)
h1 = h1[np.abs(h1) > 1e-9]
print(f"   Dasselbe g auf H^1 der vier Marken (Permutation ohne Summe): {np.sort_complex(np.round(h1, 9))}")
check("H1 der Marken gleich", np.allclose(np.sort_complex(np.round(h1, 9)), np.sort_complex(spec)))
Gamma_sign = round(np.linalg.det(marks_perm))
print(f"   Orientierter Viermarkenzustand (Produkt der vier Majoranas) unter g: Faktor {Gamma_sign} = Signum")

# ---------------------------------------------------------------------------------------------
print("\n4. E8 und die T Phase\n")


def chi(N, tau):
    nome = exp(1j * pi * tau)
    e = eta(tau)
    return sum((jtheta(a, 0, nome) / e) ** (N // 2) for a in (2, 3, 4)) / 2


tau0 = mpc(0.13, 1.21)
print("   N Majorana-Weyl, Summe über die drei geraden Spinstrukturen:")
print("    N    |S Verhältnis - 1|   T Phase / (2 pi)    T kovariant")
for N in (8, 16, 24, 32, 48):
    s_ratio = chi(N, -1 / tau0) / chi(N, tau0)
    t_ratio = chi(N, tau0 + 1) / chi(N, tau0)
    cov = abs(abs(t_ratio) - 1) < 1e-15
    ph = float(mp.arg(t_ratio) / (2 * pi)) if cov else float("nan")
    print(f"   {N:3d}    {float(abs(s_ratio - 1)):.1e}            {ph:+.6f}           {'ja' if cov else 'nein'}")
e4 = lambda tau: 1 + 240 * sum((n ** 3) * exp(2j * pi * tau * n) / (1 - exp(2j * pi * tau * n)) for n in range(1, 40))
print(f"   N = 16 ist E4/eta^8: {check('chi16 = E4/eta8', abs(chi(16, tau0) - e4(tau0) / eta(tau0) ** 8) < 1e-20)}")
j_val = chi(16, mpc(0, 1)) ** 3
print(f"   chi_E8(i) = {mp.nstr(chi(16, mpc(0, 1)).real, 12)}, chi_E8(i)^3 = {mp.nstr(j_val.real, 12)} = j(i) = 1728")
print("   S ist für jedes N exakt; die Vierteldrehung braucht nur S. T (Basiswechsel) ist für N = 16 nur bis auf"
      "\n   exp(-2 pi i/3) erfüllt. Strenge Basisunabhängigkeit ohne Zusatzstruktur verlangt N = 48 (c = 24),"
      "\n   mit Phase c/24 ganzzahlig. Die zwei konsistenten Ergänzungen sind daher:"
      "\n   (i)  eine invertierbare 3D E8 Phase (Rand c = 8), deren T^2 Hilbertraum eindimensional ist, so dass"
      "\n        die Wahl der Füllung die Vierteldrehung nicht bricht; (ii) drei identische E8 Blöcke, c = 24,"
      "\n        Zustandssumme chi^3 = j.")
omega_bar = complex(exp(-2j * pi / 3))
print(f"   Die T Phase exp(-2 pi i/3) = {omega_bar.real:+.4f}{omega_bar.imag:+.4f}i hat Argument 240 Grad;"
      f"\n   TFPTs leptonische CP Phase ist delta_PMNS = 4 pi/3 = 240 Grad (GALOIS.CP.PREDICT.01). Gleiche Zahl,"
      f"\n   ein Zusammenhang ist nicht hergeleitet.")

print("\nZusammenfassung der Prüfungen:")
for name_, cond in OK:
    print(f"   {'bestanden' if cond else 'NICHT bestanden'}: {name_}")
assert all(c for _, c in OK)
Path(__file__).with_name("z1_ergebnis.json").write_text(json.dumps({
    "checks": {n: c for n, c in OK},
    "family_spectrum_from_weil": [str(complex(x)) for x in spec],
    "c3_per_spinor_sector": float(t * per_sector),
    "e8_T_phase_over_2pi": -1 / 3,
    "strict_basis_independence_N": 48,
}, ensure_ascii=False, indent=1))
