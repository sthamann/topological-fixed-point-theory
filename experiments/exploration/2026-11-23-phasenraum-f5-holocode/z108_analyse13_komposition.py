"""Z108: Prüfung der dreizehnten externen Analyse 'Fundamentale Lösung nach H6h47' (25.11.).

Teil 1  Der Zustand aus sechs Saaten (z94): Randverteilungen, Paarkorrelationen, Zerlegung des Modulexponenten in
        Mehrkörperterme (Möbius-Inversion auf dem Booleschen Verband), [K, n_i] = 0, Varianz bezüglich der normierten Spur.
Teil 2  Komposition: Erhält der gemeinsame Modulgenerator die vollständige Dynamik jeder Naht, so ist er additiv
        (Tensorprodukt und gemeinsame CAR), also ein Produktzustand. Nichtadditiver Teil I_AB|C und bedingte gegenseitige
        Information.
Teil 3  Kompensierte Rotation zweier z94-Kopien (m = 1 und m = 3, zwölf Moden): Suche nach (A, B, C, D), Erhaltung von
        Einzelbesetzungen, Ladungen und Uhr, |K_AB| = 2 pi |sin 2 theta|, Berry-Krümmung, freier Winkel; der Übergang
        sitzt ganz im nichtadditiven Teil.
Teil 4  Rekombination: q^Delta/(1-q)^4 -> kurze Leiter + Skalar mit Delta = 3; Norm |P^2|Delta>|^2 = 32 Delta (Delta - 1)
        im Verma-Modul (d = 4); die fehlenden Nachkommen sind die aus z99.
Teil 5  Ursprung: Z(t) = Z(1/t), strikte Konvexität, einziger stationärer Punkt t = 1 (beta = 4 pi, tau = i);
        Varianzschranke Var_0(H) >= (1/4) sum m_i^2, Gleichheit genau für lineares H.
"""
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import expm, logm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(108)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def lade(z):
    return json.loads((HERE / f"{z}_ergebnis.json").read_text())


def kom(a, b):
    return a @ b - b @ a


# ------------------------------------------------------------------ Teil 1
say("Teil 1: Der Zustand aus sechs Saaten (z94)")
z94 = lade("z94")
fz = z94["zahlen"]["aussenbehauptung_gleiche_saaten"]["f(n) für N = 6"]
C6 = list(itertools.product((0, 1), repeat=6))
F6 = {c: fz.get("".join(map(str, c)), 0) for c in C6}
H6 = {c: sum(c) + 4*F6[c] for c in C6}                     # dimensionsloser Modulexponent, K = pi H + const
multi = all(Counter(H6[c] for c in C6 if c[k]) == Counter(H6[c] + 1 for c in C6 if not c[k]) for k in range(6))
mp.mp.dps = 50
yM = mp.exp(-mp.pi)
Zg = mp.fsum(yM**H6[c] for c in C6)
P6 = {c: yM**H6[c]/Zg for c in C6}
rand = max(abs(mp.fsum(v for c, v in P6.items() if c[k]) - yM/(1 + yM)) for k in range(6))
paark = []
for i, j in itertools.combinations(range(6), 2):
    ni = mp.fsum(v for c, v in P6.items() if c[i])
    nj = mp.fsum(v for c, v in P6.items() if c[j])
    paark.append(mp.fsum(v for c, v in P6.items() if c[i] and c[j]) - ni*nj)


def moebius(Hd, confs):
    return {S: sum((-1)**(sum(S) - sum(T))*Hd[T] for T in confs if all(t <= x for t, x in zip(T, S))) for S in confs}


korr = moebius({c: 4*F6[c] for c in C6}, C6)
voll = moebius(H6, C6)
n_korr = dict(sorted(Counter(sum(S) for S, v in korr.items() if v != 0 and sum(S) > 0).items()))
n_voll = dict(sorted(Counter(sum(S) for S, v in voll.items() if v != 0 and sum(S) > 0).items()))
ZAHLEN["n6_mehrkoerper"] = {"Korrektur 4f": n_korr, "voller Exponent": n_voll}
check(f"z94-Zustand y^(|n| + 4 f(n)) (51 verschobene Konfigurationen): Randverteilungen exakt (Multimengen {multi}, Fehler "
      f"{float(rand):.0e} bei 50 Stellen), {sum(abs(v) > mp.mpf(10)**-45 for v in paark)}/15 Paarkorrelationen ungleich 0; "
      f"Möbius-Zerlegung des Exponenten: Korrektur 4f hat {n_korr} nichtverschwindende Terme je Körperzahl, mit dem freien "
      f"Generator {n_voll} (sechs Einkörperterme), kein Sechskörperterm",
      multi and rand < mp.mpf(10)**-40 and sum(abs(v) > mp.mpf(10)**-45 for v in paark) == 15
      and n_korr == {1: 5, 2: 10, 3: 18, 4: 14, 5: 6} and n_voll[1] == 6 and 6 not in n_voll)
K6 = np.diag([np.pi*H6[c] for c in C6])
nops = [np.diag([float(c[k]) for c in C6]) for k in range(6)]
offdiag = np.max(np.abs(K6 - np.diag(np.diag(K6))))
var0 = float(np.var([H6[c] for c in C6]))
var0_prod = float(np.var([sum(c) for c in C6]))
ZAHLEN["n6_varianz"] = {"Var_0 Produkt": var0_prod, "Var_0 korreliert": var0}
check(f"der Modulgenerator ist diagonal in den Besetzungen: [K, n_i] = 0 für alle sechs Moden (max {max(np.max(np.abs(kom(K6, n))) for n in nops):.0e}), "
      f"keine Matrixelemente zwischen verschiedenen Besetzungen ({offdiag:.0e}): der Zustand erzeugt oder vernichtet keine "
      f"Paare. Varianz des Exponenten bezüglich der normierten Spur: Produkt {var0_prod}, korreliert {var0} (Analyse: 1,5 und 38)",
      max(np.max(np.abs(kom(K6, n))) for n in nops) == 0 and offdiag == 0 and var0_prod == 1.5 and var0 == 38.0)
befund("E: Die Rechnung der Analyse stimmt. Der Zustand aus sechs Saaten ist nichtgaußsch (Mehrkörperterme bis zur fünften "
       "Ordnung), aber rein diagonal: Er korreliert Besetzungen, kann aber keinen Übergang von einem zu zwei Paaren tragen", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Komposition mit erhaltener Einzeldynamik ist additiv")
d1, d2 = 2, 3
I1, I2 = np.eye(d1), np.eye(d2)
gens = []
for a in range(d1):
    for b in range(d1):
        E_ = np.zeros((d1, d1))
        E_[a, b] = 1
        gens.append(np.kron(E_, I2))
for a in range(d2):
    for b in range(d2):
        E_ = np.zeros((d2, d2))
        E_[a, b] = 1
        gens.append(np.kron(I1, E_))
D = d1*d2
rows = [np.kron(g.T, np.eye(D)) - np.kron(np.eye(D), g) for g in gens]     # vec([X, g]) spaltenweise
null_tp = int(np.sum(np.linalg.svd(np.vstack(rows), compute_uv=False) < 1e-10))


def jw(n):
    I_, Zp, sm = np.eye(2), np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    out = []
    for j in range(n):
        m_ = [Zp]*j + [sm] + [I_]*(n - j - 1)
        op = m_[0]
        for M in m_[1:]:
            op = np.kron(op, M)
        out.append(op)
    return out


cs = jw(4)
rows = []
for c in cs:
    for g in (c, c.T):
        rows.append(np.kron(g.T, np.eye(16)) - np.kron(np.eye(16), g))
null_car = int(np.sum(np.linalg.svd(np.vstack(rows), compute_uv=False) < 1e-10))
par = np.diag([(-1)**bin(s).count("1") for s in range(16)]).astype(float)
X1 = sum(rng.normal()*(cs[i].T @ cs[j] + cs[j].T @ cs[i]) for i in (0, 1) for j in (0, 1))
X1 = X1 + 0.7*(cs[0].T @ cs[1].T + cs[1] @ cs[0])                  # gerade, hermitesch, nur Moden 0, 1
gerade = np.allclose(kom(X1, par), 0)
fremd = max(np.max(np.abs(kom(X1, g))) for c in cs[2:] for g in (c, c.T))
check(f"Satz (Analyse Abschnitt 2): Was mit allen Erzeugern aller Teilalgebren vertauscht, ist ein Skalar. Tensorprodukt "
      f"M_2 (x) M_3: Lösungsraum von [X, A (x) 1] = [X, 1 (x) B] = 0 hat Dimension {null_tp}; gemeinsame CAR (vier Moden, zwei "
      f"Nähte, Jordan-Wigner): Dimension {null_car}. Ein gerades K_1 der ersten Naht (Parität {gerade}) vertauscht mit allen "
      f"Operatoren der zweiten (max {fremd:.0e}), also gilt [K_ges, A_i] = [K_i, A_i] genau für K_ges = sum K_i + c",
      null_tp == 1 and null_car == 1 and gerade and fremd < 1e-12)


def zufallszustand(d, rang=None):
    G = rng.normal(size=(d, rang or d)) + 1j*rng.normal(size=(d, rang or d))
    r = G @ G.conj().T
    return r/np.trace(r)


def ptrace(rho, dims, keep):
    n = len(dims)
    r = rho.reshape(dims + dims)
    for k in sorted(set(range(n)) - set(keep), reverse=True):
        r = np.trace(r, axis1=k, axis2=k + r.ndim//2)
    dk = int(np.prod([dims[k] for k in keep]))
    return r.reshape(dk, dk)


def Kmat(r):
    return -logm(r)


rk = zufallszustand(4)
KAB, KA = Kmat(rk), Kmat(ptrace(rk, [2, 2], [0]))
Atest = rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2))
abw_korr = np.max(np.abs(kom(KAB, np.kron(Atest, np.eye(2))) - np.kron(kom(KA, Atest), np.eye(2))))
rp = np.kron(zufallszustand(2), zufallszustand(2))
KABp, KAp = Kmat(rp), Kmat(ptrace(rp, [2, 2], [0]))
abw_prod = np.max(np.abs(kom(KABp, np.kron(Atest, np.eye(2))) - np.kron(kom(KAp, Atest), np.eye(2))))
r3 = zufallszustand(8)
S = lambda r: float(np.real(-np.trace(r @ logm(r))))
erw = lambda K, r: float(np.real(np.trace(r @ K)))
dims = [2, 2, 2]
rAC, rBC, rC = ptrace(r3, dims, [0, 2]), ptrace(r3, dims, [1, 2]), ptrace(r3, dims, [2])
KAC = np.einsum("acdf,be->abcdef", Kmat(rAC).reshape(2, 2, 2, 2), np.eye(2)).reshape(8, 8)    # K_AC (x) 1_B
KBC = np.kron(np.eye(2), Kmat(rBC))
KC = np.kron(np.eye(4), Kmat(rC))
Iop = KAC + KBC - KC - Kmat(r3)
cmi_op = erw(Iop, r3)
cmi_S = S(rAC) + S(rBC) - S(rC) - S(r3)
check(f"Modulflüsse: Für einen korrelierten Zustand weicht [K_AB, A (x) 1] von [K_A, A] (x) 1 ab ({abw_korr:.2f}), für einen "
      f"Produktzustand nicht ({abw_prod:.0e}): Erhaltene Einzeldynamik heißt Produktzustand. Nichtadditiver Teil "
      f"I_AB|C = K_AC + K_BC - K_C - K_ABC: <I> = {cmi_op:.6f} = S_AC + S_BC - S_C - S_ABC = {cmi_S:.6f} >= 0 (bedingte "
      "gegenseitige Information, starke Subadditivität)", abw_korr > 1e-3 and abw_prod < 1e-10 and abs(cmi_op - cmi_S) < 1e-10
      and cmi_S > 0)
befund("E: Der Satz der Analyse stimmt und ist die endliche Form eines bekannten Sachverhalts: Die Modulgruppe des Ganzen "
       "erhält die Modulgruppen der Teile genau im Produktzustand. Eine Komposition, die jede Naht als autonomes System "
       "unverändert lässt, kann keine Wechselwirkung erzeugen; die Wechselwirkung muss im nichtadditiven Teil liegen", True)

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Kompensierte Rotation zweier z94-Kopien (m = 1 und m = 3)")
qplus = lambda c: sum(c[2*j] - c[2*j + 1] for j in range(3))       # Paare (0, 1), (2, 3), (4, 5): Teilchen - Loch
info, grp = {}, defaultdict(list)
for a in C6:
    for b in C6:
        n = a + b
        E = H6[a] + 3*H6[b]                                        # Kopie m = 3: Exponent 3(|n| + 4f)
        dat = (E, sum(n), sum(a) - sum(b), qplus(a) + qplus(b))  # Q_- = sum eps nu (eps = +1 für m = 1, -1 für m = 3)
        info[n] = dat
        grp[(dat[1], dat[0], dat[2], dat[3])].append(n)
loes = []
for (exc, E, qm, qp), As in grp.items():
    if exc != 2:
        continue
    for A in As:
        for B in grp.get((4, E + 4, qm, qp), []):
            delta = tuple(b - a for a, b in zip(A, B))
            for (exc2, E2, qm2, qp2), Cs in grp.items():
                if exc2 != 6 or E2 != E:
                    continue
                for Cc in Cs:
                    Dd = tuple(c - dl for c, dl in zip(Cc, delta))
                    if all(x in (0, 1) for x in Dd) and info[Dd] == (E + 4, 4, qm2, qp2):
                        loes.append((A, B, Cc, Dd))
tab = [l_ for l_ in loes if info[l_[0]][0] == 18 and info[l_[0]][2:] == (2, 0) and info[l_[2]][2:] == (-2, 0)]
ZAHLEN["rotation_loesungen"] = {"Anzahl (A, B, C, D)": len(loes), "wie Tabelle der Analyse": len(tab)}
A, B, Cc, Dd = tab[0]
fmt = lambda n: "".join(map(str, n[:6])) + "|" + "".join(map(str, n[6:]))
ZAHLEN["rotation_tabelle"] = {k: {"Konfiguration": fmt(n), "Anregungen": info[n][1], "Exponent": info[n][0],
                                  "(Q_-, Q_+)": info[n][2:]} for k, n in zip("ABCD", (A, B, Cc, Dd))}
pAll = {n: mp.exp(-mp.pi*info[n][0]) for n in info}
Zt = mp.fsum(pAll.values())
pAll = {n: v/Zt for n, v in pAll.items()}


def gedreht(theta):
    """Diagonale nach der Drehung U (A <-> B und C <-> D um theta); Kohärenzen zwischen A und B bzw. C und D."""
    q = dict(pAll)
    c2, s2 = mp.cos(theta)**2, mp.sin(theta)**2
    for X, Y in ((A, B), (Cc, Dd)):
        q[X], q[Y] = c2*pAll[X] + s2*pAll[Y], s2*pAll[X] + c2*pAll[Y]
    return q


ok_rot, kab = True, {}
for th in (mp.pi/6, mp.mpf("0.37"), mp.pi/4, mp.mpf("1.2")):
    q = gedreht(th)
    for k in range(12):
        d_ = mp.fsum(q[n] - pAll[n] for n in info if n[k])
        ok_rot &= abs(d_) < mp.mpf(10)**-45
    for lad in (2, 3):
        for wert in set(info[n][lad] for n in info):
            ok_rot &= abs(mp.fsum(q[n] - pAll[n] for n in info if info[n][lad] == wert)) < mp.mpf(10)**-45
    kA, kB = mp.pi*info[A][0], mp.pi*info[B][0]
    U = mp.matrix([[mp.cos(th), -mp.sin(th)], [mp.sin(th), mp.cos(th)]])
    Kp = U*mp.diag([kA, kB])*U.T
    kab[float(th)] = float(abs(Kp[0, 1]))
    ok_rot &= abs(abs(Kp[0, 1]) - 2*mp.pi*abs(mp.sin(2*th))) < mp.mpf(10)**-40
uhr_ok = all((info[X][0] - info[X][2]) % 4 == 0 for X in (A, B, Cc, Dd)) and info[A][2] == info[B][2] and info[Cc][2] == info[Dd][2]
paritaet = all(info[X][1] % 2 == 0 for X in (A, B, Cc, Dd))
ZAHLEN["rotation_KAB"] = kab
check(f"Suche in zwölf Moden (Kopie m = 1 und Kopie m = 3 des z94-Zustands, Q_- = sum eps nu, Q_+ = Teilchen - Loch): "
      f"{len(loes)} kompensierte Vierergruppen, davon {len(tab)} genau wie die Tabelle der Analyse: A = {fmt(A)} (2 Anregungen, "
      f"Exponent {info[A][0]}, {info[A][2:]}), B = {fmt(B)} (4, {info[B][0]}), C = {fmt(Cc)} (6, {info[Cc][0]}, "
      f"{info[Cc][2:]}), D = {fmt(Dd)} (4, {info[Dd][0]}). Für vier Winkel bleiben alle zwölf Einzelbesetzungen und beide "
      f"Ladungsverteilungen exakt (50 Stellen), Uhr und Parität erhalten ({uhr_ok}, {paritaet}); |K_AB| = 2 pi |sin 2 theta|, "
      f"bei 30 Grad {kab[float(mp.pi/6)]:.4f} (Analyse: 5,4414)",
      len(tab) == 1 and ok_rot and uhr_ok and paritaet and abs(kab[float(mp.pi/6)] - 5.4414) < 1e-4)
th_, ph_ = sp.symbols("theta phi", real=True)
u = sp.Matrix([sp.cos(th_), sp.exp(sp.I*ph_)*sp.sin(th_)])
Ath = sp.I*(u.H*u.diff(th_))[0]
Aph = sp.I*(u.H*u.diff(ph_))[0]
Fberry = sp.simplify(sp.diff(Aph, th_) - sp.diff(Ath, ph_))
nA = [n1 for n1 in (A[:6], A[6:])]
check(f"Berry-Krümmung von |u> = cos theta |A> + e^(i phi) sin theta |B>: F = {Fberry} d theta ^ d phi (Analyse: "
      f"-sin(2 theta)); |K_AB| und F hängen an derselben Deformation sin 2 theta. Der Winkel ist frei: Jeder Winkel erfüllt "
      "alle Bedingungen (vier geprüft), also eine stetige Mehrdeutigkeit",
      sp.simplify(Fberry + sp.sin(2*th_)) == 0 and ok_rot)
zweifach = A[:6] != B[:6] and A[6:] != B[6:] and Cc[:6] != Dd[:6] and Cc[6:] != Dd[6:]
check(f"Der Übergang A <-> B ändert beide Kopien zugleich (Kopie 1: {fmt(A)[:6]} -> {fmt(B)[:6]}, Kopie 2: {fmt(A)[7:]} -> "
      f"{fmt(B)[7:]}): Die lokalen Modulgeneratoren K_1 (x) 1 und 1 (x) K_2 haben zwischen A und B kein Matrixelement, also "
      "liegt das ganze Element <A|K|B> im nichtadditiven Teil K_1 + K_2 - K_12 (gegenseitige Information): Der Übergang "
      "ist eine Eigenschaft des gemeinsamen Zustands, nicht der Teile", zweifach)
befund("E: Die Konstruktion der Analyse stimmt exakt. Unter Einzelbesetzungen, Ladungen und Uhr gibt es nichtdiagonale "
       "gemeinsame Zustände mit einem Übergang von einem zu zwei Paaren, und dieselbe Deformation gibt Berry-Krümmung. Der "
       "Winkel ist frei (stetige Mehrdeutigkeit), und die Konstruktion ist eine Zustandsdeformation zweier Kopien, keine "
       "kovariante 4D-Wechselwirkung", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: Rekombination und Normbedingung")
q, Dl = sp.symbols("q Delta")
lang1 = q/(1 - q)**4
kurz = q*(1 - q**2)/(1 - q)**4
d3 = q**3/(1 - q)**4
reihe = lambda e, n=10: [sp.series(e, q, 0, n).removeO().coeff(q, k) for k in range(1, n)]
check(f"Charaktere: q^Delta/(1-q)^4 bei Delta = 1 ist kurz + Skalar(Delta = 3) exakt; kurze Leiter {reihe(kurz)[:6]} "
      f"(z99: (k+1)^2), fehlende Nachkommen {reihe(d3)[:6]} (z99: [0, 0, 1, 4, 10, 20])",
      sp.simplify(lang1 - kurz - d3) == 0 and reihe(kurz)[:6] == [1, 4, 9, 16, 25, 36]
      and reihe(d3)[:6] == [0, 0, 1, 4, 10, 20])


def verma_norm_P2(dim):
    """<Delta| (sum_nu K_nu K_nu) (sum_mu P_mu P_mu) |Delta> im Verma-Modul eines Skalars, euklidisch:
    [D, P] = P, [K_mu, P_nu] = 2 delta D - 2 M_mu nu, [M_ab, P_c] = delta_bc P_a - delta_ac P_b, K|D> = M|D> = 0."""
    def P(mu, st):
        return {(mu,) + w: c for w, c in st.items()}

    def Dop(st):
        return {w: c*(Dl + len(w)) for w, c in st.items()}

    def M(a, b, st):
        out = defaultdict(lambda: 0)
        for w, c in st.items():
            for j, x in enumerate(w):
                if x == b:
                    out[w[:j] + (a,) + w[j + 1:]] += c
                if x == a:
                    out[w[:j] + (b,) + w[j + 1:]] -= c
        return dict(out)

    def K(nu, st):
        out = defaultdict(lambda: 0)
        for w, c in st.items():
            for j, x in enumerate(w):
                suf = {w[j + 1:]: c}
                t = defaultdict(lambda: 0)
                if x == nu:
                    for ww, cc in Dop(suf).items():
                        t[ww] += 2*cc
                for ww, cc in M(nu, x, suf).items():
                    t[ww] -= 2*cc
                for ww, cc in t.items():
                    out[w[:j] + ww] += cc
        return dict(out)

    st = defaultdict(lambda: 0)
    for mu in range(dim):
        for w, c in P(mu, P(mu, {(): 1})).items():
            st[w] += c
    tot = 0
    for nu in range(dim):
        tot += K(nu, K(nu, dict(st))).get((), 0)
    return sp.factor(sp.expand(tot))


n4 = verma_norm_P2(4)
n3 = verma_norm_P2(3)
check(f"Norm im Verma-Modul (euklidische konforme Algebra, nur Kommutatoren): |P^2|Delta>|^2 = {n4} in d = 4 (Analyse: "
      f"32 Delta (Delta - 1)), allgemein 4 d Delta (2 Delta + 2 - d) (d = 3: {n3}); Nullzustand genau bei Delta = 1 (freies "
      "Feld, Box phi = 0). Aus Box phi = lambda O mit normiertem O folgt |lambda|^2 = 32 Delta (Delta - 1), für kleine anomale "
      "Dimension gamma = Delta - 1 also gamma = |lambda|^2/32",
      sp.simplify(n4 - 32*Dl*(Dl - 1)) == 0 and sp.simplify(n3 - 4*3*Dl*(2*Dl - 1)) == 0)
befund("E/C: Der Rekombinationsweg der Analyse ist richtig und präziser als die Deformation aus z99: Eine wechselwirkende "
       "Skalarleiter entsteht, wenn die kurze Leiter (Delta = 1) mit einem unabhängigen Skalar der Dimension 3 rekombiniert "
       "(Box phi = lambda O). Die Kopplung ist dann durch Delta festgelegt, |lambda|^2 = 32 Delta (Delta - 1). Nativ fehlt der "
       "Operator O und seine Dimension", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Ursprung: Selbstdualität, Konvexität, Varianzschranke")
mp.mp.dps = 40


def logZ(t, M=401):
    return mp.pi*t/12 + 2*mp.fsum(mp.log(1 + mp.exp(-mp.pi*t*m)) for m in range(1, M, 2))


dual = max(abs(logZ(t) - logZ(1/t)) for t in (mp.mpf("0.5"), mp.mpf("0.8"), mp.mpf("1.25"), mp.mpf(2), mp.mpf(3)))
d1Z = lambda t: mp.diff(logZ, t)
d2Z = lambda t: mp.diff(logZ, t, 2)
konvex = all(d2Z(mp.mpf(t)) > 0 for t in ("0.3", "0.7", "1", "1.6", "3"))
st1 = abs(d1Z(mp.mpf(1)))
vz = [mp.sign(d1Z(mp.mpf(t))) for t in ("0.3", "0.6", "0.9", "1.1", "1.5", "3")]
beta_scan = z94["zahlen"]["original_beta_scan"]["beta = 4 pi"]
ZAHLEN["ursprung"] = {"max |log Z(t) - log Z(1/t)|": float(dual), "d log Z/dt bei t = 1": float(st1),
                      "Vorzeichen d log Z/dt": [int(v) for v in vz]}
check(f"Z(t) = e^(pi t/12) prod_(m ungerade) (1 + e^(-pi t m))^2 ist das Quadrat der Weberfunktion bei tau = i t und "
      f"erfüllt Z(t) = Z(1/t) (Abweichung {float(dual):.0e}); log Z ist streng konvex (d^2 > 0 an fünf Punkten, analytisch "
      f"als Energievarianz), also hat d log Z/dt genau eine Nullstelle, und die liegt wegen der Dualität bei t = 1 "
      f"({float(st1):.0e}; Vorzeichen {[int(v) for v in vz]}): verschwindende Energieantwort genau bei beta = 4 pi, tau = i. "
      f"Das ist P_tn aus z94 (beta = 4 pi: tau_Gesetz = {beta_scan['tau_Gesetz']}, zeitneutral {beta_scan['zeitneutral']})",
      dual < mp.mpf(10)**-30 and konvex and st1 < mp.mpf(10)**-25 and vz == [-1, -1, -1, 1, 1, 1]
      and beta_scan["zeitneutral"] is True)


def walsh(Hd, confs, N):
    return {S: float(np.mean([Hd[c]*np.prod([2*c[i] - 1 for i in range(N) if S[i]]) for c in confs])) for S in confs}


W6 = walsh(H6, C6, 6)
lin = [W6[tuple(1 if j == i else 0 for j in range(6))] for i in range(6)]
exzess = sum(v**2 for S, v in W6.items() if sum(S) >= 2)
ok_satz = True
for _ in range(20):
    N_ = 5
    cf = list(itertools.product((0, 1), repeat=N_))
    m_ = rng.integers(1, 6, size=N_)
    Hz = {c: float(np.dot(m_, c)) + rng.normal()*float(sum(c) >= 2)*float(np.prod(c[:2])) for c in cf}
    Wz = walsh(Hz, cf, N_)
    lin_z = {i: Wz[tuple(1 if j == i else 0 for j in range(N_))] for i in range(N_)}
    korrigiert = {c: Hz[c] + sum((m_[i]/2 - lin_z[i])*(2*c[i] - 1) for i in range(N_)) for c in cf}
    v_ = float(np.var(list(korrigiert.values())))
    ok_satz &= v_ >= float(np.sum(m_**2))/4 - 1e-9
    lin_only = {c: float(np.dot(m_, c)) for c in cf}
    ok_satz &= abs(float(np.var(list(lin_only.values()))) - float(np.sum(m_**2))/4) < 1e-9
ZAHLEN["varianz_walsh_n6"] = {"lineare Walsh-Koeffizienten": lin, "Exzess (|S| >= 2)": exzess,
                              "Var_0 = lin + Exzess": sum(v**2 for v in lin) + exzess}
check(f"Varianzschranke (Walsh-Zerlegung bezüglich der normierten Spur): Stimmen die Einzelbesetzungen für alle "
      f"Temperaturen, so ist der lineare Teil festgelegt (h_i = m_i/2, aus d<n_i>/ds bei s = 0), also Var_0(H) = (1/4) sum "
      f"m_i^2 + Exzess >= (1/4) sum m_i^2, Gleichheit genau für lineares H (20 Zufallsfälle). Der z94-Zustand stimmt nur bei "
      f"einer Temperatur: lineare Koeffizienten {[round(x, 3) for x in lin]}, Exzess {exzess:.1f}, zusammen "
      f"{sum(v**2 for v in lin) + exzess:.1f} = Var_0",
      ok_satz and abs(sum(v**2 for v in lin) + exzess - var0) < 1e-9)
befund("C: Der Ursprungsweg der Analyse stimmt als bedingter Satz: Verlangt A* die verschwindende Energieantwort der "
       "Zustandssumme des Gesetzes, so folgt beta = 4 pi (tau = i, P_tn) eindeutig aus Selbstdualität und Konvexität. Offen "
       "bleibt, dass die Stressbedingung der E8-Naht genau diese Zustandssumme meint. Die Varianzschranke ersetzt S(n) durch "
       "eine einzige Bedingung, minimale modulare Varianz; ob TFPT sie erzwingt, ist offen", True)

urteil = ("Die dreizehnte externe Analyse hält in allen gerechneten Punkten. (1) Der Zustand aus sechs Saaten ist "
          "nichtgaußsch, mit Mehrkörpertermen 5, 10, 18, 14, 6, 0 (Korrektur) und sechs Einkörpertermen mit dem freien "
          "Generator, aber diagonal: [K, n_i] = 0, kein Übergang von einem zu zwei Paaren; Varianz 38 gegen 1,5. (2) Eine "
          "Komposition, die jede Naht als autonomes System erhält, ist additiv (Tensorprodukt und gemeinsame CAR), also ein "
          "Produktzustand; eine Wechselwirkung muss im nichtadditiven Teil I_AB|C liegen, dessen Erwartungswert die bedingte "
          "gegenseitige Information ist. (3) Die kompensierte Rotation existiert genau wie angegeben (eine von 20 "
          "Vierergruppen in zwölf Moden): Einzelbesetzungen, Ladungen, Uhr und Parität bleiben exakt, |K_AB| = 2 pi |sin 2 "
          "theta| (5,4414 bei 30 Grad), F = -sin 2 theta; der Winkel ist frei, und der Übergang liegt ganz im nichtadditiven "
          "Teil. (4) Rekombination: q/(1-q)^4 = kurz + Skalar mit Delta = 3, und |P^2|Delta>|^2 = 32 Delta (Delta - 1) im "
          "Verma-Modul; das präzisiert z99: Gebraucht wird ein unabhängiger Operator der Dimension 3, dann ist die Kopplung "
          "durch Delta festgelegt. (5) Z(t) = Z(1/t) und strenge Konvexität geben beta = 4 pi als einzigen Punkt "
          "verschwindender Energieantwort (bedingt durch A*), und die Varianzschranke reduziert die Starrheit auf minimale "
          "modulare Varianz. Offen bleiben die Festlegung des Winkels bzw. der Kopplung, der Operator O, die Gleichsetzung "
          "der Stressbedingung mit dem Gesetz, Chiralität, G und positives Lambda.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z108_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
