"""Z139: Folgt die Ordnung der Nähte (P_ord) aus halbseitigen modularen Inklusionen?

Frage: Die Wechselwirkung des Gesamtkandidaten ist das lokale kohärente Transfergatter zwischen Nähten (z120, z121). Es braucht
die Prämisse P_ord: unterscheidbare, in einer Kette geordnete Nähte; damit ist die Kette das Destri-de-Vega-Lichtkegelgitter.
Eine halbseitige modulare Inklusion N in M (Wiesbrock) liefert kanonisch eine Translation mit positivem Erzeuger und die
geordnete Kette N_n = U(n) M U(-n). Ist P_ord eine Folge der Modultheorie?

Leitlinie 'Rekonstruktion statt Auswahl'. Firewalls: GRAV.NONCIRCULAR.01 (Richtung Zustand -> Modulgenerator -> Inklusion ->
Translation; Hislop-Longo wird aus dem Zustand nachgerechnet, nicht vorausgesetzt), DIMENSION.UPLIFT.FIREWALL.01 (die Naht ist
eindimensional; die su(2,2)-Aussagen sind Lie-Algebra-Kinematik unter z83s Identifikation, keine Dynamik). Gelöste Blöcke
(Flavour, Massen, alpha, 16, E8) bleiben unberührt, H6 wird nicht geöffnet, nichts wird befördert, keine vollständige Theorie.
Statuskürzel: E exakt, N numerisch, C bedingt (Prämisse genannt), L Lesart, O offen.

Konventionen (wie z72)
  Naht S^1, Moden e_x = e^{i x theta}/sqrt(2 pi), x in Z + 1/2 (NS-Eichung; alpha = 1/4, 3/4 ist eichäquivalent, z72).
  G_v = -i (v d_theta + v'/2); L_0 = G_1, H_1 = G_cos, H_2 = G_sin; [G_v, G_w] = -i G_(v w' - w v').
  Modular: Delta^{it} = e^{-itK}; Vakuum: K_A = 2 pi G_(v_A), v_A >= 0 das Möbiusfeld mit Nullstellen in den Endpunkten.
  Marken z^4 = 1 bei theta = k pi/2 (z83/z104). Native Gebiete (z104): Kappe O_1 = (-pi/2, pi/2), K = 2 pi H_1 = 2 pi Q^-;
  Lichtkegel V_+ = (0, pi), K = 2 pi H_2; N_+ = O_1 cap V_+ = (0, pi/2), N_- = (-pi/2, 0). Cayley y = tan(theta/2) = t.
  Halbseitig: Delta_M^{it} N Delta_M^{-it} in N für t >= 0 ('+') oder t <= 0 ('-'); P = (K_M - K_N)/(2 pi).
  Reiner Zustand am Gitter: K_A = sum_k h_k (|v_k><v_k| + |w_k><w_k|), w_k = G_(A^c A) v_k / sqrt(nu_k (1 - nu_k)); die
  unverschränkten Gittermoden des Komplements (Artefakt, im Kontinuum nach F. und M. Riesz leer) bleiben ausgespart.

Teil 1  Die native halbseitige modulare Inklusion: exakt (Möbius), aus dem Zustand (Gitter, python-flint 640 Bit), Positivität,
        Negativkontrolle (geschachtelte Bögen ohne gemeinsamen Endpunkt), Originalzustand; beide Vorzeichen und ihr
        Kommutator; zwei Arten (P_chir); su(2,2) unter z83.
Teil 2  Die Kette N_n und ihre Plätze: Inhalt (Bögen einer Naht, nicht Nähte), Typ, Kovarianz, Schrittweite, Marken, Z4.
Teil 3  Lebt das Gatter auf der Kette? Lichtkegelstruktur der z121-Kette, freie Translationen, Invarianz, Marke.
Teil 4  Urteil: Einstufung von P_ord, Raumrichtung, 3+1 und Kontinuum.
"""
import json
from pathlib import Path

import numpy as np
import sympy as sp
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(139)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def r(v, n=4):
    return float(f"{float(v):.{n}g}")


def nrm(A):
    return float(np.linalg.norm(A))


def komm(A, B):
    return A @ B - B @ A


def lie_huelle(gens, tol=1e-9, runden=6):
    """Basis der von gens erzeugten reellen Lie-Algebra (Matrizen als reelle Vektoren, Gram-Schmidt in fester Reihenfolge)."""
    def vec(X):
        return np.concatenate([np.real(X).ravel(), np.imag(X).ravel()])
    basis, ortho = [], []

    def aufnehmen(X):
        v = vec(X)
        for o in ortho:
            v = v - np.dot(o, v)*o
        n_ = np.linalg.norm(v)
        if n_ > tol*max(1.0, np.linalg.norm(vec(X))):
            ortho.append(v/n_)
            basis.append(X)
            return True
        return False
    for g in gens:
        aufnehmen(np.asarray(g, complex))
    for _ in range(runden):
        neu = False
        k = len(basis)
        for i in range(k):
            for j in range(i + 1, k):
                neu |= aufnehmen(komm(basis[i], basis[j]))
        if not neu:
            break
    return basis


def killing_rang(basis, tol=1e-8):
    """Rang der Killing-Form der von basis aufgespannten Lie-Algebra (Strukturkonstanten per Kleinste Quadrate)."""
    n = len(basis)
    A = np.array([np.concatenate([np.real(X).ravel(), np.imag(X).ravel()]) for X in basis]).T
    f = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            C = komm(basis[a], basis[b])
            f[a, b] = np.linalg.lstsq(A, np.concatenate([np.real(C).ravel(), np.imag(C).ravel()]), rcond=None)[0]
    ad = [f[a].T for a in range(n)]
    kf = np.array([[np.trace(ad[a] @ ad[b]) for b in range(n)] for a in range(n)])
    return int(np.linalg.matrix_rank(kf, tol=tol)), kf


# =================================================================== Teil 1
say("Teil 1: Die native halbseitige modulare Inklusion")
say("  (a) exakt: Möbiusfelder, Wiesbrock-Differenzen, einseitige Inklusion")
TH = sp.symbols("theta", real=True)


def v_bogen(a, b):
    return 2*sp.sin((TH - a)/2)*sp.sin((b - TH)/2)/sp.sin((b - a)/2)


def null(expr):
    return sp.simplify(sp.expand_trig(sp.expand(expr)).rewrite(sp.exp)) == 0


def klammer(v, w):
    """[G_v, G_w] = -i G_u mit u = v w' - w v'."""
    return v*sp.diff(w, TH) - w*sp.diff(v, TH)


vO, vV = sp.cos(TH), sp.sin(TH)
vN, vNm = sp.cos(TH) + sp.sin(TH) - 1, sp.cos(TH) - sp.sin(TH) - 1
felder_ok = (null(v_bogen(-sp.pi/2, sp.pi/2) - vO) and null(v_bogen(0, sp.pi) - vV) and null(v_bogen(0, sp.pi/2) - vN)
             and null(v_bogen(-sp.pi/2, 0) - vNm))
dP, dM, dM2 = vO - vN, vV - vN, vO - vNm                     # 1 - sin, 1 - cos, 1 + sin
diff_ok = null(dP - (1 - sp.cos(TH - sp.pi/2))) and null(dM - (1 - sp.cos(TH))) and null(dM2 - (1 - sp.cos(TH + sp.pi/2)))
fac = sp.cot(sp.pi/4) - sp.cot(sp.pi/2)                      # z72: cot(delta_I/2) - cot(delta_J/2)
uhr_ok = null(vO + vV - vN - 1)
lie_ok = (null(klammer(vO, dP) + dP) and null(klammer(vV, dM) - dM) and null(klammer(vO, dM2) - dM2)
          and null(klammer(dP, dM) - vN) and null(klammer(dP, dM2) - 2*vO))
v_c = v_bogen(-sp.pi/4, sp.pi/4)                             # Kontrolle: streng geschachtelt
d_c = sp.simplify(sp.expand_trig(vO - v_c))
kontr_feld = null(d_c - (1 - (sp.sqrt(2) - 1)*sp.cos(TH)))
ell_disk = sp.simplify(1 - (sp.sqrt(2) - 1)**2)              # a^2 - b^2 - c^2 > 0: elliptisch
kontr_lie = [null(klammer(vO, d_c) - e_*d_c) for e_ in (1, -1)]

# Möbiusfluss der Kappe in y = tan(theta/2): y -> (y ch + sh)/(y sh + ch), s = 2 pi t, exakt rational
y_ = sp.symbols("y", real=True)
s_ = sp.symbols("s", real=True)
fluss = (y_*sp.cosh(s_/2) + sp.sinh(s_/2))/(y_*sp.sinh(s_/2) + sp.cosh(s_/2))
feld_y = sp.simplify(sp.diff(fluss, s_).subs(s_, 0) - (1 - y_**2)/2) == 0       # cos(theta) d_theta = (1 - y^2)/2 d_y


def moeb(e, y):
    """Fluss mit e^(s/2) = e (rational): ch = (e + 1/e)/2, sh = (e - 1/e)/2."""
    ch, sh = (e + 1/e)/2, (e - 1/e)/2
    return (y*ch + sh)/(y*sh + ch)


def drin(bild, ziel):
    return ziel[0] <= bild[0] and bild[1] <= ziel[1]


einseitig = {"N+": [], "N-": [], "Kontrolle": []}
BOEGEN = {"N+": (sp.Integer(0), sp.Integer(1)), "N-": (sp.Integer(-1), sp.Integer(0)),
          "Kontrolle": (sp.Rational(-1, 2), sp.Rational(1, 2))}
for e in (sp.Integer(2), sp.Integer(3), sp.Rational(7, 5), sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(5, 7)):
    for nm in ("N+", "N-", "Kontrolle"):
        a_, b_ = BOEGEN[nm]
        einseitig[nm].append((bool(e > 1), bool(drin((moeb(e, a_), moeb(e, b_)), (a_, b_)))))
hs_plus = all(inn == vor for vor, inn in einseitig["N+"])            # genau für s >= 0
hs_minus = all(inn == (not vor) for vor, inn in einseitig["N-"])     # genau für s <= 0
hs_kontr = not any(inn for _, inn in einseitig["Kontrolle"])         # nie
ZAHLEN["exakt"] = {"v_O1": str(vO), "v_V+": str(vV), "v_N+": str(vN), "v_N-": str(vNm),
                   "(K_O1 - K_N+)/2pi": str(sp.simplify(dP)), "(K_V+ - K_N+)/2pi": str(sp.simplify(dM)),
                   "(K_O1 - K_N-)/2pi": str(sp.simplify(dM2)), "Faktor cot(pi/4) - cot(pi/2)": str(fac),
                   "Kontrolle (K_O1 - K_Ic)/2pi": "1 - (sqrt(2) - 1)*cos(theta)", "a^2 - b^2 der Kontrolle": str(ell_disk),
                   "Bild von N+ = (0, 1) für e^(s/2) = 2": [str(moeb(sp.Integer(2), sp.Integer(0))), "1"],
                   "Bild der Kontrolle (-1/2, 1/2) für e^(s/2) = 2": [str(moeb(sp.Integer(2), sp.Rational(-1, 2))),
                                                                      str(moeb(sp.Integer(2), sp.Rational(1, 2)))]}
check("[E] Native Gebiete und ihre Möbiusfelder (sympy): v_O1 = cos, v_V+ = sin, v_N+ = cos + sin - 1, v_N- = cos - sin - 1 "
      "(Nullstellen genau in den Marken). Wiesbrock-Differenzen: (K_O1 - K_N+)/2 pi = G_(1 - sin) = P^(pi/2), "
      "(K_V+ - K_N+)/2 pi = G_(1 - cos) = P^(0) = K_0, (K_O1 - K_N-)/2 pi = G_(1 + sin) = P^(-pi/2): parabolisch (doppelte "
      f"Nullstelle in der gemeinsamen Marke), Faktor cot(pi/4) - cot(pi/2) = {fac} (Einheit der Marken). Lie-Relationen "
      "[K_O1, P^(pi/2)] = +2 pi i P^(pi/2), [K_V+, P^(0)] = -2 pi i P^(0), [K_O1, P^(-pi/2)] = -2 pi i P^(-pi/2). Die Uhr ist "
      "der Moduldefekt der drei nativen Gebiete: K_O1 + K_V+ - K_(O1 cap V+) = 2 pi L_0 = K_Uhr/2",
      felder_ok and diff_ok and fac == 1 and uhr_ok and lie_ok)
check("[E] Einseitigkeit exakt (rationale Möbiusflüsse der Kappe, e^(s/2) = 2, 3, 7/5 und Kehrwerte, s = 2 pi t): N_+ = (0, 1) "
      f"wird genau für s >= 0 in sich abgebildet (Bild für e^(s/2) = 2: ({moeb(sp.Integer(2), sp.Integer(0))}, 1)), N_- = "
      "(-1, 0) genau für s <= 0: zwei halbseitige Inklusionen entgegengesetzten Vorzeichens im selben M = O_1. "
      "Negativkontrolle: der streng geschachtelte Bogen (-1/2, 1/2) ohne gemeinsamen Endpunkt verlässt sich in beiden "
      f"Richtungen (Bild ({moeb(sp.Integer(2), sp.Rational(-1, 2))}, {moeb(sp.Integer(2), sp.Rational(1, 2))})); seine Differenz "
      f"ist G_(1 - (sqrt 2 - 1) cos) (elliptisch, a^2 - b^2 = {ell_disk} > 0) und erfüllt die Lie-Relation in keiner Orientierung",
      feld_y and hs_plus and hs_minus and hs_kontr and kontr_feld and ell_disk > 0 and not any(kontr_lie))

say("  (b) aus dem Zustand: Gitter-Chiralfermion im Vakuum, Modulgeneratoren hochgenau (640 Bit)")
PREC = 640


def hp_rows(M, an, ad, cfun):
    """Korrelationsfunktion G(j - k) = (1/M) sum_x c(x) e^(2 pi i x (j - k)/M), x = k - M/2 + an/ad (wie z72)."""
    ctx.prec = PREC
    xs = [arb(k - M//2) + arb(an)/ad for k in range(M)]
    cs = [cfun(x) for x in xs]
    rows = {}
    for d in range(-(M - 1), M):
        acc = acb(0)
        for k in range(M):
            acc += cs[k]*acb(2*xs[k]*d/M).exp_pi_i()
        rows[d] = acc/M
    return rows


def c_vak(x):
    return arb(1) if float(x.mid()) < 0 else arb(0)


def c_therm(beta_durch_pi):
    def cf(x):
        return 1/(1 + (arb.pi()*beta_durch_pi*x).exp())
    return cf


def zu_np(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))


def eig_hp(rows, A):
    """Eigenzerlegung von G_AA: h_k = log((1 - nu)/nu), normierte Eigenvektoren (numpy) und die flint-Daten."""
    ctx.prec = PREC
    n = len(A)
    GA = acb_mat([[rows[j - k] for k in A] for j in A])
    E, R = GA.eig(right=True, algorithm="approx")
    hs, nus, norms = [], [], []
    V = np.zeros((n, n), complex)
    for c in range(n):
        nu = E[c].real
        nv = arb(0)
        for i in range(n):
            nv += abs(R[i, c])**2
        nv = nv.sqrt()
        hs.append(float(((1 - nu)/nu).log().mid()))
        nus.append(nu)
        norms.append(nv)
        for i in range(n):
            V[i, c] = zu_np(R[i, c]/nv)
    return np.array(hs), V, R, nus, norms


def K_rein(rows, M, A):
    """Voller Einteilchen-Modulgenerator des Bogens A im reinen Zustand (Platzbasis 0..M-1), dazu h und V auf A."""
    ctx.prec = PREC
    Ac = [j for j in range(M) if j not in set(A)]
    hs, V, R, nus, norms = eig_hp(rows, A)
    W_hp = acb_mat([[rows[j - k] for k in A] for j in Ac])*R
    W = np.zeros((len(Ac), len(A)), complex)
    for c in range(len(A)):
        sc = (nus[c]*(1 - nus[c])).sqrt()*norms[c]
        for i in range(len(Ac)):
            W[i, c] = zu_np(W_hp[i, c]/sc)
    k = np.zeros((M, M), complex)
    k[np.ix_(A, A)] = (V*hs) @ V.conj().T
    k[np.ix_(Ac, Ac)] = (W*hs) @ W.conj().T
    orth = max(float(np.abs(V.conj().T @ V - np.eye(len(A))).max()), float(np.abs(W.conj().T @ W - np.eye(len(A))).max()))
    return k, hs, V, orth


def mob_moden(xs):
    M = len(xs)
    L0 = np.diag(xs).astype(complex)
    L1 = np.zeros((M, M), complex)
    Lm1 = np.zeros((M, M), complex)
    for k in range(M):
        if k >= 1:
            L1[k - 1, k] = xs[k] - 0.5
        if k + 1 < M:
            Lm1[k + 1, k] = xs[k] + 0.5
    return L0, (L1 + Lm1)/2, 1j*(L1 - Lm1)/2


def dreh(xs, phi):
    """Drehung der Naht um phi in der Modenbasis: f(theta) -> f(theta - phi)."""
    return np.diag(np.exp(-1j*xs*phi))


GIT, LECK_DATEN = {}, {}
for M in (128, 256):
    xs = np.arange(M) - M//2 + 0.5
    thj = 2*np.pi*(np.arange(M) + 0.5)/M
    Um = np.exp(1j*np.outer(thj, xs))/np.sqrt(M)
    rows = hp_rows(M, 1, 2, c_vak)
    q = M//4
    A_O = [j % M for j in range(-q, q)]
    KOs, hO, VO, orthO = K_rein(rows, M, A_O)
    KNs, hN, VN, orthN = K_rein(rows, M, list(range(0, q)))
    KO, KN = Um.conj().T @ KOs @ Um, Um.conj().T @ KNs @ Um
    KV = dreh(xs, np.pi/2) @ KO @ dreh(xs, np.pi/2).conj().T
    KNm = dreh(xs, -np.pi/2) @ KN @ dreh(xs, -np.pi/2).conj().T
    KIc = dreh(xs, -np.pi/4) @ KN @ dreh(xs, -np.pi/4).conj().T
    L0, H1, H2 = mob_moden(xs)
    lo = np.where(np.abs(xs) < 2)[0]

    def rel(Aop, Bop, lo=lo):
        return float(np.abs((Aop - Bop)[np.ix_(lo, lo)]).max()/np.abs(Bop[np.ix_(lo, lo)]).max())
    Pp, Pm, Pm2, Dc = (KO - KN)/(2*np.pi), (KV - KN)/(2*np.pi), (KO - KNm)/(2*np.pi), (KO - KIc)/(2*np.pi)
    Pi = np.diag((xs < 0).astype(float))
    lc = np.where(np.abs(xs) < 3)[0]

    def fock_min(P, Pi=Pi, lc=lc):
        X = P @ (np.eye(len(Pi)) - 2*Pi)
        X = (X + X.conj().T)/2
        return float(np.linalg.eigvalsh(X[np.ix_(lc, lc)])[0])
    GIT[M] = {
        "hl": max(rel(KO, 2*np.pi*H1), rel(KV, 2*np.pi*H2), rel(KN, 2*np.pi*(H1 + H2 - L0))),
        "wiesbrock": max(rel(Pp, L0 - H2), rel(Pm, L0 - H1), rel(Pm2, L0 + H2)),
        "uhr": rel(KO + KV - KN, 2*np.pi*L0),
        "borchers": max(rel(komm(KO, Pp), 2j*np.pi*Pp), rel(komm(KV, Pm), -2j*np.pi*Pm), rel(komm(KO, Pm2), -2j*np.pi*Pm2)),
        "borchers_falsch": min(rel(komm(KO, Pp), -2j*np.pi*Pp), rel(komm(KV, Pm), 2j*np.pi*Pm)),
        "komm_N": rel(komm(Pp, Pm), -1j*KN/(2*np.pi)),
        "komm_M": rel(komm(Pp, Pm2), -1j*KO/np.pi),
        "komm_null": float(np.abs(komm(Pp, Pm)[np.ix_(lo, lo)]).max()/np.abs(Pp[np.ix_(lo, lo)]).max()),
        "fock_min": min(fock_min(Pp), fock_min(Pm), fock_min(Pm2)),
        "kontr_feld": rel(Dc, L0 - (np.sqrt(2) - 1)*H1),
        "kontr_luecke": fock_min(Dc),
        "kontr_lie": min(rel(komm(KO, Dc), 2j*np.pi*Dc), rel(komm(KO, Dc), -2j*np.pi*Dc)),
        "orth": max(orthO, orthN), "hmax": float(max(np.abs(hO).max(), np.abs(hN).max())),
    }
    LECK_DATEN[M] = (q, hO, VO)
g1, g2 = GIT[128], GIT[256]
ZAHLEN["gitter_vakuum"] = {str(M): {k: r(v, 3) for k, v in GIT[M].items()} for M in (128, 256)}
say(f"   M = 128: { {k: r(v, 3) for k, v in g1.items()} }")
say(f"   M = 256: { {k: r(v, 3) for k, v in g2.items()} }")
check("[N] Hislop-Longo als Resultat aus dem Zustand (Peschel, 640 Bit, niedrige Moden |x| < 2): K_O1 = 2 pi H_1, K_V+ = "
      f"2 pi H_2, K_N+ = 2 pi (H_1 + H_2 - L_0) mit rel. Abweichung {g1['hl']:.1e} (M = 128) -> {g2['hl']:.1e} (M = 256); "
      f"Wiesbrock-Differenzen = P^(pi/2), P^(0), P^(-pi/2) mit Faktor 1 ({g1['wiesbrock']:.1e} -> {g2['wiesbrock']:.1e}); "
      f"Uhr K_O1 + K_V+ - K_N+ = 2 pi L_0 ({g1['uhr']:.1e} -> {g2['uhr']:.1e}). Die Abweichung fällt mit der Gitterweite; "
      f"Rechengenauigkeit: grösstes |h| = {g2['hmax']:.0f} bei {PREC*np.log(2):.0f} verfügbaren Nepern",
      g2["hl"] < 5e-3 and g2["wiesbrock"] < 5e-3 and g2["uhr"] < 5e-3 and g2["hl"] < 0.6*g1["hl"]
      and g2["wiesbrock"] < 0.6*g1["wiesbrock"] and g2["orth"] < 1e-12 and g2["hmax"] < PREC*np.log(2) - 100)
check("[N] Borchers-Relationen aus dem Zustand (M = 256): [K_M, P] = +-2 pi i P mit den exakt vorhergesagten Vorzeichen "
      f"(rel. Defekt {g2['borchers']:.1e}; falsche Orientierung {g2['borchers_falsch']:.2f}); Fock-Positivität der drei "
      f"Translationserzeuger auf der Kompression |x| < 3 (min {g2['fock_min']:.3f} > 0). Negativkontrolle (K_O1 - K_Ic)/2 pi "
      f"mit Ic = (-pi/4, pi/4): gleich L_0 - (sqrt 2 - 1) H_1 ({g2['kontr_feld']:.1e}), positiv mit Lücke "
      f"{g2['kontr_luecke']:.3f} (vorhergesagt sqrt(2 sqrt 2 - 2)/2 = {np.sqrt(2*np.sqrt(2) - 2)/2:.3f}), Lie-Relation verfehlt "
      f"in beiden Orientierungen ({g2['kontr_lie']:.2f}): Positivität allein macht keine Translation (z72)",
      g2["borchers"] < 1e-2 and g2["borchers_falsch"] > 1.0 and g2["fock_min"] > 0.05 and g2["kontr_feld"] < 5e-3
      and abs(g2["kontr_luecke"] - np.sqrt(2*np.sqrt(2) - 2)/2) < 5e-3 and g2["kontr_lie"] > 0.3)

say("  (c) die Definition direkt: Leckage eines Pakets unter dem Modulfluss von O_1")


def leckage(hs, V, M, q, zentrum, drin_fn, sigma=0.12, zeiten=(0.1, 0.2, 0.3)):
    """max Leckage ||(e^{-i t h} f) ausserhalb||^2 für t > 0 und t < 0; f Gausspaket um zentrum im Bogen.
    Die Plätze j < 0 stehen in G mit umlaufendem Index (antiperiodisch): Eichvorzeichen -1."""
    jj = np.arange(-q, q)
    th = 2*np.pi*(jj + 0.5)/M
    innen = np.array([drin_fn(j) for j in jj])
    f = np.exp(-(th - zentrum)**2/(2*sigma**2))*np.where(jj < 0, -1.0, 1.0)
    f = f.astype(complex)/np.linalg.norm(f)
    fh = V.conj().T @ f
    vor = max(float(np.linalg.norm((V @ (np.exp(-1j*t*hs)*fh))[~innen])**2) for t in zeiten)
    rueck = max(float(np.linalg.norm((V @ (np.exp(1j*t*hs)*fh))[~innen])**2) for t in zeiten)
    return vor, rueck, float(np.linalg.norm(f[~innen])**2)


LK = {}
for M in (128, 256):
    q, hOl, VOl = LECK_DATEN[M]
    LK[f"Vakuum N+, M = {M}"] = leckage(hOl, VOl, M, q, np.pi/4, lambda j: j >= 0)
    LK[f"Vakuum N-, M = {M}"] = leckage(hOl, VOl, M, q, -np.pi/4, lambda j: j < 0)
    LK[f"Vakuum Kontrolle Ic, M = {M}"] = leckage(hOl, VOl, M, q, 0.0, lambda j, q=q: -q//2 <= j < q//2)
    rows_o = hp_rows(M, 1, 4, c_therm(4))
    hT, VT, _, _, _ = eig_hp(rows_o, [j % M for j in range(-q, q)])
    LK[f"Original N+, M = {M}"] = leckage(hT, VT, M, q, np.pi/4, lambda j: j >= 0)
    LK[f"Original N-, M = {M}"] = leckage(hT, VT, M, q, -np.pi/4, lambda j: j < 0)
ZAHLEN["leckage"] = {k: {"vorwaerts": r(v[0], 2), "rueckwaerts": r(v[1], 2), "Start": r(v[2], 2)} for k, v in LK.items()}
for k, v in LK.items():
    say(f"   {k}: vorwärts {v[0]:.1e}, rückwärts {v[1]:.1e}, Start {v[2]:.1e}")
lv, lv1 = LK["Vakuum N+, M = 256"], LK["Vakuum N+, M = 128"]
lm, lc_ = LK["Vakuum N-, M = 256"], LK["Vakuum Kontrolle Ic, M = 256"]
lo4, lo41 = LK["Original N+, M = 256"], LK["Original N+, M = 128"]
check("[N] Halbseitigkeit aus dem Zustand (Definition direkt; Gausspaket der Breite 0,12 in der Bogenmitte, |t| = 0,1; 0,2; "
      f"0,3; M = 256): im Vakuum bleibt das Paket aus N_+ für t > 0 in N_+ (Leckage {lv[0]:.0e}; bei M = 128 noch "
      f"{lv1[0]:.0e}, Gitterartefakt der gestauchten Paketbreite) und verlässt es für t < 0 ({lv[1]:.2f}); für N_- umgekehrt "
      f"({lm[1]:.0e} gegen {lm[0]:.2f}): beide Vorzeichen nativ. Negativkontrolle Ic = (-pi/4, pi/4): Leckage in beiden "
      f"Richtungen ({lc_[0]:.2f}, {lc_[1]:.2f}). Originalzustand (alpha = 1/4, beta = 4 pi): Vorwärtsleckage {lo4[0]:.1e} "
      f"(M = 128: {lo41[0]:.1e}, konvergiert), über 10^6 mal das Vakuum: keine halbseitige Inklusion (z72 bestätigt)",
      lv[0] < 1e-12 and lv[0] < lv1[0] and lv[1] > 0.5 and lm[1] < 1e-12 and lm[0] > 0.5 and min(lc_[:2]) > 0.3
      and lo4[0] > 1e-6 and lo4[0] > 1e6*max(lv[0], 1e-30) and abs(lo4[0]/lo41[0] - 1) < 0.2
      and LK["Original N-, M = 256"][1] > 1e-6 and max(v[2] for v in LK.values()) < 1e-15)

say("  (d) Positivität exakt: Spektralmass auf der diskreten Serie")


def ds_sym(h, N):
    L0s = sp.diag(*[h + n for n in range(N)])
    Lm1s = sp.zeros(N, N)
    for n in range(N - 1):
        Lm1s[n + 1, n] = sp.sqrt((n + 1)*(n + 2*h))
    return L0s, Lm1s.T, Lm1s


mom_ok, sym_ok = True, True
for hh in (sp.Rational(1, 2), sp.Integer(1)):
    Nn = 7
    L0s, L1s, Lm1s = ds_sym(hh, Nn)
    H1s, H2s = (L1s + Lm1s)/2, sp.I*(L1s - Lm1s)/2
    for Ps in (L0s - H2s, L0s + H2s, L0s - H1s):
        Pk = sp.eye(Nn)
        for k in range(1, 7):
            Pk = Pk*Ps
            mom_ok &= sp.simplify(Pk[0, 0] - sp.rf(2*hh, k)/2**k) == 0
    sym_ok &= H1s[0, 0] == 0 and sp.simplify((H1s**2)[0, 0] - hh/2) == 0 and sp.simplify((H1s**3)[0, 0]) == 0


def ds_num(hh, N):
    n = np.arange(N)
    L0n = np.diag(hh + n).astype(complex)
    Lm1n = np.diag(np.sqrt((n[:-1] + 1)*(n[:-1] + 2*hh)), -1).astype(complex)
    return L0n, (Lm1n.T + Lm1n)/2, 1j*(Lm1n.T - Lm1n)/2


spek = {}
for Nn in (200, 800):
    L0n, H1n, H2n = ds_num(0.5, Nn)
    spek[Nn] = (float(np.linalg.eigvalsh(L0n - H2n)[0]), float(np.linalg.eigvalsh(L0n - (np.sqrt(2) - 1)*H1n)[0]),
                float(np.linalg.eigvalsh(H1n)[0]))
luecke_ex = 0.5*np.sqrt(1 - (np.sqrt(2) - 1)**2)
ZAHLEN["positivitaet"] = {"min Eigenwert (P^(pi/2), elliptische Kontrolle, K_O1/2pi) bei N = 200": [r(x, 5) for x in spek[200]],
                          "bei N = 800": [r(x, 5) for x in spek[800]], "Lücke h sqrt(a^2 - b^2), h = 1/2": r(luecke_ex, 6)}
check("[E] Positivität exakt (sympy, D_1/2 und D_1): <h|P^k|h> = (2h)_k/2^k für k = 1..6 und alle drei nativen "
      "Translationserzeuger P^(pi/2), P^(-pi/2), P^(0): das Spektralmass des Grundzustands ist Gamma(2h, 1/2) auf [0, oo), "
      f"absolut stetig, ohne Lücke (kleinster Eigenwert der Trunkierung {spek[200][0]:.4f} -> {spek[800][0]:.4f}, ~ 1/N). "
      f"Kontrollen: die elliptische Differenz hat die Lücke {spek[800][1]:.4f} = h sqrt(a^2 - b^2) (Uhr, keine Translation); der "
      "Modulgenerator K_O1/2 pi = H_1 ist spektral symmetrisch (<H_1> = <H_1^3> = 0, <H_1^2> = h/2) und nach unten "
      f"unbeschränkt ({spek[200][2]:.0f} -> {spek[800][2]:.0f})",
      mom_ok and sym_ok and 3.5 < spek[200][0]/spek[800][0] < 4.5 and abs(spek[800][1] - luecke_ex) < 1e-6
      and abs(spek[200][1] - luecke_ex) < 1e-6 and spek[800][2] < 2*spek[200][2] < 0)

say("  (e) beide Vorzeichen: Kommutator der zwei Translationen")
check("[N] Die zwei nativen halbseitigen Translationen einer Naht vertauschen nicht (M = 256, aus dem Zustand): "
      f"[P^(pi/2), P^(0)] = -i K_(N+)/2 pi (rel. Defekt {g1['komm_N']:.1e} -> {g2['komm_N']:.1e}; gleiches N in O_1 und V_+) "
      f"und [P^(pi/2), P^(-pi/2)] = -i K_(O1)/pi ({g1['komm_M']:.1e} -> {g2['komm_M']:.1e}; gleiches M, N_+ und N_-); der Kommutator hat relative Grösse "
      f"{g2['komm_null']:.2f}, nicht 0. Der Kommutator der beiden Translationen ist der Modulgenerator des gemeinsamen Gebiets",
      g2["komm_N"] < 2e-2 and g2["komm_M"] < 2e-2 and g2["komm_null"] > 0.3 and g2["komm_N"] < 0.8*g1["komm_N"]
      and g2["komm_M"] < 0.8*g1["komm_M"])
# 2x2-Darstellung: l0 = s3/2, h1 = i s1/2, h2 = i s2/2 erfüllt [l0,h1] = i h2, [l0,h2] = -i h1, [h1,h2] = -i l0
S1 = np.array([[0, 1], [1, 0]], complex)
S2 = np.array([[0, -1j], [1j, 0]])
S3 = np.diag([1.0, -1.0]).astype(complex)
l0, h1, h2 = S3/2, 1j*S1/2, 1j*S2/2
konv2 = max(nrm(komm(l0, h1) - 1j*h2), nrm(komm(l0, h2) + 1j*h1), nrm(komm(h1, h2) + 1j*l0))
pP, pM, pM2 = l0 - h2, l0 - h1, l0 + h2
alg1 = lie_huelle([1j*pP, 1j*pM])                                   # reelle Form: i mal hermitescher Erzeuger
alg1b = lie_huelle([1j*pP, 1j*pM2])
kr1, _ = killing_rang(alg1)
Z2 = np.zeros((2, 2), complex)


def block(A, B):
    return np.block([[A, Z2], [Z2, B]])


# zwei Arten a, b: K = k (+) k; N_+ verkleinert Art a am rechten Ende, N_- Art b am linken Ende (P_chir)
Kc, Pa, Pb = block(2*np.pi*h1, 2*np.pi*h1), block(pP, Z2), block(Z2, pM2)
alg2 = lie_huelle([1j*Kc, 1j*Pa, 1j*Pb])
kr2, _ = killing_rang(alg2)
gew = (nrm(komm(Kc, Pa) - 2j*np.pi*Pa), nrm(komm(Kc, Pb) + 2j*np.pi*Pb), nrm(komm(Pa, Pb)))
# Kontrollen: artblind (beide Arten gleich verkleinert) und gleiches Ende für beide Arten
blind = nrm(komm(block(pP, pP), block(pM2, pM2)))
Pb_gleich = block(Z2, pP)
gew_gleich = nrm(komm(Kc, Pb_gleich) + 2j*np.pi*Pb_gleich)
ZAHLEN["zwei_vorzeichen"] = {"dim Lie-Hülle {P^(pi/2), P^(0)}": len(alg1), "dim Lie-Hülle {P^(pi/2), P^(-pi/2)}": len(alg1b),
                             "Killing-Rang (eine Naht)": kr1, "dim Lie-Hülle unter P_chir {K, P_a, P_b}": len(alg2),
                             "Killing-Rang unter P_chir": kr2, "P_chir: Defekte [K,P_a] = +2pi i P_a, [K,P_b] = -2pi i P_b, "
                             "[P_a,P_b] = 0": [r(x, 2) for x in gew], "Kontrolle artblind ||[P_+, P_-]||": r(blind, 4),
                             "Kontrolle gleiches Ende: Defekt des Gewichts -1": r(gew_gleich, 4)}
check(f"[E] Lie-Algebra (2x2-Darstellung, Konventionsfehler {konv2:.0e}): die beiden Translationen einer Naht erzeugen "
      f"sl(2,R) (Dimension {len(alg1)} bzw. {len(alg1b)}, Killing-Form vom Rang {kr1}: halbeinfach), also die Möbiusgruppe eines "
      "Lichtstrahls, nicht zwei vertauschende Lichtkegelkoordinaten. [C unter P_chir] Mit zwei Arten a, b und "
      "artabhängiger Verkleinerung (a am einen, b am anderen Ende) gilt [K, P_a] = +2 pi i P_a, [K, P_b] = -2 pi i P_b, "
      f"[P_a, P_b] = 0 (Defekte {max(gew):.0e}): die 1+1-Poincaré-Algebra (Dimension {len(alg2)}, Killing-Rang {kr2}: auflösbar). "
      f"Kontrollen: artblind vertauschen die Translationen nicht ({blind:.2f}); gleiches Ende für beide Arten gibt zweimal "
      f"dasselbe Vorzeichen (Defekt {gew_gleich:.1f}): kein +-Paar",
      konv2 < 1e-14 and len(alg1) == 3 and len(alg1b) == 3 and kr1 == 3 and len(alg2) == 3 and kr2 == 1 and max(gew) < 1e-12
      and blind > 0.5 and gew_gleich > 1.0)

say("  (f) su(2,2) unter z83s Identifikation: halbseitige Familien des Kegels und der Kappe")
I2 = np.eye(2, dtype=complex)
SIG = [I2, S1, S2, S3]
Hf = np.block([[Z2, I2], [I2, Z2]])                                  # (2,2)-Form des Twistorraums


def XP(a):
    return np.block([[Z2, 1j*sum(a[m]*SIG[m] for m in range(4))], [Z2, Z2]])


def XK(b):
    return np.block([[Z2, Z2], [1j*sum(b[m]*SIG[m] for m in range(4)), Z2]])


XD = np.diag([0.5, 0.5, -0.5, -0.5]).astype(complex)
XROT = [np.block([[1j*SIG[k]/2, Z2], [Z2, 1j*SIG[k]/2]]) for k in (1, 2, 3)]
XB3 = np.block([[S3/2, Z2], [Z2, -S3/2]])                            # Boost längs 3
EB = np.eye(4)
XE, XQ = (XP(EB[0]) + XK(EB[0]))/2, (XP(EB[0]) - XK(EB[0]))/2        # L_0 = E, H_1 = Q^-; H_2 = X_D
in_alg = max(nrm(X.conj().T @ Hf + Hf @ X) for X in [XD, XB3, XE, XQ] + XROT + [XP(EB[m]) for m in range(4)]
             + [XK(EB[m]) for m in range(4)])
radial = max(nrm(komm(XE, XQ) - XD), nrm(komm(XE, XD) + XQ), nrm(komm(XQ, XD) + XE))        # X = -i G


def kegelform(X):
    """Hermitesche Form Q_X = -i H X; X liegt im invarianten Positivkegel genau dann, wenn Q_X >= 0."""
    Qm = -1j*Hf @ X
    return np.linalg.eigvalsh((Qm + Qm.conj().T)/2)


VEKT = {"zeitartig (1, 0, 0, 0)": [1, 0, 0, 0], "zeitartig (5/4, 0, 0, 3/4)": [1.25, 0, 0, 0.75],
        "lichtartig (1, 0, 0, 1)": [1, 0, 0, 1], "lichtartig (1, 0, 0, -1)": [1, 0, 0, -1],
        "raumartig (0, 0, 0, 1)": [0, 0, 0, 1], "raumartig (3/4, 0, 0, 5/4)": [0.75, 0, 0, 1.25],
        "vergangen (-1, 0, 0, 0)": [-1, 0, 0, 0]}
kegel = {k: [r(x, 6) for x in kegelform(XP(v))] for k, v in VEKT.items()}
kegel_ok = all((min(kegelform(XP(v))) > -1e-12) == (v[0] >= np.linalg.norm(v[1:]) - 1e-12) for v in VEKT.values())
wies4 = max(nrm((XD - expm(XP(v)) @ XD @ expm(-XP(v))) - XP(v)) for v in VEKT.values())     # K_V+ - K_(V+ + a) ~ a.P
skal = max(nrm(expm(0.7*XD) @ XP(v) @ expm(-0.7*XD) - np.exp(0.7)*XP(v)) for v in VEKT.values())
abelsch = max(nrm(komm(XP(EB[m]), XP(EB[n]))) for m in range(4) for n in range(4))
alg_P = lie_huelle([XP(EB[m]) for m in range(4)])
alg_PK = lie_huelle([XP(EB[m]) for m in range(4)] + [XK(EB[m]) for m in range(4)])
# so(3)-invariante Translationen: Kern von a -> ([X_rot_k, X_P(a)])_k
Mrot = np.array([np.concatenate([np.concatenate([np.real(komm(Xr, XP(EB[m]))).ravel(), np.imag(komm(Xr, XP(EB[m]))).ravel()])
                                 for Xr in XROT]) for m in range(4)]).T
dim_inv = 4 - int(np.linalg.matrix_rank(Mrot, tol=1e-10))
ap, am = XP([1, 0, 0, 1]), XP([1, 0, 0, -1])
null_paar = (nrm(komm(ap, am)), nrm(komm(XB3, ap) - ap), nrm(komm(XB3, am) + am), nrm(komm(XD, ap) - ap), nrm(komm(XD, am) - am))
# Kappe: N_+ = O_1 cap V_+ ; Differenzen E - D und E - Q^- = K_0; Spitzenfamilie = Ad(e^(phi E)) span{P_mu}
XN = XQ + XD - XE
tip = expm(-np.pi/2*XE)
fam_tip = [tip @ XP(EB[m]) @ np.linalg.inv(tip) for m in range(4)]
fam_gegen = [np.linalg.inv(tip) @ XP(EB[m]) @ tip for m in range(4)]                # Spitze t = -1 (Kontrolle)
in_fam = (len(lie_huelle(fam_tip)), len(lie_huelle(fam_tip + [XE - XD])), len(lie_huelle(fam_gegen + [XE - XD])))
nilp = max(nrm((XE - XD) @ (XE - XD)), nrm((XE - XQ) @ (XE - XQ)))
kappe_kegel = (float(min(kegelform(XE - XD))), float(min(kegelform(XE - XQ))), [r(x, 6) for x in kegelform(XQ)])
komm4 = nrm(komm(XE - XD, XE - XQ) + XN)                             # X = -i G: [X_P, X_P'] = -X_N  <=>  [P, P'] = -i K_N/2 pi
ZAHLEN["su22"] = {"Kegelform -i H X_P(a), Eigenwerte": kegel, "dim Lie-Hülle {P_mu}": len(alg_P),
                  "dim Lie-Hülle {P_mu, K_mu}": len(alg_PK), "dim so(3)-invarianter Translationen": dim_inv,
                  "Nullpaar: [P+,P-], Boostgewichte +1/-1, Dilatationsgewichte +1/+1 (Defekte)": [r(x, 2) for x in null_paar],
                  "Kegelform von Q^- (Modulgenerator der Kappe)": kappe_kegel[2]}
check(f"[E] su(2,2) auf dem Twistorraum (Form (2,2), Erzeuger in der Algebra: {in_alg:.0e}; radiale sl(2,R) "
      f"E = L_0, Q^- = H_1, D = H_2: {radial:.0e}). Kegel V_+ (Modulgenerator = Dilatation): Für jedes a gilt "
      f"K_(V+) - K_(V+ + a) ~ a.P exakt ({wies4:.0e}) und Ad(e^(sD)) a.P = e^s a.P ({skal:.0e}); V_+ + a liegt in V_+ genau für a "
      "im abgeschlossenen Vorwärtskegel, und genau dort liegt a.P im invarianten Positivkegel (-i H X >= 0; zeitartig und "
      "lichtartig ja, raumartig und vergangen nein). Die halbseitige Familie eines Gebiets mit fester Spitze ist also die "
      f"ganze Translationsgruppe R^4 (Dimension {len(alg_P)}, abelsch {abelsch:.0e}) mit Spektralbedingung; die Familien zweier "
      f"Spitzen ({{P_mu}} und {{K_mu}}) erzeugen die ganze konforme Algebra (Dimension {len(alg_PK)})",
      in_alg < 1e-14 and radial < 1e-14 and wies4 < 1e-12 and skal < 1e-12 and kegel_ok and len(alg_P) == 4 and abelsch < 1e-14
      and len(alg_PK) == 15)
check("[E/C] Native Auswahl in su(2,2): Unter den vier Translationen ist genau eine so(3)-invariant (Dimension "
      f"{dim_inv}: die Zeit P_0); eine Raum- oder Lichtrichtung verlangt die Wahl einer Achse n in S^2. Für ein Nullpaar "
      f"a_+- = e_0 +- n vertauschen die Translationen ({null_paar[0]:.0e}), tragen unter dem Boost längs n die Gewichte +1 und -1 "
      f"(Keil: +-Paar, {max(null_paar[1:3]):.0e}) und unter der Dilatation beide +1 (Kegel: gleiches Vorzeichen, "
      f"{max(null_paar[3:]):.0e}). Native Kappe: (K_O1 - K_N+)/2 pi = E - D und (K_V+ - K_N+)/2 pi = E - Q^- = K_0 sind nilpotent "
      f"({nilp:.0e}), liegen im Positivkegel (min {min(kappe_kegel[:2]):.0e}) und E - D liegt in der Spitzenfamilie "
      f"Ad(e^(-pi E/2)) span{{P_mu}} der oberen Kappenspitze (Lie-Hülle bleibt {in_fam[1]}-dimensional; mit der Familie der "
      f"unteren Spitze dagegen {in_fam[2]}); ihr Kommutator ist der Modulgenerator von N_+ ({komm4:.0e}). Der Modulgenerator "
      f"Q^- selbst ist indefinit ({kappe_kegel[2]})",
      dim_inv == 1 and max(null_paar) < 1e-14 and nilp < 1e-14 and min(kappe_kegel[:2]) > -1e-12 and in_fam == (4, 4, 15)
      and komm4 < 1e-14
      and min(kappe_kegel[2]) < -0.1 and max(kappe_kegel[2]) > 0.1)
befund("E/N: Die halbseitige modulare Inklusion ist nativ, im Vakuumsektor der rekonstruierten Uhr: Kappe O_1 und Kegel V_+ "
       "(z104, aus E und tau) mit N_+ = O_1 cap V_+ bzw. N_- geben halbseitige Inklusionen beider Vorzeichen, exakt geometrisch "
       "und am Gitter aus dem Zustand; die Wiesbrock-Differenzen sind die parabolischen Erzeuger in den Marken mit Faktor 1, "
       "positiv mit Gamma-Spektralmass. Im Originalzustand (beta = 4 pi) gibt es sie nicht (z72: Leckage, global Takesaki)", True)
befund("E: Die zwei Vorzeichen einer Naht geben keine zwei Lichtkegelkoordinaten: Ihre Translationen vertauschen nicht, der "
       "Kommutator ist der Modulgenerator des gemeinsamen Gebiets, zusammen sl(2,R) (Möbiusgruppe eines Lichtstrahls). Zwei "
       "vertauschende lichtartige Translationen brauchen zwei unabhängige Lichtstrahlen: zwei Arten mit entgegengesetzter "
       "Zuordnung (P_chir) oder in su(2,2) (z83, C) ein Nullpaar zu einer gewählten Achse n (P_richt)", True)

# =================================================================== Teil 2
say("\nTeil 2: Die Kette N_n = U(n) M U(-n) und ihre Plätze")
say("  (a) Plätze der kanonischen Kette (M = V_+ = (0, oo), N = (1, oo), P = P_0 = -i d_y, y = tan(theta/2) = t)")
n_s = sp.symbols("n", integer=True)
ell = 2*(sp.atan(n_s + 1) - sp.atan(n_s))                           # Bogenlänge des Platzes (n, n+1)
N_s = sp.symbols("N", positive=True, integer=True)
summe = sp.limit(4*sp.atan(N_s), N_s, sp.oo)                        # sum_(n=-N)^(N-1) ell_n = 4 arctan N (Teleskop)
tele = all(sp.simplify(sum(ell.subs(n_s, k) for k in range(-Nv, Nv)) - 4*sp.atan(Nv)) == 0 for Nv in (1, 2, 5))


def slepian(l_):
    """Eigenwerte der Gram-Matrix zweier Nahtmoden mit Abstand 1 (Dublett) auf einem Bogen der Länge l."""
    return l_/(2*sp.pi) + sp.sin(l_/2)/sp.pi, l_/(2*sp.pi) - sp.sin(l_/2)/sp.pi


lam0 = slepian(sp.pi/2)
x_s = sp.symbols("x", real=True)
gram_ex = sp.integrate(sp.exp(sp.I*TH), (TH, 0, sp.pi/2))/(2*sp.pi)
gram_ok = sp.simplify(sp.Abs(gram_ex) - sp.sqrt(2)/(2*sp.pi)) == 0 and sp.simplify(lam0[0] - (sp.Rational(1, 4) + sp.sqrt(2)/(2*sp.pi))) == 0
anteile = {k: float(ell.subs(n_s, k)/(2*sp.pi)) for k in range(-4, 4)}
lam_max = {k: float(slepian(ell.subs(n_s, k))[0]) for k in range(-4, 4)}
# Gitterprobe der Gram-Matrix (M = 1024, Platz 0 und Platz 1)
Mg = 1024
thg = 2*np.pi*(np.arange(Mg) + 0.5)/Mg
gram_num = {}
for k in (0, 1):
    a_, b_ = 2*np.arctan(k), 2*np.arctan(k + 1)
    sel = (thg > a_) & (thg < b_)
    Gm = np.array([[np.sum(np.exp(1j*(xb - xa)*thg[sel]))/Mg for xb in (0.5, 1.5)] for xa in (0.5, 1.5)])
    gram_num[k] = float(np.linalg.eigvalsh(Gm)[-1])
ZAHLEN["plaetze"] = {"Anteil einer Nahtmode je Platz n = -4..3 (ell_n/2pi)": {str(k): r(v, 5) for k, v in anteile.items()},
                     "grösster Dublettanteil je Platz": {str(k): r(v, 5) for k, v in lam_max.items()},
                     "Summe der Bogenlängen": str(summe), "Gitterprobe Platz 0, 1": [r(gram_num[0], 5), r(gram_num[1], 5)]}
check("[E] Plätze der Kette sind Bögen EINER Naht: Die relativen Kommutanten N_(n+1)' cap N_n sind auf Einteilchenebene "
      "L^2(n, n+1) in y = tan(theta/2), die Bogenlängen ell_n = 2 (arctan(n+1) - arctan n) summieren sich zu "
      f"{summe} (Teleskop), d. h. alle Plätze zusammen sind eine Naht ohne den Punkt t*. Jede Nahtmode e_x liegt mit dem "
      f"Gewicht ell_n/2 pi auf Platz n (Platz 0 und -1: {anteile[0]:.2f}, Platz 1: {anteile[1]:.4f}), unabhängig von x. Das "
      "Dublett der Naht (zwei Moden im Abstand 1) hat auf dem grössten Platz die Gram-Eigenwerte 1/4 +- sqrt(2)/(2 pi) = "
      f"{float(lam0[0]):.4f}, {float(lam0[1]):.4f} (Gitterprobe {gram_num[0]:.4f}; Platz 1: {lam_max[1]:.4f} gegen {gram_num[1]:.4f} "
      "mit gerundetem Endpunkt); kein Platz trägt auch nur eine der zwei "
      "Dublettmoden zur Hälfte, die ganze Kette teilt sich ein Dublett (Spur 2)",
      tele and summe == 2*sp.pi and gram_ok and max(lam_max.values()) < 0.5 and abs(lam_max[0] - float(lam0[0])) < 1e-12
      and abs(gram_num[0] - float(lam0[0])) < 1e-5 and abs(gram_num[1] - lam_max[1]) < 1e-3
      and abs(sum(anteile[k] for k in (-1, 0)) - 0.5) < 1e-12)

say("  (b) Inhalt eines Platzes gegen Inhalt einer Naht (Gitter, Doppelgenauigkeit)")


def G_vak(M):
    xs = np.arange(M) - M//2 + 0.5
    th = 2*np.pi*(np.arange(M) + 0.5)/M
    U = np.exp(1j*np.outer(th, xs))/np.sqrt(M)
    return (U*(xs < 0)) @ U.conj().T


def entropie(G, sites):
    nu = np.linalg.eigvalsh(G[np.ix_(sites, sites)])
    nu = nu[(nu > 1e-15) & (nu < 1 - 1e-15)]
    return float(-np.sum(nu*np.log(nu) + (1 - nu)*np.log(1 - nu)))


def niveaus(G, sites, hc=20.0):
    nu = np.clip(np.linalg.eigvalsh(G[np.ix_(sites, sites)]), 1e-300, 1 - 1e-16)
    return int(np.sum(np.abs(np.log((1 - nu)/nu)) < hc))


inhalt = {}
for M in (128, 256, 512, 1024):
    Gv = G_vak(M)
    sites = np.arange(0, M//4)
    inhalt[M] = (len(sites), entropie(Gv, sites), niveaus(Gv, sites))
steig = [(inhalt[2*M][1] - inhalt[M][1])/np.log(2) for M in (128, 256, 512)]


def S_original(M, alpha=0.25, beta=4*np.pi):
    xs = np.arange(M) - M//2 + alpha
    c = 1/(1 + np.exp(np.clip(beta*xs, -700, 700)))
    c = c[(c > 1e-300) & (c < 1 - 1e-17)]
    return float(-np.sum(c*np.log(c) + (1 - c)*np.log1p(-c)))


S_or = [S_original(M) for M in (128, 1024)]
ZAHLEN["inhalt"] = {"Platz 0 (Viertelbogen): M -> (Gitterpunkte, Entropie, Modulniveaus |h| < 20)":
                    {str(M): [inhalt[M][0], r(inhalt[M][1], 6), inhalt[M][2]] for M in (128, 256, 512, 1024)},
                    "Steigung dS/d ln M": [r(x, 6) for x in steig], "Entropie des Originals (alpha = 1/4, beta = 4 pi), "
                    "M = 128, 1024": [r(x, 12) for x in S_or], "ln(2)/4": r(np.log(2)/4, 12)}
check("[N] Ein Platz ist unendlich, eine Naht endlich: Die Vakuumentropie des Platzes 0 wächst logarithmisch mit dem Gitter "
      f"(S = {inhalt[128][1]:.3f}, {inhalt[256][1]:.3f}, {inhalt[512][1]:.3f}, {inhalt[1024][1]:.3f} für M = 128..1024, Steigung "
      f"{steig[-1]:.4f} je ln M; das Gitter hat zwei Fermipunkte, das chirale Kontinuum gäbe 1/6), die Zahl der Modulniveaus "
      f"|h| < 20 ebenso ({inhalt[128][2]}, {inhalt[256][2]}, {inhalt[512][2]}, {inhalt[1024][2]}: dichtes Spektrum, hyperbolisch, "
      "beta = 2 pi). Die Naht im Originalzustand hat die gitterunabhängige Entropie "
      f"{S_or[0]:.6f} je Quarter (= ln 2/4, z68; Abweichung {abs(S_or[0] - np.log(2)/4):.0e}) und die periodische Uhr "
      "4 pi (Z + alpha). Platz (Bogen) und Naht (Kreis) sind verschiedene Objekte; nach den bekannten Sätzen Typ III_1 gegen "
      "Typ I (nicht gerechnet, gerechnet sind Entropiewachstum und Niveaudichte)",
      all(abs(s - 1/3) < 2e-3 for s in steig) and inhalt[128][2] < inhalt[256][2] < inhalt[512][2] < inhalt[1024][2]
      and abs(S_or[0] - S_or[1]) < 1e-13 and abs(S_or[0] - np.log(2)/4) < 1e-12 and S_or[0] < 0.2 < inhalt[128][1])

# Kovarianz der Kette: gegenseitige Information hängt nur vom Doppelverhältnis ab (Plätze sind gleiche Kopien)
Mk = 1024
Gk = G_vak(Mk)


def grenze(n):
    return int(round(Mk*2*np.arctan(n)/(2*np.pi)))


def platz(n):
    return np.arange(grenze(n), grenze(n + 1)) % Mk


def sehne(a, b):
    return np.sin(np.pi*abs(a - b)/Mk)


def doppelv(n, m):
    a1, b1, a2, b2 = grenze(n), grenze(n + 1), grenze(m), grenze(m + 1)
    return sehne(a1, b1)*sehne(a2, b2)/(sehne(a1, a2)*sehne(b1, b2))


def gegens_info(n, m):
    A, B = platz(n), platz(m)
    return entropie(Gk, A) + entropie(Gk, B) - entropie(Gk, np.concatenate([A, B]))


PAARE = [(-3, -1), (-2, 0), (-1, 1), (0, 2), (1, 3), (0, 3), (-1, 2), (0, 4)]
mi = {p: (gegens_info(*p), doppelv(*p)) for p in PAARE}
verh = {p: v[0]/(-np.log(1 - v[1])) for p, v in mi.items()}
ZAHLEN["gegenseitige_information"] = {str(p): {"I": r(v[0], 5), "eta (Gitter)": r(v[1], 5), "I/(-ln(1 - eta))": r(verh[p], 5)}
                                      for p, v in mi.items()}
i_nn = [mi[p][0] for p in PAARE[:5]]
check("[N] Translationskovarianz und Marke der Plätze (M = 1024): Die gegenseitige Information zweier Plätze hängt nur vom "
      "Doppelverhältnis ab, I = (1/3) ln(1/(1 - eta)) (Verhältnis "
      f"{min(verh.values()):.4f} bis {max(verh.values()):.4f} für acht Paare): Plätze im Abstand 2 haben stets eta = 1/4 "
      f"(I = {min(i_nn):.4f} bis {max(i_nn):.4f}, Streuung nur durch Gitterrundung der Endpunkte), im Abstand 3 eta = 1/9 "
      f"(I = {mi[(0, 3)][0]:.4f}), im Abstand 4 eta = 1/16 (I = {mi[(0, 4)][0]:.4f}). Die Plätze sind innerlich gleiche Kopien; "
      "unterscheidbar sind sie nur durch ihre Lage in der Kette (die Korrelation fällt mit dem Abstand: nicht austauschbar)",
      max(abs(v - 1/3) for v in verh.values()) < 2e-3 and max(i_nn) - min(i_nn) < 2e-3
      and min(i_nn) > 2*mi[(0, 3)][0] > 3*mi[(0, 4)][0] > 0)

say("  (c) kein endlicher Platz: Borchers-Relation in endlicher Dimension")
nogo = []
for d in (2, 4, 16):
    for _ in range(3):
        Ah = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
        Bh = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
        Kh, Ph = Ah + Ah.conj().T, Bh + Bh.conj().T
        Df = komm(Kh, Ph) - 2j*np.pi*Ph
        nogo.append((nrm(Df)/(2*np.pi*nrm(Ph)), abs(np.trace(Ph @ Df) + 2j*np.pi*np.trace(Ph @ Ph))/nrm(Ph)**2))
Nt = 40
L0t, H1t, H2t = ds_num(0.5, Nt)
Kt, Pt = 2*np.pi*H1t, L0t - H2t
Dt = komm(Kt, Pt) - 2j*np.pi*Pt
innen = float(np.abs(Dt[:Nt - 2, :Nt - 2]).max())
rand = nrm(Dt)/(2*np.pi*nrm(Pt))
ZAHLEN["endliche_dimension"] = {"min ||[K,P] - 2 pi i P||_F/(2 pi ||P||_F) über Zufallspaare": r(min(x[0] for x in nogo), 4),
                                "Spuridentität Tr(P D) = -2 pi i Tr P^2 (Fehler)": r(max(x[1] for x in nogo), 2),
                                "Trunkierung D_1/2, N = 40: Defekt innen": r(innen, 2), "Defekt gesamt (relativ)": r(rand, 4)}
check("[E] Kein endlicher Platz: Für hermitesche K, P in endlicher Dimension ist Tr(P ([K, P] - 2 pi i P)) = -2 pi i Tr P^2 "
      f"(Fehler {max(x[1] for x in nogo):.0e}), also ||[K, P] - 2 pi i P||_F >= 2 pi ||P||_F (kleinstes Verhältnis über neun "
      f"Zufallspaare in d = 2, 4, 16: {min(x[0] for x in nogo):.2f} >= 1): die Borchers-Relation erzwingt P = 0. Auf der "
      f"Trunkierung von D_1/2 (N = 40) gilt sie im Inneren exakt ({innen:.0e}) und der ganze Defekt {rand:.2f} sitzt am "
      "Abschneiderand. Eine Wiesbrock-Translation lebt also nie auf endlich vielen endlichen Plätzen wie (C^2)^L; die "
      "z121-Kette kann höchstens eine Diskretisierung sein",
      min(x[0] for x in nogo) >= 1 - 1e-12 and max(x[1] for x in nogo) < 1e-10 and innen < 1e-10 and rand >= 1 - 1e-12)

say("  (d) Schrittweite: frei, Einheit der Marken, Vierteldrehung")


def mat_T(s):
    return sp.Matrix([[1, s], [0, 1]])


R_q = sp.Matrix([[1, 1], [-1, 1]])/sp.sqrt(2)                       # Vierteldrehung der Naht in y: y -> (y + 1)/(1 - y)
S_h = sp.simplify(R_q*R_q)                                          # Halbdrehung y -> -1/y


def mwirk(A, y):
    if y is sp.oo:
        return sp.oo if A[1, 0] == 0 else sp.simplify(A[0, 0]/A[1, 0])
    nen = sp.simplify(A[1, 0]*y + A[1, 1])
    return sp.oo if nen == 0 else sp.simplify((A[0, 0]*y + A[0, 1])/nen)


marken_zyklus = [mwirk(R_q, y0) for y0 in (sp.Integer(0), sp.Integer(1), sp.oo, sp.Integer(-1))]     # 0 -> 1 -> oo -> -1 -> 0
lam_s = sp.symbols("lambda", positive=True)
Dil = sp.Matrix([[sp.sqrt(lam_s), 0], [0, 1/sp.sqrt(lam_s)]])       # Modulfluss von M = (0, oo): y -> lambda y
skalierung = sp.simplify(Dil*mat_T(1)*Dil.inv() - mat_T(lam_s)) == sp.zeros(2, 2)
tr_RT1 = sp.simplify((R_q*mat_T(1)).trace())
minpol = sp.Poly(sp.minimal_polynomial(tr_RT1, x_s), x_s)
ganz = minpol.LC() == 1                                            # 2 cos(pi p/q) ist stets ganzalgebraisch
tr_RT2 = sp.simplify((R_q*mat_T(2)).trace())
ord_RT2 = sp.simplify((R_q*mat_T(2))**2 + sp.eye(2)) == sp.zeros(2, 2)
ord_ST1 = sp.simplify((S_h*mat_T(1))**3 - sp.eye(2)) == sp.zeros(2, 2) or sp.simplify((S_h*mat_T(1))**3 + sp.eye(2)) == sp.zeros(2, 2)
RT1 = np.array((R_q*mat_T(1)).evalf(30), dtype=float)
Wk = np.eye(2)
naehe = []
for k in range(1, 2001):
    Wk = Wk @ RT1
    naehe.append(min(nrm(Wk - np.eye(2)), nrm(Wk + np.eye(2))))
kette_bild = [mwirk(R_q, sp.Integer(k)) for k in (2, 3, 4)]         # Bild der Kettenendpunkte unter der Vierteldrehung
vier_summe = null(sum(1 - sp.cos(TH - k*sp.pi/2) for k in range(4)) - 4)
endpunkt2 = sp.simplify(2*sp.atan(2)/(sp.pi/2))
ZAHLEN["schritt"] = {"Marken unter der Vierteldrehung (y)": [str(m) for m in marken_zyklus],
                     "Spur von R T_1": str(tr_RT1), "Minimalpolynom": str(minpol.as_expr()), "Spur von R T_2": str(tr_RT2),
                     "min Abstand (R T_1)^k zu +-1, k <= 2000": r(min(naehe), 3),
                     "Bild der Endpunkte 2, 3, 4 unter der Vierteldrehung": [str(v) for v in kette_bild],
                     "Endpunkt y = 2 in Einheiten der Marke": r(float(endpunkt2), 6)}
check("[E] Schrittweite: Der Modulfluss von M skaliert die Translation, Dil(lambda) T_1 Dil(lambda)^(-1) = T_lambda (sympy): "
      "alle Schrittweiten sind unitär äquivalent, die Inklusion allein kennt keinen Schritt (z72: keine eigene Skala). Die "
      "Marken geben eine Einheit [C, z83/z104]: N = Viertelbogen zwischen Marken, M = Halbbogen, Faktor 1, also T_1 mit den "
      f"Marken y = 0, 1, oo, -1 (Vierteldrehung: 0 -> {' -> '.join(str(m) for m in marken_zyklus)}). Aber schon der zweite "
      f"Kettenendpunkt y = 2 liegt bei {float(endpunkt2):.4f} Viertelkreisen, keiner Marke, und die Vierteldrehung bildet die "
      f"Kette nicht auf sich ab (2, 3, 4 -> {', '.join(str(v) for v in kette_bild)})",
      skalierung and marken_zyklus == [1, sp.oo, -1, 0] and abs(float(endpunkt2) - round(float(endpunkt2))) > 0.1
      and kette_bild[2] != sp.floor(kette_bild[2]))
check(f"[E] Uhr und Kette sind nicht gemeinsam diskret: R T_1 (Vierteldrehung mal Einheitstranslation) hat die Spur {tr_RT1} "
      f"mit Minimalpolynom {minpol.as_expr()} (nicht normiert, also nicht ganzalgebraisch, also kein 2 cos(pi p/q)): ein "
      "elliptisches Element unendlicher Ordnung, die erzeugte Gruppe ist nicht diskret (numerisch: (R T_1)^k kommt +-1 bis auf "
      f"{min(naehe):.3f} nahe, k <= 2000, ohne es zu erreichen). Kontrollen: mit der Halbdrehung ist S T_1 von Ordnung 3 "
      f"(Modulgruppe), mit dem Doppelschritt ist R T_2 von Ordnung 2 (Spur {tr_RT2}; Dreiecksgruppe (2, 4, oo)). Die Uhr ist das "
      "Mittel der vier nativen Translationen, L_0 = (1/4) sum_k P^(k pi/2): keine einzelne ist uhrinvariant, die Wahl von t* "
      "bricht Z4",
      abs(tr_RT1) < 2 and not ganz and 0 < min(naehe) < 0.05 and ord_RT2 and tr_RT2 == 0 and ord_ST1 and vier_summe)
befund("E: 'Geordnete Nähte' ist nicht die Kette der Inklusion. Die Kette ordnet Bögen einer einzigen Naht (in y = t: "
       "Einheitsintervalle der Minkowskizeit auf der Achse, Erzeuger P_0 = E + Q^-). Jeder Platz ist ein Bogen: unendlich "
       "viele Moden, Typ III, hyperbolischer Modulfluss bei 2 pi, kein Dublett; eine Naht ist ein Kreis: Typ I, Uhr 4 pi, "
       "endliche Entropie. Ein endlicher Platz ist ausgeschlossen (Borchers-Relation)", True)
befund("E/C: Die Beschriftung ist nicht kanonisch. Die Translation ist kontinuierlich und skalenfrei; die Marken geben den "
       "ersten Schritt (Viertelbogen in Halbbogen), danach liegen die Kettenendpunkte auf keiner Marke. Nativ diskret und "
       "geordnet ist nur der Ring der vier Viertelbögen unter der Uhr (zyklisch Z4), keine Kette Z; Uhr und "
       "Einheitstranslation erzeugen zusammen eine nichtdiskrete Gruppe", True)

# =================================================================== Teil 3
say("\nTeil 3: Lebt das Gatter auf der Kette?")
say("  (a) Lichtkegelstruktur der z121-Kette")


def swap_idx(L, x, y):
    idx = np.arange(2**L)
    bx, by = (idx >> x) & 1, (idx >> y) & 1
    return np.where(bx != by, idx ^ ((1 << x) | (1 << y)), idx)


def lage(V, L, u, start):
    """Eine Lage Gatter R(u) = (1 + i u P)/(1 + i u) auf den Bindungen (x, x+1), x = start, start + 2, ...; u = oo: Tausch."""
    a, b = (1/(1 + 1j*u), 1j*u/(1 + 1j*u)) if np.isfinite(u) else (0.0, 1.0)
    for x in range(start, L, 2):
        V = a*V + b*V[swap_idx(L, x, (x + 1) % L)]
    return V


def F_mat(L, u=1.0):
    return lage(lage(np.eye(2**L, dtype=complex), L, u, 0), L, u, 1)


def t_mat(L, u, xi):
    """Inhomogene Transfermatrix t(u) = Tr_0 prod_j (u - xi_j + i P_(0j)) als Matrix (wie z121)."""
    D = 2**L
    idx = np.arange(D)
    V = np.eye(D, dtype=complex)
    out = np.zeros_like(V)
    for a0 in (0, 1):
        X = [np.zeros_like(V), np.zeros_like(V)]
        X[a0] = V.copy()
        for j in range(L):
            w = u - xi[j]
            bj = ((idx >> j) & 1)[:, None]
            PX = []
            for a in (0, 1):
                src = idx ^ ((((idx >> j) & 1) ^ a) << j)
                PX.append(np.where(bj == 0, X[0][src], X[1][src]))
            X = [w*X[a] + 1j*PX[a] for a in (0, 1)]
        out += X[a0]
    return out


def verschiebung(L, um, nur=None):
    """Permutationsmatrix: Plätze (nur der Untergitterparität 'nur', sonst alle) wandern um 'um' Plätze."""
    D = 2**L
    idx = np.arange(D)
    img = np.zeros(D, dtype=int)
    for x in range(L):
        y = (x + um) % L if (nur is None or x % 2 == nur) else x
        img |= ((idx >> x) & 1) << y
    Pm = np.zeros((D, D))
    Pm[img, idx] = 1
    return Pm


def bis_auf_phase(A, B):
    ph = np.vdot(B.ravel(), A.ravel())
    ph = ph/abs(ph) if abs(ph) > 1e-14 else 1.0
    return nrm(A - ph*B)


def spurabstand(A, B, K=6):
    """max_k | |Tr A^k| - |Tr B^k| |/dim: null, wenn die Spektren (bis auf eine Phase) gleich sind."""
    d = A.shape[0]
    Ak, Bk, out = np.eye(d), np.eye(d), []
    for _ in range(K):
        Ak, Bk = Ak @ A, Bk @ B
        out.append(abs(abs(np.trace(Ak)) - abs(np.trace(Bk)))/d)
    return float(max(out))


LK3, FREI = {}, {}
for L in (4, 6, 8):
    D = 2**L
    Id = np.eye(D, dtype=complex)
    xi = np.array([-0.5 if j % 2 == 0 else 0.5 for j in range(L)])
    F1, Fs = F_mat(L, 1.0), F_mat(L, np.inf)
    Ue, Uo = lage(Id, L, 1.0, 0), lage(Id, L, 1.0, 1)
    UR = t_mat(L, -0.5, xi)/2**(L/4)
    UL = np.linalg.inv(t_mat(L, 0.5, xi)/2**(L/4))
    T1, T2 = verschiebung(L, 1), verschiebung(L, 2)
    LK3[L] = {"unitaer": max(nrm(UR.conj().T @ UR - Id), nrm(UL.conj().T @ UL - Id)), "[U_R, U_L]": nrm(komm(UR, UL)),
              "F = U_R U_L": bis_auf_phase(UR @ UL, F1), "T^2 = U_R U_L^-1": bis_auf_phase(UR @ np.linalg.inv(UL), T2),
              "[F, T^2]": nrm(komm(F1, T2)), "[F, T]": nrm(komm(F1, T1)), "[U_gerade, U_ungerade]": nrm(komm(Ue, Uo)),
              "U_gerade gegen U_R, U_L": min(bis_auf_phase(Ue, UR), bis_auf_phase(Ue, UL), bis_auf_phase(Uo, UR),
                                              bis_auf_phase(Uo, UL))}
    WR, WL = verschiebung(L, 2, nur=0), verschiebung(L, -2, nur=1)
    FREI[L] = {"F_Tausch = W_R W_L": nrm(WR @ WL - Fs), "T^2 = W_R W_L^-1": nrm(WR @ WL.T - T2), "[W_R, W_L]": nrm(komm(WR, WL)),
               "F_Tausch^(L/2) = 1": nrm(np.linalg.matrix_power(Fs, L//2) - Id),
               "Spurabstand F(u) zu F_Tausch": {str(u): spurabstand(F_mat(L, u), Fs) for u in (1.0, 16.0, 256.0)},
               "F(1)^(L/2) gegen Phase": min(nrm(np.linalg.matrix_power(F1, L//2) - np.exp(2j*np.pi*k/720)*Id) for k in range(720))}
ZAHLEN["lichtkegel"] = {str(L): {k: r(v, 3) for k, v in LK3[L].items()} for L in (4, 6, 8)}
ZAHLEN["frei"] = {str(L): {k: ({kk: r(vv, 4) for kk, vv in v.items()} if isinstance(v, dict) else r(v, 4))
                           for k, v in FREI[L].items()} for L in (4, 6, 8)}
mx = lambda key: max(LK3[L][key] for L in (4, 6, 8))
mn = lambda key: min(LK3[L][key] for L in (4, 6, 8))
check("[E] Die z121-Kette trägt die Lichtkegelstruktur (L = 4, 6, 8, volle Matrizen): U_R = t(-1/2)/2^(L/4) und U_L = "
      f"(t(1/2)/2^(L/4))^(-1) sind unitär ({mx('unitaer'):.0e}), vertauschen ({mx('[U_R, U_L]'):.0e}), und F = U_R U_L, "
      f"T^2 = U_R U_L^(-1) bis auf Phase ({mx('F = U_R U_L'):.0e}, {mx('T^2 = U_R U_L^-1'):.0e}): zwei vertauschende "
      f"Translationen, Zeit- und Raumschritt als Produkt und Quotient. F vertauscht mit T^2 ({mx('[F, T^2]'):.0e}), nicht mit "
      f"dem Einplatzschritt T ({mn('[F, T]'):.2f}): die Einheitszelle hat zwei Plätze (Untergitter = zwei Bewegungsrichtungen). "
      f"Die Ziegellagen selbst sind nicht das +-Paar: sie vertauschen nicht ({mn('[U_gerade, U_ungerade]'):.2f}) und sind "
      f"keine der beiden Translationen (Abstand {mn('U_gerade gegen U_R, U_L'):.2f})",
      mx("unitaer") < 1e-12 and mx("[U_R, U_L]") < 1e-12 and mx("F = U_R U_L") < 1e-12 and mx("T^2 = U_R U_L^-1") < 1e-12
      and mx("[F, T^2]") < 1e-12 and mn("[F, T]") > 0.5 and mn("[U_gerade, U_ungerade]") > 0.5
      and mn("U_gerade gegen U_R, U_L") > 0.5)
fx = lambda key: max(FREI[L][key] for L in (4, 6, 8))
sa = {u: [FREI[L]["Spurabstand F(u) zu F_Tausch"][u] for L in (4, 6, 8)] for u in ("1.0", "16.0", "256.0")}
check("[E] Lichtkegelmatrix aus den nativen Translationen: Native Translationen sind quasifrei (zweite Quantisierung eines "
      "Einteilchenflusses); auf zwei Untergittern (P_chir, Schritt 2) geben sie W_R (gerade Plätze um +2) und W_L (ungerade um "
      f"-2), und W_R W_L ist exakt die Tauschkette F(u = oo) ({fx('F_Tausch = W_R W_L'):.0e}), T^2 = W_R W_L^(-1) "
      f"({fx('T^2 = W_R W_L^-1'):.0e}), F_Tausch^(L/2) = 1 ({fx('F_Tausch^(L/2) = 1'):.0e}). Mit z121s Gatter u = 1 stimmt das "
      f"Spektrum nicht: Spurabstand max_k ||Tr F^k| - |Tr F_Tausch^k||/dim = {sa['1.0'][0]:.2f}, {sa['1.0'][1]:.2f}, "
      f"{sa['1.0'][2]:.2f} (L = 4, 6, 8), F(1)^(L/2) ist keine Phase ({min(FREI[L]['F(1)^(L/2) gegen Phase'] for L in (4, 6, 8)):.1f}); "
      f"erst für u -> oo schliesst sich der Abstand ({sa['16.0'][2]:.3f} bei u = 16, {sa['256.0'][2]:.4f} bei u = 256, L = 8). Die "
      "Translationen liefern die Geometrie (freie Rechts- und Linksläufer), nicht das Gatter",
      fx("F_Tausch = W_R W_L") < 1e-12 and fx("T^2 = W_R W_L^-1") < 1e-12 and fx("[W_R, W_L]") < 1e-12
      and fx("F_Tausch^(L/2) = 1") < 1e-12 and min(sa["1.0"]) > 0.3
      and all(sa["1.0"][i] > sa["16.0"][i] > sa["256.0"][i] for i in range(3)) and max(sa["256.0"]) < 5e-3
      and min(FREI[L]["F(1)^(L/2) gegen Phase"] for L in (4, 6, 8)) > 1.0)

say("  (b) Gatter gegen die native Translation zweier Nähte (Fockfenster wie z121)")


def jw(n):
    I_, Zp, sm = np.eye(2), np.diag([1.0, -1.0]), np.array([[0.0, 1.0], [0.0, 0.0]])
    out = []
    for j in range(n):
        m_ = [Zp]*j + [sm] + [I_]*(n - j - 1)
        op = m_[0]
        for Mx in m_[1:]:
            op = np.kron(op, Mx)
        out.append(op)
    return out


XS = [-1.75, -0.75, 0.25, 1.25]                                     # Moden einer nativen Naht (alpha = 1/4), wie z121
c8 = jw(8)
n8 = [x.T @ x for x in c8]
vak = np.array([1 if x < 0 else 0 for x in XS]*2)
Q8 = [sum(n8[4*i + k] - vak[4*i + k]*np.eye(256) for k in range(4)) for i in range(2)]
E_naht = [sum((XS[k] + 0.25)*(n8[4*i + k] - vak[4*i + k]*np.eye(256)) for k in range(4)) for i in range(2)]
Lm_naht = [sum((XS[k] + 0.75)*c8[4*i + k + 1].T @ c8[4*i + k] for k in range(3)) for i in range(2)]     # L^nat_-1 je Naht
P_naht = [E_naht[i] + (Lm_naht[i] + Lm_naht[i].T)/2 for i in range(2)]                                # dGamma(P_0) je Naht
J0 = sum(c8[k].T @ c8[4 + k] for k in range(4))
B8 = expm(-1j*np.pi/4*(J0 + J0.T))
U8 = B8.conj().T @ expm(-1j*np.pi*(Q8[0] @ Q8[0] + Q8[1] @ Q8[1])/4) @ B8
Jf = sum(c8[k].T @ c8[4 + k] for k in range(4))
Tausch8 = expm(-1j*np.pi/2*(Jf + Jf.T))                             # kanonischer Tausch der zwei Nähte (bis auf Zahlphase)
gat = {"[U, P_A + P_B]": nrm(komm(U8, P_naht[0] + P_naht[1])), "[U, P_A - P_B]": nrm(komm(U8, P_naht[0] - P_naht[1])),
       "[U, P_A]": nrm(komm(U8, P_naht[0])), "U nichttrivial": bis_auf_phase(U8, np.eye(256)),
       "Tausch: P_A -> P_B": nrm(Tausch8 @ P_naht[0] @ Tausch8.conj().T - P_naht[1]),
       "[e^(-i P_A), e^(-i P_B)]": nrm(komm(expm(-1j*P_naht[0]), expm(-1j*P_naht[1])))}
ZAHLEN["gatter_translation"] = {k: r(v, 3) for k, v in gat.items()}
check("[E] Das Gatter ist ein Skalar der nahtinternen Translation (zwei native Nähte, Moden -7/4..5/4, Fock 256; Gatter aus "
      "den Nullmoden wie z121): Es vertauscht mit der gemeinsamen Translation dGamma(P_0) beider Nähte "
      f"({gat['[U, P_A + P_B]']:.0e}), nicht mit der relativen ({gat['[U, P_A - P_B]']:.2f}) und nicht mit der einer Naht "
      f"({gat['[U, P_A]']:.2f}), und ist nicht trivial ({gat['U nichttrivial']:.1f}). Die internen Translationen zweier Nähte "
      f"vertauschen ({gat['[e^(-i P_A), e^(-i P_B)]']:.0e}) und werden vom Nahttausch aufeinander abgebildet "
      f"({gat['Tausch: P_A -> P_B']:.0e}). Die halbseitige Translation verschiebt also innerhalb jeder Naht und lässt den "
      "Nahtindex fest: Sie steht quer zur Kette der Nähte, das Gatter verbindet Nähte bei gleicher interner Lage",
      gat["[U, P_A + P_B]"] < 1e-11 and gat["[U, P_A - P_B]"] > 0.1 and gat["[U, P_A]"] > 0.1 and gat["U nichttrivial"] > 0.1
      and gat["Tausch: P_A -> P_B"] < 1e-11 and gat["[e^(-i P_A), e^(-i P_B)]"] < 1e-11)
befund("E: Das Gatter lebt nicht auf der Kette der Inklusion. Die z121-Kette hat zwar genau die Form zweier vertauschender "
       "Translationen (F = U_R U_L, T^2 = U_R U_L^(-1)), aber (i) das +-Paar sind die diagonalen Transfermatrizen, nicht die "
       "Ziegellagen, und es braucht zwei Untergitter (P_chir); (ii) die nativen Translationen sind quasifrei und geben die "
       "Tauschkette u = oo, das Gatter u = 1 bleibt P_null/P_viertel; (iii) die native Translation wirkt in jeder Naht und "
       "vertauscht mit dem Gatter nur gemeinsam; (iv) endliche Plätze tragen keine Wiesbrock-Translation", True)
befund("E/N: Unterscheidbarkeit. In einer Wiesbrock-Kette ist die Marke die Lage selbst: Die Plätze sind innerlich gleich, "
       "aber durch die Inklusion total geordnet (N_(n+1) in N_n), die Richtung ist das Vorzeichen der Halbseitigkeit "
       "(Positivität von P), und der Zustand ist nicht austauschbar (I fällt mit dem Abstand). Das wäre genau die Marke gegen "
       "den Twirl von z119, aber für Plätze (Bögen einer Naht), nicht für Nähte: Die nativen Nahtdaten tragen sie nicht "
       "(interne Translationen verschiedener Nähte vertauschen und sind unter dem Tausch gleich)", True)

# =================================================================== Teil 4
say("\nTeil 4: Bilanz und Einstufung")
BILANZ = [
    ("E", "Halbseitige modulare Inklusionen beider Vorzeichen existieren nativ im Vakuumsektor: O_1, V_+ (z104) mit "
          "N_+- = Viertelbögen zwischen Marken; Differenzen = parabolische Erzeuger in den Marken, Faktor 1, positiv "
          "(Gamma-Mass)"),
    ("N", "dasselbe aus dem Gitterzustand (640 Bit): Hislop-Longo, Wiesbrock, Borchers, Leckage einseitig; Kontrolle "
          "(geschachtelt ohne gemeinsamen Endpunkt) scheitert; Originalzustand beta = 4 pi: keine halbseitige Inklusion"),
    ("E", "Uhr = Moduldefekt: K_O1 + K_V+ - K_(O1 cap V+) = 2 pi L_0; L_0 = Mittel der vier Translationen in den Marken"),
    ("E", "die zwei Vorzeichen einer Naht: [P, P'] = -i K_gemeinsam, zusammen sl(2,R) (ein Lichtstrahl), keine "
          "Lichtkegelkoordinaten"),
    ("C", "zwei vertauschende Lichttranslationen: P_chir (zwei Arten, entgegengesetzte Enden) gibt die 1+1-Poincaré-Algebra"),
    ("E/C", "su(2,2) (z83): die halbseitige Familie einer Spitze ist R^4 mit Vorwärtskegel-Positivität (3+1-Kinematik); "
            "so(3)-invariant ist nur die Zeit; Raum- und Lichtrichtung verlangen eine Achse (P_richt)"),
    ("E", "Plätze der Kette = Bögen einer Naht (Summe 2 pi); kein Platz trägt das Dublett; Platz unendlich (Typ III), Naht "
          "endlich (Typ I); endliche Plätze ausgeschlossen"),
    ("E/C", "Schritt frei (skalenkovariant); Einheit aus den Marken nur für den ersten Schritt; Uhr und Einheitstranslation "
            "erzeugen eine nichtdiskrete Gruppe; nativ diskret ist nur der Z4-Ring der Viertelbögen"),
    ("E", "z121-Kette: F = U_R U_L, T^2 = U_R U_L^(-1), [U_R, U_L] = 0; Ziegellagen sind nicht das +-Paar"),
    ("E", "native (quasifreie) Translationen geben die Tauschkette u = oo, nicht das Gatter u = 1; das Gatter vertauscht "
          "nur mit der gemeinsamen internen Translation beider Nähte"),
    ("C", "P_ord folgt (Ordnung, Richtung, Translation, Marke) aus einer Prämisse P_platz: 'die Nähte sind die Plätze "
          "einer gemeinsamen Wiesbrock-Kette' (zweite Naht = halbseitig modular eingeschlossenes Bild der ersten)"),
    ("O", "P_platz selbst: nativ sind die Plätze Bögen, keine Nähte (Typ, Dublett, Uhr); dazu P_chir, die Wahl von t* "
          "(bricht Z4), ein diskreter Schritt und weiter P_null/P_viertel für das Gatter"),
    ("O", "3+1 und Kontinuum: die Route gibt Kinematik (Translationen mit Spektralbedingung) nur im Vakuumsektor und unter "
          "z83; keine native Raumrichtung, kein Gitter von Nähten, kein Kontinuumslimes des Gatters"),
]
for st, txt in BILANZ:
    say(f"   [{st}] {txt}")
ZAHLEN["bilanz"] = [{"status": st, "aussage": txt} for st, txt in BILANZ]
EINSTUFUNG = ("(b) bedingt: P_ord ist keine native Folge der Modulstruktur. Es reduziert sich auf die Prämisse P_platz (Nähte "
              "= Plätze einer gemeinsamen Wiesbrock-Kette, d. h. eine zweite Naht als halbseitig modular eingeschlossene "
              "Unteralgebra) samt P_chir (zwei Lichtstrahlen), der Wahl von t* und eines diskreten Schritts. P_platz ist nicht "
              "kleiner als P_ord im Sinn von schwächer, sondern strukturell: eine Inklusion statt L markierter Nähte; die "
              "native Wiesbrock-Kette erfüllt P_platz nicht (ihre Plätze sind Bögen einer Naht).")
ZAHLEN["einstufung_P_ord"] = EINSTUFUNG
say(f"   Einstufung: {EINSTUFUNG}")
befund("Einstufung P_ord: (b) bedingt auf P_platz (+ P_chir, t*, Schritt); nicht (a), weil die native Kette Bögen einer Naht "
       "ordnet und im Originalzustand gar nicht existiert; nicht als (c) gewertet, weil eine einzige Inklusion Ordnung, "
       "Richtung, Translation und Marke zugleich lieferte", True)

urteil = ("P_ord ist keine Folge der Modultheorie, sondern bedingt (b). Nativ ist viel: Im Vakuumsektor der rekonstruierten "
          "Uhr bilden die Kappe O_1 und der Kegel V_+ (z104) mit den Viertelbögen N_+- zwischen den Marken halbseitige "
          "modulare Inklusionen beider Vorzeichen, exakt (Möbius) und am Gitter aus dem Zustand (Hislop-Longo "
          f"{g2['hl']:.0e}, Borchers {g2['borchers']:.0e}, Leckage {lv[0]:.0e} einseitig gegen {lv[1]:.2f}); die "
          "Wiesbrock-Differenzen sind die parabolischen Erzeuger in den Marken mit Faktor 1, positiv mit Gamma-Mass, und die "
          "Uhr ist ihr Mittel bzw. der Moduldefekt K_O1 + K_V+ - K_(O1 cap V+) = 2 pi L_0. Im Originalzustand (beta = 4 pi) "
          f"gibt es keine solche Inklusion (Leckage {lo4[0]:.0e}; z72). Aber: (1) Die zwei Vorzeichen einer Naht geben keine "
          "Lichtkegelkoordinaten; ihre Translationen vertauschen nicht, der Kommutator ist der Modulgenerator des gemeinsamen "
          "Gebiets, zusammen sl(2,R). Zwei vertauschende Lichttranslationen verlangen zwei Arten mit entgegengesetzter "
          "Zuordnung (P_chir) oder in su(2,2) ein Nullpaar zu einer gewählten Achse. (2) Die Plätze der Kette N_n sind Bögen "
          "einer einzigen Naht (Bogenlängen summieren sich zu 2 pi), in y = t Einheitsintervalle der Minkowskizeit; jeder "
          "Platz ist unendlich (Typ III, Entropie ~ ln M, dichtes Modulspektrum), keiner trägt das Dublett (grösster Anteil "
          f"{float(lam0[0]):.3f}), während die Naht endlich ist (Entropie {S_or[0]:.4f}, Uhr 4 pi); endliche Plätze sind durch "
          "die Borchers-Relation ausgeschlossen. 'Geordnete Nähte' ist also nicht die Kette der Inklusion. (3) Der Schritt ist "
          "frei (skalenkovariant); die Marken geben nur den ersten Schritt, Uhr und Einheitstranslation erzeugen eine "
          "nichtdiskrete Gruppe, nativ diskret ist allein der Z4-Ring der vier Viertelbögen. (4) Die z121-Kette hat die Form "
          "zweier vertauschender Translationen (F = U_R U_L, T^2 = U_R U_L^(-1)), doch die nativen quasifreien Translationen "
          "geben nur die Tauschkette u = oo; das Gatter u = 1 vertauscht bloss mit der gemeinsamen internen Translation beider "
          "Nähte, die halbseitige Translation steht quer zur Nahtkette. Übrig bleibt die benannte Prämisse P_platz (Nähte = "
          "Plätze einer gemeinsamen Wiesbrock-Kette; zweite Naht als halbseitig eingeschlossene Unteralgebra) samt P_chir, t* "
          "und Schritt; unter ihr folgten Ordnung, Richtung, Translation und die Marke gegen den Twirl von z119 aus einer "
          "Inklusion. Raum: Eine Naht gibt nur einen Lichtstrahl; in su(2,2) (z83, bedingt) ist die halbseitige Familie einer "
          "Spitze die volle Translationsgruppe R^4 mit Vorwärtskegel-Positivität, aber so(3)-invariant ist nur die Zeit, eine "
          "Raumrichtung verlangt eine Achsenwahl. Für '3+1' heisst das Kinematik im Vakuumsektor ohne native Richtung und ohne "
          "Nahtgitter, für 'Kontinuum' nichts Neues: Der Limes des Gatters bleibt der von z121. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z139_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
