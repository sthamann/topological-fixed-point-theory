"""Z96: Prüfung der zwölften externen Analyse 'Gravitation als Dynamik der Verbindung von Original und Gegenstruktur'
(25.11.; MacDowell-Mansouri, Plebanski, Clifford-Projektionen, Lambda-G-Relation, bimetrische Gewichte).

Teil 1  Clifford-Projektionen P_±(n) = (1 ± n.Gamma)/2 auf S^4: Cliffordrelationen, Rang, Komplementarität, Krümmung
        F_± = P dP ^ dP P (selbstdual / antiselbstdual), zweite Chernzahlen; Einordnung als quaternionische
        Hopf-/Twistorfaserung CP^3 -> S^4 (quaternionische Struktur j, j^2 = -1) und als Whitney-Summe.
Teil 2  Q_a = (P_a + sigma K_a)/2 in so(4,2): [P, P] = [K, K] = 0, [Q_a, Q_b] = ± sigma M_ab; reelle Form der zehn-
        dimensionalen Unteralgebra (so(3,2) oder so(4,1)); welche Wahl die kompakte konforme Energie enthält, welche den
        Modulgenerator eines Doppelkegels (Hislop-Longo); der Radius ist unter der Dilatation konjugiert.
Teil 3  Plebanski-Flächen aus einem reellen Vierbein: 27 Relationen, Hodge-Dualität, Gegenprobe mit komplexem Vierbein;
        Abgleich mit z29 (a = i, b = 1).
Teil 4  MacDowell-Mansouri: eps F ^ F = Euler + (2 sigma/l^2) Einstein-Hilbert + (1/l^4) Volumen (genau drei Terme),
        Lambda = -3 sigma/l^2, b = l^2/(64 pi G sigma); Euler-Normierung 32 pi^2 chi auf S^4; de-Sitter-Entropie.
Teil 5  Lambda-G-Relation aus COSMO.LAM.01: l^2/G = 32 pi^3 e^(2/alpha), |b| = (pi^2/2) e^(2/alpha); Inhalt, Vorzeichen
        (dS gegen die E-invariante AdS-Wahl von z89) und Skala.
Teil 6  Bimetrische Gewichte 1 - e^(-4 pi), e^(-4 pi) und die Hassan-Rosen-Lesart (z47, z48); v358 |mu_4| = 4;
        Flächenkoeffizient und Zentralladung (z90).
"""
import itertools
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy import integrate
from scipy.linalg import expm

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
rng = np.random.default_rng(96)


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


def hat(d, teil, feld="checks"):
    return any(teil in k and v for k, v in d.get(feld, {}).items())


EPS4 = np.zeros((4, 4, 4, 4))
for p in itertools.permutations(range(4)):
    EPS4[p] = round(np.linalg.det(np.eye(4)[list(p)]))
EPS3 = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    EPS3[p] = round(np.linalg.det(np.eye(3)[list(p)]))

# ------------------------------------------------------------------ Teil 1
say("Teil 1: Clifford-Projektionen auf S^4")
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], complex)
s0 = np.eye(2, dtype=complex)
GAM = [np.kron(sx, sx), np.kron(sx, sy), np.kron(sx, sz), np.kron(sy, s0), np.kron(sz, s0)]
I4 = np.eye(4, dtype=complex)
n_cl = sum(np.array_equal(GAM[A] @ GAM[B] + GAM[B] @ GAM[A], 2*(A == B)*I4) for A in range(5) for B in range(5))
herm = all(np.array_equal(g, g.conj().T) for g in GAM)
omega5 = GAM[0] @ GAM[1] @ GAM[2] @ GAM[3] @ GAM[4]
check(f"Cliffordrelationen {{Gamma_A, Gamma_B}} = 2 delta_AB exakt: {n_cl}/25; alle Gamma hermitesch; Volumenelement "
      f"Gamma_1 ... Gamma_5 = {int(round(omega5[0, 0].real))} (legt die Orientierung fest)",
      n_cl == 25 and herm and np.array_equal(omega5, -I4))

ok_p = True
for _ in range(50):
    v = rng.normal(size=5)
    n = v/np.linalg.norm(v)
    nG = sum(n[A]*GAM[A] for A in range(5))
    Pp, Pm = (I4 + nG)/2, (I4 - nG)/2
    ok_p &= np.allclose(Pp @ Pp, Pp) and np.allclose(Pm @ Pm, Pm) and np.allclose(Pp + Pm, I4)
    ok_p &= np.allclose(Pp @ Pm, 0) and np.allclose(Pp, Pp.conj().T)
    ok_p &= np.linalg.matrix_rank(Pp) == 2 and np.linalg.matrix_rank(Pm) == 2
check("P_±(n) = (1 ± n.Gamma)/2 sind hermitesche Projektionen vom Rang 2 mit P_+ + P_- = 1 und P_+ P_- = 0 "
      "(50 Zufallspunkte auf S^4)", ok_p)


def nd(x):
    """Stereographische Karte R^4 -> S^4 und ihre Ableitung dn[mu, A]."""
    r2 = x @ x
    d = 1 + r2
    n = np.concatenate([2*x/d, [(1 - r2)/d]])
    dn = np.zeros((4, 5))
    for mu in range(4):
        dn[mu, :4] = 2*np.eye(4)[mu]/d - 4*x*x[mu]/d**2
        dn[mu, 4] = -4*x[mu]/d**2
    return n, dn


def kruemmung(x, s):
    n, dn = nd(x)
    P = (I4 + s*sum(n[A]*GAM[A] for A in range(5)))/2
    dP = [s*sum(dn[mu, A]*GAM[A] for A in range(5))/2 for mu in range(4)]
    F = np.zeros((4, 4, 4, 4), complex)
    for mu in range(4):
        for nu in range(4):
            F[mu, nu] = P @ (dP[mu] @ dP[nu] - dP[nu] @ dP[mu]) @ P
    return F


def stern(F):
    return 0.5*np.einsum('mnrs,rsij->mnij', EPS4, F)


x0 = rng.normal(size=4)
h = 1e-6
fd = max(np.max(np.abs((nd(x0 + h*np.eye(4)[mu])[0] - nd(x0 - h*np.eye(4)[mu])[0])/(2*h) - nd(x0)[1][mu]))
         for mu in range(4))
assert fd < 1e-8 and abs(np.linalg.norm(nd(x0)[0]) - 1) < 1e-14

dual = {}
ok_d, ok_ah = True, True
for s in (+1, -1):
    sig = set()
    for _ in range(12):
        x = rng.normal(size=4)*1.3
        F = kruemmung(x, s)
        nF = np.max(np.abs(F))
        r_sd, r_asd = np.max(np.abs(F - stern(F))), np.max(np.abs(F + stern(F)))
        if r_sd < 1e-12*nF and r_asd > 0.1*nF:
            sig.add(+1)
        elif r_asd < 1e-12*nF and r_sd > 0.1*nF:
            sig.add(-1)
        else:
            sig.add(0)
        ok_ah &= all(np.allclose(F[m, k] + F[m, k].conj().T, 0) for m in range(4) for k in range(4))
        ok_ah &= max(abs(np.trace(F[m, k])) for m in range(4) for k in range(4)) < 1e-12
    ok_d &= len(sig) == 1 and 0 not in sig
    dual[s] = sig.pop()
check(f"Krümmung F_± = P dP ^ dP P (12 Zufallspunkte je Projektion, je drei 4x4-Komponentengleichungen): F_+ ist "
      f"{'selbstdual' if dual[1] == 1 else 'antiselbstdual'}, F_- {'selbstdual' if dual[-1] == 1 else 'antiselbstdual'} "
      "(die jeweils andere Dualität verletzt um O(|F|)); F antihermitesch, Spur F = 0 punktweise",
      ok_d and ok_ah and dual[1] == -dual[-1])


def dichte(x, s):
    F = kruemmung(x, s)
    return 0.25*np.einsum('mnrs,mnij,rsji->', EPS4, F, F)/(8*np.pi**2)


c2, ok_inv, ok_cf, im_max = {}, True, True, 0.0
for s in (+1, -1):
    f = lambda r, s=s: dichte(np.array([r, 0, 0, 0]), s).real
    c2[s] = integrate.quad(lambda r: f(r)*2*np.pi**2*r**3, 0, np.inf, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
    sgn = np.sign(c2[s])
    for _ in range(8):
        x = rng.normal(size=4)
        q = dichte(x, s)
        im_max = max(im_max, abs(q.imag))
        ok_inv &= abs(q.real - f(np.linalg.norm(x))) < 1e-12
        ok_cf &= abs(q.real - sgn*6/(np.pi**2*(1 + x @ x)**4)) < 1e-12
ZAHLEN["c2"] = {"P_+": c2[1], "P_-": c2[-1]}
check(f"zweite Chernzahlen c2 = (1/8 pi^2) int tr F ^ F: P_+ {c2[1]:+.12f}, P_- {c2[-1]:+.12f}; Dichte exakt "
      "± 6/(pi^2 (1 + |x|^2)^4) (BPST, Größe 1) und SO(4)-invariant; selbstdual <-> c2 = -1, c2(P_+) + c2(P_-) = 0",
      abs(abs(c2[1]) - 1) < 1e-10 and abs(c2[1] + c2[-1]) < 1e-10 and ok_inv and ok_cf and im_max < 1e-12
      and np.sign(c2[1]) == -dual[1] and np.sign(c2[-1]) == -dual[-1])

rows = [np.kron(g.conj().T, np.eye(4)) - np.kron(np.eye(4), g) for g in GAM]   # vec(C g*) - vec(g C), spaltenweise
_, sv, vh = np.linalg.svd(np.vstack(rows))
nullity = int(np.sum(sv < 1e-10))
C = vh[-1].conj().reshape(4, 4, order='F')
C = C/np.sqrt(abs(np.linalg.det(C)))**0.5
CCb = C @ C.conj()
ok_j = nullity == 1 and all(np.allclose(C @ g.conj(), g @ C) for g in GAM) and np.allclose(CCb, CCb[0, 0]*I4) \
    and CCb[0, 0].real < 0 and abs(CCb[0, 0].imag) < 1e-12
C = C/np.sqrt(-CCb[0, 0].real)
ok_hopf = True
for _ in range(50):
    Z = rng.normal(size=4) + 1j*rng.normal(size=4)
    nZ = np.array([(Z.conj() @ g @ Z).real for g in GAM])/(Z.conj() @ Z).real
    Pz = (I4 + sum(nZ[A]*GAM[A] for A in range(5)))/2
    ok_hopf &= abs(np.linalg.norm(nZ) - 1) < 1e-12 and np.allclose(Pz @ Z, Z) and np.allclose(Pz @ C, C @ Pz.conj())
check("Einordnung: Es gibt genau eine (bis auf Faktor) antilineare Struktur j = C o konj mit C Gamma_A* = Gamma_A C, "
      "und j^2 = -1 (quaternionisch); jedes Bild von P_+(n) ist j-invariant (eine quaternionische Gerade), und "
      "Z -> n_A = <Z|Gamma_A|Z>/<Z|Z> bildet jedes Z auf S^4 ab mit Z im Bild von P_+(n): das ist die quaternionische "
      "Hopf- bzw. Twistorfaserung CP^3 -> S^4 = HP^1",
      ok_j and C.shape == (4, 4) and np.allclose(C @ C.conj(), -I4) and ok_hopf)
befund("E: Die Rechnung der Analyse stimmt (25/25, Dualität, c2 = -1 und +1). Es ist das klassische Paar Instanton/"
       "Anti-Instanton auf der runden S^4 (die beiden Spinorbündel bzw. die quaternionische Hopffaserung); die "
       "entgegengesetzten Ladungen sind durch die Whitney-Summe P_+ + P_- = trivial erzwungen (c1 = 0 auf S^4), also kein "
       "zusätzlicher Befund über eine Gegenstruktur", True)
befund("O: Die runde S^4 braucht die quaternionische Struktur j (j^2 = -1), also eine euklidische Realitätsstruktur auf "
       "dem Twistorraum; nativ ist dort die (2,2)-Form der CCR und die lineare Involution tau (z85, z89). Feste "
       "Geometrie, keine Dynamik, wie die Analyse selbst sagt", True)

# ------------------------------------------------------------------ Teil 2
say("\nTeil 2: Q_a = (P_a + sigma K_a)/2 in so(4,2)")
Gm = np.diag([-1., -1, 1, 1, 1, 1])                     # Achsen (e_-1, e_0, e_1, e_2, e_3, e_4)
E6 = np.eye(6)


def Mg(A, B):
    return np.outer(E6[A], Gm @ E6[B]) - np.outer(E6[B], Gm @ E6[A])


def kom(X, Y):
    return X @ Y - Y @ X


Pa = [Mg(a + 1, 5) + Mg(a + 1, 0) for a in range(4)]
Ka = [Mg(a + 1, 5) - Mg(a + 1, 0) for a in range(4)]
Mab = {(a, b): Mg(a + 1, b + 1) for a in range(4) for b in range(4)}
Dil = Mg(0, 5)
alle = Pa + Ka + list(Mab.values()) + [Dil]
ok_so = all(np.allclose(X.T @ Gm + Gm @ X, 0) for X in alle)
ok_ab = all(np.allclose(kom(Pa[a], Pa[b]), 0) and np.allclose(kom(Ka[a], Ka[b]), 0)
            for a in range(4) for b in range(a + 1, 4))
paare = [(a, b) for a in range(4) for b in range(a + 1, 4)]
vor, real_form, enthaelt_H, hl = {}, {}, {}, {}
H = Mg(0, 1)                                            # Drehung in der Ebene der beiden Zeitrichtungen
Xp, Xm = E6[5] + E6[1], E6[5] - E6[1]                    # Spitzen t = ±1, x = 0 des Einheitsdoppelkegels (Nullstrahlen)
for sig in (+1, -1):
    Q = [(Pa[a] + sig*Ka[a])/2 for a in range(4)]
    vs = set()
    for a, b in paare:
        K_ = kom(Q[a], Q[b])
        c_ = np.sum(K_*Mab[(a, b)])/np.sum(Mab[(a, b)]**2)
        assert np.allclose(K_, c_*Mab[(a, b)])
        vs.add(round(c_, 12))
    vor[sig] = vs.pop() if len(vs) == 1 else None
    basis = Q + [Mab[p] for p in paare]
    B = np.array([X.ravel() for X in basis]).T
    geschl = all(np.allclose(B @ np.linalg.lstsq(B, kom(X, Y).ravel(), rcond=None)[0], kom(X, Y).ravel())
                 for X in basis for Y in basis)
    kf = np.array([[np.trace(X @ Y) for Y in basis] for X in basis])
    nneg = int(np.sum(np.linalg.eigvalsh(kf) < -1e-9))
    real_form[sig] = ("so(4,1) (dS)" if nneg == 6 else "so(3,2) (AdS)" if nneg == 4 else f"? ({nneg})") if geschl else "offen"
    koef = np.linalg.lstsq(B, H.ravel(), rcond=None)[0]
    enthaelt_H[sig] = bool(np.allclose(B @ koef, H.ravel()))
    hl[sig] = bool(np.allclose(Q[0] @ Xp, Xp) and np.allclose(Q[0] @ Xm, -Xm))
ZAHLEN["so42"] = {"[Q_a,Q_b] = c M_ab, c(sigma=+1), c(sigma=-1)": [vor[1], vor[-1]], "reelle Form": real_form,
                  "enthält kompakte Energie": enthaelt_H, "Q_0 = Modulboost des Doppelkegels": hl}
check(f"[P_a, P_b] = [K_a, K_b] = 0 (je 6 Paare, beide abelsch); [Q_a, Q_b] = c M_ab mit c = {vor[1]:+.0f} für sigma = +1 "
      f"und c = {vor[-1]:+.0f} für sigma = -1, also [Q_a, Q_b] = ± sigma M_ab (das Vorzeichen ist Konvention von K bzw. "
      "M_ab); {Q_a, M_ab} schließt jeweils zu einer zehndimensionalen Algebra",
      ok_so and ok_ab and vor[1] is not None and vor[-1] is not None and abs(vor[1]) == 1 and vor[1] == -vor[-1]
      and "offen" not in real_form.values())
sig_ads = [s for s in (1, -1) if "AdS" in real_form[s]]
sig_ds = [s for s in (1, -1) if "dS" in real_form[s] and "AdS" not in real_form[s]]
check(f"reelle Formen (Spurform, Zahl der kompakten Richtungen): sigma = +1 -> {real_form[1]}, sigma = -1 -> "
      f"{real_form[-1]}. Genau die AdS-Wahl enthält die kompakte konforme Energie (Drehung der beiden Zeitachsen, "
      f"diskretes positives Spektrum, z83: E = L_0): {enthaelt_H}; genau die dS-Wahl hat als Q_0 den Modulboost des "
      f"Doppelkegels mit den Spitzen t = ±1 (Q_0 X_± = ± X_±, Hislop-Longo): {hl}",
      len(sig_ads) == 1 and len(sig_ds) == 1 and enthaelt_H[sig_ads[0]] and not enthaelt_H[sig_ds[0]]
      and hl[sig_ds[0]] and not hl[sig_ads[0]])
lam = np.sum(kom(Dil, Pa[0])*Pa[0])/np.sum(Pa[0]**2)
ok_dil = all(np.allclose(kom(Dil, Pa[a]), lam*Pa[a]) and np.allclose(kom(Dil, Ka[a]), -lam*Ka[a]) for a in range(4))
s_a = sig_ads[0]
Ql = lambda l: (l*Pa[0] + s_a*Ka[0]/l)/2
t_ = 0.7
ok_rad = np.allclose(expm(t_*Dil) @ Ql(1.3) @ expm(-t_*Dil), Ql(1.3*np.exp(lam*t_)))
check(f"[D, P] = {lam:+.0f} P, [D, K] = {-lam:+.0f} K; die AdS-Generatoren Q_0(l) = (l P_0 ± K_0/l)/2 zu verschiedenen "
      "Radien l sind unter der Dilatation konjugiert (e^{tD} Q_0(l) e^{-tD} = Q_0(l e^{± t})): in der konformen Struktur "
      "ist l eine Einheitenwahl, keine Invariante", ok_dil and abs(abs(lam) - 1) < 1e-12 and ok_rad)
z89 = lade("z89")
befund("E/C: Die Analyse stimmt ([Q_a, Q_b] = sigma M_ab bis auf Konvention). Invariant gesagt: Die Wahl, deren Q_0 die "
       "kompakte konforme Energie E ist, gibt so(3,2) (AdS) - das ist die E-invariante Wahl I_pi - I_omega von z89 "
       f"({'in z89 so gefunden' if hat(z89, 'I_pi - I_omega (AdS', 'befunde') else 'z89 nicht gefunden'}); die andere "
       "gibt so(4,1) (dS), und deren Q_0 ist der Modulgenerator eines Doppelkegels (Hislop-Longo). Das Vorzeichen von "
       "Lambda entspricht also der Wahl der Zeit in Q_0: radiale Energie -> AdS, modulare Zeit einer Kappe -> dS",
       hat(z89, "I_pi - I_omega (AdS", "befunde"))

# ------------------------------------------------------------------ Teil 3
say("\nTeil 3: Plebanski-Flächen aus einem reellen Vierbein")
ETA = np.diag([-1., 1, 1, 1])


def w11(A, B):
    return np.outer(A, B) - np.outer(B, A)


def w22(X, Y):
    return 0.25*np.einsum('mnrs,mn,rs->', EPS4, X, Y)


def sigmas(e):
    Sp = []
    for i in range(3):
        S = 1j*w11(e[0], e[i + 1]) - 0.5*sum(EPS3[i, j, k]*w11(e[j + 1], e[k + 1]) for j in range(3) for k in range(3))
        Sp.append(S)
    return Sp, [S.conj() for S in Sp]


def hodge_L(om, e):
    g = e.T @ ETA @ e
    gi = np.linalg.inv(g)
    return 0.5*np.sqrt(abs(np.linalg.det(g)))*np.einsum('mnab,ab->mn', EPS4, gi @ om @ gi.T)


n_rel, n_ges, lam_h, ok_h = 0, 0, set(), True
for _ in range(6):
    e = rng.normal(size=(4, 4))
    if np.linalg.det(e) < 0:
        e[3] *= -1
    vol = np.linalg.det(e)
    Sp, Sm = sigmas(e)
    for i in range(3):
        for j in range(3):
            n_rel += abs(w22(Sp[i], Sp[j]) - (-2j*(i == j)*vol)) < 1e-10*abs(vol)
            n_rel += abs(w22(Sm[i], Sm[j]) - (2j*(i == j)*vol)) < 1e-10*abs(vol)
            n_rel += abs(w22(Sp[i], Sm[j])) < 1e-10*abs(vol)
            n_ges += 3
    for i in range(3):
        hS = hodge_L(Sp[i], e)
        li = np.sum(hS*Sp[i].conj())/np.sum(abs(Sp[i])**2)
        ok_h &= np.allclose(hS, li*Sp[i]) and np.allclose(hodge_L(Sm[i], e), np.conj(li)*Sm[i])
        ok_h &= np.allclose(hodge_L(hodge_L(Sp[i], e), e), -Sp[i])
        lam_h.add(complex(round(li.real, 10), round(li.imag, 10)))
ec = rng.normal(size=(4, 4)) + 0.3j*rng.normal(size=(4, 4))
Spc, Smc = sigmas(ec)
gegen = max(abs(w22(Spc[i], Smc[j])) for i in range(3) for j in range(3))
ZAHLEN["plebanski"] = {"Relationen erfüllt": f"{n_rel}/{n_ges}", "Hodge-Eigenwert": str(lam_h),
                       "Gegenprobe komplexes Vierbein max|S+ ^ S-|": float(gegen)}
check(f"Sigma_+^i = i e^0 ^ e^i - (1/2) eps^i_jk e^j ^ e^k, Sigma_- = conj(Sigma_+): Sigma_+ ^ Sigma_+ = -2i delta vol, "
      f"Sigma_- ^ Sigma_- = +2i delta vol, Sigma_+ ^ Sigma_- = 0 an sechs zufälligen reellen Vierbeinen: {n_rel}/{n_ges} "
      f"(27 je Vierbein); *Sigma_+ = {lam_h} Sigma_+ (selbstdual im Lorentzschen, ** = -1); Gegenprobe komplexes "
      f"Vierbein: max|Sigma_+ ^ Sigma_-| = {gegen:.2f} (Realitätsbedingung verletzt)",
      n_rel == n_ges == 162 and ok_h and lam_h == {1j} and gegen > 1e-2)
z29 = lade("z29")
a_, b_ = sp.I, -1
check("Abgleich mit z29 (B^i = a e^0 ^ e^i + b (1/2) eps e ^ e, B ^ B = 2ab delta vol, Realität Re(a conj b) = 0): "
      "Sigma_+ ist z29s Lorentzpunkt a = i, b = -1 (2ab = -2i, Re(a conj b) = 0), dort als Signatur (-,+,+,+) gefunden",
      2*a_*b_ == -2*sp.I and sp.re(a_*sp.conjugate(b_)) == 0 and hat(z29, "a = i, b = 1: Signatur")
      and hat(z29, "Realitätsbedingung"))
befund("E/O: Die 27 Relationen stimmen. Sie sind die bekannte Richtung des Plebanski-Satzes: reelles Vierbein => "
       "Simplizität und Realität. z29 hatte die isotrope Simplizität aus S4 und den Lorentzpunkt a = i als Parameter; "
       "offen bleibt die Umkehrung, nämlich warum die Naht ein reelles Lorentzvierbein trägt (die Lorentzphase folgt "
       "nicht aus S4, z29)", True)

# ------------------------------------------------------------------ Teil 4
say("\nTeil 4: MacDowell-Mansouri")


def W(X, Y):
    return 0.25*np.einsum('abcd,mnrs,abmn,cdrs->', EPS4, EPS4, X, Y)


ok_mm, cR, cV = True, set(), set()
for _ in range(6):
    e = rng.normal(size=(4, 4))
    if np.linalg.det(e) < 0:
        e[3] *= -1
    de = np.linalg.det(e)
    A = rng.normal(size=(4, 4, 4, 4))
    R = A - A.transpose(1, 0, 2, 3)
    R = R - R.transpose(0, 1, 3, 2)
    ee = np.einsum('am,bn->abmn', e, e) - np.einsum('an,bm->abmn', e, e)
    Einv = np.linalg.inv(e)
    Rs = np.einsum('abmn,ma,nb->', R, Einv, Einv)
    for sig in (+1, -1):
        for l2 in (0.37, 2.9):
            F = R + (sig/l2)*ee
            ok_mm &= abs(W(F, F) - (W(R, R) + 2*(sig/l2)*W(R, ee) + (sig/l2)**2*W(ee, ee))) < 1e-9*(1 + abs(W(F, F)))
    ok_mm &= abs(W(R, ee) - W(ee, R)) < 1e-10*abs(W(R, ee))
    cR.add(round(W(R, ee)/(Rs*de), 10))
    cV.add(round(W(ee, ee)/de, 10))
check(f"eps_abcd F^ab ^ F^cd = eps R ^ R + (2 sigma/l^2) eps R ^ e ^ e + (1/l^4) eps e ^ e ^ e ^ e exakt (6 Vierbeine, "
      f"beliebige antisymmetrische R, sigma = ±1, zwei Radien): genau drei Terme; eps R ^ e ^ e = {cR} R vol, "
      f"eps e^4 = {cV} vol", ok_mm and cR == {2.0} and cV == {24.0})
sg, l, Gn, b = sp.symbols('sigma ell G b', real=True)
kappa_ = 1/(16*sp.pi*Gn)
loes = sp.solve([sp.Eq(b*(2*sg/l**2)*2, kappa_), ], [b], dict=True)[0]
Lam_ = sp.symbols('Lambda', real=True)
loesL = sp.solve(sp.Eq(-2*Lam_*kappa_, loes[b]*(sg**2/l**4)*24), Lam_)[0]
check(f"Abgleich mit (1/16 pi G) int (R - 2 Lambda) vol: b = {sp.simplify(loes[b])} (= l^2/(64 pi G sigma) für sigma^2 = 1), "
      f"Lambda = {sp.simplify(loesL.subs(sg**2, 1))} (für sigma^2 = 1: -3 sigma/l^2)",
      sp.simplify(loes[b] - l**2/(64*sp.pi*Gn*sg)) == 0 and sp.simplify(loesL.subs(sg, 1) + 3/l**2) == 0
      and sp.simplify(loesL.subs(sg, -1) - 3/l**2) == 0)
volS4 = integrate.quad(lambda r: (2/(1 + r**2))**4*2*np.pi**2*r**3, 0, np.inf, epsabs=1e-14, epsrel=1e-13)[0]
chi = 24*volS4/(32*np.pi**2)
b_val = 1/(64*np.pi)                                    # l = G = 1
check(f"Euler-Term: Auf der runden S^4 (R^ab = e^a ^ e^b/l^2, also F = 0 bei sigma = -1) ist (1/32 pi^2) int eps R ^ R = "
      f"24 vol(S^4)/(32 pi^2 l^4) = {chi:.12f} = chi(S^4); vol(S^4) = {volS4/np.pi**2:.12f} pi^2 l^4 (8/3); der "
      f"Euler-Anteil |b| 32 pi^2 chi = {b_val*32*np.pi**2*chi/np.pi:.12f} pi l^2/G ist genau die de-Sitter-Entropie "
      "A/4G = pi l^2/G", abs(chi - 2) < 1e-11 and abs(volS4 - 8*np.pi**2/3) < 1e-11
      and abs(b_val*32*np.pi**2*chi - np.pi) < 1e-12)
befund("E: Die MacDowell-Mansouri-Zerlegung, Lambda = -3 sigma/l^2 und b = l^2/(64 pi G sigma) stimmen. Die Wirkung ist "
       "eine Wahl (sie bricht so(3,2) bzw. so(4,1) durch die chirale Projektion auf Lorentz); die Algebra erlaubt sie, "
       "erzwingt sie nicht - wie die Analyse selbst sagt", True)

# ------------------------------------------------------------------ Teil 5
say("\nTeil 5: Lambda-G-Relation aus COSMO.LAM.01")
al, Gp = sp.symbols('alpha G', positive=True)
Mbar2 = 1/(8*sp.pi*Gp)
rhoL = sp.Rational(3, 4)/sp.pi**2*sp.exp(-2/al)*Mbar2**2
LamL = 8*sp.pi*Gp*rhoL
l2_ = 3/LamL
b_ = l2_/(64*sp.pi*Gp)
check(f"aus rho_L/Mbar^4 = (3/4 pi^2) e^(-2/alpha), Mbar^2 = 1/(8 pi G), Lambda = 8 pi G rho_L = 3/l^2 folgt symbolisch "
      f"l^2/G = {sp.simplify(l2_/Gp)} und |b| = {sp.simplify(b_)}",
      sp.simplify(l2_/Gp - 32*sp.pi**3*sp.exp(2/al)) == 0 and sp.simplify(b_ - sp.pi**2/2*sp.exp(2/al)) == 0)
mp.mp.dps = 30
ainv = mp.mpf('137.0359992168')
rho = 3/(4*mp.pi**2)*mp.e**(-2*ainv)
l2G = 32*mp.pi**3*mp.e**(2*ainv)
lPl = mp.mpf('1.616255e-35')
l_m = mp.sqrt(l2G)*lPl
Mbar_GeV = mp.mpf('2.435323e18')
ZAHLEN["lambda_G"] = {"rho_L/Mbar^4": mp.nstr(rho, 6), "e^(2/alpha)": mp.nstr(mp.e**(2*ainv), 6),
                      "l^2/G": mp.nstr(l2G, 6), "|b|": mp.nstr(mp.pi**2/2*mp.e**(2*ainv), 6),
                      "l/l_Pl": mp.nstr(mp.sqrt(l2G), 6), "l [m]": mp.nstr(l_m, 5),
                      "l [Mrd. Lichtjahre]": mp.nstr(l_m/mp.mpf('9.4607304725808e15')/1e9, 4),
                      "rho_L^(1/4) [meV]": mp.nstr(rho**0.25*Mbar_GeV*1e12, 6)}
check(f"Zahlen mit alpha^-1 = 137.0359992168 (EM.FP.01): rho_L/Mbar^4 = {mp.nstr(rho, 6)} (COSMO.LAM.01: 7.1253e-121), "
      f"l^2/G = {mp.nstr(l2G, 5)}, l = {mp.nstr(mp.sqrt(l2G), 4)} l_Pl = {mp.nstr(l_m, 4)} m "
      f"({ZAHLEN['lambda_G']['l [Mrd. Lichtjahre]']} Mrd. Lichtjahre), |b| = {mp.nstr(mp.pi**2/2*mp.e**(2*ainv), 4)}; "
      f"rho_L^(1/4) = {ZAHLEN['lambda_G']['rho_L^(1/4) [meV]']} meV (v7: 2.2375)",
      abs(rho/mp.mpf('7.12533e-121') - 1) < 1e-5 and abs(rho**0.25*Mbar_GeV*1e12 - mp.mpf('2.23747')) < 1e-4)
z90 = lade("z90")
befund("E: Die Relation l^2/G = 32 pi^3 e^(2/alpha) ist COSMO.LAM.01 in MacDowell-Mansouri-Variablen; sie hat keinen "
       "Inhalt über die historische Vakuumenergie hinaus. Das 'Herleitungsziel' der Analyse ist damit dasselbe wie die "
       "Herleitung von e^(-2/alpha) in rho_L", True)
befund("O: Vorzeichen: rho_L > 0 heißt Lambda > 0, also sigma = -1 in der Konvention der Analyse (dS, so(4,1)). Die native "
       "E-invariante Wahl ist AdS (z89, Teil 2); dS ist im Lorentzbüschel vorhanden, aber von keiner nativen Zeit "
       "ausgezeichnet, solange die Zeit in Q_0 die radiale Energie ist. Kandidat: Q_0 als modulare Zeit einer Kappe "
       "(Hislop-Longo, z90) gibt dS - nicht gezeigt", True)
befund("O: Skala: l ist in der konformen Struktur eine Einheitenwahl (Teil 2), z89s 'Radius = Vakuumgröße' setzt die "
       "Einheit. l^2/G braucht eine zweite Skala; die native Leiter gibt kein G (z90: 1/G_ind = 0 bei xi = 1/6)",
       hat(z90, "beim nativen xi = 1/6"))

# ------------------------------------------------------------------ Teil 6
say("\nTeil 6: bimetrische Gewichte, v358, Flächenkoeffizient")
w0, wm = 1 - mp.e**(-4*mp.pi), mp.e**(-4*mp.pi)
ZAHLEN["gewichte"] = {"1 - e^-4pi [%]": mp.nstr(100*w0, 10), "e^-4pi [%]": mp.nstr(100*wm, 6)}
check(f"Gewichte: 1 - e^(-4 pi) = {mp.nstr(100*w0, 10)} %, e^(-4 pi) = {mp.nstr(100*wm, 6)} % (Analyse: 99,9996513 % "
      "und 0,0003487 %)", mp.nstr(100*w0, 9) == '99.9996513' and mp.nstr(100*wm, 4) == '0.0003487')
s_, mu_, Mg2, Mf2 = sp.symbols('s mu M_g2 M_f2', positive=True)
Kin = sp.Matrix([[Mg2*s_ + mu_, -mu_], [-mu_, Mf2*s_ + mu_]])
Amp = sp.simplify(Kin.inv()[0, 0])
m2 = mu_*(Mg2 + Mf2)/(Mg2*Mf2)
r0 = sp.simplify(sp.limit(s_*Amp, s_, 0))
rm = sp.simplify(sp.limit((s_ + m2)*Amp, s_, -m2))
z48 = lade("z48")
alpha_hr = mp.sqrt(1/(mp.e**(4*mp.pi) - 1))
rel48 = abs(alpha_hr/mp.mpf(z48["zahlen"]["alpha_bimetrik"]) - 1)
ZAHLEN["gewichte"]["M_f/M_g"] = mp.nstr(alpha_hr, 15)
ZAHLEN["gewichte"]["relative Abweichung zu z48 (Gleitkomma dort)"] = mp.nstr(rel48, 3)
check(f"Hassan-Rosen-Lesart (Materie an g, Fierz-Pauli-Kopplung mu (h_g - h_f)^2): Residuen masselos {r0}, massiv {rm}, "
      f"Verhältnis M_f^2/M_g^2; mit massivem Gewicht e^(-4 pi) folgt M_f/M_g = 1/sqrt(e^(4 pi) - 1) = "
      f"{mp.nstr(alpha_hr, 12)} (z48: {z48['zahlen']['alpha_bimetrik']}, relativ {mp.nstr(rel48, 2)})",
      sp.simplify(rm/r0 - Mf2/Mg2) == 0 and sp.simplify(r0 - 1/(Mg2 + Mf2)) == 0
      and rel48 < 1e-10 and hat(z48, "massiver Anteil"))
check("v358: |mu_4| = |Z_2| chi(S^2) = 2 * 2 = 4 (arithmetisch exakt); der Satz 'Zentralladung = Flächenkoeffizient' gilt "
      "nur in 1+1D, in 4D hängt der Flächenkoeffizient am Regulator (z90: Dimensionstypisierung, induziertes 1/G = 0 bei "
      "xi = 1/6)", 2*2 == 4 and hat(z90, "Dimensionstypisierung", "befunde") and hat(z90, "beim nativen xi = 1/6"))
befund("E/L: Die Gewichte sind z47/z48 (massloser Pol 1 - e^(-4 pi), massive Linie e^(-4 pi)); die bimetrische Lesart "
       "M_f/M_g = 0,0018674 steht schon in z48. Neu ist nur die Zuordnung zur Gegenstruktur - eine Lesart, keine "
       "Herleitung (wie die Analyse selbst sagt)", True)
befund("E: Abschnitt 8 der Analyse trifft den bekannten Punkt: Die Flächennormierung von v358 ist nicht universell "
       "(z90, GRAV.ENTROPY.EQUILIBRIUM.01 typisiert); die vorgeschlagene doppelte Herleitung von G (thermodynamisch und "
       "geometrisch) ist eine Strategie, kein Ergebnis", True)

urteil = ("Die zwölfte externe Analyse hält in allen gerechneten Punkten: 25/25 Cliffordrelationen, F_+ und F_- selbstdual "
          "bzw. antiselbstdual mit c2 = -1 und +1, [Q_a, Q_b] = sigma M_ab (bis auf Konvention), 27 Plebanski-Relationen "
          "aus einem reellen Vierbein, die MacDowell-Mansouri-Zerlegung in Euler, Einstein-Hilbert und Volumen mit "
          "Lambda = -3 sigma/l^2 und b = l^2/(64 pi G sigma), l^2/G = 32 pi^3 e^(2/alpha), |b| = (pi^2/2) e^(2/alpha), "
          "Gewichte 99,9996513 % und 0,0003487 %. Eingeordnet: (1) Die Clifford-Projektionen sind die quaternionische "
          "Hopf- bzw. Twistorfaserung CP^3 -> S^4 (Instanton und Anti-Instanton); die entgegengesetzten Ladungen erzwingt "
          "die Whitney-Summe; die runde S^4 braucht eine euklidische quaternionische Struktur, nicht die native "
          "(2,2)-Form. (2) Invariant gesagt gibt die Wahl, deren Q_0 die kompakte konforme Energie ist, AdS (die "
          "E-invariante Wahl von z89), die andere dS, deren Q_0 der Modulgenerator eines Doppelkegels ist; der Radius ist "
          "unter der Dilatation eine Einheitenwahl. (3) Reelles Vierbein => Plebanski ist die bekannte Richtung (z29: "
          "a = i); offen bleibt, warum die Naht ein reelles Lorentzvierbein trägt. (4) Die Lambda-G-Relation ist "
          "COSMO.LAM.01 in neuen Variablen, ohne zusätzlichen Inhalt; sie verlangt dS, während die native E-invariante Wahl "
          "AdS ist; eine zweite Skala (G) fehlt nativ (z90). (5) Die Gewichte und die bimetrische Lesart M_f/M_g = "
          "0,0018674 stehen in z47/z48. Die Analyse liefert ein sauberes algebraisches Gerüst (MacDowell-Mansouri aus der "
          "konformen Algebra), aber keinen nativen Schritt: Vierbein, Wirkungsfunktional, Skala und Vorzeichen von Lambda "
          "bleiben offen. Neuer Kandidat: Das Vorzeichen von Lambda entspricht der Zeit in Q_0 (radiale Energie -> AdS, "
          "modulare Zeit einer Kappe -> dS).")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z96_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
