"""Z126: Folgeauftrag 'modulare Statistik': gradierte Permutationsorbifaltigkeit (P_ident mit P_CAR) über N = 2 und 4 Nähte,
P_glue unter Modulsummation und Klassifikation der modulverträglichen Statistikregeln der Naht.

Leitlinie 'Rekonstruktion statt Auswahl': Eine Statistikregel zählt als modulverträglich, wenn die Bahnergänzung ihres projizierten
Originalsektors existiert und eine Hilbertraumspur ist. Bahnergänzung heißt wie in z113: jede Einschubspur des Originalsektors wird
zu ihrer SL(2,Z)-Bahn mit trivialem Charakter ergänzt, und die Koeffizienten aller Bahnmitglieder müssen mit den
Projektorkoeffizienten übereinstimmen. H6 und die gelösten Blöcke (Flavour, Massen, alpha, die 16, E8) sind feste Eingaben und
werden nicht neu geöffnet; nichts wird befördert. Firewalls: DIMENSION.UPLIFT.FIREWALL.01 (Helizität ist Kinematik der radialen
Einpaarschicht), GRAV.NONCIRCULAR.01. Statuskürzel: E exakt, C bedingt, L Lesart, O offen.

Bestand: z113 (Rahmen: Original = zwei komplexe chirale Fermionen mit Moden Z + 1/4 und Z + 3/4, Partner antiholomorph; Uhr
R = i^(Q_a - Q_b) wirkt als T, R^2 als T^2; U_2pi = (-1)^(N - Nquer) ist die halbe Drehung der Oszillatoren, ihre verdrillten
Sektoren sind lange Nähte; N Nähte mit gemeinsamer Uhrverdrillung; mit summierter 2 pi-Drehung pro Naht nur U = +1 oder G = -1
positiv), z106 (gemeinsame CAR = Super-Tensorprodukt, Koszul-Vorzeichen = R^2-Parität, Gamma(pi) gaußsch), z119 (Gamma(pi) = +1 auf
dem Singulett, -1 auf dem Triplett), z117 (P_glue: Spinormodul neutral unter U_(pi/2) (x) R, also z97s Z4-Regel R = Z_SU^-1 auf
dem Spinormodul), z97 (Z_SU = -i^(2h), Z_SU^2 = U_2pi, G = U_2pi R^2 schließt), z82 (R = i Gamma_11).

Modell (Erweiterung von z113, exakt in den Einschüben):
- Vierteldrehung A = i^(N - Nquer) auf den Oszillatoren beider Quarter (Z_SU = -A auf der Einpaarschicht, A^2 = U_2pi), als
  Oszillatorspur e^(2 pi i o/96)/eta(w + o/4) mit o = 1, 3. Ihre Bahnen haben 96 Elemente (lange Nähte bis Länge 16).
- Permutationen: Ein k-Zyklus mit Nahteinschüben x_i hat die Spur Tr_C(K^(k-1) prod x_i q^(k L_0)) (Koszul-Identität, 1.1).
  K = R^2 des Originals (P_CAR wie z106), K = R^2 Rquer^2 (J-kovariante CAR) oder K = 1 (ungradiert).
- Gruppe = Permutationen x Nahteinschübe (Kranzprodukt) x gemeinsame Nativprojektion. Charaktere sind exakt (Phasen als Brüche).
- Gemeinsame Bahnen werden als Multimengen kanonisiert (gleiche Faktorbahn und Seite), damit eine Bahn nie doppelt zählt.
- Positivität wie z113: spinaufgelöste Fouriermoden bei Y = 0,35, 0,6 und 1,0 mit |s| <= 5, normiert auf die Nativsumme gleicher
  N; positiv heißt Minimum > -max(1e-8, 1000 x numerisches Rauschen); das Rauschen wird aus den Imaginärteilen geschätzt.

Teil 1  Gradierte Permutationsorbifaltigkeit: Koszul-Spuridentität; der gradierte N-Zyklus als gaußsches Z_N der Fouriermoden
        (Thetaidentität); Anomalien (P_CAR des Originals anomal für gerade N, J-kovariant frei, gekleideter Tausch Gamma G_1 frei);
        Positivität für N = 2, 3, 4 (pro Naht, diagonal, mit Permutationen, gekleidet); direkte S/T-Invarianz der Schlüsselsummen.
Teil 2  P_glue: Projektoralgebra; Bahnen der Vierteldrehungen (A X und U A X liegen in einer Bahn mit Koeffizient 1); keine
        modulare Ergänzung von P_glue für N = 1..4; P_A3 ergänzt sich nur zu einem Nicht-Projektor.
Teil 3  Klassifikation der Hebungen der 2 pi-Drehung nach N, Entscheidung Spin-Statistik, Folge für die Szenario-Klausel
        modular_statistics.
"""
import itertools
import json
import time
from collections import defaultdict, deque
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).parent
OUT = []
CHECKS, BEFUNDE, ZAHLEN = {}, {}, {}
T0 = time.time()
mp.mp.dps = 25
I = mp.mpc(0, 1)
PI = mp.pi
TP = [mp.mpc('0.137', '1.093'), mp.mpc('-0.211', '0.947'), mp.mpc('0.389', '1.211')]   # Testpunkte (wie z113)
NULL = mp.mpf(10)**-16
QU, HA = Fr(1, 4), Fr(1, 2)
rng = np.random.default_rng(126)


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def check(name, cond):
    CHECKS[name] = bool(cond)
    say(f"   [{'ok' if cond else 'FEHLER'}] {name}")


def befund(name, val):
    BEFUNDE[name] = bool(val)
    say(f"   [Befund: {'ja' if val else 'nein'}] {name}")


def zeit():
    return f"(Laufzeit bisher {time.time() - T0:.0f} s)"


# ------------------------------------------------------------------ Theta, eta, Modulgruppe (aus z113)
def R_(x):
    return mp.mpf(x.numerator)/x.denominator if isinstance(x, Fr) else mp.mpf(x)


def EZ(r):
    return mp.exp(2*I*PI*R_(Fr(r)))


def theta(a, b, tau):
    a, b = R_(a), R_(b)
    nmax = int(mp.sqrt((mp.mp.dps*2.31 + 12)/(PI*tau.imag))) + 3
    return mp.fsum(mp.exp(I*PI*tau*(n + a)**2 + 2*I*PI*(n + a)*b) for n in range(-nmax, nmax + 1))


def eta(tau):
    kmax = int(mp.sqrt((mp.mp.dps*2.31 + 12)/(3*PI*tau.imag))) + 3
    return mp.exp(2*I*PI*tau/24)*mp.fsum((-1)**k*mp.exp(I*PI*tau*k*(3*k - 1)) for k in range(-kmax, kmax + 1))


def act(g, tau):
    return (g[0]*tau + g[1])/(g[2]*tau + g[3])


def mul(g, h):
    return (g[0]*h[0] + g[1]*h[2], g[0]*h[1] + g[1]*h[3], g[2]*h[0] + g[3]*h[2], g[2]*h[1] + g[3]*h[3])


GEN = {"S": (0, -1, 1, 0), "T": (1, 1, 0, 1)}
EINS = (1, 0, 0, 1)


def wurzel(c, nmax=4608):
    if abs(abs(c) - 1) > 1e-12:
        return None
    r = Fr(float(mp.arg(c)/(2*PI))).limit_denominator(nmax)
    return r % 1 if abs(mp.exp(2*I*PI*R_(r)) - c) < 1e-12 else None


def werte(F, g):
    return [F(act(g, t)) for t in TP]


def proportional(v, w):
    if any(abs(x) < NULL for x in v + w):
        return mp.mpc(1) if all(abs(x) < NULL for x in v + w) else None
    c = v[0]/w[0]
    return c if all(abs(x/y - c) < mp.mpf(10)**-12*abs(c) for x, y in zip(v[1:], w[1:])) else None


def bahn(F, maxn=800):
    """SL(2,Z)-Bahn {F o g}: Vertreter, Werte, Permutationen und Phasen (Einheitswurzeln) von S und T (z113)."""
    reps, vals = [EINS], [werte(F, EINS)]
    perm, ph = {"S": {}, "T": {}}, {"S": {}, "T": {}}
    k = 0
    while k < len(reps):
        for gn, gen in GEN.items():
            g2 = mul(reps[k], gen)
            v = werte(F, g2)
            for j, w in enumerate(vals):
                c = proportional(v, w)
                if c is not None:
                    perm[gn][k], ph[gn][k] = j, wurzel(c)
                    break
            else:
                reps.append(g2)
                vals.append(v)
                perm[gn][k], ph[gn][k] = len(reps) - 1, Fr(0)
                assert len(reps) <= maxn
        k += 1
    assert all(x is not None for d in ph.values() for x in d.values())
    return {"F": F, "reps": reps, "vals": vals, "perm": perm, "ph": ph, "n": len(reps)}


def lokalisiere(B, H):
    v = werte(H, EINS)
    for j, w in enumerate(B["vals"]):
        c = proportional(v, w)
        if c is not None:
            return j, wurzel(c)
    return None


def konsistenz(n, kanten, chi):
    V = [None]*n
    V[0] = Fr(0)
    adj = defaultdict(list)
    for (i, gn, j, r) in kanten:
        adj[i].append((gn, j, r))
    stack = [0]
    while stack:
        i = stack.pop()
        for (gn, j, r) in adj[i]:
            if V[j] is None:
                V[j] = (V[i] + r - chi[gn]) % 1
                stack.append(j)
    return V, {(V[i] + r - chi[gn] - V[j]) % 1 for (i, gn, j, r) in kanten} - {Fr(0)}


def charaktere(n, kanten):
    out = []
    for k in range(24):
        for j in range(4):
            chi = {"T": Fr(k, 24), "S": Fr(j, 4)}
            V, bad = konsistenz(n, kanten, chi)
            if not bad:
                out.append((chi, V))
    return out


# ------------------------------------------------------------------ Nahtfaktoren
def L_lad(u, t, tau):
    """sum_Q q^((Q+u)^2/2) e^(2 pi i t Q): Ladungsspur eines Quarters mit Einschub e^(2 pi i t Q) (z113)."""
    return mp.exp(-2*I*PI*R_(t)*R_(u))*theta(u, t, tau)


def osz(o, w):
    """Tr(i^(o N) q^(N - 1/24)) einer chiralen Bose-Leiter: o = 0 frei, o = 2 die 2 pi-Drehung U, o = 1, 3 die Vierteldrehung."""
    return mp.exp(2*I*PI*o/96)/eta(w + mp.mpf(o)/4)


def hol(ta, tb, oa, ob, k=1):
    """holomorphe Originalspur bei k tau: Ladungseinschub e^(2 pi i (ta Q_a + tb Q_b)), Oszillatoreinschübe i^(oa N_a + ob N_b)."""
    def f(tau):
        w = k*tau
        return L_lad(-QU, ta, w)*L_lad(QU, tb, w)*osz(oa, w)*osz(ob, w)
    return f


FB, LOC = {}, {}


def spec_norm(s):
    return (Fr(s[0]) % 1, Fr(s[1]) % 1, s[2] % 4, s[3] % 4, s[4])


def locate(spec):
    """Faktorfunktion -> (Faktorbahn, Element, Phase); neue Faktorbahnen werden mit bahn() bestimmt."""
    spec = spec_norm(spec)
    if spec in LOC:
        return LOC[spec]
    f = hol(*spec)
    for nm, B in FB.items():
        if B["k"] == spec[4]:
            x = lokalisiere(B, f)
            if x is not None:
                LOC[spec] = (nm, x[0], x[1])
                return LOC[spec]
    B = bahn(f)
    nm = f"F{len(FB)}"
    B.update({"k": spec[4], "name": nm, "spec": spec})
    FB[nm] = B
    LOC[spec] = (nm, 0, Fr(0))
    return LOC[spec]


def gemeinsam(groups, start):
    """Bahn eines Produkts von Faktoren; Zustände als Multimengen je (Faktorbahn, Seite); Phasen addiert."""
    idx, states, kanten = {start: 0}, [start], []
    dq = deque([start])
    while dq:
        st = dq.popleft()
        for gn in ("S", "T"):
            nst = tuple(tuple(sorted(B["perm"][gn][s] for s in ss)) for (B, e), ss in zip(groups, st))
            phase = sum((e*B["ph"][gn][s] for (B, e), ss in zip(groups, st) for s in ss), Fr(0)) % 1
            if nst not in idx:
                idx[nst] = len(states)
                states.append(nst)
                dq.append(nst)
            kanten.append((idx[st], gn, idx[nst], phase))
    return states, idx, kanten


ORB = {}


def atom_loc(fac):
    """fac = [(spec, Exponent)], Exponent < 0: komplex konjugiert (Partner). -> (Signatur, Zustand, Phase)"""
    d = defaultdict(list)
    ph0 = Fr(0)
    for spec, e in fac:
        nm, j, r = locate(spec)
        for _ in range(abs(e)):
            d[(nm, 1 if e > 0 else -1)].append(j)
        ph0 += e*r
    keys = sorted(d)
    return tuple(keys), tuple(tuple(sorted(d[k])) for k in keys), ph0 % 1


def orbit_of(sig, state):
    for O in ORB.setdefault(sig, []):
        if state in O["idx"]:
            return O
    groups = [(FB[nm], sgn) for nm, sgn in sig]
    states, idx, kn = gemeinsam(groups, state)
    ch = charaktere(len(states), kn)
    triv = [V for c, V in ch if c["T"] == 0 and c["S"] == 0]
    O = {"sig": sig, "groups": groups, "states": states, "idx": idx, "chars": [f"T {c['T']}, S {c['S']}" for c, _ in ch],
         "V": triv[0] if triv else None, "grid": None}
    ORB[sig].append(O)
    return O


def ergaenze(atoms):
    """Bahnergänzung: atoms = [(fac, Koeffizient)] des projizierten Originalsektors. Gibt Terme (Bahn, lambda, konsistent)."""
    members = defaultdict(lambda: defaultdict(mp.mpc))
    olist = []
    for fac, coeff in atoms:
        sig, state, ph0 = atom_loc(fac)
        O = orbit_of(sig, state)
        if all(O is not x for x in olist):
            olist.append(O)
        members[id(O)][O["idx"][state]] += mp.mpc(coeff)*EZ(ph0)
    terms, ok = [], True
    for O in olist:
        mem = members[id(O)]
        if all(abs(v) < mp.mpf(10)**-18 for v in mem.values()):
            terms.append((O, mp.mpc(0), True))
            continue
        if O["V"] is None:
            terms.append((O, None, False))
            ok = False
            continue
        lam = [v/EZ(O["V"][j]) for j, v in mem.items()]
        kons = all(abs(x - lam[0]) < mp.mpf(10)**-18*max(1, abs(lam[0])) for x in lam)
        ok &= kons
        terms.append((O, lam[0] if kons else None, kons))
    return {"terms": terms, "ok": ok}


# ------------------------------------------------------------------ Gitter und Positivität
YS = (mp.mpf('0.35'), mp.mpf('0.6'), mp.mpf('1.0'))
KS = 24
GRID = [mp.mpc(mp.mpf(k)/KS, Y) for Y in YS for k in range(KS)]
GV = {}


def gridvals(B):
    if B["name"] not in GV:
        GV[B["name"]] = [[B["F"](act(g, t)) for t in GRID] for g in B["reps"]]
    return GV[B["name"]]


def orbit_grid(O):
    if O["grid"] is None:
        vals = [gridvals(B) for B, _ in O["groups"]]
        tot = [mp.mpc(0)]*len(GRID)
        for st, j in O["idx"].items():
            w = EZ(O["V"][j])
            for p in range(len(GRID)):
                x = w
                for (B, e), ss, vv in zip(O["groups"], st, vals):
                    for s in ss:
                        x *= vv[s][p] if e > 0 else mp.conj(vv[s][p])
                tot[p] += x
        O["grid"] = tot
    return O["grid"]


def Zgrid(Z):
    tot = [mp.mpc(0)]*len(GRID)
    for O, lam, kons in Z["terms"]:
        if lam is not None and abs(lam) > 1e-14:
            g = orbit_grid(O)
            tot = [a + lam*b for a, b in zip(tot, g)]
    return tot


def moden(vals):
    out = {}
    for iy, Y in enumerate(YS):
        v = vals[iy*KS:(iy + 1)*KS]
        out[Y] = {s: mp.fsum(x*mp.exp(-2*I*PI*s*k/KS) for k, x in enumerate(v))/KS for s in range(-5, 6)}
    return out


REF = {}


def positiv(Z, N):
    """(kleinster normierter Spinmode, Rauschen, positiv?) relativ zur Nativsumme gleicher N."""
    m = moden(Zgrid(Z))
    ref = REF[N]
    mn = float(min(m[Y][s].real/abs(ref[Y][0]) for Y in YS for s in m[Y]))
    noise = float(max(abs(m[Y][s].imag)/abs(ref[Y][0]) for Y in YS for s in m[Y]))
    return mn, noise, mn > -max(1e-8, 1e3*noise)


def wert(Z, tau, cache):
    """Zustandssumme an einem Punkt (für die direkte S/T-Prüfung)."""
    tot = mp.mpc(0)
    for O, lam, kons in Z["terms"]:
        if lam is None or abs(lam) < 1e-14:
            continue
        for st, j in O["idx"].items():
            x = lam*EZ(O["V"][j])
            for (B, e), ss in zip(O["groups"], st):
                for s in ss:
                    key = (B["name"], s, mp.nstr(tau, 30))
                    if key not in cache:
                        cache[key] = B["F"](act(B["reps"][s], tau))
                    x *= cache[key] if e > 0 else mp.conj(cache[key])
            tot += x
    return tot


# ------------------------------------------------------------------ Gruppen: Kranzprodukt (Permutation, Nahteinschübe)
# Nahteinschub x = (a, r) heißt A^a R^r (A = i^(N - Nquer), R = i^(Q_a - Q_b) einseitig auf dem Original);
# U = (2, 0), G = U R^2 = (2, 2), R^2 = (0, 2).
U_, G_, A_, R2_ = (2, 0), (2, 2), (1, 0), (0, 2)
Z0 = (Fr(0), Fr(0), 0, 0)
KOS = {"ungradiert": (Z0, Z0), "Koszul (Original)": ((HA, HA, 0, 0), Z0), "Koszul (J-kovariant)": ((HA, HA, 0, 0), (HA, HA, 0, 0))}


def ins_params(x):
    a, r = x
    return (Fr(r, 4), Fr(-r, 4), a % 4, a % 4), (Fr(0), Fr(0), a % 4, a % 4)


def addp(p, q):
    return ((p[0] + q[0]) % 1, (p[1] + q[1]) % 1, (p[2] + q[2]) % 4, (p[3] + q[3]) % 4)


def nat(b):
    return (Fr(b, 4), Fr(-b, 4), 0, 0)


def zyklen(perm):
    seen, out = set(), []
    for i in range(len(perm)):
        if i not in seen:
            c, j = [], i
            while j not in seen:
                seen.add(j)
                c.append(j)
                j = perm[j]
            out.append(c)
    return out


def spur_atom(perm, xs, b, kos):
    """Spur von (perm, xs) n^b im Originalsektor der N Nähte: je Zyklus der Länge k Tr_C(K^(k-1) prod x_i n^(k b) q^(k L_0))."""
    fac = []
    for c in zyklen(perm):
        L = len(c)
        h, a = Z0, Z0
        for i in c:
            hh, aa = ins_params(xs[i])
            h, a = addp(h, hh), addp(a, aa)
        for _ in range(L - 1):
            h, a = addp(h, kos[0]), addp(a, kos[1])
        h, a = addp(h, nat(b*L)), addp(a, nat(b*L))
        fac += [((*h, L), 1), ((*a, L), -1)]
    return fac


def gmul(g, h):
    (p, x), (q, y) = g, h
    n = len(p)
    return (tuple(p[q[i]] for i in range(n)),
            tuple(((x[q[i]][0] + y[i][0]) % 4, (x[q[i]][1] + y[i][1]) % 4) for i in range(n)))


def huelle(gens):
    n = len(gens[0][0])
    e = (tuple(range(n)), tuple((0, 0) for _ in range(n)))
    G = {e: ()}
    front = [e]
    while front:
        new = []
        for g in front:
            for k, s in enumerate(gens):
                h = gmul(g, s)
                if h not in G:
                    G[h] = G[g] + (k,)
                    new.append(h)
        front = new
    return G


def el(perm, xs):
    return (tuple(perm), tuple(xs))


def naht(N, k, x):
    xs = [(0, 0)]*N
    xs[k] = x
    return el(range(N), xs)


def zyk(N):
    return el(tuple((i + 1) % N for i in range(N)), [(0, 0)]*N)


def regel(gens, chis, kos="ungradiert"):
    """Orbifaltigkeit nach der von gens erzeugten Gruppe mit dem Charakter chi(gen_k) = e^(2 pi i chis[k])."""
    G = huelle(gens)
    chi = {g: sum((Fr(chis[k]) for k in w), Fr(0)) % 1 for g, w in G.items()}
    hom = all((chi[gmul(g, h)] - chi[g] - chi[h]) % 1 == 0 for g in G for h in G)
    assert hom
    atoms = [(spur_atom(g[0], g[1], b, KOS[kos]), EZ(-chi[g])) for g in G for b in range(4)]
    return ergaenze(atoms), len(G)


def bahninfo(Z):
    return [{"Bahn": len(O["states"]), "Charaktere": O["chars"], "konsistent": kons,
             "lambda": None if lam is None else mp.nstr(lam, 8)} for O, lam, kons in Z["terms"]]


POS = {}


def pos_eintrag(label, N, gens, chis, kos="ungradiert"):
    Z, ordn = regel(gens, chis, kos)
    if not Z["ok"]:
        POS[label] = {"N": N, "Gruppe": ordn, "modular": False, "Bahnen": bahninfo(Z)}
        say(f"      {label}: keine modulare Ergänzung (anomal oder inkonsistent)  {zeit()}")
        return None
    mn, noise, ok = positiv(Z, N)
    POS[label] = {"N": N, "Gruppe": ordn, "modular": True, "min": mn, "Rauschen": noise, "positiv": ok,
                  "Bahnen": bahninfo(Z)}
    say(f"      {label}: min {mn:.3g} (Rauschen {noise:.0e}) -> {'positiv' if ok else 'NEGATIV'}  {zeit()}")
    return Z


# =================================================================== Teil 1
say("Teil 1: Gradierte Permutationsorbifaltigkeit über N Nähte (P_ident mit P_CAR)")

# --- 1.1 Koszul-Spuridentität (endlicher Superraum)
p, q = 2, 2
K1 = np.diag([1.0]*p + [-1.0]*q)
par = [0]*p + [1]*q
A0 = rng.normal(size=(p, p)) + 1j*rng.normal(size=(p, p))
A1 = rng.normal(size=(q, q)) + 1j*rng.normal(size=(q, q))
Am = np.zeros((p + q, p + q), complex)
Am[:p, :p], Am[p:, p:] = A0, A1
koszul = {}
for N in (2, 3, 4):
    dim = (p + q)**N
    Cg = np.zeros((dim, dim))
    Cu = np.zeros((dim, dim))
    for idx in itertools.product(range(p + q), repeat=N):
        src = np.ravel_multi_index(idx, (p + q,)*N)
        tgt = np.ravel_multi_index((idx[-1],) + idx[:-1], (p + q,)*N)
        sgn = (-1)**(par[idx[-1]]*sum(par[i] for i in idx[:-1]))
        Cg[tgt, src] = sgn
        Cu[tgt, src] = 1
    AN = Am
    for _ in range(N - 1):
        AN = np.kron(AN, Am)
    Apow = np.linalg.matrix_power(Am, N)
    gr = np.trace(Cg @ AN)
    un = np.trace(Cu @ AN)
    koszul[N] = {"gradiert - Tr(K^(N-1) A^N)": float(abs(gr - np.trace(np.linalg.matrix_power(K1, N - 1) @ Apow))),
                 "ungradiert - Tr(A^N)": float(abs(un - np.trace(Apow))),
                 "|gradiert - ungradiert|": float(abs(gr - un))}
ZAHLEN["koszul_spuridentitaet (Superraum 2|2, zufälliges gerades A)"] = koszul
check("[E] Koszul-Spuridentität: Für den gradierten N-Zyklus Gamma(sigma) auf V^(x)N (Superraum 2|2, zufälliges gerades A) gilt "
      "Tr(Gamma(sigma) A^(x)N) = Tr(K^(N-1) A^N), für den ungradierten Tr(A^N) (Fehler "
      f"{max(max(v['gradiert - Tr(K^(N-1) A^N)'], v['ungradiert - Tr(A^N)']) for v in koszul.values()):.0e}). Gegenprobe: "
      "Für gerade N unterscheiden sich beide (|Differenz| "
      f"{koszul[2]['|gradiert - ungradiert|']:.2f}, {koszul[4]['|gradiert - ungradiert|']:.2f}), für N = 3 nicht "
      f"({koszul[3]['|gradiert - ungradiert|']:.0e}): Der Koszul-Faktor eines k-Zyklus ist K^(k-1), also die R^2-Parität für "
      "gerade Zyklen",
      all(v["gradiert - Tr(K^(N-1) A^N)"] < 1e-9 and v["ungradiert - Tr(A^N)"] < 1e-9 for v in koszul.values())
      and koszul[2]["|gradiert - ungradiert|"] > 1e-3 and koszul[4]["|gradiert - ungradiert|"] > 1e-3
      and koszul[3]["|gradiert - ungradiert|"] < 1e-9)

# --- 1.2 der gradierte N-Zyklus ist das gaußsche Z_N der Fouriermoden
fourier = {}
for N in (2, 3, 4):
    Kp = HA if N % 2 == 0 else Fr(0)
    rg, ru = [], []
    for t in TP:
        prod = mp.mpc(1)
        for k in range(N):
            prod *= hol(Fr(k, N), Fr(k, N), 0, 0)(t)
        rg.append(hol(Kp, Kp, 0, 0, N)(t)/prod)
        ru.append(hol(0, 0, 0, 0, N)(t)/prod)
    fourier[N] = {"gradiert/Fourierprodukt": mp.nstr(rg[0], 12), "|1 - Verhältnis|": float(abs(1 - rg[0])),
                  "Streuung gradiert": float(max(abs(x - rg[0]) for x in rg)),
                  "Streuung ungradiert (rel.)": float(max(abs(x - ru[0]) for x in ru)/abs(ru[0]))}
ZAHLEN["fourierzerlegung_gradierter_zyklus"] = fourier
check("[E] Der gradierte N-Zyklus ist gaußsch: Seine Spur auf N Originalen ist exakt das Produkt der Spuren der N Fouriermoden "
      "mit den Phasen e^(2 pi i k/N) auf beiden Quartern, Tr_C(K^(N-1) q^(N L_0)) = prod_k Tr(e^(2 pi i (k/N)(Q_a + Q_b)) "
      f"q^(L_0)) (Verhältnis 1, Streuung {max(v['Streuung gradiert'] for v in fourier.values()):.0e} für N = 2, 3, 4; "
      "Thetaidentität). Die langen Nähte des gradierten Zyklus sind also Fouriermoden mit Verdrillung k/N. Gegenprobe: Die "
      "ungradierte Spur ist für N = 2 und 4 kein solches Produkt (relative Streuung "
      f"{fourier[2]['Streuung ungradiert (rel.)']:.3f}, {fourier[4]['Streuung ungradiert (rel.)']:.3f}); für N = 3 sind beide "
      "Spuren gleich",
      all(v["Streuung gradiert"] < 1e-20 and v["|1 - Verhältnis|"] < 1e-20 for v in fourier.values())
      and fourier[2]["Streuung ungradiert (rel.)"] > 1e-3
      and fourier[4]["Streuung ungradiert (rel.)"] > 1e-3 and fourier[3]["Streuung ungradiert (rel.)"] < 1e-20)
say(f"   {zeit()}")

# --- 1.3 Anomalien der Permutationen
for N in (1, 2, 3, 4):
    Z, _ = regel([el(range(N), [(0, 0)]*N)], [0])
    REF[N] = moden(Zgrid(Z))
anom = {}
for N in (2, 4):
    for kos in KOS:
        Z, _ = regel([zyk(N)], [0], kos)
        anom[f"N = {N}, {kos}"] = {"modular": Z["ok"], "Bahnen": [(len(O["states"]), O["chars"]) for O, _, _ in Z["terms"]]}
for N in (2, 4):
    for nm, x in (("Gamma G_1", G_), ("Gamma U_1", U_)):
        Z, _ = regel([gmul(zyk(N), naht(N, 0, x))], [0], "Koszul (Original)")
        anom[f"N = {N}, {nm} (gekleidet)"] = {"modular": Z["ok"], "Bahnen": [(len(O["states"]), O["chars"])
                                                                              for O, _, _ in Z["terms"]]}
Z3k, _ = regel([zyk(3)], [0], "Koszul (Original)")
anom["N = 3, Koszul (Original)"] = {"modular": Z3k["ok"], "Bahnen": [(len(O["states"]), O["chars"]) for O, _, _ in Z3k["terms"]]}
# Spurfunktionen: Gamma(pi) = einseitiges R^2 bei 2 tau, Gamma(pi) G_1 = U bei 2 tau
sp_g = atom_loc(spur_atom((1, 0), ((0, 0), (0, 0)), 0, KOS["Koszul (Original)"]))
sp_r2 = atom_loc([((HA, HA, 0, 0, 2), 1), ((0, 0, 0, 0, 2), -1)])
sp_gg = atom_loc(spur_atom((1, 0), (G_, (0, 0)), 0, KOS["Koszul (Original)"]))
sp_u = atom_loc([((0, 0, 2, 2, 2), 1), ((0, 0, 2, 2, 2), -1)])
r2_n1 = orbit_of(*atom_loc([((HA, HA, 0, 0, 1), 1), ((0, 0, 0, 0, 1), -1)])[:2])
ZAHLEN["anomalie_permutationen"] = anom
ZAHLEN["anomalie_permutationen"]["Spuridentität"] = {"Gamma(pi) = einseitiges R^2 bei 2 tau": sp_g[:2] == sp_r2[:2],
                                                    "Gamma(pi) G_1 = U bei 2 tau": sp_gg[:2] == sp_u[:2],
                                                    "einseitiges R^2 (N = 1) Charaktere": r2_n1["chars"]}
check("[E, entscheidend] Anomalie der Koszul-Permutation: Mit P_CAR auf dem Original (z106) ist die Spur des gradierten Tauschs "
      "die einseitige R^2-Spur bei 2 tau (gleiche Funktion: "
      f"{sp_g[:2] == sp_r2[:2]}), und deren Bahn hat keinen Charakter, wie das einseitige R^2 einer Naht (z113: "
      f"{r2_n1['chars']}): Der gradierte Tausch (N = 2) und das gradierte Z4 (N = 4, 4-Zyklus und Doppeltausch) sind anomal "
      f"(modular: {anom['N = 2, Koszul (Original)']['modular']}, {anom['N = 4, Koszul (Original)']['modular']}). Für N = 3 "
      f"ist der Koszul-Faktor trivial und die Permutation frei ({anom['N = 3, Koszul (Original)']['modular']}). Ungradiert und "
      "mit J-kovarianter CAR (Koszul-Vorzeichen auch auf dem Partner) ist sie frei für N = 2 und 4; bei N = 2 sind beide "
      "dieselbe Bahnsumme",
      sp_g[:2] == sp_r2[:2] and r2_n1["chars"] == [] and not anom["N = 2, Koszul (Original)"]["modular"]
      and not anom["N = 4, Koszul (Original)"]["modular"] and anom["N = 3, Koszul (Original)"]["modular"]
      and all(anom[f"N = {N}, {k}"]["modular"] for N in (2, 4) for k in ("ungradiert", "Koszul (J-kovariant)"))
      and anom["N = 2, ungradiert"]["Bahnen"] == anom["N = 2, Koszul (J-kovariant)"]["Bahnen"])
check("[E] Anomalieaufhebung durch Kleidung: Der mit einem GSO-Faktor gekleidete Tausch Gamma(sigma) G_1 hat die Spur von "
      f"U = G R^2 bei N tau (gleiche Funktion für N = 2: {sp_gg[:2] == sp_u[:2]}), denn K G = R^2 U R^2 = U; er ist für N = 2 "
      f"und 4 anomaliefrei ({anom['N = 2, Gamma G_1 (gekleidet)']['modular']}, "
      f"{anom['N = 4, Gamma G_1 (gekleidet)']['modular']}), ebenso Gamma(sigma) U_1 (Spur G bei N tau; "
      f"{anom['N = 2, Gamma U_1 (gekleidet)']['modular']}, {anom['N = 4, Gamma U_1 (gekleidet)']['modular']}). Die "
      "Koszul-Permutation ist modular nur zusammen mit einer 2 pi-Drehung einer Naht, also als Tausch mit Drehung "
      "(Finkelstein-Rubinstein-Form der Spin-Statistik, Lesart L)",
      sp_gg[:2] == sp_u[:2] and all(anom[f"N = {N}, {nm} (gekleidet)"]["modular"] for N in (2, 4)
                                    for nm in ("Gamma G_1", "Gamma U_1")))
say(f"   {zeit()}")

# --- 1.4 Positivität N = 2
say("   Positivität N = 2 (gemeinsame Uhrverdrillung wie z113; Minimum normierter Spinmode):")
sw = zyk(2)
gG1, gU1 = gmul(sw, naht(2, 0, G_)), gmul(sw, naht(2, 0, U_))
pos_eintrag("N2 diag G = +1 (z113)", 2, [el(range(2), [G_]*2)], [0])
pos_eintrag("N2 diag G = -1 (z113)", 2, [el(range(2), [G_]*2)], [HA])
pos_eintrag("N2 pro Naht G = +1", 2, [naht(2, 0, G_), naht(2, 1, G_)], [0, 0])
pos_eintrag("N2 pro Naht G = -1", 2, [naht(2, 0, G_), naht(2, 1, G_)], [HA, HA])
pos_eintrag("N2 Tausch ungradiert", 2, [sw], [0])
pos_eintrag("N2 pro Naht G = +1 + Tausch ungradiert", 2, [sw, naht(2, 0, G_)], [0, 0])
pos_eintrag("N2 pro Naht G = -1 + Tausch ungradiert", 2, [sw, naht(2, 0, G_)], [0, HA])
pos_eintrag("N2 pro Naht G = +1 + Tausch Koszul (Original)", 2, [sw, naht(2, 0, G_)], [0, 0], "Koszul (Original)")
pos_eintrag("N2 diag G = +1 + Tausch ungradiert", 2, [sw, el(range(2), [G_]*2)], [0, 0])
for k in range(4):
    pos_eintrag(f"N2 gekleidet Gamma(pi) G_1 = i^{k}", 2, [gG1], [Fr(k, 4)], "Koszul (Original)")
for k in range(4):
    pos_eintrag(f"N2 gekleidet Gamma(pi) U_1 = i^{k}", 2, [gU1], [Fr(k, 4)], "Koszul (Original)")
# Kontrolle: getrennte Nativprojektion je Naht = Produkt zweier N = 1-Summen
z1n1, _ = regel([el(range(1), [(0, 0)])], [0])
ref_prod = moden([x**2 for x in Zgrid(z1n1)])
prod_pos = {}
for nm, beta in (("G = +1", 1), ("G = -1", -1)):
    zz, _ = regel([el(range(1), [G_])], [0 if beta == 1 else HA])
    m = moden([x**2 for x in Zgrid(zz)])
    prod_pos[nm] = float(min(m[Y][s].real/abs(ref_prod[Y][0]) for Y in YS for s in m[Y]))
ZAHLEN["kontrolle_getrennte_nativprojektion_N2"] = prod_pos
P2 = POS
check("[E] Reproduktion und Kontrolle: Der neue Aufbau (Kranzprodukt, Multimengen-Bahnen) gibt z113 exakt wieder: N = 2 diagonal "
      f"G = +1 min {P2['N2 diag G = +1 (z113)']['min']:.4g} (z113 4.634e-07), G = -1 {P2['N2 diag G = -1 (z113)']['min']:.6g} "
      "(z113 -77983.3). Mit getrennter Nativprojektion je Naht (Produkt zweier N = 1-Summen) bleibt G = +1 negativ "
      f"({prod_pos['G = +1']:.3g}), G = -1 positiv ({prod_pos['G = -1']:.2g}): Die gemeinsame Uhrverdrillung ist Prämisse "
      "der folgenden N = 2-Ergebnisse",
      abs(P2["N2 diag G = +1 (z113)"]["min"] - 4.633594602e-07) < 1e-12 and abs(P2["N2 diag G = -1 (z113)"]["min"]
                                                                               + 77983.2977462) < 1e-3
      and prod_pos["G = +1"] < -1 and prod_pos["G = -1"] > -1e-8)
check("[E, entscheidend] N = 2 mit summierter 2 pi-Drehung: Pro Naht ist G = +1 (P_GSO auf jeder Naht) bei gemeinsamer "
      f"Uhrverdrillung positiv (min {P2['N2 pro Naht G = +1']['min']:.3g}), ebenso G = -1 ({P2['N2 pro Naht G = -1']['min']:.3g}); "
      "z113 hatte pro Naht nur N = 1 geprüft. Mit P_ident wird es negativ oder anomal: ungradierter (= J-kovarianter) Tausch "
      f"{P2['N2 pro Naht G = +1 + Tausch ungradiert']['min']:.3g}, Koszul-Tausch (P_CAR) ohne modulare Ergänzung. Mit "
      f"Tausch bleiben pro Naht G = -1 ({P2['N2 pro Naht G = -1 + Tausch ungradiert']['min']:.2g}) und diagonal G = +1 "
      f"({P2['N2 diag G = +1 + Tausch ungradiert']['min']:.2g}) positiv",
      P2["N2 pro Naht G = +1"]["positiv"] and P2["N2 pro Naht G = -1"]["positiv"]
      and not P2["N2 pro Naht G = +1 + Tausch ungradiert"]["positiv"]
      and not P2["N2 pro Naht G = +1 + Tausch Koszul (Original)"]["modular"]
      and P2["N2 pro Naht G = -1 + Tausch ungradiert"]["positiv"] and P2["N2 diag G = +1 + Tausch ungradiert"]["positiv"])
check("[E, entscheidend] N = 2 gekleideter Koszul-Tausch: Von den vier Charakteren von <Gamma(pi) G_1> (Z4, Quadrat = G_1 G_2) "
      "ist genau Gamma(pi) G_1 = +1 positiv "
      f"(min {P2['N2 gekleidet Gamma(pi) G_1 = i^0']['min']:.2g}; i^1, i^2, i^3: "
      f"{P2['N2 gekleidet Gamma(pi) G_1 = i^1']['min']:.3g}, {P2['N2 gekleidet Gamma(pi) G_1 = i^2']['min']:.3g}, "
      f"{P2['N2 gekleidet Gamma(pi) G_1 = i^3']['min']:.3g}). Dort gilt auf allen Zuständen G_1 G_2 = +1 und auf den "
      "P_GSO-Zuständen (G_1 = G_2 = 1) Gamma(pi) = +1, also P_ident mit P_CAR (Pauli). Mit der Drehung statt GSO "
      "(Gamma(pi) U_1) ist nur -1 positiv: Identische Paare überleben dann nur mit G = -1",
      [P2[f"N2 gekleidet Gamma(pi) G_1 = i^{k}"]["positiv"] for k in range(4)] == [True, False, False, False]
      and [P2[f"N2 gekleidet Gamma(pi) U_1 = i^{k}"]["positiv"] for k in range(4)] == [False, False, True, False])
say(f"   {zeit()}")

# --- 1.5 Positivität N = 3 und N = 4
say("   Positivität N = 3 und N = 4 (zyklische Permutationen, lange Nähte der Längen 3 bzw. 2 und 4):")
for N in (3, 4):
    pos_eintrag(f"N{N} diag G = +1 (z113)", N, [el(range(N), [G_]*N)], [0])
    pos_eintrag(f"N{N} pro Naht G = +1", N, [naht(N, k, G_) for k in range(N)], [0]*N)
    pos_eintrag(f"N{N} pro Naht G = -1", N, [naht(N, k, G_) for k in range(N)], [HA]*N)
    pos_eintrag(f"N{N} Z{N} ungradiert", N, [zyk(N)], [0])
    pos_eintrag(f"N{N} Z{N} + diag G = +1", N, [zyk(N), el(range(N), [G_]*N)], [0, 0])
    pos_eintrag(f"N{N} Z{N} + pro Naht G = +1", N, [zyk(N), naht(N, 0, G_)], [0, 0])
    pos_eintrag(f"N{N} Z{N} + pro Naht G = -1", N, [zyk(N), naht(N, 0, G_)], [0, HA])
pos_eintrag("N4 Z4 J-kovariant + pro Naht G = +1", 4, [zyk(4), naht(4, 0, G_)], [0, 0], "Koszul (J-kovariant)")
E3 = gmul(zyk(3), naht(3, 0, G_))
for k in range(6):
    pos_eintrag(f"N3 gekleidet Gamma(sigma) G_1 = e^(2 pi i {k}/6)", 3, [E3], [Fr(k, 6)], "Koszul (Original)")
E4, E4u = gmul(zyk(4), naht(4, 0, G_)), gmul(zyk(4), naht(4, 0, U_))
for k in range(8):
    pos_eintrag(f"N4 gekleidet Gamma(sigma) G_1 = e^(2 pi i {k}/8)", 4, [E4], [Fr(k, 8)], "Koszul (Original)")
for k in (0, 4):
    pos_eintrag(f"N4 gekleidet Gamma(sigma) U_1 = e^(2 pi i {k}/8)", 4, [E4u], [Fr(k, 8)], "Koszul (Original)")
P = POS
check("[E] Reproduktion z113 bei N = 3, 4: diagonal G = +1 min "
      f"{P['N3 diag G = +1 (z113)']['min']:.6g} (z113 -18615205.7) und {P['N4 diag G = +1 (z113)']['min']:.6g} (z113 -82.3255)",
      abs(P["N3 diag G = +1 (z113)"]["min"] + 18615205.7377908) < 1 and abs(P["N4 diag G = +1 (z113)"]["min"]
                                                                          + 82.32551862) < 1e-5)
check("[E, entscheidend] P_ident macht G = +1 pro Naht für N = 3 und 4 nicht positiv: pro Naht allein "
      f"{P['N3 pro Naht G = +1']['min']:.3g} bzw. {P['N4 pro Naht G = +1']['min']:.3g}, mit zyklischer Permutation "
      f"{P['N3 Z3 + pro Naht G = +1']['min']:.3g} bzw. {P['N4 Z4 + pro Naht G = +1']['min']:.3g} (ungradiert; N = 3 ist "
      "gradiert dasselbe), J-kovariant bei N = 4 "
      f"{P['N4 Z4 J-kovariant + pro Naht G = +1']['min']:.3g}; die Koszul-Permutation ist bei N = 4 anomal. Dagegen ist die "
      f"vertauschte Statistik G = -1 pro Naht mit Permutation positiv ({P['N3 Z3 + pro Naht G = -1']['min']:.2g}, "
      f"{P['N4 Z4 + pro Naht G = -1']['min']:.2g})",
      not any(P[k]["positiv"] for k in ("N3 pro Naht G = +1", "N4 pro Naht G = +1", "N3 Z3 + pro Naht G = +1",
                                        "N4 Z4 + pro Naht G = +1", "N4 Z4 J-kovariant + pro Naht G = +1"))
      and P["N3 Z3 + pro Naht G = -1"]["positiv"] and P["N4 Z4 + pro Naht G = -1"]["positiv"]
      and P["N3 Z3 ungradiert"]["positiv"] and P["N4 Z4 ungradiert"]["positiv"])
pos4 = [P[f"N4 gekleidet Gamma(sigma) G_1 = e^(2 pi i {k}/8)"]["positiv"] for k in range(8)]
pos3 = [P[f"N3 gekleidet Gamma(sigma) G_1 = e^(2 pi i {k}/6)"]["positiv"] for k in range(6)]
check("[E, entscheidend] Gekleidete Koszul-Permutation bei N = 4: Von den acht Charakteren von <Gamma(sigma) G_1> (Z8, vierte "
      "Potenz = G_1 G_2 G_3 G_4) ist genau Gamma(sigma) G_1 = +1 positiv "
      f"(min {P['N4 gekleidet Gamma(sigma) G_1 = e^(2 pi i 0/8)']['min']:.3g}; positiv: {pos4}). Damit ist die diagonale "
      "GSO-Regel G_d = +1, allein bei N = 4 negativ (z113: -82.3) und mit ungradiertem Z4 ebenfalls "
      f"({P['N4 Z4 + diag G = +1']['min']:.3g}), als geeichtes Element (G_d = (Gamma G_1)^4) auf allen Zuständen erfüllt und "
      "positiv. Bei N = 3 ist dagegen +1 negativ "
      f"({P['N3 gekleidet Gamma(sigma) G_1 = e^(2 pi i 0/6)']['min']:.3g}); positiv sind dort nur die Charaktere "
      f"{[k for k in range(6) if pos3[k]]} (Sechstel), also G_d = +1 nur mit nichttrivialem Z3-Charakter ohne P_ident",
      pos4 == [True] + [False]*7 and not pos3[0] and pos3[2] and pos3[4]
      and not P["N4 Z4 + diag G = +1"]["positiv"])
say(f"   {zeit()}")

# --- 1.6 direkte S/T-Invarianz der Schlüsselsummen
st = {}
cache = {}
for label, (gens, chis, kos) in {"N2 pro Naht G = +1": ([naht(2, 0, G_), naht(2, 1, G_)], [0, 0], "ungradiert"),
                                 "N2 gekleidet Gamma(pi) G_1 = +1": ([gG1], [0], "Koszul (Original)"),
                                 "N4 gekleidet Gamma(sigma) G_1 = +1": ([E4], [0], "Koszul (Original)")}.items():
    Z, _ = regel(gens, chis, kos)
    errS, errT = [], []
    for t in TP:
        z0 = wert(Z, t, cache)
        errS.append(float(abs(wert(Z, -1/t, cache) - z0)/abs(z0)))
        errT.append(float(abs(wert(Z, t + 1, cache) - z0)/abs(z0)))
    st[label] = {"S": max(errS), "T": max(errT)}
# Gegenprobe: nur die Originalsektor-Atome ohne Bahnergänzung
Zk, _ = regel([gG1], [0], "Koszul (Original)")
roh = []
for t in TP:
    s_ = mp.mpc(0)
    for fac, c in [(spur_atom(g[0], g[1], b, KOS["Koszul (Original)"]), 1) for g in huelle([gG1]) for b in range(4)]:
        x = mp.mpc(1)
        for spec, e in fac:
            v = hol(*spec)(t)
            x *= v**e if e > 0 else mp.conj(v)**(-e)
        s_ += x
    roh.append(s_)
roh_S = []
for t in TP:
    s_ = mp.mpc(0)
    for fac, c in [(spur_atom(g[0], g[1], b, KOS["Koszul (Original)"]), 1) for g in huelle([gG1]) for b in range(4)]:
        x = mp.mpc(1)
        for spec, e in fac:
            v = hol(*spec)(-1/t)
            x *= v**e if e > 0 else mp.conj(v)**(-e)
        s_ += x
    roh_S.append(s_)
st["Gegenprobe: Originalsektor ohne Ergänzung, S"] = float(max(abs(a - b)/abs(a) for a, b in zip(roh, roh_S)))
ZAHLEN["direkte_ST_invarianz"] = st
check("[E] Direkte Gegenprobe der Modulinvarianz (unabhängig von der Bahnbuchhaltung, drei Testpunkte): Die Summen pro Naht "
      "G = +1 (N = 2) und gekleidet Gamma G_1 = +1 (N = 2, 4) sind S- und T-invariant (größter relativer Fehler "
      f"{max(max(v['S'], v['T']) for k_, v in st.items() if isinstance(v, dict)):.0e}); der Originalsektor allein ohne "
      f"Bahnergänzung ist es nicht (S-Abweichung {st['Gegenprobe: Originalsektor ohne Ergänzung, S']:.2f})",
      all(max(v["S"], v["T"]) < 1e-15 for v in st.values() if isinstance(v, dict))
      and st["Gegenprobe: Originalsektor ohne Ergänzung, S"] > 1e-2)
befund("Teil 1: Die Koszul-Permutation (P_ident mit P_CAR des Originals) ist für gerade N anomal; ihre Anomalie ist die des "
       "einseitigen R^2 und wird genau durch einen GSO-Faktor einer Naht aufgehoben (Gamma G_1 hat die Spur von U bei N tau). "
       "P_ident macht G = +1 pro Naht für N = 2, 3, 4 nicht positiv. Die gekleidete Koszul-Permutation Gamma G_1 = +1 ist bei "
       "N = 2 und 4 modulinvariant und positiv; dort gilt die Spin-Statistik als Gesamtregel G_d = +1 auf allen Zuständen, und "
       "auf den P_GSO-Zuständen ist die Austauschregel genau P_ident mit P_CAR. Bei N = 3 gelingt das nicht", True)
say(f"   {zeit()}")

# =================================================================== Teil 2
say("\nTeil 2: P_glue unter Modulsummation")

# --- 2.1 Projektoralgebra (exakt)
VALS = [1, 1j, -1, -1j]
ind_ok, quart_ok = True, True
for a_ in range(4):
    for r_ in range(4):
        A, Rv = VALS[a_], VALS[r_]
        U, Gv = A**2, A**2*Rv**2
        pg = (3 + U + Rv**2 + 3*Gv - A*Rv + A**3*Rv + A*Rv**3 - A**3*Rv**3)/8
        soll = 1 if (abs(Gv - 1) < 1e-12 and (abs(U - 1) < 1e-12 or abs(A*Rv + 1) < 1e-12)) else 0
        ind_ok &= abs(pg - soll) < 1e-12
        quart = (-A*Rv + A**3*Rv + A*Rv**3 - A**3*Rv**3)/8
        quart_ok &= abs(quart - (U - 1)*(A*Rv - A*Rv**3)/8) < 1e-12
pa3_ok = all(abs((1 - VALS[a_]*VALS[r_] + VALS[a_]**2*VALS[r_]**2 - (VALS[a_]*VALS[r_])**3)/4
                 - (1 if abs(VALS[a_]*VALS[r_] + 1) < 1e-12 else 0)) < 1e-12 for a_ in range(4) for r_ in range(4))
ZAHLEN["pglue_projektor"] = {"P_glue": "(1/8)[3 + U + R^2 + 3G - AR + A^3R + AR^3 - A^3R^3]",
                             "Viertelanteil": "(1/8)(U - 1)(AR - AR^3)", "R^2-Koeffizient": "1/8 (= U G)",
                             "P_A3": "(1/4)[1 - AR + G - A^3R^3]"}
check("[E] Projektoralgebra (exakt auf allen 16 Eigenwertpaaren von A = i^(N - Nquer) und R): Mit Z_SU = -A ist P_glue "
      "(Z_SU R = 1 auf dem Spinormodul U = -1, G = 1 auf den Bosonen) = (1/8)[3 + U + R^2 + 3G - AR + A^3R + AR^3 - A^3R^3] "
      f"(Indikator {ind_ok}), P_A3 = (1/4)[1 - AR + G - A^3R^3] ({pa3_ok}). Der Vierteldrehungsanteil von P_glue ist "
      f"(1/8)(U - 1)(AR - AR^3) ({quart_ok}), er lebt nur auf dem Spinormodul; außerdem enthält P_glue das einseitige "
      "R^2 = U G mit Koeffizient 1/8",
      ind_ok and quart_ok and pa3_ok)

# --- 2.2 Bahnen der Vierteldrehungen
viertel = {}
for N in (1, 2, 4):
    rows = {}
    for a_ in range(4):
        for r_ in range(4):
            for b in range(4):
                sig, stt, ph0 = atom_loc(spur_atom(tuple(range(N)), tuple([(a_, r_)]*N), b, KOS["ungradiert"]))
                O = orbit_of(sig, stt)
                rows[(a_, r_, b)] = (O, O["idx"][stt], ph0)
    paare, anomal = [], []
    for r_ in range(4):
        for b in range(4):
            O1, j1, p1 = rows[(1, r_, b)]
            O3, j3, p3 = rows[(3, r_, b)]
            if O1["V"] is None:
                anomal.append((r_, b))
                continue
            rel = (O3["V"][j3] - p3 - (O1["V"][j1] - p1)) % 1 if O1 is O3 else None
            paare.append((r_, b, O1 is O3, str(rel)))
    viertel[N] = {"anomal (Uhrpotenz r, nativ b)": sorted({r for r, b in anomal}),
                  "A R^r und A^3 R^r in einer Bahn, rel. Koeffizient": sorted({(r, s, c) for r, b, s, c in paare}),
                  "Bahngröße": rows[(1, 0, 0)][0]["states"].__len__(),
                  "Charaktere A": rows[(1, 0, 0)][0]["chars"]}
ZAHLEN["vierteldrehung_bahnen"] = {N: {k: str(v) for k, v in d.items()} for N, d in viertel.items()}
paarung = all(all(s and c == "0" for r, s, c in viertel[N]["A R^r und A^3 R^r in einer Bahn, rel. Koeffizient"])
              for N in (1, 2, 4))
check("[E, entscheidend] Bahnen der Vierteldrehung (96 Elemente, lange Nähte bis Länge 16): Für jede Uhrpotenz R^r und jeden "
      "Nativeinschub liegen A R^r und A^3 R^r = U A R^r in derselben SL(2,Z)-Bahn mit relativem Koeffizienten exakt 1 "
      f"(N = 1, 2, 4: {paarung}). Anomal sind bei N = 1 die ungeraden Uhrpotenzen "
      f"{viertel[1]['anomal (Uhrpotenz r, nativ b)']}, bei N = 2 und 4 keine "
      f"({viertel[2]['anomal (Uhrpotenz r, nativ b)']}, {viertel[4]['anomal (Uhrpotenz r, nativ b)']}). Folge: Jede modulare "
      "Ergänzung eines Vierteldrehungsterms ist proportional zu (1 + U) und verschwindet auf dem Spinormodul",
      paarung and viertel[1]["anomal (Uhrpotenz r, nativ b)"] == [1, 3] and viertel[2]["anomal (Uhrpotenz r, nativ b)"] == []
      and viertel[4]["anomal (Uhrpotenz r, nativ b)"] == [] and viertel[1]["Bahngröße"] == 96)

# --- 2.3 P_glue und P_A3 unter Modulsummation
PGLUE = {(0, 0): 3, (2, 0): 1, (0, 2): 1, (2, 2): 3, (1, 1): -1, (3, 1): 1, (1, 3): 1, (3, 3): -1}
PA3 = {(0, 0): 1, (1, 1): -1, (2, 2): 1, (3, 3): -1}
PGSO = {(0, 0): 1, (2, 2): 1}


def diag_regel(N, coeffs):
    return ergaenze([(spur_atom(tuple(range(N)), tuple([x]*N), b, KOS["ungradiert"]), mp.mpc(c))
                     for x, c in coeffs.items() for b in range(4)])


def ergaenzter_projektor(N, coeffs):
    """Originalsektor-Operator nach Bahnergänzung (alle Bahnmitglieder A^a R^r n^b mit erzwungenen Koeffizienten);
    kleinster Eigenwert auf den 64 Eigenwerttripeln (A, R, n), normiert wie der Projektor."""
    rows = {}
    for a_ in range(4):
        for r_ in range(4):
            for b in range(4):
                sig, stt, ph0 = atom_loc(spur_atom(tuple(range(N)), tuple([(a_, r_)]*N), b, KOS["ungradiert"]))
                O = orbit_of(sig, stt)
                rows[(a_, r_, b)] = (O, O["idx"][stt], ph0)
    c_full = defaultdict(mp.mpc)
    done = set()
    for (a_, r_), c in coeffs.items():
        O0, j0, p0 = rows[(a_, r_, 0)]
        if id(O0) in done or O0["V"] is None:
            continue
        done.add(id(O0))
        for k, (O, j, ph) in rows.items():
            if O is O0:
                c_full[k] += mp.mpc(c)*EZ(O["V"][j] - ph - (O0["V"][j0] - p0))
    norm = 4*len(coeffs)                   # Projektor (1/|Gruppe|) sum über die Gruppe, dazu die vier Nativeinschübe
    ev = []
    for al in range(4):
        for rh in range(4):
            for nu in range(4):
                ev.append(mp.fsum(c*VALS[(k[0]*al + k[1]*rh + k[2]*nu) % 4] for k, c in c_full.items())/norm)
    return c_full, ev


pglue = {}
for N in (1, 2, 3, 4):
    Z = diag_regel(N, PGLUE)
    grund = []
    for O, lam, kons in Z["terms"]:
        if O["V"] is None:
            grund.append(f"anomal (Bahn {len(O['states'])}, Charaktere {O['chars']})")
        elif not kons:
            grund.append(f"inkonsistent (Bahn {len(O['states'])})")
    pglue[N] = {"modular": Z["ok"], "Gründe": grund}
pa3 = {}
for N in (2, 4):
    cf, ev = ergaenzter_projektor(N, PA3)
    _, evg = ergaenzter_projektor(N, PGSO)
    cfm, evm = ergaenzter_projektor(N, {(0, 0): 1, (1, 1): 1, (2, 2): 1, (3, 3): 1})
    pa3[N] = {"P_A3 ergänzt: min Eigenwert": float(min(x.real for x in ev)),
              "P_A3 ergänzt: A^3R mit": mp.nstr(cf[(3, 1, 0)], 6),
              "Spiegel ergänzt: min Eigenwert": float(min(x.real for x in evm)),
              "P_GSO ergänzt: min Eigenwert": float(min(x.real for x in evg)),
              "P_GSO ergänzt: max Eigenwert": float(max(x.real for x in evg))}
pa3_n13 = {N: diag_regel(N, PA3)["ok"] for N in (1, 3)}
ZAHLEN["pglue_modular"] = pglue
ZAHLEN["pa3_ergaenzung (Originalsektor, normiert)"] = pa3
ZAHLEN["pa3_modular_ungerade_N"] = pa3_n13
check("[E, entscheidend] P_glue hat für N = 1, 2, 3, 4 keine modulare Ergänzung "
      f"({[pglue[N]['modular'] for N in (1, 2, 3, 4)]}): N = 1, 3 anomal (Vierteldrehung mit ungerader Uhrpotenz und "
      "einseitiges R^2), N = 2 einseitiges R^2 nur mit Charakter (1/2, 1/2) und die Vierteldrehungsbahnen inkonsistent, N = 4 "
      f"inkonsistent ({pglue[4]['Gründe']}): Die Vorzeichen -AR, +A^3R verlangen den relativen Koeffizienten -1, die Bahn "
      "erzwingt +1. Der definierende Anteil (U - 1)(AR - AR^3) ist genau der modular unsichtbare",
      not any(pglue[N]["modular"] for N in (1, 2, 3, 4)) and any("inkonsistent" in g for g in pglue[4]["Gründe"])
      and any("anomal" in g for g in pglue[1]["Gründe"]))
check("[E] P_A3 (die Z4-Regel ohne Einschränkung auf Spinoren, AR = -1): für ungerade N anomal "
      f"({pa3_n13}); für N = 2 und 4 ergänzt die Bahn A^3R mit demselben Koeffizienten wie AR "
      f"({pa3[2]['P_A3 ergänzt: A^3R mit']}), der Originalsektor wird (1/4)[1 + G - (1 + U)(AR + A^3R^3)] mit kleinstem "
      f"Eigenwert {pa3[2]['P_A3 ergänzt: min Eigenwert']:.3g} (N = 2) bzw. {pa3[4]['P_A3 ergänzt: min Eigenwert']:.3g} (N = 4) "
      "auf den neutralen Bosonen A = R = 1: negative Vielfachheiten schon im Originalsektor. Der Spiegel AR = +1 ebenso "
      f"({pa3[2]['Spiegel ergänzt: min Eigenwert']:.3g}). Kontrolle: P_GSO ergänzt sich zum Projektor (Eigenwerte "
      f"{pa3[2]['P_GSO ergänzt: min Eigenwert']:.3g} bis {pa3[2]['P_GSO ergänzt: max Eigenwert']:.3g})",
      not any(pa3_n13.values()) and all(pa3[N]["P_A3 ergänzt: min Eigenwert"] < -0.1 and pa3[N]["Spiegel ergänzt: min Eigenwert"]
                                        < -0.1 and abs(pa3[N]["P_GSO ergänzt: min Eigenwert"]) < 1e-12
                                        and abs(pa3[N]["P_GSO ergänzt: max Eigenwert"] - 1) < 1e-12 for N in (2, 4)))
befund("Teil 2: Die Einschränkung auf den Spinormodul umgeht z113s Hindernis nicht, sie verschärft es. Modulsummation paart jede "
       "Vierteldrehung A X mit U A X (Koeffizient 1), sieht also auf dem Spinormodul nur U und R^2, wie bei P_GSO. Der "
       "chirale Z4-Inhalt von P_glue liegt genau dort, wo die Modulsummation nichts sieht; P_glue hat für N = 1..4 keine "
       "modulare Ergänzung, P_A3 ergänzt sich nur zu einem Nicht-Projektor", True)
say(f"   {zeit()}")

# =================================================================== Teil 3
say("\nTeil 3: Klassifikation der modulverträglichen Statistikregeln der Naht")
lifts = {}
for N in (1, 2, 3, 4):
    zeile = {}
    for nm, x in (("U", U_), ("G = U R^2", G_)):
        for beta, bn in ((0, "+1"), (HA, "-1")):
            lab = f"N{N} diag {nm} = {bn}"
            if lab not in POS:
                pos_eintrag(lab, N, [el(range(N), [x]*N)], [beta])
            zeile[f"{nm} = {bn}"] = POS[lab].get("positiv") if POS[lab]["modular"] else "anomal"
    Zur = diag_regel(N, {(0, 0): 1, (2, 1): 1, (0, 2): 1, (2, 3): 1})
    zeile["U R (Z4-Hebung, Quadrat R^2)"] = "modular" if Zur["ok"] else "anomal"
    Za = diag_regel(N, {(0, 0): 1, (1, 0): 1, (2, 0): 1, (3, 0): 1})
    Za_i = diag_regel(N, {(0, 0): 1, (1, 0): -1j, (2, 0): -1, (3, 0): 1j})
    zeile["A = +1 (Vierteldrehung, nur Bosonen)"] = "modular" if Za["ok"] else "anomal"
    zeile["A = +i (Vierteldrehung auf Spinoren)"] = "modular" if Za_i["ok"] else "nicht ergänzbar"
    zeile["P_A3"] = ("anomal" if N % 2 else "Nicht-Projektor")
    zeile["P_glue"] = "keine Ergänzung"
    lifts[N] = zeile
ZAHLEN["klassifikation_hebungen_diagonal"] = lifts
for N, z in lifts.items():
    say(f"      N = {N}: " + "; ".join(f"{k}: {v}" for k, v in z.items()))
check("[E] Klassifikation der Hebungen auf N Nähten (diagonal): Die unverdrehte Hebung ist für N = 1..4 nur mit U = +1 "
      f"positiv ({[lifts[N]['U = +1'] for N in (1, 2, 3, 4)]}, U = -1: {[lifts[N]['U = -1'] for N in (1, 2, 3, 4)]}); "
      f"P_GSO diagonal G = +1: {[lifts[N]['G = U R^2 = +1'] for N in (1, 2, 3, 4)]}, G = -1: "
      f"{[lifts[N]['G = U R^2 = -1'] for N in (1, 2, 3, 4)]}; die Z4-Hebung U R (Quadrat = einseitiges R^2) ist erst bei "
      f"N = 4 modular ({[lifts[N]['U R (Z4-Hebung, Quadrat R^2)'] for N in (1, 2, 3, 4)]}); die Vierteldrehung allein ist "
      "nur mit reellem Eigenwert (A = +-1, also U = +1, keine Spinoren) ergänzbar, mit A = +-i (auf Spinoren) nie",
      [lifts[N]["U = +1"] for N in (1, 2, 3, 4)] == [True]*4 and [lifts[N]["U = -1"] for N in (1, 2, 3, 4)] == [False]*4
      and [lifts[N]["G = U R^2 = +1"] for N in (1, 2, 3, 4)] == [False, True, False, False]
      and [lifts[N]["G = U R^2 = -1"] for N in (1, 2, 3, 4)] == [True, False, False, False]
      and [lifts[N]["U R (Z4-Hebung, Quadrat R^2)"] for N in (1, 2, 3, 4)] == ["anomal"]*3 + ["modular"]
      and all(lifts[N]["A = +1 (Vierteldrehung, nur Bosonen)"] == "modular"
              and lifts[N]["A = +i (Vierteldrehung auf Spinoren)"] == "nicht ergänzbar" for N in (1, 2, 3, 4)))

TAB = {
    "unverdrehte Hebung U = +1": "positiv N = 1..4, keine Spinoren",
    "P_GSO diagonal (z113)": "positiv nur N = 2",
    "P_GSO pro Naht, gemeinsame Uhr": "positiv nur N = 2 (beide Vorzeichen); N = 3, 4 negativ",
    "P_GSO pro Naht + P_ident": "N = 2, 3, 4 negativ (ungradiert, J-kovariant) oder anomal (Koszul)",
    "vertauschte Statistik G = -1 pro Naht + P_ident": "positiv N = 2, 3, 4",
    "gekleidete Koszul-Permutation Gamma G_1 = +1": "positiv N = 2, 4 (G_d = +1 auf allen Zuständen); N = 3 negativ",
    "Z4-Hebung U R": "anomal N = 1, 2, 3",
    "P_A3": "anomal (ungerade N), Nicht-Projektor (gerade N)",
    "P_glue": "keine modulare Ergänzung N = 1..4",
}
ZAHLEN["klassifikation_gesamt"] = TAB
for k_, v_ in TAB.items():
    say(f"   [E] {k_}: {v_}")
befund("Teil 3: Spin-Statistik (2 pi-Drehung = (-1)^F) ist im Nahtrahmen modulverträglich, aber nur als Gesamtregel und nur "
       "für gerade N: bei N = 2 diagonal oder pro Naht ohne P_ident, bei N = 2 und 4 mit P_ident als gekleidete "
       "Koszul-Permutation Gamma G_1 = +1. Pro Naht zusammen mit P_ident ist sie für N = 2, 3, 4 ausgeschlossen. Ihre chirale "
       "Z4-Verfeinerung auf Spinoren (P_glue, P_A3) ist modular nie verträglich", True)
say(f"   {zeit()}")

STATUS = {
    "Koszul-Spuridentität, gaußscher Zyklus": ("E", "Tr Gamma(sigma) = Tr K^(N-1) q^(N L_0) = Produkt der Fouriermoden"),
    "Koszul-Permutation (P_CAR Original)": ("E (negativ)", "anomal für gerade N (Anomalie des einseitigen R^2)"),
    "J-kovariante CAR": ("E", "frei; bei N = 2 gleich dem ungradierten Tausch"),
    "gekleideter Tausch Gamma G_1": ("E", "anomaliefrei (Spur von U bei N tau); +1 positiv bei N = 2, 4, negativ bei N = 3"),
    "P_GSO pro Naht": ("E", "positiv nur bei N = 2 mit gemeinsamer Uhr; mit P_ident nie (N = 2, 3, 4)"),
    "Vierteldrehung": ("E", "A X und U A X in einer Bahn (Koeffizient 1): modular unsichtbar auf Spinoren"),
    "P_glue": ("E (negativ)", "keine modulare Ergänzung für N = 1..4"),
    "P_A3": ("E (negativ)", "anomal (ungerade N), ergänzt zu Nicht-Projektor (gerade N)"),
    "Spin-Statistik modular": ("C", "als Gesamtregel bei N = 2, 4 unter P_dress (gekleidete Koszul-Permutation) und "
                                    "gemeinsamer Uhrverdrillung; pro Naht mit P_ident nein"),
    "Lesart Tausch mit Drehung": ("L", "Gamma G_1 ist auf P_GSO-Zuständen der Koszul-Tausch; die Kleidung trägt eine "
                                       "2 pi-Drehung (Finkelstein-Rubinstein)"),
    "N >= 5, N = 6, 8": ("O", "nicht gerechnet"),
}
for k_, (s_, t_) in STATUS.items():
    say(f"   [{s_}] {k_}: {t_}")
ZAHLEN["status"] = STATUS
ZAHLEN["positivitaet (min normierter Spinmode, Y = 0.35, 0.6, 1.0, |s| <= 5)"] = POS
ZAHLEN["urteil_kurz"] = {
    "Teil 1": "Koszul-Permutation anomal (gerade N); gekleidet Gamma G_1 = +1 modular und positiv bei N = 2, 4; P_GSO pro "
              "Naht mit P_ident nie positiv",
    "Teil 2": "P_glue ohne modulare Ergänzung (N = 1..4); Vierteldrehung auf Spinoren modular unsichtbar",
    "Teil 3": "Spin-Statistik nur als Gesamtregel, gerade N; Spinor-Z4 nie",
    "modular_statistics": "P_GSO bedingt wahr unter P_dress (N = 2, 4, gemeinsame Uhr); P_glue, P_A3 falsch (jetzt "
                          "geprüft); unverdrehte Hebung wahr (ohne Spinoren)"}

urteil = (
    "Teil 1 (gradierte Permutationsorbifaltigkeit): Die Koszul-Permutation mit P_CAR des Originals (z106) ist gaußsch: Der "
    "gradierte N-Zyklus wirkt auf den Fouriermoden mit e^(2 pi i k/N), seine Spur ist exakt das Produkt der verdrillten "
    "Modenspuren, seine langen Nähte sind verdrillte Fouriermoden. Für gerade N ist er anomal. Seine Spur ist das einseitige "
    "R^2 bei N tau, dessen Bahn keinen Charakter hat (z113). Für N = 3 ist der Koszul-Faktor trivial. Mit J-kovarianter CAR "
    "(Partner mitgradiert) ist die Permutation frei, bei N = 2 aber identisch mit dem ungradierten Tausch. Die Anomalie wird "
    "genau durch einen GSO-Faktor einer Naht aufgehoben: Gamma(sigma) G_1 hat die Spur von U bei N tau. Ergebnis der "
    "Positivität bei gemeinsamer Uhrverdrillung (z113-Rahmen): P_GSO pro Naht ist bei N = 2 positiv (neu; z113 hatte pro Naht "
    "nur N = 1), bei N = 3, 4 negativ. Zusammen mit P_ident ist P_GSO pro Naht für N = 2, 3, 4 negativ (ungradiert, "
    "J-kovariant) oder anomal (Koszul). Die vertauschte Statistik G = -1 pro Naht bleibt mit P_ident positiv. Die gekleidete "
    "Koszul-Permutation mit Gamma G_1 = +1 ist bei N = 2 und N = 4 modulinvariant (auch direkt geprüft) und positiv; als "
    "einziger ihrer Charaktere. Dort ist G_d = (Gamma G_1)^N geeicht, also gilt die Spin-Statistik als Gesamtregel auf allen "
    "Zuständen. Auf den P_GSO-Zuständen ist die Austauschregel genau P_ident mit P_CAR (Pauli). Bei N = 4 rettet das die "
    "diagonale GSO-Regel, die allein (z113: -82,3) und mit ungradiertem Z4 negativ ist. Bei N = 3 ist +1 negativ. Teil 2 "
    "(P_glue): P_glue = (1/8)[3 + U + R^2 + 3G - AR + A^3R + AR^3 - A^3R^3] mit A = -Z_SU = i^(N - Nquer). Sein chiraler "
    "Anteil (1/8)(U - 1)(AR - AR^3) lebt nur auf dem Spinormodul. Modulsummation legt dagegen jede Vierteldrehung A X mit "
    "U A X in eine Bahn mit Koeffizient 1, also ist jede modulare Ergänzung eines Vierteldrehungsterms proportional zu "
    "(1 + U) und verschwindet auf Spinoren. Darum hat P_glue für N = 1..4 keine modulare Ergänzung (N = 1, 3 anomal, N = 2 "
    "R^2 nur mit Charakter (1/2, 1/2), N = 4 inkonsistent). P_A3 ist für ungerade N anomal und ergänzt sich für gerade N nur "
    "zu (1/4)[1 + G - (1 + U)(AR + A^3R^3)] mit Eigenwert -1/2 auf den neutralen Bosonen. Die Einschränkung auf den "
    "Spinormodul umgeht z113s Hindernis also nicht: Auf Spinoren sieht die Modulsummation nur U und R^2, wie bei P_GSO. "
    "Teil 3 (Klassifikation): Modulinvariant und positiv sind die unverdrehte Hebung (U = +1, N = 1..4, keine Spinoren), "
    "P_GSO diagonal nur bei N = 2, P_GSO pro Naht nur bei N = 2 ohne P_ident, die vertauschte Statistik pro Naht mit P_ident "
    "(N = 2, 3, 4) und die gekleidete Koszul-Permutation Gamma G_1 = +1 (N = 2, 4). Die Z4-Hebung U R ist erst bei N = 4 "
    "modular. Die Vierteldrehung ist nur auf Bosonen ergänzbar; P_A3 und P_glue nie. Spin-Statistik ist im Nahtrahmen also "
    "modulverträglich, aber nur als Gesamtregel und nur für gerade N (gerechnet N = 2, 4). Mit identischen Nähten braucht "
    "sie die Prämisse P_dress: Nähte werden über den mit einem GSO-Faktor gekleideten Koszul-Tausch komponiert, bei "
    "gemeinsamer Uhrverdrillung. Pro Naht zusammen mit P_ident ist sie ausgeschlossen. Folge für die Szenario-Klausel "
    "modular_statistics: Die Statistik P_GSO kann sie bedingt erhalten, unter P_dress mit N = 2 oder 4. P_glue und P_A3 "
    "erhalten sie nicht; das ist jetzt geprüft und negativ. Die unverdrehte Hebung behält sie, ohne Spinoren. Die Chiralität "
    "muss deshalb aus der Indexroute kommen (z114), nicht aus der Z4-Verklebung. Nichts ist befördert.")
say("\nUrteil\n   " + urteil)
say(f"\nAlle Prüfungen bestanden: {all(CHECKS.values())} ({sum(CHECKS.values())}/{len(CHECKS)})")
(HERE / "z126_ergebnis.json").write_text(json.dumps(
    {"checks": CHECKS, "befunde": BEFUNDE, "zahlen": ZAHLEN, "urteil": urteil}, indent=1, ensure_ascii=False,
    default=str) + "\n")
(HERE / "ausgabe_z126.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
