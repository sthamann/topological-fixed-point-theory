"""Z3: Prämisse A herleiten — die Naht ohne Hintergrundgeometrie.

Bisher war A gesetzt: Der orientierte Tausch der beiden Nahtzyklen ist eine Symmetrie (Fixpunkt).
Hier wird geprüft, ob A aus einer allgemeineren Forderung folgt.

Hintergrundfreiheit der Naht: Ihre Geometrie ist nicht von außen vorgegeben.
 (a) Invarianz unter allen Diffeomorphismen, auch den großen; das sind die Basiswechsel SL(2,Z).
     Daraus folgt die Summe über Spinstrukturen, also E8 mit T-Phase (Z1, Prämisse F).
 (b) Die Form ist stationär: Die Antwort des Zustands auf jede Deformation der komplexen Struktur
     verschwindet, <T_zz> = 0. Das ist die Bewegungsgleichung einer dynamischen Nahtmetrik
     (wie die Virasoro-Bedingung).

Für die chirale E8 Naht ist Z = chi = E4/eta^8, und d log chi/d tau = -(2 pi i/3) E6/E4. Geprüft wird:
 1. die Identität;
 2. die stationären Formen im ganzen Fundamentalbereich;
 3. die Energie beider Rekonstruktionen am stationären Punkt;
 4. die Stabilität (log|chi| ist harmonisch);
 5. Gegenproben mit anderem Nahtinhalt.
"""
import json
from pathlib import Path

import numpy as np
from mpmath import mp, mpc, mpf, exp, pi, qp, jtheta, findroot, diff, log, sqrt, fabs

from scipy.optimize import root as sroot

mp.dps = 20
OK = []
NN = np.arange(1, 80)


def nq(t):
    return np.exp(2j * np.pi * np.asarray(t, complex))


def neis(k, t):
    c = {2: -24, 4: 240, 6: -504}[k]
    qq = nq(t)[..., None] ** NN
    return 1 + c * np.sum(NN ** (k - 1) * qq / (1 - qq), axis=-1)


def neta(t):
    t = np.asarray(t, complex)
    return np.exp(2j * np.pi * t / 24) * np.prod(1 - nq(t)[..., None] ** np.arange(1, 200), axis=-1)


in_sl2 = lambda r: r.imag > 0.5 and abs(r.real) <= 0.5 + 1e-7 and abs(r) >= 1 - 1e-7


def zeros_in(f, domain, xr, yr, n=161):
    """Zeros of a (possibly non-holomorphic) complex function: grid minima of |f|, then 2D Newton."""
    X, Y = np.meshgrid(np.linspace(*xr, n), np.linspace(*yr, n))
    T = X + 1j * Y
    A = np.abs(f(T))
    found = []
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            if A[i, j] <= A[i - 1:i + 2, j - 1:j + 2].min():
                sol = sroot(lambda v: [f(v[0] + 1j * v[1]).real, f(v[0] + 1j * v[1]).imag], [X[i, j], Y[i, j]], tol=1e-14)
                z = sol.x[0] + 1j * sol.x[1]
                if abs(f(z)) < 1e-9 and domain(z):
                    if abs(abs(z.real) - 0.5) < 1e-6 and domain is in_sl2:
                        z = 0.5 + 1j * z.imag            # x = -1/2 and +1/2 are the same torus
                    if not any(abs(z - w) < 1e-6 for w in found):
                        found.append(z)
    return sorted((round(z.real, 6) + 0.0, round(z.imag, 6)) for z in found)


in_theta = lambda r: r.imag > 0.3 and abs(r.real) <= 1 + 1e-7 and abs(r) >= 1 - 1e-7


def check(name, cond):
    OK.append((name, bool(cond)))
    return "ja" if cond else "NEIN"


def q(tau):
    return exp(2j * pi * tau)


def eta(tau):
    return exp(1j * pi * tau / 12) * qp(q(tau))


def eis(k, tau, terms=40):
    c = {2: -24, 4: 240, 6: -504}[k]
    s = sum(n ** (k - 1) * q(tau) ** n / (1 - q(tau) ** n) for n in range(1, terms))
    return 1 + c * s


def chi(tau):
    return eis(4, tau) / eta(tau) ** 8


print("1. Identität d log chi/d tau = -(2 pi i/3) E6/E4\n")
for tau in (mpc(0.13, 1.21), mpc(-0.31, 0.97), mpc(0.4, 2.2)):
    lhs = diff(lambda t: log(chi(t)), tau)
    rhs = -(2j * pi / 3) * eis(6, tau) / eis(4, tau)
    print(f"   tau = {mp.nstr(tau, 3)}: |lhs - rhs| = {mp.nstr(abs(lhs - rhs), 3)}")
check("Identität", abs(lhs - rhs) < 1e-15)

print("\n2. Stationäre Formen: Nullstellen von E6 (und Pole von log chi bei E4 = 0) im Fundamentalbereich\n")
roots = zeros_in(lambda t: neis(6, t), in_sl2, (-0.6, 0.6), (0.8, 2.6))
roots = [r for r in roots if abs(eis(6, mpc(*r))) < 1e-5]
print(f"   Nullstellen von E6 im Fundamentalbereich: {sorted(roots)}")
check("nur tau = i stationär", roots == [(0.0, 1.0)])
rho = mpc(0.5, mp.sqrt(3) / 2)
print(f"   E4(rho) = {mp.nstr(abs(eis(4, rho)), 3)}: dort ist chi = 0, keine zulässige Form (Zustandssumme null)")
print("   Valenzformel: E6 hat Gewicht 6, also 6/12 = 1/2 Nullstelle, genau die halbe Nullstelle am")
print("   elliptischen Punkt i der Ordnung 2. Die Form ist damit ohne weitere Forderung eindeutig.")

print("\n3. Energie beider Rekonstruktionen (Casimir plus thermisch)\n")
for y in (mpf("0.8"), mpf(1), mpf("1.25")):
    # reconstruction A: space length 1, beta = y; <H> = -d/dbeta log chi(i beta)
    E_A = -diff(lambda b: log(chi(mpc(0, b))).real, y)
    E_B = -diff(lambda b: log(chi(mpc(0, b))).real, 1 / y)
    print(f"   beta/L = {mp.nstr(y, 3):>4}: <H>_A = {float(E_A):+.6f}, <H>_B (Tausch) = {float(E_B):+.6f}")
E1 = -diff(lambda b: log(chi(mpc(0, b))).real, mpf(1))
check("Energie null bei i", abs(E1) < 1e-15)
print("   Am stationären Punkt ist die Gesamtenergie exakt null: Die thermische Energie hebt die")
print("   Casimirenergie auf. Beide Rekonstruktionen sind dasselbe System (beta = L), Zeitneutralität folgt.")

print("\n4. Stabilität: F = -log|chi| ist harmonisch\n")
F = lambda x, y: -log(fabs(chi(mpc(x, y))))
h = mpf("1e-4")
Fxx = (F(h, 1) - 2 * F(0, 1) + F(-h, 1)) / h ** 2
Fyy = (F(0, 1 + h) - 2 * F(0, 1) + F(0, 1 - h)) / h ** 2
print(f"   Hesse am Punkt i: F_xx = {float(Fxx):+.5f}, F_yy = {float(Fyy):+.5f}, Spur = {float(Fxx + Fyy):+.1e}")
check("Sattel", Fxx > 0 > Fyy and abs(Fxx + Fyy) < 1e-4)
print("   Sattel: stabil gegen Verdrillung, instabil gegen Streckung (Casimir zieht zu dünnen Tori).")
print("   Allgemein: Für eine chirale Naht ist log|Z| harmonisch und hat nirgends ein lokales Minimum.")
print("   Energieminimierung kann die Form also grundsätzlich nicht wählen; Stationarität schon.")

print("\n5. Gegenproben: stationäre Formen bei anderem Nahtinhalt\n")


bos = zeros_in(lambda t: neis(2, t) - 3 / (np.pi * np.imag(t)), in_sl2, (-0.6, 0.6), (0.8, 2.6))
print(f"   freies Boson (E2* = 0): {bos}")
kk = np.arange(-12, 13)


def dlog_th3_eta(t):
    t = np.asarray(t, complex)
    e = np.exp(1j * np.pi * kk ** 2 * t[..., None])
    return np.sum(1j * np.pi * kk ** 2 * e, axis=-1) / np.sum(e, axis=-1) - (2j * np.pi / 24) * neis(2, t)


dirac = zeros_in(dlog_th3_eta, in_theta, (-1, 1), (0.5, 2.6), n=201)
print(f"   Dirac mit fester Spinstruktur (A,A), Thetagruppe: {dirac}")
e83 = zeros_in(lambda t: neis(6, t) / neis(4, t), in_sl2, (-0.6, 0.6), (0.8, 2.6))
jfn = lambda t: neis(4, t) ** 3 / neta(t) ** 24
leech = zeros_in(lambda t: neis(6, t) * neis(4, t) ** 2 / (neis(4, t) ** 3 - 720 * neta(t) ** 24), in_sl2, (-0.6, 0.6), (0.8, 2.6))
print(f"   E8 hoch 3 (Z = j): {e83}")
print(f"   Leech (Z = j - 720): {leech}")
check("E8^3 nur i", e83 == [(0.0, 1.0)])
check("Leech auch rho", any(abs(p[1] - 0.866025) < 1e-5 for p in leech))
check("Boson auch rho", any(abs(p[1] - 0.866025) < 1e-5 for p in bos))
print("   Nur wo die Zustandssumme am Sechseck verschwindet (E8, E8 hoch 3), ist i die einzige stationäre")
print("   Form. Beim Boson und beim Leech-Gitter ist auch das Sechseck stationär; dort entscheidet erst die")
print("   Spinstruktur (Sechseck hält nur PP mit Nullmoden, Z1).")

print("\nZusammenfassung")
for n, c in OK:
    print(f"   {'bestanden' if c else 'NICHT bestanden'}: {n}")
Path(__file__).with_name("z3_ergebnis.json").write_text(json.dumps({
    "checks": {n: c for n, c in OK}, "e6_zeros": sorted(roots), "boson_critical": bos,
    "dirac_AA_critical": dirac, "e8cubed_critical": e83, "leech_critical": leech,
    "hessian_at_i": [float(Fxx), float(Fyy)]}, indent=1))
assert all(c for _, c in OK)
