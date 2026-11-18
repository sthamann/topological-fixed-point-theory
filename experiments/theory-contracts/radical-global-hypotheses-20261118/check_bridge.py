#!/usr/bin/env python3
"""RADICALBRIDGE20261118: Ist c3 auf der Sphaere (Punktwert) und auf dem Torus (Energie) dieselbe Groesse?
Reproduktion: python3 -B check_bridge.py"""
import json, sys
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi; c3 = 1/(8*pi); R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)
def eis(tau):
    q = mp.e**(2j*pi*tau)
    return (1-24*mp.nsum(lambda n: n*q**n/(1-q**n), [1, mp.inf]), 1+240*mp.nsum(lambda n: n**3*q**n/(1-q**n), [1, mp.inf]),
            1-504*mp.nsum(lambda n: n**5*q**n/(1-q**n), [1, mp.inf]))
# H1  weight lemma: f(1/t) = t^k f(t)  =>  -(1/2pi) t d/dt log f at t=1 equals k/(4 pi)
th3 = lambda t: mp.jtheta(3, 0, mp.e**(-pi*t))            # sum exp(-pi t n^2), weight 1/2
E = lambda f, t: -(1/(2*pi))*t*mp.diff(lambda s: mp.log(f(s)), t)
check("H1 one lattice direction (theta3, weight 1/2): energy at t=1 is (1/2)/(4 pi) = c3", abs(E(th3, 1)-c3) < mp.mpf("1e-20"))
check("H2 rank 24 self-dual (theta3^24): energy per direction at t=1 is c3", abs(E(lambda t: th3(t)**24, 1)/24-c3) < mp.mpf("1e-20"))
# H3  exact decomposition per Cartan direction of E8:  E = 1/(8 pi Im tau) + R(tau), R = E2*/24 - E6/(24 E4), R(i)=0
def parts(tau):
    e2, e4, e6 = eis(tau); t = mp.im(tau)
    total = e2/24 - e6/(24*e4)
    return total, 1/(8*pi*t), (e2-3/(pi*t))/24 - e6/(24*e4)
ok = True
for tau in (1j, 1.7j, 0.3+0.8j):
    tot, an, rem = parts(tau); ok &= abs(tot-an-rem) < mp.mpf("1e-25")
check("H3 E8 charge energy per direction = c3/Im(tau) + modular remainder; remainder vanishes at tau=i", ok and abs(parts(1j)[2]) < mp.mpf("1e-25"))
# H4  duality sum rule  E(t) + E(1/t)/t^2 = 1/(4 pi t)  (Coulomb kernel at distance t)
Ed = lambda t: parts(1j*t)[0].real
ok = all(abs(Ed(t)+Ed(1/t)/t**2-1/(4*pi*t)) < mp.mpf("1e-22") for t in (mp.mpf("0.6"), mp.mpf("1"), mp.mpf("2.3")))
check("H4 sum rule E(t) + E(1/t)/t^2 = 1/(4 pi t)", ok)
# H5  sphere: g(r) = 1/(4 pi (1+r)), kernel between the north pole and the point at radius r on the opposite ray
g = lambda r: 1/(4*pi*(1+r))
check("H5 sphere kernel is Kelvin covariant with weight 1: g(1/r) = r g(r); g(1) = c3", all(abs(g(1/r)-r*g(r)) < mp.mpf("1e-25") for r in (mp.mpf("0.4"), mp.mpf(3))) and abs(g(1)-c3) < mp.mpf("1e-28"))
check("H6 weight fixes only the log derivative: -r g'(r) at r=1 is c3/2 = 1/(16 pi)", abs(-mp.diff(g, 1)-c3/2) < mp.mpf("1e-20"))
h = lambda r: 1/(4*pi*mp.sqrt(1+r*r))        # same weight, kernel at 90 degrees
check("H7 another weight 1 kernel (quarter marks) has a different fixed point value sqrt2 c3: the value is not fixed by weight", abs(h(mp.mpf("0.4")**-1)-mp.mpf("0.4")*h(mp.mpf("0.4"))) < mp.mpf("1e-25") and abs(h(1)/c3-mp.sqrt(2)) < mp.mpf("1e-25"))
# H8  pairing-wise dictionary  chord = 2 Im(tau): marks 1, i, -1, -i
marks = [1, 1j, -1, -1j]
chords = {"(13)(24) antipodal": abs(marks[0]-marks[2]), "(12)(34)": abs(marks[0]-marks[1]), "(14)(23)": abs(marks[0]-marks[3])}
ims = sorted([1.0, 1.0, 0.5])               # Im tau of i, 1+i, (1+i)/2  (lambda = 1/2, -1, 2)
lam = lambda tau: (mp.jtheta(2, 0, mp.e**(1j*pi*tau))/mp.jtheta(3, 0, mp.e**(1j*pi*tau)))**4
check("H8a lambda(i)=1/2, lambda(1+i)=-1, lambda((1+i)/2)=2", abs(lam(1j)-0.5) < mp.mpf("1e-20") and abs(lam(1+1j)+1) < mp.mpf("1e-20") and abs(lam((1+1j)/2)-2) < mp.mpf("1e-20"))
R["chords"] = {k: float(v) for k, v in chords.items()}; R["two_Im_tau"] = [2*x for x in ims]
match = sorted(float(v) for v in chords.values()) == sorted(2*x for x in ims)
check("H8b pairing-wise dictionary chord = 2 Im(tau) fails: chords (2, sqrt2, sqrt2) against (2, 2, 1)", not match, (sorted(R["chords"].values()), sorted(R["two_Im_tau"])))
json.dump(R, open("results_bridge.json", "w"), indent=1)
print("\nFAILS:", FAILS); sys.exit(1 if FAILS else 0)
