#!/usr/bin/env python3
"""RADICALCHARGED20261118: geladene Wurzelstroeme von (E8)_1 auf dem Torus.
Zwei unabhaengige Wege: (a) Gitterspur im D8 Chart mit Primformen,
(b) elliptische Funktionen aus Ward Identitaeten. Reproduktion: python3 -B check_charged.py"""
import json, sys
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi; R = {}; FAILS = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + str(info) if info != "" else ""))
    if not ok: FAILS.append(name)

class Torus:
    def __init__(s, tau):
        s.tau = tau; s.nome = mp.e**(1j*pi*tau); s.q = s.nome**2
        q = s.q
        s.E2 = 1-24*mp.nsum(lambda n: n*q**n/(1-q**n), [1, mp.inf])
        s.E4 = 1+240*mp.nsum(lambda n: n**3*q**n/(1-q**n), [1, mp.inf])
        s.E6 = 1-504*mp.nsum(lambda n: n**5*q**n/(1-q**n), [1, mp.inf])
        s.kappa = s.E6/(12*s.E4)
        s.t1p = mp.jtheta(1, 0, s.nome, 1)
    # (b) elliptic functions, periods 2 pi and 2 pi tau
    def zeta(s, w):
        q = s.q
        return s.E2/12*w + mp.cot(w/2)/2 + 2*mp.nsum(lambda n: q**n/(1-q**n)*mp.sin(n*w), [1, mp.inf])
    def wp(s, w):
        q = s.q
        return -s.E2/12 + 1/(4*mp.sin(w/2)**2) - 2*mp.nsum(lambda n: n*q**n/(1-q**n)*mp.cos(n*w), [1, mp.inf])
    def wpp(s, w):
        q = s.q
        return -mp.cos(w/2)/(4*mp.sin(w/2)**3) + 2*mp.nsum(lambda n: n*n*q**n/(1-q**n)*mp.sin(n*w), [1, mp.inf])
    def two_ward(s, w): return -s.wp(w) - s.kappa
    def three_ward(s, w1, w2, w3):
        a, b, c = w1-w2, w2-w3, w3-w1
        G = s.zeta(a)+s.zeta(b)+s.zeta(c)
        K = G**3 + (s.wpp(a)+s.wpp(b)+s.wpp(c))/2
        return -K/3 - s.kappa*G
    # (a) lattice trace in the D8 chart
    def prime(s, w): return 2*mp.jtheta(1, w/2, s.nome)/s.t1p
    def theta_e8(s, v):   # sum over E8 of q^{p^2/2} exp(i p.v); v has at least one zero entry
        P = lambda k: mp.fprod(mp.jtheta(k, x/2, s.nome) for x in v)
        return (P(3)+P(4)+P(2))/2
    def two_lat(s, w):
        v = [w, -w]+[0]*6
        return -s.theta_e8(v)/s.theta_e8([0]*8)/s.prime(w)**2
    def three_lat(s, w1, w2, w3):
        v = [w1-w3, w2-w1, w3-w2]+[0]*5
        return s.theta_e8(v)/s.theta_e8([0]*8)/(s.prime(w1-w2)*s.prime(w2-w3)*s.prime(w1-w3))

pts = [mp.mpf("0.7")-0.3j, mp.mpf("2.1")-0.9j, mp.mpf("1.3")+0.4j]
ok2 = True; ok3 = True; ratios = []
for tau in (1j, 1.3j, 0.2+0.9j):
    T = Torus(tau)
    for w in pts: ok2 &= abs(T.two_lat(w)-T.two_ward(w)) < mp.mpf("1e-20")
    for (a, b, c) in [(pts[0], pts[1], 0), (pts[2], pts[0], mp.mpf("0.2")), (pts[1], 0, pts[2])]:
        r = T.three_lat(a, b, c)/T.three_ward(a, b, c); ratios.append(r)
ok3 = all(abs(r-ratios[0]) < mp.mpf("1e-18") for r in ratios) and abs(abs(ratios[0])-1) < mp.mpf("1e-18")
check("G1 charged root two-point function = -wp - E6/(12 E4) (lattice trace vs Ward form, three moduli)", ok2)
check("G2 charged three-point function = c * [-(G^3 + sum wp'/2)/3 - kappa G], constant phase c, |c|=1", ok3, mp.nstr(ratios[0], 6))

Ti = Torus(1j)
check("G3 kappa = 0 at tau=i", abs(Ti.kappa) < mp.mpf("1e-25"))
varpi = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*pi))
two = Ti.two_lat(pi); three = Ti.three_ward(0, pi/2, pi)
R["two_point_antipodal"] = {"plane": -0.25, "torus": float(two.real)}
R["three_point_marks_1_i_-1"] = {"plane_abs": 0.25, "torus_abs": float(abs(three)), "torus_value": [float(three.real), float(three.imag)]}
print("   Zweipunkt am Gegenpunkt:", mp.nstr(two, 15), "  Dreipunkt an 1,i,-1:", mp.nstr(three, 15))
check("G4 charged antipodal two-point value at tau=i is -(varpi/2pi)^2", abs(two + (varpi/(2*pi))**2) < mp.mpf("1e-22"))
# plane limit of the three-point kernel
Tp = Torus(40j)
check("G5 plane limit of the three-point kernel at marks 1, i, -1 has modulus 1/4", abs(abs(Tp.three_ward(0, pi/2, pi))-mp.mpf(1)/4) < mp.mpf("1e-20"))
e1 = (varpi/(2*pi))**2
closed = (1+mp.sqrt(2))*e1**mp.mpf(1.5)
R["three_point_closed_form"] = "(1+sqrt2)*(varpi/(2 pi))^3"
check("G6 three-point modulus at tau=i is (1+sqrt2)*(varpi/2pi)^3 = 0.175454, not 1/4", abs(closed-abs(three)) < mp.mpf("1e-22"), mp.nstr(abs(three), 12))
# weight 1/2 kernel: sphere contract gives antipodal 1/(8 pi) = (1/4pi)*(plane Szego value 1/2)
c3 = 1/(8*pi)
sz = {}
for nu in (2, 3, 4):
    sz[nu] = (mp.jtheta(nu, pi/2, Ti.nome)*Ti.t1p/(mp.jtheta(nu, 0, Ti.nome)*mp.jtheta(1, pi/2, Ti.nome)*2)).real
    print(f"   Szegoe Kern am Gegenpunkt, Spinstruktur {nu}: {mp.nstr(sz[nu], 12)}  -> /(4 pi) = {mp.nstr(sz[nu]/(4*pi)/c3, 10)} c3")
R["weight_half_antipodal_over_c3"] = {str(k): float(v/(4*pi)/c3) for k, v in sz.items()}
check("G7 no spin structure at tau=i reproduces the sphere antipodal value c3 (values 0, 0.8346, 1.1803 c3)",
      all(abs(v/(4*pi)/c3-1) > 0.1 for v in sz.values()))
check("G8 the two nonzero values are theta4(i)^2/2 and theta3(i)^2/2",
      abs(sz[3]-mp.jtheta(4, 0, Ti.nome)**2/2) < mp.mpf("1e-22") and abs(sz[4]-mp.jtheta(3, 0, Ti.nome)**2/2) < mp.mpf("1e-22"))
json.dump(R, open("results_charged.json", "w"), indent=1)
print("\nFAILS:", FAILS); sys.exit(1 if FAILS else 0)
