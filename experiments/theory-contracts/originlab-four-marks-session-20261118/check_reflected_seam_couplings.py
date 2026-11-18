#!/usr/bin/env python3
"""Twelfth pasted note: two sources that agree on the whole mark circle (18 November 2026).

Left times right, with the characters of sixteen fermions
    chi_0 = (theta3^8 + theta4^8)/(2 eta^8),  chi_v = (theta3^8 - theta4^8)/(2 eta^8),  chi_s = chi_c = theta2^8/(2 eta^8):
    Z_E8 = |chi_0 + chi_s|^2                     (E8 on each side)
    Z_D8 = |chi_0|^2 + |chi_v|^2 + 2 |chi_s|^2   (diagonal coupling)

Sections
    C1  the difference: 2 |theta3/eta|^16 (Im lambda)^2
    C2  equal on the circle; 144 at the square; 7/4 for both
    C3  the tetrahedral point of the recorded family: 0 against 60.476210395...
    C4  states on the first level: 496 each, 0 against 256 of spin zero
"""
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 30
PI = mp.pi
FAIL = []


def check(name, ok, info=None):
    print(("  ok   " if ok else "  FAIL ") + name + (("   " + str(info)) if info is not None else ""))
    if not ok:
        FAIL.append(name)


def data(tau):
    N = 80
    q = lambda e: mp.exp(1j*PI*tau*e)
    t2 = mp.nsum(lambda n: q((n+mp.mpf(1)/2)**2), [-N, N])
    t3 = mp.nsum(lambda n: q(n*n), [-N, N])
    t4 = mp.nsum(lambda n: (-1)**int(n)*q(n*n), [-N, N])
    eta = mp.exp(1j*PI*tau/12)*mp.nprod(lambda n: 1-mp.exp(2j*PI*tau*n), [1, mp.inf])
    c0, cv, cs = (t3**8+t4**8)/(2*eta**8), (t3**8-t4**8)/(2*eta**8), t2**8/(2*eta**8)
    return abs(c0+cs)**2, abs(c0)**2+abs(cv)**2+2*abs(cs)**2, t2**4/t3**4, abs(t3/eta)**16


# ---------------------------------------------------------------- C1
print("C1  the difference of the two torus traces")
worst = 0
for tau in (mp.mpc("0.23", "1.07"), mp.mpc("-0.41", "0.77"), mp.mpc("0.5", "0.9")):
    ze, zd, lam, pref = data(tau)
    worst = max(worst, abs(zd-ze-2*pref*lam.imag**2)/zd)
check("Z_D8 - Z_E8 = 2 |theta3/eta|^16 (Im lambda)^2, lambda = theta2^4/theta3^4, at three moduli",
      worst < mp.mpf(10)**(-22), mp.nstr(worst, 3))

# ---------------------------------------------------------------- C2
print("C2  on the mark circle")
tau_of = lambda d: 1j*mp.ellipk(1-mp.sin(d/2)**2)/mp.ellipk(mp.sin(d/2)**2)
dev = 0
for d in (mp.mpf("0.4"), mp.atan(mp.mpf(4)/3), mp.mpf("1.2")):
    ze, zd, lam, _ = data(tau_of(d))
    dev = max(dev, abs(zd-ze)/ze, abs(lam.imag))
check("marks on the circle: lambda is real and the two traces coincide", dev < mp.mpf(10)**(-22))
ze, zd, _, _ = data(mp.mpc(0, 1))
check("at the square both are 144 = 12^2", abs(ze-144) < mp.mpf(10)**(-20) and abs(zd-144) < mp.mpf(10)**(-20))
dcw = mp.atan(mp.mpf(4)/3)


def response(d, which):
    ze_, zd_, _, _ = data(tau_of(d))
    zz = ze_ if which == "E8" else zd_
    return mp.sqrt(zz)/(16*mp.sin(d)**2)**(mp.mpf(1)/3)   # normalisation of the nineteenth round


r = {k: response(dcw, k)/response(PI/2, k) for k in ("E8", "D8")}
check("moved marks against the square, in the normalisation of the nineteenth round: 7/4 for both",
      all(abs(v-mp.mpf(7)/4) < mp.mpf(10)**(-20) for v in r.values()), {k: mp.nstr(v, 10) for k, v in r.items()})

# ---------------------------------------------------------------- C3
print("C3  the tetrahedral point")
rho = mp.exp(2j*PI/3)
ze, zd, lam, pref = data(rho)
check("there lambda = (1 + i sqrt 3)/2 up to conjugation", min(abs(lam-mp.exp(1j*PI/3)), abs(lam-mp.exp(-1j*PI/3))) < mp.mpf(10)**(-22))
check("E8 on each side: the trace vanishes", ze < mp.mpf(10)**(-40), mp.nstr(ze, 3))
check("diagonal coupling: 60.476210395...", abs(zd-mp.mpf("60.476210395")) < mp.mpf(10)**(-8), mp.nstr(zd, 14))
check("it equals (3/2) |theta3/eta|^16", abs(zd-mp.mpf(3)/2*pref) < mp.mpf(10)**(-22))
print("      in the recorded family this is the parameter value z^4 = 1, where j = 1728 (1 - z^4) = 0.")

# ---------------------------------------------------------------- C4
print("C4  the first level")
K = 9                                             # series in t = q^(1/2) up to t^8


def pmul(a, b):
    out = np.zeros(K, dtype=np.int64)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y and i+j < K:
                    out[i+j] += x*y
    return out


def ppow(a, n):
    out = np.zeros(K, dtype=np.int64)
    out[0] = 1
    for _ in range(n):
        out = pmul(out, a)
    return out


th3 = np.zeros(K, dtype=np.int64)
th4 = np.zeros(K, dtype=np.int64)
for n in range(-3, 4):
    if n*n < K:
        th3[n*n] += 1
        th4[n*n] += (-1)**n
th2r = np.zeros(K, dtype=np.int64)                # theta2 = 2 t^(1/4) (1 + t^2 + t^6 + ...)
for n in range(0, 3):
    if n*(n+1) < K:
        th2r[n*(n+1)] += 1
inv = np.zeros(K, dtype=np.int64)
inv[0] = 1
for n in range(1, K):
    f = np.zeros(K, dtype=np.int64)               # 1/(1 - t^(2n)) = sum t^(2nk)
    for k in range(0, K):
        if 2*n*k < K:
            f[2*n*k] = 1
    for _ in range(8):
        inv = pmul(inv, f)
c0 = pmul((ppow(th3, 8)+ppow(th4, 8))//2, inv)
cv = pmul((ppow(th3, 8)-ppow(th4, 8))//2, inv)
cs = pmul(ppow(th2r, 8)*128, inv)                 # times t^2 overall
check("characters: vacuum 1 + 120 q, vector 16 q^(1/2), spinor 128 q (in units of q^(-1/3))",
      c0[0] == 1 and c0[2] == 120 and cv[1] == 16 and cs[0] == 128, (int(c0[2]), int(cv[1]), int(cs[0])))
e8_level1 = 2*(int(c0[2])+int(cs[0]))            # (1,0) and (0,1)
d8_level1 = 2*int(c0[2])+int(cv[1])**2            # currents left and right, plus (1/2, 1/2)
check("first level: 248 + 248 = 496 for E8 on each side, 120 + 120 + 256 = 496 for the diagonal coupling",
      e8_level1 == 496 and d8_level1 == 496)
check("of these, spin zero: none for E8, 256 = 16 * 16 for the diagonal coupling", int(cv[1])**2 == 256)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("all checks passed")
