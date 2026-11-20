"""HOPFTWISTOR20261119: do the Hopf spheres / the twistor fibration CP1 -> CP3 -> S4
generate the recorded alpha constants, and how selective is that?"""
from common import *
from sympy import pi as PI, Rational as R, factorial, gamma, simplify, Integer
c = Checks("hopf_twistor")
c3 = 1 / (8 * PI)
volS = lambda n: 2 * PI**R(n + 1, 2) / gamma(R(n + 1, 2))          # unit S^n
volB = lambda n: PI**R(n, 2) / gamma(R(n, 2) + 1)
volCP = lambda n: PI**n / factorial(n)                              # Fubini-Study, S^(2n+1)/U(1)
w = lambda k: 1 / (2**(2 * k) * volS(2 * k - 1))
dtop = 48 * c3**4
c.ok("identity w_k = 1/(2^2k vol S^(2k-1)) = (k-1)! 2^(k-1) c3^k for all k (pure rewriting of c3 = 1/8pi)",
     all(simplify(w(k) - factorial(k - 1) * 2**(k - 1) * c3**k) == 0 for k in range(1, 10)))
c.ok("k=1 (S1, real Hopf):        w = c3", simplify(w(1) - c3) == 0)
c.ok("k=2 (S3, complex Hopf):     w = 2 c3^2, the recorded normal-form scale g = alpha/(2 c3^2)", simplify(w(2) - 2 * c3**2) == 0)
c.ok("k=4 (S7, quaternionic Hopf): w = 48 c3^4 = delta_top", simplify(w(4) - dtop) == 0)
c.ok("MISS k=3 (S5, no Hopf fibration): w = 8 c3^3, recorded cubic coefficient is 2 c3^3 (factor 4)", simplify(w(3) / (2 * c3**3)) == 4)
d2 = R(5, 4) * dtop**2
c.ok("MISS k=8 (S15, octonionic Hopf): w = 645120 c3^8, recorded second order (5/4) delta_top^2 = 2880 c3^8 (factor 224)",
     simplify(w(8) / c3**8) == 645120 and simplify(d2 / c3**8) == 2880 and simplify(w(8) / d2) == 224)
c.ok("MISS: the log coefficient (4/5) 41 c3^6 is no sphere weight (k=6 gives 3840 c3^6)", simplify(w(6) / c3**6) == 3840 and R(4, 5) * 41 != 3840)
# twistor fibration = Hopf chain
c.ok("twistor fibration: vol S7 = vol S3 * vol S4(r=1/2), vol CP3 = vol S7/2pi, fibre CP1 = S3/S1 = S2(r=1/2)",
     simplify(volS(7) - volS(3) * volS(4) / 16) == 0 and simplify(volCP(3) - volS(7) / (2 * PI)) == 0 and simplify(volS(3) / (2 * PI) - volS(2) / 4) == 0)
c.ok("so delta_top = c3 * (2 c3) * (24 c3^2): circle, twistor line, S4 base", simplify(dtop - c3 * (2 * c3) * (24 * c3**2)) == 0)
# look elsewhere: how many 'natural volumes' give delta_top with a power of two?
fam = {}
for n in range(1, 17):
    fam[f"S{n}"] = volS(n); fam[f"B{n}"] = volB(n)
    if n <= 8: fam[f"CP{n}"] = volCP(n)
hits = [(name, a) for name, v in fam.items() for a in range(0, 17) if simplify(1 / (2**a * v) - dtop) == 0]
print("  volume hits for delta_top:", hits)
c.ok("look-elsewhere: delta_top is not unique to S7 (B8 and CP4 with 2^11 give the same number)", len(hits) >= 3 and ("S7", 8) in hits, hits)
# instanton on the quaternionic Hopf bundle
al = mpf(1) / alpha_inv(lambda x: PHI_TREE + (DTOP * exp(-2 * x)) * (1 - DTOP * exp(-2 * x))**(-mpf(5) / 4))
c.ok("BPST instanton (= quaternionic Hopf bundle over S4) has weight exp(-8 pi^2/e^2) = exp(-2 pi/alpha) ~ 1e-374: numerically irrelevant, no route to the recorded factors",
     exp(-2 * pi / al) < mpf(10)**(-370))
c.ok("selectivity: the Hopf chain exists for k = 1,2,4,8 only (Adams); it singles out neither g=5 nor N_fam=3 nor 5/4", True, "5 and 3 enter only through D5+A3, not through any sphere in the chain")
c.done({"verdict": "three exact volume coincidences (k=1,2,4), three misses (k=3,6,8), non-unique; naming level, no derivation",
        "hits": [list(h) for h in hits]})
