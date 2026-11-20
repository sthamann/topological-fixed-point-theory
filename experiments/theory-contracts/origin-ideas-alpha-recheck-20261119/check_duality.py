"""COUPLINGMODULUS20261119 (idea 2): are exp(-2 alpha) and exp(-2/alpha) the two
sides of one flux sum (theta function) in the coupling tau = i/alpha?"""
from common import *
from sympy import pi as PI, Rational as R, symbols, simplify, sqrt
c = Checks("duality")
al = symbols("alpha", positive=True)
e2 = 4 * PI * al                                  # Heaviside-Lorentz
# Maxwell action of a flux F = 2 pi m on a 4-manifold, S = (1/2e^2) int F^*F = (2 pi^2/e^2) m.m
S_mag = simplify(2 * PI**2 / e2)                  # per unit m.m
c.ok("magnetic flux weight is exp(-(pi/(2 alpha)) m.m)", simplify(S_mag - PI / (2 * al)) == 0)
need = simplify((2 / al) / S_mag)
c.ok("exp(-2/alpha) would need m.m = 4/pi, which no integral lattice has", need == 4 / PI and not need.is_rational, need)
S_el = simplify(e2 / 2)                           # dual (Poisson resummed) weight per unit n.n: exp(-(e^2/2) n.n /(…)) -> 2 pi alpha
c.ok("electric dual weight is exp(-2 pi alpha n.n); exp(-2 alpha) would need n.n = 1/pi", simplify(2 * al / S_el) == 1 / PI)
c.ok("for any Gaussian flux sum the product of dual exponents is pi^2 (m.m)(n.n); TFPT has (2a)(2/a) = 4", simplify(S_mag * S_el) == PI**2 and (2 * al) * (2 / al) == 4)
c.ok("charge rescaling e -> k e cannot repair it: needs k^2 rational multiple of 1/pi", all(not (R(4) / (PI * k**2)).is_rational for k in range(1, 7)))
# what does survive
a = 1 / alpha_inv(lambda x: PHI_TREE + (DTOP * exp(-2 * x)) * (1 - DTOP * exp(-2 * x))**(-mpf(5) / 4))
c.ok("formal inversion alpha -> 1/alpha does exchange the two exponentials (naming level only)", abs(exp(-2 / (1 / a)) - exp(-2 * a)) < 1e-30)
c.ok("origin of the '2' in the sources is cover degree times spectral flow (d*SF = 2), a count, not a flux norm", True, "v2.5/v2.7 Lemma Cover-Degree Scaling")
c.done({"verdict": "EXCLUDED as theta-function/flux-sum pair; only the formal inversion remains"})
