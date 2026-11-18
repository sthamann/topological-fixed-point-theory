#!/usr/bin/env python3
"""Session check of 17 November 2026: the double cover of the four marks
(torus at tau = i) and the E8 lattice weight on it.

A. Four marks as branch points: the cover of the sphere branched at
   1, i, -1, -i is the elliptic curve with j = 1728 (tau = i).
B. Its four 2-torsion points carry the affine group AGL(2,2) = S4; the
   automorphism i fixes two of them and swaps the other two. The Moebius
   symmetry group of the four marks themselves has order 8 (not 24).
C. Any faithful three-dimensional irreducible representation of S4 gives
   four lines with squared overlap 1/9 and rows (1/9, 4/9, 4/9).
D. Lattice weight at tau = i (q = exp(-2 pi)):
   mean lattice energy per direction = 1/(8 pi) = c3, for E8 exactly;
   split over the four glue sectors of D5 + A3; Cartan two-point function
   reduces to the Weierstrass function; mean of L0 - 1/3 vanishes.
Prints JSON. mpmath 40 digits, sympy exact where stated.
"""
import json
from itertools import permutations, product
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
out = {}

# A. j invariant of y^2 = x^4 - 1 (branch points = mu4) via cross ratio
lam = sp.Rational(1, 1)
z1, z2, z3, z4 = 1, sp.I, -1, -sp.I
cr = sp.simplify((z1-z3)*(z2-z4)/((z1-z4)*(z2-z3)))
j = sp.simplify(256*(cr**2-cr+1)**3/(cr**2*(cr-1)**2))
out["A"] = {"cross_ratio": str(cr), "j": str(j)}

# B. AGL(2,2) on the four 2-torsion points; action of multiplication by i
pts = [(0, 0), (1, 0), (0, 1), (1, 1)]                       # coordinates in (1/2, i/2)
mats = [m for m in product((0, 1), repeat=4) if (m[0]*m[3]-m[1]*m[2]) % 2 == 1]
perms = set()
for m in mats:
    for t in pts:
        perms.add(tuple(pts.index(((m[0]*x+m[1]*y+t[0]) % 2, (m[2]*x+m[3]*y+t[1]) % 2)) for x, y in pts))
mult_i = tuple(pts.index(((-y) % 2, x % 2)) for x, y in pts)   # i*(x + i y)/2 = (-y + i x)/2
def moebius_symmetries():
    marks = [1, sp.I, -1, -sp.I]
    count = 0
    for img in permutations(marks):
        # Moebius map fixed by three points; check the fourth
        a, b, c, d = sp.symbols("a b c d")
        f = lambda z, co: (co[0]*z+co[1])/(co[2]*z+co[3])
        cr0 = sp.simplify((marks[0]-marks[2])*(marks[1]-marks[3])/((marks[0]-marks[3])*(marks[1]-marks[2])))
        cr1 = sp.simplify((img[0]-img[2])*(img[1]-img[3])/((img[0]-img[3])*(img[1]-img[2])))
        count += cr0 == cr1
    return count
out["B"] = {"affine_group_order": len(perms), "is_all_of_S4": len(perms) == 24,
            "multiplication_by_i": {"permutation": mult_i,
                                    "fixed": [pts[k] for k in range(4) if mult_i[k] == k]},
            "moebius_symmetries_of_the_square_marks": moebius_symmetries()}

# C. S4 as rotations of the cube
rot = []
for perm in permutations(range(3)):
    for s in product((1, -1), repeat=3):
        m = sp.zeros(3)
        for r, c in enumerate(perm): m[r, c] = s[r]
        if m.det() == 1: rot.append(m)
order = lambda m: next(k for k in range(1, 9) if m**k == sp.eye(3))
c3el = next(m for m in rot if order(m) == 3); r4 = next(m for m in rot if order(m) == 4)
w = sp.Rational(-1, 2)+sp.sqrt(3)*sp.I/2
lines = [sp.simplify(sum((sp.conjugate(z)**k*c3el**k for k in range(3)), sp.zeros(3))/3) for z in (1, w, w**2)]
T = sp.Matrix(3, 3, lambda a, b: sp.nsimplify(sp.simplify(sp.re((lines[a]*r4*lines[b]*r4.T).trace()))))
out["C"] = {"rows": [[str(x) for x in T.row(k)] for k in range(3)]}

# D. lattice weight at tau = i
def thetas(t):                      # nome exp(-pi t), tau = i t
    n = mp.e**(-mp.pi*t)
    return mp.jtheta(2, 0, n), mp.jtheta(3, 0, n), mp.jtheta(4, 0, n)
def sectors(t):
    t2, t3, t4 = thetas(t)
    D = lambda n: ((t3**n+t4**n)/2, (t3**n-t4**n)/2, t2**n/2)          # root, vector, spinor
    a, b = D(5), D(3)
    return [a[0]*b[0], a[2]*b[2], a[1]*b[1], a[2]*b[2]]                 # k = 0, 1, 2, 3
theta = lambda t: sum(sectors(t))
energy = lambda f: -mp.diff(f, 1)/(2*mp.pi)                              # q d/dq at t = 1
Th = theta(1)
mean = energy(theta)/Th
sec_w = [s/Th for s in sectors(1)]
sec_e = [energy(lambda t, k=k: sectors(t)[k])/Th for k in range(4)]
eta = lambda t: mp.e**(-mp.pi*t/12)*mp.qp(mp.e**(-2*mp.pi*t))
chi = lambda t: theta(t)/eta(t)**8
E2 = 1-24*mp.nsum(lambda n: n*mp.e**(-2*mp.pi*n)/(1-mp.e**(-2*mp.pi*n)), [1, mp.inf])
out["D"] = {"Theta_E8(i)": mp.nstr(Th, 20), "chi_E8(i)": mp.nstr(chi(1), 20),
            "mean_lattice_energy": mp.nstr(mean, 25), "1/pi": mp.nstr(1/mp.pi, 25),
            "per_direction_minus_c3": mp.nstr(mean/8-1/(8*mp.pi), 5),
            "E2(i)": mp.nstr(E2, 25), "3/pi": mp.nstr(3/mp.pi, 25),
            "mean_L0_minus_c_over_24": mp.nstr(energy(chi)/chi(1), 5),
            "sector_weights_k0..3": [mp.nstr(x, 12) for x in sec_w],
            "sector_energy_shares_k0..3": [mp.nstr(x/mean, 12) for x in sec_e]}
# Cartan two-point function: zero-mode term cancels the E2 term
zero_mode = -4*mp.pi**2*(mean/4)             # (2 pi i)^2 <p_i p_i>, <p_i^2> = 2*mean/8
e2_term = mp.pi**2/3*E2
out["D"]["cartan_two_point_constant_terms"] = {"zero_mode": mp.nstr(zero_mode, 20), "E2_term": mp.nstr(e2_term, 20),
                                               "sum": mp.nstr(zero_mode+e2_term, 5)}
# same statement for other self-dual ranks: Leech-type rank 24 needs only weight 12 and E6(i) = 0
def e8_power(t, m): return theta(t)**m
for m in (2, 3):
    out["D"][f"per_direction_minus_c3_E8^{m}"] = mp.nstr(energy(lambda t: e8_power(t, m))/e8_power(1, m)/(8*m)-1/(8*mp.pi), 5)
# a lattice that is not self-dual: Z^8 (theta3^8) and D8
t2, t3, t4 = thetas(1)
out["D"]["per_direction_Z8"] = mp.nstr(energy(lambda t: thetas(t)[1]**8)/t3**8/8, 15)
out["D"]["per_direction_D8"] = mp.nstr(energy(lambda t: (thetas(t)[1]**8+thetas(t)[2]**8)/2)/((t3**8+t4**8)/2)/8, 15)
out["D"]["c3"] = mp.nstr(1/(8*mp.pi), 15)
print(json.dumps(out, indent=1))
