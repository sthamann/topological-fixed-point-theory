"""Exact ramified spin kernel and literal archived cusp formula checks.

Analytic uniqueness of the Szego kernel is proved in GLOBAL-SPIN-KERNEL.txt.
The three theta comparisons below are diagnostics, not that proof.
"""
import json
import sympy as s
import mpmath as mp

checks = []


def check(name, value):
    if not bool(value):
        raise RuntimeError(name)
    checks.append(name)


I = s.I
eye = s.eye(2)
J = s.Matrix([[0, I], [-I, 0]])
plus = (eye + J) / 2
minus = (eye - J) / 2
check('constant orthogonal spectral projectors',
      J.H == J and J**2 == eye and plus * minus == s.zeros(2))
a = s.symbols('a', nonzero=True)
full = s.Matrix([[a + 1/a, I*(a - 1/a)],
                 [-I*(a - 1/a), a + 1/a]]) / 2
check('full two sheet algebraic kernel diagonalization',
      s.simplify(full - a * plus - minus/a) == s.zeros(2))
q, r = s.symbols('q r', positive=True)
local = s.Matrix([[1/(q-r), I/(q+r)],
                  [-I/(q+r), 1/(q-r)]])/(2*s.sqrt(q*r))
diag = (s.sqrt(q/r)*plus + s.sqrt(r/q)*minus)/(q*q-r*r)
check('ramified local Cauchy identity', s.simplify(local-diag) == s.zeros(2))
x = s.symbols('x')
ratio = (x*x-1)/(x*x+1)
eta = x/(x**4-1)
check('global logarithmic derivative', s.simplify(s.diff(ratio, x)/(4*ratio)-eta) == 0)
check('quarter clock interchanges the two character lines',
      s.simplify(ratio.subs(x, I*x)*ratio-1) == 0
      and s.simplify(I*eta.subs(x, I*x)+eta) == 0)
res = [s.residue(eta, x, z) for z in (1, I, -1, -I)]
check('four actual residues', res == [s.Rational(1,4), -s.Rational(1,4),
                                       s.Rational(1,4), -s.Rational(1,4)])
spin_meridians = [s.diag(I, -I), s.diag(-I, I)] * 2
prod = eye
for m in spin_meridians:
    prod *= m
check('global spin meridian product and central squares',
      prod == eye and all(m*m == -eye for m in spin_meridians))
P = s.Matrix([[0, 1], [1, 0]])
check('two radial resolutions have flat total holonomy one',
      P**4 == eye and (I*P)**2 == -eye)

# Literal matrices from archived TFPT45 section004, lines18-44.
w = (-1 + I*s.sqrt(3))/2
D = s.diag(1, w, w*w)
P3 = s.Matrix([[0,1,0], [0,0,1], [1,0,0]])
bad = s.simplify(D*(P3*D*P3.inv())*(P3**2*D*P3**-2))
check('archived printed fourth cusp is identity', bad == s.eye(3))

# Unchanged verification/v117_monodromy_weyl_a3.py:52-56.
M = s.Matrix([[0, -(1+I)/2, (1-I)/2],
              [-(1+I)/2, -I/2, -s.Rational(1,2)],
              [(1-I)/2, -s.Rational(1,2), I/2]])
U = s.diag(1, I, -I)
prod = s.eye(3)
for k in range(4):
    mk = U**k*M*U**(-k)
    check('v117 cusp spectrum '+str(k),
          s.expand(mk.charpoly(x).as_expr()) == x**3-1)
    prod = s.simplify(prod*mk)
check('existing v117 repairs the four cusp product',
      prod == s.eye(3) and U**4 == s.eye(3) and (M*U)**4 == s.eye(3))

# Actual theta-square modular matrices, with the scalar factors retained.
AS = s.Matrix([[0,0,1], [0,1,0], [1,0,0]])
AT = s.Matrix([[I,0,0], [0,0,1], [0,1,0]])
UT, MT = -I*AT, I*AS*AT
X = s.Matrix([[-1-I,0,0], [0,I,-I], [0,1,1]])/s.sqrt(2)
check('theta raw modular central relation', AS**2 == s.eye(3) and (AS*AT)**3 == I*s.eye(3))
check('explicit theta v117 unitary dictionary',
      s.simplify(X.H*X) == s.eye(3)
      and s.simplify(X*UT*X.H-U) == s.zeros(3)
      and s.simplify(X*MT*X.H-M) == s.zeros(3))


def group(generators):
    ident = s.eye(3)
    seen, todo = {tuple(ident)}, [ident]
    while todo:
        cur = todo.pop()
        for gen in generators:
            nxt = s.expand(cur*gen)
            key = tuple(nxt)
            if key not in seen:
                seen.add(key)
                todo.append(nxt)
                if len(seen) > 200:
                    raise RuntimeError('unexpected theta group growth')
    return seen


raw, normalized = group([AS,AT]), group([UT,MT])
check('theta raw group has 96 elements and four scalar elements',
      len(raw) == 96 and sum(tuple(a*s.eye(3)) in raw for a in [1,I,-1,-I]) == 4)
check('normalized theta group has 24 elements', len(normalized) == 24)
check('same geometric Q cannot be the v117 quarter mark',
      (-AS)**2 == s.eye(3) and (-AS).eigenvals() == {1:1, -1:2}
      and U.eigenvals() == {1:1, I:1, -I:1})

# Independent numerical diagnostics on the same square torus, in a single
# analytic chart. Decimal strings avoid introducing machine precision.
mp.mp.dps = 70
c = 2*mp.ellipk(-1)*(1+1j)
nome = mp.exp(-mp.pi)
pairs = [('0.02 0.01','0.09 0.03'), ('0.07 0.02','0.03 0.06'),
         ('0.09 -0.04','0.04 -0.02')]
errors = []
for za, wa in pairs:
    z, v = mp.mpc(*za.split()), mp.mpc(*wa.split())

    def chart(t):
        u = c*t
        xx = mp.ellipfun('sn', u, -1)
        yy = 1j*mp.ellipfun('cn', u, -1)*mp.ellipfun('dn', u, -1)
        aa = 1j*mp.sqrt(1-xx*xx)
        return xx, aa, yy/aa

    xx, aa, bb = chart(z)
    vv, cc, dd = chart(v)
    alg = -1j*c*(aa*dd+bb*cc)/(2*(xx-vv))
    theta = mp.pi*mp.jtheta(3, mp.pi*(z-v), nome)*mp.jtheta(1,0,nome,1)/(mp.jtheta(3,0,nome)*mp.jtheta(1,mp.pi*(z-v),nome))
    err = abs(alg-theta)
    if err >= mp.mpf('1e-60'):
        raise RuntimeError('same-E theta diagnostic')
    errors.append(mp.nstr(err, 8))

print(json.dumps(dict(verdict='PASS_SCOPED', exact_checks=len(checks),
                      checks=checks, theta_diagnostic_errors=errors,
                      analytic_proof='GLOBAL-SPIN-KERNEL.txt',
                      positive_equal_time_CAR_claimed=False,
                      complete_TFPT_solution=False), indent=2))
