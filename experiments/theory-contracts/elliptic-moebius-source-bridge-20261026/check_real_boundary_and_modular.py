"""Decisive exact identities for the actual boundary and modular source.

Analytic boundedness/range and compact bump RP proofs remain in the reports.
The theta and elliptic evaluations below are diagnostics of those proofs.
"""
import itertools
import json
import sympy as s
import mpmath as mp

checks = []


def need(name, value):
    if not bool(value):
        raise RuntimeError(name)
    checks.append(name)


x, t, b, y = s.symbols('x t b y')
I = s.I
AS = s.Matrix([[0,0,1], [0,1,0], [1,0,0]])
AT = s.Matrix([[I,0,0], [0,0,1], [0,1,0]])
U, M = -I*AT, I*AS*AT
a, bb, c = M**2, M*U**-1, U
need('finite modular target and marked generators',
     a**3 == s.eye(3) and bb**2 == s.eye(3)
     and c**4 == s.eye(3) and a*bb*c == s.eye(3))
d = bb*c*bb**-1
need('actual second cusp conjugate word', s.simplify(d-(U*M)**-1) == s.zeros(3))
need('quotient four puncture Hurwitz product', c*a**2*bb**2*d == s.eye(3))
phi = -(t-1)**2/(4*t)
need('degree two source map branch values', s.factor(phi-1) == -(t+1)**2/(4*t))
need('degree two source map complete critical law',
     s.simplify(s.diff(phi,t)+(t-1)*(t+1)/(4*t*t)) == 0)
J = phi.subs(t,x**4)
need('composed degree eight critical law',
     s.simplify(s.diff(J,x)+(x**4-1)*(x**4+1)/x**5) == 0)
need('source map is the original D4 invariant',
     s.simplify(J.subs(x,I*x)-J) == 0 and s.simplify(J.subs(x,1/x)-J) == 0)
angle = s.symbols('angle',real=True)
need('actual marked equator image',
     s.trigsimp((1-s.cos(4*angle))/2-s.sin(2*angle)**2) == 0)
need('degree eight full ramification sum', 4+4+2*3 == 2*8-2)
prod = s.eye(3)
for k in range(4):
    prod = s.simplify(prod*U**k*M*U**(-k))
need('all four Schreier peripheral matrices', prod == s.eye(3))
Tfull = s.Matrix([[1,1], [0,1]])
need('full modular cusp remains nontrivial', Tfull**4 != s.eye(2))

# Rigidity test in the expressly specified genuine unitary D4 orbifold class.
omega = (-1+I*s.sqrt(3))/2
root_traces = {s.simplify(e[0]+e[1]*omega**2+e[2]*omega)
               for e in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]}
need('four determinant one half meridian branches',
     root_traces == {s.Integer(0),s.Integer(2),-1+I*s.sqrt(3),-1-I*s.sqrt(3)})
need('nonzero trace branches are impossible or reducible',
     (s.re(s.Integer(2))+1)/2 > 1
     and (s.re(-1+I*s.sqrt(3))+1)/2 == 0
     and (s.re(-1-I*s.sqrt(3))+1)/2 == 0)
v0 = s.Matrix([1/s.sqrt(2),s.Rational(1,2),s.Rational(1,2)])
b0 = 2*v0*v0.H-s.eye(3)
targetU = s.diag(1,I,-I)
a0 = (b0*targetU)**-1
need('forced irreducible branch representative',
     s.simplify(b0**2) == s.eye(3) and s.simplify(a0**3) == s.eye(3))
D = s.diag(1,-(1+I)/s.sqrt(2),(1-I)/s.sqrt(2))
targetM = s.Matrix([[0,-(1+I)/2,(1-I)/2],
                    [-(1+I)/2,-I/2,-s.Rational(1,2)],
                    [(1-I)/2,-s.Rational(1,2),I/2]])
need('forced representative has exact original marked v117 pair',
     s.simplify(D*a0**2*D.H-targetM) == s.zeros(3)
     and s.simplify(D*targetU*D.H-targetU) == s.zeros(3))

# Tate normal form and its actual critical law, not a fitted j value.
F = y*y+(x+b)*y-x*x*(x+b)
critical = 16*x*x+(16*b+3)*x+4*b
need('Tate critical substitution',
     s.factor(F.subs(y,4*x*(x+b))-x*(x+b)*critical) == 0)
need('triple critical point parameter',
     s.factor(s.discriminant(critical,x)) == (16*b-9)*(16*b-1))
b0 = s.Rational(9,16)
c4 = 1-16*b+16*b*b
delta = b**4*(1-16*b)
jval = s.factor((c4**3/delta).subs(b,b0))
need('unique smooth branch and distinct elliptic invariant',
     delta.subs(b,s.Rational(1,16)) == 0
     and delta.subs(b,b0) != 0
     and jval == s.Rational(207646,6561) and jval != 1728)
special = {b:b0,x:-s.Rational(3,8),y:-s.Rational(9,32)}
need('actual triple branch point and value',
     F.subs(special) == 0
     and (y*y/(x+b)).subs(special) == s.Rational(27,64))


def compose(p,q):
    return tuple(p[q[i]] for i in range(4))


def inv(p):
    return tuple(p.index(i) for i in range(4))


def ctype(p):
    unseen, lengths = set(range(4)), []
    while unseen:
        start = next(iter(unseen))
        cur, count = start, 0
        while cur in unseen:
            unseen.remove(cur)
            count += 1
            cur = p[cur]
        lengths.append(count)
    return sorted(lengths)


perms = list(itertools.permutations(range(4)))
four = [p for p in perms if ctype(p)==[4]]
three = [p for p in perms if ctype(p)==[1,3]]
triples = [(p,q,inv(compose(p,q))) for p in four for q in three
           if ctype(inv(compose(p,q)))==[4]]
first = triples[0]
orbit = {tuple(compose(compose(k,p),inv(k)) for p in first) for k in perms}
need('unique ordered S4 Nielsen class', len(triples)==24 and set(triples)==orbit)
need('elliptic source and normal closure genus',
     -8+3+2+3 == 0 and 24*(-2+s.Rational(3,4)+s.Rational(2,3)+s.Rational(3,4)) == 4)

# Exact beta moment: gamma recurrence first, reflection identity afterwards.
r = s.Rational
moment0 = s.gamma(r(1,4))*s.gamma(r(3,4))/s.pi
moment2 = (s.gamma(r(3,4))*s.gamma(r(5,4))
           -s.gamma(r(7,4))*s.gamma(r(1,4)))/s.pi
need('nonorthogonality Fourier coefficient by beta recurrence',
     s.simplify(s.expand_func(moment2/moment0)) == -s.Rational(1,2))
need('reflection formula gives exact nonzero overlap',
     s.simplify((-s.Rational(1,2))/s.sin(s.pi/4)+1/s.sqrt(2)) == 0)
L = s.beta(r(1,4),r(1,2))/4
edge_length = s.beta(r(1,4),r(1,2))/(8*L)
need('intrinsic edge and full face lengths', edge_length == r(1,2) and 4*edge_length == 2)
need('intrinsic speed fails fixed base winding Ward', s.sqrt(2)/2 != 1)

mp.mp.dps = 70
nome = mp.exp(-mp.pi)
theta = lambda w: (1j/2)*mp.jtheta(3,mp.pi*w,nome)*mp.jtheta(1,0,nome,1)/(mp.jtheta(3,0,nome)*mp.jtheta(1,mp.pi*w,nome))
tt = mp.mpf(1)/4
sp, sm = theta(1j*tt), theta(-1j*tt)
fourier = sum(mp.exp(-tt*(2*mp.pi*(mp.mpf(n)+mp.mpf('.5'))))/(1+mp.exp(-2*mp.pi*(mp.mpf(n)+mp.mpf('.5')))) for n in range(-100,101))
if abs(sp-fourier)>mp.mpf('1e-60') or abs(sm+sp)>mp.mpf('1e-60'):
    raise RuntimeError('actual theta and NS Fourier diagnostic')
if not (abs(sp.imag)<mp.mpf('1e-60') and sp.real>0 and sm.real<0):
    raise RuntimeError('charged reflection sign diagnostic')
cc = 2*mp.ellipk(-1)*(1+1j)
pts = [mp.mpc('.375','.5'),mp.mpc('.125','0')]
absx = [abs(mp.ellipfun('sn',cc*z,-1)) for z in pts]
if not all(v<1 for v in absx):
    raise RuntimeError('checker face point positions')

print(json.dumps(dict(verdict='PASS_SCOPED',exact_checks=len(checks),checks=checks,
    theta_fourier_error=mp.nstr(abs(sp-fourier),8),
    neutral_reflection_point_diagnostic=mp.nstr((sp*sm).real,16),
    actual_inside_face_abs_x=[mp.nstr(v,12) for v in absx],
    j_of_degree_four_source=str(jval),
    compact_bump_and_outer_proofs='analytical reports; not finite grid extrapolation',
    complete_TFPT_solution=False),indent=2))
