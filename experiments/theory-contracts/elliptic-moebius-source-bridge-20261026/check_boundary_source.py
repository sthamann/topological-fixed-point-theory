"""Exact identities for the existing elliptic cut, not a P1 field simulation."""
import json
from fractions import Fraction as F
from pathlib import Path

import sympy as s

checks = []


def need(name, condition, scope):
    if not condition:
        raise RuntimeError(name)
    checks.append(dict(name=name, status="PASS", scope=scope))


k = s.Symbol("k", integer=True)
need("quarter-fold-domains", s.simplify((2*k+s.Rational(1, 2))-2*(k+s.Rational(1, 4))) == 0
     and s.simplify((2*k+s.Rational(3, 2))-2*(k+s.Rational(3, 4))) == 0,
     "Analytic mode identity for every integer k; Sobolev domains in proof")
need("half-period-characters", s.simplify(s.exp(s.I*s.pi*(2*k+s.Rational(1, 2)))-s.I) == 0
     and s.simplify(s.exp(s.I*s.pi*(2*k+s.Rational(3, 2)))+s.I) == 0,
     "Complex-linear boundary translation, not Q rotation or antiunitary OS")


def fold(n):
    block = n % 2
    return block, (n-block)//2


need("winding-with-carry", all(
    fold(n+1) == ((1, fold(n)[1]) if fold(n)[0] == 0 else (0, fold(n)[1]+1))
    for n in range(-50, 51)),
    "Finite replay of the analytic even/odd-index formula")
need("fixed-dimensionless-modulus", 2*F(1, 2)*2 == 2,
     "2 W B_Sigma = 2 pi B_cov, at W=1/2 and circumference1")
need("distinct-full-collar-moduli", F(2, 2) != F(1, 2),
     "Old full flat Pin test has modulus1, elliptic annulus has modulus1/2")

x = s.Symbol("x", positive=True)
graph = s.Matrix([[x*x, x], [x, 1]])/(1+x*x)
need("two-boundary-calderon-projector", s.simplify(graph*graph-graph) == s.zeros(2)
     and graph.T == graph and s.simplify(s.trace(graph)) == 1,
     "Exact rational identity with x=exp(W lambda); chiral product Dirac")
need("selected-face-fermi-factor", graph[1, 1] == 1/(1+x*x),
     "Face1 covariance (1+exp(2 W B))^-1; opposite face is complementary")
need("face-is-not-hardy-projector", s.simplify(graph[1, 1]**2-graph[1, 1]) != 0,
     "Every finite real mode has strictly mixed face covariance")

lam, r = s.symbols("lambda r", positive=True)
for c in (s.I, -s.I):
    t = -c*r
    p = s.Matrix([[1, s.conjugate(t)], [t, r*r]])/(1+r*r)
    b = s.diag(lam, -lam)
    comm = b*p-p*b
    need("crosscap-nonstationarity-"+str(c),
         s.simplify(comm[0, 1]-2*lam*s.conjugate(t)/(1+r*r)) == 0
         and s.simplify(comm[1, 0]+2*lam*t/(1+r*r)) == 0
         and comm != s.zeros(2),
         "Signed boundary B, finite collar r=exp(-2 L lambda)>0; no all-P1 exclusion")
    h = s.diag(0, 2*lam)
    mean = s.simplify(s.trace(h*p))
    var = s.simplify(s.trace(h*h*p)-mean*mean)
    need("crosscap-energy-variance-"+str(c),
         s.simplify(mean-2*lam*r*r/(1+r*r)) == 0
         and s.simplify(var-4*lam*lam*r*r/(1+r*r)**2) == 0,
         "Pair sector energy test, not a W8 field identification")

out = dict(status="PASS", exact_checks=len(checks), checks=checks,
           physical_gates_closed=[], complete_TFPT_solution=False,
           scope="Original elliptic cut under declared free chiral CAR interpretation; source transport remains open")
text = json.dumps(out, ensure_ascii=False, sort_keys=True, indent=2)+"\n"
Path(__file__).with_name("boundary-source-results.json").write_text(text)
print(text, end="")
