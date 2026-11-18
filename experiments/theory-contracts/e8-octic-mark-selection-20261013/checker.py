from importlib.util import spec_from_file_location, module_from_spec
from fractions import Fraction
from pathlib import Path
import json
import sys
import numpy as np
import sympy as sp

SOURCE = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4/verification/v1_e8_glue.py')
sys.path.insert(0, str(SOURCE.parent))
spec = spec_from_file_location('tfpt_v1_e8_glue', SOURCE)
mod = module_from_spec(spec)
spec.loader.exec_module(mod)
# Consume the unchanged repository root generator. Convert its binary halves to exact doubled integers.
roots2 = [tuple(int(round(2*x)) for x in beta) for beta in mod.e8_roots()]
if not (len(roots2) == 240 and all(sum(x*x for x in b) == 8 for b in roots2)):
    raise RuntimeError('Native E8 root count or norm differs')

zsyms = sp.symbols('z0:8')
def vals(z):
    # z may be any exact SymPy vector; raw homogeneous moments.
    dots = [sp.Rational(1,2)*sum(bi*zi for bi, zi in zip(beta, z)) for beta in roots2]
    return {m: sp.factor(sum(d**m for d in dots)) for m in (2,4,6,8,10,12,14)}
def grad_hess(z):
    dots2 = [sum(bi*zi for bi, zi in zip(beta, z)) for beta in roots2] # doubled-root dot
    g = [sp.factor(8*sum(sp.Rational(beta[i],2)*(sp.Rational(1,2)*d)**7 for beta,d in zip(roots2,dots2))) for i in range(8)]
    H = sp.Matrix(8,8,lambda i,j: sp.factor(56*sum(sp.Rational(beta[i]*beta[j],4)*(sp.Rational(1,2)*d)**6 for beta,d in zip(roots2,dots2))))
    return sp.Matrix(g), H

def analyze(label, coords):
    z = sp.Matrix(coords)
    z = z / sp.sqrt((z.dot(z)))
    ms = vals(z)
    g,H = grad_hess(z)
    lam = sp.simplify(z.dot(g))
    tangent = sp.simplify(g-lam*z)
    P = sp.eye(8)-z*z.T
    constrained = sp.simplify(P*(H-lam*sp.eye(8))*P)
    # E8 roots perpendicular to z; their reflections generate the Weyl stabilizer.
    zint = sp.Matrix(coords)
    ortho = [b for b in roots2 if sum(bi*zi for bi,zi in zip(b,zint)) == 0]
    eig = np.linalg.eigvalsh(np.array(constrained.evalf(), dtype=float))
    # connected components of roots in the stabilizer under nonorthogonality
    rem=set(range(len(ortho))); comps=[]
    while rem:
        stack=[rem.pop()]; comp=[]
        while stack:
            i=stack.pop(); comp.append(i)
            linked=[j for j in list(rem) if sum(a*b for a,b in zip(ortho[i],ortho[j])) != 0]
            for j in linked: rem.remove(j); stack.append(j)
        comps.append(len(comp))
    exact_eigs=constrained.eigenvals()
    return {
      'name':label,'unit_z':[str(x) for x in z],
      'moments_unit':{str(k):str(v) for k,v in ms.items()},
      'gradient_parallel_eigenvalue':str(lam),
      'gradient_tangent_norm2':str(sp.factor(tangent.dot(tangent))),
      'perpendicular_E8_roots':len(ortho),'orthogonal_root_component_sizes':sorted(comps,reverse=True),
      'projected_hessian_exact_eigenvalues':{str(k):v for k,v in exact_eigs.items()},
      'tangent_hessian_eigenvalues_numeric':[round(float(x),8) for x in eig],
      'tangent_hessian_minmax':[round(float(eig[0]),8),round(float(eig[-1]),8)]
    }

# Exact tests at obvious partition vectors. The vectors are explicit hypotheses, not source-selected inputs.
marks = [
 ('5+3 zero-sum partition (3^5,-5^3)', [3]*5+[-5]*3),
 ('3+2 zero-sum submark (2^3,-3^2,0^3)', [2]*3+[-3]*2+[0]*3),
 ('root-axis direction', [1,0,0,0,0,0,0,0]),
 ('equal-coordinate direction', [1]*8),
]
result={'source_file':str(SOURCE),'roots':len(roots2),'roots_encoding':'unchanged verification.v1_e8_glue.e8_roots(), doubled to exact integer coordinates','marks':[analyze(n,z) for n,z in marks]}

# On each permutation-fixed split R^k x R^(8-k), S8 on the unit sphere reduces to a rational
# function of t=b/a. Find exact stationary ratios of this one-dimensional fixed stratum.
t=sp.symbols('t', real=True)
strata=[]
for k in range(1,8):
    raw=[sp.Integer(1)]*k+[t]*(8-k)
    f=sp.factor(vals(raw)[8] / (sp.Integer(k)+(8-k)*t**2)**4)
    num=sp.factor(sp.together(sp.diff(f,t)).as_numer_denom()[0])
    # remove factors which only encode projective endpoints/repeated-value degeneracies
    fac=sp.factor(num)
    strata.append({'k':k,'S8_unit_rational':str(f),'stationarity_numerator_factorization':str(fac)})
result['two_level_fixed_strata']=strata
# Finite Cartan holonomy versus its infinitesimal generator. Turns t=(0^5,1/2^3).
# For beta=b/2, 4*t.beta=sum(b[5:]); the character is i**grade.
from collections import Counter
result['finite_order4_character']={
 'turns':'(0,0,0,0,0,1/2,1/2,1/2)',
 'root_grade_counts':dict(sorted(Counter(sum(b[5:])%4 for b in roots2).items())),
 'continuous_generator_orthogonal_roots':sum(sum(b[5:])==0 for b in roots2),
 'centralizer_scope':'Finite phase identity yields D5+D3 (rank8); infinitesimal orthogonality yields D5+A2 (rank7). The 5+3 split is an input to this character.'
}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
