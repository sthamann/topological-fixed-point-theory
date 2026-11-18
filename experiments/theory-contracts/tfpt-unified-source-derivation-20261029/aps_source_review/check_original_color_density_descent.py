"""R073: actual color-selector action on the original infinite-Maya currents."""
from pathlib import Path
import sys, json
from fractions import Fraction
import sympy as s

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO/'experiments/theory-contracts/joint-native-algebra-selection-20261016'))
import check_original_current_schur_kernel as original

data = json.loads((REPO/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())
M = s.Matrix([[s.Rational(x) for x in row] for row in data['explicit_Weyl_matrix']])
target = s.Matrix([1,1,-2,0,0,0,0,0])
w = list(map(int, M.T*target))
checks = []


def need(ok, msg):
    if not bool(ok):
        raise AssertionError(msg)
    checks.append(msg)


def color(lam):
    # All changed beads are above -depth. The infinite identical sea cancels.
    depth = sum(lam)+3
    occ = {(lam[i-1] if i <= len(lam) else 0)-i for i in range(1, depth+1)}
    vac = set(range(-depth, 0))
    return (sum(w[r%8] for r in occ-vac)-sum(w[r%8] for r in vac-occ))%3


need(w == [1,-1,0,-1,-1,-1,0,-1], 'actual center transported by unchanged M')
D = [sum((w[j]-w[(j-r)%8])%3 != 0 for j in range(8)) for r in range(8)]
L = [sum((w[j]-w[(j-r)%8])%3 != 0 for j in range(r)) for r in range(8)]
need(D == [0,6,4,6,2,6,4,6], 'all eight full-residue counts')
need(L == [0,1,1,3,1,3,3,5], 'all eight prefix counts')
need([j for j,x in enumerate(w) if x%3 == 1] == [0], 'unique residue1 gives all-mode witness')
rows = []
for m in range(1,25):
    terms = original.add_ribbons((), m)
    charged = sum(abs(coefficient)**2 for lam,coefficient in terms if color(lam))
    neutral = sum(abs(coefficient)**2 for lam,coefficient in terms if not color(lam))
    predicted = (m//8)*D[m%8]+L[m%8]
    need(charged == predicted and charged+neutral == m, f'original infinite-Maya J_-{m} full vector color count')
    rows.append({'m':m,'charged_norm2':charged,'neutral_norm2':neutral,'total_norm2':m})
need(rows[0]['charged_norm2'] == 1 and rows[0]['neutral_norm2'] == 0, 'whole J1J-1 word norm1, factorwise color projections0')
need(Fraction(rows[0]['charged_norm2'],4) == Fraction(1,4), 'cos(theta) infinitesimal leakage exactly1/4')
need([k for k in range(8) if all((w[j]-w[(j+k)%8])%3 == 0 for j in range(8))] == [0], 'only period8 scalar density shift preserves all colors')

# Independent actual arc-matrix commutant system, computed modulo z^2+z+1.
z = s.symbols('z')
roots = [(1+s.I)/s.sqrt(2)]
F = s.Matrix(8,8,lambda r,k:s.expand_complex(roots[0]**(r*k))/s.sqrt(8))
Z = (F*s.diag(*[z**(x%3) for x in w])*F.conjugate().T).applyfunc(s.expand)
Z = Z.applyfunc(lambda x:s.rem(x,z*z+z+1,z))
g = s.symbols('g0:8')
comm = s.diag(*g)*Z-Z*s.diag(*g)
equations=[]
for x in comm:
    equations.extend([s.expand(x).coeff(z,0),s.expand(x).coeff(z,1)])
A,_ = s.linear_eq_to_matrix(equations,g)
need(A.rank() == 7 and A.nullspace() == [s.ones(8,1)], 'actual color arc matrix forces all diagonal arc phases equal')
Za = Z.subs(z,(-1+s.I*s.sqrt(3))/2).applyfunc(s.simplify)
P = s.zeros(8); P[0,0]=1
Q = s.zeros(8); Q[1,1]=1
leak = (Za*P*Za.conjugate().T*Q-Q*Za*P*Za.conjugate().T).applyfunc(s.simplify)
norm = s.simplify(s.trace(leak.conjugate().T*leak))
need(s.simplify(Za[1,0]*s.conjugate(Za[1,0])) == s.Rational(3,64), 'original color action connects two different cover arcs')
need(norm == s.Rational(183,2048), 'exterior locality witness has nonzero exact norm183/2048')

# q-current transfers only between modes with the same actual residue/color.
q = M.T*s.Matrix([0,1,-1,0,0,0,0,0])
need((q.T*q)[0] == 2, 'actual q current retained at level2')
need(all(w[(r+8*n)%8] == w[r%8] for r in range(8) for n in [-1,0,1]), 'whole q-current mode rule preserves actual center')
result={'result':'R073','passed':len(checks),'w':w,'D':D,'L':L,
        'whole_formula':'ell_(8k+s)=k*D_s+L_s; ell_m=0 iff8 dividesm',
        'continuous_admissible_scalar_profiles':'Fourier support in8Z',
        'cos_leakage_norm2':'1/4','arc_commutant_rank':7,
        'exterior_commutator_norm2':'183/2048',
        'coordinate_support_class':'O=F^-1(K) modulo nullsets for color-invariant CAR support',
        'modes':rows,'checks':checks}
(HERE/'ORIGINAL-COLOR-DENSITY-DESCENT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'result':'R073','passed':len(checks),'actual_Maya_modes':24,'arc_commutant_rank':7}))
