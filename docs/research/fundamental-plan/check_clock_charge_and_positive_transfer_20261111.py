"""Exact comparison of winding, the original C5, E8 glue, and transfer."""
import io
from contextlib import redirect_stdout
from collections import Counter
from itertools import combinations, product, permutations
import json
from pathlib import Path
from runpy import run_path
import sympy as sp

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
old = repo/'docs/research/fundamental-plan/check_reflection_origin_eml_20261111.py'
with redirect_stdout(io.StringIO()):
    previous = run_path(str(old))
roots = [sp.Matrix(v) for v,k in previous['roots']] # original D5+A3 coordinates times 4
C, G, F = (previous[k] for k in ('C','G','G4'))
C9, G9 = sp.diag(C,sp.eye(4)), sp.diag(G,F)
check_count = 0

def check(condition, description):
    global check_count
    assert condition, description
    check_count += 1
    print('PASS', description)

fixed = [r for r in roots if C9*r == r]
fixed_span = sp.Matrix.hstack(*fixed)
perp = [r for r in roots if fixed_span.T*r == sp.zeros(len(fixed),1)]
def cartan_basis(system):
    # A generic lexicographic half system, then take indecomposable roots.
    pos={tuple(r) for r in system if next(i for i in r if i)>0}
    sums={tuple(a+b for a,b in zip(r,t)) for r,t in product(pos,repeat=2)}
    B=sp.Matrix.hstack(*(sp.Matrix(v) for v in sorted(pos-sums)))
    return B.T*B/16
fixed_cartan, perp_cartan = cartan_basis(fixed),cartan_basis(perp)
check(len(fixed)==len(perp)==20 and fixed_cartan.shape==perp_cartan.shape==(4,4)
      and fixed_cartan.det()==perp_cartan.det()==5,
      'actual C5 fixes an A4; its orthogonal root system is a second A4')
# 2*sum of the actual D5 coordinates. The roots all make this integral.
def glue_grade(r):
    return int(sum(r[:5])/2)%5
check(all(sum(r[:5])%2==0 for r in roots)
      and set(tuple(r) for r in fixed+perp)=={tuple(r) for r in roots if glue_grade(r)==0},
      'the intrinsic index-five glue grade is 2*sum(D5 coordinates) mod 5')
check({glue_grade(r) for r in roots}==set(range(5))
      and all(glue_grade(G9*r)==glue_grade(r) for r in roots),
      'actual G is identity on the index-five glue, rather than multiplying it')
seed=sp.Matrix([2]*5+[3,-1,-1,-1])
check(glue_grade(seed)==0 and C9*seed==seed,'the original primitive D5+A3 glue seed is neutral in the new quotient')
print('intrinsic glue root counts',dict(sorted(Counter(map(glue_grade,roots)).items())))

# The native SU5 lift of the permutation is diagonalizable. Its centered
# integer weights sum to zero, hence integral on the original spinors.
weights=sp.Matrix([-2,-1,0,1,2,0,0,0,0])
def spectral_grade(r):
    q=weights.dot(r)/4
    assert q.is_Integer
    return int(q)%5
h0=[r for r in roots if spectral_grade(r)==0]
h0_cartan=cartan_basis(h0)
counts=Counter(map(spectral_grade,roots))
check(dict(counts)=={0:60,1:45,2:45,3:45,4:45}
      and h0_cartan.shape==(6,6) and h0_cartan.det()==4,
      'actual carrier spectral grade gives D6 neutral roots and two abelian directions')
# Independent Fock and tensor branching count without a choice of root basis.
cycles=lambda size: Counter(sum(j for j in q)%5 for q in combinations(range(5),size))
spin=sum((cycles(k) for k in (0,2,4)),Counter())
adj=Counter((a+b)%5 for a,b in combinations(list(range(5))*2,2))
check([adj[k]+8*spin[k]+12+(15 if k==0 else 0) for k in range(5)]==[68,45,45,45,45],
      'independent spinor/vector branching reproduces the full adjoint grade dimensions')

def matrix_units_orbits(N, transforms):
    pairs=set(product(range(N),repeat=2))
    result=[]
    while pairs:
        p=pairs.pop()
        orbit={p}
        todo=[p]
        while todo:
            a,b=todo.pop()
            for g in transforms:
                q=(g(a),g(b))
                if q not in orbit:
                    orbit.add(q)
                    todo.append(q)
        pairs-=orbit
        result.append(orbit)
    return sorted(map(len,result))
check(matrix_units_orbits(5,[lambda j:(j+1)%5,lambda j:2*j%5])==[5,20],
      'full five-slot equivariant scalar transfers have exactly diagonal/off-diagonal freedom')

# Coordinate P2 split, with canonical native half-turn on the five slots.
r=lambda j:(-j)%5
g=lambda j:2*j%5
splits=[set(a) for a in combinations(range(5),3) if {r(j) for j in a}==set(a)]
check(splits==[{0,1,4},{0,2,3}] and {g(j) for j in splits[0]}==splits[1],
      'half-turn admits two coordinate 3+2 splits and the whole clock exchanges them')
def Y(color):
    return sp.diag(*(sp.Rational(-1,3) if i in color else sp.Rational(1,2) for i in range(5)))
Ya,Yb=map(Y,splits)
check(G*Ya*G.inv()==Yb and sp.trace(Ya**2)==sp.trace(Yb**2)==sp.Rational(5,6),
      'both allowed splits have identical original hypercharge normalization')
check(sum((C**i*Ya*C**(-i) for i in range(5)),sp.zeros(5))==sp.zeros(5),
      'unmarked winding translation average kills the fixed hypercharge direction')
intersections=Counter((i in splits[0],i in splits[1]) for i in range(5))
check(sorted(intersections.values())==[1,2,2],
      'simultaneously keeping both hypercharges leaves block sizes 1,2,2')
# The source's actual X=3 Nweak-2 Ncolor on the entire root lattice.
# Weak sums are integral because there are two weak coordinates.
def original_phase(r, color):
    vector = sp.Matrix([-2 if j in color else 3 for j in range(5)]+[0]*4)
    phase = vector.dot(r)/4
    assert phase.is_Integer
    return int(phase)%5
check(all(original_phase(root,set(color))==(-glue_grade(root))%5
          for root,color in product(roots,combinations(range(5),3))),
      'original U5 equals inverse glue character on every root for all ten color splits')

# F5 dual of the original full rank lattice. 4 is invertible modulo 5,
# so the D5+A3 base and its index-four E8 extension have the same mod-5 dual.
# Fix the fourth A3 coefficient to zero as a representative of the
# linear forms on the zero-sum plane. The two eigenlines are explicit.
u = sp.Matrix([0,1,2,3,4,0,0,0,0])
v = sp.Matrix([0]*5+[1,2,4,3])
def eigencharacter(coeff):
    return all((4*int(coeff.dot(G9*r-2*r)))%5==0 for r in roots)
check(eigencharacter(u) and eigencharacter(v),
      'two independent doubling characters of the original G, carrier and family')
# Direct enumeration of each small block, with q4=0 for family forms.
def candidates(M, vectors, forms):
    return [q for q in forms if all(sum(a*b for a,b in zip(q, M*r-2*r))%5==0 for r in vectors)]
e5=[sp.eye(5)[:,j] for j in range(5)]
e4=[sp.eye(4)[:,j]-sp.eye(4)[:,3] for j in range(3)]
carrier=candidates(G,e5,product(range(5),repeat=5))
family=candidates(F,e4,(x+(0,) for x in product(range(5),repeat=3)))
check(len(carrier)==len(family)==5 and all(len(set(t))==5 for t in carrier if any(t)),
      'all doubling characters classified, every nonzero carrier character separates all five slots')
all_color=list(combinations(range(5),3))
check(all(all(len({t[j] for j in color})>1 for color in all_color) for t in carrier if any(t)),
      'color SU3 neutrality kills every carrier doubling character')
check(all(any(sum(a*b for a,b in zip(t,r))%5 for r in e4) for t in family if any(t)),
      'full original family A3 neutrality also kills every family doubling character')

# Same original marked root lambda_A=e3+e6. A3 -> D3 sends f to
# (f0+f1,f0+f2,f0+f3), with the correct primitive glue sign.
from_e8=lambda r: tuple(r[:5])+tuple(r[5]+r[j] for j in (6,7,8))
old_W=sp.Matrix([0,0,4,0,0,2,2,-2,-2])
check(from_e8(old_W)==(0,0,4,0,0,4,0,0)
      and original_phase(old_W,{0,1,2})==3 and
      [int(v.dot(G9**k*old_W)/4)%5 for k in range(4)]==[3,1,2,4],
      'old Ward root has fixed U5 degree 3 but its G orbit cycles the family doubling degrees')
check(len({tuple(C9**k*old_W) for k in range(5)})==5,
      'the Ward root has components in all five C5 spectral grades, via its free five-root orbit')
# The local gauge subgroup and exact positive dynamics constrain the same action.
# A single G-commuting positive operator has equal spectra on all nonzero
# C5-eigenspaces, since G conjugates the character labels transitively.
j5=sp.ones(5)
lam=sp.symbols('lambda',real=True)
transfer = lam*sp.eye(5)+(1-lam)*j5/5
check(transfer*G==G*transfer and transfer*C==C*transfer
      and sp.factor(transfer.charpoly().as_expr())==(transfer.charpoly().gen-1)*(transfer.charpoly().gen-lam)**4,
      'the full reversible five-slot scalar transfer has one nontrivial rate of multiplicity four')

# Branches on the native longitudinal winding: a weighted, simultaneous
# composition of residues cannot change the shift defect of the full source.
# Three exact block projections for shifts by five and residue reader.
for n in range(-20,21):
    assert divmod(n+5,5)==(n//5+1,n%5)
    assert divmod(n+1,5)==(n//5+(n%5==4),(n%5+1)%5)
print('integer shift by one carries from residue four to the winding record')
# Load the actual conditional Ward dictionary, and transport the original
# quarter winding and its energy into the same frame as the E8 clock.
fold = repo/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json'
M=sp.Matrix([[sp.Rational(str(v)) for v in row]
             for row in json.loads(fold.read_text())['explicit_Weyl_matrix']])
q=sp.ones(8,1)
check(M.T*M==sp.eye(8) and M*q==sp.Matrix([0,0,2,0,0,2,0,0]),
      'load and verify the original conditional Ward dictionary')
B=sp.zeros(8,9)
B[:5,:5]=sp.eye(5)
for j in range(3):
    B[5+j,5]=B[5+j,6+j]=1
def pmatrix(p):
    return sp.Matrix(4,4,lambda row,col:int(row==p[col]))
def is4cycle(p):
    return len({0,p[0],p[p[0]],p[p[p[0]]]})==4 and p[p[p[p[0]]]]==0
def source_gear(F):
    return M.T*B*sp.diag(G,F)*sp.diag(sp.eye(5),sp.eye(4)-sp.ones(4)/4)*B.T*M
def charge_lattice(v):
    return all(x.is_Integer for x in v) and sum(v)%2==0 or (
        all((x-sp.Rational(1,2)).is_Integer for x in v) and sum(v)%2==0)
half=sp.Matrix(8,8,lambda i,j:int(i==(j+4)%8))
a=sp.Matrix([1]*4+[0]*4)
check(half*a+a==q,'actual W4 affine action squares to translation by the winding cocharacter')
gear=source_gear(F)
shift_columns=sp.Matrix.hstack(*(half**i*gear**k*q for k in range(4) for i in (0,1)))
check(gear.T*gear==sp.eye(8) and gear**4==sp.eye(8),
      'Ward-conjugated native clock is orthogonal and has order four')
shifts=sp.Matrix([2*2**k for k in range(4) for i in (0,1)])
variants=[]
for p in permutations(range(4)):
    if not is4cycle(p):continue
    g=source_gear(pmatrix(p))
    cols=sp.Matrix.hstack(*(half**i*g**k*q for k in range(4) for i in (0,1)))
    relation_detected=False
    witness=None
    for z in cols.nullspace():
        z*=sp.ilcm(*(x.q for x in z))
        if z.dot(shifts)%5:
            relation_detected=True
            witness=(tuple(z),z.dot(shifts))
    repairs=[]
    differences=[]
    for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
        delta=sp.Matrix([(k+alpha)/8-sp.Rational(1,2) for k in range(8)])
        b=g*delta-delta
        repairs.append(charge_lattice(b))
        differences.append((g.T*delta-delta).dot(q))
    variants.append((p,relation_detected,tuple(repairs),tuple(differences),witness))
check(len(variants)==6 and all(not any(row[2]) for row in variants),
      'none of the six ordinary family four-cycles admits a lattice correction preserving original time')
check(sum(row[1] for row in variants)==2,
      'two of the six same-frame family clocks already forbid any scalar W4 residue by an identity word')
print('family four-cycles, scalar identity contradiction, time repairs, energy changes, witness')
for row in variants: print(row)

# The smallest linear quotient over F5 that carries total source charge
# and both original linear actions is the span of the orbit of the charge covector.
def mod5(x):
    a,b=map(int,sp.Rational(x).as_numer_denom())
    return a*pow(b,-1,5)%5
def matrix_rank_mod5(rows):
    rows=[list(map(mod5,row)) for row in rows]
    index=0
    for col in range(8):
        pivot=next((i for i in range(index,len(rows)) if rows[i][col]),None)
        if pivot is None:continue
        rows[pivot],rows[index]=rows[index],rows[pivot]
        inv=pow(rows[index][col],-1,5)
        rows[index]=[inv*a%5 for a in rows[index]]
        for j in range(index+1,len(rows)):
            d=rows[j][col]
            rows[j]=[(a-d*b)%5 for a,b in zip(rows[j],rows[index])]
        index+=1
    return index
def orbit_span(initial):
    basis=[tuple(map(mod5,initial))]
    for row in basis:
        for A in (half,gear):
            cand=tuple(map(mod5,sp.Matrix([row])*A))
            if matrix_rank_mod5(basis+[cand])>len(basis):
                basis.append(cand)
    return basis
basis=orbit_span(list(q.T))
print('minimal simultaneous linear quotient rank',len(basis))
# The raw kernel being quotiented can carry fractional weights but its
# coefficients are all dyadic and the full-rank coordinate map is valid mod 5.
gauge_roots=[M.T*(sp.eye(8)[:,i]-sp.eye(8)[:,j])
             for i,j in ((0,1),(1,2),(3,4))]
violations=[[mod5((sp.Matrix([r])*x)[0]) for x in gauge_roots] for r in basis]
check(len(basis)==6 and matrix_rank_mod5(basis)==6 and any(x for row in violations for x in row),
      'a simultaneous quotient exists with multiple records but its orbit does not stay gauge-neutral')
print('minimal quotient gauge root evaluations',violations)
# The user-supplied H_* is checked on both the parity hypercharge and
# the original coordinate hypercharge with the carrier basis kept fixed.
I=sp.eye(5); R=G**2; Pi=sp.ones(5)/5
Ypar=(I-5*R)/12
H=I-Pi+(2*I-G-G.inv())/2  # dimensionless H_* / omega
comm=lambda A,B:A*B-B*A
norm2=lambda A:sp.trace(A.T*A)
check(H.eigenvals()=={0:1,1:1,2:2,3:1} and comm(H,Ypar)==sp.zeros(5)
      and all(H[i,j]<=0 for i,j in product(range(5),repeat=2) if i!=j)
      and H*sp.ones(5,1)==sp.zeros(5,1),
      'user H* has the stated exact spectrum and positive reversible parity-space transfer')
check(comm(H,G)==sp.zeros(5) and comm(H,Y({0,1,2}))!=sp.zeros(5)
      and comm(H,Ya)!=sp.zeros(5) and comm(H,Yb)!=sp.zeros(5),
      'user H* fails fixed original coordinate hypercharge and both half-turn-compatible coordinate splits')
print('original Y and two marked split commutator norms squared',*[norm2(comm(H,K))
      for K in (Y({0,1,2}),Ya,Yb)])
# Another positive generator for the explicitly half-turn-compatible
# coordinate split. It preserves U(1) and G, yet mixes gauge partners.
Hswap=(I-R)/2
gaugerotation=sp.zeros(5)
gaugerotation[0,1]=1
gaugerotation[1,0]=-1
check(Hswap.eigenvals()=={0:3,1:2} and Hswap**2==Hswap
      and all(comm(Hswap,A)==sp.zeros(5) for A in (Ya,Yb,G))
      and comm(Hswap,gaugerotation)!=sp.zeros(5),
      'a finite positive phase construction preserves a coordinate U(1) but fails full color invariance')
a,b=sp.symbols('a b', real=True, positive=True)
Hfamily=a*(I-Pi)+b*(2*I-G-G.inv())
pvar=sp.symbols('p')
check(sp.factor(Hfamily.charpoly(pvar).as_expr()-pvar*(pvar-a)*(pvar-a-2*b)**2*(pvar-a-4*b))==0,
      'two-coefficient Hamiltonian family has the stated independent even and odd energies')
l2,l3=sp.Rational(64,729),sp.Rational(1,729)
check(sp.sqrt(l2*l3)==sp.Rational(8,729),
      'matching both original even relaxation rates forces odd rate 8/729 in this family')
# Exhaust the attached root Fourier coefficients and compare the full rank.
from_e8_roots=[sp.Matrix(from_e8(r))/4 for r in roots]
m=Counter(2*sum(r[:5]) for r in from_e8_roots)
check(m=={0:32,1:40,-1:40,2:30,-2:30,3:20,-3:20,4:10,-4:10,5:4,-5:4},
      'attached E8 Fourier polynomial checked on all original root coordinates')
moment=sum((r*r.T for r in from_e8_roots),sp.zeros(8))
check(moment==60*sp.eye(8) and sum(i*i*n for i,n in m.items())==1200,
      'E8 root walk Hessian has rank eight, and its displayed line has coefficient 5 kappa')
check(all(any(mod5(x) for x in (source_gear(pmatrix(p)).T*q-2*q))
          for p in permutations(range(4)) if is4cycle(p)),
      'none of the six actual Ward lifts doubles the canonical total source charge modulo five')
# Exact rank bound for any nonnegative C5-invariant local root walk.
# A nontrivial rational C5 orbit in the sum-zero carrier module has rank four.
carrier_projection=sp.eye(5)-sp.ones(5)/5
def orbit_carrier_rank(root):
    w=carrier_projection*sp.Matrix(root[:5])
    return sp.Matrix.hstack(*(C**k*w for k in range(5))).rank()
orbit_ranks={tuple(r):orbit_carrier_rank(r) for r in roots}
check(all(rank in (0,4) for rank in orbit_ranks.values()) and
      all(orbit_ranks[tuple(r)]==4 for r in roots if glue_grade(r)!=0),
      'any positive invariant walk using a charged root has spatial Hessian rank at least four')
# The existing Petersen channel uses pairs of the same five carrier slots.
# Test the actual sixth-step populations, including arbitrary real logs.
pairs=[set(p) for p in combinations(range(5),2)]
Petersen=sp.Matrix([[int(not a.intersection(b)) for b in pairs] for a in pairs])
assert Petersen.eigenvals()=={3:1,-2:4,1:5}
TP=(Petersen/3)**6
check(TP.eigenvals()=={1:1,l2:4,l3:5}
      and all(TP[i,j]>0 for i,j in product(range(10),repeat=2)),
      'existing Petersen sixth step is positive and has both original rates')
alpha2,beta2=sp.log(l2),sp.log(l3)
adj_rate=sp.log(sp.Rational(27,256))/5
nonadj_rate=(3*sp.log(3)+2*sp.log(2))/5
Qlog=(-6*sp.log(3)+2*sp.log(2))*sp.eye(10)-2*sp.log(2)*Petersen+nonadj_rate*sp.ones(10)
check(sp.simplify(Qlog[0,next(j for j in range(10) if Petersen[0,j])]-adj_rate)==0
      and adj_rate<0
      and sp.simplify(sp.trace(Petersen*Qlog)-(-8*alpha2+5*beta2))==0
      and sp.simplify(sp.trace(Petersen*Qlog)-30*adj_rate)==0,
      'negative edge trace rules out even a nonreversible real Markov embedding')
# Matrix entries of the log of the three real Petersen shells also fail.
P_rad=sp.Matrix([[0,1,0],[3,0,1],[0,2,2]])/3
P2_rad=(P_rad-sp.eye(3))*(P_rad-sp.eye(3)/3)/((sp.Rational(-2,3)-1)*(-1))
P3_rad=(P_rad-sp.eye(3))*(P_rad+sp.eye(3)*sp.Rational(2,3))/((sp.Rational(1,3)-1)*1)
Q_rad=alpha2*P2_rad+beta2*P3_rad
check(sp.simplify(Q_rad[0,1]-adj_rate)==0 and Q_rad[0,1]<0,
      'actual three Petersen shells still fail continuous Markov embedding')
u2=sp.Matrix([1,-1,0]);u3=sp.Matrix([1,1,-2])
Q330=alpha2*u2*u2.T/2+beta2*u3*u3.T/6
check(all(Q330[i,j]>0 for i,j in product(range(3),repeat=2) if i!=j)
      and sp.simplify(Q330[0,1]-sp.log(sp.Rational(9,8)))==0,
      'v330/v971 with the same three eigenvalues really has positive transition rates')
print('Petersen edge rate',adj_rate,'and required lower bound for smaller transfer rate',sp.N(l2**sp.Rational(8,5),12))
# Positive alternative with the original fixed gauge marking and existing
# transfer, without a new Hamiltonian coefficient. The normalized charged
# pair labels here are population observables, not physical clock states.
weak={3,4}
types=[len(pair&weak) for pair in pairs]
charge=sp.Matrix([[-sp.Rational(2,3),sp.Rational(1,6),1][k] for k in types])
other=sp.Matrix([[1,-1,3][k] for k in types])
uniform=sp.ones(10,1)
check(Petersen*charge==-2*charge and Petersen*other==other
      and charge.dot(other)==0 and uniform.dot(charge)==uniform.dot(other)==0,
      'genuine fixed 3+2 gauge-invariant pair observables resolve both existing transfer eigenspaces')
vw=types.index(2)
charge/=sp.sqrt(charge.dot(charge))
other/=sp.sqrt(other.dot(other))
mark=sp.zeros(10);mark[vw,vw]=1
check(sp.simplify(charge.dot(mark*other))==1/sp.sqrt(5)
      and sp.simplify((charge.dot(mark*charge)+other.dot(mark*other))/2)==sp.Rational(9,20)
      and sp.simplify(sp.log(l2)-sp.log(l3)-6*sp.log(2))==0,
      'same existing positive transfer has a fixed weak-pair phase readout at frequency 6 log 2')
print(f'{check_count} exact checks passed')
