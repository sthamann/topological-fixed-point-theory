"""Reflected two-current Gram on actual charged affine E8 states."""
from contextlib import redirect_stdout
from itertools import combinations, product
from pathlib import Path
import io
import runpy
import sys

repo=Path(sys.argv[1]).resolve()
prior=repo/"experiments/theory-contracts/recursive-mirrored-source-20261103/check_seam_yukawa_reflection_20261111.py"
with redirect_stdout(io.StringIO()):
    baseline=runpy.run_path(str(prior))
data=baseline["a"]
sp=data["sp"]
charge=data["charge"]
cocycle=data["epsilon"]
vertex=data["vertex_response"]
alpha,beta,gamma=data["alpha"],data["beta"],data["complement"][0]
a,g=sp.symbols("a g",real=True)
A=tuple((root,a*c) for root,c in data["hermitian_roots"](alpha))+tuple(
    (root,g*c) for root,c in data["hermitian_roots"](gamma))
B=data["hermitian_roots"](beta)

def opposite(root):
    return tuple(-x for x in root)

def roots_of_slot(slot):
    roots=[data["matter_vertex"](slot,i)[0] for i in range(3)]
    return roots+[roots[0][:5]+(-sp.Rational(1,2),)*3]

def at_infinity(roots,positions,oracle=False):
    if any(sum(column)!=0 for column in zip(*roots)):
        return sp.S.Zero
    if oracle:
        transformed=(sp.S.One,)+tuple(x/(x+1) for x in positions)
        scale=sp.prod((x+1)**-2 for x in positions)
        return vertex(roots,transformed)*scale
    result=sp.S.One
    for i,j in combinations(range(len(roots)),2):
        result*=cocycle(roots[i],roots[j])
        if i>0:
            result*=(positions[i-1]-positions[j-1])**charge(roots[i],roots[j])
    return result

def ordered_gram(slot,bra,ket,t,s,oracle=False):
    momenta=roots_of_slot(slot)
    positions=(1/s,1/t,t,s,sp.S.Zero)
    result=sp.zeros(len(momenta))
    for i,j in product(range(len(momenta)),repeat=2):
        p,q=momenta[i],momenta[j]
        bra_sign=cocycle(opposite(p),p)
        for (r,rc),(v,vc),(x,xc),(y,yc) in product(bra[1],bra[0],ket[0],ket[1]):
            result[i,j]+=bra_sign*rc*vc*xc*yc*at_infinity(
                (opposite(p),r,v,x,y,q),positions,oracle)
    return result.applyfunc(lambda z:sp.factor(z/(t*s)**2))

words=((A,B),(B,A))
slot=data["slots"]["e"][0][0]

# The compact adjoint carries the reflected operators and state together.
# Infinity has an independent exact Möbius comparison against the native oracle.
t,s=sp.symbols("t s",positive=True)
G_AB=ordered_gram(slot,words[0],words[0],t,s)
G_BA=ordered_gram(slot,words[1],words[1],t,s)
G_cross=ordered_gram(slot,words[0],words[1],t,s)
for matrix in (G_AB,G_BA,G_cross):
    assert matrix==sp.diag(*matrix.diagonal())
    assert sp.cancel((matrix[1,1]-matrix[2,2]).subs(g,a))==0
supports=[{tuple(p_i+r_i+q_i for p_i,r_i,q_i in zip(p,r,q))
           for (r,_),(q,_) in product(A,B)} for p in roots_of_slot(slot)]
assert all(supports[i].isdisjoint(supports[j]) for i,j in combinations(range(4),2))
print("Exact charge supports of all four inputs are disjoint for both orders.")
print("The full and cross Grams are diagonal for all radii, with no 3-to-4 Schur correction.")

sample=(sp.Rational(4,5),sp.Rational(1,5))
replay={t:sample[0],s:sample[1]}
assert (ordered_gram(slot,*words,*sample,oracle=True)-G_cross.subs(replay)).applyfunc(sp.cancel)==sp.zeros(4)
all_matter=sorted({i for rows in data["slots"].values() for l,r,h in rows for i in (l,r)})
for idx in all_matter:
    for word, expected in zip(words,(G_AB,G_BA)):
        assert (ordered_gram(idx,word,word,*sample)-expected.subs(replay)).applyfunc(sp.cancel)==sp.zeros(4)
M=G_AB+G_BA
one={a:1,g:1}
scale_ab=(G_AB[2,2]/G_AB[0,0]).subs(replay).subs(one)
scale_ba=(G_BA[2,2]/G_BA[0,0]).subs(replay).subs(one)
scale_sum=sp.factor((M[2,2]/M[0,0]).subs(replay).subs(one))
print("At radii (4/5,1/5), G33/G11 by order, then sum:",
      *map(float,(scale_ab,scale_ba,scale_sum)))
for cross in (M+2*G_cross,M-2*G_cross):
    assert all(diag>0 for diag in cross.subs(replay).subs(one).diagonal())
    print("Coherent +/- G33/G11 at the same radii:",float(
          (cross[2,2]/cross[0,0]).subs(replay).subs(one)))

# With equal first-current weights, the radial model has a single
# diagonal family-contrast parameter: G2=G3. Use exact lepton ratios
# and the exact quark determinant as constraints, then hold out
# all three individual quark ratios.
mp=data["mp"]
lam,phi=data["lY"],data["phi0"]
r_e=12*phi**2/7
r_mu=8*phi/7
e0=2*lam**8/3
b0=2*lam**5/3
c0=8*lam**3/3
tau0=(c0-b0)/(1-r_mu)
needed=tau0*r_e*r_mu/e0
h3_squared=(r_mu*tau0**2-b0*c0)/lam**8
assert h3_squared>0
target_ud=baseline["target_ud"]
h1_squared=target_ud*(16+9*h3_squared)/16
expr=sp.factor((M[2,2]/M[0,0]).subs(one).subs(s,1-t))
contrast=sp.lambdify(t,expr,"mpmath")
assert contrast(mp.mpf(3)/4)<needed<contrast(mp.mpf(4)/5)
radius=mp.findroot(lambda z:contrast(z)-needed,(mp.mpf(3)/4,mp.mpf(4)/5))
assert mp.mpf(3)/4<radius<mp.mpf(4)/5
assert abs(contrast(radius)-needed)<mp.mpf("1e-25")
print("Exact lepton spectrum requires contrast",mp.nstr(needed,14),
      "and is reached on the test path at",mp.nstr(radius,14))

import numpy as np
h=data["h"]
profile={h[0]:float(mp.sqrt(h1_squared)),h[1]:1,h[2]:float(mp.sqrt(h3_squared))}
E=np.array(baseline["Eq"].subs(profile),dtype=float)
U=np.array(baseline["Uq"].subs(profile),dtype=float)
rho_c=mp.root(baseline["ed"],3)
family=np.array([1,float(needed),float(needed)])
def singulars(matrix,lengths,color_norm=1):
    end=np.sqrt(float(lam)**np.array(lengths)/family)
    return np.sort(np.linalg.svd(end[:,None]*matrix*end[None,:]/float(color_norm),compute_uv=False))

leptons=singulars(E,[8,5,3])
down=singulars(E,[7,5,2],rho_c)
up=singulars(U,[7,3,0],rho_c)
assert np.allclose(leptons[:-1]/leptons[1:],[float(r_e),float(r_mu)],rtol=1e-12)
expected=np.array([55/117,34/(47*float(phi)),3/(26*float(phi)**2)])
assert np.isclose(np.prod(up/down),np.prod(expected),rtol=1e-12)
assert all(abs((up/down)/expected-1)>.3)
print("Held-out individual quark ratios u/d, c/s, t/b:",up/down,"source:",expected)
print("Conditional e/u from this failed model:",leptons[0]/up[0])

# A bound for the full family-universal diagonal Gram model with nested
# pure word factors on the real quarter chart. Lepton and down have the
# same reduced symmetric kernel, with the first direction decoupled.
tab=data["TAB"]
e_over_d_table=lam*mp.mpf(tab["e"][1])/mp.mpf(tab["d"][1])
upper_table=lam*(mp.mpf(tab["e"][1])+data["half_printing"])/(
    mp.mpf(tab["d"][1])-data["half_printing"])
assert upper_table<lam
min_eu=mp.mpf(117)*lam*rho_c/55
assert min_eu>mp.mpf(".48")
print("Nested lengths and family-universal Gram imply e/d>=lambda*rho_c, then e/u>=",
      mp.nstr(min_eu,12),"if the exact u/d is enforced.")
print("Old table e/d =",mp.nstr(e_over_d_table,12))

# Test the simplest finite radial color-endpoint weighting at old times.
# On an invariant pair at radii u and 1/u, the normalized Gram factor is
# (1+2 cosh(2*clock_time))/3; the time identification is conditional.
times=(mp.mpf(1)/8,mp.mpf(1)/8+mp.mpf(1)/20)
color_scales=tuple((1+2*mp.cosh(2*z))/3 for z in times)
target_color=(mp.root(baseline["ed_lo"],3),mp.root(baseline["ed_hi"],3))
assert color_scales[0]<target_color[0]<target_color[1]<color_scales[1]
print("Color Gram at the two old time windows:",*[mp.nstr(x,12) for x in color_scales])

# A coherent one-current/two-current branch on the original grade-one
# charged states can mix 1 and 2. Read the effect of the actual deck
# projector before claiming a physical admixture.
import json
parent=json.loads((repo/"experiments/theory-contracts/compiler-e6-seam-readout-20260918/certificate.json").read_text())
deck=dict(zip(map(tuple,parent["seam"]["root_order"]),parent["seam"]["root_character_exponents"],strict=True))
eta=sp.symbols("eta",real=True)
def zero_mode(slot,terms):
    momenta=roots_of_slot(slot)
    answer=sp.zeros(4)
    for j,p in enumerate(momenta):
        for r,c in terms:
            if charge(r,p)==-1:
                target=tuple(x+y for x,y in zip(r,p))
                answer[momenta.index(target),j]+=c*cocycle(r,p)
    return answer
for idx in all_matter:
    A0=zero_mode(idx,A)
    B0=zero_mode(idx,B)
    coherent=B0*(sp.eye(4)+eta*A0)
    coherent_gram=sp.eye(3)+(coherent.T*coherent)[:3,:3]
    assert sp.factor(coherent_gram[0,1]**2-eta**2*a**2)==0
    D=sp.diag(*[sp.I**deck[tuple(2*x for x in p)] for p in roots_of_slot(idx)[:3]])
    assert len(set(D.diagonal()))==3
    harmonic=sp.simplify(D[0,0]/D[1,1])
    assert harmonic in (sp.I,-sp.I)
    assert data["old_spin"](roots_of_slot(idx)[0])==data["old_spin"](roots_of_slot(idx)[1])
    extracted=sum((harmonic**(-j)*D**j*coherent_gram*D**(-j)
                   for j in range(4)),sp.zeros(3))/4
    assert sp.simplify(extracted[0,1]-coherent_gram[0,1])==0
    averaged=sum((D**j*coherent_gram*D**(-j) for j in range(4)),sp.zeros(3))/4
    assert averaged==sp.diag(*coherent_gram.diagonal())
print("A marked coherent one/two-current return creates 1-2 mixing; the original C4 deck average removes it.")
