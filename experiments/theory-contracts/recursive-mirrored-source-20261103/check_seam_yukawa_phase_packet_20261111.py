"""Reflected mixed word evaluated with the existing affine lattice correlator."""
from contextlib import redirect_stdout
from itertools import product
from pathlib import Path
import io
import runpy
import sys

repo=Path(sys.argv[1]).resolve()
prior=repo/"experiments/theory-contracts/recursive-mirrored-source-20261103/check_seam_yukawa_ordered_norms_20261111.py"
with redirect_stdout(io.StringIO()):
    source=runpy.run_path(str(prior))
sp=source["sp"]
A,B=source["A"],source["B"]
t,s=source["t"],source["s"]
a,g=source["a"],source["g"]
data=source["data"]

def word_gram(slot, bra, ket, oracle=False):
    momenta=source["roots_of_slot"](slot)
    bra_marks=tuple(reversed(bra))
    reflected=tuple(1/z for z,_ in bra_marks)
    positions=reflected+tuple(z for z,_ in ket)+(sp.S.Zero,)
    scale=sp.prod(z**-2 for z,_ in bra)
    out=sp.zeros(4)
    for i,j in product(range(4),repeat=2):
        p,q=momenta[i],momenta[j]
        pbar=source["opposite"](p)
        for letters in product(*(terms for _,terms in bra_marks+ket)):
            roots=(pbar,)+tuple(root for root,_ in letters)+(q,)
            out[i,j]+=sp.prod(c for _,c in letters)*source["cocycle"](pbar,p)*source["at_infinity"](roots,positions,oracle)
    return (scale*out).applyfunc(sp.factor)

slot=source["slot"]
single=((t,B),)
AB=((t,A),(s,B))
BA=((t,B),(s,A))

K_ba=word_gram(slot,single,BA)
K_ab=word_gram(slot,single,AB)
P=word_gram(slot,single,single)
GA,GB,C=(source[k] for k in ("G_AB","G_BA","G_cross"))
mp=data["mp"]
lam=data["lY"]
rho=mp.root(source["baseline"]["ed_lo"],3)
tab_upper=source["upper_table"]
necessary=rho/tab_upper
p=1/mp.sqrt(lam)
req=mp.sqrt((necessary-p*p)*(necessary-1)/(necessary-p)**2)
assert req>mp.mpf(".878")
print("Minimum required normalized coherence in the common Gram model:",mp.nstr(req,13))

import importlib
import json
import numpy as np
from scipy.optimize import least_squares

radii=(sp.Rational(4,5),sp.Rational(1,5))
values={t:radii[0],s:radii[1],a:1,g:1}
restrict=lambda M:M.subs(values).applyfunc(sp.cancel)
assert restrict(K_ba)[0,1]==-sp.Rational(176125,2268)
assert restrict(K_ba)[1,0]==-sp.Rational(203125,2268)
for order,K in ((BA,K_ba),(AB,K_ab)):
    inner=((radii[0],B),)
    word=tuple((z.subs(values),terms) for z,terms in order)
    direct=restrict(word_gram(slot,inner,word,oracle=True))
    assert (direct-restrict(K)).applyfunc(sp.cancel)==sp.zeros(4)
for i in source["all_matter"]:
    for order,K in ((BA,K_ba),(AB,K_ab)):
        entry=restrict(word_gram(i,single,order))
        # Cocycle frame signs are retained independently; magnitudes agree.
        assert entry.applyfunc(abs)==restrict(K).applyfunc(abs)

# A true ordered two-point packet within the affine module: exact dilation
# by q = exp(-T) measured relative to the grade-one incoming state.
q=sp.Rational(4,5)
def held(word):
    return tuple((z.subs(values),terms) for z,terms in word)
def mixed_time(bra,ket):
    return restrict(q**len(ket)*word_gram(slot,held(bra),tuple((q*z,k) for z,k in held(ket))))
basis_words=(single,BA)
blocks=[[mixed_time(bra,ket) for ket in basis_words] for bra in basis_words]
timed=blocks[0][0].row_join(blocks[0][1]).col_join(blocks[1][0].row_join(blocks[1][1]))
assert timed.T==timed
assert all(timed[:i,:i].det()>0 for i in range(1,9))
raw_cross=timed[0,5]
assert raw_cross!=0
normalized_cross=raw_cross/sp.sqrt(timed[0,0]*timed[5,5])
print("Exact affine ordered W0 has rank 8; K12 at q=4/5 =", raw_cross,
      "; normalized =",normalized_cross.evalf(12))

# Use the native arithmetic consumer to check which scalar cap sees
# the computed entry, without installing an original source bundle.
runtime=repo/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime"
sys.path.insert(0,str(runtime))
adapter=importlib.import_module("original_source_adapter")
fixture=runtime.parent/"domain/arithmetic-control-bound-v2.json"
bundle=json.loads(fixture.read_text())
bundle["representation"]="original_exact_blocked_endpoint_source"
bundle.pop("fiber_records",None)
bundle["retained_words"]=[["B" if i<4 else "BA",str(i%4+1)] for i in range(8)]
bundle["time"]={"kind":"original_euclidean_ordered","id":"AFFINE_DILATION_CONTROL",
                "T":float(sp.log(1/q)),"t_U":0,"t_Y":0}
def serialize(M):
    return [[[float(sp.re(x)),float(sp.im(x))] for x in M.row(i)] for i in range(M.rows)]
bundle["gram_form"]=bundle["physical_selector_form"]=serialize(sp.eye(8))
bundle["ordered_source_forms"]={k:serialize(timed if k=="Z0" else sp.zeros(8)) for k in adapter.MARKS}
bundle["scalar_factor"]=bundle["reference_scalar"]=[1,0]
domain={"domain_id":bundle["domain_id"],"whole_source_bundle":bundle}
assembled={"domain_id":bundle["domain_id"],"native_inputs_sha256":bundle["native_inputs_sha256"]}
def cap_value(cap):
    bundle["boundary_pairing_form"]=serialize(cap)
    out=adapter.evaluate_source({},assembled,domain)
    assert out["status"]=="evaluated_arithmetic_control_only",out
    return complex(*out["Z0"])
D0=sp.diag(*[sp.I**source["deck"][tuple(2*x for x in p0)] for p0 in source["roots_of_slot"](slot)])
D=D0.row_join(sp.zeros(4)).col_join(sp.zeros(4).row_join(D0))
assert sp.simplify(D[0,0]/D[5,5]) in (sp.I,-sp.I)
avg=sum((D**j*timed*D**(-j) for j in range(4)),sp.zeros(8))/4
assert sp.simplify(avg[0,5])==0
probe=sp.zeros(8)
probe[0,5]=probe[5,0]=sp.Rational(1,2)
unmarked=sp.eye(8)
assert abs(cap_value(unmarked)-float(sp.trace(timed)))<1e-9
assert abs(cap_value(unmarked+probe)-cap_value(unmarked)-float(raw_cross))<1e-9
assert sp.trace(probe*avg)==0
print("The native cap control reads the marked entry by polarization; a deck-invariant cap is blind to it.")

# Quantitative sufficiency tests for both individual mixing orders.
def max_coherence(K,Q,where):
    b1,b2=(P[i,i].subs(where) for i in range(2))
    d1,d2=(Q[i,i].subs(where) for i in range(2))
    return (abs(K[0,1].subs(where))+abs(K[1,0].subs(where)))/(sp.sqrt(b1*d2)+sp.sqrt(b2*d1))
ba_limit=max_coherence(K_ba,GB,values)
ab_limit=max_coherence(K_ab,GA,values)
assert float(ab_limit)<float(req)<float(ba_limit)
print("Sample maximum coherence AB, BA:",float(ab_limit),float(ba_limit))

def bound_joint(where):
    ratios=[]
    for i,j in ((0,1),(1,0)):
        d=sp.Matrix([[GA[j,j],C[j,j]],[C[j,j],GB[j,j]]]).subs(where)
        k=sp.Matrix([K_ab[i,j],K_ba[i,j]]).subs(where)
        assert d.det()>0
        ratios.append(float(sp.sqrt((k.T*d.inv()*k)[0]/P[i,i].subs(where))))
    return max(ratios)
nearby={t:sp.Rational(3,4),s:sp.Rational(1,4),a:1,g:1}
assert bound_joint(nearby)<float(req)
print("At radii (3/4,1/4), arbitrary two-order coherence is bounded above by",bound_joint(nearby))

# A single original AB admixture at the declared radii can be selected
# from the exact lepton ratios. Hold all three individual quark ratios out.
arr=lambda M:np.array(M[:3,:3].subs(values),dtype=float)
gamma,phi=float(lam),float(data["phi0"])
r1=12*phi**2/7
r2=8*phi/7
ratio_target=np.array([55/117,34/(47*phi),3/(26*phi**2)])
ed=float(source["baseline"]["ed"])
color=ed**(1/3)
h=data["h"]
def spectrum(kernel,norm,lengths):
    diagonal=np.diag(gamma**(np.array(lengths)/2))
    eigenvalues,vectors=np.linalg.eigh(norm)
    inverse=(vectors*(1/np.sqrt(eigenvalues)))@vectors.conj().T
    return np.sort(np.linalg.svd(inverse@diagonal@kernel@diagonal@inverse,compute_uv=False))

def profiles(z):
    h1=np.sqrt(float(source["target_ud"])*(16+9*z*z)/16)
    mapping={h[0]:h1,h[1]:1,h[2]:z}
    return tuple(np.array(source["baseline"][key].subs(mapping),dtype=float) for key in ("Eq","Uq"))

def predict(eta,z,K,Q):
    metric=arr(P)+eta*eta*Q+eta*(K+K.T)
    E,U=profiles(z)
    e,d,u=(spectrum(k,metric,l) for k,l in ((E,(8,5,3)),(E,(7,5,2)),(U,(7,3,0))))
    return e[0]/e[1],e[1]/e[2],u/d,e[0]/d[0]*color

for sign in (1,-1):
    K,Q=arr(K_ab),arr(GA)
    def residual(v):
        eta,z=np.exp(v)
        return np.log(np.array(predict(sign*eta,z,K,Q)[:2])/np.array([r1,r2]))
    fit=least_squares(residual,np.log([.15,2.75]),ftol=1e-13,gtol=1e-13,xtol=1e-13)
    assert np.linalg.norm(fit.fun)<1e-9
    eta,z=np.exp(fit.x)
    lepton12,lepton23,quarks,ed_light=predict(sign*eta,z,K,Q)
    assert np.isclose(np.prod(quarks),np.prod(ratio_target),rtol=1e-11)
    assert all(abs(quarks/ratio_target-1)>.3)
    print("AB sign",sign,"lepton fit",lepton12,lepton23,
          "held-out u/d,c/s,t/b:",quarks,"e/u",ed_light/quarks[0])

# Exhaust the BA model's magnitude and complex phase at the same fixed
# radii by using the signed spectrum of the Hermitian generalized pencil.
from mpmath import iv
iv.dps=35
I=iv.mpf
def exact_iv(expr):
    n,dn=sp.Rational(expr).as_numer_denom()
    return I(int(n))/int(dn)
varphi=1/(6*iv.pi)+3/(256*iv.pi**4)
length=iv.sqrt(varphi*(1-varphi))
mue=12*varphi**2/7
mutau=8*varphi/7
aa,bb,cc=2*length**8/3,2*length**5/3,8*length**3/3
base_diag=[exact_iv(restrict(P)[i,i]) for i in range(3)]
word_diag=[exact_iv(restrict(GB)[i,i]) for i in range(3)]
n12=restrict(K_ba)[0,1]
n21=restrict(K_ba)[1,0]
vlow,vhigh=map(exact_iv,((n12-n21)**2,(n12+n21)**2))
unitdet=mue*mutau**2
intervals=0
for positive in range(3):
    mass=[mue*mutau,mutau,I(1)]
    mass=[z if i==positive else -z for i,z in enumerate(mass)]
    unit_trace=sum(mass)
    unit_pair=sum(mass[i]*mass[j] for i,j in ((0,1),(0,2),(1,2)))
    todo=[(0.,1.,0.,1.,0)]
    while todo:
        xl,xr,yl,yr,depth=todo.pop()
        intervals+=1
        x,y=I([xl,xr]),I([yl,yr])
        g1,g2,g3=[c0+(c1-c0)*x for c0,c1 in zip(base_diag,word_diag)]
        off_sq=x*(1-x)*(vlow+(vhigh-vlow)*y)
        delta=g1*g2-off_sq
        if delta.a>0 and g3.a>0:
            scale=((bb*g1-aa*g2)/delta-cc/g3)/unit_trace
            if scale.b<=0:
                continue
            h_sq=scale**3*unitdet*delta*g3/aa-bb*cc
            if h_sq.b<0:
                continue
            pair=(-(bb*cc+h_sq)*g1+aa*cc*g2-aa*bb*g3)/(delta*g3)
            mismatch=pair-scale**2*unit_pair
            if not (mismatch.a<=0<=mismatch.b):
                continue
        if depth>=42:
            raise RuntimeError(("unresolved BA interval",positive,xl,xr,yl,yr))
        if xr-xl>=(yr-yl)/2:
            mid=(xl+xr)/2
            todo.extend(((xl,mid,yl,yr,depth+1),(mid,xr,yl,yr,depth+1)))
        else:
            mid=(yl+yr)/2
            todo.extend(((xl,xr,yl,mid,depth+1),(xl,xr,mid,yr,depth+1)))
print("BA all complex amplitudes excluded by exact-target lepton spectrum interval test:",intervals,"boxes.")

# The word-diagonal information and block reflection positivity admit a
# continuum of cross values. This is a finite packet identifiability check,
# independent of the question of admissible full physical states.
for zeta in (0,sp.Rational(1,2),sp.I):
    cross=zeta*blocks[0][1]
    new=blocks[0][0].row_join(cross).col_join(cross.H.row_join(blocks[1][1]))
    assert all(new[:i,:i].det()>0 for i in range(1,9))
    assert new[0,5]==zeta*raw_cross
print("Fixed diagonal blocks and positivity do not fix the mixed phase or magnitude in the affine packet.")
