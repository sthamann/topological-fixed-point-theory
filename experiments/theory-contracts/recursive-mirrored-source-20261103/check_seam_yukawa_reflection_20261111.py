"""Exact tests of variable marked geometry, source inversion, and charge dressing."""
import contextlib, io, itertools, runpy, sys
from pathlib import Path
root=Path(sys.argv[1])
prior=root/"experiments/theory-contracts/recursive-mirrored-source-20261103/check_seam_ordered_yukawa_20261111.py"
with contextlib.redirect_stdout(io.StringIO()):
    a=runpy.run_path(str(prior))
sp=a["sp"]
alpha,beta=a["alpha"],a["beta"]
gamma=a["complement"][0]
h=a["h"]
T,S,d=sp.symbols("T S d",positive=True)
def coeff(roots,positions):
    if any(sum(column)!=0 for column in zip(*roots)):
        return 0,0
    i,j=len(roots)-3,len(roots)-2
    order=a["charge"](roots[i],roots[j])
    assert order>=-1
    v=sp.S.One
    rate=sp.S.Zero
    for k,l in itertools.combinations(range(len(roots)),2):
        v*=a["epsilon"](roots[k],roots[l])
        if (k,l)==(i,j):continue
        dist=positions[k]-positions[l]
        p=a["charge"](roots[k],roots[l])
        v*=dist**p
        if i in (k,l):
            rate+=p*(1 if k==i else -1)/dist
    return {-1:(v,v*rate),0:(0,v)}.get(order,(0,0))

def matrix(slot, first, second, conjugate=False, coordinates=(T,S,sp.S.Zero,sp.S.Zero,-d)):
    roots=a["matter_vertex"]
    B,P=sp.zeros(3),sp.zeros(3)
    for i,j,k in itertools.product(range(3),repeat=3):
        p,ps=roots(slot[0],i)
        q,qs=roots(slot[1],j)
        r,rs=a["conjugated_down_higgs"](k) if conjugate else roots(slot[2],k)
        for (v,vs),(w,ws) in itertools.product(first,second):
            p1,f1=coeff((v,w,p,q,r),coordinates)
            p2,f2=coeff((v,w,q,p,r),coordinates)
            assert sp.cancel(p1+p2)==0
            factor=ps*qs*rs*vs*ws*h[k]/(2*a["raw_coefficient"](slot))
            B[i,j]+=factor*(f1+f2)
            P[i,j]+=factor*(p1-p2)
    return B.applyfunc(sp.factor), P.applyfunc(sp.factor)

def rooted(p):
    return a["hermitian_roots"](p)
e=a["slots"]["e"][0]
u=a["slots"]["u"][0]
E,Ep=matrix(e,rooted(alpha),rooted(beta))
U,Up=matrix(u,rooted(gamma),rooted(beta),True)

# Write the determinants as a sum of three signed monomials before simplifying.
def diagonal_terms(block):
    assert block[0,2] == block[2,0] == 0 and block.T == block
    return tuple(sp.factor(term) for term in (
        -block[2,2]*block[0,1]**2,
        block[0,0]*block[1,1]*block[2,2],
        -block[0,0]*block[1,2]**2,
    ))

D_e = diagonal_terms(E)
D_u = diagonal_terms(U)
assert D_u[1:] == (0,0)
x = d*(T-S)/(S*(T+d))
# The first positive contribution to det(E) bounds the absolute determinant ratio.
# Parametrizing T=S+w makes the signs of the three contributions manifest.
w = sp.symbols("w",positive=True)
h1,h2,h3 = h
assert sp.factor(-D_u[0]*(1-x)**2-D_e[0]*x**2) == 0
for term, monomial in zip(D_e,(h1**2*h2,h2**3,h2*h3**2)):
    value = sp.factor((term/monomial).subs(T,S+w))
    assert value.is_nonnegative is True or sp.factor(value / (S*(S+w)+d*S-d*w)**2).is_positive is True
assert sp.factor((D_e[1]/h2**3).subs(T,S+w)).is_positive is True
assert sp.simplify(Ep-Up) == sp.zeros(3)
old = {T:3,S:2,d:1}
assert x.subs(old) == sp.Rational(1,8)
assert (x/(1-x)).subs(old)**2 == sp.Rational(1,49)
quarter = {T:2,S:sp.Rational(1,2),d:1}
points = (sp.S.One,sp.I,-sp.S.One,-sp.I)
z_h,z_m,z_s,z_t = points
assert sp.simplify((z_m-z_h)*(z_t-z_s)/((z_s-z_m)*(z_t-z_h))) == 1
assert x.subs(quarter) == 1

Eq,Uq,Pq = (X.subs(quarter) for X in (E,U,Ep))
quarter_ratio = 16*h1**2/(16*h2**2+9*h3**2)
assert sp.factor(Uq.det()+quarter_ratio*Eq.det()) == 0
assert quarter_ratio.subs(dict.fromkeys(h,1)) == sp.Rational(16,25)
for color in a["slots"]["u"]:
    other,pole = matrix(color,rooted(gamma),rooted(beta),True, (2,sp.Rational(1,2),0,0,-1))
    assert other == Uq and pole == Pq
for color in a["slots"]["d"]:
    other,pole = matrix(color,rooted(alpha),rooted(beta),False,(2,sp.Rational(1,2),0,0,-1))
    assert other == Eq and pole == Pq

# Independent nearby rational evaluation at the quarter cross ratio.
slot=e
p,ps = a["matter_vertex"](slot[0],1)
q,qs = a["matter_vertex"](slot[1],2)
hh,hs = a["matter_vertex"](slot[2],2)
near = sp.Rational(1,10**8)
raw = 0
for (r,rs),(s,ss) in itertools.product(rooted(alpha),rooted(beta)):
    pos = (2,sp.Rational(1,2),near,0,-1)
    raw += rs*ss*(a["vertex_response"]((r,s,p,q,hh),pos)+a["vertex_response"]((r,s,q,p,hh),pos))
raw *= ps*qs*hs/(2*a["raw_coefficient"](slot))
assert abs(raw-Eq[1,2].subs({h1:0,h2:0,h3:1})) < sp.Rational(1,10**6)

mp = a["mp"]
target_ud = a["reduced_ud"]
threshold = mp.sqrt(target_ud)/(1+mp.sqrt(target_ud))
print("All ordered positions: determinant bound (x/(1-x))^2; old x=1/8, threshold",mp.nstr(threshold,13))
print("Quarter marks: x=1; exact |Du/Dd|=16*h1^2/(16*h2^2+9*h3^2), unit h gives 16/25.")
print("Quarter rank three survives equal weights, and h1^2=25*R_ud/16 for h2=h3=1 reaches the exact determinant.")

# Read the actual inversion signs in their own pinned Chevalley convention.
import json
parent = json.loads((root/"experiments/theory-contracts/compiler-e6-seam-readout-20260918/certificate.json").read_text())
inv = parent["d4"]
sign = dict(zip(map(tuple,inv["inversion_root_order"]),map(int,inv["inversion_root_signs"]),strict=True))
phase = dict(zip(map(tuple,parent["seam"]["root_order"]),parent["seam"]["root_character_exponents"],strict=True))
def reflected(r):
    return r[:6]+(r[7],r[6])
def doubled(r):
    return tuple(2*x for x in r)
def neg(r):
    return tuple(-x for x in r)
assert reflected(beta) == neg(beta)
assert sign[doubled(beta)] == sign[doubled(neg(beta))] == -1
# Compact adjoint in the pinned Chevalley convention: e_r*=-e_-r.
assert a["old_spin"](beta) == a["old_spin"](reflected(beta)) == -1
for source,expected_sign in ((alpha,1),(gamma,-1)):
    partner=reflected(source)
    assert sign[doubled(source)] == sign[doubled(partner)] == expected_sign
    assert phase[doubled(source)] == 1 and phase[doubled(partner)] == 3
    assert a["old_spin"](source) == 1 and a["old_spin"](partner) == -1
    spin=sp.diag(1,-1)
    inversion=sp.Matrix([[0,expected_sign],[expected_sign,0]])
    assert spin*inversion == -inversion*spin
    projector=(sp.eye(2)+spin)/2
    inv_projector=(sp.eye(2)+inversion)/2
    assert projector*inv_projector*projector == projector/2
    assert sp.Matrix.vstack(spin-sp.eye(2),inversion-sp.eye(2)).nullspace() == []
# An apparently natural Weyl interchange of the two unit-normalized currents
# preserves the color/weak carrier but sends one marked spinor to its complement.
weyl=lambda r:r[:6]+(-r[6],-r[7])
assert weyl(alpha) == neg(gamma) and weyl(beta) == neg(beta)
family = [tuple(sp.Rational(1,2)-int(i==j) for i in range(3)) for j in range(3)]
family_complement = (sp.Rational(-1,2),)*3
assert weyl((0,)*5+family[0])[5:] == family_complement
print("Actual Chevalley inversion maps beta to its compact adjoint. It sends alpha and gamma to distinct odd-spin partners.")
print("The fixed even-spin/inversion-even intersection is zero on both first-current pairs. No mixing of alpha and gamma by the deck/reflection.")
print("The Weyl interchange of alpha and gamma changes the marked triplet projector.")

# Smooth Abelian Ward dressing at matter coincidence.
# The symmetric regular part changes by z times the common skew pole.
z=sp.symbols("z")
A=144*h2**2+81*h3**2
B=256*h1**2+16*h3**2
assert sp.factor((Eq+z*Pq).det()-2*h2*(A-B*z**2)/243)==0
assert sp.factor((Uq+z*Pq).det()-32*h1**2*h2*(4*z**2-9)/243)==0
# For a real continuation on the connected positive determinant branch,
# the more charged lepton always has a smaller determinant than down.
q2=sp.symbols("q2", positive=True)
assert sp.factor((A-B*q2)-(A-B*q2/9)) == -8*B*q2/9 or sp.factor(
    (A-B*q2)-(A-B*q2/9)+8*B*q2/9)==0

# Pure compact phase is a different hypothesis; test with the old source table
# including independent half-unit intervals on all printed coefficients.
TAB=a["TAB"]
def residual_ed(offset_lep=0,offset_down=0):
    return mp.fprod(mp.mpf(TAB[f][1])+offset_lep for f in ("e","mu","tau"))/mp.fprod(
        mp.mpf(TAB[f][1])+offset_down for f in ("d","s","b"))
rounding=a["half_printing"]
ed_lo,ed,ed_hi = residual_ed(-rounding,rounding),residual_ed(),residual_ed(rounding,-rounding)
assert ed_lo>1
def phase_weight(ed_ratio, h1_squared, h2_squared, h3_squared):
    magnitude_ratio=9*(ed_ratio-1)/(9-ed_ratio)
    quadratic=144*h2_squared+81*h3_squared
    phase_squared=magnitude_ratio*quadratic/(256*h1_squared+16*h3_squared)
    residual_ud=144*h1_squared/quadratic*(1+16*phase_squared/81)/(1+magnitude_ratio/9)
    return phase_squared,residual_ud
bounds=[phase_weight(e,1,1,1)[1] for e in (ed_lo,ed_hi)]
assert max(bounds)<target_ud
print("Real smooth Ward branch connected without determinant zeros predicts |De/Dd|<1; old table interval", tuple(mp.nstr(v,12) for v in (ed_lo,ed_hi)))
print("Compact phase with equal Higgs components gives |Du/Dd|",tuple(mp.nstr(v,12) for v in bounds),"when |De/Dd| is matched.")

# Construct two distinct profile witnesses within the compact conditional
# class. They meet both determinant tests but have different light spectra.
import numpy as np
ratios=[]
for third_squared in (0,1):
    first_squared=mp.findroot(lambda y: phase_weight(ed,y,1,third_squared)[1]-target_ud,(mp.mpf(".6"),mp.mpf("1.4")))
    phase_squared,_ = phase_weight(ed,first_squared,1,third_squared)
    vals={h1:float(mp.sqrt(first_squared)),h2:1,h3:float(mp.sqrt(third_squared))}
    pole=np.asarray(Pq.subs(vals),dtype=float)
    E_num=np.asarray(Eq.subs(vals),dtype=float)+1j*float(mp.sqrt(phase_squared))*pole
    U_num=np.asarray(Uq.subs(vals),dtype=float)+1j*float(mp.sqrt(phase_squared))*pole*2/3
    D_num=np.asarray(Eq.subs(vals),dtype=float)+1j*float(mp.sqrt(phase_squared))*pole/3
    assert abs(abs(np.linalg.det(E_num)/np.linalg.det(D_num))-float(ed))<1e-12
    assert abs(abs(np.linalg.det(U_num)/np.linalg.det(D_num))-float(target_ud))<1e-12
    light=np.linalg.svd(E_num,compute_uv=False)[-1]/np.linalg.svd(U_num,compute_uv=False)[-1]
    ratios.append(light)
    print("Compact witness h3^2 =",third_squared,"h1^2 =",mp.nstr(first_squared,12),
          "phase^2 =",mp.nstr(phase_squared,12),"light residual ratio =",light)
assert abs(ratios[0]-ratios[1])>.1

# Equal radial-time increments give geometrically spaced radii.
rho=sp.symbols("rho",positive=True)
geometric=x.subs({d:rho-1,S:rho**2-rho,T:rho**3-rho})
assert sp.cancel(geometric-rho/(rho**2+rho+1)) == 0
assert sp.cancel(sp.Rational(1,3)-rho/(rho**2+rho+1)-(rho-1)**2/(3*(rho**2+rho+1))) == 0
assert target_ud>mp.mpf(1)/4
print("Equal radial time increments force x<=1/3, hence |Du/Dd|<1/4 at equal weights.")

# If reflection is assumed to split each original lepton word length
# equally over the two endpoint factors, the compact dressed matrix
# is Hermitian. Its trace fixes the scale of any proposed spectrum.
phi, lam = a["phi0"],a["lY"]
diag_first=-mp.mpf(2)*lam**8/3
diag_middle=mp.mpf(2)*lam**5/3
diag_last=-mp.mpf(8)*lam**3/3
ratio_e_mu=mp.mpf(12)*phi**2/7
ratio_mu_tau=mp.mpf(8)*phi/7
trace=diag_first+diag_middle+diag_last
assert trace<0 and ratio_mu_tau*(1+ratio_e_mu)<1
diagonal_det=diag_first*diag_middle*diag_last
sum_diag_pairs=diag_first*diag_middle+diag_first*diag_last+diag_middle*diag_last
max_normalized_charge=9*(ed_hi-1)/(9-ed_hi)
for eigen_signs in ((1,-1,-1),(-1,1,-1)):
    eigenvalues=[s*v for s,v in zip(eigen_signs,(ratio_e_mu*ratio_mu_tau,ratio_mu_tau,1))]
    eigenvalues=[e*trace/sum(eigenvalues) for e in eigenvalues]
    pair_sum=sum(eigenvalues[i]*eigenvalues[j] for i,j in ((0,1),(0,2),(1,2)))
    off_diag_sum=sum_diag_pairs-pair_sum
    required_first=(mp.fprod(eigenvalues)-diagonal_det+diag_first*off_diag_sum)/(diag_first-diag_last)
    required_second=off_diag_sum-required_first
    if eigen_signs[0]==1:
        assert required_second < -mp.mpf("5e-5")
        print("Exact lepton ratios with positive electron eigenvalue require negative |Y23|^2:",
              mp.nstr(required_second,12))
    else:
        upper_first=lam**13*max_normalized_charge*(144+81*required_second/lam**8)/324
        assert required_second>0 and required_first>13*upper_first>0
        print("Exact lepton ratios with positive muon eigenvalue require |Y12|^2",
              mp.nstr(required_first,12),"vs upper bound",mp.nstr(upper_first,12))
        required_charge=324*required_first/(lam**13*(144+81*required_second/lam**8))
        minimum_ed=(1+required_charge)/(1+required_charge/9)
        assert minimum_ed>mp.mpf("1.9")
        print("The equal-leg compact branch would need |De/Dd| at least",mp.nstr(minimum_ed,12))
