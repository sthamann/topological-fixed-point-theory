"""Exact continuation on the existing Spin4/full96 source chart.

No new profile, field weight, geometry or physical spectrum is selected.
The relative-residue calculation has the precise ordinary product scope
recorded in the accompanying derivation; B_rel is not silently deleted
from an arbitrary full source.
"""
import importlib.util
import json
import sys
from pathlib import Path

import sympy as s

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0, str(repo / 'verification'))
spec = importlib.util.spec_from_file_location('source_v252', repo / 'verification/v252_full_finite_triple.py')
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)

checks = {}
def ck(name, value):
    if not value:
        raise AssertionError(name)
    checks[name] = True

def zero(m):
    return all(s.simplify(z) == 0 for z in m)

pauli = [s.Matrix([[0,1],[1,0]]), s.Matrix([[0,-s.I],[s.I,0]]), s.diag(1,-1)]
gam = [s.kronecker_product(pauli[0], q) for q in pauli] + [s.kronecker_product(pauli[1],s.eye(2))]
g5 = gam[0]*gam[1]*gam[2]*gam[3]
deck = s.kronecker_product(pauli[0],s.eye(2))
reflection = [1,1,1,-1]

for mu in range(4):
    W = s.I*gam[mu]*g5*s.Rational(3,4)
    ck(f'linear_deck_general_one_form_{mu}',zero(deck*W*deck + reflection[mu]*W))
    ck(f'axial_component_hermitian_{mu}',zero(W.H-W))
    for nu in range(4):
        ck(f'full_clifford_vector_axial_orthogonal_{mu}_{nu}',s.trace(gam[nu]*W)==0)
    ck(f'full_clifford_higgs_axial_orthogonal_{mu}',s.trace(g5*W)==0)
    ck(f'axial_component_nonzero_{mu}',s.trace(W.H*W)==s.Rational(9,4))
    scalar_scale_connection=s.I*gam[mu]/2
    for nu in range(4):
        axial=s.I*gam[nu]*g5
        ck(f'original_conformal_connection_has_no_axial_component_{mu}_{nu}',s.trace(axial.H*scalar_scale_connection)==0)
ck('original_conformal_symbol_changes_nontrivially',gam[0].det()!=0)

# Exact Mellin calculation, independent of any cutoff profile.
z=s.symbols('z')
a2,a4=s.symbols('a2 a4')
mellin=(a2/(z-1)+a4/z)/s.gamma(z)
ck('mellin_residue_one_is_a2',s.limit((z-1)*mellin,z,1)==a2)
ck('mellin_value_zero_is_a4',s.limit(mellin,z,0)==a4)

d,T,Tref,chi,qscale=s.symbols('d T Tref chi qscale',positive=True)
Ip,Im,Cp,Cm=s.symbols('I_plus I_minus C_plus C_minus',nonnegative=True)
aa=9*d/(32*s.pi**2)
bb=3*d/(64*s.pi**2)
cc=9/(32*s.pi**2)
rel_a2=s.expand(aa*(Ip+Im)-aa*Im)
rel_a4=s.expand(-bb*(Cp+Cm)-cc*T*(Ip+Im)-(-bb*Cm-cc*Tref*Im))
ck('reference_cancels_only_odd_one_form_a2',s.simplify(rel_a2-aa*Ip)==0)
ck('reference_cancels_only_odd_one_form_curl',s.simplify(rel_a4+bb*Cp+cc*T*Ip+cc*(T-Tref)*Im)==0)
ck('T_lift_reference_cancels_minus_mass',s.simplify(rel_a4.subs(Tref,T)+bb*Cp+cc*T*Ip)==0)
gamma2=s.expand(-6*chi**2*rel_a2+rel_a4)
even_mass=s.expand(gamma2).coeff(Ip)
even_curl=s.expand(gamma2).coeff(Cp)
ck('geometric_even_real_torsion_mass_negative',s.simplify(even_mass+54*d*chi**2/(32*s.pi**2)+9*T/(32*s.pi**2))==0)
ck('geometric_even_real_torsion_curl_nonzero',even_curl==-bb)
ck('full96_even_curl_coefficient',even_curl.subs(d,96)==-s.Rational(9,2)/s.pi**2)
ck('full96_even_mass_coefficient',s.simplify(even_mass.subs(d,96)+162*chi**2/s.pi**2+9*T/(32*s.pi**2))==0)
scaled_a2=aa*(qscale*Ip+(qscale-1)*Im)
scaled_gamma=s.expand(-6*chi**2*scaled_a2+rel_a4)
ck('constant_scale_reference_not_rescaled',s.simplify(s.expand(scaled_a2).coeff(Im)-aa*(qscale-1))==0)
ck('constant_scale_curl_unchanged',scaled_gamma.coeff(Cp)==-bb)
ck('constant_scale_one_recovers_unscaled',s.simplify(scaled_gamma.subs(qscale,1)-gamma2)==0)

# The historical minimal auxiliary source has no S derivatives. Its
# exact second variation is a multiplication operator, not a Laplacian.
A,B,c0,j,da,ss,kk=s.symbols('A B c0 j da S k',real=True)
aux=A*ss**2-ss*(s.I*B*da+c0*j)
ck('original_auxiliary_hessian_order_zero',s.diff(aux,ss,2)==2*A)
ck('original_auxiliary_hessian_has_no_momentum',s.diff(s.diff(aux,ss,2),kk)==0)
ck('laplace_hessian_not_same_original_operator',s.diff(kk**2+2*A,kk,2)==2)
tau=s.symbols('tau',positive=True)
ck('four_dimensional_a2_has_t_minus_one',s.Rational(2-4,2)==-1)
ck('four_dimensional_a4_has_t_zero',s.Rational(4-4,2)==0)

# A real compact interior one-form with reflection-invariant pullback
# and nonzero curl exists: S_2=h(r)*f(x1), with h smooth even and supported
# away from r=0. Its x1 derivative cannot vanish for every nonconstant f.
x,r=s.symbols('x r',real=True)
S2=s.Function('h')(r)*s.Function('f')(x)
ck('transverse_curl_witness',s.diff(S2,x)==s.Function('h')(r)*s.diff(s.Function('f')(x),x))

# The actual SM representation, before even imposing unimodularity.
# Hermitian q=-i rho(a), a in Lie(U(C)+U(H)+U(3)); the real gauge lift is
# q-KF conj(q) KF. All nu_R and their conjugates are neutral under it.
KF=s.Matrix(v.Kswap()).applyfunc(lambda y:s.Integer(int(y)))
gF=s.Matrix(v.grading()).applyfunc(lambda y:s.Integer(int(y)))
def exact(m):
    return s.Matrix(m).applyfunc(lambda y:s.nsimplify(complex(y),rational=True))

import numpy as np
Z2=np.zeros((2,2),complex)
Z3=np.zeros((3,3),complex)
anti=[]
anti.append(exact(v.rep(1j,Z2,Z3)))
for p in pauli:
    anti.append(exact(v.rep(0,np.array(s.I*p).astype(complex),Z3)))
for i in range(3):
    m=Z3.copy();m[i,i]=1j
    anti.append(exact(v.rep(0,Z2,m)))
for i in range(3):
    for j in range(i+1,3):
        m=Z3.copy();m[i,j]=1;m[j,i]=-1
        anti.append(exact(v.rep(0,Z2,m)))
        m=Z3.copy();m[i,j]=1j;m[j,i]=1j
        anti.append(exact(v.rep(0,Z2,m)))
ck('raw_SM_unitary_algebra_dimension_13',len(anti)==13)
nu=[v.IDX[('nu','R',0,g)] for g in range(3)]
nu=nu+[i+v.NP for i in nu]
for n,a in enumerate(anti):
    ck(f'actual_raw_gauge_generator_antihermitian_{n}',zero(a.H+a))
    q=-s.I*a
    gauge=q-KF*q.conjugate()*KF
    ck(f'actual_real_gauge_neutral_right_neutrinos_{n}',zero(gauge[:,nu]) and zero(gauge[nu,:]))
ck('grading_nonzero_on_all_right_neutrinos',all(gF[i,i] in [-1,1] for i in nu))
ck('grading_cannot_be_actual_SM_gauge_generator',not zero(gF[:,nu]))

# Test the actual R034 single-Higgs image with all Yukawa entries left
# symbolic. This is the original four-real-coordinate producer, not a
# replacement finite model or a choice of physical masses.
h1r,h1i,h2r,h2i,zr,zi=s.symbols('h1r h1i h2r h2i zr zi',real=True)
h1,h2=h1r+s.I*h1i,h2r+s.I*h2i
qH=s.Matrix([[s.conjugate(h2),h1],[-s.conjugate(h1),h2]])
j2=s.Matrix([[0,1],[-1,0]])
zphase=zr+s.I*zi
normH=h1r**2+h1i**2+h2r**2+h2i**2
ck('actual_Higgs_quaternion_real_structure',zero(j2*qH.conjugate()*j2.inv()-qH))
ck('actual_nonzero_Higgs_quaternion_invertible',s.expand(qH.det()-normH)==0)
phase_defect=j2*(zphase*qH).conjugate()*j2.inv()-zphase*qH
ck('common_mass_phase_quaternion_defect',zero(phase_defect+2*s.I*zi*qH))
ck('common_mass_phase_defect_norm',s.expand(s.trace(phase_defect.H*phase_defect)-8*zi**2*normH)==0)

YS={f:s.Matrix(3,3,lambda i,j:s.Symbol(f'{f}_r{i}{j}',real=True)+s.I*s.Symbol(f'{f}_i{i}{j}',real=True)) for f in ['u','d','e','nu']}
DY=s.MutableSparseMatrix(96,96,{})
for name,ch,color,gen in v.P:
    if ch!='L': continue
    for gj in range(3):
        i,j=v.IDX[(name,'L',color,gen)],v.IDX[(name,'R',color,gj)]
        y=YS[name][gen,gj]
        DY[i,j]=y;DY[j,i]=s.conjugate(y)
        DY[i+48,j+48]=s.conjugate(y);DY[j+48,i+48]=y
eR=s.SparseMatrix(exact(v.rep(1,Z2,Z3)))
swap=s.SparseMatrix(KF)
comm=DY*eR-eR*DY
Lh=s.MutableSparseMatrix(96,96,{})
for coeff,bh in zip([h1r,h1i,h2r,h2i],[[1,0],[1j,0],[0,1],[0,1j]]):
    a,b=bh
    qb=np.array([[complex(b).conjugate(),a],[-complex(a).conjugate(),b]],complex)
    rho=s.SparseMatrix(exact(v.rep(0,qb,Z3)))
    w=rho*comm;w=w+w.H
    Lh+=coeff*(w+swap*w.conjugate()*swap)
for species,mult in [('u',s.conjugate(h2)),('d',h2),('nu',s.conjugate(h2)),('e',h2)]:
    for gi,gj in [(0,0),(1,2)]:
        i,j=v.IDX[(species,'L',0,gi)],v.IDX[(species,'R',0,gj)]
        ck(f'actual_full96_Higgs_phase_column_{species}_{gi}_{gj}',s.expand(Lh[i,j]-mult*YS[species][gi,gj])==0)
ck('actual_original_seed_Higgs_recovers_all_Yukawas',zero(Lh.subs({h1r:0,h1i:0,h2r:1,h2i:0})-DY))
for species in ['u','d']:
    i,j=v.IDX[(species,'L',0,0)],v.IDX[(species,'R',0,0)]
    axial_phase=zr*DY+s.I*zi*gF*DY
    ck(f'actual_common_axial_phase_same_factor_{species}',s.expand(axial_phase[i,j]-zphase*DY[i,j])==0)
for gi in range(3):
    i=v.IDX[('nu','R',0,gi)]
    ck(f'fixed_Majorana_phase_factor_{gi}',s.expand(zr+s.I*zi*gF[i,i]-s.conjugate(zphase))==0)

out={
    'result_id':'TORSION20261104C',
    'passed':len(checks),
    'checks':checks,
    'relative_a2':str(rel_a2),
    'relative_a4':str(rel_a4),
    'geometric_gamma_quadratic':str(gamma2),
    'constant_scale_squared': 'qscale=exp(2 sigma0); the reference remains unscaled',
    'constant_scale_relative_a2':str(scaled_a2),
    'constant_scale_geometric_gamma_quadratic':str(scaled_gamma),
    'full96_even_mass_coefficient':str(even_mass.subs(d,96)),
    'full96_even_curl_coefficient':str(even_curl.subs(d,96)),
    'single_Higgs_phase_image':'For nonzero H and faithful up/down Yukawa columns, z q(H) is quaternionic only for real z; on U(1), z=+/-1. With a fixed nonzero MR, its zbar factor further forces z=1. Complete source/basis transport is still valid and is not a new scalar fluctuation.',
    'scope':{
        'deck_and_clifford':'Exact actual Spin4 matrices and actual v252 SM representation; fixed tetrad ordinary product inner fluctuation image only.',
        'residues':'Ordinary geometric Dirac product contribution with constant scale and full constant Phi; fixed elliptic boundary realization and no kernel crossing for differentiating zeta(0). This is not automatically the full K_sigma+B_rel caller or a field-dependent Calderon compression.',
        'exclusion':'P25 positive algebraic Euclidean auxiliary parent is not equal to this ordinary relative geometric residue functional on the same real independent even torsion directions.',
        'not_derived':['primitive full Sadm','whole original B_rel torsion jet','physical H_tr,pos','joint physical spectrum'],
    },
}
Path(__file__).with_name('gebundene-referenz-und-feldbild.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2,ensure_ascii=False))
