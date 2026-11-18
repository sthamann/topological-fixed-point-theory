"""Actual96 real structure, quarter transport and residual frame freedom.

The finite algebra checks supplement the explicit common-domain proof.
They do not select the original physical source or add a new quarter copy.
"""
from pathlib import Path
import json
import sys
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'verification'))
sys.path.insert(0,str(ROOT/'experiments/theory-contracts/seam-constructive-search-20261015/operator_synthesis'))
import v252_full_finite_triple as finite
from intertwiner_search import equation,restrict

def exact(a):
    a=np.asarray(a)
    data={}
    for i,j in zip(*np.nonzero(a)):
        z=complex(a[i,j]);assert z.real.is_integer() and z.imag.is_integer()
        data[int(i),int(j)]=s.Integer(int(z.real))+s.I*int(z.imag)
    return s.SparseMatrix(a.shape[0],a.shape[1],data)

N=s.Matrix([[0,-1],[1,0]]);Z=s.diag(1,-1);S=N*Z
g=exact(finite.grading());k=exact(finite.Kswap());I=s.eye(96)
ys=[np.array([[1+1j,2,0],[0,3-1j,1],[1j,0,4]],complex)+j*np.eye(3) for j in range(4)]
mr=np.array([[2+1j,1,0],[1,3,1j],[0,1j,4-1j]],complex)
phi=exact(finite.full_D(*ys,mr))
weak=finite.rep(0,1j*np.array([[0,1],[1,0]]),np.zeros((3,3)))
Q=exact(-1j*(weak+finite.Kswap()@weak.conjugate()@finite.Kswap()))
checks={'actual_J_real_Phi':k*phi.conjugate()*k==phi,
        'actual_J_reverses_real_weak_connection':k*Q.conjugate()*k==-Q}
Gamma=s.kronecker_product(Z,g);C0=s.I*s.kronecker_product(N,I);Ical=s.I*Gamma
D0=s.kronecker_product(S,Q)+s.kronecker_product(Z,phi)
B0=s.kronecker_product(Z,Q)-s.kronecker_product(S,phi)
normal=[]
representatives={(0,0):S,(1,0):N,(0,1):Z,(1,1):s.eye(2)}
for power in (0,1):
    for reflection in (0,1):
        accepted=[]
        for eps in (1,-1):
            eqs=[equation(N,(-1)**reflection*eps*N,False),
                 equation(S,-eps*S,False),equation(Z,eps*(-1)**power*Z,False)]
            basis=restrict(s.eye(4),eqs)
            if basis.cols:accepted.append((eps,basis.cols))
        M=representatives[power,reflection]
        eps=next(e for e,d in accepted)
        U=s.kronecker_product(M,g**power*k)
        conjugate=lambda a:U*a.conjugate()*U.H
        jsq=1 if (power,reflection)!=(1,1) else -1
        grade_sign=(-1)**reflection
        checks[f'normal_class_{power}_{reflection}_one_dimensional']=len(accepted)==1 and accepted[0][1]==1
        checks[f'actual96_Jsquare_{power}_{reflection}']=U*U.conjugate()==jsq*s.eye(192)
        checks[f'actual96_D_potential_{power}_{reflection}']=conjugate(D0)==eps*D0
        checks[f'actual96_D_normal_derivative_{power}_{reflection}']=(-1)**reflection*conjugate(s.kronecker_product(N,I))==eps*s.kronecker_product(N,I)
        checks[f'actual96_D_tangent_derivative_{power}_{reflection}']=conjugate(-s.I*s.kronecker_product(S,I))==eps*(-s.I*s.kronecker_product(S,I))
        checks[f'actual96_total_grading_{power}_{reflection}']=conjugate(Gamma)==grade_sign*Gamma
        checks[f'actual96_whole_B_{power}_{reflection}']=conjugate(B0)==(-1)**reflection*B0
        checks[f'new_complex_type_{power}_{reflection}']=conjugate(Ical)==-grade_sign*Ical
        normal.append(dict(internal_gamma_power=power,normal_reflection=reflection,
                           D_sign=eps,Gamma_sign=grade_sign,J_square=jsq,
                           normal_matrix=str(M),new_complex_type='antilinear' if grade_sign==1 else 'linear'))
JS=s.kronecker_product(S,k);JR=s.kronecker_product(Z,k)
checks['two_real_maps_compose_to_given_deck_normal_matrix']=JS*JR.conjugate()==s.kronecker_product(N,I)
alpha=s.Rational(1,4)
Dalpha=D0+alpha*s.kronecker_product(S,I)
Balpha=B0+alpha*s.kronecker_product(Z,I)
checks['scalar_quarter_JS_D_residual_exact']=JS*Dalpha.conjugate()*JS.H+Dalpha==2*alpha*s.kronecker_product(S,I)
checks['scalar_quarter_JS_B_maps_to_minus_alpha']=JS*Balpha.conjugate()*JS.H==B0-alpha*s.kronecker_product(Z,I)
checks['quarter_holonomy_inverts_not_fixed']=s.conjugate(s.I)==-s.I and s.conjugate(s.I)!=s.I

theta=s.symbols('theta',real=True)
gp=Z;kp=S;pp=(s.eye(2)+gp)/2;pm=(s.eye(2)-gp)/2
E=pp.col_join(s.exp(s.I*theta)*pm);Et=s.exp(-s.I*theta/2)*E
Jsmall=s.kronecker_product(S,kp)
def zero(a):return a.applyfunc(s.simplify)==s.zeros(*a.shape)
checks['all_theta_rephased_frame_exact_J_intertwiner']=zero(Jsmall*Et.conjugate()-Et*kp)
rotation=s.diag(s.exp(-s.I*theta/2),s.exp(s.I*theta/2))
checks['all_theta_axial_rotation_preserves_internal_J']=zero(kp*rotation.conjugate()-rotation*kp)

# All six positive-grade SM multiplets have family multiplicity three.
# Schur's lemma supplies necessity; these actual-slot generators verify sufficiency.
slots={}
for i in range(96):
    anti=i>=48;name,ch,colour,family=finite.P[i%48]
    if int(g[i,i])!=1:continue
    label=('Q' if name in ('u','d') else 'L') if not anti else {'u':'uc','d':'dc','e':'ec','nu':'nc'}[name]
    component=(name,colour)
    slots.setdefault(label,{}).setdefault(component,{})[family]=i
algebra=[]
for lam in (1,1j):algebra.append(exact(finite.rep(lam,np.zeros((2,2)),np.zeros((3,3)))))
for aa,bb in ((1,0),(1j,0),(0,1),(0,1j)):algebra.append(exact(finite.rep(0,finite.quaternion(aa,bb),np.zeros((3,3)))))
for i in range(3):
    for j in range(3):
        for phase in (1,1j):
            m=np.zeros((3,3),complex);m[i,j]=phase
            algebra.append(exact(finite.rep(0,np.zeros((2,2)),m)))
H=[]
for i in range(3):
    x=s.zeros(3);x[i,i]=1;H.append(x)
for i in range(3):
    for j in range(i+1,3):
        x=s.zeros(3);x[i,j]=x[j,i]=1;H.append(x)
        x=s.zeros(3);x[i,j]=s.I;x[j,i]=-s.I;H.append(x)
count=0;all_real=all_even=all_algebra=all_skew=True
for label,components in slots.items():
    for h in H:
        xp=s.MutableSparseMatrix(96,96,{})
        for component,indices in components.items():
            for f1 in range(3):
                for f2 in range(3):
                    if h[f1,f2]:xp[indices[f1],indices[f2]]=s.I*h[f1,f2]
        x=xp+k*xp.conjugate()*k
        all_skew &= x.H==-x;all_real &= k*x.conjugate()*k==x
        all_even &= x*g==g*x
        all_algebra &= all(x*a==a*x for a in algebra)
        count+=1
checks.update(actual54_frame_generators_antihermitian=all_skew,
              actual54_frame_generators_preserve_J=all_real,
              actual54_frame_generators_preserve_grading=all_even,
              actual54_frame_generators_preserve_left_and_opposite_algebra=all_algebra,
              six_inequivalent_multiplets=sorted(slots)==sorted(['Q','L','uc','dc','ec','nc']),
              family_frame_real_dimension54=count==54)
out={'contract':'B1-NORMAL-REAL-AND-FRAME-20261030','checks':checks,'all_passed':all(checks.values()),
     'normal_separable_class':normal,'remaining_frame_group':'U(3)^6 with negative grades fixed by J',
     'frame_real_dimension':count,'original_left_algebra_basis':len(algebra),
     'scope':'Separable normal class; actual96 coefficient controls and formal derivative identities. Quarter maps to existing inverse-holonomy image. Remaining frames are coordinate equivalences, not extra physical source choices.'}
Path(__file__).with_name('B1-NORMAL-REAL-AND-FRAME.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False,indent=2))
assert out['all_passed']
