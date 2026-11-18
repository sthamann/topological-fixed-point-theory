"""Exact controls for the old-rho, parity-conditioned odd-sector re-fielding.

This is a source/time-compatible representation candidate, not a P1 local
field construction. Its locality failure is part of the controlled result.
"""
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import json
import sympy as S


def need(value, label):
    if not value:
        raise ValueError(label)


def theta_charge(alpha, half, cutoff=128):
    state={(0,0,F(0)):1}
    for k in range(8):
        delta=(F(k)+alpha)/8-F(1,2)
        base=F(1,2) if delta<0 else F(-1,2)
        baseline=16*base*base+32*delta*base if half else 0
        options=[]
        for m in range(-4,5):
            q=F(m)+F(1,2) if half else F(m)
            e=16*q*q+32*delta*q-baseline
            if 0<=e<=cutoff:
                need(e.denominator==1,'integral quarter-energy exponent')
                options.append((int(e),m%2,q))
        nxt=defaultdict(int)
        for (energy,parity,charge),count in state.items():
            for e,p,q in options:
                if energy+e<=cutoff:
                    nxt[energy+e,(parity+p)%2,charge+q]+=count
        state=dict(nxt)
    return state


P=S.zeros(8)
for k in range(8):P[k,(k+4)%8]=1
r=S.Matrix([S.Rational(1,2)]*4+[-S.Rational(1,2)]*4)
R=S.zeros(8)
for k in range(8):R[k,(k-1)%8]=1
A=P*R
unit0=S.eye(8)[:,0]
lam=P*unit0+r
need(R*P==P*R and R*r+unit0==lam,'same primitive winding map on both parities')
need(lam.dot(lam)==2 and sum(lam)==1,'primitive winding spinor root and unit total charge')
need(A**8==S.eye(8),'primitive winding permutation order8')
need(sum((A**j*lam for j in range(8)),S.zeros(8,1))==S.ones(8,1),
     'eight primitive steps retain integer charge carry')
checks=[]
for alpha in [S.Rational(1,4),S.Rational(3,4)]:
    delta=S.Matrix([(k+alpha)/8-S.Rational(1,2) for k in range(8)])
    need(P*P==S.eye(8) and P*r==-r,'rho involution')
    need(P*delta==delta+r,'full energy linear coefficient')
    need(r.dot(r)/2+delta.dot(r)==0,'full energy constant')
    need(sum(r)==0,'total charge preserved')
    symbols=S.Matrix(S.symbols('p0:8'))
    K=lambda v:v.dot(v)/2+delta.dot(v)
    transformed=A*symbols+lam
    need(S.expand(K(transformed)-K(symbols)-(sum(symbols)+alpha)/8)==0,
         'primitive W1 exact source Ward energy increment')
    a=theta_charge(F(int(alpha.p),int(alpha.q)),False)
    b=theta_charge(F(int(alpha.p),int(alpha.q)),True)
    need(a==b,'all exact parity/energy/total-charge coefficients through128')
    checks.append({'alpha':str(alpha),'all_energy_Q_parity_coefficients_equal_through_quarter_units':128,
                   'coefficient_keys':len(a)})

# The changed current's failure of original microscopic locality has a
# bounded exact three-mode CAR witness. Modes0,1 local;mode2 remote.
def annihilator(k):
    out=S.zeros(8)
    for bits in range(8):
        if bits&(1<<k):
            out[bits^(1<<k),bits]=(-1)**((bits&((1<<k)-1)).bit_count())
    return out
a=[annihilator(k) for k in range(3)]
n=[x.T*x for x in a]
parity=S.diag(*[(-1)**bits.bit_count() for bits in range(8)])
positive=(S.eye(8)+parity)/2
negative=(S.eye(8)-parity)/2
j=n[0]-n[1]
W=S.diag(*[S.exp(S.I*S.pi*j[k,k]/2) for k in range(8)])
B=a[2]+a[2].T
changed=positive*W+negative*W.conjugate().T
need(W*B==B*W,'original local even operation commutes with remote CAR')
comm=changed*B-B*changed
need(comm!=S.zeros(8),'parity-conditioned re-fielding fails remote locality')
need(comm==parity*B*(W-W.conjugate().T),'integrated remote commutator identity')
need(set((comm.conjugate().T*comm).eigenvals())=={S.Integer(0),S.Integer(4)},
     'nonzero commutator operator norm exactly2')
need(changed.conjugate().T*changed==S.eye(8),'locality mutant is still unitary')

# Original finite Arf compiler: iota=(f1,f2,f3,a,f1+f2+f3+a),
# chi_NSR=a is n4 in the original even exterior carrier.
basis=[m for m in range(32) if m.bit_count()%2==0]
index={m:i for i,m in enumerate(basis)}
chi=S.diag(*[(-1)**((m>>3)&1) for m in basis])
root=S.zeros(16)
for m in basis:
    if m&(1<<4) and not m&(1<<3):
        sign=(-1)**((m&((1<<4)-1)).bit_count())
        t=m^(1<<4)
        sign*=(-1)**((t&((1<<3)-1)).bit_count())
        t|=1<<3
        root[index[t],index[m]]=sign
anchor_comm=chi*root-root*chi
need(anchor_comm.rank()==4,'finite anchor sign is not invariant under D5 root E45')
need(S.trace(anchor_comm.T*anchor_comm)==16,'finite anchor commutator HS square')

result={'date':'2026-10-22','result':'EXACT_TIME_Q_COMPATIBLE_ODD_REFIELDING_BUT_NOT_ORIGINAL_LOCAL_NET',
        'rho_involution_Q_and_energy_preserved':True,
        'theta_controls':checks,
        'proof_of_all_energy_equality':'rho(p)=Pp+r,Pr=-r,Pdelta=delta+r,r²/2+delta.r=0',
        'even_sector_unchanged':True,
        'primitive_W1_target_action':'q -> A q+lambda, A=P R=R^5, lambda=(1/2)^5(-1/2)^3',
        'primitive_W1_charge_increment':1,
        'primitive_W1_K_increment':'(Q+alpha)/8',
        'primitive_W1_eighth_power':'q -> q+(1)^8; charge increases8',
        'primitive_W1_is_local_spinor_field':False,
        'odd_Cartan_beta_pullback':'H_beta becomes1-H_beta for beta=e0-e4',
        'proper_locality_obstruction':'For mean-zero localized beta-current and Pbeta=-beta,Wprime=Pi_even W+Pi_odd W*. For an odd remote CAR unitary B,[Wprime,B]=parity B(W-W*),nonzero.',
        'exact_finite_CAR_commutator_norm':2,
        'finite_Arf_anchor_sign_D5_root_commutator_rank':4,
        'finite_Arf_anchor_sign_D5_root_commutator_HS_squared':16,
        'finite_Arf_scope':'Excludes literal D5-equivariant identification of the existing anchor-bit sign with a uniform central NS/R boundary label; no refutation of the finite Arf theorem or all geometric Pfaffian maps.',
        'microscopic_P1_spinor_origin_derived':False,
        'new_energy_or_source_copy_chosen':False}
Path(__file__).with_name('odd-refielding-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
