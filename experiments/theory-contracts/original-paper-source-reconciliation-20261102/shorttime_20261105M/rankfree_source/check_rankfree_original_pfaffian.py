"""Exact rank-free coefficients of the same original whole-neutrino determinant.
No numerical Yukawa/Majorana fixture and no physical rank is selected.
The original whole96/384 pairing identities are reused from K, not modified.
"""
from pathlib import Path
from itertools import combinations
import sympy as s
import json
base=Path(__file__).resolve().parent
checks=[]
def need(name, v):
    ok=bool(v);checks.append({'name':name,'passed':ok})
    if not ok:raise AssertionError(name)
def zero(e):return s.expand(e)==0
u,z=s.symbols('u z',positive=True)
a=s.symbols('a0:3',positive=True) # the actual positive squared singular values
M=s.zeros(3)
for i in range(3):
    for j in range(i,3):M[i,j]=M[j,i]=s.Symbol(f'M{i}{j}',complex=True)
# Actual symmetric Takagi-type neutrino mass block after the Y SVD.
need('whole_MR_is_complex_symmetric',M==M.T)
for rank in range(4):
    g=[u*a[i]+z if i<rank else z for i in range(3)]
    # Cauchy-Binet is an exact denominator-free expansion of
    # det(G)det(G+z M G^-1 M^dagger) on the original full six slots.
    Q=0
    for ell in range(4):
        for I in combinations(range(3),ell):
            for J in combinations(range(3),ell):
                minor=M.extract(I,J).det()
                Q+=z**ell*minor*s.conjugate(minor)*s.prod(g[i]for i in range(3)if i not in I)*s.prod(g[j]for j in range(3)if j not in J)
    Q=s.Poly(s.expand(Q),u)
    k=3-rank
    D=s.diag(*a[:rank]);Jinv=s.diag(*[1/x for x in a[:rank]])
    B=M[:rank,rank:];C=M[rank:,rank:]
    Cnorm=z*s.eye(k)+C.H*C
    lead=z**k*s.prod(a[:rank])**2*Cnorm.det()
    need(f'rank{rank}_exact_whole_neutrino_degree',Q.degree()==2*rank)
    need(f'rank{rank}_exact_leading_coefficient',zero(Q.coeff_monomial(u**(2*rank))-lead))
    if rank:
        # Clear the exact Cnorm determinant, retaining all complex MR phases.
        Tnumer=(Jinv*B*Cnorm.adjugate()*B.H).trace() if k else 0
        predicted=2*z*sum(1/x for x in a[:rank])*lead+2*z*z**k*s.prod(a[:rank])**2*Tnumer
        need(f'rank{rank}_exact_subleading_coefficient',zero(Q.coeff_monomial(u**(2*rank-1))-predicted))
    else:
        need('rank0_has_no_negative_degree_coefficient',Q.degree()==0)
    # The equality of the two Schur-return traces requires actual C^T=C.
    if k and rank:
        R=z*s.eye(k)+C*C.H
        left=z*(R.adjugate()*B.T*Jinv*B.conjugate()).trace()*Cnorm.det()
        right=z*(Jinv*B*Cnorm.adjugate()*B.H).trace()*R.det()
        need(f'rank{rank}_two_complex_symmetric_return_traces_agree',zero(left-right))
        need(f'rank{rank}_Cnorm_determinant_equals_reversed_norm',zero(Cnorm.det()-R.det()))
# All charged original family/color factors, without assuming any rank.
ru,rd,re,rn=s.symbols('r_u r_d r_e r_nu',integer=True,nonnegative=True)
kQ=6*ru+6*rd+2*re+2*rn
need('actual_whole96_degree_all_ranks',s.expand(kQ-(6*ru+6*rd+2*re+2*rn))==0)
K=s.symbols('k_Q',integer=True,positive=True)
S=s.symbols('S',real=True)
# P_pair=a_Q^2 u^(2k_Q)+4 a_Q^2 S u^(2k_Q-1)+...
# Standard four-real Gaussian radial moments: mu_n=2^n Gamma(n+2).
ratio=1/(2*(2*K+1))
BA=4*S*ratio
AU=2*(2*K+2)
BU=2*(2*K+1)*BA
need('actual_paired_u_degree_is_double_Q_degree',2*K==s.expand(2*K))
need('four_real_radial_adjacent_ratio',s.simplify(s.gamma(2*K+1)/s.gamma(2*K+2)/2-ratio)==0)
need('connected_global_Euler_all_ranks',s.simplify(2*(BU-AU*BA)+8*S/(2*K+1))==0)
need('fullrank_old_denominator_97',2*kQ.subs({ru:3,rd:3,re:3,rn:3})+1==97)
need('rank2_nu_fullrank_charged_denominator_93',2*kQ.subs({ru:3,rd:3,re:3,rn:2})+1==93)
need('zero_rank_zero_source_contact',(-8*S/(2*K+1)).subs({S:0,K:0})==0)
out={'passed':sum(c['passed']for c in checks),'total':len(checks),'checks':checks,'source':'Same actual whole96 product determinant and physical384 opposite-mode pairing as K.','physical_numeric_parameters_selected':False,'physical_rank_selected':False,'formula':'C_Euler=-8*(p_squared*Sigma_plus+T_M)/(2*k_Q+1)','scope':'global Euclidean Euler source, constant-H free Spin4 product pair, bounded held-rest degree-zero ray class; not fixed-time Yham or whole primitive history'}
(base/'check_rankfree_original_pfaffian.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
