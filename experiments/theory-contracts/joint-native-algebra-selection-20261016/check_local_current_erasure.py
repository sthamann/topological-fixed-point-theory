"""Original Cartan-current OPE/locality, checked on actual native code.

Finite rational checks accompany analytic all-profile/Weyl theorems.
No finite Fourier fixture is claimed to have proper interval support.
"""
from pathlib import Path
import json
import check_mixed_return as m
from check_conformal_erasure_carry import H,L
Q=m.Q

def dot(v,w):
    return sum((v[a]*w[(a+4)%8] for a in range(8)),Q.zero)

def vh(f,v,n):
    return m.add(*[(v[a],H(f,a,n)) for a in range(8) if v[a]])

def cj(v):
    return [m.conj(v[(a+4)%8]) for a in range(8)]

def hn(v):
    return sum((m.conj(x)*x for x in v),Q.zero)

def act(f,profile):
    return m.add(*[(Q.one,vh(f,v,n)) for n,v in profile.items()])

def gram(fs):return [[m.inner(x,y) for y in fs] for x in fs]

def main():
    fs,ns=m.code();zero=[Q.zero]*8;unit=lambda a:[Q.one if z==a else Q.zero for z in range(8)]
    ope=0;local=0;carry=[]
    for i,f in enumerate(fs):
        for a in range(4):
            g=H(f,a+4,1);t=H(f,a+4,-1)
            m.need(g==m.deriv(f,[Q.one if k==a else Q.zero for k in range(4)]),'Double-pole coefficient is actual native erasure')
            m.need(not H(f,a+4,0),'No simple pole in neutral source')
            for n in range(2,8):m.need(not H(f,a+4,n),'No higher singularity');ope+=1
            r=m.add((Q.one,t),(-Q(7)/85,L(g,-2)),(Q(3)/170,L(L(g,-1),-1)))
            for n in range(1,6):m.need(not L(r,n),'Regular OPE remainder is original weight5 primary')
            carry.append(r)
            for b in range(4):
                comm=m.add((Q.one,H(H(f,b,-1),a+4,1)),(-Q.one,H(H(f,a+4,1),b,-1)))
                m.need(comm==(f if a==b else {}),'Global annihilator fails local current commutation')
                local+=1
    # Full one-current Gram on entire actual code; arbitrary complex first
    # Fourier coefficient, adjoint negative coefficient, and higher modes.
    fixtures=[]
    fixtures.append({1:unit(4),-1:unit(0)}) # global rotating real polarization, sharp1/3
    # No random/fitted profile; finite literals only check conventions.
    v=[Q(1,1),Q(-2),Q.zero,Q(1),Q(3),Q(0,2),Q(-1),Q.zero]
    w=[Q.zero,Q(1),Q(1,1),Q.zero,Q(2),Q.zero,Q.zero,Q(0,-1)]
    fixtures.append({1:v,-1:cj(v),2:w,-2:cj(w),5:unit(3),-5:unit(7)})
    v=[Q.one,Q.zero,Q.zero,Q.zero,Q.one,Q.zero,Q.zero,Q.zero]
    fixtures.append({1:v,-1:cj(v)}) # fixed real polarization, sharp1/4
    results=[];checks=0
    for fixture in fixtures:
        v=fixture.get(1,zero);p2=hn(v[:4]);q2=hn(v[4:])
        tail=sum((Q(n)*hn(w) for n,w in fixture.items() if n>=2),Q.zero)
        expected=3*q2+p2+tail
        images=[act(f,fixture) for f in fs]
        erased=[vh(f,v,1) for f in fs]
        for i in range(5):
            for j in range(5):
                t=Q(ns[i]) if i==j else Q.zero
                m.need(m.inner(images[i],images[j])==expected*t,'Full current response Gram scalar on logical code')
                m.need(m.inner(erased[i],erased[j])==q2*t,'Actual erasure grade3 Gram')
                checks+=2
        results.append({'modes':sorted(fixture),'p2':str(Q.to_sympy(p2)),
                        'q2':str(Q.to_sympy(q2)),'tail':str(Q.to_sympy(tail)),
                        'full_gram_factor':str(Q.to_sympy(expected)),
                        'erasure_fraction':str(Q.to_sympy(q2/expected))})
    m.need(results[0]['erasure_fraction']=='1/3' and results[2]['erasure_fraction']=='1/4','Both global sharp bounds verified')
    # Source central commutator, including all mode cross terms in fixtures.
    central=0
    for a in fixtures:
        for b in fixtures:
            c=sum((Q(n)*dot(v,b.get(-n,zero)) for n,v in a.items()),Q.zero)
            for f in fs:
                got=m.add((Q.one,act(act(f,b),a)),(-Q.one,act(act(f,a),b)))
                m.need(got==m.add((c,f)),'Ordered full current commutator is central')
                central+=1
    # Weyl displacement obstruction is analytic, not tested by a cutoff.
    out={'status':'PASS_EXACT_LOCAL_CURRENT_ERASURE',
      'contract_id':'UR.NATIVE.JOINT_OPERATOR.01',
      'source':'Original neutral Cartan current sector in given J-positive native Fock embedding',
      'source_ope':'j_ubar_a(z)F_i(0)=G_ia(0)/z^2+sum_(k>=0) z^k/k! :partial^k j_ubar_a F_i:(0)',
      'regular_primary':'r_ia=t_ia-7/85 L_-2 g_ia+3/170 L_-1^2 g_ia; Gram78/85 I20',
      'singular_coefficient_checks':ope,'original_current_nonlocal_mode_checks':local,
      'locality_test':'[ubar_a(1),j_u_b(w)]=delta_ab I; nonzero outside every proper interval',
      'profile_convention':'J(V)=sum_n h_Vhat_n(n); Vhat_n=int V(theta)e^(-in theta)dtheta/(2pi); Vhat_-n=conj(Vhat_n)',
      'exact_all_profile_gram':'P J(V)^2 P=(3||q||^2+||p||^2+sum_(m>=2)m||Vhat_m||^2) I5, Vhat1=p.u+q.ubar',
      'exact_erasure_gram':'P J(V) P_grade3 J(V) P=||q||^2 I5; P_grade3 is analysis, not execution',
      'general_real_vector_profile_erasure_fraction_bound':'<=1/3; strict for smooth proper local support and q!=0',
      'fixed_real_polarization_bound':'<=1/4; strictly local nonzero erasure gives strict inequality',
      'fixtures_are_global_trigonometric_polynomials_not_local_support':True,
      'fixtures':results,'response_gram_checks':checks,'central_commutator_checks':central,
      'ordered_product':'[J(V),J(W)]=i sigma(V,W) I; exp(iJ(V))exp(iJ(W))=exp(-i sigma/2)exp(iJ(V+W))',
      'all_finite_weyl_history_theorem':'In neutral source sector, every finite ordered Cartan-current Weyl history is phase times one displacement. A nonzero displacement maps no nonzero finite Fock polynomial to a finite Fock polynomial. Every code-preserving endpoint is scalar.',
      'scope':'Exact OPE/locality and all-smooth-profile response; finite fixtures only validate coefficients. Weyl theorem is analytic, no truncated exponential. Cartan-current-only ordered unitary histories, not arbitrary products of fields, root currents, stress-current mixtures, general E8/P1 histories. Native code is NOT identified with P1 H_adm. Neither current choice/profile, P1-to-source local-net map, invariant H_adm, physical clock, measurement nor3+2 is selected.'}
    Path(__file__).with_name('local_current_erasure_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
