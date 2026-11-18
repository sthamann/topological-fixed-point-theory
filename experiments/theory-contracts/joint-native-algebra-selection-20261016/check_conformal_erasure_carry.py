"""Canonical source Virasoro action and original erasure carry.

All modes act in the original neutral Heisenberg Fock space without a cutoff.
Finite tests accompany an analytic all-word multiplicity theorem.
"""
from pathlib import Path
import json,sys
from functools import lru_cache
import sympy as s
import check_mixed_return as m
Q=m.Q


def pd(f,k):
    out={}
    for mon,c in f.items():
        for j,t in enumerate(mon):
            if t==k:
                z=mon[:j]+mon[j+1:];out[z]=out.get(z,Q.zero)+c
    return {z:c for z,c in out.items() if c}


def H(f,a,n):
    """h_a(n); complex bilinear metric pairs a and a+4."""
    if n==0:return {}
    if n<0:return m.mul({((-n-1)*8+a,):Q.one},f)
    return m.add((Q(n),pd(f,(n-1)*8+(a+4)%8)))


def L(f,n):
    out={}
    for mon,c in f.items():
        for j,k in enumerate(mon):
            mode=1+k//8
            if mode>n:
                var=(mode-n-1)*8+k%8
                z=tuple(sorted(mon[:j]+mon[j+1:]+(var,)))
                out[z]=out.get(z,Q.zero)+mode*c
    out={z:c for z,c in out.items() if c}
    if n>0:
        for a in range(4):
            for k in range(1,n):
                out=m.add((Q.one,out),(Q(k*(n-k)),pd(pd(f,(k-1)*8+a),(n-k-1)*8+a+4)))
    elif n<0:
        p=-n
        for a in range(4):
            for k in range(1,p):
                term=tuple(sorted(((k-1)*8+a,(p-k-1)*8+a+4)))
                out=m.add((Q.one,out),(Q.one,m.mul({term:Q.one},f)))
    return out


def word(f,seq):
    for n in reversed(seq):f=L(f,n)
    return f


def main():
    polys,norms=m.code();one=[{(a,):Q.one} for a in range(4)]
    erased=[pd(f,a) for f in polys for a in range(4)];tn=[n for n in norms for a in range(4)]
    vacuum={():Q.one};T=L(vacuum,-2)
    native=m._native_neutral_branch()
    m.need(Q.to_sympy(m.inner(T,T))==native['norms']['T']==4,'Original canonical stress tensor norm, c=8')
    families=[('native4',one,[1]*4,1),('native_code5',polys,norms,4),('native_erasure20',erased,tn,3)]
    primary_checks=0;gram_checks=0
    for name,fs,ns,h in families:
        for f in fs:
            m.need(L(f,0)==m.add((Q(h),f)),'Original conformal weight')
            for n in range(1,h+1):m.need(not L(f,n),'Original isotropic Fock state is primary');primary_checks+=1
        for i,f in enumerate(fs):
            for j,g in enumerate(fs):m.need(m.inner(f,g)==(Q(ns[i]) if i==j else Q.zero),'Original normalized multiplicity basis')
        level_words=[[( -1,)], [(-2,),(-1,-1)], [(-3,),(-2,-1),(-1,-1,-1)]]
        expected=[s.Matrix([[2*h]]),s.Matrix([[4*h+4,6*h],[6*h,4*h*(2*h+1)]]),
                  s.Matrix([[6*h+16,10*h,24*h],[10*h,h*(8*h+16),12*h*(3*h+1)],
                            [24*h,12*h*(3*h+1),24*h*(h+1)*(2*h+1)]])]
        for ws,ex in zip(level_words,expected):
            ds=[[word(f,w) for w in ws] for f in fs]
            for i in range(len(fs)):
                for j in range(len(fs)):
                    for a in range(len(ws)):
                        for b in range(len(ws)):
                            want=Q.from_sympy(ex[a,b])*ns[i] if i==j else Q.zero
                            m.need(m.inner(ds[i][a],ds[j][b])==want,'Virasoro descendant Gram equals common Gram times identity')
                            gram_checks+=1
    # Mode commutator on the actual code AND its conformal excursions.
    ward_checks=0
    for f in polys:
        for x in (f,L(f,-1),L(f,-2)):
            for a in range(4):
                for n in range(-3,4):
                    comm=m.add((Q.one,L(H(x,a+4,1),n)),(-Q.one,H(L(x,n),a+4,1)))
                    m.need(comm==m.add((-Q.one,H(x,a+4,n+1))),'All retained source erasure Ward coefficients')
                    ward_checks+=1
    carried=[];raw=[];full=[]
    for i,f in enumerate(polys):
        for a in range(4):
            j=4*i+a;g=erased[j];t=H(f,a+4,-1)
            m.need(not L(t,1) and L(t,2)==g,'Exact first carry highest-weight equations')
            r=m.add((Q.one,t),(-Q(7)/85,L(g,-2)),(Q(3)/170,L(L(g,-1),-1)))
            for n in range(1,6):m.need(not L(r,n),'Derived weight5 carry is primary')
            m.need(L(r,0)==m.add((Q(5),r)),'Carry conformal weight5')
            actual=H(L(f,-2),a+4,1)
            m.need(actual==m.add((Q.one,L(g,-2)),(Q.one,t)),'Unprojected original first geometric-erasure output')
            carried.append(r);raw.append(t);full.append(actual)
    for i in range(20):
        for j in range(20):
            target=Q(tn[i]) if i==j else Q.zero
            m.need(m.inner(raw[i],raw[j])==target,'Original raw carry ON20')
            m.need(m.inner(carried[i],carried[j])==Q(78)/85*target,'Derived primary Gram78/85 I20')
            m.need(m.inner(full[i],full[j])==Q(19)*target,'Full source Ward output Gram19 I20')
            m.need(not m.inner(carried[i],L(erased[j],-2)) and
                   not m.inner(carried[i],L(L(erased[j],-1),-1)),
                   'Primary remainder orthogonal to ENTIRE original20 level2 descendant block')
    # The same source coordinate quarter agrees at the original seed, but
    # not on the newly forced primary carry or general descendants.
    G=s.Matrix([[16,18],[18,84]]);coeff=G.inv()*s.Matrix([1,0])
    m.need(coeff==s.Matrix([s.Rational(7,85),-s.Rational(3,170)]),'Unique source Gram decomposition coefficients')
    # Literal existing native event lifts, also used in operator_checker.py.
    sys.path.insert(0,str(m.BASE))
    from operator_checker import native_reflection, _rays
    rays=_rays(m._lattice_data()['roots']);event_ids=[1,21,32];equivariance=0
    for event_id in event_ids:
        U=native_reflection(rays[event_id])
        @lru_cache(None)
        def monomial_image(mon):
            out={():Q.one}
            for k in mon:
                mode=k//8;a=k%8;j=a%4
                lin={}
                for row in range(4):
                    v=U[row,j] if a<4 else s.conjugate(U[row,j])
                    if v:lin[(mode*8+row+(4 if a>=4 else 0),)]=Q.from_sympy(v)
                out=m.mul(out,lin)
            return out
        def transform(f):return m.add(*[(value,monomial_image(mon)) for mon,value in f.items()])
        code_images=[transform(f) for f in polys]
        A=[[m.inner(polys[j],image)/norms[j] for image in code_images] for j in range(5)]
        for i in range(5):
            m.need(code_images[i]==m.add(*[(A[j][i],polys[j]) for j in range(5)]),'Actual native code action')
            for a in range(4):
                terms=[]
                for j in range(5):
                    for z in range(4):
                        v=A[j][i]*Q.from_sympy(s.conjugate(U[z,a]))
                        if v:terms.append((v,carried[4*j+z]))
                m.need(transform(carried[4*i+a])==m.add(*terms),'Literal native carry action rho5 tensor conjugate rho4')
                equivariance+=1
    out={
      'status':'PASS_EXACT_CANONICAL_CONFORMAL_ERASURE_CARRY',
      'contract_id':'UR.NATIVE.JOINT_OPERATOR.01','central_charge':8,
      'source_class':'ALL modes of original canonical T_E8; distinct from prior ALL neutral weight2 ZERO modes',
      'source_primary_carriers':[{'name':n,'dimension':len(fs),'weight':h} for n,fs,ns,h in families],
      'finite_primary_checks':primary_checks,'descendant_gram_checks_levels1to3':gram_checks,
      'erasure_ward_checks_including_excursions':ward_checks,
      'all_word_theorem':'For ON original primaries of weight h, descendants have Gram G_(8,h)(lambda,mu) delta_ij. Full Virasoro cyclic space is common unitary highest-weight cycle tensor multiplicity, modulo common nulls. Pure canonical stress histories cannot mix/distinguish logical code labels.',
      'erasure_ward':'[L_m,bar(u_a)(1)]=-bar(u_a)(m+1)',
      'first_higher_mode_source_output':'a_a L_-2 f_i=L_-2 g_ia+bar(u_a)(-1) f_i',
      'original_erased_weight3_level2_gram':[['16','18'],['18','84']],
      'unique_descendant_projection_coefficients':['7/85','-3/170'],
      'derived_primary_carry':'r_ia=bar(u_a)(-1) f_i -7/85 L_-2 g_ia +3/170 L_-1^2 g_ia',
      'derived_primary_carry_dimension':20,'derived_primary_carry_weight':5,
      'derived_primary_carry_gram_factor':'78/85','full_first_ward_output_gram_factor':'19',
      'literal_native_event_ray_ids':event_ids,'literal_native_carry_intertwining_checks':equivariance,
      'source_group_action':'rho5 tensor conjugate(rho4), actual native erasure20 representation; neutral, not charged E8 F20',
      'unweighted_carry_return_on_code':'312/85 I5',
      'unweighted_full_first_ward_return_on_code':'76 I5',
      'quarter_convention':'U_conf=exp(i*pi*L0/2); no physical duration inferred',
      'restricted_seed_quarter_characters':{'native4':'i','native_code5':'1','native_erasure20':'-i'},
      'new_primary_carry_quarter_characters':{'conformal':'i','native_J':'-i'},
      'scope':'Exact source-derived stress/erasure structure in given native J-positive Cartan-Fock embedding. No physical operation permission, P1 time/field assignment or3+2 selector. Pure-stress all-mode exclusion not united with prior other zero-mode exclusions; noncanonical/nonstress fields, mixed histories, other embeddings remain outside. Orthogonal primary decomposition is mathematical, never executed as projection.'}
    Path(__file__).with_name('conformal_erasure_carry_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
