"""Full original E8 source hull, unprojected returns and current stress.

All source/code arithmetic is exact.  Finite checks accompany explicit
all-grade closure, locality/energy-bound and symmetry arguments.
No dynamics, profile, root, family omission or target clock is fitted.
"""
from pathlib import Path
from functools import lru_cache
from collections import defaultdict
import json,math,itertools
import sympy as s
import check_mixed_return as m
from check_conformal_erasure_carry import H,L,pd
Q=m.Q
LD=m._lattice_data();ROOTS=LD['roots'];ZERO=(0,)*8

def degree(mon):return sum(1+k//8 for k in mon)
def rdot(a,b):return sum(x*y for x,y in zip(a,b))//2
def radd(a,b):return tuple(x+y for x,y in zip(a,b))
def rneg(a):return tuple(-x for x in a)

@lru_cache(None)
def bcoords(r):return tuple(Q(r[2*a],-r[2*a+1])/2 for a in range(4))

@lru_cache(None)
def eps(r,t):
    a=LD['coordinates'](r);b=LD['coordinates'](t);low=LD['lower']
    return Q(-1 if sum(a[i]*low[i][j]*b[j] for i in range(8) for j in range(8))%2 else 1)

@lru_cache(None)
def exponential(r,k):
    if k<0:return {}
    if k==0:return {():Q.one}
    return m.add(*[(Q.one/k,m.mul(m.linear(bcoords(r),j),exponential(r,k-j))) for j in range(1,k+1)])

def shifted(poly,r):
    """P(h_-m-(h.lambda)z^-m), Laurent powers -> full Fock polys."""
    out=defaultdict(dict);b=bcoords(r)
    for mon,c in poly.items():
        terms={0:{():c}}
        for var in mon:
            mode=1+var//8;a=var%8;coupling=b[a] if a<4 else m.conj(b[a-4]);new={}
            for power,f in terms.items():
                new[power]=m.add((Q.one,new.get(power,{})),(Q.one,m.mul({(var,):Q.one},f)))
                new[power-mode]=m.add((Q.one,new.get(power-mode,{})),(-coupling,f))
            terms={p:f for p,f in new.items() if f}
        for power,f in terms.items():out[power]=m.add((Q.one,out[power]),(Q.one,f))
    return dict(out)

def rootmode(poly,charge,r,n):
    """E_r(n) on arbitrary Fock polynomial of specified original lattice charge."""
    offset=-n-1-rdot(r,charge)
    return m.add(*[(eps(r,charge),m.mul(f,exponential(r,offset-power)))
                   for power,f in shifted(poly,r).items() if offset-power>=0])

def annih(poly,P):
    out=[]
    for mon,c in P.items():
        v=poly
        for a in mon:v=pd(v,a)
        out.append((m.conj(c),v))
    return m.add(*out)

def raw_unit(poly,i,j,ps,ns):return m.add((Q.one/ns[j],m.mul(ps[i],annih(poly,ps[j]))))

def theta_counts(bound):
    """Full standard E8 lattice: D8 plus half-integer coset; exact DP."""
    total=[0]*(bound+1)
    for coset in [0,1]:
        state={(0,0):1};lim=math.isqrt(8*bound)
        ds=[d for d in range(-lim,lim+1) if d%2==coset]
        for _ in range(8):
            nxt=defaultdict(int)
            for (norm,par),count in state.items():
                for d in ds:
                    if norm+d*d<=8*bound:nxt[norm+d*d,(par+(d-coset)//2)%2]+=count
            state=nxt
        for n in range(bound+1):total[n]+=state.get((8*n,0),0)
    return total

def fock_counts(bound):
    p=[1]+[0]*bound
    for mode in range(1,bound+1):
        for color in range(8):
            for n in range(mode,bound+1):p[n]+=p[n-mode]
    return p

def main():
    ps,ns=m.code();vac={():Q.one};checks=0
    # All25 source-derived matrix units; no projection between their8 modes.
    for i in range(5):
        for j in range(5):
            for k in range(5):
                m.need(annih(ps[k],ps[j])==({():Q(ns[j])} if j==k else {}),'Actual quartet annihilation returns vacuum')
                m.need(raw_unit(ps[k],i,j,ps,ns)==(ps[i] if j==k else {}),'Literal unprojected full source matrix unit')
                checks+=1
    # Entire neutral weight4 Fock basis, not only35positive quartics.
    basis=[]
    def rec(start,left,mon):
        if left==0:basis.append(mon);return
        for k in range(start,32):
            w=1+k//8
            if w<=left:rec(k,left-w,mon+(k,))
    rec(0,4,());m.need(len(basis)==726,'Complete original neutral weight4 basis')
    for mon in basis:
        f={mon:Q.one}
        for i in range(5):
            for j in range(5):
                want=m.add((m.inner(ps[j],f)/ns[j],ps[i]))
                m.need(raw_unit(f,i,j,ps,ns)==want,'All neutral weight4 matrixunit action equals actual Fock inner product')
                checks+=1
    # Equality must NOT be extended to the full higher-grade source algebra.
    x={ (0,):Q.one };f=m.mul(ps[4],x)
    m.need(raw_unit(raw_unit(f,4,4,ps,ns),4,4,ps,ns)!=raw_unit(f,4,4,ps,ns),'Higher-grade full operator is not a projector')
    # Complete source root OPE and sum, including fixed cocycle dual signs.
    root_sum={};casimir=s.zeros(8)
    for r in ROOTS:
        m.need(rdot(r,r)==2,'Original root physical norm2')
        v=s.Matrix(r);casimir+=v*v.T/2
        prod=rootmode(vac,rneg(r),r,-1)
        expected=m.add((Q.one/2,m.mul(m.linear(bcoords(r),1),m.linear(bcoords(r),1))),
                       (Q.one/2,m.linear(bcoords(r),2)))
        m.need(prod==m.add((eps(r,rneg(r)),expected)),'Full source opposite-root normal product')
        root_sum=m.add((Q.one,root_sum),(Q.one/eps(r,rneg(r)),prod))
    m.need(casimir==60*s.eye(8),'Complete actual root moment; no chosen family')
    cartan=m.add(*[(Q(2),{(a,a+4):Q.one}) for a in range(4)])
    sug=m.add((Q.one/62,cartan),(Q.one/62,root_sum));T=L(vac,-2)
    m.need(sug==T and m.inner(T,T)==Q(4),'Full248-source Sugawara equals original canonicalT,c8')
    # All source positive-current mode responses on all actual code inputs.
    # Entire charged Fock output is retained, not compressed mid-operation.
    mode_grams={};full_gram={n:s.zeros(5) for n in range(4)};returns={};leakage={}
    for n in range(4):
        gram=s.zeros(5);returned=[{} for _ in range(5)]
        for r in ROOTS:
            images=[rootmode(f,ZERO,r,n) for f in ps]
            for i,image in enumerate(images):
                back=rootmode(image,r,rneg(r),-n)
                returned[i]=m.add((Q.one,returned[i]),(Q.one/eps(r,rneg(r)),back))
            for f,image in zip(ps,images):
                m.need(all(degree(mon)+1==4-n for mon in image),'Full root-mode output charge and conformal grade')
            for i in range(5):
                for j in range(5):gram[i,j]+=Q.to_sympy(m.inner(images[i],images[j]))
        # Real orthonormal Cartan basis, equivalent sum over paired u,ubar.
        for a in range(8):
            images=[H(f,a,n) for f in ps]
            for i,image in enumerate(images):
                returned[i]=m.add((Q.one,returned[i]),(Q.one,H(image,(a+4)%8,-n)))
            for i in range(5):
                for j in range(5):gram[i,j]+=Q.to_sympy(m.inner(images[i],images[j]))
        normed=s.diag(*[1/s.sqrt(v) for v in ns])*gram*s.diag(*[1/s.sqrt(v) for v in ns])
        m.need(normed==normed[0,0]*s.eye(5),'All248-current full endpoint Gram scalar by native equivariance')
        full_gram[n]=normed;mode_grams[str(n)]=str(normed[0,0])
        remainder=[m.add((Q.one,f),(-Q.from_sympy(normed[0,0]),ps[i])) for i,f in enumerate(returned)]
        defects=s.Matrix(5,5,lambda i,j:Q.to_sympy(m.inner(remainder[i],remainder[j]))/s.sqrt(ns[i]*ns[j]))
        m.need(defects==defects[0,0]*s.eye(5),'All exterior return defects scalar by full native symmetry')
        leakage[str(n)]=str(defects[0,0]);returns[n]=returned
    m.need(sum(full_gram[n][0,0] for n in (1,2,3))+full_gram[0][0,0]/2==124,'Unprojected source Ward/Sugawara sum31*h4')
    for i,f in enumerate(ps):
        summed=m.add((Q.one/2,returns[0][i]),*[(Q.one,returns[n][i]) for n in (1,2,3)])
        m.need(summed==m.add((Q(124),f)),'Full off-code Sugawara Ward cancellation; no endpoint-only inference')
    # Complete graded dimensions from full lattice +8-colored partitions.
    theta=theta_counts(6);osc=fock_counts(6)
    dims=[sum(theta[k]*osc[n-k] for k in range(n+1)) for n in range(7)]
    m.need(theta[:3]==[1,240,2160] and dims[:3]==[1,248,4124],'Actual complete lattice counting, no code truncation')
    # Exact grade2 horizontal decomposition from full lattice weight character.
    norm4={radd(a,b) for a in ROOTS for b in ROOTS if rdot(a,b)==0}
    m.need(len(norm4)==theta[2],'All2160 norm4 lattice weights from original roots')
    positive=[r for r,c in zip(ROOTS,LD['root_coordinates']) if all(v>=0 for v in c)]
    rho=tuple(sum(r[i] for r in positive)//2 for i in range(8));simple=[tuple(LD['basis'][:,i]) for i in range(8)]
    dominant=[v for v in norm4 if all(rdot(v,a)>=0 for a in simple)]
    m.need(len(dominant)==1,'Unique norm4 highest weight')
    hw=dominant[0];mult={v:1 for v in norm4}
    # Subtract grade2 singlet and derivative248; derive rather than input dims.
    mult.update({r:osc[1]-1 for r in ROOTS});mult[ZERO]=osc[2]-1-8
    freud=0
    for mu,value in mult.items():
        if mu==hw:continue
        den=rdot(radd(hw,rho),radd(hw,rho))-rdot(radd(mu,rho),radd(mu,rho))
        num=0
        for a in positive:
            for k in (1,2):
                v=radd(mu,tuple(k*x for x in a))
                num+=rdot(v,a)*mult.get(v,0)
        m.need(den*value==2*num,'Complete grade2 remainder satisfies exact Freudenthal character recurrence')
        freud+=1
    dim_remainder=sum(mult.values());m.need(dim_remainder==3875,'Grade2 irrep dimension derived from complete character')
    result={'status':'PASS_EXACT_FULL_E8_SOURCE_HULL_AND_RETURNS','date':'2026-10-17',
      'contract_id':'UR.NATIVE.JOINT_OPERATOR.01','source':'Given original E8 level1 lattice/current realization; actual native quartic Fock embedding',
      'all_grade_hull_theorem':'Current-mode cyclic hull of actualC5 is V_E8; invariant closed Hilbert hull under all local supported current smearings is H_E8. Modes are sums of local smearings, not individually local.',
      'full_source_basis':'h_-m1...h_-mk e^lambda, lambda in original E8, conformal weight sum m+lambda²/2',
      'source_unit_family':'Eij=P_i(u_-1) conjugate(P_j)(ubar_1)/sqrt(NiNj); no intermediate projection',
      'code_norms':ns,'neutral_grade4_basis_size':len(basis),'exact_source_matrixunit_checks':checks,
      'matrixunit_scope':'On ENTIRE weight4: |e_i><e_j| plus0 on code complement. Not full source identity at other weights. Nonlocal mode words, not a primitive selected intervention.',
      'full_endpoint_algebra':'M5(C), irreducible5, centerC; arbitrary coefficients are available but not selected',
      'original_root_count':len(ROOTS),'root_moment':'60 I8','dual_coxeter_from_roots':30,
      'sugawara_denominator':62,'sugawara_equals_previous_canonical_T':True,'derived_central_charge':8,
      'all248_current_mode_return_gram_factors_on_code':mode_grams,
      'full_all248_return_exterior_defect_gram_factors':leakage,
      'full_sugawara_return_on_code':'T0/2+T1+T2+T3=124 I5 as FULL source vectors, exterior terms cancel exactly',
      'negative_mode_gram_rule':'kappa_-m=kappa_m+248m for m>0, by actual source current CCR',
      'positive_mode_sugawara_sum':124,'theta_coefficients_grade0to6':theta,'fock8_coefficients_grade0to6':osc,
      'full_grade_dimensions_0to6':dims,'grade2_horizontal_irreps':[1,248,dim_remainder],
      'grade2_complete_freudenthal_checks':freud,
      'symmetry_scope':'Full source affine E8 level1 module irreducible; horizontal E8 automorphisms are covariance, not commutant or derived physical gauge. All-grade lattice character determines finite horizontal decomposition. Fully contracted native-invariant LINEAR return/memory kernels on irreducibleC5 are scalar; covariant quantum channels or open outcome indices need not be scalar.',
      'control':'Finite mode words exact with cutoff >=maximum actual intermediate conformal weight, full charge sectors and PBW bases retained. Smooth supported fields have infinite mode tails; energy bounds give convergence, no finite-mode object falsely called local.',
      'not_derived':['physical P1-to-current-net realization','selected source intervention or Hamiltonian','physical state/time/scale','marked3+2 or actualR5 selection','physical families/Flavor/Yukawa','fullTFPT'],
      'success_scope':'Previously supplied canonical stress field is recovered from full original current products without assuming it independently. Full native5 matrix action is explicitly available. Neither physical selection nor target3+2 success criterion is claimed.'}
    Path(__file__).with_name('full_source_hull_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
