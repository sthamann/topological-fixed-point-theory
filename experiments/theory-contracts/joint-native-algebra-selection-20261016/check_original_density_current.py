"""Original one-channel CAR density on the full shifted native quartic response.

No root weights or Klein signs fitted: ordered sea CAR computes them directly.
Finite window is exact for the stated energy<=33 outputs (proved in text).
"""
import json
from functools import lru_cache
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
import check_mixed_return as m
from check_fold_quartic_transport import cl,series_response

Q=m.Q
def exact_rank(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ_I).rank()
LOW,HIGH=-34,34
VAC=sum(1<<(j-LOW) for j in range(LOW,0))

def clean(v):return {w:c for w,c in v.items() if c}
def add(out,v,c):
    for w,a in v.items():out[w]=out.get(w,Q.zero)+c*a
def inner(a,b):return sum((m.conj(c)*b.get(w,Q.zero) for w,c in a.items()),Q.zero)

@lru_cache(None)
def bilinear(bits,target,source):
    a,b=target-LOW,source-LOW
    if not (0<=a<=HIGH-LOW and 0<=b<=HIGH-LOW):raise RuntimeError('Unproved cutoff')
    if not (bits>>b)&1 or (bits>>a)&1:return None
    lo,hi=sorted((a,b)); mask=((1<<hi)-1)^((1<<(lo+1))-1)
    sign=(-1)**((bits&mask).bit_count())
    return bits^(1<<a)^(1<<b),sign

def density(v,shift,channel=None,weights=None):
    """sum :a_(j+shift)^* a_j:, optionally fixed original residue channel."""
    out={}
    factors={j:(Q.one if weights is None else (weights(j) if callable(weights) else weights[j%8])) for j in range(LOW,HIGH+1)}
    for bits,c in v.items():
        for j in range(LOW,HIGH+1):
            if channel is not None and j%8!=channel:continue
            factor=factors[j]
            if not factor:continue
            target=j+shift
            if not LOW<=target<=HIGH:continue
            hit=bilinear(bits,target,j)
            if hit:
                state,sign=hit;out[state]=out.get(state,Q.zero)+sign*c*factor
    return clean(out)

def energy(bits):
    return sum(j*(((bits>>(j-LOW))&1)-int(j<0)) for j in range(LOW,HIGH+1))

def charges(bits):
    return tuple(sum(((bits>>(j-LOW))&1)-int(j<0) for j in range(LOW,HIGH+1) if j%8==k) for k in range(8))

def main():
    ld=m._lattice_data();T=s.Matrix.hstack(*map(s.Matrix,m._TARGET_SIMPLE))
    O=cl(T*ld['basis'].inv()/s.sqrt(2));U=s.zeros(8,4)
    for j in range(4):U[2*j,j]=1/s.sqrt(2);U[2*j+1,j]=-s.I/s.sqrt(2)
    M=s.Matrix(json.loads((m.R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())['explicit_Weyl_matrix'])
    V=cl(M.T*O*U);Vq=[[Q.from_sympy(V[k,a]) for k in range(8)] for a in range(4)]
    ps,ns=m.code();scale=s.diag(*[s.sqrt(n) for n in ns])
    @lru_cache(None)
    def word_state(word):
        if not word:return {VAC:Q.one}
        n,a=divmod(word[-1],8);m.need(a<4,'Same original positive current coordinates')
        v=word_state(word[:-1]);out={}
        for k,c in enumerate(Vq[a]):
            if c:add(out,density(v,8*(n+1),k),c)
        return clean(out)
    def to_car(poly):
        out={}
        for word,c in poly.items():add(out,word_state(word),c)
        return clean(out)
    # Sea anomaly and first actual root, directly in original Fourier modes.
    raised=density({VAC:Q.one},1)
    m.need(len(raised)==1 and all(charges(w)==(1,0,0,0,0,0,0,-1) for w in raised),'Original lowest source current selects wrap root without fitting')
    m.need(density(raised,-1)=={VAC:Q.one},'Original current level one, fixed sign')
    cases={}
    for alpha in (s.Rational(1,4),s.Rational(3,4)):
        delta=s.Matrix([(k+alpha)/8-s.Rational(1,2) for k in range(8)])
        b=cl(V.T*delta)
        polys=[series_response(f,[Q.from_sympy(x) for x in b],4) for f in ps]
        states=[to_car(f) for f in polys]
        m.need(all(energy(w)==32 and charges(w)==(0,)*8 for v in states for w in v),'Full actual shifted response has original energy32 and neutral channel charges')
        cargram=s.Matrix(5,5,lambda i,j:Q.to_sympy(inner(states[i],states[j])))
        fockgram=s.Matrix(5,5,lambda i,j:Q.to_sympy(m.inner(polys[i],polys[j])))
        m.need(cargram==fockgram,'Independent exact original CAR embedding reproduces whole Fock Gram')
        down=[density(v,-1) for v in states];up=[density(v,1) for v in states]
        m.need(all(energy(w)==31 for v in down for w in v),'Full lowering output at source energy31')
        m.need(all(energy(w)==33 for v in up for w in v),'Full raising output at source energy33')
        gm=s.Matrix(5,5,lambda i,j:Q.to_sympy(inner(down[i],down[j])))
        gp=s.Matrix(5,5,lambda i,j:Q.to_sympy(inner(up[i],up[j])))
        m.need(gp-gm==cargram,'All entries of original current CCR retained on actual code response')
        # Strong vector CCR control, includes every external output component.
        for v,vo,vi in zip(states,up,down):
            lhs=density(vo,-1);add(lhs,density(vi,1),-Q.one);add(lhs,v,-Q.one)
            m.need(not clean(lhs),'Full vector [J1,J-1]=I, no compression')
        GN=cl(scale.inv()*gm*scale.inv());BN=cl(scale.inv()*cargram*scale.inv())
        ratio_trace=s.simplify(s.trace(BN.inv()*GN)/5)
        scalar=cl(GN-ratio_trace*BN)==s.zeros(5)
        y=s.Symbol('y');poly=s.Poly((gm-y*cargram).det(),y).monic()
        distinct=s.gcd(poly,poly.diff()).degree()==0
        m.need(distinct,'Five distinct generalized response eigenvalues, no forced degenerate 3+2')
        ward_nonzero=0;erasure_norms=[]
        for a in range(4):
            va=[m.conj(v) for v in Vq[a]]
            difference=[va[(k+1)%8]-va[k] for k in range(8)]
            for i,(input_state,jup) in enumerate(zip(states,up)):
                erased=density(input_state,-8,weights=va)
                derivative={}
                for word,c in ps[i].items():
                    if a in word:
                        w=list(word);w.remove(a);w=tuple(w)
                        derivative[w]=derivative.get(w,Q.zero)+word.count(a)*c
                expected=to_car(series_response(derivative,[Q.from_sympy(x) for x in b],3))
                m.need(erased==expected,'Original CAR erasure equals full shifted cubic derivative field')
                erasure_norms.append(str(Q.to_sympy(inner(erased,erased))/ns[i]))
                lhs=density(jup,-8,weights=va);add(lhs,density(erased,1),-Q.one)
                rhs=density(input_state,-7,weights=difference)
                m.need(clean(lhs)==rhs,'Complete actual CAR current-erasure Ward identity, no endpoint compression')
                if rhs:ward_nonzero+=1
        m.need(ward_nonzero==20,'All twenty actual current-erasure commutator outputs nonzero')
        all_modes={};total=s.zeros(5)
        for mode in range(1,33):
            lowered=down if mode==1 else [density(v,-mode) for v in states]
            m.need(all(energy(w)==32-mode for v in lowered for w in v),'Complete source mode response at exact source energy')
            Gm=gm if mode==1 else s.Matrix(5,5,lambda i,j:Q.to_sympy(inner(lowered[i],lowered[j])))
            total+=Gm
            all_modes[str(mode)]={'rank':exact_rank(Gm),'raw_gram':m.mat(Gm)}
        m.need(total==32*cargram,'Full original c1 oscillator energy identity, all32 modes retained')
        m.need(all(not density(v,-33) for v in states),'Higher annihilation modes vanish by exact positive source energy')
        # ORIGINAL c1 geometric stress, not the folded c8 subalgebra.
        def stress(v,mode):
            return density(v,-mode,weights=lambda j:Q.from_sympy(j-s.Rational(mode,2)+alpha))
        stress_down=[stress(v,1) for v in states];stress_up=[stress(v,-1) for v in states]
        da=alpha-s.Rational(1,2)
        for v,ju,ldown,lup in zip(states,up,stress_down,stress_up):
            lhs=stress(ju,1);add(lhs,density(ldown,1),-Q.one);add(lhs,v,-Q.from_sympy(da))
            m.need(not clean(lhs),'Full original c1 stress/current Ward on neutral input')
            lhs=stress(lup,1);add(lhs,stress(ldown,-1),-Q.one);add(lhs,v,-Q.from_sympy(64+da**2))
            m.need(not clean(lhs),'Full original c1 stress Virasoro identity including twisted ground energy')
        GL=s.Matrix(5,5,lambda i,j:Q.to_sympy(inner(stress_down[i],stress_down[j])))
        Bdomain=DomainMatrix.from_Matrix(cargram).convert_to(s.QQ_I)
        Binverse=Bdomain.inv()
        AJ=(Binverse*DomainMatrix.from_Matrix(gm).convert_to(s.QQ_I)).to_Matrix()
        AL=(Binverse*DomainMatrix.from_Matrix(GL).convert_to(s.QQ_I)).to_Matrix()
        commutant=s.Matrix.vstack(s.kronecker_product(s.eye(5),AJ)-s.kronecker_product(AJ.T,s.eye(5)),
                     s.kronecker_product(s.eye(5),AL)-s.kronecker_product(AL.T,s.eye(5)))
        commutant_dim=25-DomainMatrix.from_Matrix(commutant).convert_to(s.QQ_I).rank()
        m.need(commutant_dim==1,'Original current/stress response pair has only scalar common commutant')
        allcharges=sorted({charges(w) for v in down+up for w in v})
        cases[str(alpha)]={'source_input_energy':32,'source_output_energies':[31,33],
          'input_car_state_counts':[len(v) for v in states],
          'lowering_state_counts':[len(v) for v in down],'raising_state_counts':[len(v) for v in up],
          'full_input_gram':m.mat(BN),'lowering_gram':m.mat(GN),
          'raising_gram':m.mat(cl(scale.inv()*gp*scale.inv())),
          'lowering_rank':exact_rank(gm),'raising_rank':exact_rank(gp),
          'lowering_response_scalar_relative_to_input_metric':scalar,
          'generalized_lowering_characteristic_polynomial':str(poly.as_expr()),
          'generalized_response_five_distinct_eigenvalues':distinct,
          'current_erasure_full_Ward_identities':20,
          'current_erasure_nonzero_Ward_outputs':ward_nonzero,
          'full_shifted_erasure_norms':erasure_norms,
          'all_positive_density_mode_responses':all_modes,
          'full_source_oscillator_energy_sum_equals_32_input_gram':True,
          'original_c1_stress_lowering_gram':m.mat(cl(scale.inv()*GL*scale.inv())),
          'original_c1_stress_lowering_rank':exact_rank(GL),
          'source_current_stress_response_common_commutant_dimension':commutant_dim,
          'original_stress_current_Ward_verified':True,
          'original_stress_Virasoro_with_ground_shift_verified':True,
          'retained_source_charge_sectors':[list(x) for x in allcharges],
          'full_vector_CCR_verified':True,'independent_CAR_Fock_Gram_equal':True}
        print(json.dumps({'alpha':str(alpha),'car_states':[len(v) for v in states],
             'down_rank':exact_rank(gm),'down_scalar_relative_to_B4':scalar,'sectors':len(allcharges)}),flush=True)
    out={'status':'PASS_EXACT_ORIGINAL_LOCAL_DENSITY_TO_FULL_FOLD_RESPONSE','date':'2026-10-17',
      'original_current':'J_m=sum_r :a_r^* a_(r+m):; [J_m,J_n]=m delta(m+n,0)',
      'folded_raising_mode':'J_-1=sum_(k=1..7) E_(e_k-e_(k-1),0) + E_(e0-e7,-1), with ordered CAR Klein convention, not independent fitted coefficients',
      'source_erasure_Ward':'[a_a,J_-1]=sum_r (conj(V_(r+1,a))-conj(V_(r,a))) a_(r-7)^* a_r; actual full 20 outputs retained',
      'source_energy_identity':'8K=Q(Q-1)/2+alpha Q+sum_m>0 J_-m J_m; on actual neutral input grade4 sum_m1..32 Gram(J_m R_b,4)=32 B4',
      'original_geometric_stress':'Lcal_m^alpha=sum_r (r+m/2+alpha) :a_r^* a_(r+m): + ((alpha-1/2)^2/2) delta_m0. Fold c8 stress samples original modes8n only: Lhat_n=Lcal_8n/8+21/64 delta_n0. Original modes +-1 are NOT folded canonical c8 stress histories.',
      'source_current_frame':m.mat(V),'cases':cases,
      'scope':'Given original single-circle CAR source, SAME sea/fold/M/O/fixed NS dictionary. Original density is a specified observable and local U1 generator, not an already P1-selected external coupling/profile or detector. J_+-1 are global Fourier controls; all charged output sectors retained, no executed projection or closed-C5 condition. Original local scalar multipliers fold to correlated off-diagonal matrices and are a different class from independently diagonal folded Cartan profiles. No 3+2 derivation or all-field no-go.'}
    Path(__file__).with_name('original_density_current_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])

if __name__=='__main__':main()
