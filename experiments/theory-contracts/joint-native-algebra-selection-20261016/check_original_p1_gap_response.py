"""Original chi=1/min|B| and its canonical CAR spectral component.

No profile chosen to obtain ranks; exact ALL-scalar-potential theorem is in
the companion proof. This checks the already selected actual quarter gap
on the same full source states, not a postulated physical instrument.
"""
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
import check_mixed_return as m
import check_original_density_current as car
from check_original_current_schur_kernel import source_packets

def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ_I)
def positive(A):
    return A==A.H and all(s.QQ_I.to_sympy(dm(A[:j,:j]).det()).is_positive for j in range(1,A.rows+1))

def main():
    theta=s.Symbol('theta',real=True);alpha_symbol=s.Symbol('alpha',real=True)
    phi=s.Function('phi')(theta);wave=s.Function('wave')(theta)
    gauge_action=s.exp(s.I*phi)*(-s.I*s.diff(s.exp(-s.I*phi)*wave,theta)+alpha_symbol*s.exp(-s.I*phi)*wave)
    m.need(s.simplify(gauge_action-(-s.I*s.diff(wave,theta)+(alpha_symbol-s.diff(phi,theta))*wave))==0,
           'Exact original circle gauge conjugation on arbitrary smooth functions')
    # Retain original source normalization; no choice of gap mode by desired rank.
    packets,norms=source_packets();cases={}
    for alpha,states in packets.items():
        candidates=[(abs(s.Integer(n)+alpha),n) for n in range(-2,3)]
        gap,mode=min(candidates)
        m.need(gap==s.Rational(1,4),'Actual original quarter spectral gap')
        m.need(mode==(0 if alpha==s.Rational(1,4) else -1),'Actual gap vector fixed by source B')
        # Outside these original modes |n+alpha|>1/4 by integer ordering,
        # so this is the exact infinite-spectrum minimum, not a finite search.
        sign=s.sign(mode+alpha)
        def gap_excitation(bits):
            occupation=(bits>>(mode-car.LOW))&1
            return int(sign)*(occupation-int(mode<0))
        m.need(all(gap_excitation(bits) in (0,1) for v in states for bits in v),'Canonical sea-normal-ordered gap excitation is an actual full-CAR projector')
        selected=[{bits:c for bits,c in v.items() if gap_excitation(bits)} for v in states]
        complement=[{bits:c for bits,c in v.items() if not gap_excitation(bits)} for v in states]
        B=s.Matrix(5,5,lambda i,j:m.Q.to_sympy(car.inner(states[i],states[j])))
        H=s.Matrix(5,5,lambda i,j:m.Q.to_sympy(car.inner(selected[i],selected[j])))
        Hc=s.Matrix(5,5,lambda i,j:m.Q.to_sympy(car.inner(complement[i],complement[j])))
        m.need(H+Hc==B,'All actual gap occupation and complement amplitudes retained')
        ranks=[dm(H).rank(),dm(Hc).rank()]
        A=(dm(B).inv()*dm(H)).to_Matrix()
        # Full exterior defect is P Qgap(1-P)Qgap P, not a performed projection.
        defect=H-(dm(H)*dm(B).inv()*dm(H)).to_Matrix()
        scalar=A==s.trace(A)/5*s.eye(5)
        m.need(positive(H) and positive(Hc),'Actual gap and complement responses are both strictly positive rank5')
        m.need(not scalar and A*A!=A,'Canonical gap response is nonscalar and is not a sharp input projector')
        m.need(positive(defect),'Full original exterior defect is strictly positive rank5, not discarded')
        cases[str(alpha)]={'gap':str(gap),'chi':str(1/gap),'source_gap_mode':mode,
          'source_gap_eigenvalue_sign':int(sign),'input_energy':32,
          'raw_input_gram':m.mat(B),'raw_gap_excitation_gram':m.mat(H),
          'raw_complement_gram':m.mat(Hc),'response_ranks':ranks,
          'gap_and_complement_strictly_positive':[positive(H),positive(Hc)],
          'generalized_effect':m.mat(A),'effect_scalar':scalar,
          'effect_is_projector':A*A==A,'raw_exterior_defect_gram':m.mat(defect),
          'exterior_defect_rank':dm(defect).rank(),
          'exterior_defect_positive_definite':positive(defect)}
        print(json.dumps({'alpha':str(alpha),'gap_mode':mode,'gap':str(gap),'ranks':ranks,
           'effect_scalar':scalar,'projector':A*A==A,'exterior_defect_rank':dm(defect).rank()}),flush=True)
    # Exact covector chain rule lambda=1/chi: D/chi=lambda D.
    lam,dl=s.symbols('lambda dlambda');D,V=s.symbols('D V',commutative=False)
    eps=s.Symbol('eps')
    derivative=s.expand((lam+eps*dl)*(D+eps*V)).coeff(eps,1)
    m.need(derivative==lam*V+dl*D,'Gap-normalized source tangent includes original scale variation')
    result={'status':'PASS_EXACT_ORIGINAL_GAP_CAR_RESPONSE','date':'2026-10-17','cases':cases,
      'gap_normalized_source_jet':'delta(D/chi)=lambda deltaD+deltaLambda D; deltaLambda=sign(E0)<u0,deltaB u0> on simple nonzero minabs branch',
      'all_scalar_potentials_theorem':'Balpha+v(theta) is unitarily equivalent to Balpha+mean(v). Smooth periodic phi has phi_prime=mean(v)-v; g=exp(i phi). Fixed holonomy implies mean(v)=0 locally, hence gap/chi unchanged exactly. Gap vector moves as g u0.',
      'scope':'Same stipulated scalar quarter source and conditional native/reference dictionary. Qgap=sign(E0):a_mode* a_mode: is canonical CAR spectral-component operator, NOT source-derived physical event/readout rule. Tests all full original input/output amplitudes, no executed compression, mode/profile weighting, new coupling or measurement. Does not establish isospectrality of full P1 relative action unless reflection/domains/reference/eta are jointly covariant. No arbitrary-source/P1 or nonlinear-event exclusion.'}
    Path(__file__).with_name('original_p1_gap_response_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(result['status'])

if __name__=='__main__':main()
