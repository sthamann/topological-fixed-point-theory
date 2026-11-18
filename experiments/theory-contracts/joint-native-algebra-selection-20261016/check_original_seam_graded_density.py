"""Canonical ORIGINAL seam grading of the entire density response, no fitted pair."""
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
import check_mixed_return as m

def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ_I)
def positive(A):
    m.need(A==A.H,'Actual Hermitian response')
    return all(s.QQ_I.to_sympy(dm(A[:j,:j]).det()).is_positive for j in range(1,A.rows+1))

def main():
    W=Path(__file__).parent;prior=json.loads((W/'original_density_current_results.json').read_text())
    M=s.Matrix(json.loads((m.R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())['explicit_Weyl_matrix'])
    v=s.Matrix([0]*5+[s.Rational(1,2)]*3);t=M.T*v
    m.need(t==s.Matrix([s.Rational(1,2),s.Rational(1,4),0,-s.Rational(1,4)]*2),'Same original charged seam character and M')
    phases=[s.simplify(s.exp(2*s.pi*s.I*x)) for x in t]
    m.need(phases==[-1,s.I,1,-s.I]*2,'Original scalar CAR seam phase, no changed quarter holonomy')
    for degree in range(-16,17):
        m.need(all(s.simplify(phases[k]/phases[(k+degree)%8]-s.I**degree)==0 for k in range(8)), 'All seam Fourier residues, original source current/stress covariance')
    # Symbolic proof on all charge sectors, not only roots: phase difference is integer.
    p=s.symbols('p0:8',integer=True);N=s.Symbol('N',integer=True);a=s.Symbol('a');q=sum(p)
    K=N+sum(p[k]**2/2+((k+a)/8-s.Rational(1,2))*p[k] for k in range(8))
    lhs=sum(t[k]*p[k] for k in range(8))
    rhs=-2*K+(s.Rational(1,2)+a/4)*q
    residue=s.expand(lhs-rhs)
    m.need(s.expand(residue-(2*N+sum(x*x for x in p)-q+p[4]+p[5]+p[6]+p[7]))==0,'Vacuum-fixed global character/time/charge identity')
    data={}
    for alpha,c in prior['cases'].items():
        scale=s.diag(*[s.sqrt(n) for n in (96,8,8,8,1)])
        B=(scale*s.Matrix(c['full_input_gram'])*scale).applyfunc(s.simplify)
        G=[s.zeros(5) for _ in range(4)]
        for key,item in c['all_positive_density_mode_responses'].items():G[int(key)%4]+=s.Matrix(item['raw_gram'])
        m.need(sum(G,s.zeros(5))==32*B,'Complete graded response retains original energy normalization')
        m.need(all(positive(A) for A in G),'Every character-resolved response is strictly positive, not a sharp input sector')
        m.need(all(positive(32*B-A) for A in G),'Every relative effect is strictly below identity')
        even=G[0]+G[2];odd=G[1]+G[3]
        m.need(positive(even) and positive(odd),'Actual deck-even/deck-odd responses both full rank')
        inv=dm(32*B).inv();As=[(inv*dm(A)).to_Matrix() for A in G]
        m.need(sum(As,s.zeros(5))==s.eye(5),'Generalized effects sum identity with ACTUAL source metric')
        comm=s.Matrix.vstack(*[s.kronecker_product(s.eye(5),A)-s.kronecker_product(A.T,s.eye(5)) for A in As])
        dim=25-dm(comm).rank()
        m.need(dim==1,'Full canonical seam-resolved density response has scalar common commutant')
        data[alpha]={'canonical_mode_classes':[list(range(r if r else 4,33,4)) for r in range(4)],
           'class_raw_grams':[m.mat(A) for A in G],'class_ranks':[5]*4,'all_relative_effects_strict_between_zero_and_identity':True,
           'deck_even_and_odd_ranks':[5,5],'response_common_commutant_dimension':dim,
           'all_character_coarsenings_unsharp':True,'sum_equals_32_actual_input_gram':True}
    result={'status':'PASS_EXACT_ORIGINAL_SEAM_GRADED_COMPLETE_DENSITY_RESPONSE','date':'2026-10-17',
        'source_character':[str(x) for x in t],'original_CAR_phases':[str(x) for x in phases],
        'known_global_identity_reused':'t_v=exp(-4pi iK_alpha) exp(i pi(1+alpha/2)Q)=(-1)^Q W_alpha^-1; W_alpha=exp(4pi iK_alpha)exp(-i pi alpha Q/2)',
        'new_full_field_covariance':'t_v J_m t_v^-1=i^m J_m and t_v Lcal_m^alpha t_v^-1=i^m Lcal_m^alpha. Products carry total Fourier degree, all source OPE/Ward terms retained.',
        'cases':data,'scope':'Canonical original charged seam grading, entire positive-density response on existing neutral source energy32 packet, no selected profile/measurement or mode29/30. Relative effects summarize response tensors; no executed projection/instrument. Strictly positive class responses exclude only a sharp deterministic character/deck input-sector reading and common reducing3+2 spectrum, not arbitrary nonlinear histories/admissibility/P1. Known global character compensation not claimed new; odd CAR parity remains. Conditional native/reference dictionary and quotient locality retained.'}
    (W/'original_seam_graded_density_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'cases':{a:{'ranks':c['class_ranks'],'commutant':c['response_common_commutant_dimension']} for a,c in data.items()}},indent=2))

if __name__=='__main__':main()
