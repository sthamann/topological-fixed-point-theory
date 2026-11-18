"""Prevent a rank-two/mode30 coincidence from replacing the operator dictionary."""
import json,itertools
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from check_fold_quartic_transport import cl

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    W=Path(__file__).parent;data=json.loads((W/'original_density_current_results.json').read_text())
    D=s.diag(s.sqrt(96),s.sqrt(8),s.sqrt(8),s.sqrt(8),1)
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],[0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=cl(s.diag(2,s.sqrt(12),s.sqrt(12),s.sqrt(12),s.sqrt(24))*raw/s.sqrt(48));rows={}
    for alpha,c in data['cases'].items():
        G=s.Matrix(c['all_positive_density_mode_responses']['30']['raw_gram'])
        need(DomainMatrix.from_Matrix(G).convert_to(s.QQ_I).rank()==2,'Actual mode30 response rank2')
        hits=[]
        for q in range(6):
            for pair in itertools.combinations([j for j in range(6) if j!=q],2):
                v=[N[:,j]-(s.sqrt(6)-1)*N[:,q]/5 for j in pair]
                P=cl(sum([w*w.T for w in v],s.zeros(5)))
                complement=cl(D.inv()*(s.eye(5)-P)*D)
                if all(s.expand(z)==0 for z in G*complement):hits.append([q,list(pair)])
        need(not hits,'No actual native polar-flag triplet equals mode30 nullspace')
        B=cl(D*s.Matrix(c['full_input_gram'])*D)
        G29=s.Matrix(c['all_positive_density_mode_responses']['29']['raw_gram'])
        sum_rank=DomainMatrix.from_Matrix(G29+G).convert_to(s.QQ_I).rank()
        need(sum_rank==5,'Joint rank3/rank2 source responses are injective, oblique direct sum')
        bd=DomainMatrix.from_Matrix(B).convert_to(s.QQ_I).inv()
        a29=bd*DomainMatrix.from_Matrix(G29).convert_to(s.QQ_I)
        a30=bd*DomainMatrix.from_Matrix(G).convert_to(s.QQ_I)
        product=(a29*a30).to_Matrix();overlap=s.expand(s.trace(product))
        need(overlap.is_positive,'Adjacent responses are not orthogonally complementary in actual input metric')
        rows[alpha]={'mode30_rank':2,'mode30_original_polar_flag_matches':hits,
             'mode29_rank':c['all_positive_density_mode_responses']['29']['rank'],
             'rank3_rank2_generalized_overlap':str(overlap),'not_orthogonally_complementary':True,
             'joint_response_rank':sum_rank,'oblique_algebraic_direct_sum':True}
    result={'status':'PASS_EXACT_ORIGINAL_RANK_TWO_DENSITY_MODE_FLAG_CONTROL',
      'date':'2026-10-17','cases':rows,
      'scope':'Mode30 found in COMPLETE m1..32 response family, not physically selected by desired rank2. Only literal original polar-support identification and orthogonally complementary-cut shortcut tested. The adjacent rank3/rank2 responses DO define an oblique algebraic direct sum, with joint rank5. A source annihilation Fourier mode is not a unitary order30 Coxeter action. No exclusion of positive relative normalization for that supplied mode pair, no selected profile or full nonlinear-readout/P1 no-go.'}
    (W/'original_rank_two_density_mode_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'cases':{a:{'matches':c['mode30_original_polar_flag_matches'],'not_orthogonally_complementary':c['not_orthogonally_complementary'],'joint_rank':c['joint_response_rank']} for a,c in rows.items()}},indent=2))

if __name__=='__main__':main()
