"""D4-symmetric hidden history, exact positive OS source and non-Markov observation.

Only a source-level mathematical comparison, not a dilation of the
specific native TFPT pulse U=I-2DD* or its actual original 20-state R113.
"""
import sympy as sp


def controls():
    a=sp.Rational(1,10)
    b=sp.Rational(1,20)
    W=sp.Matrix(4,4,lambda i,j:(1-2*a-b) if i==j else
                b if (i-j)%4==2 else a)
    q=sp.ones(4,1)/4
    assert W==W.T and W*q==q and W*sp.ones(4,1)==sp.ones(4,1)
    assert W.eigenvals()=={1:1,sp.Rational(3,5):1,sp.Rational(7,10):2}
    assert all(v>=0 for v in W.eigenvals())
    rot=sp.Matrix(4,4,lambda i,j:int(i==(j+1)%4))
    mirror=sp.Matrix(4,4,lambda i,j:int(i==(-j)%4))
    assert rot*W==W*rot and mirror*W==W*mirror
    # Fully consistent two-level weighted history along the mark graph.
    fine={(i,k):sum(q[i]*W[i,j]*W[j,k] for j in range(4))
          for i in range(4) for k in range(4)}
    coarse={(i,k):q[i]*(W**2)[i,k] for i in range(4) for k in range(4)}
    assert fine==coarse and len(fine)==16
    # OS reflection through a central state: conditional expectation factorization.
    os_gram=sp.Rational(1,4)*W.T*W
    assert os_gram.eigenvals()=={sp.Rational(1,4):1,
                                  sp.Rational(9,100):1,
                                  sp.Rational(49,400):2}
    # Coarse observation: one particular mark vs the three other marks.
    # The full source is D4 invariant, the chosen observable is not; a fixed
    # mark is a legitimate labelled source insertion.
    collect=sp.Matrix([[1,0,0,0],[0,sp.Rational(1,3),sp.Rational(1,3),sp.Rational(1,3)]])
    embed=sp.Matrix([[1,0],[0,1],[0,1],[0,1]])
    r1=collect*W*embed
    r2=collect*W**2*embed
    defect=r2-r1**2
    assert r1==sp.Matrix([[sp.Rational(3,4),sp.Rational(1,4)],
                          [sp.Rational(1,12),sp.Rational(11,12)]])
    assert defect==sp.Matrix([[sp.Rational(1,600),-sp.Rational(1,600)],
                             [-sp.Rational(1,1800),sp.Rational(1,1800)]])
    # Failure of strong lumpability: opposite mark B->0 is b while two
    # adjacent marks B->0 give a, so hidden memory survives projection.
    assert W[1,0]==W[3,0]==a and W[2,0]==b and a!=b
    # A false claim of an autonomous generator would require r2=r1^2.
    assert r2!=r1**2 and os_gram.is_positive_definite
    return {
      'status':'EXACT_D4_HIDDEN_MARKOV_OS_DILATION_NOT_AUTONOMOUS_OBSERVED_TIME',
      'global_states':4,
      'quarter_symmetric_mark_generator_rates':{'adjacent':'1/10','opposite':'1/20'},
      'global_kernel':'D4-symmetric stochastic symmetric PSD W',
      'global_spectrum':{'1':1,'3/5':1,'7/10':2},
      'exact_weighted_two_level_boundary_pairs':len(fine),
      'OS_positive_on_full_history':True,
      'global_source_is_reversible_and_finite_block_consistent':True,
      'observed_mark':'is state 0 vs states 1,2,3',
      'observed_one_step':str(r1.tolist()),
      'observed_two_step':str(r2.tolist()),
      'observed_semigroup_defect':str(defect.tolist()),
      'observed_autonomous_markov_transfer_exists_from_one_step':False,
      'mechanism':'hidden state memory survives non-lumpable projection',
      'original_TFPT_20_state_or_R113_physical_embedding_proved':False,
      'native_no_Lindblad_log_proved_by_this_toy':False,
      'F1_solved':False,
    }

if __name__=='__main__':
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
