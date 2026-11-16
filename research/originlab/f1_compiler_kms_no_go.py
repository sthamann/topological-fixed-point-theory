"""An exact F1 source-first Lindblad/KMS obstruction for the five-point
compiler's D4 quantum history mixture.

At D4 weight x=1/2 the full M4 quantum random-unitary two-step
transfer is INVERTIBLE but is not the finite-time member of a
KMS/trace-symmetric continuously CP quantum Markov semigroup.
Proof: its UNIQUE Hilbert-Schmidt selfadjoint log violates the
conditional complete-positivity (CCP) Choi generator requirement.

This rejects a PARTICULAR mathematical history rule, not the native
TFPT physical Q0, coherent spin lift or arbitrary non-KMS dilations.
"""
from itertools import permutations
import sympy as sp


def require(ok,why):
    if not ok:
        raise ArithmeticError(why)


def permutation_matrix(p):
    return sp.eye(len(p))[:,list(p)]


def source_superoperator(cycles,weight):
    matrices=[permutation_matrix(p) for p in cycles]
    basis=[sp.eye(4)[:,i]*sp.eye(4)[j,:]
           for i in range(4) for j in range(4)]
    return sp.Matrix.hstack(*[
       (sum((weight*m*v*m.T for m in matrices),
            sp.zeros(4))).reshape(16,1)
       for v in basis])


def conditional_choi_matrix(generator):
    out=sp.zeros(16)
    eye=sp.eye(4)
    for i in range(4):
        for j in range(4):
            E=eye[:,i]*eye[j,:]
            out+=sp.kronecker_product(
                E,generator[:,4*i+j].reshape(4,4))
    return out


def controls():
    four=[p for p in permutations(range(4))
          if all(p[p[p[p[i]]]]==i for i in range(4))
          and any(p[p[i]]!=i for i in range(4))]
    r=tuple((i+1)%4 for i in range(4))
    ri=tuple(r.index(i) for i in range(4))
    geometrical=[r,ri]
    other=[p for p in four if p not in geometrical]
    require(len(four)==6 and len(other)==4,
            "original anchored D4 compiler classes changed")
    x=sp.Rational(1,2)
    M=source_superoperator(geometrical,x/2)+source_superoperator(other,(1-x)/4)
    require(M==M.T and
            M.eigenvals()=={1:2,sp.Rational(-1,2):3,
               sp.Rational(-1,4):7,sp.Rational(1,2):1,
               sp.Rational(1,4):3},
            "D4 full CP history channel eigenstructure changed")
    S=M*M
    require(S.eigenvals()=={
        1:2,sp.Rational(1,4):4,sp.Rational(1,16):10}
        and S.det()!=0 and S==S.T,
        "D4 x=1/2 quantum two-step should be strictly positive and invertible")
    # Projectors of the THREE distinct eigenvalues of S.
    lambdas=(sp.Integer(1),sp.Rational(1,4),sp.Rational(1,16))
    id16=sp.eye(16)
    projs=[]
    for lam in lambdas:
        P=sp.eye(16)
        for mu in lambdas:
            if mu!=lam:
                P=P*(S-mu*id16)/(lam-mu)
        require(P*P==P and P.H==P, "spectral projector is not exact")
        projs.append(P)
    require(sum(projs,sp.zeros(16))==id16,
            "full quantum source eigenspaces not exhaustive")
    # The ONLY trace-symmetric generator of an invertible PSD transfer
    # at one physical time is log(S), by the positive functional calculus.
    # log(1)=0, log(1/4)=-2log2, log(1/16)=-4log2.
    Lunit=(-2*projs[1]-4*projs[2])
    assert Lunit==Lunit.T
    # The GKS/Lindblad conditional CP tangent criterion: if
    # exp(tL) is CP for infinitesimal t, its Choi derivative on
    # Ω⊥ must be nonnegative. Find exact -2/3*log2 witness.
    witness=sp.Matrix([0,1,0,-1,1,0,-1,0,0,-1,0,1,-1,0,1,0])
    Omega=sp.Matrix([1 if i==j else 0 for i in range(4) for j in range(4)])
    J=conditional_choi_matrix(Lunit)
    norm=(witness.T*witness)[0]
    numerator=(witness.T*J*witness)[0]
    require((Omega.T*witness)[0]==0 and norm==8 and
            numerator==sp.Rational(-16,3) and
            numerator/norm==sp.Rational(-2,3),
            "the quantum generator CCP exact negative witness changed")
    # D4-enhanced S4 equal-six case x=1/3 is even more singular:
    averaged=source_superoperator(four,sp.Rational(1,6))
    require(averaged.det()==0 and
            (averaged*averaged).eigenvals()=={
                1:2,sp.Rational(1,9):12,0:2},
            "S4 democratic quantum history zero modes changed")
    return {
      "status":"EXACT_D4_ANCHORED_TWO_STEP_KMS_CP_NONEMBEDDABILITY",
      "physical_premise":"native five-point code, six anchored 4cycles, D4 weights x=1/2",
      "exact_full_quantum_channel_kind":"CPTP/unital 4x4 random-unitary conjugation mixture",
      "M4_one_step_eigenspace_counts":{
          "1":2,"-1/2":3,"-1/4":7,"1/2":1,"1/4":3},
      "M4_two_step_eigenspace_counts":{"1":2,"1/4":4,"1/16":10},
      "M4_two_step_transfer_invertible":True,
      "unique_trace_symmetric_log_generator":"log(S)=log(2)*(-2P_1/4-4P_1/16)",
      "CCP_Choi_tangent_witness_norm":norm,
      "CCP_Choi_tangent_witness_quadratic_over_log2":str(numerator),
      "CCP_Choi_tangent_negative_eigenvalue_over_log2":"-2/3",
      "trace_symmetric_quantum_markov_embedding_exists":False,
      "non_trace_symmetric_other_logarithms_excluded":False,
      "full_S4_democratic_quantum_M4_two_step_has_2_zero_modes":True,
      "physical_original_Padm_and_coherent_amplitudes_evaluated":False,
      "original_Q0_R113_all_fields_and_OS_time_produced":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
