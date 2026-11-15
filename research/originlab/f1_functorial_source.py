"""Finite, pre-spatial CP/OS/projectively consistent source and exact nonuniqueness.

This is a complete MATHEMATICAL source for a three-state quantum Markov
process, not the physical original TFPT Q0, full R113, 3+1D or E8 VOA.
"""
from fractions import Fraction
from itertools import permutations
import sympy as sp


def require(ok, msg):
    if not bool(ok):
        raise AssertionError(msg)


def transition(n, q):
    q = sp.Rational(q)
    require(n >= 2 and 0 <= q <= 1, 'positive contracted transfer required')
    J = sp.ones(n) / n
    return q * sp.eye(n) + (1 - q) * J


def quantum_transfer(A, q):
    n = A.rows
    return sp.expand(q) * A + (1 - q) * sp.trace(A) * sp.eye(n) / n


def rank_and_gram(q):
    n=3
    time_q=(sp.S.One,q,q*q)
    probes=[sp.diag(1,0,0),sp.diag(0,1,0),sp.diag(0,0,1)]
    def inner(A,B):return sp.trace(A.H*B)/3
    folded=[quantum_transfer(a,t) for a,t in zip(probes,time_q)]
    gram=sp.Matrix(3,3,lambda i,j:sp.factor(inner(folded[i],folded[j])))
    require(gram == gram.H and gram.is_positive_definite, 'OS finite matrix not positive definite')
    return gram


def controls():
    q=sp.Rational(1,2)
    I=sp.eye(3)
    T=transition(3,q)
    T2=transition(3,q*q)
    require(T*T==T2, 'sewing from primitive transfer failed')
    require(T*sp.ones(3,1)==sp.ones(3,1), 'normalization failed')
    require(T==T.T and all(x>0 for x in T), 'reversibility or positivity failed')
    pi=sp.ones(3,1)/3
    require(T*pi==pi, 'stationary distribution failed')
    require(T.eigenvals()=={1:1,q:2}, 'primitive transfer spectral gap changed')
    # Exact consistency under forgetting the intermediate event on ALL 9 boundary pairs.
    for left in range(3):
        for right in range(3):
            fine=sum(pi[left]*T[left,mid]*T[mid,right] for mid in range(3))
            coarse=pi[left]*T2[left,right]
            require(fine==coarse,'fine-to-coarse endpoint source pushforward failed')
    gram=rank_and_gram(q)
    # Full one-step reflected classical history cylinder Gram:
    # for F depending only on the first future state and its mirror,
    # <theta F,F>=sum_x pi_x |sum_a T(x,a) F(a)|^2.
    # Arbitrary deeper positive-half cylinders obey the same conditional
    # expectation factorization, so this is not merely pointwise positivity.
    hist_gram=T.T*sp.diag(*list(pi))*T
    require(hist_gram == T2/3 and
            hist_gram.eigenvals()=={sp.Rational(1,3):1,sp.Rational(1,12):2},
            'full history reflection gram not PSD or sewing-consistent')
    # Permutation symmetry is not sufficient to fix the transfer parameter.
    for perm in permutations(range(3)):
        W=sp.Matrix(3,3,lambda i,j:int(perm[j]==i))
        require(W*T==T*W,'state-label automorphism covariance broken')
    # Heisenberg quantum transfer has 9 explicit positive Kraus operators plus the identity.
    K0=sp.sqrt(q)*I
    K=[K0]
    for i in range(3):
        for j in range(3):
            E=sp.eye(3)[:,i]*sp.eye(3)[j,:]
            K.append(sp.sqrt((1-q)/3)*E)
    kraus_sum=sum((k.H*k for k in K),sp.zeros(3))
    require(kraus_sum==I, 'quantum channel not TP/unital')
    sample=sp.Matrix([[1,2,0],[2,3,1],[0,1,2]])
    require(sum((k.H*sample*k for k in K),sp.zeros(3))==quantum_transfer(sample,q),
            'channel Kraus rule not source-identical')
    require(quantum_transfer(quantum_transfer(sample,q),q)==quantum_transfer(sample,q*q),
            'Heisenberg CP map not a semigroup')
    # R113-compatible names only for a DECLARED arithmetic control;
    # these are never the native original TFPT Q0/Z0/ZU/ZY/ZUY data.
    def four(q):
        p=transition(3,q)
        return {'Z0':sp.S.One,'ZU':sp.Rational(1,3), 'ZY':sp.Rational(1,3),
                'ZUY':sp.Rational(1,3)*p[0,0],
                'CUY':sp.Rational(1,3)*p[0,0]-sp.Rational(1,9)}
    a,b=four(sp.Rational(1,2)),four(sp.Rational(1,4))
    require((a['ZUY'],a['CUY'],b['ZUY'],b['CUY']) ==
            (sp.Rational(2,9),sp.Rational(1,9),sp.Rational(1,6),sp.Rational(1,18)),
            'source-rate nonuniqueness control changed')
    require(all(a[k]==b[k] for k in ('Z0','ZU','ZY')) and a['CUY'] != b['CUY'],
            'same lower source marks should not select connected return')
    # Demanding the EXACT SAME numeric transfer before and after two-step blocking
    # gives q^2=q and hence only frozen q=1 or instantaneous mixing q=0.
    qx=sp.Symbol('qx',real=True)
    require(sp.solve(qx*qx-qx,qx)==[0,1], 'primitive exact scale invariant kernel not trivial')
    # GNS/OS time from an eigenvalue after the source kernel is supplied:
    # tau(q)=-log q, tau(q^2)=2*tau(q). This does not fix q or kappa.
    require(sp.simplify(-sp.log(q*q)+2*sp.log(q))==0,
            'intrinsic additive operational clock failed')
    z=sp.symbols('z',positive=True)
    heat=1+8*sp.exp(-z)
    spectral_dim=sp.simplify(-2*z*sp.diff(sp.log(heat),z))
    require(sp.limit(spectral_dim,z,0,dir='+')==0 and
            sp.limit(spectral_dim,z,sp.oo)==0,
            'finite M3 source incorrectly yields 4-dimensional continuum')
    r20=transition(20,q)
    require(r20.eigenvals()=={1:1,q:19} and
            r20*sp.ones(20,1)/20 == sp.ones(20,1)/20,
            '20-state original event stationary comparison changed')
    return {
        'status':'EXACT_CP_OS_SEWING_AND_BLOCKING_PRESPATIAL_FINITE_SOURCE',
        'primitive_source':'3-state reversible complete-graph channel Pq=qI+(1-q)J',
        'local_realization':'9-dimensional noncommutative observable algebra M3(C)',
        'uniform_stationary_state':'rho=I3/3',
        'positive_Heisenberg_channel':'Phi_q(A)=q A +(1-q)Tr(A)I3/3',
        'Kraus_operators':10,
        'semigroup_identity':'Phi_q Phi_r=Phi_qr; Pq Pr=Pqr',
        'fine_to_coarse_exact_endpoint_pairs_verified':9,
        'OS_Gram_rank':gram.rank(),
        'OS_Gram_determinant':str(sp.factor(gram.det())),
        'full_reflected_classical_history_one_step_eigenvalues':{'1/3':1,'1/12':2},
        'state_permutation_automorphism_group':'S3 (all 6 permutations)',
        'eigenvalues_P_half':{'1':1,'1/2':2},
        'control_q_half':{k:str(v) for k,v in a.items()},
        'control_q_quarter':{k:str(v) for k,v in b.items()},
        'same_lower_three_source_marks_different_connected_UY':True,
        'fixed_scale_one_transfer_equation':'q=q^2 iff q in {0,1}',
        'nontrivial_continuous_source_requires_scale_dependent_q':True,
        'emergent_dimensionless_clock':'tau(q)=-log(q), additive under sewing',
        'generator_scale_and_original_clock_not_selected':True,
        'E8_affine_vertex_operator_algebra_derived':False,
        'twenty_state_analog_stationary':'I20/20',
        'finite_source_ultraviolet_heat_trace':'1+8*exp(-kappa*t) (operator-space gap)',
        'ultraviolet_spectral_dimension_limit':0,
        'four_dimensional_continuum_derived':False,
        'original_TFPT_Q0_identified':False,
        'original_R113_computed':False,
        'emergent_spatial_geometry_and_spin2_proved':False,
        'F1_solved':False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
