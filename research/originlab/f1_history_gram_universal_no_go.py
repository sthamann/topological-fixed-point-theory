"""Exact F1 five-point compiler source tests, NOT an original Q0 evaluation.

A: every D4 orbit weight x in [0,1] fails finite-time tracial/KMS
   symmetric CP embedding for the full M4 random-unitary two-step channel.
B: six independent four-cycle histories admit a D4-invariant Hermitian
   Choi/Gram kernel with FOUR free real TP directions near I6/6.
C: among nine equal-modulus quarter-phase unitary history sums, exactly
   one has an untwisted D4-covariant rank-one history kernel, and that
   unitary is only the negative half turn.

Only the original R113 all-field source, full outer P_adm, actual OS
reflection and absolute Z0 can select a physical kernel. F1 remains OPEN.
"""
from itertools import permutations, product

import sympy as sp
from f1_compiler_kms_no_go import (
    conditional_choi_matrix, permutation_matrix, source_superoperator
)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def compose(p, q):
    return tuple(p[q[i]] for i in range(4))


def histories():
    cycles = [
        p for p in permutations(range(4))
        if all(p[p[p[p[i]]]] == i for i in range(4))
        and any(p[p[i]] != i for i in range(4))
    ]
    r = tuple((i+1) % 4 for i in range(4))
    flip = tuple((-i) % 4 for i in range(4))
    group = {tuple(range(4))}
    queue = list(group)
    for g in queue:
        for generator in (r, flip):
            result = compose(g, generator)
            if result not in group:
                group.add(result)
                queue.append(result)
    inverse = lambda g: tuple(g.index(i) for i in range(4))
    actions = []
    lookup = {c:i for i,c in enumerate(cycles)}
    for g in group:
        actions.append(tuple(lookup[compose(compose(g, c), inverse(g))]
                             for c in cycles))
    short = [r, inverse(r)]
    long = [p for p in cycles if p not in short]
    require(len(cycles)==6 and len(short)==2 and len(long)==4
            and len(group)==8, "original D4 cycle census changed")
    return cycles, short, long, actions, r


def orbit_partition(items, actions, unordered=False):
    remaining = set(items)
    orbits = []
    while remaining:
        pair = next(iter(remaining))
        orb = set()
        for action in actions:
            mapped = tuple(action[i] for i in pair)
            if unordered:
                mapped = tuple(sorted(mapped))
            orb.add(mapped)
        require(orb <= remaining or not (orb - remaining),
                "orbit action escaped its source")
        remaining -= orb
        orbits.append(orb)
    return orbits


def history_gram_gate(cycles, actions):
    ordered = orbit_partition(
        [(i,j) for i in range(6) for j in range(6)], actions)
    unordered = orbit_partition(
        [tuple(sorted((i,j))) for i in range(6) for j in range(6)],
        actions, unordered=True)
    require(len(ordered)==7 and len(unordered)==6,
            "D4 orbitals on six original histories changed")
    # Six real symmetric invariant orbitals.
    kernels = []
    for orbit in unordered:
        K = sp.zeros(6)
        for i,j in orbit:
            K[i,j]=1
            K[j,i]=1
        kernels.append(K)
    # One genuinely imaginary antisymmetric invariant orbital.
    plus = next(o for o in ordered if {(j,i) for i,j in o}!=o)
    minus = {(j,i) for i,j in plus}
    require(minus in ordered, "missing transpose orbital")
    imaginary = sp.zeros(6)
    for i,j in plus:
        imaginary[i,j] = sp.I
    for i,j in minus:
        imaginary[i,j] = -sp.I
    kernels.append(imaginary)
    require(all(K.H == K for K in kernels),
            "D4 Gram coordinates are not Hermitian")
    perms = [permutation_matrix(c) for c in cycles]
    W = sp.Matrix.hstack(*[
        sp.Matrix([U[i,j] for i in range(4) for j in range(4)])
        for U in perms])
    require(W.rank()==6, "the six original four-cycles lost independence")
    def completeness(K):
        return sum((K[c,d]*perms[d].T*perms[c]
                    for c in range(6) for d in range(6)), sp.zeros(4))
    columns = []
    for K in kernels:
        C=completeness(K)
        columns.append(sp.Matrix([C[i,j] for i in range(4) for j in range(4)]))
    rank=sp.Matrix.hstack(*columns).rank()
    K0=sp.eye(6)/6
    require(rank==3 and completeness(K0)==sp.eye(4)
            and K0.is_positive_definite,
            "D4-invariant TP Gram dimension or strict positive interior changed")
    # Since W has full column rank, J(E_K)=W K W^dagger is CP iff K >= 0.
    return {
        "D4_ordered_orbitals":len(ordered),
        "D4_real_symmetric_orbitals":len(unordered),
        "D4_Hermitian_Gram_real_dimension":len(kernels),
        "trace_preserving_constraint_rank":rank,
        "remaining_affine_positive_TP_dimensions":len(kernels)-rank,
        "six_cycle_Kraus_vectors_independent":True,
        "strictly_positive_TP_center":"K=I6/6",
        "original_physical_K_measured":False,
    }


def universal_kms_no_go(short, long):
    x=sp.symbols("x",real=True)
    Phi=source_superoperator(short,x/2)+source_superoperator(long,(1-x)/4)
    S=Phi*Phi
    require(sp.factor(Phi.det())==x**4*(x-1)**8*(3*x-1)**2/1024,
            "D4 full quantum determinant changed")
    sample=sp.Rational(2,5)
    spectrum=[
        sp.Integer(1),
        sample**2,
        (1-sample)**2/4,
        (3*sample-1)**2/4,
    ]
    S0=S.subs(x,sample)
    I16=sp.eye(16)
    projectors=[]
    for eig in spectrum:
        P=I16
        for other in spectrum:
            if eig!=other:
                P=P*(S0-other*I16)/(eig-other)
        require(P==P.T and P*P==P,
                "quantum fixed spectral subspace ceased being a projector")
        projectors.append(P)
    require([sp.trace(P) for P in projectors]==[2,4,8,2]
            and sum(projectors,sp.zeros(16))==I16,
            "D4 quantum source spectral multiplicities changed")
    recovered = (projectors[0] + x**2*projectors[1]
                 + (1-x)**2/4*projectors[2]
                 + (3*x-1)**2/4*projectors[3])
    require(S.expand()==recovered.expand(),
            "one x-independent spectral decomposition fails for all weights")
    v=sp.Matrix([0,1,0,-1,1,0,-1,0,0,-1,0,1,-1,0,1,0])
    w=sp.Matrix([1 if i in (2,7,8,13) else 0 for i in range(16)])
    omega=sp.Matrix([int(i==j) for i in range(4) for j in range(4)])
    require((omega.T*v)[0]==0 and (omega.T*w)[0]==0,
            "Choi tangent witnesses must be perpendicular to identity")
    def coefficients(vector):
        return [
            sp.cancel((vector.T*conditional_choi_matrix(P)*vector)[0]
                      / (vector.T*vector)[0])
            for P in projectors
        ]
    require(coefficients(v)==[
        sp.Rational(1,3), -1, 0, sp.Rational(2,3)],
        "first exact all-x CCP witness changed")
    require(coefficients(w)==[
        sp.Rational(1,2), 1, -2, sp.Rational(1,2)],
        "second exact all-x CCP witness changed")
    require(sp.factor(4*x**3-(3*x-1)**2)==(x-1)**2*(4*x-1),
            "first witness sign factorization changed")
    # On 0<x<=1/4, (1-x)^4 >= 81/256 and
    # 8*x^2*(1-3*x) <= 32/243 (maximum at x=2/9).
    require(sp.Rational(81,256)>sp.Rational(32,243),
            "second witness interval bound not strict")
    return {
      "status":"ALL_D4_CLASSICAL_WEIGHTS_FAIL_TRACIAL_KMS_QUANTUM_TIME",
      "M4_two_step_generic_eigenvalues":{
        "1":2,"x^2":4,"(1-x)^2/4":8,"(3*x-1)^2/4":2},
      "singular_source_weights":["0","1/3","1"],
      "first_witness":"-2 log(x) + (4/3) log(abs(3*x-1)/2), negative for 1/4<x<1 except 1/3",
      "second_witness":"2 log(x)-4 log((1-x)/2)+log(abs(3*x-1)/2), negative for 0<x<=1/4",
      "all_x_in_closed_interval_excluded_in_tracial_KMS_class":True,
      "nontracial_or_nonselfadjoint_or_dilated_sources_excluded":False,
    }


def coherent_d4_gate(cycles, actions, r):
    mats=[permutation_matrix(c) for c in cycles]
    phases=(1,-1,1j,-1j)
    valid=[]
    for tail in product(phases,repeat=5):
        weights=(1,*tail)
        a=[[0j for _ in range(4)] for _ in range(4)]
        for p,z in zip(cycles,weights):
            for col,row in enumerate(p):
                a[row][col]+=z
        unitary=all(
            sum(a[k][i].conjugate()*a[k][j] for k in range(4))
            == (4 if i==j else 0)
            for i in range(4) for j in range(4))
        if not unitary:
            continue
        # Kernel aa* invariant under the untwisted D4 conjugation action.
        covariant=all(
            weights[action[i]]*weights[action[j]].conjugate()
            == weights[i]*weights[j].conjugate()
            for action in actions for i in range(6) for j in range(6))
        valid.append((weights,covariant))
    chosen=[weights for weights,covariant in valid if covariant]
    require(len(valid)==9 and len(chosen)==1,
            "quarter-phase D4 covariant coherent clock count changed")
    chosen=chosen[0]
    U=sum((int(z.real)*mats[i] for i,z in enumerate(chosen)),sp.zeros(4))/2
    halfturn=compose(r,r)
    require(all(z.imag==0 for z in chosen)
            and U==-permutation_matrix(halfturn)
            and U**2==sp.eye(4),
            "the sole D4 covariant coherent sum is not the negative half turn")
    return {
      "relative_quarter_phase_assignments":4**5,
      "unitary_history_sums":len(valid),
      "untwisted_D4_invariant_rank_one_Grams":1,
      "only_D4_invariant_unitary":"negative half-turn permutation, U^2=I",
      "projective_native_spin_action_evaluated":False,
    }


def supplied_spectral_counterexample():
    # The spectra supplied in the user's data do not imply a quantum channel.
    # T tensor a fixed-point-free involution has EXACTLY those 16 eigenvalues
    # and squared eigenvalues, yet is not even Hermiticity-preserving on M4.
    T=(sp.ones(4)-sp.eye(4))/3
    F=permutation_matrix((1,0,3,2))
    M=sp.kronecker_product(T,F)
    vecE00=sp.eye(16)[:,0]
    out=(M*vecE00).reshape(4,4)
    require(M.eigenvals()=={
        1:2,sp.Rational(1,3):6,sp.Rational(-1,3):6,-1:2}
        and (M*M).eigenvals()=={1:4,sp.Rational(1,9):12}
        and out!=out.H,
        "provided spectral control must fail quantum Hermiticity")
    return {
      "supplied_spectrum_matched_by_T_tensor_involution":True,
      "that_spectrum_implies_CP_or_quantum_time":False,
    }


def controls():
    cycles,short,long,actions,r=histories()
    return {
        "theorem_scope":"abstract original five-point six-cycle compiler; not native all-field Q0",
        "all_weights_no_go":universal_kms_no_go(short,long),
        "history_kernel":history_gram_gate(cycles,actions),
        "coherent_symmetry":coherent_d4_gate(cycles,actions,r),
        "spectral_typing":supplied_spectral_counterexample(),
        "original_R113_full_source_bundle_available":False,
        "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
