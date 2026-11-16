"""Exact F1 projective-consistency gate for the constraint-selected P1 seam.

The previous source-first candidate used matched n-region four-sector
histories and common cyclic invariance to produce a unique GHZ_4 pure
state at each n>=2. These unique pure states are NOT a compatible
family under ordinary local partial traces. Two precise alternatives
are tested: nonlinear conditioning on a physical cap (probability
4^(-k)) and a nonlocal but deterministic CPTP copy-code decoder.

The categorical copying map is an established dagger special
commutative Frobenius structure in the supplied basis; it DOES NOT
clone unknown quantum states and cannot by itself be native R113 Q0.
All conclusions are about this exact finite toy candidate, NOT F1.
"""
import json
import sympy as sp


def require(ok, msg):
    if not ok:
        raise ArithmeticError(msg)


def embedding(d, n):
    """Columns are |j>^tensor n, j=0..d-1."""
    e=sp.zeros(d**n,d)
    for j in range(d):
        idx=sum(j*d**k for k in range(n))
        e[idx,j]=1
    require(e.T*e==sp.eye(d),"repetition-code embedding lost isometry")
    return e


def pure_density(d,n,character=0):
    e=embedding(d,n)
    phases=sp.Matrix([sp.I**(character*j) for j in range(d)])
    if d!=4 and character:
        raise ValueError("quarter characters need four sectors")
    ket=e*phases/sp.sqrt(d)
    return sp.simplify(ket*ket.H)


def mixture_density(d,n):
    e=embedding(d,n)
    return e*e.T/d


def trace_last(rho,d,n):
    size=d**(n-1)
    return sp.Matrix(size,size,lambda i,j:
                     sum(rho[d*i+k,d*j+k] for k in range(d)))


def cyclic_shift(d):
    return sp.Matrix(d,d,lambda i,j:int(i==(j+1)%d))


def cyclic_characters_gate():
    d=4
    g=sp.kronecker_product(cyclic_shift(4),cyclic_shift(4))
    match=embedding(d,2)*embedding(d,2).T
    reflection=sp.Matrix(d,d,lambda i,j:int(i==(-j)%d))
    mirror=sp.kronecker_product(reflection,reflection)
    kets=[]
    for m in range(4):
        v=embedding(d,2)*sp.Matrix([sp.I**(m*j)
                                     for j in range(d)])/2
        p=sum((sp.I**(m*k)*g**k for k in range(4)),sp.zeros(16))/4
        require(p.H==p and p*p==p and p.rank()==4
                and p*match==v*v.H and (p*match).rank()==1
                and g*v==sp.I**(-m)*v,
                "physical cyclic character projector closure changed")
        kets.append(v)
    require(all(mirror*kets[m]==kets[(-m)%4] for m in range(4)),
            "geometric reflection does not pair quarter characters")
    return {
      "Z4_character_sector_count":4,
      "matching_rank_each":4,
      "cyclic_character_projector_rank_each":4,
      "joint_rank_each":1,
      "character_phase_families":["1","i","-1","-i"],
      "reflection_pairs_characters":"m <-> -m (mod 4)",
      "neutral_character_source_physically_selected":False,
      "nontrivial_quarter_spin_source_physically_selected":False
    }


def projective_state_gate():
    rows={}
    for d in (2,3,4):
        rho2=pure_density(d,2)
        rho3=pure_density(d,3)
        mu2=mixture_density(d,2)
        mu3=mixture_density(d,3)
        reduced=trace_last(rho3,d,3)
        require(reduced==mu2 and trace_last(mu3,d,3)==mu2,
                "source quantum marginal reductions changed")
        require(rho2.rank()==1 and reduced.rank()==d
                and rho2*rho2==rho2 and sp.trace(reduced)==1,
                "pure-vs-mixed source rank obstruction changed")
        D=rho2-reduced
        ev=D.eigenvals()
        expected={sp.Rational(d-1,d):1, -sp.Rational(1,d):d-1}
        if d<4: expected[sp.Integer(0)]=d*d-d
        else: expected[sp.Integer(0)]=12
        require(ev==expected,
                "quantitative trace-distance obstruction changed")
        distance=sp.Rational(d-1,d)
        fidelity=(rho2*reduced).trace()
        require(fidelity==sp.Rational(1,d),
                "pure-vs-marginal GHZ overlap changed")
        rows[str(d)]={
            "seam_pure_rank":1,
            "trace_fine_to_coarse_rank":d,
            "trace_distance":str(distance),
            "pure_coarse_overlap":str(fidelity),
            "ordinary_partial_trace_projective_consistent":False,
            "dephased_mixtures_partial_trace_consistent":True,
        }

    # The witness is invariant under the same global Z4 shift used
    # by the constraint. It distinguishes pure vs dephased reductions
    # even within this algebraically invariant observable sector.
    rho=pure_density(4,2)
    mu=mixture_density(4,2)
    G=sp.kronecker_product(cyclic_shift(4),cyclic_shift(4))
    require((rho*G).trace()==1 and (mu*G).trace()==0,
            "cyclic-invariant observable cannot see source coherence")
    return {
        "dimensions":rows,
        "d4_source_discrepancy_trace_distance":"3/4",
        "d4_shift_invariant_joint_probe_pure":"1",
        "d4_shift_invariant_joint_probe_partial_trace":"0",
        "observable_joint_shift_commutes_with_constraint":True,
        "general_pure_marginal_extension_rule":
          "pure AB marginal => full ABC=AB tensor C (no correlated new region)",
    }


def cap_and_isometric_coarse_gate():
    d=4
    rho2=pure_density(d,2)
    rho3=pure_density(d,3)
    # A physical postselection effect must have operator norm <=1:
    # a unit-norm coherent cap recovers the pure coarse state only
    # with an unconditional success probability 1/d per deletion.
    cap=sp.kronecker_product(sp.eye(d*d),sp.ones(1,d)/sp.sqrt(d))
    coarse_raw=cap*rho3*cap.H
    require(coarse_raw==rho2/d and sp.trace(coarse_raw)==sp.Rational(1,d)
            and cap.H*cap != sp.eye(d*d*d),
            "physical cap probability and absolute normalization changed")
    # Nonlinear renormalization is part of conditioning, not a
    # deterministic quantum channel on density matrices.
    # An entirely different deterministic channel CAN decode the
    # basis-copy isometry, but it necessarily changes AB's marginal.
    delta=embedding(d,2)   # delta |j> = |jj>
    eye=sp.eye(d)
    require(delta.H*delta==eye
            and sp.kronecker_product(delta,eye)*delta
                ==sp.kronecker_product(eye,delta)*delta,
            "copy spider lacks isometry or coassociativity")
    W=sp.kronecker_product(eye,delta)   # AB -> ABC, copy B
    require(W.H*W==sp.eye(d*d) and W*rho2*W.H==rho3
            and W.H*rho3*W==rho2,
            "nonlocal reversible copy-code decoding failed")
    # CPTP decoder on the full fine Hilbert space:
    # E(X)=W^dag X W + Tr((I-WW^dag) X) tau, tau=I_16/16.
    # Both terms are CP and Tr E(X)=Tr X for all X.
    identity=sp.eye(d**3)
    tau=sp.eye(d**2)/d**2
    def decode(X):
        return W.H*X*W+sp.trace((identity-W*W.H)*X)*tau
    require(decode(rho3)==rho2
            and decode(identity)==sp.eye(d*d)
                  +sp.Rational(d**3-d**2,d*d)*sp.eye(d*d),
            "deterministic decoder not mathematically normalized")
    # The last check above is intentionally NOT unital: the decoder
    # is TP (input/output dimensions differ), not a unital channel.
    trpreserved=all(
        sp.trace(decode(sp.eye(d**3)[:,i]*sp.eye(d**3)[j,:]))
        ==int(i==j)
        for i,j in [(0,0),(15,15),(42,42),(63,63),(1,2),(2,1)])
    require(trpreserved,"nonlocal decoder must be trace preserving")

    # Raw categorical epsilon=sum_j <j| restores GHZ without norm loss
    # on this code, but epsilon^dag epsilon has max eigenvalue d > 1.
    # It is therefore NOT a valid effect on an unrestricted Hilbert space.
    epsilon=sp.ones(1,d)
    require((epsilon.H*epsilon).eigenvals()=={d:1,0:d-1}
            and sp.kronecker_product(eye,epsilon)*delta==eye,
            "Frobenius deletion is not a normalized physical effect")
    return {
        "physical_normalized_cap_single_step_success":"1/4",
        "physical_normalized_cap_k_steps_success":"4^(-k)",
        "absolute_Z0_compensation_per_step_required":True,
        "copy_map_delta_isometric_and_coassociative":True,
        "nonlocal_CPTP_decoder_on_fine_algebra_constructed":True,
        "nonlocal_decoder_preserves_old_local_AB_marginal":False,
        "unnormalized_categorical_counit_is_physical_effect":False,
        "physical_R113_block_map_identified":False,
    }


def overlap_and_clustering_gate():
    d=4
    rho=pure_density(d,2)
    mu=mixture_density(d,2)
    lam=sp.symbols("lam",real=True)
    state=sp.simplify(lam*rho+(1-lam)*mu)
    require(state.subs(lam,sp.Rational(1,2)).eigenvals()=={
        sp.Rational(5,8):1,sp.Rational(1,8):3,0:12
    },"partial environment-label overlap spectrum changed")
    require(sp.simplify((rho*state).trace())==(1+3*lam)/4,
            "incomplete seam information / coherence overlap changed")
    # In this equi-overlap family off-diagonal coherence is lam,
    # and D(psi,mixture_interpolated)=3(1-lam)/4 for 0<=lam<=1.
    z=sp.diag(1,sp.I,-1,-sp.I)
    require(sp.trace(mu*sp.kronecker_product(z,sp.eye(d)))==0
            and sp.trace(mu*sp.kronecker_product(z,z.H))==1,
            "finite local long-range label correlations disappeared")
    return {
        "environment_offdiagonal_overlap":"lambda",
        "allowed_explicit_source_interval":"0 <= lambda <= 1",
        "coarse_2state_fidelity":"(1+3lambda)/4",
        "coarse_2state_trace_distance":"3(1-lambda)/4",
        "lambda_half_spectrum":["5/8 (x1)","1/8 (x3)","0 (x12)"],
        "connected_label_correlation_at_any_distance":1,
        "gapped_short_range_unique_vacuum_with_physical_label_observables":False,
        "requires_original_R113_overlap_or_physical_observable_dictionary":True,
    }



def infinite_matched_projective_uniqueness_gate():
    """General theorem, exhibited exactly for an arbitrary 4x4 code density.

    For EVERY n>=2, projectively consistent rho_n supported on the
    diagonal code E_n C^d and invariant under S^{tensor n} must equal
    mu_n=E_n E_n^dagger / d. Marginalizing ANY supported rho_(n+1)
    kills ALL offdiagonal code coherence. Iterating one extension at
    every n forces every rho_n diagonal. Cyclic covariance forces
    its weights equal; trace-one fixes their value to 1/d.

    The symbolic arbitrary density test below checks the full
    offdiagonal loss (12 coefficients for d=4), not merely GHZ.
    """
    d=4
    M=sp.Matrix(d,d,lambda i,j:sp.symbols("h%d%d"%(i,j)))
    E2=embedding(d,2)
    E3=embedding(d,3)
    rho3=E3*M*E3.H
    coarse=trace_last(rho3,d,3)
    diagonal=sp.diag(*[M[j,j] for j in range(d)])
    require(coarse==E2*diagonal*E2.H,
            "general code-supported fine source retains forbidden coherence")
    require(all(not coarse.has(M[i,j])
                for i in range(d) for j in range(d) if i!=j),
            "some of the twelve original offdiagonal source words survived")
    shift=cyclic_shift(d)
    w=sp.symbols("w0:4",real=True)
    source=sp.diag(*w)
    difference=shift*source*shift.T-source
    A,_=sp.linear_eq_to_matrix(list(difference),w)
    require(A.rank()==3,"cyclic density symmetry no longer forces uniformity")
    equations=[difference[i,i] for i in range(d)]
    solutions=sp.solve(equations+[sum(w)-1],w, dict=True)
    require(solutions==[dict(zip(w,[sp.Rational(1,d)]*d))],
            "unique trace-one Z4 invariant diagonal source lost")
    require(trace_last(mixture_density(d,3),d,3)
            ==mixture_density(d,2),
            "unique dephased projective family does not compose")
    return {
        "theorem":"for all n>=2: matched-code support + Z4 covariance + ordinary partial-trace projectivity implies rho_n=mu_n",
        "symbolic_fine_code_density_entries":d*d,
        "offdiagonal_source_coherences_removed_per_extension":d*d-d,
        "cyclic_diagonal_constraint_rank":3,
        "normalized_unique_diagonal_weight":"1/4",
        "unique_infinite_tower_source":"mu_n=(1/4) sum_j |j>^n<j|^n",
        "coherent_GHZ_clock_survives_local_projective_limit":False,
        "physical_R113_projective_map_is_ordinary_trace_proved":False,
        "physical_original_Padm_and_Q0_selected":False,
    }


def controls():
    return {
        "status":"EXACT_P1_GHZ_PROJECTIVE_SOURCE_NO_GO_AND_REPAIRS",
        "original_candidate":"matched boundary labels plus simultaneous Z4 clock constraint",
        "ordinary_coarse_graining":projective_state_gate(),
        "infinite_projective_uniqueness":infinite_matched_projective_uniqueness_gate(),
        "time_character_selection":cyclic_characters_gate(),
        "physical_cap_vs_categorical_gluing":cap_and_isometric_coarse_gate(),
        "source_overlap_tradeoff":overlap_and_clustering_gate(),
        "native_original_R113_full_Q0_available":False,
        "actual_Padm_and_absolute_Z0_evaluated":False,
        "physical_3plus1D_locality_or_E8_current_or_gravity_derived":False,
        "F1_solved":False,
    }


if __name__=="__main__":
    print(json.dumps(controls(),indent=2,sort_keys=True))
