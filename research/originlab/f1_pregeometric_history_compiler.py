"""Source-first F1 candidate: all admissible anchored four-cycles on the
original five-point quadratic code, before time or E8 is imposed.

This is a NEW exact combinatorial/compositional mathematical candidate.
The physical ORIGINAL P1 seam anchor and the coherent amplitude rule
are NOT sourced from Q0; no 3+1D theory, physical OS time or F1 is proved.
"""
from itertools import permutations
import sympy as sp
from marked_search import phase_and_clock


def require(ok,why):
    if not ok:
        raise ArithmeticError(why)


def permutation_matrix(p):
    return sp.eye(len(p))[:,list(p)]


def anchored_clocks():
    clocks=[p for p in permutations(range(5))
            if sum(int(i==p[i]) for i in range(5))==1
            and all(p[p[p[p[i]]]]==i for i in range(5))
            and any(p[p[i]]!=i for i in range(5))]
    require(len(clocks)==30,"original five-point code has not 30 clocks")
    anchored=[p for p in clocks if p[0]==0]
    require(len(anchored)==6,"distinguished anchor no longer selects 6 cycles")
    return clocks,anchored


def quantum_random_unitary_channel(permutations4):
    """Matrix unit basis: row-major 16-dimensional Heisenberg transfer."""
    matrices=[permutation_matrix(p) for p in permutations4]
    e=[sp.eye(4)[:,i]*sp.eye(4)[j,:]
       for i in range(4) for j in range(4)]
    superop=sp.Matrix.hstack(*[
        (sum((m*v*m.T for m in matrices),sp.zeros(4))/6).reshape(16,1)
        for v in e])
    return superop


def controls():
    original=phase_and_clock({})
    base=original["clock"]
    five_q0=base["point_vectors"]
    require(len(five_q0)==5 and base["initial_clocks"]==30
            and len({tuple(p) for p in five_q0})==5,
            "originlab five-point odd quadratic code has changed")
    clocks,anchored=anchored_clocks()
    pall=sum((permutation_matrix(p) for p in clocks),sp.zeros(5))/30
    require(pall==sp.ones(5)/5 and pall.rank()==1,
            "unanchored source clock sum did not erase nontrivial sector")

    pab=sum((permutation_matrix(p) for p in anchored),sp.zeros(5))/6
    require(pab[0,0]==1 and
            pab[0,1:5]==sp.zeros(1,4) and
            pab[1:5,0]==sp.zeros(4,1),
            "source anchor is not a fixed seam point")
    Q=pab[1:5,1:5]
    I=sp.eye(4)
    J=sp.ones(4)
    P0=J/4
    require(Q==(J-I)/3 and Q.eigenvals()=={1:1,sp.Rational(-1,3):3},
            "anchored source rotation class sum altered")
    # Invariance under every relabelling of the four active points
    # proves any conjugacy-invariant positive weights on the six
    # primitive 4-cycles must be equal.
    perms4=[p for p in permutations(range(4))
            if all(p[p[p[p[i]]]]==i for i in range(4))
            and any(p[p[i]]!=i for i in range(4))]
    require(len(perms4)==6 and
            sum((permutation_matrix(p) for p in perms4),sp.zeros(4))/6==Q,
            "six anchored active rotations do not cover one S4 orbit")
    for p in permutations(range(4)):
        W=permutation_matrix(p)
        require(W*Q==Q*W,"clock class average fails four-mark covariance")
        def conjugate(q):
            return tuple(p[q[p.index(i)]] for i in range(4))
        require(set(conjugate(q) for q in perms4)==set(perms4),
                "six clocks cease to be one conjugacy class")

    B=Q*Q
    require(B==P0+(I-P0)/9 and B==B.T and
            B.eigenvals()=={1:1,sp.Rational(1,9):3} and
            B.is_positive_definite and B*sp.ones(4,1)==sp.ones(4,1),
            "paired-clock source not positive or exactly sewn")
    u=sp.symbols("u",nonnegative=True)
    semigroup=P0+sp.exp(-sp.log(9)*u)*(I-P0)
    require(semigroup.subs(u,0)==I and
            semigroup.subs(u,1)==B and
            sp.simplify(semigroup*B-(P0+sp.exp(-sp.log(9)*(u+1))*(I-P0)))==sp.zeros(4),
            "paired clock lacks nonnegative strongly continuous interpolation")
    H=sp.log(9)*(I-P0)
    require(sp.simplify(sp.exp(-H)-B)==sp.zeros(4),
            "canonical H=-log B not the whole one-source time generator")
    # A source-compatible positive energy on the active three-dimensional
    # quotient: E(f)=<f,H f>_pi, with pi uniform on four labels.
    x=sp.Matrix([1,-1,0,0])
    energy=sp.simplify((x.T*H*x)[0]/4)
    require(energy==sp.log(9)/2,
            "symmetric 3-mode quantum Dirichlet form changed")

    # Crucial COMPILER difference: incoherent average of conjugations
    # on the full M4 quantum operator algebra is a perfectly CP map,
    # but its SQUARE has zero eigenvalues, so it cannot be e^-L for
    # finite bounded generator L at a strictly positive time.
    M=quantum_random_unitary_channel(perms4)
    require(M==M.T and
            M.eigenvals()=={1:2,sp.Rational(-1,3):9,
                            sp.Rational(1,3):3,0:2} and
            (M*M).eigenvals()=={1:2,sp.Rational(1,9):12,0:2}
            and M.det()==0,
            "random-unitary mixture unexpectedly equals coherent clock")
    return {
      "status":"EXACT_ANCHORED_FIVE_POINT_HISTORY_COMPILER_CANDIDATE",
      "original_five_q0_points":five_q0,
      "all_unanchored_quarter_clocks":30,
      "uniform_unanchored_sum":"J5/5 (rank one, all nontrivial modes lost)",
      "physical_seam_anchor_selected_by_original_Q0":False,
      "one_distinguished_anchor_compatible_clocks":6,
      "anchored_one_step":"(J4-I4)/3",
      "one_step_spectrum":{"1":1,"-1/3":3},
      "minimum_even_positive_return_steps":2,
      "anchored_two_step":"P0+(I4-P0)/9",
      "two_step_spectrum":{"1":1,"1/9":3},
      "paired_two_step_OS_semigroup_compatible":True,
      "finite_source_generator":"H=log(9)*(I4-J4/4)",
      "positive_three_mode_energy_probe":"log(9)/2",
      "source_clock_unit_physically_derived":False,
      "natural_three_dimensional_nontrivial_standard_sector":True,
      "three_dimensional_sector_is_three_physical_SM_families":False,
      "quantum_random_unitary_one_step_eigenspaces":{"1":2,"-1/3":9,"1/3":3,"0":2},
      "quantum_random_unitary_two_step_eigenspaces":{"1":2,"1/9":12,"0":2},
      "incoherent_quantum_average_has_bounded_continuous_generator":False,
      "coherent_average_vs_CP_history_mixture_physically_resolved":False,
      "original_quadratic_refinement_spin_phase_selected":False,
      "physical_R113_Q0_time_and_reference_derived":False,
      "actual_pregeometric_3plus1D_gravity_and_E8_derived":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
