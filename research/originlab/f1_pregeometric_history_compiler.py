"""Source-first F1 candidate: all admissible anchored four-cycles on the
original five-point quadratic code, before time or E8 is imposed.

This is a NEW exact combinatorial/compositional mathematical candidate.
The physical ORIGINAL P1 seam anchor and the coherent amplitude rule
are NOT sourced from Q0; no 3+1D theory, physical OS time or F1 is proved.
"""
from itertools import permutations, product
import sympy as sp
from marked_search import phase_and_clock


def require(ok,why):
    if not ok:
        raise ArithmeticError(why)


def permutation_matrix(p):
    return sp.eye(len(p))[:,list(p)]


def anchored_clocks(anchor_index=4):
    clocks=[p for p in permutations(range(5))
            if sum(int(i==p[i]) for i in range(5))==1
            and all(p[p[p[p[i]]]]==i for i in range(5))
            and any(p[p[i]]!=i for i in range(5))]
    require(len(clocks)==30,"original five-point code has not 30 clocks")
    require(0<=anchor_index<5, "fifth code point anchor index out of range")
    anchored=[p for p in clocks if p[anchor_index]==anchor_index]
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


def d4_weighted_clock_gate(cycles):
    """Physically square D4 geometry has TWO inequivalent 4-cycle orbits.

    A single source parameter x allocates total weight x to the two
    geometric quarter rotations r,r^-1 and 1-x to the four remaining
    primitive 4-cycles. No native Q0 supplied this parameter.
    """
    r=tuple((i+1)%4 for i in range(4))
    s=tuple((-i)%4 for i in range(4))
    multiply=lambda p,q:tuple(p[q[i]] for i in range(4))
    group={tuple(range(4))}
    discovered=list(group)
    for p in discovered:
        for g in (r,s):
            new=multiply(p,g)
            if new not in group:
                group.add(new)
                discovered.append(new)
    require(len(group)==8, "square D4 orientation group changed")
    inverse=lambda p:tuple(p.index(i) for i in range(4))
    conjugate=lambda g,p:multiply(multiply(g,p),inverse(g))
    unclassified=set(cycles)
    orbits=[]
    while unclassified:
        c=next(iter(unclassified))
        orbit={conjugate(g,c) for g in group}
        require(orbit<=set(cycles), "D4 rotation class not stable")
        orbits.append(orbit)
        unclassified-=orbit
    require(sorted(len(orbit) for orbit in orbits)==[2,4],
            "geometric D4 should split six S4 related 4-cycles")
    shorter=next(orbit for orbit in orbits if len(orbit)==2)
    longer=next(orbit for orbit in orbits if len(orbit)==4)
    require(shorter=={r,inverse(r)}, "geometric two-cycle orbit not actual quarter rotations")
    X=sum((permutation_matrix(p) for p in shorter),sp.zeros(4))/2
    Y=sum((permutation_matrix(p) for p in longer),sp.zeros(4))/4
    x=sp.symbols("x",real=True)
    T=sp.simplify(x*X+(1-x)*Y)
    P0=sp.ones(4)/4
    v=sp.Matrix([1,-1,1,-1])
    Palt=v*v.T/4
    Pdoublet=sp.eye(4)-P0-Palt
    require(T==P0-x*Palt-(1-x)*Pdoublet/2
            and T.eigenvals()=={1:1,-x:1,(x-1)/2:2},
            "two-rate D4 clock family spectrum or irreducible decomposition changed")
    assert Palt.rank()==1 and Pdoublet.rank()==2
    squared=T*T
    require(sp.simplify(squared-P0-x*x*Palt-
                        (1-x)**2*Pdoublet/4)==sp.zeros(4),
            "D4 two-step clock is not a positive two-rate transfer")
    H=-2*sp.log(x)*Palt-2*sp.log((1-x)/2)*Pdoublet
    require(sp.simplify(H.subs(x,sp.Rational(1,3))-
                        sp.log(9)*(sp.eye(4)-P0))==sp.zeros(4),
            "S4 enhancement at equal six histories not recovered")
    for value in (sp.Rational(1,4),sp.Rational(1,2)):
        chosen=T.subs(x,value)
        require(chosen*sp.ones(4,1)==sp.ones(4,1) and
                all(z>=0 for z in chosen) and
                (chosen*chosen).is_positive_definite,
                "D4 original square symmetric countermodel does not give positive 2-step time")
    return {
        "status":"EXACT_D4_TWO_ORBIT_WEIGHTS_NOT_SELECTED_BY_SOURCE",
        "original_geometric_square_D4_group_order":8,
        "primitive_four_cycles_class_orbit_sizes":[2,4],
        "independent_source_weight_parameter":"x in (0,1)",
        "square_symmetry_preserving_one_step":"x*(r+r^-1)/2+(1-x)*(other_4_cycles)/4",
        "one_step_eigenvalues":["1","-x","(x-1)/2 (twice)"],
        "two_step_eigenvalues":["1","x^2","(1-x)^2/4 (twice)"],
        "two_step_positive_H_generator":
            "H=-2log(x)P_alternating-2log((1-x)/2)P_doublet",
        "enhanced_S4_equal_six_cycles_source_weight":"x=1/3",
        "enhanced_S4_pair_rates_equal_at_x_one_third":True,
        "distinct_D4_source_witnesses":["x=1/4","x=1/2"],
        "D4_only_fixes_unique_source_time":False,
        "physical_original_absolute_weight_x_measured":False,
        "full_S4_symmetry_before_geometric_P1_seam_derived":False
    }


def coherent_fourth_root_phase_search(cycles):
    """Enumerate ALL 4^5 relative quarter-phase assignments over the six
    anchored 4-cycles, with the first path's phase set to one.

    Test Gaussian integer column orthogonality EXACTLY before building
    any SymPy spectral matrices. This is not the original q0 phase rule.
    """
    phase=(1,-1,1j,-1j)
    retained=[]
    for relative in product(phase,repeat=5):
        weights=(1,*relative)
        amp=[[0j for _ in range(4)] for _ in range(4)]
        for p,w in zip(cycles,weights):
            for col,row in enumerate(p):
                amp[row][col]+=w
        orthogonal=all(sum(amp[k][i].conjugate()*amp[k][j]
                      for k in range(4)) == (4 if i==j else 0)
                      for i in range(4) for j in range(4))
        if orthogonal:
            A=sp.Matrix(4,4,lambda i,j:
                sp.Integer(int(amp[i][j].real))+
                sp.I*sp.Integer(int(amp[i][j].imag)))
            U=A/2
            require(U.H*U==sp.eye(4) and U**4==sp.eye(4),
                    "exact Gaussian phase solution is not unitary order <=4")
            retained.append((weights,U))
    require(len(retained)==9 and
            sum(all(w.imag==0 for w in weights)
                for weights,_ in retained)==3 and
            sum(U*U==sp.eye(4) for _,U in retained)==3,
            "anchored coherent fourth-root solutions changed")
    expected_even=(1,-1,1j,-1j)
    spectra=[{str(v):int(n) for v,n in sorted(U.eigenvals().items(),key=lambda z:str(z[0]))}
             for _,U in retained]
    require(all(U!=sp.eye(4) for _,U in retained) and
            all(U**4==sp.eye(4) for _,U in retained),
            "coherent-phase clock becomes an irreversible positive semigroup")
    return {
      "gauge_fixed_phase_conventions":"first of six 4-cycle amplitudes=1",
      "relative_phase_alphabet":["1","-1","i","-i"],
      "relative_amplitude_assignments_exhaustively_checked":4**5,
      "exact_unitary_amplitude_sum_over_two_solutions":len(retained),
      "involutive_real_solutions":3,
      "genuine_order_four_complex_solutions":6,
      "all_coherent_solutions_unitary":True,
      "all_coherent_solutions_have_U_fourth_identity":True,
      "coherent_solutions_strictly_positive_time_transfer":False,
      "coherent_phase_choice_unique":False,
      "physical_original_quadratic_refinement_selects_one_phase":False,
      "unitary_coherent_spin_clock_automatically_gives_dissipative_OS_time":False,
      "solutions": [
        {"phases":[str(w) for w in weights],
         "U":[[str(v) for v in row] for row in U.tolist()],
         "spectrum":spectra[i],
         "order":2 if U**2==sp.eye(4) else 4}
        for i,(weights,U) in enumerate(retained)],
    }


def controls():
    original=phase_and_clock({})
    base=original["clock"]
    five_q0=base["point_vectors"]
    require(len(five_q0)==5 and base["initial_clocks"]==30
            and len({tuple(p) for p in five_q0})==5,
            "originlab five-point odd quadratic code has changed")
    # The original face-label S4 permutation action fixes exactly ONE
    # vector among the five odd-quadratic null vectors, namely 1111.
    # This identifies a CANONICAL fifth code point, not yet a physical P1.
    face_group=list(permutations(range(4)))
    invariant=[i for i,v in enumerate(five_q0) if all(
        [v[g[k]] for k in range(4)]==v for g in face_group)]
    require(invariant==[4] and five_q0[4]==[1,1,1,1] and
            all([v[g[k]] for k in range(4)] in five_q0
                for v in five_q0 for g in face_group),
            "original S4 odd-quadratic code does not intrinsically anchor 1111")
    clocks,anchored=anchored_clocks(invariant[0])
    pall=sum((permutation_matrix(p) for p in clocks),sp.zeros(5))/30
    require(pall==sp.ones(5)/5 and pall.rank()==1,
            "unanchored source clock sum did not erase nontrivial sector")

    pab=sum((permutation_matrix(p) for p in anchored),sp.zeros(5))/6
    require(pab[4,4]==1 and
            pab[4,0:4]==sp.zeros(1,4) and
            pab[0:4,4]==sp.zeros(4,1),
            "S4-invariant fifth point is not fixed by the six clocks")
    Q=pab[0:4,0:4]
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

    d4_weights=d4_weighted_clock_gate(perms4)
    coherent=coherent_fourth_root_phase_search(perms4)
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
      "S4_invariant_fifth_code_point_index":4,
      "S4_invariant_fifth_code_point_bits":[1,1,1,1],
      "original_S4_action_fixes_fifth_code_point_uniquely":True,
      "physical_P1_identity_of_fifth_code_point_proved":False,
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
      "coherent_exact_relative_phase_search":coherent,
      "geometric_D4_vs_abstract_S4_weight_selection":d4_weights,
      "coherent_average_vs_CP_history_mixture_physically_resolved":False,
      "original_quadratic_refinement_spin_phase_selected":False,
      "physical_R113_Q0_time_and_reference_derived":False,
      "actual_pregeometric_3plus1D_gravity_and_E8_derived":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
