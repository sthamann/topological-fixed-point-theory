"""Exact equivariant four-mark cohomology and minimal orientation compensator.

Checks the original D=mu4 Riemann--Roch carrier (5 sections) against the
logarithmic differential family space (3 forms) and the square's C4/D4
actions. This is a finite necessary-origin control, never an original
weighted R113 source or a physical SU(3)xSU(2) representation theorem.
"""
from itertools import combinations_with_replacement

import sympy as sp


def _require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def _commutant_dimension(*actions):
    n = actions[0].rows
    values = sp.symbols(f"x0:{n*n}")
    unknown = sp.Matrix(n,n,values)
    equations = [item for action in actions
                 for item in (unknown*action-action*unknown)]
    coeff,_ = sp.linear_eq_to_matrix(equations,values)
    return n*n-coeff.rank()


def controls():
    one5 = sp.eye(5)
    r5 = sp.diag(1,sp.I,-1,-sp.I,1)
    s5 = -sp.Matrix(5,5,lambda i,j:int(i+j==4))
    c = sp.Matrix([-1,0,0,0,1])  # f4-f0 == 1
    _require(r5**4 == one5 and s5**2==one5 and
             s5*r5*s5==r5.inv() and r5*c==c and s5*c==c,
             "five-mode mu4 pullback and invariant constant are inconsistent")

    # f_(j+1)=z^(j+1)/(z^4-1), j=0,1,2.
    rc = sp.diag(sp.I,-1,-sp.I)
    sc = -sp.Matrix(3,3,lambda i,j:int(i+j==2))
    # omega_j=z^j dz/(z^4-1), j=0,1,2.
    rf = sp.diag(sp.I,-1,-sp.I)
    sf = sp.Matrix(3,3,lambda i,j:int(i+j==2))
    _require(rc==rf and sc==-sf and all(
        x**4==sp.eye(3) and y**2==sp.eye(3) and
        y*x*y==x.inv() for x,y in ((rc,sc),(rf,sf))),
        "the intrinsic reflection orientation line is missing")

    # The two fixed points of r on P1 are z=0 and z=infinity.
    # Evaluation on O(D) sends f0 -> (-1,0), f4 -> (0,1), and
    # annihilates f1,f2,f3. Multiplication by dz/z identifies
    # f_(j+1) with the original logarithmic omega_j.
    endpoint = sp.Matrix([[-1,0,0,0,0],[0,0,0,0,1]])
    inclusion = sp.Matrix(5,3,lambda i,j:int(i==j+1))
    endpoint_rotation = sp.eye(2)
    endpoint_reflection = sp.Matrix([[0,1],[1,0]])
    _require(endpoint.rank()==2 and inclusion.rank()==3 and
             endpoint*inclusion==sp.zeros(2,3) and
             r5*inclusion==inclusion*rf and
             s5*inclusion==-inclusion*sf and
             endpoint*r5==endpoint_rotation*endpoint and
             endpoint*s5==endpoint_reflection*endpoint,
             "exact D4-equivariant 3+2 fixed-point evaluation sequence fails")

    # The full D4 commutant on the five original geometric modes has
    # FOUR primitive central projections of ranks 1,1,1,2.
    # An invariant positive Hamiltonian with the constant vacuum fixed
    # therefore retains THREE independent excited energies.
    standard = [sp.eye(5)[:,k] for k in range(5)]
    zero = (standard[4]-standard[0])/sp.sqrt(2)
    weak_odd = (standard[4]+standard[0])/sp.sqrt(2)
    projections = (
        zero*zero.H,
        weak_odd*weak_odd.H,
        standard[2]*standard[2].H,
        standard[1]*standard[1].H+standard[3]*standard[3].H
    )
    _require(sum(projections,sp.zeros(5))==one5 and
             [p.rank() for p in projections]==[1,1,1,2] and
             all(p*r5==r5*p and p*s5==s5*p for p in projections),
             "D4 geometric Hamiltonian central projections changed")
    p0,pw,pm,pp = projections
    h_a = pw+2*pm+3*pp
    h_b = pw+4*pm+2*pp
    _require(h_a*zero==sp.zeros(5,1) and h_b*zero==sp.zeros(5,1) and
             h_a.eigenvals()=={0:1,1:1,2:1,3:2} and
             h_b.eigenvals()=={0:1,1:1,2:2,4:1} and
             sp.trace(h_a)==sp.trace(h_b)==9,
             "same-gap same-trace positive source witnesses changed")
    t_a = p0+pw/2+pm/4+pp/8
    t_b = p0+pw/2+pm/16+pp/4
    selected = standard[2]
    returns = [(selected.H*t*selected)[0] for t in (t_a,t_b)]
    _require(returns==[sp.Rational(1,4),sp.Rational(1,16)],
             "geometric state does not distinguish source times")

    # The SAME canonical differential eta=dz/z carries primitive
    # unit winding, is quarter-rotation invariant and reflection odd.
    z = sp.symbols("z")
    winding = sp.residue(1/z,z,0)
    _require(winding==1 and
             sp.simplify(sp.I/(sp.I*z)-1/z)==0 and
             sp.simplify((-1/z**2)/(1/z)+1/z)==0,
             "unit dlog winding and orientation response changed")
    _require(sp.simplify(1/(4*2*sp.pi)-1/(8*sp.pi))==0,
             "index-four angular normalization arithmetic changed")

    # C3 x F3^*, pairing characters +j with -j cancels clock charge.
    rt = sp.kronecker_product(rc,rf.inv().T)
    st = sp.kronecker_product(sc,sf.inv().T)
    ident = sp.eye(9)
    cyclic = (ident+rt+rt**2+rt**3)/4
    full = cyclic*(ident+st)/2
    twisted = cyclic*(ident-st)/2
    _require(all(p==p.H and p**2==p for p in (cyclic,full,twisted)),
             "averages are not orthogonal projectors")
    _require([p.rank() for p in (cyclic,full,twisted)]==[3,1,2],
             "orientation projection ranks changed")
    _require(rt*cyclic==cyclic, "cyclic compensated channels carry charge")

    # Every one of the 3 distinct nontrivial characters needs its inverse
    # character in an ancilla; any lower-dimensional ancilla misses one.
    minimum = next(n for n in range(4) if any(
        all((-j)%4 in chars for j in (1,2,3))
        for chars in combinations_with_replacement(range(4),n)))
    _require(minimum==3, "minimal compensator dimension changed")

    # The unpaired 3 modes have only diagonal clock-commuting operators:
    # no full SU3 can commute with the literal geometric quarter action.
    _require(_commutant_dimension(rc)==3 and
             _commutant_dimension(rc,sc)==2 and
             _commutant_dimension(r5)==7 and
             _commutant_dimension(r5,s5)==4,
             "geometric-clock commutant changed")

    # The compensated 3 channels, together with the 2 trivial-character
    # slots, can support U3 x U2 *if* the original source preserves
    # their distinct geometric types. Unimodularity gives the familiar
    # hypercharge diagonal, but is NOT derived here.
    y = sp.diag(*([sp.Rational(-1,3)]*3+[sp.Rational(1,2)]*2))
    _require(sp.trace(y)==0 and sp.trace(rt*cyclic)==3,
             "neutral compensated colour rank or hypercharge sum changed")

    return {
       "status":"EXACT_FINITE_SEAM_COHOMOLOGY",
       "source_theorems":["v189","v197","v228","v512","architecture shared c3 P1 P2 origin"],
       "dlog_winding_around_zero":"1 (geometric, not an original physical P1 response)",
       "dlog_orientation_under_inversion":"-1",
       "candidate_index_times_angle_normalization":"1/(8*pi) (physical identification unproved)",
       "fixed_point_evaluation_rank":2,
       "fixed_point_evaluation_kernel_dimension":3,
       "D4_orientation_twisted_three_plus_two_exact_sequence":True,
       "two_fixed_points":"0 and infinity (original quarter rotation)",
       "D4_irreducible_carrier_dimensions":[1,1,1,2],
       "D4_invariant_generator_excited_free_parameters":3,
       "matched_gap_and_trace_witness_spectra":[["0","1","2","3","3"],["0","1","2","2","4"]],
       "same_gap_trace_distinguishing_return_at_log2":["1/4","1/16"],
       "five_mode_c4_multiplicities":{"1":2,"i":1,"-1":1,"-i":1},
       "log_family_c4_multiplicities":{"1":0,"i":1,"-1":1,"-i":1},
       "functorial_orientation_twist":"C3 == F3 tensor reflection_sign on D4",
       "cyclic_neutral_compensated_rank":int(cyclic.rank()),
       "full_dihedral_invariant_rank":int(full.rank()),
       "orientation_twisted_full_invariant_rank":int(twisted.rank()),
       "minimal_distinct_character_compensator_dimension":minimum,
       "uncompensated_colour_commutant_complex_dimension":_commutant_dimension(rc),
       "uncompensated_colour_full_d4_commutant_complex_dimension":_commutant_dimension(rc,sc),
       "carrier_quarter_commutant_complex_dimension":_commutant_dimension(r5),
       "hypercharge_trace":"0 (candidate, unimodularity imposed)",
       "source_orientation_and_eta_identified":False,
       "original_weighted_Q0_produced":False,
       "original_R113_values_produced":False,
       "physical_SU3_SU2_identification_proved":False,
       "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
