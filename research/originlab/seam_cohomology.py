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
