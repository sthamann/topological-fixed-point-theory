"""Conditional original-source candidate from exchangeable Klein-twisted C* chains.

External theorem (Fidaleo/Vincenzi, arXiv:2609.09494, Sep 2026):
under the specific assumptions of Hermitian nondegenerate bicharacter
and trivial restriction on Delta_+ x Delta_+, the tractable
de-Finetti twisted families reduce to ordinary, Fermi and
Klein V4 cases up to equivalence.

This independent exact audit checks the KLEIN instance and all four
mu4-valued quadratic refinements, which include double semion phases,
while strictly avoiding an identification with the original alternating
TFPT F2^4 symplectic form or native physical R113 Q0.
"""
import json
from itertools import product

import sympy as sp


def require(ok,why):
    if not ok:
        raise ArithmeticError(why)


G=list(product((0,1),repeat=2))


def add(x,y):
    return ((x[0]+y[0])%2,(x[1]+y[1])%2)


def klein(x,y):
    return (-1)**(x[0]*y[0]+x[1]*y[1])


def alternating(x,y):
    return (-1)**(x[0]*y[1]+x[1]*y[0])


def bicharacter_gate():
    B=sp.Matrix([[klein(x,y) for y in G] for x in G])
    A=sp.Matrix([[alternating(x,y) for y in G] for x in G])
    for f in (klein,alternating):
        require(all(f(add(x,y),z)==f(x,z)*f(y,z)
                    and f(x,add(y,z))==f(x,y)*f(x,z)
                    for x,y,z in product(G,repeat=3)),
                "G bicharacter multiplicativity failed")
        require(all(f(x,y)==f(y,x) for x,y in product(G,repeat=2)),
                "G bicharacter Hermitian property changed")
    require(B*B.T==4*sp.eye(4) and A*A.T==4*sp.eye(4),
            "both G pairings must be nondegenerate orthogonal character matrices")
    delta_K=[x for x in G if klein(x,x)==1]
    delta_A=[x for x in G if alternating(x,x)==1]
    require(delta_K==[(0,0),(1,1)] and delta_A==G
            and all(klein(x,y)==1 for x,y in product(delta_K,repeat=2))
            and any(alternating(x,y)!=1 for x,y in product(delta_A,repeat=2)),
            "Klein isotropic de Finetti condition changed")
    # The native original symplectic F2^4 pairing is ALTERNATING.
    # The Klein pairing is nonalternating (self-pairing -1 on
    # e1/e2), so NO bilinear-form-preserving embedding can
    # directly identify the two.
    require(klein((1,0),(1,0))==-1
            and alternating((1,0),(1,0))==1,
            "native alternating-vs-Klein form mismatch changed")
    # The Z4 nondegenerate unitary bicharacters contain +/-i,
    # violating the Hermitian flip condition; four sectors alone
    # do not authorize an underlying cyclic Z4 exchange algebra.
    for sign in (1,-1):
        require(sp.I**sign != sp.conjugate(sp.I**sign),
                "Z4 fundamental quarter pairing unexpectedly Hermitian")
    return {
      "G":"Z2 x Z2",
      "unique_nontrivial_exchangeable_candidate_under_external_axioms":
        "Klein bicharacter modulo equivalent forms",
      "klein_pairing":"(-1)^(x1*y1+x2*y2)",
      "klein_self_isotropic_subgroup":["00","11"],
      "isotropic_order":2,
      "klein_u_on_isotropic_subgroup_trivial":True,
      "alternating_symplectic_u_on_its_isotropic_subgroup_trivial":False,
      "fourier_char_matrix_orthogonality_factor":4,
      "cyclic_Z4_primitive_quarter_bicharacter_Hermitian":False,
      "klein_directly_isometric_to_original_alternating_TFPT_pairing":False
    }


def quadratic_spin_gate():
    eps=(sp.I,-sp.I)
    solutions=[]
    for a,b in product(eps,repeat=2):
        q={(0,0):sp.Integer(1),(1,0):a,(0,1):b,(1,1):a*b}
        require(all(sp.simplify(q[add(x,y)]-q[x]*q[y]*klein(x,y))==0
                    for x,y in product(G,repeat=2)),
                "one Klein quadratic spin refinement fails polarizing identity")
        central=sp.simplify(sum(q.values())/2)
        require(central in (1,sp.I,-sp.I),
                "Klein spin Gauss-Milgram root changed")
        solutions.append({
          "q10":str(a),"q01":str(b),"q11":str(a*b),
          "spins":[str(q[x]) for x in [(0,0),(1,0),(0,1),(1,1)]],
          "Gauss_Milgram":str(central),
          "chiral_c_mod8":{1:0,sp.I:2,-sp.I:6}[central],
        })
    require(len(solutions)==4
            and sorted(s["chiral_c_mod8"] for s in solutions)==[0,0,2,6],
            "four Klein refinements and c mod8 signature changed")
    ds=[s for s in solutions if s["chiral_c_mod8"]==0]
    require(len(ds)==2 and {tuple(s["spins"]) for s in ds}
            =={("1","I","-I","1"),("1","-I","I","1")},
            "two opposite double-semion orientation refinements changed")
    return {
      "q_twisted_quadratic_refinement_count":4,
      "double_semion_c0_mod8_refinements":2,
      "reflection_exchanges_two_double_semion_orientations":True,
      "other_chiral_c_mod8_classes":[2,6],
      "quadratic_refinements":solutions,
      "native_P1_orienting_spin_character_derived":False,
      "absolute_E8_stacking_number_fixed":False,
    }


def controls():
    return {
      "status":"EXACT_KLEIN_INFINITE_SOURCE_CANDIDATE_WITH_QUARTER_SPIN_AMBIGUITY",
      "external_classification_theorem_input":
        "Fidaleo/Vincenzi 2026, Hermitian, nondegenerate and isotropy-trivial",
      "bicharacter":bicharacter_gate(),
      "quadratic_spin":quadratic_spin_gate(),
      "actual_original_R113_Padm_Q0_time_and_framing_selected":False,
      "infinite_twisted_Cstar_algebra_embedded_in_original_TFPT":False,
      "physical_3plus1D_locality_or_spin2_derived":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    print(json.dumps(controls(),indent=2,sort_keys=True))
