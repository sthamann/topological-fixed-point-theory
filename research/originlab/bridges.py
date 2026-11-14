"""Executable typed bridge checks for the A3, charged P1 and regional proposals."""
from math import gcd
import sympy as s

from language import positive
from verify import r102

SIGMA = (
    s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1),
)
DIRECTIONS = tuple(s.Matrix(v)/s.sqrt(3) for v in
                   ((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)))
I = s.eye(2)


def spinor(vector):
    return sum((component*matrix for component,matrix in zip(vector,SIGMA)),s.zeros(2))


def exact_tetrahedron():
    f = [(I+spinor(n))/4 for n in DIRECTIONS]
    e = s.symbols("e0:4",real=True)
    energy = sum(e)/4
    momentum = sum((q*n for q,n in zip(e,DIRECTIONS)),s.zeros(3,1))*s.Rational(3,4)
    dual = [3*a-I/2 for a in f]
    reconstructed = sum((q*a for q,a in zip(e,dual)),s.zeros(2))
    if (reconstructed-energy*I-spinor(momentum)).applyfunc(s.simplify) != s.zeros(2):
        raise ArithmeticError("tetrahedral reconstruction failed")
    comm = f[0]*f[1]-f[1]*f[0]
    norm2 = s.simplify((comm.H*comm)[0,0])
    if s.simplify(16*(energy**2-momentum.dot(momentum))-
                  4*(sum(e)**2-3*sum(q*q for q in e))) != 0:
        raise ArithmeticError("unexpected cone normalization")
    p = -(DIRECTIONS[0]+DIRECTIONS[1])
    sampled = [s.simplify(1+sign*n.dot(p)) for n in DIRECTIONS for sign in (-1,1)]
    if min(sampled) != s.Rational(1,3) or s.simplify(p.dot(p)) != s.Rational(4,3):
        raise ArithmeticError("invalid finite direction counterexample")
    return {"sic_sum":str(sum(f,s.zeros(2))),"dual_eigenvalues":list(map(str,dual[0].eigenvals())),
            "commutator_norm":str(s.sqrt(norm2)),"finite_direction_min":str(min(sampled)),
            "finite_direction_momentum_squared":str(s.simplify(p.dot(p))),
            "joint_cone":"sum(e) >= 0 and 3*sum(e_a**2) <= sum(e)**2"}


def translations(p):
    if len(p) != 4 or any(m.shape != p[0].shape or m != m.H for m in p):
        raise ValueError("four Hermitian generators on one carrier are required")
    commutative = all(a*b==b*a for i,a in enumerate(p) for b in p[:i])
    block = s.kronecker_product(p[0],I)
    for operator,sigma in zip(p[1:],SIGMA):
        block += s.kronecker_product(operator,sigma)
    return {"commutative_finite":commutative,"positive_block":positive(block),
            "block_spectrum":{str(e):n for e,n in block.eigenvals().items()}}


def charged_relation(omega,psi,quarter,transfer):
    identity = s.eye(quarter.rows)
    if (quarter.H*quarter!=identity or quarter**4!=identity or quarter*omega!=omega or
        quarter*psi!=s.I*psi*quarter or transfer*quarter!=quarter*transfer or
        (omega.H*omega)[0]!=1 or not positive(transfer) or not positive(identity-transfer)):
        raise ValueError("quarter, charged field, transfer and reference do not satisfy the hypotheses")
    v = psi*omega
    ordinary = s.simplify((v.H*transfer*v)[0])
    marked = s.simplify((v.H*quarter*transfer*v)[0])
    return {"ordinary":str(ordinary),"marked":str(marked),
            "relative_phase":str(s.simplify(marked/ordinary)) if ordinary else None,
            "phase_relation_holds":marked==s.I*ordinary,"absolute_Z0_derived":False}


def tensor_control():
    h = s.eye(2)
    p = [h,*SIGMA]
    bad = translations(p)
    if bad != {"commutative_finite":False,"positive_block":False,"block_spectrum":{"2":3,"-2":1}}:
        raise ArithmeticError("Pauli witness failed")
    # Explicit same-source compatible finite controls also include an
    # operator-positive noncommutative case, so positivity alone cannot select.
    weak = translations([3*h,*SIGMA])
    diagonal = translations([2*h,s.diag(1,-1),s.zeros(2),s.zeros(2)])
    return {"rotational_Pauli":bad,"operator_positive_Pauli":weak,"joint_diagonal":diagonal}


def obligations():
    # Composition is allowed only when semantic output and requested input agree.
    return [
        {"from":"A3_weight_frame","to":"tetrahedral_SIC","status":"EXACT_SCOPED"},
        {"from":"tetrahedral_SIC","to":"commuting_regional_translations",
         "status":"DIRECT_IDENTIFICATION_REFUTED","reason":"nonzero commutators"},
        {"from":"strongly_commuting_translations+transitive_rotations+positive_null_generator",
         "to":"joint_forward_cone","status":"CONDITIONAL",
         "reason":"min_n (h+n·p)=h-|p| on the joint spectrum"},
        {"from":"continuous_connected_rotations","to":"permutation_of_four_fixed_labels",
         "status":"NONTRIVIAL_ACTION_REFUTED","reason":"continuous image of a connected group in a discrete group"},
        {"from":"primitive_local_charged_source+standard_pairing","to":"R102_E8",
         "status":"CONDITIONAL","reason":"monodromy order n/gcd(b,n) plus original R102 assumptions"},
        {"from":"complete_weighted_P1","to":"charged_relative_phase_and_regional_modular_inclusions",
         "status":"SOURCE_REQUIRED","reason":"same source, state, local products and absolute reference required"},
        {"from":"uniform_physical_gap_in_spin2_sector","to":"massless_spin2",
         "status":"INCOMPATIBLE_CONDITIONALLY","reason":"spectral decay bound on vacuum-orthogonal vectors"},
        {"from":"finite_word_separation","to":"complete_original_local_algebra_generation",
         "status":"SOURCE_REQUIRED","reason":"a decoupled spectator preserves all supplied restricted readouts"},
    ]


def run(max_rank=32):
    dual = r102(3,max_rank)["solutions"]
    phase = [{"n":n,"b":b,"phase_order":n//gcd(n,b),"primitive":gcd(n,b)==1}
             for n in range(3,max_rank+1) for b in range(1,n)]
    omega,psi = s.Matrix([1,0]),s.Matrix([[0,0],[1,0]])
    charge = charged_relation(omega,psi,s.diag(1,s.I),s.diag(1,s.Rational(2,3)))
    return {"status":"EXACT_BRIDGES_AND_CONDITIONAL_OBLIGATIONS",
            "tetrahedral":exact_tetrahedron(),"regional_controls":tensor_control(),
            "primitive_phase":phase,"r102_candidates":dual,"charged_finite_control":charge,
            "spinor_current":{"weights":["5/8","3/8"],"sum":"1","dimension":16*4,
                              "e8_branch_dimensions":[45,15,64,60,64]},
            "edges":obligations(),"physical_gate_passes":[]}
