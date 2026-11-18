"""Check the existing elliptic source, its graded double and correct reflections.

This reuses the old CAR and exterior-source implementation. The symbolic
parameter z is exp(-pi/12) in the old square annulus, not a fitted coupling.
Finite NS-mode checks supplement the all-mode bounded-operator proof.
"""
import importlib.util
import json
from pathlib import Path
import sys
import sympy as s

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
OLD = REPO / 'experiments/theory-contracts/seam-collar-multifield-20261016'
sys.path.insert(0, str(OLD))
from source_sewing import source_lift
from boundary_amplitude import graph_projection

spec = importlib.util.spec_from_file_location(
    'old_half_period', REPO / 'experiments/theory-contracts/'
    'native-arc-charge-source-interface-20261026/check_half_period_gns_sewing.py')
half = importlib.util.module_from_spec(spec)
spec.loader.exec_module(half)


def equal(a, b, name):
    diff = a-b
    entries = list(diff) if isinstance(diff, s.MatrixBase) else [diff]
    if any(s.factor(x) != 0 for x in entries):
        raise RuntimeError(name)


def vec(a):
    return s.Matrix([a[i,j] for j in range(a.cols) for i in range(a.rows)])


def run():
    z = s.symbols('z', positive=True)
    # Two actual positive NS modes lambda=pi,3pi, W=1/2, beta=1.
    # z=exp(-pi/12), so Q=exp(-W H), D=exp(-W H/3).
    T = s.diag(z**6, z**18)
    Q = source_lift(T)
    equal(Q, s.diag(1,z**6,z**18,z**24), 'old exterior source lift')
    D = s.diag(1,z**2,z**6,z**8)
    equal(D**3, Q, 'same source at the three collar substeps')
    partition = s.trace(Q**2)
    data = half.double_data([Q[i,i]**2 for i in range(4)])
    a, rho, root, parity, action, word, ops, adj, P = data
    equal(root, Q/s.sqrt(partition), 'same original thermal density')
    G = graph_projection(T)
    equal(P, G, 'actual graded source double equals old graph')
    C = s.diag(z**12/(1+z**12), z**36/(1+z**36))
    X = s.diag(z**6/(1+z**12), z**18/(1+z**36))
    equal(P, (s.eye(2)-C).row_join(X).col_join(X.row_join(C)),
          'all four labelled source covariance blocks')

    def left(A): return s.kronecker_product(s.eye(4), A)
    def right(A): return s.kronecker_product(A.T, s.eye(4))
    Gamma = left(parity)*right(parity)
    Jlinear = s.zeros(16)
    for i in range(4):
        for j in range(4): Jlinear[j+4*i,i+4*j] = 1
    def modular_conjugate(A): return Jlinear*s.conjugate(A)*Jlinear
    equal(Jlinear*s.conjugate(vec(root)), vec(root), 'J fixes actual source vector')
    aR = [left(A) for A in a]
    aL = [-Gamma*right(A) for A in a]
    for j in range(2):
        equal(modular_conjugate(aR[j]), Gamma*aL[j].H,
              'graded odd partner includes the parity operator')
        nR, nL = aR[j].H*aR[j], aL[j].H*aL[j]
        equal(modular_conjugate(nR), s.eye(16)-nL,
              'J sends occupation to partner hole, not bare face swap')
        equal((aR[j]-T[j,j]*aL[j])*vec(root), s.zeros(16,1),
              'existing first source sewing null relation')
        equal((aL[j].H+T[j,j]*aR[j].H)*vec(root), s.zeros(16,1),
              'existing adjoint sewing null relation')

    # Source-bounded even unitary insertions. These are not finite matrices
    # pretending to realize a nonzero central current CCR.
    U = source_lift(s.Matrix([[0,1],[-1,0]]))
    V = source_lift(s.diag(1,s.I))
    insertions = [s.eye(4), U, V, U*V]
    for A in insertions:
        equal(A.H*A, s.eye(4), 'bounded source insertion is unitary')
        equal(A*parity, parity*A, 'neutral even source insertion')
        equal(modular_conjugate(left(A)), right(A.H), 'J acts on bounded insertions')
    for A in insertions:
        for B in insertions:
            equal(right(A.H)*right(B.H), right((A*B).H),
                  'partner Weyls are multiplicative as a conjugate representation')
            equal(left(A)*right(B.H), right(B.H)*left(A), 'two neutral sides commute')
    equal(modular_conjugate(s.I*s.eye(16)), -s.I*s.eye(16),
          'J scalar action is conjugate linear')

    # Bounded half-cut pairing from this SAME rho; no freely chosen cross metric.
    neutral_gram = s.zeros(len(insertions))
    quarter = s.diag(1,z**3,z**9,z**12)/partition**s.Rational(1,4)
    for i,A in enumerate(insertions):
        for j,B in enumerate(insertions):
            joint = (vec(root).H*right(A.H)*left(B)*vec(root))[0]
            affinity = s.trace(root*A.H*root*B)
            hs = s.trace((quarter*A*quarter).H*(quarter*B*quarter))
            equal(joint, affinity, 'actual joint bounded source functional')
            equal(affinity, hs, 'half-cut cross response is a positive HS Gram')
            neutral_gram[i,j] = s.factor(joint)

    def thermal_word(events):
        ordered = sorted(events, reverse=True, key=lambda e:e[0])
        total = s.eye(4)
        previous = 6
        for t,A in ordered:
            total = total*(D**(previous-t))*A
            previous = t
        return s.trace(total*(D**previous))/partition

    histories = [(s.eye(4),s.eye(4)), (U,V), (V,U), (U*V,V*U)]
    history_gram = s.zeros(len(histories))
    for i,(A1,A2) in enumerate(histories):
        KF = D*A2*D*A1*D
        for j,(B1,B2) in enumerate(histories):
            KG = D*B2*D*B1*D
            # Reflection is t -> -t modulo beta, NOT s -> -s.
            source_value = thermal_word([(5,A1.H),(4,A2.H),(2,B2),(1,B1)])
            gram_value = s.trace(KF.H*KG)/partition
            equal(source_value, gram_value, 'ordered reflected thermal history')
            history_gram[i,j] = s.factor(source_value)

    charge_values = []
    for j,A in enumerate(a):
        exponent = (2,6)[j]
        psi = thermal_word([(5,A.H),(1,A)])
        bar = thermal_word([(5,A),(1,A.H)])
        equal(psi,z**(2*exponent)/(1+z**(6*exponent)), 'old positive reflected psi kernel')
        equal(bar,z**(4*exponent)/(1+z**(6*exponent)), 'old positive reflected bar kernel')
        midpoint = z**(3*exponent)/(1+z**(6*exponent))
        equal(midpoint,X[j,j], 'actual Green midpoint is the joint graph cross block')
        charge_values.append(dict(psi=str(s.factor(psi)),bar=str(s.factor(bar)),midpoint=str(midpoint)))

    # Neutral bilinear currents on the actual graded double, with the same
    # source subtraction, are the negative right action of source currents.
    currents = [s.eye(2),s.diag(1,-1),s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]])]
    for A in currents:
        N = sum((A[i,j]*a[i].H*a[j] for i in range(2) for j in range(2)),s.zeros(4))
        mean = s.trace(rho*N)
        partner = sum((A[i,j]*aL[i].H*aL[j] for i in range(2) for j in range(2)),s.zeros(16))
        partner -= (s.trace(A)-mean)*s.eye(16)
        equal(partner,-right(N-mean*s.eye(4)), 'opposite current orientation is derived')
        equal(modular_conjugate(left(N-mean*s.eye(4))),-partner,
              'modular conjugation reverses centered partner current')

    # The actual old product-collar reflection identities, on a +/- mode pair.
    paired_T = s.diag(z**6,z**-6)
    paired_P = graph_projection(paired_T)
    opposite_P = graph_projection(paired_T.inv())
    R = s.Matrix([[0,1],[1,0]])
    face_swap = s.zeros(2).row_join(s.eye(2)).col_join(s.eye(2).row_join(s.zeros(2)))
    fibre_flip = s.diag(R,R)
    equal(face_swap*paired_P*face_swap,opposite_P,'bare face swap reverses B')
    equal(fibre_flip*paired_P*fibre_flip,opposite_P,'deck fibre action reverses B')
    equal(face_swap*fibre_flip*paired_P*fibre_flip*face_swap,paired_P,
          'combined inner collar reflection preserves the full graph')

    return dict(
        verdict='PASS_EXISTING_TWO_FACE_SOURCE_AND_TYPED_REFLECTION',
        actual_symbolic_source='z=exp(-pi/12), W=1/2, beta=1, NS lambda=(pi,3pi)',
        original_implementations_reused=['source_sewing.source_lift','boundary_amplitude.graph_projection','check_half_period_gns_sewing.double_data'],
        covariance=[[str(s.factor(x)) for x in row] for row in P.tolist()],
        cross=[str(X[i,i]) for i in range(2)],
        reflected_charge_kernels=charge_values,
        neutral_gram=[[str(x) for x in row] for row in neutral_gram.tolist()],
        reflected_history_gram=[[str(x) for x in row] for row in history_gram.tolist()],
        J_is_conjugate_linear_multiplicative=True,
        J_occupation_image='J n_R J = I - n_L',
        partner_current_orientation='j_L(A) = -R_(j(A)); opposite Weyl cocycle follows',
        correct_OS_coordinate_map='rho(t,s)=(-t,s); no spatial profile reversal',
        bounded_Weyl_all_mode_scope='Analytic bounded-operator HS theorem in PROOF.txt; not inferred from finite matrices',
        actual_physical_P1_second_boundary_identified=False,
        primitive_P1_measure_or_source_coupling_derived=False,
        physical_time_identified=False,
        complete_TFPT_solution=False,ledger_changed=False)


if __name__=='__main__':
    result=run()
    Path(__file__).with_name('two-face-source-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('covariance','neutral_gram','reflected_history_gram')},indent=2))
