"""Independent finite CAR checks of the analytic projector-state theorem."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as S


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def equal(a, b, message):
    if isinstance(a, S.MatrixBase):
        require(all(S.simplify(x) == 0 for x in a-b), message)
    else:
        require(S.simplify(a-b) == 0, message)


def run():
    r = S.symbols('r', positive=True)
    I = S.eye(4)
    # Basis |00>,|10>,|01>,|11>. a1 carries the occupation0 parity.
    a0 = S.zeros(4); a0[0, 1] = 1; a0[2, 3] = 1
    a1 = S.zeros(4); a1[0, 2] = 1; a1[1, 3] = -1
    a = [a0, a1]
    for i, j in itertools.product(range(2), repeat=2):
        equal(a[i]*a[j]+a[j]*a[i], S.zeros(4), 'annihilator CAR')
        equal(a[i]*a[j].T+a[j].T*a[i], I if i == j else S.zeros(4), 'star CAR')
    v = S.Matrix([0, 1, r, 0])/S.sqrt(1+r*r)
    rho = S.simplify(v*v.T)
    def cov(state):
        return S.Matrix(2, 2, lambda i, j:S.simplify(S.trace(state*a[i].T*a[j])))
    P = S.Matrix([[1, r], [r, r*r]])/(1+r*r)
    equal(cov(rho), P, 'actual two-face covariance')
    equal(P*P, P, 'full projector')
    p = r*r/(1+r*r)
    x = r/(1+r*r)
    equal(x*x, p*(1-p), 'maximal coherence')
    n0, n1 = a0.T*a0, a1.T*a1
    deph = S.diag(0, 1-p, p, 0)
    prod = S.diag(p*(1-p), (1-p)**2, p**2, p*(1-p))
    for state in [rho, deph, prod]:
        equal(S.trace(state), 1, 'normalized states')
        equal(S.trace(state*n0), 1-p, 'left occupation')
        equal(S.trace(state*n1), p, 'right occupation')
    equal(cov(deph), cov(prod), 'same mixed full two-point data')
    equal(S.trace(deph*n0*n1), 0, 'dephased four point')
    equal(S.trace(prod*n0*n1), p*(1-p), 'product four point')
    equal(S.trace(rho*n0*n1), P.det(), 'pure four point Wick')
    filled = (a0+r*a1)/S.sqrt(1+r*r)
    empty = (-r*a0+a1)/S.sqrt(1+r*r)
    Q = S.simplify(empty.T*empty+filled*filled.T)
    equal(Q*v, S.zeros(4, 1), 'unique target zero vector')
    require(Q.eigenvals() == {S.Integer(0):1, S.Integer(1):2, S.Integer(2):1}, 'counter spectrum')
    require((Q-(I-rho)).eigenvals() == {S.Integer(0):3, S.Integer(1):1}, 'gap operator inequality')
    for state in [rho, deph, prod]:
        C = cov(state)
        delta = S.trace((S.eye(2)-P)*C)+S.trace(P*(S.eye(2)-C))
        equal(S.trace(state*Q), delta, 'delta from two-point data alone')
    t = S.symbols('t', real=True)
    Ct = S.Matrix([[1-p, t*x], [t*x, p]])
    delta = S.simplify(S.trace((S.eye(2)-P)*Ct)+S.trace(P*(S.eye(2)-Ct)))
    equal(delta, 4*p*(1-p)*(1-t), 'coherence error certificate')
    # Exact noncommutative word control in the projected particle/hole basis.
    # Compare matrix Wick recursion, without assuming it in the source matrix.
    ops = [a0, a1, a0.T, a1.T]
    def expect(mat):
        return S.simplify((v.T*mat*v)[0])
    def wick(word):
        if not word:
            return S.Integer(1)
        if len(word) % 2:
            return S.Integer(0)
        return S.Add(*(S.Integer(-1)**(j-1)*expect(ops[word[0]]*ops[word[j]])*
                       wick(word[1:j]+word[j+1:]) for j in range(1, len(word))))
    checked = 0
    for word in itertools.product(range(4), repeat=4):
        mat = I
        for j in word:
            mat = mat*ops[j]
        equal(expect(mat), wick(word), 'all ordered four-field words')
        checked += 1
    require(cov(deph) != P, 'missing cross coherence detected')
    return dict(claim='UR.NATIVE.ARCCHARGE.PROJECTOR-STATE.01',
        verdict='PASS_EXACT_SCOPED_PROJECTOR_STATE_AND_MARGINAL_CONTROLS',
        actual_annulus_graph_projector=True, ordered_four_field_words_checked=checked,
        arbitrary_positive_state_analytic_theorem='PROJECTOR-STATE-RIGIDITY.txt Section2',
        additional_gaussian_premise_needed_for_projector_uniqueness=False,
        actual_lowest_ns_mode=dict(p=str(S.N(1/(1+S.exp(S.pi)),24)),
                                  x=str(S.N(1/(2*S.cosh(S.pi/2)),24))),
        same_marginal_controls_have_different_fourpoint=True,
        excitation_counter_spectrum={'0':1,'1':2,'2':1},
        coherence_defect=str(S.factor(delta)),
        infinite_source_identity_verified_by_finite_samples=False,
        primitive_p1_joint_covariance_derived=False, promotion=False)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()
    result = json.dumps(run(), indent=2, ensure_ascii=False)+'\n'
    if args.out:
        args.out.write_text(result)
    print(result, end='')
