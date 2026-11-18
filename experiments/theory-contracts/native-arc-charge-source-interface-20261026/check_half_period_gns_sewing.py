"""Exact graded HS/CAR connection and existing interaction countercontrol."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as S


def require(ok, name):
    if not ok:
        raise RuntimeError(name)


def equal(x, y, name):
    if isinstance(x, S.MatrixBase):
        require(all(S.simplify(z) == 0 for z in x-y), name)
    else:
        require(S.simplify(x-y) == 0, name)


def annihilators(n):
    mats = []
    for j in range(n):
        a = S.zeros(2**n)
        for k in range(2**n):
            if k & (1 << j):
                a[k ^ (1 << j), k] = (-1)**((k & ((1 << j)-1)).bit_count())
        mats.append(a)
    return mats


def double_data(weights):
    n = len(weights).bit_length()-1
    a = annihilators(n)
    rho = S.diag(*weights)/sum(weights)
    root = S.diag(*(S.sqrt(x/sum(weights)) for x in weights))
    parity = S.diag(*((-1)**k.bit_count() for k in range(2**n)))
    def action(op, X):
        side, j, star = op
        A = a[j].T if star else a[j]
        if side == 'R':
            return A*X
        return (1 if star else -1)*parity*X*A*parity
    def word(ops):
        X = root
        for op in reversed(ops):
            X = action(op, X)
        return S.simplify(S.trace(root*X))
    ops = [('L',j,False) for j in range(n)]+[('R',j,False) for j in range(n)]
    def adj(op):
        return op[0], op[1], not op[2]
    covariance = S.Matrix(2*n,2*n,lambda i,j:word([adj(ops[i]),ops[j]]))
    return a, rho, root, parity, action, word, ops, adj, covariance


def run():
    z = S.symbols('z', positive=True)
    equal(z/(1+z*z), 1/(z+1/z), 'half period identity')
    # One physical temporal evaluation cannot be an independent CAR copy.
    a = annihilators(1)[0]
    equal((z*a)*a.T+a.T*(z*a), z*S.eye(2), 'raw temporal anticommutator')
    require(z*S.eye(2) != S.zeros(2), 'wrong temporal identification rejected')
    free = double_data([S.Rational(8,15),S.Rational(2,15),S.Rational(4,15),S.Rational(1,15)])
    a,rho,root,parity,act,word,ops,adj,P = free
    units = []
    for i,j in itertools.product(range(4),repeat=2):
        E = S.zeros(4); E[i,j] = 1; units.append(E)
    tested = 0
    for i,j in itertools.product(range(4),repeat=2):
        for E in units:
            equal(act(ops[i],act(ops[j],E))+act(ops[j],act(ops[i],E)), S.zeros(4), 'double annihilator CAR')
            equal(act(ops[i],act(adj(ops[j]),E))+act(adj(ops[j]),act(ops[i],E)), E if i==j else S.zeros(4), 'double star CAR')
            tested += 2
    C = S.diag(S.Rational(1,5),S.Rational(1,3))
    X = S.diag(S.Rational(2,5),S.sqrt(2)/3)
    equal(P, (S.eye(2)-C).row_join(X).col_join(X.row_join(C)), 'geometric free graph covariance')
    equal(P*P, P, 'complete projector')
    # All grouped four-field matrix elements, independently from the 2pt determinant.
    grouped = 0
    for (i,j),(k,l) in itertools.product(itertools.combinations(range(4),2),repeat=2):
        actual = word([adj(ops[i]),adj(ops[j]),ops[l],ops[k]])
        equal(actual, P[i,k]*P[j,l]-P[i,l]*P[j,k], 'joint grouped four point')
        grouped += 1
    wrong = []
    for j in range(2):
        wrong.append(S.simplify(S.trace(root*a[j].T*parity*root*a[j])))
    equal(wrong[0], S.Rational(2,15), 'spectator parity wrong control0')
    equal(wrong[1], S.sqrt(2)/5, 'spectator parity wrong control1')
    require(wrong != [X[0,0],X[1,1]], 'wrong Klein factor detected')
    # Reuse precisely the old explicit quarter-clock interaction control.
    inter = double_data([S.Integer(1),S.Rational(1,4),S.Rational(1,64),S.Rational(1,65536)])
    a,rho,root,parity,act,word,ops,adj,D = inter
    c = [S.trace(rho*A.T*A) for A in a]
    m = [S.simplify(S.trace(root*A.T*root*A)) for A in a]
    equal(S.Matrix(c),S.Matrix([S.Rational(16385,82945),S.Rational(1025,82945)]),'old interaction occupations')
    equal(S.Matrix(m),S.Matrix([S.Rational(32800,82945),S.Rational(8320,82945)]),'old interaction cross')
    deficits = [S.factor(cj*(1-cj)-mj*mj) for cj,mj in zip(c,m)]
    for gap in deficits:
        equal(gap,S.Rational(14745600,82945**2),'strict saturation defect')
        require(gap>0,'nonfree source detected')
    require(D*D != D,'pure HS vector does not force projection covariance')
    # The two spectator sectors carry different modular half factors.
    factors = []
    for A in a:
        row=[]
        for i,j in itertools.product(range(4),repeat=2):
            if A[i,j]:
                row.append(S.simplify(root[i,i]/root[j,j]))
        factors.append(row)
    require(factors==[[S.Integer(2),S.Integer(32)],[S.Integer(8),S.Integer(128)]],'operator-valued modular factors')
    n0,n1 = (A.T*A for A in a)
    kappa = S.factor(S.trace(rho*n0*n1)-c[0]*c[1])
    equal(kappa,-S.Rational(16711680,82945**2),'old interaction fourpoint')
    return dict(claim='UR.NATIVE.ARCCHARGE.HALF-PERIOD-SEWING.01',
        verdict='PASS_EXACT_HALF_PERIOD_GRADED_GNS_AND_SATURATION_CONTROLS',
        half_period_is_existing_geometric_cross=True,
        raw_temporal_copy_rejected=True,
        independent_double_car_matrix_unit_checks=tested,
        grouped_joint_four_field_checks=grouped,
        free_cross=[str(X[i,i]) for i in range(2)],
        wrong_left_only_parity_cross=[str(v) for v in wrong],
        interaction_covariance=[str(v) for v in c],
        interaction_cross=[str(v) for v in m],
        interaction_saturation_deficits=[str(v) for v in deficits],
        interaction_half_modular_factors=[[str(v) for v in row] for row in factors],
        interaction_fourpoint_defect=str(kappa),
        finite_saturation_iff_proof='HALF-PERIOD-GNS-SEWING.txt Section5',
        infinite_iff_inferred_from_finite_proof=False,
        actual_primitive_p1_sewing_identified=False,promotion=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path);args=p.parse_args()
    result=json.dumps(run(),indent=2,ensure_ascii=False)+'\n'
    if args.out:args.out.write_text(result)
    print(result,end='')
