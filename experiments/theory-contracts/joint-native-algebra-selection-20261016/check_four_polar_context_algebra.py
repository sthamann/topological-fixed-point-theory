"""Exact joint algebra of the actual four polar marks; no physical instrument fit."""
from pathlib import Path
import hashlib
import itertools
import json
import sympy as s

BASE = Path(__file__).resolve().parent
cl = lambda A: A.applyfunc(s.simplify)


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


def main():
    # Literal chart in check_four_starting_points.py, whose intrinsic source
    # reconstruction and central OPE metric are separately established.
    raw = s.Matrix([[2, 2, -1, -1, -1, -1], [0, 0, -1, -1, 1, 1],
                    [0, 0, -1, 1, -1, 1], [0, 0, 1, -1, -1, 1],
                    [1, -1, 0, 0, 0, 0]])
    N = cl(s.diag(2, s.sqrt(12), s.sqrt(12), s.sqrt(12), s.sqrt(24))*raw/s.sqrt(48))
    need(N*N.T == s.eye(5), "Actual simplex normalization")
    need(cl(N.T*N) == s.eye(6)-s.ones(6)/6, "Actual six-ray Gram")
    a = (s.sqrt(6)-1)/5
    e = cl((N[:, 0]-N[:, 1])/s.sqrt(2))
    Z = e*e.T
    V = cl(s.Matrix.hstack(*[(N[:, 0]+N[:, 1]-2*a*N[:, q])/s.sqrt(2)
                            for q in range(2, 6)]))
    t = (11+4*s.sqrt(6))/25
    need(cl(V.T*V) == (1-t)*s.eye(4)+t*s.ones(4), "Complete exact overlap Gram")
    need(cl(e.T*V) == s.zeros(1, 4) and V.rank() == 4, "Actual orthogonal 1+4 vector decomposition")
    Ps = []
    for j, q in enumerate(range(2, 6)):
        A = s.Matrix.hstack(*[N[:, k]-a*N[:, q] for k in (0, 1)])
        P = cl(A*A.T)
        need(cl(P*P-P) == s.zeros(5) and s.simplify(s.trace(P)) == 2, "Actual rank-two polar mark")
        need(cl(P-Z-V[:, j]*V[:, j].T) == s.zeros(5), "Common line and four residual rays")
        Ps.append(P)
    S = sum(Ps, s.zeros(5))
    lp, lm = 1+3*t, 1-t
    spectral_Z = cl((S-lp*s.eye(5))*(S-lm*s.eye(5))/((4-lp)*(4-lm)))
    need(spectral_Z == Z, "Common line lies in generated algebra, not added by hand")
    Qs = [P-Z for P in Ps]
    units = []
    for j, k in itertools.product(range(4), repeat=2):
        A = cl(Qs[j]*Qs[k]/(1 if j == k else t))
        need(cl(A-V[:, j]*V[:, k].T) == s.zeros(5), "All residual outer products generated")
        units.append(A.reshape(25, 1))
    need(s.Matrix.hstack(*units).rank() == 16, "Full M4 and no omitted residual direction")
    commutator_norms = []
    for j, k in itertools.combinations(range(4), 2):
        C = cl(Ps[j]*Ps[k]-Ps[k]*Ps[j])
        norm = s.simplify(s.trace(C.H*C))
        need(s.simplify(norm-2*t*t*(1-t*t)) == 0 and norm != 0,
             "Every pair fails sharp joint measurability")
        commutator_norms.append(str(norm))
    rays = cl(s.sqrt(s.Rational(6, 5))*N[:, 2:6])
    ray_gram = cl(rays.T*rays)
    need(ray_gram == s.Rational(6, 5)*s.eye(4)-s.ones(4)/5,
         "Four distinct rays are nonorthogonal, not perfectly discriminable")
    # Negative control: orthogonal one-dimensional settings do commute.
    P0, P1 = s.diag(1, 0, 0, 0, 0), s.diag(0, 1, 0, 0, 0)
    need(P0*P1 == P1*P0 == s.zeros(5), "Orthogonal control is compatible")
    tensor_path = BASE/'fold_response_tensors_results.json'
    original = json.loads(tensor_path.read_text())
    expected = {'1/4': s.Rational(7882461, 4398046511104),
                '3/4': s.Rational(3201175, 13194139533312)}
    mixed = {}
    for alpha, data in original['cases'].items():
        D1 = s.Matrix(data['wick_tensors'][1])
        need(cl(D1-D1.H) == s.zeros(5), "Original full D1 is Hermitian")
        w = cl((s.eye(5)-Z)*D1*e)
        n2 = s.simplify((w.H*w)[0])
        need(n2 == expected[alpha] and n2 > 0, "Actual source correlation breaks residual central split")
        # Together with M4 this cross-block connects e to every residual
        # direction. The commutant shrinks from C Z+C(I-Z) to C I; hence M5.
        mixed[alpha] = {'cross_block_norm_squared': str(n2),
                        'response_algebra': 'M5(C)', 'commutant_dimension': 1}
    out = {'status': 'PASS_EXACT_FOUR_POLAR_CONTEXT_ALGEBRA',
           'overlap_t': str(t), 'four_ray_gram': [[str(x) for x in row] for row in ray_gram.tolist()],
           'polar_generated_algebra': 'C direct_sum M4(C)',
           'polar_commutant_dimension': 2, 'polar_algebra_dimension': 17,
           'pair_commutator_norm_squared': commutator_norms,
           'original_D1_cases': mixed,
           'source_tensor_sha256': hashlib.sha256(tensor_path.read_bytes()).hexdigest(),
           'scope': 'Actual finite quartic observation space, fixed global S={0,1}, original polar rule and NS current metric. Algebraic classification and conditional sharp-instrument constraint only. No assertion of physical implementation, full-source superselection, unique q, preferred measurement, or complete TFPT solution.'}
    (BASE/'four_polar_context_algebra_results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
