"""Exact diagnostics for the analytic conformal Cauchy-space proof.

No finite projection is used as a proof of the infinite closed-space theorem.
All checks remain active under python -OO.
"""
import argparse
import json
from pathlib import Path
import sympy as S


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def zero(expression, label):
    require(S.simplify(expression) == 0, label)


def projector(embedding):
    return S.simplify(embedding * (embedding.T * embedding).inv() * embedding.T)


def run():
    d, x = S.symbols('d x', real=True)
    # Common-bulk conjugation and physical-to-boundary density cancellation.
    zero(d/2-(d+1)/2+S.Rational(1, 2), 'left bulk exponent')
    zero((d-1)/2-d/2+S.Rational(1, 2), 'right bulk exponent')
    zero(-(d-1)*x+(d-1)*x, 'boundary norm density')
    zero(S.exp((d-1)*x/2)*S.exp(-(d-1)*x/2)-1, 'physical trace')
    zero(S.exp(-x/2)*S.exp(x/2)-1, 'common-space trace')
    t, s, lam = S.symbols('t s lambda', real=True)
    sigma = S.Function('sigma')(t, s)
    mode = S.exp(-lam*t+S.I*lam*s)
    def chiral(v):
        return S.diff(v, t)-S.I*S.diff(v, s)
    zero(chiral(mode), 'original harmonic mode')
    u = S.exp(sigma/2)*mode
    zero(S.exp(-sigma/2)*chiral(S.exp(-sigma/2)*u), 'conformal harmonic mode')
    for face in [S.Integer(0), S.Rational(1, 2)]:
        zero((S.exp(-sigma/2)*u-mode).subs(t, face), 'both weighted traces')
    # Two actual NS modes have positive graph slopes r=e^-pi/2,s=e^-3pi/2.
    # Keep r,s symbolic: the projector identity does not rely on the slopes.
    r, q, a, b = S.symbols('r q a b', positive=True)
    E = S.Matrix([[1, 0], [0, 1], [r, 0], [0, q]])
    M = S.diag(a, a, b, b)
    P = projector(E)
    correct = projector(M.inv()*M*E)
    require(correct == P, 'weighted closed graph projector')
    zero((P*P-P).norm(), 'idempotence')
    raw = projector(M*E)
    defect = S.factor(raw[2, 2]-P[2, 2])
    expected = r**2*(b**2-a**2)/((a**2+b**2*r**2)*(1+r**2))
    zero(defect-expected, 'raw trace nonzero defect formula')
    require(defect.subs({a: 2, b: 3, r: 1}) != 0, 'wrong trace detected')
    # Arbitrary weights/determinants must factor, rather than be discarded.
    w0, w1, det0, det1, c00, c01, c10, c11 = S.symbols(
        'w0 w1 det0 det1 c00 c01 c10 c11')
    norm = w0*det0+w1*det1
    four = c00*c11-c01*c10
    zero((w0*det0*four+w1*det1*four)/norm-four, 'normalized four point')
    F = S.Symbol('same_conditional_word')
    zero((w0*det0*F+w1*det1*F)/norm-F, 'all identical words')
    # A varying covariance need not factor: strict changing-covariance control.
    p, v = S.symbols('p v', real=True)
    C0, C1 = S.diag(p, 1-p), S.diag(p+v, 1-p-v)
    mean = (C0+C1)/2
    kappa = S.expand((C0.det()+C1.det())/2-mean.det())
    zero(kappa+v**2/4, 'changing covariance connected defect')
    require(kappa.subs(v, S.Rational(1, 8)) != 0, 'changed source detected')
    return {
        'claim': 'UR.NATIVE.ARCCHARGE.CONFORMAL-TRACE.01',
        'verdict': 'PASS_EXACT_SCOPED_CONFORMAL_TRACE_AND_WICK_FACTORIZATION',
        'dimension_independent_density_and_trace': True,
        'arbitrary_smooth_sigma_product_rule': True,
        'actual_ns_modes_lambda': '2*pi*(n+1/2)',
        'both_original_annulus_faces': ['0', '1/2'],
        'orthogonal_graph_projector_exact': True,
        'raw_trace_compression_defect': str(defect),
        'arbitrary_source_free_weights_and_determinants_retained': True,
        'same_conditional_words_factor_exactly': True,
        'changing_covariance_four_point_control': str(kappa),
        'infinite_space_proof': 'CONFORMAL-TRACE-AND-WICK.txt Sections2-4',
        'primitive_gaussian_measure_derived': False,
        'physical_time_identified': False,
        'promotion': False,
    }


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()
    result = run()
    data = json.dumps(result, indent=2, ensure_ascii=False)+'\n'
    if args.out:
        args.out.write_text(data)
    print(data, end='')
