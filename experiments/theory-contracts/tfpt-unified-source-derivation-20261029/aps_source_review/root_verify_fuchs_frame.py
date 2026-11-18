"""Independent substitution of the explicit v117 frame, over Q(omega)(w).

Checks all three ORIGINAL equations and the determinant.  Does not fit a
target matrix, derive row2/row1 recursively, or use the candidate solver.
"""
import json
from pathlib import Path
import sympy as S

w, om = S.symbols('w om')
K = S.QQ.frac_field(w)
mod = S.Poly(om**2 + om + 1, om, domain=K)

def r(x):
    n, d = S.fraction(S.cancel(x))
    n = S.Poly(n, om, domain=K).rem(mod)
    d = S.Poly(d, om, domain=K).rem(mod)
    return S.cancel((n*S.invert(d, mod)).rem(mod).as_expr())

s = (1+w**3)/(1-w**3)  # z^2
ell = 3*w**2/(1-w**6)  # z'/z in w coordinate
D = r(s*s-1)
H = (w+1)/(w+om)
k = r(S.diff(H, w)/(2*H))
p1 = 4*(-(1+om)*w**3+(1-om)*w**2+(om+2)*w+om)/(9*(w**3-1))
p2 = ((1+om)*w+om)/(w+1)
p3 = ((1+om)*w-om)/(w+1)
residual = [
    r(S.diff(p1,w)+k*p1-ell*s*s/D*(2*p1+S.Rational(8,9)*p2)),
    r(S.diff(p2,w)+(k+ell)*p2-ell/D*(p1+s*s*p2+S.Rational(5,9)*s*p3)),
    r(S.diff(p3,w)+(k+ell)*p3-ell/D*(s*p2+s*s*p3)),
]
P = S.Matrix([[r(p.subs(w,om**j*w)) for j in range(3)]
              for p in (p1,p2,p3)])
det_scaled = r(s*P.det()/D)
answer = {
    'original_three_equation_residuals': [str(x) for x in residual],
    's_detP_over_z4_minus_1': str(det_scaled),
    'determinant_expected_residual': str(r(det_scaled+2*(2*om+1))),
    'detG': '-27*sqrt(2)*(2*omega+1)*(z^4-1)/10',
    'det_normalized_Psi': '1-z^4',
    'method': 'exact direct substitution over Q(omega)(w), no fitted coefficients',
}
if any(x != 0 for x in residual) or r(det_scaled+2*(2*om+1)) != 0:
    raise ArithmeticError(answer)
Path(__file__).with_name('ROOT-FUCHS-FRAME-CHECK.json').write_text(
    json.dumps(answer, indent=2)+'\n')
print(json.dumps(answer, indent=2))
