# Native source structural probe; unconditional checks also run under -OO.
from pathlib import Path
import ast, itertools, json, hashlib
import numpy as np
import sympy as s
import mpmath as mp
ROOT = Path(__file__).resolve().parents[3]
src = ROOT / 'verification/v1_e8_glue.py'
node = next((n for n in ast.parse(src.read_text()).body if isinstance(n, ast.FunctionDef) and n.name == 'e8_roots'))
env = {'itertools': itertools, 'np': np}
exec(compile(ast.Module(body=[node], type_ignores=[]), str(src), 'exec'), env)
roots = [s.Matrix([s.Rational(float(x)) for x in r]) for r in env['e8_roots']()]
half = [r for r in roots if all((abs(x) == s.Rational(1, 2) for x in r))]
families = sorted(set((tuple(r[5:]) for r in half if sum((bool(x < 0) for x in r[:5])) % 2 == 1)))
T = s.Matrix.hstack(*[s.Matrix(f) for f in families])
I = s.eye(4)
J = s.ones(4)
P = I - J / 4
r = s.zeros(4)
for j in range(4):
    r[(j + 1) % 4, j] = 1
R = s.diag(0, 1, 1, 1)
u = s.ones(4, 1) / 2
H = s.diag(1, -1, 0, 0)
msrc = ROOT / 'verification/v484_seam_contact_unit.py'
run = next((n for n in ast.parse(msrc.read_text()).body if isinstance(n, ast.FunctionDef) and n.name == 'run'))
cm = next((n.value for n in run.body if isinstance(n, ast.Assign) and any((isinstance(t, ast.Name) and t.id == 'Cm' for t in n.targets))))
C = eval(compile(ast.Expression(cm), str(msrc), 'eval'), {'sp': s})
E = sum((r ** j * C * r.T ** j for j in range(4)), s.zeros(4)) / 4
x = s.symbols('x')
det = s.factor((I - x * C).det())
series = s.series(-s.log(det), x, 0, 5).removeO().expand()
residue_series = s.series(-s.log((I - x * P * C * P).det()), x, 0, 5).removeO().expand()
if not (len(roots) == 240 and len(half) == 128 and (len(families) == 4)):
    raise AssertionError('Decisive structural check failed: len(roots) == 240 and len(half) == 128 and (len(families) == 4)')
if not (T * T.T == s.eye(3) and T.T * T == P):
    raise AssertionError('Decisive structural check failed: T * T.T == s.eye(3) and T.T * T == P')
if not all((s.Matrix(f).dot(s.Matrix(f)) == s.Rational(3, 4) for f in families)):
    raise AssertionError('Decisive structural check failed: all((s.Matrix(f).dot(s.Matrix(f)) == s.Rational(3, 4) for f in families))')
if not all((sum((v * v for v in a[:5])) == s.Rational(5, 4) for a in half)):
    raise AssertionError('Decisive structural check failed: all((sum((v * v for v in a[:5])) == s.Rational(5, 4) for a in half))')
if not (P * P == P and s.trace(P) == 3 and (P * r == r * P)):
    raise AssertionError('Decisive structural check failed: P * P == P and s.trace(P) == 3 and (P * r == r * P)')
if not (R * r != r * R and P * H != H * P):
    raise AssertionError('Decisive structural check failed: R * r != r * R and P * H != H * P')
if not (E == C and E != C / 4):
    raise AssertionError('Decisive structural check failed: E == C and E != C / 4')
if not s.expand(det - (1 - 4 * x) * (1 + 2 * x) ** 2) == 0:
    raise AssertionError('Decisive structural check failed: s.expand(det - (1 - 4 * x) * (1 + 2 * x) ** 2) == 0')
if not (series.coeff(x, 4) == 72 and residue_series.coeff(x, 4) == 8):
    raise AssertionError('Decisive structural check failed: series.coeff(x, 4) == 72 and residue_series.coeff(x, 4) == 8')
deck = r * r
even = (I + deck) / 2
odd = (I - deck) / 2
if not (s.trace(P * even) == 1 and s.trace(P * odd) == 2):
    raise AssertionError('Decisive structural check failed: s.trace(P * even) == 1 and s.trace(P * odd) == 2')
carrier_monodromy = s.simplify(s.exp(-2 * s.pi * s.I * s.Rational(5, 4)))
family_monodromy = s.simplify(s.exp(-2 * s.pi * s.I * s.Rational(3, 4)))
if not (carrier_monodromy == -s.I and family_monodromy == s.I):
    raise AssertionError('Decisive structural check failed: carrier_monodromy == -s.I and family_monodromy == s.I')
if not carrier_monodromy * family_monodromy == 1:
    raise AssertionError('Decisive structural check failed: carrier_monodromy * family_monodromy == 1')
if not carrier_monodromy ** 2 == family_monodromy ** 2 == -1:
    raise AssertionError('Decisive structural check failed: carrier_monodromy ** 2 == family_monodromy ** 2 == -1')
v985 = ROOT / 'verification/v985_quillen_channel_swap.py'
top = ast.parse(v985.read_text()).body
pattern = next((n.value for n in top if isinstance(n, ast.Assign) and any((isinstance(t, ast.Name) and t.id == 'JUMP_PATTERN' for t in n.targets))))
bp = eval(compile(ast.Expression(pattern), str(v985), 'eval'), {'sp': s})
B = s.Matrix(4, 4, lambda i, j: bp[(j - i) % 4])
D = s.diag(1, -1, 1, -1)
A = D * B * D
if not (A * A == C * C and all((e >= 0 for e in A.eigenvals()))):
    raise AssertionError('Decisive structural check failed: A * A == C * C and all((e >= 0 for e in A.eigenvals()))')
if not (D * r * D == -r and D * deck == deck * D and (D != deck)):
    raise AssertionError('Decisive structural check failed: D * r * D == -r and D * deck == deck * D and (D != deck)')
if not (deck * B * deck.T == B and D * B * D != B):
    raise AssertionError('Decisive structural check failed: deck * B * deck.T == B and D * B * D != B')
Psign = D * P * D
if not ((P * C * P).rank() == 2 and (Psign * C * Psign).rank() == 3):
    raise AssertionError('Decisive structural check failed: (P * C * P).rank() == 2 and (Psign * C * Psign).rank() == 3')
mp.mp.dps = 60
c = 1 / (8 * mp.pi)
base = 1 / (6 * mp.pi)
delta = 48 * c ** 4

def solve(p, dt=delta):

    def F(a):
        q = dt * mp.exp(-2 * a)
        phi = base + q / (1 - q) ** p
        return a ** 3 - 2 * c ** 3 * a ** 2 - mp.mpf(164) / 5 * c ** 6 * mp.log(1 / phi)
    a = mp.findroot(F, mp.mpf('.0073'))
    return mp.nstr(1 / a, 40)
out = {'scope': 'Native E8 lattice and finite seam contact matrix; conformal reading conditional on the canonical lattice VOA. No physical alpha derivation.', 'native_roots': len(roots), 'half_roots': len(half), 'family_weights': [list(map(str, f)) for f in families], 'weight_gram': str(T.T * T), 'tight_frame': str(T * T.T), 'family_norm_squared': '3/4', 'carrier_norm_squared': '5/4', 'conformal_weights': {'carrier': '5/8', 'family': '3/8', 'glue_current': '1'}, 'projection': {'zero_sum_rank': 3, 'commutes_mark_rotation': True, 'coordinate_three_commutes_rotation': False, 'zero_sum_commutes_su4_generator': False, 'rank_ratio': '4/3', 'conditional_deck_r_squared_even_rank': 1, 'conditional_deck_r_squared_odd_rank': 2}, 'conditional_chiral_two_point_monodromy': {'carrier': str(carrier_monodromy), 'family': str(family_monodromy), 'full_glue': str(carrier_monodromy * family_monodromy), 'two_turns_each': '-1', 'four_turns_each': '1', 'scope': 'Canonical lattice-VOA correlator analytic continuation; not an identification with a geometric collar deck.'}, 'native_channel_swap': {'jump_matrix': str(B), 'character_twist': str(D), 'absolute_mark_matrix': str(A), 'twist_conjugates_rotation_to': '-r', 'geometric_mark_half_turn': 'r^2', 'twist_equals_permutation_half_turn': False, 'zero_sum_contact_rank': 2, 'sign_complement_contact_rank': 3, 'scope': 'The finite v985 channel match is confirmed. Interpretation as a geometric lifted deck requires an explicit fibre action.'}, 'contact': {'matrix': str(C), 'group_average_equals_input': True, 'group_average_is_one_quarter_input': False, 'determinant': str(det), 'minus_log_series': str(series), 'fourth_coefficient': 72, 'zero_sum_fourth_coefficient': 8, 'target_flat_count': 48}, 'alpha_inverse_same_assembled_closure': {'carrier_exponent_5_4': solve(mp.mpf(5) / 4), 'full_glue_exponent_2': solve(mp.mpf(2)), 'no_defect': solve(mp.mpf(5) / 4, mp.mpf(0))}, 'sources': {str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest() for f in [src, msrc, v985]}, 'verdict': 'PARTIAL', 'excluded': ['Bare unweighted four-mark logdet supplies exactly 48 c3^4 without a response/matching rule.', 'Group conjugation average of the invariant Green matrix itself divides that matrix by four.', 'Zero-sum mark projection is automatically an SU(4)-equivariant family selection.'], 'open': ['Physical inverse-norm normalization of phi_base.', 'Original puncture measure and its finite-to-continuum response dictionary.', 'Carrier-only readout, q(alpha) mapping and full U(1) determinant stationarity.', 'Raw seam to marked affine E8 current algebra.']}
print(json.dumps(out, indent=2))
