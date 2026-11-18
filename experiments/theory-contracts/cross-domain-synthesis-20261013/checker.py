#!/usr/bin/env python3
import hashlib
import importlib.util
import json
from pathlib import Path

import sympy as sp

ROOT = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
CP_CHECKER = ROOT / 'experiments/theory-contracts/marked-seam-cp-lift-20260921/checker.py'
EXPECTED_CP_SHA256 = '842c8c715d740ff0927c18c01686014d66fcbec9ea782a68c0024e929eecae98'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


actual_sha = hashlib.sha256(CP_CHECKER.read_bytes()).hexdigest()
if actual_sha != EXPECTED_CP_SHA256:
    raise RuntimeError({'expected': EXPECTED_CP_SHA256, 'actual': actual_sha})

cp = load('marked_cp_source', CP_CHECKER)
bridge = load(
    'compiler_clifford_bridge',
    ROOT / 'experiments/theory-contracts/compiler-clifford-bridge/checker.py',
)
src, _ = bridge.inherited(ROOT)
finite = bridge.finite_data(src)
q = finite['q']
zero = (0, 0, 0, 0)
words = {v: bridge.monomial(v, bridge.generators()) for v in src.W16}
o5 = [v for v in src.W16 if v != zero and q[v] == 0]
gammas = [words[v] for v in o5]

I4 = sp.eye(4)
Z4 = sp.zeros(4)
hermitian = all(g.H == g for g in gammas)
clifford = all(
    sp.simplify(gammas[i] * gammas[j] + gammas[j] * gammas[i])
    == (2 * I4 if i == j else Z4)
    for i in range(5)
    for j in range(5)
)
determinants = [sp.simplify(g.det()) for g in gammas]
det_one = determinants == [sp.Integer(1)] * 5

# For every pure q=|psi><psi|, sum_j <psi|g_j|psi>^2 = 1 follows if
# sum_j g_j tensor g_j is the identity on Sym^2(C^4).
swap = sp.zeros(16)
for i in range(4):
    for j in range(4):
        swap[4 * j + i, 4 * i + j] = 1
sym = (sp.eye(16) + swap) / 2
K = sum((sp.kronecker_product(g, g) for g in gammas), sp.zeros(16))
spin5_fierz = sp.simplify(sym * (K - sp.eye(16)) * sym) == sp.zeros(16)

# Phi(q)=7/12 q + 1/12 sum_j g_j q g_j, hence pure-state fidelity is 2/3.
pure_fidelity = sp.Rational(7, 12) + sp.Rational(1, 12)
rank3_leakage = 1 - pure_fidelity
normalized_rank3_leakage = rank3_leakage / 3
n = sp.symbols('n', integer=True, nonnegative=True)
normalized_n_step_leakage = (
    sp.Rational(1, 4)
    - sp.Rational(1, 12) * sp.Rational(1, 3) ** n
    - sp.Rational(1, 6) * sp.Rational(2, 3) ** n
)
n1 = sp.simplify(normalized_n_step_leakage.subs(n, 1))
n6 = sp.simplify(normalized_n_step_leakage.subs(n, 6))

result = {
    'source': str(CP_CHECKER.relative_to(ROOT)),
    'source_sha256': actual_sha,
    'o5_labels': [list(v) for v in o5],
    'actual_five_are_hermitian_clifford': hermitian and clifford,
    'actual_five_determinants': [str(x) for x in determinants],
    'actual_five_are_in_G_tilde_given_v783_Pauli_membership': det_one,
    'symmetric_tensor_identity': spin5_fierz,
    'identity': 'Pi_sym (sum_j g_j tensor g_j - I16) Pi_sym = 0',
    'pure_state_sum_expectation_squares': '1',
    'pure_state_phi_fidelity': str(pure_fidelity),
    'rank3_projector_leakage': str(rank3_leakage),
    'normalized_P_over_3_leakage': str(normalized_rank3_leakage),
    'normalized_n_step_leakage': '1/4 - (1/3)^n/12 - (2/3)^n/6',
    'normalized_leakage_n1': str(n1),
    'normalized_leakage_n6': str(n6),
    'correlated_lift_identity': 'For U_j=S_j tensor g_j and rho=(I16/16) tensor (P/3), every branch leaves I16/16 fixed, so the family leakage equals the bare M4 leakage at every n.',
    'conditioned_branch_scope': 'Each stored random-unitary branch transports P to P_j=g_j P g_j^dagger (and I16 tensor P to I16 tensor P_j); fixed-P leakage is caused by forgetting the branch record.',
    'scope': 'Every rank-3 projector P in the actual finite M4 average channel, and its conditional correlated G_tilde lift on maximally mixed carrier; no physical instrument selection, no no-go for recorded/conditioned branches, dilations, other states, or all TFPT.',
}
# Direct replay of the actual channel, independently of the closed n-formula.
P = sp.diag(1, 1, 1, 0)
rho = P / 3
current = rho
direct = {}
branch_covariance = True
for g in gammas:
    moved = g * P * g.H
    transformed = g * rho * g.H
    branch_covariance = branch_covariance and moved * transformed * moved == transformed
for step in range(1, 7):
    current = sp.simplify(sp.Rational(7,12)*current + sum((g*current*g.H for g in gammas), sp.zeros(4))/12)
    direct[str(step)] = str(sp.simplify(sp.trace((I4-P)*current)))
checks = {
    'native_five_Clifford': hermitian and clifford,
    'native_five_determinant_one': det_one,
    'universal_symmetric_tensor_identity': spin5_fierz,
    'actual_channel_step1': direct['1'] == str(n1),
    'actual_channel_step6': direct['6'] == str(n6),
    'transported_projector_branch_covariance': branch_covariance,
}
if not all(checks.values()):
    raise RuntimeError(checks)
result['checks'] = checks
result['direct_channel_leakage'] = direct
result['passed'] = sum(checks.values())
result['total'] = len(checks)
result['physical_gates_closed'] = []
result['complete_TFPT_solution'] = False

out = Path(__file__).with_name('results.json')
out.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
