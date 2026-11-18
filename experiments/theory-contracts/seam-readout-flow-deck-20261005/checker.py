"""Small exact tests on the existing TFPT lattice, transfer and NS lift.

The Ward calculation is conditional on the canonical free lattice current
algebra. It does not identify the physical seam measure or a 4D photon.
All checks remain active under python -OO.
"""
from pathlib import Path
from fractions import Fraction as Fr
import ast
import hashlib
import itertools
import json
import sys
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[3]
used = {}

def source(rel):
    p = ROOT / rel
    data = p.read_bytes()
    used[rel] = hashlib.sha256(data).hexdigest()
    return p, ast.parse(data.decode())

def functions(rel, names, env):
    p, tree = source(rel)
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    if {n.name for n in nodes} != set(names):
        raise RuntimeError('Original function missing: ' + rel)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(p), 'exec'), env)

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

env = {'itertools': itertools, 'np': np, 'sp': sp}
functions('verification/v1_e8_glue.py', ['e8_roots'], env)
functions('tfpt_explorer/source_charge_response.py', ['hypercharge_vector'], env)
roots = [sp.Matrix([sp.Rational(float(x)) for x in r]) for r in env['e8_roots']()]
half = [r for r in roots if all(abs(x) == sp.Rational(1, 2) for x in r)]
family = sorted(set(tuple(r[5:]) for r in half if sum(bool(x < 0) for x in r[:5]) % 2 == 1))[0]
carrier = [r for r in half if tuple(r[5:]) == family]
Y = env['hypercharge_vector']()
ynorm = Y.dot(Y)
charges = [Y.dot(p) for p in carrier]
require(len(carrier) == 16 and ynorm == sp.Rational(5, 6), 'Native marked charge carrier changed')
require(sum(q*q for q in charges) == sp.Rational(10, 3), 'Native spinor charge index changed')
require(all(sum(x*x for x in p[:5]) == sp.Rational(5, 4) for p in carrier), 'Carrier norm changed')

# Conditional lattice-VOA Ward identity; not a physical seam identification.
z, w, a, b, d = sp.symbols('z w a b d', nonzero=True)
Q, k, nu = sp.symbols('Q k nu', real=True)
Az = 1/(z-a) - 1/(z-b)
Aw = 1/(w-a) - 1/(w-b)
Z = d**(-nu)
one_z = Z * Q * Az
one_w = Z * Q * Aw
two = Z * (k/(z-w)**2 + Q**2*Az*Aw)
normalized = sp.simplify(two/Z)
connected = sp.simplify(normalized - (one_z/Z)*(one_w/Z))
require(sp.simplify(connected - k/(z-w)**2) == 0, 'Connected Ward identity failed')
require(nu not in normalized.free_symbols and nu not in connected.free_symbols, 'Common vertex amplitude did not cancel')

# Use the actual current/readout functions, including their stated limitations.
sys.path.insert(0, str(ROOT))
from tfpt_explorer.source_charge_response import exact_charge_response
from tfpt_explorer.source_flavor_correlator import exact_original_mark_port
charge_response = exact_charge_response()
mark_port = exact_original_mark_port()
for rel in ['tfpt_explorer/source_charge_response.py', 'tfpt_explorer/source_flavor_correlator.py',
            'tfpt_explorer/neutral_source_response.py']:
    source(rel)
gram = charge_response['gram']
require(gram.det() == 460 and charge_response['R_JY_JY'] == 0
        and charge_response['RY_JY_JY'] == sp.Rational(115, 36),
        'Native family/topological and charge-weighted responses changed')
Cmark, P0 = mark_port['C'], mark_port['P0']
Pfamily = sp.eye(4)-P0
require(mark_port['compressed'] == 4*P0 and Pfamily*P0 == sp.zeros(4),
        'Original scalar port/family channel split changed')
require(Pfamily*Cmark*Pfamily == Cmark-4*P0, 'Complement contact relation failed')

# Read the native v486 rule and the actual selection of its rates.
p, tree = source('verification/v486_transfer_full_rule.py')
run = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'run')
native_names = {'s', 'h', 'sol', 's0', 'h0', 'leak', 'B', 'M'}
nodes = [n for n in run.body if isinstance(n, ast.Assign)
         and any(isinstance(t, ast.Name) and t.id in native_names
                 or isinstance(t, ast.Tuple) and any(isinstance(v, ast.Name) and v.id in native_names for v in t.elts)
                 for t in n.targets)]
me = {'sp': sp}
exec(compile(ast.Module(body=nodes, type_ignores=[]), str(p), 'exec'), me)
B, M = me['B'], me['M']
swap = sp.Matrix([[0, 1], [1, 0]])
fixed = (M-sp.eye(3)).nullspace()
require(len(fixed) == 1 and set(M.eigenvals()) == {sp.Integer(1), sp.Rational(2, 3), sp.Rational(1, 3)}, 'Native fixed space changed')
require(B*swap == swap*B and B[0, 1] == sp.Rational(1, 6), 'Pair parity/hopping relation changed')
logB = (sp.log(2)/2-sp.log(3))*sp.eye(2) + sp.log(2)/2*swap
require(sp.simplify(logB.exp()-B) == sp.zeros(2), 'Pair logarithm failed')
require(sp.simplify(sum(logB.row(0))-sp.log(sp.Rational(2, 3))) == 0, 'Pair survival rate failed')

# Native SM coefficient and archived crossings, scoped to one-loop QCD.
p159, tree159 = source('verification/v159_pyrate_gauge_crosscheck.py')
native159 = p159.read_text()
require('(sp.Rational(41, 10), sp.Rational(-19, 6), sp.Integer(-7))' in native159, 'Native one-loop coefficient changed')
c3 = 1/(8*sp.pi)
phi = 1/(6*sp.pi)
beta = lambda value: -7*value**2/(2*sp.pi)
require(beta(c3) != 0 and beta(phi) != 0, 'Unexpected one-loop QCD fixed point')
arch = ROOT / '_archive/tfpt-45/technical_companion.tex'
archtext = arch.read_text()
require('Strong-coupling crossing at the seed value' in archtext
        and 'Strong-coupling crossing at the seam constant' in archtext,
        'Archived terminology changed')
used[str(arch.relative_to(ROOT))] = hashlib.sha256(arch.read_bytes()).hexdigest()

# Exact original Cl16 NS half-shift implementer; no new spin representation.
ce = {'Fr': Fr}
functions('verification/v506_seam_clock_rigidity.py',
          ['gam', 'mono_mul', 'cmul', 'cadd', 'cscale', 'shift_matrix', 'apply_lin', 'implements'], ce)
gam, mul, add, scale = (ce[n] for n in ['gam', 'cmul', 'cadd', 'cscale'])
pspin, tspin = source('verification/v506_seam_clock_rigidity.py')
native_part_b = next(n for n in tspin.body if isinstance(n, ast.FunctionDef) and n.name == 'part_b')
native_ut = next(n for n in native_part_b.body if isinstance(n, ast.Assign)
                 and any(isinstance(t, ast.Name) and t.id == 'Ut' for t in n.targets))
native_loop = next(n for n in native_part_b.body if isinstance(n, ast.For)
                   and any(isinstance(v, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'Ut' for t in v.targets) for v in n.body))
ce.update({'n': 16, 'ONE': gam()})
exec(compile(ast.Module(body=[native_ut, native_loop], type_ignores=[]), str(pspin), 'exec'), ce)
Ut = ce['Ut']
U = scale(Ut, Fr(1, 16))
parity = gam(*range(16))
D = sp.Matrix(ce['shift_matrix'](16, 8, -1))
C = sp.Matrix(ce['shift_matrix'](16, 4, -1))
require(D*D == -sp.eye(16) and C*C == D and C**4 == -sp.eye(16), 'Native NS shift relation failed')
require(mul(U, U) == parity and mul(parity, parity) == gam(), 'Spin lift square failed')
require(ce['implements'](U, D.tolist(), 16), 'Native lift does not implement NS half-shift')
require(mul(U, parity) == mul(parity, U), 'Native deck lift changed parity')
Pplus = scale(add(gam(), parity), Fr(1, 2))
Pminus = scale(add(gam(), scale(parity, Fr(-1))), Fr(1, 2))
require(not mul(mul(Pminus, U), Pplus), 'Native lift exchanges parity sectors')

out = {
 'contract': 'seam-readout-flow-deck-20261005',
 'research_id': 'UR.SEAM.READOUT_FLOW_DECK.01',
 'verdict': 'PARTIAL',
 'readout': {
   'native_Y': str(Y.T), 'Y_euclidean_norm_squared': str(ynorm),
   'T3_relative_embedding_index': '5/3',
   'one_spinor_dimension': len(carrier), 'spinor_sum_Y_squared': str(sum(q*q for q in charges)),
   'charge_multiplicities': {str(q): charges.count(q) for q in sorted(set(charges))},
   'carrier_vertex_conformal_weight': '5/8', 'carrier_vertex_two_point_exponent': '5/4',
   'full_glue_current_weight': '1', 'affine_current_two_point_exponent': '2',
   'normalized_JJ_with_fixed_vertex_pair': str(normalized),
   'connected_JJ_with_fixed_vertex_pair': str(connected),
   'common_fractional_vertex_amplitude_survives_normalization': False,
   'native_joint_response_gram': str(gram), 'native_joint_response_gram_determinant': str(gram.det()),
   'native_topological_R_current_response': str(charge_response['R_JY_JY']),
   'native_charge_weighted_RY_current_response': str(charge_response['RY_JY_JY']),
   'native_uniform_scalar_port': '4 P0',
   'zero_sum_family_projection_annihilates_uniform_scalar_port': True,
   'native_matching_epsilon': str(mark_port['matching_epsilon']),
   'matching_epsilon_status': 'Necessary target-matching equality, explicitly not a source selection of the physical contact strength.',
   'scope': 'Canonical level-one lattice current algebra, fixed charged pair and its free-boson Ward identity. A defect sector sum or a different physical seam response may differ and must be selected by the original source.'
 },
 'flow': {
   'native_transfer': str(M), 'eigenvalues': {str(k): v for k, v in M.eigenvals().items()},
   'unit_eigenspace_dimension': len(fixed), 'fixed_vectors': [str(v) for v in fixed],
   'pair_transfer': str(B), 'pair_logarithm': str(logB),
   'pair_transfer_commutes_sheet_swap': True, 'coordinate_hop': '1/6',
   'parity_changes_in_one_transfer_step': False,
   'phi_base_over_c3': str(sp.simplify(phi/c3)),
   'one_loop_QCD_beta_c3': str(beta(c3)), 'one_loop_QCD_beta_phi_base': str(beta(phi)),
   'archived_RG_values_are_crossings': True,
   'scope': 'Native finite transfer with rates selected from its target spectrum, not a derived universal RG. The QCD beta test is restricted to the declared SM one-loop coefficient; no two-loop crossing scale was reproduced.'
 },
 'deck': {
   'native_one_particle_dimension': 16, 'native_spin_lift_terms': len(U),
   'geometric_deck_order': 2, 'NS_one_particle_half_shift_square': '-I',
   'NS_one_particle_quarter_shift_fourth_power': '-I',
   'fock_lift_square': 'fermion parity', 'fock_lift_fourth_power': 'I',
   'fock_lift_commutes_fermion_parity': True,
   'P_minus_U_P_plus': '0',
   'relative_offset_delta_P_minus_is_forbidden_by_this_lift': False,
   'scope': 'The existing NS/bounding spin lift on Cl16. Geometric sheets, these Fock parity sectors, and D5/A3 chirality modules require separate explicit identifications; no identification is assumed.'
 },
 'excluded': [
   'The carrier 5/4 two-point exponent automatically remains in the normalized/connected U(1) current response of one canonical fixed charged sector.',
   'The native v486 transfer has two independent stationary modes corresponding to c3 and phi_base.',
   'The archived QCD crossings are zeros of the declared one-loop QCD beta function.',
   'The native v506 NS Fock deck lift exchanges its two fermion-parity sectors and hence fixes a relative Hamiltonian offset by that exchange.'
 ],
 'first_missing_original_choice': 'A common source-selected charged seam state, observable/defect measure and time action. These must precede interpreting the relative 4/3 factor and carrier 5/4 exponent as physical alpha data.',
 'sources': used,
 'firewall': 'Internal exact structural tests, not empirical evidence and not promotion of a physical TFPT claim.'
}
for rel, digest in used.items():
    require(hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == digest, 'Source changed during evaluation: '+rel)
print(json.dumps(out, indent=2))
