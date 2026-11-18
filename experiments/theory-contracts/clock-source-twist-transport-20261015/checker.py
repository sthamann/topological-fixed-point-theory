#!/usr/bin/env python3
"""Exact covariance/lift regressions on existing source data; no P1/P2 promotion."""
import argparse
import ast
from functools import cache
from fractions import Fraction as F
from itertools import combinations, product
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as S


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--repo', type=Path, default=HERE.parents[2])
parser.add_argument('--out', type=Path, default=HERE/'results.json')
args = parser.parse_args()
repo = args.repo.resolve()
checks = {}


def check(name, condition):
    if not bool(condition):
        raise AssertionError(name)
    checks[name] = True


manifest = json.loads((HERE/'source_manifest.json').read_text())
for item in manifest['sources']:
    digest = hashlib.sha256((repo/item['path']).read_bytes()).hexdigest()
    check('source_hash:'+item['path'], digest == item['sha256'])


def definitions(path, names, namespace):
    """Load only the cited original functions, without running their checkers."""
    tree = ast.parse((repo/path).read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    check('original_definitions:'+path, {n.name for n in nodes} == set(names))
    exec(compile(ast.Module(body=nodes, type_ignores=[]), path, 'exec'), namespace)
    return namespace


root_path = 'experiments/theory-contracts/seam-joint-ward-state-20261014/checker.py'
root_ns = definitions(root_path, {'roots', 'in_lattice'},
                      {'F': F, 'combinations': combinations, 'product': product})
roots = [S.Matrix(v) for v in root_ns['roots']('s')]
j_ns = definitions('verification/v634_st31_structure.py', {'J_vec'}, {'np': np})
J = S.Matrix.hstack(*[S.Matrix(j_ns['J_vec'](np.array(list(S.eye(8)[:,i]), dtype=object)))
                       for i in range(8)])
q = S.ones(8, 1)
s = q/2
source_path = 'experiments/theory-contracts/seam-joint-ward-state-20261014/operator_transport_check.py'
ns = definitions(source_path, {'eps_car', 'split_even_lattice', 'eps_e8'},
                 {'S': S, 'q': q, 's': s})
eps_original = ns['eps_e8']

# These are the exact original basis assignments; do not execute the source
# checker's unrelated QWZ imports or its argument parser.
tree = ast.parse((repo/source_path).read_text())
assigns = [n for n in tree.body if isinstance(n, ast.Assign)
           and any(isinstance(t, ast.Name) and t.id in {'e', 'basis'} for t in n.targets)]
bns = {'S': S}
exec(compile(ast.Module(body=assigns, type_ignores=[]), source_path, 'exec'), bns)
B = bns['basis']
Binv = B.inv()


@cache
def eps(u, v):
    return int(eps_original(S.Matrix(u), S.Matrix(v)))


@cache
def ordered_phase(n, rotated):
    """Phase of ordered generator powers, in the ACTUAL eps_e8 cocycle."""
    basis = J*B if rotated else B
    point = S.zeros(8, 1)
    phase = 1
    for i in range(8):
        k = int(n[i])
        b = basis[:, i]
        if k < 0:
            unit = eps(tuple(b), tuple(-b))  # E_b^{-1}=eps(b,-b) E_-b.
            b = -b
        else:
            unit = 1
        for _ in range(abs(k)):
            phase *= unit*eps(tuple(point), tuple(b))
            point += b
    return phase


@cache
def eta(v):
    n = Binv*S.Matrix(v)
    assert all(x.q == 1 for x in n)
    n = tuple(n)
    return ordered_phase(n, True)*ordered_phase(n, False)


check('J_orthogonal', J.T*J == S.eye(8))
check('J_square_minus_identity', J**2 == -S.eye(8))
check('J_fourth_identity', J**4 == S.eye(8))
check('original_E8_basis_unimodular', abs(B.det()) == 1)
check('original_E8_basis_even_integral', all(x.q == 1 for x in B.T*B)
      and all(int((B.T*B)[i,i]) % 2 == 0 for i in range(8)))
check('J_preserves_full_lattice_basis', all(x.q == 1 for x in Binv*J*B))
check('J_preserves_240_original_roots', {tuple(J*v) for v in roots} == {tuple(v) for v in roots})
check('source_cocycle_basis_commutator', all(
    eps(tuple(B[:,i]), tuple(B[:,j]))*eps(tuple(B[:,j]), tuple(B[:,i]))
    == (-1)**int(B[:,i].dot(B[:,j]) % 2) for i in range(8) for j in range(8)))

# The all-charge proof is in PROOF.md: this generator assignment defines a
# central-extension automorphism. Its fourth power is a lattice character,
# so checking the eight basis generators proves its fourth power everywhere.
fourth_character = []
for i in range(8):
    v = B[:,i]
    phase = 1
    for _ in range(4):
        phase *= eta(tuple(v))
        v = J*v
    fourth_character.append(phase)
check('actual_cocycle_lift_fourth_character_on_basis', fourth_character == [1]*8)
orbit_generators = [J**a*B[:,i] for a in range(4) for i in range(8)]
check('actual_cocycle_equivariance_1024_generator_pairs', all(
    eta(tuple(u+v))*eps(tuple(u),tuple(v))
    == eta(tuple(u))*eta(tuple(v))*eps(tuple(J*u),tuple(J*v))
    for u in orbit_generators for v in orbit_generators))
check('vacuum_lift_phase_is_one', eta(tuple(S.zeros(8,1))) == 1)

records = []
p = S.Matrix(S.symbols('p0:8', real=True))
r0 = S.Matrix([S.Rational(1,2)]*4+[-S.Rational(1,2)]*4)
for alpha in [S.Rational(1,4), S.Rational(3,4)]:
    delta = S.Matrix([(k+alpha)/8-S.Rational(1,2) for k in range(8)])
    ds = [J**a*delta for a in range(4)]
    qs = [J**a*q for a in range(4)]
    key = str(alpha)
    check(key+':four_distinct_twists', len({tuple(d) for d in ds}) == 4)
    check(key+':second_twist_is_minus_delta', ds[2] == -delta)
    check(key+':fixed_sector_forced_shift_not_in_Ls',
          not root_ns['in_lattice'](tuple(J*delta-delta), 's'))
    check(key+':twice_delta_not_in_Ls', not root_ns['in_lattice'](tuple(2*delta), 's'))
    check(key+':all_complex_J_affine_obstruction_identity',
          (S.eye(8)+J)*(J*delta-delta) == -2*delta)
    check(key+':four_step_translation_telescope',
          sum((ds[(a+1)%4]-ds[a] for a in range(4)),S.zeros(8,1)) == S.zeros(8,1))
    for a in range(4):
        d, next_d = ds[a], ds[(a+1)%4]
        E = lambda x,d: x.dot(x)/2+d.dot(x)
        check(key+f':sector{a}:full_quadratic_time_covariance',
              S.expand(E(J*p,next_d)-E(p,d)) == 0)
        check(key+f':sector{a}:charge_covariance_requires_rotated_q',
              S.expand(qs[(a+1)%4].dot(J*p)-qs[a].dot(p)) == 0)
        check(key+f':sector{a}:same_vacuum_and_rotated_extra_zero',
              E(S.zeros(8,1),d) == 0 and S.simplify(E(J**a*r0,d)) == 0)
        for offset in range(8):
            nu = S.Matrix([S.Rational(k+offset,8) for k in range(8)])
            tau = nu+alpha*q/8
            ta, tb = J**a*tau, J**((a+1)%4)*tau
            check(key+f':sector{a}:offset{offset}:full_uncompensated_transition',
                  all(S.expand(tb.dot(J*B[:,i])-ta.dot(B[:,i])) == 0 for i in range(8)))
    # Equality of rational exponents above retains the full character, not
    # merely its adjoint or a compensated periodic surrogate.
    check(key+':quarter_dressing_cannot_be_dropped', any(
        (alpha*q.dot(v)/8).q != 1 for v in roots))
    records.append({'alpha':str(alpha),'delta':[str(x) for x in delta],
                    'forced_fixed_sector_shift':[str(x) for x in J*delta-delta],
                    'fixed_sector_shift_in_Ls':False,'twice_delta_in_Ls':False,
                    'twist_orbit':[[str(x) for x in d] for d in ds]})

check('original_fixed_q_is_not_J_invariant', J*q != q)
check('original_fixed_q_failure_is_operator_not_constant',
      S.expand(q.dot(J*p)-q.dot(p)) != 0)
result = {
 'contract':'clock-source-twist-transport-20261015',
 'research_id':'UR.CLOCK.SOURCE_TWIST_TRANSPORT.01',
 'verdict':'PARTIAL',
 'mathematical_verdict':'PASS_SCOPED_Z4_COVARIANT_TWIST_FAMILY_AND_FIXED_SECTOR_OBSTRUCTION',
 'check_count':len(checks),'failed':0,'checks':checks,
 'actual_cocycle_fourth_character_on_basis':fourth_character,
 'cocycle_generator_pair_regressions':len(orbit_generators)**2,
 'quarter_cases':records,'source_hashes':manifest['sources'],
 'original_fixed_q_is_preserved_by_J':False,
 'covariant_q_required':'q_a=J^a q',
 'common_fock_premise':'Stipulated even E8 lattice Fock representation M(1) tensor C_eps[Ls]; no primitive-source selection',
 'ground_state_statement':'Reuse source full-lattice theorem K_delta>=0 with two ground vectors for 0<alpha<1; transport unitarily to every delta_a. The finite root checks do not prove the unbounded-lattice theorem.',
 'physical_gates_closed':[],'complete_TFPT_solution':False,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'proof_sha256':hashlib.sha256((HERE/'PROOF.md').read_bytes()).hexdigest() if (HERE/'PROOF.md').exists() else None,
 'not_claimed':['No source-selected physical twist sectors or sector-switch rate',
                'No fixed-sector J commutation with source K or fixed Q',
                'No P1/P2 selection of the E8 completion, original state or event instrument',
                'No local charged-field, odd-CAR or raw Xi dictionary',
                'No new verification/ledger/paper promotion']}
args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(f'{len(checks)}/{len(checks)} exact checks passed; failed=0')
print('Verdict: '+result['mathematical_verdict']+'; overall PARTIAL')


if __name__ == '__main__':
    pass
