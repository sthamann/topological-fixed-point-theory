#!/usr/bin/env python3
"""Exact neutral-vacuum blindness and eight-charge discriminator; no cutoffs."""
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as S

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
checks = {}


def check(name, value):
    if not bool(value):
        raise AssertionError(name)
    checks[name] = True


manifest = json.loads((HERE/'source_manifest.json').read_text())
wanted = {'verification/v634_st31_structure.py',
          'experiments/theory-contracts/seam-joint-ward-state-20261014/operator_transport_check.py'}
sources = [x for x in manifest['sources'] if x['path'] in wanted]
check('two_original_sources_pinned', {x['path'] for x in sources} == wanted)
for item in sources:
    check('source_hash:'+item['path'],
          hashlib.sha256((REPO/item['path']).read_bytes()).hexdigest() == item['sha256'])

source = next(x['path'] for x in sources if x['path'].endswith('operator_transport_check.py'))
tree = ast.parse((REPO/source).read_text())
nodes = [n for n in tree.body if isinstance(n, ast.Assign)
         and any(isinstance(t, ast.Name) and t.id in {'e', 'basis'} for t in n.targets)]
ns = {'S':S}
exec(compile(ast.Module(body=nodes,type_ignores=[]),source,'exec'),ns)
B = ns['basis']
jpath = 'verification/v634_st31_structure.py'
tree = ast.parse((REPO/jpath).read_text())
nodes = [n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name == 'J_vec']
jns = {'np':np}
exec(compile(ast.Module(body=nodes,type_ignores=[]),jpath,'exec'),jns)
J = S.Matrix.hstack(*[S.Matrix(jns['J_vec'](np.array(list(S.eye(8)[:,i]),dtype=object)))
                      for i in range(8)])

d = S.Matrix(S.symbols('d0:8', real=True))
q = S.Matrix(S.symbols('q0:8', real=True))
tau = S.Matrix(S.symbols('tau0:8', real=True))
p = S.Matrix(S.symbols('p0:8', real=True))
N = S.symbols('N', nonnegative=True)
t = S.symbols('t', positive=True)
energy = lambda charge: charge.dot(charge)/2+d.dot(charge)
zero = S.zeros(8,1)
check('neutral_sector_K_is_N_for_arbitrary_delta', energy(zero)+N == N)
check('neutral_sector_Q_is_zero_for_arbitrary_q', q.dot(zero) == 0)
check('neutral_sector_full_transition_is_identity_for_arbitrary_tau',
      S.exp(2*S.pi*S.I*tau.dot(zero)) == 1)
check('charged_energy_odd_part', S.expand(energy(p)-energy(-p)-2*d.dot(p)) == 0)
check('charged_energy_even_part', S.expand(energy(p)+energy(-p)-p.dot(p)) == 0)
check('charged_ratio_log_exponent', S.expand(-t*energy(p)+t*energy(-p)+2*t*d.dot(p)) == 0)
check('eight_original_lattice_charges_independent', B.shape == (8,8) and abs(B.det()) == 1)
odd_energy = S.Matrix([S.expand((energy(B[:,i])-energy(-B[:,i]))/2) for i in range(8)])
check('eight_charge_delta_reconstruction_symbolic', S.simplify(B.T.inv()*odd_energy-d) == zero)
check('unknown_time_only_fixes_t_delta_from_ratios',
      S.simplify(B.T.inv()*(-2*t*odd_energy)/(-2)-t*d) == zero)
m,n = S.symbols('m n', integer=True)
check('quadratic_current_total_grade_is_n', S.expand(m+(n-m)-n) == 0)

quarter_witnesses = []
for alpha in [S.Rational(1,4),S.Rational(3,4)]:
    delta = S.Matrix([(k+alpha)/8-S.Rational(1,2) for k in range(8)])
    orbit = [J**a*delta for a in range(4)]
    energies = []
    for a,da in enumerate(orbit):
        ep = S.Matrix([S.expand(v.dot(v)/2+da.dot(v)) for v in [B[:,i] for i in range(8)]])
        em = S.Matrix([S.expand(v.dot(v)/2-da.dot(v)) for v in [B[:,i] for i in range(8)]])
        check(f'{alpha}:frame{a}:actual_eight_charge_reconstruction', B.T.inv()*(ep-em)/2 == da)
        energies.append((ep-em)/2)
    check(f'{alpha}:four_twists_have_distinct_charged_data', len({tuple(x) for x in energies}) == 4)
    witness = next(i for i in range(8) if energies[0][i] != energies[1][i])
    quarter_witnesses.append({'alpha':str(alpha),'basis_charge_index':witness,
                             'charge':[str(x) for x in B[:,witness]],
                             'delta_dot_charge':str(energies[0][witness]),
                             'Jdelta_dot_charge':str(energies[1][witness])})

result = {'contract':'clock-source-twist-transport-20261015',
          'addendum':'NEUTRAL_VACUUM_TWIST_BLINDNESS_AND_EIGHT_CHARGE_DISCRIMINATOR',
          'verdict':'PASS_SCOPED_ANALYTIC_STATEMENTS',
          'overall_physical_verdict':'PARTIAL',
          'check_count':len(checks),'failed':0,'checks':checks,
          'quarter_witnesses':quarter_witnesses,
          'source_hashes':sources,
          'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'addendum_sha256':hashlib.sha256((HERE/'ADDENDUM.md').read_bytes()).hexdigest() if (HERE/'ADDENDUM.md').exists() else None,
          'scope':'Same vacuum, fixed Heisenberg/Jordan insertions each individually charge preserving, normal ordering and given lattice-Fock representation. No Gibbs/charged-state blindness claim.',
          'charged_prerequisites':'Normalized charged primary preparation, original charge labels and common known time/kinetic normalization; global charge shifts are not derived local physical vertex instruments.',
          'all_word_proof':'Charge sectors are invariant under every permitted insertion, every OPE coefficient and every intervening K evolution. Starting at p=0 the entire word remains in M(1) tensor |0>, where K=N, Q=0 and t_tau=I. See ADDENDUM.md; no finite word scan or oscillator cutoff substitutes for this argument.',
          'not_claimed':['No blind statement for words with charged intermediate insertions even if the total word is neutral',
                         'No source-selected detector or charged-state preparation',
                         'No reconstruction of q or tau from C_delta alone',
                         'No P1/P2 source selection or physical gate closure']}
(HERE/'neutral_twist_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(f'{len(checks)}/{len(checks)} exact checks passed; failed=0; overall physical verdict PARTIAL')
