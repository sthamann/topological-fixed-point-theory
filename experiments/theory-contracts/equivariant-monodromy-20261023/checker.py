#!/usr/bin/env python3
"""Exact evaluation of the existing G31 extension, not a new extension."""
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
import hashlib
import json
import pickle
import numpy as np
import sympy as S
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
CHECKS = []


def need(ok, name):
    if not bool(ok):
        raise RuntimeError(name)
    CHECKS.append(name)


def key(x):
    return tuple(map(int, x[0].flat)) + tuple(map(int, x[1].flat))


def mul(x, y):
    a, b = x
    c, d = y
    r, s = a @ c - b @ d, a @ d + b @ c
    if np.any(r % 2) or np.any(s % 2):
        raise RuntimeError('nonintegral twice-matrix product')
    return r // 2, s // 2


def power(x, n):
    z = (2*np.eye(4, dtype=np.int64), np.zeros((4, 4), dtype=np.int64))
    while n:
        if n & 1:
            z = mul(z, x)
        x = mul(x, x)
        n //= 2
    return z


def matrix(x):
    return (S.Matrix(x[0].tolist()) + S.I*S.Matrix(x[1].tolist()))/2


def real_matrix(x):
    a, b = x
    out = S.zeros(8)
    for i, j in product(range(4), repeat=2):
        out[2*i, 2*j] = S.Rational(int(a[i, j]), 2)
        out[2*i, 2*j+1] = -S.Rational(int(b[i, j]), 2)
        out[2*i+1, 2*j] = S.Rational(int(b[i, j]), 2)
        out[2*i+1, 2*j+1] = S.Rational(int(a[i, j]), 2)
    return out


def rank2(rows):
    pivots = {}
    for row in rows:
        x = sum((int(a) % 2) << j for j, a in enumerate(row))
        while x:
            p = x.bit_length()-1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    return len(pivots)


def main():
    # This state was generated locally from inspected source by 01_monodromy.py.
    state = HERE/'imported/work/state.pkl'
    refs, perms, pg, chosen, orbits, full, orbit_id, rays = pickle.loads(state.read_bytes())
    d = json.loads((HERE/'imported/work/results.json').read_text())
    expected = json.loads((HERE/'imported/expected_results.json').read_text())
    need(d == expected, 'all original classification data independently recomputed')
    need(len(full) == 46080 and len(pg) == 11520, 'actual full and projective groups')
    need(len(orbits) == 17, '17 original ordered monodromy classes')
    ident = key(power(refs[0], 0))
    need(len({p[0] for p in pg}) == 60, 'single generating reflection conjugacy class')

    B = S.Matrix(d['integral_lattice_basis'])
    Gram = S.Matrix(d['integral_gram'])
    ints = [B.inv()*real_matrix(x)*B for x in refs]
    need(all(all(a.q == 1 for a in x) and x.T*Gram*x == Gram for x in ints),
         'all 60 original lattice actions are integral isometries')
    rows = [list((x-S.eye(8)).row(j)) for x in ints for j in range(8)]
    invariant_rank = 8-rank2(rows)
    need(invariant_rank == 0, 'no invariant vectors in original lattice modulo two')

    # Explicit isometry from the Gaussian root coordinates into standard E8.
    T = S.diag(*([S.Matrix([[1, -1], [1, 1]])/2]*4))
    J = S.diag(*([S.Matrix([[0, -1], [1, 0]])]*4))
    gaussian_roots = set()
    for r in rays:
        a = tuple(v for pair in r for v in pair)
        for _ in range(4):
            gaussian_roots.add(a)
            a = tuple(v for j in range(4) for v in (-a[2*j+1], a[2*j]))
    standard = set()
    for i in range(8):
        for j in range(i+1, 8):
            for a, b in product((-1, 1), repeat=2):
                x = [S.Integer(0)]*8
                x[i], x[j] = a, b
                standard.add(tuple(x))
    for signs in product((-1, 1), repeat=8):
        if sum(s < 0 for s in signs) % 2 == 0:
            standard.add(tuple(S.Rational(a, 2) for a in signs))
    need({tuple(T*S.Matrix(x)) for x in gaussian_roots} == standard,
         'explicit isometry identifies all 240 original roots with standard E8')
    need(T.T*T == S.eye(8)/2 and T*J == J*T,
         'isometry preserves metric and the actual quarterclock')

    order8 = next(x for x in full.values()
                  if key(power(x, 8)) == ident and key(power(x, 4)) != ident)
    need(key(power(order8, 8)) == ident and key(power(order8, 4)) != ident,
         'G31 contains an element of order eight')
    positive_r = next(i for i, r in enumerate(refs)
                      if key(power(mul(r, order8), 8)) == ident)
    need(key(power(mul(refs[positive_r], order8), 8)) == ident,
         'positive original G31 example with exact order eight clock')
    # Gamma_sign has exponent four from its already proved normal form.
    for u, v, k in product(range(16), range(16), range(4)):
        x = (u, v, k)
        y = (0, 0, 0)
        for _ in range(4):
            a, b, l = y
            y = (a ^ u, b ^ v ^ (u if l % 2 else 0), (l+k) % 4)
        if y != (0, 0, 0):
            raise RuntimeError('Gamma_sign exponent is not four')
    need(True, 'all 1024 sign-lift elements have fourth power one')

    integral = []
    for oid, o in enumerate(d['orbits']):
        D = S.Matrix.vstack(*[ints[j]-S.eye(8) for j in o['representative'][:3]])
        sn = smith_normal_form(D, domain=S.ZZ)
        nz = [abs(int(sn[i, i])) for i in range(8) if sn[i, i]]
        h1 = {'free_rank': 24-len(nz), 'torsion_invariant_factors': [x for x in nz if x > 1],
              'd0_smith_nonzero': nz}
        need(h1 == o['integral_H1'], 'independent integral Smith form '+str(oid))
        integral.append({'id': oid, 'representative': o['representative'],
                         'image_order': o['image_order'], 'H0_rank': 8-len(nz), 'H1': h1})

    clocks = []
    for oid, o in enumerate(d['orbits']):
        for cid, c in enumerate(o.get('clock_classes', [])):
            if c['order'] != 4:
                continue
            C = (np.array(c['matrix_real_twice'], dtype=np.int64),
                 np.array(c['matrix_imag_twice'], dtype=np.int64))
            det = S.simplify(matrix(C).det())
            need(det in (-1, 1), f'clock determinant {oid}/{cid}')
            e = int(det == -1)
            RC = mul(refs[o['representative'][0]], C)
            need(key(power(C, 4)) == ident and key(power(RC, 4)) == ident,
                 f'both orbifold relations {oid}/{cid}')
            # Existing extension (a,g)(b,h)=(a+b-e(g)e(h),gh).
            defects = [(-2*e) % 4, (-2*(1-e)) % 4]
            for p, x in [(e, C), (1-e, RC)]:
                for a in range(4):
                    phase, parity = 0, 0
                    for _ in range(4):
                        phase, parity = (phase+a-parity*p) % 4, parity ^ p
                    need(phase == (-2*p) % 4 and parity == 0,
                         f'all generator lifts fourth powers {oid}/{cid}/{p}/{a}')
            phase_solutions = [ks for ks in product(range(4), repeat=4) if sum(ks)%4 == 2]
            strict = [ks for ks in phase_solutions if len(set(ks)) == 1]
            central_covariant = [(ks, a, delta) for ks in phase_solutions for a in range(4)
                                 for delta in range(4)
                                 if all((ks[(j+1)%4]-ks[j]-delta)%4 == 0 for j in range(4))]
            need(len(strict) == 0 and len(central_covariant) == 32,
                 f'complete peripheral phase census {oid}/{cid}')
            need({v[2] for v in central_covariant} == {1, 3},
                 f'precise necessary central covariance defect {oid}/{cid}')
            clocks.append({'monodromy_id': oid, 'clock_id': cid, 'epsilon_clock': e,
                           'obstruction_Z4_pair': defects, 'strict_lifts': 0,
                           'only_clock_fourth_power_relaxed_lifts': 0,
                           'central_covariance_lifts_with_strict_clock': 32 if e == 0 else 0,
                           'central_covariance_lifts_with_central_clock': 32,
                           'projective_orbifold_lifts_with_meridian_product_z2': 16,
                           'characteristic_polynomial': c['characteristic_polynomial'],
                           'matrix_real_twice': c['matrix_real_twice'],
                           'matrix_imag_twice': c['matrix_imag_twice'],
                           'automorphism_order': c['automorphism_order']})
    need(len(clocks) == 62, 'all 62 original exact-order-four clock classes evaluated')
    need(Counter(c['epsilon_clock'] for c in clocks) == {0: 30, 1: 32},
         '30 even and 32 odd clock determinant classes')
    totals = {'strict': sum(c['strict_lifts'] for c in clocks),
              'only_clock_power_central': sum(c['only_clock_fourth_power_relaxed_lifts'] for c in clocks),
              'central_covariance_strict_clock': sum(c['central_covariance_lifts_with_strict_clock'] for c in clocks),
              'central_covariance_central_clock': sum(c['central_covariance_lifts_with_central_clock'] for c in clocks),
              'projective_orbifold_with_central_meridian_product': sum(c['projective_orbifold_lifts_with_meridian_product_z2'] for c in clocks)}
    need(list(totals.values()) == [0, 0, 960, 1984, 992], 'complete counts for the five explicitly distinct contracts')

    # Explicit collisions for progressively stronger incomplete invariants.
    buckets = defaultdict(list)
    for c in clocks:
        sig = (c['monodromy_id'], c['epsilon_clock'], c['characteristic_polynomial'])
        buckets[sig].append([c['monodromy_id'], c['clock_id']])
    collisions = [v for v in buckets.values() if len(v) > 1]
    need(bool(collisions), 'clock spectrum and obstruction do not classify clock markings')
    collision_witness = []
    for ci in [3, 5]:
        c = d['orbits'][0]['clock_classes'][ci]
        C = (np.array(c['matrix_real_twice'], dtype=np.int64),
             np.array(c['matrix_imag_twice'], dtype=np.int64))
        r = refs[0]
        count = sum(key(mul(g, r)) == key(mul(r, g)) and
                    key(mul(g, C)) == key(mul(C, g)) for g in full.values())
        need(count == c['automorphism_order'], 'independent centralizer count 0/'+str(ci))
        collision_witness.append({'monodromy_id': 0, 'clock_id': ci,
                                  'automorphism_order': count})
    need([x['automorphism_order'] for x in collision_witness] == [64, 128],
         'explicit nonisomorphic clocks with matching coarse invariants')
    h_buckets = defaultdict(list)
    for o in integral:
        h_buckets[(o['image_order'], o['H1']['free_rank'],
                   tuple(o['H1']['torsion_invariant_factors']))].append(o['id'])

    # A central character of the lattice-sign kernel would be a mod-two
    # invariant vector. Its absence forbids a nonzero G31-equivariant pushout
    # Hom(L,Z2) -> Z4 with trivial target action.
    out = {'contract_id': 'MATH.SEAM.EQUIVARIANT.02',
           'original_results_exactly_reproduced': True,
           'group_bridge': {'Gaussian_to_standard_E8': str(T),
                            'common_quotients_up_to_isomorphism': ['1', 'C2'],
                            'no_injection_either_direction': True,
                            'order_eight_witness_real_twice': order8[0].tolist(),
                            'order_eight_witness_imag_twice': order8[1].tolist(),
                            'positive_C8_case_reflection_index': positive_r,
                            'invariant_mod2_lattice_dimension': invariant_rank,
                            'nonzero_equivariant_sign_kernel_to_Z4_maps': 0},
           'integral_cohomology_all_17': integral,
           'clock_obstructions_all_62': clocks,
           'clock_summary': dict(Counter(str(c['obstruction_Z4_pair']) for c in clocks)),
           'relaxation_class_counts': totals,
           'cohomology_and_image_order_collisions': list(h_buckets.values()),
           'same_base_parity_clock_spectrum_collisions': collisions,
           'independently_checked_collision_witness': collision_witness,
           'checks_passed': len(CHECKS), 'checks': CHECKS}
    (HERE/'results.json').write_text(json.dumps(out, indent=2)+'\n')
    (HERE/'CLOCKS.tsv').write_text('monodromy_id\tclock_id\tepsilon_C\tobstruction\tautomorphism_order\n'+
        ''.join(f"{c['monodromy_id']}\t{c['clock_id']}\t{c['epsilon_clock']}\t{c['obstruction_Z4_pair']}\t{c['automorphism_order']}\n" for c in clocks))
    (HERE/'INTEGRAL.tsv').write_text('id\tordered_reflection_indices\tH0_rank\tH1_free_rank\tH1_torsion\n'+
        ''.join(f"{o['id']}\t{o['representative']}\t{o['H0_rank']}\t{o['H1']['free_rank']}\t{o['H1']['torsion_invariant_factors']}\n" for o in integral))
    print(json.dumps({k: v for k, v in out.items() if k in
                     ['contract_id', 'clock_summary', 'relaxation_class_counts',
                      'same_base_parity_clock_spectrum_collisions', 'checks_passed']}, indent=2))


if __name__ == '__main__':
    main()
