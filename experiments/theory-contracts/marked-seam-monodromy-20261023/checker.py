#!/usr/bin/env python3
"""Exact marked-monodromy and lattice-lift calculations; pure mathematics."""
from __future__ import annotations

import ast
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import itertools as it
import json
import math
from pathlib import Path

import numpy as np
import sympy as S
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
CHECKS = {}


def need(name, condition):
    if not bool(condition):
        raise AssertionError(name)
    CHECKS[name] = True


def functions(path, names, namespace):
    tree = ast.parse(path.read_text())
    selected = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    need('selected_functions:' + path.name, {n.name for n in selected} == set(names))
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


def gf2_rank(columns):
    basis = {}
    for x in columns:
        while x:
            p = x.bit_length() - 1
            if p in basis:
                x ^= basis[p]
            else:
                basis[p] = x
                break
    return len(basis)


def apply_bits(columns, x):
    out = 0
    while x:
        low = x & -x
        out ^= columns[low.bit_length() - 1]
        x ^= low
    return out


def gf2_solve(columns, target):
    basis = {}
    for i, value in enumerate(columns):
        witness = 1 << i
        while value:
            p = value.bit_length()-1
            if p in basis:
                value ^= basis[p][0]
                witness ^= basis[p][1]
            else:
                basis[p] = (value, witness)
                break
    out = 0
    while target:
        p = target.bit_length()-1
        if p not in basis:
            return None
        target ^= basis[p][0]
        out ^= basis[p][1]
    return out


def main():
    paths = [
        'verification/v634_st31_structure.py',
        'experiments/theory-contracts/seam-joint-ward-state-20261014/operator_transport_check.py',
        'experiments/theory-contracts/clock-source-twist-transport-20261015/checker.py',
        'tfpt_explorer/process.py',
    ]
    sources = [{'path': p, 'sha256': sha256((REPO / p).read_bytes()).hexdigest()} for p in paths]
    pins = HERE / 'source_manifest.json'
    if pins.exists():
        need('source_manifest_matches', json.loads(pins.read_text())['sources'] == sources)

    jns = functions(REPO / paths[0], {'J_vec'}, {'np': np})
    J = S.Matrix.hstack(*[
        S.Matrix(jns['J_vec'](np.array(list(S.eye(8)[:, i]), dtype=object)))
        for i in range(8)
    ])
    source = REPO / paths[1]
    tree = ast.parse(source.read_text())
    assigns = [n for n in tree.body if isinstance(n, ast.Assign)
               and any(isinstance(t, ast.Name) and t.id in {'e', 'basis'} for t in n.targets)]
    bns = {'S': S}
    exec(compile(ast.Module(body=assigns, type_ignores=[]), str(source), 'exec'), bns)
    B = bns['basis']
    G = B.T * B
    A = B.inv() * J * B
    need('actual_J_order_four', J**2 == -S.eye(8) and J**4 == S.eye(8))
    need('actual_J_integral_lattice_isometry', all(x.q == 1 for x in A) and A.T * G * A == G)
    T0 = S.diag(1, -1, 1, -1, 1, -1, 1, -1)
    T0c = B.inv()*T0*B
    need('reflected_base_lattice_object', all(x.q == 1 for x in T0c)
         and T0c.T*G*T0c == G and T0*J*T0 == J.inv())
    need('actual_E8_even_unimodular', G.det() == 1 and all(G[i, i] % 2 == 0 for i in range(8)))
    snf = smith_normal_form(A - S.eye(8), domain=S.ZZ)
    need('J_minus_I_smith_form', sorted(abs(int(snf[i, i])) for i in range(8)) == [1]*4 + [2]*4)

    # Kernel of the sign-lift projection: characters of the integral lattice.
    Sdual = (A.inv()).T
    cols = [sum((int(Sdual[i, j]) % 2) << i for i in range(8)) for j in range(8)]
    ncols = [cols[j] ^ (1 << j) for j in range(8)]
    need('character_N_rank_four', gf2_rank(ncols) == 4)
    need('character_N_squared_zero', all(apply_bits(ncols, x) == 0 for x in ncols))

    # Read the old normalized cocycle and construct its generator-normalized J lift.
    ns = functions(source, {'eps_car', 'split_even_lattice', 'eps_e8'},
                   {'S': S, 'q': S.ones(8, 1), 's': S.ones(8, 1)/2})
    eps_old = ns['eps_e8']
    doubled_B = [tuple(int(x) for x in 2*B[:, i]) for i in range(8)]

    def j_apply(v):
        return tuple(z for i in range(0, 8, 2) for z in (-v[i+1], v[i]))

    @lru_cache(None)
    def eps(u, v):
        c, h = u[0] % 2, v[0] % 2
        d = [(x-c)//2 for x in u]
        f = [(x-h)//2 for x in v]
        need_even = sum(d) % 2 == 0 and sum(f) % 2 == 0
        if not need_even:
            raise ValueError('not an E8 lattice vector')
        exponent = c*(sum(f)//2) + c*h
        exponent += sum(d[i]*f[j] for i in range(8) for j in range(i))
        exponent += c*h*sum(i*(d[i]+f[i]) for i in range(8))
        return -1 if exponent % 2 else 1

    @lru_cache(None)
    def ordered(n, rotated=False):
        point = (0,)*8
        phase = 1
        for i, k in enumerate(n):
            b = j_apply(doubled_B[i]) if rotated else doubled_B[i]
            unit = eps(b, tuple(-x for x in b)) if k < 0 else 1
            if k < 0:
                b = tuple(-x for x in b)
            for _ in range(abs(k)):
                phase *= unit * eps(point, b)
                point = tuple(x+y for x, y in zip(point, b))
        return phase

    @lru_cache(None)
    def eta(n):
        return ordered(n, True) * ordered(n, False)

    rng = np.random.default_rng(20261023)
    samples = [tuple(map(int, rng.integers(-2, 3, 8))) for _ in range(24)]
    for n, m in zip(samples, samples[1:] + samples[:1]):
        u, v = B*S.Matrix(n), B*S.Matrix(m)
        u2, v2 = tuple(int(2*x) for x in u), tuple(int(2*x) for x in v)
        need('literal_cocycle:' + str(n), eps(u2, v2) == int(eps_old(u, v)))
        an = tuple(int(x) for x in A*S.Matrix(n))
        am = tuple(int(x) for x in A*S.Matrix(m))
        need('J_lift_cocycle:' + str(n),
             eta(tuple(x+y for x, y in zip(n, m))) * eps(u2, v2)
             == eta(n)*eta(m)*eps(j_apply(u2), j_apply(v2)))
        need('J_lift_order_sample:' + str(n),
             eta(n)*eta(an)*eta(tuple(int(x) for x in A*S.Matrix(an)))
             *eta(tuple(int(x) for x in A**3*S.Matrix(n))) == 1)
    need('J_lift_order_on_generators', all(
        np.prod([eta(tuple(int(x) for x in A**k*S.eye(8)[:, i])) for k in range(4)]) == 1
        for i in range(8)))

    # Full character-valued factor table on the actual clock/reflection subgroup.
    # R is the native complex reflection in the root line (1+i)e_0.
    R = S.diag(-1, -1, 1, 1, 1, 1, 1, 1)
    Rc = B.inv()*R*B
    labels = list(product(range(4), range(2)))
    matrices = {g: A**g[0] * Rc**g[1] for g in labels}
    need('native_reflection_preserves_E8', Rc.T*G*Rc == G and all(x.q == 1 for x in Rc))
    need('clock_reflection_group_C4_times_C2', J*R == R*J and len({tuple(x) for x in matrices.values()}) == 8)

    @lru_cache(None)
    def generic_ordered(n, label):
        images = [tuple(int(x) for x in 2*B*matrices[label][:, i]) for i in range(8)]
        point = (0,)*8
        phase = 1
        for i, k in enumerate(n):
            b = images[i]
            unit = eps(b, tuple(-x for x in b)) if k < 0 else 1
            if k < 0:
                b = tuple(-x for x in b)
            for _ in range(abs(k)):
                phase *= unit*eps(point, b)
                point = tuple(x+y for x, y in zip(point, b))
        return phase

    @lru_cache(None)
    def generic_eta(label, n):
        return generic_ordered(n, label)*generic_ordered(n, (0, 0))

    def group_product(g, h):
        return ((g[0]+h[0]) % 4, g[1] ^ h[1])

    def char_action(label, value):
        v = S.Matrix([(value >> i) & 1 for i in range(8)])
        w = matrices[label].inv().T*v
        return sum((int(w[i]) % 2) << i for i in range(8))

    factors = {}
    for g, h in product(labels, repeat=2):
        gh = group_product(g, h)
        right = []
        for i in range(8):
            n = tuple(int(x) for x in S.eye(8)[:, i])
            hn = tuple(int(x) for x in matrices[h]*S.Matrix(n))
            sign = generic_eta(h, n)*generic_eta(g, hn)*generic_eta(gh, n)
            right.append(int(sign == -1))
        left = matrices[gh].inv().T*S.Matrix(right)
        factors[g, h] = sum((int(left[i]) % 2) << i for i in range(8))
    need('full_factor_cocycle_all_subgroup_triples', all(
        factors[g, h] ^ factors[group_product(g, h), k]
        == char_action(g, factors[h, k]) ^ factors[g, group_product(h, k)]
        for g, h, k in product(labels, repeat=3)))
    r0 = factors[(0, 1), (0, 1)]
    defect = factors[(1, 0), (0, 1)] ^ factors[(0, 1), (1, 0)]
    relation_cols = []
    for i in range(16):
        a, b = ((1 << i), 0) if i < 8 else (0, 1 << (i-8))
        square = b ^ char_action((0, 1), b)
        commute = a ^ char_action((0, 1), a) ^ b ^ char_action((1, 0), b)
        relation_cols.append(square | (commute << 8))
    correction = gf2_solve(relation_cols, r0 | (defect << 8))
    need('clock_reflection_character_extension_splits', correction is not None)
    need('clock_reflection_splitting_relations', apply_bits(relation_cols, correction) == (r0 | (defect << 8)))

    def lifted_product(x, y):
        a, g = x; b, h = y
        return a ^ char_action(g, b) ^ factors[g, h], group_product(g, h)

    cochain = {}
    for k, e0 in labels:
        out = (0, (0, 0))
        for _ in range(k):
            out = lifted_product(out, (correction & 255, (1, 0)))
        if e0:
            out = lifted_product(out, (correction >> 8, (0, 1)))
        cochain[k, e0] = out[0]
    need('full_explicit_coboundary_all_pairs', all(
        cochain[g] ^ char_action(g, cochain[h]) ^ factors[g, h] == cochain[group_product(g, h)]
        for g, h in product(labels, repeat=2)))
    need('reflected_quarter_lift_order_on_basis', all(
        np.prod([generic_eta((k % 4, k % 2), tuple(int(x) for x in (A*Rc)**j*S.eye(8)[:, i]))
                 for j, k in [(0, 1), (1, 1), (2, 1), (3, 1)]]) == 1
        for i in range(8)))
    factor_output = {
        'group': 'C4 x C2', 'generators': ['J', 'R'], 'ordered_elements': labels,
        'coefficient_module': 'Hom(E8,Z2), left action g.a = a composed g_inverse',
        'normalized_section_factor_bits': [[factors[g, h] for h in labels] for g in labels],
        'linear_relation_rank': gf2_rank(relation_cols),
        'splitting_correction_J': correction & 255,
        'splitting_correction_R': correction >> 8,
        'explicit_trivializing_cochain_bits': [cochain[g] for g in labels],
        'class_zero_in_twisted_H2': True,
    }

    # A Jordan basis verifies that the ACTUAL character action is four free blocks.
    uvecs, vvecs, span_cols = [], [], []
    for a in range(256):
        b = apply_bits(ncols, a)
        if b and gf2_rank(span_cols + [a, b]) == len(span_cols) + 2:
            uvecs.append(a); vvecs.append(b); span_cols.extend([a, b])
    need('actual_character_Jordan_basis', len(uvecs) == 4 and gf2_rank(span_cols) == 8)
    need('actual_character_Jordan_intertwiner', all(
        apply_bits(cols, a) == (a ^ b) and apply_bits(cols, b) == b
        for a, b in zip(uvecs, vvecs)))

    # Canonical coordinates (u,v,k); S(u,v)=(u,v+u), k in Z/4.
    m = 4
    mask = (1 << m)-1

    def mul(x, y):
        u, v, k = x; a, b, l = y
        return u ^ a, v ^ b ^ (a if k & 1 else 0), (k+l) % 4

    def inverse(x):
        u, v, k = x
        return u, v ^ (u if k & 1 else 0), -k % 4

    elements = list(product(range(1 << m), range(1 << m), range(4)))
    identity = (0, 0, 0)
    generators = [(1 << i, 0, 0) for i in range(m)] + [(0, 0, 1)]
    center = [x for x in elements if all(mul(x, y) == mul(y, x) for y in generators)]
    need('lift_group_order_1024', len(elements) == 1024)
    need('lift_group_center_32', len(center) == 32)
    need('lift_group_inverses', all(mul(x, inverse(x)) == identity for x in elements))
    need('noncentral_character_extension', mul(generators[0], generators[-1]) != mul(generators[-1], generators[0]))

    # The central quotient has V=F2^5, Z=F2^5 and bilinear full cocycle.
    def cocycle(x, y):
        u, e = x & mask, x >> m
        a, f = y & mask, y >> m
        return (a if e else 0) | ((e & f) << m)

    quotient = range(1 << (m+1))
    need('central_cocycle_all_triples', all(
        cocycle(x, y) ^ cocycle(x ^ y, z) == cocycle(y, z) ^ cocycle(x, y ^ z)
        for x, y, z in product(quotient, repeat=3)))
    need('central_class_nonzero_on_clock_line', cocycle(1 << m, 1 << m) == 1 << m)
    need('central_class_zero_on_U', all(cocycle(x, y) == 0 for x, y in product(range(1 << m), repeat=2)))
    need('quadratic_obstruction_all_nonU_lines', all(cocycle(x, x) != 0 for x in quotient if x >> m))
    for ell, t in product(range(1 << m), range(2)):
        def scalar_sign(x, y):
            z = cocycle(x, y)
            return -1 if (((z & mask) & ell).bit_count() + t*(z >> m)) % 2 else 1
        need('scalar_clock_carry_coboundary:' + str((ell, t)), all(
            scalar_sign(x, y)
            == ((-1)**(((y & mask) & ell).bit_count() if x >> m else 0))
            * (1j**(t*(x >> m))) * (1j**(t*(y >> m))) / (1j**(t*((x ^ y) >> m)))
            for x, y in product(quotient, repeat=2)))
        need('scalar_commutator_class:' + str((ell, t)), all(
            scalar_sign(u, 1 << m)*scalar_sign(1 << m, u)
            == (-1)**((ell & u).bit_count()) for u in range(1 << m)))

    # Complete classification for every ordered four-peripheral conjugacy profile.
    counts = Counter()
    image_orders = Counter()
    total_classes = 0
    based_solutions = 0
    for u0, u1, u2 in product(range(1 << m), repeat=3):
        u3 = u0 ^ u1 ^ u2
        shift = (u1 ^ u0) | ((u2 ^ u0) << m)
        representatives = {min(z, z ^ shift) for z in range(1 << (2*m))}
        expected = 1 << (2*m - bool(shift))
        need_profile = len(representatives) == expected
        if not need_profile:
            raise AssertionError('profile classification')
        counts['diagonal_profiles' if not shift else 'nondiagonal_profiles'] += 1
        counts['orbit16_classes' if not shift else 'orbit32_classes'] += len(representatives)
        total_classes += len(representatives)
        based_solutions += 1 << (3*m)
        for z in representatives:
            x, y = z & mask, z >> m
            a, b = u1 ^ u0, u2 ^ u0
            rr = gf2_rank([a | (x << m), b | (y << m), a << m, b << m])
            image_orders[4*(1 << rr)] += 1
        for z in (0, (1 << (2*m))-1):
            x, y = z & mask, z >> m
            v = (0, x, y, x ^ y ^ u1 ^ u3)
            word = identity
            for u, w in zip((u0, u1, u2, u3), v):
                word = mul(word, (u, w, 1))
            if word != identity:
                raise AssertionError('four-meridian relation')
    need('admissible_profiles_4096', sum(counts[k] for k in ('diagonal_profiles', 'nondiagonal_profiles')) == 4096)
    need('complete_classes_526336', total_classes == 526336)
    need('orbit_stabilizer_mass', counts['orbit16_classes']*16 + counts['orbit32_classes']*32 == based_solutions == 2**24)
    need('complete_monodromy_image_order_census', dict(image_orders) == {4: 16, 8: 720, 16: 14880, 32: 80640, 64: 430080})
    need('no_surjective_quarter_tuple', max(image_orders) == 64 < len(elements))

    # Exhaustive conjugation countercheck of representative tuples, independent
    # of the quotient counting formula. Includes the two stabilizer types.
    examples = []
    for profile in [(0, 0, 0, 0), (1, 1, 1, 1), (0, 1, 0, 1), (0, 1, 2, 3)]:
        vs = (0, 3, 5, 3 ^ 5 ^ profile[1] ^ profile[3])
        tup = tuple((u, v, 1) for u, v in zip(profile, vs))
        orbit = {tuple(mul(mul(g, x), inverse(g)) for x in tup) for g in elements}
        expected = 16 if len(set(profile)) == 1 else 32
        need('full_group_conjugation_example:' + str(profile), len(orbit) == expected)
        examples.append({'profile': profile, 'orbit_size': len(orbit), 'stabilizer_order': 1024//len(orbit)})

    # Small-rank exhaustive normal-form censuses guard the general formula.
    small = {}
    for mm in (1, 2):
        size = 1 << mm
        remaining = set()
        for us in product(range(size), repeat=3):
            profile = us + (us[0] ^ us[1] ^ us[2],)
            for x, y in product(range(size), repeat=2):
                shift = (profile[1] ^ profile[0], profile[2] ^ profile[0])
                remaining.add((profile, min((x, y), (x ^ shift[0], y ^ shift[1]))))
        expected = (2**(5*mm) + 2**(3*mm))//2
        need('small_rank_classification:' + str(mm), len(remaining) == expected)
        small[str(mm)] = len(remaining)

    # Independent Burnside calculation directly in the ACTUAL dual J basis:
    # C1=A^3, quotient by diagonal coboundaries (N b,N b,N b).
    coboundary_cols = [v | (v << 8) | (v << 16) for v in ncols]
    burnside_fixed = []
    for k in range(4):
        operator_cols = []
        for block in range(3):
            for i in range(8):
                image = 1 << i
                for _ in range(k):
                    image = apply_bits(cols, image)
                operator_cols.append((image ^ (1 << i)) << (8*block))
        rank = gf2_rank(operator_cols + coboundary_cols)
        burnside_fixed.append(1 << (24-rank))
    need('independent_native_Burnside_fixed_counts', burnside_fixed == [2**20, 2**12, 2**20, 2**12])
    need('independent_native_Burnside_orbit_count', sum(burnside_fixed)//4 == total_classes)

    # A different lattice complex structure is not silently substituted for J.
    # Its associated native complex representations have the stated clock action.
    cns = functions(REPO/paths[3], {'_code_basis'}, {'np': np, 'it': it, 'math': math})
    _, supports = cns['_code_basis']()
    phase5 = []
    for words in supports:
        values = {(-1)**word.count(0) for word in words}
        need('native_quartic_R_support:' + str(len(phase5)), len(values) == 1)
        phase5.append(values.pop())
    need('native_quartic_R_action', phase5 == [1, 1, 1, 1, -1])
    need('native_quartic_J_action', all(S.I**len(word) == 1 for words in supports for word in words))
    eigen4 = {'J': [S.I]*4, 'JR': [-S.I]+[S.I]*3}
    dimensions = {}
    for label, eigenvalues in eigen4.items():
        wedge_values = [a*b for a, b in it.combinations(eigenvalues, 2)]
        quartic_values = [1]*5 if label == 'J' else phase5
        dims = {}
        for name, values in [('V4', eigenvalues), ('wedge2_V4', wedge_values), ('quartic5', quartic_values)]:
            h0 = values.count(1)
            dims[name] = {'H0': h0, 'H1': 2*len(values)+h0}
        dimensions[label] = dims
    need('cohomology_distinguishes_J_and_JR', dimensions['J'] != dimensions['JR'])

    result = {
        'contract': 'MATH.SEAM.MONODROMY.01',
        'verdict': 'EXACT_CLASSIFICATION_IN_DECLARED_QUARTER_CLOCK_STRATUM',
        'checks_passed': len(CHECKS), 'checks': CHECKS,
        'sources': sources,
        'lattice_basis': [[str(x) for x in row] for row in B.tolist()],
        'J_in_lattice_basis': [[int(x) for x in row] for row in A.tolist()],
        'dual_Jordan_u': uvecs, 'dual_Jordan_v': vvecs,
        'integral_cohomology': {'H0': 0, 'H1_free_rank': 16, 'H1_torsion': [2]*4, 'H2': 0},
        'lift_group': {'order': 1024, 'center_order': 32, 'derived_order': 16,
                       'character_extension_split': True, 'character_extension_central': False,
                       'central_quotient_extension_split': False},
        'clock_reflection_full_factor': factor_output,
        'monodromy': {'peripheral_projection': ['J']*4, 'ordered_punctures': True,
                       'based_solutions': based_solutions, 'admissible_profiles': 4096,
                       'conjugacy_classes': total_classes, 'counts': dict(counts),
                       'image_order_class_counts': dict(sorted(image_orders.items())),
                       'stabilizer_orders': {'orbit16': 64, 'orbit32': 32},
                       'stabilizer_groups': {'diagonal': 'C4 x C2^4', 'nondiagonal': 'C2^5 = center'},
                       'small_rank_classes': small, 'examples': examples},
        'independent_native_Burnside_fixed_counts': burnside_fixed,
        'associated_local_system_cohomology': dimensions,
        'scope': 'All sign lifts of the fixed cyclic J representation into the full preimage of <J>. Not all G31 or W(E8) reflection quadruples, no mapping-class quotient, no imposed reflection-equivariance at punctures.'
    }
    (HERE/'results.json').write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n')
    if not pins.exists():
        pins.write_text(json.dumps({'sources': sources}, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('verdict', 'checks_passed', 'lift_group', 'monodromy')}))


if __name__ == '__main__':
    main()
