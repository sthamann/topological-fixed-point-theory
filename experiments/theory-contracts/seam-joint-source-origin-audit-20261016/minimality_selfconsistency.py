"""Source-scoped minimum-rank wedge test; not a physical source construction.

The same volume tensor on Lambda^2 E x Lambda^2 E x E is restricted to
each two-block split of E=C5. No SM modules or fitted Yukawa are imported.
"""
from itertools import combinations
from pathlib import Path
from fractions import Fraction
import hashlib
import json


def wedge_sign(indices):
    if len(set(indices)) != len(indices):
        return 0
    return (-1)**sum(indices[i] > indices[j]
                     for i in range(len(indices))
                     for j in range(i+1, len(indices)))


def action(generator, basis, indices):
    """Exterior action, with columns indexed by the input basis element."""
    lookup = {x: i for i, x in enumerate(basis)}
    output = {}
    for pos, index in enumerate(indices):
        for row, col, coefficient in generator:
            if col != index:
                continue
            changed = list(indices)
            changed[pos] = row
            sign = wedge_sign(changed)
            if sign:
                target = lookup[tuple(sorted(changed))]
                output[target] = output.get(target, 0)+coefficient*sign
    return output


def generators(b, s):
    result = []
    for block in (range(b), range(b, b+s)):
        block = list(block)
        result += [[(i, j, 1)] for i in block for j in block if i != j]
        result += [[(block[i], block[i], 1),
                    (block[i+1], block[i+1], -1)]
                   for i in range(len(block)-1)]
    result.append([(i, i, -s if i < b else b) for i in range(b+s)])
    return result


def split_certificate(b):
    s = 5-b
    v = list(range(b)); w = list(range(b, 5))
    if b == 1:
        left = [(i, j) for i in v for j in w]
        right = list(combinations(w, 2))
    elif b == 2:
        left = list(combinations(v, 2))
        right = list(combinations(w, 2))
    elif b == 3:
        left = [(i, j) for i in v for j in w]
        right = list(combinations(v, 2))
    else:
        left = right = list(combinations(v, 2))
    boson = [(i,) for i in w]
    tensor = {(i, j, k): wedge_sign(alpha+beta+h)
              for i, alpha in enumerate(left)
              for j, beta in enumerate(right)
              for k, h in enumerate(boson)}
    checks = 0
    for gen in generators(b, s):
        for (i, j, k), value in tensor.items():
            residual = sum(c*tensor[ii, j, k]
                           for ii, c in action(gen, left, left[i]).items())
            residual += sum(c*tensor[i, jj, k]
                            for jj, c in action(gen, right, right[j]).items())
            residual += sum(c*tensor[i, j, kk]
                            for kk, c in action(gen, boson, boson[k]).items())
            assert residual == 0, (b, gen, i, j, k)
            checks += 1
    nonzero = sum(x != 0 for x in tensor.values())
    assert nonzero > 0
    charges = [-s]*b+[b]*s
    pieces = [(p, 2-p) for p in range(3)
              if p <= b and 2-p <= s]
    channels = [[list(x), list(y)]
                for i, x in enumerate(pieces) for y in pieces[i:]
                if x[0]+y[0] == b and x[1]+y[1]+1 == s]
    even_weights = [sum(charges[j] for j in subset)
                    for degree in (0, 2, 4)
                    for subset in combinations(range(5), degree)]
    assert len(even_weights) == 16
    assert sum(even_weights) == 0
    assert sum(q**3 for q in even_weights) == 0
    return {'negative_rank': b, 'positive_rank': s,
            'fermion_exterior_degrees': [2, 2],
            'boson_exterior_degree': 1,
            'leg_dimensions': [len(left), len(right), len(boson)],
            'nonzero_volume_tensor_entries': nonzero,
            'complexified_lie_generators': len(generators(b, s)),
            'exact_invariance_coefficients_checked': checks,
            'even_packet_dimension': len(even_weights),
            'trace_X_on_even_packet': sum(even_weights),
            'trace_X_cubed_on_even_packet': sum(q**3 for q in even_weights),
            'trace_X_squared_on_even_packet': sum(q*q for q in even_weights),
            'conditional_normalized_Y_roots': [str(-Fraction(1, b)), str(Fraction(1, s))],
            'conditional_B_times_gamma': str(Fraction(5, s*s)),
            'unordered_primitive_degree_two_leg_channels': channels,
            'primitive_channel_count': len(channels),
            'positive_minimum_carrier_degree': 5}


def smaller_rank_signatures(g):
    """Complete finite degree list in the stated holomorphic exterior class."""
    hits = []
    for b in range(1, g):
        s = g-b
        legs = [(p, q) for p in range(b+1) for q in range(s+1)
                if (p+q) >= 2 and (p+q) % 2 == 0]
        for p, q in legs:
            for r, t in legs:
                # Necessary center condition: (P,Q)=(b*m,s*m).
                if (p+r) % b == 0:
                    m = (p+r)//b
                    if m > 0 and q+t+1 == s*m:
                        hits.append([b, s, p, q, r, t, m])
    return hits


def main():
    controls = {str(g): smaller_rank_signatures(g) for g in (2, 3, 4)}
    assert all(not values for values in controls.values())
    certificates = [split_certificate(b) for b in range(1, 5)]
    assert [x['primitive_channel_count'] for x in certificates] == [1, 2, 1, 1]
    rank_balance = [[b, 5-b] for b in range(1, 5)
                    if 3*2**4 == 4*(b*b+(5-b)**2-1)]
    assert rank_balance == [[2, 3], [3, 2]]
    root = Path(__file__).resolve().parents[3]
    sources = [
        '_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex',
        '_archive/tfpt-45/source_extracts/02_carrier_source.tex',
        'tfpt_1_architecture_e8.tex',
        'tfpt_3_e8_audit_bootstrap.tex',
        'origin_theory.tex',
        'verification/v14_carrier_uniqueness.py',
    ]
    print(json.dumps({
        'contract': 'UR.SEAM.JOINT_SOURCE_ORIGIN_AUDIT.01',
        'annex': 'minimality-selfconsistency-20261016',
        'verdict': 'MINIMUM_RANK_PRIMITIVITY_DOES_NOT_SELECT_3PLUS2',
        'source_sha256': {p: hashlib.sha256((root/p).read_bytes()).hexdigest()
                          for p in sources},
        'same_tensor_for_all_splits': 'vol(alpha wedge beta wedge h)',
        'smaller_rank_necessary_center_signatures': controls,
        'minimum_rank_certificates': certificates,
        'conditional_original_occupancy_equation_solutions': rank_balance,
        'bosonic_positive_block_rank_preselected': False,
        'negative_rank_preselected': False,
        'fitted_coupling_selected': False,
        'canonical_volume_tensor_is_actual_source_vertex': False,
        'full_minimal_defect_fiber_classified': False,
        'all_original_selfconsistency_conditions_countermodelled': False,
        'physical_subalgebra_and_sector_action_selected': False,
        'full_p1_p2_origin_derived': False,
    }, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
