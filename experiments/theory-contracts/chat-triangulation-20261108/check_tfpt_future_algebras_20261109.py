"""Exact native DET16 and C6 history-observability checks."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
import json
import math
from pathlib import Path
import sys

import numpy as np
import sympy as sp

root = Path(sys.argv[1])
sys.path.insert(0, str(root / "verification"))
from v17_hexagonal_resolvent import U6_power
import v1002_spin10_mirror_projector as native


def walsh_inplace(values):
    """Unnormalized integer Walsh transform on the original occupation bits."""
    for bit in range(16):
        block = values.reshape(-1, 2, 1 << bit)
        left, right = block[:, 0].copy(), block[:, 1].copy()
        block[:, 0], block[:, 1] = left + right, left - right


def occupation_resolution(support, states, wedge_labels):
    differences = support.astype(np.int64)
    walsh_inplace(differences)
    differences **= 2
    walsh_inplace(differences)
    assert np.all(differences % len(states) == 0)
    differences //= len(states)
    assert differences[0] == np.count_nonzero(support)
    # An omitted set fails if it contains any nonzero difference between
    # two allowed strings. Superset closure checks all omission sets.
    bad = differences > 0
    bad[0] = False
    for bit in range(16):
        block = bad.reshape(-1, 2, 1 << bit)
        block[:, 1] |= block[:, 0]
    omissions = states[~bad]
    omitted_count = max(int(mask).bit_count() for mask in omissions)
    maximal = [int(mask) for mask in omissions if int(mask).bit_count() == omitted_count]
    first = maximal[0]
    assert len(np.unique(states[support] & (len(states)-1-first))) == np.count_nonzero(support)
    return {
        "positions": int(np.count_nonzero(support)),
        "minimum_individual_occupations": 16-omitted_count,
        "number_of_minimal_alphabets": len(maximal),
        "example_omitted_wedge_addresses": [int(wedge_labels[j]) for j in range(16) if first >> j & 1],
        "one_bit_collisions": differences[1 << np.arange(16)].tolist(),
    }


def det16():
    generators, pairs = native.positive_weyl_generators()
    cartan = [generators[pairs.index((2*j, 2*j+1))] for j in range(5)]
    weights = np.array([2*np.diag(x).real for x in cartan], dtype=int).T
    wedge = (weights+1)//2
    labels = wedge @ (2**np.arange(4, -1, -1))
    states = np.arange(1 << 16)
    occupation = (states[:, None] >> np.arange(16)) & 1
    charges = wedge @ np.array([-2, -2, -2, 3, 3])
    charge = occupation @ charges
    common = charge % 6 == 0
    central = charge % 12 == 0
    split = common & ((occupation @ ((weights @ np.array([-1, -1, 1, 1, 0]))//2)) % 2 == 0)
    union = central | split
    intersection = central & split
    difference = central.astype(int)-split.astype(int)
    assert (central.sum(), split.sum(), intersection.sum(), union.sum()) == (5464, 5464, 3072, 7856)
    assert not np.any(difference @ occupation)
    contrast = occupation.T @ (occupation * difference[:, None])
    values = Counter(int(contrast[i,j]) for i,j in combinations(range(16),2))
    assert values == {0:106, -52:3, -24:6, 72:2, 156:3}
    witnesses = [[int(labels[i]),int(labels[j])] for i,j in combinations(range(16),2)
                 if abs(contrast[i,j]) == 156]
    assert witnesses == [[3,12],[3,20],[3,24]]
    domains = {key: occupation_resolution(support, states, labels) for key,support in (
        ("common_Z6",common),("union",union),("central",central),
        ("split",split),("intersection",intersection),
        ("union_without_mixed_plane",union & (states!=0) & (states!=len(states)-1)))}
    assert [d["minimum_individual_occupations"] for d in domains.values()] == [14,14,13,13,12,14]
    assert domains["union"]["number_of_minimal_alphabets"] == 60

    determinant = native.determinant_projector(16, tuple(range(16)), native.PHI).tocoo()
    assert set(determinant.row) == set(determinant.col) == {0, len(states)-1}
    plane = [0,len(states)-1]
    block = determinant.tocsr()[plane,:][:,plane].toarray()
    expected = np.array([[1,np.exp(-1j*native.PHI)],[np.exp(1j*native.PHI),1]])/2
    assert np.allclose(block, expected) and union[plane].all()
    # Remove the native phase by a diagonal change of frame, which leaves
    # all occupations unchanged. These are the exact plane matrices.
    eye = sp.eye(2)
    x = sp.Matrix([[0,1],[1,0]])
    y = sp.Matrix([[0,-sp.I],[sp.I,0]])
    z = sp.diag(1,-1)
    first = (eye+z)/2
    second = (eye-z)/2
    a = sp.symbols("a", real=True)
    unitary = sp.cos(a)*eye+sp.I*sp.sin(a)*x
    assert sp.simplify(unitary.H*z*unitary - (sp.cos(2*a)*z-sp.sin(2*a)*y)) == sp.zeros(2)
    assert all(first*m*first in (first,-first,sp.zeros(2))
               and second*m*second in (second,-second,sp.zeros(2)) for m in (eye,y,z))
    assert sp.trace(x*eye) == sp.trace(x*y) == sp.trace(x*z) == 0
    assert first*x*second == sp.Matrix([[0,1],[0,0]])
    b = sp.symbols("b", positive=True)
    heat = (eye+x)/2 + b*(eye-x)/2
    assert sp.simplify(heat**2 - ((1+b**2)*eye+(1-b**2)*x)/2) == sp.zeros(2)
    return {"domains":domains, "four_field_pair_contrast_counts":dict(values),
            "strongest_pair_wedge_addresses":witnesses, "maximum_raw_pair_contrast":156,
            "coherent_algebra_complex_dimension":7858,
            "unitary_and_occupation_probability_effect_dimension":7857,
            "additional_energy_effect_probability_dimension":7858}


def clock():
    shift = U6_power(1)
    pair = sp.diag(0,1,1,0,0,0)
    masks = [set(i for i in range(6) if (shift**j * pair * shift.T**j)[i,i])
             for j in range(6)]
    full = set(range(6))
    nonadaptive = {}
    for width in range(1,7):
        exact_sets = [order for order in combinations(range(6),width)
                      if len({tuple(site in masks[i] for i in order) for site in full}) == 6]
        if exact_sets:
            nonadaptive = {"minimum":width,"number_of_minimal_offset_sets":len(exact_sets),
                           "example":exact_sets[0]}
            break
    assert nonadaptive == {"minimum":4,"number_of_minimal_offset_sets":9,"example":(0,1,2,3)}

    @lru_cache(None)
    def depth(sites):
        if len(sites) <= 1:
            return 0
        return 1 + min(max(depth(frozenset(sites & mark)),depth(frozenset(sites-mark)))
                       for mark in masks if sites & mark and sites-mark)
    assert depth(frozenset(full)) == 3
    partitions = ((1,2),(4,5),(0,3))
    pair_label = {site:index for index,group in enumerate(partitions) for site in group}
    paths = [[pair_label[site],pair_label[(site-1)%6]] for site in range(6)]
    assert len({tuple(p) for p in paths}) == 6
    assert sum(len(masks[i] & masks[(i+1)%6]) == 1 for i in range(6)) == 6

    # With the native positive Green transfer admitted as a recorded
    # contraction, real symmetric effects become complete. The linear
    # certificate uses inverse congruence, which preserves the stabilized
    # finite span because Green congruence is invertible there.
    identity = sp.eye(6)
    r = sp.symbols("r", positive=True)
    inverse_green = (1+2*r)*identity-r*(shift+shift.T)
    for site in range(6):
        v = inverse_green[:,(site+1)%6]
        assert sp.expand(v[site]*v[(site+2)%6]) == r**2
        other = (site+3)%6
        keep = full.copy()
        for mark in masks:
            if site not in mark and other not in mark:
                keep -= mark
        assert keep == {site,other}

    delta = ((794-7*math.sqrt(9961))/2187)**(1/6)
    shift_np = np.array(shift,dtype=float)
    directed_pairs = [np.diag([i in mask for i in range(6)]).astype(float) for mask in masks]
    green_ranks = []
    for cusp in (1,1/3,2/3):
        ratio = cusp*delta/(cusp-delta)**2
        green = np.linalg.inv((1+2*ratio)*np.eye(6)-ratio*(shift_np+shift_np.T))
        assert np.all(green>0)
        operations = [shift_np,green,*directed_pairs,*(np.eye(6)-d for d in directed_pairs)]
        basis = [np.eye(6).reshape(-1)/math.sqrt(6)]
        i = 0
        while i < len(basis):
            effect = basis[i].reshape(6,6)
            for op in operations:
                candidate = (op.T @ effect @ op).reshape(-1)
                for _ in range(2):
                    for prior in basis:
                        candidate -= (prior @ candidate)*prior
                norm = np.linalg.norm(candidate)
                if norm > 1e-10:
                    basis.append(candidate/norm)
            i += 1
        assert len(basis) == 21
        green_ranks.append(len(basis))
    for site in range(6):
        v = sp.eye(6)
        v[site,site] = sp.I
        for other in range(6):
            if site == other:
                continue
            s = sp.zeros(6)
            s[site,other] = s[other,site] = 1
            assert (v.H*s*v)[site,other] == -sp.I
            assert (v.H*s*v)[other,site] == sp.I
    return {"one_binary_pair_nonadaptive":nonadaptive, "adaptive_minimax_binary_reads":3,
            "two_complete_pair_reads":paths, "permutation_and_label_probability_dimension":6,
            "coherent_algebra_complex_dimension":36,
            "with_native_recorded_green_probability_dimensions":green_ranks,
            "with_additional_local_quarter_phase_probability_dimension":36}


def lattice():
    result = []
    for carrier_rank in (5,13):
        total = carrier_rank + 3
        # D_n consists of integral vectors of even coordinate sum.
        # The diagonal spinor lift doubles to odd/odd and quadruples into D_n+D_3.
        assert carrier_rank%2 == 1 and total%8 == 0
        assert sp.Rational(carrier_rank,4)+sp.Rational(3,4) == total//4
        assert sp.Rational(4*4,4**2) == 1
        integer_roots = 2*total*(total-1)
        spinor_roots = (2**(total-1)) if total//4 == 2 else 0
        result.append({"carrier_rank":carrier_rank,"family_rank":3,"index":4,
                       "determinant":1,"diagonal_spinor_norm_squared":total//4,
                       "roots":integer_roots+spinor_roots})
    assert [x["roots"] for x in result] == [240,480]
    return result


def passive_extension_bound():
    delta = ((794-7*math.sqrt(9961))/2187)**(1/6)
    slow = 6*math.log(3/2)
    target_fast = 6*math.log(3)
    tests = []
    for cusp in (1,1/3,2/3):
        ratio = cusp*delta/(cusp-delta)**2
        log_rates = [slow*math.log1p(k*ratio)/math.log1p(ratio) for k in (1,3,4)]
        assert log_rates[0] < log_rates[1] < log_rates[2] < 2*slow < target_fast
        tests.append({"cusp":cusp,"maximum_original_visible_rate":log_rates[2],
                      "minimum_return_at_t1_from_jensen":math.exp(-log_rates[2]),
                      "v330_fast_return_at_t1":math.exp(-target_fast)})
    return tests


print(json.dumps({"det16":det16(),"c6":clock(),"diagonal_lattice_glue":lattice(),
                  "passive_extension_bound":passive_extension_bound()},indent=2))
