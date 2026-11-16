"""Exact three-window order test for the original affine spatial generators."""

from itertools import product
from math import comb

import numpy as np

from reconstruction import require


def row_products(left, right):
    return (left[:, :, None]*right[:, None, :]).reshape(len(left), -1)


def three_window_jets(vectors, prime):
    """Derivatives at zero, using the Leibniz rule in the tensor algebra."""
    count = len(vectors)
    letters = [np.ones((count,1),dtype=np.int64), vectors]
    for degree in (2,3):
        letters.append(row_products(vectors,letters[degree-1]) % prime)
    prefix = letters
    for _ in range(2):
        prefix_count = len(prefix[0])
        prefix = [
            sum(comb(degree,late)*row_products(
                np.tile(letters[late],(prefix_count,1)),
                np.repeat(prefix[degree-late],count,axis=0))
                for late in range(degree+1)) % prime
            for degree in range(4)
        ]
    return prefix


def witness_rows(word_columns, prime):
    """Choose output coordinates sufficient to separate the column span."""
    reduced = word_columns.copy()
    pivots = []
    for row in range(len(reduced)):
        positions = np.flatnonzero(reduced[row])
        if not len(positions):
            continue
        index = int(positions[0])
        pivots.append(index)
        reduced[row] = reduced[row]*pow(int(reduced[row,index]),-1,prime) % prime
        reduced[row+1:] = (
            reduced[row+1:]-reduced[row+1:,index:index+1]*reduced[row]
        ) % prime
    return pivots


def direct_column_jets(generators, vectors, history, prime):
    """Independent matrix Leibniz replay of a selected history on e_0."""
    empty = np.zeros(generators[0].shape[0],dtype=np.int64)
    empty[0] = 1
    derivatives = [empty,*[np.zeros_like(empty) for _ in range(3)]]
    for index in history:
        field = sum(coefficient*generator for coefficient,generator in
                    zip(vectors[index],generators)) % prime
        powers = [derivatives]
        for _ in range(3):
            powers.append([(field@column) % prime for column in powers[-1]])
        derivatives = [
            sum(comb(degree,late)*powers[late][degree-late] for late in range(degree+1)) % prime
            for degree in range(4)
        ]
    return derivatives


def census(generators, vectors, profiles, mirrors, prime):
    """Partition all ordered pairs by jets, without enumerating the square."""
    count = len(profiles)
    histories = tuple(product(range(count),repeat=3))
    jets = three_window_jets(vectors,prime)
    column = np.zeros((generators[0].shape[0],1),dtype=np.int64)
    column[0,0] = 1
    observed, levels = [], []
    quadratic_collision = None
    for degree in (1,2,3):
        column = np.concatenate([(g@column) % prime for g in generators],axis=1)
        rows = witness_rows(column.T,prime)
        observed.append((jets[degree]@column[rows].T) % prime)
        signatures = np.hstack(observed)
        _, indices, inverse, counts = np.unique(
            signatures,axis=0,return_index=True,return_inverse=True,return_counts=True)
        levels.append({"derivative_order":degree,"matrix_rows_on_column_zero":rows,
                       "distinct_histories":len(counts),"candidate_pairs":int(counts@counts),
                       "largest_collision_class":int(counts.max())})
        if degree == 2 and counts.max() > 1:
            aliased = int(np.flatnonzero(counts>1)[0])
            first = int(indices[aliased])
            second = int(np.flatnonzero((inverse==aliased)&(np.arange(len(histories))!=first))[0])
            quadratic_collision = (first, second)
    require(len(counts)==len(histories), "the third jet did not certify all ordered histories")
    require(quadratic_collision is not None, "native negative control needs a quadratic collision")
    first,second = quadratic_collision
    history_a,history_b = histories[first],histories[second]
    direct_a,direct_b = (direct_column_jets(generators,vectors,h,prime) for h in (history_a,history_b))
    require(all(np.array_equal(direct_a[k],direct_b[k]) for k in (1,2)),
            "lower-order negative control did not replay")
    different = int(np.flatnonzero(observed[-1][first]!=observed[-1][second])[0])
    row = rows[different]
    values = [int(direct_a[3][row]),int(direct_b[3][row])]
    require(values == [int(observed[-1][i,different]) for i in quadratic_collision],
            "independent matrix replay disagrees with the cubic witness")
    return {
        "prime":prime,"left_histories":len(histories),"initial_candidate_pairs":len(histories)**2,
        "levels":levels,"unique_mirror_rule":"apply the existing profile mirror at each window in the same order",
        "equal_length_partner_count":len(histories),
        "individually_mirrored_reverse_order_partner_count":sum(i==k for i,_,k in histories),
        "quadratic_negative_control":{
            "left":[list(profiles[i]) for i in history_a],
            "other_left":[list(profiles[i]) for i in history_b],
            "matching_right":[list(profiles[mirrors[i]]) for i in history_a],
            "cubic_matrix_entry":[row,0],"cubic_residues":values},
        "scope":"Exact separation of analytic one-particle words in the declared three-window spatial menu.",
    }
