"""Exact old Ramond moments and a scoped time-sewing test against the native Pin Green source."""
from contextlib import redirect_stdout
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np
from scipy.integrate import quad
from scipy.sparse import coo_matrix

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(base / "check_tfpt_ramond_covariance_and_cover_20261110.py"))
    native = runpy.run_path(str(base / "check_tfpt_conditioned_r113_weights_20261110.py"))

Q, m = old["Q"], old["m"]
pairs = list(combinations(range(88), 2))
pair_id = {pair: index for index, pair in enumerate(pairs)}


def exact_word(state, word):
    out = {}
    for bits, coefficient in state.items():
        for i, creation in reversed(word):
            if bool(bits & (1 << i)) == creation:
                break
            if (bits & ((1 << i) - 1)).bit_count() % 2:
                coefficient = -coefficient
            bits ^= 1 << i
        else:
            out[bits] = coefficient
    return out


def inner(left, right):
    return sum((m.conj(a) * right.get(bits, Q.zero)
                for bits, a in left.items()), Q.zero)


def exact_text(value):
    return str(Q.to_sympy(value))


def q_packet(quarter, index, state, cov):
    probabilities = {}
    norm = inner(state, state)
    for bits, amplitude in state.items():
        exact = m.conj(amplitude) * amplitude / norm
        assert not exact.y
        probabilities[bits] = Fraction(str(exact.x))
    p = [Fraction(x) for x in cov["N_diagonal_by_channel_then_mode"]]
    variable = [i for i, x in enumerate(p) if 0 < x < 1]
    n = {bits.bit_count() for bits in state}
    assert len(n) == 1 and sum(p) == next(iter(n))
    pair_p = {pair: Fraction(0) for pair in combinations(variable, 2)}
    for bits, weight in probabilities.items():
        for pair in combinations([i for i in variable if bits >> i & 1], 2):
            pair_p[pair] += weight
    connected = {pair: value - p[pair[0]] * p[pair[1]]
                 for pair, value in pair_p.items()}
    wick_variance = sum(x * (1 - x) for x in p)
    assert 2 * sum(connected.values()) == -wick_variance
    for i in variable:
        assert sum(value for pair, value in connected.items() if i in pair) == -p[i]*(1-p[i])

    # All normal-ordered four-point values are a Gram of exact CAR pair
    # annihilations. The archived sparse state below is its exact factorization.
    row_id = {}
    rows, cols, coeff = [], [], []
    for bits, amplitude in state.items():
        occupied = [i for i in range(88) if bits >> i & 1]
        value = complex(float(amplitude.x), float(amplitude.y)) / np.sqrt(float(norm.x))
        for k, i in enumerate(occupied):
            for l in range(k + 1, len(occupied)):
                j = occupied[l]
                reduced = bits ^ (1 << i) ^ (1 << j)
                rows.append(row_id.setdefault(reduced, len(row_id)))
                cols.append(pair_id[i, j])
                coeff.append(value * (-1)**(k + l - 1))
    factors = coo_matrix((coeff, (rows, cols)), shape=(len(row_id), len(pairs))).tocsr()
    two_rdm = (factors.conj().T @ factors).tocoo()
    np.testing.assert_allclose(two_rdm.diagonal().sum(), next(iter(n))*(next(iter(n))-1)/2)
    # Offdiagonal pair entries are phase witnesses. Removing coherences
    # leaves every occupation moment unchanged and all these entries zero.
    offdiag = np.flatnonzero(two_rdm.row != two_rdm.col)
    largest = offdiag[np.argmax(np.abs(two_rdm.data[offdiag]))]
    i, j = pairs[two_rdm.row[largest]]
    k, l = pairs[two_rdm.col[largest]]
    word = ((i, True), (j, True), (l, False), (k, False))
    value = inner(state, exact_word(state, word)) / norm
    np.testing.assert_allclose(complex(float(value.x), float(value.y)), two_rdm.data[largest], atol=1e-13)
    assert len({i, j, k, l}) == 4
    assert sorted((i//11, j//11)) == sorted((k//11, l//11))
    assert i%11 + j%11 == k%11 + l%11

    # Four neutral halfwords on the marked old packet. Each preserves
    # charge and grade, so exact old time multiplies the whole Gram.
    a, b = min(connected, key=connected.get)
    c = next(i for i in variable if i not in (a, b)
             and 0 < sum(w for bits, w in probabilities.items()
                         if all(bits >> x & 1 for x in (a, b, i))) < pair_p[a, b])
    keys = [(), (a,), (a, b), (a, b, c)]
    moments = [[sum(w for bits, w in probabilities.items() if
                    all(bits >> x & 1 for x in set(u + v))) for v in keys] for u in keys]
    evals = np.linalg.eigvalsh(np.asarray(moments, dtype=float))
    assert evals.min() > 0
    time_points = np.array([.04, .10, .16])
    time_kernel = np.outer(np.exp(-32*time_points), np.exp(-32*time_points))
    timed = np.kron(time_kernel, np.asarray(moments, dtype=float))
    assert np.linalg.matrix_rank(timed, tol=1e-12) == 4

    filled, empty = min(((i,j) for i in range(88) for j in range(i+1,88)
                         if i//11 == j//11 and p[i] == 1 and p[j] == 0),
                        key=lambda ij: ij[1]-ij[0])
    current = ((empty,True),(filled,False))
    adjoint = ((filled,True),(empty,False))
    moved = exact_word(state, current)
    assert inner(moved, moved) == norm and not exact_word(state, adjoint)
    return {
        "quarter": quarter, "packet": index, "N": next(iter(n)),
        "wick_variance": str(wick_variance), "twice_sum_connected": str(2*sum(connected.values())),
        "connected_pairs_negative": sum(x < 0 for x in connected.values()),
        "connected_pairs_positive": sum(x > 0 for x in connected.values()),
        "connected_pairs_zero": sum(x == 0 for x in connected.values()),
        "min_connected": float(min(connected.values())), "max_connected": float(max(connected.values())),
        "nonzero_two_rdm_entries_abs_gt_1e-12": int(sum(abs(x)>1e-12 for x in two_rdm.data)),
        "phase_witness_modes_create_annihilate": [[i,j],[k,l]],
        "phase_witness_exact": exact_text(value), "phase_witness_abs": float(np.hypot(float(value.x),float(value.y))),
        "neutral_words_modes": keys, "equal_cut_exact_gram": [[str(x) for x in row] for row in moments],
        "equal_cut_gram_min_eigenvalue": float(evals.min()),
        "three_time_gram_rank": 4,
        "unmarked_KMS_even_current_counterexample": {
            "exact_empty_mode":empty,"exact_filled_mode":filled,
            "charge_neutral_old_H_energy_transfer":8*(empty-filled),
            "expect_B_dagger_B":1,"expect_B_B_dagger":0},
        "exact_state_coefficients_bithex": {format(bits,"022x"):[str(x.x),str(x.y)] for bits,x in sorted(state.items())},
        "state_norm_squared": exact_text(norm),
    }


covariances = old["cases"]
packets = []
joint = {}
quarter_states = {}
for alpha in (old["sp"].Rational(1,4), old["sp"].Rational(3,4)):
    delta = old["sp"].Matrix([(j+alpha)/8-old["sp"].Rational(1,2) for j in range(8)])
    drift = [Q.from_sympy(x) for x in old["cl"](old["V"].T*delta)]
    prepared = [old["old"]["series_response"](f,drift,4) for f in old["polynomials"]]
    states = [old["as_car"](old["vac_r"], old["old"]["vertex_zero"](polynomial))
              for polynomial in prepared]
    quarter_states[str(alpha)] = states
    batch = [q_packet(str(alpha), index, state, covariances[str(alpha)]["q_output"][index])
             for index, state in enumerate(states)]
    packets.extend(batch)
    marked = []
    keys = batch[0]["neutral_words_modes"]
    (i,j),(k,l) = batch[0]["phase_witness_modes_create_annihilate"]
    for state in states:
        marked.extend({bits:a for bits,a in state.items() if all(bits>>x&1 for x in word)}
                      for word in keys)
        marked.append(exact_word(state, ((i, True), (j, True), (l, False), (k, False))))
    gram = old["sp"].Matrix(len(marked),len(marked),
                            lambda r,c: Q.to_sympy(inner(marked[r],marked[c])))
    assert gram == gram.H
    rank = gram.to_DM().convert_to(Q).rank()
    joint[str(alpha)] = {"all_five_same_words_exact_gram_rank":rank,
                         "words_per_packet":5, "number_word_modes":keys,
                         "phase_sensitive_word_modes":[[i,j],[k,l]],
                         "all_five_same_words_exact_gram":[[str(x) for x in row] for row in gram.tolist()]}

low, high = quarter_states["1/4"], quarter_states["3/4"]
same_vectors = [i for i, (left,right) in enumerate(zip(low,high)) if left == right]
assert same_vectors == [0,3]
sp = old["sp"]
intervals = [sp.Rational(3,100),sp.Rational(4,100),sp.Rational(2,100)]
# In the same native channel and mode coordinates, H_high-H_low=Q_total/2.
# The original charged intermediate word B^* B with two removals visits
# total charges -1,-2,-1 even though the endpoints are neutral.
half_removal = ((1,False),(0,False))
assert inner(exact_word(low[0], half_removal), exact_word(low[0], half_removal)) == inner(low[0],low[0])
assert inner(exact_word(high[0], half_removal), exact_word(high[0], half_removal)) == inner(high[0],high[0])
energy_exponents = {}
for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
    omega_i,omega_j = [8*(n+alpha/8-sp.Rational(1,2)) for n in (-5,-4)]
    predicted = (32-omega_i,32-omega_i-omega_j)
    delta = [(j+alpha)/8-sp.Rational(1,2) for j in range(8)]
    # Recompute the parent energy from the charge sector and the
    # oscillator grade relative to each new charged sea.
    for bits in list(quarter_states[str(alpha)][0])[:8]:
        for actual, expected in ((bits, 32), (bits^1,predicted[0]), (bits^3,predicted[1])):
            energy = 0
            for channel in range(8):
                count = sum(actual>>old["bit"](channel,n)&1 for n in old["levels"])
                charge = sp.Rational(count)-sp.Rational(11,2)
                grade = sum(n*((actual>>old["bit"](channel,n)&1)-int(bool(n<charge)))
                            for n in old["levels"])
                energy += 8*(charge**2/2+delta[channel]*charge+grade)
            assert energy == expected
    energy_exponents[str(alpha)] = -sum(t*energy for t,energy in zip(intervals,
        (*predicted,predicted[0]), strict=True))
exponent_ratio = energy_exponents["3/4"]-energy_exponents["1/4"]
assert exponent_ratio == sp.Rational(13,200)
quarter_detection = {"exact_equal_packets_in_common_old_mode_coordinates":same_vectors,
                     "two_removed_modes":[0,1],
                     "intermediate_old_total_charges":[-1,-2,-1],
                     "time_intervals_old_H":[str(t) for t in intervals],
                     "absolute_log_returns":{q:str(x) for q,x in energy_exponents.items()},
                     "high_over_low_exact":"exp("+str(exponent_ratio)+")",
                     "high_over_low":float(sp.exp(exponent_ratio))}

# This is the first additional premise of the cross-source comparison.
# Match fractions of the two given Euclidean periods, beta_q=4*pi, beta_pin=2.
# The original q exponent is 32 in beta_q units; the matching x exponent is 64*pi.
old_energy_in_pin_x = 32 * (4*np.pi)/2
fields = []
absolute_spectra = []
for radius in (0., .55):
    phi = native["majorana"] + radius * native["inner"](np.array([0.,1.]))
    energies, basis = np.linalg.eigh(phi)
    residual = (np.linalg.norm(phi-(basis*energies)@basis.conj().T) +
                2*np.linalg.norm(basis.conj().T@basis-np.eye(len(basis)))*max(abs(energies)))
    assert residual < 1e-12
    absolute_spectra.append(abs(energies))
    d = .05
    z = np.exp(-d*old_energy_in_pin_x)
    null_polynomial_at_q = (z-np.exp(d*energies))*(z-np.exp(-d*energies))
    assert min(abs(null_polynomial_at_q)) > .99
    fields.append({"radius": radius, "normal_zero_abs_K_max": float(max(abs(energies))),
                   "distance_q_from_signed_normal_zero_spectrum": float(min(abs(np.r_[energies,-energies]-old_energy_in_pin_x))),
                   "three_time_step_pin_x": d,
                   "three_time_null_polynomial_min_at_q": float(min(abs(null_polynomial_at_q))),
                   "spectral_reconstruction_residual":float(residual),
                   "forced_mixed_block_rank_if_common_translation": 0})
shared_rate_gap = min(abs(absolute_spectra[0][:,None]-absolute_spectra[1][None,:]).ravel())
assert shared_rate_gap > .026

# Two actual metrics used for the previous Pin Higgs probe. On the NS
# double cover Vol = 2 * integral_{-1}^{1} exp(2 epsilon h(w)) dw.
def bump(w):
    x = (abs(w)-.5)/.25
    return np.exp(1-1/(1-x*x)) if abs(x)<1 else 0.


def volume(epsilon):
    return 4 + 4*quad(lambda w: np.expm1(2*epsilon*bump(w)), .25, .75, epsabs=1e-13)[0]


v0, v1 = volume(0.), volume(.1)
assert v1 > v0
print(json.dumps({"packets": packets, "joint_exact_old_q_grams":joint,
                  "quarter_detection":quarter_detection,
                  "marked_q_times_old_H": [.04,.10,.16],
                  "same_fraction_of_original_periods_comparison": {
                      "q_energy_old_H":32,"q_period":4*np.pi,"pin_period":2,
                      "q_energy_in_pin_x":old_energy_in_pin_x,"static_Higgs":fields,
                      "minimum_Higgs_atom_gap_between_static_radii":float(shared_rate_gap)},
                  "cover_volume_same_pin_metrics": {
                      "flat":v0,"epsilon_point_one":v1,
                      "shift_per_c_in_log_Z1_over_Z0_fixed_reference":-(v1-v0)}},
                 indent=2))
