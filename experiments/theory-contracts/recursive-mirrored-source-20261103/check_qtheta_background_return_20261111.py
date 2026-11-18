"""Signed mode returns with all original Q background harmonics present."""
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np


def main(repo):
    package = repo / (
        "experiments/theory-contracts/original-paper-source-reconciliation-20261102"
        "/contribution_20261106Q"
    )
    sys.path.insert(0, str(package / "runtime"))
    from assemble import Assembler

    native_path = package / "domain/native-diagnostic-v1.json"
    coarse = Assembler(native_path)
    background = {"genes": [2, 2, 2]}
    native = coarse.assemble(background)
    domain = json.loads(native_path.read_text())
    domain["modes"] = list(range(-7, 7))
    domain["field_proposal"]["slots"].append(
        {"kind": "cos", "harmonic": 8, "component": 0, "unit": "13/100"}
    )
    with TemporaryDirectory() as directory:
        path = Path(directory) / "fourier-refinement.json"
        path.write_text(json.dumps(domain))
        assembler = Assembler(path)
        profiles = {
            a: assembler.assemble({"genes": background["genes"]+[a+1]})
            for a in (-1, 0, 1)
        }

    modes = domain["modes"]
    retained = coarse.domain["modes"]
    indices = np.concatenate([96*modes.index(n) + np.arange(96) for n in retained])
    h = {a: data["h_N"] for a, data in profiles.items()}
    for op in h.values():
        assert np.array_equal(op[np.ix_(indices, indices)], native["h_N"])

    columns = {a: np.eye(len(h[0]), dtype=complex)[:, indices] for a in h}
    powers = {a: [v[indices, :]] for a, v in columns.items()}
    for _ in range(6):
        for a in columns:
            columns[a] = h[a] @ columns[a]
            powers[a].append(columns[a][indices, :])
    odd = [(a-b)/2 for a, b in zip(powers[1], powers[-1])]
    for value in odd[:4]:
        assert not np.any(value)

    def full_mode(n):
        j = modes.index(n)*96
        return slice(j, j+96)

    def kept_mode(n):
        j = retained.index(n)*96
        return slice(j, j+96)

    a = 0.13
    b = h[1][full_mode(6), full_mode(-2)]/a
    d = h[0][full_mode(4), full_mode(6)]
    # Three in-cut paths with steps +8,-2,-2,-2 ending on 0.
    # The fourth ordering visits -8 and lies outside this refinement.
    d2 = d @ d
    predicted = d2 @ d @ b + d2 @ b @ d + d @ b @ d2
    ix = kept_mode(0), kept_mode(-2)
    check_error = np.linalg.norm(odd[4][ix]/a-predicted)
    assert check_error < 1e-11
    assert np.linalg.norm(predicted) > 1e-3

    native_rep = assembler.original
    omega = np.exp(2j*np.pi/3)
    left = native_rep.rep(1, np.eye(2), omega*np.eye(3))
    center = left @ assembler.K @ left.conj() @ assembler.K
    neutral = (np.eye(96) + center + center @ center)/3
    rank_neutral = int(round(np.trace(neutral).real))
    assert rank_neutral == 24
    np.testing.assert_allclose(center @ predicted, predicted @ center, atol=1e-12)
    np.testing.assert_allclose(assembler.gamma @ predicted, predicted @ assembler.gamma, atol=1e-12)
    color_readout = neutral @ predicted @ neutral
    assert np.linalg.norm(color_readout) > 1e-3

    def real_derivative(matrix):
        return (matrix - assembler.K @ matrix @ assembler.K)/1j

    six_y = real_derivative(native_rep.rep(3j, np.zeros((2, 2)), -1j*np.eye(3)))
    two_t3 = real_derivative(native_rep.rep(0, 1j*np.diag([1, -1]), np.zeros((3, 3))))
    charge = np.diag(six_y/6 + two_t3/2).real
    electromagnetic = np.exp(1j*np.pi*charge)
    half_turn = np.concatenate([(-1.0)**n*electromagnetic for n in modes])
    for sign in (-1, 0, 1):
        np.testing.assert_allclose(
            half_turn[:, None]*h[sign]*half_turn.conj()[None, :], h[-sign], atol=1e-13
        )
    coarse_turn = half_turn[indices]
    np.testing.assert_allclose(
        coarse_turn[:, None]*odd[4]*coarse_turn.conj()[None, :], -odd[4], atol=1e-12
    )
    fixed = np.diag(neutral).real > 0.5
    np.testing.assert_allclose(electromagnetic[fixed]**2, 1, atol=1e-13)
    assert np.count_nonzero(np.abs(color_readout) > 1e-11)
    np.testing.assert_allclose(
        electromagnetic[:, None]*color_readout*electromagnetic.conj()[None, :],
        -color_readout, atol=1e-13
    )
    # The existing homogeneous collar uses length 1/8 and the fixed
    # unperturbed Cauchy graph. Test the original neutral carrier there.
    neutral_slots = np.tile(fixed, len(retained))
    neutral_ix = np.ix_(neutral_slots, neutral_slots)
    neutral_h = native["h_N"][neutral_ix]
    mark = odd[4][neutral_ix]/a
    energy, basis = np.linalg.eigh(neutral_h)
    length = 1/8
    weights = 1/(1+np.exp(-2*length*energy))
    upper, lower = np.sqrt(weights), np.sqrt(1-weights)
    gram_multiplier = np.outer(upper, upper) + np.outer(lower, lower)
    mark_eigenbasis = basis.conj().T @ mark @ basis
    on_cauchy = basis @ (gram_multiplier*mark_eigenbasis) @ basis.conj().T
    opposite_faces = basis @ (
        (np.outer(upper, upper)-np.outer(lower, lower))*mark_eigenbasis
    ) @ basis.conj().T
    assert np.linalg.norm(on_cauchy) >= np.min(gram_multiplier)*np.linalg.norm(mark)
    cauchy_turn = coarse_turn[neutral_slots]
    np.testing.assert_allclose(
        cauchy_turn[:, None]*on_cauchy*cauchy_turn.conj()[None, :], -on_cauchy, atol=1e-12
    )
    even = cauchy_turn.real > 0
    np.testing.assert_allclose(on_cauchy[np.ix_(even, even)], 0, atol=1e-12)
    np.testing.assert_allclose(on_cauchy[np.ix_(~even, ~even)], 0, atol=1e-12)

    # Same fixed collar on the full refined neutral carrier. Keep the
    # virtual intermediate state for the full return before compression.
    fine_neutral = np.tile(fixed, len(modes))
    neutral_fine_ix = np.ix_(fine_neutral, fine_neutral)
    retained_mask = np.isin(np.repeat(modes, 96), retained)[fine_neutral]
    fine_h = h[0][neutral_fine_ix]
    fine_delta = (h[1]-h[0])[neutral_fine_ix]
    open_mark = (~retained_mask)[:, None]*fine_delta*retained_mask[None, :]
    fine_e, fine_basis = np.linalg.eigh(fine_h)
    fine_A = 1/(1+np.exp(-2*length*fine_e))
    sqrt_a, sqrt_b = np.sqrt(fine_A), np.sqrt(1-fine_A)
    open_in_basis = fine_basis.conj().T @ open_mark @ fine_basis
    upper_face = sqrt_a[:, None]*open_in_basis*sqrt_a[None, :]
    lower_face = sqrt_b[:, None]*open_in_basis*sqrt_b[None, :]
    selected_open = upper_face + lower_face
    return_raw = open_in_basis.conj().T @ open_in_basis
    full_return = return_raw*(np.outer(sqrt_a, sqrt_a)+np.outer(sqrt_b, sqrt_b))
    selected_return = selected_open.conj().T @ selected_open
    missed_open = (
        -sqrt_b[:, None]*open_in_basis*sqrt_a[None, :]
        +sqrt_a[:, None]*open_in_basis*sqrt_b[None, :]
    )
    np.testing.assert_allclose(
        full_return-selected_return, missed_open.conj().T @ missed_open,
        rtol=0, atol=1e-12
    )
    return_loss = np.trace(full_return-selected_return).real
    assert return_loss > 0.1

    # A fixed scalar fibre with these three amplitudes is already separated
    # by the identity, the second return, and the fourth return.
    second = (powers[1][2]+powers[-1][2])/2 - powers[0][2]
    assert np.trace(second).real > 0
    c = color_readout
    fourth = {k: np.vdot(c, powers[k][4][ix]).real for k in (-1, 0, 1)}
    signed_fourth = (fourth[1]-fourth[-1])/2
    np.testing.assert_allclose(signed_fourth, a*np.vdot(c, c).real, rtol=0, atol=1e-12)
    assert signed_fourth > 0

    # U=H*H is the original gauge-invariant Higgs norm. Read its real
    # Fourier coefficients directly by convolution of the field profile.
    density_marks = {}
    for sign, profile in profiles.items():
        field = profile["higgs_fourier"]
        complex_field = {n: v[[0, 2]]+1j*v[[1, 3]] for n, v in field.items()}
        def cosine_density(k):
            return 2*sum(np.vdot(v, complex_field[n+k])
                         for n, v in complex_field.items() if n+k in complex_field).real
        seventh, sixteenth = cosine_density(7), cosine_density(16)
        np.testing.assert_allclose([seventh, sixteenth], [sign*a*a, (sign*a)**2/2],
                                   atol=1e-15)
        density_marks[sign] = [float(seventh), float(sixteenth)]
    print(json.dumps({
        "background": [2, 2, 2],
        "coarse_operators_identical": True,
        "odd_moment_frobenius_norms_0_through_6": [float(np.linalg.norm(v)) for v in odd],
        "odd_coarse_diagonal_traces_0_through_6": [float(np.trace(v).real) for v in odd],
        "fourth_moment_per_amplitude_on_0_minus2": float(np.linalg.norm(predicted)),
        "fourth_moment_on_color_center_fixed_block": float(np.linalg.norm(c)),
        "color_center_fixed_one_particle_rank": rank_neutral,
        "fourth_moment_path_check_error": float(check_error),
        "second_moment_even_increment_trace": float(np.trace(second).real),
        "color_contracted_fourth_moment_odd_at_a_0_13": float(signed_fourth),
        "half_turn_with_electromagnetic_rotation_conjugates_sign": True,
        "color_fixed_fourth_return_is_odd_under_half_turn": True,
        "color_fixed_fourth_return_norm_per_amplitude": float(np.linalg.norm(mark)),
        "same_face_cauchy_readout_norm_per_amplitude": float(np.linalg.norm(on_cauchy)),
        "opposite_face_cauchy_readout_norm_per_amplitude": float(np.linalg.norm(opposite_faces)),
        "minimum_same_face_cauchy_multiplier": float(np.min(gram_multiplier)),
        "gauge_invariant_higgs_density_cos_7_and_16": density_marks,
        "fixed_cauchy_open_frobenius_at_a_0_13": float(np.linalg.norm(selected_open)),
        "fixed_cauchy_full_return_trace": float(np.trace(full_return).real),
        "fixed_cauchy_reprojected_return_trace": float(np.trace(selected_return).real),
        "fixed_cauchy_discarded_return_trace": float(return_loss),
        "physical_source_weight_supplied": False,
    }, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
