#!/usr/bin/env python3
"""Operational wave-pair probe in the finite Z174 standard-form model.

No Hamiltonian is added.  The state, beta=4*pi weight, standard-form left and
right actions, modular clock, and the alternative charge weight are exactly
the finite objects used in z174.  The four source words are 00, 10, 01, 11.

The test separates three notions which are easy to conflate:

1. original--Tomita-partner matching already present in the standard vector;
2. a connected two-source response of that state;
3. selection of a diagonal Z4 sector over two complete seams by projection.
"""

from __future__ import annotations

import math

import numpy as np


TOL = 2.0e-12
PI = math.pi


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def expectation(vector: np.ndarray, operator: np.ndarray) -> complex:
    return complex(np.vdot(vector, operator @ vector))


def source_unitary(operator: np.ndarray, epsilon: float) -> np.ndarray:
    return math.cos(epsilon) * np.eye(operator.shape[0]) - 1j * math.sin(
        epsilon
    ) * operator


def hermitian_unitary(operator: np.ndarray, epsilon: float) -> np.ndarray:
    """exp(-i epsilon operator) for a finite Hermitian generator."""
    values, vectors = np.linalg.eigh(operator)
    return (vectors * np.exp(-1j * epsilon * values)) @ vectors.conj().T


def one_mode_standard_form() -> dict[str, object]:
    """One original mode and its Tomita partner, basis 00,10,01,11."""
    labels = ("00", "10", "01", "11")
    bits = ((0, 0), (1, 0), (0, 1), (1, 1))
    index = {b: j for j, b in enumerate(bits)}

    def pauli(axis: str, side: int) -> np.ndarray:
        out = np.zeros((4, 4), dtype=complex)
        for column, bit in enumerate(bits):
            target = list(bit)
            target[side] ^= 1
            if axis == "X":
                phase = 1.0
            elif axis == "Y":
                phase = 1j if bit[side] == 0 else -1j
            else:
                raise ValueError(axis)
            out[index[tuple(target)], column] = phase
        return out

    xa, ya = pauli("X", 0), pauli("Y", 0)
    xb, yb = pauli("X", 1), pauli("Y", 1)
    require(np.max(np.abs(xa @ xb - xb @ xa)) < TOL, "left/right actions do not commute")

    # z174: beta=4*pi and a quarter-energy gap.  Hence p1/p0=e^-pi.
    q = math.exp(-PI / 2.0)
    p0, p1 = 1.0 / (1.0 + q * q), q * q / (1.0 + q * q)
    omega = np.array([math.sqrt(p0), 0.0, 0.0, math.sqrt(p1)], dtype=complex)
    support_weights = np.abs(omega) ** 2
    require(abs(np.vdot(omega, omega) - 1.0) < TOL, "standard vector not normalized")
    require(
        np.max(np.abs(support_weights - np.array([p0, 0.0, 0.0, p1]))) < TOL,
        "wrong 00/10/01/11 support",
    )

    # Delta vec(Y)=vec(rho Y rho^-1), and z174's pair clock is Delta^(-i/2).
    delta = np.diag([1.0, math.exp(-PI), math.exp(PI), 1.0])
    clock = np.diag(np.exp(-0.5j * np.log(np.diag(delta))))
    expected_clock = np.diag([1.0, 1j, -1j, 1.0])
    require(np.max(np.abs(clock - expected_clock)) < TOL, "wrong modular half-time clock")
    require(np.max(np.abs(clock @ omega - omega)) < TOL, "clock moves Omega")

    # Four source histories.  Z is the vacuum persistence amplitude, and its
    # connected mixed part is the finite operational version of a pair response.
    epsilon = 0.173
    ua, ub = source_unitary(xa, epsilon), source_unitary(xb, epsilon)
    z00 = 1.0 + 0.0j
    z10 = expectation(omega, ua)
    z01 = expectation(omega, ub)
    z11 = expectation(omega, ua @ ub)
    z11_swapped = expectation(omega, ub @ ua)
    coherence = expectation(omega, xa @ xb).real
    connected = z11 - z10 * z01
    predicted = -(math.sin(epsilon) ** 2) * coherence
    require(abs(z11 - z11_swapped) < TOL, "commuting standard actions became ordered")
    require(abs(connected - predicted) < TOL, "mixed source response formula failed")
    require(abs(coherence - 1.0 / math.cosh(PI / 2.0)) < TOL, "wrong native pair weight")

    # The same diagonal KMS weights without the purification coherence do not
    # produce the pair response.  Conditioning on even support does not restore it.
    rho_product = np.diag([p0 * p0, p1 * p0, p0 * p1, p1 * p1])
    even = np.diag([1.0, 0.0, 0.0, 1.0])
    product_coherence = np.trace(rho_product @ xa @ xb)
    even_probability = float(np.trace(rho_product @ even).real)
    rho_conditioned = even @ rho_product @ even / even_probability
    conditioned_coherence = np.trace(rho_conditioned @ xa @ xb)
    require(abs(product_coherence) < TOL, "diagonal product state acquired coherence")
    require(abs(conditioned_coherence) < TOL, "support conditioning created coherence")
    require(np.max(np.abs(even @ omega - omega)) < TOL, "Omega is not already even")

    # Clock-time control.  Evolving one insertion relatively rotates its X
    # quadrature into Y.  The invariant two-quadrature norm stays fixed.
    clock_rows = []
    for modular_time in (0.0, -0.25, -0.5, -1.0):
        flow = np.diag(np.exp(1j * modular_time * np.log(np.diag(delta))))
        xb_t = flow @ xb @ flow.conj().T
        yb_t = flow @ yb @ flow.conj().T
        cx = expectation(omega, xa @ xb_t).real
        cy = expectation(omega, xa @ yb_t).real
        norm = math.hypot(cx, cy)
        require(abs(norm - coherence) < TOL, "clock changed quadrature-normalized weight")
        clock_rows.append((modular_time, cx, cy, norm))
    require(abs(clock_rows[2][1]) < TOL, "X response did not vanish at t=-1/2")
    require(abs(abs(clock_rows[2][2]) - coherence) < TOL, "clock lost the conjugate response")

    # Raw CAR fields are different from commuting standard-form left/right
    # actions.  Their creation order changes only the fermionic sign.
    a = np.zeros((4, 4), dtype=complex)
    b = np.zeros((4, 4), dtype=complex)
    a[index[(0, 0)], index[(1, 0)]] = 1.0
    a[index[(0, 1)], index[(1, 1)]] = 1.0
    b[index[(0, 0)], index[(0, 1)]] = 1.0
    b[index[(1, 0)], index[(1, 1)]] = -1.0
    ad, bd = a.conj().T, b.conj().T
    ordered = ad @ bd @ omega
    reversed_order = bd @ ad @ omega
    require(np.max(np.abs(ordered + reversed_order)) < TOL, "CAR order lacks its sign")
    require(abs(np.vdot(ordered, ordered) - p0) < TOL, "pair-creation weight changed")
    require(
        abs(np.vdot(reversed_order, reversed_order) - p0) < TOL,
        "reversed pair-creation probability changed",
    )

    return {
        "labels": labels,
        "q": q,
        "p0": p0,
        "p1": p1,
        "support_weights": support_weights,
        "coherence": coherence,
        "source_amplitudes": {
            "00": z00,
            "10": z10,
            "01": z01,
            "11": z11,
            "11_swapped": z11_swapped,
            "connected": connected,
        },
        "product_coherence": product_coherence,
        "conditioned_coherence": conditioned_coherence,
        "even_probability_product": even_probability,
        "clock_rows": clock_rows,
        "car_pair_probability": float(np.vdot(ordered, ordered).real),
    }


def two_seam_glue_control(p0: float, p1: float) -> dict[str, object]:
    """Two complete seam pairs, reduced to their 00/11 standard branches.

    This is only the occupancy branch of the finite model.  It locates which
    selection is due to the state and which is due to a diagonal Z4 projector;
    it is not the full modular-orbit sum of z174 B2--B6.
    """
    labels = ("00", "10", "01", "11")
    state = np.array(
        [p0, math.sqrt(p0 * p1), math.sqrt(p0 * p1), p1], dtype=complex
    )
    weights = np.abs(state) ** 2
    require(abs(weights.sum() - 1.0) < TOL, "two-seam state not normalized")
    require(np.all(weights > 0.0), "state itself removed a two-seam word")

    # Each occupied standard branch has clock grade one.  The diagonal glue
    # has grades 0,1,1,2 on 00,10,01,11.  X=-1 is therefore a projector choice.
    x_diag = np.diag([1.0, 1j, 1j, -1.0])
    p_x_minus = sum(np.linalg.matrix_power(x_diag, k) * ((-1.0) ** k) for k in range(4)) / 4.0
    g = x_diag @ x_diag
    p_g_plus = (np.eye(4) + g) / 2.0
    require(np.max(np.abs(p_x_minus - np.diag([0.0, 0.0, 0.0, 1.0]))) < TOL, "bad X=-1 projector")
    require(np.max(np.abs(p_x_minus @ p_g_plus - p_x_minus)) < TOL, "X=-1 lacks G=+1")
    selected_weight = float(np.vdot(state, p_x_minus @ state).real)
    require(abs(selected_weight - p1 * p1) < TOL, "wrong projector cost")

    swap = np.array(
        [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]],
        dtype=complex,
    )
    require(np.max(np.abs(swap @ state - state)) < TOL, "seam order changes source state")
    require(np.max(np.abs(swap @ p_x_minus @ swap - p_x_minus)) < TOL, "projector depends on seam order")

    return {
        "labels": labels,
        "weights": weights,
        "x_eigenvalues": np.diag(x_diag),
        "selected_weight": selected_weight,
        "normalization_after_conditioning": 1.0 / selected_weight,
    }


def even_two_mode_standard_form() -> dict[str, object]:
    """Strictly even native preparations on the K=2 Z174 standard form.

    H conserves charge and moves one fermion between the two modes.  D is an
    even but charge-changing pair field used as an honest control.  J is the
    finite standard-form modular conjugation F K.
    """
    bits = ((0, 0), (1, 0), (0, 1), (1, 1))
    index = {b: j for j, b in enumerate(bits)}
    annihilators = []
    for mode in range(2):
        c = np.zeros((4, 4), dtype=complex)
        for column, bit in enumerate(bits):
            if bit[mode] == 0:
                continue
            target = list(bit)
            target[mode] = 0
            sign = (-1) ** sum(bit[:mode])
            c[index[tuple(target)], column] = sign
        annihilators.append(c)
    c1, c2 = annihilators
    c1d, c2d = c1.conj().T, c2.conj().T
    charge = np.diag([sum(b) for b in bits]).astype(complex)
    parity = np.diag([(-1) ** sum(b) for b in bits]).astype(complex)
    hop = c1d @ c2 + c2d @ c1
    pair = c1d @ c2d + c2 @ c1
    pair_y = -1j * (c1d @ c2d - c2 @ c1)
    require(np.max(np.abs(hop - hop.conj().T)) < TOL, "hopping generator not Hermitian")
    require(np.max(np.abs(pair - pair.conj().T)) < TOL, "pair generator not Hermitian")
    require(np.max(np.abs(hop @ parity - parity @ hop)) < TOL, "hopping generator is odd")
    require(np.max(np.abs(hop @ charge - charge @ hop)) < TOL, "hopping generator changes charge")
    require(np.max(np.abs(pair @ parity - parity @ pair)) < TOL, "pair control is odd")
    require(np.max(np.abs(pair @ charge - charge @ pair)) > 1.0, "pair control conserved charge")

    r = math.exp(-PI)
    probabilities = np.array([1.0, r, r, r * r]) / (1.0 + r) ** 2
    omega = np.diag(np.sqrt(probabilities)).reshape(-1)
    identity = np.eye(4)
    flip = np.zeros((16, 16), dtype=complex)
    for i in range(4):
        for j in range(4):
            flip[j * 4 + i, i * 4 + j] = 1.0

    def lift(generator: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        left = np.kron(generator, identity)
        # For J=F K, J L J is represented by F conjugate(L) F.
        right = flip @ left.conj() @ flip
        expected_right = np.kron(identity, generator.conj())
        require(np.max(np.abs(right - expected_right)) < TOL, "wrong J-conjugate generator")
        require(np.max(np.abs(left @ right - right @ left)) < TOL, "standard sides do not commute")
        return left, right

    epsilon = 0.173

    def circuits(generator: np.ndarray) -> dict[str, object]:
        left, right = lift(generator)
        ul = hermitian_unitary(left, epsilon)
        ur_same = hermitian_unitary(right, epsilon)
        # J is antiunitary, hence J exp(-i eps H) J = exp(+i eps JHJ).
        ur_j = flip @ ul.conj() @ flip
        require(
            np.max(np.abs(ur_j - hermitian_unitary(right, -epsilon))) < TOL,
            "antiunitary J did not reverse the source parameter",
        )

        def four(ur: np.ndarray) -> dict[str, object]:
            states = {
                "00": omega,
                "10": ul @ omega,
                "01": ur @ omega,
                "11": ul @ ur @ omega,
                "11_swapped": ur @ ul @ omega,
            }
            norms = {key: float(np.vdot(value, value).real) for key, value in states.items()}
            amplitudes = {key: complex(np.vdot(omega, value)) for key, value in states.items()}
            require(max(abs(value - 1.0) for value in norms.values()) < TOL, "circuit lost normalization")
            require(
                np.max(np.abs(states["11"] - states["11_swapped"])) < TOL,
                "even standard-form circuit became order dependent",
            )
            amplitudes["connected"] = amplitudes["11"] - amplitudes["10"] * amplitudes["01"]
            return {"amplitudes": amplitudes, "norms": norms}

        variance = expectation(omega, left @ left).real
        covariance = expectation(omega, left @ right).real
        normalized = covariance / variance
        require(variance > TOL, "zero preparation variance")
        return {
            "variance": variance,
            "covariance": covariance,
            "normalized_covariance": normalized,
            "same_orientation": four(ur_same),
            "J_conjugate": four(ur_j),
            "left": left,
            "right": right,
        }

    hop_result = circuits(hop)
    pair_result = circuits(pair)
    require(abs(hop_result["normalized_covariance"] - 1.0) < TOL, "J hopping pair not perfectly matched")
    require(
        np.max(
            np.abs(
                hop_result["J_conjugate"]["amplitudes"]["11"]
                - 1.0
            )
        )
        < TOL,
        "J-conjugate hopping circuit did not return Omega",
    )
    require(
        abs(pair_result["normalized_covariance"] - 1.0 / math.cosh(PI)) < TOL,
        "wrong charge-changing even-pair normalization",
    )

    # Modular time is native state time, not an added Hamiltonian.  H commutes
    # with rho=f(Q) and is clock neutral.  D rotates into its conjugate D_y.
    delta = np.diag(
        [probabilities[i] / probabilities[j] for i in range(4) for j in range(4)]
    )
    hop_left, hop_right = lift(hop)
    pair_left, pair_right = lift(pair)
    _, pair_y_right = lift(pair_y)
    clock = []
    for modular_time in (0.0, -0.25, -0.5):
        flow = np.diag(np.exp(1j * modular_time * np.log(np.diag(delta))))
        hop_t = flow @ hop_right @ flow.conj().T
        pair_t = flow @ pair_right @ flow.conj().T
        pair_y_t = flow @ pair_y_right @ flow.conj().T
        h_value = expectation(omega, hop_left @ hop_t).real
        d_x = expectation(omega, pair_left @ pair_t).real
        d_y = expectation(omega, pair_left @ pair_y_t).real
        require(abs(h_value - hop_result["covariance"]) < TOL, "charge-neutral preparation moved in modular time")
        require(abs(math.hypot(d_x, d_y) - pair_result["covariance"]) < TOL, "even pair lost clock norm")
        clock.append((modular_time, h_value, d_x, d_y))

    return {
        "probabilities": probabilities,
        "hop": hop_result,
        "pair_control": pair_result,
        "clock": clock,
    }


def side_mode_matrix(mode: int) -> np.ndarray:
    """Pauli X on one of two occupation bits, order 00,10,01,11."""
    bits = ((0, 0), (1, 0), (0, 1), (1, 1))
    index = {b: j for j, b in enumerate(bits)}
    out = np.zeros((4, 4), dtype=complex)
    for column, bit in enumerate(bits):
        target = list(bit)
        target[mode] ^= 1
        out[index[tuple(target)], column] = 1.0
    return out


def z174_weight_control() -> dict[str, object]:
    """Compare both state weights explicitly used by z174 for K=2 modes."""
    charges = np.array([0, 1, 1, 2], dtype=float)
    weights = {
        "K=pi Q": np.exp(-PI * charges),
        "K=2pi Q^2-pi Q": np.exp(-(2.0 * PI * charges**2 - PI * charges)),
    }
    table = {}
    for name, raw in weights.items():
        p = raw / raw.sum()
        root = np.sqrt(p)
        omega = np.diag(root).reshape(-1)
        corr = np.zeros((2, 2), dtype=float)
        identity = np.eye(4)
        for j in range(2):
            left = np.kron(side_mode_matrix(j), identity)
            for k in range(2):
                right = np.kron(identity, side_mode_matrix(k).T)
                corr[j, k] = expectation(omega, left @ right).real
        n0 = np.array([0, 1, 0, 1], dtype=float)
        n1 = np.array([0, 0, 1, 1], dtype=float)
        cov_n = float(np.dot(p, n0 * n1) - np.dot(p, n0) * np.dot(p, n1))
        hop = np.zeros((4, 4), dtype=complex)
        hop[1, 2] = hop[2, 1] = 1.0
        hop_left = np.kron(hop, identity)
        hop_right = np.kron(identity, hop.T)
        hop_variance = expectation(omega, hop_left @ hop_left).real
        hop_covariance = expectation(omega, hop_left @ hop_right).real
        hop_normalized = hop_covariance / hop_variance
        require(np.max(np.abs(corr - np.diag(np.diag(corr)))) < TOL, "state stopped matching modes")
        require(abs(hop_normalized - 1.0) < TOL, "even J-pair lost normalized matching")
        table[name] = {
            "probabilities": p,
            "pair_matrix": corr,
            "number_covariance": cov_n,
            "even_hop_raw": hop_covariance,
            "even_hop_normalized": hop_normalized,
        }

    linear = table["K=pi Q"]
    quadratic = table["K=2pi Q^2-pi Q"]
    expected = 1.0 / math.cosh(PI / 2.0)
    require(np.max(np.abs(np.diag(linear["pair_matrix"]) - expected)) < TOL, "linear weight did not factorize")
    require(abs(linear["number_covariance"]) < TOL, "linear weight correlated spectators")
    require(abs(quadratic["number_covariance"]) > 1.0e-4, "alternative weight failed its control")
    require(
        np.max(np.abs(linear["pair_matrix"] - quadratic["pair_matrix"])) > 1.0e-4,
        "alternative native weight left normalization unchanged",
    )
    require(
        abs(linear["even_hop_raw"] - quadratic["even_hop_raw"]) > 1.0e-4,
        "alternative native weight left even raw response unchanged",
    )
    return table


def main() -> None:
    one = one_mode_standard_form()
    two = two_seam_glue_control(one["p0"], one["p1"])
    even = even_two_mode_standard_form()
    alternative = z174_weight_control()

    print("STANDARDFORM PAIR TEST: PASS")
    print(f"basis: {one['labels']}")
    print(
        "native one-seam weights: "
        + ", ".join(
            f"{label}={weight:.12g}"
            for label, weight in zip(one["labels"], one["support_weights"])
        )
    )
    print(f"native connected pair coefficient: sech(pi/2)={one['coherence']:.12g}")
    print("four source amplitudes:")
    for key, value in one["source_amplitudes"].items():
        print(f"  {key}: {value.real:.12g}{value.imag:+.3g}i")
    print(
        "controls: diagonal product coherence="
        f"{abs(one['product_coherence']):.1e}, even-conditioned coherence="
        f"{abs(one['conditioned_coherence']):.1e}"
    )
    print("relative modular time (t, X response, Y response, quadrature norm):")
    for row in one["clock_rows"]:
        print(f"  {row[0]:5.2f}: {row[1]: .12g}, {row[2]: .12g}, {row[3]:.12g}")
    print("two-seam standard-state weights before diagonal Z4 projection:")
    print(
        "  "
        + ", ".join(
            f"{label}={weight:.12g}" for label, weight in zip(two["labels"], two["weights"])
        )
    )
    print(
        "  X=-1/G=+1 conditioned branch weight="
        f"{two['selected_weight']:.12g}, renormalization={two['normalization_after_conditioning']:.12g}"
    )
    print("strictly even K=2 side preparation H=c1^+c2+c2^+c1:")
    print(
        f"  raw variance=covariance={even['hop']['variance']:.12g}, "
        f"normalized J-pair coefficient={even['hop']['normalized_covariance']:.12g}"
    )
    for orientation in ("same_orientation", "J_conjugate"):
        amplitudes = even["hop"][orientation]["amplitudes"]
        print(f"  {orientation} four circuits:")
        print(
            "    "
            + ", ".join(
                f"{key}={value.real:.12g}{value.imag:+.2g}i"
                for key, value in amplitudes.items()
            )
        )
    print(
        "  even charge-changing control normalized coefficient="
        f"sech(pi)={even['pair_control']['normalized_covariance']:.12g}"
    )
    print("  modular time (t, charge-neutral H, pair X, pair Y):")
    for row in even["clock"]:
        print(f"    {row[0]:5.2f}: {row[1]: .12g}, {row[2]: .12g}, {row[3]: .12g}")
    for name, row in alternative.items():
        diag = np.diag(row["pair_matrix"])
        print(
            f"native weight control {name}: matched-pair diagonal=({diag[0]:.12g}, {diag[1]:.12g}), "
            f"offdiag={row['pair_matrix'][0,1]:.1e}, even H raw/norm={row['even_hop_raw']:.12g}/"
            f"{row['even_hop_normalized']:.12g}, number covariance={row['number_covariance']:.12g}"
        )
    print(
        "VERDICT: the physically even charge-preserving preparation has an exact J-matched response but no order or "
        "clock signal; the two-seam four-word table is only an occupancy control, while the full Z4 selection belongs "
        "to z174's modular orbit sum and conditioning."
    )


if __name__ == "__main__":
    main()
