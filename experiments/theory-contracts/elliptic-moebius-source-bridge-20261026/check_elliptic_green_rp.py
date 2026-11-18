#!/usr/bin/env python3
"""Small independent checks for ELLIPTIC-GREEN-RP.txt."""

from __future__ import annotations

import hashlib
import json

import mpmath as mp
import sympy as sp


def main() -> None:
    E, t, u = sp.symbols("E t u", real=True)
    den = 1 + sp.exp(-E)
    s_plus = sp.exp(-t * E) / den
    s_minus = -sp.exp(-(t + 1) * E) / den

    jump = sp.simplify(s_plus.subs(t, 0) - s_minus.subs(t, 0)) == 1
    antiperiodic = sp.simplify(
        s_minus.subs(t, t - 1) + s_plus
    ) == 0

    a = t + u
    k_bar = sp.exp(a * E) / (1 + sp.exp(E))
    x_factor = (
        sp.exp(t * E) / sp.sqrt(1 + sp.exp(E))
        * sp.exp(u * E) / sp.sqrt(1 + sp.exp(E))
    )
    k_psi = sp.exp(-a * E) / (1 + sp.exp(-E))
    y_factor = (
        sp.exp(-t * E) / sp.sqrt(1 + sp.exp(-E))
        * sp.exp(-u * E) / sp.sqrt(1 + sp.exp(-E))
    )
    both_charge_factorizations = (
        sp.simplify(k_bar - x_factor) == 0
        and sp.simplify(k_psi - y_factor) == 0
    )

    # Exact finite-mode Gram controls.  Each mode contributes an outer
    # product, so arbitrary positive sums remain positive.
    time_points = [sp.Rational(1, 10), sp.Rational(1, 4), sp.Rational(2, 5)]
    mode_vectors = [
        sp.Matrix([sp.Rational(2) ** tt for tt in time_points]),
        sp.Matrix([sp.Rational(3) ** (-tt) for tt in time_points]),
        sp.Matrix([sp.Rational(5) ** tt for tt in time_points]),
    ]
    gram = sum((v * v.T for v in mode_vectors), sp.zeros(3))
    gram_positive = all(
        sp.simplify(gram[:k, :k].det()) > 0 for k in range(1, 4)
    )

    # Theta reverses a word, charge, time and winding.  Its real constant
    # phase is -1 per generator, so applying it twice gives +1.
    word = [
        ("psi", sp.Rational(1, 7), 2),
        ("bar", sp.Rational(1, 5), -1),
        ("psi", sp.Rational(1, 3), 0),
    ]

    def theta_word(items):
        swapped = {"psi": "bar", "bar": "psi"}
        return (-1) ** len(items), [
            (swapped[charge], -time, -winding)
            for charge, time, winding in reversed(items)
        ]

    c1, w1 = theta_word(word)
    c2, w2 = theta_word(w1)
    theta_square = c1 * c2 == 1 and w2 == word

    # Mode-label versions of B u=u(B+2pi) and Q transport.
    winding_domain = True
    q_cut_transport = True
    for n in range(-8, 9):
        lam = sp.Rational(2) * sp.pi * (sp.Rational(n) + sp.Rational(1, 2))
        lam_wound = sp.Rational(2) * sp.pi * (sp.Rational(n + 1) + sp.Rational(1, 2))
        winding_domain &= sp.simplify(lam_wound - lam - 2 * sp.pi) == 0

        q_label = -n - 1
        q_lam = sp.Rational(2) * sp.pi * (sp.Rational(q_label) + sp.Rational(1, 2))
        q_cut_transport &= sp.simplify(q_lam + lam) == 0
        q_cut_transport &= (-n - 2) == q_label - 1

    # Numerical comparison only illustrates the exact theta identity proved
    # by pole, multipliers and H0(E,alpha)=0 in the text.
    mp.mp.dps = 60
    nome = mp.e ** (-mp.pi)
    theta1_prime = mp.diff(lambda x: mp.jtheta(1, x, nome), 0)

    def spectral_kernel(tt, ss, cutoff=180):
        total = mp.mpc(0)
        for n in range(-cutoff, cutoff + 1):
            lam = 2 * mp.pi * (n + mp.mpf("0.5"))
            total += (
                mp.e ** (2j * mp.pi * (n + mp.mpf("0.5")) * ss)
                * mp.e ** (-tt * lam)
                / (1 + mp.e ** (-lam))
            )
        return total

    def theta_kernel(tt, ss):
        w = ss + 1j * tt
        return (
            1j / 2
            * theta1_prime
            * mp.jtheta(3, mp.pi * w, nome)
            / (mp.jtheta(3, 0, nome) * mp.jtheta(1, mp.pi * w, nome))
        )

    theta_samples = [(mp.mpf("0.2"), mp.mpf("0.13")),
                     (mp.mpf("0.37"), mp.mpf("-0.21"))]
    theta_match = all(
        abs(spectral_kernel(tt, ss) - theta_kernel(tt, ss)) < mp.mpf("1e-45")
        for tt, ss in theta_samples
    )

    checks = {
        "green_jump_is_identity": bool(jump),
        "temporal_NS_antiperiodicity": bool(antiperiodic),
        "both_charge_Gram_factorizations": bool(both_charge_factorizations),
        "finite_multimode_Gram_is_positive": bool(gram_positive),
        "field_reflection_squares_to_one": bool(theta_square),
        "winding_preserves_NS_labels_and_shifts_B": bool(winding_domain),
        "Q_maps_cut_spectrum_and_winding_to_inverses": bool(q_cut_transport),
        "spectral_kernel_matches_theta3_Szego_samples": bool(theta_match),
    }
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])

    payload = json.dumps(checks, sort_keys=True)
    print(f"PASS {len(checks)}/{len(checks)}")
    print(f"RESULT_SHA256 {hashlib.sha256(payload.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
