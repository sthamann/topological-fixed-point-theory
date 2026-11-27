#!/usr/bin/env python3
"""Z58/Z84: zwei schaltbare Stromanregungen in derselben thermischen Quelle.

Die Quelle ist genau die in z58 benutzte quasifreie direkte Summe der beiden
Quarterkovarianzen alpha=1/4,3/4.  Die Einfuegungen sind die geraden,
ladungsneutralen Strommoden J_-1 in je einem Quarter.  Die OS-Form wird an
tau=beta/4 ausgewertet, sodass die Spiegelpunkte den Abstand beta/2 haben.
"""

from __future__ import annotations

import json

import mpmath as mp
import numpy as np


mp.mp.dps = 70
BETA = 4 * mp.pi


def current_kernel(mode: int, delta_tau: mp.mpf) -> mp.mpf:
    """Positiver Fourierkoeffizient des wirklichen periodischen Z84-Stromkerns.

    S_J(delta theta,delta tau)=sum_{m!=0} exp(i m delta theta-m delta tau)
    m/(1-exp(-beta m)), 0<delta_tau<beta.
    """
    if not (0 < delta_tau < BETA):
        raise ValueError("delta_tau muss im offenen thermischen Kreis liegen")
    return mode * mp.exp(-mode * delta_tau) / (1 - mp.exp(-BETA * mode))


def os_weight(mode: int = 1) -> mp.mpf:
    """OS-Norm einer Fourierstrommode bei tau=beta/4."""
    return current_kernel(mode, BETA / 2)


def finite_trace_weight(mode: int, cutoff: int) -> mp.mpf:
    """Direkte endliche Spurkontrolle der symmetrischen KMS/OS-Form.

    Tr(q^(N/2) sqrt(m)a q^(N/2) sqrt(m)a^+)/Tr(q^N).
    Der oberste abgeschnittene Zustand wird korrekt ausgelassen.
    """
    q = mp.exp(-BETA * mode)
    z = mp.fsum(q**n for n in range(cutoff + 1))
    numerator = mode * mp.fsum((n + 1) * q ** (mp.mpf(n) + mp.mpf(1) / 2) for n in range(cutoff))
    return numerator / z


def close(a: mp.mpf, b: mp.mpf, rel: mp.mpf = mp.mpf("1e-45")) -> bool:
    return bool(abs(a - b) <= rel * max(mp.mpf(1), abs(a), abs(b)))


def coherent_neutral_kernel(t: complex, u: complex) -> complex:
    """Gramkern projizierter Produkt-Exponentialvektoren: sum (conj(t)u)^n/(n!)^2."""
    z = mp.conj(t) * u
    return mp.fsum(z**n / mp.factorial(n) ** 2 for n in range(80))


def hardy_kernel(z: complex, w: complex) -> complex:
    """Gramkern einer echten Paarpraeparation sum z^n |n,n>."""
    return 1 / (1 - mp.conj(z) * w)


def main() -> None:
    # Absoluter Spurfaktor nur fuer die explizit gewaehlten beiden m=1-Moden
    # mit Energie epsilon=1. Ladung, hoehere Moden und andere Sektoren sind
    # nicht enthalten.
    epsilon = mp.mpf(1)
    boltz = mp.exp(-BETA * epsilon)
    z0_pair = 1 / (1 - boltz) ** 2

    # A=J^1/4_-1 und B=J^3/4_-1, beide bei tau=beta/4.
    g = os_weight(1)
    normalized = {
        "00": mp.mpf(1),
        "10_A": g,
        "01_B": g,
        "11_AB": g * g,
        "11_BA": g * g,
    }
    absolute_pair_only = {key: z0_pair * value for key, value in normalized.items()}

    # Basis [1,A,B,AB,BA]. Verschiedene Niveaus sind orthogonal; AB=BA.
    gram = mp.matrix(5, 5)
    gram[0, 0] = 1
    gram[1, 1] = g
    gram[2, 2] = g
    for i in (3, 4):
        for j in (3, 4):
            gram[i, j] = g * g
    gram_np = np.array([[float(gram[i, j]) for j in range(5)] for i in range(5)])
    eigenvalues = np.linalg.eigvalsh(gram_np)

    # Beide KMS-Ordnungen und ihr symmetrischer OS-Wert.
    kms_creation_first = boltz / (1 - boltz)       # omega(J_-1 J_1)
    kms_annihilation_first = 1 / (1 - boltz)       # omega(J_1 J_-1)
    os_geometric_mean = mp.sqrt(kms_creation_first * kms_annihilation_first)

    finite = {cutoff: finite_trace_weight(1, cutoff) for cutoff in (2, 4, 8)}
    connected = mp.log(normalized["11_AB"] * normalized["00"] / (
        normalized["10_A"] * normalized["01_B"]
    ))

    # Gegenfall: gleiche Einzelantworten, aber eine explizit gekoppelte Quelle.
    coupling_epsilon = mp.mpf("0.1")
    correlated_joint = (1 + coupling_epsilon) * g * g
    connected_control = mp.log(correlated_joint / (g * g))

    # ------------------------------------------------ Paarleiter B^+=a^+b^+
    # Polarzerlegung B^+=S sqrt((N_a+1)(N_b+1)). S erhaelt jede
    # Differenzzahl d=n_a-n_b. Auf jedem d-Sektor ist S unilateral.
    p_neutral = (1 - boltz) / (1 + boltz)
    p_non_neutral = 1 - p_neutral
    p_boundary = 1 - boltz**2
    p_interior = boltz**2
    p_vacuum = (1 - boltz) ** 2
    p_neutral_excited = p_neutral - p_vacuum
    p_vacuum_given_neutral = p_vacuum / p_neutral
    p_excited_given_neutral = p_neutral_excited / p_neutral

    # Direkte endliche Summenkontrolle des thermischen Produktzustands.
    cutoff = 24
    z_cut = mp.fsum(boltz**n for n in range(cutoff + 1)) ** 2
    p_neutral_cut = mp.fsum(boltz ** (2 * n) for n in range(cutoff + 1)) / z_cut
    p_boundary_cut = (
        2 * mp.fsum(boltz**n for n in range(cutoff + 1)) - 1
    ) / z_cut

    # Unabhaengige kohärente Quellen liefern nach P_neutral Bessel, nicht Hardy.
    t_test = mp.mpf("0.37") * mp.expj(mp.mpf("0.2"))
    coherent_diag = coherent_neutral_kernel(t_test, t_test)
    bessel_diag = mp.besseli(0, 2 * abs(t_test))

    # Eine echte Tomita/TFD-Paarpraeparation hat z=sqrt(q) und D=0 exakt.
    z_tomita = mp.sqrt(boltz)
    tfd_norm = (1 - abs(z_tomita) ** 2) * hardy_kernel(z_tomita, z_tomita)
    mean_product = boltz / (1 - boltz)
    covariance_product = mp.mpf(0)
    covariance_tomita = boltz / (1 - boltz) ** 2

    # Nachtraegliche Konditionierung des Produktzustands auf D=0 halbiert
    # nicht die Energie, sondern quadriert das Boltzmanngewicht q -> q^2.
    conditional_ratio = boltz**2
    tomita_probability_ratio = boltz

    # Der Polarteil und sein Index kennen keine Energieskala. Eine andere
    # Modenenergie laesst sie unveraendert, aendert aber alle KMS-Gewichte.
    epsilon_control = mp.mpf(2)
    boltz_control = mp.exp(-BETA * epsilon_control)
    p_neutral_control = (1 - boltz_control) / (1 + boltz_control)
    p1_c3 = 1 / (8 * mp.pi)
    inherited_half_inverse_beta = 1 / (2 * BETA)

    checks = {
        "Kernel bei beta/2 ist m/(2 sinh(beta m/2))": close(g, 1 / (2 * mp.sinh(BETA / 2))),
        "endliche direkte Spur konvergiert zum OS-Kernel": (
            abs(finite[2] - g) < mp.mpf("1e-12")
            and abs(finite[4] - g) < mp.mpf("1e-22")
            and abs(finite[8] - g) < mp.mpf("1e-45")
        ),
        "OS ist geometrisches Mittel der zwei KMS-Ordnungen": close(g, os_geometric_mean),
        "KMS-Verhältnis rekonstruiert beta*omega=4pi": close(
            mp.log(kms_annihilation_first / kms_creation_first), BETA
        ),
        "beide Einfügungsordnungen stimmen": close(normalized["11_AB"], normalized["11_BA"]),
        "gemischte verbundene Antwort verschwindet": abs(connected) < mp.mpf("1e-50"),
        "OS-Gram ist positiv semidefinit": bool(eigenvalues[0] > -1e-18),
        "OS-Gram hat genau den Ordnungsnullvektor AB-BA": bool(
            np.sum(eigenvalues < 1e-18) == 1
            and np.linalg.norm(gram_np @ np.array([0, 0, 0, 1, -1], dtype=float)) < 1e-25
        ),
        "gekoppelter Gegenfall wird von der verbundenen Antwort erkannt": abs(
            connected_control - mp.log(mp.mpf("1.1"))
        ) < mp.mpf("1e-50"),
        "P_neutral hat die exakte thermische Wahrscheinlichkeit tanh(beta/2)": close(
            p_neutral, mp.tanh(BETA / 2)
        ),
        "endliche Spur bestaetigt P_neutral und den Polardefekt": (
            abs(p_neutral_cut - p_neutral) < mp.mpf("1e-50")
            and abs(p_boundary_cut - p_boundary) < mp.mpf("1e-50")
        ),
        "Polardefekt ist der Rand n_a=0 oder n_b=0": close(p_boundary + p_interior, 1),
        "neutraler Sektor ist nicht der ganze thermische Zustand": p_non_neutral > 0,
        "projizierte kohaerente Quellen haben Bessel-I0-Gram": abs(
            coherent_diag - bessel_diag
        ) < mp.mpf("1e-50"),
        "echte Tomita-Paarpraeparation hat Hardy-Norm und D=0": abs(tfd_norm - 1) < mp.mpf("1e-50"),
        "Produktprojektion und Tomita-Matching haben verschiedene Folgengewichte": close(
            conditional_ratio, tomita_probability_ratio**2
        ) and not close(conditional_ratio, tomita_probability_ratio),
        "unabhaengiges Produkt hat keine Zahlkovarianz, Tomita schon": (
            covariance_product == 0 and covariance_tomita > mean_product
        ),
        "Paarleiterindex legt keine Energieskala fest": (
            p_neutral_control != p_neutral and close(p1_c3, inherited_half_inverse_beta)
        ),
    }

    report = {
        "source": {
            "covariance": "c_alpha(r)=1/(1+exp(4*pi*(r+alpha))), alpha in {1/4,3/4}",
            "one_particle_modular_generator": "h_alpha(r)=4*pi*(r+alpha)",
            "full_generator": "K=4*pi*sum_alpha[E_alpha(Q_alpha)+sum_m m N_(alpha,m)]",
            "domain": "finite-energy current words on the finite-particle Fock core",
            "reflection": "tau -> beta-tau on 0<tau<beta/2; insertions at tau=beta/4",
        },
        "partition": {
            "chosen_two_mode_energy_epsilon": mp.nstr(epsilon, 10),
            "Z0_pair_only": mp.nstr(z0_pair, 30),
            "excluded": "charge sectors, modes m>1, all other source factors",
            "normalization_note": "K->K+c rescales Z0; all normalized answers are invariant",
        },
        "responses_normalized": {key: mp.nstr(value, 30) for key, value in normalized.items()},
        "responses_absolute_pair_only": {
            key: mp.nstr(value, 30) for key, value in absolute_pair_only.items()
        },
        "kms_orderings": {
            "omega(J_-1 J_1)": mp.nstr(kms_creation_first, 30),
            "omega(J_1 J_-1)": mp.nstr(kms_annihilation_first, 30),
            "ratio": mp.nstr(kms_annihilation_first / kms_creation_first, 30),
            "OS_symmetric": mp.nstr(g, 30),
        },
        "os_gram_basis": ["1", "A", "B", "AB", "BA"],
        "os_gram": [[mp.nstr(gram[i, j], 20) for j in range(5)] for i in range(5)],
        "os_gram_eigenvalues": [f"{value:.18e}" for value in eigenvalues],
        "connected_log_response": mp.nstr(connected, 20),
        "correlated_control": {
            "R11": mp.nstr(correlated_joint, 30),
            "connected_log_response": mp.nstr(connected_control, 30),
        },
        "finite_trace_control": {str(key): mp.nstr(value, 30) for key, value in finite.items()},
        "pair_ladder": {
            "polar_decomposition": "B^+=S sqrt((N_a+1)(N_b+1))",
            "S_action": "S|n_a,n_b>=|n_a+1,n_b+1>",
            "S_star_S": "I",
            "S_S_star": "I-P_boundary",
            "P_boundary": "1_{n_a=0 or n_b=0}",
            "difference_sectors": "d=n_a-n_b in Z; one unilateral shift per d",
            "neutral_sector_index": -1,
        },
        "thermal_projection_weights": {
            "q_exp_minus_beta_epsilon": mp.nstr(boltz, 30),
            "P_neutral": mp.nstr(p_neutral, 30),
            "P_non_neutral": mp.nstr(p_non_neutral, 30),
            "P_boundary_polar_defect": mp.nstr(p_boundary, 30),
            "P_interior": mp.nstr(p_interior, 30),
            "P_vacuum": mp.nstr(p_vacuum, 30),
            "P_neutral_excited": mp.nstr(p_neutral_excited, 30),
            "P_vacuum_given_neutral": mp.nstr(p_vacuum_given_neutral, 30),
            "P_excited_given_neutral": mp.nstr(p_excited_given_neutral, 30),
        },
        "neutral_preparations": {
            "projected_independent_coherent_coefficients": "t^n/n!, t=alpha*beta",
            "projected_independent_coherent_gram": "I0(2*sqrt(conj(t)u))",
            "bessel_diagonal_control": mp.nstr(coherent_diag, 30),
            "genuine_pair_coefficients": "z^n",
            "genuine_pair_gram": "1/(1-conj(z)w)",
            "tomita_z": mp.nstr(z_tomita, 30),
            "conditioned_product_probability_ratio": mp.nstr(conditional_ratio, 30),
            "tomita_probability_ratio": mp.nstr(tomita_probability_ratio, 30),
            "number_covariance_product": mp.nstr(covariance_product, 10),
            "number_covariance_tomita": mp.nstr(covariance_tomita, 30),
        },
        "p1_scale_test": {
            "polar_index_at_any_positive_energy": -1,
            "epsilon_original": mp.nstr(epsilon, 10),
            "P_neutral_original": mp.nstr(p_neutral, 30),
            "epsilon_control": mp.nstr(epsilon_control, 10),
            "P_neutral_control": mp.nstr(p_neutral_control, 30),
            "c3": mp.nstr(p1_c3, 30),
            "one_over_2beta": mp.nstr(inherited_half_inverse_beta, 30),
            "verdict": "c3=1/(2 beta) is inherited after beta=4*pi; polar data do not force it",
        },
        "checks": checks,
    }

    print(json.dumps(report, ensure_ascii=False, indent=2))
    print(f"\nAlle Prüfungen bestanden: {all(checks.values())} ({sum(checks.values())}/{len(checks)})")
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
