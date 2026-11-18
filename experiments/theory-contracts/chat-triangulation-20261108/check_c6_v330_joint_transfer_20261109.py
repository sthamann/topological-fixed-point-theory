"""Check the finite, source-linked C6 to v330 transfer comparison."""
from fractions import Fraction
from itertools import product
import json
import math
from pathlib import Path
import sys

import numpy as np
import sympy as sp

root = Path(sys.argv[1])
sys.path.insert(0, str(root / "verification"))
from v17_hexagonal_resolvent import U6_power
from v330_qgamb_admissible_measure import symmetric_transfer, correlator


def require(label, condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)


checks = []
shift = U6_power(1)
identity = sp.eye(6)
laplacians = (
    2 * identity - shift - shift.T,
    2 * identity - shift**2 - shift.T**2,
    identity - shift**3,
)
pairs = ((1, 2), (4, 5), (0, 3))
embedding = sp.zeros(6, 3)
for column, pair in enumerate(pairs):
    for row in pair:
        embedding[row, column] = 1 / sp.sqrt(2)
reflect = sp.zeros(6)
for index in range(6):
    reflect[index, (3 - index) % 6] = 1
require("Original C6 reflection quotient and original democratic vacuum",
        embedding.T * embedding == sp.eye(3)
        and embedding * embedding.T == (identity + reflect) / 2
        and embedding * sp.ones(3, 1) / sp.sqrt(3) == sp.ones(6, 1) / sp.sqrt(6))
require("Each actual shift Laplacian preserves the quotient",
        all(lap * embedding == embedding * (embedding.T * lap * embedding)
            for lap in laplacians))

slow = sp.Matrix([1, -1, 0])
fast = sp.Matrix([1, 1, -2])
require("Quotient energies of each shift Laplacian",
        [embedding.T * lap * embedding * slow for lap in laplacians]
        == [slow, 3 * slow, 2 * slow]
        and [embedding.T * lap * embedding * fast for lap in laplacians]
        == [3 * fast, 3 * fast, sp.zeros(3, 1)])

shifted = [shift**i * sp.diag(0, 1, 1, 0, 0, 0) * shift.T**i for i in range(6)]
site_projections = [shifted[i] * shifted[(i + 1) % 6] for i in range(6)]
require("Oriented clock and first coarse label recover all original sites",
        len({tuple(matrix.diagonal()) for matrix in site_projections}) == 6
        and all(matrix.trace() == 1 for matrix in site_projections)
        and sum(site_projections, sp.zeros(6)) == identity)
require("Oriented clock cannot be kept in the three-state quotient",
        (identity - embedding * embedding.T) * shift * embedding != sp.zeros(6, 3))
require("Clock alone closes the three-dimensional space at five dimensions",
        sp.Matrix.hstack(*(shift**i * embedding for i in range(6))).rank() == 5)

r1, r2, r3, es, ef = sp.symbols("r1 r2 r3 E_s E_f", real=True)
require("All reversible cyclic lifts with this quotient",
        sp.solve((r1 + 3*r2 + 2*r3 - es, 3*r1 + 3*r2 - ef), (r1, r2))
        == {r1: r3 + (ef - es) / 2, r2: (3*es - ef) / 6 - r3})

t3, omega3, lam_s, lam_f = symmetric_transfer()
require("Original v330 eigenvectors and eigenvalues",
        np.allclose(t3 @ np.array([1,-1,0]), lam_s * np.array([1,-1,0]))
        and np.allclose(t3 @ np.array([1,1,-2]), lam_f * np.array([1,1,-2]))
        and lam_s == (2/3)**6 and lam_f == (1/3)**6)
lambda_s_exact, lambda_f_exact = Fraction(64, 729), Fraction(1, 729)
require("Original v330 exact rate ratio exceeds two", lambda_f_exact < lambda_s_exact**2)
require("Exact single-site exclusion margin",
        2 * (lambda_s_exact**2 - lambda_f_exact) / 36 > Fraction(1, 4096))

sqrt_lo, sqrt_hi = Fraction(998, 10), Fraction(9981, 100)
require("Exact rational root bracket",
        sqrt_lo**2 < 9961 < sqrt_hi**2
        and Fraction(59,100)**6 < (794-7*sqrt_hi)/2187
        and (794-7*sqrt_lo)/2187 < Fraction(3,5)**6)
delta_lo, delta_hi = Fraction(59, 100), Fraction(3, 5)
cusps = (Fraction(1), Fraction(1,3), Fraction(2,3))
for y in cusps:
    require(f"Exact r>2 at y={y}",
            y * delta_lo > 2 * max(abs(y-delta_lo), abs(y-delta_hi))**2)

cusp_min_upper = cusps[1] * delta_hi / (delta_lo - cusps[1])**2
require("Exact r at y=1/3 is smaller than at both other cusps",
        all(cusp_min_upper < y * delta_lo / max(abs(y-delta_lo), abs(y-delta_hi))**2
            for y in (cusps[0], cusps[2])))

a, b = -math.log(lam_s), -math.log(lam_f)
jump_1, jump_2 = 3 * math.log(2), math.log(9/8)
require("Shortest cyclic lift has positive one- and two-step rates",
        jump_1 > 0 and jump_2 > 0 and abs(jump_1 + 3*jump_2-a) < 1e-13
        and abs(3*(jump_1+jump_2)-b) < 1e-13)

def transfer(generator, time):
    energies, vectors = np.linalg.eigh(generator)
    return (vectors * np.exp(-time * energies)) @ vectors.T


v = np.asarray(embedding, dtype=float)
omega6 = np.ones(6) / np.sqrt(6)
laps = [np.asarray(lap, dtype=float) for lap in laplacians]
h3 = v.T @ (jump_1 * laps[0] + jump_2 * laps[1]) @ v
require("Shortest cyclic lift gives original v330 transfer entry by entry",
        np.max(np.abs(transfer(h3, 1) - t3)) < 3e-15)

all_words_error = 0.0
lifts = []
alternating = np.array([(-1)**i for i in range(6)])/np.sqrt(6)
require("Opposite jump varies only the original alternating mode",
        np.allclose(laps[0]-laps[1]+laps[2], 6*np.outer(alternating, alternating)))
for opposite in (0.0, jump_2 / 2, jump_2):
    h6 = (jump_1 + opposite)*laps[0] + (jump_2-opposite)*laps[1] + opposite*laps[2]
    require(f"Native six-coordinate Markov generator, opposite={opposite:.6f}",
            np.max(h6 - np.diag(np.diag(h6))) <= 0 and np.allclose(h6.sum(axis=1), 0))
    transfers3 = [transfer(h3, time) for time in (.2, .47, .9)]
    transfers6 = [transfer(h6, time) for time in (.2, .47, .9)]
    for labels in product(range(3), repeat=4):
        word3, word6 = np.eye(3), np.eye(6)
        for k, label in enumerate(labels):
            d3 = np.diag(np.eye(3)[label])
            d6 = np.diag(np.sum(v[:, label:label+1]**2, axis=1)*2)
            if k:
                word3 = word3 @ transfers3[k-1]
                word6 = word6 @ transfers6[k-1]
            word3 = word3 @ d3
            word6 = word6 @ d6
        all_words_error = max(all_words_error,
                              abs(omega3 @ word3 @ omega3 - omega6 @ word6 @ omega6))
    site3 = 4*jump_1+6*opposite
    site_return = (2*lam_s + 2*lam_f + math.exp(-site3))/36
    probe = np.diag([1,0,0,0,0,0])
    require(f"Native single-site return resolves missing mode, opposite={opposite:.6f}",
            abs(correlator(transfer(h6,1),omega6,probe,1)-site_return) < 3e-15)
    lifts.append({"opposite_jump": opposite, "hidden_rate": site3,
                  "native_site_connected_return_t1": site_return})
require("All 243 ordered diagonal label words preserve the actual stationary response",
        all_words_error < 3e-15)

delta = ((794-7*math.sqrt(9961))/2187)**(1/6)
cusp_checks = []
for y in cusps:
    r = float(y)*delta/(float(y)-delta)**2
    E = [math.log1p(k*r) for k in (1,3,4)]
    K = np.eye(6) + r * laps[0]
    native_rates, vectors = np.linalg.eigh(K)
    native_h = (vectors * np.log(native_rates)) @ vectors.T
    target_h = jump_1*laps[0]+jump_2*laps[1]
    residual = target_h - a/E[0] * native_h
    residual_energies = np.linalg.eigvalsh(residual)
    require(f"At matched slow rate the shortest lift needs a rank-three positive correction at y={y}",
            np.allclose(residual_energies[:3], 0, atol=1e-12)
            and residual_energies[3] > 0)
    require(f"The correction is not an additive positive-jump generator at y={y}",
            np.max(residual - np.diag(np.diag(residual))) > 0)
    require(f"Actual C6 normalized Green spectrum at y={y}",
            np.allclose(np.linalg.eigvalsh(np.linalg.inv(K)),
                        sorted([1,*[math.exp(-E[0])]*2,*[math.exp(-E[1])]*2,math.exp(-E[2])])))
    scaled = [energy*a/E[0] for energy in E]
    site_return = (2*lam_s + 2*math.exp(-scaled[1]) + math.exp(-scaled[2]))/36
    require(f"Native single-site return excludes all positive cyclic lifts at y={y}",
            site_return > max(row["native_site_connected_return_t1"] for row in lifts))
    cusp_checks.append({"y": str(y), "r": r, "max_log_rate_ratio": E[-1]/E[0],
                        "even_rates_at_fixed_v330_slow": scaled[:2],
                        "site_return_at_fixed_v330_slow": site_return,
                        "shortest_lift_correction_energies": residual_energies[3:].tolist()})
rho = max(row["max_log_rate_ratio"] for row in cusp_checks)
print(json.dumps({
    "check_count": len(checks), "checks": checks, "v330_rates": [a, b],
    "v330_log_rate_ratio": b/a, "original_c6_rate_ceiling_numeric": rho,
    "minimax_log_energy_defect_lower_bound": (b-rho*a)/(1+rho),
    "minimum_positive_fast_rate_defect_with_matched_slow": b-rho*a,
    "shortest_lift_jump_rates": [jump_1, jump_2, 0.0],
    "max_error_all_words": float(all_words_error), "lift_diagnostics": lifts,
    "original_c6": cusp_checks,
}, indent=2))
