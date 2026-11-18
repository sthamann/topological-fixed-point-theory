"""Reproduce the lepton double check and the reflected Fold shell tests."""

import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 45
repo = Path(sys.argv[1]).resolve()
verification = repo / "verification"
sys.path.insert(0, str(verification))
from v3_em_alpha import make_F
from v420_fpole_qed_running import ME, MMU, MTAU, MZ
from tfpt_constants import c3, dtop, phibase, phi0

base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
alpha = mp.findroot(make_F(c3, 41), mp.mpf("0.0073"))
damped = dtop * mp.exp(-2 * alpha)
resummed = damped * (1 - damped) ** (-mp.mpf(5) / 4)
seam = phibase + resummed
family_lengths = (5, 3, 2)
coefficients = (mp.mpf(16) / 7, mp.mpf(4) / 3, mp.mpf(7) / 6)
pole = (ME, MMU, MTAU)
corner = mp.mpf(53) / 54


def ladder(seed):
    return tuple(c * seed**n for c, n in zip(coefficients, family_lengths))


def observables(masses):
    electron, muon, tau = masses
    a, b, c = map(mp.sqrt, masses)
    return {
        "mu_over_e": muon / electron,
        "tau_over_mu": tau / muon,
        "Q": sum(masses) / (a + b + c)**2,
        "K": electron * tau**2 / muon**3,
        "azimuth_rad": mp.atan2((b - a) / mp.sqrt(2), (2*c - a - b) / mp.sqrt(6)),
    }


def mobius(delta):
    m1 = (1 + delta) / (1 - delta)
    third = mp.mpf(1) / 3
    mthird = abs((third + delta) / (third - delta))
    return 1 / (m1 * mthird)**2, mp.mpf(1), m1**2


def qed_pole_from_msbar(masses):
    # First-order inverse of v420's pole -> MSbar conversion at MZ.
    normalized = tuple(m * MTAU / masses[-1] for m in masses)
    return tuple(m * (1 + alpha/mp.pi * (1 + mp.mpf(3)/4 * mp.log(MZ**2 / m**2)))
                 for m in normalized)


target = observables(pole)
seed_variants = {
    "tree": phibase, "retained_seed": phi0, "alpha_damped": phibase + damped,
    "alpha_resummed": seam, "cabibbo_squared": phi0*(1-phi0),
    "corner": corner*phi0, "corner_times_alpha_seam": corner*seam,
}
rows = {name: observables(ladder(u)) for name, u in seed_variants.items()}
assert all(abs(row["K"] - mp.mpf(21)/16) < mp.mpf("1e-37") for row in rows.values())
leading_delta = mp.mpf(3)/5
delta = leading_delta + phi0/6
delta_anchor = (mp.sqrt(MTAU/MMU) - 1)/(mp.sqrt(MTAU/MMU) + 1)
rows.update({
    "v27_mobius_leading": observables(mobius(leading_delta)),
    "v27_mobius_geometric": observables(mobius(delta)),
    "v27_mobius_alpha_seam": observables(mobius(leading_delta + seam/6)),
    "v27_mobius_empirical_anchor": observables(mobius(delta_anchor)),
    "v27_mobius_geometric_plus_qed": observables(qed_pole_from_msbar(mobius(delta))),
    "v420_qed_only": observables(qed_pole_from_msbar(ladder(phi0))),
    "corner_plus_qed": observables(qed_pole_from_msbar(ladder(corner*phi0))),
    "frozen_v420_poles": target,
})
# Determine each observable's delta separately on the same local branch.
# This is a diagnostic inverse calculation, never an input to the prediction.
delta_needed = {name: mp.findroot(lambda x: observables(mobius(x))[name] - target[name],
                                  (mp.mpf("0.604"), mp.mpf("0.612")))
                for name in ("Q", "K", "mu_over_e", "tau_over_mu")}
assert abs(delta_needed["Q"] - delta_needed["K"]) > mp.mpf("0.001")

# Hybrid diagnostic. Keep the historical geometry's unfitted tau/mu and
# impose the exact Koide cone. Select the unique 0 < sqrt(me/mmu) < 1 root.
tau_root = mp.sqrt(rows["v27_mobius_geometric"]["tau_over_mu"])
q = mp.mpf(2)/3
poly = (1-q, -2*q*(1+tau_root), 1+tau_root**2-q*(1+tau_root)**2)
root, = [r for r in mp.polyroots(poly) if 0 < r < 1]
rows["v27_tau_ratio_on_exact_koide_cone"] = observables((root**2, 1, tau_root**2))

# Old block table, using its printed formula and no recovered factors.
gamma = mp.mpf(5)/6
lam = gamma / mp.log(mp.mpf(248)/60)
planck = mp.mpf("1.22089e19")
blocks = {}
for name, n, rank, k, printed in (
    ("electroweak", 12, 2, mp.mpf(41)/32, "251.1"),
    ("hadronic", 15, 5, mp.mpf(3)/2, "0.968"),
    ("pion", 16, 5, mp.mpf(51)/32, "0.088"),
    ("seesaw", 5, 4, mp.mpf(1)/8, "1.31e15"),
    ("pq", 10, 1, mp.mpf(1)/2, "8.86e10"),
):
    phi = phi0 * mp.exp(-gamma) * (mp.mpf(60-2*n)/58)**lam
    beta = mp.mpf(8-rank)/8
    value = mp.pi*c3 * mp.exp(-beta*mp.pi*c3-k/c3) * planck * phi
    blocks[name] = {"phi_n": phi, "calculated_GeV": value,
                    "printed_GeV": mp.mpf(printed), "ratio": value/mp.mpf(printed)}

# Independent second mark on the exact original Fold shell.
shell = json.loads((base / "one-original-fold-link-results.json").read_text())
even, spinor = (set(tuple(map(sp.Rational, row["p"])) for row in shell[key])
                for key in ("even_shell", "spinor_shell"))
half = sp.Rational(1, 2)
zeta = (1 + sp.I) / sp.sqrt(2)


def reflect(p):
    return tuple(-x for x in reversed(p))


def affine_partner(p):
    return tuple(half + x for x in reflect(p))


def clock_energy(p):
    return sum(4*x*x+(k-4)*x for k, x in enumerate(p))


def harmonic(p, d):
    return sp.expand(sum(x*zeta**(d*k) for k, x in enumerate(p)))


def symmetric_cartan_square(p):
    return p[0]**2+p[7]**2


assert len(even) == len(spinor) == 7
assert {affine_partner(p) for p in spinor} == even
assert all(sum(p) == 2 and clock_energy(p) == 6 for p in even | spinor)
for p in even | spinor:
    assert affine_partner(affine_partner(p)) == p
    assert symmetric_cartan_square(p) == symmetric_cartan_square(reflect(p))
    for d in range(1, 8):
        assert sp.simplify(harmonic(affine_partner(p), d) - harmonic(reflect(p), d)) == 0
spinor_values = {symmetric_cartan_square(p) for p in spinor}
even_values = {symmetric_cartan_square(p) for p in even}
assert all(v % 1 == half for v in spinor_values)
assert all(v % 1 == 0 for v in even_values)

orbit = json.loads((base / "direct-fold-word-results.json").read_text())["source_orbit"]
middle = tuple(map(sp.Rational, orbit[1]))
assert affine_partner(middle) == (0, 0, 1, 0, 1, 0, 0, 0)
assert symmetric_cartan_square(middle) == half
assert symmetric_cartan_square(affine_partner(middle)) == 0
invariant_reader = lambda p: sp.expand(harmonic(p, 1) + harmonic(reflect(p), 1))
assert [p for p in spinor if invariant_reader(p) == invariant_reader(middle)] == [middle]

results = {
    "alpha": {"root": alpha, "inverse": 1/alpha, "phi_base": phibase,
              "topological_term": dtop, "damping_shift": damped-dtop,
              "resummation_shift": resummed-damped, "retained_seed": phi0,
              "seam_seed": seam, "one_loop_scale_alpha_over_pi": alpha/mp.pi},
    "leptons": rows,
    "v27_geometric_delta": delta,
    "v27_delta_needed_each": delta_needed,
    "corner_to_poles": {
        "delta_log_K": mp.log(target["K"]/rows["corner"]["K"]),
        "delta_log_mu_over_e": mp.log(target["mu_over_e"]/rows["corner"]["mu_over_e"]),
        "delta_log_tau_over_mu": mp.log(target["tau_over_mu"]/rows["corner"]["tau_over_mu"]),
        "delta_azimuth_rad": target["azimuth_rad"]-rows["corner"]["azimuth_rad"],
        "delta_azimuth_degrees": (target["azimuth_rad"]-rows["corner"]["azimuth_rad"])*180/mp.pi,
    },
    "historical_block_table": blocks,
    "fold": {"exact_paired_states": len(spinor), "matching_nonzero_harmonics": 7,
             "spinor_cartan_square_values": sorted(spinor_values),
             "even_cartan_square_values": sorted(even_values),
             "joint_reader_unique_target_in_given_shell": True},
}


def format_number(value):
    return mp.nstr(value, 22) if isinstance(value, mp.mpf) else str(value)


print(json.dumps(results, indent=2, default=format_number))
