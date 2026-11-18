"""S1, S2, S3 - what the three signature tests could resolve.  No real data is used here.

S1  Cubic lattice pattern on the sky.  Phenomenological branch only (no map from lattice spacing
    to arrival pattern is fixed): intensity 1 + kappa * K4(R n) with the cubic harmonic
    K4(n) = n_x^4 + n_y^4 + n_z^4 - 3/5, on top of a dipole, seen through the exposure of a
    surface detector at latitude -35.2 deg with zenith angles up to 60 deg (Sommers 2001 formula).
    Statistic: exposure-corrected sum of K4 over the events, maximised over lattice orientation.
    The null distribution includes the orientation scan.  Output: smallest kappa that is found
    with 90 % power at a global p of 0.001, for several event numbers, with and without a
    split into exploration and confirmation halves.
S2  Message in constants.  How many independent, stable sign bits are needed for a global
    p < 0.001, with decoders tried, tolerated errors, dropouts and correlated measurements.
S3  Causally separated intervention.  Block randomisation test on a phase shift: number of
    blocks for 90 % power at p < 0.001, calibration under drift, behaviour of the placebo
    assignments, and recovery of an injected local effect.
"""
import json
from math import comb
import numpy as np
from scipy.stats import norm
from scipy.spatial.transform import Rotation

rng = np.random.default_rng(0)
res = {}

# ------------------------------------------------------------------ S1
import sys, os
if '--skip-s1' in sys.argv and os.path.exists('out/s1_s2_s3_signatures.json'):
    res['S1'] = json.load(open('out/s1_s2_s3_signatures.json'))['S1']
else:
    LAT, ZMAX = np.radians(-35.2), np.radians(60.0)
    def exposure(dec):
        xi = (np.cos(ZMAX) - np.sin(LAT) * np.sin(dec)) / (np.cos(LAT) * np.cos(dec))
        am = np.where(xi > 1, 0.0, np.where(xi < -1, np.pi, np.arccos(np.clip(xi, -1, 1))))
        return np.clip(np.cos(LAT) * np.cos(dec) * np.sin(am) + am * np.sin(LAT) * np.sin(dec), 0, None)

    def K4(n):
        return (n ** 4).sum(-1) - 0.6

    # orientation scan: 200 random orientations of the lattice axes
    ROT = Rotation.random(200, random_state=1).as_matrix()
    DIP = 0.065 * np.array([np.cos(np.radians(-24)) * np.cos(np.radians(100)),
                            np.cos(np.radians(-24)) * np.sin(np.radians(100)), np.sin(np.radians(-24))])

    def sample_sky(n, kappa, Rtrue, g, dip=DIP):
        out, bound = [], 1.85 * 1.07 * (1 + 0.4 * kappa)
        while sum(len(o) for o in out) < n:
            v = g.standard_normal((3 * n, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True)
            w = exposure(np.arcsin(v[:, 2])) * (1 + v @ dip) * (1 + kappa * K4(v @ Rtrue.T))
            assert w.max() <= bound
            out.append(v[g.random(len(v)) < w / bound])
        return np.concatenate(out)[:n]

    def reference(dip):
        ref = sample_sky(400000, 0.0, np.eye(3), np.random.default_rng(5), dip)
        k = np.stack([K4(ref @ Rm.T) for Rm in ROT], 1)
        return k.mean(0), k.std(0)
    MU, SD = reference(DIP)                       # null model: exposure and the known dipole
    MU0, SD0 = reference(0 * DIP)                 # null model that forgets the dipole
    def stat(ev, rots=slice(None), mu=MU, sd=SD):
        k = np.stack([K4(ev @ Rm.T) for Rm in ROT[rots]], 1)
        return (k.mean(0) - mu[rots]) / (sd[rots] / np.sqrt(len(ev)))        # z per orientation

    # null of the scan maximum (z-scores, so nearly independent of N): effective number of orientations
    null = np.array([stat(sample_sky(20000, 0.0, np.eye(3), rng)).max() for _ in range(600)])
    m_eff = float(np.log(0.5) / norm.logcdf(np.median(null)))
    thr = float(norm.isf(1 - 0.999 ** (1 / m_eff)))
    single_thr = float(norm.isf(0.001))
    s1 = {"orientations_scanned": len(ROT), "null_scan_max_median_z": round(float(np.median(null)), 2),
          "null_scan_max_99_percent_z": round(float(np.quantile(null, 0.99)), 2),
          "effective_independent_orientations": round(m_eff, 1), "threshold_z_global_p_0.001": round(thr, 2),
          "threshold_z_without_scan": round(single_thr, 2)}
    GRID = np.geomspace(0.002, 1.0, 28)
    def power(N, kappa, split, trials=30):
        hits = 0
        for t in range(trials):
            ev = sample_sky(N, kappa, Rotation.random(random_state=1000 + t).as_matrix(), rng)
            if not split:
                hits += stat(ev).max() > thr
            else:
                best = int(stat(ev[: N // 2]).argmax())
                hits += stat(ev[N // 2:], slice(best, best + 1))[0] > single_thr
        return hits / trials
    for N in (2000, 20000, 100000):
        row = {}
        for label, split in (("full sample, orientation scanned", False), ("half to find orientation, half to confirm", True)):
            lo, hi = 0, len(GRID) - 1
            if power(N, GRID[hi], split) < 0.9:
                row[label] = None; continue
            while lo < hi:
                mid = (lo + hi) // 2
                if power(N, GRID[mid], split) >= 0.9:
                    hi = mid
                else:
                    lo = mid + 1
            row[label] = round(float(GRID[lo]), 4)
        s1[f"N={N} events: smallest kappa with 90 % power"] = row
        print("S1", N, row, flush=True)
    # an unmodelled dipole: no cubic signal injected, null model without the dipole
    fooled = [float(stat(sample_sky(N, 0.0, np.eye(3), rng), mu=MU0, sd=SD0).max()) for N in (20000, 100000) for _ in range(20)]
    s1["scan_max_z_when_the_dipole_is_left_out_of_the_null_model"] = {
        "N=20000": round(float(np.median(fooled[:20])), 2), "N=100000": round(float(np.median(fooled[20:])), 2)}
    s1["scope"] = ("Sky pattern only.  Turning kappa into a lattice spacing needs a fixed emission and propagation model "
                   "(the cited proposal, arXiv:1210.1847, argues through the energy cutoff instead).  Real events and "
                   "exposure from the public Auger release were not downloaded; the 6.5 % dipole and the exposure "
                   "formula are typed in from the literature and should be re-checked against the data release.")
    res["S1"] = s1
    json.dump(res, open("out/s1_s2_s3_signatures.json", "w"), indent=1)


# ------------------------------------------------------------------ S2
def p_match(n, e):
    return sum(comb(n, k) for k in range(e + 1)) / 2 ** n
s2 = {"bits_needed_for_global_p_below_0.001": {}}
for dec in (1, 16, 1000):
    for e in (0, 1, 2):
        n = next(n for n in range(1, 200) if dec * p_match(n, e) < 1e-3)
        s2["bits_needed_for_global_p_below_0.001"][f"decoders tried={dec}, tolerated wrong bits={e}"] = n
# power when the message is real: each bit is the sign of a residual of size r sigma
def power(n, e, r, drop_at=1.0, trials=20000):
    z = r + rng.standard_normal((trials, n))
    usable = np.abs(z) > drop_at                      # bits whose interval covers zero are dropped (preregistered)
    wrong = ((z < 0) & usable).sum(1)
    nu = usable.sum(1)
    ok = np.array([p_match(int(a), e) < 1e-3 for a in range(n + 1)])[nu] & (wrong <= e)
    return round(float(ok.mean()), 3)
s2["power_if_message_is_real"] = {f"{n} bits, residual {r} sigma, {e} wrong tolerated": power(n, e, r)
                                  for n in (12, 20, 40) for r in (1.0, 2.0, 3.0) for e in (0, 2)}
# correlated measurements: effective number of independent bits
def false_alarm(n, rho, trials=400000):
    c = rng.standard_normal((trials, 1))
    z = np.sqrt(rho) * c + np.sqrt(1 - rho) * rng.standard_normal((trials, n))
    return float(((z > 0).all(1)).mean())
s2["chance_of_matching_10_preset_equal_signs_by_accident"] = {f"correlation {rho}": float(f"{false_alarm(10, rho):.2e}") for rho in (0.0, 0.3, 0.6)}
s2["available_from_this_package"] = ("One reproduced residual (the fine-structure constant, F10).  No target string and no decoder "
                                     "exist, so the message test itself is not run.")
res["S2"] = s2

# ------------------------------------------------------------------ S3
z_a, z_b = norm.isf(0.001), norm.isf(0.10)
s3 = {"blocks_per_group_for_90_percent_power_one_sided_p_0.001":
      {f"shift = {d} standard deviations": int(np.ceil(2 * (z_a + z_b) ** 2 / d ** 2)) for d in (1.0, 0.3, 0.1, 0.03, 0.01)}}
def experiment(n_blocks, delta, drift=0.0, placebo=False, n_perm=1000, paired=False):
    bits = rng.integers(0, 2, n_blocks)
    ar = np.zeros(n_blocks)
    for i in range(1, n_blocks):
        ar[i] = 0.95 * ar[i - 1] + rng.standard_normal()
    y = rng.standard_normal(n_blocks) + drift * ar + delta * bits
    use = rng.integers(0, 2, n_blocks) if placebo else bits          # placebo: bits drawn later, unrelated
    if paired:      # neighbouring blocks compared with each other: slow drift cancels
        dy = y[1::2] - y[0::2]
        T = lambda b: (dy * (b[1::2] - b[0::2])).sum()
    else:
        T = lambda b: y[b == 1].mean() - y[b == 0].mean()
    obs = T(use)
    perm = np.array([T(rng.permutation(use)) for _ in range(n_perm)])
    return (1 + (perm >= obs).sum()) / (1 + n_perm)
def rate(n_blocks, delta, level, reps=400, **kw):
    return round(float(np.mean([experiment(n_blocks, delta, **kw) < level for _ in range(reps)])), 3)
s3["simulated_randomisation_test_400_blocks"] = {
    "false alarm rate at p<0.05, no effect, no drift": rate(400, 0.0, 0.05, reps=4000, n_perm=199),
    "false alarm rate at p<0.05, no effect, strong slow drift": rate(400, 0.0, 0.05, drift=2.0, reps=4000, n_perm=199),
    "false alarm rate at p<0.05, no effect, strong slow drift, neighbouring blocks compared": rate(400, 0.0, 0.05, drift=2.0, paired=True, reps=4000, n_perm=199),
    "detection rate at p<0.001, shift 0.5 sd": rate(400, 0.5, 0.001),
    "detection rate at p<0.001, shift 0.5 sd, strong slow drift": rate(400, 0.5, 0.001, drift=2.0),
    "detection rate at p<0.001, shift 0.5 sd, strong slow drift, neighbouring blocks compared": rate(400, 0.5, 0.001, drift=2.0, paired=True),
    "detection rate at p<0.001, placebo bits, shift 0.5 sd present": rate(400, 0.5, 0.001, placebo=True),
}
s3["note"] = "false alarm rates from 4000 repetitions (expected value 0.045 with 199 permutations), detection rates from 400"
s3["scope"] = ("Statistics of the design only.  The package names no time window and no effect size for an actual "
               "candidate, and TFPT itself predicts no such dependence, so nothing is tested about nature here.")
res["S3"] = s3
json.dump(res, open("out/s1_s2_s3_signatures.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "S1"}, indent=1))
