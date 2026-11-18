#!/usr/bin/env python3
"""Part G: do the uniform port line and the full access-shifting clock come out of the raw source?

Red-team rule: NO lambda = -1, NO square, NO uniform P0 is used as input.

Inputs only:
  * the raw antiperiodic ring, 4 stations x M modes, homogeneous half-filled vacuum C
    (hopping Hamiltonian = the v480 model choice),
  * the raw station clock Q = R_station (x) I_M (universalraum-mu4-pointing),
  * the requirement that Q acts on the SAME four access points (one per station).
The mark kernel is computed FROM THE VACUUM: density correlation S(d) = <n_0 n_d>_c,
lattice-bosonised G = Delta^-1 S (zero mode removed), K_ab = G(p_a - p_b), zero diagonal.

Outputs (computed, not assumed): mark positions, cross-ratio, [Q,K], E_+ = (I + Q^2)/2
from the clock, the rank-one port direction of E_+ K E_+ and its normalised form.
Controls: marks that are not a clock orbit; a non-logarithmic kernel (|C|^2).
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/port_from_source.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []
OUT: dict = {}


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


def vacuum(N, t2=0.0):
    H = np.zeros((N, N))
    for j in range(N):
        for r, t in ((1, 1.0), (2, t2)):
            if t == 0.0:
                continue
            k = (j + r) % N
            sgn = -1.0 if j + r >= N else 1.0          # antiperiodic
            H[j, k] += -t * sgn
            H[k, j] += -t * sgn
    e, v = np.linalg.eigh(H)
    occ = v[:, e < 0]
    return occ @ occ.T


def bosonised_green(C):
    """G with -Delta G = S, S(d) = <n_0 n_d>_c from the vacuum; zero mode removed."""
    N = len(C)
    S = -np.abs(C[0, :]) ** 2
    S[0] = C[0, 0] * (1 - C[0, 0])
    q = 2 * np.pi * np.arange(N) / N
    Sh = np.fft.fft(S).real
    Gh = np.zeros(N)
    Gh[1:] = Sh[1:] / (2 - 2 * np.cos(q[1:]))
    return np.fft.ifft(Gh).real


def station_clock(N, M):
    Rst = np.zeros((4, 4))
    for x in range(3):
        Rst[x + 1, x] = 1
    Rst[0, 3] = -1
    return np.kron(Rst, np.eye(M))


def site_image(Qfull, j):
    return int(np.nonzero(Qfull[:, j])[0][0])


def mark_kernel(G, pos):
    N = len(G)
    K = np.array([[G[(pa - pb) % N] if a != b else 0.0 for b, pb in enumerate(pos)] for a, pa in enumerate(pos)])
    return K


def cross_ratio(pos, N):
    z = [np.exp(2j * np.pi * p / N) for p in pos]
    return (z[0] - z[2]) * (z[1] - z[3]) / ((z[0] - z[3]) * (z[1] - z[2]))


Q4 = np.roll(np.eye(4), 1, axis=0)              # access shift a -> a+1 (permutation on marks)
Eplus = (np.eye(4) + Q4 @ Q4) / 2               # built from the clock, not from P0
s = np.ones(4) / 2
v = np.array([1, -1, 1, -1]) / 2

results = {}
for N in (256, 512, 1024, 2048):
    M = N // 4
    C = vacuum(N)
    G = bosonised_green(C)
    Qf = station_clock(N, M)
    # G1: access points that the raw clock maps among themselves (no positions assumed)
    base = 17                                    # arbitrary access point in station 0
    pos = [base]
    for _ in range(3):
        pos.append(site_image(Qf, pos[-1]))
    closes = site_image(Qf, pos[-1]) == base
    K = mark_kernel(G, pos)
    comm = np.abs(Q4 @ K - K @ Q4).max()
    P = Eplus @ K @ Eplus
    ev, evec = np.linalg.eigh(P)
    lam_s, lam_v = s @ K @ s, v @ K @ v
    results[N] = dict(pos=pos, closes=closes, cr=cross_ratio(pos, N), comm=comm,
                      lam_s=lam_s, lam_v=lam_v, ratio=abs(lam_v / lam_s), K=K)
r = results[2048]
step("G1.the raw clock's orbit of one access point closes after 4 steps: positions p, p+N/4, p+N/2, p+3N/4 (computed)",
     all(x["closes"] for x in results.values()) and r["pos"] == [17, 17 + 512, 17 + 1024, 17 + 1536])
step("G2.cross-ratio of the computed access points = 2 (harmonic): lambda = -1 is an OUTPUT",
     all(abs(x["cr"] - 2) < 1e-12 for x in results.values()))
step("G3.full clock covariance [Q,K] = 0 on the vacuum mark kernel (all N)",
     all(x["comm"] < 1e-12 for x in results.values()), f"max {max(x['comm'] for x in results.values()):.1e}")
ratios = [results[N]["ratio"] for N in sorted(results)]
step("G4.uniform port: E_+ K E_+ is rank one along s; alternating eigenvalue / uniform -> 0 as N grows",
     ratios[-1] < 1e-3 and all(b < a for a, b in zip(ratios, ratios[1:])),
     "lam_v/lam_s = " + ", ".join(f"{x:.2e}" for x in ratios))
Kn = r["K"] / r["K"][0, 1]
step("G5.normalised vacuum kernel -> circ(0,1,2,1) (the original C_mark of v484/v485), hence E_+ C E_+ = 4 P0 as output",
     abs(Kn[0, 2] - 2) < 2e-3 and abs(Kn[0, 3] - 1) < 1e-12,
     f"K/K01 row0 = {np.round(Kn[0], 4).tolist()}")
Cm = np.array([[0, 1, 2, 1], [1, 0, 1, 2], [2, 1, 0, 1], [1, 2, 1, 0]])
P0 = np.outer(s, s)
step("G6.with that kernel: E_+ C E_+ = 4 P0 and det(I-uC) = (1-4u)(1+2u)^2 (spectrum 4,-2,0,-2)",
     np.allclose(Eplus @ Cm @ Eplus, 4 * P0) and np.allclose(sorted(np.linalg.eigvalsh(Cm)), [-2, -2, 0, 4]))
# controls
N = 2048; C = vacuum(N); G = bosonised_green(C)
bad = [17, 17 + 380, 17 + 1024, 17 + 1404]       # one point per station, not a clock orbit
Kb = mark_kernel(G, bad)
step("G7.control: one point per station but NOT a clock orbit -> [Q,K] != 0, non-harmonic cross-ratio, port not along s",
     np.abs(Q4 @ Kb - Kb @ Q4).max() > 1e-2 and abs(cross_ratio(bad, N) - 2) > 1e-2
     and abs((v @ Kb @ v) / (s @ Kb @ s)) > 1e-2,
     f"cr={cross_ratio(bad, N).real:.4f}, lam_v/lam_s={(v @ Kb @ v) / (s @ Kb @ s):.3f}")
Kc = mark_kernel(np.abs(C[0, :]) ** 2, results[2048]["pos"])
step("G8.control: square marks but the non-logarithmic current kernel |C|^2 -> even port has rank 2 (log/Green kernel is essential)",
     abs((v @ Kc @ v) / (s @ Kc @ s)) > 0.1, f"lam_v/lam_s = {(v @ Kc @ v) / (s @ Kc @ s):.3f}")

# G9: universality -- a next-nearest-neighbour hopping (still one Fermi pair) gives the same port
rob = []
for N in (512, 2048):
    Cr = vacuum(N, t2=0.3)
    Gr = bosonised_green(Cr)
    pos = [17 + k * N // 4 for k in range(4)]
    Kr = mark_kernel(Gr, pos)
    rob.append((abs((v @ Kr @ v) / (s @ Kr @ s)), np.abs(Q4 @ Kr - Kr @ Q4).max(), Kr[0, 2] / Kr[0, 1]))
step("G9.universality: with next-nearest hopping t2=0.3 the port is again rank one along s ([Q,K]=0, ratio -> 0, K02/K01 -> 2)",
     rob[-1][0] < rob[0][0] and rob[-1][0] < 1e-3 and all(x[1] < 1e-12 for x in rob) and abs(rob[-1][2] - 2) < 5e-3,
     "; ".join(f"ratio {a:.1e}, K02/K01 {c:.4f}" for a, _, c in rob))

OUT = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
       "ratios": dict(zip(sorted(results), ratios)), "K_normalised_N2048": Kn[0].tolist()}
(HERE / "port_from_source_results.json").write_text(json.dumps(OUT, indent=2, default=str))
print(f"{OUT['passed']}/{OUT['total']} steps passed")
sys.exit(0 if OUT["passed"] == OUT["total"] else 1)
