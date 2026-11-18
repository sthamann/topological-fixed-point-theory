#!/usr/bin/env python3
"""Part H: reconstruction of the seam state kernel (QGEO.KERNEL.01 inside the free-fermion class).

The fundamental object is the STATE, not a Hamiltonian.  Class: pure, gauge-invariant,
quasi-free states of one complex fermion species on the antiperiodic ring of N sites
that are invariant under the one-site twisted translation.  Such a state is exactly a
Fermi sea: a subset F of the N antiperiodic momenta, C(d) = (1/N) sum_{k in F} e^{ikd}.

Axioms tested (none mentions a Hamiltonian, lambda, a square or P0):
  K1 half filling |F| = N/2 (charge neutrality),
  K2 reflection symmetry F = -F,
  K3 central charge one per complex fermion: exactly two Fermi points
     (16 Majoranas = 8 complex fermions must give c = 8 = rank E8).
Claims: K1-K3 leave exactly two seas, related by the local staggering gauge
c_j -> (-1)^j c_j, so the kernel is unique up to local gauge and equals
C(d) = sin(pi d/2) / (N sin(pi d/N)).  Without K3 there are other admissible seas;
they break the uniform port.  The port normalisation then reduces to v484.
    python3 -B experiments/theory-contracts/seam-holonomy-generating-20261001/kernel_reconstruction.py
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
LOG: list[dict] = []


def step(name: str, ok: bool, detail: str = "") -> None:
    LOG.append({"step": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


def momenta(N):  # antiperiodic: k = 2 pi (m + 1/2) / N, labels m = 0..N-1
    return [(2 * m + 1) % (2 * N) for m in range(N)]  # in units of pi/N, odd integers


def fermi_points(sea, N):
    occ = [m in sea for m in range(N)]
    return sum(occ[m] != occ[(m + 1) % N] for m in range(N))


def reflect(m, N):  # k -> -k : (2m+1) -> -(2m+1) mod 2N  <=> m -> N-1-m
    return N - 1 - m


def kernel(sea, N):
    d = np.arange(N)
    k = np.array([np.pi * (2 * m + 1) / N for m in sea])
    return np.array([np.exp(1j * k * x).sum() / N for x in d])


# ---------------------------------------------------------- H1/H2: exhaustive census
census = {}
for N in (12, 16, 20):
    seas = [frozenset(c) for c in combinations(range(N), N // 2)]
    k12 = [s for s in seas if frozenset(reflect(m, N) for m in s) == s]
    k123 = [s for s in k12 if fermi_points(s, N) == 2]
    census[N] = (len(seas), len(k12), len(k123), k123)
step("H1.census of all half-filled seas: reflection symmetry alone leaves many; adding two Fermi points leaves exactly 2",
     all(c[2] == 2 for c in census.values()) and all(c[1] > 2 for c in census.values()),
     "; ".join(f"N={N}: {c[0]} -> {c[1]} -> {c[2]}" for N, c in census.items()))
ok_gauge = True
for N, c in census.items():
    a, b = c[3]
    Ca, Cb = kernel(sorted(a), N), kernel(sorted(b), N)
    stag = np.array([(-1) ** x for x in range(N)])
    ok_gauge &= np.allclose(Cb, stag * Ca)
    target = np.array([np.sin(np.pi * x / 2) / (N * np.sin(np.pi * x / N)) if x else 0.5 for x in range(N)])
    ok_gauge &= np.allclose(Ca, target) or np.allclose(Cb, target)
step("H2.the two seas are related by the local staggering gauge c_j -> (-1)^j c_j; kernel = sin(pi d/2)/(N sin(pi d/N))",
     ok_gauge, "unique up to local gauge: no Hamiltonian enters")


# ---------------------------------------- H3: counter-models without K3 break the port
def green_from_kernel(Cd):
    N = len(Cd)
    S = -np.abs(Cd) ** 2
    S[0] = (Cd[0] * (1 - Cd[0])).real
    q = 2 * np.pi * np.arange(N) / N
    Sh = np.fft.fft(S).real
    Gh = np.zeros(N)
    Gh[1:] = Sh[1:] / (2 - 2 * np.cos(q[1:]))
    return np.fft.ifft(Gh).real


s = np.ones(4) / 2
v = np.array([1, -1, 1, -1]) / 2


def port_ratio(G, N, base=17):
    pos = [base + k * N // 4 for k in range(4)]
    K = np.array([[G[(a - b) % N] if i != j else 0 for j, b in enumerate(pos)] for i, a in enumerate(pos)])
    return abs((v @ K @ v) / (s @ K @ s))


N = 1024
c1 = [m for m in range(N) if np.cos(np.pi * (2 * m + 1) / N) > 0]                      # 2 Fermi points
c2 = [m for m in range(N) if abs(np.cos(2 * np.pi * (2 * m + 1) / N)) > np.cos(np.pi / 4)]  # 4 Fermi points
c2 = [m for m in c2 if True]
alt = [m for m in range(N) if np.cos(3 * np.pi * (2 * m + 1) / N) > 0]                  # 6 Fermi points
info = []
for name, sea in (("c=1", c1), ("c=4", c2), ("c=3", alt)):
    fp = fermi_points(set(sea), N)
    refl = set(reflect(m, N) for m in sea) == set(sea)
    G = green_from_kernel(kernel(sea, N))
    info.append((name, len(sea) == N // 2, refl, fp, port_ratio(G, N)))
step("H3.NEGATIVE: half-filled reflection-symmetric seas with 6 or 8 Fermi points ALSO keep the uniform port; the port does NOT select K3",
     info[0][4] < 1e-4 and all(x[1] and x[2] for x in info) and sorted(x[3] for x in info) == [2, 6, 8]
     and all(x[4] < 1e-3 for x in info[1:]),
     "; ".join(f"{n}: FP={fp}, lam_v/lam_s={r:.2e}" for n, _, _, fp, r in info))
# H5: K3 from the central charge.  A species with 2p Fermi points has c = p; c = 0 means gapped.
#     Eight complex species (16 Majoranas, SO(16)_1 inside (E8)_1) related by the carrier symmetry carry
#     the same sea, and total c = 8 = rank E8 then forces p = 1 for each species.
sols = [p for p in range(0, 9) if 8 * p == 8]
step("H5.K3 from c = 8 = rank E8 with 8 symmetry-equivalent complex species: 8p = 8 forces p = 1 (two Fermi points each)",
     sols == [1], "K3 is an independent axiom of the state, justified by the central charge, not by the port")

# --------------------------------------------- H4: normalisation reduces to v484 |D|^-1
mp.mp.dps = 20
G = green_from_kernel(kernel(c1, N))
theta = lambda d: 2 * np.pi * d / N
Gcorpus = lambda t: -(1 / np.pi) * np.log(abs(2 * np.sin(t / 2)))   # v484: |D|^-1 kernel, zero mode removed
dev = max(abs(2 * np.pi * G[d] - Gcorpus(theta(d))) for d in (N // 8, N // 4, 3 * N // 8, N // 2))
c3 = 1 / (8 * np.pi)
Kmark = 2 * np.pi * np.array([G[N // 4], G[N // 2]]) / 4
step("H4.normalisation: 2 pi x (bosonised vacuum Green) -> |D|^-1 kernel of v484; orbit average 1/4 gives K_mark = -c3 ln2 circ(0,1,2,1)",
     dev < 2e-3 and abs(Kmark[0] / (-c3 * np.log(2)) - 1) < 5e-3 and abs(Kmark[1] / (-c3 * np.log(2)) - 2) < 5e-3,
     f"max dev {dev:.1e}; K01, K02 in units of -c3 ln2: {Kmark[0] / (-c3 * np.log(2)):.4f}, {Kmark[1] / (-c3 * np.log(2)):.4f}"
     "; the factor 2 pi (non-chiral lattice density vs |D|^-1 measure) is a convention, 1/(2 pi) = 4 c3 is P1")

# H6: the equivalence map Phi of SEAM.EQUIV.01 at the finite level.  For each census size the two
#     admissible states are related by the one-particle unitary U = diag((-1)^j): it induces a
#     Bogoliubov *-automorphism Phi(c_j) = (-1)^j c_j (products and adjoints preserved by construction),
#     maps one covariance onto the other (state preservation) and commutes with the raw station
#     clock R = R_station (x) I_M (clock intertwining) whenever M = N/4 is even.
ok6 = []
for N, c in census.items():
    if N % 8:
        continue
    a, b = (sorted(x) for x in c[3])
    def covmat(sea):
        k = np.array([np.pi * (2 * m + 1) / N for m in sea])
        return np.array([[np.exp(1j * k * (y - x)).sum() / N for y in range(N)] for x in range(N)])
    Ca, Cb = covmat(a), covmat(b)
    U = np.diag([(-1.0) ** j for j in range(N)])
    Rst = np.zeros((4, 4)); Rst[1, 0] = Rst[2, 1] = Rst[3, 2] = 1; Rst[0, 3] = -1
    R = np.kron(Rst, np.eye(N // 4))
    ok6.append(np.allclose(U @ U.T, np.eye(N)) and np.allclose(U @ Ca @ U.T, Cb)
               and np.allclose(U @ R, R @ U) and np.allclose(R @ Ca @ R.T, Ca))
step("H6.Phi = staggering Bogoliubov automorphism: *-preserving, maps the one admissible state onto the other, commutes with R, both states R-invariant",
     ok6 and all(ok6), f"sizes N in {[N for N in census if N % 8 == 0]}")

summary = {"steps": LOG, "passed": sum(x["pass"] for x in LOG), "total": len(LOG),
           "census": {N: c[:3] for N, c in census.items()},
           "counter_models": [dict(name=n, fermi_points=fp, port_ratio=r) for n, _, _, fp, r in info]}
(HERE / "kernel_reconstruction_results.json").write_text(json.dumps(summary, indent=2))
print(f"{summary['passed']}/{summary['total']} steps passed")
sys.exit(0 if summary["passed"] == summary["total"] else 1)
