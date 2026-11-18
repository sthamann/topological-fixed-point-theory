"""Q1 - what can be tested on a real two-qubit system.

One TFPT register is C^4 = two qubits; the 60 source rays are the 60 two-qubit stabiliser
states.  The whole quartic structure comes from ONE measurable fact: these 60 states agree
with a uniformly random state in every moment up to the third and deviate in the fourth
(40 M4 = S + P instead of 40 M4 = (8/7) S).
  (a) frame potentials t = 1..4 against the Haar value
  (b) a hardware protocol: prepare a probe state, measure it in the 15 stabiliser bases,
      estimate W = (1/60) sum_l p_l^4.  Stabiliser probe: 1/32.  Haar average: 1/35.
      Simulated with shot noise, depolarising noise and readout error (no hardware used).
  (c) the pulse reflection U = I - 2 D D^+ on two qubits: with D = singlet it is the SWAP gate;
      as oracle/diffusion pair it is one Grover step on four items.
"""
import itertools, json
import numpy as np
import tfpt_ops as T

rng = np.random.default_rng(0)
o = T.build()
stab, pauli = o["stab"], o["pauli"]
res = {}

# (a) frame potentials
ov = np.abs(stab.conj() @ stab.T) ** 2
from math import comb
res["a_frame_potential"] = {f"t={t}": dict(sixty_states=round(float((ov ** t).mean()), 10),
                                            haar=round(1 / comb(t + 3, 3), 10)) for t in (1, 2, 3, 4)}
res["a_fourth_moment_exact"] = "1/32 = (35 + 3*5)/1600 versus Haar 1/35"

# group the 60 states into 15 orthonormal bases (states sharing one stabiliser group)
ex = np.einsum("la,Aab,lb->lA", stab.conj(), pauli, stab).real
groups = {}
for l in range(60):
    groups.setdefault(tuple(np.flatnonzero(np.abs(ex[l]) > 0.5)), []).append(l)
bases = [np.array(v) for v in groups.values()]
assert len(bases) == 15 and all(len(b) == 4 for b in bases)
assert all(np.allclose(stab[b].conj() @ stab[b].T, np.eye(4), atol=1e-9) for b in bases)

def witness_exact(phi):
    return float((np.abs(stab.conj() @ phi) ** 8).mean())

T_state = np.array([1, np.exp(1j * np.pi / 4)]) / np.sqrt(2)
probes = {"stabiliser state |00>": np.array([1, 0, 0, 0], complex),
          "Bell state": np.array([1, 0, 0, 1], complex) / np.sqrt(2),
          "T (x) T magic state": np.kron(T_state, T_state)}
# most non-stabiliser two-qubit state: minimise the witness numerically
best = None
for _ in range(300):
    v = rng.standard_normal(4) + 1j * rng.standard_normal(4); v /= np.linalg.norm(v)
    for _ in range(400):
        g = rng.standard_normal(4) + 1j * rng.standard_normal(4)
        w = v + 0.03 * g; w /= np.linalg.norm(w)
        if witness_exact(w) < witness_exact(v):
            v = w
    if best is None or witness_exact(v) < witness_exact(best):
        best = v
probes["numerically most 'magic' state"] = best
haar = [witness_exact((lambda v: v / np.linalg.norm(v))(rng.standard_normal(4) + 1j * rng.standard_normal(4)))
        for _ in range(20000)]
res["b_witness_exact"] = {k: round(witness_exact(v), 8) for k, v in probes.items()}
res["b_witness_exact"]["Haar average (20000 samples)"] = round(float(np.mean(haar)), 6)
res["b_witness_exact"]["reference 1/32, 1/35"] = [round(1 / 32, 8), round(1 / 35, 8)]

def simulate(phi, shots, depol, readout, trials=400):
    """estimate W from counts in the 15 bases; unbiased estimator of p^4 per outcome."""
    est = []
    for _ in range(trials):
        tot = 0.0
        for b in bases:
            p = (1 - depol) * np.abs(stab[b].conj() @ phi) ** 2 + depol / 4
            p = (1 - readout) * p + readout / 4          # symmetric readout confusion
            n = rng.multinomial(shots, p / p.sum()).astype(float)
            tot += (n * (n - 1) * (n - 2) * (n - 3)).sum() / (shots * (shots - 1) * (shots - 2) * (shots - 3))
        est.append(tot / 60)
    return float(np.mean(est)), float(np.std(est))

sim = {}
for depol, ro in ((0.0, 0.0), (0.01, 0.01), (0.03, 0.02)):
    for shots in (200, 1000, 5000):
        a = simulate(probes["stabiliser state |00>"], shots, depol, ro)
        b = simulate(probes["T (x) T magic state"], shots, depol, ro)
        sim[f"depol={depol} readout={ro} shots/basis={shots}"] = dict(
            stabiliser=[round(a[0], 6), round(a[1], 6)], magic=[round(b[0], 6), round(b[1], 6)],
            separation_sigma=round((a[0] - b[0]) / np.hypot(a[1], b[1]), 1))
res["b_simulated_protocol"] = sim

# (c) reflections as gates
singlet = np.array([0, 1, -1, 0]) / np.sqrt(2)
U = np.eye(4) - 2 * np.outer(singlet, singlet)
SWAP = np.eye(4)[[0, 2, 1, 3]]
res["c_reflection_on_singlet_is_SWAP"] = bool(np.allclose(U, SWAP))
s = np.ones(4) / 2
grover = {}
for marked in range(4):
    oracle = np.eye(4); oracle[marked, marked] = -1          # I - 2|m><m|
    diff = 2 * np.outer(s, s) - np.eye(4)                    # -(I - 2|s><s|)
    grover[marked] = round(float(abs((diff @ oracle @ s)[marked]) ** 2), 10)
res["c_one_grover_step_success_probability"] = grover
res["reading"] = ("The fourth-moment gap 1/32 vs 1/35 is real and measurable with a few thousand shots per basis even "
                  "at 3% noise.  It is the known statement that the two-qubit Clifford group is a 3-design and not a "
                  "4-design (Zhu, Kueng, Grassl, Gross 2016); confirming it tests quantum mechanics and the device, "
                  "not TFPT.  A TFPT-specific quantum prediction would need a number that standard theory does not give.")
json.dump(res, open("out/q1_two_qubit.json", "w"), indent=1)
print(json.dumps(res, indent=1))
