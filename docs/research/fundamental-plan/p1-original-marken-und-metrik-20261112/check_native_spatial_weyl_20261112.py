"""Bounded potential Weyl controls on the native Q spatial circle."""
import math
from pathlib import Path
import sys

import numpy as np
import sympy as sp

root = Path(sys.argv[1]).resolve()
folder = root / "experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q"
sys.path.insert(0, str(folder / "runtime"))
from assemble import Assembler

assembler = Assembler(folder / "domain/native-diagnostic-v1.json")
original = assembler.assemble({"genes": [2, 2, 2]})
gamma = original["gamma96"]
gamma_eigenvalues = np.linalg.eigvalsh(gamma)
radius = 1.0
quarter = original["alpha"]
free = np.sort(np.concatenate([np.abs((mode + quarter) * gamma_eigenvalues / radius)
                               for mode in original["modes"]]))
actual = np.sort(np.abs(np.linalg.eigvalsh(original["h_N"])))
potential = {order: -1j * gamma @ phi for order, phi in original["phi_fourier"].items()}
potential_bound = sum(float(np.linalg.norm(matrix, ord=2)) for matrix in potential.values())
max_shift = np.max(np.abs(actual - free))
assert max_shift <= potential_bound + 1e-10
print("domain_sha256", original["native_inputs_sha256"])
print("spatial_gamma_eigen_counts",
      {int(sign): int(np.sum(np.isclose(gamma_eigenvalues, sign))) for sign in (-1, 1)})
print("full96 Fourier orders", sorted(potential), "continuum_potential_bound", potential_bound,
      "native384_max_singular_shift", max_shift)


def count_free(threshold):
    return len(gamma) * max(0, math.floor(radius * threshold - quarter)
                           - math.ceil(-radius * threshold - quarter) + 1)


for threshold in (20, 80, 320, 1280):
    lower = count_free(threshold - potential_bound)
    upper = count_free(threshold + potential_bound)
    asymptotic_count = 192 * threshold
    print("Lambda", threshold, "native_continuum_spatial_count_interval", (lower, upper),
          "ratio_to_192_Lambda", (lower / asymptotic_count, upper / asymptotic_count))

# A nonconstant arclength pullback fixing the four original quarter marks.
angle, normal, epsilon = sp.symbols("angle normal epsilon", real=True)
symbolic_radius = sp.symbols("radius", positive=True)
coordinate_map = angle + epsilon * sp.sin(4 * angle) / 4
profile = symbolic_radius * sp.diff(coordinate_map, angle)
metric = sp.diag(1, profile**2)
inverse = metric.inv()
coordinates = (normal, angle)
connection = {}
for i in range(2):
    for j in range(2):
        for k in range(2):
            connection[i, j, k] = sp.simplify(sum(
                inverse[i, m] * (sp.diff(metric[m, j], coordinates[k])
                                 + sp.diff(metric[m, k], coordinates[j])
                                 - sp.diff(metric[j, k], coordinates[m])) / 2
                for m in range(2)
            ))
for j in range(4):
    mark = sp.pi * j / 2
    assert sp.simplify(coordinate_map.subs(angle, mark) - mark) == 0
for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(2):
                curvature = (
                    sp.diff(connection[i, l, j], coordinates[k])
                    - sp.diff(connection[i, k, j], coordinates[l])
                    + sum(connection[i, k, m] * connection[m, l, j]
                          - connection[i, l, m] * connection[m, k, j]
                          for m in range(2))
                )
                assert sp.simplify(curvature) == 0
print("four quarter marks fixed, all 16 two-dimensional curvature tensor entries vanish exactly")
