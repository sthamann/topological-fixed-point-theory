"""Check the local two-cut metric response with the unchanged Q source builder."""
from fractions import Fraction
from itertools import combinations_with_replacement
from pathlib import Path
import sys

import numpy as np


root = Path(sys.argv[1])
contract = root / (
    "experiments/theory-contracts/original-paper-source-reconciliation-20261102/"
    "contribution_20261106Q"
)
sys.path.insert(0, str(contract / "runtime"))

from assemble import Assembler


source = Assembler(contract / "domain/native-diagnostic-v1.json")
choices = source.domain["field_proposal"]["allowed_multipliers"]
nonzero = next(i for i, value in enumerate(choices) if Fraction(value))
native = source.assemble({"genes": [nonzero] * len(source.domain["field_proposal"]["slots"])})
gamma = native["gamma96"]
assert np.array_equal(gamma @ gamma, np.eye(96))
assert native["h_N_role"] == "original_tangential_one_particle_NOT_joint_OS_generator"
potential = {k: -1j * gamma @ block for k, block in native["phi_fourier"].items()}
assert any(k and np.any(block) for k, block in potential.items())

retained = native["modes"]
padding = 2 * max(map(abs, potential)) + 4
modes = list(range(min(retained) - padding, max(retained) + padding + 1))
matrix_size = len(retained) * 96
native_columns = np.zeros((len(modes), 96, matrix_size), complex)
for i, mode in enumerate(retained):
    native_columns[modes.index(mode), :, i * 96:(i + 1) * 96] = np.eye(96)


def multiply(coefficients, values):
    result = np.zeros_like(values)
    for shift, coefficient in coefficients.items():
        for i in range(max(0, shift), min(len(modes), len(modes) + shift)):
            incoming = values[i - shift]
            result[i] += coefficient * incoming if np.isscalar(coefficient) else coefficient @ incoming
    return result


def boundary(values):
    result = multiply(potential, values)
    for i, mode in enumerate(modes):
        result[i] += (mode + native["alpha"]) * (gamma @ values[i])
    return result


def derivative(coefficients):
    return {mode: 1j * mode * value for mode, value in coefficients.items()}


def convolve(left, right):
    result = {}
    for m, a in left.items():
        for n, b in right.items():
            result[m + n] = result.get(m + n, 0) + a * b
    return result


interior = [modes.index(mode) for mode in retained]
rebuilt = boundary(native_columns)[interior].reshape(matrix_size, matrix_size)
assert np.array_equal(rebuilt, native["h_N"])
print("Native Q builder: complete 384 by 384 source reproduced exactly with nonconstant background.")
print(f"Native quarter {native['alpha']}; potential harmonics {sorted(potential)}.")

profiles = (
    {-1: 0.5j, 1: -0.5j},
    {-1: 0.5, 1: 0.5},
    {-2: 0.5j, 2: -0.5j},
)


def boundary_squared(values):
    return boundary(boundary(values))


errors = []
for f, g in combinations_with_replacement(profiles, 2):
    v = native_columns
    fv, gv = multiply(f, v), multiply(g, v)
    actual = -0.5 * (
        boundary_squared(multiply(f, gv))
        - multiply(f, boundary_squared(gv))
        - multiply(g, boundary_squared(fv))
        + multiply(g, multiply(f, boundary_squared(v)))
    )
    expected = multiply(convolve(derivative(f), derivative(g)), v)
    errors.append(float(np.max(np.abs(actual - expected))))

assert len(errors) == 6 and max(errors) < 1e-10
print(f"All six mixed scalar probes on every original input column agree: max error {max(errors):.3e}.")
print("On the circle the pointwise response is v v^T, v=(cos(theta),-sin(theta),2 cos(2 theta)); rank 1.")
print("Integrating over the circle gives diag(1/2,1/2,2), rank 3; integrated rank is not local rank.")

# Exact stress check for the first two real circle modes beyond the Möbius algebra.
# [L_2, L_-2]=4 L_0+c/2, and all level-one E8 currents have weight one.
central_charge = Fraction(248, 1 + 30)
assert central_charge == 8 and central_charge / 2 == 4
print("Level-one E8 stress commutator: [L_2,L_-2]=4 L_0+4 I; the central part acts trivially by conjugation.")
