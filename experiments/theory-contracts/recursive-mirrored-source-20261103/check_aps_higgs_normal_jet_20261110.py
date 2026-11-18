"""Check the flat spin4 product APS heat coefficient with the existing v252 input."""
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import erfcx

repo = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(repo / "verification"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "chat-triangulation-20261108"))
from check_quarter_reference import higgs_matrix, mass_matrix

Y = higgs_matrix
M = mass_matrix
boundary_value = 0.55
phi = M + boundary_value * Y
masses, vectors = np.linalg.eigh(phi)
mass_square_derivative = 2 * masses * np.diag(vectors.conj().T @ Y @ vectors).real
A = np.trace(Y @ Y).real
B = np.trace(Y @ Y @ Y @ Y).real
C = np.trace(Y @ Y @ M @ M).real
assert np.allclose(mass_square_derivative.sum(), 2 * boundary_value * A)

# The full tangential action has two plus and two minus branches per internal eigenline.
I2 = np.eye(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1, -1])
gamma = [np.kron(sx, sx), np.kron(sx, sy), np.kron(sx, sz), np.kron(sy, I2)]
gamma5 = gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]
N = np.kron(-1j * gamma[3], np.eye(96))
momentum = np.array([0.25, 0.375, -0.125])
tangent = sum(np.kron(p * g, np.eye(96)) for p, g in zip(momentum, gamma))
tangent += np.kron(gamma5, phi)
tangent = -N @ tangent
slope = -N @ np.kron(gamma5, Y)
assert np.allclose(tangent, tangent.conj().T)
assert np.allclose(N @ tangent + tangent @ N, 0)
assert np.allclose(tangent @ slope + slope @ tangent, np.kron(np.eye(4), phi @ Y + Y @ phi))
energies = np.sqrt(momentum @ momentum + masses**2)
expected_spectrum = np.sort(np.repeat(np.r_[energies, -energies], 2))
assert np.allclose(np.linalg.eigvalsh(tangent), expected_spectrum)

def heat_pair(q):
    """Boundary correction of the D/R pair for the linear normal jet."""
    def integrand(x):
        image = np.exp(-x * x) / (2 * np.sqrt(np.pi))
        robin = q * np.exp(-x * x) * erfcx(x + q)
        return image / q - (x + 0.5 / q) * robin

    return quad(integrand, 0, np.inf, epsabs=1e-11, epsrel=1e-11)[0]

def reduced_pair(q):
    return 0.25 / q - 0.5 / np.sqrt(np.pi)

for q in [0.03125, 0.25, 0.75, 1.0, 3.0, 8.0]:
    assert np.isclose(heat_pair(q), reduced_pair(q), atol=2e-11, rtol=0)
radial = quad(lambda q: q * q * np.exp(-q * q) * heat_pair(q), 0, 9, epsabs=1e-10)[0]
assert abs(radial) < 2e-10

def fixed_heat_pair(q):
    """Same Robin kernel with a frozen free boundary projector."""
    def integrand(x):
        return -q * x * np.exp(-x * x) * erfcx(x + q)

    return quad(integrand, 0, np.inf, epsabs=1e-11, epsrel=1e-11)[0]

fixed_radial = quad(lambda q: q * q * np.exp(-q * q) * fixed_heat_pair(q),
                    0, 9, epsabs=1e-10)[0]
assert np.isclose(fixed_radial, -1 / 16, atol=2e-10, rtol=0)
free_tangent = -N @ sum(np.kron(p * g, np.eye(96)) for p, g in zip(momentum, gamma))
mass_term = tangent - free_tangent
assert np.allclose(free_tangent @ mass_term + mass_term @ free_tangent, 0)
fixed_projector = (np.eye(384) - free_tangent / np.linalg.norm(momentum)) / 2
assert np.allclose(fixed_projector @ mass_term @ fixed_projector, 0)

# At finite t use all unchanged internal masses. This is the rescaled local
# jet coefficient; higher normal Taylor coefficients are excluded analytically.
def scaled_heat_jet(t):
    def integrand(p):
        q = np.sqrt(p * p + t * masses**2)
        return p * p * np.sum(mass_square_derivative * np.exp(-q**2) * reduced_pair(q))

    return -quad(integrand, 0, 10, epsabs=2e-10, limit=160)[0] / np.pi**2

ts = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
values = [scaled_heat_jet(t) for t in ts]
assert all(abs(y) < abs(x) for x, y in zip(values, values[1:]))
assert abs(values[-1]) < 1e-3 * abs(values[0])
a, b = A / (8 * np.pi**2), B / (8 * np.pi**2)
print(f"v252 A={A:.12g} B={B:.12g} C={C:.12g}")
print(f"spin4 x 96 spectrum: 384 eigenvalues; boundary background H={boundary_value}")
print(f"quartic a={a:.12g} b={b:.12g} sharp cubic coefficient={2*np.sqrt(a*b)/3:.12g}")
print(f"moving APS radial jet by direct quadrature: {radial:.3e}")
print(f"fixed free APS radial jet by direct quadrature: {fixed_radial:.12g}")
print(f"fixed free APS net cubic coefficient: {-np.sqrt(a*b)/3:.12g}")
for t, v in zip(ts, values):
    print(f"rescaled jet t={t:.0e}: {v:+.10e}")
# Verify exact subcriticality over the dyadic input field by reducing the
# Gaussian-rational matrix modulo a prime for which -1 is a square.
prime = 65537
i_mod = 256
assert i_mod**2 % prime == prime - 1

def mod_dyadic(value):
    num, den = float(value).as_integer_ratio()
    return num * pow(den, -1, prime) % prime

Y_mod = np.array([
    [(mod_dyadic(z.real) + i_mod * mod_dyadic(z.imag)) % prime for z in row]
    for row in Y
], dtype=np.int64)
Y2_mod = Y_mod @ Y_mod % prime
A_mod = int(np.trace(Y2_mod)) % prime
B_mod = int(np.trace(Y2_mod @ Y2_mod % prime)) % prime
test = (4 * A_mod * Y2_mod - B_mod * np.eye(96, dtype=np.int64)) % prime
det_mod = 1
for j in range(96):
    pivot = next(row for row in range(j, 96) if test[row, j] != 0)
    if pivot != j:
        test[[j, pivot]] = test[[pivot, j]]
        det_mod = -det_mod
    diagonal = int(test[j, j])
    det_mod = det_mod * diagonal % prime
    factors = test[j+1:, j] * pow(diagonal, -1, prime) % prime
    test[j+1:, j:] = (test[j+1:, j:] - factors[:, None] * test[j, j:]) % prime
assert det_mod != 0
k = np.sqrt(B / A)
distance = np.min(np.abs(np.abs(np.linalg.eigvalsh(Y)) / k - 0.5))
print(f"critical-ratio check: det(4 A Y² - B I) mod {prime} = {det_mod}")
print(f"numeric distance from critical ratio 1/2: {distance:.12g}")
# Compare the full free and moving projectors, then sum over the same
# three-dimensional flat product with the original quarter shift in one axis.
I384 = np.eye(384)
moving_projector = (I384 - tangent @ np.kron(np.eye(4), (vectors / energies) @ vectors.conj().T)) / 2
assert np.allclose(moving_projector @ moving_projector, moving_projector)
expected_hs = 2 * np.sum(1 - np.linalg.norm(momentum) / energies)
assert np.isclose(np.linalg.norm(moving_projector - fixed_projector, "fro")**2, expected_hs)

ratios = []
for cutoff in [8, 16, 32, 64]:
    x = np.arange(-cutoff, cutoff + 1, dtype=np.int64)
    xx = (4 * x + 1)**2
    yy = (4 * x)**2
    p_squared = (xx[:, None, None] + yy[None, :, None] + yy[None, None, :]).ravel()
    shells, counts = np.unique(p_squared[p_squared <= (4 * cutoff)**2], return_counts=True)
    moduli = np.sqrt(shells[:, None] / 16)
    roots = np.sqrt(moduli**2 + masses[None, :]**2)
    hs = np.sum(counts * np.sum(2 * masses[None, :]**2 / (roots * (roots + moduli)), axis=1))
    limit_scale = 4 * np.pi * cutoff * np.sum(masses**2)
    ratios.append(hs / limit_scale)
    print(f"cutoff {cutoff:2d}: projector sum / leading linear term = {ratios[-1]:.9f}")
assert all(0 < x < y for x, y in zip(ratios, ratios[1:]))
assert abs(ratios[-1] - 1) < 0.06

def cap_log(q):
    return 2 * np.log1p((q - 1) / 2)

def scaled_cap_density(p, *, minimal):
    x = 1 / p**2
    e = np.sqrt(1 + x)
    deficit = x / (e * (1 + e))
    if x < 1e-4:
        first = 5 / 16 if minimal else -1 / 16
        second = -11 / 48 if minimal else 1 / 12
        return p * p * x * x * (first + x * second)
    return p * p * (cap_log(1 - deficit) + (x / 2 if minimal else deficit))

i_min = np.pi / 3 - 4 / 9
i_det2 = np.pi / 3 - 10 / 9
for minimal, exact in [(True, i_min), (False, i_det2)]:
    direct = quad(lambda p: scaled_cap_density(p, minimal=minimal), 0, np.inf,
                  epsabs=2e-10, limit=160)[0]
    assert np.isclose(direct, exact, rtol=0, atol=2e-8)
    print(f"cap minimal={minimal}: quadrature {direct:+.10f}; exact {exact:+.10f}")
third = np.sum(np.abs(np.linalg.eigvalsh(Y))**3)
assert third <= np.sqrt(A * B)
c_min, c_det2 = np.array([i_min, i_det2]) * third / (2 * np.pi**2)
assert 0 < -c_det2 < np.sqrt(a * b) / 3
print(f"v252 cubic cap log per unit area: minimal={c_min:+.10f}, det2={c_det2:+.10f}")
print(f"conditional fixed-residue-plus-det2 log coefficient: {np.sqrt(a*b)/3+c_det2:+.10f}")
def scaled_cap_product_trace(p):
    e = np.sqrt(1 + p*p)
    return p*p / (2 * e**2 * (e + p)**2)

i_seam = np.pi / 4 - 2 / 3
direct_seam = quad(scaled_cap_product_trace, 0, np.inf, epsabs=2e-10)[0]
assert np.isclose(direct_seam, i_seam, atol=2e-10)
print(f"cap det2 product correction, cubic coefficient per area: {i_seam * third / (2*np.pi**2):+.10f}")
print("28 direct checks passed, including exact modular nonvanishing")
