"""Transfer, marked resolvent and determinant on the original Pin strip."""

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq


def h(w):
    r = (abs(w) - .5) / .25
    return np.exp(1 - 1 / (1 - r*r)) if abs(r) < 1 else 0.


def density(w, eps):
    return np.exp(eps*h(w))


def flow(rhs, a, b, initial):
    sol = solve_ivp(rhs, (a, b), initial, method="DOP853",
                    rtol=2e-11, atol=2e-12, max_step=.06)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1]


def transfer(z, eps, k, cuts, marked=False):
    n = 4 if marked else 2

    def rhs(w, flat):
        rho = density(w, eps)
        a = np.zeros((n, n), dtype=complex)
        a[0, 1], a[1, 0] = -(z*rho+k), z*rho-k
        if marked:
            # Solve (D_eps-z)u=phi_0 and retain the output integral.
            a[1, 2] = a[3, 0] = np.sqrt(rho)
        return (a @ flat.reshape(n, n)).ravel()

    total = np.eye(n, dtype=complex)
    for a, b in zip(cuts, cuts[1:]):
        local = flow(rhs, a, b, np.eye(n, dtype=complex).ravel())
        total = local.reshape(n, n) @ total
    return total


def marked_green(m):
    return m[3, 2] - m[3, 0]*m[1, 2]/m[1, 0]


def characteristic_trace(eps, k):
    # M_j is the coefficient of z^j in the Cauchy transfer.
    a0 = np.array([[0., -k], [-k, 0.]])
    j = np.array([[0., -1.], [1., 0.]])

    def rhs(w, flat):
        m0, m1, m2 = flat.reshape(3, 2, 2)
        a1 = density(w, eps)*j
        return np.array([a0@m0, a0@m1+a1@m0, a0@m2+a1@m1]).ravel()

    coeffs = flow(rhs, 0., 1., np.array([np.eye(2), np.zeros((2, 2)),
                                        np.zeros((2, 2))]).ravel()).reshape(3, 2, 2)
    f0, f1, f2 = coeffs[:, 1, 0]
    return (f1/f0)**2 - 2*f2/f0


def independent_green_trace(eps, k):
    def flat(w):
        return np.array([[np.cosh(k*w), -np.sinh(k*w)],
                         [-np.sinh(k*w), np.cosh(k*w)]])

    norm_cut = np.zeros((2, 2))
    norm_cut[0, :] = flat(1.)[1, :]/flat(1.)[1, 0]
    j = np.array([[0., -1.], [1., 0.]])

    def inner(w):
        def integrand(v):
            kernel = flat(w) @ ((np.eye(2) if v <= w else 0.)-norm_cut) @ flat(-v) @ j
            return density(v, eps)*np.sum(kernel**2)
        return density(w, eps)*quad(integrand, 0, 1, points=[.25, .75, w],
                                     epsabs=5e-11, epsrel=5e-11)[0]

    return quad(inner, 0, 1, points=[.25, .75], epsabs=5e-10, epsrel=5e-10)[0]


def eigenvalue(index, eps, k, optical_length):
    def mismatch(z):
        theta = flow(lambda w, t: z*density(w, eps)-k*np.cos(2*t),
                     0., 1., np.zeros(1))[0]
        return theta-index*np.pi

    center = index*np.pi
    return brentq(mismatch, (center-abs(k)-.01)/optical_length,
                   (center+abs(k)+.01)/optical_length, xtol=1e-11)


def marked_eigenweight(z, eps, k):
    def rhs(w, v):
        yp, ym, mass, source = v
        rho = density(w, eps)
        return [-(z*rho+k)*ym, (z*rho-k)*yp,
                rho*(yp*yp+ym*ym), np.sqrt(rho)*yp]

    yp, ym, mass, source = flow(rhs, 0., 1., np.array([1., 0., 0., 0.]))
    assert abs(ym) < 2e-8
    return source**2/mass


def occupation_gain(k, hm):
    e = np.hypot(k, np.pi*np.arange(1, len(hm)+1))
    return float(np.sum(hm**2/16*((e-k)/(e+k))**2*(1-k/e)))


def main(root):
    path = root / ("experiments/theory-contracts/original-admissible-source-closure-20261025/"
                   "geometric_source_return/source_geometric_machine.py")
    spec = importlib.util.spec_from_file_location("original_pin", path)
    pin = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pin)
    k = np.pi/2
    eps = .1
    unchanged = float(np.max(abs(transfer(0., eps, k, [-1., 1.])
                                 - transfer(0., 0., k, [-1., 1.]))))
    assert unchanged == 0.
    _, d, _, _ = pin.dirac_block(128, eps, k)
    eig, vectors = np.linalg.eigh(d)
    z = 2j
    whole = transfer(z, eps, k, [0., 1.], marked=True)
    stitched = transfer(z, eps, k, [0., .37, .71, 1.], marked=True)
    sewing_error = float(np.max(abs(whole-stitched)))
    green_error = float(abs(marked_green(stitched)-np.sum(abs(vectors[128, :])**2/(eig-z))))
    assert sewing_error < 3e-9 and green_error < 3e-9

    f0 = transfer(0., 0., k, [0., 1.])[1, 0].real

    def log_factor(t, strength):
        f = transfer(1j*t, strength, k, [0., 1.])[1, 0]
        return 2*np.log(abs(f/f0))

    trace_ode = [characteristic_trace(strength, k) for strength in (0., eps)]
    trace_integral = [independent_green_trace(strength, k) for strength in (0., eps)]
    trace_error = float(np.max(abs(np.asarray(trace_ode)-trace_integral)))
    assert trace_error < 2e-8
    optical = quad(lambda w: density(w, eps), 0., 1., points=[.25, .75],
                   epsabs=1e-12)[0]
    factor = log_factor(2., eps)-log_factor(2., 0.)
    residuals = [float(log_factor(t, eps)-log_factor(t, 0.)-2*t*(optical-1))
                 for t in (8., 16., 32., 64.)]

    roots = np.array([eigenvalue(i, eps, k, optical) for i in range(-4, 5)])
    heat = float(np.sum(np.exp(-roots**2)))
    heat_direct = float(np.sum(np.exp(-eig**2)))
    heat_reference = float(np.exp(-k*k)*(1+2*np.sum(np.exp(-(np.pi*np.arange(1, 5))**2))))
    assert abs(heat-heat_direct) < 3e-9 and heat > heat_reference

    hm = np.array([2*quad(lambda w: h(w)*np.cos(np.pi*m*w), .25, .75,
                         epsabs=1e-13)[0] for m in range(1, 65)])
    coefficients = [occupation_gain(2*np.pi*(n+.25), hm)
                    + occupation_gain(2*np.pi*(n+.75), hm) for n in range(512)]
    current = float(np.dot((-1.)**np.arange(512), coefficients))
    assert current > 0 and np.all(np.diff(coefficients) < 0)
    direct_weights = []
    # Reconstruct occupied spectral mass from eigenfunction residues,
    # independently of the earlier interior matrix projection.
    for strength in (.005, -.005):
        length = quad(lambda w: density(w, strength), 0., 1., points=[.25, .75])[0]
        mass = sum(marked_eigenweight(eigenvalue(-m, strength, k, length), strength, k)
                   for m in range(1, 17))
        direct_weights.append(mass/strength**2)
    direct_gain = float(np.mean(direct_weights))
    gain_16 = occupation_gain(k, hm[:16])
    assert abs(direct_gain-gain_16) < 2e-8

    print(json.dumps({
        "static_two_boundary_transfer_max_change": unchanged,
        "spectral_cut_composition_error_at_2i": sewing_error,
        "marked_resolvent_vs_original_galerkin_error_at_2i": green_error,
        "relative_log_fredholm_factor_at_2i": factor,
        "trace_inverse_square_from_boundary": trace_ode,
        "trace_inverse_square_from_independent_green_integral": trace_integral,
        "trace_inverse_square_max_error": trace_error,
        "large_parameter_residuals_after_optical_length": residuals,
        "heat_at_1_flat": heat_reference,
        "heat_at_1_deformed": heat,
        "heat_at_1_direct_error": abs(heat-heat_direct),
        "quadratic_current_coefficient_from_boundary_eigenmodes": current,
        "lowest_positive_mode_gain_from_boundary_derivative": occupation_gain(k, hm),
        "lowest_positive_mode_gain_from_boundary_derivative_16": gain_16,
        "lowest_positive_mode_gain_from_16_residues_centered": direct_gain
    }, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
