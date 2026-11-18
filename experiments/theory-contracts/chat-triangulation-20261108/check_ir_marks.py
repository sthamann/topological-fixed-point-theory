"""Cross-check the original R032/R028 constant-H exterior chart."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

reference_path = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("r032_reference", reference_path)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)

beta = float(ref.time)
diagonal, offdiagonal = ref.grouped_higgs_coefficients()
diagonal = [(float(x), float(c)) for x, c in diagonal]
offdiagonal = [(float(x), float(y), float(c)) for x, y, c in offdiagonal if c]


def raw_prime(x):
    energy = np.sqrt(x)
    return beta * np.tanh(beta * energy / 2) / (4 * energy)


def raw_second(x):
    energy = np.sqrt(x)
    tanh = np.tanh(beta * energy / 2)
    return beta**2 * (1 - tanh**2) / (16 * x) - beta * tanh / (8 * energy**3)


def sea_prime(x):
    return -beta / (4 * np.sqrt(x))


def sea_second(x):
    return beta / (8 * x**1.5)


def cooled_prime(x):
    energy = np.sqrt(x)
    exp = np.exp(-beta * energy)
    return -beta * exp / (2 * energy * (1 + exp))


def cooled_second(x):
    energy = np.sqrt(x)
    exp = np.exp(-beta * energy)
    return beta*exp/(4*energy**3*(1+exp)) + beta**2*exp/(4*x*(1+exp)**2)


def mark_coefficient(momentum, first, second):
    q = momentum**2
    result = sum(weight * first(q + eig) for eig, weight in diagonal)
    for eig_i, eig_j, weight in offdiagonal:
        if eig_i == eig_j:
            term = second(q + eig_i)
        else:
            term = (first(q + eig_i) - first(q + eig_j)) / (eig_i - eig_j)
        result += weight * term / 2
    return result


def native_raw_y(radius, amplitude):
    gamma = ref.grading
    y = -1j * gamma @ ref.higgs_matrix * amplitude
    results = []
    for p in (2 * np.arange(4 * radius) + 1) / (4 * radius):
        h = p * gamma - 1j * gamma @ ref.mass_matrix + y
        energies, vec = np.linalg.eigh(h)
        occ = 1 / (1 + np.exp(beta * energies))
        results.append(float(np.dot(occ, np.diag(vec.conj().T @ y @ vec).real)))
    return sum(results)


def zero_kernel_by_sector():
    values, vectors = np.linalg.eigh(ref.mass_matrix @ ref.mass_matrix)
    basis = vectors[:, np.abs(values) < 1e-10]
    projector = basis @ basis.conj().T
    restricted = projector @ ref.higgs_matrix @ projector
    fields = [field[0] for field in ref.native.P] * 2
    return {name: float(np.linalg.norm(restricted[np.ix_(np.equal(fields, name), np.equal(fields, name))])**2)
            for name in ("u", "d", "e", "nu")}


def color_integral(radius, grid):
    # One actual 36-component quark particle chart, one quarter twist.
    # The entire 96-component exterior chart adds its conjugate and squares this weight.
    p = (2 * np.arange(12 * radius) + 1) / (4 * radius)
    e = np.exp(-beta * p)
    s = e / (1 + e)**2
    s_u = -beta * s * np.tanh(beta * p / 2) / (2 * p)
    yukawa_sum = float(sum(np.linalg.norm(y)**2 for y in ref.yukawa[:2]))
    angle_index = np.arange(grid)
    d = 4 * np.sin(np.pi * angle_index / grid)**2
    log_weight = np.zeros(grid)
    log_derivative = np.zeros(grid)
    for value, slope in zip(s, s_u):
        log_weight += 6 * np.log1p(-value * d)
        log_derivative -= yukawa_sum * slope * d / (1 - value * d)
    i, j = angle_index[:, None], angle_index[None, :]
    k = (-i - j) % grid
    measure = d[(i-j) % grid] * d[(i-k) % grid] * d[(j-k) % grid] / 6
    assert abs(measure.mean() - 1) < 1e-14
    term = measure * np.exp(log_weight[i] + log_weight[j] + log_weight[k])
    derivative = log_derivative[i] + log_derivative[j] + log_derivative[k]
    return float(term.mean()), float((term * derivative).sum() / term.sum()), float(-4*s_u.sum()*yukawa_sum/(6*s.sum()))


def color_center(radius):
    p = (2 * np.arange(12*radius) + 1) / (4*radius)
    e = np.exp(-beta*p)
    s = e/(1+e)**2
    sp_u = -beta*s*np.tanh(beta*p/2)/(2*p)
    ys = sum(np.linalg.norm(y)**2 for y in ref.yukawa[:2])
    log_nontrivial = 18*np.log1p(-3*s).sum()
    derivative = -9*ys*np.sum(sp_u/(1-3*s))
    weight = np.exp(log_nontrivial)
    return [radius, float(log_nontrivial), float((1+2*weight)/3),
            float(2*weight*derivative/(1+2*weight))]


def exact_slice_moments():
    # Polynomial Taylor coefficients, with no physical density assigned.
    eps, laplace = sp.symbols("epsilon s", positive=True)
    volume = sp.pi**2 * eps**2 / 2
    angular = 2*sp.pi**2
    radial_linear = sp.integrate(angular/4 * sp.Symbol("r")**5, (sp.Symbol("r"), 0, sp.sqrt(eps)))
    gaussian_linear = (sp.sqrt(sp.pi/laplace))**3 * sp.sqrt(sp.pi)/(2*laplace**sp.Rational(3, 2))
    assert sp.simplify(radial_linear/volume - eps/6) == 0
    assert sp.simplify(gaussian_linear*laplace**3 - sp.pi**2/2) == 0
    assert sp.simplify(sp.integrate(sp.pi**2*sp.Symbol("v")**2/4*sp.exp(-laplace*sp.Symbol("v")), (sp.Symbol("v"), 0, sp.oo)) - gaussian_linear) == 0
    return dict(ball_volume=str(volume), ball_x1_squared=str(radial_linear),
                laplace_x1_squared=str(gaussian_linear), leading_ratio_factor=str(radial_linear/volume))


def check_color_yukawa_cross():
    fields = ref.native.P
    indices = [i for i, field in enumerate(fields) if field[0] in ("u", "d")]
    color_signs = np.array([-1 if fields[i][2] in (0, 1) else 1 for i in indices])
    group = np.diag(color_signs)
    gamma = ref.grading[np.ix_(indices, indices)]
    yh = -1j * gamma @ ref.higgs_matrix[np.ix_(indices, indices)]
    squares = np.concatenate([np.linalg.svd(y, compute_uv=False)**2 for y in ref.yukawa[:2]])
    identity = np.eye(len(indices))
    p = 0.25

    def native_character(scale, h):
        e, basis = np.linalg.eigh(p*gamma + scale*h*yh)
        cov = (basis * (1/(1+np.exp(beta*e)))) @ basis.conj().T
        return np.linalg.det(identity - cov + cov @ group).real

    rows = []
    for h in (0.25, 0.5, 1.0, 2.0):
        e = np.sqrt(p*p + h*h*squares)
        exp = np.exp(-beta*e)
        s = exp/(1+exp)**2
        s_coupling = -beta * s * np.tanh(beta*e/2) * h*h*squares/e
        xi = float(np.prod((1-4*s)**2))
        slope = xi * np.sum(-8*s_coupling/(1-4*s))
        step = 0.00002
        direct = (native_character(1+step,h)-native_character(1-step,h))/(2*step)
        assert abs(native_character(1,h)-xi) < 2e-14
        assert abs(slope-direct) < 2e-8 and slope > 0
        rows.append([h, xi, -float(slope)/beta, -float(direct)/beta])
    return rows


def main():
    by_sector = zero_kernel_by_sector()
    coefficient = sum(by_sector.values())
    assert by_sector["nu"] < 1e-20 and min(by_sector[k] for k in ("u", "d", "e")) > 0
    rows = []
    # Keep the diagnostic momentum ceiling fixed when the circle grows.
    for r in (1, 8, 32, 128, 512, 1024):
        p = (2 * np.arange(4*r) + 1) / (4*r)
        raw = mark_coefficient(p, raw_prime, raw_second)
        sea = mark_coefficient(p, sea_prime, sea_second)
        cooled = mark_coefficient(p, cooled_prime, cooled_second)
        assert np.max(np.abs(cooled - raw - sea)) < 2e-7
        rows.append([r, float(raw.sum()/r), float(sea.sum()/r), float(cooled.sum()/r)])
    slope = (rows[-1][3] - rows[-2][3]) / np.log(2)
    assert abs(slope + beta * coefficient/2) < 2e-4
    assert abs(rows[-1][1] - rows[-2][1]) < 2e-4
    p_unit = (2 * np.arange(4) + 1) / 4
    expected = -2 * mark_coefficient(p_unit, raw_prime, raw_second).sum() / beta
    samples = [native_raw_y(1, h) / h**2 for h in (0.0005, 0.00025)]
    extrapolated = (4 * samples[1] - samples[0]) / 3
    assert abs(extrapolated - expected) < 2e-6
    color = []
    for r in (8, 16, 32, 64, 128, 256):
        low = color_integral(r, 256)
        high = color_integral(r, 512)
        assert abs(low[0]/high[0]-1) < 2e-7 and abs(low[1]/high[1]-1) < 2e-7
        color.append([r, high[0], high[1], high[2]])
    power = np.log2(color[-1][1]/color[-2][1])
    assert abs(power + 4) < 0.08
    assert abs(color[-1][2]/color[-1][3] - 1) < 0.02
    centers = [color_center(r) for r in (8, 16, 32)]
    assert centers[-1][1] < 2*centers[0][1] and abs(centers[-1][3]) < 1e-60
    print(json.dumps(dict(scope="R032 constant external H. Raw mark and sea-corrected log weight differentiated separately. SU3 tested in original quark particle chart with global Haar projection.",
        zero_kernel_weight=by_sector,
        rows_radius_raw_sea_cooled_divided_by_radius=rows,
        analytic_cooled_log_slope=-beta*coefficient/2, observed_cooled_log_slope=slope,
        native_raw_y_per_u_at_radius_1=extrapolated, analytic_raw_y_per_u_at_radius_1=expected,
        color_rows_radius_J_dlogJ_saddle_coefficient=color, color_log2_power=power,
        color_center_rows=centers, exact_slice_moments=exact_slice_moments(),
        conditional_color_yukawa_cross=check_color_yukawa_cross()), indent=2))


if __name__ == "__main__":
    main()
