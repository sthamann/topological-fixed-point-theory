"""R032 exterior reference only, at the existing v252 diagnostic inputs."""
import json

import mpmath as mp
import numpy as np
import v252_full_finite_triple as native

mp.mp.dps = 60
rng = np.random.default_rng(252)
yukawa = tuple(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)) for _ in range(4))
raw = rng.normal(size=(3, 3))
majorana = raw + raw.T
zero = np.zeros((3, 3))
mass_matrix = native.full_D(zero, zero, zero, zero, majorana)
higgs_matrix = native.full_D(*yukawa, zero)
grading = native.grading()
masses = [mp.mpf(float(m)) for m in np.linalg.svd(majorana, compute_uv=False)]
time = 4 * mp.pi


def momentum(k, radius):
    return (2 * k + 1) / (4 * radius)


def g(energy):
    return mp.log1p(mp.exp(-time * energy))


def occupation(energy):
    return 1 / (1 + mp.exp(time * energy))


def log_ratio(radius, count):
    return 2 * mp.fsum(
        g(mp.sqrt(p * p + mass * mass)) - g(p)
        for p in (momentum(k, radius) for k in range(count))
        for mass in masses
    )


def metric_score(radius, count):
    return 2 * time * mp.fsum(
        p * p / energy * occupation(energy) - p * occupation(p)
        for p in (momentum(k, radius) for k in range(count))
        for energy in (mp.sqrt(p * p + mass * mass) for mass in masses)
    )


def mode_tail(radius, count):
    first = momentum(count, radius)
    spacing = 1 / (2 * radius)
    q = mp.exp(-time * spacing)
    start = mp.exp(-time * first)
    return 6 * start / (1 - q), 6 * time * start * (first / (1 - q) + spacing * q / (1 - q)**2)


def direct_log_trace(log_radius, higgs, count=12):
    radius = np.exp(log_radius)
    phi = mass_matrix + higgs * higgs_matrix
    total = 0.0
    for k in range(count):
        p = float(momentum(k, radius))
        eigenvalues = np.linalg.eigvalsh(p * grading - 1j * grading @ phi)
        total += np.logaddexp(0, -float(time) * abs(eigenvalues)).sum()
        total -= 96 * np.logaddexp(0, -float(time) * p)
    return total


def grouped_higgs_coefficients():
    spectrum, basis = np.linalg.eigh(mass_matrix @ mass_matrix)
    groups = []
    for val in np.unique(spectrum.round(10)):
        indices = np.flatnonzero(spectrum.round(10) == val)
        eigenvalue = mp.mpf(float(np.mean(spectrum[indices]))) if val > 1e-10 else mp.mpf(0)
        groups.append((eigenvalue, indices))
    quadratic = basis.conj().T @ (higgs_matrix @ higgs_matrix) @ basis
    linear = basis.conj().T @ (mass_matrix @ higgs_matrix + higgs_matrix @ mass_matrix) @ basis
    diag = [(lam, mp.mpf(float(np.trace(quadratic[np.ix_(i, i)]).real))) for lam, i in groups]
    offdiag = [
        (lam, mu, mp.mpf(float(np.sum(np.abs(linear[np.ix_(i, j)])**2))))
        for lam, i in groups for mu, j in groups
    ]
    return diag, offdiag


def g_deriv(x, order):
    return mp.diff(lambda y: g(mp.sqrt(y)), x, order)


def divided_difference(x, y, order):
    if x == y:
        return g_deriv(x, order + 1)
    return (g_deriv(x, order) - g_deriv(y, order)) / (x - y)


def higgs_score(radius, diagonal, offdiagonal, metric=False, count=24):
    order = 2 if metric else 1
    terms = []
    for k in range(count):
        q = momentum(k, radius)**2
        terms.append(((-2 * q) if metric else 1) * (
            mp.fsum(c * g_deriv(q + lam, order) for lam, c in diagonal)
            + mp.fsum(c * divided_difference(q + lam, q + mu, order) / 2
                      for lam, mu, c in offdiagonal if c)
        ))
    return mp.fsum(terms)


def finite_higgs(log_radius, amplitude):
    return (direct_log_trace(log_radius, amplitude) + direct_log_trace(log_radius, -amplitude)
            - 2 * direct_log_trace(log_radius, 0)) / (2 * amplitude**2)


def infrared_higgs(diagonal, offdiagonal):
    spectrum, basis = np.linalg.eigh(mass_matrix @ mass_matrix)
    zero_basis = basis[:, np.abs(spectrum) < 1e-10]
    projected = zero_basis.conj().T @ higgs_matrix @ zero_basis
    kernel_weight = np.linalg.norm(projected)**2
    by_groups = float(diagonal[0][1]) - sum(float(c / mu) for lam, mu, c in offdiagonal if lam == 0 and mu > 0)
    assert abs(by_groups - kernel_weight) < 2e-10
    clock = float(time)
    predicted = -clock * kernel_weight / 2

    def prime(x):
        e = np.sqrt(x)
        f = np.exp(-clock * e)
        return -clock * f / (2 * e * (1 + f))

    def second(x):
        e = np.sqrt(x)
        f = np.exp(-clock * e)
        return clock*f/(4*e**3*(1+f)) + clock**2*f/(4*e**2*(1+f)**2)

    rows = []
    for r in (32, 64, 128, 256, 512, 1024):
        p = (2 * np.arange(32*r) + 1) / (4*r)
        q = p*p
        response = sum(float(c)*prime(q+float(lam)) for lam,c in diagonal)
        for lam,mu,c in offdiagonal:
            if c:
                x, y = q+float(lam), q+float(mu)
                response += float(c)/2 * (second(x) if lam == mu else (prime(x)-prime(y))/float(lam-mu))
        rows.append((r, float(response.sum()/r)))
    observed = (rows[-1][1] - rows[-2][1]) / np.log(2)
    assert abs(observed / predicted - 1) < 1e-5
    return dict(kernel_yukawa_squared=kernel_weight, predicted_log_slope=predicted,
                last_two_radii_log_slope=observed, response_per_radius=rows)


def main():
    one = mp.mpf(1)
    assert np.linalg.matrix_rank(mass_matrix) == 6
    assert abs(direct_log_trace(0, 0) - float(log_ratio(one, 32))) < 2e-13
    assert abs(mp.diff(lambda x: log_ratio(mp.exp(x), 32), 0) - metric_score(one, 32)) < mp.mpf('1e-54')
    density = 4 * sum(mp.quad(lambda p: g(mp.sqrt(p*p + m*m)) - g(p), [0, 1, mp.inf]) for m in masses)
    density_bessel = 4 * sum(
        mp.fsum((-1)**(j+1)*m*mp.besselk(1, j*time*m)/j for j in range(1, 22))
        - mp.pi**2/(12*time) for m in masses
    )
    assert abs(density - density_bessel) < mp.mpf('1e-55')
    diagonal, offdiagonal = grouped_higgs_coefficients()
    infrared = infrared_higgs(diagonal, offdiagonal)
    d_u = higgs_score(one, diagonal, offdiagonal)
    mixed = higgs_score(one, diagonal, offdiagonal, metric=True)
    assert abs(mp.diff(lambda x: higgs_score(mp.exp(x), diagonal, offdiagonal), 0) - mixed) < mp.mpf('1e-52')
    fd = []
    for h in (0.0005, 0.00025, 0.000125):
        step = 0.0005
        fd.append([h, finite_higgs(0, h), (finite_higgs(step, h)-finite_higgs(-step, h))/(2*step)])
    extrapolated = [(4 * fd[1][i] - fd[0][i]) / 3 for i in (1, 2)]
    assert abs(extrapolated[0]-float(d_u)) < 0.000002
    assert abs(extrapolated[1]-float(mixed)) < 0.001
    p0 = momentum(24, one)
    exp_sum = mp.exp(-time*p0)/(1-mp.exp(-time/2))
    trace_quadratic = sum(c for _, c in diagonal)
    norm_linear_squared = sum(c for _, _, c in offdiagonal)
    tail_u = exp_sum * (trace_quadratic*time/(2*p0)
        + norm_linear_squared*(time/(4*p0**3) + time**2/(4*p0**2))/2)
    tail_mixed = exp_sum * (trace_quadratic*(time/(2*p0) + time**2/2)
        + norm_linear_squared*(3*time/(4*p0**3) + 3*time**2/(4*p0**2) + time**3/(4*p0))/2)
    sn = lambda x: mp.nstr(x, 32)
    rows = []
    for scale in (1, 2, 4, 8, 16, 32):
        r = mp.mpf(scale)
        count = 32*scale
        ell = log_ratio(r, count)
        score = metric_score(r, count)
        tails = mode_tail(r, count)
        assert ell < 0 and score < 0 and max(tails) < mp.mpf('1e-80')
        rows.append(dict(radius=scale, log_ratio=sn(ell), ratio=sn(mp.exp(ell)),
            metric_score=sn(score), log_ratio_per_radius=sn(ell/r), score_per_radius=sn(score/r),
            log_tail=sn(tails[0]), metric_tail=sn(tails[1])))
    print(json.dumps(dict(masses=list(map(sn, masses)), time=sn(time), cases=rows,
        density=sn(density), independent_density_delta=sn(abs(density-density_bessel)),
        higgs_log_derivative=sn(d_u), metric_higgs_log_derivative=sn(mixed),
        finite_differences=fd, amplitude_extrapolated=extrapolated,
        higgs_mode_tails=[sn(tail_u), sn(tail_mixed)],
        native_log_ratio=direct_log_trace(0, 0), equal_T_fermion_ratio=1,
        infrared=infrared), indent=2))


if __name__ == "__main__":
    main()
