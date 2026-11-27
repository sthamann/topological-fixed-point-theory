"""Two face-current modes of the original closed annulus CAR projector.

The state assumption is explicit: the *whole* geometric projector is the
occupation covariance. The continuous boost kernel is calculated separately.
"""

from itertools import product

import mpmath as mp
import sympy as sp


def current_spectra(h, alpha, m, coherence=1):
    """Wick sums, using both blocks and no subtraction for m>0."""

    def occupation(n):
        return 1 / (1 + mp.exp(h * (n + alpha)))

    def off_diagonal(n):
        return coherence / (2 * mp.cosh(h * (n + alpha) / 2))

    low = high = cross = mp.mpf(0)
    for n in range(-65, 66):
        # Evaluate both occupation and hole probability without cancellation.
        p = occupation(n)
        hole = 1 / (1 + mp.exp(-h * (n + alpha)))
        pm = occupation(n + m)
        hm = 1 / (1 + mp.exp(-h * (n + m + alpha)))
        high += p * hm
        low += hole * pm
        cross -= off_diagonal(n) * off_diagonal(n + m)
    return low, high, cross


def telescope():
    x, y = sp.symbols("x y", positive=True)  # exp[h(n+α)/2], exp(hm/2)
    p = 1 / (1 + x**2)
    pm = 1 / (1 + (x * y) ** 2)
    cross = x / (1 + x**2) * x * y / (1 + (x * y) ** 2)
    assert sp.simplify((y - 1 / y) * cross - (p - pm)) == 0
    assert sp.simplify(p * (1 - pm) - y**2 * pm * (1 - p)) == 0


def oscillator_pair():
    for h, alpha in [
        (2 * mp.pi, mp.mpf("0.5")),
        (4 * mp.pi, mp.mpf("0.25")),
        (4 * mp.pi, mp.mpf("0.75")),
    ]:
        for m in (1, 2, 3):
            low, high, cross = current_spectra(h, alpha, m)
            scale = h * m
            same = (low + high) / 2
            assert mp.almosteq(high - low, m)
            assert mp.almosteq(same, m * mp.coth(scale / 2) / 2)
            assert mp.almosteq(cross, -m / (2 * mp.sinh(scale / 2)))
            assert mp.almosteq(cross**2, low * high)
            assert mp.almosteq(2 * mp.acosh(same / abs(cross)), scale)
            # Opposite chiralities: <a_L a_R>=cross/m, occupancies are low/m on both sides.
            t = mp.exp(-scale / 2)
            assert abs((low + 2 * t * cross + t * t * high) / m) < mp.mpf("1e-49")
            # Product control on the identical spectra: saturation and matching disappear.
            product_low, product_high, product_cross = current_spectra(
                h, alpha, m, coherence=0
            )
            assert (low, high) == (product_low, product_high) and product_cross == 0
            assert low * high > product_cross**2
            if m == 1:
                words = [
                    mp.exp(-(i * i + j * j) * same / 2 - i * j * cross)
                    for i, j in product((0, 1), repeat=2)
                ]
                assert mp.almosteq(
                    mp.log(words[-1] * words[0] / (words[1] * words[2])), -cross
                )
                print(
                    "annulus h/π =",
                    float(h / mp.pi),
                    "offset =",
                    float(alpha),
                    "m=1 sym =",
                    float(same),
                    "cross =",
                    float(cross),
                    "ratio =",
                    float(abs(cross) / same),
                    "Weyl words",
                    [float(w) for w in words],
                )
    # Charge-preserving standing quadrature (J_m+J_-m)/sqrt(2) on each face.
    # Joint Weyl formula uses the regular thermal current oscillator, not unsmeared CAR densities.


def boost_spectrum():
    # D.19 transformed to log-affine coordinates of opposite half lines.
    # ∫ exp(-ikx)/(4cosh²(x/2)) dx = Γ(1+ik)Γ(1-ik) = πk/sinh(πk).
    eps = mp.mpf("0.4")
    for k in (mp.mpf("0.5"), 1, 2, 3):
        value = mp.quad(
            lambda x, k=k: mp.cos(k * x) / (2 * mp.cosh(x / 2) ** 2),
            [0, 1, 4, 10, 30, 80, 160],
        )
        assert abs(value - mp.pi * k / mp.sinh(mp.pi * k)) < mp.mpf("1e-33")
        # From the affine measure, the same-side ordered spectrum is
        # 2πk/(1-exp(-2πk)); its even part is πk coth(πk).
        local = mp.pi * k * mp.coth(mp.pi * k)
        assert mp.almosteq(value / local, mp.sech(mp.pi * k))
        # Independent transform on a contour below the same-side pole distribution.
        # Removing the imaginary-time translation restores the physical spectra.
        ordered = []
        for freq in (k, -k):

            def integrand(x, freq=freq):
                return 2 * mp.re(
                    -mp.exp(1j * freq * x) / (4 * mp.sinh((x - 1j * eps) / 2) ** 2)
                )

            transform = mp.quad(integrand, [0, 1, 4, 10, 30, 80, 160]) * mp.exp(
                eps * freq
            )
            ordered.append(transform)
        assert abs(sum(ordered) / 2 - local) < mp.mpf("1e-32")
        assert abs(ordered[0] / ordered[1] - mp.exp(2 * mp.pi * k)) < mp.mpf("1e-26")
    k = 1
    assert not mp.almosteq(mp.sech(mp.pi * k), mp.sech(2 * mp.pi * k))
    # Every current harmonic recurs after a 2π translation on the dimensionless circle.
    # The uncompactified boost cross kernel decays and is aperiodic.
    peak = mp.mpf(1) / 4
    after = peak / mp.cosh(mp.pi) ** 2
    print(
        "Boost response after rapidity 2π divided by response at zero =",
        float(after / peak),
    )
    assert 0 < after / peak < mp.mpf("0.008")


def real_boost_clock():
    # W1 gives a strict period on its finite charge spectrum.
    for charge in range(5):
        phase = mp.exp(1j * 4 * mp.pi * charge / 2)
        assert mp.almosteq(phase, 1)
    # Under a real null boost the null tangent rescales by exp(eta).
    # In the ideal plane regime Eq. (5.29) gives an observable return test.
    eta = 4 * mp.pi
    shear_return = mp.sech(eta / 2) ** 2
    assert 0 < shear_return < mp.mpf(".000015")
    assert mp.almosteq(shear_return, 4 * mp.exp(eta) / (1 + mp.exp(eta)) ** 2)
    print(
        "At real rapidity 4π the finite W1 clock returns, the narrow shear ratio is",
        float(shear_return),
    )


def localized_packet():
    # Normalized CCR weight |f|² is uniform on [1/2,3/2] in boost frequency.
    # This is a finite-band control. Compactly supported time pulses also have finite spread.
    lo, hi = mp.mpf("0.5"), mp.mpf("1.5")

    def mean(fn):
        return mp.quad(fn, [lo, hi]) / (hi - lo)

    n = mean(lambda k: 1 / mp.expm1(2 * mp.pi * k))
    c = mean(lambda k: -1 / (2 * mp.sinh(mp.pi * k)))
    deficit = n * (1 + n) - c**2
    exact = (
        mean(lambda k: mp.tanh(mp.pi * k / 2)) * mean(lambda k: mp.coth(mp.pi * k / 2))
        - 1
    ) / 4
    assert mp.almosteq(deficit, exact) and deficit > 0
    print(
        "Finite boost packet [1/2,3/2], ordered Cauchy-Schwarz deficit =",
        float(deficit),
    )


def native_quarter_charge():
    # On a single face R=i^(Q_even-Q_odd). A density transfer n+m -> n
    # is clock-odd for odd m, clock-even for even m.
    for n in range(-4, 5):
        for m in range(1, 5):
            charge_change = (-1) ** (n % 2) - (-1) ** ((n + m) % 2)
            assert 1j**charge_change == (-1) ** m
    full = current_spectra(2 * mp.pi, mp.mpf("0.5"), 2)
    sub = [
        current_spectra(4 * mp.pi, alpha, 1) for alpha in (mp.mpf(".25"), mp.mpf(".75"))
    ]
    assert all(mp.almosteq(a, b + c) for a, b, c in zip(full, *sub))

    # Pure phase-twisted projector with identical marginal, spectrum and global translations.
    # On one face multiply the odd NS modes by i. The two orientations of odd transfers
    # have opposite phases, so their positive weights no longer sum with equal phases.
    def off(n):
        return 1 / (2 * mp.cosh(mp.pi * (n + mp.mpf(".5"))))

    raw = -sum(off(n) * off(n + 1) for n in range(-65, 66))
    twisted = -sum(
        1j ** ((n + 1) % 2 - n % 2) * off(n) * off(n + 1) for n in range(-65, 66)
    )
    low, high, cross = current_spectra(2 * mp.pi, mp.mpf(".5"), 1)
    assert mp.almosteq(raw, cross) and low * high > abs(twisted) ** 2
    print(
        "Matching deficit for the pure parity-phase control =",
        float(low * high - abs(twisted) ** 2),
    )


if __name__ == "__main__":
    mp.mp.dps = 65
    telescope()
    oscillator_pair()
    boost_spectrum()
    real_boost_clock()
    localized_packet()
    native_quarter_charge()
