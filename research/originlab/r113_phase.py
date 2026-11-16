"""Resolve charged R032 angular moments on one common conditional source."""
import json
from math import comb

from flint import acb, arb, ctx

from r113 import ball, conditional_word, pair_trace, squared_singular_values, WINDOWS
from source import NATIVE, source_refs


def affine_product(poly, c0, c1):
    return [c0 * poly[0], *[c0 * b + c1 * a for a, b in zip(poly, poly[1:])], c1 * poly[-1]]


def phase_polynomials(spectra, early, late, momenta, width):
    """Return the unmarked and late marked kernels in cos(relative phase)."""
    sectors = []
    sea = arb(0)
    pairs = []
    for p in momenta:
        for colour, spectrum in spectra:
            for squared in spectrum:
                positive, dp, energy = pair_trace(squared, early, late, p, width)
                negative, dn, _ = pair_trace(squared, -early, late, p, width)
                multiplicity = 6 if colour else 2
                sea -= width * multiplicity * energy
                pairs.extend([((positive + negative) / 2, (positive - negative) / 2,
                               (dp + dn) / 2, (dp - dn) / 2, colour)] * multiplicity)
    for twisted in (False, True):
        poly, jet = [arb(1)], [arb(0)]
        for base, slope, base_jet, slope_jet, colour in pairs:
            offset = -1 if colour and twisted else 2
            jet = [x + y for x, y in zip(affine_product(jet, base + offset, slope),
                                          affine_product(poly, base_jet, slope_jet))]
            poly = affine_product(poly, base + offset, slope)
        sectors.append((poly, jet))
    weight = sea.exp() / 3
    unmarked = [weight * (x + 2*y) for x, y in zip(sectors[0][0], sectors[1][0])]
    # The original Yukawa mark differentiates only the exterior transfer.
    marked = [-weight * (x + 2*y) / width for x, y in zip(sectors[0][1], sectors[1][1])]
    return unmarked, marked


def laurent(poly):
    return [sum(poly[d] * comb(d, (d-k)//2) / 2**d for d in range(k, len(poly), 2))
            for k in range(len(poly))]


def phase_transform(samples):
    size = len(samples)
    degree = size // 2
    return [sum(value * acb(arb(f"{2*k*j}/{size}")).exp_pi_i() for j, value in enumerate(samples))
            / size for k in range(-degree, degree+1)]


def mark_phase_scan(samples, unmarked_harmonics, marked_harmonics, phase_over_pi):
    """Read the late mark directly from unmarked scan data with the same scale."""
    size = 2 * len(unmarked_harmonics) - 1
    if len(samples) != size or len(marked_harmonics) != len(unmarked_harmonics):
        raise ValueError("require the whole uniform phase grid and matched harmonics")
    actions = []
    for unmarked, marked in zip(unmarked_harmonics, marked_harmonics):
        if unmarked.is_zero() and marked.is_zero():
            actions.append(arb(0))
        elif unmarked.contains(0):
            raise ValueError("uncertified nonzero unmarked harmonic")
        else:
            actions.append(marked / unmarked)
    degree = size // 2
    coefficients = phase_transform(samples)
    mark = sum(value * actions[abs(k)] * acb(-k*phase_over_pi).exp_pi_i()
               for k, value in enumerate(coefficients, -degree))
    return mark, coefficients, max(abs(value).upper() for value in actions)


def reduced_phase_scan(samples, unmarked, marked, phase_over_pi):
    """Reconstruct on a shorter grid and bound the missed and aliased harmonics."""
    size = len(samples)
    order = size // 2
    if not 0 < size < 2*len(unmarked)-1 or size % 2 != 1 or len(unmarked) != len(marked):
        raise ValueError("require a nonempty odd reduced grid and matched full harmonics")
    if any(not value > 0 for value in unmarked):
        raise ValueError("the reduced bound requires certified positive unmarked harmonics")
    response, coefficients, amplification = mark_phase_scan(samples, unmarked[:order+1], marked[:order+1], phase_over_pi)
    calibration = unmarked[0] - 2*sum(unmarked[size::size])
    mean = coefficients[order].real
    if not calibration > 0 or not mean >= 0:
        raise ValueError("scan mean cannot bound the total positive mass")
    error_per_mean = 2*(sum(abs(c) for c in marked[order+1:]) + amplification*sum(unmarked[order+1:])) / calibration
    return response, error_per_mean, error_per_mean * mean


def evaluate_poly(coefficients, x):
    value = arb(0)
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def replay(root):
    from verify import digest

    domain = root / NATIVE
    fixture = json.loads(domain.read_bytes())
    source_refs(root, fixture["source_originals"])
    with ctx.workprec(1536):
        spectra = [(name != "Ye", squared_singular_values(fixture["couplings"][name]))
                   for name in ("Yu", "Yd", "Ye")]
        early, late, width = map(ball, ("7/20", "11/20", "2/5"))
        momentum = ball("1/4")
        reports = []
        for momenta in ((momentum,), (momentum, ball("-3/4"))):
            plain, marked = phase_polynomials(spectra, early, late, momenta, width)
            harmonics = laurent(plain)
            degree = len(harmonics) - 1
            count = 2*degree + 1
            if not all(c > 0 for c in harmonics):
                raise ArithmeticError("native unmarked angular harmonics did not certify positive")
            # Independent direct evaluation at the two original real orientations.
            if len(momenta) == 1:
                for sign in (-1, 1):
                    k0, y = conditional_word(spectra, early, sign*late, momentum, width)
                    if not (k0.overlaps(evaluate_poly(plain, arb(sign))) and
                            (k0*y).overlaps(evaluate_poly(marked, arb(sign)))):
                        raise ArithmeticError("native real-orientation endpoint did not replay")
            # A declared unequal three-atom phase law is used solely to test the inverse.
            weights = [arb("1/5"), arb("3/10"), arb("1/2")]
            phases = [arb(0), arb("1/2"), arb(1)]  # multiples of pi
            def sample_scan(size):
                return [sum(w * evaluate_poly(plain, (p-arb(f"{2*j}/{size}")).cos_pi())
                            for w, p in zip(weights, phases)) for j in range(size)]

            samples = sample_scan(count)
            target = arb("1/7")
            marked_harmonics = laurent(marked)
            response, coefficients, amplification = mark_phase_scan(samples, harmonics, marked_harmonics, target)
            moments = [value / harmonics[abs(k)] for k, value in enumerate(coefficients, -degree)]
            expected = [sum(w * acb(k*p).exp_pi_i() for w, p in zip(weights, phases))
                        for k in range(-degree, degree+1)]
            if any(not a.overlaps(b) for a, b in zip(moments, expected)):
                raise ArithmeticError("inverse failed native harmonic round trip")
            if any(value.rad() > arb("1e-40") for value in moments):
                raise ArithmeticError("phase inversion did not achieve certified resolution")
            # No marked samples were used. Predict an unseen late response from the inverse.
            direct = sum(w * evaluate_poly(marked, (p-target).cos_pi()) for w, p in zip(weights, phases))
            if not response.real.overlaps(direct) or not response.imag.contains(0):
                raise ArithmeticError("unseen marked angular response did not replay")
            reduced_count = 25 if len(momenta) == 1 else 33
            reduced, relative_bound, absolute_bound = reduced_phase_scan(
                sample_scan(reduced_count), harmonics, marked_harmonics, target)
            if not abs(reduced - direct) < absolute_bound:
                raise ArithmeticError("reduced scan did not satisfy the certified alias bound")
            log_condition = (max(harmonics)/min(harmonics)).log() / arb(10).log()
            reports.append({"momentum_count":len(momenta), "degree":degree, "independent_phase_samples":count,
                            "log10_inverse_condition_scale":log_condition.str(10),
                            "late_mark_l2_amplification_upper":amplification.str(10),
                            "reduced_phase_samples":reduced_count,
                            "reduced_late_mark_error_over_scan_mean_upper":relative_bound.upper().str(10),
                            "recoverable_moments":[-degree,degree], "certified_test_roundtrip":True,
                            "withheld_late_mark_replayed":True})
    return {"status":"CONDITIONAL_R032_PHASE_INVERSION", "source_sha256":digest(domain),
            "source":str(WINDOWS), "reports":reports, "physical_regional_rotation_bound":False,
            "absolute_empty_reference_derived":False, "physical_source_identity_proved":False}
