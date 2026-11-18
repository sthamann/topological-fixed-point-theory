"""Two phase probes on the existing fixed-history exterior cut and native source gate."""
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import solve


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"


def run_existing(name):
    previous = base / name
    argv = sys.argv
    try:
        sys.argv = [str(previous), str(repo)]
        with contextlib.redirect_stdout(io.StringIO()) as output:
            data = runpy.run_path(str(previous))
        return data, json.loads(output.getvalue())
    finally:
        sys.argv = argv


old, _ = run_existing("check_normal_weak_curvature_bridge_20261109.py")
_, source_result = run_existing("check_original_qtheta_word_20261109.py")
gate = source_result["native_endpoint_call"]
src, native = old["source"], old["native"]
identity, swap, gamma = old["identity"], src["swap"], src["grading"]
first, cut = old["history"][0], old["history"][2]
ray = old["projectors"][2]
sigma_x = np.array([[0., 1.], [1., 0.]])
weak = ray @ sigma_x @ (np.eye(2) - ray) + (np.eye(2) - ray) @ sigma_x @ ray


def lift(g):
    lam = np.exp(.5j * np.angle(np.linalg.det(g)))
    color = np.exp(-1j * np.angle(lam) / 3) * np.eye(3)
    rep = native.rep(lam, g / lam, color)
    return rep @ swap @ rep.conj() @ swap


def mark(w):
    a = native.rep(0, 1j * w, np.zeros((3, 3)))
    return -1j * (a + swap @ a.conj() @ swap) @ gamma


def yukawa(h, species):
    legs = [native.IDX[(species, "R", 0, f)] for f in range(native.NG)]
    legs += [native.NP + j for j in legs]
    return -1j * gamma @ src["inner_form"](h) @ identity[:, legs] @ identity[legs, :]


def submasks(mask):
    part = mask
    while part:
        yield part
        part = (part - 1) & mask


def det_mark(factors):
    """Coefficient of every distinct mark in det(I + product factors).

    Each marked factor is I + t A. Products with a repeated t can be dropped
    before taking the mixed derivative.
    """
    product = {0: identity}
    n = 0
    for a, marked in factors:
        if marked:
            product.update({k | (1 << n): b @ a for k, b in product.items()})
            n += 1
        else:
            product = {k: b @ a for k, b in product.items()}
    c = identity + product.pop(0)
    phase, logabs = np.linalg.slogdet(c)
    if n == 0:
        return phase, logabs
    x = {k: solve(c, b, check_finite=False) for k, b in product.items()}
    logs = {k: np.trace(b) for k, b in x.items()}
    power = x
    for degree in range(2, n + 1):
        power = {
            k: sum(power[k ^ j] @ x[j] for j in submasks(k) if (k ^ j) in power)
            for k in x if k.bit_count() >= degree
        }
        for k, b in power.items():
            logs[k] += (-1)**(degree + 1) * np.trace(b) / degree
    result = {0: 1.}
    for k in sorted(x, key=int.bit_count):
        first_bit = k & -k
        result[k] = sum(logs[j] * result[k ^ j] for j in submasks(k) if j & first_bit)
    return phase * result[(1 << n) - 1], logabs


def exact_from_chains(bra, ket, center, sea_shift):
    factors = [(a.conj().T, marked) for a, marked in reversed(bra)]
    return sum(phase * np.exp(logabs + sea_shift)
               for phase, logabs in
               (det_mark(factors + [(z, False)] + ket) for z in center)) / len(center)


def chains(k, e, n, w):
    return [[(e, True), (k, False)], [(w, True), (k, False), (n, True)]]


def block(words, center, shift):
    return np.array([[exact_from_chains(a, b, center, shift) for b in words] for a in words])


def phase_read(block):
    a, b = block[0, 0].real, block[1, 1].real
    v0, v90 = np.array([1., 1.]), np.array([1., 1j])
    i0 = np.vdot(v0, block @ v0).real
    i90 = np.vdot(v90, block @ v90).real
    z = (i0 - a - b + 1j * (a + b - i90)) / 2
    assert abs(z - block[0, 1]) <= 1e-12 * max(a, b)
    return dict(e_only=a, weak_nu_only=b, I0=i0, I90=i90,
                cross=[z.real, z.imag], arm_sum=a+b,
                relative_phase_visibility=2*abs(z)/(a+b))


omega = np.exp(2j*np.pi/3)
c = native.rep(1, np.eye(2), omega*np.eye(3))
c = c @ swap @ c.conj() @ swap
center = [identity, c, c @ c]
assert np.linalg.norm(c @ c @ c - identity) < 1e-13

g_l = np.exp(.27j) * src["quaternion"](np.array([.6, .8j]))
g_r = np.exp(-.43j) * src["quaternion"](np.array([.8j, .6]))
l, r = lift(g_l), lift(g_r)
assert np.linalg.norm(l @ l.conj().T - identity) < 1e-13
assert np.linalg.norm(src["inner_form"](g_l @ first)
                      - l @ src["inner_form"](first) @ l.conj().T) < 1e-12
assert np.linalg.norm(src["inner_form"](g_r @ cut)
                      - r @ src["inner_form"](cut) @ r.conj().T) < 1e-12
assert np.linalg.norm(mark(g_r @ weak @ g_r.conj().T) - r @ mark(weak) @ r.conj().T) < 1e-12
assert abs(l[old["e"], old["e"]] - 1 / np.linalg.det(g_l)) < 1e-13
assert abs(l[old["nu"], old["nu"]] - 1) < 1e-13

rows = {}
for radius in (-.001, .001):
    k, sea = old["half"](radius, [identity] * len(old["links"]))
    ref, ref_sea = old["existing"]["direct_whole"]("S", radius)
    reference = exact_from_chains([], [(ref, False)], center, ref_sea)
    assert abs(reference.imag) < 1e-12 and reference.real > 0
    shift = 2 * sea - np.log(reference.real)
    e, n, w = yukawa(cut, "e"), yukawa(first, "nu"), mark(weak)
    words = chains(k, e, n, w)
    b = block(words, center, shift)
    b_frame = block(chains(r @ k @ l.conj().T, yukawa(g_r @ cut, "e"),
                           yukawa(g_l @ first, "nu"), mark(g_r @ weak @ g_r.conj().T)),
                    center, shift)
    scale = np.max(abs(b))
    assert np.max(abs(b_frame-b)) < 1e-11 * scale
    assert np.max(abs(b-b.conj().T)) < 1e-11 * scale
    assert min(np.linalg.eigvalsh(b)) >= -1e-11 * scale

    # Test the fixed-coordinate cross from the older one-particle words,
    # using both endpoint frames and the original initial columns.
    raw = src["inner_form"](first)
    raw_frame = src["inner_form"](g_l @ first)
    transported = r @ k @ l.conj().T
    col_e, col_n = old["e"], old["nu"]
    z = np.vdot(k @ raw[:, col_e], w @ k @ raw[:, col_n])
    z_frame = np.vdot(transported @ raw_frame[:, col_e],
                      mark(g_r @ weak @ g_r.conj().T) @ transported @ raw_frame[:, col_n])
    assert abs(z_frame - z/np.linalg.det(g_l)) < 1e-12
    assert abs(z_frame - z) > 1e-12

    empty = exact_from_chains([(k, False)], [(k, False)], center, shift)
    rows[str(radius)] = {
        "same_color_reference": reference.real,
        "same_empty_A_over_reference": empty.real,
        "full_exterior_probes": phase_read(b),
        "relative_independent_gauge_frame_error": np.max(abs(b_frame-b)) / scale,
        "fixed_coordinate_cross": [z.real, z.imag],
        "reframed_fixed_coordinate_cross": [z_frame.real, z_frame.imag],
        "native_full_Padm_I0": None,
        "native_full_Padm_I90": None,
    }

print(json.dumps({
    "scope": "Original R032 fixed fermionic history on its full exterior fiber with the true color-center factor only. Original weighted physical R113 endpoint is separate and uncomputed.",
    "width": sum(old["durations"]),
    "cut_time": sum(old["durations"][:3]),
    "one_body_dimension": len(identity),
    "rows": rows,
    "native_physical_endpoint": gate,
}, indent=2))
