"""Determine which original R032 matter response survives mirror pairing."""
import contextlib
import io
import json
from pathlib import Path
import runpy

import numpy as np
from scipy.linalg import eigh, solve

base = Path(__file__).resolve().parent
with contextlib.redirect_stdout(io.StringIO()):
    source = runpy.run_path(str(base / "check_r113_mixed_metric_matter_20261110.py"))

I, grade, right = (source[k] for k in ('I','grade','right'))
native = source['original']['source']['native']
swap = native.Kswap()
particle_hole = grade @ swap  # Linear part of the anti-linear Gamma_F J.
anti_right = [int(np.argmax(abs(swap[:,j]))) for j in right]
assert not set(anti_right).intersection(right)
assert np.array_equal(swap @ swap.conj(), I)
assert np.array_equal(particle_hole @ particle_hole.conj(), -I)
G, R_m, R_h, R_mh = (source[k] for k in ('G','R_m','R_h','R_mh'))
n, size = source['jet'].real, len(right)
compress = np.ix_(right, right)
C0, C_m, C_h, C_mh = tuple(M[compress] for M in (
    I-G, G@R_m@G, G@R_h@G,
    G@R_mh@G-G@R_m@G@R_h@G-G@R_h@G@R_m@G))
variance = np.real(np.array([
    np.trace(C0-C0@C0),
    np.trace(C_m-2*C0@C_m),
    np.trace(C_h-2*C0@C_h),
    np.trace(C_mh-2*C0@C_mh-2*C_m@C_h),
]))
oriented = 2*n - np.array([size,0,0,0])
paired = np.array([np.zeros(4), oriented, 2*variance]).T


def transition(momentum, sigma, h):
    R = I.copy()
    sea = 0.
    for field, width in zip(source['history'], source['durations']):
        phi = source['inner_form']((1+h)*np.asarray(field)) + source['M']
        generator = momentum*np.exp(-sigma)*grade-1j*grade@phi
        assert np.max(abs(particle_hole@generator.conj()+generator@particle_hole)) < 2e-12
        values, vectors = eigh(generator)
        R = ((vectors*np.exp(-width*values))@vectors.conj().T)@R
        sea += width*values[:len(I)//2].sum()
    return R, sea


def direct(sigma, h):
    p=source['p0']
    R_p, sea_p = transition(p, sigma, h)
    R_n, sea_n = transition(-p, sigma, h)
    # Matrix and characteristic-function forms of mirror complementation.
    assert np.max(abs(grade@R_n@grade@R_p.conj().T-I)) < 2e-12
    assert np.max(abs(swap@R_p.conj()@swap-R_n)) < 2e-12
    assert np.max(abs(R_p.T@particle_hole@R_p-particle_hole)) < 2e-12
    for z in (-.41, .21, .08+.13j):
        diag = np.ones(len(I), dtype=complex)
        diag[right] = np.exp(z)
        plus_diag = diag.conj()**-1
        test = np.linalg.det(I+diag[:,None]*R_n)
        mirrored = np.exp(size*z)*np.linalg.det(I+plus_diag[:,None]*R_p).conj()
        assert abs(test/mirrored-1) < 2e-12
    cp, cn = (solve(I+R, R)[compress] for R in (R_p,R_n))
    c_anti = solve(I+R_p, R_p)[np.ix_(anti_right,anti_right)]
    assert np.max(abs(c_anti.conj()-cn)) < 2e-12
    ep,en = np.trace(cp),np.trace(cn)
    vp,vn = np.trace(cp-cp@cp),np.trace(cn-cn@cn)
    assert np.max(abs(cp+cn.conj().T-np.eye(size))) < 2e-12
    assert abs(sea_p-sea_n) < 1e-12
    return np.real([ep+en-size,ep-en,vp+vn+(ep+en-size)**2])


central = source['central_jet'](direct)
assert np.max(abs(paired-central)) < 5e-7
assert abs(paired[-1,1]) > 1e-4 and abs(paired[-1,2]) > 1e-4
print(json.dumps({"jet_columns": ["centered mirror-even occupation", "oriented occupation",
                                  "square of centered mirror-even occupation"],
                  "jet_rows": ["value","metric","Higgs","mixed"],
                  "analytic_jet": paired.tolist(),
                  "independent_full_pair_finite_difference": central.tolist(),
                  "max_jet_difference": float(np.max(abs(paired-central)))}, indent=2))
