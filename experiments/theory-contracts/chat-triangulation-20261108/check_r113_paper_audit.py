"""Check transport covariance of the two-window R032 comparison."""
import importlib.util
import json
import sys

import numpy as np
from scipy.linalg import expm, expm_frechet

path = sys.argv.pop(1)
spec = importlib.util.spec_from_file_location("windows", path)
windows = importlib.util.module_from_spec(spec)
spec.loader.exec_module(windows)
paired = windows.paired

p, width, first, middle = 0.25, 0.8, 0.2, 0.4
tags = paired.source.native.P
positions = paired.charged[paired.connected]
W = np.diag([-1 if tags[i%len(tags)][1]=="L" else 1 for i in positions])
charges = [(1 if i<len(tags) else -1) if tags[i%len(tags)][0] in ("u","d") else 0
           for i in positions]
Z = np.diag(np.exp(2j*np.pi*np.asarray(charges)/3))
I = np.eye(len(W))
h0 = p*paired.G
y = -1j*paired.G@paired.Y
assert np.linalg.norm(W@h0@W-h0)<1e-13
assert np.linalg.norm(W@y@W+y)<1e-13
assert np.linalg.norm(W@Z-Z@W)<1e-13

outer = expm(-first*(h0+0.35*y))
changed_outer = expm(-first*(h0-0.35*y))
assert np.linalg.norm(changed_outer-W@outer@W)<1e-12
inner, inner_jet = expm_frechet(-middle*(h0-0.55*y),middle*0.55*y)
base, dbase = outer@inner@outer, outer@inner_jet@outer
fixed, dfixed = changed_outer@inner@changed_outer, changed_outer@inner_jet@changed_outer
carried, dcarried = changed_outer@W@inner@W@changed_outer, changed_outer@W@inner_jet@W@changed_outer
assert np.linalg.norm(carried-W@base@W)<1e-11
assert np.linalg.norm(dcarried-W@dbase@W)<1e-11


def read(R, jet):
    value, derivative = 0j,0j
    for charge in (I,Z,Z@Z):
        matrix = I+charge@R
        sign, logdet = np.linalg.slogdet(matrix)
        det = sign*np.exp(logdet)
        value += det
        derivative += det*np.trace(np.linalg.solve(matrix,charge@jet))
    assert abs(value.imag)<1e-12*abs(value.real)
    assert abs(derivative.imag)<1e-12*abs(value.real)
    return np.log(value.real/3),(-derivative/(middle*value)).real


a = read(base,dbase)
b = read(fixed,dfixed)
c = read(carried,dcarried)
assert np.allclose(a,c,rtol=0,atol=1e-12)
assert abs(a[1]-b[1])>1
assert abs(a[1]+3.214114733882643)<1e-11
assert abs(b[1]+13.222151965587742)<1e-11
print(json.dumps(dict(base_raw_log_and_y=a,changed_at_fixed_links_raw_log_and_y=b,
                      gauge_transported_raw_log_and_y=c,
                      transport_matrix_residual=float(np.linalg.norm(carried-W@base@W)),
                      transport_marked_matrix_residual=float(np.linalg.norm(dcarried-W@dbase@W))),indent=2))
