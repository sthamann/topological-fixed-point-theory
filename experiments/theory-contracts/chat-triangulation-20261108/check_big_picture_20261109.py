"""Cross-checks on the existing R032, R048 and v18 TFPT branches."""
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp
import v18_quark_yukawa as flavor
from tfpt_constants import phi0

mp.mp.dps = 65
repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[3]
contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
B = mp.mpf(str(json.loads((contract / "r113-temporal-cross.json").read_text())["densities"]["squared_charged_weight"]))


def memory():
    rows = []
    for k in (2, 4, 6, 8):
        times = [mp.mpf(".05")*2**i for i in range(k)]
        limiting = mp.matrix([[B/(t+s) for s in times] for t in times])
        det = mp.det(limiting)
        cauchy = B**k * mp.fprod((t-s)**2 for i,t in enumerate(times) for s in times[:i])
        cauchy /= mp.fprod(t+s for t in times for s in times)
        assert det > 0 and abs(det/cauchy - 1) < mp.mpf("1e-45")
        radius = mp.mpf(1024)
        finite = mp.matrix([[B/(2*radius*mp.sinh((t+s)/(2*radius))) for s in times] for t in times])
        error = max(abs(finite[i,j]/limiting[i,j] - 1) for i in range(k) for j in range(k))
        # Repeat with the archived finite energy band (0,1] kept fixed.
        bounded = mp.matrix([[B*(-mp.expm1(-(t+s)))/(t+s) for s in times] for t in times])
        bounded_finite = mp.matrix([[finite[i,j]*(-mp.expm1(-(times[i]+times[j])))
                                     for j in range(k)] for i in range(k)])
        assert mp.det(bounded) > 0 and mp.det(bounded_finite) > 0
        # Dimensionless normalized Gram matrices have unit diagonal.
        normalized = mp.matrix([[limiting[i,j]/mp.sqrt(limiting[i,i]*limiting[j,j])
                                 for j in range(k)] for i in range(k)])
        smallest = mp.eigsy(normalized, eigvals_only=True)[0]
        rows.append([k, float(smallest), float(error)])
    return rows


def matter_curvature():
    x, y, a, b, gxx, gxy, gyy, yi, yj, xi = sp.symbols(
        "x y a b Gxx Gxy Gyy yi yj xi", real=True)
    gaussian = sp.exp((gxx*x*x + 2*gxy*x*y + gyy*y*y)/2)
    moment = lambda i,j: sp.diff(gaussian, x, i, y, j).subs({x:0,y:0})
    # The R048 finite regulated Gaussian, including its possible linear mean.
    linear = 2*a*b*gxy
    residual = (moment(2,2)-moment(2,0)*moment(0,2))/4
    assert sp.expand(residual - gxy**2/2) == 0
    assert moment(1,2) == moment(2,1) == 0
    block = sp.Matrix([[yi*yj*gxy, yi*xi*sp.sqrt(2)*b*gxy],
                       [xi*sp.sqrt(2)*a*yj*gxy, xi**2*(linear+residual)]])
    assert sp.simplify(block.det() - xi**2*yi*yj*gxy**3/2) == 0
    return {"quadratic_remainder":str(residual.simplify()), "cross_minor":str(block.det().simplify())}


def flavor_checksum():
    e, mu, tau = map(mp.mpf, (".510999","105.66","1776.86"))
    a, b = mu/tau, e/mu
    phi_a, phi_b = mp.mpf(7)/8*a, mp.sqrt(mp.mpf(7)/12*b)
    frozen_checksum = b/a**2
    c_form = mp.mpf(21)/16
    numerical_word_checksum = (flavor.yuk("e")/flavor.yuk("mu"))/(flavor.yuk("mu")/flavor.yuk("tau"))**2
    angle = lambda seed: 45*seed/mp.pi**2
    vals = {
        "c_form_checksum":c_form,
        "frozen_pole_checksum":frozen_checksum,
        "relative_pole_gap":frozen_checksum/c_form - 1,
        "rounded_word_checksum":numerical_word_checksum,
        "seed_fixed":phi0, "seed_from_mu_over_tau":phi_a, "seed_from_e_over_mu":phi_b,
        "beta_degrees_fixed":angle(phi0),
        "beta_degrees_from_mu_over_tau":angle(phi_a), "beta_degrees_from_e_over_mu":angle(phi_b),
    }
    assert abs(phi_a-phi_b) > mp.mpf(".001")
    ze, zmu, ztau = sp.symbols("Ze Zmu Ztau", positive=True)
    assert sp.cancel((ze/zmu)/(zmu/ztau)**2) == ze*ztau**2/zmu**3
    return {k:float(v) for k,v in vals.items()}


print(json.dumps({"memory": memory(), "matter_curvature": matter_curvature(),
                  "flavor_checksum": flavor_checksum()}, indent=2))
