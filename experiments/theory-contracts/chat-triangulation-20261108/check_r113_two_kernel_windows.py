"""Original conditional R032 ordered windows and exact two-kernel identities."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm, expm_frechet

earlier_path = Path(sys.argv.pop(1))
spec = importlib.util.spec_from_file_location("paired", earlier_path)
paired = importlib.util.module_from_spec(spec)
spec.loader.exec_module(paired)


def actual_window_trace(x_first,x_middle):
    p = 0.25
    time = 0.8
    first, middle = time/4, time/2
    h0 = p*paired.G
    y = -1j*paired.G@paired.Y
    h_first = h0+x_first*y
    y_middle = x_middle*y
    h_middle = h0+y_middle
    outer = expm(-first*h_first)
    inner, jet = expm_frechet(-middle*h_middle,-middle*y_middle)
    transfer = outer@inner@outer
    derivative = outer@jet@outer
    ident = np.eye(len(transfer))

    positions = paired.charged[paired.connected]
    tags = paired.source.native.P
    charge = [(1 if i<len(tags) else -1) if tags[i%len(tags)][0] in ("u","d") else 0
              for i in positions]
    center = np.diag(np.exp(2j*np.pi*np.asarray(charge)/3))
    projectors = [ident,center,center@center]
    assert np.linalg.norm(center@transfer-transfer@center)<1e-11

    def contracted(R, dR=None):
        values, variations = [], []
        for charge_action in projectors:
            matrix = ident+charge_action@R
            sign, logdet = np.linalg.slogdet(matrix)
            det = sign*np.exp(logdet)
            values.append(det)
            if dR is not None:
                variations.append(det*np.trace(np.linalg.solve(matrix,charge_action@dR)))
        value = sum(values)/3
        assert abs(value.imag)<1e-11*value.real and value.real>0
        if dR is None:
            return float(value.real)
        variation = sum(variations)/3
        assert abs(variation.imag)<1e-10*max(1,abs(variation.real))
        return float(value.real),float(variation.real)

    n0, d0 = contracted(transfer,derivative)
    E_first = np.hypot(p,x_first*paired.y)
    E_middle = np.hypot(p,x_middle*paired.y)
    log_sea = -2*first*E_first.sum()-middle*E_middle.sum()
    log_ref_s = np.log(contracted(expm(-time*h0)))-time*p*len(paired.y)
    log_z0_S = np.log(n0)+log_sea-log_ref_s
    normalized_y = -d0/(middle*n0)
    n_U = x_first*x_first*n0
    n_Y = -d0/middle
    n_UY = x_first*x_first*n_Y
    connected = n_UY/n0-(n_U/n0)*(n_Y/n0)
    assert abs(connected)<1e-13

    def value(k, move_sea=False):
        R = outer@expm(-middle*(h_middle+k*y_middle))@outer
        log_value = np.log(contracted(R))
        if move_sea:
            log_value -= middle*np.hypot(p,(1+k)*x_middle*paired.y).sum()
        return log_value

    def slope(h, move_sea):
        return -(value(h,move_sea)-value(-h,move_sea))/(2*h*middle)

    h = 0.001
    frozen = (4*slope(h/2,False)-slope(h,False))/3
    moving = (4*slope(h/2,True)-slope(h,True))/3
    extra = (x_middle**2*paired.y**2/E_middle).sum()
    assert abs(frozen-normalized_y)<3e-9
    assert abs(moving-normalized_y-extra)<3e-9
    return dict(p=p,time=time,early_U_window=[time/16,time/8],later_Y_window=[first,first+middle],
                early_higgs=x_first,later_higgs=x_middle,
                conditional_Z3_log_Z0_over_S=float(log_z0_S),
                conditional_Z3_late_Y_mean=float(normalized_y),
                conditional_Z3_late_Y_finite_difference=float(frozen),
                conditional_Z3_late_Y_with_moving_sea=float(moving),
                extra_sea_term=float(extra),
                conditional_connected_UY=float(connected))


def algebra():
    w = sp.symbols("w0:3")
    u = sp.symbols("u0:3")
    K = sp.symbols("K0:3")
    L = sp.symbols("L0:3")
    z = sum(w[i]*K[i] for i in range(3))
    v = sum(w[i]*u[i]*K[i] for i in range(3))
    y = sum(w[i]*L[i] for i in range(3))
    m = sum(w[i]*u[i]*L[i] for i in range(3))
    doubled = sum(w[i]*w[j]*(u[i]-u[j])*(L[i]*K[j]-K[i]*L[j])
                  for i in range(3) for j in range(3))
    assert sp.expand(doubled-2*(z*m-v*y))==0
    assert sp.expand(doubled.subs(K[0],0)-2*(z*m-v*y).subs(K[0],0))==0

    # Product rule for disjoint windows when only first local weight scores remain.
    x, b = sp.symbols("x b")
    f, g, h, j, q, k, a, c = sp.symbols("f g h j q k a c")
    raw = q+x*f+b*g+x*b*h
    weight = k*(1+x*a)*(1+b*c)
    actual = sp.diff(weight*raw,x,b).subs({x:0,b:0})
    assert sp.expand(actual-k*(h+a*g+c*f+a*c*q))==0
    return dict(doubled_identity_exact=True,identity_through_zero_unmarked_fiber=True,
                disjoint_weight_first_scores_product_retained=True)


if __name__=="__main__":
    baseline = actual_window_trace(0.35,-0.55)
    same_radial_history = actual_window_trace(-0.35,-0.55)
    global_weak_flip = actual_window_trace(-0.35,0.55)
    key = "conditional_Z3_late_Y_mean"
    assert abs(baseline[key]-same_radial_history[key])>0.05
    assert abs(baseline[key]-global_weak_flip[key])<1e-12
    print(json.dumps(dict(conditional_original=baseline,
                          same_radial_history_different_early_sign=same_radial_history,
                          global_weak_flip=global_weak_flip,algebra=algebra()),indent=2))
