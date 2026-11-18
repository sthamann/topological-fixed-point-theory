"""Temporal discrimination tests on the unchanged R032 charged block."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm
from scipy.special import expit

earlier_path = Path(sys.argv.pop(1))
spec = importlib.util.spec_from_file_location("paired", earlier_path)
paired = importlib.util.module_from_spec(spec)
spec.loader.exec_module(paired)

t = paired.t
y2 = paired.y2
color_y2 = 6*sum(np.linalg.norm(y)**2 for y in paired.source.yukawa[:2])
assert 0 < color_y2 < y2


def native_probes():
    p = 0.25
    width = 0.8
    h0 = p*paired.G
    mark = -1j*paired.G@paired.Y
    identity = np.eye(len(h0))

    def raw_log(R):
        assert np.linalg.norm(R-R.conj().T)<1e-9
        return np.log1p(np.linalg.eigvalsh(R)).sum()

    base = raw_log(expm(-width*h0))

    def log_harmonic(x,n):
        frequency = 2*np.pi*n/width

        def equation(s, flat):
            f = 1 if n==0 else np.sqrt(2)*np.cos(frequency*s)
            return (-(h0+x*f*mark)@flat.reshape(h0.shape)).ravel()

        sol = solve_ivp(equation,(0,width),identity.astype(complex).ravel(),
                        rtol=2e-12,atol=2e-14,method="DOP853")
        assert sol.success
        return raw_log(sol.y[:,-1].reshape(h0.shape))

    def log_windows(x,z):
        outside = expm(-width*(h0+np.sqrt(2)*x*mark)/4)
        inside = expm(-width*(h0+np.sqrt(2)*z*mark)/2)
        return raw_log(outside@inside@outside)

    def cross(h):
        return (log_windows(h,h)+log_windows(-h,-h)
                -log_windows(h,-h)-log_windows(-h,h))/(4*h*h)

    def harmonic_jet(h,n):
        return (log_harmonic(h,n)+log_harmonic(-h,n)-2*base)/(2*h*h)

    rows = []
    for n in (0,1,2):
        exact = y2*np.tanh(width*p/2)*2*p*width/(4*p*p+(2*np.pi*n/width)**2)
        observed = (4*harmonic_jet(0.001,n)-harmonic_jet(0.002,n))/3
        assert abs(exact-observed)<3e-7
        rows.append([n,float(exact),float(observed)])
    exact_cross = y2*np.tanh(width*p/2)**2/p**2
    observed_cross = (4*cross(0.001)-cross(0.002))/3
    assert abs(exact_cross-observed_cross)<3e-7
    return dict(width=width,p=p,harmonics_n_exact_numeric=rows,
                separated_cross_exact=float(exact_cross),separated_cross_numeric=float(observed_cross))


def native_center():
    p = 0.25
    positions = paired.charged[paired.connected]
    tags = paired.source.native.P
    charges = [(1 if i<len(tags) else -1) if tags[i%len(tags)][0] in ("u","d") else 0
               for i in positions]
    Z = np.diag(np.exp(2j*np.pi*np.asarray(charges)/3))
    assert sum(q!=0 for q in charges)==72
    assert np.linalg.norm(Z@paired.Y-paired.Y@Z)<1e-12
    identity = np.eye(len(Z))

    def correction(x):
        H = p*paired.G-1j*x*paired.G@paired.Y
        R = expm(-t*H)
        sign, logtwist = np.linalg.slogdet(identity+Z@R)
        _, lograw = np.linalg.slogdet(identity+R)
        chi = sign*np.exp(logtwist-lograw)
        assert abs(chi.imag)<1e-12
        return np.log1p(2*chi.real)

    s = expit(t*p)*expit(-t*p)
    chi = (1-3*s)**36
    assert abs(correction(0)-np.log1p(2*chi))<1e-12
    base = correction(0)

    def jet(h):
        return (correction(h)+correction(-h)-2*base)/(2*h*h)

    direct = (4*jet(0.0002)-jet(0.0004))/3
    expected = (2*chi/(1+2*chi))*3*color_y2*s*np.tanh(t*p/2)*(t/(2*p))/(1-3*s)
    assert abs(direct-expected)<2e-7
    return dict(p=p,log_projected_over_raw=base,analytic_response_correction=float(expected),
                numeric_response_correction=float(direct))


def temporal_densities():
    def scores(p,n):
        frequency = 2*np.pi*n/t
        tanh = np.tanh(t*p/2)
        polar = -2*expit(-t*p)*y2*2*p*t/(4*p*p+frequency*frequency)
        subtracted = -y2*t*tanh/(2*p)*frequency**2/(4*p*p+frequency**2)
        return polar,subtracted

    def separated(p):
        tanh = np.tanh(t*p/2)
        return y2*(tanh/p)**2, -2*expit(-t*p)*y2*tanh/p**2

    def center_correction(p,n):
        s = expit(t*p)*expit(-t*p)
        lxi = 36*np.log1p(-3*s).sum()
        alpha = 2*np.exp(lxi)/(1+2*np.exp(lxi))
        jet = 2*p*t/(4*p*p+(2*np.pi*n/t)**2)
        dlogxi = 3*color_y2*np.sum(s*np.tanh(t*p/2)*jet/(1-3*s))
        return alpha*dlogxi

    harmonic_limits = {}
    for n in (1,2):
        polar = quad(lambda p: 2*scores(p,n)[0],0,np.inf,epsabs=1e-9)[0]
        difference = quad(lambda p: 2*scores(p,n)[1],0,np.inf,epsabs=1e-9)[0]
        harmonic_limits[n] = [float(polar),float(difference)]
    cross_limit = quad(lambda p: 2*separated(p)[0],0,np.inf,epsabs=1e-9)[0]
    cross_tail = quad(lambda p: 2*separated(p)[0],64,np.inf,epsabs=1e-9)[0]
    tails = [quad(lambda p: 2*scores(p,n)[1],64,np.inf,epsabs=1e-10)[0] for n in (1,2)]
    rows = []
    center = []
    for r in (1,2,4,8,32,64,128,256,512,1024):
        p = (2*np.arange(128*r)+1)/(4*r)
        polar_cross = separated(p)[1].sum()/r
        raw_cross = separated(p)[0].sum()/r+cross_tail
        polar_zero = scores(p,0)[0].sum()/r
        polar_one, subtracted_one = scores(p,1)
        polar_two, subtracted_two = scores(p,2)
        rows.append([r,*[float(v) for v in
                       (polar_zero,polar_one.sum()/r,polar_two.sum()/r,
                        subtracted_one.sum()/r+tails[0],subtracted_two.sum()/r+tails[1],
                        raw_cross,polar_cross)]])
        if r<=8:
            center.append([r,float(center_correction(p,0)),float(center_correction(p,1))])
    assert abs((rows[-1][1]-rows[-2][1])/np.log(2)+t*y2)<1e-7
    assert abs((rows[-1][-1]-rows[-2][-1])/np.log(2)+t*y2)<0.01
    assert abs(rows[-1][-2]-cross_limit)<0.003
    assert abs(rows[-1][4]-harmonic_limits[1][1])<0.001
    assert abs(rows[-1][5]-harmonic_limits[2][1])<0.001
    assert 0 < center[0][1] and 0 < center[0][2] and center[-1][1] < 1e-12
    return dict(time=t,squared_charged_weight=float(y2),squared_colored_weight=float(color_y2),
                sewn_diagonal_log_slopes=[float(-t*y2)]*3,
                polar_diagonal_log_slopes=[float(-t*y2),0,0],
                sewn_disjoint_cross_log_slope=0,polar_disjoint_cross_log_slope=float(-t*y2),
                harmonic_limits_polar_and_equal_norm_difference=harmonic_limits,
                sewn_disjoint_cross_density_limit=float(cross_limit),
                rows_r_polar0_polar1_polar2_difference1_difference2_cross_sewn_cross_polar=rows,
                center_rows_r_absolute_corrections_n0_n1=center)


def reflected_time_gram():
    radius = 1
    p = (2*np.arange(128)+1)/(4*radius)
    s = expit(t*p)*expit(-t*p)
    chi = np.exp(36*np.log1p(-3*s).sum())
    coefficient = 2*s*(y2+2*chi/(1+2*chi)*3*color_y2*s/(1-3*s))

    def kernel(tau):
        return float(np.dot(coefficient,np.cosh(p*(t-2*tau))))

    cuts = np.array([0.11,0.22])*t
    gram = np.array([[kernel(a+b) for b in cuts] for a in cuts])
    assert np.linalg.eigvalsh(gram)[0]>0
    # At one mode the determinant is the product of the two positive
    # spectral weights times the squared difference of the exponential vectors.
    a,b = cuts
    z = p[0]
    c = coefficient[0]
    first = np.array([[c*np.cosh(z*(t-2*(u+v))) for v in cuts] for u in cuts])
    exact = c*c*np.sinh(2*z*(a-b))**2
    assert abs(np.linalg.det(first)-exact)<1e-12
    # The independently recorded bare proper-time tail has the opposite sign.
    old = lambda tau: tau*np.exp(-tau*tau)
    bare = old(1)*old(2)-old(1.5)**2
    assert bare<0
    return dict(cuts_fraction=[0.11,0.22],center_projected_gram=gram.tolist(),
                center_projected_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
                original_bare_action_atom_dimensionless_determinant=float(bare))


if __name__ == "__main__":
    print(json.dumps(dict(native=native_probes(),native_center=native_center(),
                          densities=temporal_densities(),reflected_gram=reflected_time_gram()),indent=2))
