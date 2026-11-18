"""Compare the original C6 and Fuchs kernels to the printed Flavor ladder."""
from itertools import combinations, product
import json

import numpy as np
import sympy as sp

from tfpt_constants import phi0
from tfpt_explorer.flavor_path_transport import _integrate, integrate_path_transport
from v141_deck_selection import G, P, TA


def exact_ladder():
    x = sp.Symbol('x', positive=True)
    c = tuple(map(sp.Rational, ['55/117', '34/47', '3/26']))
    lengths_up, lengths_down = (7, 3, 0), (7, 5, 2)
    powers_up, powers_down = (4, 2, 0), (4, 3, 2)
    stripped = [sp.simplify(f*x**(ku-kd)*(x*(1-x))**sp.Rational(ld-lu, 2))
                for f, lu, ld, ku, kd in zip(c, lengths_up, lengths_down, powers_up, powers_down)]
    assert all(sp.cancel(f-g) == 0 for f, g in zip(stripped, [c[0], c[1]*(1-x), c[2]*(1-x)/x]))
    volume_squared = sp.factor(sp.prod(stripped)**2)
    assert sp.cancel(volume_squared-(sp.prod(c)**2*(1-x)**4/x**2)) == 0
    seed = 1/(6*sp.pi)+3/(256*sp.pi**4)
    assert abs(float(seed)-float(phi0)) < 1e-16

    # This shape assumes the Gram forms and the Yukawa map use a common
    # diagonal generation frame. Test the stated compiler D4 action in it.
    d = sp.diag(*[sp.cancel(r/f) for r, f in zip(stripped, c)])
    reflection = P.inv()*TA*P
    deck = P.inv()*G*P
    for t in (reflection, deck):
        swapped = t*d*t.inv()
        assert sp.simplify(swapped-d-sp.diag(0, (1-x)**2/x, -(1-x)**2/x)) == sp.zeros(3)
        assert sp.simplify(t*(d+swapped)*t.inv()-(d+swapped)) == sp.zeros(3)

    # Common rescaling of individual paired generations can preserve all
    # ratios and even both determinants while varying electron/up.
    cu, cd = sp.diag(sp.Rational(1,2), sp.Rational(34,41), sp.Rational(4,13)), sp.diag(
        sp.Rational(117,110), sp.Rational(47,41), sp.Rational(8,3))
    q = sp.diag(sp.Rational(11,10), sp.Rational(10,11), 1)
    assert (q*cu)*(q*cd).inv() == cu*cd.inv() == sp.diag(*c)
    assert (q*cu).det() == cu.det() and (q*cd).det() == cd.det()
    ep, mu, tau = map(sp.Rational, ['16/7','4/3','7/6'])
    e_over_u_gauge = ep*x/cu[0,0]
    assert sp.simplify((ep*x/(q*cu)[0,0])/e_over_u_gauge) == sp.Rational(10,11)

    return {
        'stripped_ratios': [str(sp.factor(r)) for r in stripped],
        'stripped_ratios_at_seed': [float(r.subs(x, seed)) for r in stripped],
        'conditional_gram_det_ratio': str(volume_squared),
        'conditional_gram_det_ratio_at_seed': float(volume_squared.subs(x, seed)),
        'reflection_odd_factor_before_halving': str((1-x)**2/x),
        'lepton_pair_ratios': [str(ep/mu*x**2), str(mu/tau*x)],
        'e_over_u_in_printed_gauge_only': float(e_over_u_gauge.subs(x,seed)),
        'e_over_u_after_determinant_preserving_counterexample': float((e_over_u_gauge*sp.Rational(10,11)).subs(x,seed)),
    }


def c6_subset_probe(target):
    """Optional numerical diagnostic for coordinate spectral triples only."""
    roots = np.exp(2j*np.pi*np.arange(6)/6)
    subsets = tuple(combinations(range(6),3))
    rows = {}
    poles = {'current_half_step': 0.5, 'v45_algebraic_pole': ((794-7*np.sqrt(9961))/2187)**(1/6)}
    for label, delta in poles.items():
        down = [float(np.prod(abs(1/3-delta*roots[list(j)])**2)) for j in subsets]
        up = [float(np.prod(abs(2/3-delta*roots[list(j)])**2)) for j in subsets]
        for power in (1, -1, 0.5, -0.5):
            error, val, i, j = min((abs(val-target),val,i,j)
                for i,j in product(range(len(subsets)), repeat=2)
                for val in [(down[i]/up[j])**power])
            rows[f'{label}_{power}'] = [error,val,subsets[i],subsets[j]]
    assert min(row[0] for row in rows.values()) > 1e-5
    return rows


def fuchs_open_probe(seed, tol):
    z = seed**.25
    into_z = _integrate(lambda t: (t*z, z), tol, tol/100)[0]
    forward = integrate_path_transport(z, orientation=1, rtol=tol, atol=tol/100)
    backward = integrate_path_transport(z, orientation=-1, rtol=tol, atol=tol/100)
    deck = np.diag([1,1j,-1j])
    m0 = np.linalg.solve(into_z, forward['monodromy'] @ into_z)
    loops = [np.linalg.matrix_power(deck,k) @ m0 @ np.linalg.matrix_power(deck,-k) for k in range(4)]
    basis = [np.diag(row) for row in np.eye(3)]
    systems = np.stack([np.concatenate([(m.conj().T@h@m-h).ravel() for m in loops]) for h in basis], axis=1)
    _, singular, vh = np.linalg.svd(np.vstack([systems.real, systems.imag]))
    h0 = np.diag(3*vh[-1]/sum(vh[-1]))
    assert min(np.diag(h0)) > 0 and singular[-2] > 0.1 and singular[-1] < 1e-8

    def transported_metric(transport):
        inverse = np.linalg.inv(transport)
        return inverse.conj().T@h0@inverse

    h_z = transported_metric(into_z)
    lower = forward['P'] @ into_z
    upper = backward['P'] @ into_z
    h_reflected = transported_metric(lower)
    metric_arc_difference = np.linalg.norm(h_reflected-transported_metric(upper))
    effective = [abs(np.linalg.det(g))*np.sqrt(np.linalg.det(h).real/np.linalg.det(h0))
                 for g,h in ((into_z,h_z),(lower,h_reflected),(upper,h_reflected))]
    raw = [np.linalg.det(g) for g in [into_z, lower, upper]]
    expected = [1-seed, -(1-seed)/seed, -(1-seed)/seed]
    assert max(abs(a-b) for a,b in zip(raw,expected)) < 1e-8
    assert max(abs(e-1) for e in effective) < 1e-8 and metric_arc_difference < 1e-8
    return {'endpoint': z, 'raw_determinants': [r.real for r in raw],
            'invariant_determinant_norms': effective, 'two_arc_metric_difference': metric_arc_difference}


if __name__ == '__main__':
    result = exact_ladder()
    result['spectral_triple_probe'] = c6_subset_probe(result['conditional_gram_det_ratio_at_seed'])
    result['original_fuchs_open_probe'] = [fuchs_open_probe(float(phi0),tol) for tol in (2e-11,5e-13)]
    print(json.dumps(result, default=float, indent=2))
