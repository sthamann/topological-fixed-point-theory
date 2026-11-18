"""Two scoped checks of original TFPT constructions."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp
from fractions import Fraction
from v128_graded_hull import build_roots

repo = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[3]
contract = repo / "experiments/theory-contracts/chat-triangulation-20261108"
spec = importlib.util.spec_from_file_location("quarter", contract / "check_quarter_reference.py")
quarter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(quarter)


def glue():
    # Derive the two labels independently from the actual root coordinates.
    roots = build_roots()
    root_counts = sp.zeros(4)
    for (d, family), recorded in roots.items():
        if d[0].denominator == 1:
            d_label = int(2*(sum(d) % 2))
        else:
            d_label = 1 + int(2*(sum(x-Fraction(1,2) for x in d) % 2))
        family_label = int((-4*family[0]) % 4)
        assert d_label == family_label == recorded
        root_counts[d_label,family_label] += 1
    assert root_counts == sp.diag(52,64,60,64)
    weights = root_counts / 240
    ones = sp.ones(4,1)
    # Uniform weighting over roots, explicitly declared, never a physical state.
    chi = sp.Matrix([1,sp.I,-1,-sp.I])
    matched = sp.simplify((chi.T * weights * chi.conjugate())[0])
    independent = sp.simplify(((chi.T*weights*ones)[0]) *
                              ((ones.T*weights*chi.conjugate())[0]))
    assert matched == 1 and independent == sp.Rational(1,900)
    return {"root_counts":[int(root_counts[i,i]) for i in range(4)],
            "lost_roots_under_independent_neutral_projection":int(240-root_counts[0,0]),
            "mixed_character_joint":str(matched), "mixed_character_independent":str(independent)}


def local_yukawa():
    D, gamma, majorana = quarter.higgs_matrix, quarter.grading, quarter.mass_matrix
    charged = np.flatnonzero(np.abs(majorana).sum(axis=0)==0)
    active = charged[np.abs(D[np.ix_(charged,charged)]).sum(axis=0)>0]
    assert len(active)==84
    gamma = gamma[np.ix_(active,active)]
    Y = -1j*gamma@D[np.ix_(active,active)]
    plus = np.diag(gamma).real > 0
    minus = ~plus
    assert not np.any(Y[np.ix_(plus,plus)]) and not np.any(Y[np.ix_(minus,minus)])
    B = np.linalg.norm(Y[np.ix_(plus,minus)])**2
    assert np.isclose(B,np.linalg.norm(Y[np.ix_(minus,plus)])**2)
    tags = quarter.native.P
    charges = np.array([(1 if i<len(tags) else -1) if tags[i%len(tags)][0] in ("u","d")
                        else 0 for i in active])
    Z = np.diag(np.exp(2j*np.pi*charges/3))
    assert np.sum(charges[plus])%3 == np.sum(charges[minus])%3 == 0
    assert np.linalg.norm(Z@Y-Y@Z)<1e-12

    # Smooth compactly supported probe in physical coordinate x in [-1,1].
    # Its integral is one, independently of the circumference 2*pi*r.
    nodes, weights = leggauss(96)
    profile = np.exp(-1/(1-nodes**2))
    profile /= weights @ profile
    f = weights*profile
    expected = B/(4*np.pi**2)

    def low_band(radius, eps):
        n = np.arange(int(np.ceil(radius*eps))+1)
        # For each n+m=s, the Fourier argument enumerates (-s,-s+2,...,s)/r.
        # Use cumulative even/odd sums to avoid an O(M^2) pair array.
        d2 = (np.cos(np.outer(n/radius,nodes)) @ f)**2
        shell = np.empty(len(n))
        shell[0::2] = 2*np.cumsum(d2[0::2])-d2[0]
        shell[1::2] = 2*np.cumsum(d2[1::2])
        weight = B*shell/(2*np.pi*radius)**2
        energies = [(n+offset)/radius for offset in (.5,1.5)]
        mass = sum(weight[E<=eps].sum() for E in energies)
        static = sum((weight[E<=eps]/E[E<=eps]).sum() for E in energies)
        return mass,static

    radius = 16384
    rows = []
    for eps in (.2,.1,.05,.025):
        mass, static = low_band(radius,eps)
        rows.append([eps,float(mass/eps**2),float(static/eps)])
    assert abs(rows[-1][1]/expected-1)<.004
    assert abs(rows[-1][2]/(2*expected)-1)<.004
    # Static energy weight per length as the *support* grows in infinite space.
    # Substitute w=(eps*ell)*t**2 to integrate the integrable log at w=0.
    q, w = leggauss(240)
    t, weights_t = (q+1)/2, w/2

    def growing_support(z):
        frequency = z*t*t
        squared = (np.cos(np.outer(frequency,nodes))@f)**2
        return 2*expected*np.dot(weights_t, -4*z*t*np.log(t)*squared)

    slope = (growing_support(80)-growing_support(40))/np.log(2)
    support_limit = (B/(2*np.pi))*(weights@(profile**2))
    assert abs(slope/support_limit-1)<2e-5
    return {"B":float(B), "local_mu_over_epsilon_squared_limit":float(expected),
            "growing_support_static_log_slope":float(slope),
            "growing_support_static_log_coefficient":float(support_limit),
            "local_static_over_epsilon_limit":float(2*expected),
            "eps_and_ratios_at_radius_16384":rows}


print(json.dumps({"E8_correlated_glue":glue(), "R032_fixed_support":local_yukawa()},indent=2))
