"""Cross-check three hypotheses with existing quarter, Hardy, and TT operators."""
import contextlib
from collections import Counter
import importlib.util
import io
from itertools import product
import json
from pathlib import Path
import runpy
import sys

import mpmath as mp
import numpy as np
from scipy.linalg import block_diag, eigh
import sympy as sp


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"


def existing(name):
    with contextlib.redirect_stdout(io.StringIO()):
        return runpy.run_path(str(base / name))


def phase_test():
    pair = existing("check_original_interference_quarter_pair_20261110.py")
    probe = pair["probe"]
    eye = pair["identity"]
    e, n, weak = pair["marks"]
    double = lambda a: block_diag(a, a)
    # -i times the actual half-translation on the first pair, i and -i.
    quarter = block_diag(eye, -eye)
    p = sp.Symbol("p", positive=True)
    f0, f1, f2 = sp.symbols("f0 f1 f2", real=True)
    momentum = sp.diag(p,-p)
    profile = sp.Matrix([[f0,f1+sp.I*f2],[f1-sp.I*f2,f0]])
    assert sp.expand((profile*momentum+momentum*profile)/2) == f0*momentum
    colors = [double(z) for z in pair["center"]]
    row = {}
    for radius in (-.001, .001):
        k0, s0 = pair["half"](radius, 1)
        k1, s1 = pair["half"](radius, -1)
        r0, sr0 = pair["half"](radius, 1, reference=True)
        r1, sr1 = pair["half"](radius, -1, reference=True)
        k, ref = block_diag(k0, k1), block_diag(r0.conj().T @ r0, r1.conj().T @ r1)
        standard = (double(e), double(n), double(weak))
        def run(marks):
            return pair["evaluate"](k, s0+s1, ref, 2*(sr0+sr1), colors, marks)
        even, re, ze = run(standard)
        odd, ro, zo = run((*standard[:2], standard[2] @ quarter))
        reversed_tag, rr, zr = run((*standard[:2], -standard[2] @ quarter))
        assert re == ro == rr and ze == zo == zr
        assert abs(even[0,1].imag) < 1e-25 and abs(odd[0,1].real) < 1e-25
        assert abs(odd[0,1].imag) > 1e-17
        assert abs(odd[0,1] + reversed_tag[0,1]) < 1e-25
        assert np.max(abs(np.diag(odd)-np.diag(reversed_tag))) < 1e-25
        row[str(radius)] = dict(even=pair["read"](even), odd=pair["read"](odd),
            reversed_odd=pair["read"](reversed_tag), same_reference=re, same_empty=ze)
    return row


def family_test():
    f = existing("check_fuchs_two_ends_and_tt_20261109.py")
    mp.mp.dps = 50
    e, m, t, r = [f[x] for x in ("e", "mu", "tau", "rel_period")]
    phase = mp.exp(-5j*mp.pi/12)
    i = mp.j
    target = [mp.mpf(x) for x in f["poles"]]
    target_norm = sum(y*y for y in target)
    s2 = sum((target[j]*target[k])**2 for j,k in ((0,1),(0,2),(1,2)))/target_norm**2
    s3 = (target[0]*target[1]*target[2])**2/target_norm**3
    x, y = f["a"], f["b"]
    first, second = f["f"], f["n"]
    mixed = sp.simplify(first.H @ second + second.H @ first)
    assert mixed == sp.I*y*y*(f["unit"](0,1)-f["unit"](1,0))
    assert first.H @ first == sp.diag(0,y*y,x*x)
    assert second.H @ second == sp.diag(y*y,0,x*x)
    a1 = e*e+m*m+t*t
    q = phase**2
    factor = 2*(e*e+m*m)-8*s2*a1
    c0 = factor*(1+r*r)-2*t*(e+m*r*r)*mp.re(q)
    c1 = 2*r*factor-4*t*r*(e+m)*q
    c2_phase = -2*t*(e*r*r+m)*q
    linear_lower_bound = c0-abs(c1)-abs(c2_phase)
    assert linear_lower_bound > mp.mpf(".0915")
    assert 1-4*s2 > 0 and (abs(r)-1)**2 > 0
    row = {}
    for label, u in zip(("1","i","-1","-i"), (1,i,-1,-i)):
        # Recovered endpoints: F=(1,r), N=-i(r,1), both with original phase.
        z1, z2 = phase*(1-i*u*r), phase*(r-i*u)
        norm = abs(z1)**2+abs(z2)**2
        a1, b1 = e*e+m*m+t*t, 2*norm
        a2 = (e*m)**2+(e*t)**2+(m*t)**2
        b2 = 2*(m*m*abs(z1)**2+e*e*abs(z2)**2-e*t*mp.re(z1**2)-m*t*mp.re(z2**2))
        c2 = norm**2
        coeff = [c2-s2*b1*b1, b2-2*s2*a1*b1, a2-s2*a1*a1]
        assert all(c > 0 for c in coeff)
        w = -i*u
        assert abs(coeff[1]-(c0+mp.re(c1*w+c2_phase*w*w))) < mp.mpf("1e-44")
        roots = mp.polyroots(coeff,maxsteps=500)
        passed = []
        for x in roots:
            if not isinstance(x,mp.mpf) or x < 0:
                continue
            d = abs(e*m*t-x*(m*z1*z1+e*z2*z2))**2/(a1+b1*x)**3
            mass = sorted(mp.svd_c(mp.matrix([[e,0,mp.sqrt(x)*z1],[0,m,mp.sqrt(x)*z2],
                              [mp.sqrt(x)*z1,mp.sqrt(x)*z2,t]]),compute_uv=False))
            assert abs(d-(mp.fprod(mass))**2/sum(y*y for y in mass)**3) < mp.mpf("1e-44")
            assert abs(d/s3-1) > mp.mpf("1e-4")
            passed.append(dict(x=float(x), determinant_over_target=float(d/s3),
                               mu_over_e=float(mass[1]/mass[0]),tau_over_mu=float(mass[2]/mass[1])))
        row[label] = {"I2_coefficients_descending":list(map(float,coeff)),
                      "I2_roots":[[float(mp.re(x)),float(mp.im(x))] for x in roots],
                      "nonnegative_candidates":passed}
    return {"target_mu_over_e":float(target[1]/target[0]), "target_tau_over_mu":float(1/target[1]),
            "F_and_N_separately_diagonal_effects":True,
            "coherent_sum_cross_effect":str(mixed),
            "continuous_unit_phase_linear_coefficient_lower_bound":float(linear_lower_bound),
            "rows":row}


def stress_test():
    path = repo / "experiments/tfpt-discovery/spin2_interacting_tt_probe.py"
    spec = importlib.util.spec_from_file_location("existing_tt",path)
    tt = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = tt
    spec.loader.exec_module(tt)
    row = {}
    for interaction in tt.INTERACTIONS:
        energy, v = eigh(tt.hamiltonian(interaction).toarray(), driver="evd")
        ground = v[:,0]
        energy -= energy[0]
        for name, momentum_tuple in tt.MOMENTA:
            momentum = np.array(momentum_tuple)
            polarizations = tt.tt_polarizations(momentum)
            stress = [tt.tt_stress_operator(momentum,interaction,p)[0] for p in polarizations]
            amplitudes = v.conj().T @ np.column_stack([op @ ground for op in stress])
            amplitudes[0,:] = 0
            hits = np.flatnonzero(np.linalg.norm(amplitudes,axis=1)>1e-10)
            start = energy[hits[0]]
            shell = np.abs(energy-start)<1e-8
            # Include the entire degenerate eigenspace. No choice of eigenvector
            # within it can change the matrix-valued spectral weight.
            gram = amplitudes[shell,:].conj().T @ amplitudes[shell,:]
            eigenvalues = np.linalg.eigvalsh(gram)
            total = amplitudes.conj().T @ amplitudes
            assert eigenvalues[0] > -1e-12
            change = np.array([[1.,1.],[1.,-1.]])/np.sqrt(2.)
            same_sized_probes = change.T @ gram @ change
            assert abs(same_sized_probes[0,0]-same_sized_probes[1,1]) < 1e-12
            normalized_cross = (same_sized_probes[0,1]/np.sqrt(
                same_sized_probes[0,0].real*same_sized_probes[1,1].real)).real
            assert abs(abs(normalized_cross)-1) < 1e-10
            row[f"{interaction}:{name}"] = dict(first_energy=float(start),
                eigenspace_dim=int(np.sum(shell)), gram=gram.real.tolist(),
                gram_imag_max=float(np.max(abs(gram.imag))), eigs=eigenvalues.tolist(),
                total_gram=total.real.tolist(), rotated_equal_weight_gram=same_sized_probes.real.tolist(),
                rotated_equal_weight_signed_correlation=float(normalized_cross))
    assert [round(row[f"{i}:k1"]["rotated_equal_weight_signed_correlation"])
            for i in tt.INTERACTIONS] == [1,-1,-1]
    return row


def glue_history_test():
    tree = existing("check_joint_glue_clock_stress_20261109.py")
    grades, rotate, roots = [tree[k] for k in ("grades","rotate","roots")]
    odd = [r for r in roots if grades[r] % 2]
    words = []
    for root in odd:
        history = []
        r = root
        for _ in range(6):
            history.append(grades[r])
            r = rotate(r)
        words.append(tuple(history))
    counts = {width:Counter(w[:width] for w in words) for width in range(1,5)}
    for width in range(1,4):
        assert counts[width] == {w:128//2**width for w in product((1,3),repeat=width)}
    assert counts[4] == {w+(w[0],):16 for w in product((1,3),repeat=3)}
    assert sum(n for w,n in counts[4].items() if w[3]==w[0]) == 128
    triples = list(product((1,3),repeat=3))
    full_count = Counter(words)
    history_future = sp.Matrix(8,8,lambda i,j: full_count[triples[i]+triples[j]])
    assert history_future == 16*sp.eye(8) and history_future.rank() == 8
    # Same root support, known quotient energy and time. Positive reweighting
    # can expose a shorter witness, but cannot change the exact return support.
    weighted = {}
    for a in (sp.Rational(1,4), sp.Rational(3,4)):
        d = tree["M"] * sp.Matrix([(j+a)/8-sp.Rational(1,2) for j in range(8)])
        energy = np.array([float(1+d.dot(sp.Matrix(r)/2)) for r in odd])
        weights = np.exp(-energy)
        weights /= weights.sum()
        joint = np.zeros((2,2))
        for w,p in zip(words,weights):
            joint[(w[0]-1)//2,(w[1]-1)//2] += p
        iid = np.outer(joint.sum(1),joint.sum(0))
        tv = np.abs(joint-iid).sum()/2
        assert tv > 1e-6
        weighted[str(a)] = dict(time=1, joint=joint.tolist(),independent_same_marginals=iid.tolist(),
                              total_variation=float(tv))
    odd_glue_minima = {n:sp.Rational(n+3,4) for n in (5,13)}
    assert odd_glue_minima == {5:2,13:4}
    return {"odd_roots":len(odd), "histories_by_length":{k:{"".join(map(str,w)):n for w,n in v.items()} for k,v in counts.items()},
            "first_uniform_history_length_distinguishing_independent":4,"native_return":1.0,
            "independent_return":0.5,"history_future_linear_rank":history_future.rank(),
            "exponential_quotient_energy_reweighting":weighted,
            "odd_glue_minimum_norm_squared_by_carrier_rank":{str(n):int(x) for n,x in odd_glue_minima.items()}}


print(json.dumps({"quarter":phase_test(),"family":family_test(),"tt":stress_test(),
                  "glue_clock":glue_history_test()},indent=2))
