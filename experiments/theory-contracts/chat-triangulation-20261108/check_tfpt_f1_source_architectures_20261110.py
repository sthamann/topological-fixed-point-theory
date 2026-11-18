"""Native tests for the R113 controller, finite-weight uniqueness and TT constraints."""
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import mpmath as mp
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
import sympy as sp

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0, str(repo/"verification"))
from tfpt_constants import Mbar
from v86_nstar_reheating import nstar_chain, starobinsky
from v494_cosmo_running_mu_record import _observables, _phi_of_N, _phi_end
from v1_e8_glue import e8_roots


def controller_and_link_representations():
    with contextlib.redirect_stdout(io.StringIO()):
        old = runpy.run_path(str(base/"check_original_interference_quarter_pair_20261110.py"))
    probe, source = old["probe"], old["probe"]["src"]
    native, swap = probe["native"], probe["swap"]
    weak, gamma = old["marks"][2], old["gamma"]
    right = [native.IDX[("nu","R",0,g)] for g in range(native.NG)]
    right += [i+native.NP for i in right]
    remaining = sorted(set(range(len(gamma)))-set(right))
    ray_blocks = []
    for h in old["old"]["history"][:3]:
        # Clifford Higgs contact between the neutral right singlets and
        # the other modes. No synthetic family operator is used.
        h_y = -1j*gamma@source["inner_form"](h)
        block = h_y[np.ix_(remaining,right)]
        sing = np.linalg.svd(block,compute_uv=False)
        assert len(sing)==6 and min(sing)>1e-8
        ray_blocks.append(float(min(sing)))
    eye = np.eye(len(gamma))
    target = eye[:,native.IDX[("nu","L",0,0)]]
    spectator0,spectator1 = eye[:,right[0]],eye[:,right[1]]
    actual_link = source["real_action"](expm(.2j*probe["weak"]))
    assert np.max(abs(actual_link-expm(.2j*weak@gamma)))<1e-13
    assert np.allclose(actual_link@spectator0,spectator0)
    assert np.allclose(actual_link@spectator1,spectator1)
    target_after = actual_link@target
    overlap = np.vdot(target,target_after)
    complementary = target_after-overlap*target
    coupling = np.linalg.norm(complementary)
    assert coupling>1e-6
    # Input is a decomposable two-particle wedge, with one right neutrino.
    # Controlled output is (r0 wedge target + r1 wedge V target)/sqrt(2).
    # In the orthogonal two-by-two support, the minor below tests Slater rank.
    control_output = np.array([[1.,0.],[overlap,coupling]])/np.sqrt(2.)
    plucker = abs(np.linalg.det(control_output))
    assert abs(plucker-coupling/2)<1e-14
    eigs = np.linalg.eigvalsh(control_output@control_output.conj().T)
    assert np.all(eigs>1e-6)
    # Native weak link representation: genuine doublet counting including
    # antiparticles. Each doublet contributes at most spin 1/2 on exterior.
    rep = native.rep(0, .5j*np.diag([1.,-1.]), np.zeros((3,3)))
    cartan = -1j*(rep+swap@rep.conj()@swap)
    values = np.linalg.eigvalsh(cartan)
    hist = {str(x):int(np.count_nonzero(values==x)) for x in (-.5,0.,.5)}
    assert hist=={"-0.5":24,"0.0":48,"0.5":24}
    return dict(neutral_right_singlet_higgs_min_singular_values=ray_blocks,
                gaussian_wedge_input_slater_rank=1,controlled_wedge_output_slater_rank=2,
                controlled_output_plucker_minor_abs=float(plucker),
                controlled_output_rank_eigenvalues=eigs.tolist(),
                one_body_weak_weights=hist,max_weak_spin_open_exterior=12,
                max_weak_spin_reflected_product=24)


def identical_low_spin_semigroups():
    # Positive central jump laws (1 ± chi_25/102)dHaar, with a common
    # diffusion a=1/1000. |chi_25|<=51. Only spin 25 distinguishes them.
    j0 = sp.Integer(25)
    eps = sp.Rational(1,102)
    assert eps*(2*j0+1)==sp.Rational(1,2)
    a = sp.Rational(1,1000)
    def generator(j, sign):
        jump_fourier = 1 if j==0 else sign*eps/(2*j+1) if j==j0 else 0
        return a*j*(j+1) + 1-jump_fourier
    for n in range(49):
        assert generator(sp.Rational(n,2),1)==generator(sp.Rational(n,2),-1)
    r_plus,r_minus = (generator(j0,sign) for sign in (1,-1))
    assert r_minus-r_plus == sp.Rational(1,2601)
    t,s = sp.symbols("t s",nonnegative=True)
    assert sp.expand(r_plus*(t+s)-r_plus*t-r_plus*s)==0
    c = sp.Symbol("cos_theta",real=True)
    assert sp.cancel((1-(4*c*c-1)/3)/(1-c))==sp.Rational(4,3)*(1+c)
    assert sp.cancel((1-(2*c*c-1))/(1-c))==2*(1+c)
    target = float(sp.log(3)/sp.log(sp.Rational(3,2)))
    # With any common lambda(g) and a(g), all measured jets through
    # spin 24 coincide identically, since the difference is in spin 25.
    return dict(difference_first_spin=25,matched_spins=[0,24],
                spin25_rates=list(map(str,(r_plus,r_minus))),
                spin25_positive_transfer_at_time1=[float(sp.exp(-r)) for r in (r_plus,r_minus)],
                density_lower_bound_relative_to_haar="1/2 for each jump law",
                u1_symmetric_jump_angle_for_target_ratio=float(np.arccos(target/2-1)),
                matched_1_and_2_metric_jets_with_same_coefficient_families=True)


def glue_reference():
    roots = [tuple(int(x) for x in 2*r) for r in e8_roots()]
    # Original D5/A3 integral charge in the fixed coordinate split.
    def grade(v):
        if v[0]%2==0:
            return 2*(sum(x//2 for x in v)%2)
        return 1+2*(sum((x-1)//2 for x in v)%2)
    assert all(grade(r[:5])==grade(r[5:]) for r in roots)
    d = [sum(grade(r[:5])==k for r in roots) for k in range(4)]
    assert d == [52,64,60,64]
    max_shift_rank = lambda sizes:sum(min(sizes[k],sizes[(k+1)%4]) for k in range(4))
    assert max_shift_rank(d)==224 and max_shift_rank([d[0]+8,*d[1:]])==240
    alpha = (1,)*8
    assert alpha in roots and grade(alpha[:5])==1
    strings = sum(sum(a*b for a,b in zip(alpha,r)) == -4 for r in roots)
    assert strings==56
    # Four grade reference. Diagonal charge states give zero trace of a shift,
    # while the coherently prepared phase eigenstate gives unit visibility.
    shift=sp.zeros(4)
    for k in range(4):
        shift[(k+1)%4,k]=1
    plus = sp.ones(4,1)/2
    assert (plus.T*shift*plus)[0]==1 and shift.trace()==0
    return dict(root_grade_dimensions=d,max_rank_of_grade_plus_one_on_240=224,
                max_rank_of_grade_plus_one_on_248=240,
                odd_root_chevalley_adjoint_rank=strings+2,
                odd_root_chevalley_nilpotence_index=3,
                prepared_four_grade_reference_visibility=1,
                charge_diagonal_four_grade_reference_visibility=0)


def geometry_and_inflation():
    # Necessary positivity and masslessness on one fixed TT component.
    m,b,g,u,v = sp.symbols("m b g u v",real=True)
    matrix = sp.Matrix([[m,b],[b,g]])
    assert sp.expand((sp.Matrix([u,v]).T*matrix*sp.Matrix([u,v]))[0]
                     - m*(u+b*v/m)**2-(g-b*b/m)*v*v)==0
    assert matrix.det().subs(g,0)==-b*b
    assert matrix.det().subs(g,b*b/m)==0
    # Heat a2 (R coefficient before the common scalar 1/2 determinant factor).
    scalar = sp.Rational(1,6)
    weyl = -sp.Rational(1,2)*4*(sp.Rational(1,6)-sp.Rational(1,4))
    vector = 4*scalar-1-2*scalar
    assert weyl/scalar==1 and vector/scalar==-4
    conformal = scalar-sp.Rational(1,6)
    assert 48*weyl+12*vector+4*conformal==0
    # The actual flat R048 Higgs exchange is exactly zero in the TT channel:
    # linear scalar curvature is k_i k_j h_ij - k² tr(h).
    q0,q3,A,B = sp.symbols("q0 q3 A B",real=True)
    momentum = sp.Matrix([q0,0,0,q3])
    h = sp.Matrix([[0,0,0,0],[0,A,B,0],[0,B,-A,0],[0,0,0,0]])
    assert sp.simplify((momentum.T*h*momentum)[0]-(momentum.dot(momentum))*h.trace())==0

    V, _, f_end = starobinsky()
    rho_end = mp.mpf(4)/3*V(f_end)*Mbar**4
    cap, _, _ = nstar_chain(4,"0.05",rho_reh_override=rho_end)
    tilt = lambda n: _observables(_phi_of_N(float(n),_phi_end()))[0]
    ns_cap_lo, ns_cap_exact = float(1-2/cap), float(tilt(cap))
    near = .9752-2*.003
    need_lo = 2/(1-near)
    need_ex = brentq(lambda n:tilt(n)-near,55,80)
    assert max(ns_cap_lo,ns_cap_exact)<near
    assert float(cap)<min(need_lo,need_ex)
    data = json.loads((base/"replica-reheating-results-20261110.json").read_text())
    table = data["reheating"] if "reheating" in data else data
    def collect(node):
        if isinstance(node,dict):
            if "n_s_exact_potential_first_slow_roll" in node and abs(node.get("N_star",0)-50.16942)<.001:
                return node
            for value in node.values():
                out = collect(value)
                if out:
                    return out
        return None
    anomaly = collect(table)
    assert anomaly
    pulls = lambda x:(.9752-x)/.003
    return dict(TT_R048_linear_scalar_curvature=0,
                conventional_SM_common_regulator_weights_per_species=[1,1,-4],
                SM_48_Weyl_12_vectors_conformal_scalar_leading_sum=0,
                positive_coupled_zero_mode_requires="g(0)=b(0)^2/m(0), with m(0)>0",
                determinant_if_bare_metric_zero_and_nonzero_constant_mixing="-b^2",
                reheating_instantaneous_N_cap=float(cap),
                ns_cap_leading=ns_cap_lo,ns_cap_exact_potential=ns_cap_exact,
                best_DR2_pull_leading=pulls(ns_cap_lo),best_DR2_pull_exact_potential=pulls(ns_cap_exact),
                N_required_for_DR2_2sigma_edge=[need_lo,need_ex],
                conformal_DR2_pull_leading=pulls(anomaly["n_s_LO"]),
                conformal_DR2_pull_exact_potential=pulls(anomaly["n_s_exact_potential_first_slow_roll"]))


print(json.dumps(dict(control=controller_and_link_representations(),
                      weight=identical_low_spin_semigroups(),
                      glue_reference=glue_reference(),
                      gravity=geometry_and_inflation()),indent=2))
