"""Run the decisive available checks for the attached F1 proposals."""
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0,str(repo/"verification"))
from tfpt_constants import phi0
from v17_hexagonal_resolvent import U6_power
from v420_fpole_qed_running import ALPHA, ME, MMU, MTAU, MZ


def effect_rank(operations):
    basis = [np.eye(6).reshape(-1)/np.sqrt(6)]
    at = 0
    while at < len(basis):
        effect = basis[at].reshape(6,6)
        for op in operations:
            moved = op.conj().T @ effect @ op
            candidate = ((moved+moved.conj().T)/2).reshape(-1).astype(complex)
            for _ in range(2):
                for q in basis:
                    candidate -= np.vdot(q,candidate).real*q
            norm = np.linalg.norm(candidate)
            if norm > 1e-10:
                basis.append(candidate/norm)
        at += 1
    return len(basis)


def clocks_and_phase():
    eye = np.eye(6)
    shift = np.array(U6_power(1),dtype=float)
    parity = np.linalg.matrix_power(shift,3)
    assert np.array_equal(parity.T,parity) and np.array_equal(parity @ parity,eye)
    # In the given signed family representation U=diag(M0,-M0),
    # its cube is exactly the two-sheet sign.
    q = 1j*parity
    root = (eye+1j*parity)/np.sqrt(2)
    assert np.linalg.norm(root@root-q)<1e-14
    k_plus, k_minus = (eye+q)/2, (eye-q)/2
    assert np.linalg.norm(k_plus-root/np.sqrt(2))<1e-14
    assert np.allclose(k_plus.conj().T@k_plus,eye/2)
    assert np.allclose(k_minus.conj().T@k_minus,eye/2)
    assert np.count_nonzero(parity[0])==1 and parity[0,3]==1
    phase = np.diag([1j,1,1,1,1,1])
    pair = np.diag([0,1,1,0,0,0])
    marks = [np.linalg.matrix_power(shift,i) @ pair @ np.linalg.matrix_power(shift,-i)
             for i in range(6)]
    base = [shift,*marks,*(eye-p for p in marks)]
    delta = ((794-7*np.sqrt(9961))/2187)**(1/6)
    a, b = 6*np.log(1.5), 6*np.log(3)
    rows = []
    for y in (1.,1/3,2/3):
        r = y*delta/(y-delta)**2
        green = np.linalg.inv((1+2*r)*eye-r*(shift+shift.T))
        ranks = [effect_rank(base), effect_rank(base+[green]), effect_rank(base+[green,q]),
                 effect_rank(base+[green,root]),effect_rank(base+[green,phase])]
        assert ranks == [6,21,21,36,36], ranks
        fastest = a*np.log1p(4*r)/np.log1p(r)
        rows.append(dict(cusp=y,ranks=ranks,fast_over_slow=float(fastest/a),
                         passive_return_lower_bound=float(np.exp(-fastest))))
    # Each real +/-1 eigenspace of the directed shift is a real line.
    # Any commuting real J preserves each line and cannot square to -1.
    assert np.linalg.matrix_rank(eye-shift)==5
    assert np.linalg.matrix_rank(eye+shift)==5
    plus = np.array([1,1j])/np.sqrt(2)
    minus = plus.conj()
    measure = np.ones((2,2))/2
    def realify(matrix):
        return np.block([[matrix.real,-matrix.imag],[matrix.imag,matrix.real]])
    def probs(operator):
        out = []
        for ket in (plus,minus):
            ket = operator @ ket
            r = np.concatenate([ket.real,ket.imag])
            measured = np.vdot(ket,measure@ket).real
            assert abs(measured-r @ realify(measure) @ r)<1e-14
            out.append(float(measured))
        return out
    assert np.allclose(probs(np.eye(2)),[.5,.5])
    assert np.allclose(probs(np.diag([1,-1j])),[1,0])
    two_flavor = np.diag([1,1j]) @ np.array([[1.,1.],[-1.,1.]])/np.sqrt(2)
    a11,a12,a21,a22 = two_flavor.ravel()
    quartic_invariant = (a11*a22*a12.conjugate()*a21.conjugate()).imag
    assert abs(quartic_invariant)<1e-14
    return dict(ranks_order=["bare","green","green_and_iU3","green_and_extra_square_root",
                            "green_and_site_local_i"],cusps=rows, v330_rate_ratio=float(b/a),
                v330_target_return=float(np.exp(-b)),
                shortest_lift_spring_ratio=float(3*np.log(2)/np.log(9/8)),
                coherent_skip_or_quarter_branch_probability=.5,
                real_shift_one_dimensional_fixed_and_alternating=True,
                two_flavor_rephasing_invariant_imaginary_part=float(quartic_invariant),
                same_quarter_probabilities_in_complex_2_and_real_4=probs(np.diag([1,-1j])))


def leptons():
    mp.mp.dps = 45
    source = (mp.mpf(16)/7*phi0**5,mp.mpf(4)/3*phi0**3,mp.mpf(7)/6*phi0**2)
    pole = (ME,MMU,MTAU)
    coords = lambda m: (mp.log(m[1]/m[0]),mp.log(m[2]/m[1]))
    k = lambda m: m[0]*m[2]**2/m[1]**3
    assert abs(k(source)-mp.mpf(21)/16) < mp.mpf("1e-37")
    da,db = [x-y for x,y in zip(coords(pole),coords(source))]
    assert da>0 and db>da
    msbar = tuple(m*(1-(ALPHA/mp.pi)*(1+mp.mpf(3)/4*mp.log(MZ**2/m**2))) for m in pole)
    source_delta = mp.log(k(pole)/k(source))
    one_loop_delta = mp.log(k(pole)/k(msbar))
    # Exact linearized log conversion. Common factor drops from K.
    charges = (1,-3,2)
    rel = [-ALPHA/mp.pi*(1+mp.mpf(3)/4*mp.log(MZ**2/m**2)) for m in pole]
    linear = -sum(q*v for q,v in zip(charges,rel))
    assert abs(linear+3*ALPHA/(2*mp.pi)*mp.log(k(pole))) < mp.mpf("1e-38")
    delta = ((794-7*mp.sqrt(9961))/2187)**(mp.mpf(1)/6)
    a, b, c = 6*mp.log(mp.mpf(3)/2), 6*mp.log(3), 12*mp.log(2)
    # v120 original addresses e+, mu-, tau-. Test the proposed *multiplicative*
    # use of the actual shortest-lift Fuchs two-sheet response on them.
    values = {}
    for y in (mp.mpf(1),mp.mpf(1)/3,mp.mpf(2)/3):
        r = y*delta/(y-delta)**2
        bn,cn = [a*mp.log1p(n*r)/mp.log1p(r) for n in (3,4)]
        plus = (1-mp.exp(-b))/(1-mp.exp(-bn))
        minus = (mp.exp(-c)-mp.exp(-a))/(mp.exp(-cn)-mp.exp(-a))
        aa,bb = coords((plus,minus,minus))
        assert aa>0 and bb==0
        values[str(float(y))] = {"M_plus":float(plus),"M_minus":float(minus),
                                  "delta_a":float(aa),"delta_b":float(bb),
                                  "delta_log_K":float(-aa)}
    return dict(K_source=float(k(source)),K_pole=float(k(pole)), delta_log_K=float(source_delta),
                fractional_K_gap=float(mp.expm1(source_delta)), required_delta_a=float(da),
                required_delta_b=float(db), one_loop_pole_minus_msbar_linear=float(linear),
                v420_unexpanded_pole_minus_msbar=float(one_loop_delta),
                fuchs_diagonal_sheet_hypothesis=values)


def haar_and_higgs():
    # The squared spinor overlap with an SU(2)-Haar transported unit
    # doublet is uniform on [0,1]. x also reads the normalized squared
    # weak endpoint vertex on the transverse chart.
    x,u,e = sp.symbols("x u epsilon",positive=True)
    v = sp.pi**2*sp.integrate(u,(u,0,e))
    v_mark = sp.pi**2*sp.integrate(u*u,(u,0,e))
    # Two positive normalized reflection-compatible weight kernels.
    weights = [sp.Integer(1), (1+x)/sp.Rational(3,2)]
    moments = []
    for weight in weights:
        assert sp.integrate(weight,(x,0,1))==1
        moments.append(sp.integrate(weight*x,(x,0,1)))
    assert moments == [sp.Rational(1,2),sp.Rational(5,9)]
    # A split link makes the overlap x=tr(P_R' P_L'), a Gram kernel
    # on Hermitian projector matrices. Adding a constant feature
    # preserves reflection positivity for these scalar test kernels.
    return dict(flat_doublet_radial_jacobian="pi^2 u du, normalized angular Haar",
                one_site_ball=str(v), same_site_two_yukawa_scale=str(v_mark),
                scalar_Haar_overlap_moments=list(map(str,moments)),
                same_normalization_both=True, no_physical_weight_assigned=True)


print(json.dumps(dict(clock=clocks_and_phase(),leptons=leptons(),
                      haar_radial=haar_and_higgs()),indent=2))
