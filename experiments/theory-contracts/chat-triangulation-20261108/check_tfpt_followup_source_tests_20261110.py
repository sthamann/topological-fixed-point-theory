"""Source-anchored control, tomography, convolution and gauge checks."""
import ast
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import block_diag, expm
import sympy as sp

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
sys.path.insert(0,str(repo/"verification"))
from v55_coxeter_cycle import _coxeter, _e8_cartan
from v973_seam_route_narrowing import _e8_simple_roots


def tomography_and_e8():
    source = repo/"experiments/tfpt-discovery/clock_fiber_product_probe.py"
    tree = ast.parse(source.read_text())
    native_defs = [n for n in tree.body if isinstance(n,ast.FunctionDef)
                   and n.name in ("J_all","sigma")]
    assert len(native_defs)==2
    native = {}
    exec(compile(ast.Module(body=native_defs,type_ignores=[]),str(source),"exec"),native)
    basis = sp.Matrix.hstack(*map(sp.Matrix,_e8_simple_roots()))
    c = basis*_coxeter(_e8_cartan())*basis.inv()
    eye = sp.eye(8)
    j,sigma = (sp.Matrix.hstack(*(sp.Matrix(native[name](tuple(eye[:,k])))
                                 for k in range(8))) for name in ("J_all","sigma"))
    x = sp.Symbol("x")
    assert j**2 == -sp.eye(8)
    assert j*sigma == sigma*j
    comm_rank = (j*c-c*j).rank()
    assert comm_rank > 0
    char = c.charpoly(x).as_expr()
    assert char == sp.cyclotomic_poly(30,x)
    assert sp.totient(30)==8 and sp.totient(60)==16
    # Explicit local-hidden state pair on two quaternionic fourers.
    j4 = sp.diag(*[sp.Matrix([[0,1],[-1,0]])]*2)
    odd = sp.diag(1,-1,0,0)
    assert j4*odd*j4.inv() == -odd
    hidden = sp.kronecker_product(odd,odd)
    j16 = sp.kronecker_product(j4,j4)
    assert j16*hidden*j16.inv() == hidden
    rho = [(sp.eye(16)+sign*hidden/2)/16 for sign in (-1,1)]
    assert all(x.trace()==1 and min(x.eigenvals())>0 for x in rho)
    # Every local even effect pairs to zero with odd, so all product
    # effects agree; the common allowed even operator sees +/-1/8.
    readings = [(x*hidden).trace() for x in rho]
    assert readings == [-sp.Rational(1,8),sp.Rational(1,8)]
    return dict(dimensions=[(16**n+(-4)**n)//2 for n in range(1,5)],
                real_five_two_copy=[225,325,625],quaternionic_four_two_copy=[36,136,256],
                explicit_hidden_joint_expectations=list(map(str,readings)),
                J_commutes_native_family_shift=True, J_commutator_with_coxeter_rank=comm_rank,
                coxeter_characteristic_polynomial=str(char),
                rational_complex_structure_commuting_with_coxeter_possible=False)


def convolution():
    c = sp.Symbol("c",real=True)
    assert sp.expand(sp.Rational(4,3)*(1+c)*(1-c)) == 1-(4*c*c-1)/3
    assert sp.Rational(8,3)<sp.log(3)/sp.log(sp.Rational(3,2))
    # f_lambda=1+lambda R_zz, hence the spin-one Fourier eigenvalue is lambda/3.
    lam = sp.Rational(1,3)
    half = sp.sqrt(3*lam)
    fourth = 3*(lam/3)**sp.Rational(1,4)
    assert half==1 and half**2/3==lam
    assert fourth==sp.sqrt(3) and fourth>1
    overlap = (1+c)/2
    mean = lambda expr: sp.integrate(expr,(c,-1,1))/2
    assert mean(overlap)==sp.Rational(1,2)
    assert mean(overlap*(1+lam*c))==sp.Rational(5,9)
    assert mean(overlap*(1+half*c))==sp.Rational(2,3)
    return dict(haar_overlap="1/2",weighted_overlap="5/9",weighted_noncentral=True,
                exact_half_step_overlap="2/3",symmetric_fourth_root_density_amplitude=str(fourth),
                reflection_symmetric_fourth_root_exists=False,
                continuous_identity_start_for_either_density=False,
                casimir_and_all_central_jump_rate_ratio_max=8/3,
                v330_rate_ratio=float(sp.log(3)/sp.log(sp.Rational(3,2))))


def det_center():
    z6 = [(a,b,k) for a in range(3) for b in range(2) for k in range(6)
          if (a-k)%3==0 and (b+k)%2==0]
    assert z6 == sorted((k%3,k%2,k) for k in range(6))
    # If a uniform central generator has phase t in SU(5), then the
    # degree-two wedge reads t² and t^5=1. Its kernel is trivial.
    assert sp.gcd(2,5)==1
    data = json.loads((base/"det16-occupation-symmetry-results-20261110.json").read_text())
    separators = data["density_wedge_comparison"]["separators"]
    partners = {b if a==3 else a for a,b,_ in separators}
    all_others = {x for x in range(32) if x.bit_count()%2==0 and x!=3}
    assert all_others-partners=={0}
    return dict(kernel_size=len(z6), kernel_generator=[1,1,1],
                blind_e_c_partner_bitmask=0, blind_partner="empty wedge, neutral singlet",
                center_preserves_each_occupation=True, conditional_flip_is_not_center=True)


def controlled_cut():
    with contextlib.redirect_stdout(io.StringIO()):
        pair = runpy.run_path(str(base/"check_original_interference_quarter_pair_20261110.py"))
    probe, old = pair["probe"], pair["old"]
    e, n, w = pair["marks"]
    gamma = pair["gamma"]
    q = w @ gamma
    source = old["source"]
    assert np.max(abs(pair["conjugate"](q)+q))<1e-13
    double = lambda x: block_diag(x,x)
    angle = .2  # declared finite response amplitude, not a fitted source value
    one = source["real_action"](expm(1j*angle*probe["weak"]))
    assert np.max(abs(one-expm(1j*angle*q)))<1e-13
    assert np.max(abs(pair["conjugate"](one)-one))<1e-13
    center = [double(z) for z in pair["center"]]
    out = {}
    for m in (-.001,.001):
        k, s = zip(*(pair["half"](m,sgn) for sgn in (1,-1)))
        r, sr = zip(*(pair["half"](m,sgn,reference=True) for sgn in (1,-1)))
        k, r = block_diag(*k), block_diag(*(a.conj().T @ a for a in r))
        e_chain = [(double(e),True),(k,False)]
        skip = [(k,False),(double(n),True)]
        active = [(double(one),False),*skip]
        local = [(double(w),True),*skip]
        qs = [(double(q),True),*skip]
        det_globals = probe["det_mark"].__globals__
        det_globals["identity"] = np.eye(len(k))
        try:
            trace = pair["trace"]
            zref = trace([],[(r,False)],center,2*sum(sr)).real
            shift = 2*sum(s)-np.log(zref)
            empty = trace([(k,False)],[(k,False)],center,shift).real
            gram = probe["block"]([e_chain,skip,active,local,qs],center,shift)
        finally:
            det_globals["identity"] = pair["identity"]
        scale = np.max(abs(gram))
        assert np.max(abs(gram-gram.conj().T))<1e-12*scale
        assert min(np.linalg.eigvalsh(gram))>-1e-12*scale
        plus,minus = (np.array([0,1,sign,0,0])/2 for sign in (1,-1))
        pplus,pminus = [np.vdot(v,gram@v).real for v in (plus,minus)]
        assert pminus>0 and pplus>0
        assert abs(pplus+pminus-(gram[1,1]+gram[2,2]).real/2)<1e-12*scale
        cross = gram[0]@minus
        assert abs(cross.imag)<1e-12*scale
        assert abs(gram[0,3].imag)<1e-12*scale and abs(gram[0,4].real)<1e-12*scale
        pack = lambda z:[float(z.real),float(z.imag)]
        out[str(m)] = dict(reference=zref, empty=empty, skip=float(gram[1,1].real),
                           active=float(gram[2,2].real),plus_port=float(pplus),
                           minus_port=float(pminus),electron_cross_minus_port=pack(cross),
                           electron_cross_clifford_weak=pack(gram[0,3]),
                           electron_cross_bare_weak_generator=pack(gram[0,4]))
    return {"angle":angle,"rows":out,"scope":"fixed R032 pair, exterior and color; no dynamic Gauss or full Padm"}


print(json.dumps(dict(tomography=tomography_and_e8(),haar=convolution(),
                      det_center=det_center(),controlled_cut=controlled_cut()),indent=2))
