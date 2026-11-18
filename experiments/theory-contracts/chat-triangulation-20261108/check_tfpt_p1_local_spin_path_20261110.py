"""Exact P1 current modes and the source caller, reusing original TFPT routines."""
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

repo = Path(sys.argv[1]).resolve()
base = repo / 'experiments/theory-contracts/chat-triangulation-20261108'
with contextlib.redirect_stdout(io.StringIO()):
    original = runpy.run_path(str(base / 'check_tfpt_family_middle_ordered_cut_20261110.py'))

sp, Q, q, r, d1, M = (original[k] for k in ('sp','Q','q','r','d1','M'))
vertex, middle, inner, grade = (original[k] for k in ('vertex','middle','inner','grade'))
old, axes = (original[k] for k in ('old','V'))
t, delay = original['t'], original['delta_t']

def scale(words, factor):
    return {mon: factor*value for mon,value in words.items() if factor*value}

def diff(left,right):
    return original['combine']([*left.items(),*((w,-x) for w,x in right.items())])

def norm(w):
    return sp.simplify(inner(w,w))

def center(p):
    target=M*p
    return sp.simplify(target[0]+target[1]-2*target[2])

out={}
for alpha in map(sp.Rational,('1/4','3/4')):
    drift=sp.Matrix([j+alpha-4 for j in range(8)])
    coeff=[Q.from_sympy(z) for z in old['cl'](axes.T*drift/8)]
    v=old['vertex_zero'](old['series_response'](original['polys'][1],coeff,4))
    B=vertex(v,d1,q)
    W=middle(v)
    assert q.dot(r)==0 and drift.dot(q)==0 and all(center(p)==0 for p in (q,2*q,r,r+q))
    assert norm(v)>0 and norm(W)>0
    assert all(grade(word)==2 for word in W)
    assert not vertex(v,q,q)
    assert vertex(B,q,r) in (W,scale(W,-Q.one))
    assert not vertex(W,q,r+q)
    BW=vertex(W,-q,r+q)
    back=vertex(BW,q,r)
    assert back in (scale(W,2*Q.one),scale(W,-2*Q.one))
    # The adjoint sign follows from the positive old Hilbert metric.
    assert norm(BW)==2*norm(W)
    sign=1 if inner(BW,B)==norm(W) else -1
    assert inner(BW,B)==sign*norm(W)
    trip=scale(BW,Q.from_sympy(sp.Rational(sign,2)))
    sing=diff(B,trip)
    assert not vertex(sing,q,r) and not vertex(sing,-q,r)
    assert inner(trip,sing)==0 and norm(trip)+norm(sing)==norm(B)
    E=32+drift.dot(r)
    modes={}
    for k in (-2,-1,0,1,2,3):
        ev=vertex(v,q,q,k)
        early=vertex(ev,d1,2*q)
        late=vertex(B,q,r,k)
        # All output words have momentum r+q and one fixed oscillator grade.
        assert all(grade(w)==2+k for w in (*early,*late))
        assert all(grade(w)==k for w in ev)
        if k:
            assert inner(W,early)==0 and inner(W,late)==0
        c_early=norm(early)/norm(v)
        c_late=norm(late)/norm(v)
        late_factor=sp.exp(-t*(E+8*k)-delay*E)
        early_factor=sp.exp(-t*(E+8*k)-delay*(32+8*k))
        # Old two-face Cauchy graph: use the already fixed window length as
        # the strip width in this separate check. The actual APS bridge is
        # not inferred from this old graph.
        strip=delay
        source_slope=sp.exp(-strip*E)
        target_slope=sp.exp(-strip*(E+8*k))
        mirror_multiplier=sp.exp(-8*k*strip)
        assert sp.simplify(mirror_multiplier*source_slope-target_slope)==0
        leakage2=(source_slope-target_slope)**2/(
            (1+source_slope**2)*(1+target_slope**2))*c_late
        assert (leakage2 == 0) == (k==0 or not late)
        # Direct application of the old orthogonal graph projector to the
        # diagonally inserted normalized input reproduces the leakage.
        v_g=sp.Matrix([1,source_slope])/sp.sqrt(1+source_slope**2)
        end_g=sp.Matrix([1,target_slope])
        p_g=(end_g*end_g.T)/(1+target_slope**2)
        assert sp.simplify(((sp.eye(2)-p_g)*v_g).dot((sp.eye(2)-p_g)*v_g)*c_late-leakage2)==0
        modes[str(k)]={
            'current_mode_index':-k,
            'old_cauchy_far_face_multiplier':str(mirror_multiplier),
            'old_cauchy_identical_face_insertion_leakage2_over_q_norm':float(leakage2.evalf(18)),
            'early_vector_zero':not early,
            'late_vector_zero':not late,
            'early_norm2_over_input':str(c_early),
            'late_norm2_over_input':str(c_late),
            'early_reflected_return_old':float((c_early*early_factor**2).evalf(18)),
            'late_reflected_return_old':float((c_late*late_factor**2).evalf(18)),
            'early_fixed_W0_mixed_pairing_zero':inner(W,early)==0,
            'late_fixed_W0_mixed_pairing_zero':inner(W,late)==0,
        }
    out[str(alpha)]={
        'reference_norm':str(norm(v)),
        'middle_momentum':list(map(str,r+q)),
        'middle_energy':str(E),
        'B_norm2':str(norm(B)),
        'family_triplet_weight_zero_norm2':str(norm(trip)),
        'family_singlet_weight_zero_norm2':str(norm(sing)),
        'middle_W_norm2':str(norm(W)),
        'P_sing_on_all_old_tested_paths':1,
        'modes':modes,
    }

# Invoke both existing endpoints on the installed, unchanged original domain.
sys.path.insert(0,str(repo/'experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime'))
with contextlib.redirect_stdout(io.StringIO()) as buffer:
    probe=runpy.run_path(str(base/'check_original_qtheta_word_20261109.py'))
caller=json.loads(buffer.getvalue())['native_endpoint_call']
print(json.dumps({'old_p1_exact':out,'original_physical_caller':caller,
                  'physical_early':None,'physical_late':None},indent=2))
