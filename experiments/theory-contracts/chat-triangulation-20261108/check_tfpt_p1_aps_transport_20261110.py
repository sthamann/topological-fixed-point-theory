"""Exact transport constraints on TFPT's original R074 two-face graph."""
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
import sympy as sp

repo = Path(sys.argv[1]).resolve()
base = repo/'experiments/theory-contracts/chat-triangulation-20261108'
with contextlib.redirect_stdout(io.StringIO()) as buffer:
    runpy.run_path(str(base/'check_original_qtheta_word_20261109.py'))
physical = json.loads(buffer.getvalue())['native_endpoint_call']

# For a fixed Fourier/channel transition, T maps (1,t) into (1,s*t).
a,b,c,d,t,s = sp.symbols('a b c d t s', real=True)
graph_constraint = sp.Poly(c+d*t-s*t*(a+b*t),t)
coefficients = graph_constraint.all_coeffs()
assert coefficients == [-b*s, d-a*s, c]
assert sp.solve(coefficients, (b,c,d)) == {b:0,c:0,d:a*s}
# The real scalar calculation represents the coefficients of a complex block.
# In the adjoint constraint replace a with its conjugate; s is strictly positive.
reverse_adjoint = sp.expand((a*s)*s*t-a*t)
assert sp.factor(reverse_adjoint) == a*t*(s-1)*(s+1)

alpha,beta = sp.Rational(1,4),sp.Rational(3,4)
# Every 8-fold transition carries an integer shift plus a quarter offset.
for left,right in ((alpha,beta),(beta,alpha),(alpha,alpha),(beta,beta)):
    for m in range(-3,4):
        for l in range(8):
            for k in range(8):
                delta=8*m+k-l+right-left
                assert (delta == 0) == (left == right and m == 0 and k == l)

def face(p):
    slope=np.exp(-2*np.pi*p)
    g=np.array([1.,slope])/np.hypot(1.,slope)
    return np.stack((g, np.array([-g[1],g[0]])),axis=1)

def cross_test(n):
    left,right=n+float(alpha),n+float(beta)
    assert (left<0)==(right<0)
    fa,fb=face(left),face(right)
    ga,gb=np.outer(fa[:,0],fa[:,0]),np.outer(fb[:,0],fb[:,0])
    d_forward=np.diag([1.,np.exp(-np.pi)])
    u=fb@fa.T
    adjoint_leak=np.linalg.norm((np.eye(2)-ga)@d_forward.T@fb[:,0])
    assert np.linalg.norm((np.eye(2)-gb)@d_forward@fa[:,0])<1e-14
    assert np.linalg.norm(u@ga@u.T-gb)<1e-14
    assert np.linalg.norm(u.T@u-np.eye(2))<1e-14
    assert abs(np.linalg.det(u)-1)<1e-14
    assert np.linalg.norm(right*u-u*left-u/2)<1e-14
    return {
        'mode':n,'tangential_source':left,'tangential_target':right,
        'same_negative_APS_sign_in_this_chart':True,
        'one_sided_multiplier_adjoint_leakage':adjoint_leak,
        'transport_angle_rad':float(np.arctan2(u[1,0],u[0,0])),
        'absolute_energy_offset':right-left,
    }

swap=np.array([[0.,1.],[1.,0.]])
for n in range(-30,31):
    pa,pb=n+float(alpha),n+float(beta)
    # Genuine paired reflection of the fixed scalar collar: reflect the
    # tangential momentum and exchange both faces.
    forward=face(pb)@face(pa).T
    reflected_reverse=face(-pb)@face(-pa).T
    assert np.linalg.norm(swap@forward@swap-reflected_reverse)<1e-14

fold=sp.Matrix(json.loads((repo/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())['explicit_Weyl_matrix'])
unit=sp.eye(8)
q=fold.T*(unit[:,7]-unit[:,6])
w=sp.Matrix([0,1,1,-1,0,0,1,0])
r=w-q
bottom=r-q
half=sp.ones(8,1)/2
def lattice_sector(momentum):
    if all(x.q==1 for x in momentum):
        return '0' if sum(momentum)%2==0 else 'v'
    offset=momentum-half
    assert all(x.q==1 for x in offset)
    return 's' if sum(offset)%2==0 else 'c'
sectors={key:lattice_sector(p) for key,p in {'q':q,'2q':2*q,'B':bottom,'r':r,'W':w}.items()}
assert sectors=={'q':'s','2q':'0','B':'0','r':'s','W':'0'}
assert (fold*w)[0]==(fold*w)[1]==(fold*w)[2]==0

rows=[cross_test(n) for n in range(-30,31)]
assert rows[30]['one_sided_multiplier_adjoint_leakage']>.2
assert abs(rows[0]['transport_angle_rad'])<1e-14 and abs(rows[-1]['transport_angle_rad'])<1e-14
assert abs(rows[30]['transport_angle_rad'])>.1
assert min(abs(n+float(alpha)) for n in range(-30,31))==.25
assert min(abs(n+float(beta)) for n in range(-30,31))==.25
# Test the concurrent original same-loop paired Pin candidate on its actual
# native 96 coefficients, then on its ordered product over three backgrounds.
with contextlib.redirect_stdout(io.StringIO()):
    native=runpy.run_path(str(base/'check_tfpt_conditioned_r113_weights_20261110.py'))
tau2=sp.Matrix([[0,-sp.I],[sp.I,0]])
z1,z2,z3,z4=sp.symbols('z1 z2 z3 z4')
x=sp.Matrix([[z1,z2],[z3,z4]])
solutions=sp.solve(list(x*tau2-tau2*x),(z3,z4))
assert solutions=={z3:-z2,z4:z1}
p0,p1=sp.diag(1,0),sp.diag(0,1)
assert (p0*sp.I*tau2*p1)**2==sp.zeros(2)
t2=np.array(tau2,dtype=complex)
identity=np.eye(96)
pin_j=np.kron(1j*t2,identity)
b_side=np.kron(np.diag([1,-1]),identity)
coeff=[np.kron(-t2,native['inner'](np.array(h)))
       for h in ((0,1),(0,1j),(1,0),(1j,0))]
generators=[np.kron(-t2,native['majorana']+native['inner'](np.array(h)))
            for h in ((0,.55),(.12j,.55),(.08,.44j))]
assert all(np.linalg.norm(a@pin_j-pin_j@a)==0 for a in coeff+generators)
assert np.linalg.norm(b_side@generators[0]-generators[0]@b_side)>1
assert np.linalg.norm(pin_j.conj().T@pin_j-np.eye(192))==0
assert np.linalg.norm(generators[0]@generators[1]-generators[1]@generators[0])>1e-3
history=np.eye(192,dtype=complex)
for k in generators:
    val,vec=np.linalg.eigh(k)
    history=(vec*np.exp(-val/20))@vec.conj().T@history
B=coeff[0]
early=history@B@pin_j@history
late=history@pin_j@B@history
defect=np.linalg.norm(early-late)
assert defect<1e-10

print(json.dumps({
    'scope':'original fixed R074 complex CAR product collar in periodic gauge; not the full physical APS domain',
    'graph_polynomial_coefficients':list(map(str,coefficients)),
    'reverse_adjoint_defect':str(reverse_adjoint),
    'one_sided_multiplier_far_face_factor_exact':'exp(-pi)',
    'comparison_scalar_tangential_gap_exact':'1/4',
    'paired_reflection_covariance_of_modewise_graph_map':True,
    'old_P1_lattice_sectors':sectors,
    'same_loop_native_Pin_label_test':{
        'label_commutant':'a I + b i tau2',
        'commuting_192d_pin_J_unitary':True,
        'whole_order_contact_defect_HS':float(defect),
        'b_side_does_not_preserve_K':True,
        'one_way_label_raising_square_zero':True,
        'geometric_pair_swap_cannot_realize_original_nilpotent_family_current':True,
    },
    'vacuum_normalized_selector_scalar_on_D8_sectors_keeps_W':True,
    'fixed_quarter_local_adjoint_preserving_bridge':False,
    'sample_modes':[r for r in rows if r['mode'] in (-2,-1,0,1,2)],
    'direct_physical_source_caller':physical,
    'P1_determinant_selection':None,
    'P1_gap_selection':None,
    'P1_original_weighted_early':None,
    'P1_original_weighted_late':None,
},indent=2))
