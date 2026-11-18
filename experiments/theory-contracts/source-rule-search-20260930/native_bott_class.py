"""Explicit normal-spin Bott class and actual E8 root trace normalization.

Computes a topological class. It does not assert that the raw source chooses
the corresponding operator family, chiral representation, or physical time.
"""
import json
from pathlib import Path
import sys
import sympy as s
from sympy.simplify.fu import fu
sys.path.insert(0,'/Users/stefanhamann/Projekte/tfpt-theoryv4')
from tfpt_explorer.current_block_geometry import _actual_current_source, _kappa

theta,beta,phi=s.symbols('theta beta phi',real=True)
a=s.cos(theta/2); b=s.sin(theta/2)
nx=s.sin(beta)*s.cos(phi); ny=s.sin(beta)*s.sin(phi); nz=s.cos(beta)
# U=(a I+i b n.sigma)(a I-i b sigma3) is periodic in theta.
q=s.Matrix([a*a+b*b*nz, a*b*nx+b*b*ny,
            a*b*ny-b*b*nx, a*b*(nz-1)])
assert fu(s.trigsimp(s.expand((q.T*q)[0]-1)))==0
assert q.subs(theta,0)==s.Matrix([1,0,0,0])
assert s.simplify(q.subs(beta,0))==s.Matrix([1,0,0,0])
assert s.simplify(q.subs(theta,2*s.pi))==s.Matrix([1,0,0,0])
# Pullback of outward volume on the unit SU2=S3.
jacobian=s.trigsimp(s.Matrix.hstack(q,q.diff(theta),q.diff(beta),q.diff(phi)).det())
expected=(1-s.cos(beta))*s.sin(beta)*s.sin(theta/2)**2/2
assert s.trigsimp(s.expand_trig(jacobian-expected))==0
volume=s.integrate(expected,(theta,0,2*s.pi),(beta,0,s.pi),(phi,0,2*s.pi))
assert volume==2*s.pi**2
degree=s.simplify(volume/(2*s.pi**2))
assert degree==1

# Full E8 adjoint restriction; NOT a fictitious 2-dimensional E8 module.
src=_actual_current_source(); alg=src['chevalley']
root=(-1,)*8; lam=s.Matrix(root)/2; H=alg.hvec_index(root)
assert _kappa(src,H,H)==2
adj_trace=sum(((lam.T*(s.Matrix(r)/2))[0])**2 for r in alg.roots)
assert adj_trace==120
fundamental_su2_trace=2
adj_index=adj_trace/fundamental_su2_trace
assert adj_index==60

result={
 'verdict':'EXACT_BOTT_TO_BASIC_E8_CLASS; ACTUAL_SOURCE_OPERATOR_IDENTIFICATION_OPEN',
 'P1_inputs':['[u_Sigma]=+1 in K1(S1)','oriented normal two-slice compactified to S2'],
 'normal_spin_structure':'unique: w2=e mod2=0 and H1(S2,Z2)=0',
 'spin_line_class':'c1(S+)=-1 for S+=K^(1/2); dual chirality has +1',
 'class_product':'[u_Sigma] external_product ([S^chirality]-[1]) generates reduced K1(S1 smash S2)=K1(S3)',
 'representative':'U(theta,n)=[I+(exp(i theta)-1)P(n)] [I+(exp(i theta)-1)P(n0)]^-1, P=(I+n.sigma)/2',
 'SU2_geometric_degree':str(degree),
 'volume_density':str(expected),
 'volume_integral':str(volume),
 'sign_guard':'Geometric degree shown with q0 I+i q.sigma and outward orientations. Chern-character/clutch conventions and normal chirality must be fixed separately; no positive physical level claimed.',
 'native_root':list(map(str,lam)),
 'basic_form_root_embedding_index':'1: kappa_E8(H,H)=2=kappa_SU2(H,H)',
 'clutching_result':'absolute basic E8 c2=1 on S4 after long-root SU2 embedding',
 'actual_adjoint_trace_H_squared':str(adj_trace),
 'actual_adjoint_index_relative_to_SU2_doublet':str(adj_index),
 'operator_guard':'A finite E8-adjoint Dirac determinant yields 60 times the basic trace normalization. The unit Bott class cannot be silently called a unit adjoint-family anomaly.',
 'normal_operator_guard':'Untwisted compact S2 spin Dirac has index0. Spin Dirac twisted by O(k) has indexk. The normal Thom/Bott class is a separate symbol class, not an already selected confined operator.',
 'clock_guard':'theta->theta+pi/2 moves the based collapsed slice. U(i z,n)=u(i,n)U(z,n)u(i,n0)^-1; this is not constant native clock conjugation.',
 'not_derived':['raw P1 response realizes this external product','Quillen/current level and vacuum selection','physical clock equivariance','per-mark c3 coupling','fundamental full solution','prime law','AI superiority']
}
Path(__file__).with_name('native_bott_class_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
