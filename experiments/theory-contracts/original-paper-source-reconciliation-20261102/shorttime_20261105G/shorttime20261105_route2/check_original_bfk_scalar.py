"""Small check using actual v974 helper and constants; does not run v974.run.
The parameter l is original spatial circle circumference, not physical strip time.
"""
import sys,importlib.util,json
from pathlib import Path
import sympy as s
root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4');sys.path.insert(0,str(root/'verification'))
p=root/'verification/v974_alpha_faces_computed.py'
spec=importlib.util.spec_from_file_location('original_v974',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
checks={}
def need(name,x):
 checks[name]=bool(x)
 if not x:raise AssertionError(name)
l=s.symbols('original_spatial_circumference_l',positive=True);a=l/4
B=m._circ([2,-1,0,-1]);R=B/a
need('actual_v974_jump_pattern_and_zero_mode',B*s.ones(4,1)==s.zeros(4,1) and B.rank()==3)
detp=s.prod(e**mult for e,mult in R.eigenvals().items() if e!=0)
need('whole_nonzero_jump_spectrum_product',s.simplify(detp-4*l/a**4)==0)
full=s.Rational(1,16)*(2*a)**4*(l/4)*detp
need('whole_original_BFK_gluing_normalization',s.simplify(full-l*l)==0)
need('original_elliptic_spatial_anchor_jump',s.simplify(R.subs(l,2*s.pi)-16*m.C3S*B)==s.zeros(4))
need('original_elliptic_spatial_anchor_determinant',s.simplify(full.subs(l,2*s.pi)-(2*s.pi)**2)==0)
out={'scope':'Native original v974 helper/constant only, no full v974 replay. Scalar second-order spatial-circle determinant with one zero mode removed and four equal interval cuts. Not primitive quantum boson/fermion strip.',
 'parameter_type':'l=spatial circle circumference, a=l/4=spatial interval length. Original anchor l=2pi; not original source semigroup time t.',
 'checked':sum(checks.values()),'total':len(checks),'checks':checks,
 'formula':'2^-4 prod(2a) (l/4) detprime(R0)=l^2; R0=(1/a)circ(2,-1,0,-1)',
 'zero_mode_Gram_factor':'l/4','physical_short_time_jet_computed':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
