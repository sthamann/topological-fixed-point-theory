"""Exact inverse using existing TFPT Grassmann helper; no model/weight choice."""
from pathlib import Path
import importlib.util,json,itertools
import sympy as s
root=Path(__file__).resolve().parents[3]
if not (root/'verification').is_dir():root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
p=root/'experiments/theory-contracts/clock-interaction-provenance/checker.py'
sp=importlib.util.spec_from_file_location('original_grassmann_helper',p);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
N=4;top=(1<<N)-1;even=[i for i in range(1<<N) if i.bit_count()%2==0]
cs=s.symbols('w0:16');w={i:cs[i] for i in even}
recovered={}
for i in even:
 j=top^i;v=m.grassmann_mul({j:s.Integer(1)},w).get(top,0)*m.sign(j,i)
 if s.expand(v-cs[i])!=0:raise RuntimeError(('inverse',i))
 recovered[i]=v
# finite log and inverse exponential, including a genuine arbitrary quartic
nil={i:cs[i]/cs[0] for i in even if i}
log=m.grassmann_add((1,nil),(-s.Rational(1,2),m.grassmann_mul(nil,nil)))
exp=m.grassmann_add((1,{0:1}),(1,log),(s.Rational(1,2),m.grassmann_mul(log,log)))
if any(s.simplify(exp.get(i,0)*cs[0]-w.get(i,0))!=0 for i in range(16)):raise RuntimeError('roundtrip')
# paired U(1) original-type phase representation. Invariant probes suffice
# to separate all invariant densities, no charged elementary access required.
charges=[1,1,-1,-1]
inv=[i for i in even if sum(charges[k] for k in range(N) if i>>k&1)==0]
B=s.Matrix([[m.grassmann_mul({a:1},{b:1}).get(top,0) for b in inv] for a in inv])
if B.rank()!=len(inv):raise RuntimeError('invariant pairing')
# Singular fiber diagnostic. w=1+a theta12+b theta34+(a*b-c)tau.
# No Pfaffian inverse enters top-insertion body recovery.
a,b,c=s.symbols('a b c');sing={0:s.Integer(1),3:a,12:b,15:a*b-c}
saturation=m.grassmann_mul({top:1},sing).get(top,0)
partition=sing[top].subs(c,a*b)
if saturation!=1 or partition!=0:raise RuntimeError('singular body')
out={'scope':'Exact generic finite exterior algebra inverse, paired invariant tests and singular-fiber control. No primitive physical top-source admissibility or numerical boson law asserted.','even_coefficients_recovered':len(recovered),'invariant_density_dimension':len(inv),'invariant_pairing_rank':B.rank(),'invariant_pairing_determinant':str(B.det()),'log_roundtrip':True,'singular_partition':str(partition),'singular_body_recovered':str(saturation),'source_density_coefficients':{str(i):str(v) for i,v in recovered.items()}}
Path(__file__).with_name('berezin-inverse-results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
