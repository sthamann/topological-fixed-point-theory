"""Exact native neutral s-orbit; arbitrarily large carries, no charge cutoff."""
import ast,json,hashlib
from pathlib import Path
import sympy as sp
ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4');source=ROOT/'experiments/theory-contracts/local-charge-transport-round20/checker.py'
t=ast.parse(source.read_text());body=next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='main').body
assign=[n for n in body if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ('eye','B','G')]
ns={'sp':sp};exec(compile(ast.Module(body=assign,type_ignores=[]),str(source),'exec'),ns);G=ns['G'];s=sp.eye(8)[:,7]
def c(a,b):return int(sum((G[i,i]//2)*a[i]*b[i] for i in range(8))+sum(G[i,j]*a[i]*b[j] for i in range(8) for j in range(i)))%2
checks=0
def ck(ok,name):
 global checks
 if not ok:raise AssertionError(name)
 checks+=1
samples=[-10**12,-19,-5,-4,-3,-2,-1,0,1,2,3,4,7,19,10**12]
def R(j,power=1):
 m,r=divmod(j,4);return 4*m+(r+power)%4
# After the native along-orbit gauge d_j=(-1)^j, W raises j and Z=W^4 raises j by4.
def E(j,a,b):
 if j%4!=b:return None
 return 4*(j//4)+a
for j in samples:
 # Exact source W phase on departure -js, arrival+js; its adjoint phase is retained.
 phase=(-1)**(c(s,-(j+1)*s)+c(s,j*s));ck(phase==-1,'original native W sign')
 ck(((-1)**(j+1))*phase*((-1)**j)==1,'coherent source-orbit gauge')
 out=j
 for _ in range(4):out=R(out)
 ck(out==j,'R4=I with original winding retained')
 ck(R(j+4)==R(j)+4,'R commutes with original carryZ')
 wrapped=R(j+(4 if j%4==3 else 0));ck(wrapped==j+1,'W=R(diag1,1,1,Z)')
 for a in range(4):
  for b in range(4):
   # R^a P0 R^-b action is exactly E_ab.
   z=R(j,-b);z=R(z,a) if z%4==0 else None
   ck(z==E(j,a,b),'native matrixunit dictionary')
   for d in range(4):
    for e in range(4):
     z=E(j,d,e);z=E(z,a,b) if z is not None else None
     expected=E(j,a,e) if b==d else None
     ck(z==expected,'matrixunit multiplication')
# Native W^4 commutes every neutral basis hopping; cocycle pairings are multiplesof4.
for p in [sp.eye(8)[:,i] for i in range(8)]:
 ck(int((p.T*G*(4*s))[0])%2==0,'all hopping words commute with original4s carry')
j=sp.Symbol('j',integer=True);D_r,u,phi2=sp.symbols('D_r u phi2',positive=True)
ck(sp.expand((((j+4)**2-j**2)-(8*j+16))*(D_r+u*phi2))==0,'same-source carry energy difference')
result={'status':'PASS','checks':checks,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'integer_samples':samples,'native_W_phase':-1,'derived_dictionary':{'grade':'r=j mod4','carry':'m=floor(j/4)','Z':'W^4; m->m+1','R':'W(I-P3)+W^-3P3; R^4=I','matrix_units':'E_ab=R^a P0 R^-b','algebra':'M4(C(S1)) on each native neutral s-orbit','energy_change_W4':'(8j+16)[D_r+u(phi_x^2+phi_y^2)]'},'scope':'Exact source-orbit operator action and algebra; no finite truncation, no identification with the earlier marked CPTP channel, no canonical origin of the already supplied grade character.'}
Path(__file__).with_name('carry-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
