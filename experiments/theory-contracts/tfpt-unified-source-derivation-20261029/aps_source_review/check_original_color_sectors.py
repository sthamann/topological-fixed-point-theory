"""R069: exact original marked color lattice/field-sector tests, not physical origin."""
from pathlib import Path
import sys,json,math
from collections import Counter,defaultdict
from fractions import Fraction as F
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
r=Path(__file__).resolve().parents[4];sys.path.insert(0,str(r/'verification'))
import v1_e8_glue as v
checks=[]
def eq(name,actual,expected):
 if actual!=expected:raise RuntimeError((name,actual,expected))
 checks.append(dict(name=name,value=str(actual)))
roots=sorted(tuple(int(2*x) for x in q) for q in v.e8_roots());R=set(roots)
eq('unchanged original root constructor',len(R),240)
def dot(a,b):return F(sum(x*y for x,y in zip(a,b)),4)
def trial(q):
 n=q[0]+q[1]-2*q[2];assert n%2==0
 return (n//2)%3
def root(i,j,si=1,sj=1):
 q=[0]*8;q[i]=2*si;q[j]=2*sj;return tuple(q)
a=root(0,1,1,-1);bb=root(1,2,1,-1)
color={q for q in R if all(x==0 for x in q[3:]) and sum(q[:3])==0}
perp={q for q in R if dot(q,a)==dot(q,bb)==0}
sector={k:{q for q in R if trial(q)==k} for k in range(3)}
eq('ACTUAL first-three color roots',len(color),6)
eq('ACTUAL color-orthogonal roots',len(perp),72)
eq('neutral root sector equals union E6 and COLOR A2',sector[0],perp|color)
eq('three root-sector dimensions',[len(sector[k]) for k in range(3)],[78,81,81])
eq('weight-one neutral incl eight original Cartans',len(sector[0])+8,86)
def basis(rs):return hermite_normal_form(sp.Matrix(sorted(rs)).T)
lb=basis(R);eb=basis(perp);ab=basis(color);nb=basis(sector[0]);prod=eb.row_join(ab)
eq('E8 original root lattice rank',lb.rank(),8)
eq('color perpendicular rank',eb.rank(),6)
eq('color root rank',ab.rank(),2)
eq('E8 determinant',sp.det(lb.T*lb/4),1)
eq('E6 determinant',sp.det(eb.T*eb/4),3)
eq('color A2 determinant',sp.det(ab.T*ab/4),3)
eq('neutral root lattice determinant',sp.det(nb.T*nb/4),9)
eq('product lattice has original neutral basis',hermite_normal_form(prod),nb)
eq('index of neutral lattice in E8',abs(sp.det(nb)/sp.det(lb)),3)
# The kernel has index3 by explicit surjective center character; equal index
# proves equality for ALL lattice vectors, not just the norm2 root census.
ward=root(2,5);higgs=root(3,5)
eq('original Ward root exists',ward in R,True)
eq('original weak mediator root exists',higgs in R,True)
eq('original Ward root triality',trial(ward),1)
eq('original weak mediator triality',trial(higgs),0)
eq('weak mediator belongs actual color centralizer',higgs in perp,True)
d5={q for q in R if all(x==0 for x in q[5:])}
eq('original D5 roots',len(d5),40)
eq('color-centralizer E6 is NOT old D5-containing E6',len(d5&perp),4)
weak={q for q in R if sum(x!=0 for x in q[:5])==1 and any(q[i] for i in [3,4]) and sum(x!=0 for x in q[5:])==1}
eq('whole original weak mediator current block has24 roots',len(weak),24)
eq('whole weak mediator block is color centralizer',weak<=perp,True)
# Project the charged roots onto the ACTUAL color plane; 3 weights x27.
weights={}
for k in [1,2]:
 counts=Counter();norms=Counter();enorms=Counter();groups=defaultdict(set)
 for q in sector[k]:
  x=[F(t,2) for t in q];avg=sum(x[:3])/3
  qc=tuple(x[i]-avg if i<3 else F(0) for i in range(8));qe=tuple(x[i]-qc[i] for i in range(8))
  counts[qc]+=1;norms[sum(t*t for t in qc)]+=1;enorms[sum(t*t for t in qe)]+=1;groups[qc].add(qe)
 eq(f'charged sector{k} actual color weight multiplicities',sorted(counts.values()),[27,27,27])
 eq(f'charged sector{k} color projection norm',dict(norms),{F(2,3):81})
 eq(f'charged sector{k} perpendicular norm',dict(enorms),{F(4,3):81})
 eq(f'charged sector{k} common27 perpendicular weights',len({frozenset(s) for s in groups.values()}),1)
 weights[k]=[dict(color=[str(t) for t in q],multiplicity=n) for q,n in sorted(counts.items())]
# Actual group action test: generate weight orbits using E8 root reflections.
def orbit(seed,subroots):
 seen={seed};todo=[seed]
 while todo:
  q=todo.pop()
  for rr in subroots:
   ip=dot(q,rr);assert ip.denominator==1
   qq=tuple(x-int(ip)*y for x,y in zip(q,rr));assert qq in R
   if qq not in seen:seen.add(qq);todo.append(qq)
 return seen
for k,seed in [(1,ward),(2,tuple(-x for x in ward))]:
 eq(f'actual E6 root action generates charged27 sector{k}',len(orbit(seed,perp)),27)
 eq(f'actual COLOR root action generates triplet sector{k}',len(orbit(seed,color)),3)
# Check all actual nonzero root brackets respect original triality. Global
# all-mode theorem follows from the additive lattice character, not a cutoff.
trans=Counter()
for q in R:
 for s in R:
  total=tuple(x+y for x,y in zip(q,s))
  if total in R:
   assert trial(total)==(trial(q)+trial(s))%3
   trans[(trial(q),trial(s),trial(total))]+=1
checks.append(dict(name='all original nonzero root bracket sector transitions',value={str(k):n for k,n in sorted(trans.items())}))
# Original charged insertion, adjoint and a three-charged balanced root word.
neg=tuple(-x for x in ward)
eq('adjoint Ward grade',trial(neg),2)
eq('Ward root normalization for original level1 vacuum coefficient',dot(ward,ward),2)
beta=root(1,5,1,-1);gamma=tuple(-x-y for x,y in zip(ward,beta))
eq('three-field balanced root word exists',all(q in R for q in [ward,beta,gamma]),True)
eq('three-field balanced root grades',[trial(q) for q in [ward,beta,gamma]],[1,1,1])
eq('three-field balanced root total',tuple(sum(t) for t in zip(ward,beta,gamma)),(0,)*8)
eq('actual root-string length for charged OPE',dot(ward,beta),-1)
eq('neutral color root is NOT a color singlet',dot(a,a),2)
# Complete lattice shells through conformal grade3, then ALL8 oscillator colors.
N=3;shell=[[0]*(N+1) for _ in range(3)]
def walk(q,left,par):
 if len(q)==8:
  if sum(q)%4==0:
   n=sum(t*t for t in q)//8
   assert sum(t*t for t in q)%8==0
   shell[trial(q)][n]+=1
  return
 bound=math.isqrt(left)
 for x in range(-bound,bound+1):
  if x%2==par:walk(q+[x],left-x*x,par)
for par in [0,1]:walk([],8*N,par)
eq('complete original E8 lattice shell totals',[sum(shell[k][n] for k in range(3)) for n in range(N+1)],[1,240,2160,6720])
osc=[1]+[0]*N
for m in range(1,N+1):
 for _ in range(8):
  for n in range(m,N+1):osc[n]+=osc[n-m]
chars=[[sum(shell[k][j]*osc[n-j] for j in range(n+1)) for n in range(N+1)] for k in range(3)]
eq('complete original E8 character grades0to3',[sum(chars[k][n] for k in range(3)) for n in range(N+1)],[1,248,4124,34752])
eq('neutral full-source character grades0to3',chars[0],[1,86,1370,11586])
eq('charged full-source characters coincide',chars[1],chars[2])
eq('charged full-source character grades0to3',chars[1],[0,81,1377,11583])
# Entire discriminant class, not just a candidate search.
def frac(x):return x-x.numerator//x.denominator
qforms=[(k,l,frac(F(2*k*k+l*l,3))) for k in range(3) for l in range(3)]
eq('all9 inherited lattice spin sectors',dict(Counter(q for k,l,q in qforms)),{F(0):5,F(1,3):2,F(2,3):2})
eq('no inherited fermionic halfinteger sector',any(q==F(1,2) for k,l,q in qforms),False)
eq('diagonal index3 original glue is isotropic',[frac(F(2*k*k+k*k,3)) for k in range(3)],[F(0)]*3)
eq('allowed integral overlattice indices from determinant9',[d for d in range(1,4) if 9%(d*d)==0],[1,3])
out=dict(verdict='SCOPED_ORIGINAL_COLOR_CENTER_FIXED_SOURCE_AND_CHARGED_FIELD_WORDS',controls=len(checks),checks=checks,lattice_shells=shell,oscillator_coefficients=osc,full_characters=chars,charged_color_weights=weights,
 scope='Existing marked E8 source and its original color center; no primitive P1/Padm/96CAR, measure, nativeJ closure or physical3+1/spin2 derivation; generic E6+A2 glue and old J obstruction reused.')
Path(__file__).with_name('ORIGINAL-COLOR-SECTORS.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(verdict=out['verdict'],controls=len(checks),full_characters=chars)))
