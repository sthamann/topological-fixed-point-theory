from pathlib import Path
from itertools import combinations, permutations
from collections import Counter
import json
import sympy as s

# Literal v809 incidence definition, with no replacement source dynamics.
duads=[frozenset(d) for d in combinations(range(6),2)]
synthemes=[frozenset(m) for m in combinations(duads,3) if len(m[0]|m[1]|m[2])==6]
synthemes.sort(key=lambda x:sorted(map(sorted,x)))
N=s.Matrix(15,15,lambda i,j:int(duads[i] in synthemes[j]))
assert all(x==3 for x in N*s.ones(15,1))
assert all(x==3 for x in N.T*s.ones(15,1))
A=s.Matrix(15,15,lambda i,j:int(not duads[i]&duads[j]))
assert N*N.T==3*s.eye(15)+A
M=N*N.T/9
assert M.eigenvals()=={s.Integer(1):1,s.Rational(4,9):9,s.Integer(0):5}
edges=[(i,15+j) for i in range(15) for j in range(15) if N[i,j]]
adj={i:set() for i in range(30)}
for i,j in edges:adj[i].add(j);adj[j].add(i)
visited={0};pending=[0]
while pending:
    i=pending.pop()
    for j in adj[i]-visited:visited.add(j);pending.append(j)
assert len(visited)==30 and len(edges)-30+1==16
cycles={}
for n in (4,6,8):
    found=set()
    def visit(path):
        if len(path)==n:
            if path[0] in adj[path[-1]]:
                rev=(path[0],*path[:0:-1])
                found.add(min(tuple(path),rev))
            return
        for j in adj[path[-1]]:
            if j>path[0] and j not in path:visit(path+[j])
    for i in range(30):visit([i])
    cycles[n]=found
assert [len(cycles[n]) for n in (4,6,8)]==[0,0,90]
ecount=Counter()
for c in cycles[8]:
    for a,b in zip(c,c[1:]+c[:1]):ecount[tuple(sorted((a,b)))]+=1
assert set(ecount.values())=={16}
base=[int(s.trace((N*N.T)**j)) for j in range(1,5)]
assert base==[45,225,1305,8865]
phase={}
i,j=next((i,j) for i in range(15) for j in range(15) if N[i,j])
for z in [-1,s.I]:
    Nz=N.copy();Nz[i,j]=z;B=Nz*Nz.conjugate().T
    phase[str(z)]=[int(s.trace(B**k)) for k in range(1,5)]
assert phase['-1']==[45,225,1305,8609]
assert phase['I']==[45,225,1305,8737]
H=s.eye(30)-s.BlockMatrix([[s.zeros(15),N],[N.T,s.zeros(15)]]).as_explicit()/3
assert H.eigenvals()=={s.Integer(0):1,s.Integer(2):1,s.Rational(1,3):9,s.Rational(5,3):9,s.Integer(1):10}
assert H[:15,15:]*H[15:,:15]==M
assert H[:15,:15]==s.eye(15) and H[15:,:15]==-N.T/3
# The original native quartic five ker(N^T) has energy one, no return coupling.
P5=(9*M-9*s.eye(15))*(9*M-4*s.eye(15))/36
assert P5**2==P5 and s.trace(P5)==5 and N.T*P5==s.zeros(15)
assert H*s.BlockMatrix([[P5],[s.zeros(15)]]).as_explicit()==s.BlockMatrix([[P5],[s.zeros(15)]]).as_explicit()
# Every complex weighted incidence has disjoint row supports on meeting duads.
q=5;inner=[d for d in duads if q not in d]
Ap=s.Matrix(10,10,lambda i,j:int(not inner[i]&inner[j]));V=Ap/3
R=s.Matrix(10,15,lambda i,j:int(inner[i]==duads[j]))
Pp=(Ap+2*s.eye(10))/3-s.ones(10)/6
Mt=M-R.T*Pp*R/3
blocked=[]
for d in inner:
    for e in inner:
        if d!=e and d&e:
            i,j=duads.index(d),duads.index(e)
            assert not any(N[i,k] and N[j,k] for k in range(15))
            assert Mt[i,j]==s.Rational(1,18)
            blocked.append((tuple(sorted(d)),tuple(sorted(e))))
assert len(blocked)==60
assert V**2==(2*s.eye(10)-Ap+s.ones(10))/9
assert all(x>0 for x in (V**2)**2)
# Character of H1 from exact cellular action. No edge signs: bipartition preserved.
def maps(g):
    dp=[duads.index(frozenset(g[i] for i in d)) for d in duads]
    sp=[synthemes.index(frozenset(frozenset(g[i] for i in d) for d in sy)) for sy in synthemes]
    return dp,sp
def chars(g):
    dp,sp=maps(g)
    ef=sum(dp[i]==i and sp[j]==j for i in range(15) for j in range(15) if N[i,j])
    return ef-sum(i==x for i,x in enumerate(dp))-sum(i==x for i,x in enumerate(sp))+1
hist=Counter(chars(g) for g in permutations(range(6)))
assert hist=={16:1,1:144,0:495,-2:80}
assert sum(n*c*c for c,n in hist.items())==720
def cyclelengths(g):
    unseen=set(range(len(g)));out=[]
    while unseen:
        a=min(unseen);b=a;count=0
        while b in unseen:unseen.remove(b);count+=1;b=g[b]
        out.append(count)
    return out
def spinchar(g):
    ls=cyclelengths(g)
    return 0 if any(n%2==0 for n in ls) else 2**(len(ls)-1)
hom=0;wedgepair=0;wedgenorm=0
for h in permutations(range(5)):
    g=(*h,5);c=chars(g);cs=spinchar(h)
    f=sum(i==h[i] for i in range(5))-1
    f2=sum(i==h[h[i]] for i in range(5))-1
    wedge=(f*f-f2)//2
    hom+=c*cs;wedgepair+=c*wedge;wedgenorm+=wedge*wedge
assert (hom,wedgepair,wedgenorm)==(120,120,120)
assert chars((1,2,0,3,4,5))==-2 and spinchar((1,2,0,3,4))==4
# Exact invariant spaces in the actual standard six-point realization.
def fixdim(gens):
    eq=s.Matrix.vstack(s.ones(1,6),*[s.eye(6)[list(g),:]-s.eye(6) for g in gens])
    return 6-eq.rank(),eq.nullspace()
id6=tuple(range(6))
def swap(i,j):
    g=list(id6);g[i],g[j]=g[j],g[i];return tuple(g)
tests={'trivial':[id6],'native_C3':[(0,2,3,1,4,5)],
       'mark_S5':[swap(i,i+1) for i in range(4)],
       'full_S6':[swap(i,i+1) for i in range(5)],
       'two_plus_four':[swap(0,1),swap(2,3),swap(3,4),swap(4,5)],
       'three_plus_three':[swap(0,1),swap(1,2),swap(3,4),swap(4,5)]}
dims={name:fixdim(gs)[0] for name,gs in tests.items()}
assert dims=={'trivial':5,'native_C3':3,'mark_S5':1,'full_S6':0,'two_plus_four':1,'three_plus_three':1}
res={'incidence_spectrum':{str(k):v for k,v in M.eigenvals().items()},
 'cycles':{k:len(v) for k,v in cycles.items()},'cycles_per_edge':16,
 'phase_moments':phase,'unphased_moments':base,'H_spectrum':{str(k):v for k,v in H.eigenvals().items()},
 'native_quartic_incidence_kernel_energy':1,'forbidden_ordered_inner_entries':len(blocked),
 'example_forbidden_transition':blocked[0], 'cycle_character_histogram':dict(hist),
 'S5_cycle_spinor_Hom_dimension':hom//120,'common_irrep_dimension':6,
 'C3_equivariant_rank_bound':12,'native_quartic_holonomy_fixed_dimensions':dims,
 'P1_source_identified':False}
Path(__file__).with_name('results.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
