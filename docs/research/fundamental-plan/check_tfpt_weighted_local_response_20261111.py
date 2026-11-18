"""Original Round20/26 signed transfers and full scalar dipole Gaussian resolvent."""
from itertools import combinations, product
from pathlib import Path
import re

import numpy as np
from scipy.integrate import quad
import sympy as s

repo=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
base=repo/'experiments/theory-contracts'
text=(base/'local-charge-transport-round20/checker.py').read_text()
assignment=re.search(r'    eye = sp.eye\(8\)\n    B = .*?\n    G = B.T \* B\n',text,re.S)
assert assignment is not None
from textwrap import dedent
state={'sp':s}
exec(dedent(assignment.group()),state)
B,G=state['B'],state['G']
basis=tuple(B[:,i] for i in range(8))
assert B.det()==G.det()==1
assert [G[i,i]/2 for i in range(8)]==[2,1,1,1,1,1,1,1]
signed={tuple(sign*v) for sign,v in product((-1,1),basis)}

# Check all actual lattice root reflections on the frozen hopping basis.
roots=set()
eye=s.eye(8)
for i,j in combinations(range(8),2):
    for a,b in product((-1,1),repeat=2):
        roots.add(tuple(a*eye[:,i]+b*eye[:,j]))
for v in product((-1,1),repeat=8):
    if v.count(-1)%2==0:
        roots.add(tuple(s.Rational(x,2) for x in v))
assert len(roots)==240 and all(all(t.is_Integer for t in B.inv()*s.Matrix(r)) for r in roots)
def reflects_basis(root):
    r=s.Matrix(root)
    return {tuple(v-v.dot(r)*r) for v in basis}<=signed
preserving={r for r in roots if reflects_basis(r)}
expected={tuple(sign*(eye[:,i]-eye[:,j])) for i,j in combinations(range(1,7),2) for sign in (-1,1)}
assert preserving==expected and len(preserving)==30
assert tuple(eye[:,3]-eye[:,4]) in preserving and tuple(eye[:,0]-eye[:,1]) not in preserving
print('Frozen hopping support: 30/240 preserving roots, exactly A5; weak e3-e4 survives as a reflection, color e0-e1 does not in the standard D5 comparison.')

# Exact charge characters and sign-sensitive one-edge action, without
# truncating a lattice or dropping cocycle phases.
zero=(0,)*8
directions=[tuple(int(i==j) for i in range(8)) for j in range(8)]
def add(n,p):
    return tuple(x+y for x,y in zip(n,p))
def neg(p):
    return tuple(-x for x in p)
def eps(n,p):
    exp=sum(n[i]*p[i]*(int(G[i,i])//2) for i in range(8))
    exp+=sum(n[i]*p[j]*int(G[i,j]) for i in range(8) for j in range(i))
    return (-1)**(exp%2)
def hop(state,p,x,y):
    out=list(state)
    out[x]=add(state[x],neg(p))
    out[y]=add(state[y],p)
    return tuple(out),eps(p,out[x])*eps(p,state[y])

vac=(zero,zero,zero)
for i,p in enumerate(directions):
    changed,phase=hop(vac,p,0,1)
    back,back_phase=hop(changed,p,1,0)
    assert phase==(-1)**int(G[i,i]/2) and back==vac and back_phase*phase==1
    assert all(changed[0][j]==-p[j] and changed[1][j]==p[j] for j in range(8))
p,r=directions[1:3]
v,c=hop(vac,p,0,1)
w,d=hop(v,r,1,2)
v2,c2=hop(vac,r,1,2)
w2,d2=hop(v2,p,0,1)
assert w==w2 and c*d==-c2*d2
print('Original cocycle retained: dipole amplitudes have unit modulus and an actual adjacent odd-Gram pair anticommutes.')

# Reuse exactly the original periodic L=3 diagnostic, with the same a, m,
# u, nu and p=e8. Here compute its zero-temperature J² coefficient rather
# than the original beta=1 static oscillator normalizer.
L=3
a=m=1.0
u=1.0/3.0
nu=0.0
sites=list(product(range(L),repeat=3))
index={x:i for i,x in enumerate(sites)}
N=len(sites)
lap=np.zeros((N,N))
edges=set()
for i,x in enumerate(sites):
    for axis in range(3):
        y=list(x)
        y[axis]=(y[axis]+1)%L
        j=index[tuple(y)]
        edge=tuple(sorted((i,j)))
        assert edge not in edges
        edges.add(edge)
        lap[i,i]+=1
        lap[j,j]+=1
        lap[i,j]-=1
        lap[j,i]-=1
assert len(edges)==3*N and np.all(np.diag(lap)==6)

# Reconstruct site stars only from the nontrivial bond commutator relation.
# The relation is available from the actual signed shift test and G mod 2.
ordered_edges=sorted(edges)
neighbors=[{j for j,f in enumerate(ordered_edges) if len(set(e)&set(f))==1}
           for e in ordered_edges]
stars=set()
for i, near in enumerate(neighbors):
    for j,k,l in combinations(sorted(near),3):
        clique={i,j,k,l}
        if all(b in neighbors[a] for a,b in combinations(clique,2)):
            stars.add(frozenset(clique | set.intersection(*(neighbors[a] for a in clique))))
assert len(stars)==N and all(len(star)==6 for star in stars)
assert all(sum(j in star for star in stars)==2 for j in range(len(ordered_edges)))
print(f'Bond-phase line graph: recovered {len(stars)} unlabeled six-bond site stars from {len(ordered_edges)} bonds.')

x,y=index[(0,0,0)],index[(1,0,0)]
dipole=np.zeros(N,dtype=int)
dipole[x],dipole[y]=-1,1
assert sum((dipole[i]-dipole[j])**2 for i,j in edges)==14
charge_energy=float(G[7,7]/2)
delta=(2/N+14*nu)*charge_energy
assert delta>0
positions=np.zeros(N)
positions[[x,y]]=1

def matrix_sqrt(mat):
    ev,U=np.linalg.eigh(mat)
    vals=np.sqrt(ev)
    return vals,U,(U*vals)@U.T

def resolvent(a,u):
    K=lap/a**2+m*m*np.eye(N)
    w0,_,W0=matrix_sqrt(K)
    w,V,_=matrix_sqrt(K+2*u*charge_energy*np.diag(positions))
    logdet0=np.log(w0).sum()
    e_diff=(w.sum()-w0.sum())/2
    def integrand(t):
        if t==0:
            return 1.0
        h=np.tanh(t*w/2)
        a_plus=W0+(V*(w*h))@V.T
        a_minus=W0+(V*(w/h))@V.T
        signp,logp=np.linalg.slogdet(a_plus)
        signm,logm=np.linalg.slogdet(a_minus)
        assert signp==signm==1
        log_prefactor=0.5*(logdet0+np.log(4*w).sum()
                            -np.log(-np.expm1(-2*t*w)).sum()-logp-logm)
        return np.exp(log_prefactor-t*(delta+e_diff))
    T=48/delta
    value,error=quad(integrand,0,T,epsabs=2e-11,epsrel=2e-11)
    tail=np.exp(-delta*T)/delta
    return value,error,tail
bare,_,_=resolvent(a,0)
assert abs(bare-1/delta)<3e-9
dressed,err,tail=resolvent(a,u)
assert 0<dressed<bare
ell=np.linalg.eigvalsh(lap)
v_a_prime=np.sum(ell/(m*m+ell/a**2)**1.5)/(2*N*a**3)
assert v_a_prime>0
mixed=-8*charge_energy*v_a_prime/delta**2
eps_u=1e-5
def metric_slope(u):
    eps_a=5e-4
    return 4*(resolvent(a+eps_a,u)[0]-resolvent(a-eps_a,u)[0])/(2*eps_a)
numeric_mixed=(metric_slope(eps_u)-metric_slope(0))/eps_u
assert abs(numeric_mixed/mixed-1)<4e-4
print(f'Original L=3 diagnostic, p=e8: gap cost={delta:.12g}; bare and dressed J² coefficients={4*bare:.12g}, {4*dressed:.12g}.')
print(f'Gaussian quadrature estimator={err:.2g}; analytic integration tail bound={tail:.2g}.')
print(f'Mixed spacing/scalar derivative of J² coefficient: analytic {mixed:.9g}, independent difference {numeric_mixed:.9g}.')
