"""Bounded implication checks, no proposed TFPT dynamics or fitted constants."""
import json
from pathlib import Path
import sympy as s
r,a,b,z,lam,t,eps,k=s.symbols('r a b z lambda t epsilon k', real=True)
I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Z=s.diag(1,-1)
results={}
def record(name,expression,target=0):
    delta=expression-target
    if isinstance(delta,s.MatrixBase):
        ok=all(s.simplify(v)==0 for v in delta)
    else: ok=s.simplify(delta)==0
    if not ok: raise RuntimeError(name+': '+str(delta))
    results[name]={'passed':True,'identity':str(expression),'target':str(target)}
# Actual original two-face graph and oriented state-amplitude pairing.
v=s.Matrix([1,r]); G=v*v.T/(1+r*r)
record('original_two_face_projection',G*G,G)
record('boundary_number_connected',0-G[0,0]*G[1,1],-r*r/(1+r*r)**2)
record('identity_transfer_entangled_boundary',G.subs(r,1),s.ones(2)/2)
# Independent Hilbert-space operator vectorization, same identity transfer.
vecI=s.Matrix([1,0,0,1])/s.sqrt(2)
rho=vecI*vecI.T
red=s.Matrix([[sum(rho[2*i+j,2*l+j] for j in range(2)) for l in range(2)] for i in range(2)])
record('identity_operator_reduced_state',red,I/2)
# Exact block decomposition of a genuinely split, deck-equivariant generator.
H=s.Matrix([[a,b],[b,a]])
U=s.Matrix([[1,1],[1,-1]])/s.sqrt(2)
record('deck_equivariance_does_not_fix_b',H*X-X*H,s.zeros(2))
record('even_odd_generator',U.T*H*U,s.diag(a+b,a-b))
record('visible_resolvent_schur',(z*I-H).inv()[0,0],1/(z-a-b*b/(z-a)))
# Positivity does not pick a coupling: e.g. a>|b| is a whole interval.
record('visible_short_time_quadratic_coefficient',(H*H)[0,0],a*a+b*b)
# Shared-source linear exponential insertions. Derivatives of log Z are
# connected algebraic cumulants. The Gaussian specialization is diagnostic only.
x,y,C=s.symbols('x y C',commutative=True)
W=C*(x+y)**2/2
record('shared_source_mixed_cumulant',s.diff(W,x,y),C)
S_eff=-W
record('shared_source_cross_bilinear_coefficient',s.diff(S_eff,x,y),-C)
# A correlated state with an exactly uncoupled bipartite dynamics.
HA=s.kronecker_product(eps*Z,I)+s.kronecker_product(I,eps*Z)
A=s.kronecker_product(X,I); B=s.kronecker_product(I,X)
record('initial_correlations_no_instant_response',A*B-B*A,s.zeros(4))
record('uncoupled_local_generator_remains_other_commutant',(HA*A-A*HA)*B-B*(HA*A-A*HA),s.zeros(4))
record('bell_covariance_without_interaction',(vecI.T*s.kronecker_product(Z,Z)*vecI)[0],1)
# An explicit coupling demonstrates the diagnostic; NOT a TFPT-selected value.
Hcoupled=HA+k*s.kronecker_product(X,X)
AZ=s.kronecker_product(Z,I); BZ=s.kronecker_product(I,Z)
record('causal_cross_commutator_derivative',s.I*((Hcoupled*AZ-AZ*Hcoupled)*BZ-BZ*(Hcoupled*AZ-AZ*Hcoupled)), -4*s.I*k*s.kronecker_product(s.Matrix([[0,-s.I],[s.I,0]]),s.Matrix([[0,-s.I],[s.I,0]])))
# Common scalar factors preserve normalized responses, carry an unfixed additive
# generator; original absolute source would distinguish them.
q0,qJ,c,t1,t2=s.symbols('q0 qJ c t1 t2',nonzero=True,real=True)
record('normalized_response_forgets_scalar',s.exp(-c*t)*qJ/(s.exp(-c*t)*q0),qJ/q0)
record('scalar_sewing_character',s.exp(-c*t1)*s.exp(-c*t2),s.exp(-c*(t1+t2)))
Path(__file__).resolve().with_name('results.json').write_text(json.dumps({'scope':'Finite exact implication checks, conditional original boundary chart; not a physical TFPT candidate, not a new fundamental solution. Shared-source cumulant hierarchy reuses R112.','checks':results},indent=2))
print(f'{len(results)} exact symbolic identities passed; no selected physical coupling, primitive source or spectrum claimed.')
print('Boundary connected occupation at r=1: -1/4; vectorized identity reduced density: I/2; split-generator eigenvalues a+b,a-b; mixed induced coefficient -C.')
