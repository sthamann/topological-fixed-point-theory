"""Exact CAR projection state -> oriented two-face amplitude.

The covariance is source data. Identification as a physical boundary CAR state
and the oriented state/map pairing are explicit premises. Wick is not supplied.
The result is an operator ray: unnormalized partition factors are not recovered.
"""
from __future__ import annotations
import sympy as s
from source_sewing import creator,source_lift


def state_from_projection(g):
    """Solve the positivity-forced annihilation conditions, without Wick input."""
    g=s.Matrix(g);n=g.rows
    if n!=g.cols or n>6 or g!=g.adjoint() or (g*g-g).applyfunc(s.simplify)!=s.zeros(n):
        raise ValueError('An exact Hermitian CAR projection of dimension <=6 is required')
    dim=1<<n
    creators=[creator(n,j) for j in range(n)]
    annihilators=[c.adjoint() for c in creators]
    constraints=[]
    for f in g.columnspace():
        constraints.append(sum((f[j]*creators[j] for j in range(n)),s.zeros(dim)))
    for f in (s.eye(n)-g).columnspace():
        constraints.append(sum((s.conjugate(f[j])*annihilators[j] for j in range(n)),s.zeros(dim)))
    null=s.Matrix.vstack(*constraints).nullspace()
    if len(null)!=1:
        raise ArithmeticError('Projection covariance did not force a unique state ray')
    omega=null[0]
    norm=(omega.adjoint()*omega)[0]
    # Recheck the actual CAR symbol convention C_ij=<c_j* c_i>.
    for i in range(n):
        for j in range(n):
            value=s.simplify((omega.adjoint()*creators[j]*annihilators[i]*omega)[0]/norm)
            if s.simplify(value-g[i,j])!=0:
                raise ArithmeticError('State covariance differs from the supplied projection')
    return omega,norm


def oriented_amplitude(g):
    """Unvectorize the boundary state using the fixed exterior input dual.

    Face order is first m input/dual coordinates, then m output coordinates:
    the original graph is (u,Ru). The complement sign moves output wedge factors
    to the end. Independent frame changes U on input and V on output give
    Q' = exterior(V) Q exterior(U*), not the reverse covariance.
    No antiunitary field map or physical quarter holonomy is inferred from it.
    """
    g=s.Matrix(g);n=g.rows
    if n%2:raise ValueError('Two equally sized labeled faces are required')
    m=n//2;dim=1<<m
    omega,norm=state_from_projection(g)
    base=omega[dim-1]
    if base==0:raise ValueError('Vacuum amplitude chart is singular; retain another operator-ray chart')
    q=s.zeros(dim)
    for row in range(dim):
        for col in range(dim):
            selected=[j for j in range(m) if col&(1<<j)]
            inversions=sum(1 for j in selected for k in range(j+1,m) if not col&(1<<k))
            index=((dim-1)^col)|(row<<m)
            q[row,col]=s.simplify((-1)**inversions*omega[index]/base)
    return q,{'boundary_modes':n,'field_modes':m,'state_ray_unique_from_CAR_projection':True,
              'oriented_state_map_pairing_assumed':True,'operator_scalar_recovered':False,
              'vacuum_coefficient_used_only_for_ray_chart':True}


def graph_projection(r):
    r=s.Matrix(r);v=s.eye(r.rows).col_join(r)
    return (v*(v.adjoint()*v).inv()*v.adjoint()).applyfunc(s.simplify)


def verify_graph_route(r):
    r=s.Matrix(r);g=graph_projection(r);q,receipt=oriented_amplitude(g)
    if q!=source_lift(r):raise ArithmeticError('Boundary-state route disagrees with full field covariance route')
    return q,receipt
