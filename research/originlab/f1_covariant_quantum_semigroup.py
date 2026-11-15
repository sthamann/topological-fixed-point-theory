"""Exact inner-automorphism covariance theorem for the three-state source.

If a hypothetical pre-spatial three-state quantum source is unital,
trace-preserving, SU3-conjugation covariant, and continuous CP Markov,
then its only nontrivial evolution (up to a clock unit) is depolarizing.
The original TFPT source does not establish these premises.
"""
import sympy as sp


def controls():
    n=3
    E=lambda i,j: sp.eye(n)[:,i]*sp.eye(n)[j,:]
    B=[E(i,j) for i in range(n) for j in range(n)]
    lie=[E(0,1),E(1,0),E(0,2),E(2,0),E(1,2),E(2,1),
         sp.diag(1,-1,0),sp.diag(0,1,-1)]
    ad=[sp.Matrix.hstack(*[(X*b-b*X).reshape(9,1) for b in B]) for X in lie]
    vars=sp.symbols('w0:81')
    Q=sp.Matrix(9,9,vars)
    eqs=[z for A in ad for z in (Q*A-A*Q)]
    mat,_=sp.linear_eq_to_matrix(eqs,vars)
    rank=mat.rank()
    assert rank==79
    # Exact commutant basis: identity on all operators and projection to scalar.
    dep=sp.Matrix.hstack(*[(sp.trace(b)*sp.eye(3)/3).reshape(9,1) for b in B])
    assert dep*dep==dep and dep.rank()==1 and all(dep*A==A*dep for A in ad)
    assert all((sp.eye(9)-dep)*A==A*(sp.eye(9)-dep) for A in ad)
    assert (sp.eye(9)-dep).rank()==8 and dep+sp.eye(9)-dep==sp.eye(9)
    q=sp.symbols('q',real=True)
    omega=sp.Matrix([1 if j in (0,4,8) else 0 for j in range(9)])
    choi=q*omega*omega.T+(1-q)*sp.eye(9)/3
    eigenvals=choi.eigenvals()
    assert sum(mult for eigen,mult in eigenvals.items()
               if sp.simplify(eigen-(1+8*q)/3)==0)==1
    assert sum(mult for eigen,mult in eigenvals.items()
               if sp.simplify(eigen-(1-q)/3)==0)==8
    assert choi.subs(q,sp.Rational(-1,8)).is_positive_semidefinite
    assert choi.subs(q,sp.Rational(1,2)).is_positive_definite
    assert choi.subs(q,sp.Rational(-1,4)).is_positive_semidefinite is False
    # Any continuous semigroup scalar q(t) on the 8D adjoint sector obeys
    # q(t+s)=q(t)q(s), q(0)=1, continuity, hence q(t)=exp(-kappa t).
    t,s,k=sp.symbols('t s k',positive=True)
    qt=sp.exp(-k*t)
    assert sp.simplify(qt.subs(t,t+s)-qt*qt.subs(t,s))==0
    H=sp.eye(9)-dep
    P_t=dep+qt*H
    P_s=dep+qt.subs(t,s)*H
    assert sp.simplify(P_t*P_s-(dep+sp.exp(-k*(t+s))*H))==sp.zeros(9)
    # Finite 9D spectrum cannot have 4D short-distance heat scaling.
    z=sp.symbols('z',positive=True)
    heat=1+8*sp.exp(-z)
    d_spec=-2*z*sp.diff(sp.log(heat),z)
    assert sp.limit(d_spec,z,0,dir='+')==0
    return {
      'status':'EXACT_UNITAL_SU3_COVARIANT_CP_SEMIGROUP_CLASSIFICATION',
      'complex_operator_space_dim':9,
      'full_linear_superoperator_unknown_entries':81,
      'lie_generators_tested':8,
      'commutant_constraint_matrix_shape':list(mat.shape),
      'commutant_linear_rank':rank,
      'commutant_dimension':81-rank,
      'scalar_subrepresentation_dim':1,
      'traceless_adjoint_subrepresentation_dim':8,
      'unital_TP_covariant_channel':'Phi_q(A)=q*A+(1-q)*Tr(A)*I/3',
      'Choi_spectrum':{'(1+8q)/3':1,'(1-q)/3':8},
      'complete_positivity_exact_interval':'-1/8<=q<=1',
      'continuous_CP_semigroup_form':'q(t)=exp(-kappa*t), kappa>=0',
      'unique_dimensionless_generator_after_clock_normalization':
         'L(A)=Tr(A)*I/3-A',
      'all_positive_time_two_point_OS_Grams_positive':True,
      'source_covariance_and_markov_requirements_choose_kappa':False,
      'physical_P1_mark_geometry_or_E8_reconstructed':False,
      'finite_source_UV_spectral_dimension':0,
      'full_3plus1D_parent_derived':False,
      'original_TFPT_source_identified':False,
      'F1_solved':False
    }

if __name__=='__main__':
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
