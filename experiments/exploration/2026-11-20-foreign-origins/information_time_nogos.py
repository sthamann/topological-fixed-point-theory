#!/usr/bin/env python3
"""Exact limited no-go probes, not general exclusions of informational or modular gravity."""
import json
from sympy import Matrix, eye, Rational, log, simplify


def test():
    p=Rational(1,2),Rational(1,3),Rational(1,6)
    # In local coordinates p1,p2 and p3=1-p1-p2, Fisher metric is
    # g_ij=sum_a (partial_i p_a)(partial_j p_a)/p_a.
    F=Matrix([[1/p[0]+1/p[2],1/p[2]], [1/p[2],1/p[1]+1/p[2]]])
    assert F==Matrix([[8,6],[6,9]]) and F.det()==36 and F[0,0]>0
    # Finite full matrix algebra M2 with U(2) invariant density rho=I/2.
    # Its modular generator log(rho) is scalar, so all commutators vanish.
    E12=Matrix([[0,1],[0,0]])
    Ktrace=log(Rational(1,2))*eye(2)
    assert Ktrace*E12-E12*Ktrace==Matrix.zeros(2)
    Knontrace=Matrix.diag(log(Rational(1,3)),log(Rational(2,3)))
    comm=Knontrace*E12-E12*Knontrace
    assert simplify(comm[0,1]-log(Rational(1,2)))==0
    return {'example_Fisher_matrix':str(F),'Fisher_first_principal_minor':8,
            'Fisher_determinant':36,'Fisher_metric_signature':'positive positive, not Lorentzian',
            'tracial_M2_modular_commutator_zero':True,
            'nontracial_M2_log_rho_commutator_offdiag':'log(1/2)',
            'limitations':'Only Fisher-as-the-physical-metric and finite fully unitarily invariant state; type-III modular flow and extra geometric structures are not excluded.'}


if __name__=='__main__':print(json.dumps(test(),indent=2,ensure_ascii=False))
