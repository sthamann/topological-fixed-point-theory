#!/usr/bin/env python3
"""Exact boundary modular clock conditional on charge-changing intertwiners.

Physical sector-preserving algebra M2⊕M3: all braid-invariant faithful
states are central => modular time is identity. If a source permits
charge-changing boundary operators (full M5), the *same* density can
have a nontrivial relative modular phase. This extension is an input,
not a derivation of physical time.
"""
from sympy import Rational, sqrt, log, simplify, diag, Matrix, eye, zeros

phi = (1+sqrt(5))/2
D2 = 1+phi**2

def density(p2,p3):
    assert simplify(p2+p3-1)==0
    return diag(*([p2/2]*2+[p3/3]*3))

def gap(p2,p3):
    return simplify(log((p3/3)/(p2/2)))

def main():
    weights = {
        "five_state_trace":(Rational(2,5),Rational(3,5)),
        "equal_sectors":(Rational(1,2),Rational(1,2)),
        "quantum_dimensions_squared":(1/D2,phi**2/D2),
        "dimension_times_quantum_dim":(2/(2+3*phi),3*phi/(2+3*phi))
    }
    Y = diag(*([Rational(1,2)]*2 + [Rational(-1,3)]*3))
    X = zeros(5)
    X[0,2] = 1
    for label,(p2,p3) in weights.items():
        rho = density(p2,p3)
        # Show log rho = a I + beta Y exactly.
        l2,l3 = log(p2/2), log(p3/3)
        beta = (l2-l3)/Rational(5,6)
        a = l2-beta/2
        assert (diag(*([l2]*2+[l3]*3))-a*eye(5)-beta*Y).applyfunc(simplify)==zeros(5)
        # On M2⊕M3, all off-sector entries are zero, and modular flow is trivial.
        # On full M5, E02 changes charge, and its commutator is nonzero if p2/2!=p3/3.
        comm = (rho*X-X*rho).applyfunc(simplify)
        print(label,"omega=",gap(p2,p3),
              "offsector_modular_action?",comm!=zeros(5))
    assert simplify(gap(1/D2,phi**2/D2)-log(1+sqrt(5)/3))==0
    print("PASS exact: omega quantum-dimension-squared = log(1+sqrt(5)/3)")
    print("FAIL uniqueness: distinct invariant states give distinct omegas.")
    print("NOTE: modular charge rotation is NOT full geometric time.")

if __name__=="__main__":
    main()
