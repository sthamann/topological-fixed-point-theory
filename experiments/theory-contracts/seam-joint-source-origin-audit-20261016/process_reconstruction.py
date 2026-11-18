"""Exact diagnostic for process-first reconstruction; no physical source fit.

Checks finite counterexamples and exhausts quartic stationary spectra.
The analytic argument and limits are in PROCESS_RECONSTRUCTION.txt.
"""
import json
import sys
from pathlib import Path
import sympy as s


def operational_checks():
    z, o, i = s.zeros(2), s.eye(2), s.I
    pauli = [o, s.Matrix([[0, 1], [1, 0]]),
             s.Matrix([[0, -i], [i, 0]]), s.diag(1, -1)]
    vec = lambda x: s.Matrix(list(x))
    f = s.Matrix.hstack(*(vec(x)/s.sqrt(2) for x in pauli))
    h = s.Matrix([[1, 1, 1, 1], [1, -1, 1, -1],
                  [1, 1, -1, -1], [1, -1, -1, 1]])/2
    e = f*h
    omega = vec(o)/s.sqrt(2)
    classical = [s.simplify(e[:, k]*e[:, k].conjugate().T) for k in range(4)]
    units = [s.Matrix(2, 2, lambda a,b: int(a == j and b == k))
             for j in range(2) for k in range(2)]
    left = [s.kronecker_product(x, o) for x in units]
    assert s.simplify(e.conjugate().T*e) == s.eye(4)
    assert all(s.simplify(x.conjugate().reshape(2,2).T-x.reshape(2,2)) == z
               for x in [e[:,k] for k in range(4)])
    assert s.Matrix.hstack(*(a*omega for a in classical)).rank() == 4
    assert s.Matrix.hstack(*(a*omega for a in left)).rank() == 4
    assert all(s.simplify(a*b-b*a) == s.zeros(4) for a in classical for b in classical)
    assert any(a*b != b*a for a in left for b in left)
    # A1 opposite action is matrix right multiplication; C4 opposite is itself.
    right = [s.kronecker_product(o, x.T) for x in units]
    assert all(a*b == b*a for a in left for b in right)
    hh = s.Matrix([[1,1],[1,-1]])/s.sqrt(2)
    rho = s.diag(1,0)
    coherent = hh*hh*rho*hh.conjugate().T*hh.conjugate().T
    transition = hh.applyfunc(lambda v: s.simplify(s.conjugate(v)*v))
    dephased = transition*transition*s.Matrix([1,0])
    assert coherent[0,0] == 1 and dephased[0] == s.Rational(1,2)
    return {"same_J_same_omega_cyclic_algebras": ["C^4", "M_2(C)_left"],
            "cyclic_dimensions": [4,4], "centers_differ": True,
            "both_order_zero": True,
            "two_coherent_pulses_vs_dephased_probability": ["1", "1/2"],
            "pulse_derived_from_P1": False}


def quartic_stationary():
    records=[]
    for k in range(1,5):
        q=s.Rational((5-k)**3+k**3,25*k*(5-k))
        records.append({"multiplicities":[k,5-k],"quartic":q})
    for m in range(1,4):
        for n in range(1,5-m):
            p=5-m-n
            if p < 1: continue
            mult=(m,n,p)
            values=(n-p,p-m,m-n)
            assert sum(a*b for a,b in zip(mult,values)) == 0
            assert sum(values) == 0
            if len(set(values)) != 3: continue
            norm=sum(a*b*b for a,b in zip(mult,values))
            q=s.Rational(sum(a*b**4 for a,b in zip(mult,values)),norm**2)
            records.append({"multiplicities":list(mult),"quartic":q})
    qmin=min(r['quartic'] for r in records)
    qmax=max(r['quartic'] for r in records)
    assert (qmin,qmax)==(s.Rational(7,30),s.Rational(13,20))
    triples=sorted(set(r['quartic'] for r in records if len(r['multiplicities'])==3))
    assert triples==[s.Rational(1,4),s.Rational(1,2)]
    # Competing even invariant polynomial controls, BOTH bounded below.
    # V_sigma=-r^2 + r^4(1+sigma*q), sigma=+/-1.
    assert 1-qmax > 0
    return {"stationary_spectra": [{**r,"quartic":str(r['quartic'])} for r in records],
            "minimum":str(qmin), "maximum":str(qmax),
            "positive_quartic_selected_split":[2,3],
            "negative_quartic_selected_split":[1,4],
            "both_polynomials_bounded_below":True,
            "field_trace_constraint_sign_derived_from_P1":False}


if __name__ == '__main__':
    result={"verdict":"SCOPED_RECONSTRUCTION_TESTS_ONLY",
            "operational":operational_checks(),"quartic":quartic_stationary(),
            "physical_algebra_selected":False,"raw_P1_32_derivation":False}
    target=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('process_reconstruction_results.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
