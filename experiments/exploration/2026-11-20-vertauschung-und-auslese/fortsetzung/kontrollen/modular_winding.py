"""KMS winding on the positive braid monoid and the endpoint meter."""
import sympy as sp
from path_os import E,A,B,inv,mul


def mod2(g):
    return tuple(e%2 for e in g[1:])


def levels(steps):
    rows=[{E}]
    for k in range(steps):
        rows.append({mul(s,g) for s in (A,B) for g in rows[-1]})
    return rows


def main():
    rows=levels(14)
    delta=mul(mul(A,B),A)
    assert delta==mul(mul(B,A),B)
    # Verify both inclusion-exclusion and the endpoint grading independently.
    for n in range(3,len(rows)):
        from_a={mul(A,g) for g in rows[n-1]}
        from_b={mul(B,g) for g in rows[n-1]}
        assert from_a & from_b == {mul(delta,g) for g in rows[n-3]}
    print("Distinct positive braids at degrees 0..14:",list(map(len,rows)))
    q=sp.symbols("q", positive=True)
    polynomial=1-2*q+q**3
    qc=(sp.sqrt(5)-1)/2
    assert sp.simplify(polynomial.subs(q,qc))==0
    assert sp.simplify(qc+qc**3)<1

    reps={mod2(g):g for row in rows[:4] for g in row}
    keys=list(reps)
    def left(s):
        out=sp.zeros(6)
        for i,k in enumerate(keys):
            out[keys.index(mod2(mul(s,reps[k]))),i]=1
        return out
    h=sp.eye(6)-q*(left(A)+left(B))+q**3*left(delta)
    rhs=sp.eye(6)[:,keys.index(mod2(E))]
    assert sp.factor(h.det() - polynomial*(1+2*q-q**3)*(1-(q+q**3)**2)**2) == 0
    probs=(polynomial*h.inv()*rhs).applyfunc(sp.cancel)
    assert sp.cancel(sum(probs)-1)==0
    assert [sp.simplify(p.subs(q,qc)) for p in probs]==[sp.Rational(1,6)]*6
    # Symbolic generating functions also reproduce direct integer counts.
    for k,p in zip(keys,probs):
        co=sp.series(p/polynomial,q,0,9).removeO().expand()
        assert all(co.coeff(q,n)==sum(mod2(g)==k for g in row)
                   for n,row in enumerate(rows[:9]))
    print("Critical fugacity:",qc,"; Gibbs endpoint probabilities:",[sp.Rational(1,6)]*6)
    print("Subcritical fugacity 1/2, in e,a,b,(ab),(ba),(aba) order:",
          [probs[keys.index(mod2(g))].subs(q,sp.Rational(1,2))
           for g in (E,A,B,mul(A,B),mul(B,A),delta)])
    print("Squared critical norm of partial hidden a^2 b^(-2) insertion:",sp.expand(qc**2))
    print("Alternative representation (a^2 a)(b^2 a)^(-1) relative norm squared:",qc)
    # Test the two off-diagonal matrix elements of the isometry equation
    # at <ba|X*X|ab>. Only a* b is supported.
    ab,ba=mul(A,B),mul(B,A)
    test=[[int(mul(i,ba)==mul(j,ab)) for j in (A,B)] for i in (A,B)]
    assert test==[[0,1],[0,0]]
    print("Primitive isometry coefficient test at <ba|X*X|ab>:",test)
    print("Stationary connected range-projection variance:",sp.simplify(qc*(1-qc)))


if __name__=="__main__":
    main()
