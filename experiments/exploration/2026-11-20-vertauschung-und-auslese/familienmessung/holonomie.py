"""Shortest SO transports of unoriented simplex normals in a fixed Euclidean frame."""
from itertools import combinations
import sympy as s
from flaggen import transfers

def main():
    eye=s.eye(7)
    v=[eye[:,i]-s.ones(7,1)/7 for i in range(7)]
    norm=s.Rational(6,7)
    def step(i,j):
        # Rodrigues in the actual 2-plane. Carries the first unit normal
        # to minus the second, the shortest projective choice.
        k=(v[i]*v[j].T-v[j]*v[i].T)/norm
        return eye+k+s.Rational(6,7)*k*k
    origin=6
    q=v[origin]*v[origin].T/norm
    tangent=eye-s.ones(7)/7-q
    weighted=s.zeros(7)
    uniform=s.zeros(7)
    lineops=s.zeros(7)
    mats=transfers()
    one=s.ones(21,1)
    def word_prob(i,j):
        a,b,c=[mats[k+1] for k in (origin,i,j)]
        return (one.T*a*b*c*a*one)[0]/(21*8**4)
    law=sum((word_prob(i,j)*2 for i,j in combinations(range(6),2)),s.S.Zero)
    for i,j in combinations(range(6),2):
        R=step(j,origin)*step(i,j)*step(origin,i)
        assert R*R==eye and R.T==R and R*v[origin]==-v[origin]
        d=(eye-R)/2-q
        assert d.rank()==1 and d*d==d and tangent*d==d
        diff=v[i]-v[j]
        plus=tangent*(v[i]+v[j])
        assert d==diff*diff.T/diff.dot(diff)
        if i==0 and j==1:
            print("Tangent axis equals difference / projected sum:",
                  d==diff*diff.T/diff.dot(diff),d==plus*plus.T/plus.dot(plus))
        uniform+=d/15
        weighted+=d*(2*word_prob(i,j)/law)
        if ((i+1)^(j+1))==origin+1:
            lineops+=d
    assert uniform==tangent/5
    assert lineops.rank()==3 and lineops**2==lineops
    assert weighted==s.Rational(7,40)*tangent+lineops/24
    from apartment import twice_projector
    assert lineops==s.Matrix(twice_projector(origin+1))[1:,1:]/2
    print("Line triangles' axis sum spectrum:",lineops.eigenvals())
    print("All closing flag triangles, exact mean tangent axis projector:",weighted.eigenvals())
    print("Fifteen equally drawn directions, exact mean:",uniform.eigenvals())
    print("Closed triangle spin adjoint = one Clifford conjugation; SO minus-rank 2")

if __name__=="__main__":
    main()
