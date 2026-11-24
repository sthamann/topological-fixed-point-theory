"""Shortest normal transport of seven calibrated channels and its lifting obstruction."""
from itertools import combinations, product
from fractions import Fraction
from flaggen import transfers
import sympy as s

def main():
    gram=s.Matrix(7,7,lambda i,j:1 if i==j else -s.Rational(1,6))
    assert gram.rank()==6 and gram*gram==s.Rational(7,6)*gram
    triples=list(combinations(range(7),3))
    for signs in product((1,-1),repeat=6):
        gauge=(1,*signs)
        entry=lambda i,j:gauge[i]*gram[i,j]*gauge[j]
        assert all(entry(i,j)*entry(j,k)*entry(k,i)==-s.Rational(1,216)
                   for i,j,k in triples)
    print("Shortest projective triangle holonomy: -1 on all 35/35 triples, all 64 gauges")
    print("Fano line triples / other triples:",sum((i+1)^(j+1)==k+1 for i,j,k in triples),28)
    # Probability of reading the same label after exactly 3 changes, with
    # the three labels distinct. The chronological fourth event closes the path.
    m=transfers()
    v=s.ones(21,1)
    def prob(a,b,c):
        return (v.T*m[a]*m[b]*m[c]*m[a]*v)[0]/(21*8**4)
    line=42*prob(1,2,3)
    nonline=168*prob(1,2,4)
    print("Four-event closed triangles, line and nonline probability:",line,nonline)
    print("Probability conditional on a triangle that it is incident:",line/(line+nonline))
    # Pf on real unit normals is fixed. Closing e^(i phi)*B changes Pf by
    # e^(2 i phi); an antipodal lifted endpoint requires an odd half-turn.
    turns=[Fraction(2*k+1,2) for k in range(-4,4)]
    assert {2*t for t in turns if abs(t)==min(map(abs,turns))}=={-1,1}

if __name__=="__main__":
    main()
