"""Stationary opposite-flag histories. All counts before normalization are integers."""
from pathlib import Path
import sys
from itertools import permutations
import sympy as s
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"fortsetzung"))
from apartment import FLAGS, POINTS

def transfers():
    mats={u:s.zeros(21) for u in POINTS}
    for i,(p,line) in enumerate(FLAGS):
        for j,(q,other) in enumerate(FLAGS):
            if p not in other and q not in line:
                r,=set(line)&set(other)
                mats[p^q^r][i,j]=1
    assert all(m==m.T for m in mats.values())
    assert sum(mats.values(),s.zeros(21))*s.ones(21,1)==8*s.ones(21,1)
    return mats

def main():
    mats=transfers()
    one=s.ones(21,1)
    def prob(word):
        v=one
        for u in reversed(word):
            v=mats[u]*v
        return sum(v)/(21*8**len(word))
    same=7*prob([1,1])
    assert same==s.Rational(1,4) and prob([1,2])==s.Rational(1,56)
    print("repeat: independently resampled 1/7, shared endpoint",same,", frozen 1")
    for a,b in [(1,1),(1,2)]:
        row=[prob([a,b,c])/prob([a,b]) for c in POINTS]
        assert sum(row)==1
        print("next label after",a,b,row)
    lines=set()
    for a,b in permutations(POINTS,2):
        history=prob([a,b])
        inferred=[c for c in POINTS if c not in (a,b) and prob([a,b,c])/history==s.Rational(1,6)]
        assert inferred==[a^b]
        assert all(prob([a,b,c])/history==s.Rational(1,8) for c in POINTS if c not in (a,b,a^b))
        lines.add(tuple(sorted((a,b,*inferred))))
    assert len(lines)==7
    tables={frozenset(tuple(sorted(p[k-1] for k in line)) for line in lines) for p in permutations(POINTS)}
    assert len(tables)==30
    print("labeled Fano geometries matching all distinct triples: 1 /",len(tables))
    basis=one
    ranks=[]
    while True:
        candidates=s.Matrix.hstack(basis,*(m*basis for m in mats.values()))
        pivots=candidates.rref()[1]
        ranks.append(len(pivots))
        if len(pivots)==basis.cols: break
        basis=candidates[:,list(pivots)]
    print("prefix/suffix span starting at length zero", [1]+ranks)
    assert ranks==[7,21,21]
    print("nonzero Gram determinant for independent columns:",(basis.T*basis).det())
    # Longest-event return has constant chance 1/8 and does not bias the
    # endpoint if the *target flag* is retained. Reset instead retains source.
    assert prob([2,1,2])/prob([2,1]) == s.Rational(1,12)
    assert prob([1,1,2])/prob([1,1]) == s.Rational(1,8)

if __name__=="__main__":
    main()
