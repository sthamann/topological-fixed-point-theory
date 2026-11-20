"""Two original trine response axes certify retention of one quaternionic normal."""
from pathlib import Path
import argparse
import sympy as s
from spintransport import native,load_bridge

def main(root):
    phi,gamma=native(root)
    br=load_bridge(root)
    src,_=br.inherited(root)
    gen=br.generators()
    F=br.monomial(src.FSIG,gen)
    A=s.I*br.monomial(src.A_BIT,gen)
    assert F*F==A*A==s.eye(4) and F*A==-A*F
    ident=next(i for i,g in enumerate(gamma) if F==g or F==-g)
    partner=next(j for j,g in enumerate(gamma) if j!=ident and
                 any(A==sign*s.I*F*g for sign in (-1,1)))
    assert phi*F.vec()==F.vec()/3 and phi*A.vec()==2*A.vec()/3
    # In a six-dimensional orthonormal frame take 0=reference normal,
    # 1=the axis x shared by F and A, and 2=the other axis y.
    # This identity works for all orientations of the other six labels.
    a,b,c=s.symbols('a b c',real=True)
    FF=(2-a*a-b*b)/3
    AA=(2-b*b-c*c)/3
    deficit=s.factor(FF-s.Rational(1,3)+s.Rational(2,3)-AA)
    assert s.expand(deficit-(1-a**2+c**2)/3)==0
    gap=s.simplify(deficit.subs({a:-s.Rational(1,6),c:0}))
    assert gap==s.Rational(35,108)
    print("Original F index and second gamma of A, zero based:",ident,partner)
    print("Native response deficit:",0,"; wrong-vertex minimum:",gap)
    print("Uniform missing label deficit:",s.Rational(1,3))
    # Original E_3=(I-F)/3 on prepared (I+F)/4. Calibrations are examples,
    # not information present in the trine source.
    for values in ((1,1,1,-1,-1,-1),(-9,-4,-1,2,5,7)):
        assert sum(values)==0
        norm=sum(x*x for x in values)
        # Each nonreference simplex vertex has coordinate -1/6 along n0
        # and sqrt(7/6)*x_i along the null-sum unit contrast x.
        probs=[s.Rational(2,9)]+[(1+s.Rational(1,36)+s.Rational(7*v*v,6*norm))/9
                                   for v in values]
        print("Original E3 probability levels with tangent contrast",values,":",len(set(probs)),probs)
    assert len(set(probs))==7

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--repo",type=Path,default=Path(__file__).resolve().parents[4])
    main(parser.parse_args().repo)
