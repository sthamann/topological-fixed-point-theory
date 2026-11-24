"""A number-preserving internal phase reference, and its sharp common-state bound."""
from argparse import ArgumentParser
from itertools import combinations
from pathlib import Path
import sympy as s
from rahmen import frame
from interferenz import simplex,simplify

def main(root):
    _,gamma,B=frame(root)
    _,dirs=simplex()
    eye=s.eye(4)
    pairs=list(combinations(range(4),2))
    wedge=lambda U:s.Matrix([[U[i,k]*U[j,l]-U[i,l]*U[j,k]
                              for k,l in pairs] for i,j in pairs])
    coefficients=lambda M:s.Matrix([M[i,j] for i,j in pairs])
    blank=coefficients(B)/s.sqrt(2)
    refs=[]
    differences=[]
    z=s.symbols('z')
    for j in range(6):
        g=sum((k*h for k,h in zip(dirs[:,j],gamma)),s.zeros(4))
        P=(eye+g)/2
        active=coefficients(simplify(P*B))
        idle=coefficients(simplify((eye-P)*B))
        assert simplify(active.H*active)==simplify(idle.H*idle)==s.ones(1)
        assert simplify(active.H*idle)==s.zeros(1)
        assert simplify((active+idle)/s.sqrt(2))==blank
        # In exterior degree four retain a fully occupied active bank
        # and a fully occupied idle bank, both with the same total charge.
        U=eye+(z-1)*P
        assert s.factor(U.det())==z*z
        refs.append(idle)
        differences.append(simplify((active-idle)/s.sqrt(2)))
        if j==0:
            U=eye+(z-1)*P
            assert simplify(wedge(U)*active-z*z*active)==s.zeros(6,1)
            assert simplify(wedge(U)*idle-idle)==s.zeros(6,1)
            w=s.symbols('w')
            assert simplify(wedge(w*U)-w*w*wedge(U))==s.zeros(6)
    refmat=s.Matrix.hstack(*refs)
    assert simplify(refmat.H*refmat)==3*s.eye(6)/5+2*s.ones(6)/5
    average=simplify(refmat*refmat.H/6)
    assert average.eigenvals()=={s.Rational(1,2):1,s.Rational(1,10):5}
    assert all(simplify((v.H*blank).conjugate()*(v.H*blank))==s.Matrix([[s.Rational(1,2)]])
               for v in refs)
    tangent_sum=sum((d for d in differences),s.zeros(6,1))
    assert simplify(tangent_sum)==s.zeros(6,1)
    tangent_average=simplify(sum((d*d.H for d in differences),s.zeros(6))/6)
    assert simplify(tangent_average-(s.eye(6)-blank*blank.H)/5)==s.zeros(6)
    print("Blind density after a common blank: cos(alpha)^2 on blank, sin(alpha)^2/5 on the other 5")
    print("Internal two-particle interferometer reads relative phase z^2 and cancels a common phase")
    print("Common unit spectator absent; best possible worst-branch overlap probability:",s.Rational(1,2))
    print("Reference Gram spectrum",s.simplify((refmat.H*refmat).trace()),":",
          {s.S(3):1,s.Rational(3,5):5})
    print("Common reference with two four-mode banks at fixed total occupation four: phase z^2 in one arm, 1 in the other")
    print("Universal reference trial saturating the bound: normalized bivector of the native B")

if __name__=="__main__":
    p=ArgumentParser()
    p.add_argument("--repo",type=Path,default=Path(__file__).resolve().parents[4])
    main(p.parse_args().repo)
