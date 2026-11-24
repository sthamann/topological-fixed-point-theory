"""Primitive positive unitary loop and positive-lift obstruction for triangles."""
from argparse import ArgumentParser
from pathlib import Path
import sympy as s
from rahmen import frame

def pf(M):
    return M[0,1]*M[2,3]-M[0,2]*M[1,3]+M[0,3]*M[1,2]

def digon(root):
    _,gamma,B=frame(root)
    eye=s.eye(4)
    g=gamma[0]
    P=(eye+g)/2
    z=s.symbols("z")
    U=eye+(z-1)*P
    Bt=U*B*U.T
    assert P*P==P and P.rank()==2 and P*B==B*P.T
    assert (Bt-(eye+(z*z-1)*P)*B).applyfunc(s.expand)==s.zeros(4)
    assert s.expand(pf(Bt)-z*z*pf(B))==0
    assert s.expand(U.det()-z*z)==0
    assert U.subs(z,-1)==-g and Bt.subs(z,-1)==B
    a=s.Rational(1,6)
    # Exact comparison: arccos(a) > pi/3 since a < 1/2; and < pi/2.
    assert 0<a<s.Rational(1,2)
    theta=s.acos(a)
    pm=(eye-g)/2
    a0,speed=s.symbols("phase_rate speed",real=True)
    Z=a0*eye+speed*g
    assert (Z*P-(a0+speed)*P).applyfunc(s.expand)==s.zeros(4)
    assert (Z*pm-(a0-speed)*pm).applyfunc(s.expand)==s.zeros(4)
    print("Digon native Pfaffian ratio:",s.expand(pf(Bt)/pf(B)),"; rank of positive generator:",P.rank())
    print("First opposite-line contact at",theta/s.pi,"of one turn")
    print("Triangle length between pi and 3*pi/2; min positive odd winding = 3")
    current=s.trace(P/2)/4/(2*s.pi)
    assert current==1/(8*s.pi)
    print("Phase-normalized stationary current density in this model:", current)
    print("Positive lift condition for either orientation: phase_rate >= abs(speed)")

def quarter(root):
    _,gamma,B=frame(root)
    I=s.eye(4)
    Y=s.Matrix([[0,-s.I],[s.I,0]])
    R=(s.eye(2)-s.I*Y)/s.sqrt(2)
    Q=s.kronecker_product(R,R)
    assert Q*B*Q.T==B and Q.det()==1
    x=s.kronecker_product(Y,Y)
    assert any(g==x or g==-x for g in gamma)
    assert Q**2==-x
    P=(I+x)/2
    # On the active two dimensional plane Q has eigenvalues i,-i;
    # outside it Q is the identity. Its two eigenlines are Kramers paired.
    J=s.expand((Q-Q.H)/(2*s.I))
    p_plus=(P+J)/2
    p_minus=(P-J)/2
    K=p_plus/2+3*p_minus/2
    assert p_plus.rank()==p_minus.rank()==1 and p_plus*p_minus==s.zeros(4)
    assert K.eigenvals()=={s.S.Zero:2,s.Rational(1,2):1,s.Rational(3,2):1}
    assert B*K.T*B.H+K==2*P
    # z = exp(i alpha/2); alpha = pi is the primitive endpoint.
    z=s.symbols("z")
    U=I-P+z*p_plus+z**3*p_minus
    bpath=U*B*U.T
    assert (bpath-(I+(z**4-1)*P)*B).applyfunc(s.expand)==s.zeros(4)
    assert U.subs(z,s.I)==Q
    assert s.expand(pf(bpath)-z**4*pf(B))==0
    print("Native signed quarter reached as a primitive endpoint; generator spectrum",K.eigenvals())
    print("Midpath quarter is excluded: its quaternionic normal is already back at the reference line.")

if __name__=="__main__":
    p=ArgumentParser()
    p.add_argument("--repo",type=Path,default=Path(__file__).resolve().parents[4])
    args=p.parse_args()
    digon(args.repo)
    quarter(args.repo)
