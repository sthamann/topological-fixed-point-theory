#!/usr/bin/env python3
"""Small exact mu8 operator-word checks; no GPU, truncation, or new state fit."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path

Q = (1,)*8
A = (1,)*4+(0,)*4
E0 = (1,)+(0,)*7
ZERO = (0,)*8


def add(a,b): return tuple(x+y for x,y in zip(a,b))
def neg(a): return tuple(-x for x in a)
def bilinear(a,b): return sum(a[i]*b[j] for i in range(8) for j in range(i))
def cyclic(n): return n[-1:]+n[:-1]
def swap(n): return n[4:]+n[:4]


def compose(*operators):
    """Usual left-to-right written operator product, applied rightmost first."""
    def result(basis):
        phase=0
        for op in reversed(operators):
            basis,p=op(basis)
            phase+=p
        return basis,phase%8
    return result


def scalar(exponent):
    return lambda basis:(basis,exponent%8)


class Model:
    def __init__(self,chi,zeta=0,nu_exp=2):
        assert chi in (1,3) and zeta in range(4) and nu_exp in range(8)
        self.chi,self.zeta,self.r=chi,zeta,nu_exp
        self.kappa=(2*nu_exp)%4

    def E(self,v):
        return lambda basis:((add(v,basis[0]),basis[1]),4*bilinear(v,basis[0])%8)

    def E_inverse(self,v):
        return compose(scalar(4*bilinear(v,v)),self.E(neg(v)))

    def X(self,basis):
        n,c=basis
        # The last cocycle factor is essential in the E_n X^c normal form.
        return (add(n,Q) if c else n,1-c),(2*self.chi*sum(n)+2*self.zeta*c+4*c*bilinear(n,Q))%8

    def X_inverse(self,basis):
        return compose(scalar(-2*self.zeta),self.E(neg(Q)),self.X)(basis)

    def C(self,basis):
        n,c=basis
        return (cyclic(n),c),(4*n[7]*sum(n[:7])+self.r*c)%8

    def C_inverse(self,basis):
        n,c=basis
        old=n[1:]+n[:1]
        return (old,c),(-4*old[7]*sum(old[:7])-self.r*c)%8

    def P(self,basis):
        n,c=basis
        return (swap(n),c),(4*sum(n[:4])*sum(n[4:])+2*self.kappa*c)%8

    def Z(self,u):
        return lambda basis:(basis,(4*sum(x*y for x,y in zip(u,basis[0]))+2*basis[1]*sum(u))%8)

    def D(self,basis): return basis,4*basis[1]
    def parity(self,basis): return basis,4*(sum(basis[0])%2)
    def sector_phase(self,r): return lambda basis:(basis,(r*basis[1])%8)
    def W1(self,basis): return compose(scalar(4),self.Z(E0),self.E(E0),self.C)(basis)
    def W4(self,basis): return compose(self.Z(A),self.E(A),self.P)(basis)


def power(op,n): return compose(*([op]*n))


def kinetic(basis,alpha):
    n,c=basis
    p=[Fraction(x)+Fraction(c,2) for x in n]
    delta=[Fraction(j,8)+alpha/8-Fraction(1,2) for j in range(8)]
    return sum(x*x/2+d*x for x,d in zip(p,delta))


def run_checks():
    charges=[ZERO,Q,tuple(1 if j<4 else -1 for j in range(8)),
             (3,-2,0,1,-4,2,5,-1)]
    charges += [tuple(sign*int(i==j) for i in range(8)) for j in range(8) for sign in (1,-1)]
    bases=[(n,c) for n in charges for c in (0,1)]
    comparisons=0
    def equal(left,right,basis):
        nonlocal comparisons
        comparisons+=1
        assert left(basis)==right(basis),(basis,left(basis),right(basis))
    identity=scalar(0)
    summaries=[]
    for chi in (1,3):
        for zeta in (0,2):
            for r in range(8):
                m=Model(chi,zeta,r)
                for basis in bases:
                    equal(compose(m.X,m.X_inverse),identity,basis)
                    equal(power(m.X,2),compose(scalar(2*zeta),m.E(Q)),basis)
                    equal(power(m.C,4),m.P,basis)
                    equal(power(m.C,8),identity,basis)
                    equal(power(m.P,2),identity,basis)
                    equal(compose(m.C,m.C_inverse),identity,basis)
                    equal(power(m.W1,4),m.W4,basis)
                    equal(power(m.W1,8),compose(m.parity,m.E(Q)),basis)
                    equal(power(m.Z(E0),2),m.D,basis)
                    equal(power(m.D,2),identity,basis)
                    equal(compose(m.D,m.X),compose(scalar(4),m.X,m.D),basis)
                    # Ad C(X)=X*(nu P_even - nu^-1 P_odd).
                    right=compose(m.X,lambda b:(b,r if b[1]==0 else (4-r)%8))
                    equal(compose(m.C,m.X,m.C_inverse),right,basis)
                    # Minimal ratio exposes the already present grade, with
                    # eigenvalues nu and -nu^-1 in the two input sectors.
                    ratio=compose(m.X_inverse,m.C,m.X,m.C_inverse)
                    equal(ratio,lambda b:(b,r if b[1]==0 else (4-r)%8),basis)
                    # Parity dressing is a source-fixing chi+ <-> chi- map.
                    other=Model(4-chi,zeta,r)
                    equal(m.X,compose(m.parity,other.X),basis)
                    for j in range(8):
                        e=tuple(int(i==j) for i in range(8))
                        equal(compose(m.X,m.E(e),m.X_inverse,m.E_inverse(e)),scalar(2*chi),basis)
                        equal(compose(m.C,m.E(e),m.C_inverse),m.E(cyclic(e)),basis)
                        equal(compose(m.P,m.E(e),m.P),m.E(swap(e)),basis)
                    equal(compose(m.P,m.X,m.P,m.X_inverse),scalar(2*m.kappa),basis)
                    # Scalar C covariance holds exactly in r=2,6.
                    direct=compose(m.C,m.X,m.C_inverse)(basis)
                    scalar_case=compose(scalar(r),m.X)(basis)
                    assert (direct==scalar_case)==(basis[1]==0 or r in (2,6))
                    for alpha in (Fraction(1,4),Fraction(3,4)):
                        n,c=basis
                        charge=sum(n)+4*c
                        image,_=m.W1(basis)
                        assert kinetic(image,alpha)-kinetic(basis,alpha)==(charge+alpha)/8
                        image4,_=m.W4(basis)
                        delta=kinetic(image4,alpha)-kinetic(basis,alpha)
                        assert delta==Fraction(charge,2)+Fraction(3,4)+alpha/2
                        expected=(6+int(4*alpha)+4*(charge%2))%8
                        assert (8*delta)%8==expected
                        comparisons+=3
                summaries.append(dict(chi=chi,zeta=zeta,nu_mu8_exp=r,kappa_mu4_exp=m.kappa,
                                      scalar_C_covariance=r in (2,6)))
    # Minimal explicit failure of the X action with the missing eps(n,q)^c.
    n=(0,1,0,0,0,0,0,0)
    correct=power(Model(1).X,2)((n,0))
    wrong=((add(n,Q),0),(4*sum(n))%8)
    assert correct!=wrong
    # The new source-invisible sector phase does not change the old source.
    base=Model(1,0,2)
    counter=Model(1,0,1)
    for basis in bases:
        equal(counter.C,compose(base.sector_phase(-1),base.C),basis)
        equal(counter.P,compose(base.D,base.P),basis)
    return dict(verdict="PASS_SCOPED_EXACT_MU8_OPERATOR_CHECKS",comparisons=comparisons,
                basis_states=len(bases),parameter_cases=len(summaries),cases=summaries,
                results={"chi_parity_dressing":True,"kappa_two_T1_counterconstruction":True,
                         "W1_fourth_W4":True,"W1_eighth_parity_Eq":True,
                         "D_equals_Ze0_squared":True,"scalar_C_requires_kappa_zero":True,
                         "D_reconstructed_from_C_X_ratio":True,
                         "new_H_or_state_parameters":False},
                scope="finite exact word checks plus accompanying all-charge algebraic proof; no physical selection")


if __name__=="__main__":
    result=run_checks()
    result["checker_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(__file__).with_name("results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="cases"},indent=2))
