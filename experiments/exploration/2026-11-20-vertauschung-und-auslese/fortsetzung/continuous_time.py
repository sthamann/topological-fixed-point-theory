"""Stationary Euclidean compressions and completely positive intermediate maps.

The moment inequalities and reduction proof are in VERDICT.md. These checks use
exact constants and two explicit CP families as controls.
"""
import sympy as sp

a = sp.Rational(1,3)
b = sp.Rational(2,3)


def pauli_weights(v,w):
    # Per-word normalized Choi weights in orbit sizes 1, 5 and 10.
    return ((1+5*v+10*w)/16, (1-3*v+2*w)/16, (1+v-2*w)/16)


def checked_bound():
    half=pauli_weights(sp.sqrt(a),sp.sqrt(b))[2]
    assert sp.ask(sp.Q.negative(half))
    one_half=(1+sp.sqrt(a))/2
    gap_half=sp.simplify(b**3/one_half**2-b**2)
    assert sp.simplify(gap_half-(28-16*sp.sqrt(3))/9)==0
    s=sp.symbols('s',positive=True)
    log_ratio=sp.log(2*b/(1+a**s))
    slope=sp.limit(log_ratio/(1-s),s,1,dir='-')
    gap_all=sp.simplify(b*sp.exp(slope)-b**2)
    assert sp.simplify(gap_all - (b*3**(-sp.Rational(1,4))-b**2))==0
    assert 0 < gap_half < gap_all < sp.Rational(1,9)
    print("Native half step normalized grade-two Choi weight:",half)
    print("Half-step-only lower bound on two-step error:",gap_half, float(gap_half))
    print("All intermediate times lower bound:",gap_all, float(gap_all))


def continuous_families():
    x=sp.symbols('x',positive=True)
    weights=tuple(sp.factor(v) for v in pauli_weights(x,(1+x)/2))
    assert all(sp.expand(u-v)==0 for u,v in zip(weights,((5*x+3)/8,(1-x)/8,0)))
    # Choose x=3^(-t). Dilation spectral measures: delta(log3) in grade 1,
    # half delta(0)+half delta(log3) in grade 2, delta(0) on identity.
    assert sp.simplify((1+a*a)/2-b*b)==sp.Rational(1,9)
    print("All-times CPTP / Euclidean Gram family, grade-two return:",(1+a*a)/2,
          "; two-step excess:",sp.Rational(1,9))

    # Five static classical memories, each selects one Lüders dephasing.
    # Per memory e^(-tL)=D_j+6^(-t)(I-D_j) is a reversible quantum semigroup.
    for t,expected in [(1,(a,b)),(2,(sp.Rational(2,9),sp.Rational(11,18)))]:
        prob=sp.Rational(1,6)**t
        rates=((1+4*prob)/5,(3+2*prob)/5)
        assert rates==expected
    print("Static five-choice reversible CP dilation, two-step response: 2/9,11/18")
    # Randomly phased periodic fresh environments retain all integer powers;
    # scalar interpolation between successive integers breaks continuous RP.
    lam=a
    spline=lambda t: 1-t+t*lam
    t=sp.Rational(1,4)
    defect=sp.simplify(spline(2*t)-spline(t)**2)
    assert defect==-sp.Rational(1,36)
    print("Fresh-choice suspension, scalar reflected 2x2 determinant:",defect)

    # Natural two-jump relaxation on the Burau block: orthogonal directions
    # to the two rank-one Hilbert-space projections are represented as Bloch axes.
    # Their squared Hilbert overlap is q/(1+q)^2; dot product is 2*overlap-1.
    q=sp.symbols('q',positive=True)
    overlap=q/(1+q)**2
    assert sp.factor(sp.Rational(1,4)-overlap)==(q-1)**2/(4*(q+1)**2)
    # Equal jump rates give the lowest middle/slow exponent ratio. From
    # |n.m|>=1/2 this is at least three, already above the native ratio.
    assert sp.Rational(3,2)**3 > 3
    print("Two local balanced/nonbalanced Hermitian jumps: decay ratio >=3;",
          "native log(3)/log(3/2)<3")


if __name__=="__main__":
    checked_bound()
    continuous_families()
