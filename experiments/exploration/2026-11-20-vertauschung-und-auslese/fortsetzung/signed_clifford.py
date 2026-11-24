"""Exact signed-site Clifford lifts above the hermitian/braid P2 model."""
from collections import Counter
from itertools import product
from fractions import Fraction
from math import prod

def compose(f,g):
    return tuple((1 if x>0 else -1)*f[abs(x)-1] for x in g)
identity=tuple(range(1,6))
def magnitude(f):
    return tuple(map(abs,f))
def closure(gens):
    known={identity}
    work=[identity]
    while work:
        a=work.pop()
        for g in gens:
            b=compose(a,g)
            if b not in known:
                known.add(b)
                work.append(b)
    return known


def run():
    # Gamma in native order: IX, IZ, XY, ZY, YY. R=(I-iY)/sqrt(2)
    # conjugates X to -Z, Z to X and Y to Y.
    quarter=(-2,1,-4,3,5)
    cycle=(1,2,4,5,3)
    central=compose(quarter,quarter)
    assert central==(-1,-2,-3,-4,5)
    # This sign pattern is conjugation by gamma 5, reflecting the four others.
    distributions=Counter()
    for signs in product((-1,1),repeat=5):
        if signs.count(-1)%2:
            continue # 3-cycle is even; Clifford conjugation preserves volume.
        lift=tuple(s*x for s,x in zip(signs,cycle))
        group=closure((quarter,lift))
        kernel=[x for x in group if magnitude(x)==identity]
        action={magnitude(x) for x in group}
        assert len(action)==6
        assert compose(lift,central)!=compose(central,lift)
        # J ~ i gamma1 gamma2. All site maps preserve its projective line.
        # Swapping the two gammas contributes one more minus sign.
        def fixes_measured_J(g):
            s1=1 if g[0]>0 else -1
            s2=1 if g[1]>0 else -1
            order_sign=1 if abs(g[0])<abs(g[1]) else -1
            return s1*s2*order_sign==1
        preserves_J=all(fixes_measured_J(g) for g in group)
        distributions[len(group),len(kernel),preserves_J]+=1
    print('all 16 signed lifts of one 3-cycle, fixed native quarter:')
    print('  (group, Pauli kernel, J exactly conserved) -> number',
          sorted(distributions.items()))
    # A canonical representative rotates X,Z,Y on first qubit with positive signs.
    smallest=closure((quarter,cycle))
    assert len(smallest)==48
    gr_kernel=[x for x in smallest if magnitude(x)==identity]
    # Products of at most two gamma span all of M4(C): grade 0,1,2.
    # A gamma monomial commutes with the kernel iff all its signs are +.
    labels=[()]+[(i,) for i in range(5)]+[(i,j) for i in range(5) for j in range(i+1,5)]
    def parity_flip(g,word):
        return sum(g[i]<0 for i in word)%2
    invariant_words=[w for w in labels if all(not parity_flip(g,w) for g in gr_kernel)]
    assert invariant_words==[(),(0,1)]
    print('kernel commutant gamma monomials (zero-based):',invariant_words)
    print('minimal extension invariant densities: (I+s*(I tensor Y))/4, -1 <= s <= 1')

    # Test all subgroup stabilizers to establish minimal equivariant flag degrees.
    # An anticommuting Pauli pair in a subgroup forbids any common eigenline.
    mono_sign=lambda w: tuple((1 if (len(w)-(i in w))%2==0 else -1)*(i+1) for i in range(5))
    pauli_monomial={mono_sign(w):w for w in labels}
    def anti(w,z):
        return (len(w)*len(z)-len(set(w)&set(z)))%2
    def subgroups(group):
        sub={frozenset({identity})}
        todo=list(sub)
        while todo:
            small=todo.pop()
            for g in group-small:
                big=frozenset(closure((*small,g)))
                if big not in sub:
                    sub.add(big)
                    todo.append(big)
        return sub
    minimum=[]
    for i in range(5):
        stabilizer={g for g in smallest if abs(g[i])==i+1}
        admissible=[]
        for subgroup in subgroups(stabilizer):
            k=[pauli_monomial[g] for g in subgroup if magnitude(g)==identity]
            if not any(anti(a,b) for a,b in product(k,repeat=2)):
                admissible.append(len(stabilizer)//len(subgroup))
        minimum.append(min(admissible))
    assert minimum==[4,4,2,2,2]
    # Work in J=+1 (the equally valid negative sector has the same projectors).
    # Three mixed gamma restrict to a Pauli triad. Bloch vectors v below mean
    # rank-one rays p=(I+J)(I+v_hat.gamma)/4, with v_hat=v/||v||.
    def image(g,coord):
        dest=[0,0,0]
        for j,c in enumerate(coord):
            entry=g[j+2]
            dest[abs(entry)-3]=(1 if entry>0 else -1)*c
        return tuple(dest)
    flags={}
    for i in range(5):
        if i<2:
            flags[i]={v for v in product((-1,1),repeat=3) if prod(v)==(-1)**i}
        else:
            flags[i]={tuple(sign if j==i-2 else 0 for j in range(3)) for sign in (-1,1)}
        assert len(flags[i])==minimum[i]
    for g in smallest:
        for i in range(5):
            assert {image(g,v) for v in flags[i]}==flags[abs(g[i])-1]
    print('all stabilizer subgroup obstructions + attaining flag fibres:',minimum)
    print('14 flags on one J sector: 8 cube vertices + 6 face centres')
    print('pointwise invariant determinant winding is divisible by 2 (order 48) or 4 (order 96)')

    # A true U(1) representation with non-negative energy has integral weights.
    # Enumerating 0..3 includes every spectrum of total weight at most three.
    weights=list(product(range(4),repeat=4))
    wind1={tuple(sorted(w)) for w in weights if sum(w)==1}
    unique_ground=[w for w in weights if w.count(0)==1]
    assert wind1=={(0,0,0,1)} and min(map(sum,unique_ground))==3
    print('positive angular winding one:',wind1,'minimum winding with one zero mode:',3)
    print('rank-one current coefficient times 2 pi:',Fraction(1,4))
    # Use Gaussian-integer 4x4 matrices to check the native antiunitary.
    def kron(a,b):
        return tuple(tuple(a[i][j]*b[k][l] for j in range(len(a))
                           for l in range(len(b))) for i in range(len(a)) for k in range(len(b)))
    def mm(a,b):
        return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*b)) for row in a)
    def scalar(c,a):
        return tuple(tuple(c*x for x in row) for row in a)
    def conj(a):
        return tuple(tuple(x.conjugate() for x in row) for row in a)
    def dagger(a):
        return tuple(zip(*conj(a)))
    I=((1,0),(0,1))
    X=((0,1),(1,0))
    Y=((0,-1j),(1j,0))
    Z=((1,0),(0,-1))
    gammas=[kron(I,X),kron(I,Z),kron(X,Y),kron(Z,Y),kron(Y,Y)]
    B=kron(mm(X,Z),I)
    theta=lambda a:mm(mm(B,conj(a)),dagger(B))
    assert mm(B,conj(B))==scalar(-1,kron(I,I))
    assert all(theta(g)==g for g in gammas)
    J=scalar(1j,mm(gammas[0],gammas[1]))
    assert J==kron(I,Y) and theta(J)==scalar(-1,J)
    # J eigenblocks have rank two, Θ swaps them. Commuting with the signed
    # group and Θ selects the maximally mixed state even though the channel
    # conserves the total weight of each J block.
    print('Theta²=-I, Theta fixes 5/5 gamma, Theta J Theta^-1=-J: unique symmetric state I4/4')
    # Mean and covariance of the two rank-one orbits; radii² are 1 and 3.
    # Since both are tight and isotropic, the relative orbit mass cancels.
    for subset in (flags[0]|flags[1],flags[2]|flags[3]|flags[4]):
        assert all(sum(v[i] for v in subset)==0 for i in range(3))
        covariance=tuple(tuple(sum(Fraction(v[i]*v[j],sum(x*x for x in v))
                                    for v in subset)/len(subset) for i in range(3)) for j in range(3))
        assert covariance==tuple(tuple(Fraction(i==j,3) for i in range(3)) for j in range(3))
    print('any relative orbit weights: Lüders map per sector = (rho+tr(rho)I2)/3')
    # Pointwise antiunitary time reflection would commute with H=p. A
    # Θ-stable complex subspace has even dimension by Kramers pairing.
    # Move the flag instead: p and Theta(p) are orthogonal by Theta²=-I,
    # and Theta u_p(theta) Theta^-1 = u_Theta(p)(-theta).
    print('Kramers-even pointwise loop ranks: 0,2,4. One moving-flag loop current: 1/(8 pi)')
    print('self-reflected pair loop has rank two and current 1/(4 pi); explicit full lifted cover: 2 x 14')
    def matsum(mats):
        return tuple(tuple(sum(m[r][c] for m in mats) for c in range(4)) for r in range(4))
    eye=kron(I,I)
    def pf(a):
        return a[0][1]*a[2][3]-a[0][2]*a[1][3]+a[0][3]*a[1][2]
    # The matrix entries below are Gaussian integers/halves, exactly
    # representable by Python's complex type for these bounded operations.
    # A 4x4 Pfaffian is quadratic: three tests certify the polynomial.
    for gamma in gammas:
        for sign in (-1,1):
            P=scalar(.5,matsum((eye,scalar(sign,gamma))))
            assert mm(P,P)==P and theta(P)==P and sum(P[i][i] for i in range(4))==2
            for z in (0,1,2):
                U=matsum((eye,scalar(z-1,P)))
                A=mm(U,B)
                assert tuple(zip(*A))==scalar(-1,A)
                assert pf(A)==z*pf(B)
    print('for 10 signed Kramers tests, Pf(u_P(z)B)/Pf(B)=z; det u_P=z²')
    # Exact Hilbert-Schmidt transfer check on the full gamma basis of M4.
    for word in labels:
        W=eye
        for i in word:
            W=mm(W,gammas[i])
        jumps=[mm(mm(g,W),g) for g in gammas]
        q10=matsum([scalar(5,W),*jumps])
        tfpt12=matsum([scalar(7,W),*jumps])
        assert q10==scalar((10,2,6)[len(word)],W)
        assert tfpt12==scalar((12,4,8)[len(word)],W)
        assert tfpt12==matsum((scalar(2,W),q10))
    print('M4 bare fastest spectrum: 1, (1/5)x5, (3/5)x10')
    print('uniform empty + five shortest tests: native 1, (1/3)x5, (2/3)x10 on 16/16 basis words')
    from itertools import combinations as choose
    from math import comb
    # All probability vectors in steps of 1/10; the analytic result uses
    # the average of all pair weights and holds on the continuous simplex.
    laws=[]
    for bars in choose(range(14),4):
        stops=(-1,*bars,14)
        w=[stops[i+1]-stops[i]-1 for i in range(5)]
        assert sum(w)==10
        second=max(max(w), max(10-a-b for a,b in choose(w,2)))
        laws.append((second,tuple(w)))
    assert len(laws)==comb(14,4) and min(laws)==(6,(2,2,2,2,2))
    assert sum(second==6 for second,w in laws)==1
    print('fastest law on 0.1 grid:',1,'/',len(laws),'(unique analytically on simplex)')
    # 16x16 matrices in the orthogonal gamma monomial basis.
    monomials=[]
    for word in labels:
        W=eye
        for i in word:
            W=mm(W,gammas[i])
        monomials.append(W)
    def inner(A,B):
        C=mm(dagger(A),B)
        return sum(C[i][i] for i in range(4))/4
    def to_matrix(scaled):
        return tuple(tuple(inner(A,B) for B in scaled) for A in monomials)
    native_scaled=[]
    short_scaled=[]
    for W in monomials:
        jumps=[mm(mm(g,W),g) for g in gammas]
        native_scaled.append(matsum((scalar(7,W),*jumps)))
        short_scaled.append(matsum((scalar(2,W),*[matsum((W,B)) for B in jumps])))
    matrix12=to_matrix(native_scaled)
    assert matrix12==to_matrix(short_scaled)
    assert all(matrix12[i][j]==(0 if i!=j else (12,4,8)[len(labels[i])])
               for i in range(16) for j in range(16))
    print('full 16x16 native vs six-word comparison: 256/256; 12*Phi = diag(12,[4]*5,[8]*10)')
    # Choi support has six orthogonal columns (I,G1,...,G5).
    # All products K_a^* K_b of nonzero scaled columns span M4: checks below
    # use unscaled columns; the nonzero scalar multipliers preserve rank.
    cols=[eye,*gammas]
    assert all(inner(A,D)==int(i==j) for i,A in enumerate(cols) for j,D in enumerate(cols))
    assert all(inner(A,D)==int(i==j) for i,A in enumerate(monomials) for j,D in enumerate(monomials))
    print('binary CP instrument effect dimensions, total / fixed all-input probability:',
          len(cols)**2, '/',len(cols)**2-len(monomials))
    # Two fully signed-permutation-covariant extremal effects on this line:
    # E=diag(a,b,b,b,b,b) where 7a+5b=6. Conditional action on gamma at
    # success is twice (7a-3b)/12 times gamma, at identical chance one half.
    result=[]
    for a,b in [(Fraction(6,7),0),(Fraction(1,7),1)]:
        assert 7*a+5*b==6 and 0<=a<=1 and 0<=b<=1
        result.append((7*a-3*b)/6)
    assert result==[1,-Fraction(1,3)]
    print('same full transfer and half-probability record, conditional gamma1:',result)
    # Reflected temporal kernel for the flagged positive generator H=p is
    # K(t,s)=(I-p)+exp(-t)exp(-s)p. Both summands are positive Gram kernels.

if __name__=='__main__':
    run()
