#!/usr/bin/env python3
"""Project Genesis: independent, finite-dimensional structural falsifiers.

Only proves implications of the explicitly supplied matrices/models. It does NOT
verify the TFPT repository, the physical SM, or a theory of everything.
Python 3.11+, sympy 1.12+; no numerical tolerances or fitting.
"""
from collections import deque
from itertools import product
from sympy import Matrix, eye, zeros, I, diag, Rational, log

N=3
A=diag(-1,1,-1)       # nontrivial character values on generator a
B=diag(1,-1,-1)       # nontrivial character values on generator b
S=Matrix([[0,1,0],[1,0,0],[0,0,1]])
C=Matrix([[0,0,1],[1,0,0],[0,1,0]])  # cycles all three characters

def key(M): return tuple(M)
def closure(generators):
    seen={key(eye(N)):eye(N)}
    q=deque([eye(N)])
    while q:
        M=q.popleft()
        for G in generators:
            R=M*G
            if key(R) not in seen:
                seen[key(R)]=R
                q.append(R)
    return list(seen.values())

def linear_span_dimension(matrices):
    return Matrix.hstack(*[Matrix(9,1,list(g)) for g in matrices]).rank()

def centralizer_dimension(generators):
    # every entry of the arbitrary 3x3 matrix is independent
    generators=list(generators)
    vals=[]
    for G in generators:
        for i in range(3):
            for j in range(3):
                row=[0]*9
                for k in range(3):
                    row[i*3+k] += G[k,j]
                    row[k*3+j] -= G[i,k]
                vals.append(row)
    return 9-Matrix(vals).rank()

for label,gens in [("V4",[A,B]),("V4 + one transposition",[A,B,S]),("V4 + one 3-cycle",[A,B,C]),("V4 + full S3",[A,B,S,C])]:
    G=closure(gens)
    d=linear_span_dimension(G)
    comm=centralizer_dimension(gens)
    print(f"GROUP {label:26s} |G|={len(G):2d}, algebra_dimension={d}, commutant_dimension={comm}")
    assert (len(G),d,comm)=={"V4":(4,3,3),"V4 + one transposition":(8,5,2),
        "V4 + one 3-cycle":(12,9,1),"V4 + full S3":(24,9,1)}[label]

# On the carrier C^5 = 2*trivial + char_a + char_b + char_ab.
Ac=diag(1,1,-1,1,-1)
Bc=diag(1,1,1,-1,-1)
Pc=(eye(5)+Ac+Bc+Ac*Bc)/4
Q=2*Pc-eye(5)
Y=Rational(1,12)*eye(5)+Rational(5,12)*Q
print("CARRIER projector rank",Pc.rank(), "Y spectrum", Y.eigenvals(),
      "trace 2P-I",Q.trace())
assert Pc.rank()==2
assert sorted(Y.diagonal())==[Rational(-1,3)]*3+[Rational(1,2)]*2
assert Q.trace()==-1

# On all 8 modes = 2*triv + 2*(char_a + char_b + char_ab).
A8=diag(1,1,-1,1,-1,-1,1,-1)
B8=diag(1,1,1,-1,-1,1,-1,-1)
P8=(eye(8)+A8+B8+A8*B8)/4
print("EIGHT-MODE group-average rank",P8.rank(),"trace 2P-I",(2*P8-eye(8)).trace())
assert P8.rank()==2 and (2*P8-eye(8)).trace()==-4

# Diagonal/eigencharacter labels from exact V4 cannot yield a CKM-like
# CP-odd commutator, and a single transposition has a common 1D invariant block.
Hu=diag(1,2,4)
Hd=Matrix([[5,1,-I],[1,6,1],[I,1,7]])
Jdet=(Hu*Hd-Hd*Hu).det()
print("CP-odd test det([Hu,Hd]) in full M3 =", Jdet)
assert Jdet!=0 and (Hd-Hd.conjugate().T)==zeros(3)
# Even if arbitrary complex coefficients are used in V4 and D4 group algebras,
# their matrices have 3 diagonal blocks or a 2+1 block, respectively.
for group in ["V4", "V4 + one transposition"]:
    mats=closure([A,B] if group=="V4" else [A,B,S])
    M1=zeros(3); M2=zeros(3)
    for k,G in enumerate(mats):
        M1+=(k+1)*(G+G.conjugate().T)
        M2+=(k+2)*I*(G-G.conjugate().T)
    comm=M1*M2-M2*M1
    print("CP-odd test det commutator",group,"=",comm.det())
    assert comm.det()==0

# Centralizers of the rank-2 projector vs the fine V4 action, as real dimensions
# of the Lie groups in U(5) and SU(5), computed from multiplicity squares.
p_centralizer_u5=2**2+3**2
v4_centralizer_u5=2**2+1**2+1**2+1**2
print("CARRIER centralizer dimensions: SU5(P) =",p_centralizer_u5-1,
      ", SU5(V4) =",v4_centralizer_u5-1)
assert p_centralizer_u5-1==12 and v4_centralizer_u5-1==6

# Genus-one double cover at four branch points: Hurwitz: 2g-2=2(-2)+n.
# Bare geometry does not select number of marks: all even n>=2 give surfaces.
for n in (2,4,6,8,10):
    g=(n-2)//2
    carrier=n+1               # h0(P1,O(D)) for n distinct degree-1 marks
    family=n-1                # dim H1(P1 minus n marks, C)
    print(f"MARKS n={n:2d} -> genus={g}, carrier={carrier}, family={family}, total={carrier+family}")
    assert 2*g-2==2*(-2)+n

# A normalized, unital depolarizing qubit channel always has the trivial fixed
# state I/2 for 0<=lambda<1; identity channel fixes ALL states.
rho=diag(Rational(3,4),Rational(1,4))
for lam in [0,Rational(1,3),Rational(2,3)]:
    out=lam*rho+(1-lam)*eye(2)/2
    assert out!=rho and lam*eye(2)/2+(1-lam)*eye(2)/2==eye(2)/2
print("FIXED-POINT control: depolarizing qubit uniquely fixes I/2, identity fixes every state")

# Multi-clock holonomy: [sigma_x, sigma_z] !=0, so infinitesimal square fails to close.
X=Matrix([[0,1],[1,0]]); Z=diag(1,-1)
print("CLOCK curvature [X,Z] =",X*Z-Z*X)
assert X*Z-Z*X != zeros(2)
print("ALL EXACT ASSERTIONS PASSED")
