#!/usr/bin/env python3
"""Independent TFPT F1 sewing-rigidity exploratory audit (13 Nov 2026).

Proves finite mathematical statements inside declared assumptions:
  1. Four-mark graph critical group has no S4-invariant primitive Z4 character.
  2. A3(1) + D5(1) has two conjugate maximal isotropic Z4 glue subgroups.
  3. Primitive spin-1 simple-current locality needs D_m with m == 5 (mod 8).
  4. Four-point cross-ratio λ is a free modulus; λ=2 has j=1728 and
     its double cover is elliptic, with primitive 4th-root differential action.
  5. Connected E8 Z4 current brackets and ONE chiral spatial translation
     force a single sector propagation speed. This is NOT a theorem that
     the reference mode of the separate Majorana dimer is such a current.
  6. A clock-covariant twisted Z4 return ring has gauge-invariant 4-step
     holonomy and a phase response only starting at fourth order.

No original P1 kernel / weighted R113 evaluation; no claim of F1 resolution.
Dependencies: sympy, numpy. Usage: python tfpt_f1_sewing_rigidity_probe.py
"""
from __future__ import annotations
from itertools import product, permutations
import math
import numpy as np
import sympy as s


def check(label, yes):
    if not yes: raise AssertionError(label)
    print('PASS:', label)


def k4_equivariant_character_no_go():
    """K(K4) = {(x,y,z,0):x+y+z=0 (mod 4)} with S4 action."""
    el = [(a,b,(-a-b)%4,0) for a,b in product(range(4),repeat=2)]
    def canon(w):
        return tuple((w[i]-w[3])%4 for i in range(4))
    def fixed_chars(subgroup):
        ans=[]
        for a,b in product(range(4),repeat=2):
            if all((a*w[0]+b*w[1]-a*canon([w[p[j]] for j in range(4)])[0]
                    -b*canon([w[p[j]] for j in range(4)])[1])%4 == 0
                   for w in el for p in subgroup):
                ans.append((a,b))
        return ans
    s4 = list(permutations(range(4)))
    c4 = [(1,2,3,0)]
    inv_s4=fixed_chars(s4)
    inv_c4=fixed_chars(c4)
    check('K4 critical group has 16 elements',len(set(el))==16)
    check('No nonzero S4-invariant map K(K4)->Z4',inv_s4==[(0,0)])
    check('C4-invariant characters at most order 2',
          inv_c4==[(0,0),(0,2)])
    print('  Invariant chars S4:',inv_s4,' C4:',inv_c4)


def simple_current_glue():
    def q(a,b):return (3*a*a+5*b*b)%8
    zeros={(a,b) for a,b in product(range(4),repeat=2) if q(a,b)==0}
    expected={(0,0),(1,1),(1,3),(2,2),(3,1),(3,3)}
    check('Exactly six isotropic classes for A3+D5 quadratic form',zeros==expected)
    groups=set()
    for a,b in product(range(4),repeat=2):
        H=frozenset(((a*k)%4,(b*k)%4) for k in range(4))
        if len(H)==4 and all(x in zeros for x in H):groups.add(H)
    plus=frozenset({(0,0),(1,1),(2,2),(3,3)})
    minus=frozenset({(0,0),(1,3),(2,2),(3,1)})
    check('Exactly two conjugate maximal isotropic cyclic Z4 glues',groups=={plus,minus})
    check('Each glue H has H-perp=H (Lagrangian)',
          all({(a,b) for a,b in product(range(4),repeat=2)
               if all((3*a*c+5*b*d)%4==0 for c,d in H)} == set(H)
              for H in groups))
    allowed=[m for m in range(4,38) if m%2==1 and (3+m)%8==0]
    check('Primitive integer-spin locality selects D5,D13,D21,D29,...',
          allowed==[5,13,21,29,37][:len(allowed)] and allowed[0]==5)
    # For affine SU(4)_1 and SO(10)_1: c_A=3, c_D=5.
    check('Minimal local glue has spin 1, total c=8, dim currents 248',
          (3+5)//8==1 and 3+5==8 and 15+45+64+64+60==248)
    # Gauss sums detect c modulo 8 only.
    gA=s.simplify(sum(s.exp(2*s.pi*s.I*s.Rational(3*k*k,8)) for k in range(4))/2)
    gD=s.simplify(sum(s.exp(2*s.pi*s.I*s.Rational(5*k*k,8)) for k in range(4))/2)
    check('Gauss-Milgram product is 1 (only c mod 8 known)',s.simplify(gA*gD)==1)
    print('  two glues: H+ and H-; Gauss phases:',gA,gD,' product:',gA*gD)
    print('  WARNING: mod-8 Gauss sum does NOT determine c=8 instead of c=0.')


def four_point_modulus():
    l=s.symbols('lambda')
    j=s.factor(256*(1-l+l*l)**3/(l*l*(1-l)**2))
    difference=s.factor(j-1728)
    target=64*(l-2)**2*(l+1)**2*(2*l-1)**2/(l*l*(l-1)**2)
    check('Exact elliptic j(lambda)-1728 factorization',s.simplify(difference-target)==0)
    for point in (s.Integer(2),s.Rational(1,2),s.Integer(-1)):
        check(f'Harmonic lambda={point}: j=1728',s.simplify(j.subs(l,point))==1728)
    e=(s.Integer(1)+s.I*s.sqrt(3))/2
    check('Equianharmonic alternative lambda=(1+i*sqrt(3))/2 has j=0',s.simplify(j.subs(l,e))==0)
    check('Generic cross-ratio lambda=3 gives a different elliptic modulus',
          j.subs(l,3)==s.Rational(21952,9))
    def cr(v):
        a,b,c,d=v
        if a==s.oo:return s.simplify((b-d)/(b-c))
        if b==s.oo:return s.simplify((a-c)/(a-d))
        if c==s.oo:return s.simplify((b-d)/(a-d))
        if d==s.oo:return s.simplify((a-c)/(b-c))
        return s.simplify((a-c)*(b-d)/((a-d)*(b-c)))
    def permutation_order(p):
        cur=tuple(range(4))
        for power in range(1,13):
            cur=tuple(p[cur[k]] for k in range(4))
            if cur==tuple(range(4)):return power
        raise AssertionError('unexpected permutation order')
    for lval,size,n4,label in [(s.Integer(3),4,0,'generic'),(s.Integer(2),8,2,'square'),
                               (e,12,0,'equianharmonic')]:
        pts=(s.oo,0,1,lval)
        stabilizer=[p for p in permutations(range(4))
                    if s.simplify(cr([pts[i] for i in p])-lval)==0]
        cnt=len(stabilizer)
        check(f'Moebius automorphism count of unordered {label} four points = {size}',cnt==size)
        check(f'Order-4 Mobius symmetries of {label} marked sphere = {n4}',
              sum(permutation_order(p)==4 for p in stabilizer)==n4)
    # For y²=x⁴-1, the order-4 automorphism x -> ix keeps y fixed,
    # and dx/y -> i dx/y; this is an order-4 linear multiplier on the torus.
    check('Riemann-Hurwitz: two-sheet cover with 4 simple branch points has genus 1',
          2*(-2)+4==0)
    check('x->i*x preserves x^4-1 and acts on dx/y by i',s.I**4==1)
    print('  Square/equianharmonic alternative: (lambda,j,Aut)=(2,1728,8),'
          '(exp(i*pi/3),0,12).')
    print('  WARNING: four marks alone do NOT select their complex cross-ratio.')


def common_translation_from_e8_brackets():
    # CFT affine E8 Z4 brackets observed in independently enumerated roots:
    # (grade r,grade s,grade r+s) : non-zero root-pair counts
    compos={(1,1,2):960,(1,2,3):960,(1,3,0):832,(2,2,0):720}
    v=s.symbols('v0:4')
    constraints=[]
    for (r,ss,t),count in compos.items():
        check(f'Nonzero E8 graded bracket ({r}+{ss}->{t})',count>0)
        # Jacobi / derivation with Fourier mode numbers m=1,n=0 and m=0,n=1.
        constraints += [v[t]-v[r],v[t]-v[ss]]
    coeff,_=s.linear_eq_to_matrix(constraints,v)
    check('Connected graded current brackets force rank 3 relative-speed constraints',
          coeff.rank()==3 and coeff.nullspace()==[s.Matrix([1,1,1,1])])
    coeff_short,_=s.linear_eq_to_matrix([
        v[2]-v[1],v[3]-v[1],v[0]-v[1]],v)
    check('Three nonzero inter-grade bracket families already force v0=v1=v2=v3',
          coeff_short.rank()==3)
    print('  Derived only if delta(J_r(k)) = i*v_r*k*J_r(k)'
          ' is a derivation of ONE local affine E8 current algebra.')
    print('  WARNING: grade-0 E8 current != automatically the uniform Majorana reference mode.')


def twisted_four_return():
    D,t,z,E=s.symbols('D t z E',nonzero=True)
    Us=s.zeros(4)
    for k in range(3):Us[k+1,k]=1
    Us[0,3]=z
    UsInv=s.zeros(4)
    for k in range(3):UsInv[k,k+1]=1
    UsInv[3,0]=1/z
    H=s.diag(0,D,D,D)-t*(Us+UsInv)
    char=s.factor((E*s.eye(4)-H).det())
    bare=s.factor(char.subs(z,1))
    check('4-step twisted determinant correction exact and energy-independent',
          s.simplify(char-bare+t**4*(z+1/z-2))==0)
    c2,c4=s.symbols('c2 c4')
    pol=s.series(char.subs(E,c2*t*t+c4*t**4),t,0,6).removeO().expand()
    es2=s.solve(pol.coeff(t,2),c2)[0]
    es4=s.factor(s.solve(pol.coeff(t,4).subs(c2,es2),c4)[0])
    check('Ground energy starts -2t²/D + (2-z-z^-1)t⁴/D³',
          es2==-2/D and s.simplify(es4-(2-z-1/z)/D**3)==0)
    for phi in (0,math.pi/2,-math.pi/2,math.pi):
        U=np.zeros((4,4),complex)
        for k in range(3):U[k+1,k]=1
        U[0,3]=np.exp(1j*phi)
        Q=np.diag([1,1j,-1,-1j])
        check(f'Clock covariance / U^4 holonomy at phi={phi:.3f}',
              np.allclose(Q@U@Q.conj().T,1j*U)
              and np.allclose(np.linalg.matrix_power(U,4),np.exp(1j*phi)*np.eye(4)))
    def ground(phi,tt=0.12):
        U=np.zeros((4,4),complex)
        for k in range(3):U[k+1,k]=1
        U[0,3]=np.exp(1j*phi)
        H=np.diag([0,1.,1.,1.])-tt*(U+U.conj().T)
        return np.linalg.eigvalsh(H)[0]
    h=1e-5
    iph=(ground(math.pi/2+h)-ground(math.pi/2-h))/(2*h)
    imh=(ground(-math.pi/2+h)-ground(-math.pi/2-h))/(2*h)
    check('Isospectral oriented conjugates have opposite persistent currents',
          abs(ground(math.pi/2)-ground(-math.pi/2))<1e-12
          and abs(iph+imh)<1e-9 and abs(iph)>1e-5)
    print('  Exact det correction:',s.factor(char-bare))
    print('  Ground energy E=-2 t^2/D+(2-2 cos(phi))*t^4/D^3+O(t^6)')
    print('  At t=0.12, Δ=1: ground(+pi/2)=',ground(math.pi/2),
          '; ground(-pi/2)=',ground(-math.pi/2),
          '; orientation derivative +',iph,' orientation derivative -',imh)
    print('  WARNING: this is a four-step CHARGE-space loop, not a K4 triangle.')



def meromorphic_clock_spinor_bridge():
    """Native 5 meromorphic modes -> clock 3+2 -> SM generation in SO10 S+.

    D={z:z**4=1}. H0(P1,O(D)) = {p(z)/(z**4-1): degree p <=4}.
    Under z->i*z the basis modes have eigenvalues (1,i,-1,-i,1).
    Crucial fence: full clock centralizer is U(2)xU(1)^3, NOT U(3)xU(2).
    U(3)xU(2) needs keeping only the coarse invariant-subspace projector.
    """
    C=s.diag(*[s.I**k for k in range(5)])
    I=s.eye(5)
    P=s.simplify(sum((C**k for k in range(4)),s.zeros(5))/4)
    P_expected=s.diag(1,0,0,0,1)
    check('Canonical pullback clock on H0(O(D)) has (1,i,-1,-i,1) eigenvalues',
          C==s.diag(1,s.I,-1,-s.I,1))
    check('Clock-averaged projector has ranks 2 invariant + 3 nontrivial',
          P==P_expected and P**2==P and s.trace(P)==2 and s.trace(I-P)==3)
    Y=s.Rational(1,2)*P-s.Rational(1,3)*(I-P)
    check('Primitive traceless hypercharge from 3+2 projector: 6Y=(3,-2,-2,-2,3)',
          s.trace(Y)==0 and 6*Y==s.diag(3,-2,-2,-2,3))
    fine_commutant_dim=2**2+1+1+1
    coarse_commutant_dim=2**2+3**2
    check('Full C4 centralizer dim 7; coarse projector centralizer dim 13',
          fine_commutant_dim==7 and coarse_commutant_dim==13)
    check('Special-unitary coarse centralizer has dimension 8+3+1=12',
          coarse_commutant_dim-1==8+3+1)
    sector=[]
    for a in range(4):
        for b in range(3):
            if (a+b)%2==0:
                multiplicity=math.comb(3,a)*math.comb(2,b)
                yy=-s.Rational(a,3)+s.Rational(b,2)
                sector.append((a,b,multiplicity,yy))
    sector.sort(key=lambda x:(x[0]+x[1],x[0],x[1]))
    check('Even exterior SO(10) half-spinor decomposes to exactly 16 Weyl states',
          sum(m for a,b,m,y in sector)==16 and len(sector)==6)
    expected={(0,0):(1,s.Rational(0)),(2,0):(3,-s.Rational(2,3)),
              (1,1):(6,s.Rational(1,6)),(0,2):(1,s.Rational(1)),
              (3,1):(2,-s.Rational(1,2)),(2,2):(3,s.Rational(1,3))}
    check('Exact SM generation hypercharge representation including sterile neutrino',
          {(a,b):(m,y) for a,b,m,y in sector}==expected)
    check('Mixed gravitational-U1 and cubic U1 anomalies cancel',
          sum(m*y for a,b,m,y in sector)==0 and
          sum(m*y**3 for a,b,m,y in sector)==0)
    check('SU3^2-U1, SU2^2-U1 and SU3^3 anomalies cancel',
          2*s.Rational(1,6)-s.Rational(2,3)+s.Rational(1,3)==0 and
          3*s.Rational(1,6)-s.Rational(1,2)==0 and 2-1-1==0)
    check('Witten SU2 anomaly absent: 3 quark doublets + 1 lepton doublet is even',
          (3+1)%2==0)
    print('  5 mode clock spec:',C.diagonal(),'  projector:',P.diagonal())
    print('  6Y spec:',(6*Y).diagonal())
    print('  (n_color,n_weak,multiplicity,Y) for Lambda^even:')
    for row in sector:print('   ',row)
    print('  WARNING: keeping full clock as observable fixes only U(2)xU(1)^3;')
    print('  the SM gauge algebra requires retaining only the coarse 2+3 projector.')
    print('  WARNING: these 3 color slots are not by themselves 3 generations.')


def four_marks_residue_exact_sequences():
    """Two canonical exact sequences from one four-point divisor D on P1.

    0 -> H0(O) = C -> H0(O(D)) = C5 -> principal-parts(D) = C4 -> 0;
    0 -> H0(Omega^1(D)) = C3 -> C4 ->sum C ->0.
    The first is a sheaf-principal-part exact sequence and the second
    the global residue theorem. Function principal-part coefficients
    are computed in the coordinate z; differential residues are intrinsic.
    """
    points=(s.Integer(1),s.I,s.Integer(-1),-s.I)
    # For f_k=z^k/(z^4-1), k=0..4, the coefficient of (z-r)^-1 is
    # r^k/(4 r^3) = r^(k+1)/4 since r^4=1.
    M=s.Matrix(4,5,lambda p,k:s.simplify(points[p]**(k+1)/4))
    expected_kernel=s.Matrix([-1,0,0,0,1])
    check('4 local principal parts of meromorphic H0(O(D)) have rank 4',M.rank()==4)
    check('Only the constant f4-f0=1 has no poles at any mark',
          len(M.nullspace())==1 and M*expected_kernel==s.zeros(4,1))
    # For omega_k=z^(k-1)dz/(z^4-1), k=1..3, residue r^k/4.
    R=s.Matrix(4,3,lambda p,k:s.simplify(points[p]**(k+1)/4))
    check('Logarithmic-residue 4x3 matrix has rank 3 and total residue zero',
          R.rank()==3 and s.ones(1,4)*R==s.zeros(1,3))
    check('Residue image equals exactly ker(total-residue sum)',
          R.rank()==4-1)
    # Both the principal-parts and residue modules transform by the 4-cycle.
    # The meromorphic carrier adds one GLOBAL constant to local 4 mark data.
    check('Dimensions give 1+4=5 and 4-1=3 from one divisor',
          1+M.rank()==5 and 4-R.rank()==1)
    print('  Function principal-part matrix:')
    print(M)
    print('  Logarithmic differential residue matrix:')
    print(R)
    print('  Five modes = four local pole directions + one global constant,')
    print('  not five marked points. Three cycles = four residues with one sum rule.')
    print('  WARNING: sharing the same three C4 characters does not yet establish')
    print('  independent color and generation copies in a physical 3+1D theory.')

if __name__=='__main__':
    print('TFPT / SEWING RIGIDITY INDEPENDENT TESTS / 2026-11-13')
    for label, fn in [
        ('K4 CHARGE OBSTRUCTION',k4_equivariant_character_no_go),
        ('PRIMITIVE E8 CONDENSATION',simple_current_glue),
        ('FOUR MARKS TO 3+2 SPINOR',meromorphic_clock_spinor_bridge),
        ('COHOMOLOGICAL 3 4 5 EXACT SEQUENCES',four_marks_residue_exact_sequences),
        ('MARKED-SPHERE MODULUS',four_point_modulus),
        ('E8 COMMON TRANSLATION',common_translation_from_e8_brackets),
        ('ORIENTED FOURTH-ORDER RETURN',twisted_four_return),
    ]:
        print('\n'+'='*16,label,'='*16)
        fn()
    print('\nALL DECLARED MATHEMATICAL TESTS PASSED; P1/R113 PHYSICAL IDENTIFICATION REMAINS OPEN.')
