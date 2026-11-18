# Original E8 current response and the native field dictionary

2026-10-16. Contract UR.NATIVE.JOINT_OPERATOR.01. Scoped mathematical
continuation; original ledger and physical gates remain unchanged/PARTIAL.
The stipulated positive E8 level-one lattice vacuum is a premise, not a
new derivation of a physical P1 source. No Standard Model module, Yukawa
coupling, detector instrument or physical time rule is inserted.

## 1. The actual mixed currents provide a nonscalar response

Use the original `_actual_current_source` Chevalley currents and compact
adjoint, not an abstract matrix algebra chosen to fit the answer. Its
twenty currents X_ia carry bar(5) tensor 4; F_ia=X_ia^dagger carry
5 tensor bar(4). The four original charged currents C_a carry 1 tensor 4.
These are weight-one currents, not established physical fermion fields.

Set B_ia,jb=[X_ia,F_jb], C_ij=(sum_a B_ia,ja)/4,
T=sum_i C_ii, S_ij=C_ij-delta_ij T/5 and Z=-4T/9. Actual brackets give
T|F=-9/4, T|C=-15/4 and hence Z|F=1, Z|C=5/3.
For a five-by-five matrix H define the actual source current

    Q(H)=-sum_ij (H-tr(H)I/5)_ij S_ji + tr(H)Z/5.

Then

    [Q(H),F_ka]=sum_j H_jk F_ja,
    [Q(H),C_a]=tr(H) C_a/3,
    kappa(Q(H),Q(K))=tr(HK)-tr(H)tr(K)/9.

`check_current_response.py` checks all 500 matrix-unit/field actions and
625 current-form matrix-unit pairs exactly in the original source.
Noncommuting H on the unresolved E+ plane therefore have an actual
noncommuting source action. This reuses the existing E8 Ward algebra;
the new result is its comparison with the actual native code/clock map.

For the existing polar direction c and projector P, put

    h_R=cc* + P tau/2.

Its eigenvalues are 1, +1/2, -1/2, 0, 0; exp(i pi h_R)=R5.
The source zero-mode exponential exp(i pi Q(h_R)_0) is an inner
automorphism, preserving the stipulated vacuum and conformal vector.
The source form is kappa(Q(h_R),Q(h_R))=25/18, identically for the WHOLE
previous CP1 marking fibre. Neither this norm nor existence of the
exponential selects the original polar point, pulse strength or time.

There is a concrete separated-point distinction. For Y=(5P-2I)/6,
the actual Ward responses <Q(Y)(-1) F(0) X(1)> have kinematic factor -1/2.
Using the original normalized source ray q and even-plane vector e1:

| Mark | q,q response | q,e1 response |
|---|---:|---:|
| Original polar mark | 0 | -sqrt(6)/12 |
| Previously defined countermark c=e0 | 1/6 | 0 |

Thus the source can distinguish these specified marks. A distinguishable
amplitude does not supply a law choosing one mark or a physical detector.

## 2. A genuine full-source order-four lift exists

Requiring R5 on F leaves three finite inner source automorphisms in this
gl5 family. Since det R5=-1, logarithms H of R5 have
tr H=pi(2k+1). Their action on C is scalar c=exp(i pi(2k+1)/3):

| tr(H)/pi modulo 6 | c | Order on full source |
|---|---|---:|
| 1 | exp(i pi/3) | 12 |
| 3 | -1 | 4 |
| 5 | exp(-i pi/3) | 12 |

The C and F currents with adjoints generate the stipulated E8 source,
so their joint action fixes the full automorphism. The order-four branch
is obtained, for example, using h_R+2cc*. Different logarithms with the
same finite endpoint do not select a unique continuous generator.
It would be incorrect to claim an unavoidable period-twelve obstruction.

## 3. What matches at the finite native group, and what does not

The checker compares ALL 60 original event charts with the actual
quartic code action using Gaussian-integer compression. The comparison
uses the original oriented determinant-path lift h=lambda r,
lambda=exp(-i pi/4). Its source five-action is

    A5=lambda^4 rho5(r)=-rho5(r).

The opposite currents F then transform as

    A5 tensor bar(h)=lambda^3 rho5(r) tensor bar(r).

This is exactly the actual W20 intertwiner action of h^tensor3. Its
central square is +i. It is a statement about the coupled finite cover
<G31,zeta8 I> intersect SU4, not the unchanged order-two word r in G31.
All phases and both representations must be transported together.

This finite match is not a current-algebra tensor dictionary. If the
four charged source currents C_a are identified with the unchanged native
four-register, every source clock branch has C^tensor3=-I20. On F20 its
action is R5 tensor I4, which is nonscalar. With the actual isometry W,

    || -W-W(R5 tensor I4) ||_HS^2=48.

Already infinitesimally, every traceless Q(H) acts trivially on C_a and
ALL their products/OPEs, but as H tensor I4 on F. Thus this literal
charged-current dictionary fails. This does not identify the original
native Cartan four-space with C4; that different realization is tested next.

## 4. Test the actual J-positive Cartan quartic realization

Use original Construction-A roots r of raw norm squared four and the
physical metric dot/2. In orthonormal Cartan coordinates alpha=r/sqrt(2),
so alpha^2=2. The ORIGINAL positive polarization is

    u_j=(h_2j - i h_(2j+1))/sqrt(2),  J u_j=i u_j.

The plus-sign convention would be the conjugate J=-i space. Keep the
minus sign when comparing the actual native charts.

The symmetric-tensor/Fock isometry sends the five actual code support
sums to homogeneous polynomials f_j(x), up to these NONZERO factors:

    sqrt(6)/24 [sum_i x_i^4],
    sqrt(2)/4 [x0^2 x1^2+x2^2 x3^2],
    sqrt(2)/4 [x0^2 x2^2+x1^2 x3^2],
    sqrt(2)/4 [x0^2 x3^2+x1^2 x2^2],
    x0 x1 x2 x3.

Their Fock Gram is exactly I5. This is the canonical oscillator embedding
of the existing code in the GIVEN source, not a derived physical
preparation or a choice of five charged root currents.

For b_j=(u_j,alpha)=(r_2j-i r_(2j+1))/2 the original lattice vertex law is

    Y(e^alpha,z) f(x) Omega
      = exp(sum_(n>0) alpha(-n) z^n/n) f(x-b/z) e^alpha.

A fixed nonzero root-field phase multiplies both sides without affecting
the argument. In the z^-1 coefficient, the PBW oscillator-partition (3)
component is EXACTLY

    (f(b)/3) alpha(-3) e^alpha.

Only the completely shifted degree-four term and the z^3 creation term
contribute to this partition. Other partitions cannot cancel it. This is
a PBW coefficient, not a coefficient in a normalized projection basis.

The five displayed quartics have no common nonzero complex zero. The last
forces a zero coordinate. The three pair quartics then force at most one
remaining nonzero coordinate, and the sum of fourth powers excludes that
last possibility. For every nonzero real root b is nonzero. Consequently
each root zero mode leaks at least one code vector into charge alpha.
Different root charges are orthogonal, so such root contributions cannot
cancel across a linear combination. Cartan zero modes annihilate the
entire neutral oscillator code. Thus the complex leakage map has rank240,
and its kernel is exactly the eight-dimensional Cartan algebra, acting
trivially on the code. `check_cartan_source_stabilizer.py` verifies the
actual 240 roots, all shift coefficients, Fock Gram and one-particle
number matrix sum_j <f_j|a_i^* a_k|f_j>=5 delta_ik exactly.

This gives a GLOBAL result. Let H be the compact inner E8 subgroup
preserving the code as a subspace of grade four. It is a closed compact
stabilizer. The leakage theorem gives Lie(H)=t8 and H^0=T8. Every element
normalizes this identity component, so its Cartan action is Weyl. The
one-particle support of the code is exactly the original H+ four-plane;
preservation of the code therefore preserves H+. Its restriction is a
code-preserving U4 operator. The already proved complete native code
stabilizer theorem forces its logical action into phase*rho5(G31).
The existing marked R5 lies outside that projective image. Therefore:

    No compact inner E8 automorphism preserving this native Cartan
    quartic code implements the marked R5 on it.

This includes arbitrary sequences of inner current-zero-mode exponentials
that return to the code at their endpoint, even if intermediate states
leave it. It is not merely an infinitesimal obstruction.

The original conformal generator does not repair this: all five code
states have exactly four mode-one oscillators and zero lattice charge,
so L0|Code=4I5. Consequently exp(-itL0) is scalar on the code, and its
conformal quarter exp(-i pi L0/2) is I5. Since [L0,Q(H)_0]=0,
adding conformal propagation to inner zero-mode pulses does not enlarge
their projective logical action. This follows from the existing conformal
grading, not from identifying L0 with physical P1 time. The current fields
themselves all have conformal weight one; their internal order-four
automorphism must likewise not be called conformal time evolution.
A cylinder Casimir shift L0-c/24 changes only the common phase.

## 5. Scope and decisive next question

Positive: the original stipulated current source has real nonscalar
responses and a full order-four inner clock; the actual coupled finite
event action matches the native twenty-action. Negative: neither the
literal charged-current tensor dictionary nor the canonical Cartan-code
inner-automorphism dictionary extends this into the required common clock.

This does NOT exclude nonlinear products of current modes, higher modes,
non-automorphism interventions, collective paths generated by such
operators, other derived composite field maps, or all P1 dynamics.
Multiplying inner-group exponentials is already excluded; introducing a
Hamiltonian polynomial in current modes is a genuinely different class.
Merely writing the desired logical matrix in the source enveloping algebra
would prove mathematical availability, not physical permission or polar
selection. The next useful test must obtain the first such operation and
its coefficients from an ORIGINAL source process, then check its action
on this actual code and its erasure carrier. Do not repeat a norm test,
an arbitrary gl5 pulse, or a dimension-based 4/5/20 field identification.
No unconditional 3+2, physical Higgs action, physical clock rate, alpha
coefficients or absolute Flavor normalization follows from this result.

## Sources and independent review

- Original TFPT modules: current_block_geometry.py, cartan_source.py,
  source_program.py, process.py; original compiler-current-product proof.
- Dong and Nagatomo, *Automorphism Groups and Twisted Modules for Lattice
  Vertex Operator Algebras*, Section2.2, p5 (vertex creation/annihilation
  formula): https://arxiv.org/pdf/math/9808088#page=5 .
- Dong and Griess, *Automorphism groups and derivation algebras of finitely
  generated vertex operator algebras*, Section2 (zero-mode derivations and
  exponentials): https://arxiv.org/pdf/math/0106051#page=2 .
- Native complete product/code stabilizer theorem and exact nonlift witness:
  MARKED-CLOCK-LIFT.md and transversal_audit.md in this same contract.
- Independent analytic review accepted source bracket signs, all three
  full-source lift orders, finite-cover phases, Ward sign, residual48,
  root leakage and global normalizer argument. It required explicit
  original J-positive minus-sign and physical root normalization, both
  used in the actual checker. No physical gate was independently verified.
