# F1: all D4 clock weights excluded; six-history Gram source remains open

**16 November 2026. Exact finite mathematics, NOT the native TFPT physical source. F1 OPEN.**

## Universal no-go for a tracial/KMS-symmetric two-step quantum clock

For the six anchored 4-cycles split into geometric rotation class C2
(two elements) and remaining class C4 (four elements), define on M4(C)

\[
\Phi_x(A)=\frac{x}{2}\sum_{c\in C_2}P_cAP_c^\dagger+
          \frac{1-x}{4}\sum_{c\in C_4}P_cAP_c^\dagger,\quad0\le x\le1.
\]

The exact determinant is x^4(x-1)^8(3x-1)^2/1024.
There exist four mutually orthogonal, x-independent projectors
P0, Pa, Pb, Pc with ranks 2,4,8,2 for which

\[
\Phi_x^2=P_0+x^2P_a+\frac{(1-x)^2}{4}P_b+
                       \frac{(3x-1)^2}{4}P_c.
\]

The cases x=0,1/3,1 are singular and cannot be a finite matrix
exponential. At all other x, the only tracially symmetric candidate
generator is L_x=log(Phi_x^2). Two explicit Omega-perpendicular Choi
test vectors (row-major matrix units) are

\[
v=(0,1,0,-1,\ 1,0,-1,0,\ 0,-1,0,1,\ -1,0,1,0),
\]
\[
w=(0,0,1,0,\ 0,0,0,1,\ 1,0,0,0,\ 0,1,0,0).
\]

Their normalized tangent expectations are

\[
f_v(x)=-2\ln x+\frac43\ln(|3x-1|/2),
\quad
f_w(x)=2\ln x-4\ln((1-x)/2)+\ln(|3x-1|/2).
\]

For 1/4<x<1, x different from 1/3, f_v is strictly negative:
4x^3-(3x-1)^2=(1-x)^2(4x-1)>0.
For 0<x<=1/4, f_w is strictly negative:
(1-x)^4 >=81/256 >32/243 >=8x^2(1-3x),
since the last cubic is maximized at x=2/9.

These exact conditional Choi violations reject **every** x as
an exact finite-time member of a tracial/KMS-symmetric completely
positive quantum Markov semigroup on full M4. They do NOT reject
non-selfadjoint logarithms, nontracial KMS states, non-Markov
memoryful sources or dilations. The classical four-state two-step
positive semigroup is a distinct source type.

## Unified primitive to investigate: positive history Gram kernel

The six P_c matrices are linearly independent. For a Hermitian 6x6
history kernel K define

\[
{\cal E}_K(\rho)=\sum_{c,d} K_{cd}P_c\rho P_d^\dagger.
\]

With W=(vec(P1),...,vec(P6)), its Choi matrix equals W K W^dagger
and rank W=6, so CP is equivalent to K>=0 in this fixed representation.
TP requires sum(K_cd P_d^dagger P_c)=I4.
Diagonal K gives random-unitary history mixing; rank-one K=aa^dagger
gives coherent history amplitudes, with their full unitary condition.

The D4 action on six histories has seven ordered orbitals and six
unordered orbitals. The space of Hermitian D4-invariant K has real
dimension 7. TP imposes three independent real linear conditions,
leaving FOUR continuous degrees of freedom. The matrix K=I6/6 is
strictly positive and trace preserving; it has a nonempty open
neighborhood of physically valid CP/TP kernels in this four-dimensional
affine set. The previous scalar x family was a proper subfamily.

Of all 1024 relative four-root phase assignments, nine normalized
coherent sums are unitary, but only ONE rank-one Gram kernel is
invariant under **untwisted** D4. Its unitary is the negative half-turn,
U=-P_(r^2), U^2=I. The genuine period-four solutions need extra
native spin/projective data or symmetry selection. Their physical
source has not been found.

## Constructive alternative requiring new physical source data

The Poisson jump semigroup exp[t gamma(Phi_x-id)] is CPTP for any
gamma>=0, and tracially symmetric for inverse-paired D4 histories.
It is NOT equal to Phi_x^2 at unit time. Hence it evades the no-go
only by changing the history composition law. Neither gamma nor
the microscopic Poisson assumption follows from the source.

The fully quantum random-unitary source has a two-dimensional fixed
algebra, spanned by I4 and J4, unlike the unique stationary direction
of the four-state classical chain. A unique vacuum needs additional
source structure.

## Important spectral type check

The user-provided eigenvalues
{1:2,1/3:6,-1/3:6,-1:2}, squared {1:4,1/9:12},
also arise from T=(J4-I4)/3 tensor a four-dimensional involution
F with eigenvalues (+1,+1,-1,-1). The resulting map on M4 is
not even Hermiticity preserving. Spectra alone do not certify a
valid quantum channel or a physical time generator.

## New constructive source test: oriented seam jump, not a fitted D4 two-step clock

The unique S4-invariant scalar projector in the active four-mark
permutation representation is P0=J4/4; put Q=I-P0. A source with
three covariant decay jumps from the standard sector Q to P0 has

\[
{\cal L}_{seam}(\rho)=\gamma\bigl[
  {\rm Tr}(Q\rho)P0-\{Q,\rho\}/2\bigr].
\]

The normalized full M4 spectrum is {0:1,-gamma/2:6,-gamma:9}.
The projected Choi tangent has eigenvalues {1:3,0:13} in units
gamma, so this is an actual conservative, completely positive
GKLS generator, not a false-positive classical clock.
It has the unique **pure** stationary density P0 and exact
finite-time sewing, for q=exp(-gamma t),

\[
{\cal T}_q(\rho)=(P0+\sqrt q Q)\rho(P0+\sqrt q Q)
      +(1-q){\rm Tr}(Q\rho)P0.
\]

Its full S4 covariance follows from P0=J4/4. The commutant of
the 24 mark permutations restricted to Q is one-dimensional;
therefore a covariant positive jump-rate matrix on the three
standard directions is proportional to I3. Thus, **inside the
explicitly assumed pure one-way S4-covariant jump class**, just
one overall decay scale gamma remains. This is a conditional
selection theorem, NOT a physical derivation of the direction,
jump operators, original time scale, native projective phases,
Q0, P_adm or physical OS positivity.

Every original four-cycle P_c leaves both P0 and Q invariant.
No linear combination or Gram-K coupling of these permutation
operators alone can transfer population from Q into P0.
Consequently this oriented dissipative source is a genuinely
additional primitive source operation, not a new value of x
inside the excluded original random-unitary family.

An exact process-memory census is also informative. After one
four-cycle step there are six word operators with Choi/Kraus
span rank 6. All 36 two-step words produce the 12 even A4
permutations, whose images span rank 10 in M4; at three steps
the odd coset also spans rank 10. Thus the minimal independent
Kraus history dimensions for positive-probability consecutive
random steps are 6,10,10. This finite ceiling is not the
original unbounded native full96 **spatial** Fourier Hankel
rank 192(k+1) and does not imply physical time.

This model has a rank-one, **nonfaithful** stationary state;
it is not a tracial/KMS-symmetric faithful physical vacuum.
It makes no claim to be the original TFPT R113 quantum source.

Reproduce with:
- research/originlab/f1_oriented_seam_source.py
- research/originlab/test_f1_oriented_seam_source.py

## Native experiment still required

Derive K_cd from the SAME original all-field R113 source, full
outer P_adm, physical OS reflection, actual spin refinement and
absolute empty reference Z0, not a synthetic positive matrix.
Then test native two-level blocking, all four original responses
Z0/ZU/ZY/ZUY, full U7 kinetic energy, family, clock and stress
on ONE shared source-generated GNS/Dirichlet form. Only after those
source checks consider native E8, 3+1D chirality and spin-2.

The actual native diagnostic remains whole_source_bundle=null
and identifies the circle as spatial_tangential_circle_NOT_OS_time.
No physical F1 completion is claimed.

Exact implementation:
- research/originlab/f1_history_gram_universal_no_go.py
- research/originlab/test_f1_history_gram_universal_no_go.py
