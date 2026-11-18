# Full Dirac lift of the relative endpoint family

Research continuation of `UR.COLLAR.FACE_CURRENT.01`, 21 October 2026.
Verdict: **conditional structural obstruction to an automatic P1 identification**.
The scalar current result in PROOF.md is retained. It must not be promoted to
the current of the entire selfadjoint Dirac module without the following test.
No general mathematical novelty or completed TFPT solution is claimed.

## 1. The original module, not an added mirror model

Use the already present product Dirac realization, from
`seam-full-collar-reconstruction-20261015/PROOF.txt:21-61`:

```
A = -i partial_theta + alpha,  alpha = 1/4 or 3/4,
B = diag(A,-A),   Jn = [[0,-I],[I,0]],   Gamma = diag(I,-I),
D = Jn (partial_r+B) = [[0,-partial_r+A],[partial_r+A,0]].
```

These are two graded components of one operator, with the original periodic
frame and the opposite-quarter relabelling retained. They are not two boundary
faces, two independent vacua, or a new set of four-dimensional particles.
The scalar transport `L=partial_r+A` is exactly its chiral block `D_plus`.
But `D_plus` maps the positive graded module to the negative one; it is not
the entire selfadjoint D on one fixed module.

For `T=exp(-ell B)`, keep BOTH physical boundary faces and both components:

```
G = [[a,c],[c,1-a]],
a=(1+exp(-2ell B))^-1, c=(2cosh(ell B))^-1.
```

On the second face the occupation is
`n(B)=diag(n(A),1-n(A))`, where `n(A)=(1+exp(2ell A))^-1`.
This is the same original graph constructor, not a fitted covariance.
Identifying it with a pure charged-CAR state remains an explicit statistical
premise, as in PROOF.md. This identification is not derived from P1.

## 2. The decisive cancellation is analytic, not a cutoff effect

For a complete orthogonal Dirac Cauchy projection C assume its actual normal
boundary unitary N satisfies `N C N*=1-C`. For the above finite collar it is
`N=diag(-Jn,Jn)`, and the identity follows from the actual graph formula.
It is the Lagrangian Cauchy property, not the geometric deck reflection and
not a definition of the original admissibility projector.

For bounded implementable generators X,Y, with `[C,X]` and `[C,Y]` Hilbert
Schmidt, use the charged-CAR cocycle

`c_C(X,Y)=Tr(C X (1-C) Y C - C Y (1-C) X C)`.

All displayed products are trace class. Cyclicity of products of Hilbert
Schmidt off-diagonal blocks gives `c_(1-C)(X,Y)=-c_C(X,Y)`.
If X and Y commute with N, unitary covariance also gives
`c_(1-C)(X,Y)=c_C(X,Y)`. Consequently **`c_C(X,Y)=0`**.
This is a sufficient condition; it does not classify all zero cocycles.

Smooth scalar multiplication on one boundary face, acting on the COMPLETE
Dirac fiber, commutes with N. On the circle its off-diagonal commutators are
Hilbert Schmidt. Therefore its full-module cocycle vanishes even on ONE face.
This cancellation is different from the cancellation BETWEEN the two faces
of a SINGLE scalar block established in PROOF.md.

In the product example, each component can be checked separately:

```
c_plus(f,g) = sum_(k>0) k (conj(f_k) g_k - f_k conj(g_k)),
c_minus(f,g) = -c_plus(f,g),  c_vector(f,g)=0.
```

The identity uses `n(-A)=1-n(A)` and the infinite telescoping argument already
proved for the scalar source. At `f=sin(4theta), g=1-cos(4theta)` these give
`-2i,+2i,0`. The commuting-profile Weyl rectangle has phase `exp(2ist)` on
the first chiral component, `exp(-2ist)` on the other, and 1 on the full
vector module. This does NOT exclude nonzero ordinary correlations, boundary
responses, mixed vector/axial currents, or global holonomies of other paths.
No finite matrix is asserted to represent the infinite CCR.

This lemma is restricted to the stated charged-CAR, Lagrangian, implementable
class. In higher boundary dimensions smooth multipliers need not have Hilbert
Schmidt commutators; the unrenormalized trace formula cannot simply be reused.
It is not a four-dimensional anomaly, mirror-decoupling, or no-physics theorem.

## 3. Which geometric unitary is actually allowed by the fixed symbols?

Retaining the original normal and tangential Clifford symbols means that a
pointwise unitary u commutes with both `Jn` and `Jn Gamma`. Their common
commutant on this two-component fiber is exactly `u=g I2`. Additional internal
multiplicity spaces would have their own commutant and are not classified here.

For `g_chi=exp(i chi(r) phi(theta))` the full lift is therefore
`U_chi=g_chi I2`, and

```
D_phi = U_chi D U_chi*
      = D - i chi_prime phi Jn - chi phi_prime Jn Gamma.
```

Both added multiplication terms together are bounded and Hermitian for smooth
real phi and chi on the compact collar. The original Clifford symbols and
grading are retained. Its Cauchy transport is
`diag(g,g) exp(-ell B)` at the second face, not a transport acting on only one
graded component. Both second-face components rotate; Section 2 applies.

The component-selective unitary `diag(g,I)` instead changes Jn and the
tangential Clifford symbol for every g != 1. It can be considered with their
full transport, or in a chiral operator/measure construction, but cannot be
called a scalar gauge perturbation of the SAME fixed full Dirac symbols.
Real structure J and the original charged representation still need their
actual bundle lift; this two-component computation does not supply them.

## 4. Fixed boundary conditions and transported ones are different sources

Let Bsa denote an original selfadjoint elliptic boundary trace subspace.
With domain `U_chi Dom(D_Bsa)`, D_phi is unitarily equivalent to D_Bsa.
If the reflection, reference, J and admissibility complex are transported
covariantly as well, the relative spectral traces and eta invariants are
unchanged. A nonzero ordinary coordinate derivative of G alone does not
establish a physical source variation. The existing covariant-source-jet
criterion removes precisely this frame motion.

Keeping Bsa FIXED is nevertheless a legitimate analytic alternative in the
displayed compact example: D_phi-D is bounded symmetric, so D_phi on the
original selfadjoint domain is selfadjoint. Conjugation by U_chi no longer
identifies the same-domain operators unless it preserves Bsa. This by itself
does not prove different spectra; a fixed-domain spectral/action response
has not been computed here. After conjugating back by U_chi*, its
boundary condition is `U_boundary* Bsa`. Thus any such change is a relative
boundary attachment, and must be specified as that physical context; one
cannot discard its variation in the reference, projector or action.

For the actual reference split, keep the geometric involution R distinct
from the Calderon polarization. At fixed R, writing `R_U=U* R U`, exactly

```
D_ref(D_phi,R)=U (D+R_U D R_U)/2 U*,
D_rel(D_phi,R)=U (D-R_U D R_U)/2 U*.
```

In particular, if U commutes with the fixed R the reference and relative
DIFFERENTIAL expressions are still conjugate; holding R fixed alone does
not make the gauge family a physical action variation. Domains remain part
of this statement. If `[U,R]!=0`, with `U_t=exp(itF)` and `K=[F,R]`,

```
delta D_ref - i[F,D_ref] = -i (K D R + R D K)/2,
delta D_rel - i[F,D_rel] = +i (K D R + R D K)/2.
```

These are the extra reference terms; deleting them is not a fixed-R
calculation. For a covariantly transported full admissibility complex its
Hodge projector instead obeys `P_t=U P_0 U*`. None of these identities
determines the original P1 physical F or the OS operator implementing it.

For the gapped APS sign subspace on a circle, a scalar multiplier and its
adjoint preserve the fixed subspace iff g is constant: preserving both Hardy
halves forces every nonzero Fourier coefficient of g to vanish. Infinitesimally,
`sin(4theta)` sends the negative mode -1 to positive mode 3 with coefficient
`1/(2i)`, and `1-cos(4theta)` does so with coefficient `-1/2`.
Fixing the four marks therefore does not preserve the original APS subspace.
This is NOT a failure of selfadjointness of the bounded perturbation on its
old domain. Nor does it rule out a physical boundary drive with that domain.

Two interior extensions chi with the same endpoints differ by a smooth unitary
equal to I on both boundary components. It preserves the fixed trace condition
and conjugates the two realizations. Therefore even the fixed-domain spectrum
depends on the endpoint attachment, not on an arbitrary interior interpolation.

There is a second, stronger trace-space distinction: the one-face endpoint
unitary `V=I direct_sum g` preserves the original scalar graph G only for g=I.
Indeed its off-diagonal block would require `c g*=c`, and c is injective.
This observation concerns the Cauchy trace projection. **It does not identify
that projection with Padm on the OS Hilbert space.**

The original physical net requires BOTH A and A* to preserve Ran(Padm),
equivalently `[A,Padm]=0` for bounded A on the full OS space. Neither a moving
Cauchy graph nor a sign-projector calculation supplies this missing physical
operator map. Compression by Padm is explicitly not the definition of the net.

## 5. The chiral route survives, but its measure is an actual missing datum

The original Gamma already defines the chiral block D_plus. One need not
invent an oscillator or a new component label to formulate that route.
However replacing the full module by its chiral determinant/quantization is
not a change of notation: it selects which field measure and line connection
are physically used. That identification is the outstanding task.

The existing source integration result already prevents a shortcut here:
`finding:joint-fermionic-source-schur-functional` shows that a DIRECT
identification of a primitive trace profile with `-log det K` requires
`x f_prime(x)+1=0`. An arbitrary declared spectral profile does not provide
that fermionic measure. A separate fermionic action is possible, but its
actual kinetic operator, measure and source coupling must then be supplied
by the original construction. This is reused prior evidence, not a new test.

For the paired adapted B, Jn pairs every eigenvalue with its negative, so
`eta(B)=0`. For the scalar A in the original periodic quarter frame,
Hurwitz zeta continuation gives `eta(A)=1-2alpha`, namely +1/2 or -1/2.
The two quantities belong to DIFFERENT operators. Even the scalar eta is
constant along `g A g*`; its value is not by itself the curvature of the
chiral implementer line. There is no license to substitute one for the other.

Original source contracts explicitly require a chiral measure and line WITH
connection: `CHIRAL4D.NOMIRROR.01` C1-C6, and
`SEAM.DETLINE.UNIFICATION.01`. Their original current manuscript passages and
the live CSV ledger still state Open. The Atlas annotation of the former has
stale source references; its status was therefore checked directly in the
current ledger, not inferred from that annotation. No four-dimensional part
of those contracts is proved by this boundary calculation.
Their four-dimensional requirements are not imposed as prerequisites for
further one-dimensional boundary work. The immediate missing map is already
the boundary source/measure/connection identification in the original datum.

## 6. Consequence for the solution search

The automatic route "full P1 Calderon graph -> same scalar one-face current
phase" fails in the stated complete Dirac/CAR realization. Restoring the
omitted graded component cancels its central phase. The previous scalar
bridge is correct on its stated scalar/chiral source, and remains conditional.

The next indispensable original map is now more specific: the P1 real/graded
source and its physical admissible measure must select the chiral block (or
supply a compensating bulk line/inflow) and identify its connection with the
original charged CAR/E8 source, while retaining the actual boundary attachment,
reference and A/A*-admissibility. A declaration of a new projector or fitted
coupling would not establish that map.

Native Hamming/E8/J, quartic4/5, erasure20, primary20 carry, all eight source
preimages, original time and winding are unchanged. No finite code-closure
condition is added. No T1-T8 goal, event-selection rule or full solution closes.

## Sources and verification

- Original full product Dirac: `seam-full-collar-reconstruction-20261015/PROOF.txt:21-61`.
- Original P1 datum, primitive action and physical net:
  `documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex:65-80,302-343`.
- Original full selector: `documents/archive--tfpt-45--source_extracts--04_qft_source/sections/002-section-closure-theorems-from-the-admissibility-complex.tex:12-123`.
- Existing frame subtraction: `seam-collar-multifield-20261016/P1-COVARIANT-INSERTION.txt:26-146`.
- Existing joint source-measure identification test: `seam-collar-multifield-20261016/PROOF.txt:271-326`.
- Original open chiral and line contracts: `documents/tfpt_research_contracts/sections/015-source.tex:235-266`, `017-source.tex:239-285`.
- [Baer and Ballmann, arXiv:1307.3021v2](https://arxiv.org/abs/1307.3021), Lemma 3.2/Remark 3.3, Theorems 4.11-4.12: adapted operator, paired spectrum, boundary domains.
- [Langmann, JMP35 (1994), arXiv:hep-th/9305014](https://arxiv.org/abs/hep-th/9305014), equations 31-34: implementable CAR currents and trace-class cocycle.

`check_dirac_lift.py` reuses the existing graph constructor and checks the
normal-complement identity, full/component cocycles, symbol commutant and
differential lift exactly. Finite matrices are controls of those analytic
identities, not a numerical proof of an infinite field theory.
