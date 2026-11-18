# Exact stabilizer of the actual quartic code

2026-10-16. Independent mathematical continuation of
`UR.NATIVE.JOINT_OPERATOR.01`; no ledger or physical-goal promotion.

**Result.** For the actual `_code_basis` code isometry
`V : C^5 -> (C^4)^{tensor 4}`, let

\[
 K=\{U\in U(4):U^{\otimes4}\operatorname{im}V=\operatorname{im}V\}.
\]

Then, using the already established native identification with the two-qubit
Clifford representation and the finite-overgroup theorem below,

\[
 \boxed{K=U(1)G_{31}.}
\]

Its logical action is therefore exactly

\[
 V^*(e^{i\theta}g)^{\otimes4}V
   =e^{4i\theta}\rho_5(g),\qquad g\in G_{31}.
\]

The actual marked quarter operator
`R5=P3(I-2F_sigma)+i P2 tau` is outside this action, **even up to an arbitrary
unit scalar on the logical code**. Consequently it cannot be obtained by
applying one and the same single-register unitary to all four native registers
while preserving the actual code. This strictly strengthens the earlier
single-native-word exclusion within this stated class; it does not exclude
collective four-register operations. Moreover, section 5 below uses the exact
pair marginals to extend this exclusion to **all product-unitary lifts**, even
with unequal local factors. Convex mixtures of native group-word channels
cannot reproduce the logical quarter channel either.

## 1. Inputs and actual source objects

The prior native contract supplies the actual matrices and identifications:

- `G31` acts on the given `C^4`, contains the exact standard Pauli group
  `P64=mu4 {X^a Z^b tensor X^c Z^d}`, and normalizes it.
- The actual code is `Fix_P64(Sym^4 C^4)` and its induced `G31/P64=S6` action
  is the real standard five-dimensional representation in the fixed process
  frame. These are actual intertwiners, not a replacement embedding.
- Native `G31` is the finite complex reflection group in this Clifford chart,
  with order 46080 and scalar center `mu4`.

Original source:
`/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/joint-native-algebra-selection-20261016/PROOF.md`
and its `operator_checker.py`.

This audit imports the actual support constructor from
`/Users/stefanhamann/Projekte/tfpt-theoryv4/tfpt_explorer/process.py::_code_basis`.
It does not compute another 46080-element group census. The five disjoint
support sets have cardinalities `(4,12,12,12,24)`, with amplitude
`1/sqrt(cardinality)` on each support. Symmetry of the tensor degree is retained.

The actual mark is read through `build_marker_data`, its recorded polar
transport, its basis rotation, and `_process_frame`, just as in the original
native operator checker. No diagonal mark is substituted in the wrong frame.

## 2. All non-scalar continuous single-register motions leak out

For Hermitian `H`, put

\[
 D(H)=\sum_{j=1}^4 I^{\otimes(j-1)}\otimes H\otimes
 I^{\otimes(4-j)},\qquad
 L(H)=(I-VV^*)D(H)V.
\]

The Lie algebra of `K` consists exactly of `iH` with `L(H)=0`. This follows
by differentiating `exp(itH)^{tensor4}=exp(it D(H))`; conversely a Hermitian
`D(H)` preserving the code leaves it invariant under its exponential.

Let `H_a`, `a=1,...,15`, be the Hermitian nonidentity two-qubit Paulis
`IX,IY,IZ,XI,XX,XY,XZ,YI,YX,YY,YZ,ZI,ZX,ZY,ZZ`.
The exact computation on the original five sparse code columns gives

\[
 V^*D(H_a)V=0,\qquad
 \operatorname{tr}(L(H_a)^*L(H_b))=32\delta_{ab}.
\]

These equalities are computed with rational/Gaussian-integer arithmetic,
restoring the exact column norm factors `1/4,1/12,1/12,1/12,1/24`; no numerical
rank tolerance is used. The first identity says the derivative is wholly
outside the code for each traceless Pauli, so the second is just the Gram
matrix of the fifteen derivatives. Its coefficient 32 also has a compact
explanation, checked independently from the same sparse supports:

\[
 \operatorname{Tr}_{34}(VV^*)=(I_{16}+\mathrm{Swap})/4.
\]

Permutation symmetry gives the same identity for every register pair; its
one-register marginal is `5 I4/4`. Thus the four equal-site terms in
`tr[V* D(H_a)D(H_b) V]` contribute `20 delta_ab`, and the twelve distinct-site
terms contribute `12 delta_ab`, since `tr(H_a H_b)=4 delta_ab`. The normalized
maximally mixed code state has two-register marginal `P_sym/10`.
For any real coefficients `h_a`,

\[
 \|L(\sum_a h_a H_a)\|_{HS}^2=32\sum_a h_a^2.
\]

Thus no nonzero traceless Hermitian generator is in the Lie algebra. Scalar
generators have `D(hI)=4hI` and zero leakage. Hence

\[
 \operatorname{Lie}K=i\mathbb R I_4.
\]

This direct rank-15 calculation removes any need to assume an irreducibility
claim about the Pauli conjugation action in this step.

## 3. Finite components, with the center treated explicitly

`K` is a closed subgroup of compact `U(4)` because it is the inverse image of
the code-subspace stabilizer under the continuous fourth tensor-power map.
It is therefore a compact Lie group. Since scalar unitaries lie in `K`, its
identity component is the scalar circle. In particular

\[
 K_1=K\cap SU(4)
\]

is compact with zero Lie algebra and consequently finite.

Let `zeta8=exp(pi i/4)` and let `X=<G31,zeta8 I4>`. Under the actual native
Clifford identification, this is precisely the finite complex Clifford group
`X2` of Nebe--Rains--Sloane. Its center is `mu8`; the source identifies `G31`
as its index-two reflection subgroup, so `|X2|=92160`. Since `G31` preserves
the code and scalar matrices preserve every subspace, `X2 <= K`.

For clarity, this step uses the prior actual Clifford identification, as well
as the literature theorem; an arbitrary abstract group with the same order
would not suffice. The identification can be conjugated into this actual
Pauli chart before applying the theorem.

Every native reflection has determinant `-1`, while
`det(zeta8 I4)=zeta8^4=-1`; thus `det(X2)={+1,-1}`. Define

\[
 F=\{U\in K:\det U\in\{+1,-1\}\}.
\]

The determinant kernel is the finite `K1`, and its image on `F` has at most
two elements. Therefore `F` is finite, and `X2 <= F <= U(4)`.

**Primary theorem used.** Nebe, Rains and Sloane prove that every finite
unitary overgroup of their `Xm` is generated by `Xm` and scalar roots of
unity: [Theorem 6.5, printed p.18 / PDF page 18](https://arxiv.org/pdf/math/0001038v2#page=18).
Their `m=2` identification and orders are in
[Remark (2), printed p.19 / PDF page 19](https://arxiv.org/pdf/math/0001038v2#page=19).

Apply that theorem to `F`: every element of `F` is a scalar times an element
of `X2`. Every element of `K` can first be multiplied by a scalar so that its
determinant is `1`, and hence lies in `F`. Thus

\[
 K=U(1)F=U(1)X_2=U(1)G_{31}.
\]

This avoids the incorrect shortcut that a projectively finite group is itself
finite. The theorem is applied only to the explicitly finite determinant
slice `F`, not to `K` with its continuous scalar center.

## 4. The actual quarter operator fails even projectively

In the fixed process frame every `rho5(g)` is a real matrix. The original
polar-transported mark is real as well. Reconstructing the original operator
gives exactly

\[
 \operatorname{tr}R_5=1,\qquad
 \chi_{R_5}(t)=(t-1)^2(t+1)(t^2+1),
\]

\[
 \operatorname{Im}R_5=P_2\tau,\qquad
 \|\operatorname{Im}R_5\|_{HS}^2=2.
\]

Suppose `R5=lambda rho5(g)`, with `|lambda|=1`. The trace of the real standard
`S6` representation is `number_of_fixed_points(g)-1`, an integer. Taking the
trace therefore forces `lambda=+1` or `-1` (and respectively the group trace
`+1` or `-1`). Either value would make `R5` real, contradicting its nonzero
imaginary part. Therefore

\[
 R_5\notin U(1)\rho_5(G_{31}).
\]

Combined with the stabilizer theorem, there is no `U in U4` with
`U^{tensor4}V=V R5`, even after adding a logical scalar phase. This is not
merely a check that the central `J` and compiler `J sigma` have the wrong
orders; it covers all identical single-register unitary code-preserving
continuations of that actual source.

## 5. All product-unitary lifts and classical native mixtures

The pair marginal identity above strengthens the initial identical-register
classification. Let `T=U1 tensor U2 tensor U3 tensor U4` preserve the actual
code. As `T` is unitary, `T Pi T*=Pi` for `Pi=VV*`. Partial trace of this
identity over the other two registers gives, for every pair `(i,j)`,

\[
 (U_i\otimes U_j)(I+\mathrm{Swap})(U_i^*\otimes U_j^*)
    =I+\mathrm{Swap}.
\]

Here all six pair marginals are the same, by symmetry of the original code;
the checker now verifies every pair directly. Subtract the identity, use
`Swap(A tensor B)=(B tensor A)Swap`, and right-multiply by `Swap`. With
`A=Ui Uj*` this becomes

\[
 A\otimes A^*=I.
\]

Since `A` is unitary, its eigenvalues `a_r` have modulus one. The displayed
identity forces `a_r conjugate(a_s)=1` for all `r,s`, so every eigenvalue is
the same. A unitary matrix with a single eigenvalue is scalar. Thus
`Ui=phase_ij Uj` for every pair. Consequently

\[
 T=e^{i\Phi}U^{\otimes4},\qquad U\in K=U(1)G_{31}.
\]

Every product-unitary code preserver therefore acts logically as a scalar
times `rho5(g)`, `g in G31`. The actual quarter is outside this class by
section 4 and its source-frame qualification below. **Unequal local unitary
factors cannot repair the quarter lift.** This is a complete classification
of the code-preserving product-unitary class, rather than only a test of an
identical-register ansatz.

There is also a direct consequence for the classical native event menu.
Suppose a convex mixture of native group-word logical channels satisfied

\[
 \sum_g p_g\,\operatorname{Ad}_{\rho_5(g)}
       =\operatorname{Ad}_{R_5},\quad p_g\ge0,\quad \sum_gp_g=1.
\]

The Choi matrix on the right is the rank-one positive operator
`|R5>><<R5|`. Each positive-weight rank-one summand
`p_g |rho5(g)>><<rho5(g)|` on the left must have its support in that same
one-dimensional space. Thus every positive-weight `rho5(g)` would have to
be a unit scalar times `R5`, contradicting section 4. Hence no such mixture
equals the quarter channel **on all logical inputs**. This argument does not
rule out agreement on selected input states or protocols with further
coherent controls, conditioning or postselection.

## 6. Scope and the next source-level issue

The excluded classes are code-preserving product unitaries, including unequal
local factors, and convex mixtures of native group-word logical channels.
Collective entangling operations, controlled native events with ancillary
systems, general LOCC, postselection, code deformation, or operations leaving
and later returning to the code are not classified by this argument.

For example, the collective unitary
`V R5 V* + (I - VV*)` always extends `R5` on the full 256-dimensional space.
Its existence is algebraic and does not select a native physical interaction.
Likewise a freely declared `R5 tensor I4` on the decoded erasure carrier does
not establish a common native tensor-source operation. The remaining bridge
is an actually source-selected collective/ancillary implementation, or a
source-derived replacement for the stipulated marking/quarter relation.
No global physical `R`, P1 principle or general coherent realization has been
excluded by this stabilizer result.

## Reproducible check

`python3 /Users/stefanhamann/Documents/Codex/2026-10-16/er/work/marked_lift/transversal_checker.py`

Observed result: `PASS_EXACT_NATIVE_TRANSVERSAL_TANGENT_AND_PHASE_AUDIT`.
It writes `transversal_checker.json` beside the script. The finite-overgroup
theorem and the prior actual Clifford identification remain explicit proof
dependencies rather than purported outputs of this finite checker.

## Independent source-frame and erasure review, 2026-10-16

The subsequent source-clock audit exposed a necessary frame qualification to
the earlier operator reconstruction. Write `C=_process_frame()` (6 by 5) and
`N=G^(1/2) W6/sqrt(48)` (5 by 6). Exact independent checks give

\[
 N=OC^T,\qquad O=-C^T P_{\mathrm{reverse}}C,\qquad \det O=1.
\]

The unchanged actual `_source_operators` four-register permutation compresses
to `rho_C((1 3 2))`, while its original `N` label is `(2 3 4)`. The original
weak transposition `(0 1)_N` is `(4 5)_C`; `q*=5` in `N` labels is `0` in `C`
labels. Therefore transporting only `Y_N` to `Y_C=O^T Y_N O` while retaining
the nominal `N` permutations in `C` coordinates does not retain that unchanged
four-register clock. This qualification must accompany the earlier phrase
"actual quarter operator" in section 4.

There is an exact repair without changing the native four-clock: use `Y_N`
with all original `N` permutations. Its quarter operator is

\[
 R_{\rm source}=O R_{\rm previous\ C}O^T.
\]

Because `O=-rho_C(reverse)` belongs to `U(1)rho5(G31)`, this conjugation
preserves the admissible logical action class `U(1)rho5(G31)`. The stabilizer
theorem and the quarter nonlift result, including arbitrary logical scalar
phase, survive exactly. The checker now verifies the actual quartic image of
the unchanged source clock and this conjugation, as well as the original
trace and imaginary-part obstruction. No unchanged `sigma4` identification is
claimed for the previous mixed-frame expression.

The corrected root `marked_lift_checker.py` was read independently; its
source-consistent frame and its four weak-transposition fibers match the
equations above. Its same-polynomial erasure conclusion has the following
analytic explanation. Let `Q=I-F_sigma5` be the color-standard two-space and
take the +1 state of a color reflection,

\[
 \rho=Q(I+s_{23})/2.
\]

`P2 Q=0`; the reflection exchanges the `omega` and `omega^2` clock eigenlines,
so `rho` is a balanced pure state on them. The logical `R5` acts as identity
on this state. On the actual retained carrier the same native clock is
`sigma20=sigma5 tensor conjugate(sigma4)`, and
`F20=(I+sigma20+sigma20^2)/3`.

The unchanged `sigma4` has eigenvalue multiplicities `(1^2,omega,omega^2)`.
In the two gauge directions with eigenvalue 1, `F20` is zero on `Q` and
`I-2F20` leaves `rho` unchanged. In each of the two remaining gauge directions,
`F20` selects one of the two color-clock eigenlines. Its reflection sends
`rho` to the same orthogonal pure state. Thus the candidate unitary's output
differs from the correctly decoded output in two gauge directions, with
nonzero eigenvalues of their difference

\[
 (+1/4,+1/4,-1/4,-1/4).
\]

Consequently the difference has rank 4, Hilbert--Schmidt squared norm `1/4`
and trace distance `1/2`. All four primitive weak-transposition lifts give
this same failure because the witness has no `P2` component: their `tau4`
terms never act on it. The candidate is nevertheless unitary and has square
`I-2 P220`, since its `P3` reflection block and `P2` involution block are
orthogonal. Commutation between each `tau4` and `sigma4` is not required for
that block argument.

This accepts precisely the scoped distinction in the corrected root checker:
the decoder-derived operation is `R5 tensor I4`, whereas substituting the
native joint clock into the same polynomial need not preserve erasure
preparation. It does not identify the decoder operation with a primitive
physical intervention, nor exclude other coherent common-source operations.

The root checker's three polar `q*` choices are tested using the same source
slot partition, positive polar factor and actual `sigma4`-centralizing
reflection lifts. Those checks establish one orbit of this specific mark
family under the native clock centralizer, rather than invariance of every
anisotropic process or physical state. No issue was found in this stated
scope. The root checker was not edited or duplicated during this review.

### Closed form of the polar mark

The root's subsequent closed form is also correct with the actual simplex
columns `n_j=N[:,j]`, whose squared norms are `5/6` and mutual inner products
are `-1/6`. Write `m=n2+n3+n4`, `q in {0,1,5}`, and

\[
 a=(1/\sqrt6-1)/5,\quad D_q=(n_j-n_q)_{j\ne q},\quad
 U_q=D_q(I+aJ).
\]

For the unit color-constant slot vector `v=sum(e_color)/sqrt3`,
`D_q v=(m-3n_q)/sqrt3`, `D_q 1=-6n_q`, and `1^T v=sqrt3`. Hence

\[
 c_q=U_qv=\frac{m-3(\sqrt6-1)n_q/5}{\sqrt3}.
\]

The image of the three color slots is the orthogonal sum of the color-standard
two-space `Q=I-F_sigma` and this unit vector `c_q` in `Fix(sigma)`. Consequently

\[
 P_2(q)=F_\sigma-c_qc_q^*.
\]

This derives the specified marked projector from the clock, simplex, fixed
point and stipulated positive polar prescription. No Standard Model charges
are needed in this identity. The positive polar prescription, choice of fixed
point and physical meaning of the mark remain declared inputs; the formula
does not select them. The derivation explains why `P2 Q=0` in the erasure
witness argument above.
