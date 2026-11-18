# Native operator audit: one actual 4/5/20 chain, its clocks, and its remaining choices

**Verdict.** There is a precise shared object, but it is representation-theoretic rather than a single clock:

\[
H_8\xrightarrow{\text{Construction A}}(E_8,J)
\longrightarrow G_{31}\curvearrowright V_4
\longrightarrow V_5=(\operatorname{Sym}^4V_4)^{P_{64}}
\longrightarrow V_{20}=V_5\otimes\overline {V_4}.
\]

All arrows above use existing constructors.  The audit corrects one important typing error: the fixed space of `P64` in the **full** tensor power `V4^tensor4` has dimension **16**; the dimension-five carrier is the fixed space inside **`Sym^4(V4)`**.  The actual sparse code lies in that symmetric subspace and is pointwise fixed by `P64`.

The chain jointly forces the block algebra `M4(C) direct-sum M5(C) direct-sum M20(C)`, its commutant `C^3`, the erasure isometry, a canonical orbit of fifteen `2+3` response splittings, and a sharply bounded normalizer.  The physical rank-two mark `P2` is **not** a member of that fifteen-projector orbit, so the orbit does not derive or explain the actual marked `3+2` split.  It also does not select the marked order-four `R`.  The carrier-five/Galois structure is stronger: the six literal copies left by the point mark form one orbit under the already marked stabilizer, hence one equivalence class rather than six independent physical choices.

## 1. Actual constructors and the corrected Pauli statement

`P64=mu4*P2` is constructed as 64 exact `4 by 4` matrices in `quartic_selection.py:49-61`.  The symmetric-power characters are computed by Newton recursion at lines 64-75 and averaged at lines 153-190.  The exact result is

\[
\dim (\operatorname{Sym}^dV_4)^{P_{64}}=(1,0,0,0,5),\qquad 0\le d\le4.
\]

Thus degree four is the first nonconstant symmetric invariant degree and has dimension five.  By contrast,

\[
\dim (V_4^{\otimes4})^{P_{64}}
=\frac1{64}\sum_{p\in P_{64}}\operatorname{tr}(p)^4=16.
\]

The actual five sparse vectors are the five permutation-orbit sums in `process.py:242-258`.  The checker verifies directly, with exact monomial phases, that every one of the 64 Paulis fixes all five vectors and that every support is closed under all four tensor-slot permutations.  No five-dimensional subspace was chosen from the 16-dimensional full-tensor fixed space after the fact.

The same existing process applies the 60 native reflections in four tensor slots and reduces them to fifteen matrices (`process.py:336-380`).  These are the standard five-dimensional `S6` transpositions.  In group terms, the exact native structure calculation finds `|G31|=46080`, its central `mu4`, the normal 64-subgroup, and

\[
G_{31}/P_{64}\cong\operatorname{Sp}_4(2)\cong S_6
\]

in `v634_st31_structure.py:238-278,561-673`.  The compiler elements satisfy `c=J sigma`, `ord(c)=12`, `sigma=c^4`, `J=c^9` at lines 724-737.

### The native `3+2` structure is an orbit, not a unique split

The original process enumerates the fifteen perfect matchings of six labels and sets

\[
\tau_m=(ab)(cd)(ef),\quad P_m=(I+\rho_5(\tau_m))/2,
\quad Y_m=(5P_m-2I)/6
\]

(`process.py:384-393,698-722`).  A fixed-point-free involution has trace `-1` on the standard `S6` module, so its multiplicities are `+1^2,-1^3`; hence every `P_m` has rank two and every `Y_m` has the claimed `2+3` spectrum.  The fifteen matchings form one canonical conjugacy orbit.  This is a positive intrinsic `3+2` family.  It does not mark one member as physical; `process.py:721` states that selection boundary explicitly.

A double transposition also yields a `3+2` decomposition, but with `+1^3,-1^2`, so `(I+tau)/2` has rank three.  The original response constructor singles out the perfect-matching class, not every involution with the same unordered dimensions.  The native quotient also carries the actual fifteen nonidentity Pauli labels.  The perfect matchings are its companion fifteen-object `S6` orbit (the classical duad/syntheme or outer-automorphism correspondence once the existing `Sp4(2)=S6` identification is fixed).  This is a real Pauli/response correspondence of families.

It is **not**, however, the actual physical mark.  The checker compares

\[
P_2=(6Y_5+2I)/5
\]

exactly with all fifteen `Q_t=(I+rho5(tau_t))/2`; all fifteen equalities are false.  There is also a basis-free proof.  `P2` commutes with the color subgroup `S3_(234)`.  If `P2=Q_t`, then `rho5(tau_t)=2P2-I` would commute with that subgroup.  Faithfulness of the standard `S6` representation would force the fixed-point-free involution `tau_t` into

\[
\operatorname{Cent}_{S_6}(S_{3,(234)})=S_{3,(015)}.
\]

That centralizer has order six and contains no fixed-point-free involution: its involutions move only two of the complementary labels.  Therefore the canonical fifteen-projector family supplies possible native `2+3` response planes, but it cannot be used as the origin of the already marked physical `P2`.  The common rank is only a shared dimension.

## 2. Actual erasure map and the 20-dimensional carrier

The code supports in `process.py:242-258` give, without a basis choice beyond the stored code basis, the sparse Schmidt isometry

\[
W:V_5\otimes\overline V_4\longrightarrow V_4^{\otimes3},
\quad W^*W=I_{20}.
\]

Its explicit row/column constructor is `hyper_common_reflection_check.py:67-77`; the actual native test is at lines 174-181.  For the existing marked color triangle, native ray IDs `(1,21,32)` compress to the logical transpositions `(23),(24),(34)` and satisfy exactly

\[
r_g^{\otimes3}W=W\bigl(C_g\otimes\overline{r_g}\bigr)
\]

for all three generating reflections (`hyper_mark_compatible_search.py:118-171`).  The associated erasure preparation is

\[
E(\rho)=\rho\otimes I_4/4,
\qquad E_{\rm phys}(\rho)=W E(\rho)W^*.
\]

This is an actual erasure/Schmidt factorization, rather than a `20=5*4` dimension match.

On `V20`, `P64` acts trivially on `V5` and irreducibly on `bar V4`; therefore its commutant is `M5 tensor I4`.  Intersecting with the irreducible standard `S6` action on the multiplicity space leaves only scalars.  Thus `V20` is irreducible.  Together with the faithful irreducible `V4` and the irreducible quotient module `V5`, Burnside's theorem gives

\[
\mathbb C[\rho_4\oplus\rho_5\oplus\rho_{20}](G_{31})
=M_4\oplus M_5\oplus M_{20},\qquad
\rho(G_{31})'=\mathbb C^3.
\]

The kernels are `ker rho4=1`, `ker rho5=P64`, `ker rho20=1`.  Central `J=iI4` acts as

\[
(iI_4,\ I_5,\ -iI_{20}).
\]

Hence the three block projectors are already polynomials in the actual central operator.  This classifies a hidden block degree of freedom; it does not choose a dissipative rate or menu of Pauli-gated operations.

## 3. Joint normalizer: proof and exact scope

The group calculation gives a normal 2-subgroup `P64` and quotient `S6`.  Since `O2(S6)=1`, every normal 2-subgroup of `G31` maps trivially to the quotient and lies in `P64`.  Therefore

\[
O_2(G_{31})=P_{64}.
\]

`O2` is characteristic, so every automorphism or unitary normalizer of `G31` preserves `P64`.  The unitary normalizer of the two-qubit Pauli group is, by definition, the two-qubit Clifford group.  Projectively it has order `11520`; the native `G31/mu4` has the same order and is already contained in it.  Consequently

\[
N_{U(4)}(G_{31})=U(1)G_{31}.
\]

For the joint representation, the dimensions `4,5,20` are distinct, so a normalizer preserves the three irreducible blocks.  Its faithful `V4` block induces an inner automorphism, by the previous equality.  Removing the corresponding common `rho(g)` leaves an operator commuting with `G31` on each irreducible block, hence a phase by Schur's lemma.  Therefore

\[
N_{U(4)\times U(5)\times U(20)}(\rho(G_{31}))
=U(1)^3\rho(G_{31}).
\]

This is the block-preserving normalizer of the actual diagonal group action.  It is not a statement about arbitrary unitaries of the full matrix algebra, nor about spacetime or gauge symmetry.

## 4. Marked `R`, compiler clock, Galois carrier, and the old sign twin

### Marked `R` is an extra coherent operator

The actual marked completion constructs

\[
R=P_3(I-2F_\sigma)+iP_2\tau,
\qquad F_\sigma=(I+\sigma+\sigma^2)/3
\]

in `marked_rr_completion.py:26-66`.  The same file explicitly types the first extra operation as a coherent polynomial/linear combination and says that neither clock nor generator is selected from the source (`:120-137`).  Its characteristic polynomial is

\[
(x-1)^2(x+1)(x^2+1),
\]

so its spectrum is `{1,1,-1,i,-i}` and its trace is `1`.  This corrects the tentative trace-three spectrum.  It has the same spectrum as an order-four element in the standard `S6` representation, so spectrum alone cannot exclude identification.

There is a short exact exclusion.  Direct symbolic commutators show that `R` centralizes the marked subgroup

\[
H=S_{\{0,1\}}\times S_{\{2,3,4\}}\cong S_2\times S_3.
\]

But `Cent_S6(H)={1,(01)}` has order two.  If `R=rho5(p)` for an actual `p in S6`, faithfulness of the standard representation would force `p` into this centralizer, hence `ord(R)<=2`, contradicting `ord(R)=4`.  Therefore `R` is not an actual native group word.  The 720-element fixed-frame census independently agrees: no matrix equals `R`; the nearest Frobenius distance is `2.6462860671`.  Its centralizer in the actual `S6` image has size 12.

The native compiler clock behaves differently.  Since `J^tensor4=I` on degree four, `c=J sigma` maps to `sigma5`, which has order three.  Thus the native compiler `c` and the marked `R` are different actual actions, even though both live in the generated algebra.

### The carrier-five Galois relation exists natively as one marked-symmetry class

The standalone cyclotomic carrier uses the rational companion `C5` and Frobenius `G` with `G C5 G^-1=C5^2`, `G^4=I` (`v419_seam_galois_carrier.py:41-44,80-93`).  The new smallest native test fixes the existing marked `S6` label 5.  On the other five labels there are six Sylow-five subgroups.  For **each of all six**, the checker constructs exact native `4 by 4` matrices satisfying

\[
C^5=I,\qquad G^4=I,\qquad GCG^{-1}=C^2.
\]

For the first subgroup one explicit pair is

\[
C={1\over2}\begin{pmatrix}
-i&1&-1&-i\\-i&-1&1&-i\\-i&-1&-1&i\\-i&1&1&i
\end{pmatrix},\quad
G={1\over2}\begin{pmatrix}
-i&-1&-i&1\\1&i&-1&i\\-1&i&1&i\\-i&1&-i&-1
\end{pmatrix}.
\]

The audit then solves, rather than chooses, an intertwiner `A` to the rational `v419` matrices:

\[
C_{\rm native}A=A C_{\rm v419},\qquad
G_{\rm native}A=A G_{\rm v419}.
\]

For every one of the six copies the joint linear system has nullity one and its nonzero solution has rank four.  Thus the intertwiner is unique up to an irrelevant nonzero scalar.  For the displayed first pair, one normalization is

\[
A=\begin{pmatrix}
1&1-2i&-1&i\\
1&-1&-i&-1+2i\\
-1+2i&-i&-1&1\\
1&i&-2-i&-i
\end{pmatrix}.
\]

So the proposed central-phase obstruction is false: honest `F20=C5 semidirect C4` subgroups exist in the native carrier, and every surviving copy is actually equivalent to the standalone cyclotomic carrier.

The six are also **not six inequivalent physical choices**.  Let

\[
H_0=S_{\{0,1\}}\times C_{3,(234)},\qquad
H=S_{\{0,1\}}\times S_{\{2,3,4\}}.
\]

On the **logical `V5` carrier**, both preserve the existing point/weak/color marks and `sigma5`; `H` also centralizes the marked `R`.  Exact conjugation on the six Sylow-five subgroups shows that `H0` acts regularly (orbit six, stabilizer one), while `H` acts transitively (orbit six, stabilizer two).  Therefore the logical marks determine one `F20` **orbit/equivalence class**.  No additional embedding has been chosen.  A literal subgroup and a normalized intertwiner matrix are gauge representatives; the remaining interpretive question is whether the marked stabilizer is quotiented as redundancy or later broken by further physical data.  This quotient-orbit result does not by itself prove that the conjugating lifts commute with the native `V4` representative of `sigma` or with a P1 physical source clock; that stronger simultaneous-lift statement was not used.

### Which constraints eliminate the old `S3` sign twin?

For the actual `(1,21,32)` triangle, restriction gives

\[
V_4\downarrow S_3=2\,\mathbf1\oplus\mathbf2,
\quad V_5\downarrow S_3=3\,\mathbf1\oplus\mathbf2,
\quad V_{20}\downarrow S_3=7\,\mathbf1\oplus\mathbf{sgn}\oplus6\,\mathbf2.
\]

The central sign projector has ranks `(0,0,1)`, so `z=1-2e_sgn` is invisible on `V4,V5` and a nontrivial involution on `V20` (`hidden_sgn_check.py:29-63`).  The actual tests sharpen its status:

- **Marked `R` does not eliminate it.** `z20` commutes exactly with `R tensor I4`.
- **The genuine six `S3` event branches preserve erasure covariance.** This is the eventwise identity already proved in `audit_source_first.py:125-138`.
- **Branchwise preservation of the actual erasure preparation eliminates it.** For an explicit rank-one state in the actual color-standard subspace,
  `Ad_z(E(rho))-E(rho)` has rank two, Hilbert-Schmidt norm squared `1/8`, and trace distance `1/4`.  Thus `z` cannot be inserted as an additional branch while keeping the same erasure map and the identity logical action.
- **Full `G31` covariance also eliminates the projector.** `V20` is irreducible, so a rank-one `S3`-central projector is not in the full commutant.
- **The primitive event menu eliminates it by definition.** `z` is a coherent group-algebra unitary, not one of the three native reflections.  Source normalization alone does not eliminate it because `z` is unitary.

The clean discriminator is therefore the full native covariance or the branchwise erasure contract, not clock commutation.

## 5. Hamming/E8/J provenance, with real maps only

The Hamming/E8 entrance uses no complex or physical-clock identification.  In doubled real coordinates, `J` is the real orthogonal block rotation

\[
(x_0,x_1,x_2,x_3,x_4,x_5,x_6,x_7)
\mapsto(-x_1,x_0,-x_3,x_2,-x_5,x_4,-x_7,x_6).
\]

It is implemented in `v689_gaussian_code_bridge.py:180-189` and `v634_st31_structure.py:90-120`.  Construction A first produces the 240 real `E8` roots.  Exact checks show that `J` and the real coordinate permutation `sigma` preserve them, commute, and satisfy `J^2=-I` (`v689:288-306`).  Only after that does one regard the real lattice as a rank-four `Z[i]` module.  Reduction modulo `(1+i)` gives `F2^4` (`v689:309-338`); the 240 roots occupy fifteen nonzero classes of sixteen roots, and each class is four of the 60 `J`-orbits (`v689:359-407`).

This is the legitimate Hamming/E8/J bridge.  Calling `J` multiplication by `i` is a repackaging of the verified real orthogonal map.  It does not identify `J` with physical time, the marked `R`, or the separate carrier-five Frobenius.

## Reproduction and pinned provenance

Run:

```text
python3 work/joint_native/operator_checker.py
```

It writes `operator_results.json` and returns `PASS_EXACT_NATIVE_OPERATOR_AUDIT`.  All structural identities are exact SymPy equalities.  Only the finite 720-way `R` equality/commutator census is evaluated in complex128 after constructing the 720 algebraic matrices exactly; its smallest nonzero residual is `2.5988895521`, far above the `1e-12` decision threshold.

Repository HEAD observed during the audit: `d44f1040d733338e1c0596433b5b2578d7c00295`.  The working tree contained unrelated changes, so the content hashes below, not HEAD alone, pin the read sources.

| Source and decisive lines | SHA-256 |
|---|---|
| `tfpt_explorer/process.py:242-258,336-393,698-722` | `e8db5fcd999c9d4c34fbee548aa2c85e1a515cf03663a62355533439ae0f1c94` |
| `tfpt_explorer/quartic_selection.py:49-75,153-205` | `122ed509b61637e9937801207cc1188f76d7112cf472c61825912c5ec5ea9296` |
| `tfpt_explorer/marker_selection.py:30-58,66-120,385-446` | `f6d9907f0e9fa242c4f50d67681c499d7bff046f0032abed37f66f332965842a` |
| `verification/v634_st31_structure.py:90-163,238-278,561-673,724-737` | `5ff53ed031ca03103187ac22edf723aabce88bee3df1085c52486ecb0044d37c` |
| `verification/v689_gaussian_code_bridge.py:95-103,180-230,288-407` | `9d08b59aa5ada89232894762398f54ed690a08d170e619f79cb80fb926b85979` |
| `verification/v419_seam_galois_carrier.py:41-44,80-93` | `cc6152ad7ea68c67725020e386e49eec32810ab25bb483c76bae99964399c2fd` |
| `hyper_common_reflection_check.py:37-77,92-181` | `cabe085e7ed5e0c3c9fd2b05b4b5a96277a62e4cb7117724f0abd41ca916e9b4` |
| `hyper_mark_compatible_search.py:52-171` | `971c389cfb4733f9ea7760b3d9432461ffc8f618f4eae431e4aebd321e5fd45c` |
| `audit_source_first.py:61-138` | `b153c3fdce33360f89c4277532d983825a2618a1d5af53b3e4e9c94338507e06` |
| `hidden_sgn_check.py:29-101` | `81e0b41e4bc23acff81a3716e5999f39f8bfb7ca83002efd7ba57ef6dd3c69fa` |
| `marked_rr_completion.py:26-66,120-137` | `f1c4e077f93e0462518f27d5abdaca3767f4c407b4a687285c0e12e3f243d6d6` |

The three October 8 source-audit scripts and the marked-`R` script are historical exact audits outside the TFPT repository; the table pins their literal content.  The native constructors they call are pinned separately above.
