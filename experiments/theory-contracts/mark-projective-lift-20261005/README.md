# SOURCE.MARK_PROJECTIVE_LIFT.01 — A4 or 2T for the oriented mark rotations?

Date 2026-10-05. Research contract under the experiments firewall: no ledger, paper or website change, no promotion.
Run: `python3 -B experiments/theory-contracts/mark-projective-lift-20261005/lift.py` (8 steps, about 2 s).

**Question.** No spin structure is inserted. Do fermionic grading and the oriented marks force the nontrivial lift
1 → Z₂ → 2T → A4 → 1, with central element (−1)^F, or the trivial linear action?

**Verdict: NOT_FORCED_IN_SOURCE; SPIN_HALF_EMERGENT_GIVEN_TWO_COMPONENT_COIN.**

**What is shown.**
- **Invariant.** The test quantity is the rephasing-invariant commutator c of the lifts of (01)(23) and (02)(13).
  - The std and regular representations give c = +1.
  - The SU(2) (2T) doublet gives c = −1.
- **Fock-implementation lemma.** Take a *genuine* orthogonal one-particle action O of A4 on Majoranas.
  - Its Clifford implementers satisfy every A4 relation up to *scalars*.
  - Their class is w₂(O) = m₃(O) mod 2, where m₃ is the multiplicity of the 3-dim irrep.
  - The central element is never (−1)^F.
  - This generalises the 20260921 lemma that U² = parity needs O² = −I.
  - Number-conserving actions (complex modes, m₃ even) are trivial.
  - The γ-only rotation of the three family Majoranas has c = −1, but on *both* parity sectors. It breaks fermion number, mixes families (Part G) and violates spin-statistics.
  - Twist zero modes give c = (−1)^N; checked for N = 1…4. The source has N = 16, so c = +1.
- **Positive control.** If the one-particle space itself carries the 2T doublet, second quantisation gives c = (−1)^F. Spin-statistics therefore comes from a projective *one-particle* action, never from the Fock implementation.
- **Seam geometry.**
  - The Möbius stabiliser of the square marks (λ = −1) is D4.
  - Only the four rotations preserve the orientation of the seam circle.
  - With A4 the seam shares only ⟨(02)(13)⟩. No 3-cycle is a seam symmetry; the A4-symmetric equianharmonic configuration is not concyclic.
  - The NS lift of the seam rotation, U(2π) = (−1)^F, is conditional on the S² spinor-section premise. Even then it reaches only ⟨(02)(13)⟩ and gives no V4 commutator.
- **Emergent spin ½ (not an IR derivation of spin-statistics).**
  - At the isotropic Weyl node of h_W, the emergent rotations act on a two-component coin necessarily through SU(2), with covariance error O(k²).
  - The 2π loop lifts to −I. This *coincides* with the fermion parity the CAR/Fock structure already has.
  - It does not derive spin-statistics from IR geometry. It also cannot be used to derive AX.WEYLCOIN, because a Weyl node presupposes the two-dimensional coin.
  - The P1 sheet doublet is no coin: all four marks act by the same swap, so its Gram matrix has rank 1.

**Consequence.** For the chiral Θ₁ (residual C3, H²(C3,U(1)) = 0) the microscopic cocycle question is void. The remaining
input is the **existence of a two-component internal-singlet one-particle coin** coupled to the three H₁ directions
(AX.WEYLCOIN). Once it exists and the IR kinetics has an isotropic isolated Weyl node, spin ½ emerges and no class bit has
to be chosen. Whether the source contains such a coin is tested in `weylcoin-origin-20261005`: it does not. This sharpens `minimal-runtime-reconstruction-20261005`, which still listed a
postulated class bit.

**Imported premises.**
- The kills from `source-oriented-polarization-20260921`.
- Part U6/K9: N = 16 species-singlet zero modes.
- The S² spinor-section premise. It is used only for the conditional NS statement.

Files: `lift.py`, `lift_results.json`, `contract_index.json`.
