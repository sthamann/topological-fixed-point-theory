/-
  GramCompactness — the finite Gram compactness theorem.
  =======================================================

  Machine-checked abstract closure of the star architecture (Modul 5
  of the GLOBAL-HANDOFF-OFFENSIVE review specification): pure linear
  algebra / functional analysis over ℝ, no number theory anywhere.

  THE THEOREM.  For matrices of a FIXED finite size:

    (a) `posSemidef_of_tendsto` — a pointwise (entrywise) limit of
        positive-semidefinite matrices is positive semidefinite:
        G_n ⪰ 0, G_n → G entrywise ⟹ G ⪰ 0.  Proof: quadratic
        forms are finite sums of products of entries, so they
        converge (`Filter.Tendsto` algebra); nonnegativity passes to
        the limit (`ge_of_tendsto'`); hermiticity passes to the limit
        by uniqueness of limits.  Stated for an arbitrary filter with
        `NeBot`, so it covers sequences (atTop) and nets alike.

    (b) `gram_posSemidef` — the Gram matrix G[i,j] = ⟪v i, v j⟫ of
        ANY finite family of vectors in ANY real inner product space
        is positive semidefinite: xᵀGx = ‖Σ_i x_i·v_i‖² ≥ 0
        (`gram_dotProduct_eq` gives the exact identity
        xᵀGx = ⟪Σ x_i v_i, Σ x_j v_j⟫).

    (c) `posSemidef_limit_of_gram` — the assembly the handoff theorem
        needs: if the Gram matrices of a sequence of vector families
        converge entrywise, the limit matrix is PSD (b) + (a).

    (d) `exists_subseq_posSemidef_limit` — the sequential-compactness
        step: an entrywise-BOUNDED sequence of PSD matrices has a
        subsequence converging entrywise to a PSD limit
        (Bolzano–Weierstraß via `tendsto_subseq_of_bounded` on the
        finite product space; the closed ball in the sup norm is
        compact since the product is proper).  Corollary
        `exists_subseq_gram_limit`: same for bounded Gram sequences.

  DESIGN.  Everything is stated over ℝ with Mathlib's
  `Matrix.PosSemidef` (over ℝ: symmetry + nonnegative quadratic
  form; `star` is trivial).  The compactness step (d) runs on the
  underlying finite product space `ι → ι → ℝ` (definitionally equal
  to `Matrix ι ι ℝ`), where Mathlib provides the proper-space and
  Bolzano–Weierstraß infrastructure; matrices carry the same (product)
  topology, so entrywise convergence transfers by `tendsto_pi_nhds`.

  HONEST SCOPE.  Formalized: (a)–(d) in full, over an arbitrary
  finite index type and an arbitrary real inner product space —
  no `sorry`, no `native_decide`, no citation gaps.  NOT formalized
  (outside the module's brief): the infinite-dimensional /
  operator-topology versions (weak-* compactness of state spaces),
  the diagonal extraction along a COUNTABLE basis of a separable
  space (here the index set is finite — the finite Gram compactness
  is exactly the module-5 specification), and any connection to the
  number-theoretic window data of the star architecture (this file
  is deliberately free of it).
-/
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Topology.MetricSpace.Sequences
import Mathlib.Tactic

set_option maxRecDepth 100000
set_option maxHeartbeats 1000000

open Filter Topology Matrix

namespace TfptCarrier.GramCompactness

variable {ι : Type*} [Fintype ι]

/-! ### (a) PSD passes to entrywise limits -/

/-- Quadratic form of a real matrix as an explicit finite double sum. -/
theorem dotProduct_mulVec_expand (M : Matrix ι ι ℝ) (x : ι → ℝ) :
    x ⬝ᵥ (M *ᵥ x) = ∑ i, ∑ j, x i * M i j * x j := by
  simp [dotProduct, Matrix.mulVec, Finset.mul_sum, mul_assoc]

/-- **(a) The PSD limit theorem**: an entrywise limit of positive
semidefinite real matrices (along ANY nontrivial filter — sequences,
nets) is positive semidefinite.  The quadratic form is a finite sum
of products of entries, hence converges; nonnegativity and symmetry
pass to the limit. -/
theorem posSemidef_of_tendsto {α : Type*} {l : Filter α} [l.NeBot]
    (f : α → Matrix ι ι ℝ) (G : Matrix ι ι ℝ)
    (hf : ∀ a, (f a).PosSemidef)
    (hlim : ∀ i j, Tendsto (fun a => f a i j) l (𝓝 (G i j))) :
    G.PosSemidef := by
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg ?_ ?_
  · -- symmetry passes to the limit by uniqueness of limits
    show G.conjTranspose = G
    ext i j
    have hsym : ∀ a, f a i j = f a j i := by
      intro a
      have h := congrFun (congrFun ((hf a).1) i) j
      simpa using h.symm
    have h1 : Tendsto (fun a => f a j i) l (𝓝 (G i j)) := by
      refine (hlim i j).congr fun a => hsym a
    simpa [Matrix.conjTranspose_apply] using
      (tendsto_nhds_unique h1 (hlim j i)).symm
  · intro x
    rw [star_trivial]
    have hquad : Tendsto (fun a => x ⬝ᵥ (f a *ᵥ x)) l (𝓝 (x ⬝ᵥ (G *ᵥ x))) := by
      simp only [dotProduct_mulVec_expand]
      exact tendsto_finset_sum _ fun i _ =>
        tendsto_finset_sum _ fun j _ =>
          ((hlim i j).const_mul (x i)).mul_const (x j)
    refine ge_of_tendsto' hquad fun a => ?_
    have h := (hf a).dotProduct_mulVec_nonneg x
    rwa [star_trivial] at h

/-! ### (b) Gram matrices are PSD -/

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- The Gram matrix of a finite family of vectors in a real inner
product space: G[i,j] = ⟪v i, v j⟫. -/
def gram (v : ι → E) : Matrix ι ι ℝ :=
  Matrix.of fun i j => inner ℝ (v i) (v j)

omit [Fintype ι] in
theorem gram_apply (v : ι → E) (i j : ι) :
    gram v i j = inner ℝ (v i) (v j) := rfl

/-- The Gram quadratic form IS an inner square:
xᵀGx = ⟪Σ_i x_i·v_i, Σ_j x_j·v_j⟫. -/
theorem gram_dotProduct_eq (v : ι → E) (x : ι → ℝ) :
    x ⬝ᵥ (gram v *ᵥ x)
      = inner ℝ (∑ i, x i • v i) (∑ j, x j • v j) := by
  rw [sum_inner]
  simp only [inner_sum, real_inner_smul_left, real_inner_smul_right,
    dotProduct, Matrix.mulVec, gram, Matrix.of_apply, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by
    ring

/-- **(b) Gram positivity**: the Gram matrix of ANY finite family in
ANY real inner product space is positive semidefinite —
xᵀGx = ‖Σ_i x_i·v_i‖² ≥ 0. -/
theorem gram_posSemidef (v : ι → E) : (gram v).PosSemidef := by
  refine Matrix.PosSemidef.of_dotProduct_mulVec_nonneg ?_ ?_
  · show (gram v).conjTranspose = gram v
    ext i j
    simpa [gram, Matrix.conjTranspose_apply] using real_inner_comm (v i) (v j)
  · intro x
    rw [star_trivial, gram_dotProduct_eq]
    exact real_inner_self_nonneg

/-! ### (c) The assembly: limits of Gram matrices are PSD -/

/-- **(c) The handoff corollary**: if the Gram matrices of a sequence
of vector families converge entrywise, the limit matrix is positive
semidefinite — [Gram_n ⪰ 0 by (b)] + [entrywise convergence] ⟹
[limit ⪰ 0 by (a)]. -/
theorem posSemidef_limit_of_gram (v : ℕ → ι → E) (G : Matrix ι ι ℝ)
    (hlim : ∀ i j,
      Tendsto (fun n => inner ℝ (v n i) (v n j)) atTop (𝓝 (G i j))) :
    G.PosSemidef :=
  posSemidef_of_tendsto (fun n => gram (v n)) G
    (fun n => gram_posSemidef (v n)) hlim

/-! ### (d) Sequential compactness: bounded PSD sequences subconverge -/

/-- **(d) Finite Gram compactness / Bolzano–Weierstraß**: an
entrywise-bounded sequence of PSD matrices of fixed finite size has
a subsequence converging entrywise to a PSD limit.  The sequence
lives in a closed ball of the (proper) finite product space, hence
subconverges; the limit is PSD by (a). -/
theorem exists_subseq_posSemidef_limit (A : ℕ → Matrix ι ι ℝ) (C : ℝ)
    (hPSD : ∀ n, (A n).PosSemidef)
    (hbound : ∀ n i j, |A n i j| ≤ C) :
    ∃ (G : Matrix ι ι ℝ) (φ : ℕ → ℕ), StrictMono φ ∧ G.PosSemidef ∧
      ∀ i j, Tendsto (fun n => A (φ n) i j) atTop (𝓝 (G i j)) := by
  -- run Bolzano–Weierstraß on the underlying finite product space
  set B : ℕ → ι → ι → ℝ := fun n i j => A n i j with hB
  have hmem : ∀ n, B n ∈ Metric.closedBall (0 : ι → ι → ℝ) (max C 0) := by
    intro n
    rw [Metric.mem_closedBall, dist_zero_right]
    refine (pi_norm_le_iff_of_nonneg (le_max_right _ _)).mpr fun i => ?_
    refine (pi_norm_le_iff_of_nonneg (le_max_right _ _)).mpr fun j => ?_
    rw [Real.norm_eq_abs]
    exact le_trans (hbound n i j) (le_max_left _ _)
  obtain ⟨g, -, φ, hφ, hconv⟩ :=
    tendsto_subseq_of_bounded Metric.isBounded_closedBall hmem
  have hij : ∀ i j, Tendsto (fun n => A (φ n) i j) atTop (𝓝 (g i j)) := by
    intro i j
    have h1 := tendsto_pi_nhds.mp hconv i
    have h2 := tendsto_pi_nhds.mp h1 j
    simpa [hB, Function.comp] using h2
  exact ⟨Matrix.of g, φ, hφ,
    posSemidef_of_tendsto (fun n => A (φ n)) (Matrix.of g)
      (fun n => hPSD (φ n)) hij, hij⟩

/-- **(d) for Gram sequences**: a sequence of vector families with
entrywise-bounded Gram matrices has a subsequence along which the
Gram matrices converge entrywise to a PSD limit (by Cauchy–Schwarz a
uniform norm bound ‖v n i‖ ≤ B gives exactly such an entry bound
with C = B²). -/
theorem exists_subseq_gram_limit (v : ℕ → ι → E) (C : ℝ)
    (hbound : ∀ n i j, |inner ℝ (v n i) (v n j)| ≤ C) :
    ∃ (G : Matrix ι ι ℝ) (φ : ℕ → ℕ), StrictMono φ ∧ G.PosSemidef ∧
      ∀ i j, Tendsto (fun n => inner ℝ (v (φ n) i) (v (φ n) j)) atTop
        (𝓝 (G i j)) :=
  exists_subseq_posSemidef_limit (fun n => gram (v n)) C
    (fun n => gram_posSemidef (v n)) hbound

end TfptCarrier.GramCompactness

#print axioms TfptCarrier.GramCompactness.exists_subseq_posSemidef_limit
