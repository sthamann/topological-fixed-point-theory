import Mathlib

/-!
# Algebraic source-reconstruction lemmas

The results below are exact identities in an abstract starred ring and over the
rational numbers.  They do not construct a Hilbert-space representation, prove
closed-range/domain statements, or derive the displayed data from TFPT P1/P2.
-/

namespace TFPT.SourceReconstruction

section Compression

variable {A : Type*} [Ring A] [StarRing A]

/-- Algebraic compression of `U` by `C`. -/
def compression (C U : A) : A := C * U * C

/-- The exact term by which the compression fails the corner-isometry equation. -/
def leakage (C U : A) : A := C * star U * (1 - C) * U * C

/--
Exact compression identity.  `C` is assumed to be a self-adjoint idempotent and
only the left-unitary relation `U⋆U=1` is used.
-/
theorem compression_star_mul_self
    (C U : A)
    (hCstar : star C = C)
    (hCid : C * C = C)
    (hU : star U * U = 1) :
    star (compression C U) * compression C U = C - leakage C U := by
  calc
    star (compression C U) * compression C U =
        C * star U * (C * C) * U * C := by
          simp only [compression, star_mul, hCstar]
          noncomm_ring
    _ = C * star U * C * U * C := by rw [hCid]
    _ = C - leakage C U := by
      simp only [leakage]
      noncomm_ring [hCid, hU]

/--
The compressed element satisfies the corner-isometry equation exactly when the
algebraic leakage term vanishes.  Calling this an isometry *on `Ran C`* requires
an external operator representation; the theorem itself is purely algebraic.
-/
theorem compression_isometry_iff_leakage_zero
    (C U : A)
    (hCstar : star C = C)
    (hCid : C * C = C)
    (hU : star U * U = 1) :
    star (compression C U) * compression C U = C ↔ leakage C U = 0 := by
  rw [compression_star_mul_self C U hCstar hCid hU]
  exact sub_eq_self

end Compression

section CharacterRecurrence

/--
A rational character with unit translation step and prescribed value at zero
is forced on every positive and negative integer.
-/
theorem integer_character_recurrence
    (b : ℤ → ℚ) (α : ℚ)
    (hstep : ∀ n : ℤ, b (n + 1) = b n + 1)
    (hzero : b 0 = α) :
    ∀ n : ℤ, b n = (n : ℚ) + α := by
  intro n
  induction n using Int.induction_on with
  | zero => simpa using hzero
  | succ n ih =>
      rw [hstep, ih]
      push_cast
      ring
  | pred n ih =>
      have h := hstep (-((n : ℤ)) - 1)
      rw [show -((n : ℤ)) - 1 + 1 = -(n : ℤ) by ring] at h
      rw [ih] at h
      push_cast at h ⊢
      linarith

/-- The recurrence stays exactly in the rational coset `1/4 + ℤ`. -/
theorem quarter_one_coset
    (b : ℤ → ℚ)
    (hstep : ∀ n : ℤ, b (n + 1) = b n + 1)
    (hzero : b 0 = 1 / 4) :
    ∀ n : ℤ, ∃ k : ℤ, b n = 1 / 4 + k := by
  intro n
  refine ⟨n, ?_⟩
  rw [integer_character_recurrence b (1 / 4) hstep hzero n]
  ring

/-- The recurrence stays exactly in the rational coset `3/4 + ℤ`. -/
theorem quarter_three_coset
    (b : ℤ → ℚ)
    (hstep : ∀ n : ℤ, b (n + 1) = b n + 1)
    (hzero : b 0 = 3 / 4) :
    ∀ n : ℤ, ∃ k : ℤ, b n = 3 / 4 + k := by
  intro n
  refine ⟨n, ?_⟩
  rw [integer_character_recurrence b (3 / 4) hstep hzero n]
  ring

end CharacterRecurrence

#print axioms compression_star_mul_self
#print axioms compression_isometry_iff_leakage_zero
#print axioms integer_character_recurrence
#print axioms quarter_one_coset
#print axioms quarter_three_coset

end TFPT.SourceReconstruction
