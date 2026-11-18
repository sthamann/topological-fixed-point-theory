import Mathlib

/-!
# Algebraic core of a two-face product-collar reconstruction

The certificate concerns an exact `2 × 2` rational graph projection and an
abstract anticommuting-involution identity.  The product-collar datum
`D = J(∂ᵣ + B)` is input context, not a conclusion of these theorems.
-/

namespace TFPT.FullCollar

section GraphProjection

abbrev V2 := Fin 2 → ℚ
abbrev M2 := Matrix (Fin 2) (Fin 2) ℚ

/-- The orthogonal graph projection associated to the rational scalar `t`. -/
def G (t : ℚ) : M2 :=
  !![1 / (1 + t ^ 2), t / (1 + t ^ 2);
     t / (1 + t ^ 2), t ^ 2 / (1 + t ^ 2)]

/-- A vector on the graph of multiplication by `t`. -/
def graphVec (t x : ℚ) : V2 := ![x, t * x]

lemma denominator_ne_zero (t : ℚ) : 1 + t ^ 2 ≠ 0 := by positivity

/-- `G(t)` is symmetric. -/
theorem G_transpose (t : ℚ) : Matrix.transpose (G t) = G t := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [G]

/-- `G(t)` is idempotent. -/
theorem G_idempotent (t : ℚ) : G t * G t = G t := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [G, Matrix.mul_apply, Fin.sum_univ_two] <;>
    field_simp [denominator_ne_zero t]

/-- Every graph vector is fixed by `G(t)`. -/
theorem G_fixes_graph (t x : ℚ) :
    Matrix.mulVec (G t) (graphVec t x) = graphVec t x := by
  funext i
  fin_cases i <;>
    simp [G, graphVec, Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;>
    field_simp [denominator_ne_zero t]

/-- Every projected vector lies explicitly on the graph of `t`. -/
theorem G_range_formula (t : ℚ) (v : V2) :
    Matrix.mulVec (G t) v = graphVec t ((v 0 + t * v 1) / (1 + t ^ 2)) := by
  funext i
  fin_cases i <;>
    simp [G, graphVec, Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;>
    field_simp [denominator_ne_zero t]

/-- The upper-left entry needed for reconstruction is nonzero. -/
theorem G₁₁_ne_zero (t : ℚ) : G t 0 0 ≠ 0 := by
  simp [G, denominator_ne_zero t]

/-- The full two-face projection recovers its graph parameter exactly. -/
theorem recover_t (t : ℚ) : G t 1 0 / G t 0 0 = t := by
  simp [G]
  field_simp [denominator_ne_zero t]

/-- Consequently the full graph-projection datum is injective. -/
theorem G_injective {t s : ℚ} (h : G t = G s) : t = s := by
  have h10 : G t 1 0 = G s 1 0 := congrFun (congrFun h 1) 0
  have h00 : G t 0 0 = G s 0 0 := congrFun (congrFun h 0) 0
  calc
    t = G t 1 0 / G t 0 0 := (recover_t t).symm
    _ = G s 1 0 / G s 0 0 := by rw [h10, h00]
    _ = s := recover_t s

/-- Diagonal-only projection data cannot distinguish `t` from `-t`. -/
def diagonalData (t : ℚ) : V2 := ![G t 0 0, G t 1 1]

theorem diagonalData_sign_blind (t : ℚ) : diagonalData t = diagonalData (-t) := by
  funext i
  fin_cases i <;> simp [diagonalData, G]

/-- The complete projection distinguishes the two signs unless `t=0`. -/
theorem G_eq_neg_iff (t : ℚ) : G t = G (-t) ↔ t = 0 := by
  constructor
  · intro h
    have ht : t = -t := G_injective h
    linarith
  · rintro rfl
    norm_num

end GraphProjection

section RetainedPolarization

variable {A : Type*} [Ring A]

/--
If the carrier projector `Q` is already retained by the positive compression
`C`, then the induced polarization `2C-1` is the identity on that projector.
This is an obstruction to reading a new `3+2` marking from positive-only
compression; `C` and `Q` are inputs here, not selected by the theorem.
-/
theorem retained_polarization_is_identity
    (C Q : A) (hCQ : C * Q = Q) :
    (2 * C - 1) * Q = Q := by
  noncomm_ring [hCQ]

/-- The complementary, negative compression vanishes on the retained projector. -/
theorem negative_restriction_zero
    (C Q : A) (hCQ : C * Q = Q) :
    (1 - C) * Q = 0 := by
  noncomm_ring [hCQ]

/-- Element form: `Q*v=v` makes the retained polarization fix `v`. -/
theorem retained_polarization_fixes_range
    (C Q v : A) (hCQ : C * Q = Q) (hQv : Q * v = v) :
    (2 * C - 1) * v = v := by
  have hCv : C * v = v := by
    calc
      C * v = C * (Q * v) := by rw [hQv]
      _ = (C * Q) * v := by rw [mul_assoc]
      _ = Q * v := by rw [hCQ]
      _ = v := hQv
  noncomm_ring [hCv]

end RetainedPolarization

section AnticommutingInvolutions

variable {A : Type*} [Ring A]

/-- Anticommuting involutions generate a square root of `-1`. -/
theorem anticommuting_involutions_square
    (R F : A)
    (hR : R * R = 1)
    (hF : F * F = 1)
    (hanti : R * F = -(F * R)) :
    (R * F) * (R * F) = -1 := by
  have hanti' : F * R = -(R * F) := by noncomm_ring [hanti]
  calc
    (R * F) * (R * F) = R * (F * R) * F := by noncomm_ring
    _ = R * (-(R * F)) * F := by rw [hanti']
    _ = -((R * R) * (F * F)) := by noncomm_ring
    _ = -1 := by rw [hR, hF, one_mul]

/-- In characteristic zero the two anticommuting involutions cannot coincide. -/
theorem anticommuting_involutions_ne
    [CharZero A]
    (R F : A)
    (hR : R * R = 1)
    (_hF : F * F = 1)
    (hanti : R * F = -(F * R)) :
    R ≠ F := by
  intro heq
  subst F
  have hone : (1 : A) = -1 := by
    calc
      1 = R * R := hR.symm
      _ = -(R * R) := hanti
      _ = -1 := by rw [hR]
  have htwo : (2 : A) = 0 := by
    have hadd : (1 : A) + 1 = -1 + 1 := congrArg (fun x : A => x + 1) hone
    calc
      (2 : A) = 1 + 1 := by norm_num
      _ = -1 + 1 := hadd
      _ = 0 := by simp
  norm_num at htwo

end AnticommutingInvolutions

#print axioms G_transpose
#print axioms G_idempotent
#print axioms G_fixes_graph
#print axioms G_range_formula
#print axioms G₁₁_ne_zero
#print axioms recover_t
#print axioms G_injective
#print axioms diagonalData_sign_blind
#print axioms G_eq_neg_iff
#print axioms retained_polarization_is_identity
#print axioms negative_restriction_zero
#print axioms retained_polarization_fixes_range
#print axioms anticommuting_involutions_square
#print axioms anticommuting_involutions_ne

end TFPT.FullCollar
