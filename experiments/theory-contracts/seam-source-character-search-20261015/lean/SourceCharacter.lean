import Mathlib

/-!
# Exact source-character identities

This file proves only rational vector identities and a lattice-saturation
statement for the displayed two-dimensional source plane.  It assumes no
physical origin, state, locality, or source-selection principle.
The identity for the compensated vector does not remove the `α` dependence
of the unmodified source character, time, or twist.
-/

namespace TFPT.SourceCharacter

abbrev V8 := Fin 8 → ℚ

/-- The source character vector at parameter `α`. -/
def d (α : ℚ) : V8 :=
  ![-1 / 2, 0, α / 4, -1 / 2, 0, (α - 1) / 4, -1 / 4, -1 / 4]

/-- The distinguished integral root `e₃ + e₆` (one-based labels). -/
def r : V8 := ![0, 0, 1, 0, 0, 1, 0, 0]

def q : V8 := 2 • r

/-- The α-independent compensated character. -/
def t : V8 := ![-1, 0, 0, -1, 0, -1 / 2, -1 / 2, -1 / 2]

/-- The fixed signed permutation from the source calculation. -/
def O (v : V8) : V8 :=
  ![v 3, -v 4, -v 5, v 0, -v 1, -v 2, v 7, v 6]

/-- This compensated combination alone is `α`-independent, coordinate by coordinate. -/
theorem compensated_character (α : ℚ) :
    2 • d α - (α / 2) • r = t := by
  funext i
  fin_cases i <;> simp [d, r, t] <;> ring

/-- The signed permutation exchanges the two source-character branches. -/
theorem O_d (α : ℚ) : O (d α) = d (1 - α) := by
  funext i
  fin_cases i <;> simp [O, d] <;> ring

/-- The distinguished doubled root is odd under `O`. -/
theorem O_q : O q = -q := by
  funext i
  fin_cases i <;> simp [O, q, r]

/-- The apparent modulo-lattice relation is in fact an exact equality. -/
theorem O_t : O t = t + (1 / 4 : ℚ) • q := by
  funext i
  fin_cases i <;> simp [O, t, q, r] <;> norm_num

/-- Explicit coordinate sum, used to avoid hiding the parity condition. -/
def coordSum (v : V8) : ℚ :=
  v 0 + v 1 + v 2 + v 3 + v 4 + v 5 + v 6 + v 7

/-- The integral, even-sum component `D₈` of the rational `E₈` model. -/
def IntegerEven (v : V8) : Prop :=
  (∀ i, ∃ z : ℤ, v i = z) ∧ ∃ k : ℤ, coordSum v = 2 * k

/-- The half-integral, even-sum component of the rational `E₈` model. -/
def HalfEven (v : V8) : Prop :=
  (∀ i, ∃ z : ℤ, v i = z + 1 / 2) ∧ ∃ k : ℤ, coordSum v = 2 * k

/-- The standard rational-coordinate model `E₈ = D₈ ∪ (D₈ + (1/2)⁸)`. -/
def E8Lattice (v : V8) : Prop := IntegerEven v ∨ HalfEven v

/-- The source plane through `d(1/4)` and `r`. -/
def sourcePlane (a b : ℚ) : V8 := a • d (1 / 4) + b • r

/-- The explicit integral point produced by parameters `s,l`. -/
def latticePoint (s l : ℤ) : V8 :=
  ![-4 * (s : ℚ), 0, (l : ℚ), -4 * (s : ℚ), 0,
    (l : ℚ) - 2 * (s : ℚ), -2 * (s : ℚ), -2 * (s : ℚ)]

theorem sourcePlane_param (s l : ℤ) :
    sourcePlane (8 * (s : ℚ)) ((l : ℚ) - (s : ℚ) / 2) = latticePoint s l := by
  funext i
  fin_cases i <;> simp [sourcePlane, latticePoint, d, r] <;> ring

/--
Because coordinates 2 and 5 (one-based) vanish on the source plane, an `E₈`
half-coset point cannot occur there.  Thus the plane's `E₈` intersection is
its integral-even component, parametrized exactly as below.
-/
theorem source_plane_saturation (a b : ℚ) :
    IntegerEven (sourcePlane a b) ↔
      ∃ s l : ℤ, a = 8 * s ∧ b = l - (s : ℚ) / 2 := by
  constructor
  · rintro ⟨hint, k, hsum⟩
    obtain ⟨z3, hz3⟩ := hint 2
    obtain ⟨z7, hz7⟩ := hint 6
    refine ⟨z3 + 3 * z7 - k, z3, ?_, ?_⟩
    · simp [sourcePlane, d, r] at hz3 hz7
      simp [coordSum, sourcePlane, d, r] at hsum
      push_cast
      linarith
    · simp [sourcePlane, d, r] at hz3 hz7
      simp [coordSum, sourcePlane, d, r] at hsum
      push_cast
      linarith
  · rintro ⟨s, l, rfl, rfl⟩
    rw [sourcePlane_param]
    constructor
    · intro i
      fin_cases i
      · exact ⟨-4 * s, by simp [latticePoint]⟩
      · exact ⟨0, by simp [latticePoint]⟩
      · exact ⟨l, by simp [latticePoint]⟩
      · exact ⟨-4 * s, by simp [latticePoint]⟩
      · exact ⟨0, by simp [latticePoint]⟩
      · exact ⟨l - 2 * s, by simp [latticePoint]⟩
      · exact ⟨-2 * s, by simp [latticePoint]⟩
      · exact ⟨-2 * s, by simp [latticePoint]⟩
    · refine ⟨l - 7 * s, ?_⟩
      push_cast
      simp [coordSum, latticePoint]
      ring

/-- The two identically zero source-plane coordinates exclude the half coset. -/
theorem sourcePlane_not_half (a b : ℚ) : ¬ HalfEven (sourcePlane a b) := by
  rintro ⟨hhalf, _⟩
  obtain ⟨z, hz⟩ := hhalf 1
  simp [sourcePlane, d, r] at hz
  have hzq : (2 : ℚ) * (z : ℚ) = -1 := by linarith
  have hzi : (2 : ℤ) * z = -1 := by exact_mod_cast hzq
  omega

/-- Full `E₈` saturation of the source plane, including exclusion of its half coset. -/
theorem source_plane_E8_saturation (a b : ℚ) :
    E8Lattice (sourcePlane a b) ↔
      ∃ s l : ℤ, a = 8 * s ∧ b = l - (s : ℚ) / 2 := by
  constructor
  · rintro (hint | hhalf)
    · exact (source_plane_saturation a b).mp hint
    · exact False.elim (sourcePlane_not_half a b hhalf)
  · intro h
    exact Or.inl ((source_plane_saturation a b).mpr h)

#print axioms compensated_character
#print axioms O_d
#print axioms O_q
#print axioms O_t
#print axioms source_plane_saturation
#print axioms sourcePlane_not_half
#print axioms source_plane_E8_saturation

end TFPT.SourceCharacter
