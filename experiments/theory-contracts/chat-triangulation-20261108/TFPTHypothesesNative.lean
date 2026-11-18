import Mathlib

namespace TFPT.HypothesesNative

/-- The central lift preserves the lepton phase and removes color order three. -/
theorem centralLiftArithmetic :
    (9 : ZMod 4) = 1 ∧ (9 : ZMod 3) = 0 ∧
    (18 : ZMod 12) = 6 ∧ (9 : ZMod 12) = 1 + 2 * 4 := by
  decide

theorem centralParentZeros (s : ℤ) (lo : -18 ≤ s) (hi : s ≤ 18) :
    0 ≤ s^2 * (s^2 - 144)^2 ∧
    (s^2 * (s^2 - 144)^2 = 0 ↔ s % 12 = 0) := by
  interval_cases s <;> norm_num

theorem colorParentZeros (t : ℤ) (lo : -6 ≤ t) (hi : t ≤ 6) :
    0 ≤ t^2 * (t^2 - 9)^2 * (36 - t^2) ∧
    (t^2 * (t^2 - 9)^2 * (36 - t^2) = 0 ↔ t % 3 = 0) := by
  interval_cases t <;> norm_num

#print axioms centralLiftArithmetic
#print axioms centralParentZeros
#print axioms colorParentZeros

end TFPT.HypothesesNative
