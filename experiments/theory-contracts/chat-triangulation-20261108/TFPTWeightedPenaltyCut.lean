import Mathlib

namespace TFPT.WeightedPenaltyCut

theorem gapRecovered (retained low x : ℚ) (hx : x ≠ 1) :
    (retained + x * low - x * (retained + low)) / (1 - x) = retained := by
  have hx' : 1 - x ≠ 0 := sub_ne_zero.mpr (Ne.symm hx)
  field_simp [hx']
  ring

theorem thetaRecovered (even odd x : ℚ) (hx : x ≠ 1) :
    (2 * (even + x * odd) - (1 + x) * (even + odd)) / (1 - x) = even - odd := by
  have hx' : 1 - x ≠ 0 := sub_ne_zero.mpr (Ne.symm hx)
  field_simp [hx']
  ring

theorem oddWeightRecovered (even odd x : ℚ) (hx : x ≠ 1) :
    ((even + odd) - (even + x * odd)) / (1 - x) = odd := by
  have hx' : 1 - x ≠ 0 := sub_ne_zero.mpr (Ne.symm hx)
  field_simp [hx']
  ring

#print axioms gapRecovered
#print axioms thetaRecovered
#print axioms oddWeightRecovered

end TFPT.WeightedPenaltyCut
