import Mathlib

namespace TFPT.SharedHiggsMarker

theorem oppositeRightPhasesForceHiggsSquare
    (h l e n : ℂ) (hunit : star h * h = 1) (hl : l ≠ 0)
    (he : e = l * star h) (hn : n = l * h)
    (hopposite : e = -n) : h * h = -1 := by
  rw [he, hn] at hopposite
  have hstar : star h = -h := by
    apply mul_left_cancel₀ hl
    simpa [mul_neg] using hopposite
  rw [hstar, neg_mul] at hunit
  linear_combination -hunit

theorem twoStepHiggsForcesSameRightPhase
    (h l e n : ℂ) (hunit : star h * h = 1) (hh : h * h = 1)
    (he : e = l * star h) (hn : n = l * h) : e = n := by
  have hh0 : h ≠ 0 := by
    intro hz
    simp [hz] at hh
  have hstar : star h = h := by
    apply mul_right_cancel₀ hh0
    exact hunit.trans hh.symm
  rw [he, hn, hstar]

theorem neutralMajoranaEvenImpliesLeftFourthOrder
    (h l : ℂ) (hh : h * h = -1) (hn : l * h = 1) : l * l = -1 := by
  have h1 : (l * l) * (h * h) = 1 := by
    rw [show (l*l)*(h*h)=(l*h)*(l*h) by ring, hn, one_mul]
  rw [hh, mul_neg, mul_one] at h1
  linear_combination -h1

noncomputable def quaternionColumn : ℂ × ℂ → ℂ × ℂ :=
  fun x => (star x.2,-star x.1)

theorem quaternionColumnIdentity (a b c d x y : ℂ)
    (h00 : star d = a * star (a*d-b*c))
    (h01 : -star c = b * star (a*d-b*c))
    (h10 : -star b = c * star (a*d-b*c))
    (h11 : star a = d * star (a*d-b*c)) :
    quaternionColumn (a*x+b*y,c*x+d*y) =
      (star (a*d-b*c)*(a*star y-b*star x),
       star (a*d-b*c)*(c*star y-d*star x)) := by
  apply Prod.ext
  · simp only [quaternionColumn, star_add, star_mul]
    rw [h00, show star c = -(b*star (a*d-b*c)) by linear_combination -h01]
    ring
  · simp only [quaternionColumn, star_add, star_mul]
    rw [h11, show star b = -(c*star (a*d-b*c)) by linear_combination -h10]
    ring

#print axioms oppositeRightPhasesForceHiggsSquare
#print axioms twoStepHiggsForcesSameRightPhase
#print axioms neutralMajoranaEvenImpliesLeftFourthOrder
#print axioms quaternionColumnIdentity

end TFPT.SharedHiggsMarker
