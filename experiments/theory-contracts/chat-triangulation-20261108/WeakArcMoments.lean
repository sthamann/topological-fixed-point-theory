import Mathlib

namespace TFPT.WeakArcMoments

theorem twoHistoryLoss
    (c₁ c₂ b₁ b₂ : ℂ) (w : ℝ)
    (h₁ : Complex.normSq c₁ + Complex.normSq b₁ = 1)
    (h₂ : Complex.normSq c₂ + Complex.normSq b₂ = 1) :
    1 - Complex.normSq ((w : ℂ) * c₁ + ((1 - w : ℝ) : ℂ) * c₂) =
      w * Complex.normSq b₁ + (1 - w) * Complex.normSq b₂ +
        w * (1 - w) * Complex.normSq (c₁ - c₂) := by
  rw [show Complex.normSq b₁ = 1 - Complex.normSq c₁ by linarith,
      show Complex.normSq b₂ = 1 - Complex.normSq c₂ by linarith]
  simp only [Complex.normSq_apply, Complex.add_re, Complex.add_im, Complex.mul_re,
    Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im, Complex.sub_re, Complex.sub_im]
  ring

theorem sameMeanConjugateButSwitchAndVarianceDiffer
    (r s : ℝ) (hu : r^2 + s^2 = 1) :
    let transverse₁ : ℂ := r
    let transverse₂ : ℂ := r
    let longitudinal₁ : ℂ := ⟨r,s⟩
    let longitudinal₂ : ℂ := ⟨r,-s⟩
    transverse₁ + transverse₂ = longitudinal₁ + longitudinal₂ ∧
      Complex.normSq (transverse₁ - transverse₂) = 0 ∧
      Complex.normSq (longitudinal₁ - longitudinal₂) / 4 = s^2 ∧
      1 - Complex.normSq transverse₁ = s^2 ∧
      1 - Complex.normSq longitudinal₁ = 0 := by
  dsimp only
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · apply Complex.ext <;> simp
  · simp
  all_goals simp only [Complex.normSq_apply, Complex.ofReal_re, Complex.ofReal_im,
    Complex.sub_re, Complex.sub_im]
  · ring
  · nlinarith
  · nlinarith

#print axioms twoHistoryLoss
#print axioms sameMeanConjugateButSwitchAndVarianceDiffer

end TFPT.WeakArcMoments
