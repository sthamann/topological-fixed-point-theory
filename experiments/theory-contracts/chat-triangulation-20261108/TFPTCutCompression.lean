import Mathlib

namespace TFPT.CutCompression

theorem projectionCompressionDefect {R : Type*} [Ring R]
    (B E : R) (hB : B*B = B) (hE : E*E = E) :
    B*E*B - (B*E*B)*(B*E*B) = B*E*(1-B)*E*B := by
  have h₁ : (B*E*B)*(B*E*B) = B*E*B*E*B := by
    calc
      _ = B*E*(B*B)*E*B := by noncomm_ring
      _ = _ := by rw [hB]
  have h₂ : B*E*E*B = B*E*B := by
    calc
      _ = B*(E*E)*B := by noncomm_ring
      _ = _ := by rw [hE]
  rw [h₁]
  calc
    _ = B*E*E*B - B*E*B*E*B :=
      congrArg (fun x => x - B*E*B*E*B) h₂.symm
    _ = _ := by noncomm_ring

theorem compressedInvolutionDefect {R : Type*} [Ring R] (E : R) :
    1 - (1 - 2*E)*(1 - 2*E) = 4*(E-E*E) := by
  noncomm_ring

#print axioms projectionCompressionDefect
#print axioms compressedInvolutionDefect

end TFPT.CutCompression
