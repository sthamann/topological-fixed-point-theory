import Mathlib
import TfptCarrier.BoundaryYukawaKernelInterface
import TfptCarrier.PositiveDescentMaster
import TfptCarrier.BWKeystone

namespace TFPT.DeterminantHiggsLift

def diagonalTest : Matrix (Fin 2) (Fin 2) ℂ := !![1, 0; 0, -1]
def swapTest : Matrix (Fin 2) (Fin 2) ℂ := !![0, 1; 1, 0]

theorem centralInvolutionDetOne (q : Matrix (Fin 2) (Fin 2) ℂ)
    (hdiag : q * diagonalTest = diagonalTest * q)
    (hswap : q * swapTest = swapTest * q)
    (hsquare : q * q = 1) : q.det = 1 := by
  have h01 := congrArg (fun m : Matrix (Fin 2) (Fin 2) ℂ => m 0 1) hdiag
  have h10 := congrArg (fun m : Matrix (Fin 2) (Fin 2) ℂ => m 1 0) hdiag
  have h00 := congrArg (fun m : Matrix (Fin 2) (Fin 2) ℂ => m 0 1) hswap
  have hsq := congrArg (fun m : Matrix (Fin 2) (Fin 2) ℂ => m 0 0) hsquare
  simp [diagonalTest, Matrix.mul_apply, Fin.sum_univ_two] at h01 h10
  simp [swapTest, Matrix.mul_apply, Fin.sum_univ_two] at h00
  simp [Matrix.mul_apply, Fin.sum_univ_two] at hsq
  have hb : q 0 1 = 0 := by linear_combination -(1 / 2 : ℂ) * h01
  have hc : q 1 0 = 0 := by linear_combination (1 / 2 : ℂ) * h10
  rw [Matrix.det_fin_two, hb, hc, ← h00]
  simpa [hb, hc] using hsq

theorem signLiftSquare (z : ℂ) (hdet : z * z = -1) :
    (z • (1 : Matrix (Fin 2) (Fin 2) ℂ)) *
      (z • (1 : Matrix (Fin 2) (Fin 2) ℂ)) = -1 := by
  rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul, one_mul, hdet, neg_one_smul]

theorem positiveSectorFraction (pos neg : ℝ) (hn : pos + neg ≠ 0) :
    (1 + (pos - neg) / (pos + neg)) / 2 = pos / (pos + neg) := by
  field_simp; ring

#print axioms centralInvolutionDetOne
#print axioms signLiftSquare
#print axioms positiveSectorFraction

#check TFPT.Carrier.BoundaryYukawaKernelInterface.BoundaryYukawaKernel.finrank_E_eq_three
#check TfptCarrier.PositiveDescentMaster.tail_bound
#check TfptCarrier.BWKeystone.qgeoFromSeam

end TFPT.DeterminantHiggsLift
