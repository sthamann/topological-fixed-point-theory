import TfptCarrier.PositiveDescentMaster
import TfptCarrier.GramCompactness
import TfptCarrier.BoundaryYukawaKernelInterface
import TfptCarrier.SeamScalingLimit
import TfptCarrier.BWKeystone
import Mathlib.Analysis.SpecialFunctions.Exp

#print TfptCarrier.PositiveDescentMaster.master_core
#print axioms TfptCarrier.PositiveDescentMaster.master_core
#print TfptCarrier.PositiveDescentMaster.tail_bound
#print axioms TfptCarrier.PositiveDescentMaster.tail_bound
#print axioms TfptCarrier.GramCompactness.exists_subseq_posSemidef_limit
#print axioms TFPT.Carrier.BoundaryYukawaKernelInterface.BoundaryYukawaKernel.nonempty_iff_finrank_eq_three
#print axioms TfptCarrier.SeamScalingLimit.seamScalingLimit
#print axioms TfptCarrier.BWKeystone.qgeoSymIsCorollary

namespace OriginalRadialTilt

/-- Here `a = sqrt(2)*v` in the existing R048 radial coordinates.
This is the ORIGINAL cubic plus quartic, not its perturbative series. -/
noncomputable def Vint (lam a x : ℝ) : ℝ := lam * x^4/4 + lam*a*x^3

theorem exact_floor_identity (lam a x : ℝ) :
    Vint lam a x + 27*lam*a^4/4 =
      lam/4 * (x+3*a)^2 * ((x-a)^2+2*a^2) := by
  unfold Vint
  ring

theorem exact_floor (lam a x : ℝ) (hlam : 0 ≤ lam) :
    -27*lam*a^4/4 ≤ Vint lam a x := by
  have hp : 0 ≤ lam/4 * (x+3*a)^2 * ((x-a)^2+2*a^2) := by positivity
  rw [← exact_floor_identity] at hp
  linarith

theorem bounded_original_tilt (lam a x : ℝ) (hlam : 0 ≤ lam) :
    Real.exp (-Vint lam a x) ≤ Real.exp (27*lam*a^4/4) := by
  apply Real.exp_le_exp.mpr
  linarith [exact_floor lam a x hlam]

theorem floor_attained (lam a : ℝ) :
    Vint lam a (-3*a) = -27*lam*a^4/4 := by
  unfold Vint
  ring

#print axioms exact_floor_identity
#print axioms exact_floor
#print axioms bounded_original_tilt
#print axioms floor_attained

end OriginalRadialTilt
