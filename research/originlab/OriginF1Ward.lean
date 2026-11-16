import Std

namespace OriginF1Ward

def coefficient (n : Nat) : Rat := ((2*n+1 : Nat) : Rat)
def spectrum (a : Rat) (n : Nat) : Rat := coefficient n * a

def Ward (f : Nat → Rat) : Prop :=
  ∀ n, coefficient n * f (n+1) = coefficient (n+1) * f n

theorem familyWard (a : Rat) : Ward (spectrum a) := by
  intro n
  simp only [spectrum, ← Rat.mul_assoc]
  rw [Rat.mul_comm (coefficient n) (coefficient (n+1))]

theorem wardDeterminesRatios (f : Nat → Rat) (ward : Ward f) :
    ∀ n, f n = spectrum (f 0) n := by
  intro n
  induction n with
  | zero =>
    have one : coefficient 0 = 1 := by decide
    simp [spectrum, one, Rat.one_mul]
  | succ n ih =>
    have hs : coefficient n * f (n+1) = coefficient n * spectrum (f 0) (n+1) := by
      rw [ward n, ih]
      simp only [spectrum, ← Rat.mul_assoc]
      rw [Rat.mul_comm (coefficient n) (coefficient (n+1))]
    have hn : coefficient n ≠ 0 := by
      have h : 0 < coefficient n := Rat.natCast_pos.mpr (by omega)
      intro eq
      rw [eq] at h
      exact (by decide : ¬ (0 : Rat) < 0) h
    have hc := congrArg (fun x : Rat => (coefficient n)⁻¹ * x) hs
    simpa only [← Rat.mul_assoc, Rat.inv_mul_cancel _ hn, Rat.one_mul] using hc

theorem wardLeavesScale : ∃ f g : Nat → Rat, Ward f ∧ Ward g ∧ f 0 ≠ g 0 :=
  ⟨spectrum 1, spectrum 2, familyWard 1, familyWard 2, by
    have one : coefficient 0 = 1 := by decide
    simpa [spectrum, one, Rat.one_mul] using (by decide : (1 : Rat) ≠ 2)⟩

end OriginF1Ward
