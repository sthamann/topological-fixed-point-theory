import OriginProofGoal

namespace OriginLabProof

def pairing (k m : Nat) : Nat → Nat → Nat
  | 0, 0 => 1
  | _+1, 0 | 0, _+1 => 0
  | a+1, b+1 => (k*(m+1)*(a+1))*pairing k m a b

def gaussian (k m : Nat) : Block → Nat
  | [] => 1
  | (a,b) :: w => pairing k m a b * gaussian k (m+1) w

theorem gaussianRules (k : Nat) : Rules k (gaussian k) := by
  constructor
  · intro m; rfl
  · intro m w; simp [gaussian, pairing]
  · intro m b w; simp [gaussian, pairing]
  · intro m a w; simp [gaussian, pairing]
  · intro m a b w; simp [gaussian, pairing, Nat.mul_assoc]

theorem sameMoments {f g : Moment} {k : Nat} (hf : Rules k f) (hg : Rules k g) :
    ∀ m w, f m w = g m w := by
  intro m w
  induction w generalizing m with
  | nil => rw [hf.vacuum, hg.vacuum]
  | cons ab w ih =>
    obtain ⟨a,b⟩ := ab
    induction a generalizing b with
    | zero =>
      cases b with
      | zero => rw [hf.skip, hg.skip]; exact ih (m+1)
      | succ b => rw [hf.leftEmpty, hg.leftEmpty]
    | succ a step =>
      cases b with
      | zero => rw [hf.rightEmpty, hg.rightEmpty]
      | succ b => rw [hf.contract, hg.contract, step b]

theorem padGaussian (k m : Nat) (w : Block) :
    gaussian k m (w ++ [(0,0)]) = gaussian k m w := by
  induction w generalizing m with
  | nil => simp [gaussian, pairing]
  | cons ab w ih =>
    obtain ⟨a,b⟩ := ab
    simp [gaussian, ih]

theorem padEnergy (m : Nat) (w : Block) :
    energy m (w ++ [(0,0)]) = energy m w := by
  induction w generalizing m with
  | nil => simp [energy]
  | cons ab w ih =>
    obtain ⟨a,b⟩ := ab
    simp [energy, ih]

theorem pairingOff (k m a b : Nat) (h : a ≠ b) : pairing k m a b = 0 := by
  induction a generalizing b with
  | zero => cases b <;> simp_all [pairing]
  | succ a ih =>
    cases b with
    | zero => rfl
    | succ b =>
      have hab : a ≠ b := fun eq => h (congrArg Nat.succ eq)
      simp [pairing, ih b hab]

theorem gaussianOff (k m : Nat) (w : Block) (h : ∃ ab ∈ w, ab.1 ≠ ab.2) :
    gaussian k m w = 0 := by
  induction w generalizing m with
  | nil => simp at h
  | cons ab w ih =>
    obtain ⟨a,b⟩ := ab
    obtain ⟨z, hz, hne⟩ := h
    rcases List.mem_cons.mp hz with eq | tail
    · subst z
      simp [gaussian, pairingOff k m a b hne]
    · simp [gaussian, ih (m+1) ⟨z, tail, hne⟩]

theorem ratMulPow (u v : Rat) (n : Nat) : (u*v)^n = u^n*v^n := by
  induction n with
  | zero => simp [Rat.pow_zero]
  | succ n ih =>
    rw [Rat.pow_succ, Rat.pow_succ, Rat.pow_succ, ih]
    calc
      (u^n*v^n)*(u*v) = u^n*(v^n*(u*v)) := Rat.mul_assoc ..
      _ = u^n*((v^n*u)*v) := by rw [Rat.mul_assoc (v^n)]
      _ = u^n*((u*v^n)*v) := by rw [Rat.mul_comm (v^n) u]
      _ = (u^n*u)*(v^n*v) := by rw [Rat.mul_assoc (u^n), Rat.mul_assoc u]

theorem closeGaussian (k : Nat) (h : k = 1) : Goal (gaussian k) := by
  subst k
  refine ⟨gaussianRules 1, ?_, padGaussian 1, padEnergy, gaussianOff 1, ?_, ?_⟩
  · intro g hg
    exact sameMoments hg (gaussianRules 1)
  · intro q hq m w
    exact Rat.mul_nonneg Rat.natCast_nonneg (Rat.pow_nonneg hq)
  · intro u v m w
    simp only [thermal, ratMulPow, Rat.mul_assoc]

theorem refuteFreeLevel : ¬ FreeLevelUnique := by
  intro h
  have equal := h 1 2 (gaussian 1) (gaussian 2) (by decide) (by decide)
    (gaussianRules 1) (gaussianRules 2) 0 [(1,1)]
  exact (by decide : (1 : Nat) ≠ 2) equal

end OriginLabProof
