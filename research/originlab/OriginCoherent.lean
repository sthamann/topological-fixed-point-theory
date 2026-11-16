import Std

namespace OriginCoherent

/-- The split of |n+1,0> into |n,1> preserves total energy. -/
theorem affineGenerator (f : Nat → Int)
    (split : ∀ n, f n + f 1 = f (n + 1) + f 0) (n : Nat) :
    f n = f 0 + (n : Int) * (f 1 - f 0) := by
  induction n with
  | zero => simp
  | succ n ih =>
    have step := split n
    rw [ih] at step
    have expand : ((n + 1 : Nat) : Int) * (f 1 - f 0)
        = (n : Int) * (f 1 - f 0) + (f 1 - f 0) := by
      simp [Int.add_mul]
    rw [expand]
    omega

/-- A fourth-root phase is stored as its exponent modulo four.
    Each even/odd pair of the oscillator carries this action. -/
abbrev PairState := Bool × Fin 4

def swap (state : PairState) : PairState := (!state.1, state.2)
def marked (state : PairState) : PairState :=
  (!state.1, state.2 + (if state.1 then 3 else 1))

def ordered (state : PairState) : PairState := (state.1, state.2 + (if state.1 then 3 else 1))
def minus (state : PairState) : PairState := (state.1, state.2 + 2)

theorem markerIdentities : ∀ s : PairState,
    swap (swap s) = s ∧ marked (marked s) = s ∧
    swap (marked s) = ordered s ∧
    swap (marked (swap (marked s))) = minus s := by
  intro ⟨odd, phase⟩
  have choices : phase.val = 0 ∨ phase.val = 1 ∨ phase.val = 2 ∨ phase.val = 3 := by
    omega
  rcases choices with h | h | h | h
  · have e : phase = 0 := Fin.ext h
    subst phase
    cases odd <;> decide
  · have e : phase = 1 := Fin.ext h
    subst phase
    cases odd <;> decide
  · have e : phase = 2 := Fin.ext h
    subst phase
    cases odd <;> decide
  · have e : phase = 3 := Fin.ext h
    subst phase
    cases odd <;> decide

/-- Equal diagonal weights on each mixer sector force a geometric sequence
    once the vacuum mass is nonzero. -/
theorem geometricWeights (p : Nat → Rat) (vacuum : p 0 ≠ 0)
    (split : ∀ n, p n * p 1 = p (n+1) * p 0) (n : Nat) :
    p n = p 0 * (p 1 / p 0)^n := by
  have step (k : Nat) : p (k+1) = p k * (p 1 / p 0) := by
    calc
      p (k+1) = p (k+1) * p 0 / p 0 := (Rat.mul_div_cancel vacuum).symm
      _ = p k * p 1 / p 0 := by rw [split k]
      _ = p k * (p 1 / p 0) := by simp only [Rat.div_def, Rat.mul_assoc]
  induction n with
  | zero => simp
  | succ n ih => rw [step n, ih, Rat.pow_succ, Rat.mul_assoc]

end OriginCoherent
