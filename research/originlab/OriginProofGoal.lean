import Std

namespace OriginLabProof

abbrev Block := List (Nat × Nat)
abbrev Moment := Nat → Block → Nat

-- The slots are independent oscillator modes m+1, m+2, ... .  The
-- exponents in each pair count the left and right occupations.
structure Rules (k : Nat) (f : Moment) : Prop where
  vacuum : ∀ m, f m [] = 1
  skip : ∀ m w, f m ((0, 0) :: w) = f (m+1) w
  leftEmpty : ∀ m b w, f m ((0, b+1) :: w) = 0
  rightEmpty : ∀ m a w, f m ((a+1, 0) :: w) = 0
  contract : ∀ m a b w, f m ((a+1, b+1) :: w) =
    (k*(m+1)*(a+1)) * f m ((a,b) :: w)

def energy (m : Nat) : Block → Nat
  | [] => 0
  | (_, b) :: w => (m+1)*b + energy (m+1) w

def thermal (f : Moment) (q : Rat) (m : Nat) (w : Block) : Rat :=
  (f m w : Rat) * q ^ energy m w

structure Goal (f : Moment) : Prop where
  equations : Rules 1 f
  unique : ∀ g, Rules 1 g → ∀ m w, g m w = f m w
  refinements : ∀ m w, f m (w ++ [(0,0)]) = f m w
  energyRefinements : ∀ m w, energy m (w ++ [(0,0)]) = energy m w
  orthogonal : ∀ m w, (∃ ab ∈ w, ab.1 ≠ ab.2) → f m w = 0
  positive : ∀ q : Rat, 0 ≤ q → ∀ m w, 0 ≤ thermal f q m w
  sew : ∀ (u v : Rat) m w, thermal f (u*v) m w =
    thermal f u m w * v ^ energy m w

def FreeLevelUnique : Prop :=
  ∀ (a b : Nat) (f g : Moment),
    0 < a → 0 < b → Rules a f → Rules b g → ∀ m w, f m w = g m w

end OriginLabProof
