import Std

-- Once the discriminant calculation gives n ∣ 4, the family range
-- leaves a unique rank parameter. The discriminant calculation remains
-- a separate mathematical premise of this lemma.
theorem originlab_rank_from_divisibility (n : Nat) (hn : 3 ≤ n) (hd : n ∣ 4) :
    n = 4 := by
  have hu : n ≤ 4 := Nat.le_of_dvd (by decide) hd
  have h3 : ¬ 3 ∣ (4 : Nat) := by decide
  by_cases h : n = 3
  · exact False.elim (h3 (h ▸ hd))
  · omega

#print axioms originlab_rank_from_divisibility
