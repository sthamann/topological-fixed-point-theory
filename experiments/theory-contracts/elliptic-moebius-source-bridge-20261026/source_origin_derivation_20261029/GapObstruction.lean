import Mathlib.Tactic

/- A source-independent spectral-set lemma. The actual quarter source supplies
   the positive spectral gap separately. This does not formalize all spectral
   theory or the physical P1 source identification. No new axiom is declared. -/
namespace TfptSourceOrigin

theorem positive_gap_blocks_half_scaling
    (s : Set ℝ) (gap : ℝ) (hpos : 0 < gap) (hmem : gap ∈ s)
    (hgap : ∀ x ∈ s, 0 < x → gap ≤ x) :
    ¬ (∀ x ∈ s, x / 2 ∈ s) := by
  grind

#print axioms positive_gap_blocks_half_scaling

end TfptSourceOrigin
