/-
  PositiveDescentMaster — the kernel-checked abstract core of
  TFPT.POSITIVE_DESCENT.MASTER.01.
  ================================================================

  THE MASTER THEOREM CANDIDATE (frozen in the sibling probe
  `experiments/tfpt-discovery/positive_descent_master_probe.py`):
  finite operator systems with unital inclusions (1), positive
  states (2), CP transitions (3), a finite symmetry group with
  character sectors (4), SUMMABLE state defects (5), a
  sector-adapted closure datum (6), and a Mosco form family (7)
  yield a unique positive limit state whose compatible sector
  projections have well-defined positive limits.

  WHAT THIS MODULE PROVES (kernel-checked, no `sorry`).  The
  load-bearing core in the finite sector-weight form: states
  resolved over the character sectors of the finite abelian
  symmetry are vectors of nonnegative sector weights in
  `Fin m → ℝ` — exactly the data BOTH instances carry (G_net:
  m = 4 μ4-sectors; prime/packet: m = 128 characters of
  C2 × F2^4 × Z4, aggregated into 12 classes).  Then:

    (a) `exists_limit` / `limit_unique` / `tail_bound` —
        hypothesis (5) telescoping: summable defects make the
        state sequence Cauchy, a unique limit exists (finite
        products of ℝ are complete), and the quantitative Cauchy
        tail `dist (p n) q ≤ ∑' k, δ (n + k)` holds on every
        rung.  (Mathlib: `cauchySeq_of_dist_le_of_summable`,
        `dist_le_tsum_of_dist_le_of_tendsto`.)

    (b) `limit_nonneg` — hypothesis (2) is closed under the
        limit: coordinate evaluation is continuous, so
        nonnegative sector weights stay nonnegative in the limit.
        The MATRIX-state version (density matrices, PSD closed
        under entrywise limits) is the GramCompactness.lean
        legacy (`posSemidef_of_tendsto`), cited not re-proved.

    (c) `sector_tendsto` / `master_core` — hypothesis (4)
        commutes with the limit: sector aggregation (any finite
        set of characters) is a finite sum of continuous
        evaluations, so sector weights converge to the sector
        weights of the limit, and each limit sector is
        nonnegative.  `mu4_char_sum` records the μ4 orthogonality
        head Σ_{j<4} i^j = 0 (the exactness token for the
        character projectors used by both instances).

  WHAT IS TYPED, NOT PROVED HERE (honest split).  The
  noncommutative upgrade (states as functionals on C*-algebras,
  GNS, normality), the CP structure of the transitions (3) — the
  probe verifies it numerically via Choi matrices / convex-update
  identities — and the closure data (6), (7).  Hypotheses (6)/(7)
  are NOT consumed by (a)–(c); they feed the later per-instance
  analytic steps.  No continuum theorem and no RH statement is
  formalized or claimed.
-/
import Mathlib.Topology.Algebra.InfiniteSum.Real
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic

namespace TfptCarrier
namespace PositiveDescentMaster

open Filter Topology

variable {m : ℕ}

/-- (a) Hypothesis (5) telescoping: if the sector-weight vectors of
consecutive states differ by a summable defect sequence, a limit
state exists (Banach completeness of the finite sector-weight
space; elementary, exactly as typed in the G_net martingale run). -/
theorem exists_limit (p : ℕ → (Fin m → ℝ)) (δ : ℕ → ℝ)
    (hδ : ∀ n, dist (p n) (p (n + 1)) ≤ δ n) (hs : Summable δ) :
    ∃ q : Fin m → ℝ, Tendsto p atTop (𝓝 q) :=
  cauchySeq_tendsto_of_complete (cauchySeq_of_dist_le_of_summable δ hδ hs)

/-- (a) Uniqueness of the limit state (Hausdorff limits). -/
theorem limit_unique {p : ℕ → (Fin m → ℝ)} {q q' : Fin m → ℝ}
    (h : Tendsto p atTop (𝓝 q)) (h' : Tendsto p atTop (𝓝 q')) :
    q = q' :=
  tendsto_nhds_unique h h'

/-- (a) The quantitative Cauchy tail on every rung: the distance to
the limit is bounded by the summed tail defects — the exact bound
the probes gate numerically. -/
theorem tail_bound (p : ℕ → (Fin m → ℝ)) (δ : ℕ → ℝ) {q : Fin m → ℝ}
    (hδ : ∀ n, dist (p n) (p (n + 1)) ≤ δ n) (hs : Summable δ)
    (hq : Tendsto p atTop (𝓝 q)) (n : ℕ) :
    dist (p n) q ≤ ∑' k, δ (n + k) :=
  dist_le_tsum_of_dist_le_of_tendsto δ hδ hs hq n

/-- (b) Hypothesis (2) closed under the limit: nonnegative sector
weights stay nonnegative (coordinate evaluation is continuous; the
matrix/PSD version is GramCompactness.posSemidef_of_tendsto). -/
theorem limit_nonneg {p : ℕ → (Fin m → ℝ)} {q : Fin m → ℝ}
    (hq : Tendsto p atTop (𝓝 q)) (hp : ∀ n i, 0 ≤ p n i) :
    ∀ i, 0 ≤ q i := by
  intro i
  have h : Tendsto (fun n => p n i) atTop (𝓝 (q i)) :=
    ((continuous_apply i).tendsto q).comp hq
  exact ge_of_tendsto' h fun n => hp n i

/-- (c) Hypothesis (4) commutes with the limit: sector aggregation
(the weight of any finite set of characters) converges to the
aggregate of the limit. -/
theorem sector_tendsto (S : Finset (Fin m)) {p : ℕ → (Fin m → ℝ)}
    {q : Fin m → ℝ} (hq : Tendsto p atTop (𝓝 q)) :
    Tendsto (fun n => ∑ i ∈ S, p n i) atTop (𝓝 (∑ i ∈ S, q i)) :=
  tendsto_finset_sum S fun i _ =>
    ((continuous_apply i).tendsto q).comp hq

/-- THE MASTER CORE, assembled: hypotheses (2) + (5) over the sector
data (4) give a limit state that is nonnegative, with every sector
aggregate converging to a nonnegative limit sector.  Uniqueness is
`limit_unique`; the quantitative tail is `tail_bound`. -/
theorem master_core (p : ℕ → (Fin m → ℝ)) (δ : ℕ → ℝ)
    (hpos : ∀ n i, 0 ≤ p n i)
    (hδ : ∀ n, dist (p n) (p (n + 1)) ≤ δ n) (hs : Summable δ) :
    ∃ q : Fin m → ℝ, Tendsto p atTop (𝓝 q) ∧ (∀ i, 0 ≤ q i) ∧
      ∀ S : Finset (Fin m),
        Tendsto (fun n => ∑ i ∈ S, p n i) atTop (𝓝 (∑ i ∈ S, q i)) ∧
          0 ≤ ∑ i ∈ S, q i := by
  obtain ⟨q, hq⟩ := exists_limit p δ hδ hs
  have hqpos := limit_nonneg hq hpos
  exact ⟨q, hq, hqpos, fun S =>
    ⟨sector_tendsto S hq, Finset.sum_nonneg fun i _ => hqpos i⟩⟩

/-- μ4 character-orthogonality head: Σ_{j<4} i^j = 0 — the exactness
token for the character projectors P_c = (1/4) Σ_j i^{-cj} U^j both
instances deploy (full projector algebra measured in the probes). -/
theorem mu4_char_sum : (∑ j : Fin 4, (Complex.I : ℂ) ^ (j : ℕ)) = 0 := by
  simp [Fin.sum_univ_four, pow_succ]

end PositiveDescentMaster
end TfptCarrier

#print TfptCarrier.PositiveDescentMaster.master_core
#print axioms TfptCarrier.PositiveDescentMaster.master_core
#print TfptCarrier.PositiveDescentMaster.tail_bound
#print axioms TfptCarrier.PositiveDescentMaster.tail_bound
