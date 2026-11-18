# SOURCE.MEDIATOR.COVARIANCE.01 — the family representation of the Gate-B mediator

Date 2026-10-05. Research contract under the experiments firewall: no promotion.
Run: `python3 -B experiments/theory-contracts/mediator-covariance-20261005/covariance.py` (7 steps, about 2 s).

**Verdict: R_PHI_DETERMINED_FOR_V1.** The mediator need not be a family singlet; the whole exchange must be invariant.

**Three-family vertices.** They are written on the pair labels S²(3), with P_{A,ff'} symmetric in f, f'.

| Vertex | Definition | Symmetry |
|---|---|---|
| V1 | decoupled, Σ_f O_{A,ff}² | Z₄³⋊S₃ and A4-monomial; not SO(3), not a rotated A4 |
| V2 | all symmetric pairs | SO(3); breaks per-family Z₄ |
| V3 | interfamily, Σ_{f≠f'} O_{A,ff'}² | — |

Z₄³⋊S₃ admits exactly one ΔN=4 family tensor, Σ_f e_f^{⊗4}. So V1 is *the* ΔN=4 structure of the decoupled symmetry.

**Pair-label representations.**

| Group | Diagonal pairs | Off-diagonal pairs |
|---|---|---|
| A4 (monomial) | 1 + 1′ + 1″ | 3 |
| SO(3) | S²(3) = 1 + 5 | (same) |
| Z₄³⋊S₃ | irreducible signed-permutation triplet | — |

**Mediator.**
- Eliminating a degenerate multiplet gives the projector onto its span.
- V1 therefore needs **one real Z₄-odd Spin(10) vector per family**, degenerate and unmixed. That is a (10, 3_signed-perm) of Z₄³⋊S₃, which equals 1 + 1′ + 1″ under A4.
- A (10,3) triplet in A4's own vector representation produces the interfamily vertex V3 instead.
- The formula Σ_f Φ_{A,f} O_{A,f} → −(λ²/Δ) Σ_f O_{A,f}² is right with the signed-permutation meaning of the index f. It is not right with the A4-vector meaning.

**Basis condition.** The v117 transport maps V1 to itself only in the A4-monomial basis (from `family-symmetry-gateb-20261005`). It does not do so in the original family-mode basis.

**Source link (representation level only).** [J(16,4) J(16,4)] in (10, S²4) contains S²(3) = 1 + 1′ + 1″ + 3 under A4. These are
per-family composites Φ_{A,f} ~ [J^f J^f]_10 of the *light* families. The omitted-line singlet S²(1) is not needed.
Still open:
- the bulk Z₄ charge;
- the gap and residue;
- the coupling λ;
- the Schrieffer–Wolff residual;
- a gapping mechanism, since the seam currents are massless.

Files: `covariance.py`, `covariance_results.json`, `contract_index.json`.
