# SOURCE.WEYLCOIN.ORIGIN.01 — is there a neutral local C² in the source?

Date 2026-10-05. Research contract under the experiments firewall: no promotion.
Run: `python3 -B experiments/theory-contracts/weylcoin-origin-20261005/coin_origin.py` (6 steps, about 15 s).

**Question.** Is there a local, family- and gauge-neutral two-dimensional factor in the existing source on which the
spatial process acts nontrivially? This asks only for an M₂(C) block of the parity-even neutral commutant. It does not ask for SU(2) or Weyl structure.

**Verdict: NO_NEUTRAL_COIN_IN_SOURCE.** The local algebra is the eight complex seam modes of U: 5 carrier modes and 3 family modes, i.e. C³²⊗C⁸.

- **Spin(10) gauge (TFPT).**
  - For family SU(4), SU(3), SO(3) and A4, every pair of irreps occurs at most once per fermion parity.
  - The parity-even neutral commutant is abelian.
  - The matter 16 carries no two-dimensional multiplicity.
- **Weakened gauge (SU(5) or SM).**
  - Neutral M₂ blocks appear, but only between the carrier-empty and carrier-full states (Λ⁰ vs Λ⁵).
  - These two states are 16 vs 16̄ under Spin(10). The mixing is an eight-fermion operator in the carrier-singlet sector and carries no matter.
- **Process action.**
  - The source process operators act on the family factor only. These are the channel-A family Majorana products and the second-quantised Γ(M_a).
  - On every neutral M₂ block they act diagonally, and their neutral part is abelian.
- **Other candidates.**
  - The sheet doublet is the parity grading itself (Part K9), and all marks act by the same swap.
  - The flavor fiber and the U(1)_N-only multiplicities are family-carrying.

**Consequence.** Together with `mark-projective-lift-20261005` this makes AX.WEYLCOIN the minimal additional microscopic
input within the declared local algebra. AX.WEYLCOIN is a two-component internal-singlet one-particle coin coupled to the three H₁ directions. Freezing it as the third axiom is the user's decision. An extended source or a different
local algebra (e.g. after a source→bulk map) is not excluded.

Files: `coin_origin.py`, `coin_origin_results.json`, `contract_index.json`.
