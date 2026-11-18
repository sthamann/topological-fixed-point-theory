# SOURCE.FAMILY.SYMMETRY.AT_GATEB.01 — which family group is exact when Gate B acts?

Date 2026-10-05. Research contract under the experiments firewall: no promotion.
Run: `python3 -B experiments/theory-contracts/family-symmetry-gateb-20261005/family_symmetry.py` (5 steps, about 3 s).

**Verdict: REAL_FAMILY_GROUP_REQUIRED.** This is conditional on Θ₁ keeping its ΔN=±4 vertex for three light families.

**What is shown.**
- **Premise.** One 16 has U(1) anomaly sums Σq = Σq³ = 16. A symmetric mirror gap therefore needs U(1) → Z₄, which is what the ΔN=4 part of Σ_A O_A² does.
- **Centre criterion.** Scalar elements of a family group act on the light families as bulk U(1) phases. Hence G ∩ U(1)·I ⊆ μ₄.

  Invariants of 3^{⊗4}, the family part of every ΔN=4 operator:

  | Group | Invariants | Allowed |
  |---|---|---|
  | U(3) | 0 | no |
  | SU(3) (centre Z₃) | 0 | no |
  | SO(3) | 3 | yes |
  | A4 | 7 | yes |
  | signed permutations B3 | 4 | yes |
  | Z₄³⋊S₃ | 1 | yes |
- **Source layers.**
  - The seam has Spin(6) = SU(4).
  - The matter selection (family number 2) leaves u(3), which contains the overall U(1).
  - In the v117 transport, the centraliser is only U(1) (Schur) and the transport is of real type.
  - In a computed basis, the *A4-monomial basis* (`a4_monomial_basis.npy`), all twelve monodromies are signed permutations.
  - In that basis the deck intertwiner is real up to a phase.

**Consequence.** At the Gate-B stage the family structure must already be real. The ordering is

seam SU(4) → matter selection U(3) → real flavor structure (transport) → Gate B → flavor/Yukawa.

The SU(3)-singlet candidate of the representation sieve (omitted line) refers to a group that cannot be exact at Gate B.
The SU(4) prohibition of the sieve is therefore not the decisive constraint. The relevant classification is under real
family groups; see `mediator-covariance-20261005`.

Files: `family_symmetry.py`, `family_symmetry_results.json`, `a4_monomial_basis.npy`, `contract_index.json`.
