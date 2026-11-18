# SOURCE.MEDIATOR.REP_SIEVE.01 — which source modes could build the (10,1) mediator at all?

Date 2026-10-05. Research contract under the experiments firewall: no promotion.
Run: `python3 -B experiments/theory-contracts/mediator-rep-sieve-20261005/sieve.py` (5 steps, about 1 s).

**Target.** A bosonic Spin(10) vector that is a family singlet. Spin-singlet is automatic, because source modes carry no spacetime spin.
**Building blocks.** The NS fermions ψ (10,1) and χ (1,6), and the E8 currents (45,1), (1,15), (10,6), (16,4), (16̄,4̄).

**Verdict.**
- **SU(4) neutrality:** forbidden at every degree.
- **SU(3), SO(3) or A4 neutrality:** minimal degree 2.

**What is shown.**
- **Centre selection rule.**
  - All 240 E8 roots satisfy q₁₀ + q₄ ≡ 0 (mod 4). This is the trivial diagonal centre already found in UR.SOURCE.OPERATOR_ORIGIN.01.
  - The NS fermions have q₁₀ + q₄ ≡ 2.
  - Hence every bosonic composite has q₁₀ + q₄ ≡ 0, while (10, SU(4)-singlet) needs 2. It is excluded at every degree; exhaustive enumeration up to degree 6 serves as the check.
- **Your candidate (16,4)(16,4)(1,6).** It has the right centre charges (2,0), but with χ it is fermionic. A bosonic (1,6) would violate the rule and does not exist in the source.
- **Spin(10) multiplicities of the 10** (numerical check):

  | Product | Multiplicity of 10 |
  |---|---|
  | S²16 | 1 |
  | Λ²16 | 0 |
  | 16⊗16̄ | 0 |
  | 10⊗10 | 0 |
  | 10⊗45 | 1 |
- **Degree-2 candidates under smaller family groups.** Counts of H-singlets:

  | Candidate | SU(3) | SO(3) | A4 |
  |---|---|---|---|
  | [J(16,4)J(16,4)] in (10, S²4) | 1 | 2 | 2 |
  | [J(10,6)J(1,15)] | 2 | 6 | 8 |

  - Under SU(3), the singlet of S²4 is the omitted-family-line component.
  - Fermion bilinears never give (10, H-singlet).

**Not decided here (next gates).** Representation theory leaves open:
- the bulk Z₄ charge, which needs the source→bulk map;
- the gap Δ and residue Z_Φ in ⟨𝓜_A(t)𝓜_B(0)⟩ ~ δ_AB Z_Φ e^{−Δt};
- the coupling λ to O_A;
- the Schrieffer–Wolff residual ‖δH‖ ≪ λ²/Δ.

The E8 currents are massless on the seam, so a heavy mode needs a gapped mechanism beyond the seam CFT.
The structural hint that the candidate lives in the omitted family line is not a claim that the mode exists.

Files: `sieve.py`, `sieve_results.json`, `contract_index.json`.
