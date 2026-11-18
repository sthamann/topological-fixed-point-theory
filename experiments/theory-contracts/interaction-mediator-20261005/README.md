# SOURCE.INTERACTION.MEDIATOR.01 — origin of −g Σ_A O_A² by eliminating a heavy mode

Date 2026-10-05. Research contract under the experiments firewall: no promotion.
Run: `python3 -B experiments/theory-contracts/interaction-mediator-20261005/mediator.py` (5 steps, about 3 s).

**Question.** Can the Gate-B vertex arise from a quasi-free light sector coupled to a heavy source mode? Does the
source contain that mode?

**Verdict: MEDIATOR_SPECIFIED_NOT_IN_ELEMENTARY_SOURCE.**

**What is shown.**
- **Mechanism.** H = H₀ + Δ Σ b†b + λ Σ (b+b†)O_A gives H_eff = H₀ − (λ²/Δ) Σ O_A² + O(1/Δ).
  - With a single channel and H₀ = 0 this is exact (completed square; error 5·10⁻¹⁴).
  - With two noncommuting channels plus hopping, the error at fixed g = λ²/Δ is 0.0088, 0.0045, 0.0023 for Δ = 10, 20, 40, i.e. it falls like 1/Δ.
- **Mediator dictionary.** This is on V = 16⊗2. The invariant symmetric bilinear form B(u,v) = ⟨z, u∧v⟩ on Λ²V (z the unique ΔN=4 invariant) has these properties:
  - It is block-diagonal in the channels.
  - It is nondegenerate on (10,1), rank 10, and on (120,3), rank 360.
  - It vanishes on (126,1).

  The resulting channel table:

  | Mediator | Generated vertex |
  |---|---|
  | complex (charge-2) | only the U(1)-invariant \|P_R\|² |
  | real 10 | \|P_10\|² + the ΔN=4 term |
  | real (120,3) | \|P_(120,3)\|² + the ΔN=4 term |
  | 126 | only complex, no ΔN=4 term |
- **Θ₁'s ray.** Σ_A O_A² corresponds to exactly **one heavy real (10,1) mode**: an SU(2)-spin singlet, Z₄-odd (Φ→−Φ under ψ→iψ), coupled on site to the filtered pair operator. The Hubbard–Stratonovich field of the running QMC has the same quantum numbers, but it is only a mathematical rewriting of the quartic, **not evidence of a physical mode**.
- **Source inventory.** The E8 roots split under D5×D3 as 40 (45,1) + 12 (1,15) + 60 (10,6) + 128 spinor roots.
  - The only bosonic Spin(10)-vector currents are (10,6). They carry a family SO(6) vector, which has no invariant vector under the family SO(3).
  - They are weight-1 seam currents, i.e. massless, not heavy.
  - So the current algebra contains no elementary heavy family-singlet real 10.

**Consequence.** The interaction question is now a concrete search for one specific mode: a heavy real Spin(10) vector,
spin singlet and Z₄-odd, with g = λ²/Δ. It could be a composite or come from an extended source. Composites are not
constructed here.

Files: `mediator.py`, `mediator_results.json`, `contract_index.json`.

**Acceptance chain for a physical mediator (all required; no step is met yet).**
1. A source operator 𝓜_A in (10,1) that is bosonic, a spin singlet and Z₄-odd. The representation sieve is in `mediator-rep-sieve-20261005`.
2. A spectral pole: ⟨𝓜_A(t)𝓜_B(0)⟩ ~ δ_AB Z_Φ e^{−Δt}, with Z_Φ > 0 and 0 < Δ < ∞, and Δ ≫ E_IR.
3. A nonzero coupling λ to O_A, i.e. its matrix element between the mediator and the filtered pair operator.
4. Schrieffer–Wolff/Feshbach: H_eff^{(2)} = −(λ²/Δ) Σ O_A² + δH with ‖δH‖ ≪ λ²/Δ. No equally large unwanted operators may appear; otherwise Θ₁ is still not selected.
