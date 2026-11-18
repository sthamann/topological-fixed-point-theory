# THETA1.UNIVERSALITY.CLASS.01 — is Gate B a property of a region in (b, θ, γ)?

Date 2026-10-05. Research contract under the experiments firewall. **Status: PREREGISTERED_NOT_RUN.** No QMC was launched.
Gate B is undecided at the frozen point, and the running Gate-B jobs are untouched.
Run the prechecks: `python3 -B experiments/theory-contracts/theta1-universality-class-20261005/prechecks.py` (9 steps, about 2 s).

**Design (frozen 2026-10-05T14:28:54).** The points are (b, θ in degrees, γ):

| Point | (b, θ, γ) |
|---|---|
| P0 (frozen Θ₁) | (0, 120, 0) |
| P_b± | (±0.25, 120, 0) |
| P_θ− | (0, 100, 0) |
| P_θ+ | (0, 140, 0) |

- **Run condition.** Gate B must first reach CONDITIONAL at P0 under TFPT.FULL.PROCESS.01 G1–G5, including the round-56 equilibration rule and SPG.BUDGET.01.
- **Run protocol.**
  - Use the same g, β, dτ, cells and estimators as at P0.
  - Run at least 2 chains per point.
  - First reproduce the free reference energies listed below.
- **Decision rule.**

  | Outcome | Condition |
  |---|---|
  | UNIVERSAL_CANDIDATE | all five points pass |
  | FINE_TUNED | P0 passes and at least one neighbour fails cleanly |
  | UNDECIDED | sign or feasibility limits |

  No larger scan and no retuning before this decision.

**Amendment 1 (2026-10-05T14:49:16, before any data).** There are two arms.
- **Arm A (fixed g).** It tests what the same fundamental parameter does in different runtime realisations.
- **Arm B (matched mirror coupling).** It sets g(b) = g₀/(1+b)⁴, so that g(b)(1+b)⁴ = g_mirror(P0):

  | Point | g(b)/g₀ |
  |---|---|
  | P_b+ | 0.4096 |
  | P_b− | 3.1605 |

- **Universality decision.** UNIVERSAL_CANDIDATE requires arm B to pass at all five points. Arm A is reported separately.
- **Velocities.** The Γ and mirror Weyl speeds are exactly 1 at every point (V_Γ = I, V_P orthogonal), so no renormalisation is needed; they are reported.
- **Feasibility.** P_b− in arm B needs 3.16·g₀ and may be sign-limited.

**Exact prechecks (9/9).**
- **Nodes.** All five points are minimally doubled: Γ (+1) and P (−1).
- **Cells.** The mirror stays at P, so the existing P-compatible cells remain valid.
- **Filter.** f(Γ) = 0 with Hessian I/2 at every point. f(+P) = 1+b and f(−P) = 1−b.
- **Normalisation.** diag(FF†) = 9/8 + b²/8, so g* = g·(9/8 + b²/8)².
  - g*/g = 1.2656 at b = 0 and 1.2833 at b = ±0.25.
  - Up to this rescaling, the one-site so(11) condition holds everywhere.
- **Mirror weight.** The mirror sees the quartic with weight (1+b)⁴: 2.441 at P_b+ and 0.316 at P_b−. Report g_mirror separately.
- **IR.** The IR is identical at every point (unit isotropic Γ cone, f(Γ) = 0 quadratic), and the IR protection is unchanged.
- **Free reference energies.** Primitive L = 4 / L = 8 cells, Γ excluded:

  | Point | L = 4 | L = 8 |
  |---|---|---|
  | P0, P_b± | 1.414214 | 0.717439 |
  | P_θ− | 1.250945 | 0.704036 |
  | P_θ+ | 1.282682 | 0.715766 |

**Note.** The points with b ≠ 0 have complex filter hoppings (c_W = −1/8 − ib/8). The production code must be checked for
complex filtered fields before any run.

Files: `prechecks.py`, `prechecks_results.json`, `contract_index.json`.
