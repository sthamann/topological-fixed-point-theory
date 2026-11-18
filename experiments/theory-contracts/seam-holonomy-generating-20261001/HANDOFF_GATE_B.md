# Handoff: TFPT Gate B (chiral matter / symmetric mirror decoupling) — QMC state of 2026-10-05

Written for a second agent that continues this work in parallel on a large CPU/GPU server.
Authoritative records: `contract_index.json` in this folder (contract UR.SOURCE.HOLONOMY_GENERATOR.01),
`../tfpt-full-process-20261002/contract_index.json` (master, gate `B_chiral_matter`), `../../next.txt` (rounds 38–57),
the report `tfpt_explorer/docs/session-2026-10-01/TFPT_Sitzungsbericht_2026-10-01.tex` (sections Runde 42–57).
Project rules: no promotion to ledger / papers / website, no git commits without the user, the theory graph only via
`python3 -B verification/build_theory_graph.py` (+ `--check`), never hand-edit generated graph files.
This folder is **untracked in git** — copy it (code + `*.jsonl`, `*.npy`, `*_results.json`, ~56 MB) to the other machine.

## 1. What we want to show

TFPT builds 3 x 16 matter, 3+1D space-time and Weyl kinematics from its source. The lattice Weyl theory has the wanted
Weyl fermion at Gamma and an unavoidable mirror (doubler). **Gate B** asks whether one fixed local Hamiltonian removes the
mirror from the low-energy physics **without breaking symmetries**, keeping Gamma massless (symmetric mass generation, SMG).

Frozen model Theta_1 (`reweight.py`, structure hash 2b23ba1ba94c9058, integrator hash **f0f2d58fada5669b** — must not change):
minimal-doubling C3 Weyl hopping h_W (nodes Gamma and P), filter F (f(Gamma) = 0, protects Gamma), interaction
H_int = -g sum_x sum_A O_A(x)^2 (ten filtered Spin(10) pairing composites), HS-decoupled, exact Pfaffian / determinant
weights, phase-quenched sampling with continuous phase tracking. Flavour twist delta = lambda (1,1,1) (Spin(10) Cartan
Q = Sigma_12) gives continuous momenta; the mirror branch is measured in the Q+ sector at P + delta ("P+d").

Gate B criterion (master contract): G1 only the desired Weyl sector light near Gamma; G2 gap outside; G3 no active sector
below the Weyl continuum; **G4 no symmetry breaking (Spin(10), Z4, translations)**; G5 V -> infinity, beta -> infinity.

## 2. Where we are (physics)

* Gamma protected (analytic sum rule; exactly free in QMC).
* Earlier "mirror light up to g ~ 0.024" was a **burn-in artefact** (8–10 burn-in sweeps). The slow mode is the overall
  HS amplitude |phi|^2 (within-chain corr with G(beta/2) -0.93 ... -0.99); local sweeps need hundreds of sweeps.
* With the validated global amplitude move: **g = 0.021, beta = 4 (16 sites) equilibrates into a strongly paired state**
  (|phi|^2 ~ 2.3 x prior, G_{P+d}(beta/2) = 0.035–0.039 vs 0.34–0.38 in short runs: mirror line strongly suppressed).
* Umbrella/WHAM at g = 0.015, beta = 6: the paired state (|phi|^2 ~ 45) is lower in free energy than the weak region in
  every WHAM variant (sign robust; magnitude not controlled because same-centre windows disagree); no barrier resolved.
* g = 0.010: weak state is the equilibrium (paired start relaxes back). g = 0.012: two stationary branches, undecided.
* **Update 2026-10-06 07:30 (G4 trend reversal):** the uniform pairing order parameter S_0/N (K = 0) grows slowly and
  steadily in every long chain (16 sites 0.2 -> 3.0 after 250-300 sweeps, still rising; 3 of 4 8-site chains condensed;
  32 sites rising). A uniform Spin(10)-vector pairing condensate is now the more likely outcome (G4 failure for this phase).
  New targeted move: `global_moves={"pcn_uniform_rho": r}` (env `EQA_GLOBAL="eps,rho,rho_uniform"`), validated (GM5).
* **Update 2026-10-07 (G4 almost decided):** umbrella in the K = 0 order parameter X = |M|^2 / N (bias kind 'uniform',
  uniform-mode scaling move `uniform_scale_eps`, env `EQA_GLOBAL="eps,rho,rho_uniform,eps_uniform"`): 16 sites minimum at
  X* ~ 16-20; 32 sites every window runs to X ~ 10-26 with strongly negative mean force -> no 1/N decay, long-range
  uniform pairing order indicated (G4 FAIL expected for the paired state; formal verdict pending stationarity).
* **G4 (decisive, open):** volume-normalised order parameters S_0/N (uniform pairing) and S_max(K != 0)/N at g = 0.021,
  beta = 4: 8 sites inconsistent between chains (one chain formed a uniform condensate), 16 sites ~1.6–2.6, 32 sites
  **still rising** (0.75–1.04 after ~110 sweeps; phase R recovering to 0.77–0.82). If 32 saturates near the 16-site value
  -> long-range order (Gate B fails for this phase); if clearly below -> symmetric paired state (SMG candidate).
* Active mirror-sector gap: the HS sector projection is **blind** in the paired state (every configuration breaks U(1) and
  translations; all 64 sector fractions = 1/64): SECTOR_TRACE_INACTIVE_IN_PAIRED_HS. Replacement preregistered:
  MIRROR.CORRELATOR.BASIS.01 (GEVP of elementary + cubic composite operators, Wick per HS configuration) — not built yet.

## 3. Code map (all in this folder; Python 3.14, numpy 2.4, scipy 1.17; BLAS threads via OPENBLAS/VECLIB_NUM_THREADS)

| file | role |
|---|---|
| `reweight.py` | **production QMC**: `Model` (cells debug4/cell8/cubic16/oh32/prim64, twist), slices, stabilised UDV chains, `greens_k`, `run()` with options `hook`, `outer` (envelope), `burn`, `monitor`, `save_phi`, `phi0` (continuation), `global_moves` (amplitude move with Jacobian s^D + pCN), `bias` (umbrella). Without options bitwise identical to the original sampler. |
| `qmc_pfaffian.py` | Pfaffian Gaussian algebra, PF1–PF7 (incl. stabilised HMC force PF5 — not used in production yet) |
| `equilibration_audit.py` | runner `run <beta> <g> <lambda> <sweeps> <seed> <tag> [start.npy] [append] [offset]` with per-sweep monitor (phase, log|W|, field low modes, G(tau) every 4, sector traces every 8; env `EQA_GLOBAL="eps,rho"`, `EQA_FULL_F=1`); analysis applies the EQUILIBRATED rule |
| `volume_test.py` | G4 volume test runner `run <cell> <beta> <g> <lambda> <sweeps> <seed> <tag> [start.npy]` (env `EQA_GLOBAL`, `VT_G_EVERY`, `VT_OFFSET` for continuation); analysis VT1–VT3 (VT3 = predeclared block rule: G4_PASS_CANDIDATE / SOLVER_LIMIT / CONTINUE) |
| `umbrella.py` | umbrella windows in |phi|^2 `run <beta> <g> <lambda> <centre> <kappa> <sweeps> <seed> <tag> [start.npy]`; WHAM (per-window kappa), binwise phase correction, Delta F with block bootstrap (UM1–UM3) |
| `global_moves.py` | validation of the amplitude / pCN moves (GM1–GM3) and benchmark (GM4) |
| `tempering.py` | parallel tempering in |W|^s (toy-validated; impractical here: ~20 replicas needed) |
| `spg_production.py` | sector-projected Green function / 64 sector traces (twisted translations x Z4, sign-exact path continuation), budget test SPG.BUDGET.01 (B1–B7), free and HF references |
| `mixing_mode.py`, `sector_sampling_feasibility.py`, `sector_envelope.py` | mixing diagnosis, projection feasibility, envelope sampler (rejected: 7x worse) |
| `groundline_residue.py`, `spectral_bounds.py`, `zero_bracket.py`, `sector_gap.py`, `symmetric_transfer.py`, `mirror_node.py`, `filter_ir.py`, `theta1_md.py` | earlier analysis / method layers (groundline candidates now BURN_IN_SUSPECT) |
| `audit.py` | dependency audit of all claims (must print AUDIT OK) |

Data: `eqa_monitor_*.jsonl`, `vt_monitor_*.jsonl`, `um_monitor_*.jsonl`, `spg*_points_*.json` (raw time series), `*_phi_*.npy`
(final HS configurations for continuation; paired-state starts: `um_phi_g015b6c46up.npy`, `eqa_phi_g015b6GB.npy`,
`eqa_phi_g021b4GB.npy`, `vt_phi_vt32{A,B}.npy`). Twist per coupling: g 0.010 -> 0.21, 0.012 -> 0.25, 0.015 -> 0.32,
0.018 -> 0.39, 0.021 -> 0.35, 0.024 -> 0.37. dtau = 0.25 everywhere.

Rules that must hold for any new number: >= 2 chains with different starts, burn-in until stationary (>= 150 sweeps with
global moves), Gamma-method / chain-level errors, phase R reported, rules fixed **before** looking at the data, retractions
recorded in the contracts. Costs on the M3 Ultra (32 cores): cubic16 beta 4 ~25–40 s/sweep, beta 6 ~40–60 s, oh32 beta 4
~6–8 min/sweep (dense complex QR / LU / SVD of 1024^2 resp. 2048^2 Nambu matrices dominate).

Running here right now (do not duplicate): 32-site G4 extension `vt32A`, `vt32B` (g 0.021, beta 4, lambda 0.35, until ~22:00).

## 4. What a big CPU/GPU server should do (in priority order)

1. **G4 on a homothetic family.** cubic16 = 4 I (16 sites) and cubic128 = 8 I (128 sites) are P-compatible and
   geometrically similar (add `"cubic128": 8*I` to `CELLS`; 6 I is not P-compatible). Plus prim64 as a third cell.
   g = 0.021, beta = 4 (and g = 0.018), >= 4 chains per cell, global moves, until stationary; report S_0/N,
   S_max(K != 0)/N, channel anisotropy, **and a density-channel structure factor** (not monitored yet). 128 sites means
   8192 x 8192 complex Nambu matrices: GPU territory.
2. **GPU port of the hot path** (CuPy or JAX, complex128): `slice_op` (SVD), `_udv_left` (pivoted QR), `logdet_split`,
   `_inv_sum`, `greens_k`. Batch many chains / umbrella windows per GPU. Keep the integrator hash: port the linear
   algebra, not the physics definitions; verify bitwise-equivalent statistics with `sector_envelope.py` ES2-style checks
   and `global_moves.py` GM3.
3. **Better updates for large V:** (a) HMC with the existing stabilised force (qmc_pfaffian PF5) for global moves of all
   slices at once; (b) site-local HS updates with Woodbury low-rank Green-function updates (classic DQMC) instead of
   whole-slice updates. Benchmark by error^2 x GPU-hour on g = 0.021, beta = 4, 16 sites against the current sampler.
4. **MIRROR.CORRELATOR.BASIS.01** (preregistered in `contract_index.json`): implement the elementary + cubic composite
   correlation matrix by Wick per configuration, GEVP with tau_0 = 0.5, beta = 4 and 6, 16 then 32 sites, in the paired
   state (start from the saved paired configurations). This is the active-sector gap test.
5. **Gamma in the paired state:** run untwisted chains from paired configurations, measure G_Gamma(tau) and the dispersion
   at the smallest momenta (desired Weyl cone must survive).
6. **Phase diagram at equilibrium:** g in {0.010, 0.012, 0.0135, 0.015, 0.018, 0.021, 0.024}, beta in {4, 6, 8}, both
   starts, global moves; umbrella / WHAM where two branches persist (more and longer windows across 38–41 to resolve or
   exclude a barrier); Delta F(g) sign to bracket the onset.
7. **Second solver (only if G4 or the gap stays undecidable):** SOLVER.FEASIBILITY.01 — symmetry-resolved Hamiltonian
   solver (neural quantum state / VMC on 32 V occupation bits on GPU, or orbital tensor networks); must first reproduce
   g = 0.010 / 0.015 / 0.021 QMC results.

Deliverables back: new `*_results.json`, monitor files, contract entries (with retractions if any), audit OK.
