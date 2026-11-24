# OriginLab

[Autonomous F1 supervisor](AUTONOMOUS-F1.md)

Independent, finite exact research campaigns for TFPT. Previous search engines
are neither imported nor executed. Historical certificates supply benchmark
data. The original-source lane calls the native R113 consumer directly.
The native domain and canonical Atlas supply original provenance.

**Routes (23 November 2026).** The original R113 endpoint call is one route
(`original_source`), not the gate of the search. Source-free generators run
first; a missing R113 packet blocks G3 of its own route only and no gate of
another origin model. Many sections below report that a lane "supplies no
R113 mark": that is a scoped statement about the original route, not a
defect of the lane. A new origin model is checked through its own candidate
contract ([F1-CANDIDATE-CONTRACT.md](F1-CANDIDATE-CONTRACT.md)).

## Process rule search, 24 November 2026

**Current certificate update after H6h27.** The paragraphs below describe the historical
search stages. The current switches preserve candidate classes, not byte-identical old
JSON after mathematical corrections. The original result files remain unchanged.
The new [source and dynamics report](../../experiments/exploration/2026-11-24-h6h26-source-dynamics/VERDICT.md)
records the full KMS condition, the distinction between current raising and the Hardy
shift, and the positive quadratic parent with weight n+1. Its full metric star has a
discrete spectrum, with three additional interlacing sectors. The same lane now computes
physical connected boundary four-point differences, exact local depth certificates,
and group-independent obstructions to the specified exact TM1 pattern.

Clauses for the analytic kernels, transfer index and HR potential class use exact
certificates. Taylor and finite spectral probes are diagnostics. Source centering is
separated from Hamilton-vacuum subtraction. The historical cost and actual pole count
are separate, and every selection reports all matches as well as minimum-cost matches.
Correctly applying mirror weight to the transfer at q=1/sqrt(2) leaves two selected
families in the old scale-free comparison. A single physical process is not claimed.
Undirected HR edges leave one quartic parameter per edge; C4 relates the arms at fixed
depth, while a common parameter across depths is an additional homogeneity rule.
The new parent comparison still calls genesis_rules.search. It does not add a second
search engine. Validation: 153 targeted tests, including the merged H6h27 tests.

`process_rules.py` turns the genesis engine (`genesis_rules.search`) on the open slots of the TFPT process instead of
on number-free structures. Four universes hold computed candidates: positive matrix kernels det(I - W^+ Z)^-nu with an
optional det Z conj det W term, spin-2 inverse propagators (Dirac spectral action with three cutoffs, the ledger's entire
ansatz, local truncations, Einstein, Schur complements of positive two-channel pencils), seam stacks of real or complex
Chern layers with equally or reversely oriented second sheet half, and pairs of residual family symmetries in S4. A fifth,
composed universe takes one candidate per slot (53 760 whole processes) and searches over nine cross-slot principles
(Hardy slice, conformal covariance, positivity everywhere, Dirac spin-2 correction, mirror sheet, real seam, E8 seam,
family mixing, family clock). Run time about three seconds.

First run: the longest minimal rule has eight principles, each load bearing, and selects one whole process: kernel
nu = 1, spin-2 kernel C = 1/x + 1/(x + 10/3), sixteen real Majorana layers with reversed second sheet, family pair
Z3 x Z4. One principle (mirror sheet) fixes two slots. The selection is a statement about the universes and the chosen
principles, not a derivation; the flavour slot is compared with data outside the lane (H6 exploration z43).

Second run (transfer extension, H6 exploration z44/z45): two new universes and two new principles. `transfer_processes`
enumerates a(z) = z^-d (1 - q z^k) on the Hardy space (40 candidates); quarter-turn covariance plus a kernel in one
grade fixes d = 1 and k in {4, 8}, but no clause fixes q (no UNIQUE verdict). `readouts` lists ten source maps of the
primitive quarter transfer; its minimal empty rules are {massless pole, connected} and {universal residue}: no rank-one
readout gives a massless, connected or universally coupled pole. The spin-2 slot gains the Marchenko-Pastur memory of
the transfer for five q; only q^2 = (8 - sqrt 55)/3 keeps a = 3/10. With `orientation_reversal` and `one_transfer`
the composed universe (67 200 processes) has a longest minimal rule of nine load-bearing principles; it selects
kernel nu = 1, the transfer memory at that q, the reversed real 16-layer seam and Z3 x Z4, and no longer needs the
mirror-sheet principle. `run(transfer=False)` restricts the current lane to those first-run candidate classes.

Third run (state rule, H6 exploration z46/z47): the number a = 3/10 behind `dirac_spin2` is a scale convention
(a = (3/10) f(0)/Int f in units of the cutoff), so it fixes no q. A new universe `q_rules` holds seven candidate states of
the transfer, with every clause computed from q. The original quarter Gaussian state (Calderon/Szego, b = 4 pi) has the
current norm omega(J_-1 J_1) = 1/(e^{4 pi} - 1); matching the ladder occupation q^2/(1-q^2) of the transfer's kernel vector
gives q = e^-2pi (clause `original_state`), and the Hardy state gives q = 0. The spin-2 slot gains the memories at e^-2pi
and 1/sqrt 2, and a stable closed form of 1/H - 1/x for small q. `composed_process.without_scale_convention` repeats the
whole-process search without `dirac_spin2` (72 576 processes). Eight principles then select kernel nu = 1, the transfer
memory in the original state, the reversed real 16-layer seam and Z3 x Z4. Without the state rule the choice falls back to
pure Einstein. The longest rules now report their tie set before the size order: n = 16 and n = 32 tie unless the E8-seam
principle is used. Flavour patterns are printed exactly in Q(sqrt 2, sqrt 3); the first runs printed (2 - sqrt 3)/6 as 2/45.
`run(state=False)` restricts the current lane to the second-run candidate classes.

Fourth run (larger family groups, H6 exploration z48/z49): `flavour_groups` enumerates residual symmetries in Delta(6 n^2),
n in {2, 3, 4, 6, 8}, with a Klein group for Majorana neutrinos and a cyclic group with three distinct eigenvalues for the
charged leptons (27 patterns). Every mixing pattern is TM2 (trimaximal middle column) at discrete angles; none meets the
rounded 3-sigma ranges, and {mixes, tm1} is empty. Delta(96), the smallest of these groups with the quarter phases i, carries
the family clock (TM2 at 15 degrees) with a Klein group, which S4 alone does not. Exact labels now reject large-denominator
fits and fall back to decimals.

Fifth run (metric stage, H6 exploration z50/z51): `metric_star(q, m)` builds a central metric with four arms of m metrics
(Einstein weights M_c and (M_c/2) q^(n+1), edge stiffness mu^2 M_child^2). Its canonical mass operator is unitarily
equivalent to the rectangular quarter transfer `hardy_rect(q, m)`, with exactly one kernel vector eta_v = M_v/M_eff, so
minimal node sources couple universally (eta_v/M_v = 1/M_eff). `readouts(metric=True)` renames the old clause `connected`
to `vacuum_subtracted` (the zero-mode-removed microscopic response of z44) and adds two computed readouts: the centred field
covariance (p^2 + mu^2 L)^-1, which is connected and massless, and the metric star with minimal node sources, which is
universal. `metric_potentials` computes Taylor data to fourth order of six Hassan-Rosen potentials with mpmath: flat
background, Fierz-Pauli mass and an undirected edge (same potential read from either metric) leave exactly V and
V +- det(I - X)/2 open, i.e. one quartic number per edge. `run(metric=False)` omits this metric stage.

```bash
.venv/bin/python -B research/originlab/process_rules.py > /tmp/process-rules.json
.venv/bin/python -B research/originlab/main.py --campaign process_rules --jobs 1 --out /tmp/originlab-process
```

A new slot is one enumerator with computed clauses, one `search(...)` call in `run`, and its clauses mapped onto the
principles in `composed_process`. Tests: [test_process_rules.py](test_process_rules.py).

## Foreign universes, 20 November 2026 (Exploration mode)

`foreign_universes.py` runs the same rule search in universes without TFPT vocabulary: all 142 connected graphs
on up to six vertices, the 88 elementary cellular automata and the 33 regular polytopes of dimension 2 to 8.
Invariants of every object go through one fixed expression grammar and are compared with the measured targets of
[`docs/research/EXPLORATION.md`](../../docs/research/EXPLORATION.md). A shuffled-invariant baseline (10 draws) decides
whether a target is hit more often than chance; with 11 targets in 3 universes, one weak lead is expected by chance.
First run: only one lead (graphs, τ/μ mass ratio, 5 objects against a baseline of 1.0), consistent with chance.

```bash
.venv/bin/python -B research/originlab/genesis_rules.py --foreign > /tmp/genesis-foreign.json
.venv/bin/python -B research/originlab/main.py --campaign genesis_foreign --jobs 1 --out /tmp/originlab-foreign
```

Add a universe as one enumerator with `clauses` and integer `invariants`, then one `search(...)` and one `score(...)`
in `run`. Tests: [test_foreign_universes.py](test_foreign_universes.py).

## Genesis rule search, 20 November 2026

The lane `genesis_rules` searches over rules instead of numbers. Each universe is a finite class of small structures,
enumerated exhaustively inside a stated bound (binary codes, two-part cuts of even self-dual lattices, lattices
generated by one vector and one clock, folds onto a sphere, sign words, relabelling groups). Every subset of a clause
pool is a rule; the smallest solutions decide EMPTY, UNIQUE or AMBIGUOUS, and each unique selection is reported with
its shortest rules and with what happens when one clause is dropped. No number of the record is an input; record
structures are named afterwards. It runs without the Atlas and belongs to the finite controls.

```bash
.venv/bin/python -B research/originlab/genesis_rules.py > /tmp/genesis-rules.json
python3 -B research/originlab/main.py --campaign genesis_rules --jobs 1 --out /tmp/originlab-genesis
```

A new universe is one enumerator and one `search(...)` call in `run`. See the
[first run, the blind hypotheses and the red-team verdicts](../../experiments/theory-contracts/genesis-20261120/README.md)
and [fourteen tests](test_genesis_rules.py). The run refuted its own best candidate: blind chirality does not survive
the self-dual completion once all three partner directions are probed.

## Structure scenarios without a model, 19 November 2026

Until now a scenario with a rule absent from the record could only come from
a model call (`--scenario-model`). The lane `structure_scenarios` generates
such scenarios deterministically. `structure_objects.py` enumerates
mathematical objects with computed features (sublattices of E8, E7, E6,
reflection groups, modular curves, elliptic surfaces, codes, anyons, boundary
theories, fermion modes and their splits, subset graphs, repeated steps,
matter content). An object may name the objects it contains by construction;
its closure enters the pool as well. `structure_scenarios.py` computes the
fingerprint from the Atlas itself: every integer and fraction that
established findings of at least three independent problems state. Objects,
pairs and triples are scored by the rarity-weighted fingerprint they carry
and compared with random fingerprints drawn from the pool, so the search
over the whole pool is part of each p-value. A layer decomposition and a
withheld test (select on two thirds, predict one third) follow. The best
scenarios are passed to `atlas_scenarios.assess` as declared hypotheses.

```bash
python3 -B research/originlab/structure_scenarios.py --fakes 1000 > /tmp/structure-scenarios.json
python3 -B research/originlab/main.py --campaign structure_scenarios --jobs 1 --out /tmp/originlab-structure
```

The supervisor runs the lane in every stage and queues its scenarios as
`atlas_scenario_case` jobs. Equal numbers are no identification; nothing is
promoted. A new object family is one more `gen_*()` function. See the
[first run and its limits](STRUCTURE-SCENARIOS-20261119.md) and
[twelve tests](test_structure_scenarios.py).

## Alternative origin scenarios, 16 November 2026

The new `f1_origin_scenarios` campaign compares five explicit hypotheses:
a spherical boundary response with radial time, changing quotient readouts
of the special Picard curve, a Lorentzian parent with E8 null reduction,
larger positive parents invisible to scalar torus data, and the newest
parallel session's spectral torus deck reader from commit `9511d244`.

The exact controls distinguish the full Picard symmetry of order 48 from
its marked cubic deck normalizer of order 12, preserve the native code
obstructions, separate the natural j=0 double torus cover from the square
j=1728 cover, and test complete quantum readout dynamics. They also expose
the scalar normalization freedom and the larger parent ambiguity. These
are alternative comparison sources; no common physical source or F1
closure is claimed. See the [full research and source audit](F1-ORIGIN-SCENARIOS-20261116.txt)
and [exact results with source hashes](F1-ORIGIN-SCENARIOS-20261116.json).

```bash
python3 -B research/originlab/main.py --finite-controls-only \
  --campaign f1_origin_scenarios --jobs 1 --out /tmp/tfpt-origin-scenarios
```

The results are included in the F1 frontier when both campaigns are requested.
The five component modules and the combined campaign use explicit failure
checks that remain active under `python -O`.
The latest torus study is included separately, with its twelve existing
tests and a clear distinction between the conical pillowcase quotient and
the smooth spherical scalar interface.
Two added numerical regressions protect the true small-time limit against
a fixed Fourier cutoff; elliptic fiber-product genera refer to normalization.
The final synchronized comparison also includes the six genus-five parent
tests from commit `8731e75`. Its theta separation uses a published theorem;
no numerical period matrix or new physical parent selection is claimed.

From the repository root with Python 3.11 or newer:

```sh
python3 -m pip install -r research/originlab/requirements.txt
python3 -B research/originlab/main.py --out research/originlab/runs/first
python3 -B research/originlab/main.py --out research/originlab/runs/first --resume --rounds 2
```

The default run executes the source-free generators first, then the
original R113 route, the calibration and model campaigns, P0 search and an
F1 frontier report. For one lane, use `--campaign`. A run records exact
witnesses, hypotheses, source hashes, the current local Atlas program
revision, bounds, versions, and seeds. A blocked source is an ordinary
scientific outcome. An exception produces an ERROR lane and a nonzero exit.
An interrupted round can be replayed without replacing its artifacts. Pass
`--watch-seconds N` with multiple rounds to run subsequent rounds when the
canonical research program or source configuration changes.

The finite controls are also runnable from a clean public checkout with
Python and SymPy alone.

```bash
python3 -B research/originlab/main.py --finite-controls-only --out /tmp/originlab-controls
```

This mode runs the finite controls for the double cell, moving Hardy carrier,
affine operators and polarization, weak fields and returns, four-mark geometry
and the family index. The regular search can also run without the local
Atlas. It reports the missing provenance and keeps available campaigns
running. Original-source campaigns require their matching, hash-checked
inputs. An unavailable original file produces a blocked result with its
path. The native R113 marks remain unset until the original weighted input
is available. Atlas publication uses the canonical local service.
The R006 benchmark can replay the compressed original inputs in
`fixtures/r006`. The checker verifies their original uncompressed hashes.

The finite OriginDSL expresses matrix processes on ordered boundary spaces.
The commands, their typing, and their denotations are in `language.py`. Its
declared binary primitives, finite dimensions, and boundary probes are inputs.
Six default grammar fragments use these inputs. Named boundaries, supplied
seams, user primitives, and joint grammar fragments can be added in a case.
These finite fragments do not constitute general
topology, arbitrary categories, a complete error correction calculus, or
general quantum field theory. Candidate matrices are deduplicated by exact
denotation in the given frames. The P0 planner generates ordered words on the
declared alphabet and ranks them by separated pairs per declared cost.
All-word equivalence for coherent models uses an exact reachable subspace.
The durable search restores informative candidates and negative certificates.

Use `--proposal-file terms.json` to submit JSON terms. An explicitly
configured `--model NAME` uses the authenticated Codex CLI. Candidates go
through the same parser and exact checks. A model is never a verifier. Pass
`--lean-toolchain leanprover/lean4:v4.29.1` to check the rank lemma and audit
its axioms. Z3 supplies the inverse integer check when installed. SymPy is the
reference arithmetic backend. See [plan](PLAN.md) for the remaining science
and backend choices.

## Equivariant P1 origin control (15 November)

Run `python3 -B research/originlab/main.py --campaign seam_cohomology --out /tmp/originlab-seam`
to audit the **actual** four-point meromorphic-section and logarithmic-differential
character representations. The quarter rotation decomposes the five
Riemann--Roch sections into two trivial and three distinct nontrivial
characters. The original three logarithmic differentials have the three
nontrivial characters, but differ by the reflection orientation sign.
A minimal three-character dual compensator creates three cyclic-neutral
channels; full untwisted D4 reflection leaves only one, while a reflection
orientation twist leaves two. The uncompensated three-dimensional sector
cannot commute with a full fundamental SU(3) action.

This is a new scoped F1 hypothesis and exact representation filter, **not**
the physical colour/family dictionary, original R113 packet or continuum
parent. It cannot be used to assign a physical source value. See
[the derivation and source conditions](P1-EQUIVARIANT-RESIDUE-ORIGIN-20261115.md)
and [independent test cases](test_seam_cohomology.py).
`--finite-controls-only` includes this control.

## Ramond spinor local-weight-one repair: canonical E8 chiral VOA, physical source still missing

The old 128 even CAR Fock states have different occupation
energies under the finite depolarizing source transfer. The standard
16-Majorana CHIRAL CFT instead assigns all 128 positive-chirality
RAMOND ground spinor fields conformal weight h=16/16=1.
Together with the 120 Neveu–Schwarz SO16 bilinear currents they
give a genuine standard abstract E8_1 weight-one current candidate.
An independently constructed D8+ even unimodular root lattice
reproduces all 240 roots, all 57,600 integral pairings and all
13,440 ordered dot-minus-one sums, with a simple-root Gram
determinant one and native (60,64,60,64) quarter phases.

Run:

    python3 -B research/originlab/main.py --out /tmp/f1-ramond \
      --campaign f1_ramond_e8_grade_bridge

This is a **mathematically sound repair of the weight-one algebraic
candidate** once a local chiral fermion CFT, spin structure and
GSO/simple-current projection are supplied. It is NOT the result of
second-quantizing the finite M3 CP transfer, and has NOT derived a
local positive original TFPT Q0, chosen P1 spin structure, the actual
R113 time, an interacting chiral 3+1D Standard Model or a graviton.
See [the exact full source and Ramond gate](F1-PREGEOMETRIC-SOURCE-CLOSURE-20261115.md)
and [five new tests](test_f1_ramond_e8_grade_bridge.py).
All six new F1 source modules have **38 locally passed regression
tests**, required in the core workflow. F1 OPEN.

## One internal positive source -> 8 CAR modes -> 16 Majoranas -> candidate E8, exact affine grade blocker

New fifth exact source-first campaign:

    python3 -B research/originlab/main.py --out /tmp/f1-source-clifford \
      --campaign f1_single_source_clifford_bridge

It starts with the SAME positive qutrit quantum algebra M3(C), not
eight arbitrarily added P1 fermionic sources. Its canonical traceless
Hilbert-Schmidt subspace V8 has complex dimension 8 and quarter-phase
multiset (2,2,2,2). The exterior/CAR functor gives eight complex
fermions, 16 Majorana generators, Fock dimension 256 and even
spinor dimension 128. All 49,152 CAR anticommutator samples have
been checked exactly, and the standard so16(120)+spinor(128)
algebraic E8 current quarter grades (60,64,60,64) follow.
This does NOT prove native physical fermionization, E8 level-one
Ward operators or a 3+1D parent.

Crucial negative result: second-quantizing the SAME source
transfer q I8 yields Gamma(q I8), with trace (1+q)^8,
supertrace (1-q)^8 and vacuum amplitude 1, as in the
original R113 exterior-empty cap contract. The even Fock
states have occupation energies 0,2,4,6,8, NOT the common
conformal weight one of E8_1 spinor root currents.
No nonzero affine rescaling of the finite number generator
remedies this. A new physical local chiral OPE / conformal
energy grading is essential. This source-first bridge
closes a previous algebraic *counting* gap but exposes the
remaining physical F1 fermionization/current-locality gap.

See [joint source theorem report](F1-PREGEOMETRIC-SOURCE-CLOSURE-20261115.md),
[exact Clifford candidate and no-go](f1_single_source_clifford_bridge.py)
and [its six regressions](test_f1_single_source_clifford_bridge.py).
Together the five new source campaigns have **33 locally passed
finite tests**, now required by the main Core Verification workflow.
No full physical F1 solution is claimed.

## F1 source from one algebra: CP, OS, exact sewing, hidden memory and conditional uniqueness

A source-first **constructive mathematical baseline** now exists. It is
NOT the original weighted TFPT source, nor the original R113 packet,
and must never be used to promote physical F1 or Atlas claims.

Run the four new finite exact campaigns:

    python3 -B research/originlab/main.py --out /tmp/f1-source-rule \
      --campaign f1_functorial_source \
      --campaign f1_covariant_quantum_semigroup \
      --campaign f1_hidden_history_dilation \
      --campaign f1_spin2_stieltjes_gate

The quantum algebra M3(C) with the explicit 10-Kraus channel
Phi_q(A)=qA+(1-q)Tr(A)I/3 has full complete positivity,
OS-positive source Gram, stationary trace state, exact time sewing,
and exact weighted coarse-graining for all 9 boundary pairs. Its
finite diagonal Markov history is reflection positive on ALL
half-cylinder observables. Operational dimensionless time is
tau=-log(q), additive under sewing. However **the source choice
q and physical clock scale are not derived**. Keeping the
same bare q at a two-step RG fixed point forces q=0 or 1,
neither of which has nontrivial ordinary decay dynamics.

Under the ADDITIONAL assumption of invariance under ALL inner
SU3 basis changes, solving 648 exact commutant equations
gives rank 79 out of 81 unknowns: a unital TP channel has
one parameter q, with exact Choi CP interval [-1/8,1].
Continuity and the Markov law then fix the semigroup up to
clock normalization as q(t)=exp(-kappa*t). The original TFPT
P1 source has not established full SU3 inner covariance, a
positive 3-state Hilbert or source star, and this finite
model has UV spectral dimension ZERO, not 3+1D.

A separate four-mark D4 full-history source provides an exact
counterexample to inferring a Markov/Lindblad generator from
an observed projection: global OS positivity, 16 exact
weighted fine-to-coarse endpoint sums and D4 symmetry all
hold, while R2-R1² has nonzero entries 1/600 and 1/1800
after projecting one mark against the other three. This is a
possible *hidden-memory mechanism*, NOT an embedding of the
native 20-state TFPT pulse.

An independently checked original 4.5 graviton source error
remains: positive Stieltjes F(z)=1+1/(z+1) yields an
additional TT pole at z=-2 of positive residue. Thus
zero-free off the physical cut does not imply massless-only
TT poles. No ghost is inferred from this example.

All statements and the physical source requirements are in
[F1-PREGEOMETRIC-SOURCE-CLOSURE-20261115.md](F1-PREGEOMETRIC-SOURCE-CLOSURE-20261115.md).
The four modules and **27 locally executed new regression tests**
are in main and in mandatory core CI. The actual original
domain retains whole_source_bundle=null and a spatial tangential
circle that is NOT physical OS time. **F1 remains OPEN.**

## Source-type reconciliation: M3 vs CAR8, E8 colour-E6 and the clock

The complete [new origin source gate](F1-QUTRIT-CAR-E8-COLOR-SOURCE-GATE-20261115.md)
separates the original qutrit M3 **operator multiplication** from the
eight-CAR one-particle Fock functor used in the newest single-source
Clifford and Ramond E8 studies. A unital M3 algebra cannot act on
the irreducible 256-dimensional CAR8 Fock; a simultaneous full
unital carrier requires at least C3 tensor exterior(C8), dim 768.
This is conditional, not a new TFPT particle count.

The exact native E8 colour roots pulled backwards through the
31-reflection source-to-native quarter clock consist of ONE D8
pair root and TWO spinor roots in their positive A2 basis. The
E6 centralizer has 72 roots, and the other 162 root currents
split into six 27-fold colour weights. ALL mixed currents carry
root conformal weight split h_SU3=1/3 plus h_E6=2/3, with
c_SU3_1=2, c_E6_1=6 and c_E8_1=8.

A 3-colour fundamental plus 5 singlet complex fermions could
produce SU3 level one. EIGHT adjoint complex fermions would
produce level six, and the natural geometric F3 triplet has
balanced C4 currents (2,2,2,2) instead of pinned native
colour's (8,0,0,0). More fundamentally the 16-Majorana
geometric/native vector C4 eigenspectra differ (4,4,4,4)
versus (10,0,6,0), excluding ANY SO16-only clock rotation
even though the full E8 current character agrees. A true
spinor/Weyl current transform and source-derived GSO are needed.

A finite 768 carrier alone also cannot represent affine E8 at
level one (nonzero central current commutator has zero finite
trace). The original full96 spatial Fourier return provides an
unbounded mode tower, not a physical original OS time or
current algebra. Fixed q at arbitrarily small microscopic step
collapses to instantaneous depolarization; physical continuum
time requires scale-dependent q(a) approaching one.

Run the exact finite bridge tests:

    python3 -B research/originlab/main.py --out /tmp/f1-color-bridge \
      --campaign f1_qutrit_car_e8_bridge

The campaign checks all 240 roots, 49,152 CAR basis
anticommutators, 55,296 qutrit/CAR tensor actions and the exact
continuum scaling. It does NOT create physical Q0 or R113,
3+1D chirality or universal gravity. F1 is OPEN.

## First priority: ORIGINAL TFPT UV measure-theorem consistency (15 November)

Run the new eight exact original-source theorem guards before interpreting
any synthetic OriginLab search candidate as physical F1 evidence:

    python3 -B research/originlab/main.py --out /tmp/f1-source-integrity \
      --campaign f1_source_theorem_audit

The ORIGINAL TFPT4.5 metrology UV fixed-point theorem as written uses
finite positive measures on unbounded [0,infinity), Lebesgue dt and
a uniform positive K(s,t)>=m>0. This makes every nonzero output
INFINITE in mass. Its published uniform-kernel Birkhoff factor is
also false: K=[[2,1],[1,2]] has true 1/3 contraction versus the
claimed 3-2sqrt(2). Additional exact controls demonstrate that
compatible compressed Dirac operators and positive completely
monotone spectral profiles do NOT guarantee projective consistency
of the full weighted measure; a pointwise positive reflection-invariant
spectral weight can be OS INDEFINITE; and four matched positive
Bernstein moments do not determine one UV profile.

A constructive integrable positive kernel repaired mathematically
does possess a unique fixed Bernstein measure, but changing its
microscopic kernel changes its profile. The actual original
source kernel must be PHYSICALLY derived, its normalization and
branch-moment invariance proved, and the original source-specific
OS and weighted fine-to-coarse pushforward checked independently.
The absent Q0 and four R113 forms cannot be fabricated.

[Original theorem correction and source-first proof requirements](F1-ORIGINAL-MEASURE-THEOREM-REPAIR-20261115.md),
[exact audit campaign](f1_source_theorem_audit.py),
[independent regression](test_f1_source_theorem_audit.py).
The candidate does not prove F1; do NOT promote the prediction ledger.

## Source candidate: original four-mark elliptic double cover spectral determinant

The **actual v512 four-point invariant** lifts to a rectangular
unit-area elliptic double cover with x=t²/(1+t²) and
y=K(1-x)/K(x). A mathematically rigorous published theorem
(M. Faulhuber, Mathematische Zeitschrift 2021,
doi:10.1007/s00209-020-02507-7) states that the regularized
scalar Laplacian determinant D(y)=y*abs(eta(i*y))^4 has a
UNIQUE MAXIMUM at y=1 among unit-area rectangular flat tori.
The exact positive eta-product derivative gives strict
log-concavity. This is a **conditional, source-adjacent**
square selector, unlike the bare chiral E8 character whose
absolute value has a minimum at the square along this family.

Original TFPT source *statistics* matter: a free bosonic path
integral supplies D^(-1/2) and prefers NOT to maximize the
square. A positive net fermion/ghost/other determinant power
D^m could select it, but the actual original all-field BRST,
Pfaffian and eta factors and their metric/area prescriptions
are uncomputed.

If the native all-field source were actually proved to
factor as |chi_E8|²*D^m (a new, unproven source premise), the
exact square log Hessian along y is 19.1570717979704
minus 1.29821132491543*m. Therefore a **local** square
maximum needs m>14.7565126188; the first INTEGER m is 15.
Neither this m nor even the factorization is an actual
physical TFPT value. A 2D elliptic double-cover determinant
is different from the original free 1D circle arc-gluing
determinant, so their opposite extremal properties are not
contradictory.

Run:

    python3 -B research/originlab/main.py --out /tmp/elliptic-source \
      --campaign p1_elliptic_determinant_source

See [source proof and original F1 gates](P1-ELLIPTIC-DETERMINANT-SOURCE-20261115.md)
and [exact regression](test_p1_elliptic_determinant_source.py).
Physical P1 Q0, one-sided map, R113 4 values, true time
and 3+1D universal gravity remain unsolved.

## F1 P1 modulus: native E8 chiral saddle and actual source-stiffness gate (15 November)

Three new **source-adjacent exact** OriginLab campaigns eliminate a
large class of false P1 solution claims:

    python3 -B research/originlab/main.py --out /tmp/originlab-modular \
      --campaign p1_original_j_witness --campaign p1_modular_source_no_go \
      --campaign p1_modulus_score

The original v512 elliptic invariant at t=tan(delta/2) is exactly
j(t)=256*(1+t²+t⁴)^3/(t⁴*(1+t²)²). Its log derivative is strictly
NEGATIVE throughout 0<t<1; j(1)=1728 and the original nonsquare
v512 witness t=1/2 yields j=148176/25. Since the standard
chiral E8_1 vacuum character is chi=E4/eta^8=j^(1/3), **the
square is a global MINIMUM of the bare chiral E8 character modulus
along this marked family**, not a maximum. At tau=i the character
is exactly 12 and its two-dimensional modulus Hessian has opposite
signs, making i a SADDLE. A purely chiral maximum-selector is
mathematically excluded.

Two explicitly constructed, time-sewing-consistent OS-positive
V4/reflection symmetric source comparison kernels have the same
time generator but select OPPOSITE modulus extrema. Thus those
generic source axioms do not determine the square. The exact
fixed-domain history identity
d²log Z/dx²=Var(S')-E(S'') at x=log t shows what the
actual weighted source must produce: net geometric stiffness must
beat fluctuation variance. In TFPT the physical moving-domain
basis measure, P_adm derivatives, eta and absolute reference also
have to be included before any positive/negative modulus claim.
A conditional factorization with |chi_E8|² would require rest
vertical log-curvature less than -19.1570718 at tau=i.

These are **not** evaluations of the original Q0, two blocked
fiber levels, same-reference R113 Z0/ZU/ZY/ZUY or 3+1D
parent. See [full source-grounded result](P1-MODULAR-E8-SOURCE-SADDLE-20261115.md),
[p1_original_j_witness.py](p1_original_j_witness.py),
[p1_modular_source_no_go.py](p1_modular_source_no_go.py),
[p1_modulus_score.py](p1_modulus_score.py) and their regression
tests. They remain in finite controls, NOT a promoted F1 theorem.

## Important correction: regional SL2 transports the marks; it is NOT an internal fixed-mode unitary

Run these two exact finite campaigns:

    python3 -B research/originlab/main.py --out /tmp/p1-moving-source \
      --campaign p1_moving_divisor_origin --campaign p1_operator_closure

The original four parabolic vector fields span the real noncompact
Möbius sl2. On a fixed four-puncture divisor its naive action
creates DOUBLE poles (at z=1 already -1/4 for the derivative of f0),
and no nonzero infinitesimal PSL2 vector field fixes four distinct
marks. Correct regional transport MUST move the divisor D(t) itself.
The exact transported meromorphic five and log-form three modes
give the standard Sym4+Sym2 complex sl2 representations. Their
canonical 3->5 evaluation inclusion is D4/quarter compatible but
NOT full sl2 compatible, and the f2 orientation compensator is
NOT a globally sl2 invariant line. Every regional Möbius motion
preserves the cross ratio, so regional symmetry alone cannot
select the original square.

The noncompact real SL2 has NO nontrivial finite-dimensional
positive unitary representation. A positive compact SU2 metric
can be chosen on the complexified polynomial modes, but the
physical source must derive the reflection, Wick/real form and
actual moving-field map. **The earlier interpretation of the
abstract positive principal SU3 Lie closure as a direct physical
consequence of the original real regional sl2 is too strong.**

Two independent conditional n=4 selection theorems are now exact:

1. If the log-form module F_n=Sym^(n-2)(C²) carries a NONZERO
   sl2-equivariant internal Lie bracket satisfying Jacobi, then
   it is perfect and semisimple, its nontrivial sl2 derivations
   are inner, and irreducibility forces dim F_n=3, n=4.
   For n=8 the candidate binary sextic third transvectant
   explicitly FAILS Jacobi with residual 5,184,000 x^6.
2. If one complete unital quantum operator algebra End(F_n)
   is built from precisely the scalar plus both original
   geometric sectors F_n and C_n, then
   (n-1)^2=1+(n-1)+(n+1), uniquely n=4. At this n,
   End(C³)=spin0+spin1+spin2, dimensions 1+3+5,
   with a positive full M3 associative algebra. Its
   traceless compact bracket is su3.

Crucial source-type exclusion: interpreting the eight geometric
modenames as **eight independent complex chiral adjoint SU3
fermions** gives affine SU3 index/level SIX, not the level-one
regular root-colour subalgebra inside native E8. The same
dimension 8 does not mean the same physical current.

See [the source-first mathematical proof and falsifiers](P1-MOVING-DIVISOR-PRIMITIVE-CLOSURE-20261115.md),
[p1_moving_divisor_origin.py](p1_moving_divisor_origin.py)
and [p1_operator_closure.py](p1_operator_closure.py).
Ten targeted regression tests were independently executed locally
and are now required by GitHub Core Verification. No original
weighted Q0, actual R113 values, physical 3+1D or graviton is
claimed. **F1 remains open.**

## Independent principal SU3 route from the original regional sl2 and 5+3 source

Two new finite/named OriginLab campaigns check a DIFFERENT fundamental
use of the SAME 5+3 original geometric modes:

    python3 -B research/originlab/main.py --out /tmp/originlab-principal-su3 \
      --campaign p1_principal_su3 --campaign p1_principal_su3_ward

The original Möbius regional span is sl2, which has spin 2
(five meromorphic modes) and spin 1 (three log-differentials)
representations. The compact SU3 adjoint under its principal SO3
decomposes as 5+3. A concrete 8x8 Lie model verifies all 28
commutators, positive compact metric and 56 Jacobi generator triples.
The original D4 representation differs by an order-two character
chi(r)=chi(s)=-1, realized by the ORIGINAL middle section
f2=z²/(z⁴-1). As a result the unmarked SU3 cubic is D4-odd,
while the formal fourth f2-mode insertion makes it D4-even.
All 112 cubic/generator comparisons are exact source-side necessary
controls; none supplies an original physical R113 mark.

Continuous SO3 covariance gives a one-dimensional space of
alternating 5x5->3 brackets. A positive ad-invariant source Gram
with norm coefficients a,b forces its coupling lambda=b/a>0
(512 checked Ward triples), and hence the compact SU3 real form
within this declared class. The native weighted Q0 has not
yet produced a,b or lambda and may not implement the assumed
regional SO3 symmetry on the original marked physical fields.

The full pinned E8 A2 root census contains 1120 regular SU3
root subalgebras, 480 with the balanced geometric quarter profile.
NONE of these 480 commutes simultaneously with native colour, weak
SU2 and hypercharge; only 10 regular A2 copies commute with all
three (and have other quarter profiles). Thus the geometric eight
modes cannot simply be renamed an independent gauge SU3 family
group with all original Standard Model charges preserved in the
regular E8 source class. A 5+3 gauge-current adjoint (weight one)
and eight independent chiral complex fermions (weight 1/2)
are competing physical assignments that dimension counting or
D4 character counts CANNOT discriminate; the native source OPE,
cubic response and Ward Gram must decide.

[Proof and source kill tests](P1-PRINCIPAL-SU3-ORIGIN-20261115.md),
[Lie/source census](p1_principal_su3.py) and
[positive Ward check](p1_principal_su3_ward.py),
plus unit tests, are in main and core CI. No original Q0/R113,
physical time or chiral interacting 3+1D/spin-2 parent follows.
F1 remains open.

## Primitive number of seam marks from first local spinor current extension

The additional OriginLab campaign p1_minimal_spinor_marks studies **all**
n>=2 distinct points of P1, without assuming four marks or E8. Its
Riemann–Roch inputs have dimensions (n+1,n-1), i.e. 2n **typed**
complex geometric modes. If the original physical source realizes
them as 2n **independent chiral complex CAR fields**, that yields
4n real Majoranas with c=2n, and a Spin(4n) halfspinor primary
with conformal weight h=n/4. A local *bosonic* spinor lattice
extension exists when n is divisible by four. Requiring a
**new weight-one spinor current** uniquely fixes n=4, hence 16
Majoranas, c=8, and the standard index-two D8+ lattice with 112
original D8 roots and 128 weight-one spinor roots (E8=240 roots,
248 currents including Cartan). For n=8 a local bosonic spinor
extension exists only at h=2 and adds no new affine weight-one
root currents. The **weight-one primitivity criterion** is essential;
bosonic locality alone does not uniquely select n=4.

Run the exact finite search:

    python3 -B research/originlab/main.py --out /tmp/originlab-four-origin \
      --campaign p1_minimal_spinor_marks

The argument is conditional on the **unproved physical CAR source
field map and the decision to require primitive weight-one spinor
currents**; it does not select the special square position of marks,
original R113 values or a local 3+1D physical TOE. A complete
proof and independent regression are in
[P1-MINIMAL-SPINOR-MARKS-20261115.md](P1-MINIMAL-SPINOR-MARKS-20261115.md)
and [test_p1_minimal_spinor_marks.py](test_p1_minimal_spinor_marks.py).

## Geometric 5+3 native E8 Weyl proof and one-sided seam source gates

The most important new finite **algebraic** F1 progress of 15 November is
the exact 31-reflection E8 Weyl conjugacy between the geometric quarter
action on the original Riemann–Roch section/differential sum (5+3)
and the repository's native D5+A3 glue C4 action. The eight typed
geometric modes have Z4 multiplicities 2,2,2,2. **If** they can be
identified with independent positive chiral CAR modes in the original
source, their sixteen real Majorana components generate a 120-current
so16 bilinear sector and a 128-state positive spinor extension, with
total exact affine E8 current phase counts 60,64,60,64. The 52 fixed
roots form D5+A3. The explicit Weyl word matches every one of the
240 root charges, not merely their multiset. It replaces the old
arbitrary 256-dimensional spectral cover with an exact standard
248-current Lie-algebra construction **conditional on the CAR/VOA source**.

Run the new focused OriginLab campaigns:

    python3 -B research/originlab/main.py --out /tmp/originlab-p1-e8 \
      --campaign p1_four_current_gate --campaign p1_e8_geometric_lift \
      --campaign p1_seam_dirac --campaign p1_geometric_gram

The genuine four-nonzero-root affine current channels have exact same
frame ratio t^(-4), which equals exp(2pi*Delta_odd) from original
v512. Ordinary current crossing swaps both operators and points,
so it does NOT force fixed-frame equality or select the square.

The exact 5+3 inclusion also defines a finite self-adjoint odd
8-dimensional seam operator with 3 positive, 3 negative and 2
endpoint zero modes. Selecting ONE endpoint gives a positive
7-dimensional transfer with unique vacuum and gap, exchanged with
its partner by geometric reflection. This is NOT the original
physical OS time or 3+1D Dirac collar; the endpoint cap is unproved.

A sharper source-level Gram test shows that the original D4
representation allows only one complex mixed coefficient z between
the five and three typed sectors after canonical marginal
normalization. The positive combined Gram has determinant
(1-|z|^2)^2: rank 8 for |z|<1, rank 6 at |z|=1 and is indefinite
above one. A naive rank-5 identification of all three logarithmic
fields with three of the five section fields violates these D4
symmetry and positive-marginal premises. **The physical original
Q0 must provide z and the actual CAR field map**; the abstract
geometric theorem does not.

The source tests and complete caveats are recorded in
[P1-FIVE-PLUS-THREE-WEYL-E8-20261115.md](P1-FIVE-PLUS-THREE-WEYL-E8-20261115.md)
and [P1-SEAM-DIRAC-GRAM-20261115.md](P1-SEAM-DIRAC-GRAM-20261115.md).
The finite tests are in test_p1_four_current_gate.py,
test_p1_e8_geometric_lift.py, test_p1_seam_dirac.py and
test_p1_geometric_gram.py, all registered in the CI workflow.
Native original weighted Q0/R113, physical locality, time,
the interacting chiral Standard Model and universal spin-2 gravity
are NOT derived. No F1 or ledger promotion follows.

## Charged E8 root triple and determinant source test

The exact P1 source search also provides the named and finite campaign
p1_charged_modulus:

    python3 -B research/originlab/main.py --out /tmp/originlab-e8-p1 --campaign p1_charged_modulus

It exhausts all 13,440 ordered pairs of E8 roots with inner product -1
and verifies that their **single free neutral chiral vertex triangle**
has amplitude modulus proportional to 1/sin(delta). Therefore it has
a minimum, not a maximum, at the original square. The algebraic
no-interior-maximum theorem extends to **any single neutral**
positive-metric free-boson four-mark vertex word. This is not a
claim about a sum of channels or the original interacting source.

A **conditional** original-source candidate multiplies that free
E8 triple by BOTH scalar determinant lines of the proposed relative
boundary polarization, producing an amplitude proportional to
sin(delta), uniquely maximized at the square. The required net
physical determinant exponent must be greater than 1/2; at 1/2 it
only gives a flat modulus. Exact withheld model ratios are
3/5, 4/5 and 12/13 at t=tan(delta/2)=1/3, 1/2 and 2/3,
respectively. These are explicitly **synthetic**, not R113 source
values. The original all-field determinant/Pfaffian/ghost/eta factors,
charge and source map remain missing. Full proof obligations are in
[P1-CHARGED-RESIDUAL-ORIGIN-20261115.md](P1-CHARGED-RESIDUAL-ORIGIN-20261115.md),
with finite regressions in [test_p1_charged_modulus.py](test_p1_charged_modulus.py).
A further exact finite two-face CAR model packages the two
determinants into one 8-Majorana antisymmetric Pfaffian. The four
original source-compatible interface directions/signs have different
determinants, and only the unit forward/forward choice cancels the
free charged root pole to a square-maximizing residual. This is **NOT**
a physical TFPT Pfaffian evaluation, and scalar positivity does not
establish OS reflection positivity.

The original [9 November R113 contract](../../experiments/theory-contracts/chat-triangulation-20261108/R113-ORIGINAL-UND-REFERENZ-20261109.txt)
provides a decisive additional source guard: on its actual exterior
empty state Gamma(R) has vacuum expectation 1. The full Fock trace
Tr Gamma(R)=det(I+R) is a **different** contraction. Consequently the
candidate paired determinant cannot be inserted into the native
zero-Fock reader unless it is *derived from the original scalar
s[b]* or the complete physical selector/cap, with the same time and
reference. The campaign checks this distinction explicitly.

These controls do not solve F1.

## Primitive P1 source falsification gate (15 November 2026)

Run `python3 -B research/originlab/main.py --campaign p1_source_falsifiers
--out /tmp/originlab-original-p1` to audit three proposed origin links
**before** any physical F1 claim. The finite campaign proves that the
3+2 Riemann–Roch evaluation sequence, the logarithmic orientation and
unit winding occur for every noncolliding antipodal four-mark modulus;
the extra geometric quarter picks the square and remains an input.

The same controls show that the apparent square preference from the
**vertex-only** Gaussian determinant cancels after integrating all four
interior arcs of the **unchanged free circle**. The exact weighted
cycle identity is `det_prime(L)*prod(2*ell_i)=64*perimeter`.
They also give a product-Dirac Cauchy-data countermodel to the automatic
identification of the Calderon projector with a local geometric deck
reflection. A two-involution **candidate**, not original source
reconstruction, relates the cross ratio to the trace of their relative
unitary: `lambda=(tr(U_rel)-2)/(tr(U_rel)+2)`, with square precisely
at trace zero in the stated two-dimensional family.

These are source-constraint and negative-control results. They supply
**none** of the original weighted `Q0,Z0,ZU,ZY,ZUY`, the physical
time Hamiltonian, full fermion/ghost factors, 3+1D locality or gravity.
See [full source audit](P1-PRIMITIVE-GLUING-SOURCE-GATE-20261115.md)
and [exact test suite](test_p1_source_falsifiers.py).
`--finite-controls-only` includes the campaign.

## R113 inputs

Without a physical packet, the native diagnostic domain produces a structured
BLOCKED result and leaves all four source marks null. A supplied packet must
bind to a native domain path and SHA256 through `original_domain_binding`.
The domain may contain the packet inline or point to a JSON file. Original
input and packet references are hash checked. Exact Gaussian rational scalars
are pairs of integers or rational strings. Two supported representations are
complete endpoint forms and a complete finite atomic measure with all four
forms per atom. The evaluator uses a null quotient, applies one outer
admissibility projection, contracts with the absolute cap, and divides by the
bound empty reference. It does not supply the physical producer or prove
continuum limits.

For supplied nested atomic levels, the optional `coarse_graining` record
checks each marked weighted fibre after the declared word embedding and
absolute reference normalization. See [the packet and native clock scope](P1-FIBRE-AND-CLOCK.md).

Run the source synthesis and inspect its decision alongside F1 with:

```sh
python3 -B research/originlab/main.py --out research/originlab/runs/r113-synthesis \
  --campaign r113 --campaign r113_synthesis --campaign f1
```

The native campaign replays the original centre-grading checker and
checks the charged word against the reviewed generators. It reports the
first required sector-changing matrix block before any outer selector.
The original affine vacuum word and the weighted R113 values have
separate fields in the report.

Run `--campaign source_bridge` to replay the absolute factor on one original
R032 history, test the full span of the native regional generators, and check
the neutral rank required by a charged affine Ward return. See
[source and session results](SOURCE-BRIDGE-2026-11-14.md) for the derivation,
scope and independent checks from the parallel sessions.

Run `--campaign marked_search` to replay the 14 November four-mark report.
The campaign enumerates both phase refinements, quarter clocks and the
four orthogonal symmetric cell sewings. It checks the integral A4 glue,
the meromorphic clock map and accessible mixed cell correlations. Pass
`--marked-search-input FILE` to apply observations and recompute exact
survivors and worst-case probe counts. Input keys `marks.observations.phase`,
`marks.observations.clock`, and `cells.observations` take mappings from probe
IDs to exact observed values. Use the point order printed in the report.
These inputs have finite control scope and do not authorize source binding.

`--campaign composed_search` joins libraries of operator fragments on one
Gram carrier and one trace contraction. Provide a case with
`--composed-search-input FILE`, or run the default D4 transfer, generated
quarter action, and odd field control. A finite case can give one `gram`
and `trace_weight`, or `boundary_choices` with alternative metrics and
trace weights in the same ordered coordinates. Each slot implements disjoint
`letters` using `choices` of exact matrices, or generates one letter
from a bounded `grammar` of OriginDSL endomorphisms. Shared `bindings`,
adjoints, `positive_steps`, and `identities` constrain the combination
as soon as their inputs are available. For example, a commutator is an
identity with two coefficient/word terms. `fit` and `withheld` observations
must name distinct ordered words. All-word classification uses an exact
reachable-space test and emits a separating word for inequivalent
candidates. An exhausted comparison budget prevents any uniqueness claim.
Use `--composed-square-resolution 64` with the built-in case to cover
33,800 combinations of exact transfer, quarter and field operators before
pruning. The default grid denominator is 8.

The original `whole_source_bundle` can supply `composed_search` with
its `source_id`, `time_id`, `reference_binding`, all four `mark_words`,
and a `case` of operator families on the native Gram quotient. Optional
`withheld_forms` carry a word and its independently supplied native form.
The evaluator fixes the Gram, cap, selector, time and absolute reference,
fits the first three R113 marks, and withholds the mixed mark and extra
forms. Candidate metadata cannot replace these observations. A result
certifies only the declared process class on its supplied alphabet;
continuum, chirality, universal gravity coupling and independent physical
validation remain separate gates.

Use `--campaign r113_synthesis --r113-synthesis-input FILE` with
`input_role: declared_finite_oracle` and `kind: inverse_fibres` to solve
for the nonnegative atomic masses of a candidate producer. Each entry
in `snapshots` supplies `fine_size`, `fine_atoms` without measure masses,
the fine scalar and reference, `time` (`id`, `t_U`, `t_Y`), and the
existing `coarse_graining` record with all four independently supplied
forms. The same atom identities, measure bindings, and reference
binding are used across snapshots. The default fit uses Z0, ZU and ZY
per coarse configuration. Optional `independent_probes` give an ID,
coefficients keyed by atom, and an expected value for a withheld linear
contrast, such as an independently derived weight ratio. An optional
`total_mass` must come from the same absolute convention.

The rational solver emits an independently replayed dual for exclusions,
maximal support and its affine dimension for feasible masses, and exact
extrema with certificates for the first varying withheld component. The
second feasibility check uses all withheld components jointly. It emits
candidate masses only after the existing forward fibre verifier accepts
every transported form separately. Proportional atom classes describe
indistinguishability on the supplied forms. The raw candidate fine Z0
and its constancy on the feasible mass set are reported before outer
contraction. Source identity still requires
proof that the atom family comes from the original operator, a bound
measure and reference, and control of regulator removal.

The native synthesis also classifies the finite Higgs histories in a
fixed flat transport chart under constant doublet U(2) and spatial
translations, with an additional table allowing reversal of the circle.
Equality of Hermitian Fourier Grams after a common translation phase
character is necessary and sufficient for that restricted equivalence.
The probe search compares reflection-even quartic Gram words with equal
mode differences, so it cannot split a class by a chart phase. The
original full96 Majorana and neutrino Higgs blocks independently replay
one sign-sensitive word and its transported chart control. This is a
spatial probe. A local gauge interpretation requires a transported
connection, and its physical weighted response requires the native
source and absolute reference.

The reflection invariant Gram measurements span the finite function
algebra on all 11 native and 22 refined spatial orbits. The reported
plans use 11 and 22 independent moments, respectively, including the
empty word. Each uses at most two probe factors evaluated in the same
fibre before integration. With the finite input kind
`orbit_source_moments`, supply `field_proposal`, `word_space_size`, and
`measurements` in the exact plan order, each with `word` and all four
`ordered_source_forms`. The inverse works entrywise on the four full
matrix forms and replays every supplied moment.

The plan also includes a factored nonnegative polynomial on selected
orbit coordinates. For any positive measure with a finite integral of
this polynomial, an integral of zero proves that the pushforward is
supported on the listed signatures. A positive integral refutes that
finite support claim. The support check concerns the chosen observables.
Microscopic distinctions within their fibres and the absolute reference
still require original data.
The test also checks constant probes such as the squared zero-mode norm,
which equals one for every supplied history. The bounded test `P/(1+P)`
has the same zeros and can be integrated against a finite positive
measure even when polynomial moments diverge. On the audited constant
radial Higgs sector the zero-mode term is `(R²−1)²`. The Calderón
volume candidates with a positive radial `R^-2 dR` tail at a fixed cutoff
therefore fail finite support on these proposals, and the polynomial
test diverges on that tail. This exclusion uses the audited replacement
volumes and does not assign a measure to the original physical source.

A factor depending only on the seam angle, positive at the tested angle,
preserves that radial obstruction. If a positive radial reweighting is
comparable to `R^-β` on the tail, finiteness of the zero-mode defect
requires `β > 3`. Even fast decay leaves a nonzero bounded defect if the
measure still has positive mass away from `R = 1` in this sector.
In the joint removal of the spatial cutoff and Gaussian regulator, the
audited replacement volumes have vanishing normalized mass on every
compact radial interval for every removal sequence covered by the proof.
The normalized bounded defect therefore tends to one, since it is at
least `(R²−1)²/(1+(R²−1)²)` and is at most one. The compact-mass bound
is established in section 3 of `TFPTCALDERONWEIGHTLIMIT20261110`.
[Atlas source record](https://www.fixpoint-theory.com/results?record=finding%3Af1-native-orbit-support-calderon-radial-exclusion-20261115).

The native R113 synthesis now also reports `native_higgs_moment_plan`.
It computes the complex Fourier coefficients of the gauge-invariant
Higgs norm in exact arithmetic on the declared native spatial proposals
and the independently audited eighth-harmonic refinement. Products of
these coefficients generate the finite evaluation algebra. The plan
provides a full-rank set of words and an exact inverse from ordered
integrated moments to class totals. On the native 27 histories the norm
algebra sees 8 classes. On the audited 81-history refinement it sees
20. For one fixed nonzero first harmonic, the constant and Fourier
harmonics 7 and 16 separate the three signs of the eighth harmonic.
These spatial observables require the original weighted fibre before
they become source moments. If fine Z0 still varies after all supplied
constraints, the inverse solver reports a separating raw matrix entry
with exact attainable bounds and dual certificates.

The same campaign checks a tetrahedral transfer. Supply `transfer.moments`
as three or four exact 4 by 4 forms on the same ordered mark vectors,
at steps zero through two or three of the same positive contraction.
Use the Gaussian rational encoding described above. The diagonal and edge
entries give the singlet and triplet moments as `d+3e` and `d-e`. For moments
`(a,b,c)` in either sector, the sharp third moment interval is
`[c²/b, c − (b−c)²/(a−b)]` when `0 < b < a`. Cauchy Schwarz with the
measures weighted by `x` and `1−x` gives these bounds. Each endpoint is
attained by a two-point measure, and their mixtures fill the interval.
The report includes those spectra and the minimum additional dimension
required by the first three moments. A positive variance in the triplet
requires three additional directions. The default control requests the
third triplet return and yields `[25/128, 31/128]`.

Set `transfer.conserved_distinct_weights` only when the four vectors carry
distinct joint Cartan weights and the transfer conserves those charges.
In that case all mixed entries must vanish and tetrahedral covariance makes
the four diagonal entries equal. The bounds allow transfer eigenvalue zero.
The form check does not supply a native word map, the physical time step,
the absolute reference, or the R113 values.

Set `transfer.symmetry` to `D4` for the square in cyclic mark order, or
`V4` for its antipodal deformation. `A4` is the tetrahedral control. The
optional `harmonic_seam` input runs the cross-ratio, transfer, and reciprocal
function searches. See the [harmonic seam search](P1-HARMONIC-SEARCH.md)
for a reproducible run and the missing physical comparison.

The audited plan has an executable benchmark campaign:

```sh
python3 -B research/originlab/main.py --out research/originlab/runs/benchmarks \
  --campaign benchmarks --lean-toolchain leanprover/lean4:v4.29.1
```

The suite runs B01 through B09 with individually scoped receipts. Its R113
baseline expects the verified missing-source outcome. A future producer
requires a new positive acceptance case. The F1 report additionally emits
the physical gate status and the upstream source blocker. `source_intertwiner`
in a supplied endpoint can test an adjoint-closed rectangular embedding
using Cartan, colour, time and charged-middle relations on the reflected
quotient. Higher-dimensional feasible spaces may remain unresolved.
To replay the execution gate with a Pawl checkout, run its `tools/run.py`
with claim `VERIFY.ORIGINLAB.BENCHMARKS.01` and `--config
research/originlab/pawl/verify.json`. The optional `ORIGINLAB_LEAN_TOOLCHAIN`
selects the same kernel check as the benchmark CLI. Pawl records execution
only. The Atlas remains the canonical authority for scientific evidence.

The report also checks the marked affine root string, its cyclic
same-colour test, and the two-time obstruction to using only the
lowest root doublet. An optional `marked_root_doublet` in an endpoint
packet accepts the two creation and adjoint forms, raising and lowering
forms, root weight, colour, hypercharge, and two positive time step
forms on the same reflected word basis. The colour form may use exact
`rational` and `sqrt3` parts. The reader reports Ward and adjoint
defect ranks, time leakage and the relative ordered returns. Absolute
values continue to use the endpoint packet's bound common reference.

`--r113-synthesis-input FILE` runs one explicit finite oracle instead.
Set `input_role` to `declared_finite_oracle` and use one of these `kind`
values:

- `real_charged_spectrum`: provide `spectrum`, `diagonal_observations`
  (one row of two values for each step starting at zero, as many rows as
  levels), and `mixed_observations` (known steps keyed by decimal integers).
  The search solves for the individual spectral weights, exhausts relative
  signs and selects the earliest separating mixed word. Positive transfer
  witnesses are replayed with the existing block Hankel reconstruction.
- `fixed_histories`: provide `histories`, a fixed positive `reference`,
  `constraints` including `Z0` (each with `coefficients` and normalized
  `value`), and `mixed_coefficients`. Coefficients are unnormalized. The
  history weights are nonnegative and sum to one. The search returns the
  exact range of the normalized mixed mark and extremizing weights, or
  reports that the declared class has no solution or exceeds the bound.
- `complex_charged_spectrum`: use the same moment fields as the real
  spectrum case, with exact complex mixed values as `[real, imaginary]`.
  Positive two-mark spectral blocks may have rank two. The search checks
  the positivity disks and constructs paired witnesses when a feasible
  seed has a certified free direction. Supply `feasible_seed` if the
  particular linear solution is outside the disks. Boundary or step
  failures can leave feasibility or uniqueness unresolved.
- `gaussian_source_jet`: provide the exact matrices `kinetic`, `U`, `Y`,
  `UY`, four scalar derivatives `scalar_weight_jet`, and `reference`.
  Both the determinant coefficients and independent replacement columns
  must agree. The mixed contact and the singular case are retained.

Oracle results are finite controls. A physical source still requires the
original field maps, ordered forms, measure, outer selector and absolute
reference on the same domain.

An optional `energy_moments` block accepts two forms `first` and `second`
for the same nonnegative shifted generator. Supply the packet's `source_id`
and `time_id` plus `generator_binding`, `ground_energy_binding`, and
`word_domain_binding`. The Gram quotient gives the exact compression
defect and the minimum number of missing directions. For nonzero defect,
the output includes two positive finite extensions with identical first and
second moments and a third-moment probe that separates them. This probe
does not determine the physical generator or the absolute reference.

## Open source components

Run `--campaign double_cell` for the exact Clifford double-cell control.
It checks the joint five-slot occupation projectors and the full permutation
action. It also checks native parity, the standard weak slot exchange, the
minimal parity reference and the joint relative-phase readout. The finite
control supplies a discriminating word for the physical source search.
It carries no original R113 value.

Run `--campaign moving_carrier` to compare repeated evolution on the fixed
regulated Hardy carrier with transport to the moved carrier. The first
comparison has a rank-four two-step defect at nonzero radius and time.
Transport between the moving carriers composes, while reproducing all fixed
carrier moments requires infinite spectral support in each nonconstant
sector. The original follow-up compares these projection choices on the
same charged word and absolute empty reference.

Run `--campaign affine_polarization` to compute the normalized connected
Cartan response of the original affine three-leg word. A single time pair
sees one outer product. Four pairs span the four matrix directions of the
two-dimensional charge plane. Its six-dimensional orthogonal complement
is an exact negative control. Restricting to polarizations orthogonal to
the middle root leaves just one direction visible to the entire word.
An optional `affine_polarization` section of the weighted endpoint packet
binds the bare and prepared word forms to the same source, time and
reference. The consumer reports the absolute bare comparison separately
from normalized Ward agreement.

Run `--campaign affine_operators` for the operator-level source. It builds
the actual mode operators of the level-one E8 lattice construction in the
fixed D8 chart and applies them to exact Fock vectors. These currents change
the charge sector and lie outside the excluded free quadratic class. The
campaign checks the level-one current relations on explicit vectors, with the
trivial cocycle as failing control. It recomputes the marked word coefficients
`min(m,n)`, the Ward response of prepared words coefficient by coefficient,
and the neutral return spectrum on the full grades zero to two. The joint
return check extends the single root line. Eight mutually orthogonal roots
give commuting returns with a rank-eight excess on neutral grade two. Two
returns commute there exactly when the roots pair to 0 or -1.

The same campaign searches for the fewest further root currents that
generate all of E8 together with the marked plane, restricted to a declared
charge class. The marked charges span one A2 with commutant E6, and their
currents generate 8 of the 248 grade-one vectors. Six further currents are
necessary by rank and sufficient in five of the eight restricted classes
and in the unrestricted control. These include colour-centre neutral D8
pair currents. Three classes are exact
negatives. Hypercharge neutral currents, and currents neutral under colour
centre and glue, generate only A1 x E7 with the marked plane. The conserved
grading is the parity of the pairing with the weak slot root `e4 - e5`. All
112 missing currents are weak doublets and all carry hypercharge. The
colour neutral D8 pair currents among them have glue class 2 and
hypercharge of modulus 1/2. The restriction applies to individual root
letters. A word can be neutral overall and still probe these currents.
The executable example splits the middle root into `e1+e4` and `e6-e4`.
Both orders preserve the end charges and give normalized mode coefficients
`1/min(m,n)` and `1/min(m,n)-1` on the tested grid. Read both orders
independently with the same weighted source, time and absolute reference.
Run `--campaign weak_return_rank` for the six late D8 weak legs. The
neutral excess return span first reaches rank six at grade two.
Products of one, two and three returns span cumulative dimensions six,
24 and 36. The phase-sensitive ordered triple distinguishes two strictly
positive densities that agree on all shorter words.

The same campaign compares those six legs with the two half-spinor
orbits of six legs each. It finds joint excess ranks 14 and 98 at neutral
grades two and three. On the 14 dimensional grade two support the spans
from ordered words of lengths up to one, two and three are 18, 108 and
196. The campaign emits a basis and an example length four holdout in
`completion_orbits`. The optional `weak_return_grade_two` packet binds
`source_id`, `time_id`, `charged_modes_and_adjoints_binding`,
`neutral_grade_and_excess_binding` and `absolute_reference_binding`.
Copy the order of `roots_doubled` and supply each emitted basis word in
`basis_words` as a `word` and whole `word_form`. A letter is the normalized
excess $(E_\alpha^*E_\alpha-I)/3$ at neutral grade two. Supply a fresh
nonzero path such as the emitted holdout in `held_out_words`, in the same
format and on the same source. The reader reconstructs the absolute linear
functional on this declared excess algebra, tests positivity, and tests
the held-out whole form. The reference and the outer selection are those
of the surrounding original packet.

Run `--campaign four_mark_geometry` for the shared four-point sphere.
The family cycle and quarter cannot both act holomorphically there.
The equianharmonic alternatives carry the family cycle holomorphically
and the quarter antiholomorphically. On the elliptic cover its square
is a half-period translation and acts trivially on the differential.
The even permutations move a reflection fixing two marks through six
distinct circles. The effective complex Möbius action also limits every
commuting real family of geometric generators to dimension two.

Run `--campaign affine_car` for the independent free D8 fermion control.
It reconstructs the pair current kernel and the weak doublet Casimir.
Its quadratic currents carry only D8 pair roots. The original marked
spinor ends and the common weighted P1 reference require their own field
map and source.

These affine vacuum coefficients leave the weighted R113 values and
physical source identity open.

For a supplied weighted endpoint, `weak_split_returns` binds `source_id`,
`time_id`, `charged_field_and_mode_binding`, `adjoint_binding` and
`absolute_reference_binding`. Supply `doubled_middle_roots` in the order
`e1+e4`, `e6-e4` and the actual `bound_cocycle_sign` (1 or −1).
Each of at least two `samples` contains mode indices `m` and `n` and
`bare_form`, `forward_form`, `reverse_form` on the common word basis.
Build each whole word before applying the single outer selector. The
reader reports the absolute reference quotients for both orders separately
from their difference and their ratios to the marked bare word.

Run `--campaign weak_source` to search the 240 roots with the middle
charge fixed and replay the surviving gauge singlet words. The report
emits a calibration rectangle and a new time pair for the original test.
The native `r113_synthesis` campaign independently checks the literal
96 dimensional Higgs blocks, their polar CAR parity and the separated
DET16 block. It also reduces the existing coupled W word, with the original
neutral pair dressing on both ends, to exact reachable time blocks of
dimensions 4 and 13. Its two interval kernel has separation rank four
for nonzero coupling. The autonomous runner schedules the affine search
and watches the native inputs for changes.

For weighted data, `weak_gauge_responses` binds `source_id`, `time_id`,
`native_three_current_binding`, `same_dressing_and_inserted_actions_binding`,
`casimir_insertions_binding`, `physical_to_radial_time_binding` and
`absolute_reference_binding`. The optional `comparison_paths` begins with
`spinor` and includes one or more of `d8`, `weak_spinor_minus` and
`weak_spinor_plus`. The default pair is `spinor` and `d8`. Each sample
supplies `left_radial_factor`, `right_radial_factor`, `held_out` and `words`.
The calibration rectangle spans at least two widths per interval, and
one held out sample uses two unseen widths. For each path supply the
`bare_form`, the four `casimir_forms` and the six `twisted_words` using
the chronological links returned by `example_networks()`. Optional
`winding_forms` use the same `original_coweight_binding` on every path.

To check a linear combination of the four fixed middle paths, add
`joint_early_insertion_binding` and the five `joint_early_forms` to the
samples in order from power zero to four. Use the diagonal weak element
with entries $(3+4i)/5,(3-4i)/5$ and the quarter insertion of the coweight
$-2e_8$ on the first interval. The reader solves for the four absolute
coefficients and reports the power four residual. These are coefficients
in the specified affine response model. A zero residual is a check on the
supplied forms. The original source and these insertions need their own
derivation to establish a physical match.

Run `--campaign weak_folded_clock` to compare the minus winding half
spinor with the marked word in both original quarter charts. Two commuting
D8 root zero modes preserve the middle current and winding. Transporting
the word co-transports the folded generator. Across all colours the six
products span the maps from the early colour dual triplet to the weak
doublet and carry $(3,2)$ with $Y=1/6$. In the literal fixed chart
the early and late spectra are separated, and the weak lower to upper
weight ratio is $\sqrt{xy}$ for the interval radial factors $x,y$.

Use `weak_clock_sectors` in a weighted endpoint packet to recover
independent complex weights for the seven eigenvalues of the same early
joint action. Bind the packet's `source_id`, `time_id` and the following
original data.

- `common_fields_and_frames_binding`
- `fixed_joint_early_action_binding`
- `same_histories_for_all_powers_binding`
- `physical_to_radial_time_binding`
- `absolute_reference_binding`

Each sample supplies `left_radial_factor`, `right_radial_factor`, a
`held_out` flag and eight entire `power_forms` ordered by powers zero to
seven. The inputs use the packet's Gaussian rational form convention.
Use a full rectangle with at least two distinct widths per interval and
one additional pair with both widths unseen. The reader solves for the
seven amplitudes at each pair using powers zero to six, checks power
seven independently, and compares the weak weights with both the fixed
fold and the W symmetry. With five widths per interval, rank five in a
resolved half spinor component rules out the fixed W word, including
separable clock factors. The campaign supplies sample widths that also
keep the literal fractional clock powers rational.

Run `--campaign p1_field_search` for the charged first variation on both
original folded charts in the specified two-face collar. The same
covariance and heat derivative can come from a charged intrinsic field
or a shared motion of the two frames. The independently anchored
hypercharge derivative separates the two contributions. The checker
reconstructs the intrinsic field and verifies all three measured matrices.

With `--p1-field-input FILE` or the optional `p1_charged_jet` field of a
weighted packet, supply `energies`, `hypercharges`, a rational
`radial_base` in $(0,1)$ and an integer `exponent_scale`. The scale
makes the energy exponents integers. Distinct hypercharge sectors must
have disjoint energies in this class. Give the whole
`two_face_covariance_tangent`, `two_face_hypercharge_tangent` and
`two_face_heat_tangent` as exact Hermitian matrices. Each matrix lists the
same energy basis on the first face followed by its copy on the second
face. The packet also binds `source_id`, `time_id`,
`two_face_frame_binding`, `independently_anchored_hypercharge_binding`,
`covariance_and_heat_tangents_binding` and `absolute_reference_binding`.
The full charge derivative includes the cross-face blocks.

Run `--campaign colour_ward` for the exact colour grading and affine
Ward relations. The campaign is part of the autonomous search and the
finite control run. It checks the 1,620 oriented cubic coefficients in
each charged sector. They belong to 270 unordered root triangles per
sector. Of 810 independent pair coefficients, the complete neutral
action fixes 809 relative coefficients. The adjoint and Jacobi identity
fix the modulus of the remaining common scale. Rephasing the two charged
sectors changes its phase by a cube. The fixed bracket admits the three
cube roots of unity as simultaneous sector rephasings.

The same report resolves the joint glue and colour grading modulo twelve
through affine grade three. It also tests the reflected form on affine
grade one before transfer and selection. Neutral Ward adjoints allow
four independent weights. The marked charged zero mode reduces them
to one common scale. See [the derivation and exact controls](P1-MARKED-AFFINE-METRIC.md).

The neutral spinor zero mode has ranks 22, 18 and 18 on the three
grade one sectors. In colour two it gives 15 transitions from spinor to
pair currents and three in the other direction. The positive neutral
mode paired with the neutral spinor returns the vacuum at level one.
Positive neutral modes annihilate the charged grade one states. The
free neutral currents and the spinor extension both give central charge
eight.

The native `r113_synthesis` comparison uses the same matrices as the
original finite diagnostic. Two successive neutral colour steps of the
affine word reach $Y=5/6$. Its native exterior comparison uses the
Majorana and Higgs matrices with an anti down spectator. Each weak
component has forward rank 27, reverse rank zero, a nonscalar Gram
matrix and rank 27 leakage of the Higgs adjoint through the spectator.
The ordinary exterior space in every degree from 1 to 95 has a nonzero
hypercharge trace on its charged colour blocks. The full exterior space
and both parity spaces also have nonzero traces. A finite $E_6$ action
with the given hypercharge inside its Cartan requires vanishing trace on
each invariant colour block. The smallest neutral completion matching
all native colour, weak and hypercharge representations has dimension
240. Its proof uses the Weyl dimension formula and joint weight counts.
The original weighted source still needs to bind the neutral currents,
the charged preparations, their time and their absolute reference.

The `weak_parity_bridge` optionally tests the source parity, the Cartan
and raising and lowering matrices for the tested root, the charged
matrix with its adjoint, and the three endpoint projectors. The source
forms are read on the reflected quotient before selection. The reader
reports the rank of the charged endpoint and separately measures the
two orderings. Both readers use the same bound time and reference as the
whole endpoint record.

- [SymPy](https://github.com/sympy/sympy) executes exact matrix and polynomial arithmetic.
- [Z3](https://github.com/Z3Prover/z3) provides inverse integer synthesis, with exact witness replay.
- [Lean](https://lean-lang.org/doc/reference/latest/Axioms/) checks the separate rank lemma and reports its dependencies.
- Exact truth tables enumerate the entire binary local rule class admitted by the declared observations.
- Typed word normalization verifies the declared compositional identities. [egglog](https://egraphs-good.github.io/egglog-tutorial/01-basics.html) can also execute the three checks when installed.
- [FLINT/Arb](https://python-flint.readthedocs.io/en/latest/general.html) encloses the free regulator tail and the conditional R032 kernels.

Pass `--finite-observations observations.json` to feed a measured finite
response back into search. Each item has `dimension`, a `word` of integer
tokens (0 candidate, 1 boundary reversal), an exact rational string `value`,
and `evidence_id`. Disagreement records an exact counterexample scoped to the
declared finite readout. These observations carry no automatic physical status.

## Autonomous scenario campaign

Install the pinned Open Source backends from requirements.txt. Run from the
TFPT root with the environment's Python:

```sh
python -B research/originlab/autonomy.py --out /path/to/campaign \
  --workers 6 --max-resolution 5 --lean-toolchain leanprover/lean4:v4.29.1
```

The durable SQLite queue runs independent subprocess workers. A failed job
gets one retry and its error remains in the run record. An interrupted
campaign resumes against its stored inputs and job keys. After the first
two sweeps the runner prioritizes topology and dynamics combinations that
share their measured return words. Increasing `--max-resolution` extends
the saved search horizon. Add `--watch-seconds 90` to keep the
process alive and restart the search on a change to the Atlas or R113 input.
The engine and source version form part of the run identity. The queue is
local to one supervisor and one machine.

The built-in scenario set compares interval, circle, figure eight, disk,
sphere, cylinder, Möbius band, torus, Klein bottle and projective plane.
It evaluates five mechanisms: a primitive local transfer, a boundary
resolvent, a flat sign local system, an edge history with vertex readout,
and a filled-face Hodge operator. Unit and alternating conductances are
explicitly supplied geometries. The transfer families are finite
mathematical candidates. Their topologies and parameters are not derived
from P1. A finite cell complex and these graph transfers do not imply
a physical spacetime, Hamiltonian or continuum field theory.

The Boolean lane enumerates all 64 tables compatible with its two uniform
observations and splits on both answers to the most discriminating probe.
Typed word normalization checks the declared compositional identities.
Arb certifies the
free periodic heat trace and a convergent regulator comparison, with
the interacting R113 remainder still open. The bridge lane checks the
tetrahedral dual basis, joint cone, charged quarter response and scoped
R102 inputs. It keeps common source identity as an explicit obligation.
An optional `ordered_bridge_forms` section of an endpoint packet can
carry `field_binding`, `charged_return`, and `quarter_marked_return`
on the same basis and absolute reference as the four R113 marks.

A completed round writes a revision-bound review outbox. `publish.py
--status` reads the canonical live Atlas. `publish.py --review FILE`
validates the exact outbox and referenced artifact hashes, then submits a
source-reviewed `record_research_update` request to the canonical service.
The reviewed envelope has `reviewer`, `outbox_sha256`, `outbox` and
`update` (the original Atlas tool arguments). The search machine never
marks a new mathematical or physical claim as reviewed by itself.

The stitch lane builds alternative derivation paths with conjunctive
premises. Facts from distinct source identities cannot be combined. Its
current frontier identifies two possible routes to the joint forward
spectrum and the original source obligations shared by both routes.

## Reconstruction and open scenario space

The reconstruction lane enumerates Pauli readouts by support, computes the
state dependent Gram rank and constructs exact finite representations from
coherent word forms. The null space must be respected by every submitted
operator. The identity is normalized, adjoints can be checked in the induced
metric, and one and two step time forms test a positive transfer and its
compression defect. Optional regional forms test unital star algebras,
conjugation of their entire spans, pairwise commutation and joint cone
positivity. A supplied regional unitary is not yet identified with an
exponential of the supplied generators. That and a continuum theorem remain
explicit physical obligations.

The regional report also computes the rank of first-order commutators
modulo the regional spans and time, and distinguishes mappings between
distinct algebra spans from automorphisms of a single span. The rank tests
movement of the declared algebras and does not establish geometric locality.

The finite all-word equality check closes a reachable subspace on the joint
alphabet. It returns equality or a distinguishing word, without a heuristic
word depth cutoff. This is operational equivalence for the matched interface.
It does not identify unobserved sectors or establish a physical equivalence
of two continuum theories. On the cyclic spaces it gives a star equivalence
when adjoints are checked in both descriptions.

The queue also exhausts labelled finite simplicial complexes and rational
positive transfers at the reported bounds. A custom manifest can describe
simplices of any finite dimension or compare two exact coherent
representations. For an executable example use
`--cases research/originlab/reconstruction-cases.json`. The built-in
campaign reports which spaces were exhausted, refined, or remain open.
There is no finite list or decision procedure for every mathematical or
physical scenario.

For original R113 data use the optional `physical_reconstruction` field
inside the same endpoint bundle. It contains `identity_index`, `forms`
and optionally `step`, `twice`, `adjoints` and `regional`. All forms use
the original word-basis pair convention. The regional record uses
`translations` in order P0, Px, Py, Pz, `algebras` keyed by region, and
`transports` with from, to and operator forms. The reconstructed basis is
the supplied Gram space before any separate outer selection. R113 remains
blocked unless the original bundle and bound empty reference are supplied.

## Reconstruction and six origin routes

See [Six origin routes](SIX-ROUTES.md) for the six checkers, finite coverage, original input schema, and campaign command.

## Coherent regional sources

The default durable search includes R104, the R032 orientation witness,
relational transport and finite process tomography. Use
`--cases research/originlab/relational-cases.json` to add exact custom jobs.
`CASE_FIELDS` in `autonomy.py` lists the accepted job names and arguments.

- `r104_variant` varies the mode count, phase in units of pi, and number
  modulus. `det_word` takes explicit ordered CAR letters, nonnegative rational
  heat intervals, and an optional single outer number projector. Its result
  is an absolute word of the separated DET factor.
- `compare_memories` takes two labelled trace models with the same interface.
  Each has a matrix `cap` and a map `operations`. `cap_constraints` takes
  operations and a list of exact observed words with values. It returns
  the affine dimension and the next test when the responses are ambiguous.
- `transport_network` accepts directed unitary edges, optional pairs of paths
  and explicit face boundaries. `transported_caps` connects two cap phases
  to the same declared Abelian network, and `charge_visibility` determines
  exactly which declared integer charges see the holonomy. `flux_inference`
  uses Z3 to exhaust a declared grid of loop phases consistent with observed
  charged words and ranks probes that distinguish the surviving candidates.
- `process_tomography` takes `dimensions` ordered as
  `[O_k,I_k,...,O_1,I_1]`, normalized tester effects per slot and word
  observations. For nontrivial input legs each tester declares its positive
  `effect` and `input_state`; the effect is bounded by
  `I_output tensor input_state`. `extension_interventions` may provide
  a larger measurement alphabet. Reconstruction and nonuniqueness are
  exact on the supplied finite tensor size of at most 16.
  `observable_alphabets` computes ordered word ranks for two alphabets
  on one carrier. The Q4 comparison recovers ranks 26 and 36.
- `stationary_symmetry` checks finite channel covariance in a charge basis
  and the fixed operator space. `oriented_shifts` checks declared products
  of one-dimensional qudit shifts under a labelled reflection.
  `gaussian_source` evaluates exact connected responses from one Hessian.
  `joint_source` checks exact connected responses and supplied linear
  Ward identities from one table of ordered marks.
- `origin_dsl` accepts a `profile` containing atoms, each with a matrix,
  input and output boundary and source. A declared seam has two boundaries,
  a matrix and source. Optional `families`, `parameters` and `probes` control
  enumeration and test ranking. The durable default search accepts
  `--proposal-file` and `--finite-observations` as revision inputs.
  `--sol-proposals` explicitly enables transmission of the default finite
  grammar to OpenAI Sol. This mode does not send custom profiles or
  native TFPT source packets.

The R113 exact fibre sum also accepts `ordered_window_kernels` in place of
the four forms for an individual atom. It declares
`kind: whole_history_fixed_bosonic_scalar_U`, a bound scalar `U`,
`original_scalar_mark_binding` and the two Gaussian rational forms
`K0` and `KY`. The contraction assembles the mixed form before applying
the common selector and reference.

`relational_path_forms` accepts two distinct ordered regional words and
their original forms on the same basis. It needs the native source and time
IDs, endpoint fields, charge and spin dictionary, transport and admissibility
bindings, and one declared connection path comparison. A supplied
`reflected_gram` is checked for positivity. The report gives both absolute
values, the ratios where defined and the ordered holonomy. With a
`unit_charge_path_comparison_binding` it checks whether the two absolute
responses satisfy the expected Abelian character relation.
`coherent_memory` can predict both regional responses as additional held-out
words. A supplied `process_tensor` binds the same source, time and intervention
dictionary and uses independently supplied observations.
`joint_source_checks` binds additional forms and supplied Ward identities
to the same original contraction, including the empty response.

See the [plan audit](PLAN-AUDIT.md) and
[R104 and R113 findings](R113-R104-RELATIONAL.md) for the checked scope.

## Four-mark selection and spectral source

The default campaigns enumerate the split and full four-mark models,
check the native quark obstruction, test a fully specified finite spectral
alternative and recompute the frozen alpha and lepton comparisons.
See [the four-mark and source results](FOUR-MARKS-AND-SOURCE.md).

- `four_mark_refinement` takes binary invertible four-dimensional actions,
  optional `pin_values` as pairs of a four-bit vector and its value, and
  optional `allowed_vectors`. It returns all compatible nondegenerate
  quadratic refinements and ranks discriminating evaluations.
- `clifford_candidate` takes four exact `operators`, with an optional
  positive `metric`, a column `state`, and `square_phases` of ±1. It tests
  the relations on the cyclic span generated by the operations and their
  adjoints. Optional `charged_loops`, `cartan`, and `transport` also enter
  that span so weight, Clifford and transport tests use the same carrier.
  The latter takes a positive
  contractive `step` and four `transported_operations`.
- `primitive_charge_compatibility` tests a `homogeneous_charge` and its
  `central_power` in an integral or cyclic charge group.
- `spectral_source_candidate` takes declared rational half-fields,
  half-masses, momentum, half-time, cutoff, and spectral strength.
  It sums all pairs, tests the reflected kernel, and reports four marks
  and their connected response. The default variant uses the 96-state
  diagnostic fixture and one real Higgs axis.

A native packet can include `four_mark_selection` with its own source and
time IDs, orientation, adjoint and state bindings, empty-word
`identity_index`, and four original `operations` as forms. `square_phases`
chooses the tested convention. With `face_code_binding`, the measured
cyclic squares select the compatible full-face refinements. Optional
`primitive_loop_forms` and `cartan_forms` use the same word basis and
require a `charge_dictionary_binding`. An independently provided
`neutral_ordered_return_form` compares the product of the four loop maps
on the quotient and after the same outer selector, cap and reference. A bound `transport`
section can additionally provide a time step and its transported marks.

The watcher hashes the selected native sources and the frozen empirical
register. A changed file starts a new epoch and the review outbox binds
these source hashes with the results.

## Source rules as purity of two arcs

```bash
python3 -B research/originlab/main.py --campaign f1_two_arc_purity --finite-controls-only --out research/originlab/runs/purity
```

`f1_two_arc_purity.py` enumerates non-negative integer multiplicities of
the three pairings of the four seam marks and reads each rule as the
purity of two opposite arcs (hypothesis H-R2, round twenty-five of
`SESSION-FOUR-MARKS-20261118.md`). Exact guards: distant arcs decouple in both
channels, the seam state is pure, the number of exchanged fermion pairs
is an integer. The lane needs no model and no source packet. It returns
`original_r113_values: None`; that the P1 form is this purity is the
hypothesis, not a result.

## Synthesis from an open source question

The first search round now also runs the [graded sewing and return
campaign](PLAN-AUDIT.md#dynamic-sewing-and-a-directed-search). It
reconstructs a positive transfer from ordered return moments, solves the
inverse two-cell coupling problem, and tests each newly selected observation
on the actual witness operators in a second queued job. Later rounds retain
the research round ID and its follow-up results in the frontier.

The custom-case manifest accepts `complete_mark_graphs`, `clifford_sewing`,
`finite_gauge_link`, `link_sewing`, `return_generator_search` and
`return_moment_reconstruction`. The link job accepts `order` (2 or 4) and an
optional list of `graphs` with `name`, `vertices` and oriented `edges`.
The two-cell input accepts exact `strength`, `reference_weight`, `gap` and
`parity_signs`. The inverse case takes `projected_return_observations` in
the four labelled even sectors, the rational complement `gap`, and
optionally `require_positive_weights` and new independently supplied
`opposite_parity_observations` with `signs` and `value`. Use
`input_role: declared_finite_oracle` for both inverse and moment cases.

The R113 packet may supply `ordered_return` with the bound source and time,
metric-adjoint insertions `first` and `last`, and one positive `step`.
The checker splits the full return into the selected, complement and mixed
terms. It also checks optional independently supplied `ordered_return_form`,
`quarter_action` with `character_power`, and `temporal_modes` on the
same time step. `ordered_return_moments` accepts independently evaluated
matrix moments and held-out answers on the same source quotient.

An optional `affine_profile_returns` section reads neutral grade-zero words
and their charged grade-one images of one full root from the same reflected
Gram before selection. `neutral_indices` and each list in `charged_indices`
use zero-based word indices in corresponding order. The first charged list
uses the unit profile. In `products`, supply real coefficients over that
same profile list for all pairs `(a,b)` with `a <= b`, in lexicographic order.
Bind `source_id`, `time_id`, `grade_and_full_root_binding`,
`profile_product_binding`, `adjoint_binding` and `absolute_reference_binding`.
An optional `neutral_transfer_form` on the ordered neutral words, with a
`transfer_binding`, measures its commutator with the recovered profile action
on the supplied quotient. A rank defect gives the minimum number of additional
neutral return directions under the stated affine Ward hypotheses. A match
on the supplied forms establishes a consistency check on those words.

`affine_return_anchor` starts with the packet's original empty word at
`identity_index` and synthesizes the map of the vacuum and its first
charged current. It checks the metric, outer selection, adjoint return,
fourth-root action and supplied independent time calibration
`expected_charged_step`. Matching a declared anchor does not supply its
physical calibration or the full source.

The `original_affine_time` job composes the original positive E8 regional
translations and evaluates the marked three-current word at two independent
widths. The [derivation](AFFINE-TIME-PROOF.txt) also gives a lower bound on
the size of a single CAR fibre carrying the full quadratic E8 action. The
[level test](AFFINE-LEVEL-PROOF.txt) checks the central term of any local
free Majorana current extending that direct construction.

A native endpoint packet can supply `marked_affine_time_grid` with the
same `source_id` and `time_id`, and `native_field_and_charge_binding`,
`physical_to_radial_time_binding` and `absolute_reference_binding`.
The original `native_circle_marks` are `["1","I","-1"]` and the
`ordered_fields` are `["J^(-p1)","J^gamma","J^p2"]`. Their
`internal_glue_classes`, `color_trialities` and `hypercharges` are
`[1,2,1]`, `[0,1,2]` and `["0","-1/3","1/3"]` in that order.
Each of `left_axis` and `right_axis` has at least two entries with
increasing positive `physical_time` and decreasing `radial_factor`
strictly between zero and one. Both fields are exact rational numbers.
Supply the corresponding rectangular `ordered_word_forms` on the original
OS word basis. The checker uses the packet's outer selection and absolute
reference for every entry and reports mismatches and both matrix ranks.
A match tests the supplied time calibration on this grid.

The `original_folded_time` job uses the original fold matrix to calculate the
two marked quarter clocks on this same word. The grid reader reports the
mixed ratio, tests whether deviations from the radial word are separate end
factors, and, for a common exponential calibration, compares both original
quarter predictions. A stationary mixture of those quarters is determined
by two entries and tested on the other entries. The clock job also checks the original gauge average
and whether the marked opposite spinor module can keep the same background.
See [the marked folded time derivation](FOLDED-TIME-PROOF.txt).
