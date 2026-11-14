# OriginLab

Independent, finite exact research campaigns for TFPT. Previous search engines
are neither imported nor executed. Historical certificates supply benchmark
data. The original-source lane calls the native R113 consumer directly.
The native domain and canonical Atlas supply original provenance.

From the repository root with Python 3.11 or newer:

```sh
python3 -m pip install -r research/originlab/requirements.txt
python3 -B research/originlab/main.py --out research/originlab/runs/first
python3 -B research/originlab/main.py --out research/originlab/runs/first --resume --rounds 2
```

The default run executes R113, the calibration and model campaigns, P0
search and an F1 frontier report. For one lane, use `--campaign`. A run records exact
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

This mode runs the double-cell, moving Hardy carrier and affine polarization
checks without a local Atlas. Original-source and full F1 campaigns require the local Atlas
program and the matching, hash-checked research inputs. A public checkout
without those inputs cannot evaluate R113.

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

- [SymPy](https://github.com/sympy/sympy) executes exact matrix and polynomial arithmetic.
- [Z3](https://github.com/Z3Prover/z3) provides inverse integer synthesis, with exact witness replay.
- [Lean](https://lean-lang.org/doc/reference/latest/Axioms/) checks the separate rank lemma and reports its dependencies.
- [cvc5](https://cvc5.github.io/papers/2022/BarbosaBBKLMMMN-TACAS22.pdf) synthesizes finite Boolean local laws.
- [egglog](https://egraphs-good.github.io/egglog-tutorial/01-basics.html) replays declared compositional identities.
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
share their measured return words. Add `--watch-seconds 90` to keep the
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

cvc5 synthesizes distinct Boolean local laws and splits on both possible
answers to an automatically chosen experiment. egglog checks algebraic
equivalences of composable process expressions. Arb certifies the
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
