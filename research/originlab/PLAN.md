# OriginLab implementation plan

Baseline, 13 November 2026. Fundamental plan revision 354, live Atlas revision
`dae6e14e742903fb1209dd9e`, F1 in progress. The Atlas added
`problem:f1-session-synthesis-20261113`. Its regional and charged-word
results still require the original weighted form and bound empty reference.
The current [plan audit](PLAN-AUDIT.md) covers the original requirements and
the three later sets of scenario tests. The first priority is the
[same-source R113 test](R113-R104-RELATIONAL.md).

## 15 November 2026 · Original square-seam re-audit and independent source-time gate

The [source re-audit](P1-ORIGINAL-REAUDIT-CROSSRATIO-TIME-20261115.md) compares the actual primitive boundary/collar theorem, operational seed, original v512 seam modulus scan, v189/v216/v499 four-mark carrier construction, v239 modular KMS, metrology and the R113 constructor trace. **Correction:** the physical square configuration has projective D4 symmetry and cross ratio -1, NOT tetrahedral A4. The recent abstract K4 graph's S4 and its A4 permutation-invariant laplacian are *not* geometrically inherited from the original square; its three graph cycles differ by the orientation sign from the three punctured-sphere H1 cycles even under the square D4 rotation. Do not treat equal dimensions as a family representation map.

The previously known v512 one-parameter V4 family has an exact free odd-covariance splitting Δ_odd=(2/pi) log cot(delta/2), zero only at the square. A new four-mode Vandermonde evaluation V_delta on these **same** four geometric points satisfies det(V_delta)=4u(u²-1)² and the normalized channel-volume identity det[(V_delta/2)(V_delta/2)*] = sin⁴(delta) = sech⁴(pi Δ_odd/2). Thus a *source-derived* isometric four-channel Fourier interface would force the square and geometric quarter clock; **no original assumption yet supplies this isometry or identifies its Fourier labels with weighted physical P1 forms**. The original non-square counterwitness remains admissible absent that new source theorem.

The original v512 square arc-weighted Laplacian (gap normalized) is unitarily intertwined with the two-cell Bell/Clifford stabilizer Hamiltonian, including its two elementary operations and vacuum (not just its spectrum). Its candidate T(ln 2) is [[9,3,3,1],[3,9,1,3],[3,1,9,3],[1,3,3,9]]/16. But D4 symmetry alone permits an independent opposite jump: H(a,b)=a(2I-X⊗I-I⊗X)+b(I-X⊗X). The full original free circular DtN two-point kernel at the same time has a *different* positive Gram spectrum; this is direct negative evidence against calling an arc-graph diagnostic the original interacting temporal source.

Combine with the newly published [native parity and fibre audit](P1-FIBRE-AND-CLOCK.md): grade-one E8 phase multiplicities (60,64,60,64) give trace(U²)=-8, while an invariant two-cell parity-reversing clock requires zero trace. Thus even the finite unitary geometrical bridge is not an E8 physical current intertwiner; extra sectors or a genuinely nontrivial source interface are necessary. Keep the two-level physical fibre pushforward, all four R113 marks and absolute empty reference as the *next executable original test*, with full selector/OS and limits. No TOE or F1 claim, no ledger promotion. [Finite regression controls](test_p1_geometry_time_bridge.py).

## 15 November 2026 · Weighted fibres and the native quarter action

The native grade-one glue action has phase multiplicities 60, 64, 60, 64. The single-cell polarization half-step lacks adjoint phase −1. Two commuting cells predict zero parity trace on a common invariant grade-one space, whereas the native trace is −8. Any parity-reversing map on that space has nullity at least eight. An enlarged 256-dimensional spectral cover is possible, while the native current and orientation intertwiner remains open. The optional R113 two-level comparison checks all supplied marked fibres before final selection. The next original input is still the actual pair of weight-bearing levels, their configuration map, time and absolute references. [Exact checks](P1-FIBRE-AND-CLOCK.md).

## 15 November 2026 · Clifford quotient, true Z4 discriminant and root bridge

[Exact independent candidate and falsification gates](P1-CLIFFORD-QUOTIENT-ROOT-BRIDGE-20261115.md) now connect two anticommuting boundary involutions to an eight-element matrix lift, the four-element adjoint quotient V4 and, *conditional on all three elementary transition types*, the K4 mark graph with three independent cycles. Mark differences then yield the A3 root lattice and its cyclic **Z4**, not V4, discriminant group. The cyclic discriminant supports a genuine order-four operator action on the E8 glue currents; combined with the audited Z3 colour triality this gives the observed Z12 additive grading. The known D5+A3 index-four E8 gluing and 240-root count replay exactly. The K4 cycle space is the sign-twisted standard S4 module, so a physical orientation lift remains necessary.

A4 symmetry of the mark permutation carrier, if physically justified, restricts its positive vacuum-preserving transfer to the one-scale tetrahedral Laplacian \`H_mark=μ(I-J/4)\`. This is **only a finite mark-channel time candidate**, not the interacting source. Even imposing the already checked seam symmetry, unique entangled ground state, fixed gap and parity leaves distinct positive two-cell spectra, with different exact returns. Neither the four quotient labels nor a candidate spectral gap derives the real physical source. The tests and note are in [test_clifford_mark_glue.py](test_clifford_mark_glue.py) and [test_seam_index_and_time.py](test_seam_index_and_time.py). The original R113 weighted \`Q0\`, time map, \`Z0/ZU/ZY/ZUY\`, geometric locality and 3+1D parent remain open. Do not infer 3 physical families or P1 c3 from the abstract cycle/index matches.

## 15 November 2026 · Seam flow and source consistency audit

Seven finite negative controls now run in [test_seam_source_obstructions.py](test_seam_source_obstructions.py). They distinguish the Calderon spectral polarization from geometric reflection in an exact collar example: P Q P = P/2, and two copies give P2 Q2 P2 = P2/4. A product of anti-commuting involutions C = R sign(B) has C^2 = -I, so its adjoint action has order at most two; it cannot alone realize the nontrivial order-four current action required by the already audited Z12 gluing grading. A separate relative parity lift U^2 = Pi has the required i and -i operator phases in a two-state model, but is not an E8 Ward/source identification.

The reflection-covariant family B(s) = diag(s,2,-2), with R swapping the fixed +/-2 pair, satisfies R B(s) R^-1 = -B(-s), has exactly one positive spectral crossing and an eta jump of two, and forces an unpaired zero mode at the seam. This is a candidate mechanism for the spectral-flow-one P1 interface, not a derivation of the geometric four marks or a Lorentzian parent. The original Calderon/deck relation on D_coll is still unverified.

The source audit further confirms that a block restriction is not a weighted conditional integration (the effective Gaussian quadratic form uses a Schur complement); pointwise positive weights need not preserve reflection positivity; and a globally uniform lower bound K(s,t) >= m > 0 on t in [0,infinity) is incompatible with a finite output measure of the stated K(s,t)dt RG map. These are counterexamples to general inference steps, not counterexamples to every possible admissible TFPT source. In the original source one must explicitly check the finite weighted fiber identity, reflected positivity and normalization on its actual configuration domain.

The latest joint Z12 marked grades 1, 2 and 9 close to zero and permit a unit vacuum word, but they still do not fix positive time: an exact three-state graded return at equal widths gives 1/8 or 1/16 for two distinct positive transfers while leaving grades, adjoints and the untimed return unchanged. This is a finite grading-only nonuniqueness result, not a refutation of the full affine E8 Ward or weighted source conditions. It rules out promoting the new grading or relative Gram uniqueness alone to the original physical Hamiltonian.

F1 remains open. The first missing original object remains the weighted null-mark reflected kernel Q0 on the original full source, together with its absolute empty reference and time. The four R113 values Z0, ZU, ZY, ZUY are not produced by these finite tests. The next decisive test must bind the original weighted collar source to the existing marked affine grade-one adjoints and Z12 charge blocks, then evaluate its same-source physical transfer and marked return before any claim of completion.

## 15 November 2026 · Reflected form in the marked affine chart

The marked charged adjoint fixes the relative reflected form on the
identified affine grade one space. Its Ward constraints leave one
common scale, which the level-one vacuum return can bind to the empty
reference. The character refined modulo twelve distinguishes the two
parts of the coarse neutral class. The next original input is still
the weighted reflected pairing with its absolute reference and physical
time insertion. Compare that pairing and the marked adjoints before the
outer selection. [Derivation and source checks](P1-MARKED-AFFINE-METRIC.md).

## 15 November 2026 · Charged P1 field search and execution

The fixed folded chart provides three early marked lines in
$(\bar 3,1)_{1/3}$ and two early minus winding weak lines in
$(1,2)_{1/2}$. The six original affine products realize
$(3,2)_{1/6}$ on these lines. The coupled W benchmark with the
original neutral pair on both ends has exact two-time rank four
for nonzero coupling. A resolved rank-five time rectangle excludes
that fixed W word.

The new two-face collar search checks all six charged components in
both original quarter clocks. On those fixed spectra an intrinsic
variation and co-moving face frames can yield identical covariance
and heat derivatives. The jointly observed ambient hypercharge
derivative removes that ambiguity in the charged class. The solver
recovers the intrinsic field and replays covariance, charge and heat
together, including a negative heat control. The physical next input
is the original absolute two-face jet, with the independently anchored
hypercharge matrix and its cross-face blocks, on the same P1 time and
weighted reference.

The Boolean lane exhausts all 64 observed-compatible rules in the
declared binary radius-one grammar. Both answers of its best probe
leave 32 rules. The autonomous runner extends its saved resolution
bound and records the original source block separately from failures.
The ordinary runner continues when the local Atlas is absent. The
public R006 replay includes compressed copies of the exact original
benchmark bytes. F1 and the weighted R113 values remain open.

## 14 November 2026: energy moments and regional motion

Fundamental plan 361, Atlas `32bb9da024c4595c2a8380f7`. The original Higgs
Fourier replay gives local kinetic coefficients `169/40000` for one moving
mode and `5117/5000` for the complete local gradient, in the same units.
The original side switch preserves the current product and gives no new
regional translation by itself. The consumer now tests regional motion
modulo algebra preserving actions and checks the compressed energy square.

For a positive Gram G on the null quotient, write D for the first energy
form and F for the second. Put E = F − D G⁻¹ D. Its rank counts the minimum
missing directions of K on the observed span. If E is nonzero, choose a
positive principal Gram E_J on independent columns and let C = E[J,:].
Then C* E_J⁻¹ C = E. In metric diag(G,E_J), the two energy forms with
blocks (D,C*,C,C D⁺ C*) and (D,C*,C,C D⁺ C* + E_J) are positive by the
Schur complement, using ker(D) ⊆ ker(C). Their visible G, D and F match
and their visible third moments differ by E. This yields a concrete
distinguishing word coordinate in the finite extension class.

Two bounded OriginDSL rounds evaluated 1,872 terms each and left 19 and
24 finite survivors. The native R113 reader reported absent original
endpoint operation and absent fibre packet both times, leaving the four
marks unset. F1 remains open. The next original check needs the reflected
weighted word form, two energy moments on the same charged and regional
operators, and the absolute empty reference. Derive an actual additional
regional action before reading the joint continuum limit or tensor sector.

## Implementation and acceptance

The double-cell campaign supplies an exact occupation and symmetry control
for the P1 to P2 transition. Bind its weak slot exchange and compensating
parity change to one weighted marked source before outer selection. Check the
adjoints, Cartan and weak actions, common time and absolute reference on that
same source. Reuse the conditional circle generator and the P1 transfer
reconstruction as time controls. A relative-phase model and a choice of
symmetric energy levels alone leave these original-source obligations open.

See [the 14 November implementation audit](PLAN-ABNAHME-2026-11-14.md)
for the full OriginLab 2.0 plan comparison, the executed benchmarks and
the distinct open software and physical proof obligations.

The 14 November finite marked campaign replays the boundary-preserving
binary transvection, the five-point A4 reduction, both integral E8 root
count partitions, the meromorphic quarter clock, and four symmetric
orthogonal cell sewings. The minimax planner needs at most three
point-image queries to distinguish the 30 clocks and one mixed query to
distinguish the four sewings with the same marginals and spectrum. These
are declared finite observation alphabets. Original physical admissibility
of the corresponding words is an independent condition. The gap scale g
is free. Both E8 constructions assume their stated integral lifts.

The exact complex R113 control has the same full equal-time matrix and
the same diagonal responses for all steps in both witnesses. Their real
mixed parts coincide too. The search chooses the imaginary part at step
one and finds opposite values of magnitude 1/4. Both positive block
spectral sequences replay through the existing Hankel verifier and two
withheld steps. Boundary feasibility remains undecided if no certified
direction or feasible seed is found. This fallback is distinct from an
exact exclusion. Source evidence changes now invalidate watched R113
inputs even when the endpoint packet itself is byte-for-byte unchanged.

The R113 synthesis campaign now exhausts real relative signs on a fixed
simple spectrum and solves the weighted response range on a bounded
rational probability polytope. The charged control has two positive
models. Both match the full zeroth moment and every individual mark
moment. The first mixed step separates them at `7/8` and `9/8`. Supplying
`9/8` leaves one of the two mixed spectra; supplying `1` excludes the
declared real sign class. The positive block Hankel reconstruction
replays both witnesses including two withheld matrix moments each.

On the unchanged native Q domain the live original checker finds zero
off the centre blocks for all ten reviewed generators and the 384
dimensional Fourier operator. In contrast, the actual neutral E8 word
has the right-to-left sector path `0 → 2 → 0 → 0`. In particular its
middle current requires the `2 → 0` block. No addition, product, adjoint
or commutator of the reviewed block-preserving generators supplies it.
The native endpoint still returns four missing values. The next native
test is to bind that middle block and its ends to the original weighted
form before outer selection and evaluate the word at two positive time
widths on the same absolute reference. The real sign control establishes
neither the complex phase nor the physical weighted form.

The 14 November Fable 5.2 review suggested cycling the three affine roots
to a horizontal matrix element between two states in colour 2. The exact
root check confirms that the states have weights p2 and p2−p1 and that
the unitary elementary string has normalized matrix element of modulus
one. The neutral root −p1 has nontrivial full Cartan weight, even though
its colour and hypercharge are zero. Its availability in the Q algebra
is a separate premise. The original charged middle is also a unit doublet
arrow on grade one, with both ends at radial grade one. The old fixed
quarter clocks split these energies by −11/2 and −9/2 in H8 units.

The affine two-time word at x,y in {1/2,1/4} has rank two and cross ratio
81/85. A closed doublet with one line for each intermediate weight gives
a product of one response per width and rank at most one. This rules out
that compression as the whole affine word. The new endpoint checker
tests the supplied horizontal root data on the reflected quotient,
including relative norms, true adjoints, colour, and time leakage. The
physical Q Cartan and mode maps, and the same weighted absolute reference,
are still open. The next native test is the colour-2 horizontal matrix
element on a source-identified Cartan pair, followed by the marked
two-time rectangle before any claim of physical closure.

The operator campaign of 14 November independently checks the closed affine
kernels with level-one mode operators and adds a class-constrained completion
search for the mark alphabet. The marked letters together with individually
hypercharge neutral letters generate A1 x E7 and preserve the weak slot
parity. All 112 odd roots are weak doublets. For the original test replace
the middle leg by `e1+e4` and `e6-e4`. Their charges add to the original
middle root, so the complete word stays neutral. The two orders have
distinct mode coefficients. Measure both on the same weighted form, time
map and absolute reference. The double-cell slot exchange inside the even
weak SU(2) supplies a separate parity test on its own carrier. The original
weighted endpoint remains open. F1 remains in progress.

For the affine weak split the mode coefficients extend to every $m,n\geq 1$.
In the forward order the lattice contraction of $(-p_1,a,b,p_2)$ has
four pairings equal to $-1$, on $(0,1),(0,3),(1,2),(2,3)$, and two
pairings equal to zero. In $|z_0|>|z_1|>|z_2|>|z_3|$, expand each
inverse difference as $z_i^{-1}\sum_{k\geq0}(z_j/z_i)^k$.
Let $u,v,w,x$ be the four nonnegative expansion indices in edge order.
The modes $(m,1-m,n-1,-n)$ require $u+v=m-1$, $u-w=m-1$ and
$x-w=n-1$. Thus $u=m-1$, $x=n-1$, $v=w=0$ is the unique term.
The forward coefficient is the cocycle times $1$. Interchanging the
middle modes subtracts the cocycle times the original middle word,
whose coefficient is $\min(m,n)$. Hence the normalized pair is
$1/\min(m,n)$ and $1/\min(m,n)-1$ for all positive mode indices in this
affine vacuum convention. The 36 independently evaluated mode pairs
check the fixed $D_8$ dictionary and phase. The weighted R113 forms
still need their own source binding.

### Weak fields, signed winding and the native source

In the fixed $D_8$ chart hold the middle root γ = e₁ + e₆ fixed.
Seek a right root r with Y = 1/2, weak weight of modulus one, zero
colour Cartan weight and vanishing colour step action. Require γ + r to
be a root. Exactly six roots qualify. They form one $D_8$ doublet and
two spinor doublets, exchanged by the reflection at e₇ − e₈. Completing
the colour orbit gives six components per alternative and three for the
marked word. All have ordered A2 pairings −1, with the cocycle-corrected
coefficient μ(m,n) = min(m,n). Six gauge step generators annihilate each
completed tensor. The Cartan actions vanish because the three charges
sum to zero.

Chronologically the marked word carries (bar3, 1) at Y = 1/3 and then
(1, 1) at Y = 0. All three alternatives carry (1, 2) at Y = 1/2 and then
(3, 2) at Y = 1/6. A fixed external weak rotation in the first interval
therefore reads a normalized doublet trace on all alternatives. Two
rotations compose as the ordered product of the interval matrices. A
fixed colour rotation reads the conjugate defining character in the first
marked interval and the defining character in the second weak interval.
For correlated dynamical links the weighted fibre integral needs its
actual gauge covariant dressing.

The original coweight −2e₈ gives signed interval windings (−1,−1),
(0,0), (−1,−1), (+1,+1) in marked, D8, minus spinor and plus spinor order.
Their glue pairs are (1,3), (2,0), (3,1), (3,1). On the early interval
compose the quarter insertion of that coweight with the weak diagonal
element with first entry z = (3 + 4i)/5. Its k-th normalized response is

$$\big((-i)^k,\;c_k,\;(-i)^k c_k,\;i^k c_k\big),\qquad
c_k=\tfrac12(z^k+\bar z^k).$$

Powers zero through three give determinant
−5376/3125 + 44928i/15625. They recover four absolute coefficients in
the specified linear response class. Power four adds an independent
check, including at a radial point with both widths unseen. The quarter
insertion is an external diagnostic of the same coweight, distinct from
the closed cocharacter loop and the derived old $U_G^4$ half phase.
The uninserted radial coefficient is the sum of min(m,n) x^m y^n over
positive m,n. Expanding min as a sum over common lower cutoffs gives
xy/[(1−x)(1−y)(1−xy)]. The campaign emits these affine targets while the
absolute whole weighted form remains an original input.

On the original 96 dimensional fibre write the two actual Higgs blocks
as A from the conjugate sterile neutrino to the lepton doublet and D
from the right down quark to the left quark doublet. Their family ranks
are three. Suppose the true affine lowering field is this Higgs
adjoint F and the new charged field is X for δ = γ + α. The affine
identity [F,[F,X]]=0 follows because δ−2α is no root. On the audited
early weight, X cannot act because the shifted gauge sector is missing,
and F squared vanishes. The projected double commutator is
−2 D* X A*. It forces the corresponding block of X and both middle
blocks of [F,X] to vanish. The same proof applies after enlargement
when the old fibre reduces the alpha actions and gauge Cartans and the
supplied charged field still has its required weight. An altered or
nonreducing field map changes this hypothesis.

For the polar normalization of that literal Higgs field, 16 disjoint
full rank families give 48 doublets. The polar Cartan has only the
weights +1 and −1. On its exterior space exp(iπ dΓ(Hα)) equals the physical
fermion parity. Every even CAR operator has even adjoint weight, whereas
the required δ field has weight +1. In a larger representation it must
exchange the sectors of K = Π exp(iπ Hα). The horizontal A2 control and
the original D8 pair action on a CAR vacuum realize that exchange when
an alpha singlet mode is present. The horizontal model has zero return
for alpha first and unit return for gamma first. These controls concern
horizontal alpha and a first D8 pair, with no extrapolation to a
physical negative mode on the weighted source.

The original separated R104 DET16 block has a unique even ground and
all odd states at energy two. Every odd insertion B from the ground
returns with exp(−2t) times its initial squared norm. The excited even
and odd dimensions are 32767 and 32768, so a time commuting odd mirror
operator has rank at most 65534 and annihilates the ground. The code
checks all sixteen Majorana ground returns and a thermal word with
unequal intervals. A time neutral unitary odd spectator on this whole
block is impossible. A coupled insertion or an energy compensation from
another part of the source requires a new original field and time map.

1. Build an independent endpoint evaluator, typed process language, exact
   search and certificate engine. Read the live Atlas for provenance. Preserve source IDs, input hashes,
   proof scope and the distinction between unavailable and false.
2. Execute the native R113 domain first. Emit the four marks and connected
   response only when supplied by the original evaluator. When blocked,
   identify the original producer, the missing operations and the proof
   obligations for the original configuration map, measure and reference.
3. Recover the R102 selection from the variable discriminant family. Replay
   the two R006 processes and generate a distinguishing word independently
   of their stored query. Check nonidentifiability of a finite positive
   transfer eigenvalue and the actual Cartan versus spatial quarter action.
4. Run typed finite process enumeration, mutation, inverse constraints and
   optional model proposals through exact checks and active word discovery.
   Keep typed counterexamples with the tested class and assumptions.
5. Make campaigns reproducible and resumable. Extend one finite class at a
   time when a surviving candidate gives a discriminating original-source
   question. Pass reviewed scientific results through the existing Atlas and
   Fundamental Plan routes.

## Scientific boundary

The acceptance registry tracks G0–G7 and requires the same source, matching
maps and original evidence across gates. The finite backend can generate and
refute hypotheses. A complete solution requires proofs of the physical
identifications and limits, followed by independent empirical validation.
Compare the original conditions in
[Osterwalder and Schrader II](https://doi.org/10.1007/BF01608978).

Use a standalone Python engine with exact arithmetic. The search and verifier
do not import or execute the previous TFPT engines.
The native replay has no dependency on a new Rust runtime or a second Atlas.
The binary synthesis lane enumerates exact truth tables. The inverse lane
uses Z3 when available. The finite campaign has an explicit
carrier and finite search bound. Any later backend must document which
additional assumptions it introduces and which question it can decide.

## Remaining implementation

The finite search runs now. Next expand the independently typed generator
families beyond common matrix primitives, implement general rewrite and
synthesis backends, and add a regulator interface with proven remainder
bounds. A production scheduler needs independent jobs, targeted retries
and reviewed publication through the canonical Atlas. The current local
runner records results and can resume or wait for changed sources. It does
not publish them automatically. The seven physical gates after G0 require
original proofs and evidence. Finite controls leave them open.

## Autonomous scenario extension, 13 November 2026

Implemented in autonomy.py: durable subprocess queue, automatic retry,
resumption, revision watching, experiment ranking and refinement of
indistinguishable scenarios. The independent scenario backends compare
ten topologies, two conductance assignments and five locations or
mechanisms of dynamics. cvc5 supplies a separate symbolic automaton
generator, egglog the declared compositional identities, and Arb a
certified continuum comparison for the free periodic regulator. The
source review outbox can be submitted to the Atlas with an actual review
and matching revision.

First scientific priority remains the missing original R113 kernel and
bound empty reference. The supplied A3 and charged P1 text adds two
decision branches for the same original source. One asks for the
physically admissible primitive charged operator and its marked return.
The other asks for compatible regional modular data and a genuinely
commuting translation representation with the joint positive spectrum.
The new finite controls exclude shortcuts in both branches. The actual
source, physical cutoff bounds, chirality, spin two and G0–G7 remain open.

## Reconstruction and operational uniqueness

First implement and test exact multipoint discrimination, state dependent
Gram reconstruction, and equality for all finite words on the same marked
interface. Validate local algebra spans, joint generators and compression
on supplied coherent data. Expand beyond the named topology list with
bounded exhaustive complex generation and arbitrary finite complex inputs.
The task can distinguish completeness of the supplied algebra, cyclic
reconstruction in a given state, and equivalence within a declared
operational interface. Continuum physical equivalence and source selection
require separate proofs.

The six-route campaign implements the tests in [Six origin routes](SIX-ROUTES.md). The first physical priority is a shared source binding for exact Petz recovery and original moment closure. Regional reconstruction, inflow, and the full scattering problem remain open beyond the finite tests.
