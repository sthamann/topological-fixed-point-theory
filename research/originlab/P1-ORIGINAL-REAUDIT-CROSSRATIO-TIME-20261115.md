# TFPT P1/F1: source re-audit of the four-mark geometry, time and full source

**Date:** 15 November 2026. **Scope:** exact finite source-matched controls and negative witnesses. **F1 is open.** None of the tests reconstructs the original weighted R113 packet or physical 3+1D parent. Preserve existing source/Atlas ownership and no ledger promotion.

## 1. Primary originals (not recent summaries)

* [Original primitive boundary collar](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex): the initial datum assumes a smooth oriented manifold with boundary and a first-order Dirac collar; the proof identifies \`P_C = (1+R_deck)/2\` without deriving the additional Calderon/deck compatibility.
* [Operational seed](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/002-section-intrinsic-closed-branch-and-canonical-rigid-obje.tex): its input already includes a local algebra, state, time evolution, reflection and an elliptic collar generator. It is not a spacetime-free generative primitive.
* [Shared c3/P1/P2 architecture](../../documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex) and [Origin Theory](../../documents/origin_theory/sections/001-section-the-whole-integer-skeleton-from-one-pair-gcar-nf.tex): the 4 punctures give \`rank H1=3\`, \`h0(O(D))=5\` by existing Gauss–Bonnet/Riemann–Roch results. Physical identification and the square modulus are not automatic.
* [Original v512 exact source](../../verification/v512_seam_tau_flag.py): a whole family \`D_delta={±1,±exp(i delta)}\`, \`0<delta<=pi/2\`, has the same free V4 action; its square is selected by an extra flag symmetry \`V4 -> D4\`. The free odd-doublet Dirichlet-to-Neumann covariance splitting is \`Delta_odd=(2/pi) log cot(delta/2)\`.
* [Original v239 thermal time](../../verification/v239_kms_thermal_time.py): Tomita–Takesaki fixes the KMS inverse temperature to one in modular units **for a supplied faithful state**. It does not independently choose the original interacting state or identify physical translations.
* [Metrology weighted source](../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex): the configuration-integrated projective consistency and OS positivity cannot be inferred from block compatibility and pointwise positive weights alone.
* [Primitive constructor audit](../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/triangulation_20261106O/caller/ACTUAL-PRIMITIVE-REGULATOR-CALLER-TRACE.txt): no executable full regulator with original bosonic/gauge/fermionic/BRST geometry, physical selector, eta, time and absolute reference was found in the traced owning paths. This is a scoped search result, not a global absence proof.

## 2. Correct projective symmetry: square D4, not tetrahedral A4

For ordered \`(1,-1,i,-i)\`, cross ratio \`(a-c)(b-d)/((a-d)(b-c))=-1\`. Of the 24 permutations, exactly eight stabilize it and none is a 3-cycle. This is a projective \`D4\` of order eight. A projective tetrahedral \`A4\` stabilizer requires \`lambda^2-lambda+1=0\`, incompatible with \`lambda=-1\`. The earlier A4-isotropic Hamiltonian was a valid abstract graph permutation control, **not** a geometrically derived square-seam Hamiltonian.

For the generic original v512 counterwitness \`D={±1, ±(3+4i)/5}\` only the four V4 permutations stabilize the cross ratio. Bare Clifford V4 or mark transitivity does **not** select the square. v512 already proved this and explicitly left flag side-indistinguishability as the remaining bit.

Furthermore the 3-dimensional graph cycle space \`H1(K4,R)\` is the \`S4\` standard representation twisted by sign; the 3-dimensional punctured-sphere homology is the untwisted standard representation. For the square rotation represented by the 4-cycle their traces are respectively \`+1\` and \`-1\`. Thus their naive identification is not even D4-equivariant. Any candidate needs an actual sign/charge/orientation twist.

## 3. Exact new relation between four-mode losslessness and the ORIGINAL odd splitting

Use the same v512 family \`D_delta={1,-1,u,-u}\` with \`u=exp(i delta)\`. Let \`V_delta[j,k]=z_j^k\`, \`k=0,1,2,3\`, be the four-mode *evaluation* matrix. Exactly,

\`det V_delta = 4u(u²-1)²\`; hence \`|det V_delta|=16 sin²(delta)\`.

Hadamard's inequality implies \`|det V|<=16\` for any four points of unit modulus. Equality requires the four points to be a rotated copy of \`mu4\`. So \`W=V/2\` is unitary if and only if \`delta=pi/2\` within the stated v512 family.

Define \`I(delta)=det(W W*)=sin⁴(delta)\`. With the original free odd-doublet splitting \`Delta_odd=(2/pi)log cot(delta/2)\`, elementary half-angle identities give the exact scalar identity

\`I(delta)=sech⁴(pi Delta_odd/2)\`.

At the old v512 counterwitness \`tan(delta/2)=1/2\`, \`Delta_odd=2log2/pi\`, \`I=256/625<1\`; at the square \`I=1\` and \`Delta_odd=0\`. This is a genuinely direct relationship between two **different diagnostics on the same geometric family**, not an equality of their physical source operators.

**Conditional new selection lemma:** If the original weighted P1 source itself demands a lossless isometric identification between exactly these first four Fourier modes and its four bound mark channels (and the physical domain/normalization agree), it selects \`mu4\` without a separately postulated flag bit. This required \`W*W=I\` **is not established by OS positivity, V4 symmetry or E8 Ward data alone**. It may be a new selection assumption; testing it is the actual origin task.

### Critical countermodel: existence of a weighted isometry does not select the square

The determinant maximum above implicitly fixes *both* the mark-channel and the four-mode Hilbert metrics to their flat Euclidean forms. That choice is not part of the original weighted R113 source. For **every** distinct four-mark configuration the normalized square evaluation matrix \`F=V/2\` is invertible. Define \`G=(F F*)^{-1}\`, which is strictly positive. Then the *weighted* isometry \`F* G F=I\` holds identically for **all** moduli; hence existence of some reflection-positive Gram metric with a lossless Fourier frame is **not** a square selector. On the actual original v512 nonsquare witness \`u=(3+4i)/5\`, the exact positive \`G\` has eigenvalues \`5/8\` (twice) and \`5/2\` (twice), and \`F*GF=I\`. The candidate would only constrain the modulus if the **source itself** supplies the channel metrics, normalization and their frame identification independently of the configuration. The existing original source constructor does not yet compute them. Thus one must not promote the finite Fekete selector into the F1 premise without the weighted physical kernel.

## 4. Original square arc Laplacian = two-cell Clifford stabilizer

v512 defines the arc-weighted cyclic four-mark Laplacian \`L_delta\` with adjacent weights \`1/delta\` and \`1/(pi-delta)\`. At the square \`delta=pi/2\`, its eigenvalues are \`0,4/pi,4/pi,8/pi\`; normalize its first gap by \`H_square=(pi/4)L_square\`.

Identify the four marks with bit labels of \`V4=Z2×Z2\`; let \`T_a=X⊗I, T_b=I⊗X\` be the two **primitive adjacent** transitions. In this order \`H_square=I-(T_a+T_b)/2\`; it has \`spec={0,1,1,2}\` and vacuum the equal four-mark superposition. Let \`S_R=X⊗X, S_eps=Z⊗Z\`; the two-cell Bell-seam Hamiltonian of the recent OriginLab is \`H_Bell=I-(S_R+S_eps)/2\`. An explicit unitary, the Bell joint-eigenbasis times the adjoint of the two-bit Fourier/Hadamard basis, sends \`T_a->S_R, T_b->S_eps\` and hence \`H_square->H_Bell\`. This tests the **intertwining actions**, not just the eigenvalues.

For \`t=ln2\`, \`exp(-t H_square)=1/16 [[9,3,3,1],[3,9,1,3],[3,1,9,3],[1,3,3,9]]\`.

**Important negative condition:** D4 invariance alone allows an *opposite* direct jump with extra rate \`b\`. The general positive D4 invariant square generator is

\`H(a,b)=a(2I-T_a-T_b)+b(I-T_a T_b)\`

with \`spec={0,4a,2a+2b,2a+2b}\`. Forbidding \`b\` means identifying only the two primitive adjacent arcs as elementary local moves. This is natural in the original circle grammar but must be shown for its physical transfer, not added as a substitute. Even then fixing the primitive time unit needs the original source.

## 5. Direct free-collar source negative control

The same v512 free massless circle with \`H=|D|\` has, after omitting the zero mode, time-separated two-point kernel

\`G_t(theta)=-(1/(2pi))log(1-2e^{-t}cos(theta)+e^{-2t})\`.

On the square at \`t=ln2\`, its finite four-mark matrix has eigenvalues

\`{log(16/15)/pi, log(3)/pi (twice), log(5/3)/pi}\`.

This is NOT the Markov-semigroup matrix \`exp(-t H_square)\`; the relative eigenvalues and the operator roles differ. It is incorrect to call the arc Laplacian the original physical time just because both share square symmetry. The full R113 source still requires all original fields, insertions, selector, unnormalized reference and limit order.

### Native affine parity obstruction forces additional grades or sectors

The parallel [native fibre/clock audit](P1-FIBRE-AND-CLOCK.md) now reports exact grade-one E8 glue-clock phase multiplicities \`(60,64,60,64)\` at eigenphases \`(1,i,-1,-i)\`. Therefore \`tr(U4²|H1)=60-64+60-64=-8\`; the affine vacuum grade zero has \`tr(U4²|H0)=1\`. If an invertible proposed seam polarization exchanger \`S\` anticommutes with \`Pi=U4²\` **and** preserves the affine energy grading \`L0\`, then on every finite grade \`tr Pi = tr(S Pi S^-1) = -tr Pi\` and hence must vanish, contradictory already at grades zero and one. In grade one any parity-flipping operator has nullity at least \`|120-128|=8\`. The existing report's minimal balanced 256-dimensional spectral *cover* at grade one remains only a clock-spectrum cover; it does not intertwine the original E8 currents or prove a physical seam.

A coherent parity pump can evade this finite-grade no-go by mixing grades: on a toy positive integer energy ladder, the self-adjoint swap \`S|2k>=|2k+1>\`, \`S|2k+1>=|2k>\` reverses parity but fails to commute with the grade operator. This makes a **dynamical sector-changing original P1 source** a concrete candidate, not an established fact. A genuine source-identified operator, full affine mode relations, physical reflected form and mixed R113 marks are necessary. Alternatively extra full graded sectors or non-invertibility can evade the premises; do not conclude that time emergence is mathematically forced by this conditional argument alone.

## 6. Modular time: a precise extension and a finite obstruction

The earlier v239 KMS witness takes the state and modular generator as inputs. A *proper half-sided modular inclusion* \`N ⊂ M\` with common cyclic-separating \`Omega\` can, under Wiesbrock/Araki–Zsidó conditions, reconstruct a positive-energy translation group from the two canonical modular groups. This is a real mathematical way to reduce primitive time assumptions, but a single modular automorphism group is insufficient, and the theorem gives a one-dimensional chiral translation, not 3+1D gravity by itself.

There is also a **finite-dimensional obstruction**: if \`N\` is finite-dimensional and \`sigma_t^M(N) ⊂ N\` for \`t>=0\`, then \`dim sigma_t^M(N)=dim N\` since \`sigma_t\` is an automorphism. Thus the inclusion is equality for each \`t\`; no proper half-sided modular compression is possible in finite dimensions. A maximally entangled two-cell Bell state additionally has tracial one-cell restriction \`I/2\` and trivial one-cell modular flow. Therefore a finite Clifford/OriginDSL model **cannot itself** be a nontrivial Wiesbrock origin of translations; only the infinite-dimensional affine/source net could, after a source-identified modular inclusion has been constructed. This is a scoped no-go for the proposed modular-inclusion mechanism, not for all finite dynamics.

Literature: Wiesbrock/Araki–Zsidó half-sided modular inclusions and positive translations; Takesaki's modular expectation theorem. These do not construct the physical original state.

## 7. One decisive original-source experiment

Do NOT create a replacement finite density or infer the original four R113 values. On the *original* source/absolute reference, calculate the unnormalized four-mark \`Q0[a,b;t]=I_orig((Theta F_a)theta_t F_b)\` on the common OS null quotient **before the outer physical selector**, together with the charge/adjoint and time actions. Three independent source equalities are needed:

1. Its four-channel Gram, after the *source-derived* Fourier/frame identification, satisfies or violates the isometry defect \`Q_frame^* Q_frame-I\`; its original odd two-point splitting must match the original v512 covariance on the same fields.
2. Its primitive interval transfer commutes with the actual square D4/deck and has no direct opposite **generator** matrix element if elementary locality is claimed. Compare the time series, not only the equal-time eigenvalues, with the exact \`9,3,3,1\` benchmark at the appropriately source-calibrated width.
3. Its actual charged E8 null-mode adjoints and mixed \`ZUY\` insertions preserve the same positive form and absolute reference, with nonzero charged returns and all intermediate channels retained.

If any fail, reject the single-rule isometric Fourier/C4/Bell candidate. If all pass, next derive the original interacting 3+1D parent, chirality, local observables and spin-2 coupling; these do not follow from dimension-three graph cycles or the chiral affine source automatically. **F1 remains unproved.**

Exact executable finite controls: [test_p1_geometry_time_bridge.py](test_p1_geometry_time_bridge.py). The seven analogous SymPy controls were independently executed before this note; these are regression scopes, not a full repository/Lean test.
