# F1: source-first five-point admissible-history compiler, phase selection and first positive time

**16 November 2026 · New exact origin-candidate mathematics and strict original-source gates, F1 OPEN.** This work does NOT begin with spacetime, the E8 lattice, a chiral VOA, 3+1D gravity, the physical time operator or a numeric fit. It begins with the ORIGINAL five-point binary odd quadratic-code and its permitted ordered quarter cycles. However the starting four source marks, the original spin refinement and the physical history amplitude are STILL taken from earlier finite TFPT mathematical premises, not derived from the physical Q0.

## Original starting point: the five q=0 points and their S4-fixed fifth point

[marked_search.py](marked_search.py) and [selection.py](selection.py) construct the original four-face S4-invariant odd quadratic refinement on F2^4. It has EXACTLY five nonzero q=0 vectors, pairwise symplectic pairing one and total vector sum zero. In the original BIT_VECTORS ordering, the fifth vector is (1,1,1,1). It is the UNIQUE fixed point of the physical four-face label-permutation S4; the remaining four vectors form a single S4 orbit. This is the original quadratic-code analogue of the user's "P1 + four marks = five points" compiler proposal. **It mathematically selects the fifth CODE POINT, not the physical identity of that point with P1 and not the four primitive marks themselves.**

An order-four permutation of five labels with exactly one fixed point is a 4-cycle on the other four. There are exactly 5*6=30 such cycles. If every five-point candidate is summed democratically before a physical seam is selected, the 5x5 classical transfer is J5/5, rank one, erasing every nontrivial direction. After the INTRINSIC S4-fixed fifth code point is held, exactly six possible ordered 4-cycles remain on the four original labels.

## Three DIFFERENT ways of combining the same six histories

### (A) Incoherent classical fixed-history averaging

The six anchored 4-cycles form one S4 conjugacy class. If the source explicitly preserves this FULL S4 relabelling symmetry and gives them equal nonnegative amplitude *probability* weights 1/6, the classical four-state transfer is

    T1 = (J4-I4)/3,
    spec(T1) = {1, -1/3, -1/3, -1/3}.

It has a negative subleading eigenvalue and is NOT itself the value at positive time of a symmetric positive generator H under a real exponential. After TWO primitive steps,

    T2=T1^2=P0+(1/9)(I4-P0), P0=J4/4,
    spec(T2)={1,1/9,1/9,1/9}.

T2 is positive definite, reversible, Markov and has the exact positive semigroup interpolation

    T(t)=P0+exp(-t log 9)(I4-P0),
    H=log(9)*(I4-P0).

There is one stationary scalar and a THREE-dimensional standard sector, but three such directions are NOT an index theorem for three physical SM chiral generations. The time unit is only "one paired source update" and requires the original physical absolute source/time dictionary.

### (B) Only the ACTUAL square geometry D4 is enforced

The six primitive anchored 4-cycles split under the square's geometric D4 into TWO inequivalent orbits of sizes 2 and 4. With relative total source weight x assigned to the two actual planar rotations r and r^-1 and weight 1-x to the other four cycles, the most general *within-orbit uniform* D4-positive transfer is

    Tx=x*(P_r+P_r^-1)/2 + (1-x)*sum_(other four) P_c/4

and its EXACT spectrum is

    spec(Tx)={1,-x,(x-1)/2,(x-1)/2}.

The square of this transfer is positive for 0<x<1; it has TWO generally different rates after paired time sewing:

    T_x^2=P0+x^2*P_alternating+((1-x)^2/4)*P_doublet,
    H_x=-2log(x)*P_alternating -2log((1-x)/2)*P_doublet.

The previous triple-degenerate rate log9 is achieved only if x=1/3. That extra equality is a demand for FULL S4 democracy of all six clocks, NOT guaranteed by the D4 square geometry and NOT supplied by original weighted Q0. The exact x=1/4 and x=1/2 witnesses both obey D4 symmetry, positivity after two steps, and full classical Markov sewing while having different time rates.

### (C) Coherent quantum history amplitudes with native-looking fourth roots

Weight all six anchored permutation operators by *coherent* fourth-root phases in {+1,-1,+i,-i}. Fix the first path phase to 1 (a common global phase has no physical action) and exhaust ALL 4^5=1024 relative assignments. Among them **exactly nine** weighted sums A=Σ_c w_c P_c obey A†A=4I, so the rescaled transfer U=A/2 is unitary. Three solutions have real phases and U²=I; six have true complex quarter phases and U⁴=I, U²≠I. All nine are exact finite-order reversible unitary clocks and none is itself a strictly positive dissipative transfer e^-H with H>=0.

Thus even when candidate spin/quarter phases are admitted, the code and unitary amplitude condition DO NOT choose a unique phase assignment, and periodic unitary motion DOES NOT automatically provide a positive OS-time heat semigroup.

Another distinct quantum operation is the honest incoherent random-unitary channel on the whole M4 algebra

    Phi(A)=(1/6)Σ_c P_c A P_c†.

Its 16-dimensional superoperator spectrum is

    {1 ×2, (-1/3) ×9, (+1/3) ×3, 0 ×2}.

After two steps its spectrum is {1 ×2, 1/9 ×12, 0 ×2}. Its two null directions persist, so Phi² cannot equal exp(-L) for ANY finite bounded, strongly continuous 16-dimensional quantum Markov generator L (every matrix exponential is invertible). This is a **no-go for directly equating the classical coherent four-state time with the full finite quantum mixed-history time**, not a prohibition on a larger physical quantum dilation or on real-time unitary quantum evolution.

### Exact nonembeddability of the full quantum two-step kernel, even away from S4 democracy

One might hope that adjusting the D4 orbit weight x removes the two zero modes of the fully democratic random-unitary channel and therefore restores a nontrivial positive quantum time. An exact test at x=1/2 shows this hope is INSUFFICIENT.

The full M4 quantum random-unitary superoperator at x=1/2 has exact spectrum

    spec(Phi) = {1 (x2), -1/2 (x3), -1/4 (x7), +1/2 (x1), +1/4 (x3)}.

Its positive symmetric paired transfer S=Phi² is INVERTIBLE:

    spec(S) = {1 (x2), 1/4 (x4), 1/16 (x10)}.

For a trace/KMS-symmetric strongly continuous quantum Markov semigroup e^{tL} with e^L=S, the selfadjoint generator L is UNIQUE by the real functional calculus, namely

    L=log(S)=log(2)*(-2 P_(1/4)-4 P_(1/16)).

This explicit exact 16x16 symmetric generator fails the necessary GKLS conditional complete positivity (CCP) tangent criterion. Its Choi matrix, restricted to the orthogonal complement of the maximally entangled identity vector Ω, has a negative exact eigenvalue -2log(2)/3. An explicit integer witness v in the 16-dimensional matrix-unit basis is

    v=(0,1,0,-1, 1,0,-1,0, 0,-1,0,1, -1,0,1,0).

It has norm squared 8, is orthogonal to Ω and obeys

    v* J(L) v / (v*v) = -2 log(2)/3 <0.

Thus no trace-symmetric quantum Markov generator can produce this full quantum transfer at finite positive time. The same original anchored *classical four-state* two-step transfer remains perfectly positive and admits a positive Markov semigroup. This demonstrates a necessary primitive SOURCE TYPING separation between classical probability sewing and quantum coherent CP sewing.

**Scope:** this excludes the unique selfadjoint log, as appropriate to KMS symmetry in the tracial finite comparison. It does NOT exclude non-selfadjoint GKLS logarithms with different branch phases, nontracial physical KMS states, non-Markov reduced dynamics, or a larger quantum source dilation. It is NOT an R113 physical original source evaluation. The declared x=1/2 is a rational held-out D4 source weight, not a physical measurement.

Exact reproducibility: [f1_compiler_kms_no_go.py](f1_compiler_kms_no_go.py) and [test_f1_compiler_kms_no_go.py](test_f1_compiler_kms_no_go.py). Both tests are registered in main OriginLab and Core Verification.

### Source-typing verdict

"Equal sum of permutations", "equal completely positive mixture of conjugations" and "sum of complex signed amplitudes" are THREE DIFFERENT OPERATIONS. The original physical q0/q1 quadratic refinement, external P_adm, original real source scalar w*s and actual original time must select one physically. Converting a classical reversible chain with spectrum {1,1/9^3} into a physical 3+1D QFT or E8 current algebra cannot be inferred from the matching number three.

## What this closes and what F1 must still derive

**Proved on the original quadratic-code geometry:** unique S4-fixed fifth code point; 30 order-four choices; exactly six clocks fixing that point; abstract rank-one erasure by five-point democracy; S4 equality gives T1, unique paired positive classical transfer T2 and H=log9; actual D4 symmetry leaves one continuous source weight x; 1024 coherent quarter-phase assignments give nine distinct finite-order unitaries, not one; the full CP random-unitary source is singular after two steps and cannot have finite bounded generator as its fixed positive-time transfer.

**Not proved physically:** why the four original marks and five-point code arise from pre-spatial microscopic source; which of the TWO S4 odd quadratic refinements actually defines the physical phase and orientation; whether the S4-fixed fifth CODE vector is physically P1; why all six histories carry the required amplitudes, phase signs and source measure; whether the full original P_adm chooses the coherent or incoherent operation; whether native Q0 has the same genuine OS time, native full96 spatial return and absolute empty reference. Native original R113 whole_source_bundle is still null, four absolute Z0/ZU/ZY/ZUY uncomputed, and 3+1D chirality/Spin2/E8 physical locality not derived. **F1 OPEN.**

## Decisive first original-source experiment

On the same actual native R113 history source, calculate the six accepted anchored history amplitudes before any averaging, using the *same* outer P_adm and absolute empty reference. Compute their D4 orbit ratio x, relative fourth-root phases and selected quadratic refinement. Then perform both the original coherent amplitude contraction and the original full CP/KMS/OS operation through two blocked levels: test whether a single positive self-adjoint H and its source Gram moments G,D,F reproduce BOTH the original measured Higgs U7 full kinetic Dirichlet form and the ordered U/Y/connected probes. If not, reject the candidate immediately. Do not substitute an invented positive Gram, an E8 partition, or full fermion trace for the original empty cap.

Reproduction: [f1_pregeometric_history_compiler.py](f1_pregeometric_history_compiler.py), [test_f1_pregeometric_history_compiler.py](test_f1_pregeometric_history_compiler.py), and the source-first companion [F1-DIRECT-SINGLE-ORIGINAL-SOURCE-20261116.md](F1-DIRECT-SINGLE-ORIGINAL-SOURCE-20261116.md).
