# TFPT F1: Charged E8 root word and a paired determinant origin test

**15 November 2026. Research result and source candidate, not a solved F1.** This document distinguishes exact source-adjacent free E8 current identities from a new *conditional* original-source hypothesis. No full weighted R113 values, Pfaffian/eta exponent, physical Hamiltonian or 3+1D parent is claimed.

## Actual source owners

- [Native marked E8 roots](affine_operators.py), [affine marked degrees](marked_source.py), [relative metric](P1-MARKED-AFFINE-METRIC.md) and [original return](R113-R104-RELATIONAL.md). The native roots are p1=p2+gamma with p1, p2, gamma all norm-square two; their marked steps carry the Z12 grades 1, 2 and 9.
- [Original boundary/collar theorem](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex), [full primitive source trace](../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/triangulation_20261106O/caller/ACTUAL-PRIMITIVE-REGULATOR-CALLER-TRACE.txt), and [weighted R113 source reader](source.py) still require the original marked field/history integration, physical selector and absolute empty reference.
- [Previous Calderon and full-Gaussian source falsifiers](P1-PRIMITIVE-GLUING-SOURCE-GATE-20261115.md) show that the square does not follow from generic V4 geometry, 3+2 sections or free endpoint Gaussian weighting.
- Standard free lattice vertex formula: a neutral chiral exponential word has a fixed cocycle phase times the pairwise product of complex coordinate differences raised to the charge inner products. This is a *chiral affine E8 source component*, not the whole TFPT interacting R113/OS source. A reference is Ginsparg, Applied Conformal Field Theory (1990).

## Exact E8 root-triangle no-go

Use the actual native E8 root charges q=(p2,gamma,-p1,0) on the original antipodal geometric marks (1, exp(i delta), -1, -exp(i delta)), with 0<delta<pi. The charges sum to zero and their three nonzero pairings are all -1. The free level-one E8 root vertex amplitude has **absolute modulus**

    |A_root(delta)| = 1/(4 sin delta).

It has a strict **minimum** at delta=pi/2 and diverges when marks collide. The fixed lattice cocycle and Z12 charge labels affect phase, not this amplitude shape.

The complete E8 census has 240 roots, each with exactly 56 other roots of pairing -1, and therefore 13,440 **ordered** root pairs leading to a root sum. All produce exactly this three-point source shape. This excludes a bare individual free E8 root triple as a mechanism selecting the square *by maximizing its own modulus*. It does NOT exclude the full interacting/charged source.

### General neutral free-word no-interior-maximum theorem

For any real Euclidean charge vectors q0,...,q3 with sum zero, the absolute free vertex correlator on the antipodal four-mark family has the shape

    W(delta)=constant*(2sin(delta/2))^A*(2cos(delta/2))^B

where A=q0·q1+q2·q3 and B=q1·q2+q3·q0. Neutrality gives

    A+B=-||q0+q2||² <=0.

Put x=sin²(delta/2)∈(0,1), so

    log W = constant + A/2 log(x)+B/2 log(1-x).

An interior stationary point requires A/B=x/(1-x)>0, hence A and B have the same sign (except degenerate zero cases). Both positive contradict A+B<=0, and both negative give a strict **minimum** by the positive second derivative of log W. Opposite-sign exponents give monotone variation. If both vanish the amplitude is flat. Thus no nonconstant single neutral Euclidean free-boson charge word has a strict interior maximum. Interference among multiple channels, descendant currents, interactions and non-positive metrics are **outside** this theorem.

## New minimal conditional candidate: paired boundary determinant and root triple

The proposed relative polarization U_t from [the primitive source audit](p1_source_falsifiers.py) is a two-dimensional unitary with t=tan(delta/2), eigenvalues u and u^(-1), where u=(1+it)/(1-it). Its two scalar secular determinants are

    Dminus=det(I-U_t)=4 t²/(1+t²),
    Dplus=det(I+U_t)=4/(1+t²),

and their product is

    Dminus*Dplus=16 t²/(1+t²)²=4sin²(delta).

If the **actual original weighted P1 source**, on precisely the same charged root-word and reference, contains one net positive power of this entire pair as a remaining **physical** defect/Pfaffian contribution, its candidate modulus is

    W_candidate = Dminus Dplus |A_root|
                = sin(delta) = 2t/(1+t²).

This has a unique positive maximum at the square. It is invariant under delta -> pi-delta, and vanishes at mark collisions. These are **synthetic model predictions**, not computed R113 source values.

Held-out exact relative weights normalized to the square, useful for the existing OriginLab marked-angle study:

| t | 1/3 | 1/2 | 2/3 | 1 |
|:--|--:|--:|--:|--:|
| W_candidate(t)/W_candidate(1) | 3/5 | 4/5 | 12/13 | 1 |

### Sharp net determinant-exponent threshold

If the **effective net power** of the paired determinant is m after all original fermion Pfaffian, boson inverse determinants, ghost, eta, regulator and admissibility factors, the declared pure shape is

    W_m(delta) proportional sin(delta)^(2m-1).

The square is a strict maximum iff m>1/2, is FLAT at m=1/2, and is a strict minimum for m<1/2. A bare half determinant cancels the native E8 pole but selects no modulus. The statement that TFPT has 16 Majoranas does not fix this effective power: their mode maps, multiplicities, boundary determinants and ghost cancellations must be derived in the **original source**. The paired determinant is nonholomorphic as an absolute amplitude; it cannot be inserted into the holomorphic affine E8 OPE without an independently constructed nonchiral/full source sector.

Importantly, an ordinary geometric reflection sending U to its inverse preserves Dplus and Dminus separately. It does **not** exchange them and cannot imply their equality. A genuine U -> -U swap would impose an additional charge/parity operation whose action must be derived on the original weighted reflected Hilbert space. The proposed **product** of the two determinant amplitudes is thus a testable source hypothesis, not an established consequence of reflection positivity.

## Mandatory original-source F1 test

On the actual finite blocked original source with its original Higgs, gauge, geometry, fermion, ghost, eta, selector and common reference, compute the charged neutral-root return and all four R113 values at the square and at t=1/3, 1/2 and 2/3. Factor out only the provably source-identical free E8 root component. Determine whether its remaining physical scalar, including its BRST/Pfaffian/eta weight, has the above pair-determinant exponent and angle dependence. Verify at TWO consecutive blocked levels, including the actual weighted fiber pushforward and the same positive-time Gram.

Failure rejects this candidate but does not refute the original E8 current algebra or P1 and P2 as conditional source axioms. Success would establish a potential square modulus selector **within a stated source class**, not the emergence of interacting 3+1D locality, the complete chiral Standard Model or universal spin-two dynamics.

Exact control and native 13,440-root census: [p1_charged_modulus.py](p1_charged_modulus.py), [test_p1_charged_modulus.py](test_p1_charged_modulus.py). The module has explicit F1=false and full-source=false gates.
