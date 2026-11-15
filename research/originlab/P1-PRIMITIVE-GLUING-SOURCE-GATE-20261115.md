# TFPT F1: return to the primitive source and full gluing (15 November 2026)

**Verdict:** New exact constraints and one source-compatible hypothesis; no solved F1, no original R113 four-value evaluation. This document audits the original owner and the most recent OriginLab results, not only their summaries.

## Source provenance

* [One-sided primitive boundary theorem](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex): input already assumes a smooth compact oriented manifold with boundary, spectral triple, Dirac product collar, real structure and grading. Its displayed identification `P_C=(I+R_deck)/2` requires additional proof as a relation on the same boundary Hilbert space; the stated proof's splitting premise is not a generic Calderon theorem.
* [Architecture shared P1/P2 origin](../../documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex), [original flag counterwitness](../../verification/v512_seam_tau_flag.py), [Riemann–Roch](../../verification/v189_riemann_roch_carrier.py) and [even Clifford](../../verification/v197_rr_carrier_clifford_d5.py).
* [Original metrology](../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex): the proposed projective geometric measure and Gaussian regulator require weighted conditional pushforwards and actual reflection positivity, not merely compatibility of blocked Dirac symbols.
* [Actual primitive source trace](../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/triangulation_20261106O/caller/ACTUAL-PRIMITIVE-REGULATOR-CALLER-TRACE.txt) and [R113 original consumer](source.py): the first missing operation is the **actual** joint weighted source `Q0[a,b;t]`, `Z0,ZU,ZY,ZUY`, with physical selector, time, all field components and absolute reference.
* External mathematical control: B. Booß-Bavnbek and K. Wojciechowski, *Elliptic Boundary Problems for Dirac Operators* (1993), Chapter 12 (Calderon is a pseudo-differential Cauchy-data projection, with the positive tangential spectral projection as principal symbol), https://thiele.ruc.dk/~booss/A-S-Index-Book/Draft/1993_Booss-Bavnbek%26Wojciechowski_.pdf ; Kirchhoff's weighted matrix-tree theorem and the continuum determinant-gluing identity.

## I. The recent 3+2 source sequence is universal on the antipodal family

Take `D_u={1,-1,u,-u}` for any `u!=0`, `u^2!=1`. Put

  `P_u(z)=(z²-1)(z²-u²)`, `f_k=z^k/P_u(z)` for k=0,...,4.

The involutions `r2:z->-z` and `s_u:z->u/z` preserve the divisor and commute. They act on `f_k` by `(-1)^k` and `u^(k-2)f_(4-k)`, respectively.

Their fixed points are 0 and infinity for `r2`. Evaluation at those points is a surjection of rank two:

  `ev_(0,infty) = [[1/u²,0,0,0,0],[0,0,0,0,1]]`.

Its kernel is span(f1,f2,f3). Multiplication by `dlog z=dz/z` identifies this kernel with all three logarithmic forms
`omega_j=z^(j-1)dz/P_u(z)` (j=1,2,3). Under `s_u`, `dlog z` changes sign, and its winding at 0 is +1. The exact 3+2 sequence, orientation and unit winding therefore hold for **every** noncolliding antipodal modulus, including the original v512 nonsquare witness `u=(3+4i)/5`.

Under the extra proposed geometric quarter `z->i z`, the divisor polynomial changes by precisely `2 z²(1+u²)`. Thus this operation is a genuine geometric automorphism **only** if `u²=-1`. The recent exact 3+2 cohomology result does not by itself establish the extra quarter. This is a necessary correction to any F1 claim built just from `3,5,2`, orientation or `dz/z`.

## II. A seemingly unique free-Gaussian square saddle cancels under full source gluing

The v512 arc-weighted four-cycle has edge lengths `ell=(delta,pi-delta,delta,pi-delta)`, conductances `1/ell`, and finite endpoint Laplacian `L=d1 diag(1/ell) d1*`.

The weighted matrix-tree theorem gives for **any four positive arc lengths**:

  `det'(L) = 4 (sum_i ell_i)/(prod_i ell_i)`.

For a circle of total circumference `2pi`, the alternating family gives

  `det'(L_delta)=8pi/[delta² (pi-delta)²]`.

If one integrates **only the four endpoint values modulo the constant mode**, their Gaussian factor is `det'(L_delta)^(-1/2) proportional to delta(pi-delta)`, apparently maximized uniquely at the square.

But the Dirichlet Laplacian on each arc of length `ell_i` has zeta-regularized determinant `2 ell_i`. The four independent interior fluctuation factors multiply to `[prod_i (2 ell_i)]^(-1/2)`. The full partition of the *unchanged free circle* therefore obeys

  `det'(L) prod_i (2 ell_i) = 64 (sum_i ell_i) = 128pi`.

The spurious square preference cancels **exactly**. It would be inconsistent for a fully glued free field to depend on arbitrary locations where we subdivided its circle. A square-selecting origin must come from a **physical puncture/defect, charged word, nontrivial holonomy, original eta/Pfaffian contribution or another original interaction**, after this base free-gluing cancellation.

This is a sharper restriction on the physical P1 origin than a free determinant saddle. The sign, magnitude and even existence of the residual interaction cannot be inferred without evaluating the original all-field source. Faddeev–Popov/ghost factors can cancel or reverse a partial determinant; never choose a sign by hand.

## III. Calderon and deck are distinct: explicit original-normal-form countermodel

In the exact product Dirac cylinder with `D=gamma_n(partial_n+B)`, `B=(-i partial_theta)*sigma3`, `gamma_n^2=-I`, and `{gamma_n,sigma3}=0`, the positive tangential spectral subspace of Fourier mode k is

  `P_C(k=+1)=diag(1,0)`,
  `P_C(k=-1)=diag(0,1)`.

It is impossible for **one frequency-independent fiberwise geometric deck involution R** to satisfy `P_C(k)=(I+R)/2` at both k=+1 and k=-1. This directly tests the type of step used in the original primitive doubling proof. It does not disprove a specially defined nonlocal/twisted boundary R on the actual TFPT source; that identification requires an actual microlocal intertwiner, not just reflection regularity.

The correct candidate primitive is therefore the **relative** pair, rather than their assumed equality:

  `U_rel=(2P_C-I) R_geom`.

It is unitary if both factors are self-adjoint involutions.

## IV. New exact *conditional* cross-ratio/relative-polarization dictionary

Let `t=tan(delta/2)>0`, `R=sigma3`, `X=sigma1`, and

  `C_t = ((1-t²)/(1+t²)) R + (2t/(1+t²)) X`,
  `P_t=(I+C_t)/2`, `U_t=C_t R`.

The projectors are orthogonal and `U_t` is unitary. Its eigenphases are `u,u^(-1)` where

  `u=(1+i t)/(1-i t)`.

For the *same* original antipodal four points `{±1,±u}`, the ordered cross-ratio is `lambda=-t²`. Exactly:

  `lambda=(tr(U_t)-2)/(tr(U_t)+2)`.

For positive t, `tr U_t=0 iff t=1 iff lambda=-1 iff U_t²=-I`. The original v512 nonsquare witness `t=1/2` produces `u=(3+4i)/5` and `tr U_t=6/5`.

This gives a remarkably small **candidate** bridge from geometric modulus to relative Cauchy/deck polarization. It is NOT a derivation of `P_t` as the actual Calderon projector of the original full TFPT Dirac operator, or a reason why its trace must be zero. Even at `U_t²=-I`, the conjugation action has square identity; the genuine E8 current order-four grading still requires the separate charged-discriminant construction already audited in OriginLab.

## V. The missing dynamical source is irreducible

The one-sided primitive datum still *assumes* a smooth manifold with boundary, the spectral triple and its Dirac collar; the operational seed includes local net, state and time. It therefore cannot be advertised as a spacetime-independent theorem deriving a Lorentzian 3+1D parent from nothing. The native source still has no complete weighted four-mark tensor or absolute empty reference. Exact finite cohomology, determinant gluing and relative-polarization matching do not supply the interacting Higgs/gauge/fermionic/BRS-T/eta source.

The measurable next original-source gates are:

1. Construct the **actual** Calderon projector and geometric reflection on the same original positive-time boundary Hilbert quotient, with spectral dependence. Evaluate their relative unitary, its trace/holonomy, and genuine reflection-positive Gram at two blocked levels.
2. Evaluate the original full weighted source at `delta` and `pi-delta` and compute the **defect-dependent residual** after free local gluing, including actual Pfaffian/ghost contributions. Test if the square is uniquely selected without imposing it or choosing a free Gram.
3. On that **same** source, compute `Z0,ZU,ZY,ZUY`, the oriented charged E8 Ward return, physical positive transfer and absolute empty reference. Compare against OriginLab's fibre pushforward and composed-search withheld words.
4. Derive locality, 3+1D and a universal stress/Spin-2 interaction; do not call the 2D geometric character lanes physical space or physical color by numerical matching.

The source-zero-mode block and native R113 are still **BLOCKED**; no amount of synthetic finite control observations closes F1. The next primitive computation must be a real evaluated matrix/field integral from the original owning source.

**Executable finite controls:** [p1_source_falsifiers.py](p1_source_falsifiers.py), [test_p1_source_falsifiers.py](test_p1_source_falsifiers.py). Their source identity is explicitly false. All displayed analytic formulas and independent numeric checks were verified separately; full repository CI is not claimed.
