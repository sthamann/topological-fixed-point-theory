# Four primitive marks from the first local weight-one spinor current extension

**15 November 2026. Exact CONDITIONAL arithmetic and root-lattice theorem, not physical F1.** It is a new alternative to the earlier Gauss-Bonnet n=4 proof and connects the original n-point Riemann-Roch geometry to the *necessity* of a nontrivial weight-one Spin extension.

## Inputs and original owners

For n>=2 distinct marks D on P1, Riemann-Roch and logarithmic differential theory give

    dim H0(P1,O(D)) = n+1
    dim H0(P1,Omega1(log D)) = n-1

independently of their moduli. This is the n-mark generalization of the exact original TFPT v189, v197 and v228 results. The existing P1 source gives n=4 through a separate **conditional** flat Z2 branch Gauss-Bonnet argument (v216); this note does not reuse that mark count.

Define the typed sum V_n=H0(O(D)) direct_sum H0(Omega1(log D)), dim_C V_n=2n.

**NEW PHYSICAL PREMISE, NEVER CLAIMED DERIVED:** the full original positive-time source realizes ALL 2n typed modes as independent positive chiral complex fermion CAR fields of conformal weight 1/2, with no Hodge/BRST identification collapsing their ranks. This maps V_n into 4n *real* Majorana fields. Their chiral free current algebra is so(4n)_1, with c=(4n)/2=2n, and the even spinor representation has conformal weight h_s=(4n)/16=n/4 and ground-state dimension 2^(2n-1).

## Exact locality and weight-one current selection theorem

The so(4n) root lattice is D_(2n), of rank 2n, discriminant 4. A representative positive half-spinor weight is s=(1/2,...,1/2) in R^(2n) with norm squared

    s² = (2n)/4 = n/2,

so the corresponding chiral conformal weight is h_s=s²/2=n/4.

The D_(2n) spinor lattice glue has index two and determinant one. The extension is an **even** lattice, hence an ordinary local BOSONIC spinor simple-current extension, exactly when the spinor norm squared is an even integer: n divisible by 4. This is equivalent to the integral conformal spin requirement h_s in Z.

If we impose not merely a local extension but the **first/primitive spinor current extension at conformal weight one**, then

    h_s=1 <=> n=4.

This is the UNIQUE n>=2 at which new nontrivial half-spinor fields can enter the **affine weight-one current algebra**. At n=8 (32 Majoranas), a bosonic spinor extension is possible at weight two, but introduces no additional weight-one currents; its current algebra remains so(32), dimension 496 (D16+ lattice). Thus n=4 is NOT selected merely by integer conformal spin; the **weight-one / current primitivity** criterion is essential and must be justified physically rather than assumed by fiat.

At n=4, 2n=8 complex fermion modes, hence 16 real Majoranas, c=8 and D8 root lattice with 112 norm-two roots. The even spinor weight has squared norm two and exactly 128 states, each a new norm-two E8 root. This gives

    112 D8 roots + 128 spinor roots = 240 E8 roots,
    240 + 8 Cartan currents = 248 current generators.

The root system is the E8 even index-two unimodular D8+ lattice, without putting E8 into the n input. The older exact E8 simple-current/glue construction remains the standard source; this is an independent **mark-count origin** under the stated fermionization rule.

## Why this is valuable but not complete F1

The result selects **four marked POINTS** from a precise minimal source-class rule, not the special square cross ratio lambda=-1. The original v512 generic antipodal configurations all remain admissible without the additional flag/quarter source condition.

The typed 5+3 spaces are geometrically related by the exact orientation-twisted evaluation sequence. Their actual physical CAR independence, positive native weighted Gram and source field weights are unproved. The new [geometric Gram gate](P1-SEAM-DIRAC-GRAM-20261115.md) restricts this to |z|<1 under normalized D4 source assumptions, and if the original physical Gamma is identified with form-degree parity and is a genuine Gram isometry, then z=0. Neither premise has been checked on the native Q0.

This note does NOT derive original spacetime 3+1, physical time, gravity, current chirality in an interacting Standard Model, real P1 c3, or actual R113 Z0/ZU/ZY/ZUY. It does not claim locality of the original source merely because a standard bosonic lattice VOA extension can be constructed mathematically.

## Reproduction and kill conditions

Exact arithmetic source: [p1_minimal_spinor_marks.py](p1_minimal_spinor_marks.py), [test_p1_minimal_spinor_marks.py](test_p1_minimal_spinor_marks.py). The census checks n=2..32, spinor weights n/4, local evenness iff n divisible4, the weight-one uniqueness n=4 and the 240 exact E8 norm-two roots. This is scoped verification of a conditional model, not native R113 evaluation.

Reject this as a PHYSICAL origin if the original full source identifies the logarithmic modes with the five meromorphic modes after BRST/Hodge reduction, assigns non-1/2 conformal field weights, or selects a spin structure that forbids the weight-one halfspinor extension. The original marked free-collar time, orientation and actual source integration must decide. **F1 remains open.**
