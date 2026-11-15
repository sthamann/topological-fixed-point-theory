# F1 original P1 square selection: chiral E8 saddle and the missing source stiffness

**15 November 2026. Exact chiral/modular no-go and positive-source counterexamples. F1 remains OPEN.**
This is the next decisive source-first negative/decision result after moving-divisor correction, the 5+3 E8 Weyl character theorem, and the original full-gluing cancellation. It does not pretend the known mathematical E8_1 VOA partition is the original physical four-mark R113 source.

## 1. Actual owning original source and independent literature

* [Original square and nonsquare RP seam v512](../../verification/v512_seam_tau_flag.py): the square P1 branch is τ=i, with moving antipodal four-mark parameter t=tan(delta/2), t=1 square. The nonsquare V4 family is already an exact counterexample to selecting the square by mark transitivity alone.
* [True R113 whole-source and bound empty-reference contract](../../experiments/theory-contracts/chat-triangulation-20261108/R113-ORIGINAL-UND-REFERENZ-20261109.txt): the source is Q[b]=s[b] Gamma(R[b]), the original zero-fermion cap reads only the absolute bosonic scalar contraction I(w s), not the synthetic full fermion trace; the ordered Z0,ZU,ZY,ZUY all remain blocked. Its source-score variation produces the mixed third cumulant and requires the actual moving selector/geometry derivatives.
* [Primitive boundary and Dirac paper](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex) and [original metrology](../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex) assume one-sided collar, field measure, RG and admissibility data; the actual fully evaluated native Q0 packet is not present in the owning consumer.
* [Corrected moving-divisor source](P1-MOVING-DIVISOR-PRIMITIVE-CLOSURE-20261115.md) and [native E8 Weyl clock](P1-FIVE-PLUS-THREE-WEYL-E8-20261115.md).
* External standard E8 VOA source: Dong and Mason, *Holomorphic Vertex Operator Algebras of Small Central Charge*, Lemma 2.1, https://citeseerx.ist.psu.edu/document?doi=b71f1d867b7e465aa9d850daccc68ca5fa08b8d2 ; E8 lattice character and modular transformations in https://www.math.ksu.edu/~gerald/voas/voa/lat-E8.html ; the Ramanujan modular derivatives in the PTEP 2021 appendix https://academic.oup.com/ptep/article/2021/8/08B104/6166659 .

## 2. The genuine chiral E8_1 partition does not maximize at the square

The unique chiral E8 level-one vacuum character is the actual E8 lattice VOA character

    chi_E8(tau)=Theta_E8(tau)/eta(tau)^8
               =E4(tau)/eta(tau)^8=j(tau)^(1/3).

For the original square τ=i, the exact modular values give E6(i)=0, j(i)=1728 and **chi_E8(i)=12**. Using the Ramanujan identities for E4, E6, eta yields

    d/dtau log(chi_E8) = -(2*pi*i/3) E6/E4,
    d2/dtau2 log(chi_E8)(i) = -(2*pi^2/3) E4(i).

Since E4(i)>0, chi has a critical point with nonzero SECOND complex derivative, not an absolute maximum. If tau=i+x+iy (real x,y small), its modulus-square logarithm has the real Hessian

    Hessian log|chi_E8|² at tau=i
      = diag(-4*pi²*E4(i)/3, +4*pi²*E4(i)/3).

Thus the square is a **strict nondegenerate saddle** of the modulus (a maximum along horizontal and minimum along vertical). This also follows from the holomorphic maximum-modulus principle for a nonconstant nonzero chiral character in an open complex neighbourhood. The chiral E8 partition is S-invariant at tau=i but modular symmetry does NOT make the modular parameter sit at its fixed point.

High-precision (60 decimal places) deterministic independent evaluation:

    |chi_E8(i)|=12,
    |chi_E8(i+0.01)| ~ 11.99425407389788 (<12),
    |chi_E8(1.01 i)| ~ 12.00569082806205 (>12).

For a *proposed*, source-derived factorization of the physical full measure

    Z_full(tau) = |chi_E8(tau)|² W_rest(tau),

the necessary vertical local-max condition at τ=i would be

    (partial_y² log W_rest)(i) < -4*pi² E4(i)/3
                                   ~ -19.1570717979704.

No actual original TFPT factorization or physical boundary-to-torus tau dictionary has been proven, so this is a *conditional quantitative source discriminator*, not a claimed TFPT prediction. The actual OS source can include additional charged, gauge, ghost, eta, gravitational and geometry contributions; an apparently nice chiral partition is not their replacement.

## 3. OS positivity, mark symmetry and exact common time sewing still DO NOT select square

A more fundamental counterexample does not rely on a chiral character. Let t>0, representing the antipodal marked family with square t=1. Let

    f(t)=4t²/(1+t²)²=sin²(delta),
    W_plus(t)=1+(1/2)f(t),
    W_minus(t)=1-(1/2)f(t).

Both W_plus and W_minus are everywhere positive, invariant under t->1/t and consistent with the same V4 mark geometry. Yet W_plus has a **strict maximum** at t=1, W_minus a **strict minimum** there. Their exact second derivatives of log W are -2/3 and +2, respectively.

For each of these two source weights take the **same** positive time generator H=diag(0,1), vector v=(1,1), and genuine time-sewing semigroup exp(-sH). The reflected Gram evaluated at times s_k=k ln 2 is

    K_ij(t)=W(t) v*exp(-(s_i+s_j)H)v
           =W(t) (1+2^(-i-j)).

This is a positive semidefinite rank-two matrix for every t, since K is a sum of two positive rank-one Grams. Its exact three-time principal minors all pass the OS positivity necessary test. The semigroup sewing law is exact. It follows that the *declared* abstract properties of positive source weights, V4/t-reflection symmetry, and an OS-positive sewn time dynamics do not determine the square selector. These are **synthetic countermodels** and do not reproduce all original field content, physical P, E8 current Ward identities, reference/eta terms or R113 source.

The two reflectively symmetric countermodels produce opposite even-modulus effective actions despite having the *same* conditional time semigroup. Adding more combinatorial model variants cannot resolve their nonidentifiability without real source input.

## 4. Concrete original-source decision: stress curvature versus fluctuations

Set x=log t, so the square is x=0 and the geometric reflection t->1/t becomes x->-x. For an actual *finite fixed-domain weighted source* of the form

    Z(x)=sum_history mu_history exp(-S_history(x))

with x-independent mu_history and no changing projector, the correct exact identity is

    d²log Z/dx² = Var_x(S'_history) - E_x[S''_history].

Reflection of the true original source can make the first derivative vanish at x=0, but **does not fix the sign of the second derivative**. The square is locally selected as a maximum iff

    E_0[S''] > Var_0(S').

A particularly sharp counterexample has two reflection-exchanged histories

    S_±(x)=c x²/2 ± a x.

Then Z=2 exp(-c x²/2)cosh(a x) is positive for every real x, symmetric under x->-x, and

    (log Z)''(0)=a²-c.

The SAME type of symmetric positive source can therefore exhibit a maximum (a=1,c=2), a minimum (a=2,c=2), or a flat quadratic response (a=c=1). The source must supply its *actual* stress mean and fluctuation covariance to decide.

In TFPT the actual original measure, geometry, external projector P_adm, boundary cap, eta phase, and reference may DEPEND on x. These additional derivatives must be included. This finite fixed-domain score formula is a necessary template only; it must NOT be used without the actual original corrections. The existing R113 score/cumulant and product-rule sections already supply the source-level warning. In particular a normalized mixed C_UY can lose an overall scalar modulus factor, whereas absolute Z0 uses the same original bound reference. Never normalize away the decisive source scalar before the modulus test.

## 5. New F1 search architecture and exact next experiment

The search should stop promoting new synthetic positive transfers to source solutions. Candidate acceptance now requires:

1. **Native absolute field source on two blocked levels:** original Q0, same geometry-dependent reference, all fields, P_adm, eta, time, true conditional fiber pushforward, four R113 forms Z0,ZU,ZY,ZUY. No freely chosen relative Gram or arbitrary exponent.
2. **Moving-mark source derivative:** evaluate the score insertion S' and stiffness S'' on the same original regulated history field measure including the moving-domain terms, P'_adm and physical source reference. Test whether the true log-partition Hessian is negative at the square, and, if claiming global selection, compare values away from the square.
3. **Affine/nonchiral decomposition if actually present:** verify whether the correct source really factors into the E8_1 character and a physical additional measure. If so, its vertical curvature must overcome +19.1570718 as above, a stringent independent quantitative necessary condition.
4. **Cubic/current algebra type:** establish physical chiral CAR/OPE weights and source Ward cubic in the original moving-divisor field map before identifying any su3 principal current with the native E8 colour algebra.

New exact OriginLab controls are [p1_modular_source_no_go.py](p1_modular_source_no_go.py) and [p1_modulus_score.py](p1_modulus_score.py) with independent regression tests. They document their own finite-versus-physical provenance and do not claim any original weighted R113 output.

**Conclusion:** source-free modular invariance, chiral E8 positivity, standard OS positivity, and ordinary sewing are individually or jointly insufficient to pick tau=i. The missing quantity is the actual *source-derived modulus stiffness relative to fluctuations*, with correct physical field/map/reference. F1 remains OPEN, and this result is a sharper proof obligation rather than a claim of solved physics.
