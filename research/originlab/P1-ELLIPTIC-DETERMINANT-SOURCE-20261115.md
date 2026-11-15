# TFPT F1: physical P1 elliptic-cover spectral determinant, exact square source and E8 competition

15 November 2026. **A new source-compatible spectral-geometry candidate, NOT a completed source derivation or F1 proof.** Its mathematical core is a known rigorous theorem for unit-area rectangular flat tori, combined with newly computed source bounds in the original TFPT four-mark t parameter.

## Original owners and independently published mathematics

* [Original P1 square and non-square modulus v512](../../verification/v512_seam_tau_flag.py): the antipodal four marked P1 branch D={1,-1,u,-u}, u=(1+i t)/(1-i t), t in (0,1], square t=1, elliptic modular j(-t²). Original nonsquare t=1/2 produces exactly j=148176/25.
* [Exact original elliptic j and E8 character bridge](p1_original_j_witness.py), [chiral E8 nonmaximum / source OS counterexamples](P1-MODULAR-E8-SOURCE-SADDLE-20261115.md), [moving divisor source and original field caveats](P1-MOVING-DIVISOR-PRIMITIVE-CLOSURE-20261115.md).
* [Original one-sided Dirac/boundary theorem](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex), [original full R113 weighted source](../../experiments/theory-contracts/chat-triangulation-20261108/R113-ORIGINAL-UND-REFERENZ-20261109.txt). These DO NOT identify the all-field TFPT geometry source with a unit-area scalar Laplacian determinant on the genus-one double cover.
* Markus Faulhuber, **Extremal determinants of Laplace–Beltrami operators for rectangular tori**, Mathematische Zeitschrift 297 (2021), 175–195, Theorem 1.1, https://doi.org/10.1007/s00209-020-02507-7 . It proves the unique maximum of the unit-area rectangular torus zeta determinant at the square, not a QFT gravitational source or a fermion/ghost exponent.

## 1. The REAL elliptic double-cover parameter from original four marks

The four punctures determine a double cover of P1 branched at those marks: a genus-one elliptic curve. Choose the standard anharmonic cross-ratio

    x=t²/(1+t²), 0<x<=1/2,

so its rectangular modular parameter tau=i*y is obtained from complete elliptic integrals

    y(t)=K(1-x)/K(x), y>=1,

where K(m)=∫_0^(pi/2) (1-m sin² theta)^(-1/2) dtheta is the complete elliptic integral of the first kind in the parameter convention. This parameter is not arbitrarily fitted to a desired modulus. Direct exact symbolic algebra shows

    j(x)=256(1-x+x²)^3/[x²(1-x)²]
        =256(1+t²+t⁴)^3/[t⁴(1+t²)²]

and hence reproduces the ACTUAL v512 full nonsquare witness as well as j(i)=1728. The sign/anharmonic conventions are explicitly fixed.

Numerical source-geometric controls at 60 decimal digits:

| original t | elliptic lambda x | rectangular y | unit-area det' Laplacian |
| --- | ---: | ---: | ---: |
| 1/3 | 1/10 | 1.598874970177602 | 0.299628724136963 |
| 1/2 | 1/5 | 1.360070638131823 | 0.327095257431999 |
| 2/3 | 4/13 | 1.201210027631647 | 0.340722699515684 |
| 1, original square | 1/2 | 1.000000000000000 | 0.348300982421419 |

## 2. A genuine spectral-geometry square selector, with the crucial field-statistics sign

For a flat torus of **unit area** and rectangular modulus tau=i*y, its zeta regularized positive Laplacian determinant (zero mode removed) is

    D(y)=det' Delta_y = y |eta(i*y)|⁴.

The Dedekind eta product gives

    log D(y)=log y - pi*y/3
             + 4 sum_(k>=1) log(1-exp(-2pi*k*y)).

It obeys D(y)=D(1/y) by modular S, and for all y>0,

    d² log D(y)/dy²
      = -1/y² - 16pi² sum_(k>=1) k² exp(-2pi*k*y)
                               / [1-exp(-2pi*k*y)]² < 0.

Therefore log D is strictly concave, its derivative vanishes at y=1 by S symmetry, and D has a **unique maximum** at the square throughout the FULL rectangular unit-area class. This also follows from Faulhuber's published rigorous theorem. It is a legitimate physical-looking source **candidate** because it is a spectral determinant of the actual elliptic double cover, not an extra free fitted numerical model.

The **QFT sign gate** is load-bearing: a standard bosonic Gaussian integral contributes D^(-1/2), which is *minimal*, not maximal, at the square. A net positive fermionic/ghost/superdeterminant power D^m with m>0 would favor it on this class. The physical TFPT BRST/ghost/Pfaffian/eta/Dirac source must compute m; guessing its sign from the existence of 16 Majoranas would be invalid. The earlier full free-circle gluing cancellation is NOT a contradiction: the present D is a global two-dimensional elliptic double-cover spectral determinant, not the four-arc subdivision determinant of one unchanged circle. The new 2D unit-area metric and field boundary dictionary are genuine extra physical source premises.

## 3. Exact quantitative competition with chiral E8 at the square

Suppose, CONDITIONALLY, the true all-field TFPT source admits a factorization at the same physical metric/time and reference:

    Z_full(i*y) = |chi_E8(i*y)|² * D(y)^m

where chi_E8=E4/eta^8 is the standard chiral E8_1 vacuum character. The bare E8 magnitude is MINIMAL along the rectangular family at y=1; it has a positive logarithmic vertical second derivative there:

    kappa_E8 = d²_y log |chi_E8(i*y)|² at y=1
             = 4pi² E4(i)/3
             ~ +19.15707179797043649.

The scalar determinant logarithmic Hessian is, exactly in terms of its convergent positive series,

    kappa_D = -1 -16pi² sum_(n>=1) n² exp(-2pi*n)
                                  /[1-exp(-2pi*n)]²
            ~ -1.29821132491543485.

The full local-square Hessian in this *declared factorization* is kappa_E8 + m*kappa_D. It becomes negative iff

    m > kappa_E8/(-kappa_D)
      ~ 14.75651261878979767523.

Thus the minimum positive INTEGER determinant exponent giving a **local maximum** along the imaginary rectangular direction in this particular multiplicative source model is exactly **m=15**. At m=14 the local square remains a minimum; at m=15 its vertical curvature becomes negative. A genuine two-complex-dimensional maximum would need independent horizontal Hessian and Hessian mixed-term checks; global maximum across all moduli still needs full comparison (and within general tori the hexagonal determinant is the global maximizer, not the square).

This 15 is NOT a prediction of 15 physical particles, 15 ghosts or a derived Majorana multiplicity; it counts a model's NET POWER of one normalized scalar Laplacian determinant. No current TFPT original physical source proof derives that factorization or its exponent. The numerical threshold is a precisely bounded F1 falsifier for any future source candidate in this restricted class, not a new verified TOE observable.

## 4. Original physical source / R113 mandatory checks

The remaining real source must decide all of:

1. Does original P1 boundary kernel really define a unit-area elliptic double-cover Laplacian (or related Dirac) on its physical BRST/OS quotient with these exact geometric mark coordinates, time and absolute reference? v512's analytic elliptic j by itself is insufficient.
2. Compute the complete original bosonic Gaussian, gauge-fixing ghost, charged E8 spinor, eta/Pfaffian and spatial/gravitational factor under moving four marks and two blocked levels. Determine the ACTUAL determinant factors, their net sign and exponent before assigning a preferred square.
3. At tau=i and t=1/3,1/2,2/3 compute physical unnormalized Z0 and ZU/ZY/ZUY with the same bound original reference. Use the original R113 ordered operator words; do not replace the empty exterior cap by a full fermion trace or assign a geometric determinant to a synthetic source.
4. Test local stiffness versus fluctuations, the full second variation with actual moving measure/selector/reference/eta, plus the global modulus comparison. Then derive actual positive time, interacting chiral 3+1D parent and universal spin-two coupling from the same source.

**F1 remains OPEN; this is an exact CONDITIONAL square selector and sharp original-source exponent test, not a completed P1 physical origin.**

Exact reproducible code: [p1_elliptic_determinant_source.py](p1_elliptic_determinant_source.py), [test_p1_elliptic_determinant_source.py](test_p1_elliptic_determinant_source.py). The native sources are not replaced or extended with fabricated numerical amplitudes.
