# F1 original theorem audit: the upstream measure and RG source cannot yet close

**15 November 2026. A rigorous *source theorem* correction, not another invented P1 modulus candidate or physical F1 solution.** This note reopens the owning original TFPT 4.5 metrology LaTeX and provides four independent exact counterexamples, plus one valid nonphysical RG repair. The original all-field Q0 and physical R113 packet remain uncomputed, as do Lorentzian 3+1D locality/chirality and spin-two gravity.

## Original owning statements

Original: [TFPT 4.5, metrology source, section 005](../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex), theorems **Unique admissible UV profile fixed point** and **Constructive geometric measure on the admissible branch**.

The first theorem explicitly declares:
- a cone of finite positive measures on the unbounded interval [0,infinity);
- the normalized transfer density T nu(dt) = [(int_0^infinity K(s,t) dnu(s))/(int_0^infinity psi(s)dnu(s))] dt;
- 0<m<=K(s,t)<=M uniformly on the entire positive quadrant, and psi(s)>0;
- a Hilbert-projective image diameter <=log(M/m), contraction <=tanh(log(M/m)/4), preservation of low-curvature Einstein/Higgs/R² moments and a unique finite fixed Bernstein measure.

The second theorem defines finite-block geometric weights proportional to exp[-Tr f*(D_V²/Lambda²)] dnu_V, and asserts projective consistency and reflection positivity from positive f*, reflection-regular compatible Dirac blocks and existing admissible structures.

Those implications require added analytic conditions. As written, they do not establish the physical source.

## A. The literal original RG map is *never finite* on its declared unbounded target

Take ANY finite positive nonzero nu with 0<d=int psi(s)dnu(s)<infinity. By the theorem's own uniform kernel lower bound,

    (T nu)(dt)/dt >= m nu([0,infinity))/d > 0

for every t>=0. The total mass of the output is at least that strictly positive constant times int_0^infinity dt, and is therefore infinity. If d is infinite or zero, the displayed normalizing quotient is undefined. Consequently **there is no nonzero finite fixed measure at all** under the stated assumptions.

This is not a subtle numerical effect or failed implementation: it is an immediate analytic contradiction. The original "unique finite UV profile" theorem is false (its nonzero finite fixed point cannot exist) as written on [0,infinity) with Lebesgue dt and global m>0.

### A mathematically correct type of repair

Provide a positive **integrable reference envelope** w(t) dt, e.g. w(t)=exp(-t), and assume

    0<m*w(t) <= K(s,t) <= M*w(t)

instead of 0<m<=K<=M. Require the correct normalizer psi(s)=int K(s,t) dt; positivity alone of psi is insufficient. A moment-matched branch slice must also be proved invariant by T, not assumed from positivity. The Hilbert-projective contraction of a strictly positive bounded-ratio kernel can then be estimated on an appropriate complete normalized cone.

An explicit *toy*, not original TFPT, repair is:

    K(s,t)=exp(-t)*(1+exp(-s-t)),
    psi(s)=1+exp(-s)/2.

For probability measures define u=int exp(-s) dnu(s). The normalized update acts on this moment as

    u_new=(1/2+u/3)/(1+u/2),     |d u_new/du|<=1/12 for u in [0,1].

Its unique positive fixed moment is

    u*=(sqrt(13)-2)/3.

The unique stationary probability density is exp(-t)*(1+u*exp(-t))/(1+u*/2), and its unique Bernstein profile is

    f*(z) = [1/(1+z)+u*/(2+z)]/[1+u*/2].

This proves that a **nonvacuous** UV profile fixed-point theorem is possible with corrected hypotheses. But different allowed positive integrable kernels produce different fixed profiles; the ORIGINAL physical kernel remains a separate, missing source law. The repaired toy has no Einstein/Higgs/R² moment assignment, no source-identical R113 and no 3+1D dynamics.

## B. The original Birkhoff projective diameter and contraction bound are numerically false

A positive 2x2 kernel with lower bound m=1 and upper bound M=2 is

    K = [[2,1],[1,2]].

Its Hilbert-projective diameter is exactly log4=2log(M/m), not <=log2. The local projective contraction at the positive identity ray is exactly (M-m)/(M+m)=1/3. The published claimed upper bound is

    tanh(log(M/m)/4)=tanh(log2/4)=3-2sqrt(2)≈0.1715729,

strictly BELOW the actually achieved contraction 1/3. Thus the stated bound fails even for a finite positive matrix.

Correct Birkhoff theorem: eta(K)<=tanh(Delta(K)/4), with Delta the ACTUAL projective diameter. The generic m<=K<=M envelope yields Delta<=2log(M/m), hence eta<=tanh(log(M/m)/2), not the claimed smaller bound; sharper bounds require additional cross-ratio or stochastic-kernel structure. See e.g. [Birkhoff–Hopf and exact positive matrix diameter](https://link.springer.com/article/10.1007/s10915-021-01524-w).

## C. Positive and compressed Dirac spectra do not force consistent finite blocked measures

Choose an exact positive completely monotone Laplace/Bernstein profile

    f(z)=1/(1+z)=int_0^infinity exp(-u) exp(-zu) du.

Use a two-value coarse field x in {0,1} and an additional fine field y in {0,1}. Put

    D_coarse(x)=[x],
    D_fine(x,y)=[[x,y],[y,0]].

Both are self-adjoint, and **the coarse Dirac is the exact compression of the fine Dirac for every y**. Compute the actual spectral actions Tr f(D²):

| x | y=0 | y=1 |
| --- | --- | --- |
| 0 | 2 | 1 |
| 1 | 3/2 | 1 |

The coarse spectral actions for x=0,1 are 1 and 1/2. Put q=exp(-1/2) in (0,1). The ratio of the fine field weights after summing over y is

    R_fine=[exp(-3/2)+exp(-1)]/[exp(-2)+exp(-1)]
          =(1+q)/(1+q²).

The coarse ratio obtained from exp[-Tr f(D_coarse²)] is R_coarse=1/q.

    R_fine-R_coarse=(q-1)/(q*(1+q²)) != 0.

Since even the relative probabilities differ, no overall normalization repairs the failure. **Commutation of blocking with the full physical spectral exponential/measure pushforward is an independent condition**, not a consequence of compatible Dirac compressions and positive f.

This counterexample is a finite source model satisfying those *particular* premises, not an instance of all original TFPT geometric/BRST restrictions. It directly refutes the logical implication that positivity plus Dirac compression by itself yields the project's all-field projective consistency.

## D. Pointwise positive spectral weights are insufficient for OS reflection positivity

On two values of a positive-time field x and its reflected negative-time field y, let reflection exchange x and y. Use the same completely monotone positive f and the scalar reflected Dirac D(x,y)=|x-y|. Then

    mu(x,y)=exp[-f(D(x,y)^2)]

is **pointwise positive** and invariant under reflection, but its OS kernel is

    G=[[q²,q],[q,q²]], q=exp(-1/2) in (0,1).

Its eigenvalues are q²+q>0 and q²-q<0, so G is NOT positive semidefinite. The determinant q²(q²-1)<0. Thus the original proof's statement that reflected OS forms stay nonnegative merely because the new spectral weight is positive is **not a valid general argument**. The full graded reflection positivity must be checked on the actual source kernels, or a stronger structure theorem for local conditional products must be provided.

This is an exact two-state counterexample for the stated weak implication, not a physical TFPT fermion covariance or a counterexample to a stronger, separately demonstrated reflection-positive Markov construction.

## E. First few fixed low-curvature moments do not select a unique Bernstein measure

On the support {0,1,2,3,4}, the distinct positive probability measures

    nu_plus  = (13,8,18,8,13)/60,
    nu_minus = (11,16,6,16,11)/60

have **identical first four moments** k=0,1,2,3 (all checked exactly), because their difference is the fourth finite difference vector (1,-4,6,-4,1)/30. Yet their Laplace profiles differ:

    f_plus(z)-f_minus(z)=(1-exp(-z))^4/30.

Hence moment matching of several Einstein/Higgs/R² coefficients cannot by itself identify the full UV profile. An actual fixed physically derived kernel CAN provide further information, but its existence/uniqueness must first be established with correct source law and a proven invariant moment slice. A free choice of a positive RG kernel cannot close the original F1.

## F. Minimal original-source repair requirements, NOT new synthetic solutions

1. **Derive an integrable positive kernel from the original microscopic one-sided source**, its exact normalizer psi(s)=int K(s,t)dt, and a well-defined source reference measure on the full noncompact domain. Prove the source's *actual* low-curvature moment slice is invariant.
2. Apply the correct Hilbert-projective diameter/cross-ratio bound. Establish a fixed normalized Bernstein measure and the source's actual uniqueness in a complete admissible cone. The explicit nonphysical kernel above is only a unit test for the repaired theorem machinery.
3. At each blocked level, prove **the full all-field weighted measure** has its own (bosonic, ghost, charged fermionic, eta, P_adm) reflection-positive source Gram and exact pushforward. Do not infer either from the positivity of f(D²) or compression of D alone.
4. Only then construct the original Q0 absolute marked form and evaluate the four physical R113 numbers Z0,ZU,ZY,ZUY on the same native source and empty reference at two blocked levels. Next obtain physically non-assumed time, an interacting local 3+1D chiral parent and universal spin-two dynamics. The original one-sided boundary T_bound ALREADY assumes a smooth manifold with boundary and a Dirac collar, so a spacetime-independent F1 rule still lies upstream.
5. Scientific status: **F1 OPEN**, no published claims or Atlas ledger promotion. These concrete errors mean earlier statements about constructive UV uniqueness, source measure compatibility and inherited OS positivity **cannot be considered proved as written**.

Exact executable reproduction: [f1_source_theorem_audit.py](f1_source_theorem_audit.py) and [test_f1_source_theorem_audit.py](test_f1_source_theorem_audit.py), integrated in finite OriginLab and main-branch Core Verification. An independent local SymPy implementation passed eight tests. The model is source-conditional/finite, does not compute the missing original Q0, and is not an alternative TOE.
