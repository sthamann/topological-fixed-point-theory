# F1: The missing moving-divisor source transport and the primitive four-mark operator theorem

**15 November 2026 · Corrected source interpretation, two new conditional four-mark theorems, exact countermodels.** This report explicitly corrects a too-strong physical reading of the earlier [principal SU3 candidate](P1-PRINCIPAL-SU3-ORIGIN-20261115.md). It does not produce the native weighted Q0, original R113 values, physical time, 3+1D spacetime or gravity. **F1 OPEN.**

## 1. Which original operations were checked?

The four original regional parabolic vector fields in [source_bridge.py](source_bridge.py) are x², (x−1)²/2, 1 and (x+1)²/2. They span precisely the three noncompact-real Möbius vector fields 1, x, x² with Lie brackets [1,x]=1, [1,x²]=2x and [x,x²]=x².

The original four-mark section and logarithmic-differential spaces in [P1-EQUIVARIANT-RESIDUE-ORIGIN-20261115.md](P1-EQUIVARIANT-RESIDUE-ORIGIN-20261115.md) are C5=H0(P1,O(D)) and F3=H0(P1,Omega¹(log D)) for D=mu4. Their bases are f_k=z^k/(z⁴−1), k=0..4, and omega_j=z^j dz/(z⁴−1), j=0..2.

The original one-sided TFPT boundary theorem [source](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex) still assumes a smooth oriented manifold with boundary and Dirac collar. The native weighted field source [source.py](source.py) still does not supply the physical same-reference R113 packet. These premises must not be silently replaced by an abstract Lie-algebra identity.

## 2. Hard correction: no infinitesimal Möbius action on one FIXED four-point divisor

The ordinary Lie derivative of f0=1/(z⁴−1) along v=1 is −4z³/(z⁴−1)², which has a double pole at z=1 with leading coefficient −1/4. C5 allows only simple poles. The same failure holds for logarithmic one-forms. More generally a holomorphic PSL2 vector field v(z)=a+bz+cz² that preserves each of four distinct points infinitesimally must vanish there; a quadratic polynomial with four zeros is identically zero. Thus the fixed D=mu4 has **no nonzero infinitesimal PSL2 stabilizer**. Its actual stabilizer is discrete D4.

The original regional sl2 cannot simply be identified as a positive internal continuous symmetry acting on the fixed C5/F3 eight-dimensional field space. The earlier principal SO3/SU3 construction is mathematically valid as a representation *after a frame and real-form choice* but was not proven to be an action of the ORIGINAL Q0 source at a fixed mark geometry.

## 3. Exact replacement: move the divisor and transport BOTH typed spaces

Let D(t)={a_j(t)} and P_t(z)=product_j(z−a_j(t)). For v(z)=a+bz+cz² let a'_j(t)=v(a_j(t)). At the original centered D=mu4, sum a_j=0, differentiation of the MOVING family gives the shared numerator correction

    sum_j (v(a_j)−v(z))/(z−a_j) = −4b−4cz.

Then the corrected transport on meromorphic sections f_k=z^k/P and logarithmic differentials omega_j=z^j dz/P has polynomial numerators

    C5: v h' + (−4b−4cz)h,
    F3: v g' + v'g + (−4b−4cz)g.

These expressions stay inside polynomials of degrees <=4 and <=2, respectively. The raw fixed-divisor Lie derivatives do not.

After the **same** scalar normalization of the dilation transport (+2 I on both typed spaces at this centered divisor), the two spaces carry the standard irreducible complex sl2 polynomial modules

    C5 ≅ Sym⁴(C²): E=partial_z, F=z²partial_z−4z, H=zpartial_z−2,
    F3 ≅ Sym²(C²): E=partial_z, F=z²partial_z−2z, H=zpartial_z−1.

Every sl2 commutator is checked exactly. The choice of scalar normalization and source-positive real structure is not thereby derived from the physical theory.

For arbitrary four points a,b,c,d, their cross-ratio obeys the exact moving-divisor Ward identity

    sum_j v(a_j) d lambda/d a_j = 0

for each v in {1,z,z²}. Hence **regional Möbius covariance CANNOT by itself select the square cross ratio**: it preserves every modulus. Additional original puncture/defect/charged source data must choose lambda=-1.

## 4. A second correction: the canonical 3->5 exact sequence is D4, not full sl2 equivariant

The original exact evaluation inclusion i:F3 -> C5 sends omega_j -> f_(j+1) using dz/z. With the standard moving-divisor sl2 representation,

    H5 i = i H3,
    E5 i − i E3 !=0,    F5 i − i F3 !=0.

The two residuals are explicit rank-three matrices; the original log-to-meromorphic inclusion preserves weights and finite quarter action but **does not** commute with full sl2 raising and lowering. This is mandatory by Schur's lemma: Sym⁴ and Sym² are inequivalent irreducibles.

Likewise the geometric compensator f2=z²/(z⁴−1) is a D4 sign character but not a full sl2-invariant line: E f2=2f1 and F f2=−2f3. There are no nontrivial one-dimensional sl2 characters. The earlier f2-marked principal SU3 cubic is a **finite D4 selection rule**; it is not an automatically continuous regional Ward invariant. Its physical source insertion must explicitly be identified.

## 5. Positive source real form is a new gate, not an automatic consequence of the real Möbius algebra

The actual regional vector fields have split-real sl2(R) coefficients. Any nontrivial finite-dimensional representation has a dilation generator H with nonzero REAL eigenvalues. If a strictly positive fixed-fiber Hermitian metric G were invariant under the full split-real action, it would require H*G+GH=0, contradicting positivity on those eigenvectors. More generally the only finite-dimensional unitary representations of SL2(R) are trivial. See the standard linear Lie-group representation notes at https://tomasz.przebinda.com/LieTheoryII_Spring_2020.pdf, Theorem 72.

A separate **compact SU2 real form** of the complexified sl2 does admit positive metrics on the polynomial modules:

    G_C5 = diag(1,1/4,1/6,1/4,1),
    G_F3 = diag(1,1/2,1).

Their normalized generators satisfy E*G=−G F and H*G=G H, and these metrics preserve the original finite square D4 action. BUT choosing this compact real form, and identifying it with the original reflected OS field pairing, are SOURCE assumptions, not theorems from the real regional vector fields. Even a compact SU2-covariant family across moving divisors need not force a fixed-fiber cross-Gram mixing z to vanish: a covariant field metric may depend on the marks and has only D4 stabilizer at the square.

An interesting, still conditional unification is to use ONE complex sl2(C) and identify its split versus compact real forms through the two different original reflection/star operations. That would distinguish geometrical transport from positive internal currents without falsely trying to make a noncompact finite state space unitary. This real-form and OS bridge remains unconstructed in the original source.

## 6. A precise four-mark theorem from internal Lie/Jacobi closure, without E8 or fermions

For general n>=2 marks, the logarithmic forms form the complex sl2 module

    F_n ≅ Sym^(n−2)(C²),  dim F_n=n−1.

**Conditional primitive-source axiom:** The physical source induces a nonzero, alternating, sl2-equivariant bilinear bracket [ , ]:wedge² F_n -> F_n satisfying Jacobi, using no additional current fields.

**Theorem:** Under this condition n=4.

Proof: The derived space [F_n,F_n] is a nonzero sl2-submodule of an irreducible module, hence is all F_n. The finite-dimensional Lie algebra F_n is perfect. Its solvable radical is characteristic and thus sl2 stable, so irreducibility forces radical=0 or all. The latter is incompatible with a perfect solvable algebra. Consequently F_n is semisimple. All its derivations are inner. Its nontrivial sl2 action therefore embeds sl2 into F_n and furnishes an invariant 3D adjoint submodule. Irreducibility implies F_n itself has dimension 3. Thus n−1=3 and n=4. Conversely Sym²(C²) is the adjoint sl2 representation with the nonzero ordinary Lie bracket.

The classic *third binary transvectant*

    (f,g)_3 = sum_{k=0}^3 (-1)^k binom(3,k)
                (partial_x^(3-k)partial_y^k f)
                (partial_x^k partial_y^(3-k) g)

gives an explicitly nonzero alternating sl2-intertwiner wedge² Sym⁴ -> Sym², the unique five-five-three structure up to scale. All 75 binary-quartic equivariance checks pass; on (x⁴,y⁴) its value is 576xy.

The first higher candidate n=8 has F_8=Sym⁶ of dimension 7 and does admit a skew equivariant third transvectant wedge² Sym⁶ -> Sym⁶, but **fails Jacobi**. With (f,g,h)=(x⁶,x⁵y,xy⁵), its Jacobi sum is the nonzero polynomial 5,184,000 x⁶. This is an exact counterexample, consistent with the general theorem.

**Scope:** Jacobi selects n=4 only IF the actual original family log-mode source carries the specified nontrivial equivariant INTERNAL bracket; the regional sl2 action must first be constructed across moving source configurations, not simply placed on a fixed D. It does not select the square modulus.

## 7. Another conditional n=4 theorem: no hidden operators in one quantum algebra

Let the family differential space F_n be a candidate quantum state Hilbert space, with all observables in its unital associative full matrix algebra End(F_n). Require that the original two geometric readouts together with the scalar identity are COMPLETE, with no extra operator directions:

    End(F_n) = C*I direct_sum F_n direct_sum C_n

as typed representation spaces. Dimensions require

    (n−1)^2=1+(n−1)+(n+1) = 2n+1
    => n(n−4)=0 => n=4 for n>=2.

At n=4, the usual spin-one representation of compact SU2 has

    End(C³) = spin0 (dim1) + spin1 (dim3) + spin2 (dim5).

There is an explicit positive unital full 3x3 matrix algebra. Its skew-Hermitian traceless commutator is compact su3; all 81 matrix products close, and the Hermitian Hilbert–Schmidt metric is positive. This is a stronger *associative* source completion than matching the number eight.

It still requires the actual original physical source to identify the log differentials with one Hilbert 3 and the meromorphic sections with its remaining spin-two observable operators (with the correct D4 sign-line twist). It is not an automatic consequence of Riemann–Roch.

## 8. Important obstruction to conflating physical fermions and those su3 operators

If all eight typed geometric modes are INDEPENDENT complex chiral CAR fermions transforming in the adjoint SU3 representation, their bilinear SU3 current algebra has level 6 in the standard convention where the fundamental representation has index one. The exact 64 generator-pair trace identities give

    Tr_ad(T_a T_b) = 6 Tr_3(T_a T_b).

A regular A2 root subalgebra of the native E8 level-one currents instead has level 1 (root length and embedding index one). Thus the *adjoint eight-complex-fermion* source cannot simply be relabelled the same pinned native SU3 level-one colour current. Extra coset/higher-index embeddings or a different representation assignment could change the source interpretation; this is a scoped exclusion only.

See standard free-fermion current-algebra level/index relation in the Cambridge nonperturbative-field-theory chapter https://www.cambridge.org/core/services/aop-cambridge-core/content/view/4558A6AAC013339E351F7C39D979B907/9781009401647c3_39-60.pdf/theories_invariant_under_affine_current_algebras.pdf .

## 9. What the actual F1 search must now compute

The new architecture should distinguish three **typed** operations rather than identifying them:
1. Mark-divisor changes D -> gD with physical source configuration functor, moving-field transport and same-reference positive Gram. A fixed D only has discrete D4; no SL2 noncompact unitary fixed-carrier representation.
2. A source star/reflection selecting a compatible positive compact internal real form, physical operator/OPE type (CAR fermion h=1/2, compact su3 currents h=1, or both in distinct layers) and, if appropriate, the original marked 3->5 operator and third-transvectant vertex.
3. Actual weighted original four-mark Q0, Z0/ZU/ZY/ZUY, the absolute empty cap, original physical time and cross-ratio dependence, before regulator limits, preserving the original full-field BRST/fermion/gauge/boson integration at TWO block scales.

Do not infer actual source values from the finite models. The existing OriginLab composed search should use these as **necessary type/no-go gates**: fixed-divisor sl2 and positive unitary split sl2 on finite modules must be rejected; categorical moving transport is permitted, but it MUST carry actual original source data. An unmarked cubic and an f2-marked cubic are physically distinct operations.

**F1 remains open.** The underlying one-sided TFPT input already assumes smooth boundary and Dirac collar; this note has not generated spacetime or universal spin-2 dynamics.

### Reproduction

Executable controls: [p1_moving_divisor_origin.py](p1_moving_divisor_origin.py), [p1_operator_closure.py](p1_operator_closure.py) and independent unit tests of both. The same exact identities were independently run in a standalone SymPy audit before they were committed.
