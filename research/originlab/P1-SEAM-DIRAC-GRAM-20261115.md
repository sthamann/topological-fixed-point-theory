# P1 geometric source candidate: an 8-mode odd seam operator and exact 8-vs-6 CAR Gram gate

15 November 2026. This is a finite **constructive and falsifying** continuation of the exact fourmark 5+3 evaluation sequence, not an original physical P1 or R113 construction. The exact root-lattice E8 Weyl theorem is in [P1-FIVE-PLUS-THREE-WEYL-E8-20261115.md](P1-FIVE-PLUS-THREE-WEYL-E8-20261115.md).

## Original source geometry

C5=H0(P1,O(mu4)) has basis f0,...,f4 and F3=H0(P1,Omega1(log mu4)) has basis omega1,...,omega3. The exact D4-equivariant inclusion i:F3 -> C5 tensor orientation maps omega_j -> f_j and its cokernel is the two endpoint evaluation channels at z=0 and infinity. Rotation r:z->iz gives degrees C5=(0,1,2,3,0), F3=(1,2,3); inversion s:z->1/z acts as s5=-rev5, sF=+rev3, so s5*i=-i*sF.

## One natural finite *odd operator* from the same exact sequence

On the typed sum C5 direct_sum F3 define

    D = [[0_5,i],[i^*,0_3]]
    Gamma = diag(I5,-I3)
    H(m)=D+m*Gamma.

With the declared standard positive D4-invariant Gram and the injection i as isometry,

    spec(D) = { -1 (x3), 0 (x2), +1 (x3) }.
    spec(H(m)) = { -sqrt(1+m²) (x3), m (x2), +sqrt(1+m²) (x3) }.
    det H(m) = -m²(1+m²)^3.

The two zero modes are precisely the geometric quarter fixed-point evaluation channels (0,infinity). The original geometric D4 rotation commutes with D, while full inversion S anticommutes with D, and obeys

    S H(m) S = -H(-m).

So the doubled finite operator has two signed zero-mode spectral crossings at m=0, i.e. an **index 2**. Restricting to a specified single endpoint gives one crossing, i.e. index 1, and the geometric inversion exchanges the two one-sided choices. The one-sided endpoint projection is an *additional physical cap assumption*, not a theorem about the actual primitive P1 kernel.

The positive Hodge-type generator D² on the full 8D typed space has spec {0(x2),1(x6)}; on either one-sided 7D complement it has {0(x1),1(x6)}. Hence a **finite** candidate semigroup has a unique vacuum and gap one:

    T_side(t) = P_vac + exp(-t)*(I-P_vac).

The reflection exchanges two such finite semigroups. This establishes a mathematically explicit "side selection then positive dynamics" model, but NOT the actual original OS time, its units, common physical state, spacetime locality, or a 3+1D parent. It is a finite 8x8 matrix, not the original first-order elliptic Dirac collar or its Calderon projector.

The geometric r and s have determinant -1 on each separate C5 and F3 sector, but determinant +1 on their eight-dimensional typed direct sum. Thus the combined D4 representation is special-unitary and can be represented on its exterior fermion space without an independent determinant-phase obstruction. This is a *geometric determinant-character cancellation*, NOT a proven physical gauge or chiral anomaly cancellation.

## An exact physical CAR-rank gate from the original D4 geometry

Assume, CONDITIONALLY, that the original physical source produces positive normalized Hermitian forms on C5 and F3, with the same geometric linear D4 actions. A mixed source Gram block M:C3 -> C5 must satisfy

    R5^* M RF=M,   S5^* M SF=M.

Solving all 15 complex coefficients shows the mixed intertwiner space has **complex dimension one**, represented in the original bases by

    M[1,0]=z, M[3,2]=-z, all other entries zero.

Crucially the middle character (-1) has no permitted mixed pairing, because its reflection parity disagrees between C5 and F3. The full normalized 8x8 Gram is

    G(z)=[[I5,M],[M^*,I3]],   det G=(1-|z|²)^2.

Its eigenvalues are 1 (x4), 1+|z| (x2), and 1-|z| (x2). Therefore:

* |z|<1: strictly positive rank 8. Eight independent complex CAR one-particle modes are *compatible* with the source symmetry.
* |z|=1: positive semidefinite rank 6. A paired null quotient removes two complex modes.
* |z|>1: indefinite. This violates the assumed positive Hermitian source Gram.

A 5-dimensional quotient obtained by identifying *all three* logarithmic modes with the three original meromorphic kernel modes is **impossible under these D4 and positive-marginal assumptions**: only the 2D common irreducible D4 representation may mix to saturation. This is stronger than saying that 5+3=8, but STILL does not prove |z|<1 in the physical source. There is one complex parameter not fixed by symmetry alone.

An extra full SU2-invariant spherical structure would force the cross block zero by Schur's lemma because C5=Sym^4(C²) (spin 2) and F3=Sym^2(C²) (spin 1) are inequivalent irreducible SU2 modules. However the marked D=mu4 geometry breaks the full SU2 action. Using SU2 to set z=0 is an additional parent-source hypothesis, not already proven by v512.

## Stronger result: equivariant source normalization and one grading Ward gate

On the whole strictly positive range |z|<1, the unknown mixing coefficient is **not an abstract D4/fermion CAR character parameter**. Let M(z) be the 5x3 cross block, and define the positive 3x3 matrix

    K = sqrt(I3-M* M) = diag(h,1,h),
    h=sqrt(1-|z|²)>0.

The explicit 8x8 triangular frame transformation

    T = [[I5,-M K^-1],[0,K^-1]]

obeys exactly

    T* G(z) T = I8
    [T,r8]=[T,s8]=0.

Here * denotes Hermitian adjoint. It follows that **all** |z|<1 yield the same geometric unitary D4 and quarter-clock representations after a D4-equivariant positive Gram normalization. An exact rational and genuinely complex original-symmetry compatible witness and symbolic Schur identity are checked by the executable test. Thus the **source rank condition** |z|<1 suffices for the same abstract CAR representation and the exact algebraic E8 character match, even if the original source does not fix the complex value of z. Physical time, charged sources and absolute reference STILL depend on their actual weighted forms and cannot be transformed away without checking all original observables.

There is a sharper optional **one-line source Ward gate**. The primitive spectral-triple datum already contains some grading Gamma, but its action on these geometric C5 and F3 labels has **not** been identified. If the original physical source produces the form-degree action

    Gamma_form = diag(+I5,-I3)

as a genuine *isometry* of its same reflected positive Gram, then

    Gamma_form* G(z) Gamma_form = G(z)
    implies M(z)=0 and z=0.

The actual geometric odd operator D anticommutes with Gamma_form and the D4 actions commute with Gamma_form. So this is algebraically self-consistent. **It is not automatic:** the original spectral triple's Gamma may act differently on the physical CAR source and Euclidean time reflection may exchange chirality. One must demonstrate Gamma_orig ∘ field_map = field_map ∘ Gamma_form AND source Ward invariance, not assume them. If both are proven, the full 8-mode positive Gram follows without a free mixing modulus.

A weaker demand that the constructed odd D be self-adjoint in G(z) fixes only Im(z)=0, not z=0: the exact skew Hermiticity residual has diagonal entries ±(z-z*). Hence source Hermiticity cannot replace the missing grading Ward identity.

### BRST/Hodge double-use obstruction

The finite D has six nonzero paired modes and two zero modes. If one identifies its Hodge kernel projector **itself** with the entire physical admissibility projector and insists that the E8 current fields all survive within that kernel, then the resulting physical state space has rank 2, or rank 1 after side selection. Such a finite zero-mode space cannot support 248 independent E8 current operators; the two interpretations would contradict each other. The geometric D should instead be treated as an **auxiliary seam interface**, or the actual physical projector must have a different graded/source action. This is a necessary distinction, not a general no-go for BRST constructions with nonzero excitations and infinite-dimensional Hilbert space.

## Scope and mandatory original tests

The map i is a meromorphic multiplication by dz/z, not the exterior derivative. The derivative of a general f_k has double poles at D, so it does not preserve the logarithmic differential space without a new projection/connection. Do not call D the original de Rham or APS Dirac operator. Also 5 and 3 holomorphic section degrees give different geometric conformal representations; their realization as eight chiral CAR fermions all with weight 1/2 requires a real source spin/twist dictionary.

The decisive physical checks are:

1. Obtain the actual positive Q0 source Gram and its physical eight one-particle field map. Verify the original D4 transformations, the mixed source coefficient z and whether |z|<1, =1, or fails positivity.
2. Derive the physical one-sided boundary cap, native Calderon/Dirac operator and its action on the two fixed-point channels. Confirm/refute spectral flow one rather than inserting an arbitrary endpoint projection.
3. Prove chiral CAR locality, the Spin16 simple-current E8 extension, the explicit 31-reflection Weyl clock map and actual E8 Ward current operators on the same source.
4. Compute actual R113 Z0,ZU,ZY,ZUY and their original weighted two-level pushforward and absolute reference; obtain real time and extend to a local, interacting 3+1D chiral Standard Model/Spin-2 parent.

Exact finite code: [p1_seam_dirac.py](p1_seam_dirac.py) and [p1_geometric_gram.py](p1_geometric_gram.py), with regression tests. **F1 remains open; no original R113 values assigned.**
