# TFPT F1 original regional sl2 -> 5+3 principal SU3: exact new gates

15 November 2026. This is a new exact *algebraic* origin candidate and a
new E8-native exclusion. F1 is NOT solved. The physical cubic source,
CAR/fermion interpretation, true R113 all-field Q0, time, geometric
3+1D and universal gravitation remain open. Never promote the ledger.

## 1. Actual original inputs, not a new arbitrary 8D module

- [Original c3/P1/P2 architecture](../../documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex): the *same* square divisor on P1 produces 5 meromorphic modes and 3 logarithmic one-forms. Both types are exact but their physical identification remains conditional.
- [Original regional Lie span](source_bridge.py): four native Möbius/parabolic polynomial vector fields span exactly {1,x,x²} with nontrivial Lie brackets [1,x]=1, [1,x²]=2x, [x,x²]=x², hence sl2 over C, NOT three commuting independent local translations. The physical source action on the 5 and 3 modes is still unproved.
- [Exact five/three source and positive Gram](P1-SEAM-DIRAC-GRAM-20261115.md) and [native E8 Weyl clock](P1-FIVE-PLUS-THREE-WEYL-E8-20261115.md). They prove mathematical characters/clock conjugacy but not a physical field OPE.
- Standard Kostant principal sl2 theorem: a simple rank-two Lie algebra's adjoint representation decomposes into two sl2 irreducibles of dimensions 2e1+1 and 2e2+1, with e_i its exponents. The possibilities A2: (3,5), B2: (3,7), G2: (3,11). Source: Bushek and Kumar (2014), https://kumar.math.unc.edu/papers/kumar71.pdf . The compact real su3 decomposition is also the classical SO3 Cartan representation on real symmetric traceless 3x3 matrices (3+5).

## 2. Exact principal SO3 -> compact SU3 realization

Use real matrices A in so3 (the 3 real skew-symmetric 3x3 matrices) and S in Sym0(R³) (the 5 real symmetric traceless 3x3 matrices).

    su3 = so3 ⊕ i Sym0(R³), dim_R 3+5=8.

Under the SO3 adjoint action, the first part is spin one, the second spin two; their complexifications are exactly Sym²(C²) and Sym⁴(C²), the same representation types as the original logarithmic-forms and meromorphic-section spaces after a chosen source-compatible line-bundle trivialization. The compact Lie bracket has types

    [3,3] subset 3,
    [3,5] subset 5,
    [5,5] subset 3.

All 28 independent commutators close, span all eight generators, satisfy all 56 independent Jacobi generator triples, and have positive invariant metric -Tr(XY) with eigenvalues 2 (multiplicity seven) and 6 (multiplicity one) in the explicit nonorthonormal basis. This is a genuine nonabelian *mathematical* completion, unlike dimension-only readings.

**Caution:** the fixed four-point divisor has only D4, not the full SO3/PSL2 as automorphisms. One must derive how the original regional sl2 acts on the actual physical section and differential source fields across source configurations, instead of asserting an unbroken continuous compact group.

## 3. Full D4 sign/clock twist, with an internal original line

Let r:z->i z and s:z->1/z. In the original geometric basis f_k=z^k/(z^4-1), omega_j=z^(j-1)dz/(z^4-1):

    r f_k = i^k f_k,       s f_k = -f_(4-k),
    r omega_j = i^j omega_j, s omega_j=+omega_(4-j).

The original 5- and 3-dimensional D4 representations are NOT literally the principal SO3 adjoint representation. They differ by the **same** 1D D4 character

    chi(r)=-1, chi(s)=-1.

This is proved with explicit invertible equivariant intertwiners T5 and T3; their determinants are 32768 i and -512 i. Both r and s acts as 90-degree and 180-degree SO3 rotations in the twisted representation.

The necessary character is already present: the single original geometric central meromorphic section

    f2(z)=z²/(z⁴-1)

has r f2 = -f2 and s f2 = -f2. It has zeros of order two at the geometric quarter fixed points z=0 and infinity. Thus its ONE-DIMENSIONAL REPRESENTATION supplies the exact twist without introducing a new abstract character. This does NOT give a globally invertible multiplication by f2 across those zeros; a source field/cap and its singularity/line-bundle treatment must be checked.

In the original *untwisted* 8D geometric representation a standard SU3 bracket is not a D4-equivariant ordinary scalar bracket: it is chi-twisted. Passing to the typed representation V⊗<f2> makes it an ordinary D4-equivariant SU3 bracket. Whether the physical source allows this tensor-line insertion is a new source question, not derived by character arithmetic alone.

## 4. Source positivity and unique nonabelian completion up to coupling

The SO3-intertwining space

    Hom_SO3( wedge²(Sym0(R³)), so3 )

has dimension EXACTLY ONE (30 unknown coefficients, linear constraint rank 29). Therefore, after fixing the natural so3 bracket and the standard spin-two action, the only remaining alternating spin2-spin2 bracket is a scalar lambda times the standard commutator.

Define for real 3x3 skew A,B and real symmetric traceless P,Q:

    [(A,P),(B,Q)]_lambda =
       ([A,B]-lambda[P,Q], [A,Q]+[P,B]).

The family satisfies Jacobi for every lambda. The positive source-adjoint metric, if the ORIGINAL weighted OS/Ward source supplies it, has sector norms

    <A,B>=a*(-Tr AB),     <P,Q>=b*Tr(PQ),   a>0,b>0.

Already one explicit [P,Q],A Ward component (and independently all 8^3=512 basis coefficient identities) gives

    4*(a*lambda-b)=0,
    lambda=b/a>0.

Thus under the PRECISE additional premises of full principal SO3 action, a nontrivial compact ad-invariant **positive** Ward form, and this one-parameter bracket class, the Lie algebra must be the compact su3 (up to the source-defined overall coupling). lambda=0 gives a semidirect so3 ⋉ R5 without a positive ad-invariant form compatible with the nontrivial spin-two action, and lambda<0 the split sl3(R), likewise noncompact.

This is a potential **single scalar original-source decision**. A physical proof must bind the cubic [5,5,3] Ward response lambda and positive two-point sector norms a,b from the SAME actual Q0 and time/reference, and check the Ward equality. Their values are currently missing. Positive pointwise weights without the full OS Hermitian form do not suffice.

## 5. Why a FOURTH oriented insertion is mathematically required for an invariant cubic

The compact SU3 Cartan cubic

    kappa(X,Y,Z) = -Tr([X,Y] Z)

is nonzero; a chosen [5,5,3] generator triple yields kappa=4. But the *actual geometric* 8D D4 action is chi times the SU3 adjoint automorphism. Therefore its cubic transforms as chi cubed = chi, i.e.

    kappa(rX,rY,rZ)=-kappa(X,Y,Z),
    kappa(sX,sY,sZ)=-kappa(X,Y,Z).

All 56 independent generator triples and both D4 generators were checked (112 exact identities). An unmarked cubic scalar expectation in a fully D4 invariant source must vanish. Multiplication by the source-compatible chi character from f2 cancels the sign:

    f2 * kappa(X,Y,Z)  is invariant under r and s.

This is a new **four-insertion source hypothesis**: three nonabelian current fields plus ONE compensating geometric f2-mode insertion can produce a D4-even marked scalar response. f2 is a MODE of the original four punctures, not a fifth puncture. The physical original R113 U/Y insertions are NOT automatically these fields. Without a proved native field map this is a selection rule, not a source amplitude.

## 6. Count and kill 1,120 regular E8 SU3 embeddings

Enumerate every native D8+ E8 norm-two root triangle, identifying identical A2 six-root sets. There are EXACTLY 1,120 regular A2 root subsystems. Their 8-generator native glue C4 distributions are:

| C4 profile (grades 0,1,2,3) | Count |
|---|---:|
| (2,2,2,2), matching geometric V8 | 480 |
| (4,2,0,2) | 416 |
| (4,0,4,0) | 180 |
| (8,0,0,0), fully fixed | 44 |

Of 480 quarter-balanced A2 copies, 48 commute with the repository's native pinned SU3 color roots. But among ALL E8 regular A2 subsystems commuting with the **entire** pinned colour SU3, weak SU2 and physical hypercharge Cartan, only 10 survive: six with profile (4,2,0,2), four fully fixed. **None** has (2,2,2,2).

Therefore the *geometric* principal SU3 algebra CANNOT be an extra independent SU3 family gauge algebra both (i) acting on exactly the proposed 5+3 representation with balanced quarter action AND (ii) commuting with the entire native SM gauge algebra, **if one restricts to regular E8 root A2 subalgebras and preserves the exact pinned charge dictionary**. This is a useful falsification, not a no-go for nonregular affine embeddings, cosets, or dynamical sector-changing operators.

Also, the repository's pinned physical color SU3 is fully quarter-neutral and has profile (8,0,0,0); it is NOT automatically the same as the balanced geometric principal SU3. One must not identify these two eight-dimensional spaces by their dimension alone.

## 7. Two competing physical completions share the same geometry

- **Fermionic route:** treat C5⊕F3 as eight independent *chiral complex CAR fields* of conformal weight 1/2. Their sixteen Majoranas may give the standard so16+128 spinor E8 level-one completion (the exact 31-reflection Weyl compatibility is proved on all 240 roots). Physical independence, source CAR weights, chiral spin structure, E8 OPE and OS remain unproved.
- **Principal-current route:** treat the same 5+3 types as an 8-generator compact SU3 current algebra of conformal weight 1 under a source-derived principal sl2 action and f2-twisted cubic. Its positive Ward structure is exact in the declared representation. The actual local current OPE and pinned physical color embedding remain unproved.

These are NOT physically identified by D4 characters or 8-dimensional arithmetic. The original source needs an actual two-point short-distance scaling and the **three-plus-f2** ordered cubic response. A fermion field's leading OPE pole scales as (z-w)^(-1), a current two-point as (z-w)^(-2); these are different fields and may of course coexist in one theory. Do not force an exclusive choice if the original theory contains both primitive fermions and composite currents.

The crucial original F1 input remains the real weighted Q0, absolute empty reference, original field/BRST map, time, Z0/ZU/ZY/ZUY and 3+1D chirality/gravity. Do NOT promote an algebraic candidate to physical F1.

## Executable exact evidence

[p1_principal_su3.py](p1_principal_su3.py), [p1_principal_su3_ward.py](p1_principal_su3_ward.py), with unit regression tests. The 8-generator Lie closure, explicit D4 intertwiners, 1120-root centralizer census and positive Ward relations were also reproduced independently in a local standalone SymPy audit. CI is asynchronous and no green status has been assumed.

**F1 OPEN.**
