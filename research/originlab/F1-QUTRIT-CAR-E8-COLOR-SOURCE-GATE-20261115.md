# TFPT F1: source-type resolution of the original qutrit, CAR8, E8 colour and infinite time tower

15 November 2026. **Exact representation and root-lattice theorem with necessary physical source gates, NOT a solved F1, new standalone TOE, or original R113 values.** This integrates the current main-branch [single M3 source to Clifford comparison](f1_single_source_clifford_bridge.py), [Ramond E8 conformal grading](f1_ramond_e8_grade_bridge.py), [31-step geometric Weyl witness](p1_e8_geometric_lift.py), [native spatial infinite Fourier return](F1-NATIVE-SPATIAL-RETURN-20261115.md), and the exact native [R113 source owner](source.py).

## 1. One 3-state operator algebra is not one irreducible 8-CAR Fock algebra

The geometric "complete original operator algebra" candidate is the full complex M3(C) with Hilbert dimension 3 and traceless adjoint vector space dimension 8. If this EIGHT-DIMENSIONAL VECTOR SPACE is USED AS the one-particle space of eight new complex CAR modes, its finite fermionic Fock space is exterior(C8), dimension 2^8=256, with algebra

    CAR(C8)=M256(C).

This is a well-defined mathematical **Fock FUNCTOR**; it does not embed the product law of M3 inside the CAR algebra. The source identification must keep the two types distinct.

**Exact unital finite-matrix obstruction:** Every finite-dimensional nonzero unital *-representation of M3(C) is a direct sum of its defining 3-dimensional representation, so a unital representation on C256 would require 3 dividing 256, which is false. Even the two parity blocks M128(C) of the CAR even algebra cannot carry a unital M3(C) representation individually, since 3 does not divide 128. A rank-3 code compression P inside C256 cannot preserve all CAR products: choose A=|outside><inside|; then C(A* A)=|inside><inside| while C(A*)C(A)=0.

If BOTH full algebras have to remain simultaneously present in ONE finite Hilbert space with shared identity, its dimension must be divisible by 3 and 256, so the minimum is

    lcm(3,256)=768.

There is an exact minimal compatible realization H=C3 tensor exterior(C8), on which M3 tensor I and I tensor M256 commute and jointly generate M768. The source code independently verifies ALL 49,152 CAR anticommutator actions on every 256-mode occupation mask and ALL 55,296 combinations of qutrit matrix units and CAR annihilator on its 768 tensor basis states. No physical evidence in TFPT presently demands these two full *associative* algebras act concurrently. 768 is a conditional minimum carrier, **not** a new prediction of the original number of events or particles.

This also refines the earlier claim "M3 -> 8 traceless operators -> 16 Majoranas": the vector-space/Fock functor is valid, but one cannot preserve the original full associative M3 multiplication on the original irreducible Fock Hilbert space without an additional factor or an explicitly nonunital code with lost products.

## 2. Which SU3 representation gives the native affine LEVEL ONE?

For three complex chiral free fermions transforming as ONE fundamental SU3 triplet (plus five colour singlet fermions), their bilinear color currents have SU3 affine level k=1. This is the standard fundamental Dynkin index normalization.

By contrast, if the eight one-particle CAR modes inherit the eight-dimensional **adjoint** SU3 action of the traceless M3, their fermionic bilinears have SU3 level k=6:

    Tr_ad(T_a T_b)=6 Tr_3(T_a T_b).

The previous M3 single-source Clifford construction produces an eight-dimensional vector space with **adjoint SU3 transformation** unless a different source representation assignment is supplied. The native pinned regular SU3 root subalgebra inside E8_1 is level one. Therefore that naive adjoint-fermion current cannot be identified with pinned physical colour at level one.

A physically compatible **candidate field assignment** is instead

    V8 = F3 (one fundamental 3 of colour) plus C5 (five colour singlets),

provided the original source actually supplies eight independent positive chiral weight-1/2 CAR modes in those representation types. The five meromorphic modes are not thereby declared the five Standard Model fields or five families. This assignment yields SU3_1 rather than SU3_6 but still does not solve the quarter-clock or the actual source.

## 3. Geometric quarter clock and native colour clock: the mismatch is exact

In the original 5+3 basis the geometric quarter weights are (0,1,2,3,0) on C5 and (1,2,3) on F3. In the grouped eight-dimensional lattice chart used by the existing **explicit** E8 Weyl word, the torus is

    q_geom=(0,0,1,1,2,2,3,3).

The actual THREE F3 axes in this grouped chart are indices 3,5,7 (zero-based), not 5,6,7. Their quarter eigenphases are i,-1,-i.

Under the naïve SU3 fundamental triplet bilinear current construction, the six A2 roots and its two Cartans have C4 grading

    (2,2,2,2).

The NATIVE pinned colour SU3 roots (on lattice axes 5,6,7) instead have quarter grading

    (8,0,0,0).

Multiplying the entire F3 triplet by the original f2 orientation sign character changes its three **fundamental** phases, but leaves their pairwise differences and hence their SU3 ROOT-current C4 charges invariant. Thus a uniform orientation twist cannot simply transform the natural F3 colour currents into the native pinned **quarter-neutral** ones.

**Stronger group-level no-go, independent of the Weyl witness:** On the REAL 16-Majorana vector carrier the geometric and native C4 actions have different eigenphase multiplicities,

    geometric: (4,4,4,4),
    native:    (10,0,6,0).

No SO16 free-fermion rotation can conjugate these. On all 248 E8 currents, however, BOTH C4 actions have the same (60,64,60,64) multiplicities and are E8-Weyl conjugate, as independently proved by the 31-reflection witness. The full nontrivial Spin16 spinor extension is therefore genuinely required for that CLOCK conjugacy, not merely a convenient change of fermion basis.

## 4. Exact backtransport of the NATIVE SU3 color and E6 commutant

Pin the two simple native color A2 roots, in doubled D8+ lattice coordinates, at

    beta1=(0,0,0,0,0,2,-2,0),
    beta2=(0,0,0,0,0,0,2,-2).

Apply the INVERSE of the 31-reflection Weyl witness, keeping the actual SOURCE E8 lattice and vector-phase conventions. The complete positive A2 triplet goes to

    w^-1(beta1)       = (0,0,0,0,-2,-2,0,0)       [D8 pair root],
    w^-1(beta2)       = (1,1,1,1,1,1,-1,-1)       [spinor root],
    w^-1(beta1+beta2) = (1,1,1,1,-1,-1,-1,-1)     [spinor root].

The six source roots therefore consist of 2 D8 bilinear and 4 spinor roots under this explicit Weyl convention. The root preimages all remain GEOMETRIC-C4 neutral, as required by the native current clock. This exact mixed realization shows why the native physical color cannot be identified with just the original three free-fermion bilinears THROUGH THIS PARTICULAR witnessed E8 clock map.

Enumerating all 240 E8 roots and their pairings with beta1,beta2 gives:
- 6 SU3 root currents and 2 Cartans -> SU3 adjoint dimension 8
- 72 orthogonal E6 roots and 6 Cartans -> E6 adjoint dimension 78
- 162 other roots, split into 6 SU3 weights with **27 roots per weight**

so precisely

    E8 adjoint 248=(8,1)+(1,78)+(3,27)+(3bar,27bar).

This is the standard E8 -> SU3xE6 branching, here replayed on the EXACT native root dictionary and transported to the original geometric clock by a full explicit Weyl witness. The three weights of a color triplet are COLOR states, not yet THREE fermion generations.

Even the conformal weights are fixed by the root geometry, not numerically guessed. The inverse A2 Cartan metric is (1/3)[[2,1],[1,2]]. For EVERY ONE of the 162 mixed roots, orthogonal projection yields

    ||alpha_SU3||²=2/3    => h_SU3=1/3,
    ||alpha_E6||²=4/3     => h_E6=2/3,
    h_total=1.

The central charges at level one similarly satisfy

    c(SU3_1)=8/(1+3)=2,
    c(E6_1)=78/(1+12)=6,
    c(E8_1)=248/(1+30)=8.

Three free complex color fermions have central charge 3: after SU3_1 removes 2 the remaining c=1 plus five complex singlets c=5 agrees numerically with the E6_1 residual c=6. The central charge alone does not identify an E6 chiral coset or physical coupling; the 72-root centralizer plus the standard E8 VOA simple-current OPE would also be required.

## 5. A finite qutrit/768 host can NEVER be the full affine E8 level-one CFT

Any nonzero-level affine current Cartan relation has the form

    [H_1,H_-1] = k I

on a correctly normalized chiral mode. In finite dimension, the trace of the left side is zero, whereas the right is k*dimH nonzero for k=1. Therefore no finite-dimensional Hilbert space, 256 OR 768, can carry the full positive-level affine E8 current algebra.

The known Ramond spinor correction of the parallel [grade source](f1_ramond_e8_grade_bridge.py) supplies the required CONFORMAL weight h=1 for all 128 spinor root fields, whereas the simple finite Markov second quantization Gamma(q I8) assigns the even-occupation energies 0,2,4,6,8 and does not supply Virasoro L0. The original physical source must generate an INFINITE CHIRAL MODE TOWER, Ramond/GSO spin fields, correct current locality/OPE and a compatible positive OS source. Finite CAR zero-mode counts and a 768-state Hilbert factor are insufficient.

A new original input from the independent [native full96 spatial Fourier source](F1-NATIVE-SPATIAL-RETURN-20261115.md) DOES supply an unbounded spatial mode tower whose compressed return Hankel rank grows like 192(k+1); it cannot be replaced by one finite transfer. But its coordinate is expressly spatial, not original OS time, and its 96x96 couplings have not been tied to the E8_1 current central Schwinger term or the Ramond spinor insertion. It is a concrete candidate location for the missing INFINITE structure, not the completed affine source.

## 6. A second independent obstruction: an unscaled microscopic Markov clock has no continuum time

The existing CP source Phi_q with fixed 0<q<1 and one step of microscopic duration a gives the n-step damping q^n. At a fixed macroscopic time t>0 with n≈t/a and a->0,

    q^(t/a) -> 0,

so the resulting time family would be complete depolarization at EVERY positive t, but identity at t=0. This is discontinuous at zero, so cannot be a strongly continuous, nontrivial continuum quantum time evolution. The primitive microscopic amplitude MUST scale toward one, e.g.

    q(a)=exp(-rate*a) -> 1,
    q(a)^(t/a)=exp(-rate*t).

The rate is still a free source quantity. In a critical gapless E8 CFT, the physical mode gap must vanish appropriately; simply choosing rate positive gives a finite gapped depolarizing process, not a chiral or Lorentzian 3+1D theory.

## 7. The concrete new F1 acceptance program, not more synthetic amplitudes

The original R113 native domain still declares whole_source_bundle=null and its Fourier circle coordinate spatial_tangential_circle_NOT_OS_time. A genuine origin must therefore execute **on the native original full-field Q0 source and its bound empty reference**:

1. Identify an original positive-reflected **physical representation map** of the C5 plus F3 geometric modes into a FULL chiral one-particle CAR carrier, with their physical SU3 colour representation, grade, and no double counting; distinguish the qutrit M3 algebra (if genuinely retained) from the one-particle vector space and its Fock functor.
2. Construct the ACTUAL E8 affine spinor/NS and Ramond sectors with weight-one OPEs and k=1 Schwinger term from the native INFINITE Fourier/spatial field source, with the 31-Weyl geometric-to-native clock and native colour/E6 root dictionary. Verify all 162 mixed h=1/3+2/3 current channels on original field words.
3. Establish full interacting Osterwalder-Schrader reflection positivity and an emergent strongly continuous time with source-derived q(a), **not** the source's tangential spatial circle used as time. Produce the original Q0/Z0/ZU/ZY/ZUY and exact positive two-level fibre pushforward with the same selector and absolute reference.
4. Derive a local 3+1-dimensional chiral Standard Model and universal spin-two gravity, including the correct infinite-mode/lattice scaling and gauge anomaly constraints, rather than relying on a finite qutrit clock.

The two source-compatible routes (qutrit operator algebra and 8-CAR chiral mode algebra) should not be conflated or discarded; the theorem shows exactly which algebraic functor and physical representations would be needed for a genuine shared source.

**All statements are exact FINITE representation or E8 lattice results unless explicitly labelled SOURCE-PHYSICAL UNPROVEN. F1 IS OPEN.** No original R113 numbers have been manufactured and no Atlas/prediction ledger promotion follows.

### Exact independent reproduction

- [qutrit/CAR8 E8 source audit](f1_qutrit_car_e8_bridge.py)
- [targeted regression](test_f1_qutrit_car_e8_bridge.py)
- [31 Weyl geometric-to-native proof](p1_e8_geometric_lift.py)
- [real native spatial tail, not OS time](F1-NATIVE-SPATIAL-RETURN-20261115.md)
- [original source consumer](source.py)
