# F1 root-source research: exact CP/OS/sewing/renormalization from a pre-spatial finite algebra

**15 November 2026. New constructive mathematical source with three exact theorem gates and one original TT counterexample. This is NOT a full TFPT F1 / 3+1D TOE proof.** The objective is to derive a single internally consistent quantum source *before* treating P1, E8, spacetime or the original R113 values as outputs. These tests must never be submitted as a source-reviewed original TFPT packet.

## 0. Actual original source and why the research must start upstream

The owning original TFPT 4.5 source is [archive boundary kernel](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex), which **assumes** a smooth compact manifold with boundary and a first-order Dirac collar in its primitive datum. The [archive metrology theorem](../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex) has an impossible uniformly positive lower kernel bound on unbounded Lebesgue target; see [exact RG/source erratum](F1-ORIGINAL-MEASURE-THEOREM-REPAIR-20261115.md).

The actual [native source domain](../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/domain/native-diagnostic-v1.json) still explicitly contains whole_source_bundle=null and coordinate_role=spatial_tangential_circle_NOT_OS_time. Its owner [source.py](source.py) requires one physical whole-source finite blocked field packet, true original weighted forms, a single physical projector and an absolute empty reference. Neither the present finite algebra nor its synthetic four control marks can replace that packet. Its normalized Z0, ZU, ZY, ZUY remain physical NULL/uncomputed.

## 1. A COMPLETE mathematical finite source from ONE pre-spatial algebra

Take the basic quantum observable algebra A=M3(C), with trace state omega(A)=Tr(A)/3. It is a 9-complex-dimensional full quantum operator algebra, not a position space. With a single real amplitude q in (0,1), define the Heisenberg quantum channel

    Phi_q(A) = q*A + (1-q)*Tr(A)*I3/3.

For every q in [0,1] this has an explicit Kraus representation

    K0=sqrt(q)*I,
    K_ij=sqrt((1-q)/3)*|i><j| for i,j=0,1,2.

The ten Kraus operators satisfy sum K* K=I. Thus Phi_q is **completely positive**, unital and trace preserving in the dual. It is self-adjoint in the Hilbert-Schmidt trace inner product, with spectrum 1 on the scalar identity and q on the eight traceless directions.

The time-composition and exact state algebra obey

    Phi_q Phi_r = Phi_(qr),
    P_q P_r = P_(qr),  P_q=q*I3+(1-q)*J3, J3=11^T/3.

The finite *classical event sector* obtained by restricting A to the three diagonal projectors is a reversible Markov process with stationary state pi=(1/3,1/3,1/3). Assign each history s0,...,sn its full finite probability

    Pr(s0,...,sn)=pi(s0)*prod P_(q_k)(s_k,s_(k+1)).

Sum over an intermediate event exactly gives the coarser P_(q1*q2). The code checks all nine endpoint pairs at q1=q2=1/2 and proves the identity symbolically. This is a real two-level weighted source, not an approximate equality of Dirac spectral weights.

### Exact Osterwalder–Schrader reflection positivity

For arbitrary finite collections of positive-time quantum observables A_i and corresponding semigroup parameters q_i, the OS one-word Gram is

    G_ij=Tr[(Phi_(q_i)A_i)^* (Phi_(q_j)A_j)]/3

because Phi_(q_iq_j)=Phi_(q_i)Phi_(q_j). It is a Hilbert-space Gram, hence PSD for all choices. A 3x3 mixed indicator test has full rank three and exact positive determinant 49/15552.

For the full classical history measure, take reflection through a central state x0. Conditioned on x0, the positive and negative half histories are independent and have identical transition laws, so for ANY positive-time cylinder function F,

    E[(theta F)^* F]
      =sum_(x0) pi(x0)*|E[F | x0]|² >= 0.

This is a genuine full classical-history OS proof for the constructed finite process, not just positive pointwise weights. The exact first-time reflected cylinder Gram is (1/3)P_q² with eigenvalues 1/3, q²/3, q²/3. For q=1/2 the exact eigenvalues are 1/3 and 1/12 twice.

The finite model is BOTH exactly OS positive and projectively consistent, fixing the gaps found in the earlier arbitrary spectral-action source. These properties hold ONLY for this explicitly given source, not for the original TFPT Q0.

## 2. A dimensionless clock follows the composition; its physical scale does NOT

For any nontrivial eigenvalue q in (0,1), the multiplicative sewing character defines an additive operational clock

    tau(q)=-log(q),  tau(qr)=tau(q)+tau(r).

The clock coordinate is **derived after selecting the microscopic transfer amplitude**, not inserted as a one-particle spacetime coordinate. However the value of q for one physical step and the conversion to a dimensionful elapsed time remain free and are not fixed by positivity.

A critical RG selection no-go follows immediately: demanding literal identical transfer probabilities before and after blocking two steps without renormalizing the clock means P_q=P_q², i.e. q=q², with only q=0 (instantaneous mixing) or q=1 (no evolution). Nontrivial dynamics requires a scale-dependent q (or time rescaling), not an exact unchanged bare one-step kernel. This sharpens the original "time must be born through gluing" hypothesis: a physically useful source must specify a time/scale functor, not merely a fixed matrix and a blocking map.

### Two different source kernels already produce different connected responses

The same arithmetic four mark *control* with U and Y both indicators for the same state at two adjacent steps yields

    Z0=1,  ZU=ZY=1/3,
    ZUY=(1+2q)/9,
    C_UY=2q/9.

At q=1/2: ZUY=2/9 and C_UY=1/9.
At q=1/4: ZUY=1/6 and C_UY=1/18.

Both are strict positive, stationary, permutation-symmetric, exactly OS and sewn sources with the same first three control values, yet their connected U/Y responses differ. All are **synthetic** control values, NOT the four actual R113 Z0/ZU/ZY/ZUY. This exposes the minimal observational information required to identify a dynamical rate even in this tiny declared model.

## 3. A stronger UNIQUENESS theorem for the whole quantum internal source

Add the explicit new premise: the original quantum source is basis-free under **all inner SU3 automorphisms** A->UAU*. For a general linear map M3->M3 there are 81 free complex coefficients. Construct eight exact traceless sl3 generator commutators on the whole 9D matrix space and solve all 648 linear commutant equations. The coefficient matrix has exact rank 79; the commutant has dimension TWO, spanned by projection onto scalars E(A)=Tr(A)I/3 and its traceless complement I-E. This is the explicit finite Schur lemma for 1+adjoint8.

A unital, trace-preserving SU3-covariant *-preserving channel is therefore exactly Phi_q. Its full 9x9 Choi matrix has eigenvalues

    (1+8q)/3  (multiplicity 1),
    (1-q)/3   (multiplicity 8),

hence complete positivity iff -1/8<=q<=1. A continuous positive-time CP semigroup starting at the identity requires q(t+s)=q(t)q(s), q(0)=1 and q(t)>=0 continuously, hence

    q(t)=exp(-kappa*t),  kappa>=0.

After measuring elapsed time in units of the unknown kappa, the dimensionless generator is **unique**:

    L(A)=Tr(A)I/3-A.

Thus **a pre-spatial, finite quantum dynamics is unique modulo overall clock scale in this exact declared covariance/source class**. It is a real mathematical theorem, not an arbitrary fit. But full SU3 basis-free covariance, the physical algebra M3 and its trace state are not derived from the original TFPT source. The resulting finite generator also cannot spontaneously generate the required particle spectrum, symmetry breaking, E8 VOA, local 3+1D kinematics or a universal spin-two graviton.

## 4. D4 symmetric full hidden history explains a nonautonomous *observed* clock

To avoid conflating an observed reduced transfer with the full source, take a second finite EXACT model on four mark-labelled states, with adjacent transition probability a=1/10 and opposite transition probability b=1/20:

    W = [
      [3/4,1/10,1/20,1/10],
      [1/10,3/4,1/10,1/20],
      [1/20,1/10,3/4,1/10],
      [1/10,1/20,1/10,3/4]
    ].

This W commutes with the geometric quarter rotation and inversion D4 permutations, is reversible, stochastic and positive as a transfer operator. Its exact eigenvalues are (1,3/5,7/10,7/10). All 16 weighted fine-to-coarse endpoint sums of two transitions reproduce W², and the global reflected-history Gram has strictly positive eigenvalues 1/4,9/100,49/400,49/400. The full source has the same general classical conditional-reflection OS proof as above.

Now observe ONLY whether the system is at one particular mark or at one of the other three. The reduced one-step observed transfer is

    R1=[[3/4,1/4],[1/12,11/12]].

The actual two-step observation is not R1²:

    R2-R1²=[[1/600,-1/600],[-1/1800,1/1800]] != 0.

The reason is exact non-lumpability: from the hidden "other three" block, two states return to the marked state with rate a and the opposite state with rate b. The full process has legitimate dynamics; the observed projection has memory and is not governed by the naive one-step Markov transfer. This is a **constructive possible mechanism** behind the previous TFPT reduced-channel no-homogeneous-Lindblad diagnostic, NOT a dilation or physical evaluation of the actual TFPT U=I-2DD*, q2=(25/1568)(1-g)² or original 20-to-150-to-400 channel data. The original hidden degrees and source field map still have to be identified.

## 5. Independent original gravity failure: positive Stieltjes does NOT guarantee only one TT pole

The original [TFPT4.5 metrology source](../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex) asserts from positive Stieltjes form factors F(z)=a_F+int dmu_F(t)/(t+z) that the only TT propagator pole of [zF(z)]^-1 is at z=0 and thus no additional physical spin-two poles exist. The inference from zero-free on C minus the negative real axis to no extra poles ON the physical negative real axis is invalid.

Already the exact allowed positive Stieltjes form factor

    F(z)=1+1/(z+1)

gives TT propagator

    [zF(z)]^-1 = (z+1)/[z(z+2)]
               = 1/(2z)+1/[2(z+2)].

There is an ADDITIONAL massive spin-two pole at z=-2. Its residue is positive, so this example does not prove a ghost, but it refutes the claimed massless-only TT spectrum. In general F(z)=a+mu/(z+m), a,m,mu>0, has a negative-axis zero at z=-(m+mu/a), yielding another pole with positive residue mu/[a(a*m+mu)]. The original physical source must control its zero/branch structure to prove the actual graviton spectrum. Stieltjes positivity alone does not.

## 6. A single internal quantum algebra can supply the CAR16 / spinor128 *algebraic* bridge, but not the affine E8 conformal grade

The previous P1 5+3, 16 Majorana and E8 statements were separate candidate identifications. Starting NOW from the exact same **positive M3 Hilbert–Schmidt source**, the traceless subspace

    V8 = {A in M3(C):Tr A=0}

has exactly eight complex dimensions, an explicit positive orthogonal Hermitian basis (metric diagonal 2 repeated seven times, 6 once), and inherits the principal quarter rotation by conjugation with the real SO3 90-degree matrix. The actual 8D matrix has FOUR eigenspaces of multiplicity (2,2,2,2), exactly the previously audited P1 geometric C4 character. The original full geometric D4 intertwiners remain a separate finite twist and must not be conflated with this abstract eigenspace census.

Apply the **canonical fermionic exterior/Clifford functor** to the Hilbert space V8. It has eight complex CAR operators and therefore 16 real Majorana generators. The Fock space has dimension 2^8=256 and its even halfspinor contains 128 states. We tested **all 49,152 source CAR anticommutators** on all 256 occupation masks. The 16-dimensional Majorana vector has phase classes (4,4,4,4), its 120 bilinear so16 generators have classes (28,32,28,32), and the even 128 spinor states (32,32,32,32); the combined 248 vector-space classes (60,64,60,64) agree with the known standard E8=so16+spinor finite Lie algebra construction and native current character. This is the first direct route *from the same positive M3 source* to the 16 Majorana algebra, rather than separately postulating eight independent physical fermion fields from the 5+3 geometry.

The fermionic functor itself is a **CHOICE OF QUANTIZATION**: it is mathematically canonical once CAR is selected but is not forced by the original TFPT source. Nor does counting the E8 character prove that the 248 root current Lie brackets or level-one vertex OPEs occur on the physical reflected Hilbert space.

### A sharp additional obstruction: the lifted microscopic clock does not give E8 weight one

Second quantizing the ORIGINAL chosen M3 depolarizing transfer on the eight traceless modes yields the positive exterior transfer Gamma(q I8). Its action on a k-particle state is **q^k** and its common vacuum matrix element is exactly **1**. The ordinary trace and supertrace differ sharply:

    Tr_Fock Gamma(q I8) = (1+q)^8,
    Str_Fock Gamma(q I8) = (1-q)^8.

For q=1/2 these are 6561/256 and 1/256, versus the native R113 EMPTY Fock expectation 1. No full trace is licensed by the original zero-fermion cap, consistent with the prior R113 source contract.

The unique normalized source clock on Fock is the occupation operator N, so the 128 even Fock states lie in **different time-energy grades** 0,2,4,6,8. In contrast, the 128 spinor ROOT CURRENTS of the actual affine E8 at level one all have conformal weight ONE. No nonzero affine rescaling H=kappa N+c can assign equal weight one to all even Fock states; the only algebraic solution is kappa=0, c=1, which makes the source time evolution TRIVIAL on that sector. Thus naive second quantization of the uniquely covariant M3 channel does **NOT** produce the local affine E8 current dynamics, even though its CAR dimensions and quarter character match exactly.

The missing primitive is now more precise: a source-derived **chiral conformal energy grading / local operator product and spinor extension** on a continuum one-particle mode space, rather than another 248-dimension match. Finite CAR Fock is only an abstract zero-mode module. An honest E8_1 construction needs infinitely many chiral fermion oscillators, current algebra locality/OPE, spin structure/GSO selection and a physical one-sided reflected source that realizes those fields. This is precisely what the original one-sided TFPT Dirac/BRST source has not supplied.

Executable exact bridge: [f1_single_source_clifford_bridge.py](f1_single_source_clifford_bridge.py) and [six independent tests](test_f1_single_source_clifford_bridge.py). **F1 remains open.**

## 7. What has REALLY been closed; strict gates still separating TFPT from F1

Established in concrete finite mathematics, with reproducible exact tests:
- A pre-spatial M3 quantum source with explicit Kraus complete positivity, unitality, stationary state and a single common microscopic rule.
- Exact two-level weighted source pushforward, classical full OS and quantum one-word OS, with derived additive dimensionless transfer clock.
- Conditional uniqueness of its internal CP semigroup modulo an overall time unit under FULL inner SU3 covariance.
- A D4 full-history source whose projection is non-Markov yet remains an observable of a completely consistent OS and sewn history process.
- An independent source-paper TT spectral extra-pole counterexample.

Not established:
- Why the original universe should use M3, the full SU3 automorphism source covariance, the Markov step, or the chosen D4 hidden rate pattern; these are new toy premises, not consequences of the source-free origin.
- The actual **original** original_blocked_fiber_packet, Q0 with native physical projector, boson/fermion/ghost/effective action and absolute reference, with all original source marks Z0,ZU,ZY,ZUY.
- A source-selected original square P1 modulus, charged root E8 Ward operators, local unitary interaction and 3+1D continuum, Standard Model chirality, physical clock/scale, graviton/spin-two spectrum and cosmology.
- The finite M3 heat trace is 1+8exp(-kappa*t) and tends to 9 as t->0, so its spectral dimension tends to ZERO in the ultraviolet; a genuine 3+1D continuum requires a nontrivial family of growing source algebras and a proved continuum limit, not another finite algebra identification.

**F1 remains OPEN.** The next original theorem must construct one family of physical positive reflected source functionals across ALL scales, with an exact weighted conditional expectation/BRST null-quotient functor and a unique microscopic rule. It must then realize the source-specific E8 and gravity observables on that same source. Without this last identification, neither mathematical CP uniqueness nor a finite source example is a theory of everything.

Exact code: [f1_functorial_source.py](f1_functorial_source.py), [f1_covariant_quantum_semigroup.py](f1_covariant_quantum_semigroup.py), [f1_hidden_history_dilation.py](f1_hidden_history_dilation.py), [f1_spin2_stieltjes_gate.py](f1_spin2_stieltjes_gate.py) and their corresponding independent tests. These have been executed locally before being committed to main and hooked into OriginLab CI.
