# Native marking, internal carry and the exact response memory

14 October 2026. Research ID `UR.SOURCE.NATIVE_MARK_MEMORY.01`.
Overall physical contract **PARTIAL**. The statements below are source-specific mathematical results; no physical T1–T8 closure follows.

## 1. Target and source dictionary

The requested target is a common original source selecting marking, state and event law. Work with the exact original Round20–26 neutral scalar/integer-charge source already reconstructed in `native-source-reconstruction-20261014`, not a new bath, collision or reset model. N=L^3, L>=3, a,m>0, u,nu,J>=0; the nonzero-memory assertion requires hopping J>0. Its original H, positive transfer and faithful Gibbs rho_beta are retained, including the hopping constant48JN. Beta>0 remains prescribed.

Three distinct objects must not be identified by their number of labels:

* The original grade k:L->Z4 gives four commuting charge projectors M_a. Their algebra is C4, not M4(C).
* The original D5/D3 splitting gives a rank5 projector in the eight-dimensional current index space. It is not a projector onto five quantum states.
* The geometric four punctures, an earlier Clifford M4 channel and a physically occurring measuring instrument need separate operator maps. None is identified with M_a merely because there are four labels.

The physical E8 basis is B=[2e1,e1+ei for i=2..7,(1/2)^8]. A site coefficient vector n corresponds to Bn, with Gram G=B^T B and the unchanged triangular sign cocycle c. The existing grade is

    k(q)=2 sum_(i=1)^5 q_i mod4.

Its coefficient-basis increments are (0,0,0,0,0,2,2,1). This grade is inherited from the earlier carrier construction. Its microscopic physical selection is not silently reclassified as a consequence of H.

## 2. Recover the marked subspaces from the original grade

The known carrier-module theorem supplies an intrinsic reconstruction once norm and k are specified. The norm2 roots with k=0 are the40 D5 roots and12 D3 roots. Their second moment is

    M0=sum_(alpha root,k(alpha)=0) alpha alpha^T =16P5+8P3,
    P5=(M0-8I)/8, P3=I-P5.

This eliminates a separate arbitrary coordinate projector; it does not eliminate the premise k. The calculation follows from the D_m roots +/-e_i+/-e_j: off-diagonal moments cancel and the diagonal coefficient is4(m-1). The checker repeats it on all240 original roots for each of the15 source-symmetry-related grades described below. It is a prior theorem, not a newly discovered marking rule.

## 3. The carry is an internal source coordinate

Round16 already proves the exact arithmetic decomposition, including every original neutral direction. Its integer lift of k is

    j(q)=4 sum_(i=6)^8 q_i-2 sum_(i=1)^5 q_i,
    j(s)=1, s=(1/2)^8, j mod4=k.

With K=ker j, every q has a unique representation

    q=K u+(4m+r)s, u in Z7, m in Z, r=0,1,2,3.

Thus ell2(L) is explicitly unitarily equivalent to ell2(Z7) tensor ell2(Z) tensor C4, with the original cocycle phases retained. The complement coordinates are a recording convention, not a physically preferred subsystem. Integer carry m is not an appended spectator.

The same construction applies to the neutral charge configuration lattice C0 at N sites. Fix adjacent departure x and arrival y. The primitive neutral vector v=(-s at x,+s at y) has j_y(v)=1, so

    C0=ker j_y direct_sum Z v.

Every configuration is z+(4m+r)v with z in ker j_y. This directly gives a four-grade factor and the entire original remaining charge memory, plus the unchanged scalar Hilbert space. A grade-preserving complement choice only changes this coordinate representation.

On every bilateral orbit of the native W=W_xy,s, recursively determined unit phases remove its one-dimensional hopping sign. On the orbit starting at zero the original W phase is -1 at every j; d_j=(-1)^j makes W raise j with coefficient+1. This rephasing does not delete the cocycle of other directions.

## 4. Derive an internal M4 and keep its winding

Let P_r be the four grade projectors on that orbit and Z=W^4. Define

    R=W(I-P3)+W^(-3)P3.

In the original (m,r) coordinates, R increments r modulo4 while holding m fixed. Therefore R is unitary, R^4=I, [R,Z]=0, and

    W=R[(I-P3)+ZP3].

The exact native matrix units are

    E_ab=R^a P0 R^(-b), a,b=0..3,
    E_ab E_cd=delta_(b,c) E_ad, E_ab*=E_ba.

Consequently a constant M4(C) factor is generated inside the same source-orbit algebra. It is not a quotient that discards winding. The full algebra is

    C*(W,P_r) = M4(C(S1)),

and its von Neumann closure in the bilateral representation is M4(L_infinity(S1)). Indeed Fourier transformation of the original carry makes Z multiplication by z on the unit circle; W is the four-step cyclic matrix with wrap entry z. The E_ab generate all constant matrix units, and Z has the full unit-circle spectrum. This proves both inclusions and faithfulness.

No selected Bloch phase or finite periodic carry is introduced. Replacing Z by a fixed phase is a representation of this smaller algebra, not an established reduction of the original time evolution. The R construction uses a cut after grade3; shifting that cut changes the coordinate frame. Nor is this abstract M4 already identified with the old Clifford CPTP channel.

## 5. The same source couples grade and carry

For the original adjacent v, D(v)=2/N+14nu and, at n=jv,

    V(jv,phi)=j^2[D(v)+u(phi_x^2+phi_y^2)]+V_scalar(phi).

The14 consists of4 from the common bond and10 from the other incident bonds. Thus

    V((j+4)v,phi)-V(jv,phi)
      =(8j+16)[D(v)+u(phi_x^2+phi_y^2)].

Every original hopping commutes with W^4: its cocycle commutator with a4s shift has even pairing. K also commutes with W^4. On the common core,

    [H,W^4]=W^4[V(n+4v)-V(n)],

which is nonzero, already at n=0. The center of the grade/transport subalgebra is not central in the full dynamical source. A fixed-z fiber therefore does not carry the full H as an invariant direct-integral fiber. Fourfold wrap information affects subsequent source responses.

More strongly, use the same normalized scalar Schwartz vector and z=0 in the full factorization from section3. The product-basis diagonal energies for j=0,1,4,5 obey

    E5+E0-E4-E1=8[D(v)+u<phi_x^2+phi_y^2>]>0.

The scalar baseline and48JN cancel; every nonconstant hopping word has zero diagonal charge matrix element. Any tensor-sum H_grade tensor I + I tensor H_memory would have zero mixed difference. Hence this derived grade/memory factorization is dynamically coupled.

The faithful full Gibbs state cannot be rho_grade tensor rho_memory in this factorization. If it were, spectral logarithms of the positive injective densities would give

    H=-beta^(-1)(log rho_grade tensor I + I tensor log rho_memory)
      -beta^(-1)log Z_beta I,

a tensor-sum, contradicting the diagonal witness. This proves that the actual source equilibrium already correlates marks and their memory. It does not forbid every other mathematical factorization or derive a particular measurement apparatus.

## 6. Exact memory for the four observed grades

Let K_rho be the completion of bounded operators in <X,Y>_rho=Tr(rho X*Y). It is isometric to Hilbert-Schmidt space through X->X rho^(1/2). Source Heisenberg evolution alpha_t(X)=e^(itH)Xe^(-itH) extends to a strongly continuous unitary group e^(itL), with self-adjoint Liouvillian L and the commutator action L X=[H,X] on its domain.

At site y put S=span{M0,M1,M2,M3}. All four are nonzero even in the total-neutral sector. Faithfulness gives weights w_a=Tr(rho M_a)>0, so S has dimension4. Its orthogonal GNS projection is

    P X=sum_a M_a Tr(rho M_a X)/w_a, Q=I-P.

The unbounded diagonal and scalar terms commute with M_a. The only commutator is [Hhop,M_a], which is bounded. Thus S is in Dom L; the finite-rank maps A=PLP and B=QLP are bounded. The opposite off-diagonal block is B*. Removing B+B* from L is a bounded self-adjoint perturbation, leaving a self-adjoint block diagonal operator. In particular L_Q=QLQ on the Q subspace is self-adjoint with its induced domain.

Duhamel's formula for these bounded off-diagonal couplings gives an exact closed Volterra equation for the uninserted compressed response C(t)=P e^(itL)P:

    C'(t)=iA C(t)-integral_0^t K(t-s)C(s)ds,
    K(t)=B* e^(it L_Q)B, C(0)=I_S.

All ingredients are determined by the same source H, rho and original grade; there is no independently chosen bath state or forgetting rate. This is the source-specific projection-memory construction, not new general Mori mathematics.

## 7. All three nonconstant grade modes need memory

The exact native e8=s neutral hop sends all-zero charges to (-s,+s). Its original W sign is-1; the -J W hopping term therefore has matrix coefficient+J. At the arrival site it changes grade0 to1, and

    <(-s,+s)|[H,M0(y)]|(0,0)>=+J.

For any diagonal f(k), [H,f(k)] is off-diagonal in charges. If B f=0 then Lf lies in S and is diagonal, so every off-diagonal entry vanishes. The same e8 hop connects all four grades, forcing f(a+1)=f(a). Thus, for J>0,

    ker B=span{I}, rank B=3,
    K(0)=B*B is positive semidefinite, nonzero, rank3.

It is NOT strictly positive on all four dimensions: its constant direction is zero. The finite compressed matrix has the second-derivative identity

    C''(0)=-A^2-B*B.

Only S subset Dom L is required: matrix-element differentiation uses <Lu,Lv> rather than assuming every vector is in Dom L^2. An autonomous group e^(itA) would instead have second derivative -A^2. Therefore the original four-grade compression is not that memoryless group. This statement does not rule out every time-dependent or nonunitary phenomenological fit.

The minimal unitary realization of this two-time memory kernel is the closed subspace generated by e^(itL_Q) B S for all real t. It is internal to the source GNS space. No claim that its dimension is finite, or that this two-time subspace suffices for arbitrary interventions, is made.

## 8. Preserve complete inserted histories

For a complete grade-history response, use the larger source-selected GNS closure generated from I by e^(itL), all its inverses, and left multiplication by the M_a. Left multiplication is bounded and has adjoint left multiplication by M_a itself. The resulting closed subspace is reducing for these operations. It contains every prescribed ordered history vector

    C_h=M_(a_n)(t_n)...M_(a_1)(t_1).

Its Gram matrix reproduces the full source decoherence functional

    D(h,h')=Tr(C_h rho C_h'*)=<C_h',C_h>_rho.

This is the minimal closed representation retaining these history vectors and actions; the usual Gram/GNS uniqueness is up to an intertwining unitary. It retains interference and any information required by the stated probe family. The unitary kernel in section6 describes only its uninserted two-time compression. Classical event probabilities additionally require a specified measuring instrument or decoherence condition. Their physical occurrence is not selected by forming a Gram matrix.

## 9. A grade measurement cannot preserve the equilibrium by simple forgetting

There is no normal rho-preserving conditional expectation from the full bounded source algebra onto the four-grade abelian algebra when J>0. Bimodularity would imply E(M_a X M_b)=0 for a!=b. State preservation and trace duality would force M_b rho M_a=0, so [rho,M_a]=0. Because exp(-beta H)/Z has an injective scalar spectral function, this would force [H,M_a]=0, contradicting section7.

Conversely, a grade-block-diagonal rho would give a CP expectation by the coefficient formula in section6. The present source does not satisfy that condition. Thus GNS orthogonal projection is an exact response calculation, not a physical CP erasure. Nonselective sharp grade measurement changes rho.

## 10. What the unmarked H does and does not select

Permuting physical coordinates2..7 permutes the six equal-weight root hopping directions, fixes the long and half-root directions, and preserves G. Its triangular cocycle changes by the exact inversion coboundary. Generator rephases cancel on neutral bonds, so the same H, Q_delta and rho are invariant under this coordinate S6.

For each four-element S subset{2,...,7}, set k_S(q)=2(q1+sum_(i in S)q_i) mod4. These15 distinct grades have the same source symmetry orbit and each reconstructs a5+3 projector by section2. A unique intrinsically equivariant rule applied to this unmarked H alone cannot choose a named member of that orbit.

This is NOT15 different physical worlds and NOT a no-go for the whole TFPT source. Choices related by simultaneously transporting all data are equivalent. Fixed original k or actual clock/field data may supply a relative alignment. The ambient coordinate S6 is also not the separate Sp(4,2) compiler S6. The calculation prevents silently treating the existing k as if H had already selected it.

## 11. The actual joint clocks require coherent transported marking

Use unchanged orthogonal C,J from `compiler-integral-triality-20260918/certificate.json`, in the same eight-component current basis as P5. Four Hermitian projectors

    P5, C P5 C^T, J P5 J^T, C^2 P5 (C^2)^T

generate M8(C). The unconstrained64-variable commutant equations have successive ranks30,52,58,63 modulo101. I is in their exact rational kernel, giving exact rank63 and scalar complex commutant. The finite-dimensional self-adjoint algebra decomposition proves fullness: a direct sum of matrix blocks with multiplicities has scalar commutant only when it is the single8x8 block of multiplicity1. Also

    ||[P5,C P5 C^T]||_HS^2=95/128 !=0.

These transported marks are not a common commuting partition. A scalar trace-preserving invariant average retains5I/8, whereas coherent products retain the full index algebra. This makes the required operator information precise. It does not identify eight current components with eight physical states or infer that C,J are symmetries of the Round20 H. No Clifford/CP channel or common physical clock time is inferred from this index-space theorem.

## 12. Consequence and remaining first premise

The source's original discrete grade, the internal matrix factor, its integer carry, their equilibrium correlations, and the exact response-memory law now have one explicit dictionary. Complete prescribed grade histories have an internal GNS realization; no external register or independent reset is needed for their mathematical construction. Actual joint clock access requires coherently transported index observables, not just a fixed average.

The earliest unsupplied physical identification is still the selection and realization of the original character k and its geometric/field/instrument dictionary from the primitive TFPT seam. It is already present in the carrier sources and can reconstruct P5, but the unmarked source H alone does not select it. The operator theorems above neither derive a primitive physical kernel nor choose beta, geometry, apparatus, allowed interventions, full interacting3+1D physics or gravity. The contract remains PARTIAL. These boundaries are consequences of the explicit source maps, not reasons to replace the common-origin target with an arbitrary finite simulator.

## 13. Evidence and references

Owning originals: charged-lift-round16; carrier-module-conjugation; local-charge-transport-round20; hopping-source-domain-round20; scalar-charge-energy-round20; history-sum-control-round26; native-source-reconstruction-20261014; compiler-integral-triality-20260918; charged-sector-readout-clock-20261013. Exact source hashes are registered in source_manifest.json.

`checker.py` checks original lattice grades, all240 roots, the15 moments,720 cocycle-corrected source isometries and the native commutator witness. `carry_check.py` tests the exact bilateral integer action and matrix units, including negative and trillion-scale carries. `clock_probe.py` supplies an exact rational commutant certificate. These checks do not simulate a truncated Gibbs state or numerically compute the infinite memory kernel. Analytical proofs are independently reviewed in MEMORY-REVIEW.md and MARK-REVIEW.md.

The projection-memory method is standard: Hazime Mori, Transport, Collective Motion, and Brownian Motion, Progress of Theoretical Physics33 (1965),423–455, https://doi.org/10.1143/PTP.33.423. The source-specific applications here do not constitute a new general theorem of statistical mechanics.
