# Round20–26: what the native Hamiltonian remembers about a 5+3 marking

## Result

The specified unmarked Round20–26 Hamiltonian and its Gibbs state do **not** select one of the 15 (5+3) characters (k_S). Under the proposed ambient-coordinate (S_6), all 15 markings are carried into one another while the source Hamiltonian, its transfer, and its normalized thermal state are preserved. This is a source-specific covariance statement, not a no-go for the complete TFPT: the full construction also supplies fixed (J,\sigma,C), and a chosen (k), and those extra data are not part of this Round20–26 Hamiltonian. In particular, the already-derived recovery of (P_5) from a **given** (k) remains valid.

## Source ingredients and the action

The Round20–26 source takes the E8 Gram matrix (G=B^TB), frozen bilinear cocycle

\[
c(p,n)=\sum_a p_an_aG_{aa}/2+\sum_{a>b}p_an_bG_{ab}\pmod 2,
\]

and the eight integral basis hopping directions. Its diagonal charge/scalar terms depend on (e(n)=n^TGn/2), and its hopping is the equal-coupling sum over all eight directions; see `local-charge-transport-round20/PROOF.md` §§2–3 and `native-source-reconstruction-20261014/PROOF.md` §§1–2. The latter proves that this exact declared source has a self-adjoint (H), a trace-class Gibbs state for each prescribed \(\beta>0\), and the associated transfer/time construction. It expressly does not derive a physical marking from P1/P2.

Use the physical basis columns stated for this check:

\[
b_0=2e_0,\quad b_i=e_0+e_i\ (1\le i\le6),\quad b_7=(1/2)^8.
\]

An ambient permutation \(\tau\in S_6\) of coordinates 1 through 6 fixes (b_0,b_7) and permutes the six (b_i). In the lattice coefficient basis it is a permutation (P_\tau) of entries 1 through 6 fixing 0 and 7. Hence (P_\tau^TGP_\tau=G). It preserves every diagonal term (e,D,A), fixes the scalar-field coordinates, and permutes the eight hopping generators (the two fixed generators stay fixed). Since the hopping has the same coefficient (J) in all eight directions, (V_\tau H V_\tau^*=H).

The cocycle is not strictly fixed in the chosen triangular convention, but its change is a coboundary. Let \(I(\tau)\) be the pairs \(i<j\) whose order is reversed by \(\tau\), and set

\[
f_\tau(x)=\sum_{\{i,j\}\in I(\tau)}x_ix_j\pmod2.
\]

Only the off-diagonal entries (G_{ij}=1\) among the six root columns change order, so direct expansion gives

\[
c(P_\tau x,P_\tau y)+c(x,y)
=\sum_{\{i,j\}\in I(\tau)}(x_i y_j+x_j y_i)
=f_\tau(x)+f_\tau(y)+f_\tau(x+y)\pmod2.
\]

Thus the lattice-basis lift \(\widetilde P_\tau|q\rangle=(-1)^{f_\tau(q)}|P_\tau q\rangle\) intertwines the cocycle translations up to the generator phase \((-1)^{f_\tau(p)}\). For a neutral bond (W_{xy,p}=T_{y,p}T_{x,p}^*\), that phase occurs once in each factor and cancels. Therefore \(W_{xy,p}\mapsto W_{xy,P_\tau p}\); summing over all eight (p) leaves (H_{hop}) invariant. This establishes the claimed symmetry with the actual frozen cocycle, rather than assuming that the triangular cocycle itself is permutation invariant.

Consequently the normalized Gibbs state \(\rho_\beta=e^{-\beta H}/\mathrm{Tr}(e^{-\beta H})\) is invariant under the lifted action. The same applies to the source transfer built from these covariant terms. This means no covariant, unique construction using only this unmarked source and its thermal state can return one preferred member of the 15-element orbit without adding a non-invariant input.

## The 15 marks are actually distinct

For \(|S|=4\), define (k_S(q)=2(q_0+\sum_{i\in S}q_i)\pmod4\). Under (P_\tau), (k_S\circ P_\tau=k_{\tau^{-1}S}), so the orbit comprises all \(\binom64=15\) choices. The marks are distinct on native root directions: the charge along (b_i) has grade 2 exactly when (i\in S\), and grade 0 otherwise. The fixed long-root and half-sum directions do not distinguish the choices. Thus this is not merely a relabelling invisible on the charge lattice.

For each (S), the corresponding neutral-root second moment is the rotated version of the familiar (16P_5+8P_3); the source Hamiltonian sees only the common E8 norm and the complete equal-weight direction sum, so that moment is not an observable selector in this source. The stabilizer of a chosen (S) in this ambient action is (S_4\times S_2); specifying the mark itself breaks the orbit symmetry.

## Where the full TFPT can break the covariance

The conclusion stops at the Round20–26 unmarked source. The full TFPT carries additional fixed structure:

- `charged-cocycle-lift/README.md` §1 defines (J) as rotation of each coordinate pair and \(\sigma\) as the cycle of the first three pairs. These are actual E8 lattice isometries whose coherent lifts require phase corrections (§4); they are not derived from the scalar (H) above.
- `carrier-module-conjugation/README.md` §1 proves that a *supplied* (k(q)=2\sum_{i=1}^5q_i\pmod4) recovers (D_5\oplus D_3) and (M_0=16P_5+8P_3). It explicitly says this (k) is an input, not derived there. Its §2 tests a particular (K) that preserves (k,\sigma) while reversing (J), showing that preserving the carrier/mark does not by itself preserve oriented (J).
- The (\sigma)-selected Arf refinement in `verification/v774_arf_spinor_compiler.py` belongs to the compiler quotient (V=L/(1+J)L\cong\mathbb F_2^4\). Its reported (S_6\cong Sp(4,2)) action is the action on the six Arf-1 refinements, not automatically the ambient-coordinate (S_6) used here. The two permutation groups must not be identified by matching their names or orders.

Since the fixed full-source data can be moved by the candidate \(\tau\), the Round20 covariance alone neither proves that the full TFPT has 15 equivalent sectors nor closes the selection problem. To decide which elements survive, one must transport the actual (J,\sigma,C,k) in one common lattice basis and test their intertwining equations. The supplied (k) already reduces the 15-way ambiguity to its stabilizer; the remaining mark is then an explicit source input, not something inferred from the unmarked thermal response.

## Narrow conclusion

The first (5+3) split is not encoded in the Round20–26 scalar energy plus equal eight-channel hopping and Gibbs state: those data are covariant across the 15 (k_S). The full marked TFPT may distinguish them through its extra fixed structures, and existing carrier results explain what follows **once (k) is given**. The exact unresolved identification is whether those full-source structures are themselves derived compatibly from the same Round20 source, not whether E8 has a unique unmarked (5+3) subspace.

## Joint clock probe on the eight-dimensional index space

I checked the actual C,J vector matrices from compiler-integral-triality-20260918/certificate.json (SHA-256 `a342f865bec164ae6dd54e9e7b7c7cc991efcf7a939dda2ecab586b4a6c6cfe3`) against the fixed coordinate projector

    P5 = diag(1,1,1,1,1,0,0,0).

The four projectors are P5, C P5 C^T, J P5 J^T, and C^2 P5 (C^2)^T. I imposed no symmetry restriction on an unknown 8x8 matrix X: for each projector Q, the linear constraints are XQ-QX=0, with 64 complex scalar unknowns. The successive constraint ranks modulo 101 are 30, 52, 58, and 63. An exact rational nullspace computation gives a one-dimensional kernel, spanned by I8, so the exact rational rank is 63. Independently,

    ||[P5, C P5 C^T]||_HS^2 = 95/128.

Thus even the first actual clock conjugation moves the fixed projector, while the four-projector joint commutant is scalar. The unital complex *-algebra generated by these four Hermitian projectors therefore has scalar commutant and equals M8(C). Its coherent operator record is not exhausted by the scalar invariant part: averaging a trace-5 projector over the compact closure of the generated orthogonal clock group gives the unique invariant coefficient 5 I8/8, which erases non-scalar distinctions that products/words in the projectors retain.

This is only an eight-dimensional **index-space** statement. These matrices are not eight physical states and do not establish physical time. Nor does the calculation prove C or J is a symmetry of the Round20–26 Hamiltonian; the earlier S6 covariance result concerns that unmarked Hamiltonian and must remain separate from the fixed-clock action. The reproducible exact check is clock_probe.py, with matrices, ranks, exact kernel, norm, and certificate fingerprint in clock-result.json.
