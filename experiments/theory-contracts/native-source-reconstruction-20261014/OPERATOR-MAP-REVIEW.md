# Operator-Map Review: E8 charge shifts, finite Clifford quotients, and marked time

## Findings

**The source-defined E8 shifts do not yet have a demonstrated identification with the four-bit Clifford/compiler module.** There is already a precise obstruction to identifying the *entire charged commutator form* with that four-bit quotient: the E8 form is nondegenerate on eight binary coordinates, while the quotient has only four. A form-preserving linear map from the former into the latter would have to be injective, which is impossible by dimension. The existing charged-cocycle analysis confirms this at the level of actual roots: two roots identified by the compiler quotient have different cocycle commutators against a third root. This rules out that quotient as a lossless encoding of all charged commutator signs. It does not rule out a selected subgroup, added coordinates, or a different map with extra structure.

The distinction is larger for the local operators: (L/2L) captures a finite sign algebra, whereas (T_{x,p}) acts on the full infinite charge lattice ℓ²(E8). Its integer carry and nonzero shifts (T_{x,4n_s}≠ I) survive. The finite Clifford representation therefore cannot, by itself, be the full local charge-transport representation. The intrinsic E8-to-integral-Clifford construction is positive mathematical progress, but its own proof explicitly leaves the map to physical charges, group actions, and state rules to be supplied.

## Three-site neutral-transfer consequence

The existing local-transport proof gives, for distinct sites (x,y,z),

\[
W_{xy,p}W_{yz,r}=(-1)^{p^TGr}W_{yz,r}W_{xy,p},
\]

and same-bond transports commute with one another. Take (p=e_i,r=e_j), and set (A_i=W_{xy,e_i}), (B_j=W_{yz,e_j}). Since the E8 Gram matrix (G\bmod 2) is alternating and nondegenerate (rank eight), the (A_i) form a commuting family, the (B_j) form a commuting family, and their cross-commutation sign matrix has rank eight.

In any finite-dimensional unitary representation of these 16 bond generators, simultaneous diagonalization gives a common eigenvector of all (A_i). Acting on it with the 256 products (B_1^{v_1}\cdots B_8^{v_8}), (v\in\mathbb F_2^8), flips the joint (A)-eigenvalue character by (Gv). Nondegeneracy makes all 256 characters distinct, hence their eigenspaces are orthogonal. Therefore the representation dimension is at least (2^8=256). The bound is sharp for the *finite sign-commutation algebra* (eight Pauli (Z_i)'s paired with eight independent Pauli (X)-products).

This is a deduction from the source's commutation relation, not a theorem stated in that source. It does **not** exhibit a 256-dimensional faithful representation of the original infinite lattice shifts: quotienting away integer carry can collapse (W_{xy,4e_i}\ne I). The original Hilbert space remains infinite-dimensional. The lower bound says that an exact finite representation retaining these three-site sign relations cannot be the 16-dimensional four-bit module.

More generally, for the neutral mod-2 profile space
\[
K=\{(n_1,\ldots,n_N)\in(\mathbb F_2^8)^N:\sum_x n_x=0\}
\]
with pairing \(\Omega(n,m)=\sum_x n_x^T Gm_x\), its rank is (8(N-1)) for odd (N) and (8(N-2)) for even (N). To see the parity factor, the scalar site pairing is the standard dot product restricted to \(\mathbf1^\perp\): its radical is zero for odd (N), and is \(\langle\mathbf1\rangle\) for even (N). Thus for (N=3) the neutral form has rank 16 and its irreducible projective representation has dimension (2^{16/2}=256), agreeing with the bond-generator argument. This rank statement assumes the displayed inherited pairing; it does not add a physical sector-selection rule.

## Marked-time no-go scope

The September 22 contract settles a different, specific identification. For its original M4 channel \(\Phi\), every positive integer power \(\Phi^n\) has no time-independent GKSL logarithm on the same M4, with no covariance assumption on intermediate times. The proof handles all logarithm branches. It does not exclude the original GNS unitary time, time-dependent GKSL evolution, or a larger environment, and it does not settle the TFPT physical time law.

## Additional source-specific corollary: parity does not determine this hopping response

For the declared finite periodic cubic source with (L\ge3), choose a positive edge (x\to y), (p=e_8), (e(p)=1), and (r_x=-p,r_y=p), zero elsewhere. The original charge sector is (n_k=-2kr). Every (n_k) is neutral and has the same onsite (L/2L) parity (zero at every site); (W=W_{xy,p}) shifts it to (n_k+r). This uses the full integer-charge model, not a parity quotient.

With (e(n)=n^TGn/2) and the six-neighbor cubic penalty (D(n)=N^{-1}\sum_xe(n_x)+\nu\sum_{\langle ab\rangle}e(n_a-n_b)),
\[
D(r)=\frac{2e(p)}N+14\nu e(p),\qquad
D(n_k+r)-D(n_k)=(1-4k)D(r).
\]
The coefficient (14) is (4) from the shared bond (xy), where the charge difference is (2p), plus (5+5) other bonds incident to (x,y). For the scalar-charge coupling (V_{sc}=(\lambda/N)\sum_a e(n_a)\phi_a^2), the same quadratic calculation gives
\[
V_{sc}(n_k+r)-V_{sc}(n_k)
=(1-4k)\frac{\lambda}{N}e(p)(\phi_x^2+\phi_y^2).
\]
Thus the prompt's (u) is (\lambda/N) in the Round20 source convention. The continuous part of the Hamiltonian commutes with (W). Conjugating any hopping term by (W) leaves its nonzero integer translation vector unchanged (only its phase changes), so its diagonal charge matrix element in the basis vector at charge profile `n_k` is zero; the common 2I term is unchanged. On the finite-charge scalar Schwartz core this yields
\[
R'_k(0)=i(1-4k)\left[\frac{2e(p)}N+14\nu e(p)
+\frac{\lambda}{N}e(p)\langle\psi,(\phi_x^2+\phi_y^2)\psi\rangle\right],
\]
for (R_k(t)=\langle n_k,\psi|W^*e^{itH}We^{-itH}|n_k,\psi\rangle). Since (e(p)=1), (N<\infty), and \(\lambda,\nu\ge0\), the bracket is strictly positive. Distinct integers (k) therefore give distinct initial slopes—and hence distinct short-time future responses—despite identical onsite parity and identical scalar wavefunction.

This is a source-specific response corollary for the declared Round20 hopping Hamiltonian. It establishes that onsite parity alone is not a future-complete state descriptor under this (W)-dressed readout. It does not exclude a selected partial readout, further state projection, or a different physical dynamics; nor does it promote the declared hopping coupling to a derived TFPT law. The unbounded-energy derivative is evaluated only on the stated common core, where the source proves the domain conditions.

## Source trail

1. [`local-charge-transport-round20/PROOF.md`](</Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/local-charge-transport-round20/PROOF.md:48>) — E8 cocycle and local shifts, lines 48–90; bond transport and positivity, 92–131; commutation relations used above, especially equation (8), lines 133–150. This source establishes the relations and their infinite-lattice scope, not the derived 256-dimensional lower bound.
2. [`charged-cocycle-lift/README.md`](</Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/charged-cocycle-lift/README.md:1>) — eight-coordinate (L/2L) charged sign refinement and 16-dimensional irreducible sign representation; lines 1–35; actual root counterexample to the four-bit compiler quotient and the nondegenerate form argument, lines 37–100. It explicitly distinguishes this finite sign algebra from complete integer charges/fields.
3. [`source-continuous-marked-time-search-20260922/PROOF.md`](</Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/source-continuous-marked-time-search-20260922/PROOF.md:27>) — all positive integer powers and all logarithm branches, lines 27–39 and 66–124; retained alternatives and limits, lines 178 onward. [`compiler-triality-rigidity-20260918/PROOF.txt`](</Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/compiler-triality-rigidity-20260918/PROOF.txt:105>) separately records the intrinsic E8-to-Clifford result and explicitly states the remaining physical identification at lines 118–126.
4. [`hopping-source-domain-round20/PROOF.md`](</Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/hopping-source-domain-round20/PROOF.md:3>) — exact common source Hamiltonian and (D,V,L_J), lines 9–23; domain/core, lines 48–59. [`scalar-charge-energy-round20/PROOF.md`](</Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/scalar-charge-energy-round20/PROOF.md:1>) — projective shift calculus and quadratic energy increments, lines 1–32; scalar coupling and operator-domain scope, lines 35–70 and 109–130. The derivative formula above follows directly from these definitions and the local charge-transport relation (8); it is a corollary, not a separately reported source theorem.

The decisive next mathematical check for a proposed finite map is not matching dimensions: verify that its induced map preserves the actual E8 mod-2 commutator form on the chosen charge subspace and also accounts for integer carry. The three-site neutral subspace already rules out representing all these bond-sign degrees of freedom in the four-bit module without extra structure.
