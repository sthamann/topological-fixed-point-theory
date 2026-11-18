# Native source reconstruction and the information needed by neutral transport

14 October 2026. Research contract, **PARTIAL**. No T1–T8 physical gate is closed.

The target is a common origin of marking, state and event law. This work uses the *existing declared* finite-volume g=0 charge/scalar source from Round20–26 without replacing its hopping words, cocycle signs, scalar modes, or charge integers. It proves a constructive source-to-time-and-response connection at fixed spatial lattice, and a necessary information bound for an exact finite representation of its local neutral transport. It does not derive this source's spatial lattice, couplings, temperature or physical marking from P1/P2.

## 1. Source, assumptions, and prior results

Fix the Round26 finite periodic spatial lattice, N=L^3, L>=3, a,m>0, u,nu,J>=0. The neutral charge configurations are C0={n=(n_x): n_x in Z^8, sum_x n_x=0}. Work on

    H0 = l2(C0) tensor L2(R^N).

The actual integral even unimodular Gram matrix G and the frozen bilinear lattice cocycle c are those of `local-charge-transport-round20`. Let

    e(n)=n^T G n/2,
    D(n)=sum_x e(n_x)/N + nu sum_edges e(n_x-n_y),
    A(n)=-Delta_a + m^2 I + 2u diag(e(n_x)),
    K=sum_x pi_x^2/2,
    V(n,phi)=D(n)+phi^T A(n)phi/2.

The native two-site unitaries are W_xy,p=T_y,p T_x,p*, p in the eight integral basis directions. Their phases are retained. The hopping term is

    Hhop=J sum_(positive edges,p)(2I-W_xy,p-W_xy,p*) >= 0.

It is bounded at fixed N, with norm <=96JN. The constant 48JN is retained. There is no conversion of eight lattice coordinates into eight particle species or spacetime dimensions.

Round26 supplies the positive trace-class transfer

    Q_delta = e^(-delta V/2) e^(-delta (K+Hhop)) e^(-delta V/2), delta>0.

K and Hhop act on different factors. Round25 supplies the scalar determinants and their nonzero time-correlated response. Round26 proves absolute convergence of the *signed* Euclidean history sum and a positive total partition function. These are prior results, not newly discovered here. A signed history expansion is not a classical event distribution.

The earlier `charged-cocycle-lift` already establishes the minimal onsite sign quotient L/2L and its 16-dimensional irreducible representation; the four-bit Gaussian compiler quotient is not a quotient of the complete charged sign algebra. The older `source-continuous-marked-time-search-20260922` already excludes every homogeneous GKLS logarithm of the marked M4 channel and its positive integer powers. Our previous 14 October rank argument is an alternative scoped proof, not the first project proof of that broader exclusion.

## 2. The existing source has a common self-adjoint Hamiltonian

At each charge configuration, A(n)>=m^2 I and the scalar operator is a finite-dimensional confining oscillator. Define the direct sum

    Hdiag = direct_sum_n [K + phi^T A(n)phi/2 + D(n)].

This is positive self-adjoint. Finite charge support tensored with scalar Schwartz functions is a core. The original full operator

    H = Hdiag + Hhop

is positive self-adjoint on Dom(Hdiag), by bounded self-adjoint perturbation. The same core applies. Equivalently it is the densely defined positive form sum of V and K+Hhop. No new oscillator, pulse, bath supply, coupling or reset is introduced in this step.

G is positive definite and D(n)>=||n||^2/(62N), while every oscillator frequency is at least m. Thus Hdiag has compact resolvent. Hhop preserves compact resolvent. At every prescribed beta>0 the thermal state exists, with

    0 < Z(beta)=Tr exp(-beta H)
      <= Tr exp(-beta Hdiag)
      <= Theta(beta) [2 sinh(beta m/2)]^(-N) < infinity,

where Theta(beta)=sum_(n in C0)exp(-beta D(n)). The middle inequality is an eigenvalue/min-max trace bound, not an invalid assertion that the exponential is operator monotone. Consequently

    rho_beta = exp(-beta H)/Z(beta)

is faithful, invariant under U(t)=exp(-itH), and a Gibbs/KMS state for the same H. This is a thermal family conditional on beta; it is not a proof that the physical source chooses beta or a unique pure ground state. Frustrated cocycle signs are not bypassed by a Perron–Frobenius assumption.

## 3. Reconstruct time from the actual transfer

### 3.1 A fixed regulator

All three positive factors in Q_delta are injective contractions: their spectral functions e^(-delta x) are positive at every finite x. For any psi !=0, e^(-delta V/2)psi !=0 and the quadratic form against e^(-delta(K+Hhop)) is strictly positive. Hence ker Q_delta=0. Its spectrum may approach zero; no bounded inverse is assumed.

Functional calculus defines the densely defined positive self-adjoint operator

    H_delta = -log(Q_delta)/delta,
    U_delta(t) = exp(-it H_delta) = Q_delta^(it/delta).

The domain consists of vectors with square-integrable log(lambda)/delta in their Q_delta spectral measure. This is an all-real-time unitary group of the *specified regulated operator*, not a logarithm of a reduced quantum channel. If beta=T delta with integer T>=1, its thermal density is precisely

    rho_(beta,delta)=Q_delta^T / Tr(Q_delta^T).

Thus the original Euclidean partition state and this regulated time are generated by the same operator. At fixed delta one must not silently replace H_delta by H.

### 3.2 Remove only the temporal regulator

Put S_delta=(I-Q_delta)/delta. It is bounded positive self-adjoint for each delta. On the common finite-charge/Schwartz core,

    S_delta psi -> (V+K+Hhop)psi=H psi.

This follows by telescoping the three semigroup factors and applying their strong derivatives to core vectors. The standard common-core convergence theorem for self-adjoint operators yields S_delta -> H in strong resolvent sense.

The logarithmic generator has the same limit. For z>0, put x=delta h>=0 and b=delta z. Spectral calculus gives

    0 <= (S_delta+z)^(-1) - (H_delta+z)^(-1) <= delta I.

Indeed the scalar difference is

    delta [1/(1-exp(-x)+b) - 1/(x+b)].

It is nonnegative and bounded by delta: adding b only reduces the difference, and

    0 <= 1/(1-exp(-x)) - 1/x <= 1

for x>0, using exp(x)>=1+x. The expression at x=0 is interpreted by continuity with b>0. Therefore H_delta -> H in strong resolvent sense. Self-adjoint functional calculus then gives

    U_delta(t)psi -> exp(-itH)psi

strongly, uniformly for t in each compact real interval.

This closes the mathematical transfer-to-real-time connection for the original declared source at fixed finite spatial volume. It is not a spatial continuum limit, Lorentz theorem, selection of the primitive source, or a uniform efficient numerical error bound. No operator-norm convergence or trace-norm convergence of rho_(beta,delta) is inferred solely from strong resolvent convergence.

### 3.3 The thermal states converge in trace norm

This stronger conclusion uses the original confining Gaussian source, not strong resolvent convergence alone. Remove the hopping exponential from the middle factor to obtain Q_(delta,0). Since K and Hhop commute and Hhop>=0,

    0<Q_delta<=Q_(delta,0).

The logarithm is operator monotone. Apply it first to Q_delta+epsilon I and Q_(delta,0)+epsilon I; their inverses are now bounded. Taking the monotone closed-form limit epsilon down to zero gives

    H_delta >= H_(delta,0) = -log(Q_(delta,0))/delta.

For a fixed charge n, diagonalize A(n) orthogonally. In each scalar eigenmode of frequency omega, direct comparison with the Mehler kernel gives

    alpha = 2 asinh(delta omega/2),
    -delta^(-1)log[e^(-delta omega^2 q^2/4) e^(-delta p^2/2)
                   e^(-delta omega^2 q^2/4)]
       = (a_delta p^2+b_delta q^2)/2,
    a_delta=alpha/sinh(alpha),  b_delta=alpha sinh(alpha)/delta^2.

No additive constant is omitted: both kernels have prefactor (2 pi delta)^(-1/2), cross coefficient 1/delta and diagonal coefficient cosh(alpha)/(2 delta). The scalar D(n) is retained separately and exactly. Using asinh(y)>=y/sqrt(1+y^2),

    a_delta >= 1/(1+delta^2 omega^2/4),   b_delta >= omega^2.

Fix delta_0>0 and define a single comparison operator, independent of delta:

    c_n = 1/(1+delta_0^2 ||A(n)||/4),
    B = direct_sum_n [D(n)+c_n K+m^2 |phi|^2/2].

For every 0<delta<=delta_0, H_delta>=B as closed forms. Also H>=B. The n-block of B has N oscillator frequencies m sqrt(c_n), so for every gamma>0,

    Tr exp(-gamma B)
       = sum_n exp(-gamma D(n))
         [2 sinh(gamma m sqrt(c_n)/2)]^(-N) < infinity.

Indeed ||A(n)|| grows at most quadratically in the charge coordinates; c_n^(-N/2) grows only polynomially. The coercive quadratic D(n) dominates this growth. This proves compactness and Gibbs trace summability even though c_n approaches zero at large charge.

By min-max and closed-form order, R_delta=exp(-beta H_delta) obeys

    Tr R_delta <= Tr exp(-beta B),
    Tr(B^(1/2) R_delta B^(1/2))
       <= Tr(H_delta exp(-beta H_delta))
       <= 2/(e beta) Tr exp(-(beta/2) B) = C_beta.

For the fixed finite-rank spectral projection P_R=1_[0,R](B),

    Tr[(I-P_R)R_delta] <= C_beta/R.

The same estimates hold for R=exp(-beta H). Strong semigroup convergence gives convergence in each fixed finite-dimensional P_R block. Uniform tails then give trace-norm convergence: for any positive trace-class S, writing t=Tr[(I-P_R)S],

    ||S-P_R S P_R||_1 <= t+2 sqrt(Tr(S)t).

This follows by factoring the off-diagonal blocks through S^(1/2) and applying the Hilbert-Schmidt Cauchy-Schwarz inequality. First take delta to zero at fixed R, then R to infinity. Consequently

    ||exp(-beta H_delta)-exp(-beta H)||_1 -> 0,
    Z_(beta,delta) -> Z_beta > 0,
    ||rho_(beta,delta)-rho_beta||_1 -> 0.

The proof works for every fixed beta>0 as delta down to zero; in particular it applies to the original integer slicing delta=beta/T. It also establishes that Q_delta^(beta/delta) is trace class for every beta>0, beyond integer powers. No thermodynamic or spatial-continuum limit is taken.

## 4. One source produces consistent ordered responses

For a prescribed finite family of bounded insertions, the same H and rho_beta define the ordered quantum response. For example, with finite charge-cylinder partitions {P_a} and ordered times, let

    C_h=P_(a_r)(t_r) ... P_(a_1)(t_1),
    P_a(t)=U(t)* P_a U(t),
    D(h,h')=Tr(C_h rho_beta C_h'*).

D is a positive semidefinite Gram matrix: for any complex coefficients z_h,

    sum_h,h' conjugate(z_h) D(h,h') z_h'
      = Tr(B rho_beta B*) >=0, B=sum_h conjugate(z_h) C_h.

Equivalent conjugation conventions give the same positivity statement. For complete finite partitions at each time, sum_h C_h=I, so sum_h,h' D(h,h')=1. For an actually specified sequence of ideal projective instruments, sum_h D(h,h)=1 by repeated completeness and trace cyclicity. The diagonal is its sequence-of-measurement probability law. Coherent coarse graining sums both indices; deleting off-diagonal D is not generally valid for unmeasured histories.

The projector partitions are chosen probes of already existing charge operators, not a derivation of which measurements occur physically. Likewise replacing projectors by a particular six-Kraus instrument is an additional choice until an operator dictionary derives it. The construction supplies all responses to stated probes without introducing stochastic resets.

For a fixed finite number of bounded insertions, locally uniform strong convergence of U_delta and its adjoint gives strong convergence of each history operator and its adjoint. Section 3.3 also gives convergence of the actual regulated states in trace norm. Hence

    Tr(C_(h,delta) rho_(beta,delta) C_(h',delta)*)
      -> Tr(C_h rho_beta C_h'*).

All fixed finite bounded ordered responses and the entire finite decoherence matrix converge using each regulator's own Gibbs state. This supplies a common state/time/response limit from the existing source; it does not select beta, the probes or the physical occurrence of their outcomes.

## 5. What a lossless local finite map would have to retain

### 5.1 Two adjacent neutral bonds

Choose distinct neighboring x,y,z and set A_i=W_xy,e_i, B_j=W_yz,e_j. The original exact word action gives

    A_i A_j=A_j A_i,
    B_i B_j=B_j B_i,
    A_i B_j=(-1)^(G_ij) B_j A_i.

These identities preserve total charge and apply inside the neutral sector; they are not the onsite charged operators misread as neutral. G mod2 is invertible because det G=1.

Suppose these relations admit a nonzero finite-dimensional unitary representation on C^d. Pick a joint eigenvector v of all A_i, with nonzero unit-modulus eigenvalues a_i. For b in F2^8, let B^b=product_j B_j^(b_j). Then B^b v is a joint eigenvector with character

    (a_i (-1)^((G b)_i))_i.

All 256 characters are distinct because G is invertible mod2. The corresponding vectors are nonzero and mutually orthogonal. Therefore d>=256. Decomposing the joint spectrum into these free orbits also shows 256 divides d.

The bound is sharp for the abstract pairwise commutator algebra: on 256 binary basis states t, take A_i=X_i and B_j=Z^(G_column_j). These are only a finite representation of the sign relations. They do not replace the unbounded integer translations, their nontrivial even carries, electric energies or original Hamiltonian.

Accordingly no unitary word intertwiner into M4 or M16 can preserve this entire pair of adjacent neutral bond algebras. A selected partial readout or a nonhomomorphic CP compression remains possible. If used, its lost information and memory must be kept explicit. This is not a lower bound on a single particle's internal dimension.

### 5.2 General neutral sign rank

At N connected sites, the neutral mod2 vectors satisfy sum_x n_x=0. Parametrize by n_1,...,n_(N-1), with n_N=sum n_x. Their commutator form is

    Omega = (I_(N-1)+11^T) tensor (G mod2).

For odd N the first factor is nonsingular; for even N its radical has dimension one. Hence

    rank Omega = 8(N-1)  if N is odd,
                 8(N-2)  if N is even.

The associated minimal finite irreducible sign representation has dimension 2^(rank(Omega)/2), after fixing the center. Connected bond transfers generate the neutral module. The argument concerns sign information only; the full integer charge/state algebra remains infinite-dimensional. The N=3 witness above already suffices to reject a lossless M4 identification.

### 5.3 Equal sign labels do not determine future native responses

The existing carry obstruction has a direct dynamical witness for this source. Choose adjacent x,y, p=e8 with e(p)=1, and r=(-p at x,+p at y). Let n_k=-2kr, k any integer. These configurations are all neutral and have the same zero onsite L/2L label. Use the same normalized scalar Schwartz vector psi for all k.

The native W=W_xy,p translates n to n+r. Quadraticity gives

    V(n_k+r,phi)-V(n_k,phi)
       = (1-4k)[D(r)+u(phi_x^2+phi_y^2)],
    D(r)=2/N+14 nu.

On the periodic cubic lattice with L>=3, the common bond contributes e(2p)=4; the other ten incident bonds contribute 10. This explains 14 without changing the original neighbor convention. The scalar free part cancels. K commutes with W. Every nonconstant hopping word, including its conjugate by W, shifts integer charge nontrivially; thus its diagonal matrix element in n_k vanishes. The constant hopping contribution cancels between W*Hhop W and Hhop.

Consequently the actual ordered native response

    R_k(t)=<n_k,psi | W* exp(itH) W exp(-itH) | n_k,psi>

has the exact derivative on the common core

    R_k'(0)=i(1-4k)[2/N+14 nu+u <phi_x^2+phi_y^2>_psi].

The bracket is strictly positive. Therefore infinitely many source inputs with exactly the same parity data have different future responses already at first order in time. This proves that onsite parity alone is not a future-equivalent state description for these native probes. It is a source-specific consequence of existing integer carries, not a new universal no-go for every compressed readout. In particular a projection to selected states or a restricted observable set changes the question and needs its own preservation proof.

## 6. Consequence for the research route

The existing source supports a common Hamiltonian, a thermal family and all ordered response functionals. This is a constructive common-object connection rather than combining different model states. Its time regulator can be removed at fixed spatial lattice, including trace-norm convergence of its thermal states and all fixed finite bounded responses. At the same time the full local neutral algebra provably cannot be identified with the small marked channel as a lossless word representation.

The sign quotient also fails dynamically: equal parity inputs can give different native future responses, as section 5.3 proves. The useful route is therefore to identify a **specified physical subsystem or response algebra** inside the complete source, preserving the ordered responses relevant to the marks. The target is not a quotient justified solely by matching endpoint values. The actual selection of that algebra, its marking and physical state must follow from the original source principles; they are not fixed by the formulas above. The declared original source itself contains the finite lattice, eight arithmetic hopping directions, and coupling parameters. Those are inherited inputs, not newly derived P1/P2 consequences.

## 7. Verification and references

`neutral_transport_check.py` reads the original Gram assignments and tests actual unbounded integer-profile word actions with exact signs, the neutral ranks, and a sharp finite representation of the pairwise relations. Its finite checks accompany, but do not replace, the proofs for arbitrary profiles and arbitrary N. The three-site dimension proof was independently reviewed. `reconstruction_check.py` independently compares the original static neutral Gaussian block integral with the exact reconstructed oscillator spectrum, including the normalization. This is a finite numerical diagnostic, not a full moving signed-charge history simulation. Both reconstruction and the uniform Gibbs-tail proof were independently reviewed in `TRANSFER-REVIEW.md`.

Original project sources: `local-charge-transport-round20/PROOF.md`, `history-sum-control-round26/PROOF.md`, `moving-history-memory-round25/PROOF.md`, `charged-cocycle-lift/README.md`, and `source-continuous-marked-time-search-20260922`.

The functional-calculus/common-core convergence facts are standard self-adjoint operator theory, not claimed as new general mathematics. A primary mathematical review of the underlying self-adjoint/form-sum and semigroup framework is Barry Simon, *Tosio Kato's work on non-relativistic quantum mechanics*, https://arxiv.org/abs/1711.00528. The distinction between process responses and reduced channel powers is treated by Pollock et al., https://arxiv.org/abs/1801.09811.
