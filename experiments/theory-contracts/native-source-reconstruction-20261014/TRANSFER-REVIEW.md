# Transfer-log reconstruction on the native finite lattice

**Independent review, 2026-10-14.** This note reviews only the neutral `g=0`
scalar/charge model on the fixed finite periodic spatial lattice used in
Round 20 and Round 26. It does not add a new Hamiltonian or a new source law.

## Verdict

The proposed fixed-regulator construction is mathematically valid, with two
qualifications.

1. For every fixed `delta > 0`, the Round-26 transfer operator determines a
   unique nonnegative self-adjoint effective Hamiltonian

   \[
   H_\delta=-\delta^{-1}\log Q_\delta,
   \qquad Q_\delta=e^{-\delta H_\delta}.
   \]

   It also determines an exact faithful Gibbs state at `beta = T delta`, an
   exact strongly continuous unitary group, and a strongly positive
   decoherence functional for any *prescribed* finite family of native charge
   projectors. This is a conditional source-to-state/time/response result at
   the stated regulator.

2. The effective Hamiltonian is not merely an unrelated replacement for the
   Round-20 Hamiltonian. With

   \[
   H=K_\phi\dotplus V+H_{\rm hop}
   \]

   denoting the existing nonnegative closed-form sum on the neutral reducing
   subspace, Kato's form-product theorem gives

   \[
   Q_{t/n}^{n}\xrightarrow[n\to\infty]{\rm strong}e^{-tH}.
   \]

   The associated logarithmic generators satisfy

   \[
   H_\delta\xrightarrow[\delta\downarrow0]{\rm strong\ resolvent}H,
   \qquad
   e^{-itH_\delta}\xrightarrow[\delta\downarrow0]{\rm strong}e^{-itH},
   \]

   with the unitary convergence locally uniform in real time. This is a
   mathematical real-time reconstruction at fixed finite spatial lattice.
   It is not a spatial-continuum, Lorentzian-causality, or physical-clock
   selection result.

   The original `H` also has a trace-class Gibbs operator for every
   `beta > 0`. Consequently its own faithful Gibbs state and native bounded
   history responses exist without taking a limit of the regulated density
   matrices.

The strong limits alone would not imply trace-norm convergence. However, the
original Round-26 Gaussian and charge coercivity provide the missing uniform
trace tail. The addendum below proves

\[
\|e^{-\beta H_\delta}-e^{-\beta H}\|_1\longrightarrow0,
\qquad
\|\rho_{\beta,\delta}-\rho_\beta\|_1\longrightarrow0
\]

for every fixed `beta>0` along admissible slices `delta=beta/T -> 0`.
Consequently finite bounded ordered responses converge while using their
respective regulated and limiting Gibbs states.

## Owning source statements

The exact Round-26 transfer is in
`source:0d7cc9da475ee4945fe8`,
`history-sum-control-round26/PROOF.md`:

- lines 3--8 fix the finite lattice, neutral sector, and Euclidean regulator;
- lines 12--31 define `K_phi`, `V`, the fully normalized positive hopping
  Hamiltonian, `Q_delta`, and its trace-class property;
- lines 52--84 prove absolute convergence of the signed word expansion and a
  positive partition-function lower bound while retaining the sign problem;
- lines 119--128 restrict the available approximation estimate to a bounded
  one-time charge observable and expressly decline arbitrary signed
  multi-time insertions;
- the final paragraph leaves real-time reconstruction open in that contract.

The matching original positive dynamics is not inferred from dimensions. It
is explicitly present in the Round-20 sources:

- `source:e39426b76acca89849ef`,
  `hopping-source-domain-round20/PROOF.md`, lines 17--31 defines
  `X = X_d + L_J` and retains the hopping constant; lines 48--60 prove
  self-adjointness, positivity, and the common core using the bounded
  perturbation theorem.
- `source:18d340aeed16115daafa`,
  `scalar-charge-energy-round20/PROOF.md`, lines 38--42 identifies the same
  scalar/charge Hamiltonian
  `H_sc = H_m + D + (lambda/N) sum e(n_x) phi_x^2 + L_J`; lines 94--137 give
  its closed positive local forms and qualified real-time evolution.
- `source:5959236bb77ca96b5633`,
  `local-charge-transport-round20/PROOF.md`, lines 72--90 constructs the
  native local charge spectral projectors and equal-time site net; lines
  94--108 proves `H_hop >= 0`, boundedness, and retention of the constant;
  lines 183--187 shows that total-charge spectral projectors strongly commute
  with the full Hamiltonian and all-time group.

Thus the neutral total-charge sector used in Round 26 is a reducing subspace
for every factor and for `H`; the product formula may be applied after this
restriction.

## Fixed-regulator proof

Write

\[
A_\delta=e^{-\delta V/2},\qquad
B_\delta=e^{-\delta(K_\phi+H_{\rm hop})},\qquad
Q_\delta=A_\delta B_\delta A_\delta.
\]

Both exponentials are positive injective contractions. Therefore

\[
0\le Q_\delta\le A_\delta^2\le I,
\qquad
\langle\psi,Q_\delta\psi\rangle
=\|B_\delta^{1/2}A_\delta\psi\|^2.
\]

The last identity shows `ker Q_delta = {0}`. Round 26 separately proves that
`Q_delta` is trace class. Hence it is a positive compact injective operator
with nonzero eigenvalues `q_j in (0,1]` tending to zero. Spectral calculus
then defines

\[
H_\delta\psi
=-\delta^{-1}\sum_j(\log q_j)
  \langle e_j,\psi\rangle e_j,
\]

on the domain where the squared weighted sum converges. This operator is
self-adjoint, nonnegative, and has compact resolvent. It is the unique
nonnegative self-adjoint operator whose `delta` heat step is `Q_delta`.

For `beta=T delta`,

\[
Q_\delta^T=e^{-\beta H_\delta},\qquad
Z_\beta=\operatorname{Tr}e^{-\beta H_\delta}\in(0,\infty).
\]

Because `Q_delta` is injective, the normalized density matrix

\[
\rho_{\beta,\delta}=e^{-\beta H_\delta}/Z_\beta
\]

has support `I`; it is therefore a faithful normal state on the bounded
operators of the neutral Hilbert space. The functional calculus expression

\[
U_\delta(t)=Q_\delta^{it/\delta}=e^{-itH_\delta}
\]

is unambiguous despite `0` belonging to the spectrum: the spectral
projection at the eigenvalue zero vanishes. It is a strongly continuous
unitary group. With
`alpha_t(A)=e^{itH_delta} A e^{-itH_delta}`,
`rho_{beta,delta}` is the standard `beta`-KMS Gibbs state.

The additive energy level must be retained. In particular, Round 20 and
Round 26 both keep the hopping constant `kappa=48JN`. Replacing `Q_delta`
by `c Q_delta` shifts `H_delta` by
`-(log c)/delta` times the identity. This leaves normalized Gibbs states
and Heisenberg dynamics unchanged but changes the partition function and
free energy. Such a rescaling is not licensed by the source.

## Identification with the existing Round-20 Hamiltonian

Let `A=V` and `B=K_phi+H_hop` on the neutral subspace. Both are nonnegative
self-adjoint, and their closed forms have a dense intersection. The
Round-20 positive polynomial Schrödinger construction identifies their
closed form sum with the declared scalar/charge Hamiltonian `H`. The
hopping contribution is bounded and does not change the form domain.

Kato's product theorem for closed nonnegative forms, applied to the
symmetric product, yields for every `t >= 0`

\[
\left(e^{-tV/(2n)}e^{-tB/n}e^{-tV/(2n)}\right)^n
=Q_{t/n}^{n}\longrightarrow e^{-tH}
\]

strongly, locally uniformly in `t` on compact nonnegative intervals. This
is exactly the unsliced heat semigroup of the already declared `H`; it is
not a newly fitted Hamiltonian.

There is also a direct generator proof which avoids extracting the
logarithmic limit from the product theorem. Define the bounded nonnegative
self-adjoint operator

\[
S_\delta=\frac{I-Q_\delta}{\delta}.
\]

On the Round-20 operator core of finite charge support and continuous
Schwartz vectors, the three semigroup factors are strongly differentiable
at zero. The elementary decomposition

\[
Q_\delta-I
=A_\delta B_\delta(A_\delta-I)
 +A_\delta(B_\delta-I)+(A_\delta-I),
\quad
A_\delta=e^{-\delta V/2},\quad
B_\delta=e^{-\delta(K_\phi+H_{\rm hop})},
\]

and contractivity of the factors give
`S_delta psi -> H psi` on that core. Since `H` is self-adjoint, the standard
core-convergence theorem for uniformly lower-bounded self-adjoint operators
gives `S_delta -> H` in strong resolvent sense.

It remains to compare `S_delta` with the requested logarithm. They are
functions of the same `Q_delta`. If a spectral value is written as
`q=e^{-x}`, `x>=0`, then for every `z>0`, with `a=delta z`,

\[
0\le (S_\delta+z)^{-1}-(H_\delta+z)^{-1}
\le \delta I.
\]

Indeed, the scalar difference divided by `delta` is

\[
\frac1{1-e^{-x}+a}-\frac1{x+a}
\le \frac1{1-e^{-x}}-\frac1x\le1.
\]

The first inequality follows by adding `a` to both positive denominators;
the right-hand expression is read by continuous extension at `x=0`, and the
final inequality is equivalent to `e^x-1>=x`. Hence the two resolvents differ
in operator norm by at most `delta`, and
`H_delta -> H` strongly resolvent. The self-adjoint
Trotter--Kato/Stone theorem then gives
`e^{-itH_delta} -> e^{-itH}` strongly, uniformly for `t` in compact real
intervals.

This direct argument also shows exactly which domain assertion is used: the
finite-charge Schwartz space must be an operator core for the original `H`,
not merely a form core. Round 20 states precisely this for the positive
polynomial Schrödinger direct sum and preserves it under bounded hopping.

The original Gibbs operator exists independently of the regulated Gibbs
limit. Let `H_diag=K_phi+V`. The form inequality `H>=H_diag`, compact
resolvent, and the min--max principle give

\[
\operatorname{Tr}e^{-\beta H}
\le \operatorname{Tr}e^{-\beta H_{\rm diag}}
\le \Theta(\beta)
   \left[2\sinh\!\left(\frac{\beta m}{2}\right)\right]^{-N}<\infty.
\]

Here each fixed-charge scalar stiffness has all `N` oscillator frequencies
at least `m`, and summing `e^{-beta D(n)}` over the neutral charges is bounded
by the Round-26 `Theta(beta)`. Thus
`rho_beta=e^{-beta H}/Tr(e^{-beta H})` is a faithful Gibbs/KMS state of the
original fixed-lattice Hamiltonian for every `beta>0`.

## Addendum: uniform Gibbs tail and trace-norm convergence

The original transfer structure supplies the extra compactness that strong
resolvent convergence alone lacks. Remove only the hopping exponential from
the middle factor and define

\[
Q_{\delta,0}
=e^{-\delta V/2}(e^{-\delta K_\phi}\otimes I)e^{-\delta V/2}.
\]

Because `K_phi` and `H_hop` act on different factors and
`0<=e^{-delta H_hop}<=I`, congruence by `e^{-delta V/2}` gives

\[
0<Q_\delta\le Q_{\delta,0}.
\]

The operators are injective but are not boundedly invertible, so applying
the logarithm informally would leave a domain gap. Regularize instead:
`Q_delta+epsilon I <= Q_delta,0+epsilon I`. Operator monotonicity of `log`
gives the reversed inequality for `-log`; monotone form convergence as
`epsilon downarrow 0` yields

\[
H_\delta\ge H_{\delta,0}
:=-\delta^{-1}\log Q_{\delta,0}
\]

as a closed quadratic-form inequality.

The logarithm on the right is exactly computable block by block. In a fixed
charge block, `D(n)` is scalar and the orthogonal eigenbasis of `A(n)` leaves
`K_phi` invariant. For one eigenfrequency `omega`, put

\[
\alpha=2\operatorname{arsinh}(\delta\omega/2).
\]

Direct comparison with the Mehler kernel gives

\[
e^{-\delta\omega^2q^2/4}e^{-\delta p^2/2}
e^{-\delta\omega^2q^2/4}
=\exp\!\left[-\frac\delta2
  (a_\delta p^2+b_\delta q^2)\right],
\]

where

\[
a_\delta=\frac{\alpha}{\sinh\alpha},
\qquad
b_\delta=\frac{\alpha\sinh\alpha}{\delta^2}.
\]

There is no missing scalar normalization: both kernels have prefactor
`(2 pi delta)^(-1/2)`, cross coefficient `1/delta`, and diagonal coefficient
`cosh(alpha)/(2 delta)`, with
`cosh(alpha)=1+delta^2 omega^2/2`. Thus the logarithm is the displayed
quadratic Hamiltonian, not that Hamiltonian plus an undetected constant.

For `y=delta omega/2`, the elementary inequality
`arsinh(y)>=y/sqrt(1+y^2)` gives

\[
a_\delta\ge\frac1{1+\delta^2\omega^2/4},
\qquad b_\delta\ge\omega^2.
\]

Fix `delta_0>0`. For `0<delta<=delta_0`, define

\[
c_n=\frac1{1+\delta_0^2\|A(n)\|/4},
\qquad
B=\bigoplus_n\left[
D(n)+c_nK_\phi+\frac{m^2}{2}|\phi|^2\right]
\]

on the neutral charge sum. Every eigenvalue of `A(n)` is at least `m^2`,
so the preceding mode inequalities prove the uniform form bound

\[
H_\delta\ge H_{\delta,0}\ge B.
\]

The comparison operator has compact resolvent and is Gibbs for every
positive inverse temperature. Its `n`-block consists of `N` oscillators of
frequency `m sqrt(c_n)`, hence

\[
\operatorname{Tr}e^{-\gamma B}
=\sum_{n\ {\rm neutral}} e^{-\gamma D(n)}
 \left[2\sinh\!\left(\frac{\gamma m\sqrt{c_n}}2\right)
 \right]^{-N}<\infty.
\]

Indeed, `A(n)` grows at most quadratically in the finite-dimensional charge
coordinates. Therefore `c_n^{-N/2}` has only polynomial growth, whereas the
Round-26 coercivity `D(n)>=a_c||n||^2` supplies Gaussian decay. This checks
the potentially dangerous large-charge blocks; no uniform positive lower
bound on `c_n` is being assumed.

The min--max principle now gives, uniformly for `delta<=delta_0`,

\[
\operatorname{Tr}e^{-\gamma H_\delta}
\le\operatorname{Tr}e^{-\gamma B}<\infty.
\]

Form order also gives a uniform first-moment bound in the fixed compact
operator `B`:

\[
\begin{aligned}
\operatorname{Tr}\!\left(B^{1/2}e^{-\beta H_\delta}B^{1/2}\right)
&\le \operatorname{Tr}(H_\delta e^{-\beta H_\delta})\\
&\le \frac{2}{e\beta}
       \operatorname{Tr}e^{-(\beta/2)H_\delta}\\
&\le \frac{2}{e\beta}
       \operatorname{Tr}e^{-(\beta/2)B}=:C_\beta.
\end{aligned}
\]

Let `P_R=1_[0,R](B)`, a finite-rank projection. Since
`I-P_R<=B/R`,

\[
\sup_{0<\delta\le\delta_0}
\operatorname{Tr}\!\left[(I-P_R)e^{-\beta H_\delta}\right]
\le C_\beta/R.
\]

Strong-resolvent convergence already proved above implies
`e^{-beta H_delta}->e^{-beta H}` strongly. On each fixed finite-dimensional
`P_R` block the traces therefore converge, while the displayed estimate
controls the regulated tails. The same bound applies to the limit because
the original `H` also satisfies `H>=B`. It follows that

\[
\operatorname{Tr}e^{-\beta H_\delta}
\longrightarrow \operatorname{Tr}e^{-\beta H}.
\]

For positive trace-class operators, strong convergence together with
convergence of traces implies trace-norm convergence. Hence

\[
\|e^{-\beta H_\delta}-e^{-\beta H}\|_1\to0,
\qquad
\|\rho_{\beta,\delta}-\rho_\beta\|_1\to0.
\]

This argument uses the actual Mehler normalization, the original coercive
`D(n)`, and the original positive hopping. It introduces no new physical
source or auxiliary model.

This theorem sharpens the Round-26 limitation. Round 26 itself did not
perform real-time reconstruction. Combining its transfer with the already
proved Round-20 form Hamiltonian does reconstruct mathematical real-time
dynamics in the `delta -> 0` limit at fixed finite spatial lattice. Nothing
in this argument selects that lattice theory from microscopic TFPT, takes a
spatial or volume continuum, establishes Lorentz symmetry, or identifies
the parameter as physical clock time.

For reference, Kato's closed-form product theorem is summarized as Theorem
18.2 in Barry Simon's review of Kato's work. Strong product convergence is
the general conclusion. Trace-norm convergence of Gibbs products requires
additional compactness or trace-ideal hypotheses; the source-specific
fixed-`B` tail above supplies them here. For general results see, for example, Ichinose and
Tamura, *Error bound in trace norm for Trotter--Kato product formula of
Gibbs semigroups*, Asymptotic Analysis 17 (1998), 239--266,
doi:10.3233/ASY-1998-310.

## Native response and decoherence functional

Choose, rather than derive, a finite sequence of projective decompositions
`{P_a^(k)}` from the diagonal local charge spectral projectors (or finite
coarse grainings of them) in the neutral site net. At prescribed times set

\[
P_a^{(k)}(t_k)=e^{it_kH_\delta}P_a^{(k)}e^{-it_kH_\delta},
\qquad
C_h=P_{a_n}^{(n)}(t_n)\cdots P_{a_1}^{(1)}(t_1),
\]

and

\[
D_\delta(h,h')
=\operatorname{Tr}\!\left(C_h\rho_{\beta,\delta}C_{h'}^*\right).
\]

This functional is Hermitian, normalized for exhaustive histories, additive
under coarse graining, and strongly positive. Indeed, for any finite family
of coefficients `z_h`,

\[
\sum_{h,h'}\overline{z_h}D_\delta(h,h')z_{h'}
=\operatorname{Tr}(B\rho_{\beta,\delta}B^*)\ge0,
\qquad B=\sum_h\overline{z_h}C_h.
\]

In particular `D_delta(h,h) >= 0`. Off-diagonal terms generally remain, so
these diagonal values form classical additive probabilities only after a
consistency/decoherence condition or a specified measuring instrument. The
word “positive” here means strong matrix positivity of the decoherence
functional, not that every history-word coefficient is a probability and
not that an instrument has been selected.

The Euclidean word weights of Round 26 must remain distinct from these
real-time history values. The former contain the original cocycle signs;
Round 26 proves their *absolute convergence* by an unsigned comparison and
explicitly retains a potentially enormous cancellation budget. No termwise
identification with positive event probabilities follows.

The addendum permits the stronger limiting response statement with each
regulator's own Gibbs state. Locally uniform strong convergence of the
unitary groups implies strong convergence of every finite product of bounded
inserted projectors, including the adjoints. Together with
`rho_{beta,delta}->rho_beta` in trace norm, trace-class duality gives

\[
\operatorname{Tr}(C_{h,\delta}\rho_{\beta,\delta}C_{h',\delta}^*)
\longrightarrow
\operatorname{Tr}(C_h\rho_\beta C_{h'}^*).
\]

This proves convergence of all fixed finite native bounded ordered
responses, including the full finite decoherence matrix. It still does not
select the projectors, observation times, instrument, or classical record.

## Remaining boundary

The combined sources support the following scoped chain:

`declared Round-20 finite-lattice source/Hamiltonian`
`-> positive Round-26 transfer`
`-> exact regulated Gibbs state and unitary time`
`-> strongly positive responses for prescribed native projectors`
`-> strong-resolvent recovery of the Round-20 Hamiltonian as delta -> 0`
`-> trace-norm Gibbs-state and finite bounded-response convergence`.

They do not derive:

- the Round-20 coupling, hopping channels, neutral-sector choice, Euclidean
  slicing, inverse temperature, projector schedule, apparatus, or record;
- a classical probability measure for interfering histories;
- a spatial/thermodynamic continuum, Lorentzian causality, a relativistic
  stress tensor, a physical-time selector, or full TFPT/Standard-Model/
  gravity reconstruction.

The proof of trace-norm convergence depends essentially on the newly
identified fixed-`B` spectral tail. A mere uniform bound on the partition
functions would still be insufficient, because finite-rank mass could escape
to higher modes.

## Atlas search note

Targeted Atlas searches for the exact combinations
`H_delta = -delta^-1 log Q_delta`, `strong resolvent`, `KMS`, `decoherence
functional`, and the symmetric Trotter identification did not return an
existing contract proving this chain. They returned Round 26, unrelated
marked-channel Lindblad work, and older finite GNS/OS scripts. Accordingly,
the fixed-regulator spectral-calculus proof and the Round-20/Round-26
Trotter identification above are a new synthesis of existing sources, not a
claim already present in either contract.

One Atlas retrieval call failed because `research_context` requires an
explicit record `id`; the attempted query-only call returned MCP error
`-32602: missing string argument id`. A subsequent `search_knowledge` plus
`get_dossier(contract:hopping-source-domain-round20)` correctly resolved the
owning proof as `source:e39426b76acca89849ef`. This was a caller-contract
mistake, not an empty or corrupt Atlas response.
