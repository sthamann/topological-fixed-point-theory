# Fixed Peierls phases do not remove the Euclidean charged-pair obstruction

15 October 2026. Addendum to `UR.NATIVE.CHARGED_CONJUGATION.01`.
**Exact conditional antiunitary statement; overall physical verdict PARTIAL.**
No source or phase is selected from P1/P2. The original proof/checker and
source manifest remain unchanged.

## 1. The tempting repair and its exact remaining symmetry

Keep the original real bilinear-cocycle onsite fields and bond unitaries
from `local-charge-transport-round20`. Replace only each fixed hopping term
by its Peierls-phase form

    h_(b,theta)=J_b(2I-exp(i theta_b)W_b-exp(-i theta_b)W_b^*),
    H_theta=H_diag+sum_b h_(b,theta),  J_b>=0.

The theta_b may vary by oriented link/channel and may have nonzero flux.
They are fixed, not time-dependent. The diagonal/scalar source Hamiltonian
is the same original real C-even operator from the parent proof.

Let K be complex conjugation in the original integer-charge and scalar
position representation. The cocycle signs are real, so K T_p K=T_p and
K W_b K=W_b. The scalar momentum reverses under K, but its original square
and the real potential are unchanged. The linear charge conjugation C is
real in that basis, hence commutes with K. Define the antiunitary

    A=C K,  A²=I.

The two separate transformations obey

    K H_theta K=H_-theta,
    C H_-theta C=H_theta,
    A H_theta A^-1=H_theta,
    A T_p A^-1=T_-p.

For the hopping this is termwise: antiunitarity conjugates the coefficient
exp(i theta_b), while C sends W_b to W_b^*. Their two effects cancel.
No gauge transformation, gauge-trivial flux or spatial reflection is
assumed. C and K individually can genuinely fail at a quarter phase; A
still survives. The checker includes those failures as negative controls.

The onsite adjoint identity remains
`T_p^*=(-1)^e(p)T_-p`. Its actual sign is retained in the regression checks;
antiunitary conjugation also preserves adjoints.

## 2. Positive operator and invariant Gibbs starts

Because W_b is the original unitary,

    2I-exp(i theta_b)W_b-exp(-i theta_b)W_b^*
      =(I-exp(i theta_b)W_b)^*(I-exp(i theta_b)W_b)>=0.

The hopping has norm at most `4 sum_b J_b`; for the original uniform J this
is the original `96 J N` bound. It retains the unchanged diagonal/scalar
domain. In particular H_theta is the same type
of nonnegative selfadjoint source operator; its A covariance extends
from the core to its spectral calculus. The earlier coercivity and
trace-class finite-volume Gibbs construction remain applicable.

For clarity the original total-zero-sector Gibbs argument extends to the
full Round20 field space as follows. At fixed finite N let
`c=lambda_min(G)/(2N)>0`. Then

    D(n)>=c sum_x |n_x|²,   A(n)>=m² I,
    Tr exp(-beta H_diag)
      <= [2 sinh(beta m/2)]^(-N)
         [sum_(k in Z) exp(-beta c k²)]^(8N) < infinity.

The first bound is the actual quadratic coercivity; the second bounds the
N oscillator frequencies below by m. No neutrality constraint is needed
for this full integer-profile theta sum. The direct-sum diagonal operator
has compact resolvent. Its bounded positive hopping perturbation preserves
compact resolvent, and eigenvalue min-max comparison with H_diag bounds
`Tr exp(-beta H_theta)` by the displayed finite trace. This supplies the
full-space Gibbs statement; restriction to C0 also remains trace class.

An antiunitary conjugates a real spectral function without changing it.
Thus A exp(-beta H_theta) A^-1=exp(-beta H_theta). Normalized full Gibbs
densities are A invariant. Since A Q_tot A^-1=-Q_tot, its total-zero
projector P0 is invariant too, and

    rho_(theta,beta,0)=P0 exp(-beta H_theta)P0/Z0

is invariant. This is an A-invariance statement: a generic quarter-phase
Gibbs state need not remain C invariant by itself.

## 3. All-time Euclidean opposite-charge returns remain equal

For any A-invariant density rho define the same response as in the parent
contract,

    G_(rho,p)^theta(t)=Tr[rho T_p^* exp(-tH_theta) T_p], t>=0.

The operator `T_p^* exp(-tH_theta) T_p` is positive, so this finite trace is
real and nonnegative (strictly positive for nonzero rho). Antiunitary
conjugation of a trace-class operator conjugates its trace. Therefore

    G_(rho,-p)^theta(t)=conjugate(G_(rho,p)^theta(t))
                      =G_(rho,p)^theta(t).

This holds for every t, every original integer p and arbitrary fixed link
phases/flux, at the stated A-invariant starts. It is not inferred from a
few sampled times or a zero-moment comparison. Charged onsite insertion
and removal still use the full field space and exterior charge-sector
propagation when rho is supported in C0.

It follows that simply putting a fixed quarter Peierls phase into the
original W hopping cannot reproduce the asymmetric Euclidean pair ratio
of a single fixed chiral delta at those starts. This excludes that proposed
repair in its actual original source class. It does not remove the
positive local T fields, balanced W dynamics, or their phase/flux-sensitive
responses in other channels.

## 4. Precisely what changes the premise

The result does not assert that all real-time words, all prepared states,
all response instruments or every Hamiltonian extension are identical.
The argument used the positive Euclidean return to remove the complex
conjugation introduced by A. It must not be transferred to an arbitrary
complex multitime amplitude without checking its transformation.

Time-dependent phases, an A-breaking start or boundary prescription, an
additional A-odd real source term, or a new field dictionary need their
own source derivation. For example the uniform chemical term eta.Q_tot
with real eta is A odd. The parent proof already shows that this term can
alter the charged-corner ratio while vanishing inside C0. That is a
conditional discriminator, not a selected physical repair.

No new candidate source world, state mixture, sector-switch process,
physical event rate or P1/P2 completion is introduced here. The next
original-source question is whether its actual quarter seam supplies the
required A-asymmetric marked operator/state/time data, rather than merely
a complex number multiplying otherwise symmetric bond hopping.

## 5. Exact validation and authority

Run `python3 experiments/theory-contracts/native-charged-conjugation-20261015/peierls_checker.py`.
The saved run passed **16/16 exact groups**, zero failures. It reads the
actual original Round20 Gram, cocycle, onsite translation and bond-transfer
functions without executing its historical prerequisite suite. Arbitrary
symbolic c+i s verifies the coefficient/conjugation identity; zero, +/-
quarter and half phases provide finite controls. Three unbounded-charge
profiles are regressions of the monomial identities, not a cutoff Hilbert
space or a physical three-site model.

The all-charge/all-time theorem is the analytic argument above. No Gibbs
matrix diagonalization, heat-kernel simulation, finite moment inference,
GPU experiment or physical-gate promotion is claimed. `peierls_results.json`
records the original source pin and checker/addendum hashes.
