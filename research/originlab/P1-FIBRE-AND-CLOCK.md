# Marked fibre integration and the native quarter action

## Two supplied source levels

Set `coarse_graining` in an `original_exact_finite_fiber_sum` endpoint
bundle to compare two levels before the endpoint cap and selector.
The fine level is the existing `fiber_records` list. Supply:

- `time_id` and nonempty `configuration_map_binding`, `word_embedding_binding`
  and `reference_binding` fields identifying the common source operations.
- `coarse_fiber_records` in the same atom format, and `fine_to_coarse`,
  a map covering all fine atoms and all declared coarse atoms.
- `word_embedding`, an injective matrix from the coarse word space to
  the fine word space, with each entry in the existing exact pair encoding.
- `coarse_scalar_factor` and `coarse_reference_scalar` in the same pair
  encoding as the fine level. Both reference scalars must be supplied.

The consumer forms the weighted pushforward with that map and matrix.
It checks the independently supplied ratio of overall scalar and empty
reference at each level. Every coarse atom must match separately for
each of Z0, ZU, ZY and ZUY. On failure the result includes the first
coarse configuration, the ordered mark and the residual rank. It
emits no candidate source values from inconsistent input. In the
negative control two mixed defects cancel in the aggregate but fail
this per-configuration check.

The current native diagnostic has no original atom packet or two-level
weighted forms, so it returns `BLOCKED`. Agreement of supplied forms
checks the stated finite identity. Physical fibre provenance, the
OS transfer and the regulator limit require further evidence.

## Quarter lift on native currents

For anticommuting involutions R and ε, set C = Rε and
W = (I+C)/√2. W rotates the two polarizations. After a scalar
rephasing, its eigenvalues are 1 and −i. Conjugation on any such
single-cell operator space has no phase −1. The native affine
grade-one glue multiplicities are (60, 64, 60, 64) in phases
(1, i, −1, −i), including the eight Cartan directions. Thus this
specified single-cell conjugation cannot intertwine the native currents.

Two commuting copies have all four adjoint phases. Suppose their
first-cell polarization ε₁ preserves the identified native grade-one
space. Conjugation by ε₁ negates the square of the tensor clock,
because ε₁ anticommutes with its own C₁ and commutes with the
second-cell C₂. The square must then have trace zero on that space.
The native glue square has trace 120 − 128 = −8. The corresponding
native parity traces in affine grades 0 to 3 are 1, −8, 28 and −64.

This obstruction holds in every affine grade. Put q = exp(2πiτ),
and let θ₂, θ₃, θ₄ be the standard theta constants at τ. The
integer D₈ class contributes (θ₃⁸+θ₄⁸)/2 and the half-integer
spinor class contributes θ₂⁸/2 with opposite glue parity. The eight
neutral oscillators contribute ∏(1−qⁿ)⁻⁸. With
θ₂⁴=θ₃⁴−θ₄⁴ and the [Jacobi products](https://dlmf.nist.gov/20.5),
the resulting trace series without the vacuum energy factor is

```text
Σ trace(U² | grade n) qⁿ = ∏ (1 − q^(2j−1))^8, j ≥ 1.
```

Every term at grade n has sign (−1)ⁿ, because it uses only odd parts.
There is at least one term in each grade, including grade two from
two colours of the part one. Integer theta convolution independently
matched the product through grade 32. A degree-preserving balancing
extension needs additional directions in every grade.

Every parity-reversing linear map on the original grade-one space
has nullity at least eight. If only clock equivariance of the
embedding into a balanced enlarged space is required, the minimum
dimension for the spectral counts is 256. A direct sum of 32 tensor
copies and 32 with a common sign shift gives 64 states in each phase.
The embedding drops four states in each even phase. This gives a
spectral cover. The original currents, adjoints and orientation on
the compressed space still require an intertwiner. The original
collar may instead have a non-invertible or sector-changing
polarization, to be checked on its actual domain.

The local tensor clock exchanges the two commuting Bell stabilizers
R⊗R and ε⊗ε. It preserves their common positive vacuum line and
equates the two single excitation energies. The general simultaneous
diagonal Hamiltonian with that vacuum and clock has energies
(0, a, a, b). Only the additive stabilizer ansatz gives b=2a.
Thus the double excitation provides a held-out dynamical probe.

Run `--campaign colour_ward --campaign r113 --campaign r113_synthesis --campaign f1`
from the standard OriginLab entry point. The clock result is under
`marked_source.clock_lift` in the colour report, and the supplied
fibre comparison or its missing inputs are under `fibre_pushforward`.
