# Native spatial return at every Fourier depth

The original full96 diagnostic proposal `[2,1,1]` has constant Higgs
component `(0,0,1,0)` and first harmonic `(13/100) cos(x)` in component
zero. All other harmonics vanish. The original Fourier assembler uses
`h[n,m] = (n+alpha) Gamma δ[n,m] − i Gamma Phi[n−m]`. Keep modes `−1,0`
and extend that same rule to all integers. The constant Majorana term
adds to the diagonal blocks. Both off-diagonal blocks come from the
first harmonic. The original spatial coordinate has no assigned OS time.

## Exact rank and first missing return

Write `B = h[n,n+1]` and `d = 96`. The first harmonic is a real unit
quaternion direction multiplied by `13/200` on each Fourier edge.
Its left weak action is unitary. The chiral blocks of the native lift
therefore repeat the singular values of the four Yukawa matrices and
of their conjugates, scaled by `13/200`. Multiplication by `−i Gamma`
preserves these singular values. Reading the source decimals as exact
rationals gives four rank three family matrices. Thus `B` has rank `d`.
In particular, this argument does not depend on a floating point
singular value cutoff.

Let `P` embed modes `−1,0`. The span of `[P,hP,...,h^kP]` is exactly the
sum of modes `−1−k,...,k`, so its rank is `192(k+1)`. To see this, the
two new end modes of `h^k P` can only be reached by moving outwards at
every step. Their coefficients from the respective end of `P` are the
invertible products `B^k` and `(B*)^k`. Appending these blocks to the
previous span proves the induction. The block Hankel matrix with entries
`P* h^(i+j) P` is the Gram matrix of these columns and has the same rank.
No fixed finite matrix model can reproduce all the spatial return moments
on this input. In particular a finite auxiliary matrix self-energy cannot
reproduce its entire compressed spectral response.

Cutting after `m` extra modes on each side reproduces all return powers
through `2m+1`. At `2m+2`, the always-outward and then always-inward
paths are the only omitted returns. The difference is a positive
definite matrix on the two retained modes, with blocks
`(B*)^(m+1) B^(m+1)` and `B^(m+1) (B*)^(m+1)`. If a positive average
of these spatial moments contains this configuration with nonzero
atomic mass, its Hankel forms retain the unbounded rank. This last
statement concerns that spatial moment average and assumes the moments
exist. The actual R113 source support and observable dictionary are open.

## Computing the return

The finite Schur recursion eliminates both Fourier tails using only
`96×96` solves and then inverts the retained `192×192` block. If `A_n`
is the diagonal block, each inward update is
`Σ_n = C_n (zI − A_(n±1) − Σ_(n±1))^(-1) C_n*`, with `C_n` equal to
`B` or `B*` on the corresponding side. The Fourier diagonal is
self-adjoint on its natural domain. The bounded Hermitian remainder
preserves self-adjointness, and the finitely supported vectors remain a
core. The cut resolvents converge strongly off the real axis. For a
finitely supported core vector `u`, once the cut contains `u` and its
neighbours, the cut resolvent on `(z−h)u` equals `u`. Density of these
vectors and the uniform `1/|Im z|` bound give the convergence. On the
fixed retained modes this also gives norm convergence.

The rational family Gram traces give `||B|| < 5/16` by bounding the
operator norm of each Yukawa block by its Frobenius norm. At `z=2i`, the
Neumann series in the off-diagonal hopping around the diagonal
resolvent has ratio at most `r=5/16`. The infinite and cut series
agree on the retained pair for all terms with fewer than `2(m+1)`
hops. Bounding both remaining geometric tails gives the explicit error
upper bound `r^(2m+2)/(1−r)`. This bounds the exact truncated resolvent
before floating point evaluation. The implementation also reproduces
the direct four-mode inversion through the original assembler.

The same elimination gives the relative spatial determinant by multiplying
the Schur pivots and dividing by the determinants of the diagonal reference
blocks. Its infinite-cut limit is the regularized determinant
`det₂(I − (z−D)^(-1) E)`, where `D` is the block diagonal and `E` the hopping.
Indeed the diagonal resolvent is Hilbert–Schmidt since the block norms
are `O(1/|n|)` at large `|n|`. The hopping is bounded, and the finite
compressions converge in Hilbert–Schmidt norm. The finite trace of
`(z−D)^(-1) E` vanishes because the diagonal blocks are zero. The usual
finite determinant ratios therefore equal their `det₂` counterparts.
Continuity of `det₂` gives the limit. This determinant uses the spatial
zero-hopping reference. The absolute empty R113 reference is a separate
input.

Run the native `source_bridge` campaign to replay the determinant,
neighbourhood, cutoff and resolvent checks. This computes a spatial
part of the first required source operation. The boson and geometry
basis integration, selector, absolute reference, physical time and
four weighted R113 forms remain open.
