# C11 algebraic state expectation

This continuation is isolated from the pinned production certificate.

`StateExpectation.lean` defines the concrete two-component normal form
`a0 + a1 X`, its crossed multiplication, its star operation, and the
degree-zero projection `E0`.  From `X a X^-1 = phi(a)`, `X^2 = u`, and the
corresponding coherence relations it proves

```text
E0(y* y) = a0* a0 + phi^-1(a1* a1).
```

Lean also proves that the displayed multiplication is associative and that the
displayed star is involutive and reverses multiplication.  It then proves that
`omega0 o E0` is positive on these finite words whenever `omega0` is positive
on the base algebra.

Before the abstract crossed construction, the file also instantiates the
actual ordered eight-channel cocycle on basis monomials `i^p E_n`.  Lean proves
its cocycle associativity and both unitary inverse identities for every
`n : Z^8`.

The concrete base algebra of finitely supported sums and its candidate
automorphism are connected to the general theorem in
[`SOURCE_BRIDGE.md`](SOURCE_BRIDGE.md).  That note checks every `CrossedData`
field from the production identities and extends the basis identities by
finite bilinearity.  The repository does not yet contain a kernel-level
Mathlib `Finsupp` instantiation of `CrossedData`; the receipt therefore reports
this as an analytical finite-linearity bridge rather than an end-to-end Lean
instantiation.

This result is algebraic.  It does not construct a C-star completion, select
`omega0`, prove Gibbs/KMS properties, or derive a physical state from P1/P2.
Dual-Z2 invariance is not needed for positivity.  A separate uniqueness theorem
shows that it is exactly what forces an additive extension to vanish on the odd
coefficient, hence to equal `omega0 o E0`.

Time invariance is also separate.  The file proves a conditional theorem: if a
time evolution preserves the grading and its even restriction leaves `omega0`
invariant, then the induced state is invariant.  C11 positivity alone does not
establish either premise for the physical `K_alpha` evolution.  In particular,
the span of twisted lattice monomials is not automatically stable under the
charge-dependent quadratic `K_alpha` flow; a representation containing the
Cartan/number operators and an invariant vacuum must be supplied separately.

