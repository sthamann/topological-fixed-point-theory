# Source bridge for C11

This note supplies the concrete mathematical bridge between the fixed source
relations and the abstract `CrossedData A` used by `StateExpectation.lean`.
It is deliberately explicit about what Lean checks and what is transferred by
the standard finite-linearity argument.

## 1. The concrete base star algebra

Let

```text
A = C_c[Z^8]epsilon
```

be the complex vector space of finitely supported functions on `Z^8`, written
as finite sums `a = sum_n a_n E_n`.  On basis elements define

```text
E_n E_m = epsilon(n,m) E_(n+m),
epsilon(n,m) = (-1)^(sum_(i>j) n_i m_j).
```

Extend the product bilinearly.  Every coefficient in a product is a finite
sum because both inputs have finite support.  The ordered exponent is exactly
`epsExp` in the Lean file.  Lean proves its cocycle identity for all
`n,m,l : Z^8`, hence basis multiplication is associative; bilinearity then
gives associativity on every finite sum.

Define

```text
E_n^* = epsilon(-n,n)^(-1) E_(-n)
```

and extend conjugate-linearly.  Since `epsilon` takes values in `{+1,-1}` and
`epsilon(-n,n)=epsilon(n,-n)`, Lean's two `basis_unitary` theorems prove
`E_n^* E_n = E_n E_n^* = 1` for every integer charge.  The cocycle identity
then gives `(ab)^*=b^*a^*` on basis products, and finite conjugate-linearity
extends it to all of `A`.

## 2. The candidate automorphism and the square element

For any accepted candidate set

```text
chi(n) = i^L(n),       L(n) = sum_j chi_exp[j] n_j mod 4,
u = i^zeta E_q,        q = (1,1,1,1,1,1,1,1).
```

Define `phi(E_n)=chi(n)E_n` and extend complex-linearly.  Its inverse is the
multiplier `chi(n)^(-1)`, so `phi` is bijective.  The following table maps each
`CrossedData` field to an already certified source identity or its immediate
basis consequence.

| `CrossedData` field | Concrete verification |
|---|---|
| `phi` multiplicative | Production theorem `chi_add`; on `E_nE_m` both sides have coefficient `chi(n+m)=chi(n)chi(m)`. Extend bilinearly. |
| `phi_star` | `chi(-n)=chi(n)^(-1)=conj(chi(n))`, so both sides on `E_n` equal `conj(chi(n)) E_n^*`. Extend conjugate-linearly. |
| `uInv` | Set `uInv=u^*=i^(-zeta)E_q^*`. The Lean basis-unitarity identities give `uInv*u=u*uInv=1`; the scalar is a unit-modulus mu4 phase. |
| `phi_sq` | Production theorem `conjugation_square_all_charges` proves `chi(n)^2=epsilon(q,n)/epsilon(n,q)` for every `n`. This is exactly `phi^2(E_n)=u E_n uInv`; the scalar `i^zeta` cancels in the adjoint action. Extend linearly. |
| `phi_u` | Production theorem `chi_q` gives `chi(q)=1`, hence `phi(u)=u`. |
| `u_unitary_left/right` | Follows from the preceding `uInv=u^*` calculation and Lean's two basis-unitarity theorems. |

No candidate branch is fixed in this argument.  It uses only the identities
certified for both uniform orientations `chi_exp=1` and `chi_exp=3`.  The
arbitrary `zeta` is a unitary scalar gauge, and `kappa=0,2` does not enter C11.
Thus the bridge applies to all four orientation-by-P-lift classes and every
accepted zeta representative.

## 3. The crossed algebra and expectation

With these concrete data, the algebraic crossed extension is the vector space

```text
B_alg = A direct_sum A X,
```

with the multiplication and star defined in `StateExpectation.lean`.  The
field-by-field verification above provides an actual `CrossedData A`
mathematically.  The general Lean theorems therefore apply and give

```text
E0(y^* y) = a0^* a0 + phi^(-1)(a1^* a1)
```

for every algebraic finite word `y=a0+a1 X`.  Consequently `omega0 o E0` is
positive whenever `omega0` is positive on `A`.

## 4. Exact formal boundary

Lean checks the universal coefficient identities, basis associativity and
unitarity, the full abstract crossed multiplication/star laws, the SOS formula,
positivity, and dual-invariant uniqueness.  This continuation does **not**
construct `C_c[Z^8]epsilon` as a Mathlib `Finsupp` ring and instantiate
`CrossedData` with it inside the kernel.  The passage from basis identities to
finite sums is the ordinary finite bilinear/conjugate-linear extension given
above; it is a mathematical proof step, not a hidden Lean axiom.

It also does not prove a C-star completion.  Unconditional `K_alpha` state
invariance is separate: the quadratic charge flow is not automatically an
endomorphism of the bare span of the `E_n`.  A grading-preserving regular
representation with the required Cartan/number operators and an invariant
base vacuum must be supplied before the conditional time-invariance theorem
can be applied.

