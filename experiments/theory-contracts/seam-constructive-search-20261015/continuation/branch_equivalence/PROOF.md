# Exact branch equivalence and a constructive cyclic lift

## Result and scope

There are two different questions. In the four-key search presentation, the
uniform characters `chi=1,3` and the permutation signs `kappa=0,2` give four
normal-form solutions. In the actual operator envelope containing the already
given source parity `Pi=exp(i*pi*Q)`, the two characters are related by the
complex-linear, source-fixing change `X_plus = Pi X_minus`. It preserves the
supplied charge, time, quarter character, canonical state restriction and
`U_s R_c=W4`. Thus there are **two marked-P lift classes**, rather than four
classes modulo every source-compatible choice of the unmarked defect generator.
If the physical operation X itself is fixed, its two different word responses
remain distinguishable. Changing its definition is then not an equivalence of
that additional operational marking.

The original unit winding does not yet remove the remaining kappa sign.
A completely explicit, unitary, charge- and time-compatible extension has
`kappa=2`, `C^8=I`, `W1^4=W4` and `W1^8=Pi E_q`. It closes using only the
canonical defect grade D. The inference `kappa=0` is correct **if one additionally
requires the scalar field relation `C X C^-1=nu X`**. That scalar relation does
not follow just from the original integer-sector T1 action. In the constructive
counterexample the conjugate is a fixed linear combination of X and XD.

These statements concern the declared global crossed-word and charge-space
construction. They neither select its physical source from P1/P2 nor construct
a local field net. No new Hamiltonian, state weight, beta or fitted rate is used.
The analytic identities below hold for every integer charge. `checks.py` checks
their exact mu8 arithmetic on a small signed-charge set; it is not a Lean proof
or a finite replacement for those analytic arguments.

## 1. Fixed source and a unitary all-charge realization

Write

```
B(n,m) = sum_(i>j) n_i*m_j,  eps(n,m)=(-1)^B(n,m),
q=(1,...,1), s=q/2, a=(1,1,1,1,0,0,0,0),
C(n)=(n7,n0,...,n6), P=C^4,
gamma(n)=(-1)^(n7*(n0+...+n6)), eta(n)=(-1)^(A(n)*Bblock(n)).
```

Use the algebraic orthonormal charge basis `|n,c>`, n in Z^8 and c=0,1,
with charge p=n+c*s. This realizes the finite-word algebra on
`l2(Z^8 x Z2)`; there is no finite cutoff or boundary wraparound. For
`chi=1,3`, `zeta in Z4` and `nu in mu8`, define

```
E_m |n,c> = eps(m,n) |m+n,c>,
X |n,c> = i^(chi*Q(n)+zeta*c) eps(n,q)^c |n+c*q,1-c>,
C_nu |n,c> = gamma(n) nu^c |Cn,c>,
P_kappa |n,c> = eta(n) i^(kappa*c) |Pn,c>,
Z_u |n,c> = (-1)^(u.n) i^(c*Q(u)) |n,c>,
D |n,c> = (-1)^c |n,c>,
Pi |n,c> = (-1)^Q(n) |n,c>.
```

All are unitary monomial operators: their charge maps are bijections and their
coefficients have unit modulus. Their inverses therefore define the star
operation exactly. Bilinearity gives `E_n E_m=eps(n,m) E_(n+m)`.
The character condition `2*chi=2 mod4` gives

```
2*chi*Q(n)+2*B(n,q) = 2*B(q,n) mod4,
```

and `chi*Q(q)=0 mod4`. Consequently `X^2=i^zeta E_q` and
`X E_n X^-1=i^(chi*Q(n)) E_n` for all n.

The factor `eps(n,q)^c` in the displayed X action is essential. Omitting it
does not represent the declared normal form `E_n X^c`; already n=e1 makes
`X^2=i^zeta E_q` fail. One can move cocycle factors to other formulas by a
consistent basis change, but cannot delete this one while retaining the
displayed left-regular E action unchanged.

Append the original eight-channel current oscillator Fock space if required.
C and P permute its channels, E and X perform the stated charge transport,
and total oscillator level N is unchanged. Every zero-mode identity below
therefore holds on every finite current descendant in that tensor extension.

## 2. The canonical grade closes both cyclic lifts

The source cocycle identity for C is

```
gamma(n+m)/(gamma(n)*gamma(m)) = eps(Cn,Cm)/eps(n,m).
```

It follows by expanding the two terms crossing the cyclic cut. Thus C_nu
implements exactly the original cyclic source action on every E_n, for every
nu. In particular `gamma(q)=-1`.

For the fourth power, the gamma exponent is

```
sum_(j=4..7) n_j*(Q(n)-n_j)
 = Bblock(n)*Q(n)-sum_(j=4..7) n_j^2
 = A(n)*Bblock(n) mod2.
```

For the eighth power its exponent is `Q(n)^2-sum_j n_j^2`, which is even.
Therefore

```
C_nu^4 |n,c> = eta(n) nu^(4c) |Pn,c>,
C_nu^8 = I,
C_nu^4 = P_kappa  iff  nu^4=i^kappa.
```

Both kappa signs are realized. For kappa=0 take nu=i. For kappa=2 take
`nu=exp(i*pi/4)`.

The exact conjugation of X is

```
C_nu X C_nu^-1 = X*(nu*P_even - nu^-1*P_odd)
 = ((nu-nu^-1)/2)*X + ((nu+nu^-1)/2)*X*D,
P_even=(I+D)/2, P_odd=(I-D)/2.
```

To verify it, compare `C_nu X` and `X C_nu` on each input grade. The ratios
are nu in grade zero and `-nu^-1` in grade one. The minus sign is precisely
`gamma(q)=-1`, combined with the source cocycle identity. It holds for both
chi orientations and every zeta. The conjugate is scalar exactly when
`nu=-nu^-1`, equivalently `nu^2=-1`. This implies `nu^4=1` and kappa=0.

For the concrete kappa=2 choice the formula is particularly small:

```
C X C^-1 = (i*X + X*D)/sqrt(2).
```

This also reconstructs D rather than postulating an unrelated ancilla.
With `R=X^-1 C X C^-1`,

```
R=a*I+b*D,  a=(nu-nu^-1)/2, b=(nu+nu^-1)/2.
```

Whenever b is nonzero, closure under the existing C/X word forces
`D=(R-a*I)/b`. In the kappa=2 example this is `D=sqrt(2)*R-i*I`.
The resulting D has eigenvalues +1 and -1, is self-adjoint, commutes with
the old integer E_n, and anticommutes with X. It is exactly the already
present integer/half-integer transport grade. It adds no independent state
coordinate to the explicit charge representation. Calling it a physically
accessible measurement or a derived autonomous memory law would require
additional operational evidence.

Hence the scalar-C calculation correctly excludes kappa=2 in the **smaller
scalar-eigenfield closure**. The counterconstruction realizes the same original
integer C action in the canonical grade-closed operator envelope.

This closure uses no new state or time parameter. D is already the canonical
dual-Z2 grade implementer, with `D E_n=E_n D`, `D X=-X D`, `D^*=D`, `D^2=I`.
The span `A + A D` is closed under multiplication and star. Equivalently,
with a base lift C_i and `R_t=P_even+t*P_odd`,

```
C_nu=R_(nu/i) C_i,   P_(kappa=2)=D*P_(kappa=0).
```

The relative sector lift phase nu is explicit; its selection is not derived.
In the stipulated half-charge torus representation D is even already
`Z_e0^2`, so R_t is a polynomial in an operator used in W1.

There is an important global-group qualification: `Z_e0^2=I` holds on the
old integer sector, whereas `Z_e0^2=D` on this half-charge extension. An
extension of the original torus Lie action or its covering action is not an
injective embedding of an operator algebra already quotiented by the global
identity `Z_e0^2=I`. If that latter identity is demanded globally, none of these
half-charge lifts with `Z_e0 X Z_e0^-1=iX` is permitted. The construction does
not conceal this change of representation sector/global torus cover.

## 3. Actual W1, W4, carry, time and quarter character

Define the original source-normalized words

```
W1=-Z_e0 E_e0 C_nu,   W4=Z_a E_a P_kappa.
```

Their charge maps are respectively `p -> Cp+e0` and `p -> Pp+a`.
Writing Q=Q(n), the W1 coefficient is

```
(-1)^(n7*(Q-n7+1)) * i^c * nu^c.
```

Its fourth-power sign reduces to `(-1)^(A*Bblock+n5+n7)`; this is also the
sign in `Z_a E_a P_kappa`. Since `i^(4c)=1` and `nu^4=i^kappa`,
`W1^4=W4` as complete operators, not just as charge maps. Similarly
`W1^8=Pi E_q`. Thus the cyclic carry is not reduced modulo four and the
kappa=2 construction retains the eighth-power source identity as well.

For `K_alpha(p,N)=N+|p|^2/2+delta_alpha.p`,
`delta_alpha[j]=(j+alpha)/8-1/2`, direct expansion gives

```
K_alpha(Cp+e0,N)-K_alpha(p,N) = (Q(p)+alpha)/8,
K_alpha(Pp+a,N)-K_alpha(p,N) = Q(p)/2+3/4+alpha/2.
```

The first equality uses `delta.(Cp-p)=Q(p)/8-p7`; its p7 cancels the
quadratic shift. The second equality follows by summing four consecutive
integer charge increments. The coefficients nu do not change charge or
energy. For alpha=1/4,3/4 the complete W4 character is exactly

```
exp(2*pi*i*(3/4+alpha/2)) * (-1)^Q(p).
```

Neither the quarter phase nor forward-time convention has been removed.

## 4. Character equivalence in the supplied parity envelope

The actual source already supplies `Pi=exp(i*pi*Q)`, with

```
Pi^2=I,   Pi E_n Pi^-1=(-1)^Q(n) E_n,
[Pi,X]=[Pi,P]=[Pi,C]=[Pi,K]=0.
```

Therefore `Y=Pi X_minus` has

```
Y E_n Y^-1=i^Q(n) E_n,  Y^2=i^zeta E_q,
P Y P^-1=i^kappa Y.
```

It is exactly X_plus in the explicit realization. This defines an invertible,
complex-linear, source-fixing star-isomorphism, fixes K and the full clock
operators, and preserves E0 because it preserves the defect grade. In
`U_s=i^-zeta X`, `R_c=Z_a E_(-b) X P`, both factors acquire Pi, so their
product remains W4. The old scalar-section gauge `X->tX` could not do this;
the existing source parity can.

State preservation is explicit. For `a=a0+a1*X_plus`,

```
F(a)=a0+a1*Pi*X_minus,
E0_minus(F(a))=a0=E0_plus(a),
omega_minus(F(a))=omega0(a0)=omega_plus(a).
```

Thus `U_F[a]_plus=[F(a)]_minus` is an isometry of the GNS pre-Hilbert
spaces: `omega_minus(F(a)^*F(b))=omega_plus(a^*b)`. F is invertible, so it
extends to a unitary, sends Omega to Omega, and intertwines every fixed
source operator. It intertwines X_plus with the **dressed** Pi X_minus,
not with a separately operationally marked X_minus. This proves the scoped
global crossed/zero-mode equivalence with canonical state, not a local-field
intertwiner, an OPE identification, or a full physical equivalence. Those
stronger maps have not been constructed here.

By contrast, dressing with a charged source word E_v cannot achieve an
equivalence with the same full charge and time marking: its transport becomes
`s+v`, its square has charge `q+2v`, and the kinetic difference contains the
nonconstant term `v.p`. Keeping all those identities forces v=0.

For the declared single-defect charge module, an isomorphism fixing the
source E_n, full charge labels, K and the named P cannot flip kappa. In a
charge-diagonal multiplier presentation its image has form `f(n,c) X`.
Intertwining every E_m forces f to be constant in n when chi is unchanged,
or a constant times `(-1)^Q(n)` when chi is flipped. P preserves both grade
and Q, so these multipliers commute with P and leave its X sign unchanged.
The same conclusion holds in the canonical D envelope. This is a statement
about this marked, one-copy crossed module, not a classification of arbitrary
additional multiplicity spaces or unrestricted C-star completions.

Relabelling `P -> D P` changes kappa but changes the named extended source
implementer. It is not a marked-P isomorphism. Once that lift marking is
forgotten, both kappa choices live constructively in the same D-closed
operator algebra. The old integer source alone does not distinguish them.

Complex conjugation also flips chi, sends nu to its conjugate and zeta to
minus zeta, while preserving the real kappa signs. It is antiunitary: it
conjugates quarter phases and reverses `exp(-itK)` to `exp(+itK)`. It is
therefore not an equivalence fixing the oriented complex clock. The Pi
dressing above achieves the character equivalence without this reversal.

## 5. Minimal return words and the canonical state

The shortest mixed closed words exposing the two presentation labels are

```
L_j = X E_ej X^-1 E_ej^-1 = i^chi I,
M   = P X P^-1 X^-1      = i^kappa I.
```

They each use four factors; the two-path comparisons `X E_ej` versus
`E_ej X`, and `P X` versus `X P`, are the corresponding open-path version.
For any normalized `omega_tilde=omega0 o E0`,

```
omega_tilde(L_j)=+i or -i,
omega_tilde(M)=+1 or -1.
```

The canonical conditional expectation does not erase completed-loop phases.
It erases a single odd-grade insertion, `omega_tilde(aX)=0`. Positivity is
independent of the choice of omega0: for `y=a0+a1X`,

```
E0(y^*y)=a0^*a0 + phi^-1(a1^*a1) >=0.
```

In the canonical GNS representation dual invariance implements D by
`D[a]=[theta(a)]`, fixes Omega, and gives a parameter-free positive state on
the D envelope: `(a+bD)Omega=(a+b)Omega`. Thus adding the canonical grade
carrier requires no new fitted state. Whether that global carrier is a
physically available measurement is a separate source/instrument question.

The L_j response is not invariant under the allowed redefinition `X->Pi X`:
it flips sign because Pi anticommutes with the odd source E_ej. It therefore
distinguishes two **named operations**, not the two isomorphic unmarked
defect-generator presentations. The source's fixed complex quarter clock
can calibrate its imaginary quadrature if a concrete X operation is supplied;
neither that calibration nor omega0 selects which X is physical.

M remains a genuine marked-P criterion. Its two exact real responses survive
for any supplied omega0. To use this as a source decision one must supply
the actual action of the named P on the half-sector operation X, or derive
the scalar cyclic-field covariance from the original local field structure.
The original integer winding, its charge carry, its Ward identity and its
quarter character alone do not supply this missing mixed response.

There is a concrete conditional route to that missing criterion. If the
original source supplies an injective covariant state-field map Y, a fixed
vacuum Omega, and a one-dimensional marked half-charge primary x=X Omega,
then `C Omega=Omega` and `Cx=nu*x` imply
`C Y(x) C^-1=Y(Cx)=nu Y(x)`. To apply this to the global unitary X one must
also prove its identification with the appropriate charge part of that
specific field, including the oscillator factors and domain. Under those
additional identifications scalar C covariance and hence kappa=0 follow.
The grade-closed counterconstruction does not supply such a local state-field
identification: XD and X have the same vacuum-created vector but different
actions on the half sector. An injective single-field identification cannot
silently treat both as the same field. The original source explicitly calls
its W_n a global implementer rather than an already localized vertex field;
thus this field bridge cannot be imported solely from its charge map.

## Sources and reproducibility

Original source: `experiments/theory-contracts/seam-joint-ward-state-20261014/PROOF.txt`,
lines 752-861 and 1003-1087, SHA256
`7c0cd92dfc5eca88260baae0af9edad66f958a27ce0d15d7965414b3971f8da0`.
These give the actual T1/W4 action, time, parity, ordered CAR cocycle, the
stipulated extension and the explicitly unselected original P-lift signs.
The cyclic gamma cochain and the D-closed counterconstruction are derived
here; they are not quoted from those source lines as pre-existing proofs.

Run `python3 checks.py` in this directory. Results are written to `results.json`.
No existing engine, objective, prior result, ledger or paper is modified.
