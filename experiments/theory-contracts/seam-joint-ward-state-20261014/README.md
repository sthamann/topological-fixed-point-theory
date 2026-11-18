# Joint seam Ward, chirality, state and time contract

This theory contract asks whether the same folded source completion can keep
the source charge/clock class, determinant-only packet response, original fold
time law, and a unique vacuum.  In the declared class of even unimodular
lattice extensions containing the original `D8`, the answer is no.  The full
lattice time calculation selects `Lc`.  In `Ls`, however, `lambda=q/2` is
already a norm-two root.  A genuine E8 isometry keeps it a root and its
projection pair `(1,1)` forces the mixed `(10,6)` class.  Thus `q=2 lambda`
and the primitive carrier/family norms follow within the same transported
`q` frame.  Together with the inherited clock/source-charge projection norms,
this selects `Ls` and the packet with exactly two lift entries `-1`.

The result is deliberately scoped.  It assumes the unchanged energy

```text
K(p,N) = |p|^2/2 + delta.p + N,
delta_k = (k+alpha)/8 - 1/2,
```

the actual `D8` currents, source cocharacter, retained genuinely eight-periodic
U(1)-dressed clock with inherited square `C=g^2`, and the original
determinant-only packet formula.  It does not derive those structures from raw
P1/P2, derive an E8 field extension, select physical time or state, or close a
physical gate.  The separate `Ls` lift with projection pair `(8,0)` is retained
as an independent `C` comparison, not as the square of the same unchanged `g`.

The new source-operator audit identifies an actual operator before making any
E8 comparison.  Multiplication by `exp(i theta)` on the original chiral CAR
one-particle space has fourth power

```text
B4(x) = exp(i x/2) P
```

on folded source-fermion `V16_source`, where `P` exchanges the two four-component channel blocks.  It is
CAR preserving and boundary equivariant, shifts the folded derivative by
`1/2`, and its Fock implementer `W4` obeys, on the common finite-energy core,
`[Q,W4]=4 W4` and
`[K_alpha,W4]=W4(Q/2+3/4+alpha/2)`.  This is a genuine charge-four chiral
spectral flow under the retained source premises.  Those premises are not
selected by raw P1/P2, and the global implementer is not thereby a local
vertex or a native field dictionary.  Here `theta` is the chosen source
circle.  It is not yet the P1 seam-core circle; reading source winding four as
core winding two needs an additional double-cover identification.  This
`V16_source` is also distinct from the target
`V16=(10,1)+(1,6)` in the packet/E8 dictionary.  The existing `M8` label map
is not a Clifford/CAR automorphism, so the operator transfer between those
spaces remains unproved.

On charge labels the fourth winding has the exact affine form
`T=tau_s rho`, with `rho(p)=Pp+r0`.  Here `rho` preserves `K` and `Q`,
exchanges `0` with `r0` and `D8` with `D8+s`, and separates the half-root
translation from the reflected residue label.  This affine formula alone does not supply a local involution, a field cocycle,
a physical memory system, or the common operator identity.  The adjusted
global even-sector identity proved below uses the additional stipulated E8
and explicit cocycle premises.  In particular `W4`
reflects `beta=e0-e4` and permutes oscillators, so it is not the same operator
as pure E8 half-translation.  The unmodified geometric segment transition is
`t_alpha=exp(2 pi i alpha Q/8) g_0`; `g_0` is the compensated periodic
clock.  Its `W4` character retains both the quarter phase and `(-1)^Q`.

Inside the separately stipulated unitary integrable E8 level-one
representation `Ls`, the affine factor now has a genuine global operator
lift.  Four mutually orthogonal affine `su2` triples construct a unitary `N`
with `N^2=-I`; hence `R=iN` satisfies

```text
R^2=I,   [Q,R]=[K,R]=0.
```

Its adjoint action implements `rho` on all Heisenberg modes and root currents,
including the affine central term, and preserves the quarter transition.  The
lift is global.  It is not a local memory-field insertion, a vacuum-preserving
VOA automorphism, or an E8 selection from raw P1/P2.  This fixed Tits/current
lift `R` need not equal the adjusted lift below.

On the original **even `D8` Fock sector**, relative to the same stipulated E8
completion, the common global factorization is now explicit.  The ordered
source cocycle extends from `D8` to `Ls`; the torus factor
`Z_a=exp(i pi a.H_0)` cancels the actual CAR Klein signs.  For the displayed
allowed convention,

```text
R_c = Z_a E_r0 P_hat,   U_s = -E_s,
W4 = U_s R_c,           R_c^2 = I,
[Q,R_c] = [K,R_c] = 0.
```

The equality includes every oscillator and the affine central term.  It does
not make `R_c` identical to the fixed Tits `R`, and the original even source
does not canonically select between the two allowed `kappa` lifts.  Odd vector
CAR fields still require spin/defect-module transport or the complement
construction; no ordinary common local odd-CAR/E8 algebra, local memory
insertion, native-field map, raw-P1 selection, or event law follows.

The finite v988 half-string is a useful negative control.  It leaves the full
Wilson loop unchanged and is charge neutral; the opposite endpoint cancels
the apparent seam flip.  Its scalar sign comparison is not a monodromy test,
and it cannot be normalized into charge-four transport.  The native `f/b`
operators are still not identified with either source operator.

Run the exact checker with:

```bash
cd experiments/theory-contracts/seam-joint-ward-state-20261014
python3 checker.py
```

Expected output:

```text
178/178 checks passed
67fb9789a0c623c5a22a7f94ae2ce32cdfedee2a1bf94ee6bd7107d10e57d903
```

The checker uses exact `Fraction` arithmetic for the lattice and energy
statements and exact SymPy matrices for the double-cover and two-mode Fock
checks.  The finite root census is only a cross-check; the unbounded-lattice
minimum follows from the coordinatewise convexity and parity argument in
`PROOF.txt`.

Inside the stipulated multiplicity-one `Ls` representation, the full
normal-ordered Cartan Ward values select the initial normal ground density
`|0><0|`.  This does not make the Hamiltonian ground eigenspace
one-dimensional: `q`, clock, and deck remain blind to the second zero vector,
and the physical origin of the state rule is still open.

The event-operator audit also corrects three overclaims.  An allowed original
Cartan term can restrict to `diag(0,1)` on the ground doublet;
`cos(4 theta)` folds to an alternating channel operator rather than a uniform
`q` identity; and a stationary pure vacuum under autonomous selfadjoint
dynamics forbids a direct off-diagonal doublet mixer.  These facts do not give
a physical event-law no-go.  Actual `A_E`/`Phi` coefficients and virtual
effects through the excited spectrum remain open.

The retained fold-time cocharacter also has a non-scalar family projection.
For every dressed-clock offset, `delta_F` has squared norm `41/256` or
`33/256` and six-weight spectrum `{+/-x,+/-1/4,+/-1/4}`, with
`x=3/16` or `1/16`.  This is a result inside the stipulated E8/A3 frame.
Its field action is `rho_R(delta_F)`, while `ad(delta_F)h=0` on the
geometric Cartan residue labels.  It therefore gives neither their time
evolution nor a native W insertion or necessary bilocal vertex.  The raw
reference's `Y_I=0` observation does not select its vertex, and compressed
one-leg actions exist on the antisymmetric pair space.

The existing native-clock dictionary also gives a positive, explicitly
conditional bridge of weight labels.  With its supplied finite metric `N`, marked `M3/U3`
correspondences, and `L=HX/2`, the six minimal marks become plus/minus three
orthonormal geometric `A3` Cartan-residue axes.  This does not identify the
residue vector with a parallel `End(F3)` operator, the neutral `D5` memory root,
or native `U3` with source `g^2`.  In the old one-`A(h)` class, changing among
these axes retains the single-`C` squared spectrum `(1,1,4)`; common-axis real
clock words remain CP-even there.  The rank-three label map is not a
complex-linear, CCR-preserving isometry of the native six-boson module; an
actual source-selected 6-to-3 operator intertwiner is still missing.  No general mass, CP, multi-vertex or TFPT
no-go follows.

Conversely, if one additionally identifies the determinant lift with the
separate irreducible triplet plus singlet and requires it to be parallel, its
family norm is `3 t^2 (2r+1)^2`, hence neither `4` nor `2`.  Only this extra
parallel determinant-only comparison class is excluded.

Run the exact projection/packet addendum with:

```bash
python3 projection_check.py --repo ../../..
```

Expected output: `149/149 checks passed`.

Run the exact source-operator addendum with:

```bash
python3 operator_transport_check.py --repo ../../..
```

Expected output: `208/208 checks passed`.  The normal and optimized
certificates are byte-identical.  The checker calls the actual original
`verification/v988_psi_lambda_reduction.py` constructor for the finite-source
control and does not substitute a new Hamiltonian.

Files:

- `PROOF.txt`: scoped mathematical proof and source reconciliation.
- `checker.py`, `results.json`: independent exact checks and certificate.
- `projection_check.py`, `projection-results.json`,
  `projection-results-optimized.json`: 149 exact checks of the
  intrinsic mixed-root source direction, fixed-root components, projection
  norms, integral packet class, separate parallel flavor comparison, and
  constructive Cartan/native-clock bridge;
  normal and optimized certificates are byte-identical.
- `REVIEW.txt`: independent review of the full-lattice and current-algebra
  steps.
- `operator_transport_check.py`, `operator-transport-results.json`,
  `operator-transport-results-optimized.json`: 208 exact checks of the chiral
  winding, affine factorization, global neutral E8 current-algebra lift,
  common global even-sector cocycle/Klein factorization in the stipulated
  `Ls`, CAR action and original v988 negative control; normal and optimized
  certificates are byte-identical.
- `OPERATOR-REVIEW.txt`: independent source/operator and scope review.
- `source_manifest.json`: pinned local source evidence.
- `claude_reference/`: reproduced upstream source files at commits `603e5d91`,
  `b14ed366`, `f3e91a87`, and the corrected event-scope reference
  `6b9ba728`; it intentionally contains no nested `contract_index.json`.
  The latter reproduces `176/176` locally with byte-identical results.  The
  `c8aa5fa8` time/family reference reproduces `86/86` byte-identically
  (`83` rational/exact checks and `3` NumPy-tolerance controls), while its
  bilocal-necessity inference is not imported.  The
  `80e302e5` native-vertex note is retained as a raw, unexecuted 18-check
  PARTIAL scope report: its mixed `(10,6)` weight comparison is useful, but
  noncommutation with scalar `WW^dagger=8 I_60` is not a valid test.  None of
  these reference contracts' metadata is imported.

Verdict: `EXCLUDED_MARKED_FOLD_COMPLETION` only within the historical
unique-vacuum/unchanged-`K`/determinant-only comparison class.  It does not
exclude the positive global even-sector factorization above.
This remains an experiment-only contract with no ledger, verification, paper,
website, Lean, or full-TFPT status upgrade.
