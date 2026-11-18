# Original seam, mixed weight-one current, and the proposed selector

14 October 2026. `UR.SEAM.MIXED_CURRENT_ORIGIN.01`. Verdict **PARTIAL**.

The requested result is a local mixed spinor weight-one field derived from
the original TFPT seam, preserving its source, state, charged products and
time. This document does not claim that result. It completes the identification
test of the pasted proposal against the actual available microscopic source
and constructs the charge-sensitive test that its scalar pair test lacks.
No source copies, hopping coefficients, temperatures, stress tensors or
observables have been selected to fit the desired answer.

## 1. Source boundary and what is actually given

The Atlas identifies `SEAM.SIMPLECURRENT.GENERATOR.01` and
`SEAM.MMST.TYPEIII.CHARGED.01` as open at microscopic charged-field realization.
The actual one-copy QWZ source used by the charged-limit contracts is
`verification/v1033_charged_disorder.py`, mass one and width eight. Its
transitively pinned adapter is `microscopic-charged-car-limit/checker.py`.
The present checker imports and executes that adapter and the original
`strip_at_momentum`; it does not rewrite the strip or import an eight-copy
target Hamiltonian. The established charged CAR scaling theorem and the
global crossing argument are inherited from the cited proofs, not proven by
29 finite checks.

The original operational definition is also important. In archived
`01_boundary_kernel_source.tex:480–520`, the seed already consists of a local
net, time, reflection, state, seam class and elliptic collar operator.
Its almost-commutative factorization at lines 84–91 is explicitly a
realization choice. The original geometric construction of
`V16=E_R+Fix(Lambda²F)` and its Spin(16) lift is a real positive result;
it is not the identification of this source operator's coefficient bundle
with V16. The synthesis, lines 4425–4470, states exactly the required
symbol/field map and preserves that distinction.

Thus “original seam” must retain two meanings here: the proposed primitive
TFPT source, still specified through structural conditions, and the explicit
one-copy microscopic representative used by existing limit proofs. Excluding
a construction in the latter is not an exclusion of the former.

## 2. Exact check of the unchanged representative

For its 16-component transverse Bloch matrix,

    h(p)=H0−sin(p)X+(1−cos(p))Z,
    H0³=H0, rank H0=14,
    {H0,X}={Z,X}=0.

The zero space has the two normalized vectors
`rt=(0,...,0,1,1)/sqrt(2)` and `rb=(1,-1,0,...,0)/sqrt(2)`.
On it `h'(0)=diag(−1,+1)`. Also

    h(p)²=sin²(p) I+[H0+(1−cos(p))Z]².

A real zero requires `sin p=0`; the original matrix `H0+2Z` is invertible.
There are no other Fermi crossings. The established top-edge limit has one
complex chiral channel; the bottom channel has opposite orientation.
The eight transverse rows are not eight co-propagating species.

For every transverse bilinear B, its leading top-edge current coefficient
is `b(B)=rt*B rt`. The Gram on all 256 matrix units is `b b*`, rank one.
This previously established gate is independently re-executed here using
exact Gaussian dyadic coefficients. It concerns the established CAR subnet,
not an asserted exhaustion theorem for every microscopic scaling sequence.

## 3. All quarter-twist fusions in that channel, not just bilinears

Consider the same one-channel free U(1) Heisenberg/Fock conformal structure.
Integer-charge CAR fields have integral charge. Even if quarter-twist
intertwiners of charge ±1/4 are admitted provisionally, every finite OPE
fusion has charge `q=m/4`, `m in Z`. In the corresponding untwisted
Heisenberg oscillator module its states have weights

    h=m²/32+N,  N in Z_{≥0}.

The integer oscillator levels for actual source backgrounds are in
`source-half-sector-bridge/README.md`; for these provisional quarter-twist
modules the assertion follows from the free-boson Heisenberg module class.
It is not a proof that an actual microscopic local quarter-twist field exists.

Weight one gives `m²=32(1−N)`. For N≥2 the right side is negative; for N=0
32 is not a square; for N=1, m=0. Thus the only weight-one state in this
fusion class is the neutral oscillator current. This is an all-fusions
argument, not an unbounded conclusion from scanning a finite list.
The check enumerates precisely the exhaustive range |m|≤5, N≤1 implied
by h≤1.

In particular two quarter twists give q=1/2 and h=1/8. Eight independent
such channels give h=1 for `(1/2)^8`; eight repetitions of an operator in
one channel add its charge and are not eight independent channels.
Adding descendants cannot lower h. This excludes the suggested fusion
as a derivation of the requested current inside the known one-channel
Fock/OPE class. It does not exclude other source realizations, other scaling
fields, a selected extension, or the full TFPT theory.

All weights in this section are conformal L0 weights. Actual quarter-holonomy
source time is `H_a=L0+(a−1/2)Q`; at a=1/4 it is `L0−Q/4`.
The backgrounds r=1 and r=3 are not extra species and no intersector field
convergence follows from having separate convergent background theories.

## 4. Even eight CAR copies do not contain the spinor field internally

This is a conditional check of the *already used* eight-channel target,
not a new source proposal. The even NS Fock module of 2n real fermions has
exactly `choose(2n,2)` states at weight one: two distinct creation modes
at level 1/2. Four creation modes have weight at least two, and a level
3/2 mode cannot contribute to weight one. For n=1 the dimension is one;
for n=8 it is 120, the D8 current algebra.

Equivalently the even eight-complex-fermion charge lattice is
`D8={z in Z^8: sum z even}`. The mixed spinor
`lambda=(1/2)^8` is outside this lattice. CAR products and their OPEs do
not introduce its half-integral Cartan charge. A weight-one field acting
on the same NS vacuum Hilbert space must create a weight-one vector there;
no such spinor vector is present. This does not prohibit a local extension
on a larger vacuum Hilbert space, nor invalidate the localized spinor DHR
sector construction. A sector representation is different from a field
creating that sector from the original vacuum.

The extension `D8 union (D8+lambda)` adds the 128 spinor ground states and
gives 120+128=248 E8 currents. The base D5+D3 currents plus an actual
spinor root and its adjoint generate all 240 E8 roots by nonzero lattice
VOA root brackets. The checker closes the exact root sums. A mixed vector
root instead closes only the 112 D8 roots. These positive algebraic
constructions use the actual field extension as a premise.

The standard distinction is visible directly in Boyle Smith et al.,
arXiv:2303.16917, section 5.1, equations 5.3–5.6: the E8 vacuum character
contains the D8 vacuum and one spinor character, whereas the NS free-fermion
character contains the vacuum and vector character. Boeckenhauer,
arXiv:hep-th/9507047, sections 4–7, constructs the distinct representations
and localized sector morphisms. These are established results, not new TFPT
origin theorems.

## 5. Why the proposed pair test cannot identify the field

In the normalized target, both a neutral abelian current j and a charged
root current V_lambda have two-point function `(w−u)^−2`.
Conformal mapping gives the corresponding circle kernel, with its standard
coordinate factors. Equality of this scalar kernel does not determine the
charge or its source operator. On z=0, w=1, u=3 both two-point values are 1/4.

The charge-sensitive Ward identity is

    <H_i(z) V_lambda(w) V_{−lambda}(u)>
      =lambda_i (w−u)^−2 [(z−w)^−1−(z−u)^−1].

For lambda=(1/2)^8 all eight values at those points are −1/12.
For the neutral abelian j the corresponding three-point is zero.
The residues measure the *actual eight charges*, not just their norm.
They must be determined using the same source fields, Cartan action,
state and common time; inserting target charges is not a source measurement.

There is also a higher-word counterexample: on 0,1,3,7 the alternating
charged four-point is 81/784, while the abelian-current Wick four-point is
the different rational sum stored in results.json. Thus even the specific
two-point equality has a direct distinguishing test.

More generally, the D8 net and its E8 extension have identical vacuum
correlators and L0 evolution on the *shared D8 vacuum subnet*: its cyclic
subspace is invariant and the extension restricts to its original vacuum
representation. No amount of those existing-subnet correlators alone
decides whether the spinor extension is physically included. This statement
does not identify full Hilbert spaces, charged words or original source time.
An E0 average must act on full words: a supplied charged pair can survive
even when its individual fields do not. We do not identify the original
glue character with the old finite 16-dimensional v993 clock.

Ward residues are necessary identification evidence, not sufficient proof
of microscopic existence, relative locality, energy bounds or convergence.

## 6. The pasted stronger selector, fully audited within its lattice class

Let `R_d=A_{d−1}+D_{d+1}`, d≥2. A-class k has minimum
`k(d−k)/d`, D-vector class minimum one, D-spinor class minimum (d+1)/4.
Within a fixed simply laced root-lattice coset the norm differs from its
minimum by a nonnegative even integer. A new nontrivial norm-two root
outside R_d therefore requires the stated minima to sum to two.
“Mixed” here means nonzero charge on both factors, not a neutral linear
combination of their Cartan generators.

Spinor equation: `4k(d−k)=d(7−d)`; d>7 is excluded immediately.
The solutions are (4,1),(4,3),(7,0). Only the first two are mixed.
Vector equation: `k(d−k)=d`; taking k≤d/2 gives
`d=k+1+1/(k−1)`, hence only (4,2).
Pure A equation: `k(d−k)=2d`; k≤d/2 implies k≤4, and direct solution
gives (8,4),(9,3),(9,6). The pure exceptions yield respectively
A6+E8, E7+D9 and E8+D10. They are new root directions, not a claim that
the already present base roots disappear for other d.

For d=4 the spinor glue has order four, so det16/index4²=1;
the generated root lattice is E8. Vector glue has order two and determinant
four, giving D8. “One spinor produces E8” must mean adjoining the field
and its adjoint to the *base R_d current algebra*, not that powers of one
single operator generate an entire nonabelian algebra without the base.

The 240 roots have glue counts 52/64/60/64; the full 248 weight-one currents
include eight neutral Cartan fields and have counts 60/64/60/64.
This corrects the pasted root/current terminology. The stronger selector
is exact; applying it to a scalar pair kernel remains invalid.

There is no extra spinor-versus-vector bit after the *actual original glue*
grade 1 or 3 has already been established on the proposed field: every
half-root is odd-grade, every integer D8 root even-grade. The vector
extension lies in grade 2. This inference cannot use the unrelated old
finite clock without proving its operator dictionary. Moreover, an
eight-copy test with the stipulated 5+3 block and that glue class already
implements the d=4 target; it can validate a realization but cannot prove
that the primitive seam selected d=4. Selection and target validation must
not be interchanged.

## 7. Constructive fusion retained, without changing its premise

The original A3→D3 isometry sends the family fundamental weight to
`(1/2,1/2,1/2)`. Together with the D5 spinor weight it gives lambda.
In the stipulated eight-channel free-boson system put q=lambda/2.
Then q²=1/2, h(q)=1/4, and

    V_q(z)V_q(w)=c(q,q)(z−w)^(1/2)[V_lambda(w)+higher terms], c(q,q)≠0.

This is the valid leading charged intertwiner coefficient in that system.
It is not the action of two finite clock buttons, not a CAR-vacuum
observable, and not a proof that the actual one-copy seam supplies it.
Lambda is relatively local to the even D8 lattice but has half-integral
pairing with an elementary CAR charge e_i. The bosonic spinor extension
and the unprojected CAR field algebra cannot both be called the same
ordinary local net.

## 8. Decision and verification

The original executable source check passes. The proposed derivation from
its known one-channel quarter-twist fusion class fails by an exhaustive
weight/charge argument. The pasted algebraic selector succeeds, after
correcting its root count and the meaning of “one generator”; its pair
identification fails by an explicit charge-sensitive counterexample.
Both conclusions are more specific than merely marking the bridge open.

The remaining origin problem has not been solved here. The existing source
is not its solution, and the abstract existence of the geometric V16 bundle
does not specify the actual operator/field map. The decisive usable target
is the original coefficient-bundle and Cartan-Ward action, with the actual
charged fusion and common-time realization. A matched two-point formula,
added copies or source-subnet neutral reconstruction cannot replace it.
No T1–T8 gate, full TFPT solution, Lean theorem, ledger status or paper
promotion is claimed.

29 checks pass, with exact arithmetic after validating the original source's
Gaussian dyadic entries. Analytic continuum and sector facts have the
explicit upstream/literature premises above. An independent reviewer
checked the Fock-state count, subnet ambiguity, Ward test and all-fusions
argument; the review is stored separately. This is not independent
verification of an entire continuum or physical theory.
