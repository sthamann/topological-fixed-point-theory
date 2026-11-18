# Recovering the marked V16 action from the original packet

14 October 2026. `UR.SEAM.V16_OPERATOR_PULLBACK.01`. **PARTIAL**.

This continuation addresses the actual original seam, not a new eight-copy
QWZ model. It constructs an inverse for the connection/action already carried
by the original geometric spinor/family packet, and integrates the later
Möbius-clock folding proposal in [FOLD-AUDIT.md](FOLD-AUDIT.md).
Neither result is presented as a complete raw-seam field derivation.

## Original data and the positive construction

The original bundle result is

    V16=E_R + Fix(J6:Lambda²F -> Lambda²F),
    sum m_j=-2, beta_j=m_j+1/2.

Its local family weights are beta_j times the determinant connection of E.
The half-determinant contributions cancel in the existing rank64 complex
packet

    W_m=Lambda^even E tensor sum_j (det E)^(m_j).

For a supplied Hermitian U5 connection A_E, t=tr A_E, the existing exact
packet connection is

    A_W=rho_even(A_E) tensor I4 + I16 tensor diag(m_j) t.

This is equation (1) of `source-boundary-lift-selection-20260922/PROOF.md`.
The data include the chosen lift m; neither A_E nor m is selected by the
present continuation. A further supplied SU4 family connection is transported
by its exterior-square representation in the same way.

The realification of W has dimension128. With a coherent chirality and
orientation dictionary it is the half-spin representation of the geometric
V16. This bundle-representation identification is denoted Upsilon here.
It is NOT the requested raw multiplicity map

    Xi: Hom_Cl(S,E_DSigma) -> (V16)_C.

They have different sources, targets and roles. A change of representation
from a128-dimensional spinor packet to a16-dimensional vector representation
is not a unitary identification of their physical field Hilbert spaces.

## A fixed inverse, including every mixed action

Use Euclidean Hermitian Clifford generators {gamma_a,gamma_b}=2 delta_ab,
and Sigma_ab=gamma_a gamma_b/2 for a<b. On either complex chiral128 module,

    tr(Sigma_ab Sigma_cd)=-32 delta_(ab),(cd).

For equal pairs Sigma²=-I/4. A distinct pair product has degree2 or4 and
zero chiral trace. Thus the Lie-algebra action is injective and, for an
actual packet connection in its image,

    omega_ab=-tr(Sigma_ab Omega_W)/32,
    Omega_V=sum_(a<b) omega_ab L_ab,
    L_ab=e_a e_b^T-e_b e_a^T.

No new rotation coefficient is fitted. The inverse recovers it from the
packet action. Before using it as an identification, compute

    R=Omega_W-sum_(a<b) omega_ab Sigma_ab.

R must vanish. Otherwise the formula is only an orthogonal projection and
all of R remains source data; it may not be discarded. The checker retains
an explicit skewadjoint quartic-Clifford term with zero projected coefficients
and nonzero residual.

Apply the formula in the spin frame supplied by Upsilon. If the frame varies,
the transformed connection must include its derivative. For new fields
psi_new=Upsilon psi_old,

    Omega_new=Upsilon Omega_old Upsilon^-1-(dUpsilon)Upsilon^-1.

A pointwise conjugation without this term is not a connection transport.
No raw field frame Upsilon is inferred merely from equal dimensions.

## Exact application to the original lift formula

Choose the three real family planes corresponding to the pairs (i,4) and
their complements in Lambda²F. Their charges under the determinant line are

    a_i=m_i+m_4+1, i=1,2,3.

The source-packet inverse gives exactly

    Omega_V=rho_R(i A_E) + rho_6(i A_F),
    A_F=diag(beta_j) tr A_E,

on E_R plus the real family6. This is the same induced connection, now
reconstructed backwards from the packet rather than chosen on V16.

The checker independently builds the actual formula on the64 complex packet,
completes its realification using charge conjugation, and applies the trace
inverse. It compares the result with the independently realified U5 matrix
and the determinant charges a_i. All25 real Hermitian U5 basis directions
are checked for three existing valid lifts,75 cases in total. These matrices
are representation basis probes, not selected physical couplings or new
source models. Linearity and the displayed formulas give the general result.
In this K_f-reduced formula the connection has no60 cross components.
That fact does NOT forbid mixed currents: a block-diagonal background
connection and the existence of a mixed symmetry-current operator are
different statements.

## The original C/X word has a concrete V16 action

The existing phase-aligned current implementation is executed, not replaced:
`reconstruct_marked_fermion_source(0)` reads the original Chevalley data and
returns

    [C_0,X_(0,1)]=E_(-e1-e6), coefficient +1.

It also explicitly guards against interpreting that algebraic source as an
already derived microscopic P1/P2 field. Preserve that guard.

In the compatible spin representation let A=c_0 c_5 be the corresponding
negative-root operator. Its compact real span has the two generators

    A-A*=Sigma_0,10-Sigma_1,11,
    i(A+A*)=-Sigma_0,11-Sigma_1,10.

The trace inverse therefore yields the explicit cross actions

    L_0,10-L_1,11,  -L_0,11-L_1,10

on V16. The checker constructs the annihilation/creation matrices directly
and independently recovers these coefficients. This fixes a normalized
root-vector convention for this marked representation; it is not a new
microscopic implementation or a complete recheck of all lattice cocycles.
The same inverse checks every one of the60 real10-by6 mixed directions.

This supplies a positive operator dictionary. Once the actual source action
is known in this packet/frame, its vector action is calculable, including
the cross terms. The calculation does not assert that the raw source has
already supplied those fields or the corresponding covariance Ward identity.

## Global transport and actual raw differential operator

The half-spin representation is locally faithful at the Lie-algebra level
but has a central group kernel. With chirality fixed, the volume element
or its negative acts as +I on that spinor and as -I on V16. Thus spinorial
holonomy alone does not determine a vector transition's central sign.
The specified original E/F/K_f lift resolves this within its given branch;
it may not be replaced by a freely selected spin structure.

The raw boundary operator has geometric symbol, coefficient bundle,
connection, zero-order term and a domain. The local inverse just constructed
is a connection/action map. It does not determine the raw Clifford symbol,
identify its coefficient bundle with V16, or replace its original Phi.
It does not transport state/time between inequivalent field representations.

The original archive's primitive projector/3+2 literal-collapse issue was
already recorded in `source-local-line-20260920`; it is not a fresh result
of this continuation. The valid internal carrier grade after retaining the
P2/branch data is `(2X-I)/5`, not another restriction of the seam-even
projector. We do not repeat that excluded route as a new origin attempt.

## Verification and precise completion state

The trace inverse,75 source-connection cases,60 mixed actions, original C/X
witness and residual control pass32,010 exact checks. Phases and entries are
Gaussian integers or dyadic halves of bounded magnitude, represented exactly
by the checker; no fitted tolerances or asymptotic extrapolation are used.
An independent review checked the32 normalization, frame derivative,
realification/central distinction and the limits of the connection claim.

The fold's explicitly transported source charge supplies one additional
connection case. In its fixed E8/ordered-plane dictionary the original
lift equations select m=(0,-1,-1,0), and the original packet action
reconstructs exactly -2 L_4,5-2 L_10,11 on V16. This joins the two maps;
it does not derive the upstream quotient/E8 field dictionary. Including
this charge-matching case, the final operator checker has32,133 checks.

The accompanying fold audit adds a unitary source-operator pushforward and
an explicit periodic Weyl/triality conjugacy. Together these are constructive
bridges toward the requested source action. They do not establish the
original source's quotient locality, Pin lift, E8 extension, or raw Xi.
No complete TFPT/T1–T8 solution, raw local mixed current, Lean proof,
verification promotion or ledger change is claimed.
