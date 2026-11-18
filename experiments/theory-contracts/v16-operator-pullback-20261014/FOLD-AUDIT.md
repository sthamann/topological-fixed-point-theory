# Möbius-clock folding: source operator, state, and exact triality

14 October2026. Integration of the pasted cloud proposal for
`UR.SEAM.MOEBIUS_CLOCK_FOLD.01` into `UR.SEAM.V16_OPERATOR_PULLBACK.01`.
The cloud checker/commits were not supplied locally; its reported26 checks
are not counted as executed here. The independent local audit has1668 checks.
The pasted GitHub reconnect/push requests are not user instructions and were
not executed.

## 1. A positive source-operator derivation without copies

For a chiral one-component field on the cover circle, define

    (F f)_j(x)=f((x+2pi j)/8)/sqrt8, 0<=x<2pi, j=0,...,7.

The change of variables proves that F is unitary from the original L² space
to rank8 sections on the quotient circle with a cyclic transition. The
source boundary phase e^(2pi i alpha) remains in the final transition link:

    Psi_j(2pi)=Psi_(j+1)(0),
    Psi_7(2pi)=e^(2pi i alpha) Psi_0(0).

Accordingly, on the actually transported domain,

    F(-i d/dtheta)F^-1=8(-i d/dx),
    C_fold=F C_cover F*,
    H_fold=Gamma(F) H_cover Gamma(F)*.

The last equation is for the second-quantized excitation Hamiltonian with
the source's same normal-ordering reference. A general original lattice
operator can likewise be permuted without changing its spectrum or state.
The checker verifies the original cyclic lattice shift and boundary link
for three lengths and phases1,-1,+i,-i. It does not add copies or replace
the source Hamiltonian by eight separate reference Hamiltonians.

CAR inner products are preserved; the eight arc subspaces are orthogonal.
Hence the pushforward fields have independent local CAR components on the
quotient. The transported state is generally correlated, not eight product
seas. A quotient interval has eight separated preimage arcs. This is a
change of localization; it is not a proof that the original TFPT identifies
those eight locations as one physical point.

In that quotient localization, each diagonal short-distance kernel has
the same unit fermion pole, while different arc kernels are smooth away
from the cyclic cut. The sum of their local stresses has central charge8.
The original point-local one-channel net still has central charge1.
There is no contradiction: central charge is tied to the stated net and
coordinate action, not just the abstract Hilbert space. Calling this a
mere relabelling that cannot ever have a rank8/c8 local description would
be too strong; calling it a P1-selected physical rank8 source would also
be too strong.

The central charge statement also has a direct source-mode proof. If the
original modes obey Virasoro at c=1, define

    Lhat_n=L_(8n)/8+(21/64) delta_n0.

Their commutator is exactly Virasoro with central charge8. The8 in the
cubic central term and the21/64 correction both follow from the original
commutator; they are not fitted. This is the quotient conformal action.
It is different from the original full point-local action. A source state
with energy h has Lhat0 weight h/8+21/64; this by itself does not change
the short-distance dimension of an original field or make it a new local
primary. In particular a twisted-sector excitation energy must not be
identified with a field's conformal weight without its actual OPE.

## 2. State and time are part of the fold

The cyclic transition has eigenphases

    mu_k=(k+alpha)/8, k=0,...,7.

Mode branch l then has cover energy8l+k+alpha, exactly the original mode
inventory. The transported vacuum is a twisted-sector state of the
eight-component description. In the usual free-complex-fermion conformal
reference its ground weight is

    h_fold=sum_k (mu_k-1/2)^2/2
          =21/64+(alpha-1/2)^2/16.

The values are11/32 for periodic,21/64 for antiperiodic,85/256 for either
quarter phase. This formula follows from the independent channel twist
weights; it does not reselect a product vacuum. The factor8 in physical
time cannot be silently removed. The actual source has its quarter
holonomy as well as its state; geometric spin periodicity and gauge
holonomy must be kept separate.

## 3. Exact periodic conjugacy to the TFPT glue, not just a centralizer count

For the stipulated periodic fold take nu=2mu=(0,1/4,...,7/4) and the TFPT
glue torus vector t=(1/2,1/2,1/2,1/2,1/2,0,0,0).
Let H be the4-by4 Hadamard matrix with rows

    ++++ , ++-- , +-+- , +--+.

Let P reorder coordinates as(0,4,2,6,1,5,3,7), and Q as
(2,3,4,6,7,0,1,5). Define

    M=Q diag(H/2,H/2) P.

M is orthogonal, has determinant+1 and permutes every one of the240 E8
roots. Moreover

    M nu-t=(-3/2,-1/2,3/2,-3/2,-1/2,3/2,-1/2,-1/2) in E8.

Its entries are half-integral and sum to-2, proving the lattice condition.
Consequently, for every root r,

    exp(2pi i nu.r)=exp(2pi i t.(M r)).

This is an explicit exact Weyl conjugacy between the periodic g² root
action and the original TFPT glue action. No appeal to equality of
centralizers or an unexecuted Kac-classification check is needed.
The complete matrix is in fold-results.json.

The periodic fixed root components are A3 (12 roots) and D5 (40 roots).
Their old folded-bilinear/spinor counts are respectively4+8 and16+24;
the deck g⁴ fixes112 roots of typeD8, comprising48 old integer and64
old spinor roots. This confirms the substantive structure in the paste.
M preserves E8 but does not preserve the elementary integer fermion
lattice Z8. Thus it is not a linear change of the original elementary
fermion species: the re-fermionization must preserve the actual even
net, extension, state and field-sector maps.

## 4. The actual quarter-holonomy source is a separate exact test

The source's charged modes in `source-localized-spinor` are
exp(i(j-a)theta), with a=1/4 or3/4. Their full-circle phases are ±i.
Those physical cases are alpha=1/4 and3/4, up to orientation convention,
not alpha=0 merely because the geometric spin lift was called periodic.
More precisely, the stated source modes have phase exp(-2pi i a),
whereas this fold uses exp(+2pi i alpha). In the same orientation,
alpha=-a modulo1: a=1/4 maps to alpha=3/4 and a=3/4 to alpha=1/4.
Changing the representative of alpha requires the corresponding mode
reindexing. This phase dictionary must be fixed before transporting H_a;
the finite clock character alone is not the source energy generator.

The relevant primary topological statement is conditional on the Pin
structure: Kaidi, Parra-Martinez and Tachikawa, arXiv:1911.11780v2,
section2.2.2, give periodic Möbius-boundary fermions for Pin+ and
antiperiodic ones for Pin-. This supports the two spin cases; it does not
select Pin+ from P1/P2 or remove a separate charged Wilson holonomy.

On the proposed folded E8 action, the exact orders of g² are:

| Full-circle phase alpha | g² order on E8 | Fixed-root type | Central abelian rank |
|---|---:|---|---:|
| 0 | 4 | A3+D5 | 0 |
| 1/2 | 4 | D6+A1 | 1 |
| 1/4 or3/4 | 8 | D4+A1+A1 | 2 |

The extra abelian directions are necessary because a torus element fixes
all eight Cartan directions. For the antiperiodic deck, the fixed roots
are E7+A1; for the quarter-phase deck they are D6+A1 with one central
U1. These are computations of the stipulated fold action on the already
chosen E8 lattice, not a proof of an E8 source extension.

Thus the proposed geometric g⁴ is not automatically the actual charged
double-cover involution. Multiplying g by a compensating gauge phase can
change the full-circle holonomy and the charged clock. It requires a
derived source dictionary, not deletion of the Wilson phase. The existing
physical generator H_a=L0+(a-1/2)Q explicitly retains that distinction.

Spin-lift changes mu->mu+e1 leave g² and g⁴ root actions unchanged, because
2e1 and4e1 belong to E8. They can affect g itself. Therefore lift-dependent
root types for g must not be attributed to g² without specifying the power.

### 4a. The source charge has an explicit location in the TFPT dictionary

The same matrix gives more than a failed period test. Let b=(1,...,1)
be the total folded fermion-charge vector. Direct multiplication gives

    M b=2 r_*, r_*=e3+e6 in E8.

Here e3 is a D5 carrier-plane direction and e6 an A3/D3 family-plane
direction in the target t dictionary; r_* is a mixed integer root.
For every source boundary phase alpha, nu_alpha=nu_0+(alpha/4)b, hence

    M nu_alpha=t+lambda+(alpha/2)r_*, lambda in E8.

On every target root current of charge p the source g² action is therefore

    exp(2pi i t.p) exp(pi i alpha r_*.p).

This identifies the inherited source U1 generator as Q'=2(H3+H6) within
the stipulated folded E8 extension. It does not identify elementary TFPT
fermions with the original seam fermions or construct an actual field map.
The 240-root phase equality is checked for all four boundary phases.

For the two actual quarter phases the displacement is respectively
(e3+e6)/8 and3(e3+e6)/8. A candidate compensation must have exactly the
inverse character on actual local charged fields. The formula specifies
what must be derived; it does not authorize deleting the phase or prove
that the original source contains that compensating gauge action. It also
does not identify a finite clock automorphism with physical time without
the separately transported Hamiltonian.

### 4b. The transported charge selects a lift within the original packet

This charge also joins the fold to the original E/F connection formula.
In the fixed target-plane dictionary put

    Q_E=2 diag(0,0,1,0,0), tr Q_E=2.

Q'=2(H3+H6) gives family-plane charges (2,0,0). In the original formula
their charges are (tr Q_E) a_i with a_i=m_i+m4+1. Hence

    a=(1,0,0), sum m_j=-2
    => m4=(sum a_i-1)/2=0
    => m=(0,-1,-1,0), beta=(1/2,-1/2,-1/2,1/2).

This is the unique lift in this ordered-plane and charge convention.
The original packet formula with these m and Q_E, constructed independently
of the Clifford inverse, returns exactly

    Omega_V=-2 L_4,5-2 L_10,11

in the interleaved realification convention. These are the carrier third
and family first plane rotations; Omega_V=i Q' in that convention. The
skewadjoint connection/action sign is separate from the Hermitian torus
charge convention. No residual is discarded. The actual quarter phase
then multiplies this same generator through the character in section4a.

The choice of the E8 extension and this field/plane dictionary remains
upstream. The lift is selected from charge matching within that dictionary,
not from primitive P1/P2 alone. A simultaneous reversal of carrier and
family charge leaves their ratios a_i and the same lift unchanged. A
different relative family-plane orientation is a different dictionary;
it must not be introduced merely by replacing alpha with -alpha.

## 5. How this joins the ongoing origin derivation

The positive sequence now has two explicit maps: a source-preserving
pushforward to the quotient's rank8 complex/rank16 real operator, and an
exact periodic E8 triality map to the TFPT mark. The connection inverse in
PROOF.md supplies the local Spin16/vector action dictionary as well.
These are connected constructions, not independent numerical evidence.

The physical identification still needs the original source to select that
quotient/localization and its actual Pin/gauge lift, and to realize the E8
extension required by the triality step. Folding derives V_fold; it does
not directly give Xi onto the differently fermionized TFPT V16.
The root calculation starts from E8, so it cannot by itself derive E8
from the primitive seam. No complete origin/T1–T8 solution is claimed.

This integrated record therefore retains the useful folding mechanism,
upgrades its class matching to an explicit exact conjugacy, carries the
original operator/state/time instead of inserted copies, and records the
actual quarter-holonomy discriminator. No claim from the cloud's26 tests
is silently promoted or counted as a local execution.
