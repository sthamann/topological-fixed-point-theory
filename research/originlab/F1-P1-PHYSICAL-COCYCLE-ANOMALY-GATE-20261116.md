# An independent P1 source route: local gauge gluing, nontrivial cocycle and semionic quarter phases

**16 November 2026. Exact conditional finite and topological source tests.
This is a new external physically motivated possibility, NOT native TFPT F1.**

## Independent established literature

- Kitaev toric code: local four-edge star A_v and plaquette B_p
  gauge/flux checks on a square lattice; exact commuting projector
  topological code. https://www.nature.com/articles/s41467-024-45584-6
- Hu, Wan, Wu, Phys Rev B 87, 125114 (2013), twisted quantum double
  lattice models classified by normalized G-group 3-cocycles:
  https://doi.org/10.1103/PhysRevB.87.125114
- Bullivant, Hu, Wan, Phys Rev B 96, 165138 (2017),
  boundaries classified by subgroup K and boundary 2-cochain
  subject to a Frobenius compatibility constraint:
  https://doi.org/10.1103/PhysRevB.96.165138
- Lan, Kong, Wen, 2016, 2+1D bosonic topological orders and
  E8 invertible stacking shifting c by 8 without new anyons:
  https://academic.oup.com/nsr/article/3/1/68/2460241

## Route 1: P1 as physical Gauss vertex, not chosen singularity

The physical square-lattice toric code uses edge qubits with
A_v=product_(e incident v) X_e and B_p=product_(e in boundary p) Z_e.
Star and plaquette checks overlap in zero or two edges and commute;
an odd overlap would anticommute and spoil simultaneous admissibility.

We checked the exact local 4-qubit overlap normal form
A=X tensor X tensor X tensor X,
B=Z tensor Z tensor I tensor I.
Their projectors (I+A)/2 and (I+B)/2 commute, but their
joint +1 sector has rank FOUR in dimension sixteen, not one.
The four marks / one Gauss vertex / one flux check is a possible
local gauge-source grammar, not yet a uniquely selected P1 vacuum.
Additional bulk constraints and choice of topology/boundary are
required. A square spatial lattice is assumed, not derived, and
the four incident edges are not yet the original TFPT four faces.

## Route 2: complex phase from physical associativity obstruction

For G=Z2, write a,b,c in {0,1} with group addition mod two and
choose the nontrivial normalized 3-cocycle

    omega(a,b,c)=(-1)^(a*b*c).

All sixteen quadruple cocycle identities were checked exactly.
For normalized two-cochains beta with beta(0,g)=beta(g,0)=1
and arbitrary nonzero beta(1,1)=z, the three-coboundary obeys

    (delta beta)(1,1,1)=+1,
    omega(1,1,1)=-1.

Hence no such two-cochain on the FULL subgroup K=Z2
trivializes omega|_K: the simplest fully matched gapped
boundary of the specific Bullivant-Hu-Wan subgroup/2-cochain
construction is obstructed. The trivial K={0} restriction
is unobstructed. This does not forbid all boundaries or
interface degrees of freedom in the double-semion phase.
Indeed physical boundary structure is constrained by anomaly
cancellation, rather than any arbitrary microscopic projection
P_adm. The theorem is conditional on the TOPICAL GAUGE model G=Z2
and cannot be inserted as native TFPT source without deriving G,
omega and the actual outer selector.

## Route 3: derive semionic fourth roots, then discover an E8 ambiguity

The standard double-semion Abelian K-matrix is
K_DS=diag(2,-2), with determinant -4. Its charge labels
(0,0),(1,0),(0,1),(1,1) have topological twists

    theta(l)=exp(pi*i*l^T*K_DS^(-1)*l)
            =1, i, -i, 1, respectively.

Thus a physically established topological associativity twist
produces fourth roots of unity WITHOUT arbitrary complex
path weights. This is a possible physical source model of
TFPT-like quarter phases, but the four anyon TYPES are not
yet the four TFPT marked points, and their braiding is 2+1D.

The normalized Gauss-Milgram sum is

    [1+i-i+1]/2=1=exp(2*pi*i*c_minus/8).

Thus topological twists determine c_minus only modulo 8.
The E8 positive even unimodular integral Cartan K_E8 has
rank 8, determinant 1 and positive signature 8. Stacking
K_DS direct sum K_E8 leaves |det|=4 and the same four
anyon twist labels but shifts c_minus: 0 -> 8.

This proves an exact original-source ambiguity:

   (DS, four quarter spins, c=0) and
   (DS tensor E8, identical four spins, c=8)

have identical anyon types, fusion/twist data in this
comparison. The absolute vacuum modular T factor
exp(-2*pi*i*c/24) differs by exp(-2*pi*i/3).
The *relative* quarter spins cannot choose the
native E8 gravitational framing/anomaly. Original
Ref/Eta/framing and absolute R113 Z0 must select
a physical source, not a name match to E8.

This 2+1D topological order and its E8 layer do NOT
derive a chiral interacting 3+1D Standard Model,
local universal spin-2, or pregeometric P1. Stacking
E8 is an extra physical choice.

## Integration with new projective P1 result

The previous candidate "match + Z4 global clock" produced
pure GHZ_n at each block, but was inconsistent under
physical marginal trace from n to n-1. Topological lattice
models avoid the demand for pure *local* states:
one global state defines generally mixed local restrictions.
True boundary gluing uses an explicitly specified
topological/categorical bulk-to-boundary operation, not
ordinary marginal trace relabeled as an amplitude cap.

The local gauge plus twisted-cochain approach supplies an
established template for F1 source typing: original
bulk algebra/group, 3-cocycle or associator, restricted
boundary data, and absolute framing must be defined
together; only THEN q0/P_adm, physical clock and the
original R113 joint reflected Gram can be evaluated.
Neither the topology nor symmetry determines the
source's absolute scalar or continuum time.

## The decisive native-source falsifier

Can the **original** four-mark q0 and the actual original
outer P_adm be obtained from one primitive algebraic
boundary admissibility condition, and can the *same*
boson/fermion/eta/R113 source detect its physical
complex relative phases, exact absolute vacuum phase,
and weighted two-cutoff pushforward?

If no, reject a topological-cocycle route as physical F1.
If yes, one STILL must reconstruct genuine 3+1D
locality, interacting chiral matter and physical spin2
from this same source. No F1 completion is claimed.

Reproduce: f1_p1_cocycle_boundary_gate.py and
test_f1_p1_cocycle_boundary_gate.py.
