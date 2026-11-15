# TFPT F1: an equivariant 5/3 residue bridge and its orientation/gauge obstruction

**15 November 2026. Research candidate, not a physical F1 theorem.** This connects the actual original four-point Riemann--Roch data to the existing P1 orientation problem. All stated finite character identities have independent exact checks in [seam_cohomology.py](seam_cohomology.py) and [test_seam_cohomology.py](test_seam_cohomology.py). The original weighted R113 packet, physical source state, physical time and 3+1D parent remain unavailable.

## Original owners consulted

* [Shared c3/P1/P2 architecture](../../documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex).
* [Riemann--Roch carrier](../../verification/v189_riemann_roch_carrier.py), [even Clifford closure](../../verification/v197_rr_carrier_clifford_d5.py), [one index gate](../../verification/v228_rr_index_gate.py) and [the original square/flag census](../../verification/v512_seam_tau_flag.py).
* [Original one-sided boundary theorem](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex) and [operational seed](../../documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/002-section-intrinsic-closed-branch-and-canonical-rigid-obje.tex).
* [R113 source consumer](source.py), [full producer trace](../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/triangulation_20261106O/caller/ACTUAL-PRIMITIVE-REGULATOR-CALLER-TRACE.txt), [new joint fragment search](composed_search.py) and [native affine metric](P1-MARKED-AFFINE-METRIC.md).

## The strengthened exact 5 = 3+2 mechanism

For the original square divisor D = mu4, put

  C5 = H0(P1,O(D)), F3 = H0(P1,Omega1(log D)).

In the fixed affine frame, a basis for C5 is f_k(z) = z^k/(z^4-1), k=0..4.
The holomorphic square rotation r(z)=iz gives r* f_k = i^k f_k and the
reflection s(z)=1/z gives s* f_k = -f_(4-k). Here f4-f0 = 1 is the
constant global section.

A basis for F3 is omega_j = z^(j-1) dz/(z^4-1), j=1..3.
Its pullbacks are r* omega_j = i^j omega_j and
s* omega_j = +omega_(4-j).

It follows as a statement about **representation isotypic components**:

  C5|C4 = 2 chi0 + chi1 + chi2 + chi3
  F3|C4 = chi1 + chi2 + chi3.

So the geometric quarter rotation distinguishes the nontrivial 3-dimensional
part of the **same** five-mode carrier from its 2-dimensional fixed part,
rather than only giving the numerical dimensions (5,3). The comparison
C3 vs F3 is equivariant for the rotations but differs by the reflection
orientation character: s_C = -s_F. The transformation of dz supplies
that sign in the actual source geometry; an arbitrary graph orientation
sign need not be inserted.

This does NOT establish that physical color and weak SU(3)xSU(2) act
on these isotypic components. It is a new typed candidate for the
previously physical/conditional 3+2 identification.

## Negative gauge test and minimal compensating family channel

On the nontrivial C3, rotation has three **distinct** eigenvalues
(i,-1,-i). Its complex commutant has dimension 3, the diagonal
algebra; it cannot contain the full End(C3) or commute with a
fundamental irreducible SU(3) action. On the full C5 the rotation
commutant has complex dimension 2²+1²+1²+1²=7, and full D4 reduces it to 4.
This excludes a naive identification in which the physical color gauge
SU(3) and geometric quarter rotation must commute on the **uncompensated**
three character lines.

Proposed compensation: take the dual of the original F3 as a minimal
character ancilla. On C3 tensor F3* the simultaneous rotation has an
invariant projection

  E_C4 = (1 + r + r² + r³)/4,

with rank **3**, and rotation acts trivially on its image. Every
nontrivial C3 character requires the corresponding inverse ancilla
character to produce an invariant channel, so dimension >=3 is a sharp
minimality bound if all three physical source lines must remain visible.
The original F3* attains it.

This is insufficient to establish three independent *families* and
three *colors*. The neutralized rank-3 channel by itself has only three
states, whereas a full triplet with three independent generations
requires nine states in the corresponding tensor product. A separate
physical family action, or a larger original source representation,
must be derived; do not use the same three states twice.

If the source additionally preserves the geometric three-versus-two
type projector, the compensated C3 and the two C4-invariant carrier
slots admit a U(3)xU(2) unitary automorphism group. Requiring determinant
one permits S(U3xU2) with traceless Y = diag(-1/3,-1/3,-1/3,1/2,1/2).
Those *group and unimodularity selections* are conditions, not yet a
derived original physical gauge dictionary.

## The orientation is load-bearing, not cosmetic

With the original geometric reflection on both factors,

  rank E_D4 on C3 tensor F3* = 1.

If the compensating factor is twisted by one reflection-orientation
sign the corresponding rank becomes 2. Only restricting to the
orientation-preserving C4 retains all three channels.

Thus a full linear D4-invariant **scalar** pairing cannot simultaneously
identify all three character channels in this exact minimal tensor
construction. A physical solution may use an orientation-charged
operator, antiunitary reflection, extra sectors, or more complicated
graded action; the finite counterexample does **not** rule out TFPT.

This is a concrete origin hypothesis: the missing P1 orientation/Pfaffian
or eta factor could supply the charge conjugation and oriented twist
needed to reconcile three carriers, three families and the current algebra.
The hypothesis is falsified if the original source requires the literal
untwisted D4-invariant scalar projection **and** all three independent
charged channels on this same minimal carrier.

## Original F1 gates, before promotion

1. Obtain actual weighted Q0[a,b;t] and the bound absolute empty reference
   from the **original** finite BRST/Higgs/gauge/fermion/geometry source,
   before the outer selector and with its physical time dictionary.
2. Identify the actual rotation and reflection operations on the
   reflected OS null quotient. Test whether their action matches the
   character matrices above **on the original field basis**; check
   the correct Calderon/deck rather than asserting P_C=(1+R)/2.
3. Evaluate the charged/neutral orientation word, its genuine adjoint
   and physical gauge Lie algebra. Inspect the native E8 120/128 grade-one
   parity obstruction and any grade-mixing source action.
4. At two actual blocking levels, test conditional fiber integration,
   all four marked forms, reflection positivity, the same reference,
   and the ordered mixed R113 value before regulator limits.
5. Prove independent color and family multiplicities, chiral
   3+1D locality and spin-2 coupling from the **same** interacting source.

**No new physical R113 Z0, ZU, ZY or ZUY is asserted. F1 is open.**

## A relevant new OS warning (13 November 2026)

Handrack and Salmhofer, *Reflection positivity in free fermionic theories*,
Letters in Mathematical Physics 116, article 161 (2026),
https://doi.org/10.1007/s11005-026-02191-4:
for a free flat Dirac **covariance** given by a real rational function,
reflection positivity holds iff its poles are real/simple with
nonnegative residues; they also exhibit exponential **covariance**
regulators violating RP. This does not automatically disprove TFPT's
spectral **action weight** exp[-Tr f_star(D²/Lambda²)]; those are different
objects. But it forbids replacing the weighted original OS/fermion
condition by unproved pointwise positivity. Check the actual supplied
propagator and covariance at every block level.
