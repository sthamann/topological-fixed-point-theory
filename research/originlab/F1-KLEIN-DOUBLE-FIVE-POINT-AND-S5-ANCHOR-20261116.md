# F1: Klein group doubled to native odd five-point code, but full S5 forbids canonical P1

**16 November 2026. New exact mathematical bridge and a stronger source-first
symmetry no-go. ORIGINAL TFPT all-field R113 Q0, absolute Z0, OS time,
P_adm and F1 are NOT solved.**

## 1. Four-sector source -> sixteen-dimensional phase space, no E8/spacetime axiom

Start with the optional algebraic source G=F2^2 (four elements),
motivated conditionally by the 2026 Fidaleo-Vincenzi classification
of certain infinite twisted exchangeable algebras,
arXiv:2609.09494. The original TFPT primitive has not yet
been proven to belong to that class.

Introduce its Pontryagin dual G*=Hom(G,Z2), and the canonical
Heisenberg/doubled phase space

    V=G + G* = F2^4,
    B((g,chi),(h,psi))=chi(h)+psi(g) mod2.

This B is bilinear, alternating and nondegenerate. It has the
correct mathematical SOURCE TYPE for the existing TFPT
four-coordinate symplectic quadratic code. Importantly, the
Klein theorem's *nonalternating* 2D bicharacter u_K is NOT
identified with B; one obtains B only after dual/doubling.

## 2. Additional odd quadratic/spin polarization -> exactly five code points

Choose ONE further (UNDERIVED) polarizing axis e1 in G and its
dual e1* in G*. With coordinates v=(a,b,c,d), use

    q(v)=ac+bd+a+c mod2.

Its polar form is precisely B. The Arf invariant is odd
(6 zeros including 0000; 10 ones).

The five NONZERO q=0 vectors, ordered lexicographically, are

    0001 0100 0111 1101 1111.

They obey B(vi,vj)=1 for all i!=j, their total xor is 0000,
and their first four vectors form a basis (integral det -1).
Thus one recovers exactly a finite five-point symplectic simplex
from a phase-space double of FOUR source sectors, without
postulating the E8 algebra or the original four face labels.

This is isomorphic, as odd quadratic F2^4 space, to the native
TFPT marked_search five-null-point code; an explicit
physical identification of their elementary operators and
orientation is still required. The use of q above is an
EXTRA spin-polarization source choice, not implied by u_K.

## 3. NEW IMPORTANT ORIGIN NO-GO: full O-minus(4,2) = S5 has no fixed fifth point

The q-preserving symplectic automorphism group is

    O^-(4,2) ≅ S5, order 5!=120.

Proof in this concrete model: every permutation of the five
q-null nonzero vectors extends uniquely to an F2-linear
automorphism of V because their only linear relation is
v1+...+v5=0. All pairwise symplectic pairings are one, so
each extension preserves B; its q preservation follows because
the 16 vectors consist of 0, five null singletons and
ten q=1 pair sums. Conversely every q-isometry must
permute the five null points, hence the full isometry
group is exactly S5.

The executable test constructs all 120 maps explicitly,
checks q on all 16 vectors and B on all 256 vector pairs
for each map, and verifies one-point stabilizer order 24.

**Consequently not one of the five possible anchors is
intrinsic under the UNMARKED original q symmetry.**
Selecting a fifth/P1 code point reduces the symmetry to its
stabilizer S4 (order 24). The previous observation
"1111 is fixed by permutations of the four pre-given marks"
is true only AFTER the S4 four-mark action is taken as
physical input; it cannot replace the S5->S4 selection theorem.

This materially changes the F1 source question. The original
physical Q0/P_adm/OS history must supply a *source-derived
symmetry-breaking observable/sector/defect* that selects a
point AND an orienting spin refinement, without merely
declaring the four remaining original marks elementary.
Otherwise the original origin problem has moved into the
choice of a subgroup rather than been solved.

## 4. Honest compatibility verdict

This route is a new conditional mathematical reconstruction:
4 source event sectors -> G x G* = 16 phase-space labels ->
one odd quadratic q -> five null points -> S5 ->
S4 stabilizer of chosen anchor.

It still does NOT imply:
- the 5th point equals the physical P1 quantum seam,
- which source dynamically breaks S5 to S4,
- which of the four Klein semionic quadratic refinements and
  native two odd F2^4 refinements is physically selected,
- six anchored D4 2+4 weights, original quarter phase/amplitude,
- ordinary trace versus cap versus coding decoder for
  true blocked quantum histories,
- original R113 weighted Q0 and all absolute Z probes,
- an infinite physical modular half-sided inclusion and OS time,
- 3+1D chirality, local E8_1 currents, physical graviton spin2.

The odd code, its automorphism group, and its stabilizer are
necessary mathematical consistency candidates, not a TOE.

Exact origin module:
research/originlab/f1_klein_phase_space_double.py
and test_f1_klein_phase_space_double.py.
