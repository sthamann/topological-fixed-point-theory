# Reproduzierte Programmausgaben

Python 3.14.3, SymPy 1.14.0. Alle acht Programme beendeten sich erfolgreich.

## associativity.py

```text
arity 1 bracketings 1 branching 0
arity 2 bracketings 1 branching 0
arity 3 bracketings 2 branching 0
arity 4 bracketings 5 branching 1
permutation 123 component 5 arrows 5 exact_signature True
permutation 132 component 10 arrows 10 exact_signature False
permutation 213 component 10 arrows 10 exact_signature False
permutation 231 component 10 arrows 10 exact_signature False
permutation 312 component 4 arrows 4 exact_signature False
permutation 321 component 4 arrows 4 exact_signature False
three-variable right-bracketed control hits 1 / 6
all different three-variable binary right sides: 1 / 11
  (((ab)c)d) extremal -> ((a(bc))d),((ab)(cd))
  ((a(bc))d) transit -> (a((bc)d))
  ((ab)(cd)) transit -> (a(b(cd)))
  (a((bc)d)) transit -> (a(b(cd)))
  (a(b(cd))) extremal -> 
cycle directions: degree split 20 / 32; and arm lengths 2,3 10 / 32
arity 5: pentagons 6 vertices 14 vertices with inconsistent inherited roles 8
  (((a(bc))d)e) [('first four', 'transit'), ('fused bc', 'extremal')]
  (((ab)(cd))e) [('first four', 'transit'), ('fused ab', 'transit'), ('fused cd', 'extremal')]
  (((ab)c)(de)) [('fused ab', 'transit'), ('fused de', 'extremal')]
  ((a(b(cd)))e) [('first four', 'extremal'), ('fused cd', 'transit')]
  ((ab)(c(de))) [('fused ab', 'extremal'), ('fused de', 'transit')]
  (a(((bc)d)e)) [('last four', 'extremal'), ('fused bc', 'transit')]
  (a((bc)(de))) [('last four', 'transit'), ('fused bc', 'extremal'), ('fused de', 'transit')]
  (a(b((cd)e))) [('last four', 'transit'), ('fused cd', 'extremal')]
all ten local two-sets: conflict counts [7, 7, 7, 7, 8, 8, 8, 8, 9, 9]
possible off-diagonal local matrix entries (unordered) 0 / 10
size of forced diagonal equality class 5 / 5
smallest basis-state lift retaining the original tree and its local role 22
arity 4 trees 5 mixed 0 minimal tree/role states 5
arity 5 trees 14 mixed 8 minimal tree/role states 22
arity 6 trees 42 mixed 35 minimal tree/role states 77
arity 7 trees 132 mixed 124 minimal tree/role states 256
arity 8 trees 429 mixed 420 minimal tree/role states 849
arity 9 trees 1430 mixed 1420 minimal tree/role states 2850
```

## braid_origin.py

```text
Nonidentity endpoints, odd/even: 3 2
Characters (identity, swap, 3-cycle), endpoints vs meromorphic: [[5, 1, 2], [5, 3, 2]]
Five pole-mode determinant / Vandermonde: 1
Elliptic j at t=-1, 3, (1+i sqrt(3))/2: [1728, 21952/9, 0]
Hidden commutator: endpoint=id, winding=0, torus= [[13, 8], [8, 5]] trace= 18
Orders of hidden commutator on torsion: {2: 1, 3: 4, 4: 1, 5: 10, 6: 4, 7: 8}
Identical-cost pure exchange commutators on 3-torsion: {1: 8, 4: 8}
Order four in SL(2,F3): 6/24
Braid E[2] eigendirections over minimal splitting field F4: 3 + 2
Exact h=Weil pairing after Frobenius; endpoint-equivariant bijections: 2 / 120
Ordered Vandermonde fiber: Y²=1−4X³, j=0
Shortest winding-zero, endpoint-neutral torus witnesses: 8/108, length 4, trace 6
All hidden words on 3-torsion: Q8; element orders {1: 1, 2: 1, 4: 6}
Fixed real Hermitian forms (a,d,Re z,Im z): {(0, 0, 0, y)} primitive Eisenstein determinant = -3
Invariant positive semidefinite forms: zero only
Moving positive Hodge metric: determinant 1, covariant under both braid generators
Torsion permutation trace, reflected Gram on a,a^2: [[1/3, 1], [1, 1/3]] eigenvalues= [-2/3, 4/3]
Reflection Gram on positive words <=3: rank 6 at r=1/2, 2/3, 1
Same double exchange response: [1/4, 4/9, 1]
```

## hermitian_origin.py

```text
null rays : count (nondegenerate sesquilinear): [(0, 60), (1, 90), (3, 30)]
null rays : count (all Hermitian): [(1, 5), (3, 10), (5, 1)]
null rays : count (nondegenerate Hermitian): [(3, 10)]
distinct marked triples / realisations each: 10 {1}
diagonal form null rays: [(1, 1), (1, 2), (1, 3)]
q, total, null, nonnull: 2 5 3 2
q, total, null, nonnull: 3 10 4 6
q, total, null, nonnull: 4 17 5 12
q, total, null, nonnull: 5 26 6 20
native quadratic/marked agreement: 16/16 labels, 5/5 rays; Arf = 1
same-cost fixed-coordinate quadratic formulas: 1/18; GL truth table + mark: 12 / 20160 (quadratic alone: 120 )
integral polarization eigenvalues: [-3, -1, -1, 1] discriminant: -3
finite intervals [0,nI]: [(1, 4), (2, 15), (4, 103), (8, 1305), (16, 20041)]
integral SU(2) site image: 6 involutions: 3 specified native permutation: (1, 0, 3, 2, 4)
primitive lifted traces on ordered five rays: [1, 1, 2, 2, 2]
|U(2,2)|, |SU(2,2)|, |projective site image|: 18 6 6
site character on 1, transposition, 3-cycle: (5,1,2); equivariant self-bijections: 2
native involution cycle type: 3 / 6; exact named action: 1 / 6
normalized RP controls at one tick: 3/8 5/16
```

## kontrollen/bell_carrier.py

```text
Local-product dimension minus uniform composite dimension: n**2*r*(n - 1)**2*(r - 2)/4
Scalar dimension, local-product count, joint count: [(1, 9, 10), (2, 16, 16), (4, 36, 28)]
Real locally indistinguishable states; global witness: [-1, 1]
Bell outcome difference lattice: [[2, -1, 0], [-1, 2, -1], [0, -1, 2]] ; determinant = 4 ; roots = 12
Axial quarter turn fixes all four outcome projectors.
Fixed-degree labelled graph baseline: 30240 / 11180820
Dodecahedral control: Petersen edges = 15 ; independent four-sets = 5 ; rank of product with labelled Atlas Gram = 0
Both fields: CHSH = 2*sqrt(2) ; squared table = Matrix([[1/2, 1/2], [1/2, 1/2]])
Ranks of sign-blind two-trial responses (real bit, complex bit): {2: 2, 3: 5}
Axial generator characteristic polynomial: t*(t**2 + 1)*(t**2 + 4)
Half-turn eigenspaces: {1: 3, -1: 2}
Nonclassical sign games: 8 / 16
Game x two actually realized Bell fields: 8 / 32
After conditioning on nonclassical games (R,C): 8 / 16
Game x three candidate scalar algebras with tomography filter: 8 / 48
```

## kontrollen/operation_exchange.py

```text
choices for one representative per product class: [2, 8, 3]
involutive exchanges by rank(I-B): {0: 1, 1: 1, 3: 2, 4: 2, 6: 3, 7: 3}
involutive Yang-Baxter exchanges by rank(I-B): {0: 1, 6: 2, 7: 2}
nondegenerate involutive exchanges: 0
signed involutive Yang-Baxter exchanges by rank(I-B): {0: 1, 6: 4, 7: 8, 9: 8, 10: 4, 15: 4, 16: 8, 18: 8, 19: 4, 25: 1}
Smith diagonals for all U(1) phase extensions: {(1, 1, 1, 1, 1, 2): 1, (1, 1, 1, 2, 2, 0): 2, (1, 1, 2, 2, 2, 0): 2}
all U(1) components / rank-ten components: 26 2
alternating dimensions 0..6 for real rank-ten representatives: {(1, 5, 10, 20, 40, 80, 160): 4}
Hurwitz cycle lengths: [1, 1, 1, 1, 1, 2, 3, 3, 4, 4, 4]
ranks for Hurwitz and ordinary flip: 14 10
ordinary flip, products changed: 18 / 25
basis dimensions of (occupied set, ordered product), and oriented quotient: [1, 5, 19, 30, 15, 3] [1, 5, 10, 10, 5, 1]
worst-case distinguishable memory states by event number: [1, 1, 2, 3, 3, 3]
optimal squared scale of the coherent erasure Kraus operator: 1/3
native and class-averaged X,Z factors: ('2/3', '1/3') ('1/2', '1/2')
native distance from class-averaged label distribution: 1/9
```

## kontrollen/operation_origin.py

```text
order  reduced Latin squares  labelled groups  group types  noncommuting types
1 1 1 1 0
2 1 1 1 0
3 1 1 1 0
4 4 4 2 0
5 56 6 1 0
6 9408 80 2 1
aut size / nonidentity orbit sizes: 6 [3, 2]
3 and 2, order 6: 1 / 2
3 and 2, all types through 6: 1 / 8
true-set sizes among eight one-product equalities: {0: 4, 3: 1, 5: 3}
predicates selecting exactly the involutions: ['xx=e']
distinct semantic predicates: 3
```

## kontrollen/operational_rotations.py

```text
SO(d) invariant precessions for d=2..12: [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
Centered half-turn span: 5 ; marked ranks: (3, 2)
Centered primitive rotation spans for q=2..12: [5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
Five-dimensional span: 1/11 order symbols; 1/45 primitive angles.
Labeled sign-involution control (3+2): 6/8
Squared principal cosines: {x, (2*x-1)^2}; x=cos(theta)^2
Picard/weak common orbit at {1/2,1/2}: impossible (exact gcd=1)
Finite polyhedral quadratic character norms: [3, 2, 1] (A4,S4,A5)
Maximal commuting pulse contexts: 5 ; pair orders: {2: 30, 3: 60, 5: 120}
Shortest noncommuting recurrence: 3 ; fixed/moved: (2, 3) ; unordered marks: 10
Dodecahedral test against Q4 Gram: sum=2I, product=0; reversed triangle signs: 120
Petersen-preserving labelings: 120 ; matching after signed unrestricted relabeling: 720
```

## kontrollen/qec_origin.py

```text
minimum n, attained d: 5 3
correctable erasure sets by lost-site count: {0: 1, 1: 5, 2: 10, 3: 0, 4: 0, 5: 0}
complete distinct syndrome table: 2592 / 782595 = 96/28985
incidence automorphisms / site permutations: 360 / 120
Frobenius site involutions / site involutions: 10 / 25
Frobenius fixed triples / unconditionally invariant triples: 10 / 0
```

## Unabhängiger Lauf des Originalmoduls

Am geprüften Repository Stand a1273187b82488d34d2a474470152a41d42d4055 wurde zusätzlich in research/originlab die Funktion f1_binary_clifford_origin.controls() ausgeführt.

```text
Original controls: 51 ; status: EXACT_SCOPED_CONTROLS_PASSED_F1_OPEN
Original two windows: [2, 3] ; original quarter: [1, 0, 3, 2, 4]
```

