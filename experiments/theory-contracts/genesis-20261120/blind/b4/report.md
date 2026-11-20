# Blind agent b4 — four rules, two exhaustive tests

Roles: self-compilers / fixed points; selection between competing consistency conditions.
Worked blind (own mathematics + Python only). Labels: ED = exact derivation, DA = derivation under
named extra assumption, NO = numerical observation, HA = heuristic analogy, UC = untested conjecture.

## H1 — The self-interpreting magma (ontology: one binary operation, nothing else)

1. Primitives: a finite set X and one operation `*`. No program/data/machine split: every element is
   simultaneously an argument and a transformation (y -> x*y).
2. Rule: **L : X -> End(X,*), x -> (y -> x*y), is a bijection.** The elements are exactly the
   structure-preserving self-maps of the structure they form ("X = [X -> X]" internal to X's own category).
3. Why: X ≅ Set(X,X) forces |X|=1 (n = n^n); X ≅ Hom_lin(X,X) forces dim 1 (n = n^2). The only way a
   finite thing can equal its own function space is if "function" means "morphism of the structure
   itself". Alternatives compared: all n^(n^2) magmas; Aut instead of End; anti-endomorphisms.
4. Forced consequences (ED): (a) left self-distributivity x*(y*z)=(x*y)*(x*z) (each L_x is an
   endomorphism); (b) a left identity e exists (id ∈ End) — a "universal machine" element;
   (c) for every idempotent a, const_a is an endomorphism, hence some element k_a with k_a*y=a for all y:
   every stable value has an eraser; e is idempotent so an eraser-to-e always exists for n>1;
   (d) the reversible variant (L: X ≅ Aut(X)) has no solution for 2<=n<=4 (NO, by enumeration; no proof).
   (e) n=2 solution is unique: e*y=y, k*y=e. This is Boolean implication with e=True
   (True->y = y, False->y = True). ED by hand, confirmed by enumeration.
5. Emergence: states = elements; dynamics = iterated self-application (UC). Irreversibility is built in
   (erasers are forced, so no solution is a quasigroup): an arrow of time before any time (HA).
   QM, space, gauge structure: UC — nothing derived. Linearised version (algebra A with
   A ≅ End_alg(A)) not tested.
6. Falsifier: if solutions exist for every large n with no rigidity (count growing freely), the rule
   selects nothing. If no infinite/limit object with nontrivial symmetry exists, it cannot reach physics.

## H2 — Exact self-distinction (ontology: distinctions only)

1. Primitives: a finite set V of distinctions; two distinctions combine by "differ in exactly one"
   (so V is an F2-vector space — this is the definition of a distinction, not an extra field choice).
2. Rule: **the system distinguishes exactly its own distinctions and nothing distinguishes itself:**
   a bilinear B: V×V -> F2 with x -> B(x,·) a bijection V -> V*, and B(x,x)=0.
3. Why: compared with all 2^(n^2) bilinear forms. Degenerate = some distinction is invisible (not
   "exactly"); B(x,x)=1 = a distinction that separates itself from the null distinction (Russell-type
   self-application). The surviving class is a single GL-orbit.
4. Consequences (ED, classical theorems): dim V = 2m is even; B is unique up to isomorphism
   (symplectic); B is not a coboundary-free object: the 2-cocycle β with β(x,y)-β(y,x)=B(x,y) defines a
   central extension. Over centre Z2 there are two (real/quaternionic) extraspecial groups; over centre
   Z4 exactly one, with a unique faithful irreducible representation, of dimension 2^m over **C**, on which
   V acts as the m-qubit Pauli group; non-commutation = B. Complex numbers appear as the minimal centre
   that makes the extension unique (ED for the group theory).
5. States = vectors of that irrep: DA (assumption: "state" := irreducible carrier of the extension).
   Born rule: DA (Gleason-type, needs the full projective unitary group, whereas the rule only gives the
   Clifford normaliser Sp(2m,2) — stabiliser QM, classically simulable; universality is missing).
   Space/time/gauge group: UC.
6. Falsifier: if no principled extension from Clifford to universal dynamics exists without an added
   non-Clifford primitive, the rule gives only a classically simulable shadow of QM.

## H3 — The description fixed point (ontology: one infinite word)

1. Primitives: a word s over a set of positive integers; D = run-length description (the shortest
   lossless "say what you see" code).
2. Rule: **D(s) = s.** The world is the word that is its own description.
3. Why: compared with all description alphabets {a,b}. Alphabet {1} gives only the trivial word "1";
   {1,2} is the minimal nontrivial alphabet, and the fixed point is then unique up to the first symbol
   (Kolakoski word) — construction is forced symbol by symbol (ED).
4. Consequences: (ED) not eventually periodic: a period p containing r runs would give D(s)=s a period
   r<p, contradicting minimality. (ED) D^-1 is an inflation, so s carries an exact discrete
   self-similarity with no input scale. (UC, famous open problem) letter density 1/2, inflation
   factor 3/2.
5. States = finite factors; "time" = inflation depth; hierarchy of scales from powers of the inflation
   factor: all HA/UC. No route to Hilbert space identified.
6. Falsifier: non-existence of letter frequencies (then no invariant measure, no statistical physics);
   or proof that the factor complexity/dynamics is isomorphic to a trivial odometer.

## H4 — Timeless admissible histories (ontology: whole histories, no step rule)

1. Primitives: an alphabet A and one ternary relation R ⊂ A^3, the only admissibility constraint,
   imposed at every vertex of a trivalent network of A-labelled edges.
2. Rule: **(P) any two legs determine the third; (S) R is invariant under all permutations of its legs;
   (X) gluing two constraints through a hidden leg gives a 4-leg relation independent of the pairing:
   (a*b)*c = (a*c)*b.** No leg is "time", no vertex is "earlier", and the constraint does not depend
   on how a history is cut into pieces.
3. Why: (P) is the weakest condition making every direction a Cauchy direction (differs from
   't Hooft/Wolfram: no update rule and no preferred foliation; time is a gauge choice). (S),(X) are
   the statements that decomposition of a history is unobservable. Compared with: all Latin squares.
4. Consequence (ED, 3 lines): (S) gives commutativity and x*(x*y)=y; (X) with commutativity gives
   (b*c)*a=(b*a)*c=(a*b)*c, i.e. associativity; so (A,*) is an abelian group of exponent 2:
   **A = F2^k, R = {x+y+z=0}.** Admissible histories = F2^k-valued 1-cocycles of the network: they
   depend only on its first homology (a Z2^k gauge theory at zero coupling); local degrees of freedom
   vanish — a topological, background-free theory.
5. The rule lands on H2's primitive (F2-spaces) from an independent starting point (ED that both are
   F2; the pairing B is NOT produced by H4). States = cocycle classes; QM via H2: DA. Lorentzian
   3+1: UC — a topological theory has no light cones; a metric would need (X) to be broken or weighted.
   Three families/SU(3)xSU(2)xU(1): UC, nothing found.
6. Falsifier: the result itself is a partial falsification: (X) kills all local propagation, so the
   rule as stated cannot give gravity with local degrees of freedom.

## Computation 1 (H1) — `h1_self_interpreter.py`, `h1_n5.py`

Exhaustive over all magmas (the rule implies left-distributivity, used only as sound pruning).

| n | magmas in class | LD | LD+faithful | RULE (L: X≅End) labelled / iso | Aut-version | anti-End control |
|---|---|---|---|---|---|---|
| 1 | 1 | 1 | 1 | 1 / 1 | 1 | 1 |
| 2 | 16 | 9 | 5 | 2 / **1** | 0 | 1 |
| 3 | 19683 | 224 | 65 | 6 / **1** | 0 | 1 |
| 4 | 4^16 | 14067 | 2207 | 48 / **2** | 0 | 205 |
| 5 | left-identity class 5^20, complete | – | 7987 | 0 / **0** | not completed | – |
| 6 | incomplete (45 s timeout) | – | – | ≥ 2 iso classes | 0 so far | – |

Tables: n=2 `((0,1),(0,0))` (implication); n=3 `((0,0,0),(0,1,2),(1,1,1))`; n=4
`((0,0,0,0),(0,1,2,3),(0,1,1,3),(1,1,1,1))` and `((0,1,2,3),(0,0,0,0),(0,1,0,1),(0,2,0,2))`.

Reading: the rule is extremely selective (1–2 classes out of up to 4.3·10^9 magmas) and not vacuous,
but it does **not** single out one object: solutions exist at n=1,2,3,4,6, none at n=5, and the n=6
search is incomplete. The Aut-version (reversible self-interpreter) has no solution for 2≤n≤4 (n=5 not completed for
this variant): on this evidence a self-interpreter is irreversible. The anti-End control is satisfiable too (205 at n=4, but the
control as coded only requires anti-endomorphisms ⊂ rows, a weaker condition — not a clean control).
No physical structure is visible in the solutions; no numbers were matched to anything.

## Computation 2 (H4) — `h4_histories.py`

All Latin squares n≤5; n=6..9 by backtracking with the rule as pruning.

| n | all | S only | X only | assoc | RULE S+X (labelled / iso) |
|---|---|---|---|---|---|
| 1 | 1 | 1 | 1 | 1 | 1 / 1 |
| 2 | 2 | 2 | 2 | 2 | 2 / 1 |
| 3 | 12 | 3 | 6 | 3 | 0 |
| 4 | 576 | 16 | 96 | 16 | 4 / 1 (Klein four-group) |
| 5 | 161280 | 30 | 720 | 30 | 0 |
| 6 | – | 480 | – | – | 0 |
| 7 | – | 1290 | – | – | 0 |
| 8 | – | ≥21794 (timeout) | – | – | 240 (first table verified: identity, associative, exponent 2) |
| 9 | – | timeout | – | – | 0 |

Reading: the rule singles out exactly n = 1,2,4,8 with a unique class at n≤4 (240 = 8!/|GL(3,2)| = 8!/168
is the labelled count of one class at n=8, consistent with uniqueness, but iso classes at n=8 were
not separately canonicalised). Each control alone is permissive (S: Steiner-type systems at 3,5,6,7;
X: 7 and 10 classes at n=4,5, including subtraction). Caveat against over-reading: this confirms a
theorem provable in three lines; "1,2,4,8" here means 2^k, not the division-algebra dimensions —
n=16 would also pass.

## Honest overall verdict

- H4 is the only rule that singles out a unique family, and it is forced into a topological theory
  (negative for gravity/locality as stated).
- H1 is selective and yields implication at n=2 and forced irreversibility, but no uniqueness.
- H2 gives complex stabiliser QM exactly, stopping short of universality and the Born rule.
- None of the four reaches 3+1 dimensions, the gauge group, or three families; all such claims are UC.

Files: report.md, h1_self_interpreter.py, h1_output.txt, h1_n5.py, h1_n5_output.txt,
h4_histories.py, h4_output.txt (this directory).
