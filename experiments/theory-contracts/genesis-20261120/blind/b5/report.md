# Blind disruptor report b5 — five ontologies, five one-line rules, two computations

Worked blind (no repository, graph, web). Labels: ED = exact derivation, DA = derivation under named
extra assumption, NO = numerical observation, HA = heuristic analogy, UC = untested conjecture.
Files: `autology.py`/`autology_output.txt`, `parallax.py`/`parallax_output.txt` (same directory).

## H1 Autology — "a thing is the set of its own changes"
1. Primitives: a finite-dimensional real vector space V and one bilinear form B on it (a bare notion of
   "overlap"; not assumed symmetric or nondegenerate).
2. Rule: **Lie(Aut(V,B)) ≅ V as Aut(V,B)-modules**, where g = {X : XᵀB + BX = 0}.
3. Why: compared with the class "any tensor structure and its stabiliser". Usually object and symmetry
   are two things; the rule is the unique fixed point of "replace an object by its symmetries".
   Without it one must postulate both a space and a group.
4. Forced: for symmetric nondegenerate B, dim so(n) = n(n−1)/2 = n ⇒ n = 3 (ED); the iso is the
   cross product / Hodge map. Over C: so(3,C) ≅ C³, and so(3,C) ≅ so(3,1) as real Lie algebras (ED).
5. Space: 3 dimensions, ED inside the class of nondegenerate symmetric forms; NO for general B (below).
   Lorentz algebra so(3,1): DA (extra assumption: V carries a complex structure; then V = C³ is the space
   of self-dual 2-forms, not Minkowski space; R^{3,1} appears only as Hermitian forms on the spinor of
   sl(2,C) — second extra assumption). Gravity: HA (so(3,C)-valued self-dual variables à la Plebanski).
   QM, Born rule, gauge group, families: UC — nothing derived. Signature not fixed: so(3) and so(2,1)
   both satisfy the rule.
6. Falsifier: any (V,B) with dim ≠ 3 satisfying the rule. **Found**: see computation — a degenerate
   rank-1 form on R² (g = aff(1)) satisfies it. The rule therefore needs "B nondegenerate" as an
   additional axiom, or it must accept a 2-dimensional "scale + translation" world as a second solution.

## H2 Parallax — "reality is the mismatch of mutually incomplete projections"
1. Primitives: a finite set X; projections = partitions of X. No points with properties, only views.
2. Rule: **every two projections are jointly complete and non-redundant (X → P_i × P_j is a
   bijection), and the family is maximal.**
3. Why: alternatives are "one complete projection" (classical state, no mismatch) and "arbitrary
   overlapping projections" (no constraint). Pairwise-exact complementarity is the extremal case in which
   no single view carries any information about another, yet any two exhaust reality.
4. Forced: |X| = n², every projection has n classes of size n; k projections ⇔ k−2 mutually orthogonal
   Latin squares; k ≤ n+1 (ED: a point lies in one class of each projection; the k classes through it share no second point, so k(n−1) ≤ n²−1).
   k = n+1 ⇔ affine plane of order n.
5. Mathematics as emergent: for n ≤ 5 the maximal family forces the coordinate loop to be the
   elementary abelian group, i.e. the additive group of F_n (NO, computed; known theorem). QM: the n+1
   projections correspond to n+1 mutually unbiased bases of Cⁿ, translations give a Weyl–Heisenberg
   group (HA; close to Wootters' discrete phase space — essential difference: here the projections are the
   ontology and the Hilbert space is derived, not the reverse). Space, time, gauge, gravity: UC.
6. Falsifier: (a) an n where k = n+1 is reached without a field: known to exist at n = 9
   (non-Desarguesian planes; quoted from memory, not computed) — so "field emerges" is false beyond small
   n; (b) the rule selects no particular n (all prime powers) — it singles out a family, not a world.

## H3 Biprojection — "observer and observed are the two marginals of one relation"
1. Primitives: one real square array M (entry = how much i registers j); no states, no observers apart
   from row/column roles.
2. Rule: **M Mᵀ = Mᵀ M** (the overlap structure seen from the row role equals that seen from the
   column role).
3. Why: alternatives: M symmetric (roles identical, no arrow), M arbitrary (roles unrelated). Normality
   is the largest class in which the two projections cannot be told apart.
4. Forced (ED): a real normal matrix is orthogonally block-diagonal with 1×1 blocks and 2×2 blocks
   [[a,−b],[b,a]] = a+ib: complex numbers appear exactly where M ≠ Mᵀ. Polar factors M = UP commute
   (reversible part and weight part are independent), and spec(MᵀM) = |λ|²: co-registration counts are
   squared moduli of the eigen-amplitudes.
5. Complex amplitudes: ED (as linear algebra). Born rule: HA (|λ|² is a spectrum identity, not a
   probability rule). Unitary time = powers of U: DA (assumes iteration is time). Space, gauge, gravity: UC.
6. Falsifier: a single normal M generates a commutative *-algebra, hence no contextuality. If no
   canonical second relation forced by the rule fails to commute with M, the model is classical. I expect
   this to kill it unless composition (M₁⊗M₂ vs M₁M₂) forces non-commutativity.

## H4 Kernel — "existence is the description that survives its own negation"
1. Primitives: a set X of candidate statements and one relation R ("refutes").
2. Rule: **K = X ∖ R[K]** (what exists is exactly what the existing does not refute).
3. Why: alternatives K ⊆ X∖R[K] (merely consistent; many solutions) and K ⊇ X∖R[K] (merely complete).
   Equality is the only self-referential closure. (Near Spencer-Brown/Kauffman re-entry; difference: a
   relation on many statements, with solution-counting as the physics.)
4. Forced (ED): R acyclic ⇒ unique K, obtained by backward induction, which defines a well-founded
   order (a proto-time). A directed odd cycle ⇒ no K. In the linear relaxation x = 1 − Ax, x ∈ [0,1]^X,
   the cycle gives x_i = 1 − x_{i−1}: odd length forces x ≡ 1/2 uniquely (irreducible chance from
   paradox); even length leaves a one-parameter family (x,1−x,…) (a free modulus from consistent loops).
   Unique solution iff −1 ∉ spec(A).
5. Time as evaluation order: ED for acyclic R. Probability: DA (linear relaxation). Complex QM / Born
   rule: UC, and I expect failure. Space, matter: UC.
6. Falsifier: the relaxed solution set is a polytope; if all its correlations obey CHSH ≤ 2 the model
   is a hidden-variable theory and is dead as a route to QM.

## H5 Perfect actuality — "every possible world is the unique minimal corruption of one actual world"
1. Primitives: an alphabet Q, words Qⁿ ("possible worlds"), counting of differing places.
2. Rule: **the balls of radius e around the actual set C ⊂ Qⁿ tile Qⁿ exactly.**
3. Why: alternatives: packings (some possible worlds belong to no actual one) or coverings (ambiguous).
   Exact tiling is the only case where "possible" and "actual" determine each other with no remainder.
4. Forced (ED): |C|·Σ_{i≤e} C(n,i)(q−1)^i = qⁿ; for q=2, e=1: n = 2^r − 1. Classification (known
   theorem, from memory): only trivial codes, Hamming codes and the two Golay codes.
5. Law of physics = decoding map; time = decoding steps (HA). Everything else UC.
6. Falsifier: automorphism groups are finite (GL(r,2), M₂₃); if no scaling limit of Hamming codes
   carries a continuous Lorentz-type symmetry, there is no macroscopic limit. Any numerical match of
   7, 23, 24 with physics would be numerology.

## Computation 1 (H1) — `autology.py`, 16 s
Class tested: every bilinear form B with entries in {−1,0,1} on Rⁿ, n = 1,2,3 (exhaustive: 2, 80,
19682 forms); entries {0,1} for n = 4 (exhaustive, 65535); random samples for n = 4 ({−1,0,1}, 60000)
and n = 5 (20000). For each: stabiliser algebra, test dim g = n, then solve for intertwiners V → g and
test invertibility. Floating-point SVD with tolerance 1e-8, not exact rational arithmetic.

| n | forms | dim g = n | V ≅ g | which |
|---|---|---|---|---|
| 1 | 2 | 0 | 0 | — |
| 2 | 80 | 8 | **8** | all rank-1 symmetric (degenerate), g = aff(1) |
| 3 | 19682 | 696 | **516** | exactly the nondegenerate symmetric forms: 2 definite, 514 of signature (2,1)/(1,2) |
| 3 | | 180 | 0 | nonsymmetric rank-3 forms with g ≅ sl(2,R): dimension matches, module does not |
| 4 exh. {0,1} | 65535 | 4732 | 0 | |
| 4 sample | 60000 | 1474 | 0 | |
| 5 sample | 20000 | 50 | 0 | |

Named controls (Part B): so(p,q) for n ≤ 7 — only so(3), so(2,1) pass; sp(2m,R), sl, gl, g₂ fail on
dimension; **u(2) and u(1,1) on R⁴ and gl(1) on R match in dimension but have Hom_g(V,g) = 0** (clean
negative controls: dimension counting alone is not what selects n=3); so(3,C) on R⁶ passes with
Killing signature (3,0,3), the signature of so(3,1); so(2,C), so(4,C) fail.
Part C (different class: structure = the bracket alone, criterion "smallest nontrivial irrep is the
adjoint"), exact Weyl formula for all 31 simple types of rank ≤ 8: only **E8** (248 = 248).

Verdict: within single-bilinear-form structures the rule singles out dimension 3 with a nondegenerate
symmetric form — plus one unexpected degenerate solution in dimension 2. It does not select the
signature. n = 4,5 are only partly covered (restricted entries / sampling), so "no solutions for n ≥ 4"
is a numerical observation, not a theorem. The E8 line is a separate criterion and must not be read as
a physics match (E8 has no complex representations, so no chiral families directly).

## Computation 2 (H2) — `parallax.py`, 1.2 s
Exhaustive over all reduced Latin squares of order 2..6 (1, 1, 4, 56, 9408; counts agree with the known
values), all transversals, all orthogonal mates, largest mutually orthogonal set.

| n | reduced squares | with a mate | max projections k | n+1 |
|---|---|---|---|---|
| 2 | 1 | 0 | 3 | 3 |
| 3 | 1 | 1 | 4 | 4 |
| 4 | 4 | 1 | 5 | 5 |
| 5 | 56 | 6 | 6 | 6 |
| 6 | 9408 | **0** | **3** | 7 |

Which squares reach the maximum: n=4 only the Z₂×Z₂ table (the three Z₄ tables have 0 transversals);
n=5 only the 6 tables of Z₅ (the 50 non-associative loops have 3 transversals each and no mate).
n=6: 7308 squares have transversals, none decomposes into six — the negative control fails as it should.

Verdict: the rule does select something — it rejects n = 6 and, for n ≤ 5, rejects every coordinate
structure except the additive group of a field. But this reproduces classical theorems (Tarry; uniqueness
of small affine planes); it does not pick one n, and field-emergence is known to break at n = 9.
The parallel "n = 6 admits only 3 projections / only 3 mutually unbiased bases are known in C⁶" is an
analogy, not a derivation.

## Ranking
H1 is the only rule that reaches a physical target (dimension 3, and so(3,1) under one named
assumption) from a class-wide exhaustive check, at the price of an extra nondegeneracy axiom.
H2 gives genuine emergence of a number field from non-numerical primitives but no physics yet.
H3 is the cheapest source of complex amplitudes and probably classical. H4 and H5 are mainly
ontological provocations; I expect H4's falsifier to succeed.
