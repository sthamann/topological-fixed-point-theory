# Blind agent b3 — four ontologies, two decisive computations

Roles: topology/geometry, time/causality. Worked blind (own mathematics + Python only).
Labels: ED = exact derivation, DA = derivation under named extra assumption, NO = numerical observation, HA = heuristic analogy, UC = untested conjecture.

---

## H1. "Twice" — what exists is composable doing; the only law counts returns

1. **Primitives.** A set A of "acts" that can be superposed (real vector space) and composed (associative bilinear product). No metric, no involution, no order.
2. **Rule.** Q(x) = tr(x∘x): the return count of doing x twice (trace in the regular representation). The sign of Q *is* the geometry. A must be simple (no mutually independent sectors).
3. **Why.** tr(x) is linear (no geometry); tr(x∘x) is the lowest-order return invariant that is a quadratic form and needs no extra structure (no transpose, no star). Compared class: all finite-dimensional real associative algebras; controls: non-simple, non-associative, another quadratic invariant, group algebras.
4. **Forced consequences.** By Wedderburn, simple A = M_n(K), K = R, C, H. Symmetric/antisymmetric (hermitian/antihermitian) parts are Q-orthogonal with opposite signs, giving
   M_n(R): (n(n+1)/2, n(n−1)/2); M_n(C): (n², n²); M_n(H): (n(2n−1), n(2n+1)).
   Exactly one direction of the minority sign ⇔ A ∈ {C (1,1), M_2(R) (3,1), H (1,3)}. Hence *simple + noncommutative + a single odd direction ⇒ dimension 4, Lorentzian* (ED, verified below). For traceless x in M_2(R), Q(x) = 0 ⇔ x² = 0: "light = acts whose repetition is nothing". The odd direction in M_2(R) is J with J² = −1, a cyclic clock of order 4.
5. **Emergence.**
   - Time/space: 4 dimensions, one odd sign — ED within the class. Which real form (clock-like time J in M_2(R), or identity-like time in H with Aut(H) = SO(3) as rotations) is undecided.
   - Lorentz group: O(3,1) is the isometry group of Q, but only PGL_2(R) ≅ SO(2,1) (resp. SO(3)) preserves the product. Full Lorentz symmetry needs the product to be forgotten at large scale — UC.
   - QM: both real forms complexify to M_2(C); complex amplitudes as their common envelope — HA.
   - Gauge group, matter, gravity: nothing derived — UC.
6. **Falsifier.** An equally structure-free quadratic return invariant with another index, or a proof that a non-simple class is equally minimal (control B shows simplicity carries real weight).

## H2. "Three clocks" — what exists is finite cyclic clocks; the only law is joint return

1. **Primitives.** k cyclic acts x_1..x_k of pairwise coprime finite orders p_i (no shared sub-rhythm). No space, no common time.
2. **Rule.** x_1 x_2 ⋯ x_k = 1 (running all clocks in turn returns), and the generated group is finite and nontrivial (the world closes).
3. **Why.** One relation is the minimum that couples the clocks at all. Compared class: all k and all order tuples; controls: non-coprime tuples, k = 2, k ≥ 4.
4. **Forced consequences.** k = 2: x^p = x^q = 1 ⇒ trivial. k = 3: von Dyck group, finite iff 1/p+1/q+1/r > 1; with pairwise coprimality only (2,3,5): order 60, ≅ A5. k ≥ 4: Σ(1−1/p_i) ≥ 2 ⇒ infinite. So the number of clocks (3) and their orders (2,3,5) are both outputs (ED; theorem, checked below). A5 is perfect; its minimal faithful representation is real of degree 3; its double cover 2I lies in the unit quaternions S³; S³/2I is a homology 3-sphere; the McKay graph of 2I is affine E8.
5. **Emergence.**
   - Space: three incompatible cyclic orders whose synchronisation group first acts faithfully in real dimension 3; compact space S³/2I — DA (assumption: space = minimal faithful real carrier of the synchronisation group).
   - Time: none fundamental; observed time = a section order on the group. No Lorentz signature — UC.
   - QM: spinorial double cover (2I irreps with Frobenius–Schur −1) — HA.
   - Gauge/matter: E8 via McKay contains SU(3)×SU(2)×U(1), but E8 has only real representations, so chirality is not delivered — HA with a known obstruction.
   - Gravity: UC.
6. **Falsifier.** A coprime tuple other than (2,3,5) with finite nontrivial group (none exists). Empirically, the rule demands dodecahedral-space topology, which would have to show up in large-angle sky correlations.

## H3. "Solved both ways" — what exists is one unordered constraint; time is an elimination order

1. **Primitives.** A finite set V of variables over a field and one sparse linear relation system A·φ = 0. No order on V.
2. **Rule.** Only boundary data fix φ (past and future together). A *time* is any total preorder of V in which A is block-triangular with invertible diagonal blocks (forward substitution). Different admissible orders are incompatible times; observed time is their common coarsening.
3. **Why.** Triangularisability is the one order-theoretic property of a relation invariant under relabelling. Compared class: all sparsity patterns; "elliptic" ones admit no such order, "hyperbolic" ones a cone of them.
4. **Forced consequences.** For translation-invariant second-order A on Z^d with symbol a quadratic form q, forward substitution along a covector n works iff q is definite on a complement with sign opposite to q(n), i.e. iff q has index 1 (Gårding hyperbolicity). Existence of *any* time forces Lorentzian signature; admissible times form the timelike cone; relativity of simultaneity = non-uniqueness of elimination order (DA: translation invariance, second order).
5. **Emergence.** Signature index 1 — DA. Dimension 3 — UC. QM as a sum over boundary-compatible solutions with characters as amplitudes — UC. Closed causal loops = irreducible non-triangular blocks, allowed, perceived as "events" — HA. Gravity = variable sparsity pattern — UC.
6. **Falsifier.** If generic sparse constraint systems have a common coarsening with two or more independent ordered directions, "one time" fails. Differs from causal sets (order is derived and plural, not primitive) and Page–Wootters (no Hilbert space or clock subsystem presupposed).

## H4. "Two pairings" — what exists is perfect matchings; the only law is loop counting

1. **Primitives.** 2n points and perfect matchings of them (nothing is alone).
2. **Rule.** ⟨a|b⟩ = d^{loops(a∪b)}: overlay two matchings, each closed return contributes a factor d; the form on formal sums must be positive semidefinite for all n.
3. **Why.** A matching is the smallest self-returning structure (fixed-point-free involution); overlaying two is the smallest comparison. Compared class: planar (order-respecting) versus arbitrary matchings.
4. **Forced consequences.** Planar (Temperley–Lieb): positivity for all n iff d ≥ 2 or d = 2cos(π/k) (Jones). Arbitrary (Brauer): positivity for all n forces d a positive integer. So a linear order of the points (planarity as proto-time) permits a quantised family of loop values; no order forces a classical integer dimension. Known theorems; the hypothesis is their physical reading.
5. **Emergence.** Hilbert space from positivity — ED given the rule. Born rule as two-sided overlay, hence quadratic — HA. Complex numbers from loop orientation — UC. Time = planarity order — HA. Space, gauge group, gravity — UC.
6. **Falsifier.** If no d yields a tensor-factorising family of spaces containing a qubit with full SU(2), the rule cannot reach standard QM.

---

## Computation 1 (H1): `sig_enum.py` → `sig_enum.out`

Exhaustive over simple algebras M_n(K), K ∈ {R,C,H}, n ≤ 8 (24 algebras, dim up to 256). Signature of tr(x²) from explicit bases and, for dim ≤ 64, independently from the regular representation (all 17 checks OK).

- Index-1 simple algebras: **M_2(R) (3,1), C (1,1), H (1,3)**; none above dimension 4. Matches the closed forms.
- Control B (drop simplicity; all semisimple algebras of dim ≤ 8): index 1 occurs in *every* dimension 2..8 (R^k × C gives (k+1,1)); in dim 4 also R×R×C (3,1).
- Control C (drop associativity; Cayley–Dickson): (1,1), (1,3), (1,7), (1,15) — no selection.
- Control D (other invariant e₂ on M_n(R)): (2,2), (4,5), (7,9), (11,14), (16,20) — never index 1. The choice of form is an input.
- Control E (group algebras, 14 groups of order ≤ 8): index 1 for Z3 (2,1), Z4 (3,1), S3 (5,1), D4 (7,1) — no selection.

**Verdict.** The rule singles out dimension 4 with one odd direction only inside "simple, associative, noncommutative, form tr(x²)". It does not choose between (3,1) and (1,3), and selects nothing in the control classes. It is a theorem about quaternion algebras rather than a numerical coincidence, but the selecting power sits in the word "simple".

## Computation 2 (H2): `clocks_tc.py` → `clocks_tc.out`

Own Todd–Coxeter coset enumeration (cap 30000 cosets; sanity Z5 → 5, S3 → 6, (2,3,7) → cap).

- All 120 triples 2 ≤ p ≤ q ≤ r ≤ 9: finite exactly for (2,2,n), (2,3,3), (2,3,4), (2,3,5); 0 mismatches with the spherical criterion and order formula. Of 18 pairwise-coprime triples **only (2,3,5), order 60**, is finite.
- Controls: k = 2 coprime → always order 1. Five coprime 4-tuples, one 5-tuple and (2,2,2,2): cap exceeded (enumeration cannot prove infinity; consistent with the Fuchsian/Euclidean theorem). Non-coprime k = 3 gives ten other finite groups, so coprimality does the selecting.
- Forced group built as quaternions: 2I has 120 elements, 9 classes, is perfect, irreps 1,2,2,3,3,4,4,5,6; McKay graph is a tree with arms (1,2,5) = affine E8 and (2−A)·dims = 0. Quotient: 60 elements, perfect, irreps 1,3,3,4,5, all real; 120 pairs of orders (2,3) with product of order 5, each generating all 60.

**Verdict.** The rule singles out the number of clocks (3), their orders (2,3,5), a group (A5/2I) and a carrier dimension (3). Caveats: finiteness is an extra axiom; "3" is partly built in, because finite three-clock groups are exactly rotation groups of the 2-sphere; no time direction, signature or chiral matter follows.

## Overall assessment

Neither computed rule reaches the targets. H1 gives 3+1 as an index but has no dynamics; H2 gives a compact 3-dimensional space and E8 but no time. They are complementary (H1's odd direction is itself an order-4 clock), yet no derivation joins them. H3 gives signature-from-solvability under assumptions; H4 gives Hilbert-space positivity. No hypothesis yields SU(3)×SU(2)×U(1) with three chiral families.
