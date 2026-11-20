# Blind agent b2 — four rules, two computations

Worked blind (own mathematics + Python only). Labels: ED = exact derivation, DA = derivation under named extra assumption, NO = numerical observation, HA = heuristic analogy, UC = untested conjecture.

## H1 — "Every syndrome is a generator" (ontology: distinctions with a size)

1. **Primitives.** A finite set of N alternatives; formal combinations with coefficients in a division algebra K; a size function ‖x‖_p. Nothing else.
2. **Rule.** *The number of independent non-trivial outcome statistics equals the number of independent continuous reversible moves that act non-trivially on those statistics*: `dim(g) − dim(g_trivial) = dim span{|(Ux)_i|^p} − 1`, for every N ≥ 2.
   Coding reading: each thing you can learn (syndrome) is exactly one way to move reversibly (recovery), no surplus on either side.
3. **Why.** Class compared: all (K, p) with K ∈ {R, C, H}, p ∈ [1,∞). It is a single-system condition: no composites, no local tomography, no purification (this is the essential difference from Hardy/Chiribella-type reconstructions; it is close to the Alfsen–Shultz "dynamical correspondence", but stated as a count over the (K,p) class with an effectiveness clause).
4. **Forced consequence.** For p ≠ 2 the isometry group is monomial (Lamperti), so effective generators = 0 < N−1. For p = 2: R gives N(N−1)/2 vs N(N+1)/2−1; H gives N(2N+1) vs N(2N−1)−1; C gives N²−1 = N²−1. Only (C, 2) passes. Note the *raw* count g = r also holds for C at every p (N phases vs N statistics) — the phases act trivially, which is why the effectiveness clause is needed.
5. States = rays of C^N: **ED** (within the class). Born rule |x_i|²: **DA** (statistics are the p-th powers of the size function — assumed form). Hilbert tensor product for composites: **UC**. Space, time, matter, gravity: not addressed (**UC**).
6. **Falsifier.** Any (K,p) other than (C,2) satisfying the count for some N ≥ 2; or a physically motivated enlargement of the class (octonions, non-ℓ_p gauges) where a second solution appears.

## H2 — "Chiral and invisible to all probes of order ≤ 3" (ontology: syndromes of n binary checks)

1. **Primitives.** n commuting two-valued labels (checks/logical operators with eigenvalues ±1); a state is a sign string s ∈ {±1}^n; a multiset of states with a global flip s → −s (conjugation).
2. **Rule.** *Matter is the shortest sign code that is not flip-symmetric yet unbiased against every product of at most three checks*: an odd integer function x on {±1}^n with Σ_s x(s) s_i s_j s_k = 0 for all i, j, k (i=j gives the linear condition for free).
3. **Why.** Class: all odd integer-valued functions on {±1}^n for all n (not only linear codes). Order 3 is the triangle-anomaly order; this is an input equivalent to D = 4 (not derived).
4. **Forced consequence.** An odd function has only odd-degree Walsh coefficients; the rule kills degrees 1 and 3; so n ≤ 4 has no solution and at n = 5 the only survivor is the degree-5 character: x = ± parity. The unique solution is the even-weight code [5,4,2]: **16 states** (or its 16 conjugates). In coding language: the dual code has an odd word but none of weight 1 or 3. These 16 sign strings are the weights of the Spin(10) spinor.
5. Matter content "16 chiral states per family": **ED** given binary charges and order 3. Gauge group SU(3)×SU(2)×U(1): **UC** (requires an unexplained 3+2 split of the five checks; only the Cartan charges are derived). Three families: **not obtained**. QM, space, time, gravity: not addressed.
6. **Falsifier.** A smaller chiral solution inside the class (excluded exhaustively below); or the binary restriction being unjustifiable — the control shows it is essential.

## H3 — "Exactly one of a part and its complement is recoverable" (ontology: a recoverability relation on a finite set)

1. **Primitives.** A finite set V of parts and a family Δ of subsets ("erasable/recoverable sets"), closed under taking subsets — which any erasure-correcting code supplies automatically.
2. **Rule.** *For every A ⊆ V exactly one of A, V∖A is in Δ* (Δ equals its Alexander dual), *and Δ is locally spherical* (every link is a sphere: "every recoverable set is surrounded uniformly").
3. **Why.** Class: all down-closed set families (= all possible correctable-erasure patterns of any code, classical or quantum). Self-complementarity is the sharpest no-cloning statement (never both, never neither). Difference from holographic complementary recovery: no bulk/boundary, no tensor network, no entropy; the code's correctable family *is* the space.
4. **Forced consequence (recalled theorem, Brehm–Kühnel / Arnoux–Marin, not re-verified here).** A combinatorial d-manifold satisfying complementarity has |V| = 3d/2 + 3 and d ∈ {0, 2, 4, 8, 16}; examples RP²₆, CP²₉ (unique), and 8-/16-dimensional "projective-plane-like" manifolds. So dimension is not free; d = 4 forces a 9-part object with the homology of CP² — a complex structure appears uninvited.
5. Space with d = 4 as one of four options: **DA** (local sphericity). Lorentzian signature, time: **UC**. Existence of a quantum code on nine qutrits whose correctable sets are exactly the faces of CP²₉: **UC** (concrete, checkable). Matter: χ(CP²)=3 ↔ three families is **HA** only, flagged as numerology. Gravity as variation of Δ under code deformation: **HA**.
6. **Falsifier.** Proof that no code (stabilizer or otherwise) on 6 parts realises RP²₆ or on 9 parts realises CP²₉ as its correctable family.

## H4 — "Five pairings" (ontology: a set with fixed-point-free involutions; pure incidence)

1. **Primitives.** A set X and c fixed-point-free involutions (c perfect matchings).
2. **Rule.** *c is the least number of pairings whose mutual 2-cycles carry a chiral invisible code in the sense of H2*; cells = elements, gluing = pairings (a (c)-coloured graph is a (c−1)-dimensional pseudo-manifold).
3. **Why.** Class: c-coloured graphs for all c. Differs from coloured tensor models by having no sum over graphs and by tying c to the matter code rather than choosing it.
4. **Forced consequence.** The minimal H2 code has n = 5; taking X = {±1}^5 with the five single-flip involutions gives the 5-cube graph, whose encoded complex is the boundary of the 5-dimensional cross-polytope: a 4-sphere with 32 simplices, bipartite into 16 + 16 by parity = orientation. The H2 states are exactly one orientation class. In general n_min = m+2 with m = D/2+1 gives n = D/2+3, and D/2+3 = D+1 only for D = 4.
5. Space dimension 4: **DA** (anomaly order m = D/2+1 imported) and the D = 4 coincidence is one equation tested against four candidates — weak, flagged as possible numerology. Matter = orientation class of cells: **HA**. Gravity as deficit of bicoloured cycle lengths (Regge-like): **UC**. QM: none.
6. **Falsifier.** The selection disappears if anomaly order or alphabet changes (it does for the alphabet, see control).

**Unification reading (role 2).** H1 and H3 share one cause: *recoverability*. QM is the unique single-system arena where learnable = reversible; geometry is the pattern of which parts are recoverable. This link is **HA**, not derived.

## Computation 1 (H1): `c1_field_norm.py` → `c1_output.txt`

Numerical nullspaces (SVD, tolerance 1e-9, random sampling, seed 7) for K ∈ {R,C,H}, N ∈ {2,3,4}, p ∈ {1, 1.5, 2, 3, 4}: 45 cases. g = isometry Lie algebra dimension, t = generators trivial on all statistics, r = dimension of span of statistics.

| K, p=2 | N=2 (eff, r−1) | N=3 | N=4 |
|---|---|---|---|
| R | 1, 2 | 3, 5 | 6, 9 |
| C | 3, 3 | 8, 8 | 15, 15 |
| H | 10, 5 | 21, 14 | 36, 27 |

All p ≠ 2: eff = 0, r−1 = N−1 (R: g=0; C: g=t=N; H: g=t=3N). **Rule satisfied only for (C, p=2), all three N.** The raw count g = r is satisfied by C at *every* p (15 of 15 complex rows) — so without the effectiveness clause the rule does not select the 2-norm. Caveats: finite p-list, N ≤ 4, numerical rank; the result reproduces known theorems (Lamperti; dim u(N) = dim Herm) rather than discovering one; "statistics = |·|^p" is assumed.

## Computation 2 (H2): `c2_chiral_code.py`, `c2b_ternary_control.py` → `c2_output.txt`, `c2b_output.txt`

Exhaustive meet-in-the-middle over all x ∈ {−1,0,1}^(antipodal pairs):

- Binary, cubic rule: n = 1,2,3,4 → 0 solutions (rational nullity 0, so none with any multiplicity); n = 5 → nullity 1, exactly 2 solutions, each 16 states of one parity. n = 6 → nullity 6, exactly 72 solutions (60 with 16 states, 12 with 32).
- Order scan: m=1 → n_min 3 (4 states); m=3 → 5 (16); m=5 → 7 (64). D/2+3 = D+1 only at D=4 among D = 4, 8, 12, 16.
- Controls: alphabets {−1,0,1} (n ≤ 3), {±1,±3}, {±1,±2}, {0,±1,±2} (n ≤ 2): 0 solutions. **But** {−1,0,1} at n = 4 has rational nullity 20; the 15 weights of SU(5) 5̄+10 in Dynkin labels are a verified solution (15 < 16), while exhaustive search finds none with ≤ 8 states. Sizes 9–14 at n=4 were not searched.

**Verdict.** Inside the binary class the rule singles out exactly one object (the 16) at minimal length — exact, not numerology, but it is the known Spin(10) spinor seen as a code. The uniqueness is destroyed by allowing zero charges (15-state solution with fewer checks), so "every check acts on every state" is a load-bearing, unexplained assumption. The order-3 input presupposes D = 4. No gauge group, no three families.

## What was not achieved
No rule here yields Lorentzian signature, three families, SU(3)×SU(2)×U(1), gravity dynamics or any dimensionless constant. H3's key consequence rests on a recalled theorem and an untested code-existence conjecture; that (search for a 6-part code with correctable family RP²₆) is the most informative next computation.
