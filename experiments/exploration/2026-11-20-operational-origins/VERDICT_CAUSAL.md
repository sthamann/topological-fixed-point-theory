# Operational dynamics and a Lorentz response cone

## Principle

Start with the simplest nonclassical homogeneous, self-dual statistical cone of rank two, without fixing its dimension. Its normalized state space is a Euclidean ball of dimension d. Require isotropy under all orientation-preserving rotations and require that differences of observables themselves label infinitesimal reversible generators through a nonzero alternating bilinear bracket. Require covariance of that bracket under all rotations.

This does not initially contain P1, P2, four marks, five sites, E8, or a spacetime dimension. It does contain substantial operational assumptions. In particular, the dynamics correspondence is an assumption, not a consequence of positivity alone.

## Exact result and baseline

For d >= 2, the relevant space is Hom_SO(d)(Λ² R^d, R^d). Its dimension is 1 for d=3 and 0 for every other d. The script explicitly checks d=2,...,8 with rational/integer linear algebra.

| d | alternating coefficients | remaining after π rotations | invariant dimension |
|---|---:|---:|---:|
|2|2|0|0|
|3|9|3|1|
|4|24|0|0|
|5|50|0|0|
|6|90|0|0|
|7|147|0|0|
|8|224|0|0|

The scan is a structural baseline, not a chance probability or a look-elsewhere corrected significance. The comparison grammar is the same rank-two spin-factor cone plus the same equivariant bracket requirement, with only d varied. Every d has a Lorentz-type cone of signature d+1; only d=3 has the specified bracket. Thus the discrimination is caused by the extra dynamics assumption, not by fitting a numerical observable.

A proof beyond the finite scan is short. Under a π rotation in an arbitrary coordinate plane, an invariant coefficient c_ij^k acquires the product of the three coordinate signs. Its index multiplicities must therefore have the same parity in every coordinate. The total index multiplicity is three, so all parities cannot vanish; all must be odd. This requires d<=3. For d=2 the sum of two odd counts is even and cannot be three. At d=3 the surviving tensors have three distinct indices; invariance under infinitesimal rotations leaves one coefficient, the alternating epsilon tensor. At d=1 an alternating bracket is zero. This establishes the result for all d, not merely d<=8.

The surviving bracket is κ times the cross product. It satisfies Jacobi. Together with the spin-factor Jordan product it has the Pauli algebra realization on Herm_2(C). The scalar κ remains free; it fixes the conversion between an observable and a time generator. Its sign is an orientation choice. If the symmetry demand is O(d), including reflections, rather than SO(d), even the d=3 nonzero bracket is excluded. Equivalently, a nonzero bracket itself introduces orientation.

The exact determinant is

    det(t I + x σx + y σy + z σz) = t² - x² - y² - z².

Positive semidefinite matrices are the forward Lorentz response cone. A rational SL2(C) element diag(2,1/2) produces the exact transformation

    (t,x,y,z) -> ((17t+15z)/8, x, y, (15t+17z)/8).

The rational sample only verifies the identity; its coefficients are not TFPT predictions or selected fits.

## What follows, and what does not

The result supplies a conditional origin of a 3+1 signature for a space of possible responses. It does not supply physical events, spatial locality, a Lorentzian manifold, a continuum field integral, Newton's constant, area entropy, P1, P2, or a selected state.

Congruences H -> B H B† for nonunitary B are not deterministic trace-preserving quantum evolution. They are positive filters before normalization. A Lorentz-shaped state cone is consequently not sufficient to identify physical Lorentz boosts.

Real, complex and quaternionic rank-two matrix cones respectively give 2+1, 3+1 and 5+1 signatures; the rank-two octonionic cone gives 9+1. Positivity alone does not select the complex case. The simple local-tomography counts for two ordinary real or complex bits are 3²=9 versus dim Herm_4(R)=10 and 4²=16 versus dim Herm_4(C)=16. The quaternionic numbers 6²=36 versus dim Herm_4(H)=28 are only a naive same-field dimension comparison, not a claimed valid composition law. General composition requires stronger hypotheses.

## Sharp finite-dimensional obstruction to modular boosts

Let K and P be finite-dimensional Hermitian matrices, c a nonzero real scalar, and use the full Hilbert-Schmidt norm. Cyclicity gives

    Tr(P i[K,P]) = 0,

hence the exact identity

    ||i[K,P] - cP||_HS² = ||i[K,P]||_HS² + c² ||P||_HS².

Therefore an exact dilation relation i[K,P]=cP forces P=0. Moreover, the relative residual is at least one:

    ||i[K,P] - cP||_HS / (|c| ||P||_HS) >= 1.

This is a theorem for every finite matrix size. It is checked symbolically for arbitrary Hermitian 2x2 matrices in the script. It forbids a finite-clock construction from approximating the generator of a nontrivial boost/translation relation globally in this norm, even as matrix size increases. It does not forbid approximation in a selected low-energy sector, a different weak topology, unbounded infinite-dimensional limits, nonunitary finite Lorentz representations, or non-Hermitian generators. Those changes reopen the route and must be stated.

For a faithful finite density matrix rho, modular flow is sigma_s(A)=rho^(is) A rho^(-is). On rho=diag(p,1-p), E01 acquires the phase (p/(1-p))^(is). The state parameter p and the conversion of modular parameter to physical time remain inputs. Requiring rho to be invariant under every unitary selects rho=I/n, for which modular flow is trivial. A maximally symmetric state plus a finite algebra therefore does not generate nontrivial modular time by itself.

This gives a concrete next target: construct a sequence of growing local algebras together with a state and inclusion structure, and test the boost relation on an explicitly specified observable sector. A list of finite internal clocks is not this structure.

## Verdict

KANDIDAT, Stufe 1, for the conditional response-cone dimension selector only. The finite-dimensional identification of the same clock with an exact physical boost is VERWORFEN under the Hermitian/full-Hilbert-Schmidt assumptions. A continuum limit with physical spacetime remains OFFEN.

First free structure: a rank-two homogeneous self-dual cone, followed by full SO(d) isotropy and an equivariant alternating dynamics bracket.

Free continuous choices remaining: κ, state rho, physical time scale, spacetime metric scale if a physical embedding is later supplied. Orientation sign is also unselected.

Strongest objection: the bracket requirement could merely encode the special algebraic feature of an ordinary qubit. Until it follows from a deeper compositional or variational rule, the result relocates an assumption. It nevertheless makes that assumption exact and falsifiable over the full family of spin factors.

Novelty: the operational reconstruction of complex quantum theory and connections between bit geometry and spatial dimension are established research areas. This calculation is not a claim of a new mathematical theorem. Its TFPT-specific novelty must be judged against the repository/Atlas. The explicit full-norm boost obstruction is an audit criterion, not a TOE.

## Primary sources

1. Howard Barnum and Alexander Wilce, Local tomography and the Jordan structure of quantum theory, arXiv:1202.4513, https://arxiv.org/abs/1202.4513 . The theorem includes Jordan systems, suitable locally tomographic composites and an existing complex qubit; it does not prove the complex qubit from nothing.
2. Markus P. Müller and Lluis Masanes, Three-dimensionality of space and the quantum bit: an information-theoretic approach, arXiv:1206.0630, https://arxiv.org/abs/1206.0630 . Their theorem assumes spatial dimension d and minimal direction-information systems with continuous reversible interactions. It is a conditional reconstruction, not emergence of an entire spacetime manifold.
3. Erik M. Alfsen and Frederic W. Shultz, Orientation in operator algebras, PNAS 95 (1998), https://pmc.ncbi.nlm.nih.gov/articles/PMC22571/ . The order/Jordan structure does not fix an associative product without orientation; opposite algebras have the same order.
4. Ricardo Correa da Silva and Gandalf Lechner, Modular Structure and Inclusions of Twisted Araki-Woods Algebras, CMP 402 (2023), https://link.springer.com/article/10.1007/s00220-023-04773-y . Discusses half-sided modular inclusions and the Borchers dilation/translation relation, with its operator-algebraic assumptions.

All new equations and bounds stated above are proved directly or checked by causal_origin.py; sources are contextual support, not substituted evidence for a claimed TFPT derivation.

## Repository novelty audit received and checked

The 20 November catalog has prior projects `compiler-cone-object-audit` (Versuchsprojekte row 107) and `compiler-process-reconstruction-20260918` (row 127). The Lorentz-cone dictionary is therefore already known in TFPT and must not be presented as new. Atlasbefunde row 31, `finding:process-composition-and-cyclicity-do-not-select-algebra`, explicitly records that classical and quantum models satisfy generic composition/causality/local tomography, so those premises alone cannot select the algebra. The present selector adds rank-two nonclassical cone structure and a nonzero covariant vector dynamics bracket. No matching SO(d)-equivariant bracket test was located by the catalog audit; absence of a keyword hit is not a proof of novelty. The old Lorentz dictionary should not be counted as a new success, only the conditional selector and the scope-specific finite boost bound are additional checked results.
