# Independent audit of the QEC syndrome bridge

## Finite group chain

The implementation in `syndrome_spread.py` correctly enumerates GL(4,2) and the spread stabilizer, and checks the derived series by explicit permutation multiplication. The orders 20160 -> 360 -> 180 -> 60 and faithful A5 image are sound within the syndrome-incidence model.

An independent conceptual identification is F2^4 = F4^2. A Desarguesian spread consists of the five one-dimensional F4 subspaces. Its semilinear automorphism group is GammaL(2,4), order 2(16-1)(16-4)=360. The quotient on lines is PGammaL(2,4)=S5 and the pointwise line kernel is F4^*=C3. The canonical perfect core is SL(2,4)=A5. This is an abstract syndrome symmetry, not automatically its physical logical action.

The independent exact matrix test `qec_bridge_audit.py` constructs the code projector as a 32x32 rational matrix and checks every one of the 120 physical tensor permutations including stabilizer phases. Exactly ten preserve it. Thus the pure physical permutation group is D5, whereas S5 is the larger syndrome-incidence action.

The additional test `qec_lift_check.py` constructs the physical symplectic lift of every one of the 360 syndrome automorphisms. Each has a unique site permutation plus local one-qubit symplectic transformations; all preserve the stabilizer subspace. The full induced logical symplectic image has order 6, Sp(2,2)=S3. The entire perfect A5 core acts as the identity on logical Pauli labels.

This does not mean the full physical action is trivial. It means that the derived A5 cannot simply be identified with an icosahedral rotation action on the encoded logical qubit. Any consistent faithful projective A5 action on a qubit has image A5 in PU(2), while logical Clifford actions have image contained in S4. No injective A5 -> S4 exists. Pauli correction choices and phase lifts must be distinguished from this verified symplectic statement.

## Central extension and E8

For the perfect group A5, its universal central extension is canonically 2.A5, isomorphic to SL(2,5), once taking the universal central extension is the declared construction. The two minimal faithful complex irreducible representations have dimension two and are exchanged by an outer automorphism/Galois conjugation. Their McKay graphs are isomorphic affine E8. Thus choosing between those two does not change the abstract graph.

This can be stated as a mathematical association from the syndrome geometry, conditional on using universal projective completion followed by the minimal representation McKay construction. It is not yet a physical origin theorem. Neither the perfect-core extraction nor the central extension alone proves that source processes carry that projective symmetry, that its central involution is TFPT sheet parity, or that the McKay graph is the physical E8 current algebra. In particular, the tested logical Clifford action does not supply the needed faithful qubit representation.

The first extra bridge choice is the decision that the universal projective lift of an abstract syndrome-incidence symmetry has the physical role of a source symmetry. The minimal representation/McKay operation is another declared mathematical construction, even when its output has no continuous fit parameters. Therefore this is OFFEN as a source bridge and at most a conditional Stufe-1 structural correspondence, never Stufe 2 or 3.

## Why 63 is not a derived memory count

For an isometric encoding V:C²->(C²)^tensor5 and a three-site set A, consider

    F_A(X) = Tr_(A complement)(V X V†),  X in M2.

The domain has complex dimension four, so rank(F_A)<=4. Because every three-share set exactly recovers the logical qubit, F_A is injective and rank(F_A)=4. Its Hermitian trace-zero image has real dimension three. This is established exactly for all ten triples by `qec_bridge_audit.py`.

For pure logical input, every such reduced density matrix has spectrum 1/4 (multiplicity four), 0 (multiplicity four). The maximally mixed logical state maps to I8/8. Equivalently, up to a decoder unitary the family is rho_logical tensor I4/4.

The full three-qubit matrix algebra has dimension 64 and 63 traceless directions. But the encoded family occupies only four/three of those directions. Consequently 63=4^3-1 is a full-subsystem dimension coincidence, not a count of the physically reachable memory directions of the encoded qubit. The TFPT process count 68=5+63 is not reproduced by this code.

Arbitrary physical errors, coherent syndrome excitations, measurement records or additional dynamics can enlarge the reachable state/process space. That is an explicit change of the model; it needs its own reachability and minimal-process-rank calculation. The rank-four result only concerns the fixed encoded one-qubit operator family and must not be generalized to all protocols.

## Exact outputs

- phase-sensitive physical tensor permutation count: 10
- spread symplectic lifts: 360
- full logical symplectic image: 6
- A5 perfect-core logical symplectic image: identity only
- all ten three-share operator ranks: 4
- all ten three-share traceless ranks: 3
- pure logical three-share spectra: {1/4:4, 0:4}

## Useful primary references for the extension/McKay qualification

- Quantum codes and irreducible products of characters, Designs, Codes and Cryptography (2025), https://link.springer.com/article/10.1007/s10623-025-01599-8 . Discusses the two degree-two Galois-conjugate irreps of SL(2,5) and explicitly notes their equivalence for choosing a fundamental SU(2) embedding. It does not establish this five-qubit syndrome bridge.
- Pierre-Philippe Dechant, The birth of E8 out of the spinors of the icosahedron, Proc. R. Soc. A (2016), https://pmc.ncbi.nlm.nih.gov/articles/PMC4786034/ . Explains the binary icosahedral McKay graph and the separate spinor/icosian constructions.

These standard group correspondences support the conditional association; the physical obstruction and the rank test above are proved directly.
