# Gemeinsame Herkunft der beiden TFPT-Normierungen

Research ID: UR.SEAM.TWO_NORMALIZATIONS.ORIGIN.01  
Date: 2026-10-05  
Verdict: PARTIAL (exact internal geometry; physical response selection open)

Question: Does one original seam operation supply c3 = 1/(8π), phi_base = 1/(6π), the deck action and the actual alpha response, without inserting its target coefficients?

Success criterion: Use the native source operators and their functional to obtain the absolute modular response, the conditional three-channel response, the puncture term 48 c3^4 and the complete U(1) determinant stationarity. Matching dimensions, names or coefficient identities is insufficient.

Decisive exclusions: Do not replace the original source by a new oscillator or guessed E8 Hamiltonian. Stop a direct identification if its group action, operator rank or cumulant coefficient differs. This does not exclude a different, explicitly derived continuum or sector readout of the same source.

## Result

The native 240-root generator supplies the tetrahedral A3 weight frame: T^T T = I4 - J4/4 and TT^T = I3. Its weights have squared norm 3/4; the matching D5 spinor weights have squared norm 5/4. Their sum is the E8 root norm 2. The four-to-three rank ratio and the two glue norm representatives therefore have a common internal source.

With beta_angle = 2π and a normalized response trace, response per direction is 1/(4 beta_angle) or, after conditioning and renormalizing a rank-three response, 1/(3 beta_angle). This explains the ratio 4/3 conditionally. The choice of response functional, absolute normalization and selected physical three-channel subspace remains an origin obligation. The current Gram I4 is a current correlation, not a vacuum probability distribution.

The zero-sum projection commutes with the literal four-mark cyclic permutation. It is not equivariant under the entire SU(4) fundamental action. A coordinate three-family cut and the zero-sum H1 space are distinct embeddings; any physical identification must supply its operator dictionary. No exclusion of a separately derived SU(3) selection is claimed.

In the canonical lattice-VOA reading, h_D = 5/8, h_A = 3/8 and h_glue = 1. The carrier chiral two-point factor has exponent -5/4, the family factor -3/4, and the complete local E8 current has exponent -2. This derives a candidate meaning of the alpha exponent from the native weights; it does not derive carrier-only detection, separation 1-q, q(alpha) = 48 c3^4 exp(-2 alpha), or the full Ward functional.

For those conditional chiral correlators, analytic continuation by one turn gives phases -i and +i, cancelling in the full product. After two turns each factor still has phase -1; four turns restore each. These are CFT analytic monodromies, not the geometric collar deck or physical spacetime spin.

The original v484 contact matrix C has determinant (1-4u)(1+2u)^2. Its unweighted -logdet fourth coefficient is 72, yielding 72 (epsilon c3 ln2)^4 with v484's insertion variable. It is not the target 48 c3^4. Projection onto the naive zero-sum space gives coefficient 8 and rank 2. This excludes only reading the puncture coefficient directly from that bare unweighted determinant; continuum matching, admissible-state weights and a physical response derivative are not tested or disproved here.

v985's finite Fourier-channel swap is reproduced exactly. Multiplication by D = diag(1,-1,1,-1) maps the jump matrix to |C| and conjugates the mark rotation r to -r. The permutation half-turn S = r^2 leaves the jump matrix unchanged. D commutes with S but D != S. Calling both operations Z2 does not supply their identification. A geometric deck with this fibre lift is possible only with an explicit fibre action. The established finite v985 result is confirmed, not reopened.

## Scope and sources

checker.py extracts the unchanged e8_roots function from v1, the actual contact matrix assignment from v484 and the jump pattern from v985. It adds no replacement kernel, Hamiltonian or physical source state. Checks remain active under Python -OO. Source hashes in results.json describe this replay only; they are not external evidence.

Original sources: verification/v1_e8_glue.py; v154_simple_current_theorem.py; v341_alpha_quillen.py; v470_alpha_inflow_level.py; v484_seam_contact_unit.py; v485_contact_diagonal_closed.py; v813_p1_index_kms.py; v985_quillen_channel_swap.py. Original affine-current input: compiler-vacuum-current-cubic-20260918/PROOF.txt, UR.COMPILER.VACUUM_CURRENT_CUBIC.06. The PBW checker from that earlier contract was not re-run here; this probe does not claim a new full current-algebra proof.

Related ledger obligations: AX.P1.01; P1.INDEX.KMS.01; HYP.PHI0.PUNCTURE.01; ALPHA.QUILLEN.EXACT.01; SEAM.DETLINE.UNIFICATION.01; QGEO.KERNEL.01. Their status is not changed.

Literature: Lam–Yamada–Yamauchi, arXiv:math/0403010 (E8 index-four A3+D5 sublattice; its separate sqrt(2) VOA construction is not imported); Zuber, A short introduction to conformal field theory, equation (80) (primary two-point exponent); Longo–Tanimoto, arXiv:1608.08903 (rotational KMS states with specified generator and beta, not a derivation of the TFPT source).

## Replay

python3 -B experiments/theory-contracts/seam-two-normalizations-origin-20261005/checker.py

The same assembled alpha closure gives inverse alpha = 137.0359992168407125 for exponent 5/4; substituting the complete-glue exponent 2 gives 137.0360022863583280. This ablation diagnoses a meaningful readout selection. Neither number is a new physical derivation.

## Next origin step

Derive the original charged puncture operator and source functional, including its measure and fibre/deck lift. Evaluate that single functional on the family-conditioned response and the carrier correlator, and establish the finite-to-continuum determinant identification. Do not tune the bare determinant to 48 or insert the desired response before evaluating it.

The earlier compiler-vacuum-current contract locates the upstream obligation: raw seam state plus actual local charged operators -> the marked affine E8 level-one vacuum current algebra. This work gives an explicit common weight frame and response target for that obligation, not its completion.

## Firewall

Research-only theory contract. No empirical scorecard, verification module, paper, ledger, compiler/runtime implementation or physical gate is promoted. Existing experiments/next.txt has other work in progress and is left intact; this standalone contract with its own index carries the next research step. The live theory graph is stale because of these new files and another changing active contract; it is not regenerated here.
