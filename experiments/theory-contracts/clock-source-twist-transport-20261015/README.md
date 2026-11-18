# Actual J: covariant source twist family, not a fixed-source symmetry

Research ID `UR.CLOCK.SOURCE_TWIST_TRANSPORT.01`. Overall **PARTIAL**.

The actual Gaussian E8 J cannot have an affine charge-lattice lift preserving
the fixed source K: it would require 2 delta in Ls, which fails exactly. The
original fixed Q is not J-invariant either.

The positive result retains the full tuple instead. In the already stipulated
E8 Fock representation, the actual ordered CAR/E8 cocycle admits a unitary
order-four lift U_J. It transports K_delta, Q_q and the full uncompensated
quarter transition to the frames delta_a=J^a delta, q_a=J^a q. This preserves
time laws and states *between the corresponding frames*. Keeping q fixed
would fail. Positivity and the twofold ground degeneracy are transported;
the existing conditional full Cartan Ward rule selects the same vacuum
density without removing the second ground vector.

The four frames are a covariance orbit of supplied source data. They are not
derived physical sectors, a switching law, or a primitive-source selection.
No new Hamiltonian, toy source, physical event instrument or E8 completion
is selected here. The actual lift phase convention is stated in `PROOF.md`.

Run from the repository root:

```sh
python3 experiments/theory-contracts/clock-source-twist-transport-20261015/checker.py
```

The saved run passed **123/123 exact groups**, including source pins, the
order-four character on all eight lattice basis generators, 1,024 cocycle
pair regressions, both quarter cases, four frames and all eight clock
offsets. All-charge statements have analytic proofs; these finite regressions
do not replace them. The checker reads the original definitions without
running unrelated original checkers. Python, NumPy and SymPy are required.

This is an experiment-only theory contract. No physical T1–T8 closure,
verification/ledger/paper promotion, primitive quotient locality, odd-CAR
field dictionary, raw Xi or detector access is claimed. Original hashes are
saved in `source_manifest.json`; detailed checks are in `results.json`.

## Addendum: which responses identify the twist?

`ADDENDUM.md` proves that all vacuum words whose individual insertions preserve all eight Cartan charges are twist-blind, even with complete quadratic tomography. Eight independent original charged primary directions (each ±) reconstruct delta when their normalized return responses and common time scale are supplied. `neutral_twist_checker.py` passed23/23 exact checks, recorded in `neutral_twist_results.json`. Global charge shifts are not derived local physical instruments. This adds no physical T1–T8 closure.
