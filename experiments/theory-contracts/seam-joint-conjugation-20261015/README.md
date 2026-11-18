# Joint quarter conjugation and covariant marking

**Result:** exact transport inside the supplied source and stipulated E8 dictionary;
no complete derivation of that source from P1/P2.

The original quarter-mode particle-hole map fixes the folded action
`O(x)=(x4,-x5,-x6,x1,-x2,-x3,x8,x7)`. It transports charge, vacuum and positive
time between the actual quarter domains. It also transports the carrier/family
marking. No time-preserving orthogonal correction restores the old fixed marking:
the carrier projections have squared norms `129/256` and `137/256`.

Positive completion: transport the original Z4 character with O. The existing
neutral-root formula then reconstructs the transported 5+3 projector without a
new choice. The second marking follows from the first; their physical origin
is not selected by this equivalence.

The proof also corrects two tempting false bridges: complex-bilinear cross-sheet
spinor pairing does not imply a sheet-odd Hermitian Calderon kernel, and a
complex half-spinor 16 is a state/sector representation, not the 16-real-fermion
one-particle coefficient space of the c=8 realization.

See [the derivation](PROOF.txt), [exact checks](results.json) and
[original sources](source_manifest.json). Reproduce with:

```sh
python3 experiments/theory-contracts/seam-joint-conjugation-20261015/checker.py
```

Dependencies: NumPy and SymPy. Checks use exact rational lattice data and small
integer/Gaussian-integer matrices with exact equality. The analytic proofs are
not inferred from the finite counts. No H200 scan is needed for the invariant.

**Firewall:** experiment-only; no ledger, paper, website or T1–T8 promotion.
Verdict `EXACT_COVARIANT_QUARTER_TRANSPORT`; `raw_P1P2_derivation=false`.
