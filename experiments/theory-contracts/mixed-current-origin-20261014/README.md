# Original seam: mixed-current origin audit

14 October 2026. `UR.SEAM.MIXED_CURRENT_ORIGIN.01`. **PARTIAL**.

The original source and the pasted mixed-current proposal are checked against
actual charge, state, conformal weight and field-sector identification.
[PROOF.md](PROOF.md) gives the analytic scope; [results.json](results.json)
records 29 exact checks. [REVIEW.md](REVIEW.md) is the independent review.

```sh
python3 experiments/theory-contracts/mixed-current-origin-20261014/checker.py --repo . --out /tmp/tfpt-mixed-origin-results.json
```

Requires the repository's existing NumPy and SymPy environment. The checker
uses the original transitive source adapter and does not add source copies.
The all-fusions exclusion concerns the explicitly specified one-channel
U(1) Fock module class; target extension and Ward formulas remain conditional.
No full scaling-limit exhaustion, complete TFPT solution, physical gate
closure, ledger change or verification/paper promotion is claimed.
