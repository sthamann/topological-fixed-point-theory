# Constructive V16 action and audited Mobius clock fold

Research ID: `UR.SEAM.V16_OPERATOR_PULLBACK.01`. Verdict: **PARTIAL**.

The local action inverse and source-preserving fold are proved in their stated scopes. The periodic clock has an explicit exact Weyl conjugacy to the TFPT glue. The actual quarter-phase source is audited separately. Neither construction selects primitive physical locality or supplies raw Xi.

Read [PROOF.md](PROOF.md), [FOLD-AUDIT.md](FOLD-AUDIT.md), and [REVIEW.md](REVIEW.md). The independent reruns passed 32,133 operator checks and 1,668 fold checks.

From the repository root:

```sh
python3 experiments/theory-contracts/v16-operator-pullback-20261014/checker.py --repo . --out /tmp/tfpt-v16-pullback.json
python3 experiments/theory-contracts/v16-operator-pullback-20261014/fold_check.py --out /tmp/tfpt-fold-audit.json
```

The source-charge-matched lift case, 75 original U5 basis/lift cases, 60 mixed Lie actions, residual control, original C/X witness, finite cyclic source links, twisted-sector weights, Virasoro formula, all240 E8 roots, and nonperiodic holonomies are checked. Basis probes do not select physical couplings. Reported cloud26 tests were not executed locally.

No verification, ledger, paper or complete T1-T8 promotion.
