# Native sector calibration

The existing conditional v1033 source retains a four-sector CAR algebra
and its original disorder operator D. The internal algebra leaves four
energy offsets undetermined; the D response fixes their three relative
modes, leaving one common scalar. Four native D-translated preparations
provide an explicit reconstruction protocol. D² alone loses the parity
component of this calibration information.

Run `python3 checker.py` for the 256-dimensional source check; read
`PROOF.txt` for the all-parameter proof and its preparation/readout scope.
`source_manifest.json` identifies the original inputs.

Verdict: **PARTIAL**. This is a project-specific application of standard
finite-dimensional algebra, not a new physical model or general theorem.
The Fock lift, preparation and physical time are still conditional.
Absolute geometric normalization, P1/P2 selection and T1–T8 remain open.
Experiment only; no empirical scorecard or verification-ledger promotion.
