# Native charge conjugation and the charged response test

Research ID `UR.NATIVE.CHARGED_CONJUGATION.01`. Overall **PARTIAL**.

The original full Round20 onsite field source has a linear charge conjugation
C with C T_p C=T_-p and C W_p C=W_p^*. Its actual ordered cocycle signs are
retained: T_p^*=(-1)^e(p)T_-p. The original quadratic/scalar/hopping
Hamiltonian is C invariant, as are its full and total-zero-sector Gibbs starts.
Opposite-charge adjoint-paired responses are therefore exactly equal for
all times. This excludes equality of all such *marked* responses with a
single fixed nonzero chiral quarter twist at the specified starts.

The full field net exists before the later C0 restriction. Balanced bonds
stay in C0; onsite fields leave it and return through charged corners. A
uniform eta.Q_tot is zero on C0 but remains visible in those insertion
responses, with exact +/- ratio exp(-2t eta.p). Adding or selecting eta is
not derived here. A symmetric twist mixture is not inserted as a solution.

Read `PROOF.md` for the all-charge/all-time argument and the precise state,
marking and source hypotheses. Run:

```sh
python3 experiments/theory-contracts/native-charged-conjugation-20261015/checker.py
```

The saved run passed **24/24 exact groups**, zero failures. `results.json`
and `source_manifest.json` record the checks and original source hashes.
No numerical heat-kernel/Gibbs simulation, primitive P1/P2 selection,
local field identification, unrestricted TFPT no-go or physical gate
closure is claimed. This is an experiment-only theory contract.
