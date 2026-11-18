# Collar source and multi-field sewing

Research contract `UR.SEAM.COLLAR_MULTIFIELD.01`, 2026-10-16.

Run with the positive-kernel search dependencies:

```sh
python collar_control.py --search-controls --out results/collar-control.json
python charged_pair_control.py --search-controls --out results/charged-pair-control.json
python source_sewing.py --out results/source-sewing.json
python chiral_tail_control.py --out results/chiral-tail-control.json
python joint_source.py --out results/joint-source.json
python kernel_functional.py --completion ../seam-positive-kernel-search-20261015/runs/source-clustering/min-trace.json --out results/kernel-functional-min-trace.json
python kernel_functional.py --completion ../seam-positive-kernel-search-20261015/runs/source-clustering/max-trace.json --out results/kernel-functional-max-trace.json
python e8_boundary_response.py --output results/e8-boundary-response.json
python -m unittest discover -p 'test_*.py'
```

The complete argument, source premises and physical limits are in `PROOF.txt`.
`source_sewing.py` reconstructs the exterior lift from full field covariance.
`boundary_amplitude.py` independently reconstructs the same operator ray from
an actual two-face CAR projection state and a fixed oriented state/map pairing.
An unnormalized scalar partition function is not obtained from normalized data.
`collar_control.py` tests the known quarter-domain collar and separately tests
the actual thermal GNS transfer; it does not substitute a vacuum contraction.
Its equal-mode control explicitly duplicates the lowest source mode. The
distinct `charged_pair_control.py` uses the actual lowest particle/hole energies
1/4 and 3/4, giving response factors 1/4 and 1/64 and joint response 1/256.
Its exact reconstruction proves feasibility independently of the SCS
`optimal_inaccurate` status, which is retained in the result.
In the v2 collar report the optional graph-parent matrices are explicitly
separated from the source amplitude; the input-dual and output faces are labeled.
Section 8 of the proof supplies explicit operator-norm and trace-tail bounds for
the complete known chiral quarter source in its declared sea-relative CAR
representation. A finite cutoff is not used as a continuum proof.

The existing unknown-kernel search now defaults to feasibility. Its contraction
and thermal transfer classes are explicit; a pure-state completion failure is
inconclusive. The 2026-10-15 min/max-trace runs remain historical calibration.
It now also accepts an explicit reference eigenvalue lambda, without confusing
conditional state normalization with an absolute source partition function.

Section10 and `joint_source.py` integrate a supplied quadratic fermion source
jointly: the source covariance, effective retained-mode action and determinant
factor come from the same Berezin integration. The original quarter circle
operator has zeta determinant1-i and an explicitly checked Green/eta response.
These are conditional quantizations of given operators. They do not identify
the archived comparison part D_rel with the full physical kinetic operator.
An arbitrary spectral profile is not silently converted into -logdet; the
necessary profile variation, measure, reference and phase conventions are
checked separately. The concrete Schur control is a mode partition, not a
new geometric boundary-pairing proof.

Section12 and `kernel_functional.py` now reconstruct the complete non-Gaussian
coherent field kernel and finite link action from a supplied parity-even CAR
matrix, including all higher vertices. The identity/composition requirement
uniquely fixes its finite CAR gluing functional. This is standard coherent-state
reconstruction implemented on the existing source-search path, not a derivation
of the physical CAR dictionary from P1/P2. Matrix logarithm and Grassmann link
logarithm are kept distinct, as are the empty-vacuum entry and the reference
eigenvalue. The interacting two-mode control retains the actual quarter-time
equation, the prescribed sea scalar and the central loop phase.

The adapter checks the existing diagonal five-mode min/max completion matrices
and exports every kernel/action coefficient. Both actions have degree10. The
min-trace result simplifies to `T=(I+3|31><31|)/4`, with link action `-B-3M`;
this is a historical search-objective artifact, not a selected physical vertex.
The original internal action `W=i^N` is preserved separately from normalized
`R=i^(N-N_Omega)`. Neither proves a physical quarter-time equation for these
search outputs, which supply no physical collar length. Full details and
independently reviewed proof limits are in section12 of `PROOF.txt`.

`E8-BOUNDARY-DICTIONARY.txt` retains the origin paper's reverse-bootstrap result within its
discrete E8-closure framework. `e8_boundary_response.py` supplies an explicit
dictionary of all248 currents of the stipulated local E8 level-one lattice
vacuum representation. Root vertices, their actual cocycle, all Cartan Ward
insertions and internal quarter grades give every separated-point current
source coefficient. The nonzero three-point and connected four-point probes
prevent a Gaussian-current shortcut. The original charged quarter-time
controls remain distinct from this internal grading. These are standard
lattice-current identities made executable here, not a new physical origin
theorem. The physical insertion map into the raw P1 spectral kernel,
original field-dependent operators, absolute amplitude and contact data
are not derived by the numerical reverse-bootstrap checks.

**PARTIAL for physical TFPT.** Exact source-relative constructions and positive
controls do not identify the physical field algebra or select a source from P1/P2.
The original action permits interaction vertices; the general OS sewing rule is
retained alongside the conditional free determinant rule. Experiment-only work;
no ledger, verification, paper or physical-gate promotion.

Existing results reused: `CAR.QFREE.01`, `GATE.METRIC.17`, `GATE.METRIC.18`,
`QFT.NETRP.01`, and `seam-full-collar-reconstruction-20261015`.

Primary background: [Dereziński, representations of CCR/CAR, Theorem 18 and
Section 5](https://www.fuw.edu.pl/~derezins/derezinski.pdf), and
[Bravyi, Lagrangian representation for fermionic linear optics](https://arxiv.org/abs/quant-ph/0404180).
