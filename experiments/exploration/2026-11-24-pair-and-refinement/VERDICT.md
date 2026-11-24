# Pair dynamics and edge refinement (2026-11-24)

Exploratory results after H6h28, source baseline 5537c03fe. These are independently executed model checks, not a proof of the physical source identification or of T1–T8. The two scripts were rerun during integration: 43/43 and 87/87 grouped checks passed.

- `tfpt_h6h29_checks.py`: local Spin(10) pair Hamiltonian, conditional Spin(11) spectrum and moments, strong projection counterexample, native family-clock compression and area-variable potential. The local gap is not a four-dimensional mirror gap.
- `tfpt_refinement_check.py`: classical diagonal tetrad refinement tests, geodesic alternative, native Clifford pair contractions and chiral stress spectral moments. The new geodesic potential is outside the HR class; no nonlinear constraint or quantum positivity result is asserted.
- Matching pair tensor and CFT responses do not supply the physical four-dimensional source map.
- The continuum chiral regulator route remains an external research candidate, not a completed continuum proof.

Keep the ledger and frozen predictions unchanged. Results of the exploration are conditional; no bench promotion.
