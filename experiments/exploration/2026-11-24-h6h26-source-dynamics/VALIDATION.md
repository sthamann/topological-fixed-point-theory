# Validation of the H6h27 continuation

## Observed results

153 targeted tests passed in 146.224 seconds. This includes the existing process and
genesis tests, the merged H6h27 metric tests, the analytic certificates, the weighted
metric parent, the physical source and amplitude audit, the boundary four-point
response, local depth controls and the full five-parameter HR classification.

The complete lane was then run in separate processes with PYTHONHASHSEED=1 and 77.
The two result files are byte-identical. The canonical result is lane_results.json.

SHA256:
2b9ebdf44bb1749799ca3d7ef413897d76af31f9aaf3e2d996c2ce2647e89d01

This validation does not claim that the entire repository suite or all legacy
exploration scripts were rerun. The original z48 through z51 sources were preserved
and their hashes are included in the result. The earlier H6h26 certificate checkpoint
has its separate 116-test log; the current result supersedes that count.

No ledger promotion, Atlas update or assertion of full interacting quantum
positivity follows from passing these checks.

## Reproduce the targeted tests

From research/originlab:

~~~bash
python3 -m unittest \
  test_process_rules test_genesis_rules \
  test_process_matrix_certificates test_process_spin2_certificates \
  test_process_transfer_certificates test_process_spectral_no_go \
  test_process_star_metric test_process_weighted_star \
  test_process_gaussian_parent test_process_source_bimetric \
  test_process_flavour_obstructions test_process_source_dynamics \
  test_process_boundary_fourpoint test_process_depth_control \
  test_process_metric_potential_certificates -v
~~~

## Reproduce the complete lane

From the repository root:

~~~bash
env PYTHONHASHSEED=1 python3 \
  experiments/exploration/2026-11-24-h6h26-source-dynamics/run_audit.py \
  --lane --output /tmp/tfpt-h6h27-seed1.json
~~~

Repeat with PYTHONHASHSEED=77 and a different output path. Compare the files byte
for byte. Omitting --lane runs only the new certificate collection.

The runner serializes with sorted keys and rejects NaN. It records all process
module hashes and the four historical input-script hashes. Numerical values are
explicit diagnostics where used; exact identities decide the new certificates.

## Environment

| Component | Version |
| --- | --- |
| Python | 3.12.14 |
| NumPy | 2.3.5 |
| SciPy | 1.17.0 |
| SymPy | 1.14.0 |
| mpmath | 1.3.0 |

The physical model is conditional. The JSON fields physical_solution and
global_process_unique are false by design. A finite enumeration with a unique
minimum is not a theorem about every possible process.

## Independent reviews

The full weighted star was independently checked in all four sectors. The
boundary four-point flux identity, physical source normalization, closed flat-star
profile, weighted Li3 infrared factor and determinant polarization were reviewed
independently of their implementation. The proofs are stored next to this file;
the tests include counterexamples and are not substitutes for the general proofs.
