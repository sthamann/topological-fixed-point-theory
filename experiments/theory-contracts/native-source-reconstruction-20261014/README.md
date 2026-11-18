# Native source: a joint state, time and response reconstruction

Research ID: `UR.SOURCE.NATIVE_TIME_RESPONSE.01`. Overall verdict **PARTIAL**; the mathematical finite-volume reconstruction is proved conditionally on the existing source. No physical T1–T8 gate is closed.

Read `PROOF.md` for the complete argument, `TRANSFER-REVIEW.md` for independent domain/limit review and `OPERATOR-MAP-REVIEW.md` for independent native-operator review.

The original Round20–26 Hamiltonian and transfer give a single state/time/response family. The temporal regulator can be removed with Gibbs trace-norm convergence. Native adjacent transports also prove why a lossless four-dimensional identification fails and why onsite parity cannot predict all future responses.

Run `python3 neutral_transport_check.py` and `python3 reconstruction_check.py` with SymPy and NumPy. The first reads only the original Gram assignments and checks exact native integer word actions and energy differences. The second checks the normalized static-charge Gaussian block. Neither executes the entire moving signed-charge history sum. The analytical proofs cover its declared Hamiltonian.

The source, couplings, finite spatial lattice, beta and physical probes remain stated premises. This result does not select them from P1/P2.
