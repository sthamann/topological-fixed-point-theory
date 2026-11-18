# Joint seam readout, flow and deck investigation

Research ID: `UR.SEAM.READOUT_FLOW_DECK.01`. Verdict: **PARTIAL**.

This continues `seam-two-normalizations-origin-20261005` and answers the three
requested follow-ups together. It tests the existing native objects rather
than introducing a replacement dynamics or fitting a new response.

## Decisive results

1. The canonical D5 spinor block has weight 5/8 and two-point exponent 5/4;
   the full glued E8 root current and the affine U(1) current have weight 1.
   In a fixed charged-pair Ward calculation its common fractional amplitude
   cancels from normalized JJ and from connected JJ. The original charge
   response also separates R from R_Y: their Gram determinant is 460,
   R_JY_JY=0 and RY_JY_JY=115/36. A common source is possible, but a common
   scalar response is not established. The actual mark port retains 4 P0;
   the zero-sum family projection annihilates P0. The contact strength
   epsilon=12 c3^3 exp(-2 alpha)/ln(2) is only an explicitly declared
   target-matching equality in the original implementation.
2. Native v486 has spectrum 1, 2/3, 1/3 and a one-dimensional unit
   eigenspace. Its pair transfer hops in sheet coordinates but commutes
   with the sheet swap, so it does not exchange parity eigenspaces.
   This excludes two independent stationary unit modes in this linear
   transfer, not every possible nonlinear or conditioned dynamics.
   The archival RG protocol explicitly describes crossings. With the
   native one-loop QCD coefficient -7 the beta function is nonzero at
   c3 and phi_base (and at the corrected positive phi0). No two-loop
   crossing scale is newly reproduced here.
3. On the declared bounding/NS Cl16 spin structure, the exact native
   half-shift implementer satisfies U^2=fermion parity and U^4=1.
   It commutes with parity and P_minus U P_plus=0. Consequently this
   symmetry cannot remove an otherwise allowed relative sector offset
   by exchanging sectors. This is a genuine projective spin lift, but
   geometric sheets, Fock parity and D5/A3 chirality still need distinct
   explicit dictionaries. No compiler/process separation follows.

The conditional 4/3 relation from the prior contract remains available:
phi_base=(4/3)c3. It is consistent with two readouts of a common object;
the source selection of the physical inverse-norm readout remains open.

## Run

From the repository: `python3 -B experiments/theory-contracts/seam-readout-flow-deck-20261005/checker.py`.
Dependencies are the existing Python environment's NumPy and SymPy. The
checker extracts unchanged native lattice, transfer and Clifford functions,
calls the actual joint charge and mark-port readouts, and performs exact
symbolic Ward calculations. It checks source consistency during evaluation.
Its checks remain active under `-OO`. Results are in `results.json`.

## First missing source choice

One physically selected charged seam functional must specify its state,
defect/observable measure and time action. The existing conditional
constructions must be preserved: a charged sector sum or another defect
response can retain fractional powers even where a fixed-sector connected
current response cannot. The present calculation does not exclude those
possibilities. It identifies the precise readout choice needing derivation.

The active neighbouring source-holonomy contract has substantial conditional
source constructions. Its evolving summaries are not promoted here to a
P1-derived physical source or coupling. The actual explorer functions used
here explicitly retain their physical state/action/matching limitations.

## Firewall and integration

These are internal structural tests, not empirical evidence. No verification
module, paper, ledger or scorecard is modified. No new load-bearing claim is
promoted. A local contract index is supplied; the shared theory graph is
already stale and actively edited elsewhere, so it is not regenerated here.
The user-facing research note is recorded in this contract and the output
report, leaving the concurrently edited `experiments/next.txt` intact.

For the standard current weights and Ward/OPE conventions, see J.-B. Zuber,
*An Introduction to Conformal Field Theory*, equations 93, 99–100 and the
two-point discussion: https://www.lpthe.jussieu.fr/~zuber/MesPapiers/zako.pdf.
