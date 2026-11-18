# Boundary, seam and clocks — part 8

[Topic entry](../boundary-clocks.md) · [Registry entry](../../INDEX.md)

Exact ledger statuses remain at the linked original rows and in registry.json. Citations are source links, not proof certification.

## 2026-07-21 (II) (PG.06b FULL: the complete NICER Vela archive reduced — the first X-ray [formula in source] leg of the recovery-comb search, verdict data_limited; the old PG.06b ``pipeline proven'' fold corrected on record; experiments level only — no verification module, no gate closes, no marker moves)

[documents/changelog/sections/007-source.tex:310-342](../../../../documents/changelog/sections/007-source.tex#L310) · `section:e2d9ada460279647968e` · subsection

experiments/pulsar-glitch-recovery PG.06b FULL executed ( tfpt-pulsar vela-full). The full NICER Vela-pulsar (PSR B0833-45) L2 archive is reduced, with the manifest frozen before the bulk download ( manifest_full_reduction.json: all 665 ObsIDs, preregistered 2019/2021/2024 glitch windows, 20-GB abort gate): 665/665 obs downloaded (2.6\,GB, raw gitignored) [formula in source] 491 usable [formula in source] 240 per-obs

Explicit source references:

- [experiments/README.md](../../../../experiments/README.md#L1) · explicit_file_citation
- [experiments/next.txt](../../../../experiments/next.txt#L1) · explicit_file_citation
- [experiments/pulsar-glitch-recovery/data/nicer_vela/manifest_full_reduction.json](../../../../experiments/pulsar-glitch-recovery/data/nicer_vela/manifest_full_reduction.json#L1) · explicit_file_citation

## 2026-07-20 (VIII) (E8 mass-ladder blind-recovery bed EXECUTED: verdict validated at the detector level — the Cartan–PF ladder pipeline blindly recovers the published Zamolodchikov-E8 spectra of BaCo [formula in source] V [formula in source] O [formula in source] (8/8) and CoNb [formula in source] O [formula in source] (6/8 in-window) and rejects every control class; experiments level only — confirms Zamolodchikov E8, not TFPT; no verification module, no gate closes, no marker moves)

[documents/changelog/sections/008-source.tex:1-30](../../../../documents/changelog/sections/008-source.tex#L1) · `section:967792ae6a0e25788ac8` · subsection

experiments/e8-ladder-bed/ executed (preregistered hypotheses/e8_ladder_bed_v1.yaml). The blind ladder detector (assignment DP [formula in source] scale fit [formula in source] 5000-ladder Monte-Carlo null; the same detector class TFPT uses on search surfaces, e.g.\ hfqpo-ladder) on the only public known-answer mass-ladder dataset: BaCo [formula in source] V [formula in source] O [formula in source] (Zou [formula in

Explicit source references:

- [verification/v1_e8_glue.py](../../../../verification/v1_e8_glue.py#L1) · explicit_script_identifier
- [experiments/next.txt](../../../../experiments/next.txt#L1) · explicit_file_citation
- [experiments/README.md](../../../../experiments/README.md#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v1) E8 glue arithmetic ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L46)

## 2026-07-20 (VII) ( v493 — WP2 of CELEST.SEAM.01 executed: the clock-invariant deformation [formula in source] is a PURE SCALE with the [formula in source] pillowcase frozen, the [formula in source] clock IS the Picard–Lefschetz/Coxeter monodromy of the seam deformation, and the BHS deformed-algebra pattern transfers [formula in source] with no sector leak; verdict B, contract updated to WP2 done, no gate closes, no marker moves, SEAM.EQUIV.01 and the v216 residual stay )

[documents/changelog/sections/008-source.tex:31-72](../../../../documents/changelog/sections/008-source.tex#L31) · `section:7750a727b522e861fdbc` · subsection

v493_celestial_wp2_clock_deformation.py ( CELEST.WP2.01). WP2 of the celestial contract, sympy exact (47 checks), promoted from the experiments/tfpt-discovery/ probe celestial_seam_wp2_clock_deformation_probe.py. Geometry: [formula in source] forces [formula in source] sharply and two-sided ( [formula in source] ); [formula in source] is smooth iff [formula in source] ( [formula in source] ); the branch points [formu

Explicit source references:

- [verification/v491_p2_partition_corollary.py](../../../../verification/v491_p2_partition_corollary.py#L1) · explicit_script_identifier
- [GATE.WOLFRAM.02](../../../../verification/status_ledger.csv#L168) · explicit_claim_identifier · original status at linked ledger row
- [verification/wolfram/tfpt_readouts_extension.wl](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1) · explicit_file_citation
- [verification/v493_celestial_wp2_clock_deformation.py](../../../../verification/v493_celestial_wp2_clock_deformation.py#L1) · explicit_file_citation
- [CELEST.SEAM.01](../../../../verification/status_ledger.csv#L556) · explicit_claim_identifier · original status at linked ledger row
- [verification/v168_qgeo_rigidity.py](../../../../verification/v168_qgeo_rigidity.py#L1) · explicit_script_identifier
- [experiments/tfpt-discovery/celestial_seam_wp2_clock_deformation_probe.py](../../../../experiments/tfpt-discovery/celestial_seam_wp2_clock_deformation_probe.py#L1) · explicit_file_citation
- [verification/v492_celestial_z4_orbifold.py](../../../../verification/v492_celestial_z4_orbifold.py#L1) · explicit_script_identifier
- [verification/v214_seam_pillowcase.py](../../../../verification/v214_seam_pillowcase.py#L1) · explicit_script_identifier
- [CELEST.WP2.01](../../../../verification/status_ledger.csv#L557) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [(* ==== v491 round: P2.PARTITION.01 -- g_car = 5 as the partition corollary of the four marks (the sharpening of v53). Exact content: 4 into exactly 3 positive integer parts is uniquely {1,1,2} => e = (4,5,2) => g_car = e2 = 5, |Z2| = e3 = 2, N_fam = e2 - e3 = 3 as corollaries; every assumption load-bearing (positivity / sum / part count / integrality negative controls); v53 needs all three targets, the partition only e1 = 4; sum rule 4(N-1)/2 = 4 iff N = 3; mirrors v491. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3117)
- [experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean#L1)
- [(* ---- (v168) QGEO rigidity: mu4 square cross-ratio, mu4 characters = A3 exponents ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1453)
- [(* ---- (v214) QGEO.PILLOW.01: pillowcase reduction -- cross-ratio 2 => j=1728 => order-4 CM ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1727)
- [(* ==== v493 round: CELEST.WP2.01 -- the clock-invariant deformation XY = Z^4 + a0 of the A3 seam orbifold (WP2 of CELEST.SEAM.01). Exact content mirrored: clock-invariance selects the 1-parameter slice sharply (P(iZ) = P(Z) forces a3 = a2 = a1 = 0; prod(Z - i^k w) = Z^4 - w^4); disc = 256 a0^3; binary-quartic invariants I = 12 a0, J = 0 identically => j = 1728 frozen (tau = i pillowcase), cross-ratio 2; negative controls a1 (j = 0) and a2-instance (j = 1556068/81); clock on H2 = Coxeter element of W(A3) (char x^3+x^2+x+1, eigenvalues {i,-1,-i}, no invariant cycle, chi_1-Fourier eigenline, period point t(i-1)(1,i,i^2)); versal weights (1,2,3,0) = regular mu4 rep; S-algebra corrections linear in a0 at the Z4 wrap (6/16 pairs) with the EH k = 2 analogue; mirrors v493. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3148)

## 2026-07-20 (V) ( v492 [formula in source] the new research contract CELEST.SEAM.01 — the celestial-holographic route: the [formula in source] [formula in source] -glue is the flat [formula in source] monodromy of the equivariant celestial chiral algebra on the [formula in source] ALE space, with the exact correction clock [formula in source] deck (spin bridge [formula in source] , [formula in source] ); verdict B, no gate closes, no marker moves, SEAM.EQUIV.01 stays )

[documents/changelog/sections/008-source.tex:98-130](../../../../documents/changelog/sections/008-source.tex#L98) · `section:0bcac4df826bc0d57d23` · subsection

v492_celestial_z4_orbifold.py ( CELEST.WP1.01). WP1 of the new contract, sympy exact (36 checks): the v128 [formula in source] -glue grading of [formula in source] ( [formula in source] , additive on all [formula in source] root pairs) is inner — [formula in source] reads the glue class mod [formula in source] on all [formula in source] roots, so the glue is [formula in source] , an order- [formula in source] inner t

Explicit source references:

- [CONTRACT.F.01](../../../../verification/status_ledger.csv#L302) · explicit_claim_identifier · original status at linked ledger row
- [verification/v492_celestial_z4_orbifold.py](../../../../verification/v492_celestial_z4_orbifold.py#L1) · explicit_script_identifier
- [CELEST.WP1.01](../../../../verification/status_ledger.csv#L555) · explicit_claim_identifier · original status at linked ledger row
- [verification/v92_glue_uniqueness.py](../../../../verification/v92_glue_uniqueness.py#L1) · explicit_script_identifier
- [verification/v128_graded_hull.py](../../../../verification/v128_graded_hull.py#L1) · explicit_script_identifier
- [verification/v377_seam_s3_e8character.py](../../../../verification/v377_seam_s3_e8character.py#L1) · explicit_script_identifier
- [CELEST.SEAM.01](../../../../verification/status_ledger.csv#L556) · explicit_claim_identifier · original status at linked ledger row
- [verification/v168_qgeo_rigidity.py](../../../../verification/v168_qgeo_rigidity.py#L1) · explicit_script_identifier
- [verification/v125_glue_qsystem.py](../../../../verification/v125_glue_qsystem.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v214_seam_pillowcase.py](../../../../verification/v214_seam_pillowcase.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)
- [(* ---- (v92) glue uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L151)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ---- (v125) glue Q-system ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L872)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v168) QGEO rigidity: mu4 square cross-ratio, mu4 characters = A3 exponents ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1453)
- [experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean#L1)
- [(* ---- (v214) QGEO.PILLOW.01: pillowcase reduction -- cross-ratio 2 => j=1728 => order-4 CM ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1727)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [(* ---- (v128) graded hull ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L918)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)

## 2026-07-20 (IV) (evidence scorecard: five new typed rows — the [formula in source] running prediction of record, the [formula in source] -distortion [formula in source] branch discriminator, and three parked analog beds (bosonic-E8 thermal Hall, E8 mass-ladder recovery, Efimov DSI tuning); experiments level only — no verification module, no gate closes, no marker moves)

[documents/changelog/sections/008-source.tex:131-182](../../../../documents/changelog/sections/008-source.tex#L131) · `section:7b5da22330a89ad066dc` · subsection

[formula in source] running (prediction of record, consistent). The Starobinsky branch fixes the running of the spectral index: [formula in source] at [formula in source] (exact slow roll [formula in source] ; band [formula in source] : [formula in source] ) and [formula in source] (leading order; NOT [formula in source] ) — the identity [formula in source] holds only at leading order ( [formula in source] 7\% off ex

Explicit source references:

- [experiments/tfpt-discovery/alpha_s_running_probe.py](../../../../experiments/tfpt-discovery/alpha_s_running_probe.py#L1) · explicit_file_citation
- [experiments/build_evidence_scorecard.py](../../../../experiments/build_evidence_scorecard.py#L1) · explicit_file_citation
- [experiments/tfpt-discovery/efimov_s0_massratio_probe.py](../../../../experiments/tfpt-discovery/efimov_s0_massratio_probe.py#L1) · explicit_file_citation
- [verification/v451_seam_edge_cardy_tower.py](../../../../verification/v451_seam_edge_cardy_tower.py#L1) · explicit_script_identifier
- [verification/v450_seam_edge_entanglement.py](../../../../verification/v450_seam_edge_entanglement.py#L1) · explicit_script_identifier
- [experiments/tfpt-discovery/mu_distortion_probe.py](../../../../experiments/tfpt-discovery/mu_distortion_probe.py#L1) · explicit_file_citation
- [verification/v489_seam_modular_commutator.py](../../../../verification/v489_seam_modular_commutator.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [(* ==== v450/v451 round: SEAM.EQUIV.EDGE.* -- the edge central charge c_-=8 from entanglement (c=1/2 per Majorana) and the Ising minimal-model Kac weights (Cardy tower); mirrors v450_seam_edge_entanglement.py / v451_seam_edge_cardy_tower.py (the numerical fits are Python-only; the c_-=8 assembly and the Kac weights are exact and mirrored here). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2793)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)

## 2026-07-20 (II) ( v490 — the [formula in source] spin-structure parity census: [formula in source] witnessed on the lattice, the 16-fold-way bridge [formula in source] on measured numbers, and the honest TEE non-discriminator fence; no gate closes, no marker moves)

[documents/changelog/sections/008-source.tex:208-228](../../../../documents/changelog/sections/008-source.tex#L208) · `section:67f92d5b57079948ad55` · subsection

v490_seam_parity_census.py ( SEAM.MU.01). The [formula in source] ground-state fermion parity of the v367 [formula in source] collar, resolved per spin structure [formula in source] and determined TWICE — Read–Green TRIM product AND real-space Majorana Pfaffian (Schur), relative to the trivial [formula in source] phase, with the twisted real-space spectra matching the Bloch bands to [formula in source] . \ per layer

Explicit source references:

- [SEAM.MU.01](../../../../verification/status_ledger.csv#L552) · explicit_claim_identifier · original status at linked ledger row
- [verification/v490_seam_parity_census.py](../../../../verification/v490_seam_parity_census.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v378_seam_s3_modular.py](../../../../verification/v378_seam_s3_modular.py#L1) · explicit_script_identifier
- [verification/v489_seam_modular_commutator.py](../../../../verification/v489_seam_modular_commutator.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)

## 2026-07-20 ( v489 — the chiral central charge [formula in source] measured DIRECTLY from the bulk ground state via the Kim–Shehab–Kim modular commutator: the first ground-state witness for premise (i) of SEAM.EQUIV.01, independent of the Chern route ( v367) and the KLM count ( v378); no gate closes, no marker moves)

[documents/changelog/sections/008-source.tex:229-250](../../../../documents/changelog/sections/008-source.tex#L229) · `section:fa69e504a169a4dd0514` · subsection

v489_seam_modular_commutator.py ( SEAM.CCC.01). The Kim–Shehab–Kim modular commutator [formula in source] (PRL 128 (2022); arXiv:2110.06932/ 2110.10400) is evaluated on the EXACT BdG ground state of the v367 [formula in source] seam collar model, entirely from the Peschel Majorana covariance ( [formula in source] , [formula in source] ), on a disk about the torus centre cut into three [formula in source] wedges. \ gr

Explicit source references:

- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v489_seam_modular_commutator.py](../../../../verification/v489_seam_modular_commutator.py#L1) · explicit_file_citation
- [SEAM.CCC.01](../../../../verification/status_ledger.csv#L551) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v378_seam_s3_modular.py](../../../../verification/v378_seam_s3_modular.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)

## 2026-07-15 (II) ( v487 — the lazy split FORCED by the resummed clock: the walk's one-step spectrum is the complete v124 clock ladder below the wall, deck parity [formula in source] rung index, and the two comb tones are one bend apart — the v486 \ ``split selection'' closes onto established structure; no gate closes, no marker moves)

[documents/changelog/sections/008-source.tex:291-317](../../../../documents/changelog/sections/008-source.tex#L291) · `section:3d01f167d55724948536` · subsection

v487_transfer_clock_rungs.py ( HYP.REWRITE.03). v486 left the split selection open (`` [formula in source] vs uniform [formula in source] is selected by the physical spectrum, not yet derived''; named candidate: the vague ``sheet-symmetric choice''). Closed: \ RUNG READING — the lazy rule's one-step spectrum is EXACTLY the complete resummed-clock ladder below the wall ( v124: [formula in source] , pole at [formula in

Explicit source references:

- [verification/v487_transfer_clock_rungs.py](../../../../verification/v487_transfer_clock_rungs.py#L1) · explicit_script_identifier
- [HYP.REWRITE.03](../../../../verification/status_ledger.csv#L548) · explicit_claim_identifier · original status at linked ledger row
- [verification/v107_quantum_clock_target.py](../../../../verification/v107_quantum_clock_target.py#L1) · explicit_script_identifier
- [experiments/theory-contracts/lazy01_clock_rung_forcing.py](../../../../experiments/theory-contracts/lazy01_clock_rung_forcing.py#L1) · explicit_file_citation
- [verification/v124_resummed_clock.py](../../../../verification/v124_resummed_clock.py#L1) · explicit_script_identifier
- [verification/v486_transfer_full_rule.py](../../../../verification/v486_transfer_full_rule.py#L1) · explicit_script_identifier
- [verification/v327_hypergraph_rewrite.py](../../../../verification/v327_hypergraph_rewrite.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v327 NON-GRAPH-SPECTRAL: 2/3 is NOT a root of the affine-E8 charpoly (sharpens v312) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2232)
- [(* v327 HYP.REWRITE.01: the minimal 3-channel/1-absorbing rewrite M has spectrum {0,2/3,1} *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2227)
- [(* ---- (v107) quantum-clock target ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L472)
- [(* ---- (v124) resummed clock ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L857)
- [(* ==== v325/v327: the keystone pillowcase orbifold + the cusp-weight rewrite (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2221)

## 2026-07-15 ( v486 — generator hunt: the FULL transfer spectrum [formula in source] from ONE forced local rule — the unique lazy [formula in source] -pair walk [formula in source] ; the previously generator-less [formula in source] mode gets its origin; v327's uniform rule superseded as full-spectrum generator — no gate closes, no marker moves)

[documents/changelog/sections/008-source.tex:318-342](../../../../documents/changelog/sections/008-source.tex#L318) · `section:4d7b24e6cd695cd3409f` · subsection

v486_transfer_full_rule.py ( HYP.REWRITE.02). The verified seam transfer spectrum is [formula in source] ( v54/ v221), but the only local generator in the suite ( v327) produces [formula in source] — its zero mode persists under every power, so no iterate reaches the third mode. \ UNIQUENESS: for the symmetric [formula in source] -channel rule [formula in source] (one absorbing family channel [formula in source] [for

Explicit source references:

- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [verification/v355_e8_unmapped_bandwidth.py](../../../../verification/v355_e8_unmapped_bandwidth.py#L1) · explicit_script_identifier
- [experiments/tfpt-discovery/generator_hunt_full_transfer_rule.py](../../../../experiments/tfpt-discovery/generator_hunt_full_transfer_rule.py#L1) · explicit_file_citation
- [verification/v486_transfer_full_rule.py](../../../../verification/v486_transfer_full_rule.py#L1) · explicit_file_citation
- [verification/v327_hypergraph_rewrite.py](../../../../verification/v327_hypergraph_rewrite.py#L1) · explicit_script_identifier
- [verification/v54_seam_horizon_keystones.py](../../../../verification/v54_seam_horizon_keystones.py#L1) · explicit_script_identifier
- [verification/v354_e8_reverse_audit.py](../../../../verification/v354_e8_reverse_audit.py#L1) · explicit_script_identifier
- [HYP.REWRITE.02](../../../../verification/status_ledger.csv#L547) · explicit_claim_identifier · original status at linked ledger row
- [verification/v221_seam_qecc.py](../../../../verification/v221_seam_qecc.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v355 E8.UNMAPPED.BANDWIDTH.01: the unmapped region's forced content is COLLECTIVE ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2363)
- [(* ---- (v54) seam=horizon keystones: triply-forced 8 + shared transfer operator ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L235)
- [(* ==== v354 E8.REVERSE.AUDIT.01: 3 of 8 E8 Casimir degrees are load-bearing, 5 unmapped ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2353)
- [(* v327 NON-GRAPH-SPECTRAL: 2/3 is NOT a root of the affine-E8 charpoly (sharpens v312) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2232)
- [(* ==== v325/v327: the keystone pillowcase orbifold + the cusp-weight rewrite (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2221)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [(* v327 HYP.REWRITE.01: the minimal 3-channel/1-absorbing rewrite M has spectrum {0,2/3,1} *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2227)

## 2026-07-14 (IV) ( v485 — the merged analytic target's computable side COMPLETED: the diagonal self-energy is renormalised to ZERO exactly at the KMS seam scale, the mark contact calculus resums in closed form (BFK route), and the [formula in source] / [formula in source] multiplicities are ONE state set under two response weights — the residual collapses onto the SEAM.EQUIV.01 face alone; no gate closes, no marker moves)

[documents/changelog/sections/008-source.tex:343-369](../../../../documents/changelog/sections/008-source.tex#L343) · `section:d1b557c6c16fcdfd5e58` · subsection

v485_contact_diagonal_closed.py ( SEAM.CONTACT.UNIT.02). The v484 merge left one analytic step: the diagonal (self-energy) zeta-renormalisation at the marks [formula in source] the multiplicity matching. Settled at the computable level: \ DIAGONAL ZERO AT THE KMS SCALE — the renormalised coincident-point Green function of [formula in source] on circumference [formula in source] is [formula in source] EXACTLY (symboli

Explicit source references:

- [SEAM.CONTACT.UNIT.02](../../../../verification/status_ledger.csv#L546) · explicit_claim_identifier · original status at linked ledger row
- [verification/v239_kms_thermal_time.py](../../../../verification/v239_kms_thermal_time.py#L1) · explicit_script_identifier
- [verification/v485_contact_diagonal_closed.py](../../../../verification/v485_contact_diagonal_closed.py#L1) · explicit_file_citation
- [verification/v151_bfk_split.py](../../../../verification/v151_bfk_split.py#L1) · explicit_script_identifier
- [verification/v484_seam_contact_unit.py](../../../../verification/v484_seam_contact_unit.py#L1) · explicit_script_identifier
- [ALPHA.QUILLEN.EXACT.01](../../../../verification/status_ledger.csv#L445) · explicit_claim_identifier · original status at linked ledger row
- [HYP.PHI0.PUNCTURE.01](../../../../verification/status_ledger.csv#L464) · explicit_claim_identifier · original status at linked ledger row
- [verification/v470_alpha_inflow_level.py](../../../../verification/v470_alpha_inflow_level.py#L1) · explicit_script_identifier
- [verification/v382_alpha_quillen_exact.py](../../../../verification/v382_alpha_quillen_exact.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v151) BFK split ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1342)
- [(* ==== v470 round: the alpha^3 level / embedding-index rigidity -- the exact k_Y = 5/3 arithmetic (tr Y^2 / tr T3^2 over the 5bar), the b1 conversion (3/5)(41/6) = 41/10, and the carrier decomposition 41/6 = 20/3 + 1/6; mirrors v470 (the FHS Chern and the alpha^-1 root re-verification are numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2943)

## 2026-07-14 (III) ( v484 — the two open analytic targets MERGE: the shared `` [formula in source] per boundary insertion'' rule of ALPHA.QUILLEN.EXACT.01 and HYP.PHI0.PUNCTURE.01 is the KMS seam unit [formula in source] with [formula in source] , derived on the seam circle for the finite cycle sector — no gate closes, no marker moves)

[documents/changelog/sections/009-source.tex:1-26](../../../../documents/changelog/sections/009-source.tex#L1) · `section:25a3a7ca87e0691a61b8` · subsection

v484_seam_contact_unit.py ( SEAM.CONTACT.UNIT.01). The v483 round left the two analytic targets sharing one rule (``weight [formula in source] per boundary contact insertion'': the [formula in source] ladder [formula in source] , the [formula in source] puncture [formula in source] ). This module shows the shared rule is NOT a new unknown: \ the unit chain [formula in source] — one bare boundary insertion divided by

Explicit source references:

- [verification/v483_pillowcase_exact_traces.py](../../../../verification/v483_pillowcase_exact_traces.py#L1) · explicit_script_identifier
- [verification/v484_seam_contact_unit.py](../../../../verification/v484_seam_contact_unit.py#L1) · explicit_script_identifier
- [SEAM.CONTACT.UNIT.01](../../../../verification/status_ledger.csv#L545) · explicit_claim_identifier · original status at linked ledger row
- [verification/v354_e8_reverse_audit.py](../../../../verification/v354_e8_reverse_audit.py#L1) · explicit_script_identifier
- [verification/v355_e8_unmapped_bandwidth.py](../../../../verification/v355_e8_unmapped_bandwidth.py#L1) · explicit_script_identifier
- [verification/v58_seam_horizon_chain.py](../../../../verification/v58_seam_horizon_chain.py#L1) · explicit_script_identifier
- [ALPHA.QUILLEN.EXACT.01](../../../../verification/status_ledger.csv#L445) · explicit_claim_identifier · original status at linked ledger row
- [verification/v239_kms_thermal_time.py](../../../../verification/v239_kms_thermal_time.py#L1) · explicit_script_identifier
- [verification/v73_k_c3_half.py](../../../../verification/v73_k_c3_half.py#L1) · explicit_script_identifier
- [verification/v60_lambda_metrology_branch.py](../../../../verification/v60_lambda_metrology_branch.py#L1) · explicit_script_identifier
- [HYP.PHI0.PUNCTURE.01](../../../../verification/status_ledger.csv#L464) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v60) Lambda branch selection: prefactor 3/(4pi^2), branch mis-scale 64pi^3/3 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L283)
- [(* ---- (v73) k=c3/2 forced: variational factor x Gauss-Bonnet topology; Fursaev-Solodukhin S=A/4 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L353)
- [(* ==== v354 E8.REVERSE.AUDIT.01: 3 of 8 E8 Casimir degrees are load-bearing, 5 unmapped ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2353)
- [(* ==== v355 E8.UNMAPPED.BANDWIDTH.01: the unmapped region's forced content is COLLECTIVE ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2363)
- [(* ---- (v58) seam-horizon chain: one-sided Gauss-Bonnet, KMS, S_BH, RT in seam units ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L270)

## 2026-07-14 ( v479– v481 — the bird's-eye ``overlooked tools'' round: R3's graph [formula in source] geometry bridge discharged to Kronheimer 1989; the four [formula in source] marks realised as the exactly solvable four-interval multilocal modular geometry with [formula in source] manifest; and the absolute [formula in source] -scale collapsed to a one-parameter, PS-window-bracketed candidate under the carrier normalisation [formula in source] — no gate closes, no marker moves)

[documents/changelog/sections/009-source.tex:131-192](../../../../documents/changelog/sections/009-source.tex#L131) · `section:e99d6098550f62e277de` · subsection

v479_kronheimer_quiver_bridge.py ( SEAM.KRONHEIMER.01). The R3 residual of the [formula in source] keystone ( v344: ``the rewrite attractor graph IS the du Val singularity'') is an ESTABLISHED THEOREM — Kronheimer's ALE hyper-K\"ahler quotient construction [formula in source] Torelli classification (J. Diff. Geom. 29 (1989) 665 and 685) — with all finite input data already compiler outputs, machine-checked : the Kac

Explicit source references:

- [verification/v467_thirdgen_mixing_pattern.py](../../../../verification/v467_thirdgen_mixing_pattern.py#L1) · explicit_script_identifier
- [verification/v468_dm2_ratio_jarlskog.py](../../../../verification/v468_dm2_ratio_jarlskog.py#L1) · explicit_script_identifier
- [verification/v478_entropic_continuum_legs.py](../../../../verification/v478_entropic_continuum_legs.py#L1) · explicit_script_identifier
- [verification/v354_e8_reverse_audit.py](../../../../verification/v354_e8_reverse_audit.py#L1) · explicit_script_identifier
- [SEAM.DETK.01](../../../../verification/status_ledger.csv#L110) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v479_kronheimer_quiver_bridge.py](../../../../verification/v479_kronheimer_quiver_bridge.py#L1) · explicit_file_citation
- [FLAV.NUSCALE.02](../../../../verification/status_ledger.csv#L542) · explicit_claim_identifier · original status at linked ledger row
- [verification/v272_nu_mass_scale.py](../../../../verification/v272_nu_mass_scale.py#L1) · explicit_script_identifier
- [verification/v346_seam_geometric_bridge.py](../../../../verification/v346_seam_geometric_bridge.py#L1) · explicit_script_identifier
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.MULTILOCAL.01](../../../../verification/status_ledger.csv#L541) · explicit_claim_identifier · original status at linked ledger row
- [FLAV.NUSCALE.01](../../../../verification/status_ledger.csv#L365) · explicit_claim_identifier · original status at linked ledger row
- [verification/v249_ps_unification.py](../../../../verification/v249_ps_unification.py#L1) · explicit_script_identifier
- [experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean#L1) · explicit_file_citation
- [verification/v312_hypergraph_substrate.py](../../../../verification/v312_hypergraph_substrate.py#L1) · explicit_script_identifier
- [verification/v344_detk_synthesis.py](../../../../verification/v344_detk_synthesis.py#L1) · explicit_script_identifier
- [verification/v481_seesaw_carrier_ladder.py](../../../../verification/v481_seesaw_carrier_ladder.py#L1) · explicit_script_identifier
- [SEAM.KRONHEIMER.01](../../../../verification/status_ledger.csv#L540) · explicit_claim_identifier · original status at linked ledger row
- [experiments/theory-contracts/phi0_donnelly_equivariant_probe.py](../../../../experiments/theory-contracts/phi0_donnelly_equivariant_probe.py#L1) · explicit_file_citation
- [HYP.PHI0.PUNCTURE.01](../../../../verification/status_ledger.csv#L464) · explicit_claim_identifier · original status at linked ledger row
- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v355_e8_unmapped_bandwidth.py](../../../../verification/v355_e8_unmapped_bandwidth.py#L1) · explicit_script_identifier
- [verification/v480_multilocal_four_interval.py](../../../../verification/v480_multilocal_four_interval.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [(* ==== v344 SEAM.DETK.01: the ADE |det Cartan| = |H_1(link)| sequence, only E8 -> 1 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2272)
- [(* ==== v479 round: the Kronheimer quiver bridge -- the exact content: the Kac marks are the null vector of the affine-E8 Cartan matrix AND the Perron vector (A delta = 2 delta) with multiset {1,2,2,3,3,4,4,5,6} = 2I irrep dims; the hyper-Kaehler quotient dimension count dim_R M - 4 dim_R G = 4 with sum_edges d_i d_j = sum d_i^2 = 120 = |mu4| h(E8), h(E8) = 30 = 2*3*5; and the finite-E8 Cartan is unimodular with Smith normal form = identity (the Poincare-sphere link, the det K = 1 face); mirrors v479. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3089)
- [experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean#L1)
- [(* ==== v355 E8.UNMAPPED.BANDWIDTH.01: the unmapped region's forced content is COLLECTIVE ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2363)
- [(* ==== v354 E8.REVERSE.AUDIT.01: 3 of 8 E8 Casimir degrees are load-bearing, 5 unmapped ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2353)

## 2026-07-13 ( v478 — first computable steps on the two remaining legs of the entropic-action bridge: the state-side modular data flows to the CHM/Bisognano–Wichmann geometric form (Calabrese–Cardy [formula in source] , CHM parabola, even bands exactly zero) — the bridge's modular side meets TFPT's Einstein-derivation input; and the measure condition reduced to one exact KMS time [formula in source] with the [formula in source] near-miss explicitly declined; both legs stay , no status moves)

[documents/changelog/sections/009-source.tex:193-215](../../../../documents/changelog/sections/009-source.tex#L193) · `section:671808877239add5f598` · subsection

v478_entropic_continuum_legs.py ( GRAV.ENTROPIC.LEGS.01). Leg A (AP2 continuum, exhibited): on the exactly solvable critical chain (infinite-volume sinc kernel, [formula in source] -digit precision) the compressed interval state has [formula in source] with [formula in source] (monotone [formula in source] the Dirac [formula in source] , the v450 reading; gapped control [formula in source] , area law); the Peschel mo

Explicit source references:

- [verification/v478_entropic_continuum_legs.py](../../../../verification/v478_entropic_continuum_legs.py#L1) · explicit_file_citation
- [verification/v450_seam_edge_entanglement.py](../../../../verification/v450_seam_edge_entanglement.py#L1) · explicit_script_identifier
- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [GRAV.ENTROPIC.LEGS.01](../../../../verification/status_ledger.csv#L539) · explicit_claim_identifier · original status at linked ledger row
- [verification/v355_e8_unmapped_bandwidth.py](../../../../verification/v355_e8_unmapped_bandwidth.py#L1) · explicit_script_identifier
- [verification/v354_e8_reverse_audit.py](../../../../verification/v354_e8_reverse_audit.py#L1) · explicit_script_identifier
- [verification/v476_entropic_compression_toy.py](../../../../verification/v476_entropic_compression_toy.py#L1) · explicit_script_identifier
- [verification/v477_entropic_scaleflow.py](../../../../verification/v477_entropic_scaleflow.py#L1) · explicit_script_identifier
- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v477 (i): Frullani + the moment dictionary at flat weight *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3071)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [(* ==== v355 E8.UNMAPPED.BANDWIDTH.01: the unmapped region's forced content is COLLECTIVE ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2363)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [(* ==== v450/v451 round: SEAM.EQUIV.EDGE.* -- the edge central charge c_-=8 from entanglement (c=1/2 per Majorana) and the Ising minimal-model Kac weights (Cardy tower); mirrors v450_seam_edge_entanglement.py / v451_seam_edge_cardy_tower.py (the numerical fits are Python-only; the c_-=8 assembly and the Kac weights are exact and mirrored here). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2793)
- [(* ==== v354 E8.REVERSE.AUDIT.01: 3 of 8 E8 Casimir degrees are load-bearing, 5 unmapped ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2353)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [(* ==== v477 round: the scale-flow representation of the entropic action -- the exact content: the Frullani identity Integrate[(e^{-a t}-e^{-t})/t] = -Log[a], the moment dictionary (flat weight = unit moments = the v475 raw coefficients 3 b R + (17/24) b^2 R^2), the ONE moment condition mu2/mu1^2 = (72/17)(8 pi)^9 with the exact closure (4608 pi^2/17)/((72/17)(8 pi)^9) = c3^7, and the v36 consistency f0 = 6 (4 pi)^2/c3^7 <=> M^2 = c3^7 Mbar^2; mirrors v477 (the chi-form quadrature is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3062)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [(* v477 (ii): the one moment condition + closure + v36 consistency *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3081)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)

## 2026-07-13 ( v476 — the AP2 compression conjecture of the entropic-action bridge made well-posed on an exactly solvable model: the literal operator-side reading is ill-posed on a pure bulk, the state-side reading is forced, the mismatch is exactly quadratic in the cross-cut correlations and gap-suppressed; AP2 stays , no status moves)

[documents/changelog/sections/009-source.tex:239-262](../../../../documents/changelog/sections/009-source.tex#L239) · `section:67117a94f5952e2ebab9` · subsection

v476_entropic_compression_toy.py ( GRAV.ENTROPIC.COMPRESS.01). The last open work package of the v473 round, exercised on a gapped staggered-mass chain ( [formula in source] ). The seam situation is reproduced (pure bulk: the ground-state covariance is a projection to [formula in source] ; thermal boundary: [formula in source] , Peschel [formula in source] finite, [formula in source] — v155 in the bridge context). Th

Explicit source references:

- [GRAV.ENTROPIC.COMPRESS.01](../../../../verification/status_ledger.csv#L537) · explicit_claim_identifier · original status at linked ledger row
- [verification/v476_entropic_compression_toy.py](../../../../verification/v476_entropic_compression_toy.py#L1) · explicit_file_citation
- [verification/v337_decoupling_theorem.py](../../../../verification/v337_decoupling_theorem.py#L1) · explicit_script_identifier
- [verification/v155_quasifree_boundary.py](../../../../verification/v155_quasifree_boundary.py#L1) · explicit_script_identifier
- [verification/v473_entropic_action_bridge.py](../../../../verification/v473_entropic_action_bridge.py#L1) · explicit_script_identifier
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v473 (iii): the R^2 kill-test gap *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2985)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ==== v337 DECOUPLING.THEOREM.01: the ambient back-reaction scale identity ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2241)
- [(* ==== v473 round: the entropic-action bridge (Bianconi, PRD 111, 066001 (2025); arXiv:2408.14391) -- the exact algebraic content: the carrier Hodge count 1+5+10 = 16 = 2^(g_car-1) (vs the 4d fiber 11), the coefficient pin beta'_B = c3/6 = 1/(48 pi) from 3 beta' = 1/(16 pi), the entropy potential s(x) = x-1-ln x (nonnegative, quadratic minimum), the Lambda-channel target Tr Q^2 = 32 c3^4 = 1/(128 pi^4) = (2/3) delta_top, and the R^2 kill-test gap 3 c3^-9 = 3 (8 pi)^9; mirrors v473 (the Frullani quadrature is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2958)
- [(* v473 (ii): entropy potential + Lambda-channel target *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2975)
- [(* v473 (i): carrier Hodge count + coefficient pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2968)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)

## 2026-07-03 (surfacing round — three already-computed unique readouts promoted to the public prediction surfaces: the axion–photon coefficient [formula in source] , the free-seam Higgs band [formula in source] – [formula in source] \,GeV, and the never-run BH HFQPO [formula in source] ladder-tooth search; no new vN, no status moves)

[documents/changelog/sections/009-source.tex:346-371](../../../../documents/changelog/sections/009-source.tex#L346) · `section:1d709d6b4fbbfe522d4a` · subsection

Uniqueness scan follow-up (no new results — presentation only). A structure scan for ``what does TFPT predict that nobody else pins?'' found three concrete, falsifiable readouts that were fully computed and ledgered but missing from the public prediction surfaces; they are now surfaced, with their existing typing unchanged. (i) Axion–photon coefficient: the determinant-line anomaly coefficient [formula in source] , [

Explicit source references:

- [GRAV.ASYMP.01](../../../../verification/status_ledger.csv#L292) · explicit_claim_identifier · original status at linked ledger row
- [HIGGS.FREESEAM.01](../../../../verification/status_ledger.csv#L251) · explicit_claim_identifier · original status at linked ledger row
- [verification/v207_asymptotic_safety.py](../../../../verification/v207_asymptotic_safety.py#L1) · explicit_script_identifier
- [verification/v166_higgs_free_seam.py](../../../../verification/v166_higgs_free_seam.py#L1) · explicit_script_identifier

## 2026-07-03 ( v472 — the determinant line over the [formula in source] -twist moduli of the collar model carries curvature [formula in source] the inflow level: the v470 bridge lemma exhibited at the finite/model level; ALPHA.QUILLEN.EXACT.01 stays . Plus: bird's-eye rounds — electron-sector/Cabibbo whitespot probes, Calder\'on UV-universality, parameter-free flat-budget closure (sandbox); the flat-budget closure enters the evidence scorecard)

[documents/changelog/sections/010-source.tex:1-55](../../../../documents/changelog/sections/010-source.tex#L1) · `section:2327b1b5b72d8d0ad7fe` · subsection

v472_quillen_detline_moduli.py ( ALPHA.QUILLEN.DETLINE.01). ALPHA.QUILLEN.EXACT.01 names the [formula in source] level as the curvature of the determinant line over the [formula in source] seam moduli (Quillen 1985 / Bismut–Freed CMP 106 (1986) curvature; Dai–Freed JMP 35 (1994) section); v470 computed the bulk Chern invariant only over the Bloch Brillouin zone (the translation-invariant shortcut). v472 computes the

Explicit source references:

- [verification/v472_quillen_detline_moduli.py](../../../../verification/v472_quillen_detline_moduli.py#L1) · explicit_file_citation
- [verification/v470_alpha_inflow_level.py](../../../../verification/v470_alpha_inflow_level.py#L1) · explicit_script_identifier
- [verification/v471_seam_horizon_replica.py](../../../../verification/v471_seam_horizon_replica.py#L1) · explicit_script_identifier
- [experiments/tfpt-discovery/quillen_detline_twist_moduli.py](../../../../experiments/tfpt-discovery/quillen_detline_twist_moduli.py#L1) · explicit_file_citation
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [experiments/tfpt-discovery/electron_sector_cabibbo_probe.py](../../../../experiments/tfpt-discovery/electron_sector_cabibbo_probe.py#L1) · explicit_file_citation
- [experiments/tfpt-discovery/kernel_calderon_uv_universality.py](../../../../experiments/tfpt-discovery/kernel_calderon_uv_universality.py#L1) · explicit_file_citation
- [verification/v468_dm2_ratio_jarlskog.py](../../../../verification/v468_dm2_ratio_jarlskog.py#L1) · explicit_script_identifier
- [experiments/tfpt-discovery/flat_budget_closure.py](../../../../experiments/tfpt-discovery/flat_budget_closure.py#L1) · explicit_file_citation
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [QGEO.KERNEL.01](../../../../verification/status_ledger.csv#L257) · explicit_claim_identifier · original status at linked ledger row
- [ALPHA.QUILLEN.DETLINE.01](../../../../verification/status_ledger.csv#L533) · explicit_claim_identifier · original status at linked ledger row
- [verification/v204_muon_seam_g2.py](../../../../verification/v204_muon_seam_g2.py#L1) · explicit_script_identifier
- [ALPHA.QUILLEN.EXACT.01](../../../../verification/status_ledger.csv#L445) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.THEOREM.01](../../../../verification/status_ledger.csv#L131) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ---- (v204) FR.MUONG2.01: the muon seam-vertex value a_mu = 45/(524288 pi^9) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1696)
- [(* ==== v470 round: the alpha^3 level / embedding-index rigidity -- the exact k_Y = 5/3 arithmetic (tr Y^2 / tr T3^2 over the 5bar), the b1 conversion (3/5)(41/6) = 41/10, and the carrier decomposition 41/6 = 20/3 + 1/6; mirrors v470 (the FHS Chern and the alpha^-1 root re-verification are numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2943)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)

## 2026-07-02 ( v471 — the Seam–Horizon replica chain exercised numerically on the discretized collar: the kernel-identification premise of SEAM.THEOREM.01 discharged at the finite/model level; residual retypes to the continuum MMST leg [formula in source] the anchor; gate stays . Plus: the GW ringdown-echo channel hardened — 10 loudest events incl.\ GW250114, redshift-corrected templates, stacked search and a 12-variant signature battery, all null)

[documents/changelog/sections/010-source.tex:56-94](../../../../documents/changelog/sections/010-source.tex#L56) · `section:02788e96f36ad336a0ac` · subsection

v471_seam_horizon_replica.py ( SEAM.EHMODEL.04). The residual isolated by v150/ v151/ v152 — ``the RP seam kernel is this BFK Calder\'on datum'' — is exercised numerically on a discretized seam collar with the seam's own kernel and real replica sheets: the matched conical deficit has slope [formula in source] on deficit cones ( [formula in source] – [formula in source] ) and replica sheets [formula in source] ( [form

Explicit source references:

- [SEAM.EHMODEL.04](../../../../verification/status_ledger.csv#L532) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [experiments/gw-ringdown-echo/scripts/fetch_strain_4096.py](../../../../experiments/gw-ringdown-echo/scripts/fetch_strain_4096.py#L1) · explicit_file_citation
- [verification/v152_norm_is_anchor.py](../../../../verification/v152_norm_is_anchor.py#L1) · explicit_script_identifier
- [verification/v150_replica_eh_model.py](../../../../verification/v150_replica_eh_model.py#L1) · explicit_script_identifier
- [experiments/theory-contracts/seam_horizon_replica.py](../../../../experiments/theory-contracts/seam_horizon_replica.py#L1) · explicit_file_citation
- [verification/v151_bfk_split.py](../../../../verification/v151_bfk_split.py#L1) · explicit_script_identifier
- [SEAM.THEOREM.01](../../../../verification/status_ledger.csv#L131) · explicit_claim_identifier · original status at linked ledger row
- [verification/v471_seam_horizon_replica.py](../../../../verification/v471_seam_horizon_replica.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v302_seam_gap.py](../../../../verification/v302_seam_gap.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ---- (v150) replica EH model ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1322)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ---- (v151) BFK split ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1342)
- [(* ---- (v152) R3 normalisation = the anchor ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1356)

## 2026-07-02 ( v469/ v470 — the closure-route round on the two named targets: the SEAM.EQUIV.01 extension leg re-founded on 1995–2001 peer-reviewed crossed-product theorems [formula in source] the realisation axiom reduced to invariant level; the [formula in source] level [formula in source] the computed bulk Chern invariant [formula in source] the seam [formula in source] -normalisation [formula in source] the embedding index [formula in source] ; both targets stay )

[documents/changelog/sections/010-source.tex:108-175](../../../../documents/changelog/sections/010-source.tex#L108) · `section:ce5afee193d0ce1ea07e` · subsection

v469_seam_crossedproduct_route.py ( SEAM.EQUIV.CROSSEDPRODUCT.01). The [formula in source] -spinor extension leg of SEAM.EQUIV.01 no longer rests on 2024/2025 preprints: \ the [formula in source] spinor has [formula in source] , so its statistics phase is [formula in source] — exactly the Longo–Rehren locality criterion — and the [formula in source] simple-current crossed product [formula in source] is a local index-

Explicit source references:

- [verification/v463_seam_c8_holomorphic_uniqueness.py](../../../../verification/v463_seam_c8_holomorphic_uniqueness.py#L1) · explicit_script_identifier
- [ALPHA.QUILLEN.EXACT.01](../../../../verification/status_ledger.csv#L445) · explicit_claim_identifier · original status at linked ledger row
- [FORM.SEAM.RESIDUAL.01](../../../../verification/status_ledger.csv#L523) · explicit_claim_identifier · original status at linked ledger row
- [verification/v435_alpha_quillen_chernlevel.py](../../../../verification/v435_alpha_quillen_chernlevel.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamResidualAxiom.klm_holomorphic](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L198) · explicit_lean_symbol
- [verification/v433_alpha_quillen_heatkernel.py](../../../../verification/v433_alpha_quillen_heatkernel.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamResidualAxiom.crossedproduct_route_theorem](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L260) · explicit_lean_symbol
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [TfptCarrier.SeamResidualAxiom.crossedproduct_arithmetic](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L208) · explicit_lean_symbol
- [verification/v160_seam_gaussianity_from_pf.py](../../../../verification/v160_seam_gaussianity_from_pf.py#L1) · explicit_script_identifier
- [verification/v125_glue_qsystem.py](../../../../verification/v125_glue_qsystem.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamResidualAxiom.collar_invariants](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L258) · explicit_lean_symbol
- [verification/v469_seam_crossedproduct_route.py](../../../../verification/v469_seam_crossedproduct_route.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamResidualAxiom.lr_locality_integer](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L188) · explicit_lean_symbol
- [verification/v382_alpha_quillen_exact.py](../../../../verification/v382_alpha_quillen_exact.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamResidualAxiom.sixteenfold_pin](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L205) · explicit_lean_symbol
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [verification/v459_seam_lattice_voa_route.py](../../../../verification/v459_seam_lattice_voa_route.py#L1) · explicit_script_identifier
- [verification/v434_alpha_quillen_betafunction.py](../../../../verification/v434_alpha_quillen_betafunction.py#L1) · explicit_script_identifier
- [verification/v456_seam_chirality_from_c3.py](../../../../verification/v456_seam_chirality_from_c3.py#L1) · explicit_script_identifier
- [experiments/next.txt](../../../../experiments/next.txt#L1) · explicit_file_citation
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1) · explicit_file_citation
- [verification/v461_seam_strict_locality.py](../../../../verification/v461_seam_strict_locality.py#L1) · explicit_script_identifier
- [verification/v155_quasifree_boundary.py](../../../../verification/v155_quasifree_boundary.py#L1) · explicit_script_identifier
- [verification/v470_alpha_inflow_level.py](../../../../verification/v470_alpha_inflow_level.py#L1) · explicit_script_identifier
- [verification/v302_seam_gap.py](../../../../verification/v302_seam_gap.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamResidualAxiom.glue_locality_integers](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L195) · explicit_lean_symbol
- [ALPHA.QUILLEN.INFLOW.01](../../../../verification/status_ledger.csv#L531) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.CROSSEDPRODUCT.01](../../../../verification/status_ledger.csv#L530) · explicit_claim_identifier · original status at linked ledger row
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v454/v456/v457/v459 round (G-block): the edge current-algebra/Sugawara c=8, the S3-from-P1 one-sided arithmetic, the (E8)_1 character/sector kill-test, and the lattice-VOA second route; mirrors v454/v456/v457/v459. The Casimir + Kac-Moody fits (v454), the Chern flip (v456) and the Tomita-Takesaki tower (v455) are numerical, Python-only; the exact algebra is mirrored here. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2828)
- [(* v463 (i): c=dim/(1+h^v)=8 has three simple solutions; holomorphy forces 248=dim E8 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2910)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* v461: tied obstruction integers + the reflection forcing (strict-local => C=-C => C=0) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2881)
- [(* v459: E8 even unimodular + the 240 roots split 112+128; 248=8+240=120+128 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2866)
- [(* ==== v469 round: the net-level crossed-product certification -- the LR locality integers (h_s(SO(16)_1) = 16/16 = 1 in Z; the mu4 glue h(J^k) = {1,1,1} all integer), the KLM mu-index arithmetic (4/2^2 = 1 and 16/4^2 = 1, both holomorphic) and the Kitaev 16-fold-way pin (nu = 16 x |C| = 16 = 0 mod 16, c_- = nu/2 = 8 = g_car+N_fam); mirrors v469 (the FHS Chern itself is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2924)
- [(* ---- (v125) glue Q-system ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L872)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* v469 (i): the LR locality integers of both extension routes *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2930)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [(* v456: c3's '8' is the one-sided count; reflection C->-C so two-sided forces C=0 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2841)
- [(* ==== v470 round: the alpha^3 level / embedding-index rigidity -- the exact k_Y = 5/3 arithmetic (tr Y^2 / tr T3^2 over the 5bar), the b1 conversion (3/5)(41/6) = 41/10, and the carrier decomposition 41/6 = 20/3 + 1/6; mirrors v470 (the FHS Chern and the alpha^-1 root re-verification are numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2943)
- [(* ==== v463 round: the c=8 holomorphic-uniqueness pin -- c=8 has THREE level-1 candidates (A8/80, D8/120, E8/248), but holomorphy forces dim V_1 = E4/eta^8 q^1 coeff = 248 => E8 uniquely; mirrors v463 (the c8 candidates + the forced 248 + the E4/eta^8 tower). v464 (one-particle realization rigidity) is numerical (symbol convergence), Python-only. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2905)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* v463 (ii): the (E8)_1 tower = E4/eta^8 coefficients {1,248,4124,34752,213126} *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2916)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* v469 (ii): KLM mu-index + the 16-fold-way pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2937)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v461/v462 round: the Kapustin-Fidkowski strict-locality obstruction (integer logic; the Wilson-loop winding/Chern themselves are numerical, Python-only) and the 128-spinor extension at character level -- the Jacobi/E8 theta identity theta2^8+theta3^8+theta4^8=2E4 (so the (E8)_1 character is the SO(16)_1 vacuum+spinor sum) and the 120+128=248 sector + finite 16-Majorana Fock counts; mirrors v461/v462. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2875)

## 2026-06-30 ( v463/ v464 + Lean + Wolfram — the [formula in source] holomorphic-uniqueness pin and the one-particle realisation rigidity: the SEAM.EQUIV.01 residual reduced to entirely certification on both halves)

[documents/changelog/sections/010-source.tex:272-308](../../../../documents/changelog/sections/010-source.tex#L272) · `section:013b22f7d203ce0c7da1` · subsection

The bird's-eye reframing. The recent G-block modules all sharpened the target of the scaling limit; this round attacks the two halves that actually remain. The continuum existence is already citable ( v458 MMST [formula in source] v459 AGT/AMT); what was implicit was (B) why a [formula in source] holomorphic limit must be [formula in source] , and (A) the realisation input [formula in source] that the abstract seam I

Explicit source references:

- [SEAM.EQUIV.ONEPARTICLE.01](../../../../verification/status_ledger.csv#L527) · explicit_claim_identifier · original status at linked ledger row
- [verification/v463_seam_c8_holomorphic_uniqueness.py](../../../../verification/v463_seam_c8_holomorphic_uniqueness.py#L1) · explicit_file_citation
- [verification/v464_seam_oneparticle_rigidity.py](../../../../verification/v464_seam_oneparticle_rigidity.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.UNIQ.01](../../../../verification/status_ledger.csv#L526) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [TfptCarrier.SeamResidualAxiom.holomorphy_selects_e8](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L179) · explicit_lean_symbol
- [verification/v161_one_particle_reduction.py](../../../../verification/v161_one_particle_reduction.py#L1) · explicit_script_identifier
- [verification/v458_seam_mmst_citation_audit.py](../../../../verification/v458_seam_mmst_citation_audit.py#L1) · explicit_script_identifier
- [verification/v450_seam_edge_entanglement.py](../../../../verification/v450_seam_edge_entanglement.py#L1) · explicit_script_identifier
- [verification/v160_seam_gaussianity_from_pf.py](../../../../verification/v160_seam_gaussianity_from_pf.py#L1) · explicit_script_identifier
- [verification/v113_quasifree_kernel.py](../../../../verification/v113_quasifree_kernel.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamResidualAxiom.c8_three_candidates](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L174) · explicit_lean_symbol
- [TfptCarrier.SeamResidualAxiom.residual_arithmetic](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L286) · explicit_lean_symbol
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1) · explicit_file_citation
- [verification/v459_seam_lattice_voa_route.py](../../../../verification/v459_seam_lattice_voa_route.py#L1) · explicit_script_identifier
- [verification/v426_seam_lto_rp_kms.py](../../../../verification/v426_seam_lto_rp_kms.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [(* ==== v463 round: the c=8 holomorphic-uniqueness pin -- c=8 has THREE level-1 candidates (A8/80, D8/120, E8/248), but holomorphy forces dim V_1 = E4/eta^8 q^1 coeff = 248 => E8 uniquely; mirrors v463 (the c8 candidates + the forced 248 + the E4/eta^8 tower). v464 (one-particle realization rigidity) is numerical (symbol convergence), Python-only. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2905)
- [(* ---- (v113) quasi-free kernel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L608)
- [(* v463 (ii): the (E8)_1 tower = E4/eta^8 coefficients {1,248,4124,34752,213126} *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2916)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* v459: E8 even unimodular + the 240 roots split 112+128; 248=8+240=120+128 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2866)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ==== v450/v451 round: SEAM.EQUIV.EDGE.* -- the edge central charge c_-=8 from entanglement (c=1/2 per Majorana) and the Ising minimal-model Kac weights (Cardy tower); mirrors v450_seam_edge_entanglement.py / v451_seam_edge_cardy_tower.py (the numerical fits are Python-only; the c_-=8 assembly and the Kac weights are exact and mirrored here). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2793)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [(* ==== v454/v456/v457/v459 round (G-block): the edge current-algebra/Sugawara c=8, the S3-from-P1 one-sided arithmetic, the (E8)_1 character/sector kill-test, and the lattice-VOA second route; mirrors v454/v456/v457/v459. The Casimir + Kac-Moody fits (v454), the Chern flip (v456) and the Tomita-Takesaki tower (v455) are numerical, Python-only; the exact algebra is mirrored here. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2828)
- [(* v463 (i): c=dim/(1+h^v)=8 has three simple solutions; holomorphy forces 248=dim E8 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2910)

## 2026-06-30 ( v461/ v462 + Lean axiom merge — the strict-locality obstruction made explicit, the 128-spinor extension exhibited at character level, and the SEAM.EQUIV.01 residual reduced to ONE realisation axiom [formula in source] ONE combined cited theorem)

[documents/changelog/sections/010-source.tex:309-347](../../../../documents/changelog/sections/010-source.tex#L309) · `section:9d4f17be173ce81f109e` · subsection

v461_seam_strict_locality.py ( SEAM.S3.LOCALITY.01). Sharpens v460's verdict from a cited statement to an exhibited topological obstruction on the same v367 [formula in source] collar. \ CHERN INTEGER: the FHS Chern is [formula in source] (chiral [formula in source] ), [formula in source] (trivial [formula in source] ), [formula in source] (opposite mass). \ WILSON-LOOP / WANNIER-CENTRE WINDING [formula in source] CH

Explicit source references:

- [SEAM.EQUIV.SPINOR.01](../../../../verification/status_ledger.csv#L525) · explicit_claim_identifier · original status at linked ledger row
- [verification/v462_seam_spinor_continuum.py](../../../../verification/v462_seam_spinor_continuum.py#L1) · explicit_script_identifier
- [FORM.SEAM.RESIDUAL.01](../../../../verification/status_ledger.csv#L523) · explicit_claim_identifier · original status at linked ledger row
- [verification/v460_seam_s3_lto.py](../../../../verification/v460_seam_s3_lto.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamResidualAxiom.seam_realisation_theorem](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L226) · explicit_lean_symbol
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v424_seam_lto_rp_reduction.py](../../../../verification/v424_seam_lto_rp_reduction.py#L1) · explicit_script_identifier
- [verification/v461_seam_strict_locality.py](../../../../verification/v461_seam_strict_locality.py#L1) · explicit_file_citation
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1) · explicit_file_citation
- [SEAM.S3.LOCALITY.01](../../../../verification/status_ledger.csv#L524) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v458_seam_mmst_citation_audit.py](../../../../verification/v458_seam_mmst_citation_audit.py#L1) · explicit_script_identifier
- [verification/v459_seam_lattice_voa_route.py](../../../../verification/v459_seam_lattice_voa_route.py#L1) · explicit_script_identifier
- [verification/wolfram/tfpt_readouts_extension.wl](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1) · explicit_file_citation
- [verification/v456_seam_chirality_from_c3.py](../../../../verification/v456_seam_chirality_from_c3.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamResidualAxiom.fock_sectors](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L153) · explicit_lean_symbol
- [TfptCarrier.SeamResidualAxiom.strict_locality_forbidden](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L138) · explicit_lean_symbol
- [TfptCarrier.SeamResidualAxiom.residual_arithmetic](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L286) · explicit_lean_symbol

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* v462 (ii): (E8)_1 weight-1 = SO(16)_1 vacuum 120 + spinor 128; finite 16-Majorana Fock *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2896)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [(* v461: tied obstruction integers + the reflection forcing (strict-local => C=-C => C=0) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2881)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* v456: c3's '8' is the one-sided count; reflection C->-C so two-sided forces C=0 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2841)
- [(* v462 (i): the Jacobi/E8 theta identity theta2^8+theta3^8+theta4^8 = 2 E4 (numerical at q=1/50) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2887)
- [(* v459: E8 even unimodular + the 240 roots split 112+128; 248=8+240=120+128 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2866)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [(* ==== v454/v456/v457/v459 round (G-block): the edge current-algebra/Sugawara c=8, the S3-from-P1 one-sided arithmetic, the (E8)_1 character/sector kill-test, and the lattice-VOA second route; mirrors v454/v456/v457/v459. The Casimir + Kac-Moody fits (v454), the Chern flip (v456) and the Tomita-Takesaki tower (v455) are numerical, Python-only; the exact algebra is mirrored here. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2828)
- [(* ==== v461/v462 round: the Kapustin-Fidkowski strict-locality obstruction (integer logic; the Wilson-loop winding/Chern themselves are numerical, Python-only) and the 128-spinor extension at character level -- the Jacobi/E8 theta identity theta2^8+theta3^8+theta4^8=2E4 (so the (E8)_1 character is the SO(16)_1 vacuum+spinor sum) and the 120+128=248 sector + finite 16-Majorana Fock counts; mirrors v461/v462. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2875)

## 2026-06-30 ( v460 — the explicit flat-band commuting-projector (LTO) realisation of the S3 input: advancing the SEAM.EQUIV.01 sub-step (i) from abstract existence to an explicit gap-protected quasi-local projector net carrying the [formula in source] reflection and the [formula in source] clock)

[documents/changelog/sections/010-source.tex:348-371](../../../../documents/changelog/sections/010-source.tex#L348) · `section:586aafa793fa108ed39d` · subsection

v460_seam_s3_lto.py ( SEAM.S3.LTO.01). NPW26 ( arXiv:2605.10693) prove the intrinsic Bisognano–Wichmann lemma for commuting-projector (local-topological-order) models with a [formula in source] reflection; v424 reduced the keystone BW residual to that theorem modulo two open sub-steps, of which (ii) was settled by v426 and the continuum lift of BW by the uniform-in- [formula in source] tower v455. This module attacks

Explicit source references:

- [SEAM.S3.LTO.01](../../../../verification/status_ledger.csv#L522) · explicit_claim_identifier · original status at linked ledger row
- [verification/v424_seam_lto_rp_reduction.py](../../../../verification/v424_seam_lto_rp_reduction.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v460_seam_s3_lto.py](../../../../verification/v460_seam_s3_lto.py#L1) · explicit_script_identifier
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/v426_seam_lto_rp_kms.py](../../../../verification/v426_seam_lto_rp_kms.py#L1) · explicit_script_identifier
- [verification/v455_seam_bw_uniform.py](../../../../verification/v455_seam_bw_uniform.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)

## 2026-06-30 ( v454– v459 + Lean + Wolfram — the post-F ``G-block'': the lattice Sugawara current algebra, a uniform-in- [formula in source] Tomita–Takesaki tower, chirality forced by [formula in source] , an [formula in source] character kill test, the exact MMST citation audit and the complementary lattice-VOA route)

[documents/changelog/sections/011-source.tex:1-52](../../../../documents/changelog/sections/011-source.tex#L1) · `section:7ca01292e7152e165e80` · subsection

v454_seam_edge_virasoro.py ( SEAM.EQUIV.EDGE.VIRASORO.01). Tests the cited MMST ``Koo–Saleur lattice generators [formula in source] continuum Virasoro with the correct [formula in source] '' fact on our own free-fermion edge (G2). \ the Affleck–Cardy Casimir energy fits [formula in source] ( [formula in source] per Majorana, [formula in source] ); \ the lattice [formula in source] current closes its affine level exac

Explicit source references:

- [verification/v458_seam_mmst_citation_audit.py](../../../../verification/v458_seam_mmst_citation_audit.py#L1) · explicit_script_identifier
- [verification/v454_seam_edge_virasoro.py](../../../../verification/v454_seam_edge_virasoro.py#L1) · explicit_file_citation
- [SEAM.EQUIV.LATTICEVOA.01](../../../../verification/status_ledger.csv#L521) · explicit_claim_identifier · original status at linked ledger row
- [verification/v459_seam_lattice_voa_route.py](../../../../verification/v459_seam_lattice_voa_route.py#L1) · explicit_file_citation
- [verification/v424_seam_lto_rp_reduction.py](../../../../verification/v424_seam_lto_rp_reduction.py#L1) · explicit_script_identifier
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1) · explicit_file_citation
- [verification/v456_seam_chirality_from_c3.py](../../../../verification/v456_seam_chirality_from_c3.py#L1) · explicit_file_citation
- [verification/v426_seam_lto_rp_kms.py](../../../../verification/v426_seam_lto_rp_kms.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.KILLTEST.01](../../../../verification/status_ledger.csv#L519) · explicit_claim_identifier · original status at linked ledger row
- [verification/v457_seam_character_killtest.py](../../../../verification/v457_seam_character_killtest.py#L1) · explicit_file_citation
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/wolfram/tfpt_readouts_extension.wl](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1) · explicit_file_citation
- [SEAM.EQUIV.MMST.AUDIT.01](../../../../verification/status_ledger.csv#L520) · explicit_claim_identifier · original status at linked ledger row
- [FORM.SEAM.RESIDUAL.01](../../../../verification/status_ledger.csv#L523) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.EDGE.VIRASORO.01](../../../../verification/status_ledger.csv#L516) · explicit_claim_identifier · original status at linked ledger row
- [verification/v455_seam_bw_uniform.py](../../../../verification/v455_seam_bw_uniform.py#L1) · explicit_file_citation
- [SEAM.EQUIV.BW.UNIFORM.01](../../../../verification/status_ledger.csv#L517) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* v459: E8 even unimodular + the 240 roots split 112+128; 248=8+240=120+128 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2866)
- [(* v454: level-1 Sugawara c = dim G/(1+h^v) = 8 for both SO(16)_1 and (E8)_1 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2835)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* v457: the (E8)_1 character tower {1,248,4124,34752} and the kill-test discriminators *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2848)
- [(* v456: c3's '8' is the one-sided count; reflection C->-C so two-sided forces C=0 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2841)
- [(* ==== v454/v456/v457/v459 round (G-block): the edge current-algebra/Sugawara c=8, the S3-from-P1 one-sided arithmetic, the (E8)_1 character/sector kill-test, and the lattice-VOA second route; mirrors v454/v456/v457/v459. The Casimir + Kac-Moody fits (v454), the Chern flip (v456) and the Tomita-Takesaki tower (v455) are numerical, Python-only; the exact algebra is mirrored here. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2828)

## 2026-06-30 ( v449– v453 + Lean + Wolfram — the edge CFT named end-to-end: uniform-in- [formula in source] convergence, three independent reads of [formula in source] , the Ising operator tower, the [formula in source] modular data, and [formula in source] -symmetry derived from the four marks)

[documents/changelog/sections/011-source.tex:53-104](../../../../documents/changelog/sections/011-source.tex#L53) · `section:12f49e008971e532539c` · subsection

v449_seam_mmst_uniform.py ( SEAM.EQUIV.MMST.UNIFORM.01). Strengthens the cited MMST hypothesis from convergence at a single cross-ratio ( v444/ v448) to the uniformity the theorem actually asserts, on the same v367/ v368 collar. \ cross-ratio A [formula in source] , the deviation shrinking monotonically; \ an independent cross-ratio B [formula in source] its own value [formula in source] , so convergence is not a tun

Explicit source references:

- [GATE.WOLFRAM.02](../../../../verification/status_ledger.csv#L168) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/wolfram/tfpt_readouts_extension.wl](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1) · explicit_file_citation
- [verification/v447_seam_edge_chern.py](../../../../verification/v447_seam_edge_chern.py#L1) · explicit_script_identifier
- [verification/v450_seam_edge_entanglement.py](../../../../verification/v450_seam_edge_entanglement.py#L1) · explicit_file_citation
- [SEAM.RIGIDITY.MU4FROMMARKS.01](../../../../verification/status_ledger.csv#L514) · explicit_claim_identifier · original status at linked ledger row
- [verification/v448_seam_edge_fourpoint.py](../../../../verification/v448_seam_edge_fourpoint.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.EDGE.ENTANGLEMENT.01](../../../../verification/status_ledger.csv#L511) · explicit_claim_identifier · original status at linked ledger row
- [verification/v452_seam_e8_modular.py](../../../../verification/v452_seam_e8_modular.py#L1) · explicit_script_identifier
- [FORM.SEAM.EDGE.CHERN.01](../../../../verification/status_ledger.csv#L515) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v446_seam_clock_invariance.py](../../../../verification/v446_seam_clock_invariance.py#L1) · explicit_script_identifier
- [verification/v449_seam_mmst_uniform.py](../../../../verification/v449_seam_mmst_uniform.py#L1) · explicit_script_identifier
- [verification/v156_seam_net_construction.py](../../../../verification/v156_seam_net_construction.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.EDGE.CARDY.01](../../../../verification/status_ledger.csv#L512) · explicit_claim_identifier · original status at linked ledger row
- [verification/v451_seam_edge_cardy_tower.py](../../../../verification/v451_seam_edge_cardy_tower.py#L1) · explicit_script_identifier
- [verification/v444_seam_edge_scaling.py](../../../../verification/v444_seam_edge_scaling.py#L1) · explicit_script_identifier
- [QGEO.REALIZE.01](../../../../verification/status_ledger.csv#L253) · explicit_claim_identifier · original status at linked ledger row
- [verification/v168_qgeo_rigidity.py](../../../../verification/v168_qgeo_rigidity.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.E8.MODULAR.01](../../../../verification/status_ledger.csv#L513) · explicit_claim_identifier · original status at linked ledger row
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1) · explicit_file_citation
- [verification/v453_seam_mu4_from_marks.py](../../../../verification/v453_seam_mu4_from_marks.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.MMST.UNIFORM.01](../../../../verification/status_ledger.csv#L510) · explicit_claim_identifier · original status at linked ledger row
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamEdgeChern.edge_kernel_facts](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L100) · explicit_lean_symbol
- [verification/v445_seam_rigidity_forcing.py](../../../../verification/v445_seam_rigidity_forcing.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v452 round: SEAM.EQUIV.E8.MODULAR.01 -- the (E8)_1 torus modular data (one primary, S-invariant, T-phase e^{-2pi i/3}, holomorphic c=8); mirrors v452_seam_e8_modular.py. EisensteinE is left unevaluated by some engines, so E4 is built from its q-series; the modular T-branch e^{i pi/12} is carried correctly by the native DedekindEta. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2770)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* v445 (iii): the order is the discriminator -- order-2 commutant strictly larger; marks=4 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2759)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [(* ==== v450/v451 round: SEAM.EQUIV.EDGE.* -- the edge central charge c_-=8 from entanglement (c=1/2 per Majorana) and the Ising minimal-model Kac weights (Cardy tower); mirrors v450_seam_edge_entanglement.py / v451_seam_edge_cardy_tower.py (the numerical fits are Python-only; the c_-=8 assembly and the Kac weights are exact and mirrored here). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2793)
- [(* v445 (ii): exact commutant dimension = sum n_s^2, a PROPER subspace of N^2 (N=16) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2753)
- [(* v445 (i): entrywise forcing -- [rho,E_ab]=(i^a-i^b)E_ab=0 <=> a==b mod 4 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2747)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v168) QGEO rigidity: mu4 square cross-ratio, mu4 characters = A3 exponents ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1453)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v453 round: SEAM.RIGIDITY.MU4FROMMARKS.01 -- the mu4-symmetry of the seam RP data (QGEO.SYM.01) derived from the four marks; mirrors v453_seam_mu4_from_marks.py. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2809)
- [(* ==== v445 round: SEAM.RIGIDITY.FORCING.01 -- commuting with the order-4 clock FORCES mu4-block-diagonality (the converse of v442; mirrors v445_seam_rigidity_forcing.py) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2745)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ---- (v156) seam-net construction: DtN = |k| + the (E8)_1 character ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1391)
- [experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean#L1)

## 2026-06-30 ( v446/ v447/ v448 + Lean + Wolfram — deriving clock-invariance from [formula in source] -symmetry, an independent topological reading of the edge [formula in source] , and pinning the one external MMST fact at four-point order)

[documents/changelog/sections/011-source.tex:105-149](../../../../documents/changelog/sections/011-source.tex#L105) · `section:e41d583f251c87e99261` · subsection

v446_seam_clock_invariance.py ( SEAM.RIGIDITY.CLOCKINV.01). Closes v445's single residual one notch further: the operator condition [formula in source] (clock-invariance) is shown to be a consequence of a structural symmetry, not an independent assumption. \ the OS canonical transfer [formula in source] (fixed by RP [formula in source] reflection, v54) of a [formula in source] -invariant RP covariance satisfies [form

Explicit source references:

- [verification/v447_seam_edge_chern.py](../../../../verification/v447_seam_edge_chern.py#L1) · explicit_script_identifier
- [SEAM.RIGIDITY.CLOCKINV.01](../../../../verification/status_ledger.csv#L506) · explicit_claim_identifier · original status at linked ledger row
- [verification/wolfram/tfpt_readouts_extension.wl](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1) · explicit_file_citation
- [SEAM.EQUIV.EDGE.FOURPOINT.01](../../../../verification/status_ledger.csv#L508) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [TfptCarrier.SeamRigidityForcing.forcing_kernel_facts](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean#L84) · explicit_lean_symbol
- [SEAM.EQUIV.EDGE.CHERN.01](../../../../verification/status_ledger.csv#L507) · explicit_claim_identifier · original status at linked ledger row
- [verification/v448_seam_edge_fourpoint.py](../../../../verification/v448_seam_edge_fourpoint.py#L1) · explicit_script_identifier
- [verification/v444_seam_edge_scaling.py](../../../../verification/v444_seam_edge_scaling.py#L1) · explicit_script_identifier
- [verification/v445_seam_rigidity_forcing.py](../../../../verification/v445_seam_rigidity_forcing.py#L1) · explicit_script_identifier
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean#L1) · explicit_file_citation
- [FORM.SEAM.RIGIDITY.FORCING.01](../../../../verification/status_ledger.csv#L509) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v446_seam_clock_invariance.py](../../../../verification/v446_seam_clock_invariance.py#L1) · explicit_file_citation
- [verification/v442_seam_rigidity_uniform.py](../../../../verification/v442_seam_rigidity_uniform.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [GATE.WOLFRAM.02](../../../../verification/status_ledger.csv#L168) · explicit_claim_identifier · original status at linked ledger row
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_script_identifier
- [verification/v54_seam_horizon_keystones.py](../../../../verification/v54_seam_horizon_keystones.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v445 (ii): exact commutant dimension = sum n_s^2, a PROPER subspace of N^2 (N=16) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2753)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ==== v445 round: SEAM.RIGIDITY.FORCING.01 -- commuting with the order-4 clock FORCES mu4-block-diagonality (the converse of v442; mirrors v445_seam_rigidity_forcing.py) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2745)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* v445 (i): entrywise forcing -- [rho,E_ab]=(i^a-i^b)E_ab=0 <=> a==b mod 4 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2747)
- [(* v445 (iii): the order is the discriminator -- order-2 commutant strictly larger; marks=4 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2759)
- [(* ---- (v54) seam=horizon keystones: triply-forced 8 + shared transfer operator ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L235)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)

## 2026-06-30 ( v445 — the rigidity FORCING converse: commuting with the order- [formula in source] clock is exactly [formula in source] -block-diagonality, upgrading the rigidity step from ``permit'' to ``force'')

[documents/changelog/sections/011-source.tex:150-174](../../../../documents/changelog/sections/011-source.tex#L150) · `section:04b7cfb5801c917940c7` · subsection

v445_seam_rigidity_forcing.py ( SEAM.RIGIDITY.FORCING.01). The rigidity route ( v398/ v442) showed a [formula in source] -block-diagonal [formula in source] commutes with the order- [formula in source] clock [formula in source] . This module proves the converse and the exact structure, isolating the single remaining physical premise. \ FORCING (entrywise iff): [formula in source] iff [formula in source] whenever [for

Explicit source references:

- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_script_identifier
- [verification/v89_carrier_index_lemma.py](../../../../verification/v89_carrier_index_lemma.py#L1) · explicit_script_identifier
- [verification/v54_seam_horizon_keystones.py](../../../../verification/v54_seam_horizon_keystones.py#L1) · explicit_script_identifier
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_script_identifier
- [verification/v154_simple_current_theorem.py](../../../../verification/v154_simple_current_theorem.py#L1) · explicit_script_identifier
- [SEAM.RIGIDITY.FORCING.01](../../../../verification/status_ledger.csv#L504) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v445_seam_rigidity_forcing.py](../../../../verification/v445_seam_rigidity_forcing.py#L1) · explicit_script_identifier
- [verification/v398_seam_state_rigidity.py](../../../../verification/v398_seam_state_rigidity.py#L1) · explicit_script_identifier
- [SEAM.RIGIDITY.01](../../../../verification/status_ledger.csv#L468) · explicit_claim_identifier · original status at linked ledger row
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v422_galois_net_bridge.py](../../../../verification/v422_galois_net_bridge.py#L1) · explicit_script_identifier
- [verification/v351_continuum_realisation_sharpened.py](../../../../verification/v351_continuum_realisation_sharpened.py#L1) · explicit_script_identifier
- [verification/v442_seam_rigidity_uniform.py](../../../../verification/v442_seam_rigidity_uniform.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ---- (v54) seam=horizon keystones: triply-forced 8 + shared transfer operator ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L235)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [(* v445 (ii): exact commutant dimension = sum n_s^2, a PROPER subspace of N^2 (N=16) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2753)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [(* ==== v445 round: SEAM.RIGIDITY.FORCING.01 -- commuting with the order-4 clock FORCES mu4-block-diagonality (the converse of v442; mirrors v445_seam_rigidity_forcing.py) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2745)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [(* ==== v422 round: the Galois<->Net bridge (mu4 = Gal(Q(zeta5)) = (E8)_1 glue) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2583)
- [(* ---- (v154) Simple-Current Extension Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1381)
- [(* v445 (iii): the order is the discriminator -- order-2 commutant strictly larger; marks=4 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2759)
- [(* ==== v351 SEAM.EQUIV.CONTINUUM.02: the order-4 clock + Coxeter-match pin E8 over SO(16) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2335)
- [(* ---- (v89) carrier index lemma ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L118)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [(* v445 (i): entrywise forcing -- [rho,E_ab]=(i^a-i^b)E_ab=0 <=> a==b mod 4 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2747)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)

## 2026-06-30 (Lean kernel-hardening of two of the seven SEAM.EQUIV.01 routes — the Borchers/Wiesbrock standard-pair relations and the applicability-ledger count, machine-verified with clean axioms; rigor without new mathematics)

[documents/changelog/sections/011-source.tex:175-202](../../../../documents/changelog/sections/011-source.tex#L175) · `section:fcc9c4aa50365441ef12` · subsection

Context. Two of the seven convergent attacks on the keystone residual are turned from numerical witnesses into kernel-checked facts, in the audit-contract style of SeamScalingLimit.lean / BWKeystone.lean: the load-bearing algebra/bookkeeping is proved by the Lean kernel, while the cited theorems and the one open continuum input remain named axioms. lake build OK; \#print axioms clean. Tracked as two new FORM.* ledger

Explicit source references:

- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1) · explicit_file_citation
- [TfptCarrier.SeamApplicabilityLedger.internal_count](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L77) · explicit_lean_symbol
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1) · explicit_file_citation
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1) · explicit_file_citation
- [TfptCarrier.SeamApplicabilityLedger.mmst_range](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L104) · explicit_lean_symbol
- [TfptCarrier.SeamApplicabilityLedger.external_count](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L80) · explicit_lean_symbol
- [TfptCarrier.SeamApplicabilityLedger.total_count](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L72) · explicit_lean_symbol
- [FORM.SEAM.APPLICABILITY.01](../../../../verification/status_ledger.csv#L505) · explicit_claim_identifier · original status at linked ledger row
- [FORM.SEAM.BW.HSMI.01](../../../../verification/status_ledger.csv#L503) · explicit_claim_identifier · original status at linked ledger row
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1) · explicit_file_citation
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1) · explicit_file_citation
- [verification/v441_seam_applicability_ledger.py](../../../../verification/v441_seam_applicability_ledger.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamStandardPair.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamStandardPair.lean#L1) · explicit_file_citation
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1) · explicit_file_citation
- [verification/v438_seam_hsmi_borchers.py](../../../../verification/v438_seam_hsmi_borchers.py#L1) · explicit_script_identifier
- [TfptCarrier.SeamApplicabilityLedger.npw_outside_bucket](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L107) · explicit_lean_symbol

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamStandardPair.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamStandardPair.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)

## 2026-06-30 (the chiral-edge scaling motor v444 — the edge half of the one SEAM.EQUIV.01 continuum residual: the lattice edge already carries the conformal data of a chiral [formula in source] Majorana CFT and converges under refinement, grounding the EDGE the way v439 grounded the BULK; the keystone stays )

[documents/changelog/sections/011-source.tex:203-225](../../../../documents/changelog/sections/011-source.tex#L203) · `section:cccf91ab17849ec5feae` · subsection

Context. v439 grounded the bulk half of the cited MMST existence (uniform gap, exponential clustering, finite light cone) and explicitly narrowed the one open analytic step to the chiral massless edge scaling limit. v444 attacks exactly that residual on the SAME explicit v367/ v368 [formula in source] collar. vN-registered, runs in run_all.py (suite ends ALL CHECKS PASSED), veri-cited in the body of tfpt_research_con

Explicit source references:

- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.CONTINUUM.EDGE.01](../../../../verification/status_ledger.csv#L502) · explicit_claim_identifier · original status at linked ledger row
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/v439_seam_lieb_robinson.py](../../../../verification/v439_seam_lieb_robinson.py#L1) · explicit_script_identifier
- [verification/v444_seam_edge_scaling.py](../../../../verification/v444_seam_edge_scaling.py#L1) · explicit_file_citation
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)

## 2026-06-30 (six convergent attacks on the one SEAM.EQUIV.01 continuum residual — intrinsic Bisognano–Wichmann, a Lieb–Robinson continuum motor, the [formula in source] -homotopy, an applicability audit, uniform rigidity, and falsifiable kill tests; none closes the keystone, together they reduce it to a single analytic fact and make the claim falsifiable)

[documents/changelog/sections/011-source.tex:226-270](../../../../documents/changelog/sections/011-source.tex#L226) · `section:c6b6da016e7f30bbb7c1` · subsection

Context. The keystone SEAM.EQUIV.01 is ``closed modulo a cited theorem'': its one open input is the continuum existence of the seam's chiral massless scaling limit. Six new modules attack that residual from independent directions; each is honest about not closing it. All six are vN-registered, run in run_all.py (suite ends ALL CHECKS PASSED), and are veri-cited in the body of tfpt_research_contracts (the [formula in

Explicit source references:

- [verification/v90_conical_defect_chain.py](../../../../verification/v90_conical_defect_chain.py#L1) · explicit_script_identifier
- [verification/v439_seam_lieb_robinson.py](../../../../verification/v439_seam_lieb_robinson.py#L1) · explicit_file_citation
- [SEAM.E8.DISCRIMINATOR.01](../../../../verification/status_ledger.csv#L501) · explicit_claim_identifier · original status at linked ledger row
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [verification/v440_seam_lto_rp_beta.py](../../../../verification/v440_seam_lto_rp_beta.py#L1) · explicit_file_citation
- [SEAM.EQUIV.BW.HSMI.01](../../../../verification/status_ledger.csv#L496) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.LTORP.BETA.01](../../../../verification/status_ledger.csv#L498) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.APPLICABILITY.01](../../../../verification/status_ledger.csv#L499) · explicit_claim_identifier · original status at linked ledger row
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v426_seam_lto_rp_kms.py](../../../../verification/v426_seam_lto_rp_kms.py#L1) · explicit_script_identifier
- [verification/v442_seam_rigidity_uniform.py](../../../../verification/v442_seam_rigidity_uniform.py#L1) · explicit_file_citation
- [verification/v441_seam_applicability_ledger.py](../../../../verification/v441_seam_applicability_ledger.py#L1) · explicit_file_citation
- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [verification/v398_seam_state_rigidity.py](../../../../verification/v398_seam_state_rigidity.py#L1) · explicit_script_identifier
- [verification/v443_seam_e8_discriminator.py](../../../../verification/v443_seam_e8_discriminator.py#L1) · explicit_file_citation
- [SEAM.RIGIDITY.UNIFORM.01](../../../../verification/status_ledger.csv#L500) · explicit_claim_identifier · original status at linked ledger row
- [verification/v347_seam_closure_modes.py](../../../../verification/v347_seam_closure_modes.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v438_seam_hsmi_borchers.py](../../../../verification/v438_seam_hsmi_borchers.py#L1) · explicit_file_citation
- [SEAM.EQUIV.CONTINUUM.LR.01](../../../../verification/status_ledger.csv#L497) · explicit_claim_identifier · original status at linked ledger row
- [verification/v424_seam_lto_rp_reduction.py](../../../../verification/v424_seam_lto_rp_reduction.py#L1) · explicit_script_identifier
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/v215_seam_deck_killtest.py](../../../../verification/v215_seam_deck_killtest.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v347 SEAM.EQUIV.CLOSURE.01: the flat seam base vs the spherical icosahedral base ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2298)
- [(* ---- (v90) conical defect chain ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L126)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamStandardPair.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamStandardPair.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean#L1)

## 2026-06-29 (discoverability — an SEO/GEO pass on the public website: a generated /llms.txt, explicit AI-crawler robots.txt, query-targeted metadata, and the first real backlinks from the README and the Zenodo deposit)

[documents/changelog/sections/011-source.tex:299-342](../../../../documents/changelog/sections/011-source.tex#L299) · `section:44fb2c4f56d2211a1e7f` · subsection

New /llms.txt context file for AI / answer engines (ChatGPT, Perplexity, Claude, Gemini), generated at build time from the same data layer as the site ( papers.ts, predictions.ts, version.ts) so it never drifts. It carries the one-line summary, the ``claims / does not claim'' split, the document set with PDF links, the headline predictions, and an explicit disambiguation from the unrelated Brouwer–Lefschetz ``topolog

Explicit source references:

- [verification/v437_e8_degree_joint.py](../../../../verification/v437_e8_degree_joint.py#L1) · explicit_script_identifier
- [README.md](../../../../README.md#L1) · explicit_file_citation
- [verification/v431_e8_degree_ladder.py](../../../../verification/v431_e8_degree_ladder.py#L1) · explicit_script_identifier
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation
- [verification/v355_e8_unmapped_bandwidth.py](../../../../verification/v355_e8_unmapped_bandwidth.py#L1) · explicit_script_identifier
- [verification/wolfram/tfpt_readouts_extension.wl](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1) · explicit_file_citation
- [verification/v66_e8_casimir_degrees.py](../../../../verification/v66_e8_casimir_degrees.py#L1) · explicit_script_identifier
- [E8.DEGREE.JOINT.01](../../../../verification/status_ledger.csv#L99) · explicit_claim_identifier · original status at linked ledger row
- [verification/v6_bootstrap.py](../../../../verification/v6_bootstrap.py#L1) · explicit_script_identifier
- [verification/v354_e8_reverse_audit.py](../../../../verification/v354_e8_reverse_audit.py#L1) · explicit_script_identifier
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v431 (1): two-family split forced by 30 = 2*3*5 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2664)
- [(* v431 (4): 18 = 6*N_fam spine member, not det-ladder holdout *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2690)
- [(* ==== v355 E8.UNMAPPED.BANDWIDTH.01: the unmapped region's forced content is COLLECTIVE ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2363)
- [(* v437 (3): the three mapped degrees are the canonical invariants *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2738)
- [(* v431 (2): 6k family = 6 x spine {2,3,4,5} = v91 anchor a=(1,1,2) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2669)
- [(* ==== v437 round: E8.DEGREE.JOINT.01 -- (g_car,N_fam)=unique roots of x^2-8x+15; E8 fixed by its degrees ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2715)
- [(* ==== v431 round: E8.DEGREE.LADDER.01 -- unmapped degrees = 6*spine U det-ladder ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2650)
- [(* v431 (6): all five 'unmapped' degrees homed; v355 functorial line upheld [P] *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2707)
- [(* v431 (5): clean {0,2} split special to E8 among exceptionals *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2694)
- [(* ---- (v66) compiler atoms = E8 Casimir invariant degrees ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L308)
- [(* v437 (2): the joint quadratic -- (g_car,N_fam) unique roots, coeffs (rank, h/2) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2730)
- [(* v431 (3): 6k+2 family = {2} U det-ladder {8,14,20} on winding line 6s+8 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2679)
- [(* ==== v354 E8.REVERSE.AUDIT.01: 3 of 8 E8 Casimir degrees are load-bearing, 5 unmapped ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2353)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean#L1)
- [(* v437 (1): all of E8 fixed by its degree multiset *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2725)

## 2026-06-27 ( v431, v432 — the ``unmapped'' [formula in source] degrees are forced spine [formula in source] flavor structure; the unconditional over-determination floor)

[documents/changelog/sections/012-source.tex:97-132](../../../../documents/changelog/sections/012-source.tex#L97) · `section:b37031f13445096fcccf` · subsection

New module v431 ( E8.DEGREE.LADDER.01) — the reverse-audit ``unmapped'' [formula in source] Casimir degrees are NOT diffuse overhead but a forced two-family decomposition. The structural complement of the reverse audit ( v354/ v355) and its sheet/deck complement ( v430): [formula in source] . \ the two-family split is forced by [formula in source] — the exponents are the [formula in source] totatives of [formula in s

Explicit source references:

- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [verification/v135_det_surface.py](../../../../verification/v135_det_surface.py#L1) · explicit_script_identifier
- [verification/v100_numerology_null_mc.py](../../../../verification/v100_numerology_null_mc.py#L1) · explicit_script_identifier
- [ALPHA.QUILLEN.EXACT.01](../../../../verification/status_ledger.csv#L445) · explicit_claim_identifier · original status at linked ledger row
- [verification/v431_e8_degree_ladder.py](../../../../verification/v431_e8_degree_ladder.py#L1) · explicit_script_identifier
- [verification/v432_overdet_floor.py](../../../../verification/v432_overdet_floor.py#L1) · explicit_script_identifier
- [verification/wolfram/tfpt_readouts_extension.wl](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1) · explicit_file_citation
- [verification/v91_spine_tetrahedron.py](../../../../verification/v91_spine_tetrahedron.py#L1) · explicit_script_identifier
- [verification/v427_overdet_witness_map.py](../../../../verification/v427_overdet_witness_map.py#L1) · explicit_script_identifier
- [verification/v428_overdet_reclass.py](../../../../verification/v428_overdet_reclass.py#L1) · explicit_script_identifier
- [verification/v430_other_side_reverse_audit.py](../../../../verification/v430_other_side_reverse_audit.py#L1) · explicit_script_identifier
- [E8.DEGREE.LADDER.01](../../../../verification/status_ledger.csv#L98) · explicit_claim_identifier · original status at linked ledger row
- [verification/v382_alpha_quillen_exact.py](../../../../verification/v382_alpha_quillen_exact.py#L1) · explicit_script_identifier
- [OVERDET.FLOOR.01](../../../../verification/status_ledger.csv#L100) · explicit_claim_identifier · original status at linked ledger row
- [verification/v187_ftransfer_laws.py](../../../../verification/v187_ftransfer_laws.py#L1) · explicit_script_identifier
- [verification/v354_e8_reverse_audit.py](../../../../verification/v354_e8_reverse_audit.py#L1) · explicit_script_identifier
- [verification/v355_e8_unmapped_bandwidth.py](../../../../verification/v355_e8_unmapped_bandwidth.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v431 (5): clean {0,2} split special to E8 among exceptionals *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2694)
- [(* v431 (2): 6k family = 6 x spine {2,3,4,5} = v91 anchor a=(1,1,2) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2669)
- [(* v431 (1): two-family split forced by 30 = 2*3*5 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2664)
- [(* v430 (2): the two sheets are the 128-spinor (matched, collective) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2639)
- [(* ==== v355 E8.UNMAPPED.BANDWIDTH.01: the unmapped region's forced content is COLLECTIVE ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2363)
- [(* v431 (3): 6k+2 family = {2} U det-ladder {8,14,20} on winding line 6s+8 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2679)
- [(* ==== v431 round: E8.DEGREE.LADDER.01 -- unmapped degrees = 6*spine U det-ladder ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2650)
- [(* ---- (v91) spine tetrahedron ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L137)
- [(* v430 (1): the deck is the matched degree-2 invariant, not one of the unmapped degrees *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2634)
- [(* ==== v430 round: the 'other side' (double-cover deck / conjugate sheet S^-) vs the unmapped E8 degrees ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2625)
- [(* v431 (4): 18 = 6*N_fam spine member, not det-ladder holdout *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2690)
- [(* v430 (3): the sheet/deck set is forced-disjoint from the unmapped region *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2644)
- [(* ==== v354 E8.REVERSE.AUDIT.01: 3 of 8 E8 Casimir degrees are load-bearing, 5 unmapped ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2353)
- [(* ---- (v135) determinant surface ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1027)
- [(* v431 (6): all five 'unmapped' degrees homed; v355 functorial line upheld [P] *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2707)

## 2026-06-27 ( v430 — the seam's ``other side'' is forced-disjoint from [formula in source] 's unmapped Casimir degrees)

[documents/changelog/sections/012-source.tex:133-161](../../../../documents/changelog/sections/012-source.tex#L133) · `section:f9ed8e942e9df8db94c7` · subsection

New module v430 ( E8.OTHERSIDE.AUDIT.01) — the sheet/deck complement of the reverse audit ( v354/ v355). The reverse audit found five [formula in source] Casimir degrees [formula in source] with no primary readout (``hull overhead''); the double cover has an other side (the one-sided [formula in source] collar / conjugate half-spinor [formula in source] ), so the sharp adversarial follow-up is ``is the leftover where

Explicit source references:

- [verification/v354_e8_reverse_audit.py](../../../../verification/v354_e8_reverse_audit.py#L1) · explicit_script_identifier
- [verification/v227_degree_exponent_channel_split.py](../../../../verification/v227_degree_exponent_channel_split.py#L1) · explicit_script_identifier
- [verification/v355_e8_unmapped_bandwidth.py](../../../../verification/v355_e8_unmapped_bandwidth.py#L1) · explicit_script_identifier
- [verification/v58_seam_horizon_chain.py](../../../../verification/v58_seam_horizon_chain.py#L1) · explicit_script_identifier
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [verification/v430_other_side_reverse_audit.py](../../../../verification/v430_other_side_reverse_audit.py#L1) · explicit_script_identifier
- [E8.OTHERSIDE.AUDIT.01](../../../../verification/status_ledger.csv#L97) · explicit_claim_identifier · original status at linked ledger row
- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v355 E8.UNMAPPED.BANDWIDTH.01: the unmapped region's forced content is COLLECTIVE ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2363)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [(* v430 (3): the sheet/deck set is forced-disjoint from the unmapped region *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2644)
- [(* ==== v354 E8.REVERSE.AUDIT.01: 3 of 8 E8 Casimir degrees are load-bearing, 5 unmapped ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2353)
- [(* v227 degree/exponent channel split: 248 = 120 + 128 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1814)
- [(* ---- (v58) seam-horizon chain: one-sided Gauss-Bonnet, KMS, S_BH, RT in seam units ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L270)
- [(* v430 (1): the deck is the matched degree-2 invariant, not one of the unmapped degrees *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2634)
- [(* ==== v430 round: the 'other side' (double-cover deck / conjugate sheet S^-) vs the unmapped E8 degrees ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2625)
- [(* v430 (2): the two sheets are the 128-spinor (matched, collective) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2639)

## 2026-06-25 (External-works references appendix — web-verified)

[documents/changelog/sections/012-source.tex:189-206](../../../../documents/changelog/sections/012-source.tex#L189) · `section:0b1a0250b55a22849072` · subsection

New shared appendix tex-artefacts/tfpt_references.tex (``External works built upon''), input into the six math-heavy documents ( tfpt_1, tfpt_2, tfpt_4, tfpt_horizon_readouts, origin_theory, tfpt_research_contracts). The active papers named the foundational results they build on by author only, with no reference list; the appendix supplies web-verified locators (journal, volume, year / arXiv) grouped by topic: modula

Explicit source references:

- [verification/make_manifest.py](../../../../verification/make_manifest.py#L1) · explicit_file_citation

## 2026-06-25 (External AQFT citations [formula in source] the four-fold [formula in source] identity sharpened)

[documents/changelog/sections/012-source.tex:207-227](../../../../documents/changelog/sections/012-source.tex#L207) · `section:a9b1f3ad9bb703cf14c6` · subsection

Two recent Adamo–Giorgetti–Tanimoto theorems cited where the keystone builds on them. The open SEAM.EQUIV.01 continuum/realisation legs ( v336) now cite, alongside MMST and Adamo–Moriwaki–Tanimoto ( arXiv:2407.18222, already in use): (a) the construction and classification of two-dimensional conformal nets from even lattices ( arXiv:2506.01008, 2025) — ``the lattice is a conformal net'' is a theorem of the same opera

Explicit source references:

- [SEAM.EQUIV.CONTINUUM.01](../../../../verification/status_ledger.csv#L438) · explicit_claim_identifier · original status at linked ledger row
- [verification/v215_seam_deck_killtest.py](../../../../verification/v215_seam_deck_killtest.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v424_seam_lto_rp_reduction.py](../../../../verification/v424_seam_lto_rp_reduction.py#L1) · explicit_script_identifier
- [verification/v422_galois_net_bridge.py](../../../../verification/v422_galois_net_bridge.py#L1) · explicit_script_identifier
- [verification/v428_overdet_reclass.py](../../../../verification/v428_overdet_reclass.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.GALOIS.NET.01](../../../../verification/status_ledger.csv#L489) · explicit_claim_identifier · original status at linked ledger row
- [verification/v429_axion_pentagon_phi.py](../../../../verification/v429_axion_pentagon_phi.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v422 round: the Galois<->Net bridge (mu4 = Gal(Q(zeta5)) = (E8)_1 glue) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2583)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* v429 (3): the golden character is unique to g_car=5 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2618)
- [(* ==== v429 round: theta_i = the carrier-pentagon interior angle = golden ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2603)
- [(* v429 (1): theta_i = interior angle of the regular g_car-gon (pentagon) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2607)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [(* v429 (2): cos(theta_i) is golden *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2612)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)

## 2026-06-25 ( v428 — honest self-correction of the over-determination map)

[documents/changelog/sections/012-source.tex:228-248](../../../../documents/changelog/sections/012-source.tex#L228) · `section:3ff315fdd4e19bf322ac` · subsection

v428 ( OVERDET.WITNESS.RECLASS.01) — applying v427's own multiply-vs-compress test to v427 itself. A pattern search and the right objection exposed that calling the seven arithmetic witnesses ``disjoint grammars that multiply'' overstates independence. \ by the Brieskorn classification ( v236) the [formula in source] singularity is the ONE generator behind the order- [formula in source] clock (Milnor [formula in sour

Explicit source references:

- [verification/v313_golden_atoms_spectral.py](../../../../verification/v313_golden_atoms_spectral.py#L1) · explicit_script_identifier
- [verification/v428_overdet_reclass.py](../../../../verification/v428_overdet_reclass.py#L1) · explicit_script_identifier
- [verification/v236_brieskorn_capstone.py](../../../../verification/v236_brieskorn_capstone.py#L1) · explicit_script_identifier
- [verification/v427_overdet_witness_map.py](../../../../verification/v427_overdet_witness_map.py#L1) · explicit_script_identifier
- [verification/v305_witness_independence.py](../../../../verification/v305_witness_independence.py#L1) · explicit_script_identifier
- [OVERDET.WITNESS.MAP.01](../../../../verification/status_ledger.csv#L494) · explicit_claim_identifier · original status at linked ledger row
- [OVERDET.WITNESS.RECLASS.01](../../../../verification/status_ledger.csv#L495) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v236: the (2,3,5) Brieskorn singularity generates the skeleton ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1975)
- [(* v313 GOLD.ATOMS.01 -- the golden ratio is the g_car=5 spectral signature *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2185)
- [(* ==== v313-v320: the cyclotomic / Galois arithmetic arc (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2183)

## 2026-06-25 (Three closure-strategy modules v425– v427 [formula in source] the Safeguards paper)

[documents/changelog/sections/012-source.tex:268-304](../../../../documents/changelog/sections/012-source.tex#L268) · `section:f6d50ad015506eb0de3c` · subsection

v425 ( DYN.TRANSFER.UNIVERSAL.01) — the single-flow reduction of [formula in source] . The four frontier transfers do not each need a separate external solver: their dynamics is the ONE native flow (the seam modular/recovery semigroup, v238/ v240, gap [formula in source] ) restricted to each sector. \ the v99 Koide generator linearised at [formula in source] has rate [formula in source] , time- [formula in source] co

Explicit source references:

- [verification/v383_dynamics_universal.py](../../../../verification/v383_dynamics_universal.py#L1) · explicit_script_identifier
- [verification/v427_overdet_witness_map.py](../../../../verification/v427_overdet_witness_map.py#L1) · explicit_script_identifier
- [verification/v238_modular_lindblad_dynamics.py](../../../../verification/v238_modular_lindblad_dynamics.py#L1) · explicit_script_identifier
- [verification/v99_koide_flow_time.py](../../../../verification/v99_koide_flow_time.py#L1) · explicit_script_identifier
- [verification/make_docs_map.py](../../../../verification/make_docs_map.py#L1) · explicit_file_citation
- [verification/v305_witness_independence.py](../../../../verification/v305_witness_independence.py#L1) · explicit_script_identifier
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/make_manifest.py](../../../../verification/make_manifest.py#L1) · explicit_file_citation
- [verification/v425_dyn_transfer_universal.py](../../../../verification/v425_dyn_transfer_universal.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.LTORP.KMS.01](../../../../verification/status_ledger.csv#L493) · explicit_claim_identifier · original status at linked ledger row
- [verification/v426_seam_lto_rp_kms.py](../../../../verification/v426_seam_lto_rp_kms.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v424_seam_lto_rp_reduction.py](../../../../verification/v424_seam_lto_rp_reduction.py#L1) · explicit_script_identifier
- [OVERDET.WITNESS.MAP.01](../../../../verification/status_ledger.csv#L494) · explicit_claim_identifier · original status at linked ledger row
- [DYN.TRANSFER.UNIVERSAL.01](../../../../verification/status_ledger.csv#L492) · explicit_claim_identifier · original status at linked ledger row
- [verification/v320_galois_cp_relation.py](../../../../verification/v320_galois_cp_relation.py#L1) · explicit_script_identifier
- [verification/v240_gns_os_reconstruction.py](../../../../verification/v240_gns_os_reconstruction.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ---- (v99) Koide flow time ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L304)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)
- [(* v320 GALOIS.CP.PREDICT.01 -- the falsifiable CP lock *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2216)

## 2026-06-25 (The (LTO-RP) reduction of the keystone SEAM.EQUIV.01: v424)

[documents/changelog/sections/012-source.tex:305-328](../../../../documents/changelog/sections/012-source.tex#L305) · `section:debae20192b75c7a987d` · subsection

New module v424 ( SEAM.EQUIV.LTORP.01). An honest MMST-style reduction of the one open keystone — SEAM.EQUIV.01's intrinsic Bisognano–Wichmann lemma (the v308/ v309 residual) — to the recent commuting-projector theorem of Naaijkens–Penneys–Wallick ( arXiv:2605.10693, [NPW26]), with the central mapping corrected. \ the [formula in source] clock is a modular symmetry ( [formula in source] , v309/ v199), not the flow —

Explicit source references:

- [verification/v89_carrier_index_lemma.py](../../../../verification/v89_carrier_index_lemma.py#L1) · explicit_script_identifier
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v308_seam_equiv_chain.py](../../../../verification/v308_seam_equiv_chain.py#L1) · explicit_script_identifier
- [verification/v309_modular_bedrock.py](../../../../verification/v309_modular_bedrock.py#L1) · explicit_script_identifier
- [verification/v199_seam_state_invariance.py](../../../../verification/v199_seam_state_invariance.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.LTORP.01](../../../../verification/status_ledger.csv#L491) · explicit_claim_identifier · original status at linked ledger row
- [verification/v424_seam_lto_rp_reduction.py](../../../../verification/v424_seam_lto_rp_reduction.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ---- (v89) carrier index lemma ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L118)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [(* ---- (v199) QGEO.STATE.01: [rho,H]=0 <=> H is mu4-character-block-diagonal ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1649)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)

## 2026-06-25 (The Galois [formula in source] Net bridge: v422)

[documents/changelog/sections/012-source.tex:352-372](../../../../documents/changelog/sections/012-source.tex#L352) · `section:4defd98603f63d26dd3f` · subsection

New module v422 ( SEAM.GALOIS.NET.01). The seam [formula in source] is the SAME cyclic [formula in source] as the [formula in source] simple-current glue — not a mere order- [formula in source] coincidence — bridging the Galois gearbox of v419 to [formula in source] / SEAM.EQUIV.01. \ [formula in source] (one Smith invariant factor [formula in source] [formula in source] cyclic, [formula in source] ); [formula in sou

Explicit source references:

- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v419_seam_galois_carrier.py](../../../../verification/v419_seam_galois_carrier.py#L1) · explicit_script_identifier
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [verification/v89_carrier_index_lemma.py](../../../../verification/v89_carrier_index_lemma.py#L1) · explicit_script_identifier
- [verification/v422_galois_net_bridge.py](../../../../verification/v422_galois_net_bridge.py#L1) · explicit_script_identifier
- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation
- [SEAM.GALOIS.NET.01](../../../../verification/status_ledger.csv#L489) · explicit_claim_identifier · original status at linked ledger row
- [verification/v154_simple_current_theorem.py](../../../../verification/v154_simple_current_theorem.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v419 (1): 30 squarefree => no order-4 on the cyclic clock *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2564)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [(* v419 (3): the seam mu4 = Gal(Q(zeta5)) via the explicit Frobenius G *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2576)
- [(* v419 (2): (Z/30)^x = mu4 x Z2 = (Z/5)^x x (Z/3)^x *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2569)
- [(* ==== v422 round: the Galois<->Net bridge (mu4 = Gal(Q(zeta5)) = (E8)_1 glue) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2583)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [(* ---- (v154) Simple-Current Extension Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1381)
- [(* ==== v419 round: the seam mu4 = Gal(Q(zeta5)) = (Z/5)^x ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2557)
- [(* ---- (v89) carrier index lemma ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L118)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)

## 2026-06-25 (F_pole external: the Koide source [formula in source] pole transfer is one-loop EM: v420)

[documents/changelog/sections/013-source.tex:1-19](../../../../documents/changelog/sections/013-source.tex#L1) · `section:8ac4310a6a4cf6c449f1` · subsection

New module v420 ( FR.POLE.QED.01). The honest external-physics follow-up to v93 and the clock No-Go: the Koide source [formula in source] pole transfer is a one-loop electromagnetic effect, not an algebraic step. \ Koide [formula in source] is rescaling-invariant ( [formula in source] ), so flavor-universal running leaves it fixed and the source [formula in source] pole shift is purely flavor-split; \ the tree defici

Explicit source references:

- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation
- [verification/v187_ftransfer_laws.py](../../../../verification/v187_ftransfer_laws.py#L1) · explicit_script_identifier
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [FR.POLE.QED.01](../../../../verification/status_ledger.csv#L488) · explicit_claim_identifier · original status at linked ledger row
- [verification/v327_hypergraph_rewrite.py](../../../../verification/v327_hypergraph_rewrite.py#L1) · explicit_script_identifier
- [verification/v420_fpole_qed_running.py](../../../../verification/v420_fpole_qed_running.py#L1) · explicit_script_identifier
- [verification/v93_koide_relaxation_toy.py](../../../../verification/v93_koide_relaxation_toy.py#L1) · explicit_script_identifier
- [verification/v82_koide_attractor_splitting.py](../../../../verification/v82_koide_attractor_splitting.py#L1) · explicit_script_identifier
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/v371_ftransfer_pole.py](../../../../verification/v371_ftransfer_pole.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v327 HYP.REWRITE.01: the minimal 3-channel/1-absorbing rewrite M has spectrum {0,2/3,1} *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2227)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)
- [(* ==== v325/v327: the keystone pillowcase orbifold + the cusp-weight rewrite (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2221)
- [(* ---- (v93) Koide relaxation toy ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L166)
- [(* v327 NON-GRAPH-SPECTRAL: 2/3 is NOT a root of the affine-E8 charpoly (sharpens v312) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2232)

## 2026-06-25 (The seam [formula in source] is the carrier's Galois group: v419 [formula in source] figure)

[documents/changelog/sections/013-source.tex:20-44](../../../../documents/changelog/sections/013-source.tex#L20) · `section:36246a833c23944e824a` · subsection

New module v419 ( SEAM.GALOIS.01). A positive resolution of the cyclic/Galois asymmetry ( v409, RES.COXETER.SYMMETRY.01): why the seam [formula in source] sits on the Galois side, not on the cyclic [formula in source] dial. \ [formula in source] is squarefree, so [formula in source] has element orders [formula in source] — no order [formula in source] , so the seam [formula in source] cannot be a rotation hand and is

Explicit source references:

- [verification/v418_cyclotomic_norm_triple.py](../../../../verification/v418_cyclotomic_norm_triple.py#L1) · explicit_script_identifier
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [RES.COXETER.SYMMETRY.01](../../../../verification/status_ledger.csv#L465) · explicit_claim_identifier · original status at linked ledger row
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [SEAM.GALOIS.01](../../../../verification/status_ledger.csv#L487) · explicit_claim_identifier · original status at linked ledger row
- [verification/v316_galois_readout.py](../../../../verification/v316_galois_readout.py#L1) · explicit_script_identifier
- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation
- [verification/make_figures.py](../../../../verification/make_figures.py#L1) · explicit_file_citation
- [verification/v409_coxeter_prime2_lemma.py](../../../../verification/v409_coxeter_prime2_lemma.py#L1) · explicit_script_identifier
- [verification/v419_seam_galois_carrier.py](../../../../verification/v419_seam_galois_carrier.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v316 GALOIS.READOUT.01 -- CP phases are the family-factor cyclotomic data *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2206)
- [(* v418 (1): the carrier-5 clock C5 + the golden ratio *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2540)
- [(* v419 (1): 30 squarefree => no order-4 on the cyclic clock *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2564)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [(* ==== v419 round: the seam mu4 = Gal(Q(zeta5)) = (Z/5)^x ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2557)
- [(* v418 (3): the triple + negative controls *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2549)
- [(* v419 (3): the seam mu4 = Gal(Q(zeta5)) via the explicit Frobenius G *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2576)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ==== v418 round: the cyclotomic norm triple (7,13,55) over the atom-rings ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2533)
- [(* v418 (2): the carrier norm = the quark numerator *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2545)
- [(* v419 (2): (Z/30)^x = mu4 x Z2 = (Z/5)^x x (Z/3)^x *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2569)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean#L1)

## 2026-06-25 (The cyclotomic norm triple [formula in source] the carrier-5 clock: v418)

[documents/changelog/sections/013-source.tex:45-67](../../../../documents/changelog/sections/013-source.tex#L45) · `section:4061f274c55499fb2ed1` · subsection

New module v418 ( ARITH.NORMTRIPLE.01). The carrier split [formula in source] read in the three atom-rings gives [formula in source] (scalaron, [formula in source] , quark numerator), and the carrier- [formula in source] clock missing from v417 is found as the [formula in source] [formula in source] -companion [formula in source] . \ [formula in source] , [formula in source] , [formula in source] (the four primitive

Explicit source references:

- [verification/v418_cyclotomic_norm_triple.py](../../../../verification/v418_cyclotomic_norm_triple.py#L1) · explicit_script_identifier
- [verification/v222_cm_norm_duality.py](../../../../verification/v222_cm_norm_duality.py#L1) · explicit_script_identifier
- [verification/v417_eisenstein_cp_operator.py](../../../../verification/v417_eisenstein_cp_operator.py#L1) · explicit_script_identifier
- [ARITH.NORMTRIPLE.01](../../../../verification/status_ledger.csv#L486) · explicit_claim_identifier · original status at linked ledger row
- [verification/v14_carrier_uniqueness.py](../../../../verification/v14_carrier_uniqueness.py#L1) · explicit_script_identifier
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [verification/v416_atom_trichotomy.py](../../../../verification/v416_atom_trichotomy.py#L1) · explicit_script_identifier
- [verification/v411_ud_ratio_vpower.py](../../../../verification/v411_ud_ratio_vpower.py#L1) · explicit_script_identifier
- [verification/v313_golden_atoms_spectral.py](../../../../verification/v313_golden_atoms_spectral.py#L1) · explicit_script_identifier
- [verification/v410_sheet_generator_binary.py](../../../../verification/v410_sheet_generator_binary.py#L1) · explicit_script_identifier
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/v349_raw_seam_golden_test.py](../../../../verification/v349_raw_seam_golden_test.py#L1) · explicit_script_identifier
- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v349 SEAM.EQUIV.GOLDEN.01: the raw carrier is NOT golden (D5 h=8 -> sqrt2); 5|h only on the output side ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2315)
- [(* v222 CM-norm duality: 41 (square), 7 (hex), 13 (square) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1794)
- [(* v410: closed form V^n, binary spine, bilinear families, theta cross-link *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2403)
- [(* v417 (2): the CP clock W = -P^2 and both CP phases *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2518)
- [(* v418 (2): the carrier norm = the quark numerator *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2545)
- [(* ==== v417 round: the Eisenstein/CP operator + the order-30 clock ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2502)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean#L1)
- [(* v416 (2): the Z[omega] (hex/flavor) trichotomy + ramified prime = own atom *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2490)
- [(* v418 (1): the carrier-5 clock C5 + the golden ratio *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2540)
- [(* v417 (1): the family rotation = Eisenstein deck + the hex norm 7 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2509)
- [(* v411: u/d ratio as a pure V-power readout *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2418)
- [(* v417 (3): one order-30 Coxeter clock; seam mu4 is the Galois side *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2525)
- [(* ==== v313-v320: the cyclotomic / Galois arithmetic arc (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2183)
- [(* ==== v418 round: the cyclotomic norm triple (7,13,55) over the atom-rings ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2533)
- [(* v313 GOLD.ATOMS.01 -- the golden ratio is the g_car=5 spectral signature *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2185)
- [(* ==== v410-v414 round: the sheet generator V as a binary internal compiler ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2394)
- [(* v418 (3): the triple + negative controls *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2549)
- [(* v416 (1): the Z[i] (square/seam) trichotomy *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2485)
- [(* v416 (3): the three facets ramify over one atom each (v403) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2496)

## 2026-06-25 (The Eisenstein/CP operator [formula in source] the order-30 clock: v417)

[documents/changelog/sections/013-source.tex:68-90](../../../../documents/changelog/sections/013-source.tex#L68) · `section:b3159de77fcabb12cc0c` · subsection

New module v417 ( DIAMOND.EISEN.01). The hexagonal dual of v415: the [formula in source] flavor side put on operators and unified with the seam into one order- [formula in source] Coxeter clock. \ the family rotation [formula in source] is the order- [formula in source] Eisenstein deck ( [formula in source] , so [formula in source] on the plane [formula in source] ), the order- [formula in source] dual of [formula in

Explicit source references:

- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/v417_eisenstein_cp_operator.py](../../../../verification/v417_eisenstein_cp_operator.py#L1) · explicit_script_identifier
- [verification/v320_galois_cp_relation.py](../../../../verification/v320_galois_cp_relation.py#L1) · explicit_script_identifier
- [verification/v313_golden_atoms_spectral.py](../../../../verification/v313_golden_atoms_spectral.py#L1) · explicit_script_identifier
- [verification/v316_galois_readout.py](../../../../verification/v316_galois_readout.py#L1) · explicit_script_identifier
- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation
- [verification/v415_gaussian_operator.py](../../../../verification/v415_gaussian_operator.py#L1) · explicit_script_identifier
- [DIAMOND.EISEN.01](../../../../verification/status_ledger.csv#L485) · explicit_claim_identifier · original status at linked ledger row
- [verification/v349_raw_seam_golden_test.py](../../../../verification/v349_raw_seam_golden_test.py#L1) · explicit_script_identifier
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [verification/v233_cp_triality_phase.py](../../../../verification/v233_cp_triality_phase.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v233: CP = the universal family/triality phase, sheet-split (mu6 = mu3 x mu2) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1916)
- [(* v415 (4): the intrinsic order-4 mu4 deck D = -I + J - J^2 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2477)
- [(* v313 GOLD.ATOMS.01 -- the golden ratio is the g_car=5 spectral signature *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2185)
- [(* v320 GALOIS.CP.PREDICT.01 -- the falsifiable CP lock *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2216)
- [(* ==== v417 round: the Eisenstein/CP operator + the order-30 clock ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2502)
- [(* v415 (3): the square-norm dictionary as operator identities *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2470)
- [(* v415 (2): the Gaussian integers 3+2i, 5+4i are eigenvalues *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2464)
- [(* v417 (3): one order-30 Coxeter clock; seam mu4 is the Galois side *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2525)
- [(* ==== v415-v416 round: Gaussian operator + atom trichotomy ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2448)
- [(* v415 (1): J = [U,V]/3 is an integer mu4 quarter-turn complex structure *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2457)
- [(* v417 (2): the CP clock W = -P^2 and both CP phases *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2518)
- [(* ==== v349 SEAM.EQUIV.GOLDEN.01: the raw carrier is NOT golden (D5 h=8 -> sqrt2); 5|h only on the output side ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2315)
- [(* v316 GALOIS.READOUT.01 -- CP phases are the family-factor cyclotomic data *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2206)
- [(* ==== v313-v320: the cyclotomic / Galois arithmetic arc (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2183)
- [(* v417 (1): the family rotation = Eisenstein deck + the hex norm 7 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2509)

## 2026-06-25 (The Gaussian operator [formula in source] the atom trichotomy: v415, v416)

[documents/changelog/sections/013-source.tex:91-121](../../../../documents/changelog/sections/013-source.tex#L91) · `section:bb2a02a767e560376342` · subsection

New module v415 ( DIAMOND.GAUSS.01). The square CM-norm dictionary of v222/ v230 is operator-realised, and the [formula in source] deck becomes an intrinsic diamond object. \ the integer commutator [formula in source] has [formula in source] , [formula in source] on the rank- [formula in source] image and [formula in source] — a [formula in source] quarter-turn born from the binary [formula in source] and ternary [fo

Explicit source references:

- [ARITH.TRICHOTOMY.01](../../../../verification/status_ledger.csv#L484) · explicit_claim_identifier · original status at linked ledger row
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation
- [verification/v410_sheet_generator_binary.py](../../../../verification/v410_sheet_generator_binary.py#L1) · explicit_script_identifier
- [verification/script_registry.csv](../../../../verification/script_registry.csv#L1) · explicit_file_citation
- [verification/v222_cm_norm_duality.py](../../../../verification/v222_cm_norm_duality.py#L1) · explicit_script_identifier
- [verification/v403_facet_compositum.py](../../../../verification/v403_facet_compositum.py#L1) · explicit_script_identifier
- [verification/v95_centered_diamond.py](../../../../verification/v95_centered_diamond.py#L1) · explicit_script_identifier
- [verification/wolfram/tfpt_readouts_extension.wl](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1) · explicit_file_citation
- [verification/v140_canonical_map.py](../../../../verification/v140_canonical_map.py#L1) · explicit_script_identifier
- [verification/v412_sheet_source_corner_J.py](../../../../verification/v412_sheet_source_corner_J.py#L1) · explicit_script_identifier
- [DIAMOND.GAUSS.01](../../../../verification/status_ledger.csv#L483) · explicit_claim_identifier · original status at linked ledger row
- [verification/v416_atom_trichotomy.py](../../../../verification/v416_atom_trichotomy.py#L1) · explicit_script_identifier
- [GATE.WOLFRAM.02](../../../../verification/status_ledger.csv#L168) · explicit_claim_identifier · original status at linked ledger row
- [verification/v415_gaussian_operator.py](../../../../verification/v415_gaussian_operator.py#L1) · explicit_script_identifier
- [verification/v230_center_budget_norms.py](../../../../verification/v230_center_budget_norms.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v415 (3): the square-norm dictionary as operator identities *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2470)
- [(* v416 (2): the Z[omega] (hex/flavor) trichotomy + ramified prime = own atom *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2490)
- [(* ---- (v140) canonical map ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1123)
- [(* v415 (2): the Gaussian integers 3+2i, 5+4i are eigenvalues *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2464)
- [(* v415 (4): the intrinsic order-4 mu4 deck D = -I + J - J^2 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2477)
- [(* v415 (1): J = [U,V]/3 is an integer mu4 quarter-turn complex structure *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2457)
- [(* v230 center budget (7,11,13) = three local norms *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1839)
- [(* ---- (v95) centered flavor diamond ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L205)
- [(* v222 CM-norm duality: 41 (square), 7 (hex), 13 (square) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1794)
- [(* ==== v410-v414 round: the sheet generator V as a binary internal compiler ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2394)
- [(* v410: closed form V^n, binary spine, bilinear families, theta cross-link *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2403)
- [(* ==== v415-v416 round: Gaussian operator + atom trichotomy ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2448)
- [(* v416 (3): the three facets ramify over one atom each (v403) *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2496)
- [(* v412: the sheet source corner J = M(1,-2) on the Z2 wall *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2424)
- [(* v416 (1): the Z[i] (square/seam) trichotomy *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2485)

## 2026-06-24 (Two new tracked targets + axion-branch consistency: v408, v409)

[documents/changelog/sections/013-source.tex:159-192](../../../../documents/changelog/sections/013-source.tex#L159) · `section:cb81efd922f7b6f4776f` · subsection

New module v408 ( HYP.PHI0.PUNCTURE.01). Names the [formula in source] puncture term [formula in source] as the order- [formula in source] local boundary heat-kernel (Seeley–DeWitt) contact term of the [formula in source] puncture divisor, summed over [formula in source] admissible states — the analytic complement of the icosahedral tree term ( v396). \ the exact split [formula in source] , the puncture count, the cl

Explicit source references:

- [verification/v408_phi0_puncture_heatkernel.py](../../../../verification/v408_phi0_puncture_heatkernel.py#L1) · explicit_script_identifier
- [verification/v409_coxeter_prime2_lemma.py](../../../../verification/v409_coxeter_prime2_lemma.py#L1) · explicit_script_identifier
- [ARCH.QUAD.01](../../../../verification/status_ledger.csv#L189) · explicit_claim_identifier · original status at linked ledger row
- [verification/v185_axion_relic_solver.py](../../../../verification/v185_axion_relic_solver.py#L1) · explicit_script_identifier
- [HYP.PHI0.PUNCTURE.01](../../../../verification/status_ledger.csv#L464) · explicit_claim_identifier · original status at linked ledger row
- [verification/v396_phi0_icosahedral.py](../../../../verification/v396_phi0_icosahedral.py#L1) · explicit_script_identifier
- [FORM.COXETER.PRIME2.01](../../../../verification/status_ledger.csv#L466) · explicit_claim_identifier · original status at linked ledger row
- [RES.COXETER.SYMMETRY.01](../../../../verification/status_ledger.csv#L465) · explicit_claim_identifier · original status at linked ledger row
- [experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean#L1) · explicit_file_citation
- [verification/v394_coxeter_atoms.py](../../../../verification/v394_coxeter_atoms.py#L1) · explicit_script_identifier
- [verification/v407_dn_pairings_omega.py](../../../../verification/v407_dn_pairings_omega.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)

## 2026-06-24 (The [formula in source] dual keystone + the flavor residuals it homes: v405– v407)

[documents/changelog/sections/013-source.tex:202-208](../../../../documents/changelog/sections/013-source.tex#L202) · `section:20ef03b5c6c62b042300` · subsection

/ / \ v405_seam_equiv_omega.py ( SEAM.EQUIV.02): the dual keystone at [formula in source] – the family/flavor sector as the order- [formula in source] Eisenstein [formula in source] face of [formula in source] , dual to the seam [formula in source] at [formula in source] ( v404). A named contract (genre of v286/ v398), NOT a closure. \ [formula in source] (the family CM point degenerates [formula in source] ) vs [for

Explicit source references:

- [FLAV.SELECTOR.CLOSE.01](../../../../verification/status_ledger.csv#L477) · explicit_claim_identifier · original status at linked ledger row
- [verification/v407_dn_pairings_omega.py](../../../../verification/v407_dn_pairings_omega.py#L1) · explicit_file_citation
- [FLAV.PMNS.CLOSE.01](../../../../verification/status_ledger.csv#L476) · explicit_claim_identifier · original status at linked ledger row
- [verification/v286_seam_equivalence_contract.py](../../../../verification/v286_seam_equivalence_contract.py#L1) · explicit_script_identifier
- [verification/v270_pmns_jarlskog_assembly.py](../../../../verification/v270_pmns_jarlskog_assembly.py#L1) · explicit_script_identifier
- [verification/v272_nu_mass_scale.py](../../../../verification/v272_nu_mass_scale.py#L1) · explicit_script_identifier
- [verification/v140_canonical_map.py](../../../../verification/v140_canonical_map.py#L1) · explicit_script_identifier
- [verification/v136_dual_normal_selector.py](../../../../verification/v136_dual_normal_selector.py#L1) · explicit_script_identifier
- [verification/v404_e8_char_cm_duality.py](../../../../verification/v404_e8_char_cm_duality.py#L1) · explicit_script_identifier
- [verification/v184_etaB_anchored_boltzmann.py](../../../../verification/v184_etaB_anchored_boltzmann.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.02](../../../../verification/status_ledger.csv#L475) · explicit_claim_identifier · original status at linked ledger row
- [verification/v9_neutrino_texture.py](../../../../verification/v9_neutrino_texture.py#L1) · explicit_script_identifier
- [verification/v406_pmns_phase_origin.py](../../../../verification/v406_pmns_phase_origin.py#L1) · explicit_file_citation
- [verification/v405_seam_equiv_omega.py](../../../../verification/v405_seam_equiv_omega.py#L1) · explicit_script_identifier
- [verification/v398_seam_state_rigidity.py](../../../../verification/v398_seam_state_rigidity.py#L1) · explicit_script_identifier
- [verification/v268_theta13_carrier_trace.py](../../../../verification/v268_theta13_carrier_trace.py#L1) · explicit_script_identifier
- [verification/v153_no_unit_theorem.py](../../../../verification/v153_no_unit_theorem.py#L1) · explicit_script_identifier
- [verification/v139_selector_triangle.py](../../../../verification/v139_selector_triangle.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v153) No-Unit Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1369)
- [(* ---- (v140) canonical map ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1123)
- [(* ---- (v139) selector triangle ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1110)
- [experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean#L1)
- [(* ==== v268: reactor-angle exponent = carrier hypercharge trace (FLAV.TH13.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2070)
- [(* ---- (v136) dual-normal selector ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1051)

## 2026-06-24 (The arithmetic hull + the [formula in source] -character CM duality: v403, v404)

[documents/changelog/sections/013-source.tex:209-214](../../../../documents/changelog/sections/013-source.tex#L209) · `section:8e3e6d08d3c102a09298` · subsection

/ \ v403_facet_compositum.py ( ARITH.HULL.01): the three single-prime number-field facets of the order- [formula in source] clock ( v390/ v394, used so far only one at a time) compose to one field [formula in source] . \ the squarefree parts [formula in source] are independent in [formula in source] , so [formula in source] with Galois group [formula in source] (order [formula in source] ); \ exactly [formula in sour

Explicit source references:

- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v404_e8_char_cm_duality.py](../../../../verification/v404_e8_char_cm_duality.py#L1) · explicit_script_identifier
- [ARITH.HULL.01](../../../../verification/status_ledger.csv#L473) · explicit_claim_identifier · original status at linked ledger row
- [verification/v267_qgeo_rigidity_minimal_axiom.py](../../../../verification/v267_qgeo_rigidity_minimal_axiom.py#L1) · explicit_script_identifier
- [verification/v390_prime2_facet.py](../../../../verification/v390_prime2_facet.py#L1) · explicit_script_identifier
- [verification/v403_facet_compositum.py](../../../../verification/v403_facet_compositum.py#L1) · explicit_script_identifier
- [verification/v220_cp_hexagonal_modulus.py](../../../../verification/v220_cp_hexagonal_modulus.py#L1) · explicit_script_identifier
- [E8.CMDUAL.01](../../../../verification/status_ledger.csv#L474) · explicit_claim_identifier · original status at linked ledger row
- [verification/v282_e8_tau_i_unification.py](../../../../verification/v282_e8_tau_i_unification.py#L1) · explicit_script_identifier
- [verification/v214_seam_pillowcase.py](../../../../verification/v214_seam_pillowcase.py#L1) · explicit_script_identifier
- [verification/v394_coxeter_atoms.py](../../../../verification/v394_coxeter_atoms.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v267: QGEO rigidity / minimal-axiom form of QGEO.SYM.01 (QGEO.SYM.02) ==== v262 (alpha_s RG), v263 (seesaw numpy), v264 (DtN FFT), v265 (RG+text guard), v266 (proton RG) are numerical -> Python-only; v267's exact cross-ratio/j core is mirrored here (the DtN/FFT part of v267 stays Python-only). *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2052)
- [(* ---- (v214) QGEO.PILLOW.01: pillowcase reduction -- cross-ratio 2 => j=1728 => order-4 CM ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1727)
- [(* ==== v282: E8-at-tau=i unification, chi_E8(i)=12 (QGAMB.UNIFY.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2171)

## 2026-06-24 (The five TOE/QFT closure contracts v398– v402 + a ledger rescope)

[documents/changelog/sections/013-source.tex:215-224](../../../../documents/changelog/sections/013-source.tex#L215) · `section:75d473d613aba16a08b8` · subsection

/ \ v398_seam_state_rigidity.py ( SEAM.RIGIDITY.01, Paper A): the Seam State Rigidity Theorem as ONE named target, consolidating the scattered state-invariance reductions ( v177/ v199/ v308/ v309/ v335) and hardening the finite-block evidence to a multi-mode sweep. \ RP-definability hinge (OS canonical transfer [formula in source] from RP [formula in source] reflection alone, v54, quasi-free v155, no [formula in sour

Explicit source references:

- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v386_grav_amplitude.py](../../../../verification/v386_grav_amplitude.py#L1) · explicit_script_identifier
- [verification/v274_scale_overdetermination.py](../../../../verification/v274_scale_overdetermination.py#L1) · explicit_script_identifier
- [METROLOGY.CLOSURE.01](../../../../verification/status_ledger.csv#L471) · explicit_claim_identifier · original status at linked ledger row
- [verification/v309_modular_bedrock.py](../../../../verification/v309_modular_bedrock.py#L1) · explicit_script_identifier
- [verification/v36_spectral_action_g2.py](../../../../verification/v36_spectral_action_g2.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v384_residual_certification.py](../../../../verification/v384_residual_certification.py#L1) · explicit_script_identifier
- [verification/v399_gravity_complete.py](../../../../verification/v399_gravity_complete.py#L1) · explicit_file_citation
- [verification/v155_quasifree_boundary.py](../../../../verification/v155_quasifree_boundary.py#L1) · explicit_script_identifier
- [GRAV.NONLINEAR.01](../../../../verification/status_ledger.csv#L75) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.RIGIDITY.01](../../../../verification/status_ledger.csv#L468) · explicit_claim_identifier · original status at linked ledger row
- [verification/v398_seam_state_rigidity.py](../../../../verification/v398_seam_state_rigidity.py#L1) · explicit_script_identifier
- [verification/v402_ftransfer_suite.py](../../../../verification/v402_ftransfer_suite.py#L1) · explicit_file_citation
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier
- [verification/v187_ftransfer_laws.py](../../../../verification/v187_ftransfer_laws.py#L1) · explicit_script_identifier
- [verification/v199_seam_state_invariance.py](../../../../verification/v199_seam_state_invariance.py#L1) · explicit_script_identifier
- [GRAV.NONLOCAL.01](../../../../verification/status_ledger.csv#L470) · explicit_claim_identifier · original status at linked ledger row
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_script_identifier
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v54_seam_horizon_keystones.py](../../../../verification/v54_seam_horizon_keystones.py#L1) · explicit_script_identifier
- [verification/v364_vgeo_sharpen.py](../../../../verification/v364_vgeo_sharpen.py#L1) · explicit_script_identifier
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v253_kappa_scalaron_ps.py](../../../../verification/v253_kappa_scalaron_ps.py#L1) · explicit_script_identifier
- [verification/v153_no_unit_theorem.py](../../../../verification/v153_no_unit_theorem.py#L1) · explicit_script_identifier
- [GATE.METRIC.01](../../../../verification/status_ledger.csv#L55) · explicit_claim_identifier · original status at linked ledger row
- [verification/v177_seam_marking_kernel.py](../../../../verification/v177_seam_marking_kernel.py#L1) · explicit_script_identifier
- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_script_identifier
- [FTRANSFER.SUITE.01](../../../../verification/status_ledger.csv#L472) · explicit_claim_identifier · original status at linked ledger row
- [verification/v401_metrology_closure.py](../../../../verification/v401_metrology_closure.py#L1) · explicit_file_citation
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [GRAVITY.COMPLETE.01](../../../../verification/status_ledger.csv#L469) · explicit_claim_identifier · original status at linked ledger row
- [verification/v308_seam_equiv_chain.py](../../../../verification/v308_seam_equiv_chain.py#L1) · explicit_script_identifier
- [verification/v400_grav_nonlocal_action.py](../../../../verification/v400_grav_nonlocal_action.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)
- [(* ---- (v199) QGEO.STATE.01: [rho,H]=0 <=> H is mu4-character-block-diagonal ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1649)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v36) G2: spectral-action Seeley-DeWitt coefficients -> R + R^2 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L113)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* ---- (v153) No-Unit Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1369)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CohomologyGrading.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CohomologyGrading.lean#L1)
- [(* ---- (v54) seam=horizon keystones: triply-forced 8 + shared transfer operator ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L235)

## 2026-06-24 (A falsifiable clock probe of the frontier: v397 how external is the ``external physics''?)

[documents/changelog/sections/013-source.tex:236-249](../../../../documents/changelog/sections/013-source.tex#L236) · `section:cb09371605fac20095ba` · subsection

v397 ( FR.CLOCK.PROBE.01) — do the ``external'' frontier scheme-numbers land on the [formula in source] clock? A falsifiable, anti-numerology probe (strict No-Free-Pattern, v100: ``on the clock'' means an exact match, not a wide-band fit). \ the [formula in source] decuple [formula in source] is ON the clock ( v212); \ the Koide rate [formula in source] is ON the clock ( v303); \ the axion spine angle [formula in sou

Explicit source references:

- [verification/v373_ftransfer_relic.py](../../../../verification/v373_ftransfer_relic.py#L1) · explicit_script_identifier
- [verification/v374_ftransfer_qcd.py](../../../../verification/v374_ftransfer_qcd.py#L1) · explicit_script_identifier
- [verification/v303_ftransfer_dynamics.py](../../../../verification/v303_ftransfer_dynamics.py#L1) · explicit_script_identifier
- [verification/v397_external_clock_probe.py](../../../../verification/v397_external_clock_probe.py#L1) · explicit_script_identifier
- [verification/v211_axion_spine_angle.py](../../../../verification/v211_axion_spine_angle.py#L1) · explicit_script_identifier
- [FR.CLOCK.PROBE.01](../../../../verification/status_ledger.csv#L467) · explicit_claim_identifier · original status at linked ledger row
- [verification/v100_numerology_null_mc.py](../../../../verification/v100_numerology_null_mc.py#L1) · explicit_script_identifier
- [verification/v187_ftransfer_laws.py](../../../../verification/v187_ftransfer_laws.py#L1) · explicit_script_identifier
- [verification/v212_leptogenesis_decuple.py](../../../../verification/v212_leptogenesis_decuple.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)

## 2026-06-24 (Attacking the last open lemma: v392 the scaling-limit central charge converges)

[documents/changelog/sections/013-source.tex:266-279](../../../../documents/changelog/sections/013-source.tex#L266) · `section:a553d00a39600c2d24cd` · subsection

v392 ( SEAM.S3.SCALINGLIMIT.01) — an honest attack on SEAM.EQUIV.01's one open lemma (continuum existence, v336). It does not close v336 (the abstract existence is the cited MMST theorem); it strengthens the case with the existence-relevant observable v376 did not test: that the finite-size data converge. \ the Calabrese–Cardy entanglement slope [formula in source] at [formula in source] is a monotone Cauchy-like seq

Explicit source references:

- [verification/v377_seam_s3_e8character.py](../../../../verification/v377_seam_s3_e8character.py#L1) · explicit_script_identifier
- [verification/v378_seam_s3_modular.py](../../../../verification/v378_seam_s3_modular.py#L1) · explicit_script_identifier
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [SEAM.S3.SCALINGLIMIT.01](../../../../verification/status_ledger.csv#L459) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v392_seam_s3_scalinglimit.py](../../../../verification/v392_seam_s3_scalinglimit.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)

## 2026-06-24 (Website: an animated ``Universal Gap Lab'' visualising v383– v391)

[documents/changelog/sections/013-source.tex:296-308](../../../../documents/changelog/sections/013-source.tex#L296) · `section:6b87e1b142e887eedca7` · subsection

Website — the UniversalGapLab animation on /verification. A new client component (Framer Motion, no new dependency) with three linked, module-grounded views: (A) every sector as the gap contraction [formula in source] racing to its attractor at [formula in source] (seam-gapped, v303/ v387/ v388) or the golden [formula in source] ( v312), with a gapless [formula in source] track that never forgets its start (the Lean

Explicit source references:

- [verification/v386_grav_amplitude.py](../../../../verification/v386_grav_amplitude.py#L1) · explicit_script_identifier
- [verification/v387_corrections_gap.py](../../../../verification/v387_corrections_gap.py#L1) · explicit_script_identifier
- [verification/v391_alpha_quillen_progress.py](../../../../verification/v391_alpha_quillen_progress.py#L1) · explicit_script_identifier
- [FORM.SPECTRALGAP.01](../../../../verification/status_ledger.csv#L454) · explicit_claim_identifier · original status at linked ledger row
- [verification/v384_residual_certification.py](../../../../verification/v384_residual_certification.py#L1) · explicit_script_identifier
- [verification/v312_hypergraph_substrate.py](../../../../verification/v312_hypergraph_substrate.py#L1) · explicit_script_identifier
- [verification/v303_ftransfer_dynamics.py](../../../../verification/v303_ftransfer_dynamics.py#L1) · explicit_script_identifier
- [verification/v388_corrections_budget.py](../../../../verification/v388_corrections_budget.py#L1) · explicit_script_identifier
- [verification/v383_dynamics_universal.py](../../../../verification/v383_dynamics_universal.py#L1) · explicit_script_identifier
- [verification/v389_grav_loop_finiteness.py](../../../../verification/v389_grav_loop_finiteness.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)

## 2026-06-24 (Completing the Coxeter clock: v390 the prime-2 (Gaussian / [formula in source] -seam) facet of [formula in source] )

[documents/changelog/sections/013-source.tex:341-356](../../../../documents/changelog/sections/013-source.tex#L341) · `section:f1054498b80cd7df6676` · subsection

v390 ( COXETER.PRIME2.01) — the prime-2 facet of the order-30 Coxeter clock. v383/ v387 named two number-field facets of the clock [formula in source] — the golden [formula in source] (prime-5, carrier) and [formula in source] (prime-3, family). This names the prime-2 facet and identifies it with structure already in the theory: the Gaussian field [formula in source] , i.e.\ the CM point [formula in source] and the o

Explicit source references:

- [verification/v387_corrections_gap.py](../../../../verification/v387_corrections_gap.py#L1) · explicit_script_identifier
- [verification/v390_prime2_facet.py](../../../../verification/v390_prime2_facet.py#L1) · explicit_script_identifier
- [verification/v288_full_l2_subprincipal_z4.py](../../../../verification/v288_full_l2_subprincipal_z4.py#L1) · explicit_script_identifier
- [verification/v383_dynamics_universal.py](../../../../verification/v383_dynamics_universal.py#L1) · explicit_script_identifier
- [verification/v333_cm_selection.py](../../../../verification/v333_cm_selection.py#L1) · explicit_script_identifier
- [verification/v282_e8_tau_i_unification.py](../../../../verification/v282_e8_tau_i_unification.py#L1) · explicit_script_identifier
- [COXETER.PRIME2.01](../../../../verification/status_ledger.csv#L453) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v282: E8-at-tau=i unification, chi_E8(i)=12 (QGAMB.UNIFY.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2171)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)

## 2026-06-23 (Track 2 of the forward plan v2 — the AMBIENT REDUNDANCY statement ( v369): the holographic route to Gravity-complete that sidesteps building the general Euclidean QG measure)

[documents/changelog/sections/014-source.tex:107-125](../../../../documents/changelog/sections/014-source.tex#L107) · `section:acb731fecb30c7953915` · subsection

v369 ( QGAMB.REDUNDANCY.01) — Track 2 (forward plan v2). Instead of constructing the non-perturbative ambient measure (the conformal-factor swamp, v332/ v334), assemble the established discriminators that prove the ambient sector carries no physical degrees of freedom relative to the holomorphic boundary algebra, so the QG measure (gate C7/ QG.AMB.01) is a non-fundamental certification object (like a coordinate choic

Explicit source references:

- [verification/v329_os_gap_reduction.py](../../../../verification/v329_os_gap_reduction.py#L1) · explicit_script_identifier
- [verification/v334_conformal_resolution.py](../../../../verification/v334_conformal_resolution.py#L1) · explicit_script_identifier
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v337_decoupling_theorem.py](../../../../verification/v337_decoupling_theorem.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_script_identifier
- [verification/v332_qgamb_metric_sector.py](../../../../verification/v332_qgamb_metric_sector.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [verification/v366_mmst_seam_collar.py](../../../../verification/v366_mmst_seam_collar.py#L1) · explicit_script_identifier
- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier
- [QGAMB.REDUNDANCY.01](../../../../verification/status_ledger.csv#L94) · explicit_claim_identifier · original status at linked ledger row
- [verification/v2_carrier_pascal.py](../../../../verification/v2_carrier_pascal.py#L1) · explicit_script_identifier
- [verification/v221_seam_qecc.py](../../../../verification/v221_seam_qecc.py#L1) · explicit_script_identifier
- [verification/v309_modular_bedrock.py](../../../../verification/v309_modular_bedrock.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v2) carrier / Pascal ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L51)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/ArfSpinorCompiler.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/ArfSpinorCompiler.lean#L1)
- [(* ==== v337 DECOUPLING.THEOREM.01: the ambient back-reaction scale identity ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2241)

## 2026-06-23 (Status propagation: the v365– v380 [formula in source] FORM.SEAM.MMST.01 status carried into the remaining papers, the README and the Zenodo description)

[documents/changelog/sections/014-source.tex:126-139](../../../../documents/changelog/sections/014-source.tex#L126) · `section:43e53494baf54171d70d` · subsection

Deep-sync of the post-S3 status. The canonical surfaces ( tfpt_status, introduction) were already current; this pass propagated the same honest wording into origin_theory, tfpt_1– tfpt_5, tfpt_research_contracts, the root README.md, the Zenodo description and website/lib/papers.ts: SEAM.EQUIV.01 is closed modulo a cited theorem (lattice v367/ v368 [formula in source] the [formula in source] stack v376– v379 [formula

Explicit source references:

- [verification/v371_ftransfer_pole.py](../../../../verification/v371_ftransfer_pole.py#L1) · explicit_script_identifier
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v375_observatory_registry.py](../../../../verification/v375_observatory_registry.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [verification/v365_qg_oneloop_saddle.py](../../../../verification/v365_qg_oneloop_saddle.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [README.md](../../../../README.md#L1) · explicit_file_citation
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [verification/v374_ftransfer_qcd.py](../../../../verification/v374_ftransfer_qcd.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)

## 2026-06-23 ( v380 — the KMS Entire Hessian: the Stelle ghost is EXACTLY the Seeley–DeWitt truncation, and resummation pushes it to infinity)

[documents/changelog/sections/014-source.tex:140-153](../../../../documents/changelog/sections/014-source.tex#L140) · `section:6db706756a796987fc17` · subsection

v380 ( GRAV.KMS.HESSIAN.01). Upgrades v304/ v370 from the assumption ``the resummed graviton form factor is entire'' to a derived statement. \ the seam KMS cutoff [formula in source] ( v259) gives the dressed spin-2 propagator [formula in source] with [formula in source] entire and nowhere zero (only the [formula in source] pole). \ each finite Seeley–DeWitt truncation is a Taylor partial sum of [formula in source] a

Explicit source references:

- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_script_identifier
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier
- [GRAV.KMS.HESSIAN.01](../../../../verification/status_ledger.csv#L93) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)

## 2026-06-23 (Lean: the S3 continuum leg formalised as ``closed modulo a cited theorem'' — FORM.SEAM.MMST.01)

[documents/changelog/sections/014-source.tex:154-167](../../../../documents/changelog/sections/014-source.tex#L154) · `section:ded29520f6b18df4c062` · subsection

FORM.SEAM.MMST.01 — Lean 4 formalisation. The S3 continuum leg is now machine-formalised in a new module ( SeamScalingLimit.lean), lake build OK and \#print axioms clean, in the audit-contract style of SeamEquivChain.lean. The MMST applicability hypotheses for the seam collar are provable by kernel decide (no axioms): [formula in source] , [formula in source] , the range [formula in source] , [formula in source] , an

Explicit source references:

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1) · explicit_file_citation
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1) · explicit_file_citation
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)

## 2026-06-23 (S3 closure stack — the target of the seam scaling limit pinned at EVERY computable level: central charge [formula in source] from the lattice ( v376), the [formula in source] character ( v377), the genus-1 torus count [formula in source] ( v378), and reflection positivity ( v379); only the abstract continuum existence (cited MMST) stays external)

[documents/changelog/sections/014-source.tex:168-192](../../../../documents/changelog/sections/014-source.tex#L168) · `section:b450ecffaba882770fc8` · subsection

v376 ( SEAM.S3.CENTRALCHARGE.01). The Calabrese–Cardy finite-size entanglement-entropy scaling of the critical free-fermion content gives [formula in source] per complex mode (Peschel correlation-matrix EE; the [formula in source] sizes avoid the half-filling zero-mode), so the 16-Majorana ( [formula in source] complex) collar has [formula in source] , chiral ( [formula in source] , v367/ v368). Numerical [formula in

Explicit source references:

- [SEAM.S3.E8CHARACTER.01](../../../../verification/status_ledger.csv#L89) · explicit_claim_identifier · original status at linked ledger row
- [verification/v377_seam_s3_e8character.py](../../../../verification/v377_seam_s3_e8character.py#L1) · explicit_script_identifier
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [SEAM.S3.MODULAR.01](../../../../verification/status_ledger.csv#L90) · explicit_claim_identifier · original status at linked ledger row
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [SEAM.S3.CENTRALCHARGE.01](../../../../verification/status_ledger.csv#L88) · explicit_claim_identifier · original status at linked ledger row
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v378_seam_s3_modular.py](../../../../verification/v378_seam_s3_modular.py#L1) · explicit_script_identifier
- [SEAM.S3.RP.01](../../../../verification/status_ledger.csv#L91) · explicit_claim_identifier · original status at linked ledger row
- [verification/v344_detk_synthesis.py](../../../../verification/v344_detk_synthesis.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v344 SEAM.DETK.01: the ADE |det Cartan| = |H_1(link)| sequence, only E8 -> 1 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2272)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)

## 2026-06-23 (Track 1 — the S3 input made CONCRETE: an explicit gapped chiral-Majorana ( [formula in source] ) lattice model with numerical Chern [formula in source] realises the [formula in source] edge ( v367), and a strip shows the chiral edge exists by anomaly inflow ( v368))

[documents/changelog/sections/014-source.tex:224-242](../../../../documents/changelog/sections/014-source.tex#L224) · `section:1a6d16256ec9d6226d97` · subsection

v367 ( SEAM.S3.LATTICE.01) — Track 1. The abstract S3 residual of v356/ v366 (``the collar realised as a genuine lattice chiral free-fermion invertible phase'') is turned into a RUNNABLE model. \ the [formula in source] Bloch model is gapped at [formula in source] with a numerically computed Fukui–Hatsugai–Suzuki Chern number [formula in source] ( [formula in source] trivially at [formula in source] ); \ [formula in

Explicit source references:

- [SEAM.S3.LATTICE.01](../../../../verification/status_ledger.csv#L81) · explicit_claim_identifier · original status at linked ledger row
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v356_continuum_mmst_applicability.py](../../../../verification/v356_continuum_mmst_applicability.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v366_mmst_seam_collar.py](../../../../verification/v366_mmst_seam_collar.py#L1) · explicit_script_identifier
- [verification/v154_simple_current_theorem.py](../../../../verification/v154_simple_current_theorem.py#L1) · explicit_script_identifier
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v351_continuum_realisation_sharpened.py](../../../../verification/v351_continuum_realisation_sharpened.py#L1) · explicit_script_identifier
- [SEAM.S3.INFLOW.01](../../../../verification/status_ledger.csv#L82) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v154) Simple-Current Extension Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1381)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v351 SEAM.EQUIV.CONTINUUM.02: the order-4 clock + Coxeter-match pin E8 over SO(16) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2335)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)

## 2026-06-23 (Directions 3 \& 7 — the two hard functional-analysis residuals ADVANCED: the QG measure as one-loop fluctuations around the parameter-free saddle ( v365), and the seam collar verified IN MMST's free-fermion class so the scaling limit is [formula in source] ( v366))

[documents/changelog/sections/014-source.tex:258-280](../../../../documents/changelog/sections/014-source.tex#L258) · `section:4a157afd200102fcaaee` · subsection

v365 ( QGAMB.SADDLE.01) — Direction 3. With the parameter-free saddle ( v358/ v359) the QG measure (gate C7/ QG.AMB.01) organises as a one-loop Gaussian fluctuation determinant. \ the fluctuation stiffness is [formula in source] (no free dial, so the one-loop problem is now well-posed); \ the fluctuation operator is gapped [formula in source] ( v76/ v337) [formula in source] mode-by-mode convergence, finite model log

Explicit source references:

- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_script_identifier
- [verification/v351_continuum_realisation_sharpened.py](../../../../verification/v351_continuum_realisation_sharpened.py#L1) · explicit_script_identifier
- [verification/v334_conformal_resolution.py](../../../../verification/v334_conformal_resolution.py#L1) · explicit_script_identifier
- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v160_seam_gaussianity_from_pf.py](../../../../verification/v160_seam_gaussianity_from_pf.py#L1) · explicit_script_identifier
- [verification/v356_continuum_mmst_applicability.py](../../../../verification/v356_continuum_mmst_applicability.py#L1) · explicit_script_identifier
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [verification/v297_route_a_literature_stack.py](../../../../verification/v297_route_a_literature_stack.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v155_quasifree_boundary.py](../../../../verification/v155_quasifree_boundary.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v337_decoupling_theorem.py](../../../../verification/v337_decoupling_theorem.py#L1) · explicit_script_identifier
- [verification/v302_seam_gap.py](../../../../verification/v302_seam_gap.py#L1) · explicit_script_identifier
- [verification/v330_qgamb_admissible_measure.py](../../../../verification/v330_qgamb_admissible_measure.py#L1) · explicit_script_identifier
- [QGAMB.SADDLE.01](../../../../verification/status_ledger.csv#L78) · explicit_claim_identifier · original status at linked ledger row
- [verification/v364_vgeo_sharpen.py](../../../../verification/v364_vgeo_sharpen.py#L1) · explicit_script_identifier
- [SEAM.MMST.INCLASS.01](../../../../verification/status_ledger.csv#L79) · explicit_claim_identifier · original status at linked ledger row
- [verification/v366_mmst_seam_collar.py](../../../../verification/v366_mmst_seam_collar.py#L1) · explicit_script_identifier
- [verification/v365_qg_oneloop_saddle.py](../../../../verification/v365_qg_oneloop_saddle.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v337 DECOUPLING.THEOREM.01: the ambient back-reaction scale identity ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2241)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [(* ==== v351 SEAM.EQUIV.CONTINUUM.02: the order-4 clock + Coxeter-match pin E8 over SO(16) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2335)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)

## 2026-06-23 (Direction 2, the disciplined result: the gap-induced GR correction is the known [formula in source] scalaron — already falsifiable via [formula in source] ; a new forced near-term deviation does NOT survive the discriminator and is declined: v360)

[documents/changelog/sections/014-source.tex:331-347](../../../../documents/changelog/sections/014-source.tex#L331) · `section:b9e3cb43b2bc303303e6` · subsection

v360 ( GRAV.GAPCORR.01) — Direction 2, run honestly. Does the seam transfer gap [formula in source] yield a NEW, calculable, FORCED deviation from GR (a near-term kill test) beyond the known [formula in source] scalaron? With the v354/ v355 discriminator the answer is no, and this is recorded rather than manufacturing a coefficient. \ the leading higher-curvature term IS the [formula in source] scalaron (spectral-act

Explicit source references:

- [verification/v354_e8_reverse_audit.py](../../../../verification/v354_e8_reverse_audit.py#L1) · explicit_script_identifier
- [verification/v355_e8_unmapped_bandwidth.py](../../../../verification/v355_e8_unmapped_bandwidth.py#L1) · explicit_script_identifier
- [verification/v36_spectral_action_g2.py](../../../../verification/v36_spectral_action_g2.py#L1) · explicit_script_identifier
- [verification/v360_grav_gap_corrections.py](../../../../verification/v360_grav_gap_corrections.py#L1) · explicit_script_identifier
- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_script_identifier
- [GRAV.GAPCORR.01](../../../../verification/status_ledger.csv#L74) · explicit_claim_identifier · original status at linked ledger row
- [verification/v28_gravity_fR.py](../../../../verification/v28_gravity_fR.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v36) G2: spectral-action Seeley-DeWitt coefficients -> R + R^2 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L113)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)
- [(* ==== v354 E8.REVERSE.AUDIT.01: 3 of 8 E8 Casimir degrees are load-bearing, 5 unmapped ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2353)
- [(* ==== v355 E8.UNMAPPED.BANDWIDTH.01: the unmapped region's forced content is COLLECTIVE ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2363)

## 2026-06-23 (the two frontier walls attacked directly: Direction A assembles the continuum chain of SEAM.EQUIV.01 theorem-by-theorem and ties its residual to the seam one-sidedness; Direction B was found ALREADY closed by v40 (harmonic metric [formula in source] finite linear algebra, not a PDE): v356)

[documents/changelog/sections/015-source.tex:17-40](../../../../documents/changelog/sections/015-source.tex#L17) · `section:5f8b87dd447560390c4e` · subsection

v356 ( SEAM.EQUIV.CONTINUUM.03) — Direction A, the continuum chain assembled. After the reverse-audit/bandwidth work showed no hidden discrete lever remains, the only fundamental progress is to assemble the continuum chain. Link-by-link: S1 the gapped ( [formula in source] , v302) quasi-free [formula in source] -Majorana collar ( v155/ v160) [formula in source] S2 invertible (gapped free fermions have no intrinsic to

Explicit source references:

- [verification/v31_R_dictionary.py](../../../../verification/v31_R_dictionary.py#L1) · explicit_script_identifier
- [verification/v302_seam_gap.py](../../../../verification/v302_seam_gap.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.CONTINUUM.03](../../../../verification/status_ledger.csv#L77) · explicit_claim_identifier · original status at linked ledger row
- [verification/v356_continuum_mmst_applicability.py](../../../../verification/v356_continuum_mmst_applicability.py#L1) · explicit_script_identifier
- [verification/v155_quasifree_boundary.py](../../../../verification/v155_quasifree_boundary.py#L1) · explicit_script_identifier
- [verification/v160_seam_gaussianity_from_pf.py](../../../../verification/v160_seam_gaussianity_from_pf.py#L1) · explicit_script_identifier
- [verification/v297_route_a_literature_stack.py](../../../../verification/v297_route_a_literature_stack.py#L1) · explicit_script_identifier
- [verification/v351_continuum_realisation_sharpened.py](../../../../verification/v351_continuum_realisation_sharpened.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v301_route_a_invertible.py](../../../../verification/v301_route_a_invertible.py#L1) · explicit_script_identifier
- [verification/v58_seam_horizon_chain.py](../../../../verification/v58_seam_horizon_chain.py#L1) · explicit_script_identifier
- [verification/v40_harmonic_metric.py](../../../../verification/v40_harmonic_metric.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v351 SEAM.EQUIV.CONTINUUM.02: the order-4 clock + Coxeter-match pin E8 over SO(16) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2335)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ---- (v58) seam-horizon chain: one-sided Gauss-Bonnet, KMS, S_BH, RT in seam units ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L270)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)

## 2026-06-22 (two consequences: the [formula in source] which-net ambiguity resolved by the order-4 clock, and BOTH axioms shown bootstrap-forced so the only irreducibles are the framework [formula in source] [formula in source] : v351, v352)

[documents/changelog/sections/015-source.tex:100-119](../../../../documents/changelog/sections/015-source.tex#L100) · `section:a5034e3be2279bec9f75` · subsection

v351 ( SEAM.EQUIV.CONTINUUM.02) — the which-net ambiguity falls. The long-standing `` [formula in source] ambiguity'' ( [formula in source] det 1 vs [formula in source] det 4, both [formula in source] , flagged open in v277/ v344) is RESOLVED, non-circularly, by the seam's ORDER-4 [formula in source] clock: the index-4 simple-current extension of [formula in source] is [formula in source] ( v154), the index-2 would b

Explicit source references:

- [verification/v154_simple_current_theorem.py](../../../../verification/v154_simple_current_theorem.py#L1) · explicit_script_identifier
- [verification/v6_bootstrap.py](../../../../verification/v6_bootstrap.py#L1) · explicit_script_identifier
- [verification/v344_detk_synthesis.py](../../../../verification/v344_detk_synthesis.py#L1) · explicit_script_identifier
- [verification/v54_seam_horizon_keystones.py](../../../../verification/v54_seam_horizon_keystones.py#L1) · explicit_script_identifier
- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [verification/v352_framework_irreducible.py](../../../../verification/v352_framework_irreducible.py#L1) · explicit_script_identifier
- [verification/v351_continuum_realisation_sharpened.py](../../../../verification/v351_continuum_realisation_sharpened.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.CONTINUUM.02](../../../../verification/status_ledger.csv#L71) · explicit_claim_identifier · original status at linked ledger row
- [verification/v83_e8net_holomorphic_uniqueness.py](../../../../verification/v83_e8net_holomorphic_uniqueness.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [TFPT.IRREDUCIBLE.01](../../../../verification/status_ledger.csv#L103) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v344 SEAM.DETK.01: the ADE |det Cartan| = |H_1(link)| sequence, only E8 -> 1 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2272)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ---- (v154) Simple-Current Extension Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1381)
- [(* ==== v352 TFPT.IRREDUCIBLE.01: both axioms bootstrap-forced; irreducibles = framework + pi ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2345)
- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean#L1)
- [(* ---- (v54) seam=horizon keystones: triply-forced 8 + shared transfer operator ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L235)
- [(* ==== v351 SEAM.EQUIV.CONTINUUM.02: the order-4 clock + Coxeter-match pin E8 over SO(16) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2335)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)

## 2026-06-22 (correction: the inputs are bootstrap-forced fixed points, not free axioms – [formula in source] IS the icosahedral [formula in source] of [formula in source] , and the golden ratio is emergent from the [formula in source] -glue, not external: v350)

[documents/changelog/sections/015-source.tex:120-134](../../../../documents/changelog/sections/015-source.tex#L120) · `section:2ef143ddeb6e7890a438` · subsection

v350 ( SEAM.EQUIV.BOOTSTRAP.01) — correcting v349. An honest correction prompted by the right objection that the inputs are not axioms but are back-determined by the theory. \ the inputs ARE bootstrap-forced ( v6): [formula in source] is over-determined three ways – rank-fill ( [formula in source] ), Coxeter-match ( [formula in source] max prime of [formula in source] ), Pascal ( [formula in source] ). \ so [formula

Explicit source references:

- [verification/v349_raw_seam_golden_test.py](../../../../verification/v349_raw_seam_golden_test.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v350_bootstrap_inputs_correction.py](../../../../verification/v350_bootstrap_inputs_correction.py#L1) · explicit_script_identifier
- [verification/v154_simple_current_theorem.py](../../../../verification/v154_simple_current_theorem.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.BOOTSTRAP.01](../../../../verification/status_ledger.csv#L104) · explicit_claim_identifier · original status at linked ledger row
- [verification/v6_bootstrap.py](../../../../verification/v6_bootstrap.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ---- (v154) Simple-Current Extension Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1381)
- [(* ==== v350 SEAM.EQUIV.BOOTSTRAP.01: g_car=5 IS the icosahedral 5 of h(E8)=30 (the v349 correction) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2324)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v349 SEAM.EQUIV.GOLDEN.01: the raw carrier is NOT golden (D5 h=8 -> sqrt2); 5|h only on the output side ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2315)

## 2026-06-22 (the honest test – does the raw seam carry [formula in source] ? – answer NO: the bare carrier gives [formula in source] , golden is the icosahedral input; the keystone reduces to ``is [formula in source] a pentagon or a count?'': v349)

[documents/changelog/sections/015-source.tex:135-151](../../../../documents/changelog/sections/015-source.tex#L135) · `section:fac429be371842092df9` · subsection

v349 ( SEAM.EQUIV.GOLDEN.01). The honest, computed test of the v348 question, and a clarifying NEGATIVE: the raw seam does not carry the golden ratio. \ the golden ratio [formula in source] is the signature of a [formula in source] -fold rotation ( [formula in source] ); the carrier [formula in source] has [formula in source] (eigenvalue [formula in source] ) and [formula in source] has [formula in source] – no 5-fol

Explicit source references:

- [verification/v349_raw_seam_golden_test.py](../../../../verification/v349_raw_seam_golden_test.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.GOLDEN.01](../../../../verification/status_ledger.csv#L105) · explicit_claim_identifier · original status at linked ledger row
- [verification/v348_seam_rigidity_route.py](../../../../verification/v348_seam_rigidity_route.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v349 SEAM.EQUIV.GOLDEN.01: the raw carrier is NOT golden (D5 h=8 -> sqrt2); 5|h only on the output side ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2315)
- [(* ==== v348 SEAM.EQUIV.RIGID.01: phi bridges crystallographic mu4 to non-crystallographic 2I ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2307)

## 2026-06-22 (Route B carried out: rigidity [formula in source] the icosian-ring identity reduce the whole keystone to ONE question – does the raw seam carry the golden ratio? – the simplest form of L2: v348)

[documents/changelog/sections/015-source.tex:152-167](../../../../documents/changelog/sections/015-source.tex#L152) · `section:414a84d23afa6566ca5e` · subsection

v348 ( SEAM.EQUIV.RIGID.01). Executes the rigidity route (mode B of v347) and lands on a strikingly simple statement. \ rigidity closes the crystallographic part: the order-4 clock has the unique normal form [formula in source] (Nielsen/Kerckhoff, v180), Troyanov gives the unique flat pillowcase metric ( v284). \ the downstream is a LATTICE identity, not just a graph: by Conway–Sloane the ring of icosians ( [formula

Explicit source references:

- [verification/v284_route_i_rp_uniqueness.py](../../../../verification/v284_route_i_rp_uniqueness.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.RIGID.01](../../../../verification/status_ledger.csv#L106) · explicit_claim_identifier · original status at linked ledger row
- [verification/v347_seam_closure_modes.py](../../../../verification/v347_seam_closure_modes.py#L1) · explicit_script_identifier
- [verification/v312_hypergraph_substrate.py](../../../../verification/v312_hypergraph_substrate.py#L1) · explicit_script_identifier
- [verification/v348_seam_rigidity_route.py](../../../../verification/v348_seam_rigidity_route.py#L1) · explicit_script_identifier
- [verification/v180_clock_is_mobius.py](../../../../verification/v180_clock_is_mobius.py#L1) · explicit_script_identifier
- [verification/v313_golden_atoms_spectral.py](../../../../verification/v313_golden_atoms_spectral.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v313-v320: the cyclotomic / Galois arithmetic arc (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2183)
- [(* ==== v348 SEAM.EQUIV.RIGID.01: phi bridges crystallographic mu4 to non-crystallographic 2I ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2307)
- [(* ==== v347 SEAM.EQUIV.CLOSURE.01: the flat seam base vs the spherical icosahedral base ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2298)
- [experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean#L1)
- [(* v313 GOLD.ATOMS.01 -- the golden ratio is the g_car=5 spectral signature *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2185)

## 2026-06-22 (can the one open arrow be fully solved? the closure-mode classification: the flat [formula in source] spherical base locus, four modes, no free theorem – only rigidity or axiom [formula in source] : v347)

[documents/changelog/sections/015-source.tex:168-182](../../../../documents/changelog/sections/015-source.tex#L168) · `section:9800e2e20cb3dc2e280a` · subsection

v347 ( SEAM.EQUIV.CLOSURE.01). A direct, honest answer to ``is there any other way to fully solve the one open arrow L2?''. \ the precise locus: the seam gives the flat pillowcase base [formula in source] ( [formula in source] , v216) while the [formula in source] structure is the Seifert base over the spherical [formula in source] ( [formula in source] ), so L2 bridges two geometric types – why it needs the carrier

Explicit source references:

- [verification/v347_seam_closure_modes.py](../../../../verification/v347_seam_closure_modes.py#L1) · explicit_script_identifier
- [verification/v331_necessity_of_H.py](../../../../verification/v331_necessity_of_H.py#L1) · explicit_script_identifier
- [verification/v246_unification_data_crosscheck.py](../../../../verification/v246_unification_data_crosscheck.py#L1) · explicit_script_identifier
- [verification/v300_flataway_rigid.py](../../../../verification/v300_flataway_rigid.py#L1) · explicit_script_identifier
- [verification/v284_route_i_rp_uniqueness.py](../../../../verification/v284_route_i_rp_uniqueness.py#L1) · explicit_script_identifier
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.CLOSURE.01](../../../../verification/status_ledger.csv#L107) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [(* ==== v347 SEAM.EQUIV.CLOSURE.01: the flat seam base vs the spherical icosahedral base ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2298)

## 2026-06-22 (the geometric bridge made explicit: the keystone is ONE arrow – raw-seam [formula in source] is a six-link chain, five theorems [formula in source] one open orbifold realisation, the group forced by the [formula in source] atoms: v346)

[documents/changelog/sections/015-source.tex:183-197](../../../../documents/changelog/sections/015-source.tex#L183) · `section:00c96659d1c56ecb746e` · subsection

v346 ( SEAM.EQUIV.GEOM.01). The honest capstone of the SEAM.EQUIV.01 investigation: it lays out the full path raw-seam [formula in source] as six links and shows exactly one is open. / \ L1 raw seam [formula in source] [formula in source] four marks (Gauss–Bonnet, v216); \ L2 [formula in source] [formula in source] the [formula in source] atoms [formula in source] the [formula in source] orbifold [formula in source]

Explicit source references:

- [verification/v219_icosahedral_mckay.py](../../../../verification/v219_icosahedral_mckay.py#L1) · explicit_script_identifier
- [verification/v346_seam_geometric_bridge.py](../../../../verification/v346_seam_geometric_bridge.py#L1) · explicit_script_identifier
- [verification/v345_hypergraph_homotopy.py](../../../../verification/v345_hypergraph_homotopy.py#L1) · explicit_script_identifier
- [verification/v232_e8_kleinian_seam.py](../../../../verification/v232_e8_kleinian_seam.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.GEOM.01](../../../../verification/status_ledger.csv#L108) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v345 SEAM.DETK.02: the plumbing-link H_1 = coker(Cartan) via Smith normal form ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2285)
- [(* ==== v219-v230 round: icosahedral McKay, CM-norm duality, structural finds ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1781)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [(* v219 McKay: affine E8 marks = 2I irrep degrees, sums 30 and 120 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1783)
- [(* ==== v232: the seam as the E8 Kleinian singularity (du Val resolution + link) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1895)

## 2026-06-22 (executing route R3: the [formula in source] -rewrite attractor link computed to be the Poincar\'e homology sphere ( [formula in source] ), the combinatorial core complete, only the geometric bridge open: v345)

[documents/changelog/sections/015-source.tex:198-213](../../../../documents/changelog/sections/015-source.tex#L198) · `section:f45ff82a01e81b33f3c1` · subsection

v345 ( SEAM.DETK.02). Carries out the combinatorial core of the hypergraph-homotopy route (R3, the most promising fresh angle from v344) as far as it computably goes. By plumbing topology (Brieskorn/Hirzebruch) the link of a tree plumbing with the ADE Cartan intersection form has [formula in source] . \ the Smith normal form of the [formula in source] Cartan matrix is the identity, so [formula in source] – the [formu

Explicit source references:

- [SEAM.DETK.02](../../../../verification/status_ledger.csv#L109) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v344_detk_synthesis.py](../../../../verification/v344_detk_synthesis.py#L1) · explicit_script_identifier
- [verification/v345_hypergraph_homotopy.py](../../../../verification/v345_hypergraph_homotopy.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v344 SEAM.DETK.01: the ADE |det Cartan| = |H_1(link)| sequence, only E8 -> 1 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2272)
- [(* ==== v345 SEAM.DETK.02: the plumbing-link H_1 = coker(Cartan) via Smith normal form ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2285)

## 2026-06-22 (the one keystone bit [formula in source] mapped: six equivalent faces, the genus-1 obstruction, three access routes, the hypergraph-homotopy route the most promising fresh angle: v344)

[documents/changelog/sections/015-source.tex:214-229](../../../../documents/changelog/sections/015-source.tex#L214) · `section:ee9a2e2d9ccdc3799662` · subsection

v344 ( SEAM.DETK.01). A bird's-eye, full-spectrum map of SEAM.EQUIV.01's only open residual, the holomorphy bit [formula in source] . \ it is ONE statement with six equivalent faces (holomorphy / homology-sphere link [formula in source] / one [formula in source] -dim irrep / [formula in source] Cartan [formula in source] / full [formula in source] condensation / [formula in source] CM), unified by [formula in source]

Explicit source references:

- [verification/v87_bulk_uniqueness_reduction.py](../../../../verification/v87_bulk_uniqueness_reduction.py#L1) · explicit_script_identifier
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_script_identifier
- [SEAM.DETK.01](../../../../verification/status_ledger.csv#L110) · explicit_claim_identifier · original status at linked ledger row
- [verification/v312_hypergraph_substrate.py](../../../../verification/v312_hypergraph_substrate.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v232_e8_kleinian_seam.py](../../../../verification/v232_e8_kleinian_seam.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v344_detk_synthesis.py](../../../../verification/v344_detk_synthesis.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ==== v344 SEAM.DETK.01: the ADE |det Cartan| = |H_1(link)| sequence, only E8 -> 1 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2272)
- [(* ---- (v87) SO(16)_1 modular data and its six invariants ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L98)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v232: the seam as the E8 Kleinian singularity (du Val resolution + link) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1895)

## 2026-06-22 (the four black-hole-birth completion routes investigated: none closes all or is 100\% provable alone; all four converge on SEAM.EQUIV.01 [formula in source] the decoupled QG contour: v343)

[documents/changelog/sections/015-source.tex:230-245](../../../../documents/changelog/sections/015-source.tex#L230) · `section:4dcf513f2ae74de8645b` · subsection

v343 ( FOUR.ROUTES.01). A direct, honest answer to ``which black-hole-birth route closes ALL problems and is 100\% provable'': none individually — but the four converge on the same two residuals as v335. \ (A) the finite causal diamond reframes QG.AMB to a finite thermal trace (de Sitter [formula in source] finite, [formula in source] ) but inherits the GHP contour ( v334, ); does not touch the seam. \ (B) Carlip–Car

Explicit source references:

- [FOUR.ROUTES.01](../../../../verification/status_ledger.csv#L111) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v334_conformal_resolution.py](../../../../verification/v334_conformal_resolution.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v239_kms_thermal_time.py](../../../../verification/v239_kms_thermal_time.py#L1) · explicit_script_identifier
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_script_identifier
- [verification/v343_four_routes_analysis.py](../../../../verification/v343_four_routes_analysis.py#L1) · explicit_script_identifier
- [verification/v56_unique_attractor.py](../../../../verification/v56_unique_attractor.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v56) unique attractor: gapped transport, Coxeter planes ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L257)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)

## 2026-06-22 (the heat-kernel origin of the [formula in source] terms: the Calder\'on quadratic is the Gilkey gauge-curvature coefficient; the cubic Maxwell moment stays the EM.WARD.01 residual: v342)

[documents/changelog/sections/015-source.tex:246-259](../../../../documents/changelog/sections/015-source.tex#L246) · `section:ed458665aa45a0bcd44b` · subsection

v342 ( EM.WARD.02) — the heat-kernel sharpening of EM.WARD.01. Derives the ORIGIN and STRUCTURE of the [formula in source] determinant-line terms from textbook Seeley–DeWitt / Gilkey coefficients, advancing the EM-Ward origin from \ (conjectured) to \ for the term structure; no new gate. \ [formula in source] is the boundary Gauss–Bonnet coefficient ( [formula in source] , [formula in source] , v216). \ the [formula

Explicit source references:

- [EM.WARD.01](../../../../verification/status_ledger.csv#L70) · explicit_claim_identifier · original status at linked ledger row
- [verification/v3_em_alpha.py](../../../../verification/v3_em_alpha.py#L1) · explicit_script_identifier
- [verification/v342_em_ward_heatkernel.py](../../../../verification/v342_em_ward_heatkernel.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [EM.WARD.02](../../../../verification/status_ledger.csv#L112) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v342 EM.WARD.02: the heat-kernel origin of the F_U(1) terms ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2260)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [(* ---- (v3) electromagnetic closure: alpha^-1 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L37)

## 2026-06-22 ( [formula in source] as the stationarity of a [formula in source] determinant line: every ingredient an index / heat-kernel coefficient / discriminant form, the open EM-Ward obligation named: v341)

[documents/changelog/sections/015-source.tex:260-278](../../../../documents/changelog/sections/015-source.tex#L260) · `section:656d00cdb6c9c469cdbe` · subsection

v341 ( ALPHA.QUILLEN.01). Answers the external review's sharpest valid point — ``derive the form of [formula in source] , not more decimals''. The EM closure [formula in source] is rewritten as a stationarity of the boundary [formula in source] determinant line: [formula in source] . \ every ingredient is a NAMED object, not a knob: [formula in source] is a [formula in source] index ( [formula in source] , [formula i

Explicit source references:

- [verification/v341_alpha_quillen.py](../../../../verification/v341_alpha_quillen.py#L1) · explicit_script_identifier
- [verification/v48_em_ward.py](../../../../verification/v48_em_ward.py#L1) · explicit_script_identifier
- [EM.WARD.01](../../../../verification/status_ledger.csv#L70) · explicit_claim_identifier · original status at linked ledger row
- [ALPHA.QUILLEN.01](../../../../verification/status_ledger.csv#L443) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v341 ALPHA.QUILLEN.01: the seam-form discriminant identities ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2249)
- [(* ---- (v48) EM boundary Ward decomposition (Theorem C) ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L191)

## 2026-06-22 (the first-principles boundary, mapped: can [formula in source] be computed? can [formula in source] be first-principles? both NO, with the exact reason: v339)

[documents/changelog/sections/015-source.tex:295-310](../../../../documents/changelog/sections/015-source.tex#L295) · `section:a850fdd649c27ed35b89` · subsection

v339 ( FIRSTPRINCIPLES.BOUNDARY.01). An honest boundary audit answering two head-on questions and fabricating nothing. \ The Planck mass cannot be computed — [formula in source] is forbidden by theorem (No-Unit, v153): a mass is mass-dimension [formula in source] while every TFPT datum is dimension [formula in source] , so no dimensionless map outputs an absolute scale. It is not a gap but over-determined (gravity [f

Explicit source references:

- [FIRSTPRINCIPLES.BOUNDARY.01](../../../../verification/status_ledger.csv#L441) · explicit_claim_identifier · original status at linked ledger row
- [verification/v35_quark_qcd_boundary.py](../../../../verification/v35_quark_qcd_boundary.py#L1) · explicit_script_identifier
- [verification/v339_firstprinciples_boundary.py](../../../../verification/v339_firstprinciples_boundary.py#L1) · explicit_script_identifier
- [verification/v153_no_unit_theorem.py](../../../../verification/v153_no_unit_theorem.py#L1) · explicit_script_identifier
- [verification/v246_unification_data_crosscheck.py](../../../../verification/v246_unification_data_crosscheck.py#L1) · explicit_script_identifier
- [verification/v262_fqcd_mp_me.py](../../../../verification/v262_fqcd_mp_me.py#L1) · explicit_script_identifier
- [verification/v274_scale_overdetermination.py](../../../../verification/v274_scale_overdetermination.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v153) No-Unit Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1369)

## 2026-06-22 (next-steps round: the one open lemma literature-anchored, the Decoupling Theorem, the holomorphy discriminator hardened in Lean, and the [formula in source] budget: v336, v337, v338, FORM.CARTAN.DET.01)

[documents/changelog/sections/015-source.tex:311-341](../../../../documents/changelog/sections/015-source.tex#L311) · `section:cfd99b18e189abbc1f87` · subsection

v336 ( SEAM.EQUIV.CONTINUUM.01) — the one open lemma, literature-anchored. SEAM.EQUIV.01's only open input (the continuum OS reconstruction of the gapped quasi-free collar) is split into a CITABLE part and a narrow open part. \ continuum leg: Osborne–Stottmeister ( CFT from Lattice Fermions, doi:10.1007/s00220-022-04521-8, CMP 398 (2023) 219) give, building on the MMST (Morinelli–Morsella–Stottmeister–Tanimoto) opera

Explicit source references:

- [verification/v328_theta13_pressure.py](../../../../verification/v328_theta13_pressure.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.CONTINUUM.01](../../../../verification/status_ledger.csv#L438) · explicit_claim_identifier · original status at linked ledger row
- [DECOUPLING.THEOREM.01](../../../../verification/status_ledger.csv#L439) · explicit_claim_identifier · original status at linked ledger row
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1) · explicit_file_citation
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_script_identifier
- [FORM.CARTAN.DET.01](../../../../verification/status_ledger.csv#L431) · explicit_claim_identifier · original status at linked ledger row
- [verification/v338_theta13_budget.py](../../../../verification/v338_theta13_budget.py#L1) · explicit_script_identifier
- [verification/v311_gap_decoupled_measure.py](../../../../verification/v311_gap_decoupled_measure.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v330_qgamb_admissible_measure.py](../../../../verification/v330_qgamb_admissible_measure.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v337_decoupling_theorem.py](../../../../verification/v337_decoupling_theorem.py#L1) · explicit_script_identifier
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [THETA13.BUDGET.01](../../../../verification/status_ledger.csv#L440) · explicit_claim_identifier · original status at linked ledger row
- [verification/v332_qgamb_metric_sector.py](../../../../verification/v332_qgamb_metric_sector.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ==== v337 DECOUPLING.THEOREM.01: the ambient back-reaction scale identity ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2241)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)

## 2026-06-22 (keystone unification + status reframe: ONE open theorem, ONE decoupled foreign problem: v335, FORM.QGEO.BW.01 qgeoSymIsCorollary)

[documents/changelog/sections/015-source.tex:342-373](../../../../documents/changelog/sections/015-source.tex#L342) · `section:2cf3ee9ae0a89a313d79` · subsection

v335 ( SEAM.EQUIV.UNIFY.01) — the keystone unification. A bird's-eye consolidation of the v308/v323/v325/v329/v333 arc that collapses the ``two open structural items'' narrative. \ QGEO.SYM.01 is a COROLLARY of SEAM.EQUIV.01 (no extra premise): a chiral conformal net has, by axiom, a M\"obius-covariant vacuum (the unique invariant positive-energy vector); rotations [formula in source] are the compact subgroup of the

Explicit source references:

- [FORM.QGEO.BW.01](../../../../verification/status_ledger.csv#L428) · explicit_claim_identifier · original status at linked ledger row
- [verification/v332_qgamb_metric_sector.py](../../../../verification/v332_qgamb_metric_sector.py#L1) · explicit_script_identifier
- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier
- [verification/v325_pillowcase_keystone.py](../../../../verification/v325_pillowcase_keystone.py#L1) · explicit_script_identifier
- [verification/v330_qgamb_admissible_measure.py](../../../../verification/v330_qgamb_admissible_measure.py#L1) · explicit_script_identifier
- [verification/v333_cm_selection.py](../../../../verification/v333_cm_selection.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.UNIFY.01](../../../../verification/status_ledger.csv#L437) · explicit_claim_identifier · original status at linked ledger row
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [experiments/next.txt](../../../../experiments/next.txt#L1) · explicit_file_citation
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_script_identifier
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [verification/v329_os_gap_reduction.py](../../../../verification/v329_os_gap_reduction.py#L1) · explicit_script_identifier
- [verification/v308_seam_equiv_chain.py](../../../../verification/v308_seam_equiv_chain.py#L1) · explicit_script_identifier
- [verification/v311_gap_decoupled_measure.py](../../../../verification/v311_gap_decoupled_measure.py#L1) · explicit_script_identifier
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1) · explicit_file_citation
- [TfptCarrier.BWKeystone.conformal_net_vacuum_rotation_invariant](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L87) · explicit_lean_symbol

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [(* ==== v325/v327: the keystone pillowcase orbifold + the cusp-weight rewrite (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2221)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* v325 QGEO.KEYSTONE.01: the flat pillowcase orbifold S^2(2,2,2,2) has chi_orb = 0 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2223)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)

## 2026-06-22 (closing round: the CM-selection mechanism (the last keystone residual), the conformal-mode resolution, and the Conway–Sloane even-unimodular side in Lean: v333, v334, FORM.CARTAN.DET.01)

[documents/changelog/sections/016-source.tex:1-22](../../../../documents/changelog/sections/016-source.tex#L1) · `section:ee232865bdd19669d0cd` · subsection

v333_cm_selection.py ( QGEO.CMSELECT.01, / ). The CM-selection mechanism for the last keystone residual: [formula in source] (the square) is the unique order-4 modulus — [formula in source] has [formula in source] and fixes [formula in source] , and an order-4 deck selects [formula in source] (vs the order-3 hexagon [formula in source] ) — AND the square is the Fekete attractor of symmetric mark-equilibration (four p

Explicit source references:

- [verification/v334_conformal_resolution.py](../../../../verification/v334_conformal_resolution.py#L1) · explicit_file_citation
- [QGEO.CMSELECT.01](../../../../verification/status_ledger.csv#L435) · explicit_claim_identifier · original status at linked ledger row
- [FORM.CARTAN.DET.01](../../../../verification/status_ledger.csv#L431) · explicit_claim_identifier · original status at linked ledger row
- [QGAMB.CONFORMAL.01](../../../../verification/status_ledger.csv#L436) · explicit_claim_identifier · original status at linked ledger row
- [verification/v333_cm_selection.py](../../../../verification/v333_cm_selection.py#L1) · explicit_script_identifier
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1) · explicit_file_citation
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row

## 2026-06-22 (endgame round: the necessity of [formula in source] reduced to the order-4 CM selection, the QG metric sector characterized as the conformal-mode problem, [formula in source] and the [formula in source] commutation grounded in Lean: v331, v332, FORM.MU4COMM.01)

[documents/changelog/sections/016-source.tex:23-48](../../../../documents/changelog/sections/016-source.tex#L23) · `section:8611c18776ca96eaa8ee` · subsection

v331_necessity_of_H.py ( QGEO.NECESS.01, / ). The necessity companion to v329: a concrete Steklov (DtN) computation shows the seam [formula in source] commutes with the [formula in source] clock iff the four marks sit at the square ( [formula in source] ) — moving a mark off the square breaks [formula in source] and [formula in source] . The square has cross-ratio [formula in source] ( [formula in source] ), the uniq

Explicit source references:

- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1) · explicit_file_citation
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v331_necessity_of_H.py](../../../../verification/v331_necessity_of_H.py#L1) · explicit_script_identifier
- [QGEO.NECESS.01](../../../../verification/status_ledger.csv#L432) · explicit_claim_identifier · original status at linked ledger row
- [verification/v332_qgamb_metric_sector.py](../../../../verification/v332_qgamb_metric_sector.py#L1) · explicit_file_citation
- [verification/v201_seam_subprincipal_marks.py](../../../../verification/v201_seam_subprincipal_marks.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [QGAMB.METRIC.01](../../../../verification/status_ledger.csv#L433) · explicit_claim_identifier · original status at linked ledger row
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [verification/v330_qgamb_admissible_measure.py](../../../../verification/v330_qgamb_admissible_measure.py#L1) · explicit_script_identifier
- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1) · explicit_file_citation
- [verification/v267_qgeo_rigidity_minimal_axiom.py](../../../../verification/v267_qgeo_rigidity_minimal_axiom.py#L1) · explicit_script_identifier
- [verification/v329_os_gap_reduction.py](../../../../verification/v329_os_gap_reduction.py#L1) · explicit_script_identifier
- [verification/v214_seam_pillowcase.py](../../../../verification/v214_seam_pillowcase.py#L1) · explicit_script_identifier
- [FORM.MU4COMM.01](../../../../verification/status_ledger.csv#L434) · explicit_claim_identifier · original status at linked ledger row
- [FORM.CARTAN.DET.01](../../../../verification/status_ledger.csv#L431) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [(* ---- (v201) QGEO.SUBPRIN.01: a mu4-mark sum is Z4-invariant => sub-principal block-diagonal ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1666)
- [(* ---- (v214) QGEO.PILLOW.01: pillowcase reduction -- cross-ratio 2 => j=1728 => order-4 CM ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1727)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)
- [(* ==== v267: QGEO rigidity / minimal-axiom form of QGEO.SYM.01 (QGEO.SYM.02) ==== v262 (alpha_s RG), v263 (seesaw numpy), v264 (DtN FFT), v265 (RG+text guard), v266 (proton RG) are numerical -> Python-only; v267's exact cross-ratio/j core is mirrored here (the DtN/FFT part of v267 stays Python-only). *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2052)

## 2026-06-22 (consistency + de-accretion editorial pass: current-state alignment across papers, website and ledger; no new modules)

[documents/changelog/sections/016-source.tex:109-138](../../../../documents/changelog/sections/016-source.tex#L109) · `section:ac6f5be914309816a43d` · subsection

Keystone naming reconciled to SEAM.EQUIV.01. Following v323/ QGEO.BW.01 (the deck premise is downstream of the seam being a chiral net), the single named open keystone is now SEAM.EQUIV.01 everywhere, with QGEO.SYM.01 stated as its downstream conformal-deck face. The ledger row QGEO.SYM.01 was updated (it no longer calls itself ``the single structural residual of the whole theory''; it now depends on SEAM.EQUIV.01/ Q

Explicit source references:

- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier
- [verification/v213_ftransfer_functor.py](../../../../verification/v213_ftransfer_functor.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [FORM.SEAMEQUIV.01](../../../../verification/status_ledger.csv#L423) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v324_hypergraph_fiber.py](../../../../verification/v324_hypergraph_fiber.py#L1) · explicit_script_identifier
- [QGEO.BW.01](../../../../verification/status_ledger.csv#L421) · explicit_claim_identifier · original status at linked ledger row
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)

## 2026-06-22 (Lean: the SEAM.EQUIV.01 composition chain formalised — the keystone residual machine-pinned to the named cited steps: FORM.SEAMEQUIV.01)

[documents/changelog/sections/016-source.tex:229-244](../../../../documents/changelog/sections/016-source.tex#L229) · `section:d2dce761d475a2f11d08` · subsection

SeamEquivChain.lean ( FORM.SEAMEQUIV.01, Lean 4, / ). A Lean 4 mirror of the v308 closing chain: the four nodes [formula in source] assembled as a well-typed composition theorem seamEquivChain, with the three cited literature steps (Kitaev free-fermion invertibility; M\"uger / Kawahigashi–Longo–M\"uger holomorphy; Conway–Sloane [formula in source] uniqueness) as named axioms and the carrier recovery gap (the OS step)

Explicit source references:

- [TfptCarrier.SeamEquivChain.recovery_gap](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L64) · explicit_lean_symbol
- [FORM.SEAMEQUIV.01](../../../../verification/status_ledger.csv#L423) · explicit_claim_identifier · original status at linked ledger row
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1) · explicit_file_citation
- [verification/v308_seam_equiv_chain.py](../../../../verification/v308_seam_equiv_chain.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)

## 2026-06-22 (the translation clock: the discrete [formula in source] dynamic bridge is the order- [formula in source] Coxeter element, a clock with two coprime hands [formula in source] , read law-inclusive ( [formula in source] ) or live-only ( [formula in source] ): v319)

[documents/changelog/sections/016-source.tex:260-283](../../../../documents/changelog/sections/016-source.tex#L260) · `section:cc398ba7e3e94c7addf2` · subsection

v319_translation_clock.py ( TRANSLATE.CLOCK.01, / ). A synthesis of the v313– v318 arc with v124/ v55, answering the question ``is the static [formula in source] dynamic translation a clock running [formula in source] or [formula in source] ?''. \ the only object holding both the static carrier ( [formula in source] ) and the dynamic family (rational) data is the order- [formula in source] Coxeter element ( v314); [f

Explicit source references:

- [verification/v318_arithmetic_capstone.py](../../../../verification/v318_arithmetic_capstone.py#L1) · explicit_script_identifier
- [verification/v319_translation_clock.py](../../../../verification/v319_translation_clock.py#L1) · explicit_script_identifier
- [verification/v315_coxeter_coupling.py](../../../../verification/v315_coxeter_coupling.py#L1) · explicit_script_identifier
- [verification/v316_galois_readout.py](../../../../verification/v316_galois_readout.py#L1) · explicit_script_identifier
- [TRANSLATE.CLOCK.01](../../../../verification/status_ledger.csv#L417) · explicit_claim_identifier · original status at linked ledger row
- [verification/v314_rate_translation.py](../../../../verification/v314_rate_translation.py#L1) · explicit_script_identifier
- [verification/v56_unique_attractor.py](../../../../verification/v56_unique_attractor.py#L1) · explicit_script_identifier
- [verification/v313_golden_atoms_spectral.py](../../../../verification/v313_golden_atoms_spectral.py#L1) · explicit_script_identifier
- [verification/v124_resummed_clock.py](../../../../verification/v124_resummed_clock.py#L1) · explicit_script_identifier
- [verification/v55_coxeter_cycle.py](../../../../verification/v55_coxeter_cycle.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v316 GALOIS.READOUT.01 -- CP phases are the family-factor cyclotomic data *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2206)
- [(* ==== v313-v320: the cyclotomic / Galois arithmetic arc (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2183)
- [(* v318 ARITH.CAPSTONE.01 -- the magnitude seed reduces to a pure pi-function *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2213)
- [(* v315 COX.COUPLE.01 -- the order-30 Coxeter sectors couple as a cyclotomic field *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2198)
- [(* ---- (v124) resummed clock ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L857)
- [(* v314 RATE.TRANSLATE.01 -- the number-field split between discrete kernel and dynamic rates *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2195)
- [(* ---- (v55) cyclic element in the hull: E8 Coxeter order 30 + one-alpha cross-sector ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L246)
- [(* v313 GOLD.ATOMS.01 -- the golden ratio is the g_car=5 spectral signature *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2185)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)
- [(* ---- (v56) unique attractor: gapped transport, Coxeter planes ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L257)

## 2026-06-21 (are the 3 generations a [formula in source] /Galois orbit? Yes ( [formula in source] ), Galois-refined [formula in source] with the fixed generation [formula in source] the recovery attractor: v317)

[documents/changelog/sections/016-source.tex:300-314](../../../../documents/changelog/sections/016-source.tex#L300) · `section:fa6e96662b31faa3de47` · subsection

v317_galois_family.py ( GALOIS.FAMILY.01, / ). The v316 follow-up to the generation structure. \ the three cusp weights [formula in source] (with [formula in source] , the [formula in source] exponents) are the cube roots of unity [formula in source] under [formula in source] , and the cyclic family symmetry [formula in source] is a transitive [formula in source] orbit (the family IS the cube-root orbit); \ the trans

Explicit source references:

- [verification/v317_galois_family.py](../../../../verification/v317_galois_family.py#L1) · explicit_file_citation
- [verification/v316_galois_readout.py](../../../../verification/v316_galois_readout.py#L1) · explicit_script_identifier
- [GALOIS.FAMILY.01](../../../../verification/status_ledger.csv#L414) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v317 GALOIS.FAMILY.01 -- the three generations are the mu3 orbit, Galois-refined 1+2 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2210)
- [(* v316 GALOIS.READOUT.01 -- CP phases are the family-factor cyclotomic data *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2206)

## 2026-06-21 (do the discrete-kernel and dynamic rates translate? An exact number-field split [formula in source] vs [formula in source] : v314)

[documents/changelog/sections/016-source.tex:348-364](../../../../documents/changelog/sections/016-source.tex#L348) · `section:a444fd35a837e55b2dcd` · subsection

v314_rate_translation.py ( RATE.TRANSLATE.01, / ). Answers, exactly, whether the discrete kernel (the affine- [formula in source] network, carrying the golden ratio) and the dynamic part (the seam transfer, carrying the recovery rate [formula in source] ) translate. \ FIELD SPLIT: the adjacency spectrum is the rational part [formula in source] plus the golden part [formula in source] , while the transfer spectrum [fo

Explicit source references:

- [verification/v314_rate_translation.py](../../../../verification/v314_rate_translation.py#L1) · explicit_script_identifier
- [RATE.TRANSLATE.01](../../../../verification/status_ledger.csv#L416) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v314 RATE.TRANSLATE.01 -- the number-field split between discrete kernel and dynamic rates *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2195)

## 2026-06-21 (TOE-roadmap solver round: the SEAM.EQUIV.01 keystone certificate, the modular bedrock, the carrier [formula in source] SM closure, the gap-decoupled measure, and the sharpened hypergraph substrate: v308– v312)

[documents/changelog/sections/017-source.tex:19-55](../../../../documents/changelog/sections/017-source.tex#L19) · `section:796200681c870006983a` · subsection

v308_seam_equiv_chain.py ( SEAM.EQUIV.CHAIN.01, / ). TOE attack point 1: the keystone ``the raw seam-Calder\'on net is the holomorphic [formula in source] [formula in source] net'' assembled into ONE well-typed logical chain [formula in source] , with the exact discriminators [formula in source] vs [formula in source] , the carrier tower [formula in source] and [formula in source] all machine-checked; the input gap i

Explicit source references:

- [QGEO.MODHAM.01](../../../../verification/status_ledger.csv#L407) · explicit_claim_identifier · original status at linked ledger row
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [HYP.INJECT.01](../../../../verification/status_ledger.csv#L410) · explicit_claim_identifier · original status at linked ledger row
- [verification/v308_seam_equiv_chain.py](../../../../verification/v308_seam_equiv_chain.py#L1) · explicit_script_identifier
- [verification/v312_hypergraph_substrate.py](../../../../verification/v312_hypergraph_substrate.py#L1) · explicit_script_identifier
- [QGAMB.CLUSTER.01](../../../../verification/status_ledger.csv#L409) · explicit_claim_identifier · original status at linked ledger row
- [verification/v309_modular_bedrock.py](../../../../verification/v309_modular_bedrock.py#L1) · explicit_file_citation
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v311_gap_decoupled_measure.py](../../../../verification/v311_gap_decoupled_measure.py#L1) · explicit_file_citation
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.CHAIN.01](../../../../verification/status_ledger.csv#L406) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v214_seam_pillowcase.py](../../../../verification/v214_seam_pillowcase.py#L1) · explicit_script_identifier
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [CAR.SM.ANOM.01](../../../../verification/status_ledger.csv#L408) · explicit_claim_identifier · original status at linked ledger row
- [verification/v310_carrier_sm_anomaly.py](../../../../verification/v310_carrier_sm_anomaly.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [(* ---- (v214) QGEO.PILLOW.01: pillowcase reduction -- cross-ratio 2 => j=1728 => order-4 CM ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1727)
- [experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/KronheimerMarks.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [experiments/lean4-carrier-rigidity/TfptCarrier/ArfSpinorCompiler.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/ArfSpinorCompiler.lean#L1)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)

## 2026-06-19 (an infinite-derivative narrowing of the Stelle-ghost attack, the double-pushout figure, and the ``constancy of [formula in source] '' framing of SEAM.EQUIV.01: v304)

[documents/changelog/sections/017-source.tex:92-115](../../../../documents/changelog/sections/017-source.tex#L92) · `section:355a87dfaea9de9cb7cc` · subsection

v304_idg_nonlocal_ghost.py ( QGAMB.IDG.01, / ). A conditional, parameter-free narrowing of the red-team rt_F/ v278 Stelle-ghost attack on the perturbative gravity sector. \ the local [formula in source] Weyl [formula in source] four-derivative propagator [formula in source] carries a negative-residue massive spin-2 mode (the Stelle ghost) only as a truncation artefact; \ for an entire graviton form factor [formula in

Explicit source references:

- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v278_lsz_bridge_unitarity.py](../../../../verification/v278_lsz_bridge_unitarity.py#L1) · explicit_script_identifier
- [verification/v92_glue_uniqueness.py](../../../../verification/v92_glue_uniqueness.py#L1) · explicit_script_identifier
- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [verification/v234_seam_holomorphy_selection.py](../../../../verification/v234_seam_holomorphy_selection.py#L1) · explicit_script_identifier
- [verification/v267_qgeo_rigidity_minimal_axiom.py](../../../../verification/v267_qgeo_rigidity_minimal_axiom.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_script_identifier
- [verification/v1_e8_glue.py](../../../../verification/v1_e8_glue.py#L1) · explicit_script_identifier
- [QGAMB.IDG.01](../../../../verification/status_ledger.csv#L400) · explicit_claim_identifier · original status at linked ledger row
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [(* ---- (v92) glue uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L151)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* ==== v267: QGEO rigidity / minimal-axiom form of QGEO.SYM.01 (QGEO.SYM.02) ==== v262 (alpha_s RG), v263 (seesaw numpy), v264 (DtN FFT), v265 (RG+text guard), v266 (proton RG) are numerical -> Python-only; v267's exact cross-ratio/j core is mirrored here (the DtN/FFT part of v267 stays Python-only). *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2052)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [(* ==== v234: Seam-Holomorphy selection -- one condition, three faces, E8 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1930)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [(* ---- (v1) E8 glue arithmetic ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L46)
- [(* ==== v278: S_pert -> S_phys LSZ bridge, one-loop optical theorem (QFT4D.SPERT.04) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2142)
- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)

## 2026-06-19 (the discrete [formula in source] dynamic lens on F_transfer: one shape, one seam-rate, v303)

[documents/changelog/sections/017-source.tex:116-132](../../../../documents/changelog/sections/017-source.tex#L116) · `section:0968f1b6845750a87eeb` · subsection

v303_ftransfer_dynamics.py asks whether the main-branch mechanism — a gapped, positivity-preserving relaxation to a unique attractor (the local-averaging update [formula in source] the [formula in source] marks, gap [formula in source] ) — is also the shape of the four F_transfer instances, or an unrelated bolt-on. It is the same shape: \ F_pole (Koide) is the seam dynamics (M\"obius multiplier [formula in source] ,

Explicit source references:

- [verification/v187_ftransfer_laws.py](../../../../verification/v187_ftransfer_laws.py#L1) · explicit_script_identifier
- [verification/v303_ftransfer_dynamics.py](../../../../verification/v303_ftransfer_dynamics.py#L1) · explicit_script_identifier
- [verification/v82_koide_attractor_splitting.py](../../../../verification/v82_koide_attractor_splitting.py#L1) · explicit_script_identifier
- [verification/v54_seam_horizon_keystones.py](../../../../verification/v54_seam_horizon_keystones.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v54) seam=horizon keystones: triply-forced 8 + shared transfer operator ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L235)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)

## 2026-06-19 (the last lever: the seam mass gap is the derived Recovery gap [formula in source] , v302)

[documents/changelog/sections/017-source.tex:133-151](../../../../documents/changelog/sections/017-source.tex#L133) · `section:5b9fc9c550f79525abcd` · subsection

v302_seam_gap.py closes the chain on the single statement v301 had left: ``the quasi-free seam bulk is gapped''. That gap is not a free input — it is the established Recovery gap of the seam transfer operator, derived from the carrier. \ the frozen transfer spectrum is [formula in source] ( v160/ v162), a simple Perron eigenvalue [formula in source] separated from the rest; \ the gap is [formula in source] (multiplic

Explicit source references:

- [verification/v302_seam_gap.py](../../../../verification/v302_seam_gap.py#L1) · explicit_file_citation
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [verification/v162_seam_transport_identification.py](../../../../verification/v162_seam_transport_identification.py#L1) · explicit_script_identifier
- [verification/v301_route_a_invertible.py](../../../../verification/v301_route_a_invertible.py#L1) · explicit_script_identifier
- [verification/v160_seam_gaussianity_from_pf.py](../../../../verification/v160_seam_gaussianity_from_pf.py#L1) · explicit_script_identifier
- [verification/v300_flataway_rigid.py](../../../../verification/v300_flataway_rigid.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v297_route_a_literature_stack.py](../../../../verification/v297_route_a_literature_stack.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)

## 2026-06-19 (Route A's last hypothesis discharged: the gapped quasi-free bulk is invertible by the free-fermion classification, v301)

[documents/changelog/sections/017-source.tex:152-170](../../../../documents/changelog/sections/017-source.tex#L152) · `section:41a8556f07331f818583` · subsection

v301_route_a_invertible.py attacks the single remaining unconditional gap of SEAM.EQUIV.01. Since v300 collapsed Route B into Route A's rationality, that gap is Route A's open hypothesis ( v297 LIT-A): ``the raw quasi-free seam bulk is invertible (SRE)''. A gapped free-fermion (Gaussian) [formula in source] D bulk cannot carry intrinsic topological order — free fermions give the integer (IQHE/Chern) phases, all inver

Explicit source references:

- [verification/v154_simple_current_theorem.py](../../../../verification/v154_simple_current_theorem.py#L1) · explicit_script_identifier
- [verification/v301_route_a_invertible.py](../../../../verification/v301_route_a_invertible.py#L1) · explicit_file_citation
- [verification/v297_route_a_literature_stack.py](../../../../verification/v297_route_a_literature_stack.py#L1) · explicit_script_identifier
- [verification/v160_seam_gaussianity_from_pf.py](../../../../verification/v160_seam_gaussianity_from_pf.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v148_fock_census.py](../../../../verification/v148_fock_census.py#L1) · explicit_script_identifier
- [verification/v300_flataway_rigid.py](../../../../verification/v300_flataway_rigid.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v148) NS/R sector census (corrected) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1281)
- [(* ---- (v154) Simple-Current Extension Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1381)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)

## 2026-06-19 (the Flat-Away attack: a hard discrete obstruction [formula in source] the pin derived from [formula in source] , v300)

[documents/changelog/sections/017-source.tex:171-189](../../../../documents/changelog/sections/017-source.tex#L171) · `section:2f80e058daba5055768b` · subsection

v300_flataway_rigid.py attacks the single external fact the three Flat-Away routes ( v291/ v293/ v295) left open, on two fronts. Harden: the flat seam was only a soft minimum (``the spectral mismatch grows with [formula in source] ''); here the obstruction becomes a discrete on/off invariant. \ the flat fingerprint — [formula in source] is the integer ladder with every level [formula in source] of multiplicity exactl

Explicit source references:

- [verification/v300_flataway_rigid.py](../../../../verification/v300_flataway_rigid.py#L1) · explicit_file_citation
- [verification/v295_flataway_a2_exact.py](../../../../verification/v295_flataway_a2_exact.py#L1) · explicit_script_identifier
- [verification/v293_flataway_spectral_hessian.py](../../../../verification/v293_flataway_spectral_hessian.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v291_flataway_contract.py](../../../../verification/v291_flataway_contract.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)

## 2026-06-18 (Flat-Away attacked on all three routes, v292– v294)

[documents/changelog/sections/017-source.tex:250-269](../../../../documents/changelog/sections/017-source.tex#L250) · `section:db81095caf4bfe7ff937` · subsection

v292_flataway_heat_reduction.py ( FLATAWAY.HEAT.01): the heat route, made precise. The heat-trace deviation [formula in source] is a positive-definite quadratic form in the off-mark curvature ( [formula in source] for every [formula in source] mode, positive on all combinations), leading invariant [formula in source] — so Flat-Away reduces to the single fact ``the seam fixes the [formula in source] heat coefficient''

Explicit source references:

- [verification/v293_flataway_spectral_hessian.py](../../../../verification/v293_flataway_spectral_hessian.py#L1) · explicit_file_citation
- [FLATAWAY.HEAT.01](../../../../verification/status_ledger.csv#L387) · explicit_claim_identifier · original status at linked ledger row
- [FLATAWAY.ENERGY.01](../../../../verification/status_ledger.csv#L389) · explicit_claim_identifier · original status at linked ledger row
- [verification/v292_flataway_heat_reduction.py](../../../../verification/v292_flataway_heat_reduction.py#L1) · explicit_file_citation
- [FLATAWAY.SPEC.01](../../../../verification/status_ledger.csv#L388) · explicit_claim_identifier · original status at linked ledger row
- [verification/v294_flataway_rp_energy.py](../../../../verification/v294_flataway_rp_energy.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)

## 2026-06-18 (Route-B red-team: the smooth [formula in source] adversary + Flat-Away as a mini-theorem, v290– v291)

[documents/changelog/sections/017-source.tex:270-285](../../../../documents/changelog/sections/017-source.tex#L270) · `section:1dd51adeb586c882e7c4` · subsection

v290_z4_smooth_curvature_adversary.py ( SEAM.ADVERSARY.01) is the honest red-team that [formula in source] block-diagonality is not mark-locality: a smooth [formula in source] off-mark curvature [formula in source] \ passes the v288 commutator ( [formula in source] , same as mark-local) yet \ shifts the DtN spectrum ( [formula in source] ) and the heat trace ( [formula in source] ). So the commutator is necessary, no

Explicit source references:

- [verification/v291_flataway_contract.py](../../../../verification/v291_flataway_contract.py#L1) · explicit_script_identifier
- [verification/v290_z4_smooth_curvature_adversary.py](../../../../verification/v290_z4_smooth_curvature_adversary.py#L1) · explicit_script_identifier
- [verification/v288_full_l2_subprincipal_z4.py](../../../../verification/v288_full_l2_subprincipal_z4.py#L1) · explicit_script_identifier
- [FLATAWAY.RP.01](../../../../verification/status_ledger.csv#L386) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.ADVERSARY.01](../../../../verification/status_ledger.csv#L385) · explicit_claim_identifier · original status at linked ledger row

## 2026-06-18 (Route B sharpened: raw mark-locality decomposed + reviewer firewalls, v289 / v287– v288)

[documents/changelog/sections/017-source.tex:286-303](../../../../documents/changelog/sections/017-source.tex#L286) · `section:f3f0a4f56db397726de9` · subsection

v289_marklocal_raw.py ( MARKLOCAL.RAW.01) decomposes the one remaining Route-B residual (`why is the raw seam sub-principal term mark-local?') into a 5-lemma chain and isolates the single open analytic lemma: \ L1 conic parametrix (b-calculus/Melrose, imported); \ L2 Flat-Away — RP [formula in source] gap [formula in source] 4 marks [formula in source] curvature vanishes away from the marks; \ L3 orbit-source ( v195/

Explicit source references:

- [verification/v288_full_l2_subprincipal_z4.py](../../../../verification/v288_full_l2_subprincipal_z4.py#L1) · explicit_script_identifier
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [verification/v289_marklocal_raw.py](../../../../verification/v289_marklocal_raw.py#L1) · explicit_file_citation
- [verification/v287_free_rp_bulk_to_holomorphic_boundary.py](../../../../verification/v287_free_rp_bulk_to_holomorphic_boundary.py#L1) · explicit_script_identifier
- [MARKLOCAL.RAW.01](../../../../verification/status_ledger.csv#L384) · explicit_claim_identifier · original status at linked ledger row
- [verification/v195_marks_lefschetz_character.py](../../../../verification/v195_marks_lefschetz_character.py#L1) · explicit_script_identifier
- [verification/v264_raw_seam_marklocal.py](../../../../verification/v264_raw_seam_marklocal.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v195) QGEO.MARKS.02: Lefschetz/character forcing of the free mu4 orbit ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1608)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)

## 2026-06-18 (closing sprint: the Seam Equivalence Theorem named + both routes attacked, v286– v288)

[documents/changelog/sections/017-source.tex:304-324](../../../../documents/changelog/sections/017-source.tex#L304) · `section:98e28332378f4d77a163` · subsection

v286_seam_equivalence_contract.py ( SEAM.EQUIV.01) names the one open theorem — the raw RP seam state is the holomorphic [formula in source] net at [formula in source] — and concentrates the whole open structure on it: four objects, six arrows (2 closed, 1 conditional v282, 1 open Gral), an import firewall proving the two routes ( v276/ v280/ v284 vs v277/ v281/ v285) never import each other (no ` [formula in source]

Explicit source references:

- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [verification/v276_qgeo_flat_closes_commutator.py](../../../../verification/v276_qgeo_flat_closes_commutator.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v287_free_rp_bulk_to_holomorphic_boundary.py](../../../../verification/v287_free_rp_bulk_to_holomorphic_boundary.py#L1) · explicit_file_citation
- [verification/v288_full_l2_subprincipal_z4.py](../../../../verification/v288_full_l2_subprincipal_z4.py#L1) · explicit_file_citation
- [verification/v201_seam_subprincipal_marks.py](../../../../verification/v201_seam_subprincipal_marks.py#L1) · explicit_script_identifier
- [verification/v285_route_ii_seam_condensation.py](../../../../verification/v285_route_ii_seam_condensation.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.B01](../../../../verification/status_ledger.csv#L383) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.A01](../../../../verification/status_ledger.csv#L382) · explicit_claim_identifier · original status at linked ledger row
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_script_identifier
- [verification/seam_equivalence.json](../../../../verification/seam_equivalence.json#L1) · explicit_file_citation
- [verification/v284_route_i_rp_uniqueness.py](../../../../verification/v284_route_i_rp_uniqueness.py#L1) · explicit_script_identifier
- [verification/v286_seam_equivalence_contract.py](../../../../verification/v286_seam_equivalence_contract.py#L1) · explicit_script_identifier
- [verification/v280_pillowcase_steklov.py](../../../../verification/v280_pillowcase_steklov.py#L1) · explicit_script_identifier
- [verification/v282_e8_tau_i_unification.py](../../../../verification/v282_e8_tau_i_unification.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)
- [(* ---- (v201) QGEO.SUBPRIN.01: a mu4-mark sum is Z4-invariant => sub-principal block-diagonal ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1666)
- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [(* ==== v282: E8-at-tau=i unification, chi_E8(i)=12 (QGAMB.UNIFY.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2171)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)

## 2026-06-18 (self-investigation round: pillowcase DtN, modular data, the E8-at- [formula in source] unification, live scorecard, v280– v283)

[documents/changelog/sections/017-source.tex:341-359](../../../../documents/changelog/sections/017-source.tex#L341) · `section:5573768c07b5b0367bc7` · subsection

v280_pillowcase_steklov.py ( QGEO.STEKLOV.01) runs the QGEO chain on the actual flat [formula in source] pillowcase orbifold: the order-4 deck permutes the 4 cone points (2 fixed [formula in source] 1 swap) and the geometric chain [formula in source] is numerically exact — route-(i) evidence for QGEO.OBLIG.01; the premise stays . v281_holomorphy_modular_data.py ( QGAMB.MODULAR.01): \#anyons [formula in source] (E8 [f

Explicit source references:

- [verification/v280_pillowcase_steklov.py](../../../../verification/v280_pillowcase_steklov.py#L1) · explicit_script_identifier
- [verification/v283_live_killtest_scorecard.py](../../../../verification/v283_live_killtest_scorecard.py#L1) · explicit_file_citation
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_file_citation
- [QGEO.OBLIG.01](../../../../verification/status_ledger.csv#L373) · explicit_claim_identifier · original status at linked ledger row
- [QGAMB.MODULAR.01](../../../../verification/status_ledger.csv#L376) · explicit_claim_identifier · original status at linked ledger row
- [PRED.LIVE.01](../../../../verification/status_ledger.csv#L378) · explicit_claim_identifier · original status at linked ledger row
- [QGAMB.UNIFY.01](../../../../verification/status_ledger.csv#L377) · explicit_claim_identifier · original status at linked ledger row
- [verification/v282_e8_tau_i_unification.py](../../../../verification/v282_e8_tau_i_unification.py#L1) · explicit_script_identifier
- [QGEO.STEKLOV.01](../../../../verification/status_ledger.csv#L375) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [(* ==== v282: E8-at-tau=i unification, chi_E8(i)=12 (QGAMB.UNIFY.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2171)

## 2026-06-18 (the Gr\"al: QGEO.SYM.01 as one precise lemma [formula in source] the all-orders step Lean-formalised, v279/ FORM.QGEO.04)

[documents/changelog/sections/017-source.tex:360-377](../../../../documents/changelog/sections/017-source.tex#L360) · `section:47adde543bdad0ea081f` · subsection

v279_qgeo_obligation_lemma.py ( QGEO.OBLIG.01) writes the bedrock as ONE precise constructive-QFT lemma. Given the premise the raw seam carries the flat [formula in source] pillowcase metric (Troyanov-unique, [formula in source] ), then [formula in source] mark-locality [formula in source] . \ the reduction DAG from `free RP seam' has exactly one open leaf — `the raw seam state is the flat state'; \ the geometric sid

Explicit source references:

- [TFPT.Carrier.SeamDeckClosure.flat_all_orders_clock](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L202) · explicit_lean_symbol
- [verification/v276_qgeo_flat_closes_commutator.py](../../../../verification/v276_qgeo_flat_closes_commutator.py#L1) · explicit_script_identifier
- [FORM.QGEO.04](../../../../verification/status_ledger.csv#L374) · explicit_claim_identifier · original status at linked ledger row
- [verification/v279_qgeo_obligation_lemma.py](../../../../verification/v279_qgeo_obligation_lemma.py#L1) · explicit_file_citation
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1) · explicit_file_citation
- [QGEO.OBLIG.01](../../../../verification/status_ledger.csv#L373) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [TFPT.Carrier.SeamDeckClosure.diagonal_commutesClock](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L182) · explicit_lean_symbol

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)

## 2026-06-18 (QG.AMB.01 Tier-B made concrete: the seam-Calder\'on [formula in source] match, v277)

[documents/changelog/sections/018-source.tex:15-26](../../../../documents/changelog/sections/018-source.tex#L15) · `section:f1a0b6f25c4c201d4774` · subsection

v277_seam_calderon_e8_match.py ( QGAMB.TIERB.01) reduces Tier-B to ONE bit. \ the full [formula in source] matching target is computed and matches: [formula in source] , [formula in source] , the holomorphic character [formula in source] with a single primary and exactly [formula in source] weight-1 currents ( [formula in source] roots). \ the discriminator is holomorphy — the [formula in source] counterexample [form

Explicit source references:

- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [QGAMB.TIERB.01](../../../../verification/status_ledger.csv#L371) · explicit_claim_identifier · original status at linked ledger row
- [verification/v276_qgeo_flat_closes_commutator.py](../../../../verification/v276_qgeo_flat_closes_commutator.py#L1) · explicit_script_identifier
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)

## 2026-06-18 (QGEO.SYM.01 narrowed: the flat pillowcase closes the commutator to all orders, v276)

[documents/changelog/sections/018-source.tex:27-40](../../../../documents/changelog/sections/018-source.tex#L27) · `section:25e03edabfb4b9161f7f` · subsection

v276_qgeo_flat_closes_commutator.py ( QGEO.SYM.03) collapses the bedrock to ONE geometric postulate. The state-level residual the whole bedrock rests on is the operator equation [formula in source] . \ if [formula in source] and the order-4 deck [formula in source] ( [formula in source] ) is a flat- [formula in source] isometry, then [formula in source] exactly, to all orders (verified on the discretised square torus

Explicit source references:

- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v201_seam_subprincipal_marks.py](../../../../verification/v201_seam_subprincipal_marks.py#L1) · explicit_script_identifier
- [verification/v276_qgeo_flat_closes_commutator.py](../../../../verification/v276_qgeo_flat_closes_commutator.py#L1) · explicit_file_citation
- [verification/v198_modular_commutator_reduction.py](../../../../verification/v198_modular_commutator_reduction.py#L1) · explicit_script_identifier
- [QGEO.SYM.03](../../../../verification/status_ledger.csv#L370) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v198) QGEO.MODULAR.01: principal symbol |k| commutes with the clock rho EXACTLY ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1639)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)
- [(* ---- (v201) QGEO.SUBPRIN.01: a mu4-mark sum is Z4-invariant => sub-principal block-diagonal ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1666)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)

## 2026-06-18 (The bedrock in its minimal-axiom form: a rigidity theorem for QGEO.SYM.01, v267)

[documents/changelog/sections/018-source.tex:160-175](../../../../documents/changelog/sections/018-source.tex#L160) · `section:d22159cf08f411c39c54` · subsection

v267_qgeo_rigidity_minimal_axiom.py ( QGEO.SYM.02) sharpens the one bedrock postulate to its irreducible kernel (it does not close it). \ the single bare statement ``the seam admits the carrier order- [formula in source] clock as a conformal symmetry permuting its four marks'' ( [formula in source] ) forces everything: an order- [formula in source] M\"obius orbit [formula in source] has cross-ratio [formula in source

Explicit source references:

- [verification/v267_qgeo_rigidity_minimal_axiom.py](../../../../verification/v267_qgeo_rigidity_minimal_axiom.py#L1) · explicit_file_citation
- [QGEO.SYM.02](../../../../verification/status_ledger.csv#L360) · explicit_claim_identifier · original status at linked ledger row
- [verification/v198_modular_commutator_reduction.py](../../../../verification/v198_modular_commutator_reduction.py#L1) · explicit_script_identifier
- [verification/v214_seam_pillowcase.py](../../../../verification/v214_seam_pillowcase.py#L1) · explicit_script_identifier
- [verification/v201_seam_subprincipal_marks.py](../../../../verification/v201_seam_subprincipal_marks.py#L1) · explicit_script_identifier
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v264_raw_seam_marklocal.py](../../../../verification/v264_raw_seam_marklocal.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ==== v267: QGEO rigidity / minimal-axiom form of QGEO.SYM.01 (QGEO.SYM.02) ==== v262 (alpha_s RG), v263 (seesaw numpy), v264 (DtN FFT), v265 (RG+text guard), v266 (proton RG) are numerical -> Python-only; v267's exact cross-ratio/j core is mirrored here (the DtN/FFT part of v267 stays Python-only). *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2052)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [(* ---- (v214) QGEO.PILLOW.01: pillowcase reduction -- cross-ratio 2 => j=1728 => order-4 CM ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1727)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)
- [(* ---- (v198) QGEO.MODULAR.01: principal symbol |k| commutes with the clock rho EXACTLY ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1639)
- [(* ---- (v201) QGEO.SUBPRIN.01: a mu4-mark sum is Z4-invariant => sub-principal block-diagonal ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1666)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)

## 2026-06-18 (The 4D-QFT fork frozen as a branch policy: boundary-only canonical, Pati–Salam a falsifiable UV branch, SM-only GUT killed, v265)

[documents/changelog/sections/018-source.tex:191-208](../../../../documents/changelog/sections/018-source.tex#L191) · `section:0e1e7c50a8be56e30768` · subsection

v265_qft4d_fork_freeze.py ( QFT4D.FORK.01) turns the 4D-QFT fork from an open ambiguity into a frozen decision tree ( qft4d_fork_freeze.json), no new physics. A (canonical): boundary-only is the default 4D reading — a 4D-GUT is not claimed by default ( [formula in source] is the audit hull). \ SM-only 4D-GUT unification is killed (1-loop: at the [formula in source] scale [formula in source] off by [formula in source]

Explicit source references:

- [verification/v253_kappa_scalaron_ps.py](../../../../verification/v253_kappa_scalaron_ps.py#L1) · explicit_script_identifier
- [verification/v265_qft4d_fork_freeze.py](../../../../verification/v265_qft4d_fork_freeze.py#L1) · explicit_script_identifier
- [verification/v249_ps_unification.py](../../../../verification/v249_ps_unification.py#L1) · explicit_script_identifier
- [QFT4D.BOUNDARY.DEFAULT.01](../../../../verification/status_ledger.csv#L352) · explicit_claim_identifier · original status at linked ledger row
- [verification/v247_e8_branching_no126.py](../../../../verification/v247_e8_branching_no126.py#L1) · explicit_script_identifier
- [PS.THRESHOLD.01](../../../../verification/status_ledger.csv#L355) · explicit_claim_identifier · original status at linked ledger row
- [PS.UVBRANCH.01](../../../../verification/status_ledger.csv#L354) · explicit_claim_identifier · original status at linked ledger row
- [PS.PROTON.KILL.01](../../../../verification/status_ledger.csv#L356) · explicit_claim_identifier · original status at linked ledger row
- [verification/qft4d_fork_freeze.json](../../../../verification/qft4d_fork_freeze.json#L1) · explicit_file_citation
- [verification/v246_unification_data_crosscheck.py](../../../../verification/v246_unification_data_crosscheck.py#L1) · explicit_script_identifier
- [QFT4D.FORK.01](../../../../verification/status_ledger.csv#L351) · explicit_claim_identifier · original status at linked ledger row
- [verification/v248_ps_algebra.py](../../../../verification/v248_ps_algebra.py#L1) · explicit_script_identifier
- [QFT4D.NO_OVERCLAIM.01](../../../../verification/status_ledger.csv#L357) · explicit_claim_identifier · original status at linked ledger row
- [QFT4D.SMONLYGUT.01](../../../../verification/status_ledger.csv#L353) · explicit_claim_identifier · original status at linked ledger row

## 2026-06-18 (Three open-residual advances: the [formula in source] solver, [formula in source] from [formula in source] , and the raw-seam [formula in source] mark-locality reduction, v262– v264)

[documents/changelog/sections/018-source.tex:209-235](../../../../documents/changelog/sections/018-source.tex#L209) · `section:4a3638ac1b45d771412e` · subsection

v262_fqcd_mp_me.py ( FR.MPME.02) implements the one contract-only [formula in source] pipeline. A real two-loop [formula in source] solver runs [formula in source] with [formula in source] thresholds (carrier [formula in source] , ) to [formula in source] and [formula in source] \,GeV; with the lattice bridge [formula in source] (external [formula in source] ) and the closed [formula in source] it reproduces [formula

Explicit source references:

- [verification/v262_fqcd_mp_me.py](../../../../verification/v262_fqcd_mp_me.py#L1) · explicit_script_identifier
- [FR.MPME.01](../../../../verification/status_ledger.csv#L49) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v198_modular_commutator_reduction.py](../../../../verification/v198_modular_commutator_reduction.py#L1) · explicit_script_identifier
- [verification/v264_raw_seam_marklocal.py](../../../../verification/v264_raw_seam_marklocal.py#L1) · explicit_file_citation
- [FR.MPME.02](../../../../verification/status_ledger.csv#L348) · explicit_claim_identifier · original status at linked ledger row
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [verification/v231_cp_mu6_phases.py](../../../../verification/v231_cp_mu6_phases.py#L1) · explicit_script_identifier
- [QGEO.ENERGY.03](../../../../verification/status_ledger.csv#L350) · explicit_claim_identifier · original status at linked ledger row
- [verification/v233_cp_triality_phase.py](../../../../verification/v233_cp_triality_phase.py#L1) · explicit_script_identifier
- [verification/v201_seam_subprincipal_marks.py](../../../../verification/v201_seam_subprincipal_marks.py#L1) · explicit_script_identifier
- [verification/v263_mnu_from_df.py](../../../../verification/v263_mnu_from_df.py#L1) · explicit_file_citation
- [FLAV.PMNS.02](../../../../verification/status_ledger.csv#L349) · explicit_claim_identifier · original status at linked ledger row
- [verification/v9_neutrino_texture.py](../../../../verification/v9_neutrino_texture.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v231: both CP phases are mu6 powers of one hexagonal unit, split by the sheet ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1875)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [(* ---- (v198) QGEO.MODULAR.01: principal symbol |k| commutes with the clock rho EXACTLY ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1639)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ==== v233: CP = the universal family/triality phase, sheet-split (mu6 = mu3 x mu2) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1916)
- [(* ---- (v201) QGEO.SUBPRIN.01: a mu4-mark sum is Z4-invariant => sub-principal block-diagonal ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1666)

## 2026-06-18 (Paper consistency pass: the closure ledger and the status/QFT sections updated for the Modular Spectral Closure, with an explicit ``what is not closed'' line)

[documents/changelog/sections/018-source.tex:236-254](../../../../documents/changelog/sections/018-source.tex#L236) · `section:817ba5dfc4c0aadf9e21` · subsection

Five paper sections that predated the v238– v261 round are brought up to date (additive, no status marker moved). The introduction closure-ledger table and the ``five questions'' item now carry a boundary QFT (Modular Spectral Closure) row — one relative object [formula in source] closed modulo QGEO.SYM.01 ( v258– v261) — plus an explicit list of what the closure does not reach: the full measured pole masses, the abs

Explicit source references:

- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_script_identifier
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_script_identifier
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v238_modular_lindblad_dynamics.py](../../../../verification/v238_modular_lindblad_dynamics.py#L1) · explicit_script_identifier
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_script_identifier
- [verification/v255_spectral_action_expansion.py](../../../../verification/v255_spectral_action_expansion.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)

## 2026-06-18 (Consistency pass: the Modular Spectral Closure surfaced structurally across the website and the intro status card)

[documents/changelog/sections/018-source.tex:255-267](../../../../documents/changelog/sections/018-source.tex#L255) · `section:fe4822aaeee66083bf56` · subsection

The v258– v261 closure is now reflected on every ``what is open'' surface, not just the research-contracts paper. The homepage open-gates section and compiler-pipeline diagram, the hostile-referee FAQ, the glossary [formula in source] entry, the reviewer page and the verification residual-chain caption now state that the entire emergent-QFT layer (Dirac, cutoff, gauging, glue, orientability) collapses onto the same p

Explicit source references:

- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_script_identifier
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_script_identifier
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row

## 2026-06-17 (The Modular Spectral Closure: Dirac as covariance induction, the modular cutoff, the Kummer/K3 unification, and the assembly certificate, v258– v261)

[documents/changelog/sections/018-source.tex:268-310](../../../../documents/changelog/sections/018-source.tex#L268) · `section:ca7f6aced60b17ae6709` · subsection

v258_dirac_covariance_induction.py ( PS.DIRAC.03) closes the ``posited Yukawa vs derived operator'' gap of v250/ v252. The finite Dirac operator [formula in source] is the modular/covariance induction of the boundary seam quasi-free (KMS, [formula in source] ) state: a Gaussian state is fixed by its covariance [formula in source] (a positive contraction) with one-particle modular Hamiltonian [formula in source] (Arak

Explicit source references:

- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_file_citation
- [verification/v239_kms_thermal_time.py](../../../../verification/v239_kms_thermal_time.py#L1) · explicit_script_identifier
- [PS.SPECACT.02](../../../../verification/status_ledger.csv#L345) · explicit_claim_identifier · original status at linked ledger row
- [verification/v214_seam_pillowcase.py](../../../../verification/v214_seam_pillowcase.py#L1) · explicit_script_identifier
- [PS.DIRAC.03](../../../../verification/status_ledger.csv#L344) · explicit_claim_identifier · original status at linked ledger row
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_script_identifier
- [verification/v1_e8_glue.py](../../../../verification/v1_e8_glue.py#L1) · explicit_script_identifier
- [verification/v252_full_finite_triple.py](../../../../verification/v252_full_finite_triple.py#L1) · explicit_script_identifier
- [verification/v232_e8_kleinian_seam.py](../../../../verification/v232_e8_kleinian_seam.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v18_quark_yukawa.py](../../../../verification/v18_quark_yukawa.py#L1) · explicit_script_identifier
- [verification/v197_rr_carrier_clifford_d5.py](../../../../verification/v197_rr_carrier_clifford_d5.py#L1) · explicit_script_identifier
- [verification/v253_kappa_scalaron_ps.py](../../../../verification/v253_kappa_scalaron_ps.py#L1) · explicit_script_identifier
- [verification/v58_seam_horizon_chain.py](../../../../verification/v58_seam_horizon_chain.py#L1) · explicit_script_identifier
- [verification/v238_modular_lindblad_dynamics.py](../../../../verification/v238_modular_lindblad_dynamics.py#L1) · explicit_script_identifier
- [verification/v255_spectral_action_expansion.py](../../../../verification/v255_spectral_action_expansion.py#L1) · explicit_script_identifier
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_file_citation
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_script_identifier
- [QFT.MSC.01](../../../../verification/status_ledger.csv#L347) · explicit_claim_identifier · original status at linked ledger row
- [verification/v236_brieskorn_capstone.py](../../../../verification/v236_brieskorn_capstone.py#L1) · explicit_script_identifier
- [verification/v250_dirac_triple_contract.py](../../../../verification/v250_dirac_triple_contract.py#L1) · explicit_script_identifier
- [ARCH.K3.01](../../../../verification/status_ledger.csv#L346) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v197) ARCH.RRCAR.02: even Clifford of the 5-dim carrier = the D5 half-spinor ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1631)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [(* ---- (v58) seam-horizon chain: one-sided Gauss-Bonnet, KMS, S_BH, RT in seam units ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L270)
- [(* ==== v232: the seam as the E8 Kleinian singularity (du Val resolution + link) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1895)
- [(* ==== v236: the (2,3,5) Brieskorn singularity generates the skeleton ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1975)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [(* ---- (v1) E8 glue arithmetic ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L46)
- [(* ---- (v214) QGEO.PILLOW.01: pillowcase reduction -- cross-ratio 2 => j=1728 => order-4 CM ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1727)
