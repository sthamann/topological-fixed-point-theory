# Archives and dated records — part 25

[Topic entry](../history.md) · [Registry entry](../../INDEX.md)

Exact ledger statuses remain at the linked original rows and in registry.json. Citations are source links, not proof certification.

## 2026-06-12 (tfpt_1 / tfpt_2 style alignment + box reduction)

[documents/changelog/sections/023-source.tex:221-241](../../../../documents/changelog/sections/023-source.tex#L221) · `section:5cb9ec06306c576c2af2` · subsection

Shared header/footer everywhere. tfpt_1_architecture_e8 dropped its custom footmarkers footer (which after the marker collapse showed four redundant ) and now inherits the shared tfpt_style header/footer with the version stamp, identical to every other document. Fewer boxes, embedded in text. The narrative/commentary callouts of both documents are converted to run-in bold paragraphs (and the boxed one-line formulas t

Explicit source references:

- [COSMO.NSTAR.01](../../../../verification/status_ledger.csv#L160) · explicit_claim_identifier · original status at linked ledger row
- [verification/v86_nstar_reheating.py](../../../../verification/v86_nstar_reheating.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v86) reheating arithmetic (scale chain; pivot solve stays Python) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L84)

## 2026-06-12 (marker simplification to 4 classes + status reconciliation)

[documents/changelog/sections/023-source.tex:242-263](../../../../documents/changelog/sections/023-source.tex#L242) · `section:4f6387455280982ea775` · subsection

Four-class display markers. The reader-facing marker system is collapsed from eleven symbols to four: \ exact / machine-proven (identity, Lie/lattice, formalised, numerical FP), \ conditional (physical / bridge / readout under named hypotheses), \ open / axiom, and \ kill test. tex-artefacts/tfpt_style.tex defines mE/ mC/ mO/ mX and keeps the legacy tI tX names as aliases that render their class. All marker call-site

Explicit source references:

- [COSMO.NSTAR.01](../../../../verification/status_ledger.csv#L160) · explicit_claim_identifier · original status at linked ledger row
- [verification/v86_nstar_reheating.py](../../../../verification/v86_nstar_reheating.py#L1) · explicit_script_identifier
- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v86) reheating arithmetic (scale chain; pivot solve stays Python) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L84)

## 2026-06-12 (version system, predictions/K5, readability pass B3)

[documents/changelog/sections/023-source.tex:264-302](../../../../documents/changelog/sections/023-source.tex#L264) · `section:9b652b6d86057eb55f7d` · subsection

Authors, auto-date and auto-version on every paper. New single source tex-artefacts/version.tex defines TFPTauthors (Stefan Hamann, Alessandro Rizzo), a curated TFPTversion, and a title/footer date stamp ( today [formula in source] version [formula in source] revision). build.sh writes the auto build revision (git commit count) into the gitignored tex-artefacts/version-auto.tex, so every document's title page and pag

Explicit source references:

- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier
- [verification/v1_e8_glue.py](../../../../verification/v1_e8_glue.py#L1) · explicit_script_identifier
- [verification/v20_lepton_c_derivation.py](../../../../verification/v20_lepton_c_derivation.py#L1) · explicit_script_identifier
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier
- [verification/v18_quark_yukawa.py](../../../../verification/v18_quark_yukawa.py#L1) · explicit_script_identifier
- [verification/v23_anchor_generator.py](../../../../verification/v23_anchor_generator.py#L1) · explicit_script_identifier
- [verification/v14_carrier_uniqueness.py](../../../../verification/v14_carrier_uniqueness.py#L1) · explicit_script_identifier
- [verification/v6_bootstrap.py](../../../../verification/v6_bootstrap.py#L1) · explicit_script_identifier
- [verification/v108_pascal_ladder.py](../../../../verification/v108_pascal_ladder.py#L1) · explicit_script_identifier
- [verification/v113_quasifree_kernel.py](../../../../verification/v113_quasifree_kernel.py#L1) · explicit_script_identifier
- [verification/v46_grand_mass_volume.py](../../../../verification/v46_grand_mass_volume.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v20/v21) lepton c's: delta=1/2 resolvent + product rule ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L78)
- [(* ---- (v46) Grand Mass Volume: absolute mass-scaling = det ~ u^25 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L178)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean#L1)
- [(* ---- (v1) E8 glue arithmetic ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L46)
- [(* ---- (v108) Pascal ladder ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L486)
- [experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean#L1)
- [(* ---- (v23) anchor generator a=(1,1,2) ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L96)
- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)
- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [(* ---- (v113) quasi-free kernel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L608)

## 2026-06-12 (tex-artefacts refactor + introduction visual pass B2)

[documents/changelog/sections/023-source.tex:303-330](../../../../documents/changelog/sections/023-source.tex#L303) · `section:826ddf772de04049b71a` · subsection

Shared imports moved to tex-artefacts/. The shared, imported-only TeX fragments ( tfpt_style.tex, tfpt_figures.tex, tfpt_docset.tex, tfpt_status.tex) now live in a dedicated tex-artefacts/ folder; all nine active documents input them via tex-artefacts/<name> (and tfpt_style loads tex-artefacts/tfpt_figures). changelog.tex stays at the repo root because it is also a standalone build.sh notes deliverable. The make_mani

Explicit source references:

- [README.md](../../../../README.md#L1) · explicit_file_citation
- [verification/make_manifest.py](../../../../verification/make_manifest.py#L1) · explicit_file_citation

## 2026-06-12 (visual pass: box reduction B1)

[documents/changelog/sections/023-source.tex:331-340](../../../../documents/changelog/sections/023-source.tex#L331) · `section:37c689285c5fc6a40772` · subsection

Box reduction (B1). The [formula in source] non-load-bearing aside boxes (audit / recorded / interpretation: auditbox/ audbox/ intbox/ resbox/ roadbox) across tfpt_1, tfpt_2, tfpt_3 and origin_theory are converted to inline bold run-in paragraphs — the reviewer's prescription that an audit fingerprint must not carry the same visual volume as a theorem. The loud box families are now only keybox (claim), winbox (result

## 2026-06-12 (visual pass: shared figures, badges, marker typology)

[documents/changelog/sections/023-source.tex:341-355](../../../../documents/changelog/sections/023-source.tex#L341) · `section:58cade7e0ae49570c7f2` · subsection

Shared figure library tfpt_figures.tex ( input by tfpt_style): reusable TikZ macros for the architecture claim stack ( TFPTclaimstack[zone], B3), the anchor cockpit [formula in source] (B4), the [formula in source] glue diagram (B5), the [formula in source] motif map (B6) and the [formula in source] action lattice (B7). The claim stack and a marker badge bar are now shown near the front of every document, each highli

## 2026-06-12 (review pass on list.txt)

[documents/changelog/sections/023-source.tex:356-381](../../../../documents/changelog/sections/023-source.tex#L356) · `section:b7fc94b9c310d988cbec` · subsection

Claim hygiene (review section A). The introduction's one-sentence claim is softened from ``derives all constants'' to ``constructs a discrete compiler; physical readouts run through named, status-typed bridges''; the ledger zombie in FLAV.QRATIO.01 is removed (the statement no longer says ``stays [P]'' — ratios are \ closed on the derived stratum, only the absolute scale stays , with an explicit forbidden-wording not

Explicit source references:

- [FLAV.QRATIO.01](../../../../verification/status_ledger.csv#L25) · explicit_claim_identifier · original status at linked ledger row

## 2026-06-12 (later)

[documents/changelog/sections/024-source.tex:1-11](../../../../documents/changelog/sections/024-source.tex#L1) · `section:e2e448f834c8b18c1b82` · subsection

Shared style and central status artifact. Two new single-source fragments: tfpt_style.tex (the canonical palette, status markers, marker key, veri reference, the four tcolorbox families with back-compat aliases, and one running header/footer) is now input by all eight documents — the per-paper duplicated preambles are gone; and tfpt_status.tex (``TFPT status at a glance'') gives one canonical overall-status card impo

## 2026-06-12

[documents/changelog/sections/024-source.tex:12-75](../../../../documents/changelog/sections/024-source.tex#L12) · `section:a7ddfbed7669a08bf2ee` · subsection

Margin theorem (v122). The established frozen selectors ( [formula in source] annihilator [formula in source] , [formula in source] , [formula in source] ) pin [formula in source] uniquely among hexagon matrices ( [formula in source] ); the atom margins become theorems — H2 introduces no new residual class. Inventory update (v123). Residual table re-pinned post v110–v122: exactly five structural classes (R1 clock, R2

Explicit source references:

- [verification/v138_vw_firewall.py](../../../../verification/v138_vw_firewall.py#L1) · explicit_script_identifier
- [verification/v140_canonical_map.py](../../../../verification/v140_canonical_map.py#L1) · explicit_script_identifier
- [verification/v126_clock_wall_bridge.py](../../../../verification/v126_clock_wall_bridge.py#L1) · explicit_script_identifier
- [verification/v134_dual_anchor.py](../../../../verification/v134_dual_anchor.py#L1) · explicit_script_identifier
- [verification/v128_graded_hull.py](../../../../verification/v128_graded_hull.py#L1) · explicit_script_identifier
- [verification/v129_entropy_power_law.py](../../../../verification/v129_entropy_power_law.py#L1) · explicit_script_identifier
- [verification/v133_zeta_budget.py](../../../../verification/v133_zeta_budget.py#L1) · explicit_script_identifier
- [verification/v127_ring_resummation.py](../../../../verification/v127_ring_resummation.py#L1) · explicit_script_identifier
- [verification/v122_margin_theorem.py](../../../../verification/v122_margin_theorem.py#L1) · explicit_script_identifier
- [verification/v131_measure_is_area.py](../../../../verification/v131_measure_is_area.py#L1) · explicit_script_identifier
- [verification/v132_detratio_anomaly.py](../../../../verification/v132_detratio_anomaly.py#L1) · explicit_script_identifier
- [verification/v139_selector_triangle.py](../../../../verification/v139_selector_triangle.py#L1) · explicit_script_identifier
- [verification/v124_resummed_clock.py](../../../../verification/v124_resummed_clock.py#L1) · explicit_script_identifier
- [verification/v136_dual_normal_selector.py](../../../../verification/v136_dual_normal_selector.py#L1) · explicit_script_identifier
- [verification/v110_calderon_sheet.py](../../../../verification/v110_calderon_sheet.py#L1) · explicit_script_identifier
- [verification/v125_glue_qsystem.py](../../../../verification/v125_glue_qsystem.py#L1) · explicit_script_identifier
- [verification/v123_inventory_update.py](../../../../verification/v123_inventory_update.py#L1) · explicit_script_identifier
- [verification/v135_det_surface.py](../../../../verification/v135_det_surface.py#L1) · explicit_script_identifier
- [verification/v130_born_square.py](../../../../verification/v130_born_square.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v135) determinant surface ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1027)
- [(* ---- (v133) zeta budget ---- Numerical continuations are Python-only (mpmath); the exact arithmetic of both routes is mirrored, plus a numeric spot-check in WL. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L995)
- [(* ---- (v131) measure is area ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L969)
- [(* ---- (v132) det-ratio anomaly ---- The numerical heat-trace continuation is Python-only (mpmath); the exact arithmetic is mirrored here, plus an independent exact check of the constant term via the full heat-trace expansion. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L980)
- [(* ---- (v136) dual-normal selector ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1051)
- [(* ---- (v128) graded hull ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L918)
- [(* ---- (v129) entropy power law ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L949)
- [(* ---- (v126) clock-wall bridge ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L891)
- [(* ---- (v124) resummed clock ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L857)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v139) selector triangle ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1110)
- [(* ---- (v140) canonical map ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1123)
- [(* ---- (v125) glue Q-system ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L872)
- [(* ---- (v130) Born square ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L961)
- [(* ---- (v123) inventory update: ledger-bookkeeping module (CSV parsing of status_ledger.csv, typing contract) -- Python-only by nature, like the statistical v100; no algebraic content to mirror. ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L853)
- [(* ---- (v134) dual anchor ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1011)
- [(* ---- (v110) Calderon-sheet selection ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L525)
- [(* ---- (v127) ring resummation ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L904)
- [(* ---- (v138) VW firewall ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1103)
- [(* ---- (v122) margin theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L835)

## 2026-06-11

[documents/changelog/sections/024-source.tex:76-108](../../../../documents/changelog/sections/024-source.tex#L76) · `section:9575cef792be9f73f7ff` · subsection

External-review integration (v83) and release hygiene. Manifest –check mode and the release rule; holomorphy pins [formula in source] (the unique even unimodular rank-8 lattice theory). Sheet diamond and centered form (v94/v95). All four flavor operators are points of one two-parameter surface [formula in source] ; centered cross [formula in source] , [formula in source] with the cofactor seam normal [formula in sour

Explicit source references:

- [verification/v118_hexagon_family_dictionary.py](../../../../verification/v118_hexagon_family_dictionary.py#L1) · explicit_script_identifier
- [verification/v94_sheet_diamond.py](../../../../verification/v94_sheet_diamond.py#L1) · explicit_script_identifier
- [verification/v95_centered_diamond.py](../../../../verification/v95_centered_diamond.py#L1) · explicit_script_identifier
- [verification/v109_sheet_pairing.py](../../../../verification/v109_sheet_pairing.py#L1) · explicit_script_identifier
- [verification/v101_horizon_anchor.py](../../../../verification/v101_horizon_anchor.py#L1) · explicit_script_identifier
- [verification/v113_quasifree_kernel.py](../../../../verification/v113_quasifree_kernel.py#L1) · explicit_script_identifier
- [verification/v83_e8net_holomorphic_uniqueness.py](../../../../verification/v83_e8net_holomorphic_uniqueness.py#L1) · explicit_script_identifier
- [verification/v104_nariai_clock.py](../../../../verification/v104_nariai_clock.py#L1) · explicit_script_identifier
- [verification/v108_pascal_ladder.py](../../../../verification/v108_pascal_ladder.py#L1) · explicit_script_identifier
- [verification/v103_trisection_normal_form.py](../../../../verification/v103_trisection_normal_form.py#L1) · explicit_script_identifier
- [verification/v121_address_pinning.py](../../../../verification/v121_address_pinning.py#L1) · explicit_script_identifier
- [verification/v114_torsion_delta.py](../../../../verification/v114_torsion_delta.py#L1) · explicit_script_identifier
- [verification/v119_review_validation_2.py](../../../../verification/v119_review_validation_2.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v108) Pascal ladder ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L486)
- [(* ---- (v113) quasi-free kernel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L608)
- [(* ---- (v103) trisection normal form ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L393)
- [(* ---- (v119) second review validation ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L781)
- [(* ---- (v101) horizon anchor: SdS in seam units ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L334)
- [(* ---- (v94) sheet diamond & winding line ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L180)
- [(* ---- (v104) the classical Nariai clock ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L421)
- [(* ---- (v114) torsion normal form + delta = 1/2 theorem ---- The [N] branch census (scipy Nelder-Mead sampling of the D4-fixed flat cusp locus) is Python-only by convention; the exact statements are mirrored here. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L640)
- [experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean#L1)
- [(* ---- (v118) hexagon-family dictionary ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L750)
- [(* ---- (v95) centered flavor diamond ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L205)
- [(* ---- (v121) address pinning ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L815)
- [(* ---- (v109) sheet pairing ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L500)

## 2026-06-10

[documents/changelog/sections/024-source.tex:109-118](../../../../documents/changelog/sections/024-source.tex#L109) · `section:6b03fdb0bab2125a7ce4` · subsection

B/C/D round closed. Koide attractor toy model (v93), glue uniqueness (v92), the Wolfram extension as the independent second verification path, Lean hardening; papers and website synced with the round. Content rounds v96–v100. Sheet conjugation bridge, discriminant dictionary, Koide flow time, and the look-elsewhere-corrected numerology null test (v100: joint null probability [formula in source] , conditional on the d

Explicit source references:

- [verification/v92_glue_uniqueness.py](../../../../verification/v92_glue_uniqueness.py#L1) · explicit_script_identifier
- [verification/v96_branch_kernel_selection.py](../../../../verification/v96_branch_kernel_selection.py#L1) · explicit_script_identifier
- [verification/v93_koide_relaxation_toy.py](../../../../verification/v93_koide_relaxation_toy.py#L1) · explicit_script_identifier
- [verification/v100_numerology_null_mc.py](../../../../verification/v100_numerology_null_mc.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v92) glue uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L151)
- [(* ---- (v96) branch-kernel selection (P1) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L221)
- [(* ---- (v93) Koide relaxation toy ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L166)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)

## 2026-06-09

[documents/changelog/sections/024-source.tex:119-126](../../../../documents/changelog/sections/024-source.tex#L119) · `section:c1b46441dd277d84f5eb` · subsection

Blind registry frozen (v84). Machine-enforced prediction registry ( predictions_frozen.json): exactly one [formula in source] number of record before JUNO; conditional entries carry their hypotheses by name. Research-contracts document added to the set.

Explicit source references:

- [verification/predictions_frozen.json](../../../../verification/predictions_frozen.json#L1) · explicit_file_citation
- [verification/v84_frozen_registry.py](../../../../verification/v84_frozen_registry.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v84) frozen registry decimals recompute from the axioms ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L39)

## 2026-06-07 / 2026-06-08

[documents/changelog/sections/024-source.tex:127-136](../../../../documents/changelog/sections/024-source.tex#L127) · `section:b7a05826756c1fc3ef1e` · subsection

Compiler-closure consolidation. The TFPT development is consolidated into the current document set with the verification suite (modules v1–v82: [formula in source] glue, carrier Pascal, EM fixed point, flavor matrix, cascade, bootstrap, gravity/cosmology, transport, masses, registry and gate audits) as the single proof layer; manifests ( manifest.sha256, lean_manifest.sha256) as content identity; website mirrors the

Explicit source references:

- [verification/v1_e8_glue.py](../../../../verification/v1_e8_glue.py#L1) · explicit_script_identifier
- [verification/v82_koide_attractor_splitting.py](../../../../verification/v82_koide_attractor_splitting.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)
- [(* ---- (v1) E8 glue arithmetic ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L46)

## 2026-04-27

[documents/changelog/sections/024-source.tex:137-146](../../../../documents/changelog/sections/024-source.tex#L137) · `section:38453ef5bba2b057b4e6` · subsection

Initial commit. TFPT papers and website baseline: electromagnetic closure and flavor transport (UFE bridge, birefringence seed, dyonic intercept), PDF-download tracking, accessibility and metadata passes.

## Changelog

[documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/004-section-interface-table-tfpt-paper-2-texorpdfstring-left.tex:320-323](../../../../documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/004-section-interface-table-tfpt-paper-2-texorpdfstring-left.tex#L320) · `section:3f2aaf619c0723316625` · section

No opening prose excerpt; see the source.

## 2. Prüfung des nachgereichten Forschungspapers

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/001-section-tfpt-universalraum-gemeinsame-schatten-und-ein.tex:60-103](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/001-section-tfpt-universalraum-gemeinsame-schatten-und-ein.tex#L60) · `section:8fc907113cbfc5c7a92b` · subsection

Geprüft wurden das 33-seitige PDF, sein vollständiger Markdowntext, Prüfprotokoll, Quellenmanifest und eigene Ergebnisdatei sowie das zusätzlich vorhandene Reproduktionsarchiv. Die vollständige Synthese wurde gelesen; nicht jede der 74 historischen Quellen wurde neu bewiesen. Die relevanten PDF-Seiten 16 bis 20 und 22 wurden zusätzlich bildlich geprüft.

## Historischer Teil: vollständige Hauptfassung v1.6.7

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/009-subsection-3-die-neue-pr-zise.tex:280-287](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/009-subsection-3-die-neue-pr-zise.tex#L280) · `section:bc4cfb01673779ef28a9` · section

Die folgende Fassung bleibt wortgetreu erhalten. Den aktuellen Stand liefern die vorangestellten Kapitel v1.6.8. Frühere Aussagen, Fragen und Vorhaben werden dadurch nicht rückwirkend neu bewiesen. Die neue Synthese ist separat vollständig im eingefrorenen Prüfpaket enthalten.

## 11. Reproduktion und Veröffentlichung

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/011-subsubsection-4-3-zehn-vollst-ndige-energiemomente-der-g.tex:353-376](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/011-subsubsection-4-3-zehn-vollst-ndige-energiemomente-der-g.tex#L353) · `section:83de894249e2d334406d` · subsection

Das Prüfpaket enthält sechs neue Python-Prüfer, die eingefrorenen Quellen, die alte und neu wiederholte Kontraktion, Ergebnisdateien, Quellenmanifest und die vollständige v1.6.6 als historische Basis. replay.py prüft alle Pins, beide Python-Modi und die C++-Kontraktion. Die ausgelieferte ZIP wird nach dem Entpacken erneut geprüft.

## Historischer Anhang: vollständige Fassung v1.6.6

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/013-subsubsection-b5-2-warum-das-noch-kein-holografischer.tex:284-292](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/013-subsubsection-b5-2-warum-das-noch-kein-holografischer.tex#L284) · `section:ec9e31bbdf525854ea3d` · section

Die folgende vollständige Fassung wird unverändert bewahrt. Maßgeblich für den aktuellen Status sind die vorangestellten Kapitel v1.6.7. Ältere Vorhaben, Modellannahmen und Zahlen erhalten dadurch keinen neuen Beweisstatus; die kleinen Energieverschärfungen und die neuen Singulett- und Bilinearergebnisse stehen oben.

## Historischer vollständiger Herleitungsstand v1.6.5 einschließlich v1.6.4

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/020-subsection-10-ein-korrigierter-wirklich-entscheidbarer.tex:230-238](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/020-subsection-10-ein-korrigierter-wirklich-entscheidbarer.tex#L230) · `section:b4c5f2a28f2af8f92d4e` · section

Der folgende Text wird unverändert bewahrt. Neue Gesamtstatusaussagen und Reichweitenkorrekturen stehen in der vorangestellten Revision v1.6.6. Historische Arbeitsaufträge oder Verfügbarkeitsannahmen werden dadurch nicht erneut zu aktuellen Beweisen erklärt.

## Historischer Beweisanhang: vollständiger nativer Herleitungsstand v1.6.4

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/023-subsection-1-bersicht-bernehmen-korrigieren-oder.tex:301-311](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/023-subsection-1-bersicht-bernehmen-korrigieren-oder.tex#L301) · `section:2a26f6d442de76283ca7` · section

Dieser Anhang bewahrt die frühere Herleitung vollständig. Seine damaligen Statussätze werden durch die aktuelle Konsolidierung und den Nachtrag oben ergänzt; er ist kein ungeprüft übernommener neuer Gesamtstatus. Die externe laufende Symmetrieergänzung wurde separat geprüft und ist nicht unbemerkt in diesen eingefrorenen Anhang eingegangen.

## 11. Reproduktion und Veröffentlichung

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/003-subsection-10-die-drei-n-chsten-entscheidenden.tex:34-57](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/003-subsection-10-die-drei-n-chsten-entscheidenden.tex#L34) · `section:e0eb8570f3877ef0f462` · subsection

Das Prüfpaket enthält sechs neue Python-Prüfer, die eingefrorenen Quellen, die alte und neu wiederholte Kontraktion, Ergebnisdateien, Quellenmanifest und die vollständige v1.6.6 als historische Basis. replay.py prüft alle Pins, beide Python-Modi und die C++-Kontraktion. Die ausgelieferte ZIP wird nach dem Entpacken erneut geprüft.

## Historischer Anhang: vollständige Fassung v1.6.6

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/004-subsection-b5-holografie-und-schatten-die-konkrete-nativ.tex:322-330](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/004-subsection-b5-holografie-und-schatten-die-konkrete-nativ.tex#L322) · `section:de75af879906bbb81aef` · section

Die folgende vollständige Fassung wird unverändert bewahrt. Maßgeblich für den aktuellen Status sind die vorangestellten Kapitel v1.6.7. Ältere Vorhaben, Modellannahmen und Zahlen erhalten dadurch keinen neuen Beweisstatus; die kleinen Energieverschärfungen und die neuen Singulett- und Bilinearergebnisse stehen oben.

## Historischer vollständiger Herleitungsstand v1.6.5 einschließlich v1.6.4

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/011-subsection-7-operationsgeometrie-ben-tigt-mehr-als.tex:298-306](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/011-subsection-7-operationsgeometrie-ben-tigt-mehr-als.tex#L298) · `section:0bfc1a1996ea7551fd8a` · section

Der folgende Text wird unverändert bewahrt. Neue Gesamtstatusaussagen und Reichweitenkorrekturen stehen in der vorangestellten Revision v1.6.6. Historische Arbeitsaufträge oder Verfügbarkeitsannahmen werden dadurch nicht erneut zu aktuellen Beweisen erklärt.

## Historischer Beweisanhang: vollständiger nativer Herleitungsstand v1.6.4

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/014-subsection-11-reproduktion-und-ehrliche.tex:342-352](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/014-subsection-11-reproduktion-und-ehrliche.tex#L342) · `section:af18f417563cb829f222` · section

Dieser Anhang bewahrt die frühere Herleitung vollständig. Seine damaligen Statussätze werden durch die aktuelle Konsolidierung und den Nachtrag oben ergänzt; er ist kein ungeprüft übernommener neuer Gesamtstatus. Die externe laufende Symmetrieergänzung wurde separat geprüft und ist nicht unbemerkt in diesen eingefrorenen Anhang eingegangen.

## Historischer Teil: vollständige Hauptfassung v1.6.8

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/006-subsubsection-was-auf-dem-lichtkegel.tex:266-275](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/006-subsubsection-was-auf-dem-lichtkegel.tex#L266) · `section:b0ef97d528a5180e126b` · section

Die gesamte bisherige Hauptfassung folgt wortgetreu. Aktuell maßgeblich sind die vorangestellten Kapitel v1.6.9. Sie vereinfachen den früheren Z4-Versuch durch Nb, präzisieren die Clock-Komposition und erweitern den Lorentzadapter. Alte Vorhaben oder Bedingungen gelten dadurch nicht als neu bewiesen. Die frühere vollständige Reproduktion ist als unverändertes Archiv eingefroren.

## 2. Prüfung des nachgereichten Forschungspapers

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/006-subsubsection-was-auf-dem-lichtkegel.tex:329-372](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/006-subsubsection-was-auf-dem-lichtkegel.tex#L329) · `section:0336d590fe1b12c9e284` · subsection

Geprüft wurden das 33-seitige PDF, sein vollständiger Markdowntext, Prüfprotokoll, Quellenmanifest und eigene Ergebnisdatei sowie das zusätzlich vorhandene Reproduktionsarchiv. Die vollständige Synthese wurde gelesen; nicht jede der 74 historischen Quellen wurde neu bewiesen. Die relevanten PDF-Seiten 16 bis 20 und 22 wurden zusätzlich bildlich geprüft.

## Historischer Teil: vollständige Hauptfassung v1.6.7

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/015-subsection-50-gegen-die-lteren-gepinnten-zertifikate.tex:74-81](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/015-subsection-50-gegen-die-lteren-gepinnten-zertifikate.tex#L74) · `section:db2201637a8a61196517` · section

Die folgende Fassung bleibt wortgetreu erhalten. Den aktuellen Stand liefern die vorangestellten Kapitel v1.6.8. Frühere Aussagen, Fragen und Vorhaben werden dadurch nicht rückwirkend neu bewiesen. Die neue Synthese ist separat vollständig im eingefrorenen Prüfpaket enthalten.

## 11. Reproduktion und Veröffentlichung

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/017-subsection-8-die-fundamentale-frage-was-fehlt.tex:144-167](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/017-subsection-8-die-fundamentale-frage-was-fehlt.tex#L144) · `section:bcf8a968d07f46b36fdc` · subsection

Das Prüfpaket enthält sechs neue Python-Prüfer, die eingefrorenen Quellen, die alte und neu wiederholte Kontraktion, Ergebnisdateien, Quellenmanifest und die vollständige v1.6.6 als historische Basis. replay.py prüft alle Pins, beide Python-Modi und die C++-Kontraktion. Die ausgelieferte ZIP wird nach dem Entpacken erneut geprüft.

## Historischer Anhang: vollständige Fassung v1.6.6

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/019-subsection-b7-abgleich-des-zuletzt-eingesandten-berichts.tex:86-94](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/019-subsection-b7-abgleich-des-zuletzt-eingesandten-berichts.tex#L86) · `section:c699c91fe3777a13f172` · section

Die folgende vollständige Fassung wird unverändert bewahrt. Maßgeblich für den aktuellen Status sind die vorangestellten Kapitel v1.6.7. Ältere Vorhaben, Modellannahmen und Zahlen erhalten dadurch keinen neuen Beweisstatus; die kleinen Energieverschärfungen und die neuen Singulett- und Bilinearergebnisse stehen oben.

## Historischer vollständiger Herleitungsstand v1.6.5 einschließlich v1.6.4

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/026-subsection-3-die-richtige.tex:84-92](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/026-subsection-3-die-richtige.tex#L84) · `section:c6834d91a743322cbb61` · section

Der folgende Text wird unverändert bewahrt. Neue Gesamtstatusaussagen und Reichweitenkorrekturen stehen in der vorangestellten Revision v1.6.6. Historische Arbeitsaufträge oder Verfügbarkeitsannahmen werden dadurch nicht erneut zu aktuellen Beweisen erklärt.

## Historischer Beweisanhang: vollständiger nativer Herleitungsstand v1.6.4

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/029-subsubsection-5-1-anomalien-label-anomalien.tex:139-149](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/029-subsubsection-5-1-anomalien-label-anomalien.tex#L139) · `section:c9b69df333f6d33064c0` · section

Dieser Anhang bewahrt die frühere Herleitung vollständig. Seine damaligen Statussätze werden durch die aktuelle Konsolidierung und den Nachtrag oben ergänzt; er ist kein ungeprüft übernommener neuer Gesamtstatus. Die externe laufende Symmetrieergänzung wurde separat geprüft und ist nicht unbemerkt in diesen eingefrorenen Anhang eingegangen.

## 13. Reproduktion und Lesereihenfolge

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/002-subsection-7-derselbe-grundzustand-besitzt-eine-kontroll.tex:349-368](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/002-subsection-7-derselbe-grundzustand-besitzt-eine-kontroll.tex#L349) · `section:1d121fa86acdab3e383f` · subsection

verify_spin_lift.py prüft native Tensorpins, die vollständige benötigte Casimir-Grundlage, die N3-Reduktion, Vorzeichenklassen, unabhängige CAR-Kontrollen, exakte Momente, rationale Restkonstanten sowie die Operations- und Archivpräzisierungen. replay.py wiederholt diese Rechnung normal und optimiert und führt den historischen v480-Zeugen separat aus. Der vollständige analytische Grundzustands-/Antwortbeweis steht in

Explicit source references:

- [verification/v480_multilocal_four_interval.py](../../../../verification/v480_multilocal_four_interval.py#L1) · explicit_script_identifier

## Historischer Bestand: vollständige Hauptfassung v1.6.9

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex:1-11](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex#L1) · `section:78361d4ce34a5362eac3` · section

Die vorherige Hauptfassung folgt vollständig und unverändert. Ihr offener allgemeiner Quellen-/Grundzustandsstatus wird durch v1.6.10 nur für den vorn ausdrücklich definierten Viererzyklus und Kopplungsbereich ergänzt. Aussagen über N=64, g/Delta=1/20, Lorentzfelder und größere Systeme bleiben in ihrem ursprünglichen Vertrag. Die neue N=256-Konstruktion darf nicht rückwirkend in diese alten Sätze eingesetzt werden.

## TFPT / Universalraum: Zusammensetzung aus der gemeinsamen Antwort

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex:12-17](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex#L12) · `section:f87b3bfe1382abc86cc3` · section

No opening prose excerpt; see the source.

## 1. Ergebnis: Die fehlende Verbindung wird berechenbar, aber noch nicht ursprünglich ausgewählt

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex:18-74](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex#L18) · `section:88d218fc42e29b9cc352` · subsection

Die drei beauftragten Arbeitslinien wurden parallel ausgeführt: ursprüngliche Komposition, ein kompletter kleiner Ablauf und ein früher relativistischer Feldtest. Hinzu kommt eine unabhängige globale CAR-Rechnung im Hauptstrang. Sie ergeben keinen Abschluss der TOE, aber vier konkrete Fortschritte:

## 2. Einfaches Gesamtbild

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex:75-130](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex#L75) · `section:f37af42956109f5032e6` · subsection

Man kann sich den ursprünglichen Tensor wie eine Kupplung vorstellen: Zwei Fermionen können gemeinsam in ein Boson übergehen und zurück. Diese Kupplung ist intern wesentlich besser bestimmt, als es die offenen Grundsatzfragen vermuten lassen. Die offene Frage lautet aber noch: Welche Kupplungen greifen auf dieselben Ressourcen zu, und welche sind unabhängig?

## 3. Quellenweg: ursprünglicher Clock ist keine Kopiermaschine

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex:131-200](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex#L131) · `section:024d57f8a44059142cc6` · subsection

Die alten Quellen wurden bis zu ihren tatsächlichen Trägern zurückverfolgt. Die sieben frühen Wörter wirken auf demselben C3; ihre 49 Produkte schließen dort eine siebendimensionale Algebra. Die 60 ursprünglichen Gaußwurzelstrahlen wirken auf demselben C4. Der spätere W-Konstruktor liefert einen gemeinsamen 64-Fermionträger und 60 Bosonkanäle. Keiner dieser Schritte darf allein durch neue Indexnamen zu Raumorten oder

## 4. Stärkerer Vergleich: selbst gleiche gesamte Zweiladungsspektren reichen nicht

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex:201-259](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex#L201) · `section:0bc9fb4657b46991d593` · subsection

Eine unabhängige Familie verwendet nur eine gemeinsame Bosonbank und einen zusätzlichen symmetrischen normierten Tensor C auf L Fermionkopien. Ihre globalen Paarmatrizen sind [formula in source] . In dieser Familie bleiben das globale Paargramm [formula in source] und das gesamte N=2-Hamiltonspektrum unverändert. Die lokalen Einzelkopienamplituden bleiben dagegen nicht unverändert; diese Annahme ist ausdrücklich vers

## 5. Vereinfachter vollständiger Prozess mit einer echten Aufzeichnung

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex:260-327](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex#L260) · `section:db5a452fec41f5bba106` · subsection

Der native Dreizustandsraum [formula in source] aus v1.6.8 bleibt unverändert. Neu ist die mittlere Operation:

## 6. Relativistischer Anschluss: ein positiver Adapter und ein negativer Kinetiktest

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex:328-382](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/004-section-historischer-bestand-vollst-ndige-hauptfassung.tex#L328) · `section:04b70e9419497cecdfc6` · subsection

Für ein gemischthändiges Komposit [formula in source] ist ohne zusätzliche Lorentzstruktur keine lineare Auslese in einen linken Weylspinor möglich. Mit genau einem Impulsfaktor P gibt es dagegen einen bis auf Normierung eindeutigen solchen Operator C_P. Für jedes P≠0 hat er Rang zwei und einen vierdimensionalen Kern.

## Historischer Teil: vollständige Hauptfassung v1.6.8

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/009-subsubsection-was-auf-dem-lichtkegel.tex:266-275](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/009-subsubsection-was-auf-dem-lichtkegel.tex#L266) · `section:d2f13c0dabf80d971481` · section

Die gesamte bisherige Hauptfassung folgt wortgetreu. Aktuell maßgeblich sind die vorangestellten Kapitel v1.6.9. Sie vereinfachen den früheren Z4-Versuch durch Nb, präzisieren die Clock-Komposition und erweitern den Lorentzadapter. Alte Vorhaben oder Bedingungen gelten dadurch nicht als neu bewiesen. Die frühere vollständige Reproduktion ist als unverändertes Archiv eingefroren.

## 2. Prüfung des nachgereichten Forschungspapers

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/009-subsubsection-was-auf-dem-lichtkegel.tex:329-372](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/009-subsubsection-was-auf-dem-lichtkegel.tex#L329) · `section:db270f9986c1acd1c4fe` · subsection

Geprüft wurden das 33-seitige PDF, sein vollständiger Markdowntext, Prüfprotokoll, Quellenmanifest und eigene Ergebnisdatei sowie das zusätzlich vorhandene Reproduktionsarchiv. Die vollständige Synthese wurde gelesen; nicht jede der 74 historischen Quellen wurde neu bewiesen. Die relevanten PDF-Seiten 16 bis 20 und 22 wurden zusätzlich bildlich geprüft.

## Historischer Teil: vollständige Hauptfassung v1.6.7

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/018-subsection-50-gegen-die-lteren-gepinnten-zertifikate.tex:74-81](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/018-subsection-50-gegen-die-lteren-gepinnten-zertifikate.tex#L74) · `section:c10f556505c11202e2ba` · section

Die folgende Fassung bleibt wortgetreu erhalten. Den aktuellen Stand liefern die vorangestellten Kapitel v1.6.8. Frühere Aussagen, Fragen und Vorhaben werden dadurch nicht rückwirkend neu bewiesen. Die neue Synthese ist separat vollständig im eingefrorenen Prüfpaket enthalten.

## 11. Reproduktion und Veröffentlichung

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/020-subsection-8-die-fundamentale-frage-was-fehlt.tex:144-167](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/020-subsection-8-die-fundamentale-frage-was-fehlt.tex#L144) · `section:bdfff177be89c539256f` · subsection

Das Prüfpaket enthält sechs neue Python-Prüfer, die eingefrorenen Quellen, die alte und neu wiederholte Kontraktion, Ergebnisdateien, Quellenmanifest und die vollständige v1.6.6 als historische Basis. replay.py prüft alle Pins, beide Python-Modi und die C++-Kontraktion. Die ausgelieferte ZIP wird nach dem Entpacken erneut geprüft.

## Historischer Anhang: vollständige Fassung v1.6.6

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/022-subsection-b7-abgleich-des-zuletzt-eingesandten-berichts.tex:86-94](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/022-subsection-b7-abgleich-des-zuletzt-eingesandten-berichts.tex#L86) · `section:6c6e891c0496b6a29275` · section

Die folgende vollständige Fassung wird unverändert bewahrt. Maßgeblich für den aktuellen Status sind die vorangestellten Kapitel v1.6.7. Ältere Vorhaben, Modellannahmen und Zahlen erhalten dadurch keinen neuen Beweisstatus; die kleinen Energieverschärfungen und die neuen Singulett- und Bilinearergebnisse stehen oben.

## Historischer vollständiger Herleitungsstand v1.6.5 einschließlich v1.6.4

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/029-subsection-3-die-richtige.tex:84-92](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/029-subsection-3-die-richtige.tex#L84) · `section:d2a156e75004b899c2d6` · section

Der folgende Text wird unverändert bewahrt. Neue Gesamtstatusaussagen und Reichweitenkorrekturen stehen in der vorangestellten Revision v1.6.6. Historische Arbeitsaufträge oder Verfügbarkeitsannahmen werden dadurch nicht erneut zu aktuellen Beweisen erklärt.

## Historischer Beweisanhang: vollständiger nativer Herleitungsstand v1.6.4

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/032-subsubsection-5-1-anomalien-label-anomalien.tex:139-149](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/032-subsubsection-5-1-anomalien-label-anomalien.tex#L139) · `section:eb47e1c679e877d523c3` · section

Dieser Anhang bewahrt die frühere Herleitung vollständig. Seine damaligen Statussätze werden durch die aktuelle Konsolidierung und den Nachtrag oben ergänzt; er ist kein ungeprüft übernommener neuer Gesamtstatus. Die externe laufende Symmetrieergänzung wurde separat geprüft und ist nicht unbemerkt in diesen eingefrorenen Anhang eingegangen.

## The companion documents

[documents/introduction/sections/001-section-what-this-is-about-the-jump-in-construction-prin.tex:270-302](../../../../documents/introduction/sections/001-section-what-this-is-about-the-jump-in-construction-prin.tex#L270) · `section:8ae7510979d5941b1583` · subsection

The whole development is consolidated into eight companion documents (plus this introduction), organised in layers that build on one another; this introduction is the entry point and reading guide. Each document collects the material of its layer into self-contained parts:

Explicit source references:

- [verification/v101_horizon_anchor.py](../../../../verification/v101_horizon_anchor.py#L1) · explicit_script_identifier
- [verification/v138_vw_firewall.py](../../../../verification/v138_vw_firewall.py#L1) · explicit_script_identifier
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier
- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [(* ---- (v101) horizon anchor: SdS in seam units ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L334)
- [(* ---- (v138) VW firewall ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1103)
- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)

## Honest status

[documents/origin_theory/sections/006-section-honest-status-label-sec-status.tex:1-332](../../../../documents/origin_theory/sections/006-section-honest-status-label-sec-status.tex#L1) · `section:888628823bf71f709698` · section

Source reconciliation, 5 November 2026. P1 is the structured RP seam kernel; [formula in source] is its scalar normalization. The shared four-mark/anchor reductions, analytic finite monodromy and the two finite W-bridge premise closures are consolidated in the Architecture chapter ``Shared origin'' ( REC20261105C3). The later charged MMST result retains the named microscopic local-algebra exhaustion condition for the

Explicit source references:

- [verification/v180_clock_is_mobius.py](../../../../verification/v180_clock_is_mobius.py#L1) · explicit_script_identifier
- [verification/v54_seam_horizon_keystones.py](../../../../verification/v54_seam_horizon_keystones.py#L1) · explicit_file_citation
- [verification/v175_net_existence_full_cone.py](../../../../verification/v175_net_existence_full_cone.py#L1) · explicit_script_identifier
- [verification/v56_unique_attractor.py](../../../../verification/v56_unique_attractor.py#L1) · explicit_file_citation
- [SEAM.BIT.FREEDOM.01](../../../../verification/status_ledger.csv#L572) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1034_parent_mirror.py](../../../../verification/v1034_parent_mirror.py#L1) · explicit_file_citation
- [verification/v123_inventory_update.py](../../../../verification/v123_inventory_update.py#L1) · explicit_file_citation
- [verification/v499_p2_weights_deligne_bg.py](../../../../verification/v499_p2_weights_deligne_bg.py#L1) · explicit_script_identifier
- [SEAM.CLOCK.RIGIDITY.01](../../../../verification/status_ledger.csv#L570) · explicit_claim_identifier · original status at linked ledger row
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_file_citation
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v330_qgamb_admissible_measure.py](../../../../verification/v330_qgamb_admissible_measure.py#L1) · explicit_script_identifier
- [SEAM.BIT.ORIGIN.01](../../../../verification/status_ledger.csv#L571) · explicit_claim_identifier · original status at linked ledger row
- [verification/v325_pillowcase_keystone.py](../../../../verification/v325_pillowcase_keystone.py#L1) · explicit_file_citation
- [verification/v143_graded_frobenius.py](../../../../verification/v143_graded_frobenius.py#L1) · explicit_script_identifier
- [verification/v308_seam_equiv_chain.py](../../../../verification/v308_seam_equiv_chain.py#L1) · explicit_script_identifier
- [verification/v1035_matter_coupling.py](../../../../verification/v1035_matter_coupling.py#L1) · explicit_file_citation
- [verification/v55_coxeter_cycle.py](../../../../verification/v55_coxeter_cycle.py#L1) · explicit_file_citation
- [verification/v275_qgamb_roadmap.py](../../../../verification/v275_qgamb_roadmap.py#L1) · explicit_script_identifier
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [verification/v469_seam_crossedproduct_route.py](../../../../verification/v469_seam_crossedproduct_route.py#L1) · explicit_script_identifier
- [verification/v510_seam_bit_freedom.py](../../../../verification/v510_seam_bit_freedom.py#L1) · explicit_file_citation
- [verification/v329_os_gap_reduction.py](../../../../verification/v329_os_gap_reduction.py#L1) · explicit_file_citation
- [verification/v529_seam_interacting_toy_fk.py](../../../../verification/v529_seam_interacting_toy_fk.py#L1) · explicit_script_identifier
- [verification/v210_mark_local_dtn.py](../../../../verification/v210_mark_local_dtn.py#L1) · explicit_script_identifier
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_file_citation
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v280_pillowcase_steklov.py](../../../../verification/v280_pillowcase_steklov.py#L1) · explicit_script_identifier
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_file_citation
- [verification/v201_seam_subprincipal_marks.py](../../../../verification/v201_seam_subprincipal_marks.py#L1) · explicit_script_identifier
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_file_citation
- [verification/v267_qgeo_rigidity_minimal_axiom.py](../../../../verification/v267_qgeo_rigidity_minimal_axiom.py#L1) · explicit_script_identifier
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_file_citation
- [verification/v57_horizon_crosslinks.py](../../../../verification/v57_horizon_crosslinks.py#L1) · explicit_file_citation
- [FORM.QGEO.01](../../../../verification/status_ledger.csv#L296) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1031_quantum_tensor_curvature.py](../../../../verification/v1031_quantum_tensor_curvature.py#L1) · explicit_file_citation
- [verification/v111_quadratic_transport.py](../../../../verification/v111_quadratic_transport.py#L1) · explicit_file_citation
- [verification/v176_seam_collar_realisation.py](../../../../verification/v176_seam_collar_realisation.py#L1) · explicit_script_identifier
- [verification/v181_clock_is_conformal_symmetry.py](../../../../verification/v181_clock_is_conformal_symmetry.py#L1) · explicit_script_identifier
- [verification/v302_seam_gap.py](../../../../verification/v302_seam_gap.py#L1) · explicit_script_identifier
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_file_citation
- [verification/v179_conformal_realisation.py](../../../../verification/v179_conformal_realisation.py#L1) · explicit_script_identifier
- [verification/v507_seam_bit_origin.py](../../../../verification/v507_seam_bit_origin.py#L1) · explicit_file_citation
- [verification/v113_quasifree_kernel.py](../../../../verification/v113_quasifree_kernel.py#L1) · explicit_file_citation
- [verification/v112_selfcounting_channel.py](../../../../verification/v112_selfcounting_channel.py#L1) · explicit_file_citation
- [verification/v1030_alg2_sector_frame.py](../../../../verification/v1030_alg2_sector_frame.py#L1) · explicit_file_citation
- [verification/v178_marks_kernel_attempt.py](../../../../verification/v178_marks_kernel_attempt.py#L1) · explicit_script_identifier
- [verification/v198_modular_commutator_reduction.py](../../../../verification/v198_modular_commutator_reduction.py#L1) · explicit_script_identifier
- [verification/v326_ftransfer_suite.py](../../../../verification/v326_ftransfer_suite.py#L1) · explicit_file_citation
- [verification/v195_marks_lefschetz_character.py](../../../../verification/v195_marks_lefschetz_character.py#L1) · explicit_script_identifier
- [verification/v108_pascal_ladder.py](../../../../verification/v108_pascal_ladder.py#L1) · explicit_script_identifier
- [verification/v214_seam_pillowcase.py](../../../../verification/v214_seam_pillowcase.py#L1) · explicit_script_identifier
- [verification/v303_ftransfer_dynamics.py](../../../../verification/v303_ftransfer_dynamics.py#L1) · explicit_script_identifier
- [verification/v528_seam_bit_twist_class_definition.py](../../../../verification/v528_seam_bit_twist_class_definition.py#L1) · explicit_script_identifier
- [verification/v1033_charged_disorder.py](../../../../verification/v1033_charged_disorder.py#L1) · explicit_file_citation
- [verification/v194_raw_seam_dtn_rp.py](../../../../verification/v194_raw_seam_dtn_rp.py#L1) · explicit_script_identifier
- [FORM.QGEO.BW.01](../../../../verification/status_ledger.csv#L428) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1032_covariant_curvature.py](../../../../verification/v1032_covariant_curvature.py#L1) · explicit_file_citation
- [verification/v493_celestial_wp2_clock_deformation.py](../../../../verification/v493_celestial_wp2_clock_deformation.py#L1) · explicit_script_identifier
- [verification/v124_resummed_clock.py](../../../../verification/v124_resummed_clock.py#L1) · explicit_script_identifier
- [verification/v521_seam_bit_rp_blind.py](../../../../verification/v521_seam_bit_rp_blind.py#L1) · explicit_script_identifier
- [verification/v480_multilocal_four_interval.py](../../../../verification/v480_multilocal_four_interval.py#L1) · explicit_script_identifier
- [verification/v165_premise_a_closure.py](../../../../verification/v165_premise_a_closure.py#L1) · explicit_script_identifier
- [verification/v110_calderon_sheet.py](../../../../verification/v110_calderon_sheet.py#L1) · explicit_file_citation
- [FORM.CARTAN.DET.01](../../../../verification/status_ledger.csv#L431) · explicit_claim_identifier · original status at linked ledger row
- [verification/v177_seam_marking_kernel.py](../../../../verification/v177_seam_marking_kernel.py#L1) · explicit_script_identifier
- [verification/v202_rare_kaon.py](../../../../verification/v202_rare_kaon.py#L1) · explicit_script_identifier
- [verification/v331_necessity_of_H.py](../../../../verification/v331_necessity_of_H.py#L1) · explicit_file_citation
- [SEAM.EQUIV.TWISTOR.01](../../../../verification/status_ledger.csv#L584) · explicit_claim_identifier · original status at linked ledger row
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v187_ftransfer_laws.py](../../../../verification/v187_ftransfer_laws.py#L1) · explicit_script_identifier
- [verification/v167_inventory_post_closure.py](../../../../verification/v167_inventory_post_closure.py#L1) · explicit_file_citation
- [verification/v83_e8net_holomorphic_uniqueness.py](../../../../verification/v83_e8net_holomorphic_uniqueness.py#L1) · explicit_script_identifier
- [verification/v160_seam_gaussianity_from_pf.py](../../../../verification/v160_seam_gaussianity_from_pf.py#L1) · explicit_script_identifier
- [verification/v333_cm_selection.py](../../../../verification/v333_cm_selection.py#L1) · explicit_file_citation
- [verification/v53_compiler_core.py](../../../../verification/v53_compiler_core.py#L1) · explicit_file_citation
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [TFPT.TOE.COMPLETE.01](../../../../verification/status_ledger.csv#L764) · explicit_claim_identifier · original status at linked ledger row
- [verification/v512_seam_tau_flag.py](../../../../verification/v512_seam_tau_flag.py#L1) · explicit_script_identifier
- [verification/v332_qgamb_metric_sector.py](../../../../verification/v332_qgamb_metric_sector.py#L1) · explicit_file_citation
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v133_zeta_budget.py](../../../../verification/v133_zeta_budget.py#L1) · explicit_script_identifier
- [FORM.MU4COMM.01](../../../../verification/status_ledger.csv#L434) · explicit_claim_identifier · original status at linked ledger row
- [verification/v525_seam_bit_twist_blind.py](../../../../verification/v525_seam_bit_twist_blind.py#L1) · explicit_script_identifier
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier
- [SEAM.BIT.TWISTCLASS.01](../../../../verification/status_ledger.csv#L595) · explicit_claim_identifier · original status at linked ledger row
- [FORM.SEAMEQUIV.01](../../../../verification/status_ledger.csv#L423) · explicit_claim_identifier · original status at linked ledger row
- [verification/v334_conformal_resolution.py](../../../../verification/v334_conformal_resolution.py#L1) · explicit_file_citation
- [verification/v503_qgeo_emergence_light.py](../../../../verification/v503_qgeo_emergence_light.py#L1) · explicit_file_citation
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v471_seam_horizon_replica.py](../../../../verification/v471_seam_horizon_replica.py#L1) · explicit_script_identifier
- [QGEO.ISO.01](../../../../verification/status_ledger.csv#L260) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1026_telb_round4_closure.py](../../../../verification/v1026_telb_round4_closure.py#L1) · explicit_file_citation
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v105_residual_inventory.py](../../../../verification/v105_residual_inventory.py#L1) · explicit_file_citation
- [SEAM.STRADDLE.CONE.01](../../../../verification/status_ledger.csv#L597) · explicit_claim_identifier · original status at linked ledger row
- [FORM.QGEO.02](../../../../verification/status_ledger.csv#L297) · explicit_claim_identifier · original status at linked ledger row
- [verification/v506_seam_clock_rigidity.py](../../../../verification/v506_seam_clock_rigidity.py#L1) · explicit_file_citation
- [verification/v534_seam_straddle_cone.py](../../../../verification/v534_seam_straddle_cone.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v493 round: CELEST.WP2.01 -- the clock-invariant deformation XY = Z^4 + a0 of the A3 seam orbifold (WP2 of CELEST.SEAM.01). Exact content mirrored: clock-invariance selects the 1-parameter slice sharply (P(iZ) = P(Z) forces a3 = a2 = a1 = 0; prod(Z - i^k w) = Z^4 - w^4); disc = 256 a0^3; binary-quartic invariants I = 12 a0, J = 0 identically => j = 1728 frozen (tau = i pillowcase), cross-ratio 2; negative controls a1 (j = 0) and a2-instance (j = 1556068/81); clock on H2 = Coxeter element of W(A3) (char x^3+x^2+x+1, eigenvalues {i,-1,-i}, no invariant cycle, chi_1-Fourier eigenline, period point t(i-1)(1,i,i^2)); versal weights (1,2,3,0) = regular mu4 rep; S-algebra corrections linear in a0 at the Z4 wrap (6/16 pairs) with the EH k = 2 analogue; mirrors v493. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3148)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [(* ---- (v55) cyclic element in the hull: E8 Coxeter order 30 + one-alpha cross-sector ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L246)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [(* ---- (v110) Calderon-sheet selection ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L525)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)
- [(* v1034: exact source parent envelope. The source's 3/5 is a dispersion amplitude; a standard symmetric NN bond has half of it. *)](../../../../verification/wolfram/tfpt_round7.wl#L161)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean#L1)
- [(* v325 QGEO.KEYSTONE.01: the flat pillowcase orbifold S^2(2,2,2,2) has chi_orb = 0 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2223)
- [(* ==== v512 round: SEAM.TAU.FLAG.01 -- "the flag-transitivity equivalence web": the tau = i attack on the free RP seam circle. Test family M(delta) = {+-1, +-e^(i delta)}; bare mark-transitivity is delta-blind (pair-exchanging V4), FLAG transitivity <=> delta = pi/2 <=> tau = i. Six exact mirrors: (i) the V4/D4 census over five exact members incl. flag/arc orbit counts, (ii) the solveset equivalences + family identities, (iii) the odd-split closed form + arc Laplacian + K3 indicator, (iv) the twist mode identity + AB delta-blindness + j-rationality control, (v) the H4 centraliser rigidity (Schur), and (vi) the negative controls (hexagonal/silver/Z16 census/n = 3). Mirrors v512. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5676)
- [(* ==== v325/v327: the keystone pillowcase orbifold + the cusp-weight rewrite (exact) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2221)
- [(* ---- (v56) unique attractor: gapped transport, Coxeter planes ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L257)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [(* v469 (i): the LR locality integers of both extension routes *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2930)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ---- (v112) self-counting channel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L583)
- [(* ---- (v201) QGEO.SUBPRIN.01: a mu4-mark sum is Z4-invariant => sub-principal block-diagonal ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1666)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [(* ==== v267: QGEO rigidity / minimal-axiom form of QGEO.SYM.01 (QGEO.SYM.02) ==== v262 (alpha_s RG), v263 (seesaw numpy), v264 (DtN FFT), v265 (RG+text guard), v266 (proton RG) are numerical -> Python-only; v267's exact cross-ratio/j core is mirrored here (the DtN/FFT part of v267 stays Python-only). *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2052)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [(* v1032: a separate Lorentz-covariant curvature calculation at one exact nonzero null momentum. Pi is the de Donder trace-reversed numerator; positivity is inferred only after the curvature sandwich. *)](../../../../verification/wolfram/tfpt_round7.wl#L92)
- [experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean#L1)
- [(* ==== v528 round: SEAM.BIT.TWISTCLASS.01 -- bit = twist-class choice. Exact mirrors: (i) cardinality lemma; (ii) 224-axis flip census; (iii) O-values. Pfaffian spectra Python-only. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8214)
- [(* v1031: exact rational-fibre curvature projector and positivity. *)](../../../../verification/wolfram/tfpt_round7.wl#L80)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [(* v1033: exact Z4 deck algebra and the load-bearing weight mismatch. *)](../../../../verification/wolfram/tfpt_round7.wl#L142)
- [(* ==== v503 round: QGEO.EMERGE.LIGHT.01 -- the v215 deferred emergence test (lever 7) executed light + the residual convergence QGEO-R1 == P2-R1. Exact/numeric content mirrored: (i) the CM-POINT IDENTITIES -- E2(i) = 3/Pi and E2(rho) = 2 Sqrt[3]/Pi (30-digit q-series vs the classical closed forms), hence E2*(i) = E2*(rho) = 0 EXACTLY by arithmetic; j(i) = 1728, j(rho) = 0 (KleinInvariantJ); DedekindEta[I] = Gamma[1/4]/(2 Pi^(3/4)) exact; (ii) NO DYNAMICAL SELECTION -- logF = log(tau2 |eta|^4) is S-invariant (30 digits at a generic point), the HEXAGONAL point wins (logF(rho) > logF(i)) and tau = i is a saddle (max along it, min along tau1 + i, sampled); (iii) the PILLOWCASE HALVING -- exact N = 8 discrete-torus spectrum (algebraic eigenvalues), Z2 pairing with exactly 4 self-paired cone modes, dim(even) = (N^2+4)/2 = 34, and det'_even^2 = det'_torus x prod_self EXACTLY; (iv) CLOCK => SQUARE + COXETER -- P(iZ) = P(Z) forces a3 = a2 = a1 = 0 (SolveAlways), j = 6912 I^3/(4 I^3 - J^2) = 1728 identically for Z^4 + a0, the order-4 element on H^1 has char poly lam^3 + lam^2 + lam + 1 (Coxeter of W(A3), order 4 = |mu4| = N_fam + 1), and Z6 = Aut(hexagonal) has element orders {1,2,3,6} -- no order 4. The v280-style lattice DtN curves (mod-4 off-character weights, sqrtm/expm) are numerical and stay Python-only (v503 B3-B5). Mirrors v503. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4054)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [(* ==== v499 round: P2.TYPING.01 -- the P2 weight-typing postulate hardened: the anchor a = (1,1,2) as the Birkhoff-Grothendieck splitting exponents of the Deligne canonical extension E of the flavor connection (types the v491 check-9 residual). Exact content mirrored: (i) partition enumeration under the full typing (BG dictionary h^0 = 0 <=> a_i >= 1, exhaustive window) + the naive- typing negative control (h^1(E') = 3 forces deg -6, three splittings, e2 never 5) + the h^1 discrepancy (h^1(anchor bundle) = 1 != 3); (ii) the Deligne residue-trace degree (cusp exponents {0,1/3,2/3} = the unique lift of spec(lam^3 - 1) to [0,1), trace 1 per mark, deg E = -4, par-deg 0, witness A0* charpoly); (iii) the Schur-Horn permutohedron integer points {(2,2,0),(2,1,1)} with h^0 = 0 killing {2,2,0}; (iv) the Biswas Z3-cover degree arithmetic (genus 2, eigensheaves O(-2)^2, chi cross-check, regular-rep weights, the decomposable model {0,2,2} with h^0 = 1 = the unstable companion; deck double genus 1, j = 1728); (v) the n-mark formula e2 = (n-2)(n+1)/2 = 5 only at n = 4. The sympy symbolic parts (H^1 residue basis, monodromy group closure |<U,M0>| = 24, the 324-case stability skeleton) stay Python-side (v499). Mirrors v499. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3621)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* ==== v469 round: the net-level crossed-product certification -- the LR locality integers (h_s(SO(16)_1) = 16/16 = 1 in Z; the mu4 glue h(J^k) = {1,1,1} all integer), the KLM mu-index arithmetic (4/2^2 = 1 and 16/4^2 = 1, both holomorphic) and the Kitaev 16-fold-way pin (nu = 16 x |C| = 16 = 0 mod 16, c_- = nu/2 = 8 = g_car+N_fam); mirrors v469 (the FHS Chern itself is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2924)
- [(* v1035: physical-position Fourier convention D -> i k. Retain I: c0=-s.q and c_i=i (B p)_i. The source Ward equations then propagate both affine constraints exactly. *)](../../../../verification/wolfram/tfpt_round7.wl#L182)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [(* ==== v507 round: SEAM.BIT.ORIGIN.01 -- the tautology attack on the v506 alignment bit, executed and refuted (Part A torus origin + Part B fermion dichotomy). Exact content mirrored: (i) the HALF-PERIOD DESCENT -- the fixed points of w -> -w on C/(Z + tau Z) are the 4 half-periods (2w in Lambda <=> coefficients in {0, 1/2}; the 1/3-control fails), a translation t_c descends iff 2c in Lambda (exactly the 3 nonzero half-periods), and the descended map sigma_i(x) = e_i + (e_i-e_j)(e_i-e_k)/(x-e_i) satisfies the exact curve-lift identity prod(sigma_i - e_l) = psi_i^2 prod(x - e_l) with RATIONAL psi_i, closing to Klein V4 symbolically; (ii) the CORE SOLVE -- the pair cross-ratio condition CR = -1 solves to lambda = -1 (sigma_1), 2 (sigma_2), 1/2 (sigma_3) exactly, union = the harmonic orbit with j = 1728, against j(3) = 21952/9 and j(hex) = 0; (iii) the STABILISER COUNTS -- the complete 24-triple enumeration gives (|Stab|, #involutions, #free) = (4,3,3) generic / (8,5,3) harmonic / (12,3,3) hexagonal, and in ALL three types the free involutions are EXACTLY the three sigma_c (the deck CLASS is forced); (iv) the CM FIXED POINT -- mult-by-i mod 2 fixes exactly the half-period (1,1) = c* = (1+tau)/2 while mult-by-rho 3-cycles all three, the lemniscatic clock ghat: x -> (x-1)/(x+1) needs the curve-lift factor psi = 2i/(x+1)^2 (the CM unit; t^2+1 irreducible over Q) and squares to the central -1/x; (v) the Cl(16) DICHOTOMY -- on explicit 256x256 Jordan-Wigner gammas the central implementer squares to 256 gamma_1...gamma_16 = 256 (-1)^F (nonsplit) while the EDGE implementer Utilde_e = P_3 (g1+g5)(g2+g4)(g6-g16)(g7-g15)(g8-g14) (g9-g13)(g10-g12) implements the NS edge reflection and squares to the SCALAR +2^7 (split), and in the 32-element NS dihedral lift group the central deck has EXACTLY 2 square roots (the +-quarter shifts) vs 0 for the edge deck; (vi) the SILVER BIJECTION CENSUS -- all 8 mark bijections {+-1, +-(3+2 sqrt 2)} -> {0, 1, -1, inf} carry the deck z -> -z to the EDGE sigma_2/sigma_3 (never the central sigma_1), while all 8 mu4 bijections carry it to the CENTRAL sigma_1. The v493 circularity-audit direction and the 28/21 arrangement census stay Python-side (v507). Mirrors v507. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4546)
- [(* v469 (ii): KLM mu-index + the 16-fold-way pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2937)
- [(* ==== v510 round: SEAM.BIT.FREEDOM.01 -- the freedom attack on the v506/v507 alignment bit: the POSITION half of the bit is topology. Exact content mirrored: (i) the DIHEDRAL CENSUS -- on the refined 2N-point seam circle D_16 has 17 involutions with EXACTLY ONE free (the antipodal shift by 8) and all 16 reflections have exactly 2 circle fixed points (8 site-axis + 8 edge-axis; robustness N = 8, 12; odd control N = 9: zero free), plus the quotient topology chi(C/antipode) = 0 (circle) vs chi(C/reflection) = 1 (interval, 2 boundary corners); (ii) the AUT COMPLETENESS -- Aut(C_8) = D_8 by brute force over all 8! = 40320 permutations and Aut(C_16) = D_16 by degree-2 propagation (all 32 candidates), so the census covers ALL circle symmetries; (iii) the Cl(16) DICHOTOMY CENSUS -- on explicit 256x256 Jordan-Wigner gammas the free deck's NS implementer squares to 256 gamma_1...gamma_16 = 256 (-1)^F (nonsplit) while the implementers of ALL 16 reflections exist (one of the two spin lifts +-R) and square to SCALARS {2^7 (site-axis), 2^8 (edge-axis)} (split): NONSPLIT <=> FREE over the complete census; root dichotomy 2 (the +-quarter shifts) vs 0 total over all 16 reflections; (iv) the REAL-STRUCTURE TABLE -- the deck-commuting antiholomorphic involutions solve symbolically to exactly two families (M diagonal: mu conj z, fix line through the poles; M antidiagonal: mu/conj z, fix circle |z| = sqrt mu), and the unique mark-fixing real structure per configuration is equatorial/FREE for mu4 and generic-unit, real-axis/NOT-free for silver and generic-real, and NONEXISTENT for the hexagonal frame; (v) HARMONIC + FREE => mu4 -- the deck-pair cross-ratio (1-b)^2/(1+b)^2 = -1 solves exactly to b = +-i, with the lambda table (2, 1/2, 3/4, 5/4) real and the hexagonal lambda complex; (vi) the HONEST COUNTERWITNESSES -- the discrete witness {0,1,8,9} (shift-invariant, crossing, but not quarter-shift-invariant) and the continuous witness {+-1, +-(3+4i)/5} (pair cross-ratio -1/4, j = 148176/25 != 1728, clock does not preserve): the Fock dichotomy measures the CLASS, not the modulus. The Moebius/Klein covering models and the pillowcase branch control stay Python-side (v510). Mirrors v510. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5155)
- [(* ---- (v143) graded Frobenius ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1176)
- [(* ==== v525 round: SEAM.BIT.TWISTBLIND.01 -- "the twist-state kill": the NINTH side-blind test on the v512 scoreboard (8 -> 9); the free-plus-twist class is exhausted. Three exact mirrors: (i) the integer incidence frame + the sigma o r = +-sigma congruence mechanism (swap axes: sigma o r = sigma => the twisted Gram is the diagonal unitary congruence D M D = spectrally invisible; fixing axes at pi/2: sigma o r = -sigma = the eta flip); (ii) the pi/2-existence solveset of the mark-fixing axes + the exact counterwitness census ((3+4i)/5 unimodular but NEVER a root of unity up to order 64 -- no lattice mark-compatible reflection at any dyadic refinement); (iii) the delta <-> pi - delta mirror equivariance (the lattice reflection maps the m-mark set EXACTLY to the (8-m)-mark set -- the exact mechanism behind the Casimir mirror symmetry ln|omega(D)|(m) = ln|omega(D)|(8-m) and behind pi/2 as the symmetric point of every invariant readout). The Pfaffian insertion spectra, the Bogoliubov 2-of-16 pattern census, the defect energy and the 40-digit inertia certificates are numerical (mpmath), Python-only, flagged here. Mirrors v525. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8059)
- [(* ==== v521 round: SEAM.BIT.RPBLIND.01 -- "free OS positivity does not see delta": the anti-linear census solvesets on M(delta), the counterwitness census, the V4 -> D4 closure (face #13), the cut/mark incidence, the continuum axis elimination + Cauchy-Stieltjes factorisation, the clock-tower clause and the hex/silver precondition controls. The 40-digit lattice Gram inertia spectra (N = 16 / N = 32) are Python-only by the suite convention. Mirrors v521. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L7480)
- [(* ---- (v214) QGEO.PILLOW.01: pillowcase reduction -- cross-ratio 2 => j=1728 => order-4 CM ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1727)
- [(* ---- (v54) seam=horizon keystones: triply-forced 8 + shared transfer operator ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L235)
- [(* ---- (v111) quadratic transport ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L543)
- [(* ---- (v113) quasi-free kernel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L608)
- [(* ==== v534 round: SEAM.STRADDLE.CONE.01 -- the RP cone on the interacting FK quartets (the first dynamical selection of the alignment bit). Exact mirrors: the stabilizer census of the mark family, the pi/2 partition + straddle incidence, the reflection-closure signs of the straddling quartets, and the equivariant Fix rays (uniform on m=2..6, NS-wrap-twisted on m=1,7). The 256-dim Fock spectra, the PT Gram corrections and the finite-g survival scan are Python-only (numpy), flagged per GATE.WOLFRAM.02. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8402)
- [(* ---- (v198) QGEO.MODULAR.01: principal symbol |k| commutes with the clock rho EXACTLY ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1639)
- [(* v1030: actual connected-sector words span M4; a disconnected mutant leaves one extra commutant degree. *)](../../../../verification/wolfram/tfpt_round4.wl#L163)
- [(* ---- (v57) horizon cross-links: Jacobson 8pi, Hod ln(3), Hawking |W(D5)| ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L264)
- [(* ---- (v53) compiler core: everything from (g_car,N_fam)=(5,3) ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L223)
- [(* ---- (v105) residual inventory ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L436)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v529 round: SEAM.INT.FKTOY.01 -- interacting FK toy / Kill-Test-2 shadow. Exact mirrors: (i) Z8 integer identities; (ii) straddle combinatorics; (iii) FK binomial ladder multiplicities. 256-dim spectra Python-only. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8246)
- [(* ---- (v123) inventory update: ledger-bookkeeping module (CSV parsing of status_ledger.csv, typing contract) -- Python-only by nature, like the statistical v100; no algebraic content to mirror. ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L853)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CohomologyGrading.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CohomologyGrading.lean#L1)
- [(* ---- (v195) QGEO.MARKS.02: Lefschetz/character forcing of the free mu4 orbit ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1608)
- [(* ---- (v108) Pascal ladder ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L486)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [(* ==== v506 round: SEAM.CLOCK.RIGIDITY.01 -- clock rigidity, Route B of the v503 classification (Part A Moebius + Part B Fock). Exact content mirrored: (i) the COMPLETE SQUARE-ROOT SOLVE -- g^2 = deck in PGL2(C) has exactly two solutions z -> +-iz (tr != 0 forces diagonal with (a/d)^2 = -1; tr = 0 forces det = 0), both automatically of projective order 4; (ii) the D4 CENTER + ORDER-4 COUNT + the COUNTEREXAMPLE -- the mu4 stabiliser D4 = {i^k z, i^k/z} has exactly 2 order-4 elements (the clocks), both squaring to the deck, center = {id, deck}; the silver harmonic set {+-1, +-(3+2 sqrt 2)} has j = 1728 and full D4 but its 4-cycle squares to the CROSSING involution v/z, not the deck (deck non-central: 0 roots -- harmonicity alone insufficient); (iii) the Cl(16) IDENTITY -- on explicit 256x256 Jordan-Wigner gamma matrices, Utilde = prod_j (1 - g_j g_{j+8}) satisfies Utilde^2 = 256 g_1...g_16 EXACTLY (U^2 = (-1)^F, the forced nonsplit Z4), the R implementer squares to +256 (split control), and the volume element anticommutes with every gamma; (iv) the N-ARITHMETIC -- S4 element orders {1,2,3,4} (faithfulness kills n >= 5: max cyclic image order < n for n = 5..8), free orbits force n | 4, the cusp side kills n = 6 (deg(lam^5-1) = 5 != 3 = deg(lam^3-1)), transfer spectrum {1, 64/729, 1/729} with gap (2/3)^6, and the generic 4-set stabiliser V4 = {+-z, +-2/z} has exponent 2 (not even the deck has a root). The full 24-triple stabiliser enumeration and the A4/hexagonal profile stay Python-side (v506 Part A). Mirrors v506. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4395)
- [(* ---- (v175) net existence + full-cone RP: E8 Cartan even unimodular, character order 5, Gamma(t) functoriality integers ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1523)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v133) zeta budget ---- Numerical continuations are Python-only (mpmath); the exact arithmetic of both routes is mirrored, plus a numeric spot-check in WL. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L995)
- [(* ---- (v124) resummed clock ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L857)

## Front matter / support

[documents/tfpt_1_architecture_e8/sections/004-source.tex:1-246](../../../../documents/tfpt_1_architecture_e8/sections/004-source.tex#L1) · `section:f22c82858ff073c133d4` · frontmatter

equilibrium kernel witnesses carry the full [formula in source] mixing at zero entropy production — no non-equilibrium state is needed, and whether the deployed pure-i wiring is compiler-forced became the registered successor contract ( v908_seam_equilibrium_wiring.py, SEAM.STATE.WIRING.SELECTOR.01); the \ premise stays unmoved. Dated addendum (2026-08-12, rounds 69–71): the successor contract is closed, and the answ

Explicit source references:

- [verification/v845_cfin_normal_form.py](../../../../verification/v845_cfin_normal_form.py#L1) · explicit_file_citation
- [verification/v491_p2_partition_corollary.py](../../../../verification/v491_p2_partition_corollary.py#L1) · explicit_file_citation
- [verification/v776_boundary_hamming_uniqueness.py](../../../../verification/v776_boundary_hamming_uniqueness.py#L1) · explicit_script_identifier
- [QGEO.REALIZE.01](../../../../verification/status_ledger.csv#L253) · explicit_claim_identifier · original status at linked ledger row
- [verification/v109_sheet_pairing.py](../../../../verification/v109_sheet_pairing.py#L1) · explicit_script_identifier
- [verification/v113_quasifree_kernel.py](../../../../verification/v113_quasifree_kernel.py#L1) · explicit_file_citation
- [SEAM.STRADDLE.CONE.01](../../../../verification/status_ledger.csv#L597) · explicit_claim_identifier · original status at linked ledger row
- [verification/v112_selfcounting_channel.py](../../../../verification/v112_selfcounting_channel.py#L1) · explicit_file_citation
- [verification/v499_p2_weights_deligne_bg.py](../../../../verification/v499_p2_weights_deligne_bg.py#L1) · explicit_script_identifier
- [verification/v908_seam_equilibrium_wiring.py](../../../../verification/v908_seam_equilibrium_wiring.py#L1) · explicit_file_citation
- [verification/v529_seam_interacting_toy_fk.py](../../../../verification/v529_seam_interacting_toy_fk.py#L1) · explicit_script_identifier
- [AX.P2.01](../../../../verification/status_ledger.csv#L3) · explicit_claim_identifier · original status at linked ledger row
- [verification/v525_seam_bit_twist_blind.py](../../../../verification/v525_seam_bit_twist_blind.py#L1) · explicit_script_identifier
- [verification/v110_calderon_sheet.py](../../../../verification/v110_calderon_sheet.py#L1) · explicit_file_citation
- [verification/v510_seam_bit_freedom.py](../../../../verification/v510_seam_bit_freedom.py#L1) · explicit_file_citation
- [verification/v799_seam_code_typeii.py](../../../../verification/v799_seam_code_typeii.py#L1) · explicit_script_identifier
- [verification/v506_seam_clock_rigidity.py](../../../../verification/v506_seam_clock_rigidity.py#L1) · explicit_file_citation
- [SEAM.TAU.FLAG.01](../../../../verification/status_ledger.csv#L576) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.STATE.WIRING.SELECTOR.01](../../../../verification/status_ledger.csv#L736) · explicit_claim_identifier · original status at linked ledger row
- [verification/v106_review_validation.py](../../../../verification/v106_review_validation.py#L1) · explicit_file_citation
- [SEAM.BIT.TWISTCLASS.01](../../../../verification/status_ledger.csv#L595) · explicit_claim_identifier · original status at linked ledger row
- [verification/v911_wiring_freedom.py](../../../../verification/v911_wiring_freedom.py#L1) · explicit_file_citation
- [verification/v2_carrier_pascal.py](../../../../verification/v2_carrier_pascal.py#L1) · explicit_script_identifier
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_file_citation
- [verification/v108_pascal_ladder.py](../../../../verification/v108_pascal_ladder.py#L1) · explicit_script_identifier
- [verification/v507_seam_bit_origin.py](../../../../verification/v507_seam_bit_origin.py#L1) · explicit_file_citation
- [verification/v111_quadratic_transport.py](../../../../verification/v111_quadratic_transport.py#L1) · explicit_file_citation
- [verification/v534_seam_straddle_cone.py](../../../../verification/v534_seam_straddle_cone.py#L1) · explicit_file_citation
- [verification/v521_seam_bit_rp_blind.py](../../../../verification/v521_seam_bit_rp_blind.py#L1) · explicit_script_identifier
- [verification/v512_seam_tau_flag.py](../../../../verification/v512_seam_tau_flag.py#L1) · explicit_file_citation
- [verification/v528_seam_bit_twist_class_definition.py](../../../../verification/v528_seam_bit_twist_class_definition.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v510 round: SEAM.BIT.FREEDOM.01 -- the freedom attack on the v506/v507 alignment bit: the POSITION half of the bit is topology. Exact content mirrored: (i) the DIHEDRAL CENSUS -- on the refined 2N-point seam circle D_16 has 17 involutions with EXACTLY ONE free (the antipodal shift by 8) and all 16 reflections have exactly 2 circle fixed points (8 site-axis + 8 edge-axis; robustness N = 8, 12; odd control N = 9: zero free), plus the quotient topology chi(C/antipode) = 0 (circle) vs chi(C/reflection) = 1 (interval, 2 boundary corners); (ii) the AUT COMPLETENESS -- Aut(C_8) = D_8 by brute force over all 8! = 40320 permutations and Aut(C_16) = D_16 by degree-2 propagation (all 32 candidates), so the census covers ALL circle symmetries; (iii) the Cl(16) DICHOTOMY CENSUS -- on explicit 256x256 Jordan-Wigner gammas the free deck's NS implementer squares to 256 gamma_1...gamma_16 = 256 (-1)^F (nonsplit) while the implementers of ALL 16 reflections exist (one of the two spin lifts +-R) and square to SCALARS {2^7 (site-axis), 2^8 (edge-axis)} (split): NONSPLIT <=> FREE over the complete census; root dichotomy 2 (the +-quarter shifts) vs 0 total over all 16 reflections; (iv) the REAL-STRUCTURE TABLE -- the deck-commuting antiholomorphic involutions solve symbolically to exactly two families (M diagonal: mu conj z, fix line through the poles; M antidiagonal: mu/conj z, fix circle |z| = sqrt mu), and the unique mark-fixing real structure per configuration is equatorial/FREE for mu4 and generic-unit, real-axis/NOT-free for silver and generic-real, and NONEXISTENT for the hexagonal frame; (v) HARMONIC + FREE => mu4 -- the deck-pair cross-ratio (1-b)^2/(1+b)^2 = -1 solves exactly to b = +-i, with the lambda table (2, 1/2, 3/4, 5/4) real and the hexagonal lambda complex; (vi) the HONEST COUNTERWITNESSES -- the discrete witness {0,1,8,9} (shift-invariant, crossing, but not quarter-shift-invariant) and the continuous witness {+-1, +-(3+4i)/5} (pair cross-ratio -1/4, j = 148176/25 != 1728, clock does not preserve): the Fock dichotomy measures the CLASS, not the modulus. The Moebius/Klein covering models and the pillowcase branch control stay Python-side (v510). Mirrors v510. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5155)
- [(* ---- (v109) sheet pairing ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L500)
- [(* ==== v499 round: P2.TYPING.01 -- the P2 weight-typing postulate hardened: the anchor a = (1,1,2) as the Birkhoff-Grothendieck splitting exponents of the Deligne canonical extension E of the flavor connection (types the v491 check-9 residual). Exact content mirrored: (i) partition enumeration under the full typing (BG dictionary h^0 = 0 <=> a_i >= 1, exhaustive window) + the naive- typing negative control (h^1(E') = 3 forces deg -6, three splittings, e2 never 5) + the h^1 discrepancy (h^1(anchor bundle) = 1 != 3); (ii) the Deligne residue-trace degree (cusp exponents {0,1/3,2/3} = the unique lift of spec(lam^3 - 1) to [0,1), trace 1 per mark, deg E = -4, par-deg 0, witness A0* charpoly); (iii) the Schur-Horn permutohedron integer points {(2,2,0),(2,1,1)} with h^0 = 0 killing {2,2,0}; (iv) the Biswas Z3-cover degree arithmetic (genus 2, eigensheaves O(-2)^2, chi cross-check, regular-rep weights, the decomposable model {0,2,2} with h^0 = 1 = the unstable companion; deck double genus 1, j = 1728); (v) the n-mark formula e2 = (n-2)(n+1)/2 = 5 only at n = 4. The sympy symbolic parts (H^1 residue basis, monodromy group closure |<U,M0>| = 24, the 324-case stability skeleton) stay Python-side (v499). Mirrors v499. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3621)
- [(* ==== v491 round: P2.PARTITION.01 -- g_car = 5 as the partition corollary of the four marks (the sharpening of v53). Exact content: 4 into exactly 3 positive integer parts is uniquely {1,1,2} => e = (4,5,2) => g_car = e2 = 5, |Z2| = e3 = 2, N_fam = e2 - e3 = 3 as corollaries; every assumption load-bearing (positivity / sum / part count / integrality negative controls); v53 needs all three targets, the partition only e1 = 4; sum rule 4(N-1)/2 = 4 iff N = 3; mirrors v491. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3117)
- [experiments/lean4-carrier-rigidity/TfptCarrier/ArfSpinorCompiler.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/ArfSpinorCompiler.lean#L1)
- [(* ==== v507 round: SEAM.BIT.ORIGIN.01 -- the tautology attack on the v506 alignment bit, executed and refuted (Part A torus origin + Part B fermion dichotomy). Exact content mirrored: (i) the HALF-PERIOD DESCENT -- the fixed points of w -> -w on C/(Z + tau Z) are the 4 half-periods (2w in Lambda <=> coefficients in {0, 1/2}; the 1/3-control fails), a translation t_c descends iff 2c in Lambda (exactly the 3 nonzero half-periods), and the descended map sigma_i(x) = e_i + (e_i-e_j)(e_i-e_k)/(x-e_i) satisfies the exact curve-lift identity prod(sigma_i - e_l) = psi_i^2 prod(x - e_l) with RATIONAL psi_i, closing to Klein V4 symbolically; (ii) the CORE SOLVE -- the pair cross-ratio condition CR = -1 solves to lambda = -1 (sigma_1), 2 (sigma_2), 1/2 (sigma_3) exactly, union = the harmonic orbit with j = 1728, against j(3) = 21952/9 and j(hex) = 0; (iii) the STABILISER COUNTS -- the complete 24-triple enumeration gives (|Stab|, #involutions, #free) = (4,3,3) generic / (8,5,3) harmonic / (12,3,3) hexagonal, and in ALL three types the free involutions are EXACTLY the three sigma_c (the deck CLASS is forced); (iv) the CM FIXED POINT -- mult-by-i mod 2 fixes exactly the half-period (1,1) = c* = (1+tau)/2 while mult-by-rho 3-cycles all three, the lemniscatic clock ghat: x -> (x-1)/(x+1) needs the curve-lift factor psi = 2i/(x+1)^2 (the CM unit; t^2+1 irreducible over Q) and squares to the central -1/x; (v) the Cl(16) DICHOTOMY -- on explicit 256x256 Jordan-Wigner gammas the central implementer squares to 256 gamma_1...gamma_16 = 256 (-1)^F (nonsplit) while the EDGE implementer Utilde_e = P_3 (g1+g5)(g2+g4)(g6-g16)(g7-g15)(g8-g14) (g9-g13)(g10-g12) implements the NS edge reflection and squares to the SCALAR +2^7 (split), and in the 32-element NS dihedral lift group the central deck has EXACTLY 2 square roots (the +-quarter shifts) vs 0 for the edge deck; (vi) the SILVER BIJECTION CENSUS -- all 8 mark bijections {+-1, +-(3+2 sqrt 2)} -> {0, 1, -1, inf} carry the deck z -> -z to the EDGE sigma_2/sigma_3 (never the central sigma_1), while all 8 mu4 bijections carry it to the CENTRAL sigma_1. The v493 circularity-audit direction and the 28/21 arrangement census stay Python-side (v507). Mirrors v507. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4546)
- [(* ---- (v2) carrier / Pascal ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L51)
- [(* ==== v525 round: SEAM.BIT.TWISTBLIND.01 -- "the twist-state kill": the NINTH side-blind test on the v512 scoreboard (8 -> 9); the free-plus-twist class is exhausted. Three exact mirrors: (i) the integer incidence frame + the sigma o r = +-sigma congruence mechanism (swap axes: sigma o r = sigma => the twisted Gram is the diagonal unitary congruence D M D = spectrally invisible; fixing axes at pi/2: sigma o r = -sigma = the eta flip); (ii) the pi/2-existence solveset of the mark-fixing axes + the exact counterwitness census ((3+4i)/5 unimodular but NEVER a root of unity up to order 64 -- no lattice mark-compatible reflection at any dyadic refinement); (iii) the delta <-> pi - delta mirror equivariance (the lattice reflection maps the m-mark set EXACTLY to the (8-m)-mark set -- the exact mechanism behind the Casimir mirror symmetry ln|omega(D)|(m) = ln|omega(D)|(8-m) and behind pi/2 as the symmetric point of every invariant readout). The Pfaffian insertion spectra, the Bogoliubov 2-of-16 pattern census, the defect energy and the 40-digit inertia certificates are numerical (mpmath), Python-only, flagged here. Mirrors v525. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8059)
- [(* ---- (v108) Pascal ladder ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L486)
- [(* ==== v521 round: SEAM.BIT.RPBLIND.01 -- "free OS positivity does not see delta": the anti-linear census solvesets on M(delta), the counterwitness census, the V4 -> D4 closure (face #13), the cut/mark incidence, the continuum axis elimination + Cauchy-Stieltjes factorisation, the clock-tower clause and the hex/silver precondition controls. The 40-digit lattice Gram inertia spectra (N = 16 / N = 32) are Python-only by the suite convention. Mirrors v521. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L7480)
- [(* ==== v528 round: SEAM.BIT.TWISTCLASS.01 -- bit = twist-class choice. Exact mirrors: (i) cardinality lemma; (ii) 224-axis flip census; (iii) O-values. Pfaffian spectra Python-only. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8214)
- [(* ---- (v113) quasi-free kernel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L608)
- [(* ==== v529 round: SEAM.INT.FKTOY.01 -- interacting FK toy / Kill-Test-2 shadow. Exact mirrors: (i) Z8 integer identities; (ii) straddle combinatorics; (iii) FK binomial ladder multiplicities. 256-dim spectra Python-only. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8246)
- [(* ==== v534 round: SEAM.STRADDLE.CONE.01 -- the RP cone on the interacting FK quartets (the first dynamical selection of the alignment bit). Exact mirrors: the stabilizer census of the mark family, the pi/2 partition + straddle incidence, the reflection-closure signs of the straddling quartets, and the equivariant Fix rays (uniform on m=2..6, NS-wrap-twisted on m=1,7). The 256-dim Fock spectra, the PT Gram corrections and the finite-g survival scan are Python-only (numpy), flagged per GATE.WOLFRAM.02. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8402)
- [(* ---- (v106) review validation ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L451)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean#L1)
- [(* ---- (v112) self-counting channel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L583)
- [(* ---- (v111) quadratic transport ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L543)
- [(* ---- (v110) Calderon-sheet selection ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L525)
- [experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean#L1)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [(* ==== v512 round: SEAM.TAU.FLAG.01 -- "the flag-transitivity equivalence web": the tau = i attack on the free RP seam circle. Test family M(delta) = {+-1, +-e^(i delta)}; bare mark-transitivity is delta-blind (pair-exchanging V4), FLAG transitivity <=> delta = pi/2 <=> tau = i. Six exact mirrors: (i) the V4/D4 census over five exact members incl. flag/arc orbit counts, (ii) the solveset equivalences + family identities, (iii) the odd-split closed form + arc Laplacian + K3 indicator, (iv) the twist mode identity + AB delta-blindness + j-rationality control, (v) the H4 centraliser rigidity (Schur), and (vi) the negative controls (hexagonal/silver/Z16 census/n = 3). Mirrors v512. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5676)
- [(* ==== v506 round: SEAM.CLOCK.RIGIDITY.01 -- clock rigidity, Route B of the v503 classification (Part A Moebius + Part B Fock). Exact content mirrored: (i) the COMPLETE SQUARE-ROOT SOLVE -- g^2 = deck in PGL2(C) has exactly two solutions z -> +-iz (tr != 0 forces diagonal with (a/d)^2 = -1; tr = 0 forces det = 0), both automatically of projective order 4; (ii) the D4 CENTER + ORDER-4 COUNT + the COUNTEREXAMPLE -- the mu4 stabiliser D4 = {i^k z, i^k/z} has exactly 2 order-4 elements (the clocks), both squaring to the deck, center = {id, deck}; the silver harmonic set {+-1, +-(3+2 sqrt 2)} has j = 1728 and full D4 but its 4-cycle squares to the CROSSING involution v/z, not the deck (deck non-central: 0 roots -- harmonicity alone insufficient); (iii) the Cl(16) IDENTITY -- on explicit 256x256 Jordan-Wigner gamma matrices, Utilde = prod_j (1 - g_j g_{j+8}) satisfies Utilde^2 = 256 g_1...g_16 EXACTLY (U^2 = (-1)^F, the forced nonsplit Z4), the R implementer squares to +256 (split control), and the volume element anticommutes with every gamma; (iv) the N-ARITHMETIC -- S4 element orders {1,2,3,4} (faithfulness kills n >= 5: max cyclic image order < n for n = 5..8), free orbits force n | 4, the cusp side kills n = 6 (deg(lam^5-1) = 5 != 3 = deg(lam^3-1)), transfer spectrum {1, 64/729, 1/729} with gap (2/3)^6, and the generic 4-set stabiliser V4 = {+-z, +-2/z} has exponent 2 (not even the deck has a root). The full 24-triple stabiliser enumeration and the A4/hexagonal profile stay Python-side (v506 Part A). Mirrors v506. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4395)

## The [formula in source] compiler: [formula in source] as a closure, not an input

[documents/tfpt_1_architecture_e8/sections/005-source.tex:13-114](../../../../documents/tfpt_1_architecture_e8/sections/005-source.tex#L13) · `section:eea2436d67f35ca41c00` · section

This section upgrades the [formula in source] atlas (Section lattice construction. The whole content is consolidated and verified in the standalone companion note ``The Even Carrier Code and the [formula in source] Coxeter Compiler.'' The thesis: the carrier code is the [formula in source] half-spinor, the family geometry is an [formula in source] , and [formula in source] is the standard unimodular glue closure of [

Explicit source references:

- [GNET.RAMIFIED.01](../../../../verification/status_ledger.csv#L667) · explicit_claim_identifier · original status at linked ledger row
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v722_phys_ramified_ns_r.py](../../../../verification/v722_phys_ramified_ns_r.py#L1) · explicit_file_citation
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v1_e8_glue.py](../../../../verification/v1_e8_glue.py#L1) · explicit_script_identifier
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v726_phys_car_pp_index.py](../../../../verification/v726_phys_car_pp_index.py#L1) · explicit_file_citation
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [verification/v378_seam_s3_modular.py](../../../../verification/v378_seam_s3_modular.py#L1) · explicit_script_identifier
- [verification/v779_gnet_gns_arf.py](../../../../verification/v779_gnet_gns_arf.py#L1) · explicit_file_citation
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v469_seam_crossedproduct_route.py](../../../../verification/v469_seam_crossedproduct_route.py#L1) · explicit_file_citation
- [SEAM.EQUIV.TWISTOR.01](../../../../verification/status_ledger.csv#L584) · explicit_claim_identifier · original status at linked ledger row
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [verification/v802_gnet_martingale.py](../../../../verification/v802_gnet_martingale.py#L1) · explicit_file_citation
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_script_identifier
- [verification/v92_glue_uniqueness.py](../../../../verification/v92_glue_uniqueness.py#L1) · explicit_script_identifier
- [verification/v234_seam_holomorphy_selection.py](../../../../verification/v234_seam_holomorphy_selection.py#L1) · explicit_script_identifier
- [GATE.METRIC.08](../../../../verification/status_ledger.csv#L210) · explicit_claim_identifier · original status at linked ledger row
- [verification/v489_seam_modular_commutator.py](../../../../verification/v489_seam_modular_commutator.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [(* ---- (v92) glue uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L151)
- [(* ==== v469 round: the net-level crossed-product certification -- the LR locality integers (h_s(SO(16)_1) = 16/16 = 1 in Z; the mu4 glue h(J^k) = {1,1,1} all integer), the KLM mu-index arithmetic (4/2^2 = 1 and 16/4^2 = 1, both holomorphic) and the Kitaev 16-fold-way pin (nu = 16 x |C| = 16 = 0 mod 16, c_- = nu/2 = 8 = g_car+N_fam); mirrors v469 (the FHS Chern itself is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2924)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/ProjectiveHamming.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/ProjectiveHamming.lean#L1)
- [(* ---- (v1) E8 glue arithmetic ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L46)
- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)
- [(* ==== v234: Seam-Holomorphy selection -- one condition, three faces, E8 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1930)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* v469 (ii): KLM mu-index + the 16-fold-way pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2937)
- [(* v469 (i): the LR locality integers of both extension routes *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2930)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)

## Consolidated picture

[documents/tfpt_1_architecture_e8/sections/009-source.tex:272-327](../../../../documents/tfpt_1_architecture_e8/sections/009-source.tex#L272) · `section:718f04ace2423814ada5` · section

A small grammar generates many sectors: the Standard-Model packet, flavor, the electroweak hierarchy, the cosmological constant, the metrology, and the inflation/CMB interface. The strongest new statements are exact: the decuple bridge ties [formula in source] to the tenth power of the stripped electroweak scale with [formula in source] ; the Hubble scale is the square root of [formula in source] ; and the Einstein n

Explicit source references:

- [verification/v531_coxeter_tensor_stage_a.py](../../../../verification/v531_coxeter_tensor_stage_a.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v531 round: COX.STAGEA.TENSOR.01 -- the Stage A tensor clock. Exact mirrors: order 30, chi = Phi30, Ramanujan trace tomography, honest Phi15 scope, coprimality controls. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8301)

## The muon anomalous magnetic moment — a seam vertex readout

[documents/tfpt_4_frontier/sections/002-source.tex:281-314](../../../../documents/tfpt_4_frontier/sections/002-source.tex#L281) · `section:8e607432a187ac5887cb` · section

A \ downstream readout (archive integration), not a compiler power. The carrier algebra carries a second-order topological defect beyond the one that fixes [formula in source] : with [formula in source] and the carrier compression quotient [formula in source] ( v51_boundary_half_step.py), [formula in source] . Projected through the seam-loop phase [formula in source] (the same [formula in source] unit that normalises

Explicit source references:

- [FR.MUONG2.01](../../../../verification/status_ledger.csv#L289) · explicit_claim_identifier · original status at linked ledger row
- [verification/v51_boundary_half_step.py](../../../../verification/v51_boundary_half_step.py#L1) · explicit_file_citation
- [verification/v204_muon_seam_g2.py](../../../../verification/v204_muon_seam_g2.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v204) FR.MUONG2.01: the muon seam-vertex value a_mu = 45/(524288 pi^9) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1696)
- [(* ---- (v51) delta=1/2 glue-norm origin + (g_car,N_fam)/|mu4| unification ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L206)

## Prüfung des neuen Gesamtvorschlags: Die Clock als Auslese mit Gedächtnis

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/013-subsection-die-vier-urspr-nglichen-marken-w-hlen-eine-ec.tex:250-264](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/013-subsection-die-vier-urspr-nglichen-marken-w-hlen-eine-ec.tex#L250) · `section:f30b6ecba1ba7b4c8a5f` · section

Der neue Vorschlag enthält eine tragfähige zusätzliche Verbindung. Die drei Replica-Paarungen liefern exakt normierte positive Gewichte. Zusammen mit den bereits vorhandenen nativen Reflexionen entsteht daraus eine konkrete Operation, deren Erwartungswerttransfer die gesamte historische Matrix [formula in source] ergibt. Dafür muss kein Clockeigenwert eingesetzt werden. Dieselbe Rechnung zeigt, welche Information bei

## Der ganze Quantenkanal und die gemeinsame Ereignisregel

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/015-subsection-der-uniforme-sewing-test-wird-zu-einem-exakt-.tex:74-120](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/015-subsection-der-uniforme-sewing-test-wird-zu-einem-exakt-.tex#L74) · `section:7300c871735b76c1b717` · subsection

Auf der tatsächlichen markierten Zweierebene gilt in der Pauli-Reihenfolge [formula in source] Die ursprüngliche Trine erfüllt [formula in source] . Am historischen, weiterhin bedingt ausgewählten Punkt [formula in source] ist Der partiell transponierte normierte Choioperator hat den Eigenwert [formula in source] . Dieser Kanal kann also Verschränkung mit einer Referenz erhalten. Die dreidimensionale vollständige Dep

## Verifikation, Lean und Experimente: Der tatsächliche Sitzungsstand

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/017-subsection-ein-gemeinsamer-algorithmus-statt-f-nf-neuer-.tex:35-38](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/017-subsection-ein-gemeinsamer-algorithmus-statt-f-nf-neuer-.tex#L35) · `section:0a82087fd3ce6765d02d` · section

Die historischen Ergebnisse dieses Abschnitts wurden aus den gespeicherten Originalprotokollen der Sitzung gelesen. Die Fortsetzung hat die neuen Herkunfts- und Realisierungsrechnungen gezielt frisch ausgeführt. Der gesamte tausendteilige Bestand wurde dafür nicht erneut ausgeführt; die folgenden Prüfstände bleiben voneinander unterscheidbar.

## Historischer Algorithmusvergleich

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/018-subsection-ausf-hrbare-ergebnisse-der-sitzung.tex:87-91](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/018-subsection-ausf-hrbare-ergebnisse-der-sitzung.tex#L87) · `section:93eeab96199c517de8ad` · subsection

Der im Anwendungsteil verwendete Vergleich ist am datierten Original nachgelesen: /Users/stefanhamann/Documents/Codex/2026-09-25/und-verifiziere-und-beweise-und-teste/BEURTEILUNG.md. Er umfasst den reparierten endlichen Causal-Compiler, AALpy-Vergleiche, unabhängige Modellprüfungen, SQLite-Traces und Gegenfälle. Er belegt weder eine Hylæan-Integration noch eine Überlegenheit gegenüber modernen neuronalen Weltmodellen

## 1.2 Quellenhierarchie und Stichtag

[documents/universal_room--new-2--TFPT_Universalraum_Forschungspaper_2026-09-15--TFPT_Universalraum_Forschungspaper_2026-09-15/sections/001-section-zusammenfassung-label-zusammenfassung.tex:108-167](../../../../documents/universal_room--new-2--TFPT_Universalraum_Forschungspaper_2026-09-15--TFPT_Universalraum_Forschungspaper_2026-09-15/sections/001-section-zusammenfassung-label-zusammenfassung.tex#L108) · `section:6b396e90964f8288a741` · subsection

Der Ausgangskorpus umfasst 74 eingefrorene Dokumentdateien einschließlich des nachgereichten Texts. Darunter befinden sich elf PDFs mit zusammen 830 Seiten; ihre Texte wurden extrahiert. Die zentralen und neu geänderten Abschnitte wurden gezielt gelesen, repräsentative Ausgangsseiten visuell kontrolliert. Die Zahl der inventarisierten Seiten bedeutet ausdrücklich keinen vollständigen erneuten Beweisaudit jeder histor

## 1.2 Quellenhierarchie und Stichtag

[documents/universal_room--new-3--TFPT_Universalraum_Forschungspaper_2026-09-15--TFPT_Universalraum_Forschungspaper_2026-09-15/sections/001-section-zusammenfassung-label-zusammenfassung.tex:108-167](../../../../documents/universal_room--new-3--TFPT_Universalraum_Forschungspaper_2026-09-15--TFPT_Universalraum_Forschungspaper_2026-09-15/sections/001-section-zusammenfassung-label-zusammenfassung.tex#L108) · `section:0b9455e2f6dcab42783c` · subsection

Der Ausgangskorpus umfasst 74 eingefrorene Dokumentdateien einschließlich des nachgereichten Texts. Darunter befinden sich elf PDFs mit zusammen 830 Seiten; ihre Texte wurden extrahiert. Die zentralen und neu geänderten Abschnitte wurden gezielt gelesen, repräsentative Ausgangsseiten visuell kontrolliert. Die Zahl der inventarisierten Seiten bedeutet ausdrücklich keinen vollständigen erneuten Beweisaudit jeder histor

## T2: uniform relaxed TEL-B theorem.

[experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/tfpt/tex-artefacts/toe_round4_proofs.tex:9-22](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/tfpt/tex-artefacts/toe_round4_proofs.tex#L9) · `section:83c8266dcbf844e64c6b` · paragraph

SEAM.MMST.TELB.NORM.01. For the frozen [formula in source] , [formula in source] collar, the original remainder obeys [formula in source] The interval inputs are v1022 (TV and potential), v1025 (C1/C2a), and v1026_ telb_ round4_ closure.py (CF, DG, signed CROSS, validated Fourier core and all-volume aliases). This is not a proof of the historical sharper [formula in source] or Bound-B target [formula in source] . In

Explicit source references:

- [verification/v1026_telb_round4_closure.py](../../../../verification/v1026_telb_round4_closure.py#L1) · explicit_file_citation
- [SEAM.MMST.TELB.NORM.01](../../../../verification/status_ledger.csv#L780) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1025_telb_c1_c2a_certificates.py](../../../../verification/v1025_telb_c1_c2a_certificates.py#L1) · explicit_script_identifier
- [verification/v1022_telb_round3_certificates.py](../../../../verification/v1022_telb_round3_certificates.py#L1) · explicit_script_identifier

## Validated Fourier core and the missing finite-grid aliases.

[experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/tfpt/tex-artefacts/toe_round4_proofs.tex:175-242](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/tfpt/tex-artefacts/toe_round4_proofs.tex#L175) · `section:bce0cf0522b937838636` · paragraph

Use the sigma- [formula in source] basis [formula in source] , [formula in source] , [formula in source] . The jump is [formula in source] (zero-based indices), and [formula in source] for [formula in source] . At the 512 midpoint nodes, 256-bit Acb validated eigensystems strictly separate eight occupied and eight unoccupied eigenvalues. The outward matrix DFT [formula in source] obeys the exact midpoint alias formul

Explicit source references:

- [verification/v1022_telb_round3_certificates.py](../../../../verification/v1022_telb_round3_certificates.py#L1) · explicit_script_identifier

## Research consolidation: source fields, dynamics and the remaining bridge (2026-09-09)

[experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tex-artefacts/toe_research_20260909.tex:2-15](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tex-artefacts/toe_research_20260909.tex#L2) · `section:18b3031aad48b3897a19` · section

tocsectionResearch consolidation: source fields and dynamics (2026-09-09) Status and provenance. The research archive after Round 27 is now integrated, including the Round 28–50 development chain and the subsequent named branches. These are written, model-specific arguments accompanied by executable controls, not newly formalized or independently reviewed theorems. All physical T1–T8 gates and TFPT.TOE.COMPLETE.01 re

Explicit source references:

- [experiments/theory-contracts/RESEARCH_2026-09-09.md](../../../../experiments/theory-contracts/RESEARCH_2026-09-09.md#L1) · explicit_file_citation
- [TFPT.TOE.COMPLETE.01](../../../../verification/status_ledger.csv#L764) · explicit_claim_identifier · original status at linked ledger row

## 2. Neutral multipoint limits, without changing the source

[experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tex-artefacts/toe_research_20260909.tex:97-114](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tex-artefacts/toe_research_20260909.tex#L97) · `section:65daa31841ba9073fdfb` · subsection

For the same strip and filled projector, the alternating four-word [formula in source] cancels the outer history ramps. Normal ordering is the actual finite-sea trace subtraction. Positive heat regularization, Jensen and Cauchy–Schwarz give uniform [formula in source] bounds for [formula in source] on the full four-endpoint torus, including collisions. They promote the separated-point limit to an [formula in source]

## Non-RH research through Round 27: a full-reference control theorem

[experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tex-artefacts/toe_round27_full_reference.tex:2-17](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tex-artefacts/toe_round27_full_reference.tex#L2) · `section:294c00be71c5846e45ed` · subsection

The previously local research contracts through Round 27 are now shipped with their proof sources, source-pinned checkers, negative controls and saved validation records under experiments/theory-contracts/. They extend the positive matter/constraint constructions to an interacting scalar/charge model, local charge transport, exact scalar elimination and controlled signed history sums. Additional off-shell fields, sel

Explicit source references:

- [TFPT.TOE.COMPLETE.01](../../../../verification/status_ledger.csv#L764) · explicit_claim_identifier · original status at linked ledger row

## The [formula in source] compiler: [formula in source] as a closure, not an input

[experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tfpt_1_architecture_e8.tex:1809-1910](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tfpt_1_architecture_e8.tex#L1809) · `section:cb1bdab881371b610e8f` · section

This section upgrades the [formula in source] atlas (Section lattice construction. The whole content is consolidated and verified in the standalone companion note ``The Even Carrier Code and the [formula in source] Coxeter Compiler.'' The thesis: the carrier code is the [formula in source] half-spinor, the family geometry is an [formula in source] , and [formula in source] is the standard unimodular glue closure of [

Explicit source references:

- [verification/v726_phys_car_pp_index.py](../../../../verification/v726_phys_car_pp_index.py#L1) · explicit_file_citation
- [verification/v802_gnet_martingale.py](../../../../verification/v802_gnet_martingale.py#L1) · explicit_file_citation
- [GATE.METRIC.08](../../../../verification/status_ledger.csv#L210) · explicit_claim_identifier · original status at linked ledger row
- [verification/v469_seam_crossedproduct_route.py](../../../../verification/v469_seam_crossedproduct_route.py#L1) · explicit_file_citation
- [verification/v92_glue_uniqueness.py](../../../../verification/v92_glue_uniqueness.py#L1) · explicit_script_identifier
- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [verification/v1_e8_glue.py](../../../../verification/v1_e8_glue.py#L1) · explicit_script_identifier
- [verification/v722_phys_ramified_ns_r.py](../../../../verification/v722_phys_ramified_ns_r.py#L1) · explicit_file_citation
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [verification/v234_seam_holomorphy_selection.py](../../../../verification/v234_seam_holomorphy_selection.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v779_gnet_gns_arf.py](../../../../verification/v779_gnet_gns_arf.py#L1) · explicit_file_citation
- [verification/v489_seam_modular_commutator.py](../../../../verification/v489_seam_modular_commutator.py#L1) · explicit_file_citation
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [GNET.RAMIFIED.01](../../../../verification/status_ledger.csv#L667) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.TWISTOR.01](../../../../verification/status_ledger.csv#L584) · explicit_claim_identifier · original status at linked ledger row
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_script_identifier
- [verification/v378_seam_s3_modular.py](../../../../verification/v378_seam_s3_modular.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [(* ==== v469 round: the net-level crossed-product certification -- the LR locality integers (h_s(SO(16)_1) = 16/16 = 1 in Z; the mu4 glue h(J^k) = {1,1,1} all integer), the KLM mu-index arithmetic (4/2^2 = 1 and 16/4^2 = 1, both holomorphic) and the Kitaev 16-fold-way pin (nu = 16 x |C| = 16 = 0 mod 16, c_- = nu/2 = 8 = g_car+N_fam); mirrors v469 (the FHS Chern itself is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2924)
- [(* ---- (v92) glue uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L151)
- [(* v469 (i): the LR locality integers of both extension routes *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2930)
- [experiments/lean4-carrier-rigidity/TfptCarrier/ProjectiveHamming.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/ProjectiveHamming.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)
- [(* ==== v234: Seam-Holomorphy selection -- one condition, three faces, E8 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1930)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [(* v469 (ii): KLM mu-index + the 16-fold-way pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2937)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v1) E8 glue arithmetic ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L46)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)

## Consolidated picture

[experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tfpt_1_architecture_e8.tex:4820-4875](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/transformationen/source-snapshots/tfpt/tfpt_1_architecture_e8.tex#L4820) · `section:9e65ba01694e50f9c888` · section

A small grammar generates many sectors: the Standard-Model packet, flavor, the electroweak hierarchy, the cosmological constant, the metrology, and the inflation/CMB interface. The strongest new statements are exact: the decuple bridge ties [formula in source] to the tenth power of the stripped electroweak scale with [formula in source] ; the Hubble scale is the square root of [formula in source] ; and the Einstein n

Explicit source references:

- [verification/v531_coxeter_tensor_stage_a.py](../../../../verification/v531_coxeter_tensor_stage_a.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v531 round: COX.STAGEA.TENSOR.01 -- the Stage A tensor clock. Exact mirrors: order 30, chi = Phi30, Ramanujan trace tomography, honest Phi15 scope, coprimality controls. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8301)

## Die Ausgangsannahmen enthalten mehr als zwei Zahlen

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/02_algebra.tex:3-12](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/02_algebra.tex#L3) · `section:185051c6a24906465080` · section

Die historische Kurzfassung nennt [formula in source] und einen fünfteiligen Träger. Vollständig gelesen setzt P1 einen orientierten Rand beziehungsweise eine Naht, einen geeigneten positiven Kern und eine Normierung voraus. P2 setzt Art und Struktur des Trägers voraus; die spätere [formula in source] -Lesart wird mit Farb- und schwachen Richtungen verbunden. Orientierung, Positivität, komplexe Struktur und Markierun

## Warum ein Logarithmus noch keine Sekunde ist

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/06_zeit_echo.tex:58-62](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/06_zeit_echo.tex#L58) · `section:205f9626334b2508e7c4` · section

Aus [formula in source] folgt [formula in source] . Die historischen Größen [formula in source] und [formula in source] beschreiben deshalb sechsfach wiederholte Faktoren [formula in source] und [formula in source] . Eine physikalische Rate ist erst [formula in source] für einen festgelegten Schritt [formula in source] .

## Flavor und Massenschemata

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/09_physik.tex:120-145](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/09_physik.tex#L120) · `section:4dc71a57d9f707c818de` · section

Der historische Quellmassenansatz lautet Die sektorabhängigen [formula in source] sind nicht automatisch eins. Die Restklassenmatrix mit Zeilen [formula in source] , [formula in source] , [formula in source] hat Determinante acht und charakteristisches Polynom [formula in source] . Das sind interne Daten, keine direkte Massenmessung.

## Die übrigen Zweige bleiben Teil des Bildes

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/09_physik.tex:166-181](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/09_physik.tex#L166) · `section:1ce347ba9067944a60fa` · section

L41mmL46mmL63mm Zweig& Historischer Quellwert& Offener Anteil\\ Starke CP-Phase& [formula in source] im Pfaffian-Zweig&Positives physisches Maß und nichtperturbativer Transfer.\\ Elektrische Dipolmomente&Führend null oder stark unterdrückt&Keine pauschale Null sämtlicher SM-Beiträge.\\ Axion& [formula in source] GeV; [formula in source] eV&Externe Planckskala und QCD-Massenrelation.\\ Neutrinomassen&Summe um [formula

## Zuse und Wolfram: eine Verwandtschaft ohne Rangordnung

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/10_ausblick.tex:2-15](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/10_ausblick.tex#L2) · `section:6f65165d9f307d71c65f` · section

Der eingefügte Gesprächstext vergleicht den Universalraum mit Konrad Zuse und Stephen Wolfram. Zuses Werk heißt Rechnender Raum, nicht „denkender Raum“. In seinem Buch von 1969 untersucht er digitale Beschreibungen physikalischer Felder und Teilchen und die Rolle von Automaten. Das ist ein historischer Bezugspunkt für die Vorstellung von Physik als regelgeleitetem Prozess. https://link.springer.com/book/10.1007/978-3

## Evidenz nicht durch Addition aufwerten

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/20_release_v13.tex:19-30](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/baseline-v1.3/chapters/20_release_v13.tex#L19) · `section:73540c32f1ed46362317` · section

L37mmL30mmL83mm Paket&Umfang&Prüfgrenze\\ Basisbuch v1.1&210 Bedingungen&Historische Basisprüfung; nicht erneut vollständig ausgeführt.\\ Schatten v1.1&128 Bedingungen&Historische Erweiterung; alte RH-Aktualitätsgrenze bleibt sichtbar.\\ Sechs-Quellen-Audit&100 eigene&78 exakt/Integrität, 22 numerisch; in dieser Ausgabe redaktionell integriert.\\ Gemeinsame Ausführung&231 eigene exakte&Für diese Ausgabe erneut ausgef

## Die Ausgangsannahmen enthalten mehr als zwei Zahlen

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/02_algebra.tex:3-12](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/02_algebra.tex#L3) · `section:0b66b511703173d3aaa8` · section

Die historische Kurzfassung nennt [formula in source] und einen fünfteiligen Träger. Vollständig gelesen setzt P1 einen orientierten Rand beziehungsweise eine Naht, einen geeigneten positiven Kern und eine Normierung voraus. P2 setzt Art und Struktur des Trägers voraus; die spätere [formula in source] -Lesart wird mit Farb- und schwachen Richtungen verbunden. Orientierung, Positivität, komplexe Struktur und Markierun

## Warum ein Logarithmus noch keine Sekunde ist

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/06_zeit_echo.tex:58-62](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/06_zeit_echo.tex#L58) · `section:4f2301304656920e0b55` · section

Aus [formula in source] folgt [formula in source] . Die historischen Größen [formula in source] und [formula in source] beschreiben deshalb sechsfach wiederholte Faktoren [formula in source] und [formula in source] . Eine physikalische Rate ist erst [formula in source] für einen festgelegten Schritt [formula in source] .

## Flavor und Massenschemata

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/09_physik.tex:123-148](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/09_physik.tex#L123) · `section:597060f5ce663fc8590c` · section

Der historische Quellmassenansatz lautet Die sektorabhängigen [formula in source] sind nicht automatisch eins. Die Restklassenmatrix mit Zeilen [formula in source] , [formula in source] , [formula in source] hat Determinante acht und charakteristisches Polynom [formula in source] . Das sind interne Daten, keine direkte Massenmessung.

## Die übrigen Zweige bleiben Teil des Bildes

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/09_physik.tex:169-184](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/09_physik.tex#L169) · `section:0cb56bce05f74c1c714b` · section

L41mmL46mmL63mm Zweig& Historischer Quellwert& Offener Anteil\\ Starke CP-Phase& [formula in source] im Pfaffian-Zweig&Positives physisches Maß und nichtperturbativer Transfer.\\ Elektrische Dipolmomente&Führend null oder stark unterdrückt&Keine pauschale Null sämtlicher SM-Beiträge.\\ Axion& [formula in source] GeV; [formula in source] eV&Externe Planckskala und QCD-Massenrelation.\\ Neutrinomassen&Summe um [formula

## Zuse und Wolfram: eine Verwandtschaft ohne Rangordnung

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/10_ausblick.tex:5-18](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/10_ausblick.tex#L5) · `section:423d4d8b3fa434f09d43` · section

Der eingefügte Gesprächstext vergleicht den Universalraum mit Konrad Zuse und Stephen Wolfram. Zuses Werk heißt Rechnender Raum, nicht „denkender Raum“. In seinem Buch von 1969 untersucht er digitale Beschreibungen physikalischer Felder und Teilchen und die Rolle von Automaten. Das ist ein historischer Bezugspunkt für die Vorstellung von Physik als regelgeleitetem Prozess. https://link.springer.com/book/10.1007/978-3

## Evidenz nicht durch Addition aufwerten

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/20_release_v13.tex:19-30](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/chapters/20_release_v13.tex#L19) · `section:77c15945b2d84e94ad8b` · section

L37mmL30mmL83mm Paket&Umfang&Prüfgrenze\\ Basisbuch v1.1&210 Bedingungen&Historische Basisprüfung; nicht erneut vollständig ausgeführt.\\ Schatten v1.1&128 Bedingungen&Historische Erweiterung; alte RH-Aktualitätsgrenze bleibt sichtbar.\\ Sechs-Quellen-Audit&100 eigene&78 exakt/Integrität, 22 numerisch; in dieser Ausgabe redaktionell integriert.\\ Gemeinsame Ausführung&231 eigene exakte&Für Ausgabe v1.3 erneut ausgefü

## Umfang und Versionsgrenze

[experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/main.tex:141-180](../../../../experiments/theory-contracts/universalraum-five-source-frontier-20260914/main-v1.4/main.tex#L141) · `section:ce2c74b64eb41b5c557a` · section

Die folgenden Prüfberichte dokumentieren die historische Basisfassung v1.1. Sie behaupten keinen erneuten Lauf ihrer Prüfungen in v1.3 oder v1.4. Evidenz der Version 1.3 steht in Anhang

## Die Ausgangsannahmen enthalten mehr als zwei Zahlen

[experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/02_algebra.tex:3-12](../../../../experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/02_algebra.tex#L3) · `section:ad25b2a892a5fc64ddfa` · section

Die historische Kurzfassung nennt [formula in source] und einen fünfteiligen Träger. Vollständig gelesen setzt P1 einen orientierten Rand beziehungsweise eine Naht, einen geeigneten positiven Kern und eine Normierung voraus. P2 setzt Art und Struktur des Trägers voraus; die spätere [formula in source] -Lesart wird mit Farb- und schwachen Richtungen verbunden. Orientierung, Positivität, komplexe Struktur und Markierun

## Warum ein Logarithmus noch keine Sekunde ist

[experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/06_zeit_echo.tex:58-62](../../../../experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/06_zeit_echo.tex#L58) · `section:11e557e98cb7d0b998fb` · section

Aus [formula in source] folgt [formula in source] . Die historischen Größen [formula in source] und [formula in source] beschreiben deshalb sechsfach wiederholte Faktoren [formula in source] und [formula in source] . Eine physikalische Rate ist erst [formula in source] für einen festgelegten Schritt [formula in source] .

## Flavor und Massenschemata

[experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/09_physik.tex:123-148](../../../../experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/09_physik.tex#L123) · `section:444fd13d59ba42ba8139` · section

Der historische Quellmassenansatz lautet Die sektorabhängigen [formula in source] sind nicht automatisch eins. Die Restklassenmatrix mit Zeilen [formula in source] , [formula in source] , [formula in source] hat Determinante acht und charakteristisches Polynom [formula in source] . Das sind interne Daten, keine direkte Massenmessung.

## Die übrigen Zweige bleiben Teil des Bildes

[experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/09_physik.tex:169-184](../../../../experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/09_physik.tex#L169) · `section:2b64a8efd0f8afb9064e` · section

L41mmL46mmL63mm Zweig& Historischer Quellwert& Offener Anteil\\ Starke CP-Phase& [formula in source] im Pfaffian-Zweig&Positives physisches Maß und nichtperturbativer Transfer.\\ Elektrische Dipolmomente&Führend null oder stark unterdrückt&Keine pauschale Null sämtlicher SM-Beiträge.\\ Axion& [formula in source] GeV; [formula in source] eV&Externe Planckskala und QCD-Massenrelation.\\ Neutrinomassen&Summe um [formula

## Zuse und Wolfram: eine Verwandtschaft ohne Rangordnung

[experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/10_ausblick.tex:5-18](../../../../experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/10_ausblick.tex#L5) · `section:902ddfee46d5493bbc68` · section

Der eingefügte Gesprächstext vergleicht den Universalraum mit Konrad Zuse und Stephen Wolfram. Zuses Werk heißt Rechnender Raum, nicht „denkender Raum“. In seinem Buch von 1969 untersucht er digitale Beschreibungen physikalischer Felder und Teilchen und die Rolle von Automaten. Das ist ein historischer Bezugspunkt für die Vorstellung von Physik als regelgeleitetem Prozess. https://link.springer.com/book/10.1007/978-3

## Evidenz nicht durch Addition aufwerten

[experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/20_release_v13.tex:19-30](../../../../experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/20_release_v13.tex#L19) · `section:88dcee359b10187b7775` · section

L37mmL30mmL83mm Paket&Umfang&Prüfgrenze\\ Basisbuch v1.1&210 Bedingungen&Historische Basisprüfung; nicht erneut vollständig ausgeführt.\\ Schatten v1.1&128 Bedingungen&Historische Erweiterung; alte RH-Aktualitätsgrenze bleibt sichtbar.\\ Sechs-Quellen-Audit&100 eigene&78 exakt/Integrität, 22 numerisch; in dieser Ausgabe redaktionell integriert.\\ Gemeinsame Ausführung&231 eigene exakte&Für Ausgabe v1.3 erneut ausgefü

## Eine autonome Uhr ist konstruierbar, aber sie wählt kein Programm

[experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/25_revision_v15.tex:258-268](../../../../experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/chapters/25_revision_v15.tex#L258) · `section:64d5bad2cc812343b92e` · section

Für ein vorgegebenes Programm [formula in source] setze [formula in source] , [formula in source] und [formula in source] . Der offene History-Hamiltonoperator mit Links [formula in source] und Laplacediagonalen erfüllt exakt [formula in source] Jedes Programm hat also dasselbe Spektrum [formula in source] , wiederholt mit der Datendimension. Feste Zwischenzeit-Auslesungen unterscheiden die Programme, ihre Eigenwerte

## Umfang und Versionsgrenze

[experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/main.tex:148-188](../../../../experiments/theory-contracts/universalraum-native-closure-20260914/main-v1.5/main.tex#L148) · `section:2ccb31885ee34c689030` · section

Die folgenden Prüfberichte dokumentieren die historische Basisfassung v1.1. Sie behaupten keinen erneuten Lauf ihrer Prüfungen in v1.3 oder v1.4. Evidenz der Version 1.3 steht in Anhang

## The [formula in source] compiler: [formula in source] as a closure, not an input

[experiments/theory-contracts/universalraum-native-exterior-reset-20260915/sources/tfpt_1_architecture_e8.tex:1809-1910](../../../../experiments/theory-contracts/universalraum-native-exterior-reset-20260915/sources/tfpt_1_architecture_e8.tex#L1809) · `section:d3f2b7d7ec6acb95794a` · section

This section upgrades the [formula in source] atlas (Section lattice construction. The whole content is consolidated and verified in the standalone companion note ``The Even Carrier Code and the [formula in source] Coxeter Compiler.'' The thesis: the carrier code is the [formula in source] half-spinor, the family geometry is an [formula in source] , and [formula in source] is the standard unimodular glue closure of [

Explicit source references:

- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [verification/v489_seam_modular_commutator.py](../../../../verification/v489_seam_modular_commutator.py#L1) · explicit_file_citation
- [SEAM.EQUIV.TWISTOR.01](../../../../verification/status_ledger.csv#L584) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v92_glue_uniqueness.py](../../../../verification/v92_glue_uniqueness.py#L1) · explicit_script_identifier
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [GNET.RAMIFIED.01](../../../../verification/status_ledger.csv#L667) · explicit_claim_identifier · original status at linked ledger row
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_script_identifier
- [verification/v726_phys_car_pp_index.py](../../../../verification/v726_phys_car_pp_index.py#L1) · explicit_file_citation
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [GATE.METRIC.08](../../../../verification/status_ledger.csv#L210) · explicit_claim_identifier · original status at linked ledger row
- [verification/v378_seam_s3_modular.py](../../../../verification/v378_seam_s3_modular.py#L1) · explicit_script_identifier
- [verification/v1_e8_glue.py](../../../../verification/v1_e8_glue.py#L1) · explicit_script_identifier
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_script_identifier
- [verification/v722_phys_ramified_ns_r.py](../../../../verification/v722_phys_ramified_ns_r.py#L1) · explicit_file_citation
- [verification/v234_seam_holomorphy_selection.py](../../../../verification/v234_seam_holomorphy_selection.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v802_gnet_martingale.py](../../../../verification/v802_gnet_martingale.py#L1) · explicit_file_citation
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v469_seam_crossedproduct_route.py](../../../../verification/v469_seam_crossedproduct_route.py#L1) · explicit_file_citation
- [verification/v779_gnet_gns_arf.py](../../../../verification/v779_gnet_gns_arf.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/ProjectiveHamming.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/ProjectiveHamming.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ==== v469 round: the net-level crossed-product certification -- the LR locality integers (h_s(SO(16)_1) = 16/16 = 1 in Z; the mu4 glue h(J^k) = {1,1,1} all integer), the KLM mu-index arithmetic (4/2^2 = 1 and 16/4^2 = 1, both holomorphic) and the Kitaev 16-fold-way pin (nu = 16 x |C| = 16 = 0 mod 16, c_- = nu/2 = 8 = g_car+N_fam); mirrors v469 (the FHS Chern itself is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2924)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ==== v234: Seam-Holomorphy selection -- one condition, three faces, E8 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1930)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [(* v469 (ii): KLM mu-index + the 16-fold-way pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2937)
- [(* ---- (v92) glue uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L151)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [(* ---- (v1) E8 glue arithmetic ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L46)
- [(* v469 (i): the LR locality integers of both extension routes *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2930)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)

## Consolidated picture

[experiments/theory-contracts/universalraum-native-exterior-reset-20260915/sources/tfpt_1_architecture_e8.tex:4820-4875](../../../../experiments/theory-contracts/universalraum-native-exterior-reset-20260915/sources/tfpt_1_architecture_e8.tex#L4820) · `section:72f0339c6c28a48d15a7` · section

A small grammar generates many sectors: the Standard-Model packet, flavor, the electroweak hierarchy, the cosmological constant, the metrology, and the inflation/CMB interface. The strongest new statements are exact: the decuple bridge ties [formula in source] to the tenth power of the stripped electroweak scale with [formula in source] ; the Hubble scale is the square root of [formula in source] ; and the Einstein n

Explicit source references:

- [verification/v531_coxeter_tensor_stage_a.py](../../../../verification/v531_coxeter_tensor_stage_a.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v531 round: COX.STAGEA.TENSOR.01 -- the Stage A tensor clock. Exact mirrors: order 30, chi = Phi30, Ramanujan trace tomography, honest Phi15 scope, coprimality controls. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8301)

## Die Ausgangsannahmen enthalten mehr als zwei Zahlen

[experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/02_algebra.tex:3-12](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/02_algebra.tex#L3) · `section:3bbb9db694838bd28e8d` · section

Die historische Kurzfassung nennt [formula in source] und einen fünfteiligen Träger. Vollständig gelesen setzt P1 einen orientierten Rand beziehungsweise eine Naht, einen geeigneten positiven Kern und eine Normierung voraus. P2 setzt Art und Struktur des Trägers voraus; die spätere [formula in source] -Lesart wird mit Farb- und schwachen Richtungen verbunden. Orientierung, Positivität, komplexe Struktur und Markierun

## Warum ein Logarithmus noch keine Sekunde ist

[experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/06_zeit_echo.tex:58-62](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/06_zeit_echo.tex#L58) · `section:9653c027c434993816ba` · section

Aus [formula in source] folgt [formula in source] . Die historischen Größen [formula in source] und [formula in source] beschreiben deshalb sechsfach wiederholte Faktoren [formula in source] und [formula in source] . Eine physikalische Rate ist erst [formula in source] für einen festgelegten Schritt [formula in source] .

## Flavor und Massenschemata

[experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/09_physik.tex:120-145](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/09_physik.tex#L120) · `section:74663ace7c3a41db7d1a` · section

Der historische Quellmassenansatz lautet Die sektorabhängigen [formula in source] sind nicht automatisch eins. Die Restklassenmatrix mit Zeilen [formula in source] , [formula in source] , [formula in source] hat Determinante acht und charakteristisches Polynom [formula in source] . Das sind interne Daten, keine direkte Massenmessung.

## Die übrigen Zweige bleiben Teil des Bildes

[experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/09_physik.tex:166-181](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/09_physik.tex#L166) · `section:2c0528ffec71ceb9090c` · section

L41mmL46mmL63mm Zweig& Historischer Quellwert& Offener Anteil\\ Starke CP-Phase& [formula in source] im Pfaffian-Zweig&Positives physisches Maß und nichtperturbativer Transfer.\\ Elektrische Dipolmomente&Führend null oder stark unterdrückt&Keine pauschale Null sämtlicher SM-Beiträge.\\ Axion& [formula in source] GeV; [formula in source] eV&Externe Planckskala und QCD-Massenrelation.\\ Neutrinomassen&Summe um [formula

## Zuse und Wolfram: eine Verwandtschaft ohne Rangordnung

[experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/10_ausblick.tex:2-15](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/10_ausblick.tex#L2) · `section:3922126d83272bf28332` · section

Der eingefügte Gesprächstext vergleicht den Universalraum mit Konrad Zuse und Stephen Wolfram. Zuses Werk heißt Rechnender Raum, nicht „denkender Raum“. In seinem Buch von 1969 untersucht er digitale Beschreibungen physikalischer Felder und Teilchen und die Rolle von Automaten. Das ist ein historischer Bezugspunkt für die Vorstellung von Physik als regelgeleitetem Prozess. https://link.springer.com/book/10.1007/978-3

## Evidenz nicht durch Addition aufwerten

[experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/20_release_v13.tex:19-30](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/chapters/20_release_v13.tex#L19) · `section:14e7bc55aae4b4730ed6` · section

L37mmL30mmL83mm Paket&Umfang&Prüfgrenze\\ Basisbuch v1.1&210 Bedingungen&Historische Basisprüfung; nicht erneut vollständig ausgeführt.\\ Schatten v1.1&128 Bedingungen&Historische Erweiterung; alte RH-Aktualitätsgrenze bleibt sichtbar.\\ Sechs-Quellen-Audit&100 eigene&78 exakt/Integrität, 22 numerisch; in dieser Ausgabe redaktionell integriert.\\ Gemeinsame Ausführung&231 eigene exakte&Für diese Ausgabe erneut ausgef

## Umfang und Versionsgrenze

[experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/main.tex:136-172](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/main.tex#L136) · `section:74e57e61eb59c02c2b82` · section

Die folgenden Prüfberichte dokumentieren die historische Basisfassung v1.1. Sie behaupten keinen erneuten Lauf ihrer Prüfungen in v1.3. Neue Evidenz, Quellen und Abdeckung stehen in Anhang

## 1. Was neu bestätigt oder korrigiert wurde

[experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/update.tex:75-87](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/main-v1.3/update.tex#L75) · `section:ba1676965a08414756db` · section

L39mmL61mmL50mm Bereich&Änderung&Beweisgrenze\\ History und Austausch&Geordnete innere Records zerstören [formula in source] und erzeugen [formula in source] . Die Reparatur erhält die kohärenten Austauschwege.&Korrekte lokale Ausführung, noch keine native Auswahl.\\ Zweizellenmodell&Eindeutiger Grundzustand und Gap [formula in source] für alle endlichen [formula in source] .&Zwei vollständige Tetramer mit einer Brüc

## Die Ausgangsannahmen enthalten mehr als zwei Zahlen

[experiments/theory-contracts/universalraum-paired-release-20260914/source-v1.1/TFPT_Universalraum_LaTeX_2026-09-14/chapters/02_algebra.tex:3-12](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/source-v1.1/TFPT_Universalraum_LaTeX_2026-09-14/chapters/02_algebra.tex#L3) · `section:1a2d5258f1b335d55b2b` · section

Die historische Kurzfassung nennt [formula in source] und einen fünfteiligen Träger. Vollständig gelesen setzt P1 einen orientierten Rand beziehungsweise eine Naht, einen geeigneten positiven Kern und eine Normierung voraus. P2 setzt Art und Struktur des Trägers voraus; die spätere [formula in source] -Lesart wird mit Farb- und schwachen Richtungen verbunden. Orientierung, Positivität, komplexe Struktur und Markierun

## Warum ein Logarithmus noch keine Sekunde ist

[experiments/theory-contracts/universalraum-paired-release-20260914/source-v1.1/TFPT_Universalraum_LaTeX_2026-09-14/chapters/06_zeit_echo.tex:58-62](../../../../experiments/theory-contracts/universalraum-paired-release-20260914/source-v1.1/TFPT_Universalraum_LaTeX_2026-09-14/chapters/06_zeit_echo.tex#L58) · `section:0cac91ff575e48441883` · section

Aus [formula in source] folgt [formula in source] . Die historischen Größen [formula in source] und [formula in source] beschreiben deshalb sechsfach wiederholte Faktoren [formula in source] und [formula in source] . Eine physikalische Rate ist erst [formula in source] für einen festgelegten Schritt [formula in source] .
