# Quantum fields, dynamics and continuum — part 6

[Topic entry](../quantum-continuum.md) · [Registry entry](../../INDEX.md)

Exact ledger statuses remain at the linked original rows and in registry.json. Citations are source links, not proof certification.

## Chain A — trace identity.

[documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/003-subsection-layer-6b-yukawa-stage-b-algebraic-content.tex:256-281](../../../../documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/003-subsection-layer-6b-yukawa-stage-b-algebraic-content.tex#L256) · `section:35a557d0d4bb41900f4f` · paragraph

[formula in source] together with the rank discharges [formula in source] and the determinant-normalised coefficient choice [formula in source] give, via TraceProjection, [formula in source] The Yukawa input is now an alternating trilinear form with two explicit Prop predicates. The contraction iso and the rank-3 conclusion are both derived theorems.

## Yukawa: Stage D resolves the v5 critique.

[documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/004-section-interface-table-tfpt-paper-2-texorpdfstring-left.tex:114-134](../../../../documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/004-section-interface-table-tfpt-paper-2-texorpdfstring-left.tex#L114) · `section:5b5deee14edc9e9ff37f` · paragraph

The v5 critique of Stage C was that it bundled the contraction iso as data, so it was equivalent (in finite dimensions) to [formula in source] and did not "explain" the rank. v6 resolved this by introducing Stage D; v7 renamed its operational predicates to keep semantic-load words out of the Lean interface:

Explicit source references:

- [verification/v7_gravity_cosmo.py](../../../../verification/v7_gravity_cosmo.py#L1) · explicit_script_identifier
- [verification/v5_e8_cascade.py](../../../../verification/v5_e8_cascade.py#L1) · explicit_script_identifier
- [verification/v6_bootstrap.py](../../../../verification/v6_bootstrap.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/YukawaTrilinearForm.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/YukawaTrilinearForm.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CalderonInterface.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CalderonInterface.lean#L1)

## Why Nontrivial\, [formula in source] is necessary.

[documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/004-section-interface-table-tfpt-paper-2-texorpdfstring-left.tex:135-146](../../../../documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/004-section-interface-table-tfpt-paper-2-texorpdfstring-left.tex#L135) · `section:c0b758d352b686c50e93` · paragraph

Without the non-triviality hypothesis the rank-3 slogan would be false: for [formula in source] the spaces [formula in source] and [formula in source] are both zero, the contraction [formula in source] is trivially bijective, both operational predicates hold vacuously, but [formula in source] . The transparency theorem finrank_E_eq_zero_or_three (v9) records the unconditional disjunction [formula in source] ; the ran

Explicit source references:

- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1) · explicit_file_citation
- [verification/v9_neutrino_texture.py](../../../../verification/v9_neutrino_texture.py#L1) · explicit_script_identifier
- [TFPT.Carrier.YukawaTrilinearForm.YukawaTrilinearForm.finrank_E_eq_zero_or_three](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/YukawaTrilinearForm.lean#L221) · explicit_lean_symbol

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean#L1)

## v11 — review-pass closure and typed targets for the remaining upstream theorems (the present revision)

[documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/005-subsection-v11-review-pass-closure-and-typed-targets-for.tex:1-94](../../../../documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/005-subsection-v11-review-pass-closure-and-typed-targets-for.tex#L1) · `section:307d97065ae218f24093` · subsection

Relative to v10, this version is a closure-and-naming pass. It addresses the v10 reviewer's three actionable findings:

Explicit source references:

- [TFPT.Carrier.BoundaryYukawaKernelInterface.BoundaryYukawaKernel.finrank_E_eq_three](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryYukawaKernelInterface.lean#L99) · explicit_lean_symbol
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamWindingInterface.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamWindingInterface.lean#L1) · explicit_file_citation
- [verification/v11_unique_KQ.py](../../../../verification/v11_unique_KQ.py#L1) · explicit_script_identifier
- [verification/v9_neutrino_texture.py](../../../../verification/v9_neutrino_texture.py#L1) · explicit_script_identifier
- [verification/v10_projection_involution.py](../../../../verification/v10_projection_involution.py#L1) · explicit_script_identifier
- [experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryYukawaKernelInterface.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryYukawaKernelInterface.lean#L1) · explicit_file_citation
- [TFPT.Carrier.BoundaryPolarization.BoundaryPolarization.ofCarrierInvolution](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryPolarization.lean#L184) · explicit_lean_symbol

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BoundaryYukawaKernelInterfaceCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BoundaryYukawaKernelInterfaceCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamWindingInterface.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamWindingInterface.lean#L1)
- [(* ---- (v10) K/Q/Sigma determinant ladder ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L62)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryPolarization.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryPolarization.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/HiggsSchemeCohomologyShadow.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/HiggsSchemeCohomologyShadow.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryYukawaKernelInterface.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryYukawaKernelInterface.lean#L1)

## Next Steps

[documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/006-subsection-v6-yukawa-stage-d-and-calder-n-interface.tex:219-283](../../../../documents/experiments--lean4-carrier-rigidity--note_carrier_rigidity_lean4/sections/006-subsection-v6-yukawa-stage-d-and-calder-n-interface.tex#L219) · `section:be1ea56966a7d77ab977` · section

CI workflow. Wire audit.sh into a GitHub Actions workflow that refuses any merge with AUDIT:\ FAIL.

Explicit source references:

- [verification/v11_unique_KQ.py](../../../../verification/v11_unique_KQ.py#L1) · explicit_script_identifier
- [verification/v10_projection_involution.py](../../../../verification/v10_projection_involution.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BoundaryYukawaKernelInterfaceCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BoundaryYukawaKernelInterfaceCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryPolarization.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryPolarization.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/HiggsSchemeCohomologyShadow.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/HiggsSchemeCohomologyShadow.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryYukawaKernelInterface.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BoundaryYukawaKernelInterface.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamWindingInterface.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamWindingInterface.lean#L1)
- [(* ---- (v10) K/Q/Sigma determinant ladder ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L62)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean#L1)

## 2.2 Welche Gleichheit ist tatsächlich beobachtbar?

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/001-section-tfpt-universalraum-gemeinsame-schatten-und-ein.tex:135-162](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/001-section-tfpt-universalraum-gemeinsame-schatten-und-ein.tex#L135) · `section:a28d5bf4b965ffd38a72` · subsubsection

Ein markierter komplexer Wortkern ist stärker als eine Liste gewöhnlicher Quantenkanäle. Beispielsweise liefern [formula in source] und [formula in source] für jeden Zustand denselben Kanal:

## Ergebnis und Leseschlüssel

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/009-subsection-3-die-neue-pr-zise.tex:294-340](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/009-subsection-3-die-neue-pr-zise.tex#L294) · `section:0a809facfc4dd471f5b2` · subsubsection

Die zuletzt ergänzte Vogelperspektive steht im Teil B dieser Fassung: Gesucht wird der kleinste gemeinsame beobachtbare Prozess hinter den verschiedenen Ansichten. Ein nativer Rang-60-Kodierer ist bereits vorhanden. Seine Informationsverluste, die minimalen symmetrieverträglichen Dynamiken und der Unterschied zwischen Korrelation und kausalem Eingriff wurden neu berechnet. Das verschiebt die fundamentale Priorität vo

## B2. Noethers Frage: Welche Unterschiede kann die Dynamik überhaupt verändern?

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/012-subsubsection-schluss-in-einem-satz-label-schluss-in-ein.tex:68-120](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/012-subsubsection-schluss-in-einem-satz-label-schluss-in-ein.tex#L68) · `section:241cc65d5f3d1f2d1c6c` · subsection

Bei ungebrochener innerer Gruppe G zerfällt ein Zustandsraum als

## B3. Korrelation ist noch kein Eingriff: ein exakt gelöstes Gegenbeispiel

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/012-subsubsection-schluss-in-einem-satz-label-schluss-in-ein.tex:121-173](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/012-subsubsection-schluss-in-einem-satz-label-schluss-in-ein.tex#L121) · `section:cac38694d999ad339b62` · subsection

Die Kreuzkorrelation startet bei null und verändert sich mit der Zeit. Es gibt dennoch keinerlei Signal von A nach B: Alle A-lokalen Operationen kommutieren mit jedem zeitentwickelten B-lokalen Observablen. Für jedes spurtreue lokale Quanteninstrument bleibt daher die unbedingte B-Statistik unverändert. Eine Konditionierung auf ein Messergebnis in A wäre etwas anderes.

## B4.2 Eine weitere Vereinfachung: neue Dynamik oder andere Teilchenvariablen?

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/012-subsubsection-schluss-in-einem-satz-label-schluss-in-ein.tex:282-330](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/012-subsubsection-schluss-in-einem-satz-label-schluss-in-ein.tex#L282) · `section:b8942f44eecb86a5a53f` · subsubsection

Manche scheinbaren Erweiterungen beschreiben nur andere Bosonvariablen. Auf der reellen, phasengleichen Teilfamilie setze man

## 1. Welche Quellen zusammengeführt werden

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/020-subsection-10-ein-korrigierter-wirklich-entscheidbarer.tex:309-354](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/020-subsection-10-ein-korrigierter-wirklich-entscheidbarer.tex#L309) · `section:4c9884700b1e237edec5` · subsection

Die Benutzeranlage addbf22d-4764-4685-bdfe-f422bfc563d0/pasted-text.txt, eingefroren unter sources/user_attachment.txt. Die externe Untersuchung Universalraum_Beweisversuch_2026-09-15.md: Spektralzählung, Primzahl-Dynamik, Phasencocycle und kontrollierte Polnäherung. Die externe Untersuchung Universalraum_Urspruenglicher_Austausch_2026-09-15.md: interne Paarumwandlung, zusätzlicher Rotor-Link und älteres Round37-Beis

Explicit source references:

- [experiments/theory-contracts/universalraum-two-bank-transfer-20260915/sources/user_attachment.txt](../../../../experiments/theory-contracts/universalraum-two-bank-transfer-20260915/sources/user_attachment.txt#L1) · explicit_file_citation
- [experiments/theory-contracts/universalraum-two-bank-transfer-20260915/LATE_AUDIT.md](../../../../experiments/theory-contracts/universalraum-two-bank-transfer-20260915/LATE_AUDIT.md#L1) · explicit_file_citation
- [experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/user-attachments/39228cf8-f71d-43e6-9a05-4d49d38d721b/pasted-text.txt](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/user-attachments/39228cf8-f71d-43e6-9a05-4d49d38d721b/pasted-text.txt#L1) · explicit_file_citation

## 6.2 Primzahlen können eine gewählte Dynamik organisieren, erzwingen sie aber nicht

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/022-subsection-5-warum-die-herkunftsl-cke-nicht-durch-l-nger.tex:71-101](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/022-subsection-5-warum-die-herkunftsl-cke-nicht-durch-l-nger.tex#L71) · `section:e3b031accbaec6f14145` · subsubsection

Auf [formula in source] kann man [formula in source] definieren. Die additive Primfaktorbesetzung erklärt dann die Eulerstruktur. Der Definitionsbereich lautet [formula in source] . Dass quantenstatistische Systeme mit Zeta-Zustandssumme existieren, ist bereits ein Ergebnis von https://alainconnes.org/wp-content/uploads/bostconnesscan.pdfBost und Connes. Das ist ein relevanter Anschluss, nicht die fehlende native Her

## 8.1 Konkrete Präparationsgrenze

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/025-subsection-8-was-der-operationssatz-jetzt-tats-chlich.tex:58-73](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/025-subsection-8-was-der-operationssatz-jetzt-tats-chlich.tex#L58) · `section:f334ddf0d171d0767f17` · subsubsection

Alle genannten Fockkontrollen erhalten N. Ein Wort aus ihnen und jeder einzeln zahlenerhaltende ausgewählte Messzweig bleiben im Ausgangssektor. Sie können aus dem leeren N=0-Zustand nicht den N=64-Grundzustand erzeugen. Diese Aussage betrifft zahlenerhaltende Krauszweige; eine bloß U(1)-kovariante offene Dynamik kann sehr wohl geladene Krausoperatoren haben.

## Das Verfügbarkeits-Theorem (Parallelsession, symmetry_availability.py)Das Verfügbarkeits-Theorem (Parallelsession, )

[documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/027-subsubsection-die-isotypie-zerlegung-von-n-3-beide-spure.tex:31-52](../../../../documents/experiments--theory-contracts--universalraum-common-process-20260915--tmp--pdfs--TFPT_Universalraum_Gemeinsamer_Prozess_Konsolidierung_2026-09-15_v1.6.8/sections/027-subsubsection-die-isotypie-zerlegung-von-n-3-beide-spure.tex#L31) · `section:79e6c1b533d8873f2f4a` · subsubsection

Die gewährte Dynamik (Wörter aus X, N_b) liegt in der Verflechteralgebra (Dim Σmᵢ² = 16); die 60 Lie-Erzeuger liegen nicht darin - keine Komposition gewährter Dynamik erzeugt sie. Die Symmetrie ist ein eigenständiges Primitiv: gewährt oder nicht gewährt. Der gesamte verbleibende symmetrieverträgliche Spielraum ist exakt 2 Dimensionen breit = die Fähigkeit, die drei dunklen Blöcke getrennt anzusprechen (= beide Casimi

## Ausgangslage

[documents/experiments--theory-contracts--universalraum-fable-kernfragen-20260910--fable-runde3--paper--tfpt_anschluss_zellen_seam_2026-09-14/sections/001-section-evidenzkonvention-und-abgrenzung.tex:20-26](../../../../documents/experiments--theory-contracts--universalraum-fable-kernfragen-20260910--fable-runde3--paper--tfpt_anschluss_zellen_seam_2026-09-14/sections/001-section-evidenzkonvention-und-abgrenzung.tex#L20) · `section:5926380330792d2bbe2b` · section

Das Manuskript (Version 1.1) konstruiert auf den 60 Quellstrahlen des Compilers eine endliche Kette mit klassischem und quantenmechanischem Schatten, eine Vor-Messung mit Registern und ein ausdrücklich zusätzlich gesetztes Vierträger-Austauschmodell [formula in source] mit eindeutigem antisymmetrischem Grundzustand [formula in source] , Gap [formula in source] , dem Spektrum [formula in source] , [formula in source]

## Kopplungsvereinigung mit dem quelleigenen Pati–Salam-Inhalt

[documents/experiments--theory-contracts--universalraum-fable-kernfragen-20260910--fable-runde3--paper--tfpt_anschluss_zellen_seam_2026-09-14/sections/001-section-evidenzkonvention-und-abgrenzung.tex:129-143](../../../../documents/experiments--theory-contracts--universalraum-fable-kernfragen-20260910--fable-runde3--paper--tfpt_anschluss_zellen_seam_2026-09-14/sections/001-section-evidenzkonvention-und-abgrenzung.tex#L129) · `section:ba9522cc4af6db9d84e8` · section

QFT4D.RGTEST.01 trägt [X] für den reinen SM-Lauf. Mit dem Inhalt aus PS.ALGEBRA.01 – Eichgruppe [formula in source] , Skalare nur aus [formula in source] , kein 126 – ist die Zweistufe [formula in source] bei 1-Loop ein lineares System und reproduziert v249 ziffergenau. Bei 2-Loop (gauge-only, Machacek–Vaughn; validiert an [formula in source] , [formula in source] , [formula in source] ): @lcccccc@ Inhalt oberhalb [f

Explicit source references:

- [verification/v249_ps_unification.py](../../../../verification/v249_ps_unification.py#L1) · explicit_script_identifier
- [PS.ALGEBRA.01](../../../../verification/status_ledger.csv#L334) · explicit_claim_identifier · original status at linked ledger row
- [QFT4D.RGTEST.01](../../../../verification/status_ledger.csv#L332) · explicit_claim_identifier · original status at linked ledger row

## Ergebnis und Leseschlüssel

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/001-section-tfpt-universalraum-die-richtige-einfache.tex:13-59](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/001-section-tfpt-universalraum-die-richtige-einfache.tex#L13) · `section:47b40fbb6906a60b7264` · subsubsection

Die zuletzt ergänzte Vogelperspektive steht im Teil B dieser Fassung: Gesucht wird der kleinste gemeinsame beobachtbare Prozess hinter den verschiedenen Ansichten. Ein nativer Rang-60-Kodierer ist bereits vorhanden. Seine Informationsverluste, die minimalen symmetrieverträglichen Dynamiken und der Unterschied zwischen Korrelation und kausalem Eingriff wurden neu berechnet. Das verschiebt die fundamentale Priorität vo

## B2. Noethers Frage: Welche Unterschiede kann die Dynamik überhaupt verändern?

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/003-subsection-10-die-drei-n-chsten-entscheidenden.tex:125-177](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/003-subsection-10-die-drei-n-chsten-entscheidenden.tex#L125) · `section:bb57f2d79e3d7a50a419` · subsection

Bei ungebrochener innerer Gruppe G zerfällt ein Zustandsraum als

## B3. Korrelation ist noch kein Eingriff: ein exakt gelöstes Gegenbeispiel

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/003-subsection-10-die-drei-n-chsten-entscheidenden.tex:178-230](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/003-subsection-10-die-drei-n-chsten-entscheidenden.tex#L178) · `section:cea93c1ca326f1224317` · subsection

Die Kreuzkorrelation startet bei null und verändert sich mit der Zeit. Es gibt dennoch keinerlei Signal von A nach B: Alle A-lokalen Operationen kommutieren mit jedem zeitentwickelten B-lokalen Observablen. Für jedes spurtreue lokale Quanteninstrument bleibt daher die unbedingte B-Statistik unverändert. Eine Konditionierung auf ein Messergebnis in A wäre etwas anderes.

## B4.2 Eine weitere Vereinfachung: neue Dynamik oder andere Teilchenvariablen?

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/003-subsection-10-die-drei-n-chsten-entscheidenden.tex:339-387](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/003-subsection-10-die-drei-n-chsten-entscheidenden.tex#L339) · `section:a16af0d89acc08f425a9` · subsubsection

Manche scheinbaren Erweiterungen beschreiben nur andere Bosonvariablen. Auf der reellen, phasengleichen Teilfamilie setze man

## 1. Welche Quellen zusammengeführt werden

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/012-subsection-1-welche-quellen-zusammengef-hrt.tex:1-46](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/012-subsection-1-welche-quellen-zusammengef-hrt.tex#L1) · `section:584bd08c65b5915fd0cf` · subsection

Die Benutzeranlage addbf22d-4764-4685-bdfe-f422bfc563d0/pasted-text.txt, eingefroren unter sources/user_attachment.txt. Die externe Untersuchung Universalraum_Beweisversuch_2026-09-15.md: Spektralzählung, Primzahl-Dynamik, Phasencocycle und kontrollierte Polnäherung. Die externe Untersuchung Universalraum_Urspruenglicher_Austausch_2026-09-15.md: interne Paarumwandlung, zusätzlicher Rotor-Link und älteres Round37-Beis

Explicit source references:

- [experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/user-attachments/39228cf8-f71d-43e6-9a05-4d49d38d721b/pasted-text.txt](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/user-attachments/39228cf8-f71d-43e6-9a05-4d49d38d721b/pasted-text.txt#L1) · explicit_file_citation
- [experiments/theory-contracts/universalraum-two-bank-transfer-20260915/LATE_AUDIT.md](../../../../experiments/theory-contracts/universalraum-two-bank-transfer-20260915/LATE_AUDIT.md#L1) · explicit_file_citation
- [experiments/theory-contracts/universalraum-two-bank-transfer-20260915/sources/user_attachment.txt](../../../../experiments/theory-contracts/universalraum-two-bank-transfer-20260915/sources/user_attachment.txt#L1) · explicit_file_citation

## 6.2 Primzahlen können eine gewählte Dynamik organisieren, erzwingen sie aber nicht

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/013-subsubsection-4-5-dieselbe-urspr-ngliche-fermionantwort-.tex:106-136](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/013-subsubsection-4-5-dieselbe-urspr-ngliche-fermionantwort-.tex#L106) · `section:96445a07838ee8ca073c` · subsubsection

Auf [formula in source] kann man [formula in source] definieren. Die additive Primfaktorbesetzung erklärt dann die Eulerstruktur. Der Definitionsbereich lautet [formula in source] . Dass quantenstatistische Systeme mit Zeta-Zustandssumme existieren, ist bereits ein Ergebnis von https://alainconnes.org/wp-content/uploads/bostconnesscan.pdfBost und Connes. Das ist ein relevanter Anschluss, nicht die fehlende native Her

## 8.1 Konkrete Präparationsgrenze

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/016-subsection-7-einfache-organisation-aber-keine-falsche.tex:116-131](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/016-subsection-7-einfache-organisation-aber-keine-falsche.tex#L116) · `section:a053c7de09481b95f1f5` · subsubsection

Alle genannten Fockkontrollen erhalten N. Ein Wort aus ihnen und jeder einzeln zahlenerhaltende ausgewählte Messzweig bleiben im Ausgangssektor. Sie können aus dem leeren N=0-Zustand nicht den N=64-Grundzustand erzeugen. Diese Aussage betrifft zahlenerhaltende Krauszweige; eine bloß U(1)-kovariante offene Dynamik kann sehr wohl geladene Krausoperatoren haben.

## Das Verfügbarkeits-Theorem (Parallelsession, symmetry_availability.py)Das Verfügbarkeits-Theorem (Parallelsession, )

[documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/018-subsubsection-die-isotypie-zerlegung-von-n-3-beide-spure.tex:31-52](../../../../documents/experiments--theory-contracts--universalraum-singlet-observable-20260915--tmp--pdfs--TFPT_Universalraum_Richtige_Kette_Konsolidierung_2026-09-15_v1.6.7/sections/018-subsubsection-die-isotypie-zerlegung-von-n-3-beide-spure.tex#L31) · `section:a04f6efe04777d89b399` · subsubsection

Die gewährte Dynamik (Wörter aus X, N_b) liegt in der Verflechteralgebra (Dim Σmᵢ² = 16); die 60 Lie-Erzeuger liegen nicht darin - keine Komposition gewährter Dynamik erzeugt sie. Die Symmetrie ist ein eigenständiges Primitiv: gewährt oder nicht gewährt. Der gesamte verbleibende symmetrieverträgliche Spielraum ist exakt 2 Dimensionen breit = die Fähigkeit, die drei dunklen Blöcke getrennt anzusprechen (= beide Casimi

## 4. Gleiche vollständige Zweiladungsspektren, unterschiedliche globale Dynamik

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/002-subsection-7-was-die-vier-rechnungen-gemeinsam-haben-und.tex:309-346](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/002-subsection-7-was-die-vier-rechnungen-gemeinsam-haben-und.tex#L309) · `section:e0abf53daa3d1aa5f478` · subsection

wobei J_L die Matrix aus lauter Einsen bezeichnet. Beide sind symmetrisch, normiert, invariant unter Permutation der Kopien und stimmen für L=1 mit dem nativen Modell überein. Sie liefern jedoch verschiedene spätere Physik. Permutationssymmetrie ist dabei keine geometrische Lokalität.

## 2.2 Welche Gleichheit ist tatsächlich beobachtbar?

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/007-subsubsection-2-1-tats-chliche-korrektur-die-erste-zeile.tex:32-59](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/007-subsubsection-2-1-tats-chliche-korrektur-die-erste-zeile.tex#L32) · `section:999196bb84979b38dde5` · subsubsection

Ein markierter komplexer Wortkern ist stärker als eine Liste gewöhnlicher Quantenkanäle. Beispielsweise liefern [formula in source] und [formula in source] für jeden Zustand denselben Kanal:

## Ergebnis und Leseschlüssel

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/015-subsection-50-gegen-die-lteren-gepinnten-zertifikate.tex:88-134](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/015-subsection-50-gegen-die-lteren-gepinnten-zertifikate.tex#L88) · `section:956a0de4bed7881a86c1` · subsubsection

Die zuletzt ergänzte Vogelperspektive steht im Teil B dieser Fassung: Gesucht wird der kleinste gemeinsame beobachtbare Prozess hinter den verschiedenen Ansichten. Ein nativer Rang-60-Kodierer ist bereits vorhanden. Seine Informationsverluste, die minimalen symmetrieverträglichen Dynamiken und der Unterschied zwischen Korrelation und kausalem Eingriff wurden neu berechnet. Das verschiebt die fundamentale Priorität vo

## B2. Noethers Frage: Welche Unterschiede kann die Dynamik überhaupt verändern?

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/017-subsection-8-die-fundamentale-frage-was-fehlt.tex:235-287](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/017-subsection-8-die-fundamentale-frage-was-fehlt.tex#L235) · `section:60481da614058c6b744a` · subsection

Bei ungebrochener innerer Gruppe G zerfällt ein Zustandsraum als

## B3. Korrelation ist noch kein Eingriff: ein exakt gelöstes Gegenbeispiel

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/017-subsection-8-die-fundamentale-frage-was-fehlt.tex:288-340](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/017-subsection-8-die-fundamentale-frage-was-fehlt.tex#L288) · `section:efae809ab8e6f2b08bb8` · subsection

Die Kreuzkorrelation startet bei null und verändert sich mit der Zeit. Es gibt dennoch keinerlei Signal von A nach B: Alle A-lokalen Operationen kommutieren mit jedem zeitentwickelten B-lokalen Observablen. Für jedes spurtreue lokale Quanteninstrument bleibt daher die unbedingte B-Statistik unverändert. Eine Konditionierung auf ein Messergebnis in A wäre etwas anderes.

## B4.2 Eine weitere Vereinfachung: neue Dynamik oder andere Teilchenvariablen?

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/018-subsubsection-b4-1-zwei-zus-tzliche-terme-am-unver-ndert.tex:85-133](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/018-subsubsection-b4-1-zwei-zus-tzliche-terme-am-unver-ndert.tex#L85) · `section:fa23964ac61e246531cc` · subsubsection

Manche scheinbaren Erweiterungen beschreiben nur andere Bosonvariablen. Auf der reellen, phasengleichen Teilfamilie setze man

## 1. Welche Quellen zusammengeführt werden

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/026-subsection-3-die-richtige.tex:163-208](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/026-subsection-3-die-richtige.tex#L163) · `section:b99a78128c4df5f41a04` · subsection

Die Benutzeranlage addbf22d-4764-4685-bdfe-f422bfc563d0/pasted-text.txt, eingefroren unter sources/user_attachment.txt. Die externe Untersuchung Universalraum_Beweisversuch_2026-09-15.md: Spektralzählung, Primzahl-Dynamik, Phasencocycle und kontrollierte Polnäherung. Die externe Untersuchung Universalraum_Urspruenglicher_Austausch_2026-09-15.md: interne Paarumwandlung, zusätzlicher Rotor-Link und älteres Round37-Beis

Explicit source references:

- [experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/user-attachments/39228cf8-f71d-43e6-9a05-4d49d38d721b/pasted-text.txt](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/user-attachments/39228cf8-f71d-43e6-9a05-4d49d38d721b/pasted-text.txt#L1) · explicit_file_citation
- [experiments/theory-contracts/universalraum-two-bank-transfer-20260915/LATE_AUDIT.md](../../../../experiments/theory-contracts/universalraum-two-bank-transfer-20260915/LATE_AUDIT.md#L1) · explicit_file_citation
- [experiments/theory-contracts/universalraum-two-bank-transfer-20260915/sources/user_attachment.txt](../../../../experiments/theory-contracts/universalraum-two-bank-transfer-20260915/sources/user_attachment.txt#L1) · explicit_file_citation

## 6.2 Primzahlen können eine gewählte Dynamik organisieren, erzwingen sie aber nicht

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/027-subsubsection-4-2-die-l-cke-zum-gesamten-brigen.tex:270-300](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/027-subsubsection-4-2-die-l-cke-zum-gesamten-brigen.tex#L270) · `section:cc4fcf2ba9ea3af45a10` · subsubsection

Auf [formula in source] kann man [formula in source] definieren. Die additive Primfaktorbesetzung erklärt dann die Eulerstruktur. Der Definitionsbereich lautet [formula in source] . Dass quantenstatistische Systeme mit Zeta-Zustandssumme existieren, ist bereits ein Ergebnis von https://alainconnes.org/wp-content/uploads/bostconnesscan.pdfBost und Connes. Das ist ein relevanter Anschluss, nicht die fehlende native Her

## 8.1 Konkrete Präparationsgrenze

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/030-subsection-5-dieselbe-antwort-auf-demselben.tex:290-305](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/030-subsection-5-dieselbe-antwort-auf-demselben.tex#L290) · `section:9e8b66ab27abb9f9adf5` · subsubsection

Alle genannten Fockkontrollen erhalten N. Ein Wort aus ihnen und jeder einzeln zahlenerhaltende ausgewählte Messzweig bleiben im Ausgangssektor. Sie können aus dem leeren N=0-Zustand nicht den N=64-Grundzustand erzeugen. Diese Aussage betrifft zahlenerhaltende Krauszweige; eine bloß U(1)-kovariante offene Dynamik kann sehr wohl geladene Krausoperatoren haben.

## Das Verfügbarkeits-Theorem (Parallelsession, symmetry_availability.py)Das Verfügbarkeits-Theorem (Parallelsession, )

[documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/032-subsection-0-gesamtlage-zwei-spuren-label-gesamtlage-zwe.tex:132-153](../../../../documents/experiments--theory-contracts--universalraum-source-composition-20260915--tmp--pdfs--TFPT_Universalraum_Quellkomposition_Konsolidierung_2026-09-15_v1.6.9/sections/032-subsection-0-gesamtlage-zwei-spuren-label-gesamtlage-zwe.tex#L132) · `section:07e6cbc4eee55cbfbfe0` · subsubsection

Die gewährte Dynamik (Wörter aus X, N_b) liegt in der Verflechteralgebra (Dim Σmᵢ² = 16); die 60 Lie-Erzeuger liegen nicht darin - keine Komposition gewährter Dynamik erzeugt sie. Die Symmetrie ist ein eigenständiges Primitiv: gewährt oder nicht gewährt. Der gesamte verbleibende symmetrieverträgliche Spielraum ist exakt 2 Dimensionen breit = die Fähigkeit, die drei dunklen Blöcke getrennt anzusprechen (= beide Casimi

## 4. Gleiche vollständige Zweiladungsspektren, unterschiedliche globale Dynamik

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/005-subsection-7-was-die-vier-rechnungen-gemeinsam-haben-und.tex:309-346](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/005-subsection-7-was-die-vier-rechnungen-gemeinsam-haben-und.tex#L309) · `section:6e6f40763fa9b9a72e59` · subsection

wobei J_L die Matrix aus lauter Einsen bezeichnet. Beide sind symmetrisch, normiert, invariant unter Permutation der Kopien und stimmen für L=1 mit dem nativen Modell überein. Sie liefern jedoch verschiedene spätere Physik. Permutationssymmetrie ist dabei keine geometrische Lokalität.

## 2.2 Welche Gleichheit ist tatsächlich beobachtbar?

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/010-subsubsection-2-1-tats-chliche-korrektur-die-erste-zeile.tex:32-59](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/010-subsubsection-2-1-tats-chliche-korrektur-die-erste-zeile.tex#L32) · `section:767848613f2baca7dbab` · subsubsection

Ein markierter komplexer Wortkern ist stärker als eine Liste gewöhnlicher Quantenkanäle. Beispielsweise liefern [formula in source] und [formula in source] für jeden Zustand denselben Kanal:

## Ergebnis und Leseschlüssel

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/018-subsection-50-gegen-die-lteren-gepinnten-zertifikate.tex:88-134](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/018-subsection-50-gegen-die-lteren-gepinnten-zertifikate.tex#L88) · `section:392a638f98c5ff1d697a` · subsubsection

Die zuletzt ergänzte Vogelperspektive steht im Teil B dieser Fassung: Gesucht wird der kleinste gemeinsame beobachtbare Prozess hinter den verschiedenen Ansichten. Ein nativer Rang-60-Kodierer ist bereits vorhanden. Seine Informationsverluste, die minimalen symmetrieverträglichen Dynamiken und der Unterschied zwischen Korrelation und kausalem Eingriff wurden neu berechnet. Das verschiebt die fundamentale Priorität vo

## B2. Noethers Frage: Welche Unterschiede kann die Dynamik überhaupt verändern?

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/020-subsection-8-die-fundamentale-frage-was-fehlt.tex:235-287](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/020-subsection-8-die-fundamentale-frage-was-fehlt.tex#L235) · `section:030e0693c77684ad1988` · subsection

Bei ungebrochener innerer Gruppe G zerfällt ein Zustandsraum als

## B3. Korrelation ist noch kein Eingriff: ein exakt gelöstes Gegenbeispiel

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/020-subsection-8-die-fundamentale-frage-was-fehlt.tex:288-340](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/020-subsection-8-die-fundamentale-frage-was-fehlt.tex#L288) · `section:bc975c23a1ebffb4e56e` · subsection

Die Kreuzkorrelation startet bei null und verändert sich mit der Zeit. Es gibt dennoch keinerlei Signal von A nach B: Alle A-lokalen Operationen kommutieren mit jedem zeitentwickelten B-lokalen Observablen. Für jedes spurtreue lokale Quanteninstrument bleibt daher die unbedingte B-Statistik unverändert. Eine Konditionierung auf ein Messergebnis in A wäre etwas anderes.

## B4.2 Eine weitere Vereinfachung: neue Dynamik oder andere Teilchenvariablen?

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/021-subsubsection-b4-1-zwei-zus-tzliche-terme-am-unver-ndert.tex:85-133](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/021-subsubsection-b4-1-zwei-zus-tzliche-terme-am-unver-ndert.tex#L85) · `section:cdad006b6f1b4120a8f6` · subsubsection

Manche scheinbaren Erweiterungen beschreiben nur andere Bosonvariablen. Auf der reellen, phasengleichen Teilfamilie setze man

## 1. Welche Quellen zusammengeführt werden

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/029-subsection-3-die-richtige.tex:163-208](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/029-subsection-3-die-richtige.tex#L163) · `section:b6d20e6827aeff2134ad` · subsection

Die Benutzeranlage addbf22d-4764-4685-bdfe-f422bfc563d0/pasted-text.txt, eingefroren unter sources/user_attachment.txt. Die externe Untersuchung Universalraum_Beweisversuch_2026-09-15.md: Spektralzählung, Primzahl-Dynamik, Phasencocycle und kontrollierte Polnäherung. Die externe Untersuchung Universalraum_Urspruenglicher_Austausch_2026-09-15.md: interne Paarumwandlung, zusätzlicher Rotor-Link und älteres Round37-Beis

Explicit source references:

- [experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/user-attachments/39228cf8-f71d-43e6-9a05-4d49d38d721b/pasted-text.txt](../../../../experiments/theory-contracts/universalraum-fable-kernfragen-20260910/context/fortsetzung/source-snapshots/user-attachments/39228cf8-f71d-43e6-9a05-4d49d38d721b/pasted-text.txt#L1) · explicit_file_citation
- [experiments/theory-contracts/universalraum-two-bank-transfer-20260915/sources/user_attachment.txt](../../../../experiments/theory-contracts/universalraum-two-bank-transfer-20260915/sources/user_attachment.txt#L1) · explicit_file_citation
- [experiments/theory-contracts/universalraum-two-bank-transfer-20260915/LATE_AUDIT.md](../../../../experiments/theory-contracts/universalraum-two-bank-transfer-20260915/LATE_AUDIT.md#L1) · explicit_file_citation

## 6.2 Primzahlen können eine gewählte Dynamik organisieren, erzwingen sie aber nicht

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/030-subsubsection-4-2-die-l-cke-zum-gesamten-brigen.tex:270-300](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/030-subsubsection-4-2-die-l-cke-zum-gesamten-brigen.tex#L270) · `section:a6cb9954a55bedb598a0` · subsubsection

Auf [formula in source] kann man [formula in source] definieren. Die additive Primfaktorbesetzung erklärt dann die Eulerstruktur. Der Definitionsbereich lautet [formula in source] . Dass quantenstatistische Systeme mit Zeta-Zustandssumme existieren, ist bereits ein Ergebnis von https://alainconnes.org/wp-content/uploads/bostconnesscan.pdfBost und Connes. Das ist ein relevanter Anschluss, nicht die fehlende native Her

## 8.1 Konkrete Präparationsgrenze

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/033-subsection-5-dieselbe-antwort-auf-demselben.tex:290-305](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/033-subsection-5-dieselbe-antwort-auf-demselben.tex#L290) · `section:a842aa80d4afcd93581c` · subsubsection

Alle genannten Fockkontrollen erhalten N. Ein Wort aus ihnen und jeder einzeln zahlenerhaltende ausgewählte Messzweig bleiben im Ausgangssektor. Sie können aus dem leeren N=0-Zustand nicht den N=64-Grundzustand erzeugen. Diese Aussage betrifft zahlenerhaltende Krauszweige; eine bloß U(1)-kovariante offene Dynamik kann sehr wohl geladene Krausoperatoren haben.

## Das Verfügbarkeits-Theorem (Parallelsession, symmetry_availability.py)Das Verfügbarkeits-Theorem (Parallelsession, )

[documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/035-subsection-0-gesamtlage-zwei-spuren-label-gesamtlage-zwe.tex:132-153](../../../../documents/experiments--theory-contracts--universalraum-spin-lift-basics-20260915--tmp--pdfs--TFPT_Universalraum_Spinlift_Hauptdokument_2026-09-15_v1.6.10/sections/035-subsection-0-gesamtlage-zwei-spuren-label-gesamtlage-zwe.tex#L132) · `section:9251f4ff96eb011613de` · subsubsection

Die gewährte Dynamik (Wörter aus X, N_b) liegt in der Verflechteralgebra (Dim Σmᵢ² = 16); die 60 Lie-Erzeuger liegen nicht darin - keine Komposition gewährter Dynamik erzeugt sie. Die Symmetrie ist ein eigenständiges Primitiv: gewährt oder nicht gewährt. Der gesamte verbleibende symmetrieverträgliche Spielraum ist exakt 2 Dimensionen breit = die Fähigkeit, die drei dunklen Blöcke getrennt anzusprechen (= beide Casimi

## The compact status formula: three honest layers

[documents/introduction/sections/002-section-the-dependency-dag-and-the-proof-ledger.tex:60-94](../../../../documents/introduction/sections/002-section-the-dependency-dag-and-the-proof-ledger.tex#L60) · `section:da803316f5620018a265` · subsection

[formula in source] TFPT v \ closes the discrete compiler; the discrete Standard-Model structure (charges, families, flavor operators) is derived under the axioms, while [formula in source] D existence, chiral measure, mirror gap, couplings and neutrino dynamics remain open. It closes the EM fixed point, the admissible gapped IR sector and the [formula in source] spectral-action shadow. It does not yet certify a stri

Explicit source references:

- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v68_seeley_dewitt_residual.py](../../../../verification/v68_seeley_dewitt_residual.py#L1) · explicit_script_identifier
- [verification/v56_unique_attractor.py](../../../../verification/v56_unique_attractor.py#L1) · explicit_file_citation
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [verification/v14_carrier_uniqueness.py](../../../../verification/v14_carrier_uniqueness.py#L1) · explicit_file_citation
- [verification/v60_lambda_metrology_branch.py](../../../../verification/v60_lambda_metrology_branch.py#L1) · explicit_file_citation
- [verification/v77_e8_conformal_net.py](../../../../verification/v77_e8_conformal_net.py#L1) · explicit_script_identifier
- [verification/v6_bootstrap.py](../../../../verification/v6_bootstrap.py#L1) · explicit_file_citation
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v75_upoint_to_vgeo.py](../../../../verification/v75_upoint_to_vgeo.py#L1) · explicit_script_identifier
- [verification/v78_vgeo_floor.py](../../../../verification/v78_vgeo_floor.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v77) G6 via E8 level-1 net: level-1 central charges + conformal embedding c-sum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L377)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ---- (v68) Seeley-DeWitt a2 for the carrier; 1/G UV-sensitive => k=c3/2 is a normalization ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L322)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)
- [(* ---- (v78) v_geo floor: S_dS*rho_Lambda = 1/(128 c3^4) = 32 pi^4 (one scale pinned by one measurement) ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L382)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ---- (v75) Gate 1: U_point -> v_geo. lepton product 32/9; (ratios,product)<=>(individuals) bijection ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L366)
- [(* ---- (v60) Lambda branch selection: prefactor 3/(4pi^2), branch mis-scale 64pi^3/3 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L283)
- [(* ---- (v56) unique attractor: gapped transport, Coxeter planes ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L257)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean#L1)

## The hard reduction: the live residual as typed interfaces

[documents/introduction/sections/002-section-the-dependency-dag-and-the-proof-ledger.tex:95-387](../../../../documents/introduction/sections/002-section-the-dependency-dag-and-the-proof-ledger.tex#L95) · `section:81a6dc002dc46cf7b6c5` · subsection

After the compiler closure, the entire residual is the pair of accounting formulae (the compiler rest, and the strict-physical-TOE rest of the 2026-08-27/28 contract wave; tfpt_research_contracts): [formula in source] [formula in source] The first is the compiler rest — one dimensionful scale anchor, one metric-sector inclusion, and a set of deliberately typed QFT/cosmology interfaces. The second is the strict-physic

Explicit source references:

- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_script_identifier
- [verification/v373_ftransfer_relic.py](../../../../verification/v373_ftransfer_relic.py#L1) · explicit_script_identifier
- [CONTRACT.F.01](../../../../verification/status_ledger.csv#L302) · explicit_claim_identifier · original status at linked ledger row
- [verification/v232_e8_kleinian_seam.py](../../../../verification/v232_e8_kleinian_seam.py#L1) · explicit_file_citation
- [AX.P1.01](../../../../verification/status_ledger.csv#L2) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1006_mmst_lemma_battery.py](../../../../verification/v1006_mmst_lemma_battery.py#L1) · explicit_file_citation
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v146_moebius_d4.py](../../../../verification/v146_moebius_d4.py#L1) · explicit_script_identifier
- [verification/v372_ftransfer_boltzmann.py](../../../../verification/v372_ftransfer_boltzmann.py#L1) · explicit_script_identifier
- [verification/v160_seam_gaussianity_from_pf.py](../../../../verification/v160_seam_gaussianity_from_pf.py#L1) · explicit_script_identifier
- [verification/v75_upoint_to_vgeo.py](../../../../verification/v75_upoint_to_vgeo.py#L1) · explicit_script_identifier
- [verification/v76_gmetric_reduction.py](../../../../verification/v76_gmetric_reduction.py#L1) · explicit_script_identifier
- [GRAV.SPIN2.EMERGENCE.01](../../../../verification/status_ledger.csv#L762) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1034_parent_mirror.py](../../../../verification/v1034_parent_mirror.py#L1) · explicit_file_citation
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_script_identifier
- [SEAM.MMST.TYPEIII.CHARGED.01](../../../../verification/status_ledger.csv#L767) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1032_covariant_curvature.py](../../../../verification/v1032_covariant_curvature.py#L1) · explicit_file_citation
- [verification/v1008_master_assembly_scaffold.py](../../../../verification/v1008_master_assembly_scaffold.py#L1) · explicit_file_citation
- [verification/v1011_3p1d_ladder_battery.py](../../../../verification/v1011_3p1d_ladder_battery.py#L1) · explicit_file_citation
- [verification/v234_seam_holomorphy_selection.py](../../../../verification/v234_seam_holomorphy_selection.py#L1) · explicit_file_citation
- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v137_qplus_cohomology.py](../../../../verification/v137_qplus_cohomology.py#L1) · explicit_script_identifier
- [verification/v1010_simplicity_bridge_census.py](../../../../verification/v1010_simplicity_bridge_census.py#L1) · explicit_file_citation
- [verification/v1028_gauss_d4_cap.py](../../../../verification/v1028_gauss_d4_cap.py#L1) · explicit_file_citation
- [verification/v175_net_existence_full_cone.py](../../../../verification/v175_net_existence_full_cone.py#L1) · explicit_script_identifier
- [verification/v213_ftransfer_functor.py](../../../../verification/v213_ftransfer_functor.py#L1) · explicit_script_identifier
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [verification/v60_lambda_metrology_branch.py](../../../../verification/v60_lambda_metrology_branch.py#L1) · explicit_script_identifier
- [FLAV.NU.TEXTURE.MECHANISM.01](../../../../verification/status_ledger.csv#L768) · explicit_claim_identifier · original status at linked ledger row
- [verification/v999_weak_collision_hp.py](../../../../verification/v999_weak_collision_hp.py#L1) · explicit_file_citation
- [verification/v1030_alg2_sector_frame.py](../../../../verification/v1030_alg2_sector_frame.py#L1) · explicit_file_citation
- [verification/v1001_nu_pentagon_chain.py](../../../../verification/v1001_nu_pentagon_chain.py#L1) · explicit_file_citation
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v1029_local_tensor_constraints.py](../../../../verification/v1029_local_tensor_constraints.py#L1) · explicit_file_citation
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_script_identifier
- [verification/v6_bootstrap.py](../../../../verification/v6_bootstrap.py#L1) · explicit_script_identifier
- [verification/v147_clock_gaussian.py](../../../../verification/v147_clock_gaussian.py#L1) · explicit_script_identifier
- [verification/v14_carrier_uniqueness.py](../../../../verification/v14_carrier_uniqueness.py#L1) · explicit_file_citation
- [verification/v308_seam_equiv_chain.py](../../../../verification/v308_seam_equiv_chain.py#L1) · explicit_script_identifier
- [E8.DIRECTED.READOUT.01](../../../../verification/status_ledger.csv#L771) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1019_coxeter_euler_completion.py](../../../../verification/v1019_coxeter_euler_completion.py#L1) · explicit_file_citation
- [verification/v725_phys_vgeo_torsor_audit.py](../../../../verification/v725_phys_vgeo_torsor_audit.py#L1) · explicit_script_identifier
- [verification/v153_no_unit_theorem.py](../../../../verification/v153_no_unit_theorem.py#L1) · explicit_script_identifier
- [verification/v112_selfcounting_channel.py](../../../../verification/v112_selfcounting_channel.py#L1) · explicit_script_identifier
- [verification/v995_kms_subnet_rigidity.py](../../../../verification/v995_kms_subnet_rigidity.py#L1) · explicit_file_citation
- [verification/v303_ftransfer_dynamics.py](../../../../verification/v303_ftransfer_dynamics.py#L1) · explicit_script_identifier
- [verification/v1031_quantum_tensor_curvature.py](../../../../verification/v1031_quantum_tensor_curvature.py#L1) · explicit_file_citation
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [verification/v113_quasifree_kernel.py](../../../../verification/v113_quasifree_kernel.py#L1) · explicit_script_identifier
- [verification/v1009_rho0_minimizer.py](../../../../verification/v1009_rho0_minimizer.py#L1) · explicit_file_citation
- [verification/v219_icosahedral_mckay.py](../../../../verification/v219_icosahedral_mckay.py#L1) · explicit_file_citation
- [verification/v141_deck_selection.py](../../../../verification/v141_deck_selection.py#L1) · explicit_script_identifier
- [GRAV.WEINBERG.WITTEN.01](../../../../verification/status_ledger.csv#L770) · explicit_claim_identifier · original status at linked ledger row
- [verification/v86_nstar_reheating.py](../../../../verification/v86_nstar_reheating.py#L1) · explicit_script_identifier
- [verification/v1013_thermodynamic_dynamics.py](../../../../verification/v1013_thermodynamic_dynamics.py#L1) · explicit_file_citation
- [verification/v624_external_lattice_audit.py](../../../../verification/v624_external_lattice_audit.py#L1) · explicit_script_identifier
- [verification/v1004_milnor_frobenius_gray.py](../../../../verification/v1004_milnor_frobenius_gray.py#L1) · explicit_file_citation
- [verification/v89_carrier_index_lemma.py](../../../../verification/v89_carrier_index_lemma.py#L1) · explicit_script_identifier
- [CHIRAL4D.MIRROR.DET16.01](../../../../verification/status_ledger.csv#L766) · explicit_claim_identifier · original status at linked ledger row
- [verification/v111_quadratic_transport.py](../../../../verification/v111_quadratic_transport.py#L1) · explicit_script_identifier
- [verification/v1012_toe_gate_witnesses.py](../../../../verification/v1012_toe_gate_witnesses.py#L1) · explicit_file_citation
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v154_simple_current_theorem.py](../../../../verification/v154_simple_current_theorem.py#L1) · explicit_script_identifier
- [verification/v994_mmst_criteria_battery.py](../../../../verification/v994_mmst_criteria_battery.py#L1) · explicit_file_citation
- [verification/v996_gauge_grammar_kill.py](../../../../verification/v996_gauge_grammar_kill.py#L1) · explicit_file_citation
- [verification/v528_seam_bit_twist_class_definition.py](../../../../verification/v528_seam_bit_twist_class_definition.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.TWISTOR.01](../../../../verification/status_ledger.csv#L584) · explicit_claim_identifier · original status at linked ledger row
- [verification/v998_seam_modular_closure.py](../../../../verification/v998_seam_modular_closure.py#L1) · explicit_file_citation
- [verification/v1026_telb_round4_closure.py](../../../../verification/v1026_telb_round4_closure.py#L1) · explicit_file_citation
- [verification/v1033_charged_disorder.py](../../../../verification/v1033_charged_disorder.py#L1) · explicit_file_citation
- [AX.P2.01](../../../../verification/status_ledger.csv#L3) · explicit_claim_identifier · original status at linked ledger row
- [verification/v144_detratio_family.py](../../../../verification/v144_detratio_family.py#L1) · explicit_script_identifier
- [verification/v1020_spectator_sector_obstruction.py](../../../../verification/v1020_spectator_sector_obstruction.py#L1) · explicit_script_identifier
- [verification/v1027_signed_det_car_wall.py](../../../../verification/v1027_signed_det_car_wall.py#L1) · explicit_file_citation
- [verification/v993_minimal_defect_selector.py](../../../../verification/v993_minimal_defect_selector.py#L1) · explicit_file_citation
- [verification/v108_pascal_ladder.py](../../../../verification/v108_pascal_ladder.py#L1) · explicit_script_identifier
- [verification/v1005_relative_pencil_bridges.py](../../../../verification/v1005_relative_pencil_bridges.py#L1) · explicit_file_citation
- [E8.COXETER.EULER.COMPLETION.01](../../../../verification/status_ledger.csv#L772) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1035_matter_coupling.py](../../../../verification/v1035_matter_coupling.py#L1) · explicit_file_citation
- [GATE.METRIC.06](../../../../verification/status_ledger.csv#L163) · explicit_claim_identifier · original status at linked ledger row
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_file_citation
- [verification/v23_anchor_generator.py](../../../../verification/v23_anchor_generator.py#L1) · explicit_file_citation
- [TFPT.TOE.COMPLETE.01](../../../../verification/status_ledger.csv#L764) · explicit_claim_identifier · original status at linked ledger row
- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier
- [verification/v136_dual_normal_selector.py](../../../../verification/v136_dual_normal_selector.py#L1) · explicit_script_identifier
- [verification/v36_spectral_action_g2.py](../../../../verification/v36_spectral_action_g2.py#L1) · explicit_file_citation
- [verification/v371_ftransfer_pole.py](../../../../verification/v371_ftransfer_pole.py#L1) · explicit_script_identifier
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [verification/v133_zeta_budget.py](../../../../verification/v133_zeta_budget.py#L1) · explicit_script_identifier
- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier
- [verification/v1016_state_gap_batteries.py](../../../../verification/v1016_state_gap_batteries.py#L1) · explicit_file_citation
- [verification/v1015_axiom_core_closure.py](../../../../verification/v1015_axiom_core_closure.py#L1) · explicit_file_citation
- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier
- [verification/v138_vw_firewall.py](../../../../verification/v138_vw_firewall.py#L1) · explicit_script_identifier
- [verification/v374_ftransfer_qcd.py](../../../../verification/v374_ftransfer_qcd.py#L1) · explicit_script_identifier
- [verification/v149_cusp_normal.py](../../../../verification/v149_cusp_normal.py#L1) · explicit_script_identifier
- [verification/v510_seam_bit_freedom.py](../../../../verification/v510_seam_bit_freedom.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [GAUGE.DETLINE.FIXPOINT.01](../../../../verification/status_ledger.csv#L761) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1007_det16_stability_battery.py](../../../../verification/v1007_det16_stability_battery.py#L1) · explicit_file_citation
- [verification/v1018_e8_directed_readout.py](../../../../verification/v1018_e8_directed_readout.py#L1) · explicit_file_citation
- [verification/v997_spin2_smg_witnesses.py](../../../../verification/v997_spin2_smg_witnesses.py#L1) · explicit_file_citation
- [verification/v165_premise_a_closure.py](../../../../verification/v165_premise_a_closure.py#L1) · explicit_script_identifier
- [verification/v56_unique_attractor.py](../../../../verification/v56_unique_attractor.py#L1) · explicit_script_identifier
- [verification/v1000_spin2_sk_witnesses.py](../../../../verification/v1000_spin2_sk_witnesses.py#L1) · explicit_file_citation
- [FTRANSFER.SK.RHO0.01](../../../../verification/status_ledger.csv#L763) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.MILNOR.LOCALRING.01](../../../../verification/status_ledger.csv#L769) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1014_bridge_refinements.py](../../../../verification/v1014_bridge_refinements.py#L1) · explicit_file_citation
- [verification/v124_resummed_clock.py](../../../../verification/v124_resummed_clock.py#L1) · explicit_script_identifier
- [verification/v512_seam_tau_flag.py](../../../../verification/v512_seam_tau_flag.py#L1) · explicit_script_identifier
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v235_seam_chern_simons.py](../../../../verification/v235_seam_chern_simons.py#L1) · explicit_file_citation
- [verification/v145_pairing_atoms.py](../../../../verification/v145_pairing_atoms.py#L1) · explicit_script_identifier
- [verification/v110_calderon_sheet.py](../../../../verification/v110_calderon_sheet.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* v1034: exact source parent envelope. The source's 3/5 is a dispersion amplitude; a standard symmetric NN bond has half of it. *)](../../../../verification/wolfram/tfpt_round7.wl#L161)
- [(* ---- (v108) Pascal ladder ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L486)
- [(* v1028: exact ordinary D4 cap on the selected 12D even pair space. *)](../../../../verification/wolfram/tfpt_round4.wl#L113)
- [(* v219 McKay: affine E8 marks = 2I irrep degrees, sums 30 and 120 *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1783)
- [(* v1035: physical-position Fourier convention D -> i k. Retain I: c0=-s.q and c_i=i (B p)_i. The source Ward equations then propagate both affine constraints exactly. *)](../../../../verification/wolfram/tfpt_round7.wl#L182)
- [(* ---- (v133) zeta budget ---- Numerical continuations are Python-only (mpmath); the exact arithmetic of both routes is mirrored, plus a numeric spot-check in WL. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L995)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ---- (v111) quadratic transport ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L543)
- [experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/PascalLadder.lean#L1)
- [(* ==== v512 round: SEAM.TAU.FLAG.01 -- "the flag-transitivity equivalence web": the tau = i attack on the free RP seam circle. Test family M(delta) = {+-1, +-e^(i delta)}; bare mark-transitivity is delta-blind (pair-exchanging V4), FLAG transitivity <=> delta = pi/2 <=> tau = i. Six exact mirrors: (i) the V4/D4 census over five exact members incl. flag/arc orbit counts, (ii) the solveset equivalences + family identities, (iii) the odd-split closed form + arc Laplacian + K3 indicator, (iv) the twist mode identity + AB delta-blindness + j-rationality control, (v) the H4 centraliser rigidity (Schur), and (vi) the negative controls (hexagonal/silver/Z16 census/n = 3). Mirrors v512. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5676)
- [(* v1032: a separate Lorentz-covariant curvature calculation at one exact nonzero null momentum. Pi is the de Donder trace-reversed numerator; positivity is inferred only after the curvature sandwich. *)](../../../../verification/wolfram/tfpt_round7.wl#L92)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [(* ==== v235: the closing step in abelian Chern-Simons -- holomorphic <=> det K = 1 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1951)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [(* ==== v232: the seam as the E8 Kleinian singularity (du Val resolution + link) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1895)
- [(* ==== v528 round: SEAM.BIT.TWISTCLASS.01 -- bit = twist-class choice. Exact mirrors: (i) cardinality lemma; (ii) 224-axis flip census; (iii) O-values. Pfaffian spectra Python-only. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8214)
- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [(* ---- (v136) dual-normal selector ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1051)
- [(* v1027: exact signed-wall determinant and unsquared low branch. *)](../../../../verification/wolfram/tfpt_round4.wl#L29)
- [(* ---- (v153) No-Unit Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1369)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [(* ==== v219-v230 round: icosahedral McKay, CM-norm duality, structural finds ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1781)
- [(* ---- (v141) deck selection ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1135)
- [(* ---- (v113) quasi-free kernel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L608)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)
- [(* ---- (v994) SEAM.EQUIV.MMST.01 C3: Theta_E8 / eta^8 = {1,248,4124,34752} ---- Python-only: C1 entropy, C2 QWZ spectrum, C4 two-point numeric, C5 sector recombination beyond det, C6 cited, C7 open, GSO/N=8 must-fail combinatorics. *)](../../../../verification/wolfram/tfpt_readouts.wl#L681)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [(* ---- (v145) pairing atoms ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1227)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CohomologyGrading.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CohomologyGrading.lean#L1)
- [(* v1033: exact Z4 deck algebra and the load-bearing weight mismatch. *)](../../../../verification/wolfram/tfpt_round7.wl#L142)
- [(* ---- (v60) Lambda branch selection: prefactor 3/(4pi^2), branch mis-scale 64pi^3/3 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L283)
- [(* ---- (v110) Calderon-sheet selection ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L525)
- [(* ---- (v147) clock Gaussian model ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1265)
- [(* ---- (v86) reheating arithmetic (scale chain; pivot solve stays Python) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L84)
- [experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean#L1)
- [(* ---- (v993) AX.P2.01 census: SNF disc groups, unique (D5,A3), [lambda]^2=[v], CE dims ---- Python-only: random CE axiom sampling on M_d^{Z4} / crossed product. *)](../../../../verification/wolfram/tfpt_readouts.wl#L558)
- [(* ---- (v144) det-ratio family cancellation ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1208)
- [experiments/lean4-carrier-rigidity/TfptCarrier/LorentzCongruence.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/LorentzCongruence.lean#L1)
- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ---- (v146) Moebius D4 realisation ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1240)
- [(* ---- (v23) anchor generator a=(1,1,2) ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L96)
- [(* ---- (v89) carrier index lemma ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L118)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [(* v1029: exact generic local tensor identities and the k=0 instability. *)](../../../../verification/wolfram/tfpt_round4.wl#L136)
- [(* ==== v234: Seam-Holomorphy selection -- one condition, three faces, E8 ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1930)
- [(* ---- (v56) unique attractor: gapped transport, Coxeter planes ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L257)
- [(* v1031: exact rational-fibre curvature projector and positivity. *)](../../../../verification/wolfram/tfpt_round7.wl#L80)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ---- (v76) Gate 2: decoupling margin Delta_eff = 6log(3/2) - 31/(4pi^2) > 0 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L372)
- [(* ---- (v75) Gate 1: U_point -> v_geo. lepton product 32/9; (ratios,product)<=>(individuals) bijection ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L366)
- [(* ---- (v149) cusp normal ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1304)
- [(* ---- (v137) Q+ cohomology grading ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1087)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CarrierRankUniqueness.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [(* ---- (v124) resummed clock ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L857)
- [(* ---- (v36) G2: spectral-action Seeley-DeWitt coefficients -> R + R^2 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L113)
- [(* ---- (v112) self-counting channel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L583)
- [(* ==== v510 round: SEAM.BIT.FREEDOM.01 -- the freedom attack on the v506/v507 alignment bit: the POSITION half of the bit is topology. Exact content mirrored: (i) the DIHEDRAL CENSUS -- on the refined 2N-point seam circle D_16 has 17 involutions with EXACTLY ONE free (the antipodal shift by 8) and all 16 reflections have exactly 2 circle fixed points (8 site-axis + 8 edge-axis; robustness N = 8, 12; odd control N = 9: zero free), plus the quotient topology chi(C/antipode) = 0 (circle) vs chi(C/reflection) = 1 (interval, 2 boundary corners); (ii) the AUT COMPLETENESS -- Aut(C_8) = D_8 by brute force over all 8! = 40320 permutations and Aut(C_16) = D_16 by degree-2 propagation (all 32 candidates), so the census covers ALL circle symmetries; (iii) the Cl(16) DICHOTOMY CENSUS -- on explicit 256x256 Jordan-Wigner gammas the free deck's NS implementer squares to 256 gamma_1...gamma_16 = 256 (-1)^F (nonsplit) while the implementers of ALL 16 reflections exist (one of the two spin lifts +-R) and square to SCALARS {2^7 (site-axis), 2^8 (edge-axis)} (split): NONSPLIT <=> FREE over the complete census; root dichotomy 2 (the +-quarter shifts) vs 0 total over all 16 reflections; (iv) the REAL-STRUCTURE TABLE -- the deck-commuting antiholomorphic involutions solve symbolically to exactly two families (M diagonal: mu conj z, fix line through the poles; M antidiagonal: mu/conj z, fix circle |z| = sqrt mu), and the unique mark-fixing real structure per configuration is equatorial/FREE for mu4 and generic-unit, real-axis/NOT-free for silver and generic-real, and NONEXISTENT for the hexagonal frame; (v) HARMONIC + FREE => mu4 -- the deck-pair cross-ratio (1-b)^2/(1+b)^2 = -1 solves exactly to b = +-i, with the lambda table (2, 1/2, 3/4, 5/4) real and the hexagonal lambda complex; (vi) the HONEST COUNTERWITNESSES -- the discrete witness {0,1,8,9} (shift-invariant, crossing, but not quarter-shift-invariant) and the continuous witness {+-1, +-(3+4i)/5} (pair cross-ratio -1/4, j = 148176/25 != 1728, clock does not preserve): the Fock dichotomy measures the CLASS, not the modulus. The Moebius/Klein covering models and the pillowcase branch control stay Python-side (v510). Mirrors v510. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5155)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v138) VW firewall ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1103)
- [(* v1030: actual connected-sector words span M4; a disconnected mutant leaves one extra commutant degree. *)](../../../../verification/wolfram/tfpt_round4.wl#L163)
- [(* ---- (v997) GRAV.SPIN2.EMERGENCE.01: Barnes-Rivers ranks (5,3,1,1) ---- Python-only: Weyl^2 ghost residues, random-p numeric, SMG N=8 Fidkowski-Kitaev toy (spectrum / Spin(7) / Schmidt / edge). *)](../../../../verification/wolfram/tfpt_readouts.wl#L700)
- [(* ---- (v175) net existence + full-cone RP: E8 Cartan even unimodular, character order 5, Gamma(t) functoriality integers ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1523)
- [(* ---- (v154) Simple-Current Extension Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1381)
- [(* v1028: exact characteristic-five reduction of the full 192-mode one-site character. P^12 == P^2 P(z^5)^2 mod 5 is evaluated by coefficient dictionaries, never by floating group quadrature. *)](../../../../verification/wolfram/tfpt_round4.wl#L41)

## Cosmology is the Pascal action grammar of the carrier

[documents/introduction/sections/003-section-the-compiler-why-30-2-cdot3-cdot5.tex:139-175](../../../../documents/introduction/sections/003-section-the-compiler-why-30-2-cdot3-cdot5.tex#L139) · `section:f098dd9c071ad5279d43` · subsection

With [formula in source] the cosmological readouts are a power tower on the five-slot carrier graph [formula in source] : [formula in source] EW is the unit [formula in source] , Hubble is the product over the [formula in source] vertices, [formula in source] is the product over the [formula in source] edges — the same Pascal triple that reads one generation, the action ladder and the [formula in source] root count [

Explicit source references:

- [verification/v20_lepton_c_derivation.py](../../../../verification/v20_lepton_c_derivation.py#L1) · explicit_file_citation
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier
- [verification/v46_grand_mass_volume.py](../../../../verification/v46_grand_mass_volume.py#L1) · explicit_file_citation
- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier
- [verification/v18_quark_yukawa.py](../../../../verification/v18_quark_yukawa.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)
- [(* ---- (v20/v21) lepton c's: delta=1/2 resolvent + product rule ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L78)
- [(* ---- (v46) Grand Mass Volume: absolute mass-scaling = det ~ u^25 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L178)

## What this means for the theory-of-everything question

[documents/introduction/sections/003-section-the-compiler-why-30-2-cdot3-cdot5.tex:217-282](../../../../documents/introduction/sections/003-section-the-compiler-why-30-2-cdot3-cdot5.tex#L217) · `section:72b43011ac4fb932c47a` · section

0.46 Flavor residual — ratios closed; only the absolute scale is the wall-selection. Mehta–Seshadri supplies the equivalence-class framework (stable parabolic [formula in source] unitary representation); the algebraic selector [formula in source] , the [formula in source] spectrum and the splitting type are fixed (now derived: [formula in source] lattice, [formula in source] from [formula in source] , v69/ v70; since

Explicit source references:

- [verification/v137_qplus_cohomology.py](../../../../verification/v137_qplus_cohomology.py#L1) · explicit_script_identifier
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v60_lambda_metrology_branch.py](../../../../verification/v60_lambda_metrology_branch.py#L1) · explicit_script_identifier
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_script_identifier
- [verification/v136_dual_normal_selector.py](../../../../verification/v136_dual_normal_selector.py#L1) · explicit_script_identifier
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [verification/v364_vgeo_sharpen.py](../../../../verification/v364_vgeo_sharpen.py#L1) · explicit_script_identifier
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v374_ftransfer_qcd.py](../../../../verification/v374_ftransfer_qcd.py#L1) · explicit_script_identifier
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v473_entropic_action_bridge.py](../../../../verification/v473_entropic_action_bridge.py#L1) · explicit_script_identifier
- [verification/v70_q_integer_lift.py](../../../../verification/v70_q_integer_lift.py#L1) · explicit_script_identifier
- [verification/v69_d4_q_geometry.py](../../../../verification/v69_d4_q_geometry.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_script_identifier
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_script_identifier
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier
- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier
- [verification/v139_selector_triangle.py](../../../../verification/v139_selector_triangle.py#L1) · explicit_script_identifier
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_script_identifier
- [verification/v360_grav_gap_corrections.py](../../../../verification/v360_grav_gap_corrections.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v70) integer-lift of Q: Z4 puncture action unimodular, det Q = N_fam ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L337)
- [(* v473 (ii): entropy potential + Lambda-channel target *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2975)
- [(* ---- (v60) Lambda branch selection: prefactor 3/(4pi^2), branch mis-scale 64pi^3/3 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L283)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)
- [(* ==== v473 round: the entropic-action bridge (Bianconi, PRD 111, 066001 (2025); arXiv:2408.14391) -- the exact algebraic content: the carrier Hodge count 1+5+10 = 16 = 2^(g_car-1) (vs the 4d fiber 11), the coefficient pin beta'_B = c3/6 = 1/(48 pi) from 3 beta' = 1/(16 pi), the entropy potential s(x) = x-1-ln x (nonnegative, quadratic minimum), the Lambda-channel target Tr Q^2 = 32 c3^4 = 1/(128 pi^4) = (2/3) delta_top, and the R^2 kill-test gap 3 c3^-9 = 3 (8 pi)^9; mirrors v473 (the Frullani quadrature is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2958)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ---- (v137) Q+ cohomology grading ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1087)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* v473 (i): carrier Hodge count + coefficient pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2968)
- [(* ---- (v136) dual-normal selector ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1051)
- [(* v473 (iii): the R^2 kill-test gap *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2985)
- [(* ---- (v69) D4-equivariant Q-geometry: Q_+ from cusp weights, Q_- E-block coupling ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L330)
- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)
- [(* ---- (v139) selector triangle ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1110)

## The closure ledger: what is done, and what genuinely remains

[documents/introduction/sections/003-section-the-compiler-why-30-2-cdot3-cdot5.tex:283-377](../../../../documents/introduction/sections/003-section-the-compiler-why-30-2-cdot3-cdot5.tex#L283) · `section:213a1f15f7373957fc2f` · subsection

A complete theory must close the discrete structure, the dimensionless constants, the mass spectrum, gravity, and the quantum measure. The current state, graded honestly: 6pt @p4.6cmcY@ TOE piece & Status & If not complete: what closes it \\ gauge group, [formula in source] , hypercharge, [formula in source] & & — complete \\ [formula in source] , group closure & & — complete \\ bootstrap [formula in source] , `` [fo

Explicit source references:

- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v330_qgamb_admissible_measure.py](../../../../verification/v330_qgamb_admissible_measure.py#L1) · explicit_script_identifier
- [verification/v986_nu_scalaron_rung.py](../../../../verification/v986_nu_scalaron_rung.py#L1) · explicit_script_identifier
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v274_scale_overdetermination.py](../../../../verification/v274_scale_overdetermination.py#L1) · explicit_script_identifier
- [verification/v124_resummed_clock.py](../../../../verification/v124_resummed_clock.py#L1) · explicit_script_identifier
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_script_identifier
- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_script_identifier
- [verification/v265_qft4d_fork_freeze.py](../../../../verification/v265_qft4d_fork_freeze.py#L1) · explicit_script_identifier
- [verification/v122_margin_theorem.py](../../../../verification/v122_margin_theorem.py#L1) · explicit_script_identifier
- [verification/v308_seam_equiv_chain.py](../../../../verification/v308_seam_equiv_chain.py#L1) · explicit_script_identifier
- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_script_identifier
- [verification/v282_e8_tau_i_unification.py](../../../../verification/v282_e8_tau_i_unification.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v481_seesaw_carrier_ladder.py](../../../../verification/v481_seesaw_carrier_ladder.py#L1) · explicit_script_identifier
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v468_dm2_ratio_jarlskog.py](../../../../verification/v468_dm2_ratio_jarlskog.py#L1) · explicit_script_identifier
- [verification/v276_qgeo_flat_closes_commutator.py](../../../../verification/v276_qgeo_flat_closes_commutator.py#L1) · explicit_script_identifier
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier
- [EM.FP.01](../../../../verification/status_ledger.csv#L11) · explicit_claim_identifier · original status at linked ledger row
- [verification/v133_zeta_budget.py](../../../../verification/v133_zeta_budget.py#L1) · explicit_script_identifier
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_script_identifier
- [verification/v332_qgamb_metric_sector.py](../../../../verification/v332_qgamb_metric_sector.py#L1) · explicit_script_identifier
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v374_ftransfer_qcd.py](../../../../verification/v374_ftransfer_qcd.py#L1) · explicit_script_identifier
- [verification/v118_hexagon_family_dictionary.py](../../../../verification/v118_hexagon_family_dictionary.py#L1) · explicit_script_identifier
- [verification/v272_nu_mass_scale.py](../../../../verification/v272_nu_mass_scale.py#L1) · explicit_script_identifier
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_script_identifier
- [verification/v139_selector_triangle.py](../../../../verification/v139_selector_triangle.py#L1) · explicit_script_identifier
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_script_identifier
- [CELEST.SEAM.01](../../../../verification/status_ledger.csv#L556) · explicit_claim_identifier · original status at linked ledger row
- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [verification/v334_conformal_resolution.py](../../../../verification/v334_conformal_resolution.py#L1) · explicit_script_identifier
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [verification/v136_dual_normal_selector.py](../../../../verification/v136_dual_normal_selector.py#L1) · explicit_script_identifier
- [verification/v469_seam_crossedproduct_route.py](../../../../verification/v469_seam_crossedproduct_route.py#L1) · explicit_script_identifier
- [verification/v246_unification_data_crosscheck.py](../../../../verification/v246_unification_data_crosscheck.py#L1) · explicit_script_identifier
- [verification/v142_frame_integrality.py](../../../../verification/v142_frame_integrality.py#L1) · explicit_script_identifier
- [verification/v278_lsz_bridge_unitarity.py](../../../../verification/v278_lsz_bridge_unitarity.py#L1) · explicit_script_identifier
- [verification/v273_eg_gauge_running.py](../../../../verification/v273_eg_gauge_running.py#L1) · explicit_script_identifier
- [verification/v270_pmns_jarlskog_assembly.py](../../../../verification/v270_pmns_jarlskog_assembly.py#L1) · explicit_script_identifier
- [verification/v467_thirdgen_mixing_pattern.py](../../../../verification/v467_thirdgen_mixing_pattern.py#L1) · explicit_script_identifier
- [verification/v144_detratio_family.py](../../../../verification/v144_detratio_family.py#L1) · explicit_script_identifier
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier
- [verification/v371_ftransfer_pole.py](../../../../verification/v371_ftransfer_pole.py#L1) · explicit_script_identifier
- [verification/v271_eg_oneloop_quartic.py](../../../../verification/v271_eg_oneloop_quartic.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v375_observatory_registry.py](../../../../verification/v375_observatory_registry.py#L1) · explicit_script_identifier
- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier
- [verification/v86_nstar_reheating.py](../../../../verification/v86_nstar_reheating.py#L1) · explicit_script_identifier
- [FR.MPME.01](../../../../verification/status_ledger.csv#L49) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.TWISTOR.01](../../../../verification/status_ledger.csv#L584) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v138_vw_firewall.py](../../../../verification/v138_vw_firewall.py#L1) · explicit_script_identifier
- [verification/v269_spert_paqft_skeleton.py](../../../../verification/v269_spert_paqft_skeleton.py#L1) · explicit_script_identifier
- [verification/v145_pairing_atoms.py](../../../../verification/v145_pairing_atoms.py#L1) · explicit_script_identifier
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* v469 (ii): KLM mu-index + the 16-fold-way pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2937)
- [(* ---- (v133) zeta budget ---- Numerical continuations are Python-only (mpmath); the exact arithmetic of both routes is mirrored, plus a numeric spot-check in WL. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L995)
- [(* ---- (v142) frame integrality ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1157)
- [(* ==== v282: E8-at-tau=i unification, chi_E8(i)=12 (QGAMB.UNIFY.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2171)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ---- (v138) VW firewall ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1103)
- [(* ---- (v145) pairing atoms ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1227)
- [(* v469 (i): the LR locality integers of both extension routes *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2930)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ==== v469 round: the net-level crossed-product certification -- the LR locality integers (h_s(SO(16)_1) = 16/16 = 1 in Z; the mu4 glue h(J^k) = {1,1,1} all integer), the KLM mu-index arithmetic (4/2^2 = 1 and 16/4^2 = 1, both holomorphic) and the Kitaev 16-fold-way pin (nu = 16 x |C| = 16 = 0 mod 16, c_- = nu/2 = 8 = g_car+N_fam); mirrors v469 (the FHS Chern itself is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2924)
- [(* ---- (v118) hexagon-family dictionary ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L750)
- [(* ---- (v144) det-ratio family cancellation ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1208)
- [(* ---- (v86) reheating arithmetic (scale chain; pivot solve stays Python) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L84)
- [(* ==== v278: S_pert -> S_phys LSZ bridge, one-loop optical theorem (QFT4D.SPERT.04) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2142)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)
- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [(* ==== v271: concrete Epstein-Glaser one-loop quartic renormalization (QFT4D.SPERT.02) ==== v269 (S_pert skeleton) + v270 (PMNS Jarlskog, numerical) are Python-only; v271's exact EG core (scaling degree, counterterm count, loop factor) is mirrored here. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2083)
- [(* ---- (v139) selector triangle ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1110)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)
- [(* ---- (v136) dual-normal selector ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1051)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* ==== v273: EG one-loop gauge self-energy -> (b1,b2,b3) (QFT4D.SPERT.03) ==== v272 (nu-scale, numerical) + v274 (over-determination, numerical) + v275 (QG roadmap) are Python-only; v273's exact one-loop beta coefficients from the content are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2103)
- [(* ---- (v124) resummed clock ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L857)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEquivChain.lean#L1)
- [(* ---- (v122) margin theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L835)

## The four-bit quotient and the [formula in source] 15 x 16 census

[documents/note_e8_gaussian_code/sections/001-section-the-three-objects-label-sec-setup.tex:136-226](../../../../documents/note_e8_gaussian_code/sections/001-section-the-three-objects-label-sec-setup.tex#L136) · `section:080a8e1bc2a51ba0adbf` · section

All checks of this section are from gaussian_code_bridge_probe.py ( [formula in source] checks, [formula in source] failures, verdict gaussian-code-bridge-exact; exact integer and Fraction arithmetic throughout, sympy only for the Smith normal form), with the finite certificates re-proved by the Lean kernel in GaussianCodeBridge.lean.

Explicit source references:

- [TfptCarrier.GaussianCodeBridge.Q_unimodular](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L271) · explicit_lean_symbol
- [TfptCarrier.GaussianCodeBridge.index_sixteen](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L284) · explicit_lean_symbol
- [TfptCarrier.GaussianCodeBridge.mu4_acts_trivially_on_labels](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L880) · explicit_lean_symbol
- [TfptCarrier.GaussianCodeBridge.census_equidistribution](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L832) · explicit_lean_symbol
- [TfptCarrier.GaussianCodeBridge.census_mu4_stable](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L890) · explicit_lean_symbol
- [TfptCarrier.GaussianCodeBridge.census_zero_class_empty](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L827) · explicit_lean_symbol
- [experiments/tfpt-discovery/gaussian_code_bridge_probe.py](../../../../experiments/tfpt-discovery/gaussian_code_bridge_probe.py#L1) · explicit_file_citation
- [TfptCarrier.GaussianCodeBridge.smith_converse_ingredient](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L280) · explicit_lean_symbol
- [TfptCarrier.GaussianCodeBridge.snf_exponent_two](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L291) · explicit_lean_symbol
- [experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L1) · explicit_file_citation

## The two probes

[documents/note_e8_gaussian_code/sections/002-section-explicit-character-basis-of-the-eight-coxeter-ch.tex:255-294](../../../../documents/note_e8_gaussian_code/sections/002-section-explicit-character-basis-of-the-eight-coxeter-ch.tex#L255) · `section:0ee9c6150b29e6e55585` · subsection

@lccl@ probe & checks & verdict & arithmetic\\ gaussian_code_bridge_probe.py & [formula in source] & gaussian-code-bridge-exact & exact [formula in source] / Fraction; SNF via sympy\\ quartic_half_probe.py & [formula in source] & quartic-half-alive & exact [formula in source] / [formula in source] ; no floats\\ Both probes cross-validate their own machinery: the class map is computed twice (Hermite-normal-form labels

Explicit source references:

- [verification/v798_seam_clifford_modular_s.py](../../../../verification/v798_seam_clifford_modular_s.py#L1) · explicit_script_identifier
- [verification/v783_two_qubit_clifford.py](../../../../verification/v783_two_qubit_clifford.py#L1) · explicit_script_identifier
- [verification/v803_e8_jetcode_affine_nsr.py](../../../../verification/v803_e8_jetcode_affine_nsr.py#L1) · explicit_script_identifier
- [experiments/tfpt-discovery/quartic_half_probe.py](../../../../experiments/tfpt-discovery/quartic_half_probe.py#L1) · explicit_file_citation
- [experiments/tfpt-discovery/gaussian_code_bridge_probe.py](../../../../experiments/tfpt-discovery/gaussian_code_bridge_probe.py#L1) · explicit_file_citation

## Front matter / support

[documents/origin_theory/sections/005-section-consequences-and-the-whole-story-mc-label-sec-co.tex:1-37](../../../../documents/origin_theory/sections/005-section-consequences-and-the-whole-story-mc-label-sec-co.tex#L1) · `section:7a3a0ecd6a4836d0730c` · frontmatter

The coefficient [formula in source] is internally forced ( v73_k_c3_half.py). The dimensionless coefficient statement splits into two cutoff-independent pieces, neither a free normalisation: [formula in source] where (a) the factor [formula in source] is the universal action [formula in source] EOM factor ( [formula in source] action vs [formula in source] field equation, forced by [formula in source] ); and (b) [for

Explicit source references:

- [SEAM.THEOREM.01](../../../../verification/status_ledger.csv#L131) · explicit_claim_identifier · original status at linked ledger row
- [verification/v59_area_law_evidence.py](../../../../verification/v59_area_law_evidence.py#L1) · explicit_file_citation
- [verification/v68_seeley_dewitt_residual.py](../../../../verification/v68_seeley_dewitt_residual.py#L1) · explicit_script_identifier
- [verification/v73_k_c3_half.py](../../../../verification/v73_k_c3_half.py#L1) · explicit_file_citation
- [verification/v471_seam_horizon_replica.py](../../../../verification/v471_seam_horizon_replica.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v59) area-law evidence: exact seam-unit parts (the EE area-law is numerical, Python-only) ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L277)
- [(* ---- (v68) Seeley-DeWitt a2 for the carrier; 1/G UV-sensitive => k=c3/2 is a normalization ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L322)
- [(* ---- (v73) k=c3/2 forced: variational factor x Gauss-Bonnet topology; Fursaev-Solodukhin S=A/4 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L353)

## Grammar III — the action grammar (exponentials of [formula in source] )

[documents/tfpt_1_architecture_e8/sections/006-source.tex:168-228](../../../../documents/tfpt_1_architecture_e8/sections/006-source.tex#L168) · `section:4393abb43a4d948d7898` · section

The relation in plain words. The fine-structure constant is not just ``a number near [formula in source] '': it is the one exponential engine behind three wildly separated scales. Take the same [formula in source] and divide the carrier into it: [formula in source] of it sets the electroweak scale, [formula in source] it sets the Hubble scale, [formula in source] it sets the cosmological constant. The three action ch

Explicit source references:

- [verification/make_figures.py](../../../../verification/make_figures.py#L1) · explicit_file_citation

## The [formula in source] orders of the cosmological constant

[documents/tfpt_1_architecture_e8/sections/007-source.tex:315-339](../../../../documents/tfpt_1_architecture_e8/sections/007-source.tex#L315) · `section:d499db91c855166290b4` · subsection

The `` [formula in source] '' is the action [formula in source] ; the `` [formula in source] '' is the seam defect [formula in source] .

Explicit source references:

- [verification/v7_gravity_cosmo.py](../../../../verification/v7_gravity_cosmo.py#L1) · explicit_file_citation

## The ladder is binomial; [formula in source] singles out [formula in source]

[documents/tfpt_1_architecture_e8/sections/007-source.tex:340-360](../../../../documents/tfpt_1_architecture_e8/sections/007-source.tex#L340) · `section:b2783b0f143d528bcc39` · subsection

In units of the electroweak action [formula in source] , the three rungs are the first three binomial coefficients of the five-slot carrier: The first two are automatic ( [formula in source] , [formula in source] ). The non-trivial content is the cosmological-constant charge: since [formula in source] . So the ``double overlap'' charge [formula in source] equals the number of carrier-slot pairs [formula in source] on

Explicit source references:

- [verification/v2_carrier_pascal.py](../../../../verification/v2_carrier_pascal.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v2) carrier / Pascal ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L51)
- [experiments/lean4-carrier-rigidity/TfptCarrier/ArfSpinorCompiler.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/ArfSpinorCompiler.lean#L1)

## Discriminant action theorem: the ladder reads the hypercharge discriminant

[documents/tfpt_1_architecture_e8/sections/007-source.tex:361-376](../../../../documents/tfpt_1_architecture_e8/sections/007-source.tex#L361) · `section:3f07c167284f01e9956c` · subsection

The [formula in source] in the action ladder is not merely a rank: it is the square root of the discriminant of the carrier (hypercharge) polynomial. For a general split [formula in source] the carrier polynomial is [formula in source] with discriminant Hence the action ladder is, exactly, The electroweak hierarchy therefore reads the discriminant of the hypercharge polynomial, tying the Standard-Model hypercharge de

## A two-dimensional action lattice

[documents/tfpt_1_architecture_e8/sections/008-source.tex:1-28](../../../../documents/tfpt_1_architecture_e8/sections/008-source.tex#L1) · `section:4bd3b166d51412804107` · subsection

All scales sit on the integer lattice spanned by two generators — the gauge unit [formula in source] and the flavor unit [formula in source] — as [formula in source] : @llll@ scale & [formula in source] & [formula in source] & axis \\ [formula in source] & [formula in source] & [formula in source] & gauge \\ [formula in source] & [formula in source] & [formula in source] & gauge \\ [formula in source] & [formula in s

## Level 2.7 — the carrier Coxeter theorem (the strongest new bridge)

[documents/tfpt_1_architecture_e8/sections/008-source.tex:241-261](../../../../documents/tfpt_1_architecture_e8/sections/008-source.tex#L241) · `section:98d7438dbd43ec4fd921` · subsection

Use the universal root-system identity [formula in source] (Coxeter number [formula in source] ). With the carrier values [formula in source] and [formula in source] , The Coxeter number of [formula in source] is exactly three families times the cosmological pair action: each family contributes one [formula in source] -decuple ( [formula in source] ) to the Coxeter period. Equivalently the action ladder is a Coxeter

## Closures that the architecture rests on

[documents/tfpt_1_architecture_e8/sections/009-source.tex:135-179](../../../../documents/tfpt_1_architecture_e8/sections/009-source.tex#L135) · `section:ea5c1a7e1cb184a3e1e8` · section

Machine checks for this section: v36_ spectral_ action_ g2.py v4_ flavor_ matrix.py

Explicit source references:

- [verification/v36_spectral_action_g2.py](../../../../verification/v36_spectral_action_g2.py#L1) · explicit_file_citation
- [verification/v4_flavor_matrix.py](../../../../verification/v4_flavor_matrix.py#L1) · explicit_file_citation
- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)
- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [experiments/lean4-carrier-rigidity/TfptCarrier/YukawaPrimitive.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/YukawaPrimitive.lean#L1)
- [(* ---- (v36) G2: spectral-action Seeley-DeWitt coefficients -> R + R^2 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L113)

## The Coxeter–cyclotomic compiler and the explicit [formula in source] construction

[documents/tfpt_1_architecture_e8/sections/009-source.tex:328-365](../../../../documents/tfpt_1_architecture_e8/sections/009-source.tex#L328) · `section:b41c336f017dcb22d331` · part

This short companion compresses the TFPT bridge layer to a single discrete object. The claim is that the Standard-Model family, the cosmological action ladder, the [formula in source] counting data, the flavor budget, inflation and the cosmological constant are not independent coincidences but projections of one structure: the even-weight code [formula in source] on [formula in source] slots, together with its [formu

## Theorem 2 — Exterior–Action Duality

[documents/tfpt_1_architecture_e8/sections/010-source.tex:35-93](../../../../documents/tfpt_1_architecture_e8/sections/010-source.tex#L35) · `section:52f95907d935acc58bce` · section

tocsectionTheorem 2 On [formula in source] the flip [formula in source] is even, so [formula in source] . The three relations [formula in source] form a symmetric association scheme (the even Hamming scheme on five bits) with The eigenvalues (computed directly) are [formula in source] and [formula in source] , with [formula in source] . The two graphs are named objects: The Standard-Model family is the valency side [

Explicit source references:

- [verification/v79_review_identities.py](../../../../verification/v79_review_identities.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v79) Inverse Anchor Theorem: 1^T M^-1 1 = 1/atom; a^T M^-1 a = 1 for R,K,L ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L391)
- [(* ---- (v79) hypercharge Lucas-Binet D_n + the quark-ratio reading ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L386)

## Theorem 3 — the Coxeter closure equation selects [formula in source]

[documents/tfpt_1_architecture_e8/sections/010-source.tex:94-106](../../../../documents/tfpt_1_architecture_e8/sections/010-source.tex#L94) · `section:4b3c5060d4ff767121af` · section

tocsectionTheorem 3 Requiring the even family to equal the action truncation through pair degree, [formula in source] , i.e. the unique positive-integer solution. At [formula in source] the two sides equal [formula in source] , so carrier exhaustion [formula in source] and the Coxeter period [formula in source] coincide only at five slots: [formula in source] and [formula in source] .

## Master statement and status

[documents/tfpt_1_architecture_e8/sections/011-source.tex:317-366](../../../../documents/tfpt_1_architecture_e8/sections/011-source.tex#L317) · `section:2e4cf2c588758ec271e5` · section

tocsectionMaster statement Five carrier slots give [formula in source] ,\ four family corners give [formula in source] ,\ and [formula in source] gives the [formula in source] atlas. Equivalently: [formula in source] , [formula in source] , [formula in source] ; [formula in source] , [formula in source] , the [formula in source] Coxeter structure and the action ladder read the same five-slot pair-decuple. The Standar

## Shared origin of the seam normalization and carrier data

[documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex:2-10](../../../../documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex#L2) · `section:ee3f21621def60705393` · section

What is being explained. P1 supplies a primitive reflection-positive kernel on an oriented seam with unit winding. The number [formula in source] is its normalization. P2 declares the five-slot [formula in source] readout; its algebraic consequences and its conditional source reductions are already available. Calling both entries ``axioms'' records these interfaces and does not erase those reductions.

## The common explanation and its boundary.

[documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex:84-94](../../../../documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex#L84) · `section:415785fac2ef8712799f` · paragraph

The marked seam supplies the geometric setting; its cycles, sections, parabolic anchor and primitive normalization explain the shared [formula in source] data in their specified branches. P1, its scalar normalization and P2 are three levels of this construction. The operational original seed already contains the local net, state, reflection, time flow and collar generator. Existing conditional OS and common-action re

## Front matter / support

[documents/tfpt_2_standard_model/sections/001-part-the-standard-model-in-one-master-formula.tex:1-24](../../../../documents/tfpt_2_standard_model/sections/001-part-the-standard-model-in-one-master-formula.tex#L1) · `section:3b5318f0d65e3c120bab` · frontmatter

The Standard Model as one [formula in source] -ladder mass formula (all masses, CKM, PMNS, neutrinos), the flavor matrix from parabolic transport on [formula in source] , the five worked closures ( [formula in source] , quark [formula in source] , the explicit mass gap, Starobinsky [formula in source] , the H2 splitting), and gravity/quantum gravity as the seam response (the Seam Response Grammar [formula in source]

## (1) Flavor residual — ratios closed combinatorially; only the absolute scale is the wall-selection [formula in source] /

[documents/tfpt_2_standard_model/sections/006-section-neutrinos.tex:348-377](../../../../documents/tfpt_2_standard_model/sections/006-section-neutrinos.tex#L348) · `section:52b4cce0c7c3cdd6eeca` · subsection

Mehta–Seshadri supplies the equivalence-class framework — stable parabolic bundles [formula in source] irreducible unitary representations — but it does not by itself give ``parabolic-stable [formula in source] transport-minimal'' or the hexagonal word-lengths; that is the TFPT-specific step. The family connection [formula in source] provides the representation, its weights are the cusp exponents, the stable degree i

Explicit source references:

- [verification/v31_R_dictionary.py](../../../../verification/v31_R_dictionary.py#L1) · explicit_file_citation
- [verification/v34_h2_bridge_attempt.py](../../../../verification/v34_h2_bridge_attempt.py#L1) · explicit_file_citation
- [verification/v69_d4_q_geometry.py](../../../../verification/v69_d4_q_geometry.py#L1) · explicit_script_identifier
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier
- [verification/v33_explicit_flat_bundle.py](../../../../verification/v33_explicit_flat_bundle.py#L1) · explicit_file_citation
- [verification/v30_d4_character_variety.py](../../../../verification/v30_d4_character_variety.py#L1) · explicit_file_citation
- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier
- [verification/v70_q_integer_lift.py](../../../../verification/v70_q_integer_lift.py#L1) · explicit_script_identifier
- [verification/v27_wall_representative.py](../../../../verification/v27_wall_representative.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v69) D4-equivariant Q-geometry: Q_+ from cusp weights, Q_- E-block coupling ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L330)
- [(* ---- (v70) integer-lift of Q: Z4 puncture action unimodular, det Q = N_fam ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L337)
- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)

## (2) Gravity — the classical covariant field equation carries no free dimensionless Newton coupling under the named QFT/entanglement premises ; the [formula in source] inflation layer up to one truncation /

[documents/tfpt_2_standard_model/sections/007-subsection-2-gravity-the-classical-covariant-field-equat.tex:1-50](../../../../documents/tfpt_2_standard_model/sections/007-subsection-2-gravity-the-classical-covariant-field-equat.tex#L1) · `section:b6f6bed3bead3afd28d6` · subsection

The metric-sector field equation itself carries no free dimensionless Newton coupling under those premises: fixed- volume entanglement equilibrium (Jacobson 2015) forces the full covariant Einstein equation [formula in source] with both coefficients fixed ( [formula in source] , [formula in source] from [formula in source] ; v358/ v359). What this subsection details is the complementary spectral-action / inflation la

Explicit source references:

- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_script_identifier
- [verification/v36_spectral_action_g2.py](../../../../verification/v36_spectral_action_g2.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v36) G2: spectral-action Seeley-DeWitt coefficients -> R + R^2 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L113)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)

## (3) Strong CP — decoupled from the mass-gap problem /

[documents/tfpt_2_standard_model/sections/007-subsection-2-gravity-the-classical-covariant-field-equat.tex:51-91](../../../../documents/tfpt_2_standard_model/sections/007-subsection-2-gravity-the-classical-covariant-field-equat.tex#L51) · `section:e563b224f9d38e2a2c72` · subsection

The key reduction: [formula in source] needs neither the full OS package nor the mass gap. It follows from three purely structural facts plus reflection positivity on the admissible branch:

Explicit source references:

- [verification/v173_pfaffian_cp_car.py](../../../../verification/v173_pfaffian_cp_car.py#L1) · explicit_file_citation

## (4) Constructive QFT closure of the admissible physical sector /

[documents/tfpt_2_standard_model/sections/007-subsection-2-gravity-the-classical-covariant-field-equat.tex:92-179](../../../../documents/tfpt_2_standard_model/sections/007-subsection-2-gravity-the-classical-covariant-field-equat.tex#L92) · `section:a29d0e884b184665e5ef` · subsection

This is the constructive QFT closure of the admissible physical sector — not the unprojected ambient theory: the unprojected ambient metric measure is not constructed here and remains the [formula in source] frontier item ( ambient QG). On the admissible sector, the input previously listed as ``open'' — the parent-domination / mass-gap lemma from (3) — is not analytically open: it is computed by the explicit discrete

Explicit source references:

- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [QFT4D.EG.ALLORDER.01](../../../../verification/status_ledger.csv#L444) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier
- [QFT4D.OS.RECON.01](../../../../verification/status_ledger.csv#L745) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_file_citation
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_file_citation
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v374_ftransfer_qcd.py](../../../../verification/v374_ftransfer_qcd.py#L1) · explicit_script_identifier
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v252_full_finite_triple.py](../../../../verification/v252_full_finite_triple.py#L1) · explicit_file_citation
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_file_citation
- [verification/v381_qft4d_eg_allorder.py](../../../../verification/v381_qft4d_eg_allorder.py#L1) · explicit_file_citation
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)

## Conclusion

[documents/tfpt_2_standard_model/sections/007-subsection-2-gravity-the-classical-covariant-field-equat.tex:279-304](../../../../documents/tfpt_2_standard_model/sections/007-subsection-2-gravity-the-classical-covariant-field-equat.tex#L279) · `section:be245dd6a6eb1776e773` · section

tocsectionConclusion The simplification is substantial: the discrete Standard-Model structure (charges, families, flavor operators) is derived under the axioms; [formula in source] D existence, chiral measure, mirror gap, couplings and neutrino dynamics are open. The [formula in source] -ladder formula [formula in source] with one seed [formula in source] and one fixed residue matrix whose spectral invariants ( [form

Explicit source references:

- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_script_identifier
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier
- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)

## Exact scaffolding (proved)

[documents/tfpt_2_standard_model/sections/008-section-the-setup-parabolic-data-on-the-family-curve.tex:50-86](../../../../documents/tfpt_2_standard_model/sections/008-section-the-setup-parabolic-data-on-the-family-curve.tex#L50) · `section:c491ce4d24c1e7428721` · section

For a unitary (flat) family local system, [formula in source] is parabolic-polystable of parabolic degree zero (Mehta–Seshadri). With full flags and weights [formula in source] at each of the four punctures, The balanced (stable) underlying splitting is therefore [formula in source] .

## The explicit mass gap closes dynamics and clustering

[documents/tfpt_2_standard_model/sections/010-section-the-explicit-mass-gap-closes-dynamics-and-cluste.tex:1-160](../../../../documents/tfpt_2_standard_model/sections/010-section-the-explicit-mass-gap-closes-dynamics-and-cluste.tex#L1) · `section:cbcb5b486b59ced607c8` · section

unique limit state and exponential clustering] The admissible [formula in source] transport transfer operator has the closed spectrum [formula in source] (the cusp cubic). Hence [formula in source] The contraction ratio [formula in source] satisfies the Dobrushin uniqueness condition for a finite-range (range- [formula in source] , cyclic) transfer operator; therefore the thermodynamic limit has a unique Gibbs/vacuum

Explicit source references:

- [verification/v798_seam_clifford_modular_s.py](../../../../verification/v798_seam_clifford_modular_s.py#L1) · explicit_file_citation
- [verification/v783_two_qubit_clifford.py](../../../../verification/v783_two_qubit_clifford.py#L1) · explicit_file_citation
- [verification/v802_gnet_martingale.py](../../../../verification/v802_gnet_martingale.py#L1) · explicit_file_citation
- [GNET.LOCALNET.01](../../../../verification/status_ledger.csv#L675) · explicit_claim_identifier · original status at linked ledger row
- [verification/v746_phys_gnet_local_functor.py](../../../../verification/v746_phys_gnet_local_functor.py#L1) · explicit_file_citation
- [verification/v779_gnet_gns_arf.py](../../../../verification/v779_gnet_gns_arf.py#L1) · explicit_file_citation

## The seam elliptic complex is Fredholm (APS)

[documents/tfpt_2_standard_model/sections/010-section-the-explicit-mass-gap-closes-dynamics-and-cluste.tex:161-183](../../../../documents/tfpt_2_standard_model/sections/010-section-the-explicit-mass-gap-closes-dynamics-and-cluste.tex#L161) · `section:0a3610a9e32985012149` · section

Compactify the oriented two-dimensional normal slice [formula in source] to [formula in source] . The geometric BRST/Hodge differential, restricted to the seam-even sector with Atiyah–Patodi–Singer boundary conditions on [formula in source] , is an elliptic operator on a compact manifold; hence it is Fredholm, with finite-dimensional kernel and cokernel and an integer index. The compact bosonic Higgs index selects [f

## The [formula in source] truncation is [formula in source] -bounded

[documents/tfpt_2_standard_model/sections/010-section-the-explicit-mass-gap-closes-dynamics-and-cluste.tex:184-272](../../../../documents/tfpt_2_standard_model/sections/010-section-the-explicit-mass-gap-closes-dynamics-and-cluste.tex#L184) · `section:ee5c68ec39d878e448bd` · section

On the [formula in source] attractor the curvature during observable inflation is [formula in source] . Higher curvature terms [formula in source] , [formula in source] , enter with [formula in source] , so their contribution to the observables is [formula in source] . With [formula in source] this is [formula in source] ( [formula in source] ), and the leading predictions are [formula in source] matching Planck. The

Explicit source references:

- [verification/v86_nstar_reheating.py](../../../../verification/v86_nstar_reheating.py#L1) · explicit_file_citation
- [verification/v494_cosmo_running_mu_record.py](../../../../verification/v494_cosmo_running_mu_record.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v86) reheating arithmetic (scale chain; pivot solve stays Python) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L84)

## The ambient closure: one spectral-action measure, IR-complete by decoupling

[documents/tfpt_2_standard_model/sections/010-section-the-explicit-mass-gap-closes-dynamics-and-cluste.tex:273-408](../../../../documents/tfpt_2_standard_model/sections/010-section-the-explicit-mass-gap-closes-dynamics-and-cluste.tex#L273) · `section:9d6c7d63f5df5a4df579` · section

The two formerly-hardest items — the ambient (unprojected) QFT measure beyond the selector [formula in source] , and full quantum gravity beyond the [formula in source] scalaron — are not two separate problems. They are one object: the relative spectral-action measure [formula in source] with [formula in source] the matter/gauge Dirac operator and the gravitational data in the same trace. Three facts reduce both open

Explicit source references:

- [verification/v28_gravity_fR.py](../../../../verification/v28_gravity_fR.py#L1) · explicit_file_citation
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v255_spectral_action_expansion.py](../../../../verification/v255_spectral_action_expansion.py#L1) · explicit_file_citation
- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_script_identifier
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_file_citation
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_file_citation
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_file_citation
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)

## Status: what is now written out, and what is honestly left

[documents/tfpt_2_standard_model/sections/011-subsection-seam-response-grammar-gravity-is-the-inverse-.tex:307-360](../../../../documents/tfpt_2_standard_model/sections/011-subsection-seam-response-grammar-gravity-is-the-inverse-.tex#L307) · `section:00b05ea692eac25e3d3e` · section

tocsectionStatus @p5.6cmcY@ Item & Status & Result \\ Dynamics / clustering (Lemma Gravity Fredholm (Lemma [formula in source] truncation (Lemma Starobinsky [formula in source] / [formula in source] & / & [formula in source] GeV exact ; [formula in source] a \ band over [formula in source] \\ [formula in source] typing (Prop. [formula in source] (Prop. quark [formula in source] -rationals & / & hierarchy + source rat

Explicit source references:

- [verification/v723_phys_modular_clock.py](../../../../verification/v723_phys_modular_clock.py#L1) · explicit_file_citation
- [verification/v19_monodromy_moduli.py](../../../../verification/v19_monodromy_moduli.py#L1) · explicit_file_citation
- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/v375_observatory_registry.py](../../../../verification/v375_observatory_registry.py#L1) · explicit_script_identifier
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier
- [verification/v371_ftransfer_pole.py](../../../../verification/v371_ftransfer_pole.py#L1) · explicit_script_identifier
- [verification/v724_phys_t3b_modular_flows.py](../../../../verification/v724_phys_t3b_modular_flows.py#L1) · explicit_file_citation
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier

## Front matter / support

[documents/tfpt_3_e8_audit_bootstrap/sections/004-source.tex:1-373](../../../../documents/tfpt_3_e8_audit_bootstrap/sections/004-source.tex#L1) · `section:3d2a5b53f3345969fcbc` · frontmatter

The interacting Klein–Landau semigroup ( WOIT.GAMMA.SEMI.01). The gamma net-dynamics slice ( [formula in source] checks, [formula in source] \,s). At [formula in source] the v524 Klein–Landau pattern reproduces exactly at [formula in source] on the admissible cut: all seven steps [formula in source] Hermitian (domain dims [formula in source] ), even steps PD, odd steps [formula in source] indefinite with exactly zero

Explicit source references:

- [QGEO.ORBCAS.01](../../../../verification/status_ledger.csv#L642) · explicit_claim_identifier · original status at linked ledger row
- [verification/v605_translator.py](../../../../verification/v605_translator.py#L1) · explicit_script_identifier
- [verification/v645_klein_rp.py](../../../../verification/v645_klein_rp.py#L1) · explicit_file_citation
- [verification/v624_external_lattice_audit.py](../../../../verification/v624_external_lattice_audit.py#L1) · explicit_file_citation
- [verification/v623_covered_seam.py](../../../../verification/v623_covered_seam.py#L1) · explicit_file_citation
- [verification/v599_real_structure.py](../../../../verification/v599_real_structure.py#L1) · explicit_script_identifier
- [verification/v425_dyn_transfer_universal.py](../../../../verification/v425_dyn_transfer_universal.py#L1) · explicit_script_identifier
- [verification/v534_seam_straddle_cone.py](../../../../verification/v534_seam_straddle_cone.py#L1) · explicit_script_identifier
- [verification/v639_twist_ope.py](../../../../verification/v639_twist_ope.py#L1) · explicit_file_citation
- [verification/v629_root_incidence.py](../../../../verification/v629_root_incidence.py#L1) · explicit_file_citation
- [CONTRACT.F.01](../../../../verification/status_ledger.csv#L302) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.SEAMID.01](../../../../verification/status_ledger.csv#L638) · explicit_claim_identifier · original status at linked ledger row
- [verification/v622_seam_identification.py](../../../../verification/v622_seam_identification.py#L1) · explicit_file_citation
- [E8.RM13REV.01](../../../../verification/status_ledger.csv#L650) · explicit_claim_identifier · original status at linked ledger row
- [verification/v628_orbifold_casimir.py](../../../../verification/v628_orbifold_casimir.py#L1) · explicit_file_citation
- [verification/v632_ftransfer_pgl2.py](../../../../verification/v632_ftransfer_pgl2.py#L1) · explicit_file_citation
- [verification/v597_mu4_cover_model.py](../../../../verification/v597_mu4_cover_model.py#L1) · explicit_script_identifier
- [verification/v533_ftransfer_disc7_norm.py](../../../../verification/v533_ftransfer_disc7_norm.py#L1) · explicit_script_identifier
- [E8.ST31.01](../../../../verification/status_ledger.csv#L646) · explicit_claim_identifier · original status at linked ledger row
- [FTR.PGL2.01](../../../../verification/status_ledger.csv#L644) · explicit_claim_identifier · original status at linked ledger row
- [QGEO.KLEINRP.01](../../../../verification/status_ledger.csv#L649) · explicit_claim_identifier · original status at linked ledger row
- [verification/v183_koide_f_corner_transfer.py](../../../../verification/v183_koide_f_corner_transfer.py#L1) · explicit_script_identifier
- [WOIT.GAMMA.SEMI.01](../../../../verification/status_ledger.csv#L637) · explicit_claim_identifier · original status at linked ledger row
- [verification/v805_e8_syndrome_algebra.py](../../../../verification/v805_e8_syndrome_algebra.py#L1) · explicit_file_citation
- [verification/v633_orbit60_quotient.py](../../../../verification/v633_orbit60_quotient.py#L1) · explicit_file_citation
- [verification/v638_code_semantics.py](../../../../verification/v638_code_semantics.py#L1) · explicit_file_citation
- [experiments/lean4-carrier-rigidity/TfptCarrier/G31Orders.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/G31Orders.lean#L1) · explicit_file_citation
- [QGEO.TWISTOPE.01](../../../../verification/status_ledger.csv#L648) · explicit_claim_identifier · original status at linked ledger row
- [verification/v626_e8_code.py](../../../../verification/v626_e8_code.py#L1) · explicit_file_citation
- [verification/v510_seam_bit_freedom.py](../../../../verification/v510_seam_bit_freedom.py#L1) · explicit_script_identifier
- [E8.CODE.01](../../../../verification/status_ledger.csv#L641) · explicit_claim_identifier · original status at linked ledger row
- [verification/v621_interacting_semigroup.py](../../../../verification/v621_interacting_semigroup.py#L1) · explicit_file_citation
- [E8.INCIDENCE.01](../../../../verification/status_ledger.csv#L643) · explicit_claim_identifier · original status at linked ledger row
- [verification/v602_duality_forms.py](../../../../verification/v602_duality_forms.py#L1) · explicit_script_identifier
- [verification/v634_st31_structure.py](../../../../verification/v634_st31_structure.py#L1) · explicit_file_citation
- [verification/v575_cr_continuous_uniflow.py](../../../../verification/v575_cr_continuous_uniflow.py#L1) · explicit_script_identifier
- [AX.P2.01](../../../../verification/status_ledger.csv#L3) · explicit_claim_identifier · original status at linked ledger row
- [verification/v614_transport_deck.py](../../../../verification/v614_transport_deck.py#L1) · explicit_script_identifier
- [verification/v524_woit_beta2_os_quotient.py](../../../../verification/v524_woit_beta2_os_quotient.py#L1) · explicit_script_identifier
- [QGEO.COVERSEAM.01](../../../../verification/status_ledger.csv#L639) · explicit_claim_identifier · original status at linked ledger row
- [EXTREV.LATTICE.01](../../../../verification/status_ledger.csv#L640) · explicit_claim_identifier · original status at linked ledger row
- [verification/v578_native_jets.py](../../../../verification/v578_native_jets.py#L1) · explicit_script_identifier
- [verification/v519_woit_theta_rp_free.py](../../../../verification/v519_woit_theta_rp_free.py#L1) · explicit_script_identifier
- [E8.CODESEM.01](../../../../verification/status_ledger.csv#L647) · explicit_claim_identifier · original status at linked ledger row
- [verification/v617_seam_cover_forcing.py](../../../../verification/v617_seam_cover_forcing.py#L1) · explicit_script_identifier
- [verification/v646_rm13_reverse.py](../../../../verification/v646_rm13_reverse.py#L1) · explicit_file_citation
- [E8.ZIQUOTIENT.01](../../../../verification/status_ledger.csv#L645) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v519 round: WOIT.THETA.FREE.01 -- the WOIT-alpha milestone of the OS twistor bridge: the real structure exists (two-family classification, mu4 torsor, Z8/C^4/Cl(16) levels) and the exact free-RP preconditions + site-cut degeneration + continuum control. The 40-digit positive- definiteness inertia certificates of the bond-cut Grams (R3/R5) are Python-only by the suite convention (mpmath eigh at dps 40). Mirrors v519. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L7078)
- [(* ==== v533 round: FTR.DISC7.NORM.01 -- the disc -7 norm line. Exact mirrors: norm identity, translation, norms (2,4,14,32), constant discriminant -7, vertex, negative controls. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8373)
- [(* ==== v524 round: WOIT.BETA2.OS.01 -- "the OS quotient of the free system made explicit" (verdict SUCCESS, [C]-typed per contract precision (iii)). Four exact mirrors: (i) the euclidean-rotation preconditions (alpha_1^N = (-1)^F wrap arithmetic at N = 16 and N = 8) + the exact PD certificate identity sin^2(3pi/8) - sin(pi/8) sin(5pi/8) = 1/2; (ii) the Klein-Landau one-particle transfer forms tau_k on the shrinking domains at N = 16, k = 1..7: exactly Hermitian, odd k with exactly ZERO diagonal (the C(even) = 0 checkerboard), even k = 2j obeying the exact square identity tau_2j(a, b) = <alpha_j a, alpha_j b> entrywise; (iii) the clock identities (C(9) = C(7), C(11) = C(5), the transparency product formula sin(7pi/16) - sin(5pi/16) = 2 cos(3pi/8) sin(pi/16) > 0) + the N = 8 compressed clock spectrum EXACTLY {1, sqrt(2)-1} = {1, 1/delta_S} via the generalized eigenvalues of (tau_2, G) on the one-particle domain, with sin(pi/8)/sin(3pi/8) = sqrt(2)-1; (iv) the pre-declared KMS witnesses (C(1)/C(3) = 1 + sqrt(2) = delta_S at N = 8; -(C5-C7)((C3-C7)+(C5-C7)) < 0 exactly at N = 16, i.e. det(G - tau_4) < 0 on the odd one-particle clock block). The 37x37/16x16 H_phys Grams, the 40-digit inertia certificates, the spectral projections and the reconstructed rotation group are numerical (mpmath eighe), Python-only, flagged here. Mirrors v524. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L7966)
- [experiments/lean4-carrier-rigidity/TfptCarrier/HammingCode.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/HammingCode.lean#L1)
- [(* ==== v510 round: SEAM.BIT.FREEDOM.01 -- the freedom attack on the v506/v507 alignment bit: the POSITION half of the bit is topology. Exact content mirrored: (i) the DIHEDRAL CENSUS -- on the refined 2N-point seam circle D_16 has 17 involutions with EXACTLY ONE free (the antipodal shift by 8) and all 16 reflections have exactly 2 circle fixed points (8 site-axis + 8 edge-axis; robustness N = 8, 12; odd control N = 9: zero free), plus the quotient topology chi(C/antipode) = 0 (circle) vs chi(C/reflection) = 1 (interval, 2 boundary corners); (ii) the AUT COMPLETENESS -- Aut(C_8) = D_8 by brute force over all 8! = 40320 permutations and Aut(C_16) = D_16 by degree-2 propagation (all 32 candidates), so the census covers ALL circle symmetries; (iii) the Cl(16) DICHOTOMY CENSUS -- on explicit 256x256 Jordan-Wigner gammas the free deck's NS implementer squares to 256 gamma_1...gamma_16 = 256 (-1)^F (nonsplit) while the implementers of ALL 16 reflections exist (one of the two spin lifts +-R) and square to SCALARS {2^7 (site-axis), 2^8 (edge-axis)} (split): NONSPLIT <=> FREE over the complete census; root dichotomy 2 (the +-quarter shifts) vs 0 total over all 16 reflections; (iv) the REAL-STRUCTURE TABLE -- the deck-commuting antiholomorphic involutions solve symbolically to exactly two families (M diagonal: mu conj z, fix line through the poles; M antidiagonal: mu/conj z, fix circle |z| = sqrt mu), and the unique mark-fixing real structure per configuration is equatorial/FREE for mu4 and generic-unit, real-axis/NOT-free for silver and generic-real, and NONEXISTENT for the hexagonal frame; (v) HARMONIC + FREE => mu4 -- the deck-pair cross-ratio (1-b)^2/(1+b)^2 = -1 solves exactly to b = +-i, with the lambda table (2, 1/2, 3/4, 5/4) real and the hexagonal lambda complex; (vi) the HONEST COUNTERWITNESSES -- the discrete witness {0,1,8,9} (shift-invariant, crossing, but not quarter-shift-invariant) and the continuous witness {+-1, +-(3+4i)/5} (pair cross-ratio -1/4, j = 148176/25 != 1728, clock does not preserve): the Fock dichotomy measures the CLASS, not the modulus. The Moebius/Klein covering models and the pillowcase branch control stay Python-side (v510). Mirrors v510. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5155)
- [experiments/lean4-carrier-rigidity/TfptCarrier/OrbitPacket.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/OrbitPacket.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CoverEmbedding.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CoverEmbedding.lean#L1)
- [(* ==== v534 round: SEAM.STRADDLE.CONE.01 -- the RP cone on the interacting FK quartets (the first dynamical selection of the alignment bit). Exact mirrors: the stabilizer census of the mark family, the pi/2 partition + straddle incidence, the reflection-closure signs of the straddling quartets, and the equivariant Fix rays (uniform on m=2..6, NS-wrap-twisted on m=1,7). The 256-dim Fock spectra, the PT Gram corrections and the finite-g survival scan are Python-only (numpy), flagged per GATE.WOLFRAM.02. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8402)
- [(* ---- (v183) Koide 53/54 = a^T(R+Q)1 / (2 1^T R a): missing Sheet-Diamond corner ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1547)
- [experiments/lean4-carrier-rigidity/TfptCarrier/LorentzCongruence.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/LorentzCongruence.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/G31WordOrders.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/G31WordOrders.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SquareParity.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SquareParity.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/G31Orders.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/G31Orders.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GaussianCodeBridge.lean#L1)

## The celestial and twistor continuum route

[documents/tfpt_3_e8_audit_bootstrap/sections/011-section-the-celestial-and-twistor-continuum-route-label-.tex:1-363](../../../../documents/tfpt_3_e8_audit_bootstrap/sections/011-section-the-celestial-and-twistor-continuum-route-label-.tex#L1) · `section:b63f8484d81bde8026ec` · section

The previous section identified the seam's boundary [formula in source] -sphere with the celestial sphere of null directions ( [formula in source] , clock [formula in source] order- [formula in source] M\"obius map v180_clock_is_mobius.py). Taking that identification seriously and reading the seam through the celestial-holography lens (the celestial-chiral-algebra construction of Bittleston–Homans–Sharma on ALE space

Explicit source references:

- [verification/v449_seam_mmst_uniform.py](../../../../verification/v449_seam_mmst_uniform.py#L1) · explicit_script_identifier
- [verification/v504_celestial_wp5d_klm_completeness.py](../../../../verification/v504_celestial_wp5d_klm_completeness.py#L1) · explicit_file_citation
- [WOIT.OS.TWISTOR.01](../../../../verification/status_ledger.csv#L585) · explicit_claim_identifier · original status at linked ledger row
- [verification/v511_celestial_wp5e_delta2_tensor_ledger.py](../../../../verification/v511_celestial_wp5e_delta2_tensor_ledger.py#L1) · explicit_file_citation
- [verification/v499_p2_weights_deligne_bg.py](../../../../verification/v499_p2_weights_deligne_bg.py#L1) · explicit_script_identifier
- [verification/v180_clock_is_mobius.py](../../../../verification/v180_clock_is_mobius.py#L1) · explicit_file_citation
- [SEAM.CLOCK.RIGIDITY.01](../../../../verification/status_ledger.csv#L570) · explicit_claim_identifier · original status at linked ledger row
- [verification/v503_qgeo_emergence_light.py](../../../../verification/v503_qgeo_emergence_light.py#L1) · explicit_file_citation
- [verification/v490_seam_parity_census.py](../../../../verification/v490_seam_parity_census.py#L1) · explicit_script_identifier
- [CELEST.DTERM.NONDERIV.01](../../../../verification/status_ledger.csv#L577) · explicit_claim_identifier · original status at linked ledger row
- [verification/v500_celestial_wp5c_gns_limit_state.py](../../../../verification/v500_celestial_wp5c_gns_limit_state.py#L1) · explicit_file_citation
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [verification/v497_celestial_wp5a_null_ideal.py](../../../../verification/v497_celestial_wp5a_null_ideal.py#L1) · explicit_file_citation
- [verification/v510_seam_bit_freedom.py](../../../../verification/v510_seam_bit_freedom.py#L1) · explicit_file_citation
- [verification/v214_seam_pillowcase.py](../../../../verification/v214_seam_pillowcase.py#L1) · explicit_script_identifier
- [verification/v498_celestial_wp5b_singular_vector.py](../../../../verification/v498_celestial_wp5b_singular_vector.py#L1) · explicit_file_citation
- [verification/v495_celestial_wp3_gs_alignment.py](../../../../verification/v495_celestial_wp3_gs_alignment.py#L1) · explicit_file_citation
- [QGEO.REALIZE.01](../../../../verification/status_ledger.csv#L253) · explicit_claim_identifier · original status at linked ledger row
- [verification/v168_qgeo_rigidity.py](../../../../verification/v168_qgeo_rigidity.py#L1) · explicit_script_identifier
- [verification/v509_celestial_wp5e_eps2_level_from_flux.py](../../../../verification/v509_celestial_wp5e_eps2_level_from_flux.py#L1) · explicit_file_citation
- [verification/v502_celestial_wp5e_prefactor_level.py](../../../../verification/v502_celestial_wp5e_prefactor_level.py#L1) · explicit_file_citation
- [verification/v505_celestial_wp5e_beta_equivariant_ledger.py](../../../../verification/v505_celestial_wp5e_beta_equivariant_ledger.py#L1) · explicit_file_citation
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [verification/v506_seam_clock_rigidity.py](../../../../verification/v506_seam_clock_rigidity.py#L1) · explicit_file_citation
- [verification/v507_seam_bit_origin.py](../../../../verification/v507_seam_bit_origin.py#L1) · explicit_file_citation
- [verification/v518_celestial_wp5e_delta1_derived_kill.py](../../../../verification/v518_celestial_wp5e_delta1_derived_kill.py#L1) · explicit_script_identifier
- [SEAM.TAU.FLAG.01](../../../../verification/status_ledger.csv#L576) · explicit_claim_identifier · original status at linked ledger row
- [verification/v377_seam_s3_e8character.py](../../../../verification/v377_seam_s3_e8character.py#L1) · explicit_script_identifier
- [verification/v508_celestial_wp5e_gamma_pairing.py](../../../../verification/v508_celestial_wp5e_gamma_pairing.py#L1) · explicit_file_citation
- [verification/v512_seam_tau_flag.py](../../../../verification/v512_seam_tau_flag.py#L1) · explicit_file_citation
- [verification/v215_seam_deck_killtest.py](../../../../verification/v215_seam_deck_killtest.py#L1) · explicit_script_identifier
- [verification/v496_celestial_wp4_character_shadow.py](../../../../verification/v496_celestial_wp4_character_shadow.py#L1) · explicit_file_citation
- [verification/v513_celestial_dterm_nonderivation.py](../../../../verification/v513_celestial_dterm_nonderivation.py#L1) · explicit_file_citation
- [verification/v514_celestial_wp5e_eps1_axion_slot.py](../../../../verification/v514_celestial_wp5e_eps1_axion_slot.py#L1) · explicit_file_citation
- [SEAM.BIT.FREEDOM.01](../../../../verification/status_ledger.csv#L572) · explicit_claim_identifier · original status at linked ledger row
- [verification/v493_celestial_wp2_clock_deformation.py](../../../../verification/v493_celestial_wp2_clock_deformation.py#L1) · explicit_file_citation
- [verification/v492_celestial_z4_orbifold.py](../../../../verification/v492_celestial_z4_orbifold.py#L1) · explicit_file_citation
- [verification/v128_graded_hull.py](../../../../verification/v128_graded_hull.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [CELEST.SEAM.01](../../../../verification/status_ledger.csv#L556) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.BIT.ORIGIN.01](../../../../verification/status_ledger.csv#L571) · explicit_claim_identifier · original status at linked ledger row
- [verification/v501_celestial_wp5d_two_interval_index.py](../../../../verification/v501_celestial_wp5d_two_interval_index.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [(* ==== v511 round: CELEST.WP5E.DELTA2.01 -- WP5e-delta-2 "the full-tensor ledger" (the named escape route of the v508 exchange no-go, executed: the collapse is CONFIRMED full-tensorially, and one cubic door opens). Exact content mirrored: (i) the INNERNESS COUNT -- 240 glue roots with class split (52, 64, 60, 64), h = (2,2,2,2,2;0^4) and h' = (0^5;1,1,1,-3) reading the class mod 4 on ALL 240 roots, the class-0 roots splitting 40 pure D5 + 12 pure A3 with NO mixed root (g0 = d5 (+) a3 semisimple, NO u(1); root span rank 8), and class(-alpha) = -class(alpha); (ii) the BILINEAR HOM TABLE -- dim Hom_{g0}(g_j (x) g_j') by exact Kostant alternating sums on integer-scaled minuscule weight orbits: nonzero ONLY at j' = -j with dims (2, 1, 1, 1), all 12 charged pairs ZERO (the inner-torus theorem at arity 2); (iii) the TRILINEAR HOM TABLE -- all 15 non-neutral sector triples have Hom = 0 (arity 3) and the five neutral multisets carry dims (3, 2, 2, 1, 1); (iv) the UNIQUE SYMMETRIC SURVIVOR -- Inv(S^3 adj) = (0, 1) for (d5, a3) (so(10) has NO cubic Casimir, su(4) has exactly one), Inv(L^3 adj) = (1, 1), Inv(S^2 adj) = (1, 1), and on g0: Inv(S^3 g0) = 1, Inv(L^3 g0) = 2 (of the three g0-trilinears exactly ONE is totally symmetric: the su(4) d-symbol), with the dual machinery guard 4 x 4 = 0 vs 4 x 4bar = 1; (v) Q_dd BY TWO ROUTES -- explicit su(4) matrices: bracket vertices f(x,x,B) = 0, root-direction contractions d(x,x,E_ij) = 0, tadpole sum G^{-1}_ab d(x,B_a,B_b) = 0, and the exchange quartic Q_dd = sum d(x,x,B_a)(60 G)^{-1}_ab d(x,x,B_b) = (0, 0, -1/240, 0, 1/60) = (1/60)(T3 - P3/4) agreeing exactly with the hand route |x^2 - (tr x^2/4) 1|^2/60, so Phi_T3(Q_dd) = 1/60 != 0; plus Tr X^3 = 0 for antisymmetric X (no T5 analogue); (vi) the PSI CERTIFICATE + THE RELAXED SOLUTION -- E13 = (16, -96, 144, 0, 0), E22 = (16, 32, 16, 0, 0), E00 = 225 E22, A_fix = (9, -30, -15, 0, 32); M = [E13, E22, E00, Q_dd] has rank 3 with rank([M | A_fix]) = 4 (still unsolvable), psi = Phi_P - Phi_T3/4 annihilates all four channels with psi(A_fix) = 72 - 8 = 64, and the relaxed system {u, v, w, Q_dd} solves EXACTLY: A_fix = -u + 8v + 2w + 1920 Q_dd (rank 4, c_d = 1920 = |W(D5)| = 8 x 240); (vii) the CONTROLS -- SO(16)/D8: bilinears (adj adj, adj spin, spin spin) = (1, 0, 1), trilinears (1, 0, 1, 0), Inv(S^3 adj_d8) = 0 and Inv(L^3 adj_d8) = 1 (NO symmetric cubic, no odd sectors, no A3 block -- the delta-2 channel cannot even be posed); false g0 = d5+a2+u(1): bilinear entries (5, 2, 2) != (2, 1, 1) and Inv(S^3 g0_false) = 6 != 1; false sector content (16_s, 4bar): both pairings die (0, 0) vs the true value 1. Mirrors v511. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5380)
- [(* ---- (v214) QGEO.PILLOW.01: pillowcase reduction -- cross-ratio 2 => j=1728 => order-4 CM ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1727)
- [(* ==== v495 round: CELEST.WP3.01 -- the Green-Schwarz axion coefficient lambda_g (Costello) vs c3 = 1/(8 pi) (WP3 of CELEST.SEAM.01, alignment test only). Exact content mirrored: the Deligne closed form lambda~^2 = 10 h_vee^2/(dim+2) = h_vee + 6 for all 8 algebras on Costello's list, with the dimension formula dim = 2(5h_vee-6)(h_vee+1)/(h_vee+6) and the fundamental-index values (h_vee+6)/i^2; the E8 alignment (kappa/c3)^2 = lambda^2/3 = 12 = |mu4| N_fam resp. 1/300 with kappa/c3 itself irrational; the so8 factor-2 slip (exact 3, printed 3/2); the look-elsewhere pass counts (smooth 8/8, lambda-integrality 2/8 = {sl3, e8}, g_car | h_vee 1/8 = {e8}, phi(h_vee) = rank 4/8, phi(h) = rank 5/8) and the e8 exponents = totatives of 30; mirrors v495. The Okubo root-sum identity itself is polynomial (verified in Python from explicit root systems); its arithmetic consequence lambda~^2 = h_vee + 6 is mirrored here. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3210)
- [(* ==== v493 round: CELEST.WP2.01 -- the clock-invariant deformation XY = Z^4 + a0 of the A3 seam orbifold (WP2 of CELEST.SEAM.01). Exact content mirrored: clock-invariance selects the 1-parameter slice sharply (P(iZ) = P(Z) forces a3 = a2 = a1 = 0; prod(Z - i^k w) = Z^4 - w^4); disc = 256 a0^3; binary-quartic invariants I = 12 a0, J = 0 identically => j = 1728 frozen (tau = i pillowcase), cross-ratio 2; negative controls a1 (j = 0) and a2-instance (j = 1556068/81); clock on H2 = Coxeter element of W(A3) (char x^3+x^2+x+1, eigenvalues {i,-1,-i}, no invariant cycle, chi_1-Fourier eigenline, period point t(i-1)(1,i,i^2)); versal weights (1,2,3,0) = regular mu4 rep; S-algebra corrections linear in a0 at the Z4 wrap (6/16 pairs) with the EH k = 2 analogue; mirrors v493. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3148)
- [(* ==== v509 round: CELEST.WP5E.EPS2.01 -- WP5e-epsilon-2 "the CPS level-from-flux dial" (verdict B; k = 1 from flux quantisation on the lockstep spheres, Costello-Paquette-Sharma Burns holography). Exact content mirrored: (i) the CPS PERIOD INTEGRALS -- the S^3 period of the brane backreaction 3-form = (2 pi i)^2 (pullback density -sin 2 theta) and the exceptional-sphere Kahler flux = 2 pi N exactly (flat part dies at t -> 0); (ii) the LEVEL ~ FLUX contraction 2N (symplectic bosons, N = 1, 2, 3) with the Dynkin normalisation T(so8 vector) = 1, T(adj) = 6; (iii) the PAIRING-MATRIX ENUMERATION -- row candidates (8, 6, 8), 48 unimodular solutions of P C^-1 P^T = C^-1, with effectivity P >= 0 exactly 2 = {identity, A3 diagram flip}; (iv) the FLUXES (64, 60, 64) = dim g_i with total 188 = 248 - 60, flip-invariant, one quantum per charged current; (v) the SECTOR-COUNTER dial #primaries((E8)_k) = (1, 3, 5, 10, 15, 27) for k = 1..6 from the affine marks -- exactly one sector iff k = 1; (vi) the LOCKSTEP THEOREM -- clock invariance forces a3 = a2 = a1 = 0, periods t(i-1)(1, i, i^2) with squared moduli 2t^2 equal on all three spheres, monodromy = pure phase i, det(A-1) = -4, and the negative probe Z^4 - Z (disc -27): 0/24 orderings lockstep; (vii) the NEGATIVE CONTROLS -- SO(16) fluxes (0, 60, 0) with ch2 defect -30 and condensation stuck at 4, A2 form det 3, glue diagonal 23/24 not integer, mu = 12 not a perfect square, Coxeter order 3 != 4. Mirrors v509. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4981)
- [(* ==== v508 round: CELEST.WP5E.GAMMA.01 -- WP5e-gamma "the sphere-axion pairing check", an honest rigid NEGATIVE result (the exchange route for pairing the three twisted sphere axions with the rigid 32 T3 residual of the v505 Atiyah-Bott ledger is killed with a certificate). Exact content mirrored: (i) the VERTEX-SPACE DIMS -- the W(D5) x W(A3)-invariant quadratics on the 8-dim glue Cartan have bidegree dims (1,0,1), total 2 = span{s5, s3}, and the invariant quartics have bidegree dims (2,0,1,0,2), total 5 = span{P1,P2,P3,T5,T3} (Weyl nullspace arithmetic, generator by generator); (ii) the TWISTED QUADRATICS + the K^(0) = -15 K^(2) collapse (the two even-sector quadratics are PARALLEL); (iii) the PRODUCT THEOREM -- (a s5 + b s3)(c s5 + d s3) = ac P1 + (ad+bc) P2 + bd P3 identically (zero T5/T3 content of every exchange image); (iv) the RANK-3 CERTIFICATE -- channels E13 = (16,-96,144,0,0), E22 = (16,32,16,0,0), E00 = 225 E22, exchange matrix rank 2 with kernel 0, rank([M | A_fix]) = 3, annihilator certificate (Phi_T5, Phi_T3, Phi_P)(A_fix) = (0, 32, 72), and the naturalness images (25/4,-11/2,97/4,0,0) / (12,-40,76,0,0) with certificate (0,0); (v) the SO(16) AB SUM (15,-18,-15,-4,40) -- T5 does not even cancel, odd classes empty; (vi) the D8 NATIVE QUARTIC 16 T8 + 12 s8^2 with Killing 28 s8 -- no T3 structure. Mirrors v508. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4785)
- [experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/MobiusUniformisation.lean#L1)
- [(* ==== v499 round: P2.TYPING.01 -- the P2 weight-typing postulate hardened: the anchor a = (1,1,2) as the Birkhoff-Grothendieck splitting exponents of the Deligne canonical extension E of the flavor connection (types the v491 check-9 residual). Exact content mirrored: (i) partition enumeration under the full typing (BG dictionary h^0 = 0 <=> a_i >= 1, exhaustive window) + the naive- typing negative control (h^1(E') = 3 forces deg -6, three splittings, e2 never 5) + the h^1 discrepancy (h^1(anchor bundle) = 1 != 3); (ii) the Deligne residue-trace degree (cusp exponents {0,1/3,2/3} = the unique lift of spec(lam^3 - 1) to [0,1), trace 1 per mark, deg E = -4, par-deg 0, witness A0* charpoly); (iii) the Schur-Horn permutohedron integer points {(2,2,0),(2,1,1)} with h^0 = 0 killing {2,2,0}; (iv) the Biswas Z3-cover degree arithmetic (genus 2, eigensheaves O(-2)^2, chi cross-check, regular-rep weights, the decomposable model {0,2,2} with h^0 = 1 = the unstable companion; deck double genus 1, j = 1728); (v) the n-mark formula e2 = (n-2)(n+1)/2 = 5 only at n = 4. The sympy symbolic parts (H^1 residue basis, monodromy group closure |<U,M0>| = 24, the 324-case stability skeleton) stay Python-side (v499). Mirrors v499. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3621)
- [(* ==== v497 round: CELEST.WP5A.01 -- WP5a of CELEST.SEAM.01, "the null ideal from the limit" (the v496 boundary-limit shadow made precise). Exact content mirrored: (i) the chi_w stabilisation family (threshold w = n+1, w = 2 = the chiral jet grading, loop Fock 897266 at u^4); (ii) the null ideal DERIVED via the Weyl dimension formula in the standard E8 frame (Sym^2(248) = 27000 + 3875 + 1, 27000 = h_vee^3, Casimirs/Dynkin indices, level-1 integrability, quotient 31124 - 27000 = 4124); (iii) the two-route identity (mu4 theta-split sector sum at q^2 = (1036,1024,1040,1024) = 4124, q^3 = 34752, glue diagonal h = (0,1,1,1)); (iv) the negative controls (SO(16)_1: 5304 + 1820 + 135 + 1, 5304 != 14^3, quotient 2076 = Theta_D8/eta^8, block weights (0,1/2,1,1); false periodisation; swapped glue). The Freudenthal/peeling residual-zero derivation itself is Python-side (v497 S3); the Weyl-dimension and orbit-independent counting consequences are what is mirrored here. Mirrors v497. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3313)
- [(* ==== v501 round: CELEST.WP5D.01 -- WP5d-alpha of CELEST.SEAM.01, "the two-interval index from the lattice" (the local-net question, alpha stage). Exact content mirrored: (i) the Cartan determinants det(D5) = det(A3) = 4, product 16 = mu(carrier), det(D8) = 4 = mu(SO(16)_1), det(E8) = 1 = mu((E8)_1) (v378 replication, explicit matrices); (ii) the KLM/Longo-Rehren quotients mu/|L|^2 -- 16/4^2 = 1 (one-step mu4 glue, |L| = |mu4| = 4) and 16/2^2 = 4 -> 4/2^2 = 1 (via SO(16)_1), both routes exact; (iii) the Sigma d_i^2 cross-checks (Ising 1+1+2 = 4, SO(16)_1 four abelian sectors = 4, (E8)_1 one sector = 1) + the index reading arithmetic ([F:F_even] = e^(ln 2) = 2, mu_gauged = e^(2 ln 2) = 4, per-layer gauging 16 ln 2 = ln 2^16, budget ln 4 = ln mu(SO(16)_1)); (iv) the 16-fold-way condensability theta_v = e^(2 pi i nu/16) = 1 EXACTLY at nu = 2 c_- = 16 (c_- = 8 = g_car + N_fam, 16 layers = 2^(g_car - 1)) vs the rivals nu = 1, 2, 8. The lattice curves (Peschel covariance, Renyi-2 offsets, the ln 2 plateau, the duality-asymmetry witness, all ED validations) are numerical and stay Python-only (v501 S0-S4, S6). Mirrors v501. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3831)
- [(* ==== v500 round: CELEST.WP5C.01 -- WP5c of CELEST.SEAM.01, "the GNS limit state" (the state whose GNS kernel contains the v497/v498 null ideal). Exact content mirrored: (i) the LEVEL-2 BLOCK CENSUS from the standard E8 root data -- 31124 = 248 + 30876 monomials in 9361 weight blocks with dim profile {164:1, 37:240, 7:2160, 1:6960}, orbit census by lattice norm (8,6,4,2,0) -> (240, 6720, 2160, 240, 1) with constant block dims (1,1,7,37,164); (ii) the RANK/KERNEL ARITHMETIC of the full Gram -- rank table per orbit (0,0,1,8,44) (the Python-side Gram result) gives total rank 4124 = 1 + 248 + 3875 = chi_2 (dim V(omega) = 3875 recomputed by Weyl from the dominant norm-4 weight) and kernel per orbit (1,1,6,29,120) summing to 27000 = 30^3 = dim V(2 theta) (Weyl); (iii) the CLOCK-CLASS SPLIT at the state level -- the same census in the glue frame with class(lambda) = <lambda,h> mod 4 turns the rank table into the split (1036,1024,1040,1024) = Theta_Cj/eta^8 at q^2, and level 1 gives (52+8,64,60,64) = (60,64,60,64); (iv) the THRESHOLD + LEVEL-DIAL + D8 arithmetic -- x-adic order w*r >= N+1 for all r >= 1 iff w >= N+1 (sharp at w = N via r = 1), level-1 root-mode norm = k (no 248 layer at k = 0), <s|s> = 2k(k-1) (0 at k = 1, +4 at k = 2), D8 census 7380 = 120 + 7260 in 2705 blocks with quotient 7380 - 5304 = 2076. The exact Gram computation itself (9361 blocks over Fractions, PSD inertia, the highest-weight norms (3844, 49, 2, 0), the contravariance and anti-involution machine checks on 61504 pairs) is Python-side (v500 S1-S4). Mirrors v500. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3710)
- [(* ==== v498 round: CELEST.WP5B.01 -- WP5b of CELEST.SEAM.01, "the singular vector as an operator" (the deleting object of the v497 null ideal made explicit). Exact content mirrored: (i) the case classification of the J^a_1 condition over the 248-element basis (190/57/1, exactly one case sees the level k) from the standard E8 root data + the theta-grading (1,56,134,56,1); (ii) the LEVEL DIAL derived from the affine sl2 commutation -- F^theta_1 (E^theta_{-1})^n|0> = n(k-n+1) (E_{-1})^{n-1}|0>, square coefficient 2(k-1) = -2(1-k) with values (-2,0,2) at k = 0,1,2, Shapovalov 2k(k-1) (0 at k = 1, 4 at k = 2), cube singular exactly at k = 2; (iii) the glue frame -- 240 roots with classes (52,64,60,64), inner grading <r,h> = class mod 4 on all 240, theta_glue = argmax with <theta,h> = 5 => j(theta) = 1, clock phase i^(2j) = -1, cross-table rows (1,56,126,56,1); (iv) multiplicity 1 of weight 2 theta in the level-2 Fock (U(g)|s> is THE 27000 = 30^3), lattice reading (2t,2t)/2 = 4 > 2, and the comark discriminator (E8 comarks all >= 2 -- vacuum the only level-1 primary; D8 has three comark-1 weights with h = (1/2,1,1)); (v) the twisted-slot tension (two sector-C_1 modes never sum to 8 quarters, minimum 6 = q^(3/2)) + the per-period 248. The full Chevalley/PBW singularity run (J^a_1|s> = 0 on all 248 basis elements, engine unit tests on all 61504 pairs) is Python-side (v498 S1-S3); the counting/dial/grading consequences are what is mirrored here. Mirrors v498. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3464)
- [(* ==== v503 round: QGEO.EMERGE.LIGHT.01 -- the v215 deferred emergence test (lever 7) executed light + the residual convergence QGEO-R1 == P2-R1. Exact/numeric content mirrored: (i) the CM-POINT IDENTITIES -- E2(i) = 3/Pi and E2(rho) = 2 Sqrt[3]/Pi (30-digit q-series vs the classical closed forms), hence E2*(i) = E2*(rho) = 0 EXACTLY by arithmetic; j(i) = 1728, j(rho) = 0 (KleinInvariantJ); DedekindEta[I] = Gamma[1/4]/(2 Pi^(3/4)) exact; (ii) NO DYNAMICAL SELECTION -- logF = log(tau2 |eta|^4) is S-invariant (30 digits at a generic point), the HEXAGONAL point wins (logF(rho) > logF(i)) and tau = i is a saddle (max along it, min along tau1 + i, sampled); (iii) the PILLOWCASE HALVING -- exact N = 8 discrete-torus spectrum (algebraic eigenvalues), Z2 pairing with exactly 4 self-paired cone modes, dim(even) = (N^2+4)/2 = 34, and det'_even^2 = det'_torus x prod_self EXACTLY; (iv) CLOCK => SQUARE + COXETER -- P(iZ) = P(Z) forces a3 = a2 = a1 = 0 (SolveAlways), j = 6912 I^3/(4 I^3 - J^2) = 1728 identically for Z^4 + a0, the order-4 element on H^1 has char poly lam^3 + lam^2 + lam + 1 (Coxeter of W(A3), order 4 = |mu4| = N_fam + 1), and Z6 = Aut(hexagonal) has element orders {1,2,3,6} -- no order 4. The v280-style lattice DtN curves (mod-4 off-character weights, sqrtm/expm) are numerical and stay Python-only (v503 B3-B5). Mirrors v503. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4054)
- [(* ---- (v128) graded hull ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L918)
- [(* ==== v510 round: SEAM.BIT.FREEDOM.01 -- the freedom attack on the v506/v507 alignment bit: the POSITION half of the bit is topology. Exact content mirrored: (i) the DIHEDRAL CENSUS -- on the refined 2N-point seam circle D_16 has 17 involutions with EXACTLY ONE free (the antipodal shift by 8) and all 16 reflections have exactly 2 circle fixed points (8 site-axis + 8 edge-axis; robustness N = 8, 12; odd control N = 9: zero free), plus the quotient topology chi(C/antipode) = 0 (circle) vs chi(C/reflection) = 1 (interval, 2 boundary corners); (ii) the AUT COMPLETENESS -- Aut(C_8) = D_8 by brute force over all 8! = 40320 permutations and Aut(C_16) = D_16 by degree-2 propagation (all 32 candidates), so the census covers ALL circle symmetries; (iii) the Cl(16) DICHOTOMY CENSUS -- on explicit 256x256 Jordan-Wigner gammas the free deck's NS implementer squares to 256 gamma_1...gamma_16 = 256 (-1)^F (nonsplit) while the implementers of ALL 16 reflections exist (one of the two spin lifts +-R) and square to SCALARS {2^7 (site-axis), 2^8 (edge-axis)} (split): NONSPLIT <=> FREE over the complete census; root dichotomy 2 (the +-quarter shifts) vs 0 total over all 16 reflections; (iv) the REAL-STRUCTURE TABLE -- the deck-commuting antiholomorphic involutions solve symbolically to exactly two families (M diagonal: mu conj z, fix line through the poles; M antidiagonal: mu/conj z, fix circle |z| = sqrt mu), and the unique mark-fixing real structure per configuration is equatorial/FREE for mu4 and generic-unit, real-axis/NOT-free for silver and generic-real, and NONEXISTENT for the hexagonal frame; (v) HARMONIC + FREE => mu4 -- the deck-pair cross-ratio (1-b)^2/(1+b)^2 = -1 solves exactly to b = +-i, with the lambda table (2, 1/2, 3/4, 5/4) real and the hexagonal lambda complex; (vi) the HONEST COUNTERWITNESSES -- the discrete witness {0,1,8,9} (shift-invariant, crossing, but not quarter-shift-invariant) and the continuous witness {+-1, +-(3+4i)/5} (pair cross-ratio -1/4, j = 148176/25 != 1728, clock does not preserve): the Fock dichotomy measures the CLASS, not the modulus. The Moebius/Klein covering models and the pillowcase branch control stay Python-side (v510). Mirrors v510. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5155)
- [(* ==== v507 round: SEAM.BIT.ORIGIN.01 -- the tautology attack on the v506 alignment bit, executed and refuted (Part A torus origin + Part B fermion dichotomy). Exact content mirrored: (i) the HALF-PERIOD DESCENT -- the fixed points of w -> -w on C/(Z + tau Z) are the 4 half-periods (2w in Lambda <=> coefficients in {0, 1/2}; the 1/3-control fails), a translation t_c descends iff 2c in Lambda (exactly the 3 nonzero half-periods), and the descended map sigma_i(x) = e_i + (e_i-e_j)(e_i-e_k)/(x-e_i) satisfies the exact curve-lift identity prod(sigma_i - e_l) = psi_i^2 prod(x - e_l) with RATIONAL psi_i, closing to Klein V4 symbolically; (ii) the CORE SOLVE -- the pair cross-ratio condition CR = -1 solves to lambda = -1 (sigma_1), 2 (sigma_2), 1/2 (sigma_3) exactly, union = the harmonic orbit with j = 1728, against j(3) = 21952/9 and j(hex) = 0; (iii) the STABILISER COUNTS -- the complete 24-triple enumeration gives (|Stab|, #involutions, #free) = (4,3,3) generic / (8,5,3) harmonic / (12,3,3) hexagonal, and in ALL three types the free involutions are EXACTLY the three sigma_c (the deck CLASS is forced); (iv) the CM FIXED POINT -- mult-by-i mod 2 fixes exactly the half-period (1,1) = c* = (1+tau)/2 while mult-by-rho 3-cycles all three, the lemniscatic clock ghat: x -> (x-1)/(x+1) needs the curve-lift factor psi = 2i/(x+1)^2 (the CM unit; t^2+1 irreducible over Q) and squares to the central -1/x; (v) the Cl(16) DICHOTOMY -- on explicit 256x256 Jordan-Wigner gammas the central implementer squares to 256 gamma_1...gamma_16 = 256 (-1)^F (nonsplit) while the EDGE implementer Utilde_e = P_3 (g1+g5)(g2+g4)(g6-g16)(g7-g15)(g8-g14) (g9-g13)(g10-g12) implements the NS edge reflection and squares to the SCALAR +2^7 (split), and in the 32-element NS dihedral lift group the central deck has EXACTLY 2 square roots (the +-quarter shifts) vs 0 for the edge deck; (vi) the SILVER BIJECTION CENSUS -- all 8 mark bijections {+-1, +-(3+2 sqrt 2)} -> {0, 1, -1, inf} carry the deck z -> -z to the EDGE sigma_2/sigma_3 (never the central sigma_1), while all 8 mu4 bijections carry it to the CENTRAL sigma_1. The v493 circularity-audit direction and the 28/21 arrangement census stay Python-side (v507). Mirrors v507. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4546)
- [(* ---- (v168) QGEO rigidity: mu4 square cross-ratio, mu4 characters = A3 exponents ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1453)
- [(* ==== v496 round: CELEST.WP4.01 -- the (E8)_1 character E4/eta^8 vs the glue-equivariant E8[C^2] S-algebra (WP4 of CELEST.SEAM.01, the type mismatch head-on). Exact content mirrored: the character coefficients through q^6 with the q j = E4^3/eta^24 cross-check; the equivariant Hilbert series (60+128t+60t^2)/(1-t^2)^2 with the quadratic cumulative count 31 s^2 + 92 s + 60 and the SU(2)_+ slice; the two-sided jet-Fock mismatch f = (1,128,8436,381056) with strictly increasing f_n/chi_n (n = 1..12); the level-2 null ideal 27000 = 30^3 = h_vee^3; the boundary bridges (mu4 period sum 248, 248 not a tower dimension, loop Fock 897266); mirrors v496. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3258)
- [(* ==== v512 round: SEAM.TAU.FLAG.01 -- "the flag-transitivity equivalence web": the tau = i attack on the free RP seam circle. Test family M(delta) = {+-1, +-e^(i delta)}; bare mark-transitivity is delta-blind (pair-exchanging V4), FLAG transitivity <=> delta = pi/2 <=> tau = i. Six exact mirrors: (i) the V4/D4 census over five exact members incl. flag/arc orbit counts, (ii) the solveset equivalences + family identities, (iii) the odd-split closed form + arc Laplacian + K3 indicator, (iv) the twist mode identity + AB delta-blindness + j-rationality control, (v) the H4 centraliser rigidity (Schur), and (vi) the negative controls (hexagonal/silver/Z16 census/n = 3). Mirrors v512. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5676)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [(* ==== v513 round: CELEST.DTERM.NONDERIV.01 -- "the c_d negative certificate": can geometry derive c_d = 1920? Anti-numerology battery on the v511 [C] fence 1920 = |W(D5)|. Six exact mirrors: (i) the propagator anchor + T3-slot pinning, (ii) the five-way convention table + the convention-stable factorisation 32 x 60, (iii) the flux enumeration + look-elsewhere ledger with control targets, (iv) the K2 quantisation lattices, (v) the K3 provenance clash + Weyl data, and (vi) the negative controls. Mirrors v513. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L5860)
- [(* ==== v514 round: CELEST.WP5E.EPS1.01 -- "the O(-2) bulk-axion slot": WP5e-epsilon-1 executed (verdict B). Eight exact mirrors: (i) the equivariant Penrose ledger (pushforward blocks = wave-operator kernels, block by block, all d <= 6, all 4 characters), (ii) the character series P_m with the surviving bulk-axion zero mode and the twisted minimal content, (iii) the Molien hypersurface + the diag(i,i) Veronese control, (iv) the Okubo residue (6<x,x>)^2 on the 240 glue roots + the so8 irrationality control, (v) the measure cancellation excluding 3 and 12 + the flux/sector dial, (vi) the GH centre dichotomy re-deriving the v493 family, (vii) the period lockstep + the Coxeter clock from geometry, and (viii) the multipole selection rule + the EH log ledger. Mirrors v514. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L6057)
- [(* ==== v505 round: CELEST.WP5E.BETA.01 -- WP5e-beta "the equivariant anomaly ledger on twistor space" (the Atiyah-Bott/Lefschetz skeleton of the one-loop inflow on C^2/Z4 with the E8 mu4-glue monodromy; the CFT side is v502). Exact content mirrored: (i) the AB DENOMINATORS + CHARACTERS -- det(1 - g_j^{-1}) = (2,4,2), Dedekind sum 5/4 = (|Z4|^2-1)/12, equivariant characters (248,0,-8,0) by TWO routes (graded dims / explicit root sum), invariant average 60 = carrier, Frobenius 61568; (ii) the INDEX BRIDGE -- f(m) = (1/4) sum_j (i^{jm}-1)/det_j = ch2(T_m) = -(C^{-1})_mm/2 exactly, glue defect -78 by both routes with classes (-24,-30,-24), bridge identity sum_m h_A3(m) = 5/4; (iii) OKUBO PER SECTOR -- the invariant quartic sum over the 240 glue roots is EXACTLY 36 <x,x>^2 (36 = lambda~^2_e8), the twisted quartics are NOT proportional to squares (T5/T3 content), and the AB-weighted sum leaves the rigid residual 32 T3; (iv) the LEVEL DIALS -- lattice current count 240/0/0/0 for k = 1..4, embedding residual 47k^2+219k-266 = (0,360,814,1362), det(A-1) = -4 (one scale => one level); (v) the AXION-SLOT WEIGHTS -- a0 fills the graviton slot O(2) (Hamiltonian weight arithmetic X = -2 != +2 = Y, mismatch 4 = |mu4|), the clock weights on (a3,a2,a1,a0) are the regular mu4 representation, and the Coxeter eigencharacters {i,-1,-i} = the three twisted sector characters. Mirrors v505. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4245)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ==== v504 round: CELEST.WP5DB.01 -- WP5d-beta "split + strong additivity from the lattice", the two remaining KLM legs of complete rationality (the mu-index is v501). Exact content mirrored: (i) the GF2 SPAN THEOREM -- Majorana monomials multiply as labels in GF(2)^n (c_S c_T = +- c_{S xor T}), so span dims are exact set counts: with the shared boundary Majorana Even(A) v Even(B) = Even(A u B) EXACTLY (64/64, 256/256, 1024/1024), disjoint EXACTLY half (64/128, 256/512, 1024/2048) = index 2; (ii) the PIMSNER-POPA IDENTITY -- E(a) = (a + PaP)/2 gives E(a) - a/2 = PaP/2 symbolically on a full 4x4 matrix, and the attainment x = 1 + h (h = c_1 odd, h^2 = 1, PhP = -h, explicit gamma matrices) has E(x*x) - lam x*x = 2(1 - lam) - 2 lam h with spec(h) = {+-1}: nonnegative iff lam <= 1/2 = 1/[F:F_even] EXACTLY; (iii) the TWO-INTERVAL GROUP -- E4 over {1,PA,PB,PAB} on the uniform 4-sector vector gives E4(vv+) = I/4 exactly, so lam_E4 = 1/4 = 1/mu, and the index consistency (1/lam_PP)^2 = 4 = 1/lam_E4 = mu(SO(16)_1) = det Cartan(D8); (iv) the LAMBDA^2 COMPOUND -- the Wick cross-Gram of even bilinears is the second compound Minors[C, 2]: Cauchy-Binet functoriality Minors[A.B, 2] = Minors[A, 2].Minors[B, 2] on explicit integer matrices, Minors[diag(s), 2] = diag(s_a s_b) symbolically (singular values = products, the same ladder at doubled rate), and Minors[-C, 2] == Minors[C, 2] IDENTICALLY (the parity flip C -> -C is invisible to the even bilinears -- the exact orbifold split-inheritance mechanism); (v) the U(1) PP CONTROL -- charge dephasing on m modes has m + 1 sectors and the uniform-spread vector gives E_U(vv+) = P_diag/(m + 1), so lam = 1/(m + 1) -> 0 EXACTLY (m = 2, 3): infinite index. The lattice curves (deficit tables, Ising exponent p = 0.2444, elliptic-nome ladder, Klich-Levitov slope, ED validations) are numerical and stay Python-only (v504 S0/S2/S3/S5). Mirrors v504. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4133)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* ==== v518 round: CELEST.WP5E.DELTA1.01 -- "the delta-1 chain decided: the derived chiral measure kills the delta-1 route" (delta-1b/1c/1d consolidated; an honest, decided NEGATIVE result with the v516 tension stated). Six exact mirrors: (i) the 16-component Weil system of the D5 (+) A3 discriminant module Z4 x Z4 with q = (5x^2 + 3y^2)/8 -- Gauss sums 2 zeta_8^5 x 2 zeta_8^3 = 4 = sqrt(16) (signature 0 mod 8 = rank E8), S symmetric + unitary, S^2 = C, (ST)^3 = S^2, S^4 = 1, T^8 = 1, tensor split S = S_D5 (x) S_A3 with factor signatures e(-5/8), e(-3/8); (ii) the two invariant Lagrangians (diagonal H and anti-diagonal H', both S- and T-fixed = two E8 gluings), 6 isotropic classes, P_H S P_H = (1/4) J and 3/4 of every S-image weight on the 12 shifted classes (why the naive 4-character rule fails); (iii) the KOBOUNDARY TEST on the recorded exact transport phases (t_p, s_p) of the dressed gauge blocks M = G[a,b] eta^{-8}: the SL(2,Z) relation defects c_{S^4}, c_{(ST)^3 S^-2}, c_{[S^2,T]} are (1, 1, 1) EXACTLY on all 15 sector pairs -- the mu4 multiplier system is a CHARACTER of the orbit stabilisers, not a genuine 2-cocycle; (iv) the Gamma_1(4) CHARACTER FORMULA lambda(gamma) = i^{2B + C/4} on the recorded loop matrices (all in Gamma_1(4), exact integer check), with all four equivalent monomial solutions agreeing; (v) the T-FIX MECHANISM: v502 vacuum energies E_b exact, the f1 f3 dressing shifts every T-fixed node by e(2E + 2E) = e(-1/6) and (-1) x e(-1/6) = e(1/3) = chi_4(T) UNIFORMLY, while f1f2f3 demands {chi_3, chi_4} and f2 demands {chi_5, chi_6} (inconsistent: only the twisted fibre block supplies one character); (vi) the misassignment controls (wrong form (x^2+y^2)/4: |Gauss|^2 = 64 != 16 and S^2 != C exactly; wrong signature sigma = e(-1/4): S^2 = e(-1/2) C exactly). The Harvey-Moore integrals, the SVD contraction solves, the theta covariance certificates and the block transport measurements are numerical (mpmath 30+ digits / numpy) and stay Python-only, flagged here. Mirrors v518. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L6877)
- [(* ==== v502 round: CELEST.WP5E.ALPHA.01 -- WP5e-alpha of CELEST.SEAM.01, "prefactor and level bookkeeping" (the CFT-side anchor of the twistor uplift). Exact content mirrored: (i) the HURWITZ-ZETA vacuum energies -- zeta(-1,theta) = -B2(theta)/2 at rational twists, E_b(theta) = -1/24 + theta(1-theta)/4, Ising shift E_R - E_NS = 1/16 per Majorana, the untwisted vacuum 8 x (-1/24) = -1/3 = -c/24 with the Sugawara route 248/(1+30) = 8, and the eta^8 exponent 8/24 = 1/3; (ii) DISCRIMINANT FORM = CASIMIR ENERGY -- coset minima h_D5 = (0,5/8,1/2,5/8), h_A3 = (0,3/8,1/2,3/8), glue diagonal (0,1,1,1) INTEGER (exact lattice enumeration), spectral flow j^2 (h,h)/32 with (h,h) = 20, (h',h') = 12, sum 32 (mod-1 identity), free-fermion route n/16 = 5/8, 3/8, 1 on 10 + 6 = 16 = 2^(g_car-1) Majoranas, and the deck failure 3/16 != 3/8; (iii) the k = 1 FIXING EQUATIONS -- simple-current closed forms 5k/8, 3k/8, glue-diagonal integrality for ALL k = 1..8 (the honest 'fixes nothing' finding), current condition h(J) = k = 1 unique, conformal embedding 47k^2 + 219k - 266 = (k-1)(47k+266) with unique positive root 1, D8 route 128k(1-k) = 0, central charge 248k/(k+30) = 8 <=> 240(k-1) = 0, Shapovalov values 2k(k-1) = (0,4,12,24), and the Weyl-dim count 31124 - 27000 = 4124 (dims 248/3875/27000 recomputed); (iv) the CHARACTER SUMS -- the four glue-coset thetas sum to E4 = (1,240,2160, 6720), the sector characters (leading coefficients (1,64,60,64), currents (60,64,60,64)) sum to E4/eta^8 = q^{-1/3}(1,248,4124,34752), and the D8 control Theta_D8 + Theta_s = E4 with h_v = 1/2 non-integer. The sympy symbolic manipulations and the PBW engine are Python-side (v502 S1/S5). Mirrors v502. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L3894)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ==== v506 round: SEAM.CLOCK.RIGIDITY.01 -- clock rigidity, Route B of the v503 classification (Part A Moebius + Part B Fock). Exact content mirrored: (i) the COMPLETE SQUARE-ROOT SOLVE -- g^2 = deck in PGL2(C) has exactly two solutions z -> +-iz (tr != 0 forces diagonal with (a/d)^2 = -1; tr = 0 forces det = 0), both automatically of projective order 4; (ii) the D4 CENTER + ORDER-4 COUNT + the COUNTEREXAMPLE -- the mu4 stabiliser D4 = {i^k z, i^k/z} has exactly 2 order-4 elements (the clocks), both squaring to the deck, center = {id, deck}; the silver harmonic set {+-1, +-(3+2 sqrt 2)} has j = 1728 and full D4 but its 4-cycle squares to the CROSSING involution v/z, not the deck (deck non-central: 0 roots -- harmonicity alone insufficient); (iii) the Cl(16) IDENTITY -- on explicit 256x256 Jordan-Wigner gamma matrices, Utilde = prod_j (1 - g_j g_{j+8}) satisfies Utilde^2 = 256 g_1...g_16 EXACTLY (U^2 = (-1)^F, the forced nonsplit Z4), the R implementer squares to +256 (split control), and the volume element anticommutes with every gamma; (iv) the N-ARITHMETIC -- S4 element orders {1,2,3,4} (faithfulness kills n >= 5: max cyclic image order < n for n = 5..8), free orbits force n | 4, the cusp side kills n = 6 (deg(lam^5-1) = 5 != 3 = deg(lam^3-1)), transfer spectrum {1, 64/729, 1/729} with gap (2/3)^6, and the generic 4-set stabiliser V4 = {+-z, +-2/z} has exponent 2 (not even the deck has a root). The full 24-triple stabiliser enumeration and the A4/hexagonal profile stay Python-side (v506 Part A). Mirrors v506. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4395)

## -1.2cm TFPT — Frontier items\\[2mm] [formula in source] , [formula in source] , the Koide relation, dark matter and quantum gravity —\\ the honest status, not forced onto the ladder

[documents/tfpt_4_frontier/preamble.tex:1-37](../../../../documents/tfpt_4_frontier/preamble.tex#L1) · `section:e04339a135a123052daa` · frontmatter

-1.2cm TFPT — Frontier items\\[2mm] [formula in source] , [formula in source] , the Koide relation, dark matter and quantum gravity —\\ the honest status, not forced onto the ladder

## Front matter / support

[documents/tfpt_4_frontier/sections/001-source.tex:1-90](../../../../documents/tfpt_4_frontier/sections/001-source.tex#L1) · `section:fbc55d19f88e095d54eb` · frontmatter

The honest frontier: which physics has a genuine TFPT handle and which does not. For each of [formula in source] , [formula in source] , the Koide relation, dark matter and full quantum gravity, this note states the genuine structural handle, the precision with which it currently lands, and — crucially — what is not a clean compiler power and is deliberately not forced onto the ladder. This document is the status aut

Explicit source references:

- [verification/status_ledger.csv](../../../../verification/status_ledger.csv#L1) · explicit_file_citation
- [CONTRACT.F.01](../../../../verification/status_ledger.csv#L302) · explicit_claim_identifier · original status at linked ledger row
- [verification/v375_observatory_registry.py](../../../../verification/v375_observatory_registry.py#L1) · explicit_script_identifier
- [verification/v211_axion_spine_angle.py](../../../../verification/v211_axion_spine_angle.py#L1) · explicit_script_identifier
- [FTRANSFER.SUITE.01](../../../../verification/status_ledger.csv#L472) · explicit_claim_identifier · original status at linked ledger row
- [verification/v187_ftransfer_laws.py](../../../../verification/v187_ftransfer_laws.py#L1) · explicit_file_citation
- [verification/v402_ftransfer_suite.py](../../../../verification/v402_ftransfer_suite.py#L1) · explicit_file_citation
- [verification/v35_quark_qcd_boundary.py](../../../../verification/v35_quark_qcd_boundary.py#L1) · explicit_file_citation
- [verification/v792_ftransfer_kms_coset.py](../../../../verification/v792_ftransfer_kms_coset.py#L1) · explicit_file_citation
- [verification/v303_ftransfer_dynamics.py](../../../../verification/v303_ftransfer_dynamics.py#L1) · explicit_script_identifier
- [verification/v188_frontier_wording_guard.py](../../../../verification/v188_frontier_wording_guard.py#L1) · explicit_file_citation
- [verification/v212_leptogenesis_decuple.py](../../../../verification/v212_leptogenesis_decuple.py#L1) · explicit_script_identifier
- [verification/v371_ftransfer_pole.py](../../../../verification/v371_ftransfer_pole.py#L1) · explicit_script_identifier
- [verification/v213_ftransfer_functor.py](../../../../verification/v213_ftransfer_functor.py#L1) · explicit_file_citation
- [verification/v777_ftransfer_clock_jets.py](../../../../verification/v777_ftransfer_clock_jets.py#L1) · explicit_file_citation
- [verification/v397_external_clock_probe.py](../../../../verification/v397_external_clock_probe.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)

## Route A — from the closed baryon fraction.

[documents/tfpt_4_frontier/sections/001-source.tex:95-109](../../../../documents/tfpt_4_frontier/sections/001-source.tex#L95) · `section:84d26a09031108cede24` · paragraph

The seed already fixes the baryon fraction (the [formula in source] decoder): [formula in source] matching the observed [formula in source] . The baryon-to-photon ratio is then the standard cosmological readout [formula in source] . With the TFPT Hubble value ( [formula in source] , from [formula in source] ), [formula in source] So [formula in source] is a downstream cosmological readout from the closed [formula in

## The Higgs quartic — near-criticality from the free seam

[documents/tfpt_4_frontier/sections/001-source.tex:262-300](../../../../documents/tfpt_4_frontier/sections/001-source.tex#L262) · `section:cc7b2d0df421c394451c` · section

The seam UV is the free chiral [formula in source] fixed point ( v158_fixed_point_stable.py): Gaussian, with no relevant and no marginal interactions. The Higgs self-coupling [formula in source] is the one marginal scalar coupling, so a free seam forces it to the Gaussian point there: [formula in source] the Shaposhnikov–Wetterich double-criticality condition — here derived from the free seam rather than postulated.

Explicit source references:

- [verification/v166_higgs_free_seam.py](../../../../verification/v166_higgs_free_seam.py#L1) · explicit_file_citation
- [verification/v164_qcd_scale.py](../../../../verification/v164_qcd_scale.py#L1) · explicit_file_citation
- [verification/v159_pyrate_gauge_crosscheck.py](../../../../verification/v159_pyrate_gauge_crosscheck.py#L1) · explicit_file_citation
- [verification/v158_fixed_point_stable.py](../../../../verification/v158_fixed_point_stable.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v158) free chiral c=8 fixed point is stable ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1422)
- [(* ---- (v159) PyR@TE gauge cross-check: carrier content -> SM 1-loop b_i ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1433)

## Full quantum gravity — induced from the boundary kernel

[documents/tfpt_4_frontier/sections/003-source.tex:1-376](../../../../documents/tfpt_4_frontier/sections/003-source.tex#L1) · `section:f5274c8105ff62e99dcd` · section

Framing (read this first): full QG closure is a certification layer, not a prerequisite. TFPT does not need the complete quantum-gravity measure to test its Standard-Model and cosmology readouts — those run through [formula in source] and [formula in source] and are falsifiable now (JUNO, CMB-S4, EDM). [formula in source] certifies the metric sector; its absence narrows what can be claimed, not what can be tested. \\

Explicit source references:

- [verification/v473_entropic_action_bridge.py](../../../../verification/v473_entropic_action_bridge.py#L1) · explicit_script_identifier
- [verification/v993_minimal_defect_selector.py](../../../../verification/v993_minimal_defect_selector.py#L1) · explicit_file_citation
- [verification/v1031_quantum_tensor_curvature.py](../../../../verification/v1031_quantum_tensor_curvature.py#L1) · explicit_file_citation
- [verification/v3_em_alpha.py](../../../../verification/v3_em_alpha.py#L1) · explicit_script_identifier
- [verification/v471_seam_horizon_replica.py](../../../../verification/v471_seam_horizon_replica.py#L1) · explicit_script_identifier
- [verification/v1008_master_assembly_scaffold.py](../../../../verification/v1008_master_assembly_scaffold.py#L1) · explicit_file_citation
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [GRAV.WEINBERG.WITTEN.01](../../../../verification/status_ledger.csv#L770) · explicit_claim_identifier · original status at linked ledger row
- [verification/v9_neutrino_texture.py](../../../../verification/v9_neutrino_texture.py#L1) · explicit_script_identifier
- [verification/v367_seam_s3_lattice.py](../../../../verification/v367_seam_s3_lattice.py#L1) · explicit_script_identifier
- [verification/v994_mmst_criteria_battery.py](../../../../verification/v994_mmst_criteria_battery.py#L1) · explicit_file_citation
- [verification/v1018_e8_directed_readout.py](../../../../verification/v1018_e8_directed_readout.py#L1) · explicit_file_citation
- [QFT4D.OS.RECON.01](../../../../verification/status_ledger.csv#L745) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1006_mmst_lemma_battery.py](../../../../verification/v1006_mmst_lemma_battery.py#L1) · explicit_file_citation
- [verification/v1009_rho0_minimizer.py](../../../../verification/v1009_rho0_minimizer.py#L1) · explicit_file_citation
- [FTRANSFER.SK.RHO0.01](../../../../verification/status_ledger.csv#L763) · explicit_claim_identifier · original status at linked ledger row
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v998_seam_modular_closure.py](../../../../verification/v998_seam_modular_closure.py#L1) · explicit_file_citation
- [experiments/theory-contracts/RESEARCH_2026-09-09.md](../../../../experiments/theory-contracts/RESEARCH_2026-09-09.md#L1) · explicit_file_citation
- [CHIRAL4D.NOMIRROR.01](../../../../verification/status_ledger.csv#L746) · explicit_claim_identifier · original status at linked ledger row
- [verification/v277_seam_calderon_e8_match.py](../../../../verification/v277_seam_calderon_e8_match.py#L1) · explicit_script_identifier
- [verification/v379_seam_s3_rp.py](../../../../verification/v379_seam_s3_rp.py#L1) · explicit_script_identifier
- [verification/v999_weak_collision_hp.py](../../../../verification/v999_weak_collision_hp.py#L1) · explicit_file_citation
- [QFT4D.LATTICE.FUNDAMENTAL.01](../../../../verification/status_ledger.csv#L760) · explicit_claim_identifier · original status at linked ledger row
- [verification/v275_qgamb_roadmap.py](../../../../verification/v275_qgamb_roadmap.py#L1) · explicit_file_citation
- [verification/v624_external_lattice_audit.py](../../../../verification/v624_external_lattice_audit.py#L1) · explicit_script_identifier
- [FORM.SEAM.MMST.01](../../../../verification/status_ledger.csv#L92) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1030_alg2_sector_frame.py](../../../../verification/v1030_alg2_sector_frame.py#L1) · explicit_file_citation
- [E8.DIRECTED.READOUT.01](../../../../verification/status_ledger.csv#L771) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1026_telb_round4_closure.py](../../../../verification/v1026_telb_round4_closure.py#L1) · explicit_file_citation
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v996_gauge_grammar_kill.py](../../../../verification/v996_gauge_grammar_kill.py#L1) · explicit_file_citation
- [verification/v205_xi_threequarter.py](../../../../verification/v205_xi_threequarter.py#L1) · explicit_file_citation
- [verification/v1032_covariant_curvature.py](../../../../verification/v1032_covariant_curvature.py#L1) · explicit_file_citation
- [verification/v1010_simplicity_bridge_census.py](../../../../verification/v1010_simplicity_bridge_census.py#L1) · explicit_file_citation
- [SEAM.MMST.TYPEIII.CHARGED.01](../../../../verification/status_ledger.csv#L767) · explicit_claim_identifier · original status at linked ledger row
- [verification/v368_seam_s3_inflow.py](../../../../verification/v368_seam_s3_inflow.py#L1) · explicit_script_identifier
- [verification/v269_spert_paqft_skeleton.py](../../../../verification/v269_spert_paqft_skeleton.py#L1) · explicit_script_identifier
- [verification/v354_e8_reverse_audit.py](../../../../verification/v354_e8_reverse_audit.py#L1) · explicit_script_identifier
- [verification/v989_lattice_gate_battery.py](../../../../verification/v989_lattice_gate_battery.py#L1) · explicit_file_citation
- [verification/v1015_axiom_core_closure.py](../../../../verification/v1015_axiom_core_closure.py#L1) · explicit_file_citation
- [verification/v1033_charged_disorder.py](../../../../verification/v1033_charged_disorder.py#L1) · explicit_file_citation
- [verification/v1029_local_tensor_constraints.py](../../../../verification/v1029_local_tensor_constraints.py#L1) · explicit_file_citation
- [verification/v423_no_ambient_dependency.py](../../../../verification/v423_no_ambient_dependency.py#L1) · explicit_file_citation
- [verification/v480_multilocal_four_interval.py](../../../../verification/v480_multilocal_four_interval.py#L1) · explicit_script_identifier
- [verification/v1000_spin2_sk_witnesses.py](../../../../verification/v1000_spin2_sk_witnesses.py#L1) · explicit_file_citation
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_file_citation
- [verification/v355_e8_unmapped_bandwidth.py](../../../../verification/v355_e8_unmapped_bandwidth.py#L1) · explicit_script_identifier
- [verification/v36_spectral_action_g2.py](../../../../verification/v36_spectral_action_g2.py#L1) · explicit_file_citation
- [verification/v273_eg_gauge_running.py](../../../../verification/v273_eg_gauge_running.py#L1) · explicit_script_identifier
- [verification/v7_gravity_cosmo.py](../../../../verification/v7_gravity_cosmo.py#L1) · explicit_script_identifier
- [GRAV.SPIN2.EMERGENCE.01](../../../../verification/status_ledger.csv#L762) · explicit_claim_identifier · original status at linked ledger row
- [verification/v358_grav_entropy_equilibrium.py](../../../../verification/v358_grav_entropy_equilibrium.py#L1) · explicit_script_identifier
- [verification/v268_theta13_carrier_trace.py](../../../../verification/v268_theta13_carrier_trace.py#L1) · explicit_script_identifier
- [verification/v1013_thermodynamic_dynamics.py](../../../../verification/v1013_thermodynamic_dynamics.py#L1) · explicit_file_citation
- [verification/v334_conformal_resolution.py](../../../../verification/v334_conformal_resolution.py#L1) · explicit_script_identifier
- [E8.COXETER.EULER.COMPLETION.01](../../../../verification/status_ledger.csv#L772) · explicit_claim_identifier · original status at linked ledger row
- [verification/v60_lambda_metrology_branch.py](../../../../verification/v60_lambda_metrology_branch.py#L1) · explicit_file_citation
- [verification/v1012_toe_gate_witnesses.py](../../../../verification/v1012_toe_gate_witnesses.py#L1) · explicit_file_citation
- [verification/v332_qgamb_metric_sector.py](../../../../verification/v332_qgamb_metric_sector.py#L1) · explicit_script_identifier
- [verification/v1028_gauss_d4_cap.py](../../../../verification/v1028_gauss_d4_cap.py#L1) · explicit_file_citation
- [SEAM.MILNOR.LOCALRING.01](../../../../verification/status_ledger.csv#L769) · explicit_claim_identifier · original status at linked ledger row
- [verification/v995_kms_subnet_rigidity.py](../../../../verification/v995_kms_subnet_rigidity.py#L1) · explicit_file_citation
- [verification/v1020_spectator_sector_obstruction.py](../../../../verification/v1020_spectator_sector_obstruction.py#L1) · explicit_script_identifier
- [verification/v360_grav_gap_corrections.py](../../../../verification/v360_grav_gap_corrections.py#L1) · explicit_file_citation
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_script_identifier
- [verification/v1014_bridge_refinements.py](../../../../verification/v1014_bridge_refinements.py#L1) · explicit_file_citation
- [verification/v23_anchor_generator.py](../../../../verification/v23_anchor_generator.py#L1) · explicit_script_identifier
- [AX.P2.01](../../../../verification/status_ledger.csv#L3) · explicit_claim_identifier · original status at linked ledger row
- [verification/v323_bw_geometric_modular.py](../../../../verification/v323_bw_geometric_modular.py#L1) · explicit_script_identifier
- [verification/v1001_nu_pentagon_chain.py](../../../../verification/v1001_nu_pentagon_chain.py#L1) · explicit_file_citation
- [FLAV.NUSCALE.06](../../../../verification/status_ledger.csv#L765) · explicit_claim_identifier · original status at linked ledger row
- [verification/v271_eg_oneloop_quartic.py](../../../../verification/v271_eg_oneloop_quartic.py#L1) · explicit_script_identifier
- [verification/v1011_3p1d_ladder_battery.py](../../../../verification/v1011_3p1d_ladder_battery.py#L1) · explicit_file_citation
- [verification/v1002_spin10_mirror_projector.py](../../../../verification/v1002_spin10_mirror_projector.py#L1) · explicit_file_citation
- [verification/v376_seam_s3_centralcharge.py](../../../../verification/v376_seam_s3_centralcharge.py#L1) · explicit_script_identifier
- [verification/v206_h0_lambda_branch.py](../../../../verification/v206_h0_lambda_branch.py#L1) · explicit_file_citation
- [verification/v997_spin2_smg_witnesses.py](../../../../verification/v997_spin2_smg_witnesses.py#L1) · explicit_file_citation
- [TFPT4D.LATTICE.ACTION.01](../../../../verification/status_ledger.csv#L757) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1003_nu_structural_nulls.py](../../../../verification/v1003_nu_structural_nulls.py#L1) · explicit_file_citation
- [verification/v152_norm_is_anchor.py](../../../../verification/v152_norm_is_anchor.py#L1) · explicit_file_citation
- [verification/v337_decoupling_theorem.py](../../../../verification/v337_decoupling_theorem.py#L1) · explicit_file_citation
- [TFPT.TOE.COMPLETE.01](../../../../verification/status_ledger.csv#L764) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1019_coxeter_euler_completion.py](../../../../verification/v1019_coxeter_euler_completion.py#L1) · explicit_file_citation
- [verification/v1035_matter_coupling.py](../../../../verification/v1035_matter_coupling.py#L1) · explicit_file_citation
- [verification/v1004_milnor_frobenius_gray.py](../../../../verification/v1004_milnor_frobenius_gray.py#L1) · explicit_file_citation
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [SEAM.EQUIV.MMST.01](../../../../verification/status_ledger.csv#L583) · explicit_claim_identifier · original status at linked ledger row
- [verification/v281_holomorphy_modular_data.py](../../../../verification/v281_holomorphy_modular_data.py#L1) · explicit_script_identifier
- [verification/v1005_relative_pencil_bridges.py](../../../../verification/v1005_relative_pencil_bridges.py#L1) · explicit_file_citation
- [QGAMB.NODEP.01](../../../../verification/status_ledger.csv#L490) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1034_parent_mirror.py](../../../../verification/v1034_parent_mirror.py#L1) · explicit_file_citation
- [CHIRAL4D.MIRROR.DET16.01](../../../../verification/status_ledger.csv#L766) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1027_signed_det_car_wall.py](../../../../verification/v1027_signed_det_car_wall.py#L1) · explicit_file_citation
- [verification/v255_spectral_action_expansion.py](../../../../verification/v255_spectral_action_expansion.py#L1) · explicit_file_citation
- [QGAMB.ROADMAP.01](../../../../verification/status_ledger.csv#L368) · explicit_claim_identifier · original status at linked ledger row
- [verification/v77_e8_conformal_net.py](../../../../verification/v77_e8_conformal_net.py#L1) · explicit_script_identifier
- [verification/v1016_state_gap_batteries.py](../../../../verification/v1016_state_gap_batteries.py#L1) · explicit_file_citation
- [verification/v58_seam_horizon_chain.py](../../../../verification/v58_seam_horizon_chain.py#L1) · explicit_script_identifier
- [verification/v216_marks_gauss_bonnet.py](../../../../verification/v216_marks_gauss_bonnet.py#L1) · explicit_script_identifier
- [FLAV.NU.TEXTURE.MECHANISM.01](../../../../verification/status_ledger.csv#L768) · explicit_claim_identifier · original status at linked ledger row
- [verification/v1007_det16_stability_battery.py](../../../../verification/v1007_det16_stability_battery.py#L1) · explicit_file_citation
- [GAUGE.DETLINE.FIXPOINT.01](../../../../verification/status_ledger.csv#L761) · explicit_claim_identifier · original status at linked ledger row
- [verification/v361_grav_backreaction.py](../../../../verification/v361_grav_backreaction.py#L1) · explicit_file_citation
- [verification/v359_grav_nonlinear_einstein.py](../../../../verification/v359_grav_nonlinear_einstein.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v58) seam-horizon chain: one-sided Gauss-Bonnet, KMS, S_BH, RT in seam units ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L270)
- [(* ==== v355 E8.UNMAPPED.BANDWIDTH.01: the unmapped region's forced content is COLLECTIVE ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2363)
- [(* ---- (v60) Lambda branch selection: prefactor 3/(4pi^2), branch mis-scale 64pi^3/3 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L283)
- [(* ==== v271: concrete Epstein-Glaser one-loop quartic renormalization (QFT4D.SPERT.02) ==== v269 (S_pert skeleton) + v270 (PMNS Jarlskog, numerical) are Python-only; v271's exact EG core (scaling degree, counterterm count, loop factor) is mirrored here. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2083)
- [(* v1028: exact ordinary D4 cap on the selected 12D even pair space. *)](../../../../verification/wolfram/tfpt_round4.wl#L113)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* ==== v268: reactor-angle exponent = carrier hypercharge trace (FLAV.TH13.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2070)
- [(* ---- (v994) SEAM.EQUIV.MMST.01 C3: Theta_E8 / eta^8 = {1,248,4124,34752} ---- Python-only: C1 entropy, C2 QWZ spectrum, C4 two-point numeric, C5 sector recombination beyond det, C6 cited, C7 open, GSO/N=8 must-fail combinatorics. *)](../../../../verification/wolfram/tfpt_readouts.wl#L681)
- [(* v1030: actual connected-sector words span M4; a disconnected mutant leaves one extra commutant degree. *)](../../../../verification/wolfram/tfpt_round4.wl#L163)
- [(* ---- (v216) QGEO.MARKS.03: the four marks from Gauss-Bonnet + Euclidean-orbifold uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1741)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [(* v1027: exact signed-wall determinant and unsquared low branch. *)](../../../../verification/wolfram/tfpt_round4.wl#L29)
- [(* ==== v281: anyon-condensation tower + Gauss-Milgram (QGAMB.MODULAR.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2154)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/BWKeystoneCurrent.lean#L1)
- [(* v473 (i): carrier Hodge count + coefficient pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2968)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/YukawaStageDExistence.lean#L1)
- [(* ==== v354 E8.REVERSE.AUDIT.01: 3 of 8 E8 Casimir degrees are load-bearing, 5 unmapped ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2353)
- [experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/BWKeystone.lean#L1)
- [(* v473 (iii): the R^2 kill-test gap *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2985)
- [(* ==== v358 GRAV.ENTROPY.EQUILIBRIUM.01: parameter-free linearised Einstein + thermo=geo ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2372)
- [experiments/lean4-carrier-rigidity/TfptCarrier/LorentzCongruence.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/LorentzCongruence.lean#L1)
- [(* ---- (v3) electromagnetic closure: alpha^-1 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L37)
- [(* v1029: exact generic local tensor identities and the k=0 instability. *)](../../../../verification/wolfram/tfpt_round4.wl#L136)
- [(* ==== v337 DECOUPLING.THEOREM.01: the ambient back-reaction scale identity ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2241)
- [(* ---- (v77) G6 via E8 level-1 net: level-1 central charges + conformal embedding c-sum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L377)
- [(* v1028: exact characteristic-five reduction of the full 192-mode one-site character. P^12 == P^2 P(z^5)^2 mod 5 is evaluated by coefficient dictionaries, never by floating group quadrature. *)](../../../../verification/wolfram/tfpt_round4.wl#L41)
- [(* v1034: exact source parent envelope. The source's 3/5 is a dispersion amplitude; a standard symmetric NN bond has half of it. *)](../../../../verification/wolfram/tfpt_round7.wl#L161)
- [(* ---- (v23) anchor generator a=(1,1,2) ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L96)
- [(* ==== v473 round: the entropic-action bridge (Bianconi, PRD 111, 066001 (2025); arXiv:2408.14391) -- the exact algebraic content: the carrier Hodge count 1+5+10 = 16 = 2^(g_car-1) (vs the 4d fiber 11), the coefficient pin beta'_B = c3/6 = 1/(48 pi) from 3 beta' = 1/(16 pi), the entropy potential s(x) = x-1-ln x (nonnegative, quadratic minimum), the Lambda-channel target Tr Q^2 = 32 c3^4 = 1/(128 pi^4) = (2/3) delta_top, and the R^2 kill-test gap 3 c3^-9 = 3 (8 pi)^9; mirrors v473 (the Frullani quadrature is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2958)
- [(* v1035: physical-position Fourier convention D -> i k. Retain I: c0=-s.q and c_i=i (B p)_i. The source Ward equations then propagate both affine constraints exactly. *)](../../../../verification/wolfram/tfpt_round7.wl#L182)
- [(* ---- (v205) GRAV.XI.01: xi = c3/phi_tree = 3/4, the independent gravitational 3/4 ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1709)
- [(* v473 (ii): entropy potential + Lambda-channel target *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2975)
- [(* ---- (v152) R3 normalisation = the anchor ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1356)
- [(* v1031: exact rational-fibre curvature projector and positivity. *)](../../../../verification/wolfram/tfpt_round7.wl#L80)
- [experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/Mu4Commutation.lean#L1)
- [(* ---- (v36) G2: spectral-action Seeley-DeWitt coefficients -> R + R^2 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L113)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AnchorLadder.lean#L1)
- [(* ==== v273: EG one-loop gauge self-energy -> (b1,b2,b3) (QFT4D.SPERT.03) ==== v272 (nu-scale, numerical) + v274 (over-determination, numerical) + v275 (QG roadmap) are Python-only; v273's exact one-loop beta coefficients from the content are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2103)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [(* ---- (v993) AX.P2.01 census: SNF disc groups, unique (D5,A3), [lambda]^2=[v], CE dims ---- Python-only: random CE axiom sampling on M_d^{Z4} / crossed product. *)](../../../../verification/wolfram/tfpt_readouts.wl#L558)
- [(* v1033: exact Z4 deck algebra and the load-bearing weight mismatch. *)](../../../../verification/wolfram/tfpt_round7.wl#L142)
- [(* ==== v359 GRAV.NONLINEAR.01: full covariant Einstein eq, both coefficients parameter-free ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2383)
- [(* ==== v277: seam-Calderon -> (E8)_1 matching certificate (QGAMB.TIERB.01) ==== v276 (flat-pillowcase commutator, numerical) is Python-only; v277's exact lattice discriminator (det E8 vs det D8) + the (E8)_1 character 248 are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2120)
- [(* ---- (v997) GRAV.SPIN2.EMERGENCE.01: Barnes-Rivers ranks (5,3,1,1) ---- Python-only: Weyl^2 ghost residues, random-p numeric, SMG N=8 Fidkowski-Kitaev toy (spectrum / Spin(7) / Schmidt / edge). *)](../../../../verification/wolfram/tfpt_readouts.wl#L700)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CartanDeterminants.lean#L1)
- [(* v1032: a separate Lorentz-covariant curvature calculation at one exact nonzero null momentum. Pi is the de Donder trace-reversed numerator; positivity is inferred only after the curvature sandwich. *)](../../../../verification/wolfram/tfpt_round7.wl#L92)

## Appendix: computational verification

[documents/tfpt_4_frontier/sections/004-source.tex:200-215](../../../../documents/tfpt_4_frontier/sections/004-source.tex#L200) · `section:0c7b2f5d7f77f90fe85c` · section

These are the frontier items: most carry \ (open) or \ (conditional) and are therefore not part of the pass/fail suite. The one quantitative handle that is reproduced is [formula in source] (whence the [formula in source] order) in verification/v7_gravity_cosmo.py. The [formula in source] structure and the gap-decoupling of the physical QG sector are machine-checked ( v28_gravity_fR.py, v36_spectral_action_g2.py); th

Explicit source references:

- [verification/v28_gravity_fR.py](../../../../verification/v28_gravity_fR.py#L1) · explicit_file_citation
- [verification/run_all.py](../../../../verification/run_all.py#L1) · explicit_file_citation
- [verification/README.md](../../../../verification/README.md#L1) · explicit_file_citation
- [verification/v7_gravity_cosmo.py](../../../../verification/v7_gravity_cosmo.py#L1) · explicit_file_citation
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v374_ftransfer_qcd.py](../../../../verification/v374_ftransfer_qcd.py#L1) · explicit_script_identifier
- [verification/v36_spectral_action_g2.py](../../../../verification/v36_spectral_action_g2.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v36) G2: spectral-action Seeley-DeWitt coefficients -> R + R^2 ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L113)

## -1.2cm TFPT — Red Team / Stress Test layer\\[2mm] Adversarial targets A–E, the QFT round (F) and the seam round (G): try to break the reductions, not confirm them

[documents/tfpt_5_redteam/preamble.tex:1-37](../../../../documents/tfpt_5_redteam/preamble.tex#L1) · `section:57628e792bbc3f4e2f1c` · frontmatter

-1.2cm TFPT — Red Team / Stress Test layer\\[2mm] Adversarial targets A–E, the QFT round (F) and the seam round (G): try to break the reductions, not confirm them

## Target F — the perturbative 4D-QFT [formula in source] scale round (v269–v275)

[documents/tfpt_5_redteam/sections/001-section-method.tex:186-266](../../../../documents/tfpt_5_redteam/sections/001-section-method.tex#L186) · `section:c00f76407fefd9782fd7` · section

Minimal claim. The round published in release 5.2 — the Epstein–Glaser [formula in source] layer ( v269/ v271/ v273), the PMNS Jarlskog CP strength ( v270), the absolute [formula in source] -scale typing ( v272), the scale over-determination ( v274) and the QG.AMB.01 roadmap ( v275). redteam/rt_F_qft4d.py

Explicit source references:

- [verification/v304_idg_nonlocal_ghost.py](../../../../verification/v304_idg_nonlocal_ghost.py#L1) · explicit_script_identifier
- [verification/redteam/rt_F_qft4d.py](../../../../verification/redteam/rt_F_qft4d.py#L1) · explicit_file_citation
- [GRAV.KMS.HESSIAN.01](../../../../verification/status_ledger.csv#L93) · explicit_claim_identifier · original status at linked ledger row
- [verification/v335_seam_equiv_unify.py](../../../../verification/v335_seam_equiv_unify.py#L1) · explicit_script_identifier
- [QFT4D.OS.RECON.01](../../../../verification/status_ledger.csv#L745) · explicit_claim_identifier · original status at linked ledger row
- [verification/v270_pmns_jarlskog_assembly.py](../../../../verification/v270_pmns_jarlskog_assembly.py#L1) · explicit_script_identifier
- [verification/v481_seesaw_carrier_ladder.py](../../../../verification/v481_seesaw_carrier_ladder.py#L1) · explicit_script_identifier
- [verification/v385_uvbranch_killtest.py](../../../../verification/v385_uvbranch_killtest.py#L1) · explicit_file_citation
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_script_identifier
- [LAMBDA.BRANCH.01](../../../../verification/status_ledger.csv#L129) · explicit_claim_identifier · original status at linked ledger row
- [verification/v265_qft4d_fork_freeze.py](../../../../verification/v265_qft4d_fork_freeze.py#L1) · explicit_script_identifier
- [verification/v381_qft4d_eg_allorder.py](../../../../verification/v381_qft4d_eg_allorder.py#L1) · explicit_script_identifier
- [verification/v482_seesaw_rung_decision.py](../../../../verification/v482_seesaw_rung_decision.py#L1) · explicit_script_identifier
- [verification/v272_nu_mass_scale.py](../../../../verification/v272_nu_mass_scale.py#L1) · explicit_script_identifier
- [verification/v274_scale_overdetermination.py](../../../../verification/v274_scale_overdetermination.py#L1) · explicit_script_identifier
- [QFT4D.EG.ALLORDER.01](../../../../verification/status_ledger.csv#L444) · explicit_claim_identifier · original status at linked ledger row
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v370_grav_spin2_unitarity.py](../../../../verification/v370_grav_spin2_unitarity.py#L1) · explicit_script_identifier
- [verification/v271_eg_oneloop_quartic.py](../../../../verification/v271_eg_oneloop_quartic.py#L1) · explicit_script_identifier
- [verification/v275_qgamb_roadmap.py](../../../../verification/v275_qgamb_roadmap.py#L1) · explicit_script_identifier
- [GRAV.SPIN2.UNITARITY.01](../../../../verification/status_ledger.csv#L80) · explicit_claim_identifier · original status at linked ledger row
- [verification/v269_spert_paqft_skeleton.py](../../../../verification/v269_spert_paqft_skeleton.py#L1) · explicit_script_identifier
- [verification/v380_grav_kms_hessian.py](../../../../verification/v380_grav_kms_hessian.py#L1) · explicit_script_identifier
- [QGAMB.IDG.01](../../../../verification/status_ledger.csv#L400) · explicit_claim_identifier · original status at linked ledger row
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row
- [verification/v334_conformal_resolution.py](../../../../verification/v334_conformal_resolution.py#L1) · explicit_script_identifier
- [verification/v273_eg_gauge_running.py](../../../../verification/v273_eg_gauge_running.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [(* ==== v271: concrete Epstein-Glaser one-loop quartic renormalization (QFT4D.SPERT.02) ==== v269 (S_pert skeleton) + v270 (PMNS Jarlskog, numerical) are Python-only; v271's exact EG core (scaling degree, counterterm count, loop factor) is mirrored here. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2083)
- [(* ==== v273: EG one-loop gauge self-energy -> (b1,b2,b3) (QFT4D.SPERT.03) ==== v272 (nu-scale, numerical) + v274 (over-determination, numerical) + v275 (QG roadmap) are Python-only; v273's exact one-loop beta coefficients from the content are mirrored. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2103)

## Was der funktionierende Rest über den Ursprung aussagt

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/001-section-korrektur-der-suchrichtung-am-29-september.tex:93-107](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/001-section-korrektur-der-suchrichtung-am-29-september.tex#L93) · `section:7fe2b9bd33d71d6cdb49` · subsection

Die Trennung ``alles Erfolgreiche ist invariant, alles Dynamische scheitert'' ist am Original zu grob. Vorhanden sind exakte Transportiteration, positive Markov-Ausführung, unitäre Wortoperationen, KMS-Rekonstruktion sowie in dieser Sitzung analytische lokale Quellgrenzwerte. Ihre gemeinsame Herkunft und physische Identifikation bleibt zu zeigen. Das ist eine andere Aufgabe als die Existenz irgendeiner Dynamik.

Explicit source references:

- [experiments/theory-contracts/compiler-origin-audit-20260913/label_lift.py](../../../../experiments/theory-contracts/compiler-origin-audit-20260913/label_lift.py#L1) · explicit_file_citation
- [experiments/theory-contracts/seam-state-dynamics-injectivity-20260918/README.md](../../../../experiments/theory-contracts/seam-state-dynamics-injectivity-20260918/README.md#L1) · explicit_file_citation
- [experiments/theory-contracts/universalraum-inversion-20260914/checker.py](../../../../experiments/theory-contracts/universalraum-inversion-20260914/checker.py#L1) · explicit_file_citation

## Wie dieses Dokument zu lesen ist

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/001-section-korrektur-der-suchrichtung-am-29-september.tex:282-302](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/001-section-korrektur-der-suchrichtung-am-29-september.tex#L282) · `section:adf159777928dff103cb` · section

tocsectionWie dieses Dokument zu lesen ist Der Auftrag lautet, die vollständige Verbindung zu finden. Deshalb ist der Maßstab dieses Textes eine gemeinsame Herleitung mit denselben Feldern, demselben Zustand und derselben Dynamik. Eine Sammlung passender Zahlen, ein funktionierendes Hilfsmodell oder eine große Zahl grüner Prüfungen reicht dafür nicht.

## Chirale Diagnostik: bestätigte Literatur, ausdrücklich bedingte Übertragung

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/009-subsection-zustand-und-wirkung-aus-derselben-r-ckwirkung.tex:329-361](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/009-subsection-zustand-und-wirkung-aus-derselben-r-ckwirkung.tex#L329) · `section:a4621408bc230f19ea9a` · subsection

Für einen vollrangigen exakten Quanten-Markovzustand gilt Der letzte Schritt folgt durch zyklische Spur und [formula in source] . Eine nichtverschwindende solche Antwort darf daher nicht durch vollständige Entfernung der betreffenden räumlichen bedingten Korrelationen ausgelöscht werden. Dies verbietet keine zeitliche Markov- oder GKLS-Dynamik.

## Der stärkste neue Selektor verbindet Orientierung und gesamten Randinhalt

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/010-subsection-der-st-rkste-neue-selektor-verbindet-orientie.tex:1-25](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/010-subsection-der-st-rkste-neue-selektor-verbindet-orientie.tex#L1) · `section:98a4a2e5a2b1c9c2845e` · subsection

Für dieselbe geeignete Bulk-/Randquelle sei der chirale Zusammenhang [formula in source] begründet. Für die Vakuum-Intervallentropie ihrer Rand-CFT auf der Linie gilt bei festem Regulator- und Lorentzrahmen castroentropy Dann ergibt sich Die Schranke folgt aus Unitarität, [formula in source] . Ihre Sättigung schließt gegenläufige gaplose CFT-Beiträge aus. Insbesondere Mit tatsächlich vorhandenen unitären affinen E8-S

## Existenzsätze erhalten, ihre Reichweite nicht ausweiten

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/011-subsection-bei-gr-erer-tiefe-reichen-vorzeichen-nicht-au.tex:140-144](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/011-subsection-bei-gr-erer-tiefe-reichen-vorzeichen-nicht-au.tex#L140) · `section:1b7b8367dc194d85c806` · subsection

Die ältere thermodynamische Dynamikarbeit enthält bereits Lieb–Robinson-Schranken, Zeitentwicklung, Gauss-Kompatibilität sowie Grundzustands-/KMS-Existenz für die dort spezifizierte lokale Hamiltonklasse mit endlichen Faktoren auf [formula in source] . Diese Ergebnisse werden hier nicht erneut als ungelöste reine Existenzfragen behandelt. Sie begründen aber auch nicht von selbst die ausgewählte räumliche [formula in

## Warum eine endliche Projektion Gedächtnis verliert

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/012-section-die-quelle-bewegt-sich-wirklich-eine-exakte-zeit.tex:28-40](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/012-section-die-quelle-bewegt-sich-wirklich-eine-exakte-zeit.tex#L28) · `section:db796bce86ab27e87184` · subsection

Für den Projektor [formula in source] auf die anfänglichen Grad-eins-Ströme folgt Bei [formula in source] ist diese Amplitude [formula in source] . Die ununterbrochene Bewegung bis [formula in source] hat Amplitude [formula in source] , während zweimaliges Projizieren [formula in source] ergibt. Die scheinbare endliche Dynamik komponiert somit nicht wie der volle Prozess. Sie hat Anteile entfernt, die später zurückwi

## Der ganze Quantenkanal und die gemeinsame Ereignisregel

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/015-subsection-der-uniforme-sewing-test-wird-zu-einem-exakt-.tex:74-120](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/015-subsection-der-uniforme-sewing-test-wird-zu-einem-exakt-.tex#L74) · `section:7300c871735b76c1b717` · subsection

Auf der tatsächlichen markierten Zweierebene gilt in der Pauli-Reihenfolge [formula in source] Die ursprüngliche Trine erfüllt [formula in source] . Am historischen, weiterhin bedingt ausgewählten Punkt [formula in source] ist Der partiell transponierte normierte Choioperator hat den Eigenwert [formula in source] . Dieser Kanal kann also Verschränkung mit einer Referenz erhalten. Die dreidimensionale vollständige Dep

## Baustein 2: Ein aktionsabhängiger Gedächtniscompiler

[documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/016-subsection-die-zur-ckgewonnene-darstellung-ist-der-urspr.tex:296-328](../../../../documents/tfpt_explorer--docs--session-2026-09-28--TFPT_Gesamtsynthese_2026-09-28/sections/016-subsection-die-zur-ckgewonnene-darstellung-ist-der-urspr.tex#L296) · `section:11da13aac90475d2e0cd` · subsection

Für bekannte endliche lineare Dynamiken sind die Bedingungen konkret. Ein Encoder [formula in source] erhält sämtliche Antworten, wenn für jede erlaubte Aktion [formula in source] Dann folgt durch Induktion für jedes endliche Wort. Dies ist ein wirklicher Vertrag über die Zukunft, keine Ähnlichkeit nur des momentanen Zustands.

## Front matter / support

[documents/tfpt_horizon_readouts/sections/001-section-the-seam-constant-is-the-universal-horizon-tempe.tex:1-63](../../../../documents/tfpt_horizon_readouts/sections/001-section-the-seam-constant-is-the-universal-horizon-tempe.tex#L1) · `section:0d5c0120adefeb817d59` · frontmatter

A change of bookkeeping, not new gravitational physics: if gravity is the geometry-channel readout of the seam, then all horizons read the same boundary constant [formula in source] . This note collects the readouts — Hawking, de Sitter and Unruh temperature, black-hole thermodynamics, the Page time, scrambling, the Nariai bound, [formula in source] and cosmic birefringence — in seam units, with two genuine compiler

Explicit source references:

- [verification/v526_seam_thermal_kms_nariai_bridge.py](../../../../verification/v526_seam_thermal_kms_nariai_bridge.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ==== v526 round: SEAM.THERMAL.KMS.01 -- the thermal seam / third leg of c3. Exact mirrors: (i) Q-polynomial / measured beta = N at N = 8; (ii) beta_angle = 2 pi = 1/(4 c3) chain + Hawking/Nariai identities; (iii) nu char-poly / closed modular spectrum at N = 8. Tomita 40-digit residuals and the entropy tower are Python-only. Mirrors v526. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8141)

## Schwarzschild thermodynamics in four [formula in source] -lines

[documents/tfpt_horizon_readouts/sections/001-section-the-seam-constant-is-the-universal-horizon-tempe.tex:74-122](../../../../documents/tfpt_horizon_readouts/sections/001-section-the-seam-constant-is-the-universal-horizon-tempe.tex#L74) · `section:6a9749089feb7aba080e` · section

For a non-rotating, uncharged hole in Planck units ( [formula in source] ): [formula in source] Two of the constants are compiler numbers, not coincidences: the Hawking temperature factor is [formula in source] — the same P1 boundary normaliser; the radiated-power denominator is [formula in source] , the Weyl-group order of the [formula in source] carrier. Scope of [formula in source] : this is the idealised massless

Explicit source references:

- [verification/v28_gravity_fR.py](../../../../verification/v28_gravity_fR.py#L1) · explicit_file_citation
- [verification/v208_bh_thermodynamics.py](../../../../verification/v208_bh_thermodynamics.py#L1) · explicit_file_citation
- [verification/v198_modular_commutator_reduction.py](../../../../verification/v198_modular_commutator_reduction.py#L1) · explicit_file_citation
- [SEAM.THERMAL.KMS.01](../../../../verification/status_ledger.csv#L593) · explicit_claim_identifier · original status at linked ledger row
- [verification/v526_seam_thermal_kms_nariai_bridge.py](../../../../verification/v526_seam_thermal_kms_nariai_bridge.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamDeckClosure.lean#L1)
- [(* ==== v526 round: SEAM.THERMAL.KMS.01 -- the thermal seam / third leg of c3. Exact mirrors: (i) Q-polynomial / measured beta = N at N = 8; (ii) beta_angle = 2 pi = 1/(4 c3) chain + Hawking/Nariai identities; (iii) nu char-poly / closed modular spectrum at N = 8. Tomita 40-digit residuals and the entropy tower are Python-only. Mirrors v526. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8141)
- [(* ---- (v208) HOR.BHTHERMO.01: scalaron Wald factor + modular 2 pi = 1/(4 c3) ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1717)
- [(* ---- (v198) QGEO.MODULAR.01: principal symbol |k| commutes with the clock rho EXACTLY ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1639)

## De Sitter and the Nariai bound from the [formula in source] closure

[documents/tfpt_horizon_readouts/sections/001-section-the-seam-constant-is-the-universal-horizon-tempe.tex:153-168](../../../../documents/tfpt_horizon_readouts/sections/001-section-the-seam-constant-is-the-universal-horizon-tempe.tex#L153) · `section:dfadaead644b49d4cd45` · section

With the action-grammar Hubble/ [formula in source] readouts ( [formula in source] , [formula in source] ), the asymptotic de Sitter horizon gives [formula in source] This is the area-law horizon capacity of the universe from [formula in source] alone. The maximal black hole in de Sitter (Nariai) is [formula in source] the absolute gravitational-capacity limit, not an observed structure.

## Appendix: computational verification

[documents/tfpt_horizon_readouts/sections/002-section-turnaround-radius-gravitational-waves-light.tex:342-353](../../../../documents/tfpt_horizon_readouts/sections/002-section-turnaround-radius-gravitational-waves-light.tex#L342) · `section:f14cd99bec2b5b92a146` · section

The seam-unit identities of this note (the thermal grammar [formula in source] , the power denominator [formula in source] , the de Sitter entropy form, the Page fraction, the scrambling [formula in source] , and [formula in source] ) are re-derived and machine-checked in verification/v8_horizon.py. Run python verification/v8_horizon.py (needs mpmath); it imports the same [formula in source] primitive as the rest of

Explicit source references:

- [verification/v8_horizon.py](../../../../verification/v8_horizon.py#L1) · explicit_file_citation

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [experiments/lean4-carrier-rigidity/TfptCarrier/HiggsTopForm.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/HiggsTopForm.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CalderonProjector.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CalderonProjector.lean#L1)

## Front matter / support

[documents/tfpt_research_contracts/sections/001-source.tex:1-160](../../../../documents/tfpt_research_contracts/sections/001-source.tex#L1) · `section:8d1e955857da6655fadd` · frontmatter

After the compiler closure the entire residual is [formula in source] This note turns the two genuine research gates into contracts: a numbered chain of lemmas, the single theorem that closes each gate, and — for every step — whether it is machine-certifiable today. [formula in source] (Koide, [formula in source] , axion relic scale, [formula in source] ) is not a gate: those are deliberately typed QFT/cosmology inte

Explicit source references:

- [verification/v519_woit_theta_rp_free.py](../../../../verification/v519_woit_theta_rp_free.py#L1) · explicit_script_identifier
- [verification/v384_residual_certification.py](../../../../verification/v384_residual_certification.py#L1) · explicit_file_citation
- [RES.COXETER.SYMMETRY.01](../../../../verification/status_ledger.csv#L465) · explicit_claim_identifier · original status at linked ledger row
- [verification/v508_celestial_wp5e_gamma_pairing.py](../../../../verification/v508_celestial_wp5e_gamma_pairing.py#L1) · explicit_script_identifier
- [verification/v374_ftransfer_qcd.py](../../../../verification/v374_ftransfer_qcd.py#L1) · explicit_script_identifier
- [verification/v89_carrier_index_lemma.py](../../../../verification/v89_carrier_index_lemma.py#L1) · explicit_script_identifier
- [verification/v529_seam_interacting_toy_fk.py](../../../../verification/v529_seam_interacting_toy_fk.py#L1) · explicit_script_identifier
- [QGEO.SYM.01](../../../../verification/status_ledger.csv#L261) · explicit_claim_identifier · original status at linked ledger row
- [verification/v518_celestial_wp5e_delta1_derived_kill.py](../../../../verification/v518_celestial_wp5e_delta1_derived_kill.py#L1) · explicit_script_identifier
- [verification/v383_dynamics_universal.py](../../../../verification/v383_dynamics_universal.py#L1) · explicit_script_identifier
- [verification/v113_quasifree_kernel.py](../../../../verification/v113_quasifree_kernel.py#L1) · explicit_file_citation
- [verification/v517_celestial_wp5e_m3_a0_uplift.py](../../../../verification/v517_celestial_wp5e_m3_a0_uplift.py#L1) · explicit_script_identifier
- [verification/v258_dirac_covariance_induction.py](../../../../verification/v258_dirac_covariance_induction.py#L1) · explicit_script_identifier
- [verification/v533_ftransfer_disc7_norm.py](../../../../verification/v533_ftransfer_disc7_norm.py#L1) · explicit_file_citation
- [GATE.METRIC.05](../../../../verification/status_ledger.csv#L161) · explicit_claim_identifier · original status at linked ledger row
- [verification/v224_diamond_ftransfer_path.py](../../../../verification/v224_diamond_ftransfer_path.py#L1) · explicit_file_citation
- [FORM.COXETER.PRIME2.01](../../../../verification/status_ledger.csv#L466) · explicit_claim_identifier · original status at linked ledger row
- [verification/v75_upoint_to_vgeo.py](../../../../verification/v75_upoint_to_vgeo.py#L1) · explicit_script_identifier
- [RESIDUAL.CERTIFICATION.01](../../../../verification/status_ledger.csv#L447) · explicit_claim_identifier · original status at linked ledger row
- [verification/v261_modular_spectral_closure.py](../../../../verification/v261_modular_spectral_closure.py#L1) · explicit_script_identifier
- [verification/v71_simple_r_bridge.py](../../../../verification/v71_simple_r_bridge.py#L1) · explicit_script_identifier
- [verification/v336_continuum_limit.py](../../../../verification/v336_continuum_limit.py#L1) · explicit_script_identifier
- [TFPT.TOE.COMPLETE.01](../../../../verification/status_ledger.csv#L764) · explicit_claim_identifier · original status at linked ledger row
- [verification/v515_celestial_wp5e_m1_omega_n.py](../../../../verification/v515_celestial_wp5e_m1_omega_n.py#L1) · explicit_script_identifier
- [verification/v49_readout_rigidity.py](../../../../verification/v49_readout_rigidity.py#L1) · explicit_script_identifier
- [verification/v382_alpha_quillen_exact.py](../../../../verification/v382_alpha_quillen_exact.py#L1) · explicit_script_identifier
- [verification/v565_woit_beta3_pt_duality.py](../../../../verification/v565_woit_beta3_pt_duality.py#L1) · explicit_script_identifier
- [verification/v524_woit_beta2_os_quotient.py](../../../../verification/v524_woit_beta2_os_quotient.py#L1) · explicit_script_identifier
- [verification/v259_modular_cutoff_kappa.py](../../../../verification/v259_modular_cutoff_kappa.py#L1) · explicit_script_identifier
- [verification/v520_celestial_wp5e_measure_decision.py](../../../../verification/v520_celestial_wp5e_measure_decision.py#L1) · explicit_script_identifier
- [verification/v522_woit_beta1_gso_gauge.py](../../../../verification/v522_woit_beta1_gso_gauge.py#L1) · explicit_script_identifier
- [verification/v92_glue_uniqueness.py](../../../../verification/v92_glue_uniqueness.py#L1) · explicit_script_identifier
- [verification/v472_quillen_detline_moduli.py](../../../../verification/v472_quillen_detline_moduli.py#L1) · explicit_script_identifier
- [verification/v153_no_unit_theorem.py](../../../../verification/v153_no_unit_theorem.py#L1) · explicit_script_identifier
- [CONTRACT.F.01](../../../../verification/status_ledger.csv#L302) · explicit_claim_identifier · original status at linked ledger row
- [verification/v523_celestial_wp5e_wm_derived.py](../../../../verification/v523_celestial_wp5e_wm_derived.py#L1) · explicit_script_identifier
- [verification/v394_coxeter_atoms.py](../../../../verification/v394_coxeter_atoms.py#L1) · explicit_file_citation
- [ALPHA.QUILLEN.EXACT.01](../../../../verification/status_ledger.csv#L445) · explicit_claim_identifier · original status at linked ledger row
- [verification/v371_ftransfer_pole.py](../../../../verification/v371_ftransfer_pole.py#L1) · explicit_script_identifier
- [verification/v632_ftransfer_pgl2.py](../../../../verification/v632_ftransfer_pgl2.py#L1) · explicit_file_citation
- [verification/v187_ftransfer_laws.py](../../../../verification/v187_ftransfer_laws.py#L1) · explicit_script_identifier
- [verification/v68_seeley_dewitt_residual.py](../../../../verification/v68_seeley_dewitt_residual.py#L1) · explicit_script_identifier
- [CELEST.SEAM.01](../../../../verification/status_ledger.csv#L556) · explicit_claim_identifier · original status at linked ledger row
- [verification/v409_coxeter_prime2_lemma.py](../../../../verification/v409_coxeter_prime2_lemma.py#L1) · explicit_file_citation
- [WOIT.OS.TWISTOR.01](../../../../verification/status_ledger.csv#L585) · explicit_claim_identifier · original status at linked ledger row
- [verification/v469_seam_crossedproduct_route.py](../../../../verification/v469_seam_crossedproduct_route.py#L1) · explicit_script_identifier
- [SEAM.EQUIV.01](../../../../verification/status_ledger.csv#L381) · explicit_claim_identifier · original status at linked ledger row
- [verification/v973_seam_route_narrowing.py](../../../../verification/v973_seam_route_narrowing.py#L1) · explicit_file_citation
- [verification/v390_prime2_facet.py](../../../../verification/v390_prime2_facet.py#L1) · explicit_script_identifier
- [verification/v364_vgeo_sharpen.py](../../../../verification/v364_vgeo_sharpen.py#L1) · explicit_script_identifier
- [verification/v369_qgamb_redundancy.py](../../../../verification/v369_qgamb_redundancy.py#L1) · explicit_script_identifier
- [verification/v534_seam_straddle_cone.py](../../../../verification/v534_seam_straddle_cone.py#L1) · explicit_script_identifier
- [verification/v972_seam_interaction_front.py](../../../../verification/v972_seam_interaction_front.py#L1) · explicit_file_citation
- [experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean#L1) · explicit_file_citation
- [verification/v260_k3_kummer_unification.py](../../../../verification/v260_k3_kummer_unification.py#L1) · explicit_script_identifier
- [verification/v470_alpha_inflow_level.py](../../../../verification/v470_alpha_inflow_level.py#L1) · explicit_script_identifier
- [QG.AMB.01](../../../../verification/status_ledger.csv#L47) · explicit_claim_identifier · original status at linked ledger row

Navigation joins (both sources name the same Python identifier; statement equivalence is not inferred):

- [(* ---- (v113) quasi-free kernel ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L608)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditContract.lean#L1)
- [(* ==== v522 round: WOIT.BETA1.GSO.01 -- "the clock is time-like, GSO is the gauge datum": the Z8/NS tower + Ramond control, the invariant- dimension chain, the time-like reflection census, the exact Hermiticity witness of the clock insertion, the GSO 2-torsion census and the family-A tower-normalisation failure. The 40-digit GSO-fixed RP inertia certificates ((29,0,0)/(8,0,0)/(7,9,6)/(17,12,0)) are Python-only by the suite convention. Mirrors v522. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L7610)
- [(* ---- (v49) Readout Rigidity: c_u/c_d constant on the selector stratum ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L202)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AuditCheck.lean#L1)
- [(* ---- (v153) No-Unit Theorem ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1369)
- [(* ---- (v68) Seeley-DeWitt a2 for the carrier; 1/G UV-sensitive => k=c3/2 is a normalization ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L322)
- [(* ---- (v92) glue uniqueness ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L151)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SpectralGapAttractor.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/CoxeterPrime2.lean#L1)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean#L1)
- [experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean](../../../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105J/shorttime20261105J_lean/SeamScalingLimitCurrent.lean#L1)
- [(* ==== v534 round: SEAM.STRADDLE.CONE.01 -- the RP cone on the interacting FK quartets (the first dynamical selection of the alignment bit). Exact mirrors: the stabilizer census of the mark family, the pi/2 partition + straddle incidence, the reflection-closure signs of the straddling quartets, and the equivariant Fix rays (uniform on m=2..6, NS-wrap-twisted on m=1,7). The 256-dim Fock spectra, the PT Gram corrections and the finite-g survival scan are Python-only (numpy), flagged per GATE.WOLFRAM.02. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8402)
- [(* ==== v524 round: WOIT.BETA2.OS.01 -- "the OS quotient of the free system made explicit" (verdict SUCCESS, [C]-typed per contract precision (iii)). Four exact mirrors: (i) the euclidean-rotation preconditions (alpha_1^N = (-1)^F wrap arithmetic at N = 16 and N = 8) + the exact PD certificate identity sin^2(3pi/8) - sin(pi/8) sin(5pi/8) = 1/2; (ii) the Klein-Landau one-particle transfer forms tau_k on the shrinking domains at N = 16, k = 1..7: exactly Hermitian, odd k with exactly ZERO diagonal (the C(even) = 0 checkerboard), even k = 2j obeying the exact square identity tau_2j(a, b) = <alpha_j a, alpha_j b> entrywise; (iii) the clock identities (C(9) = C(7), C(11) = C(5), the transparency product formula sin(7pi/16) - sin(5pi/16) = 2 cos(3pi/8) sin(pi/16) > 0) + the N = 8 compressed clock spectrum EXACTLY {1, sqrt(2)-1} = {1, 1/delta_S} via the generalized eigenvalues of (tau_2, G) on the one-particle domain, with sin(pi/8)/sin(3pi/8) = sqrt(2)-1; (iv) the pre-declared KMS witnesses (C(1)/C(3) = 1 + sqrt(2) = delta_S at N = 8; -(C5-C7)((C3-C7)+(C5-C7)) < 0 exactly at N = 16, i.e. det(G - tau_4) < 0 on the odd one-particle clock block). The 37x37/16x16 H_phys Grams, the 40-digit inertia certificates, the spectral projections and the reconstructed rotation group are numerical (mpmath eighe), Python-only, flagged here. Mirrors v524. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L7966)
- [(* v469 (ii): KLM mu-index + the 16-fold-way pin *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2937)
- [(* ==== v518 round: CELEST.WP5E.DELTA1.01 -- "the delta-1 chain decided: the derived chiral measure kills the delta-1 route" (delta-1b/1c/1d consolidated; an honest, decided NEGATIVE result with the v516 tension stated). Six exact mirrors: (i) the 16-component Weil system of the D5 (+) A3 discriminant module Z4 x Z4 with q = (5x^2 + 3y^2)/8 -- Gauss sums 2 zeta_8^5 x 2 zeta_8^3 = 4 = sqrt(16) (signature 0 mod 8 = rank E8), S symmetric + unitary, S^2 = C, (ST)^3 = S^2, S^4 = 1, T^8 = 1, tensor split S = S_D5 (x) S_A3 with factor signatures e(-5/8), e(-3/8); (ii) the two invariant Lagrangians (diagonal H and anti-diagonal H', both S- and T-fixed = two E8 gluings), 6 isotropic classes, P_H S P_H = (1/4) J and 3/4 of every S-image weight on the 12 shifted classes (why the naive 4-character rule fails); (iii) the KOBOUNDARY TEST on the recorded exact transport phases (t_p, s_p) of the dressed gauge blocks M = G[a,b] eta^{-8}: the SL(2,Z) relation defects c_{S^4}, c_{(ST)^3 S^-2}, c_{[S^2,T]} are (1, 1, 1) EXACTLY on all 15 sector pairs -- the mu4 multiplier system is a CHARACTER of the orbit stabilisers, not a genuine 2-cocycle; (iv) the Gamma_1(4) CHARACTER FORMULA lambda(gamma) = i^{2B + C/4} on the recorded loop matrices (all in Gamma_1(4), exact integer check), with all four equivalent monomial solutions agreeing; (v) the T-FIX MECHANISM: v502 vacuum energies E_b exact, the f1 f3 dressing shifts every T-fixed node by e(2E + 2E) = e(-1/6) and (-1) x e(-1/6) = e(1/3) = chi_4(T) UNIFORMLY, while f1f2f3 demands {chi_3, chi_4} and f2 demands {chi_5, chi_6} (inconsistent: only the twisted fibre block supplies one character); (vi) the misassignment controls (wrong form (x^2+y^2)/4: |Gauss|^2 = 64 != 16 and S^2 != C exactly; wrong signature sigma = e(-1/4): S^2 = e(-1/2) C exactly). The Harvey-Moore integrals, the SVD contraction solves, the theta covariance certificates and the block transport measurements are numerical (mpmath 30+ digits / numpy) and stay Python-only, flagged here. Mirrors v518. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L6877)
- [(* ---- (v75) Gate 1: U_point -> v_geo. lepton product 32/9; (ratios,product)<=>(individuals) bijection ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L366)
- [(* ==== v533 round: FTR.DISC7.NORM.01 -- the disc -7 norm line. Exact mirrors: norm identity, translation, norms (2,4,14,32), constant discriminant -7, vertex, negative controls. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8373)
- [(* ---- (v89) carrier index lemma ---- *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L118)
- [(* ==== v260: one Kummer/K3 carries seam + carrier-16 + E8 (ARCH.K3.01) ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2033)
- [experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean#L1)
- [(* ==== v470 round: the alpha^3 level / embedding-index rigidity -- the exact k_Y = 5/3 arithmetic (tr Y^2 / tr T3^2 over the 5bar), the b1 conversion (3/5)(41/6) = 41/10, and the carrier decomposition 41/6 = 20/3 + 1/6; mirrors v470 (the FHS Chern and the alpha^-1 root re-verification are numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2943)
- [(* ==== v517 round: CELEST.WP5E.M3.01 -- "the a0 uplift": milestone M3 of the CELEST.SEAM.01 back-reaction programme executed (SUCCESS on the preregistered v514 S8.3 criterion). Six exact mirrors: (i) the GH multipole ledger (selection rule m = 0 mod 4, (4, +-4) amplitude (35/16) sin^4 th) + the twistor family at general orbit phase (e1 = e3 = 0, a2 = 4 z0 zb0 zeta^2, a0(zeta), v515 slice), (ii) the kernel bridge (null coordinate: any kernel harmonic; residue V-ledger identity 1/(2r); source flux -4 pi per centre), (iii) the asymptotic kernel log (coefficient 4 = |mu4|, first seam correction a0/eta^4, exact m-grading, GLT tower p_{4k} = 4(-a0)^k with the n = 0 mod 4 selection rule, no log x power terms), (iv) the exceptional-locus log (chi(0) = log a0 = 4 log t0 + i(4 phi0 + pi), modulus response 1 at the locus vs power-law at infinity), (v) the period response (linear in t0, d log Pi/d log a0 = 1/4 = 1/|mu4|, integrated monodromy i = one Coxeter clock step, centre permutation, a0-rigid node support, K4 discriminants ~ t0^2), and (vi) the negative controls ((4,0) clock-invariant, Z2/EH reads 2 on every dial, k = 3/5 orbits move the coefficient, the forbidden family fails both dials). The contour-transform typing itself is the [C] GLT dictionary; the BM kernel anchor is already mirrored in the v515 round. Mirrors v517. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L6679)
- [(* ---- (v71) SIMPLE R-bridge: selector stratum derived => ratios are integer Plucker, no monodromy ---- *)](../../../../verification/wolfram/tfpt_readouts.wl#L343)
- [(* ==== v529 round: SEAM.INT.FKTOY.01 -- interacting FK toy / Kill-Test-2 shadow. Exact mirrors: (i) Z8 integer identities; (ii) straddle combinatorics; (iii) FK binomial ladder multiplicities. 256-dim spectra Python-only. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L8246)
- [(* ==== v523 round: CELEST.WP5E.WM.01 -- "the constructive w_m derivation": the residual [O] of v516/v520 closed at the fibre-zero-mode level (verdict ERFOLG; the typed premises TP-REG/TP-Q/TP-NUM/TP-CH are the [C] fence, a typing, not mirrorable). Five exact mirrors: (i) the equivariant mode ledger = the closed form 1/((1-q i^j)(1-q i^{-j})) at every order n <= 40 with the Abel value (1/2, 1/4, 1/2) = 1/det_j, Dedekind sum 5/4 and product 16; (ii) the zeta-determinant route (reflection formula symbolic, det_zeta = 4 sin^2(pi j/4) = det_j, Quillen split det Delta = det_j^2, the unique REAL POSITIVE section, the diag(i,i) control breaking realness and the conjugation pairing); (iii) the derived-weight contact vectors + per-sector Okubo squares + total 45<x,x>^2 + psi certificates on the 240 glue roots + the Z2/EH anchor; (iv) the consistency chain w = 4h = -4 ch2 with the T3 budget c = 4 and the two wrong-weight negative controls (1/det^2 leftovers (2, 12) + Quillen contradiction; 1/|1-i^j| irrational with T5/T3 leftovers in Q(sqrt2)); (v) the k = 3, 5 cyclotomic wandering controls w_m = m(k-m)/2 = k h_m. The (C,2) Cesaro truncation means, the Hurwitz/Lerch 40-digit certificates and the delta-1f stripped block limits are numerical (mpmath), Python-only, flagged here. Mirrors v523. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L7793)
- [(* ==== v259: the modular/KMS spectral-action cutoff fixes f_2/f_0 = 1 (PS.SPECACT.02) ==== v258 (Dirac as covariance induction) is numerical matrix-log linear algebra -> Python-only. *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2014)
- [(* ==== v469 round: the net-level crossed-product certification -- the LR locality integers (h_s(SO(16)_1) = 16/16 = 1 in Z; the mu4 glue h(J^k) = {1,1,1} all integer), the KLM mu-index arithmetic (4/2^2 = 1 and 16/4^2 = 1, both holomorphic) and the Kitaev 16-fold-way pin (nu = 16 x |C| = 16 = 0 mod 16, c_- = nu/2 = 8 = g_car+N_fam); mirrors v469 (the FHS Chern itself is numerical, Python-only). ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2924)
- [experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean#L1)
- [(* v224 diamond F_transfer path: sheet axis curved, winding axis flat, Plucker steps *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L1849)
- [(* ==== v520 round: CELEST.WP5E.MEASURE.01 -- "the WP5e BCOV measure question decided at probe level" (delta-1e + delta-1f consolidated): the exact Weil data, the invariant-subspace forcing, the character arithmetic behind the two constructive sources (F-independence + Quillen pairing), the completion/contact/psi and Z2-anchor targets, the level-1 4-design of BOTH Lagrangian gluings, and the exact negative controls. The Harvey-Moore integrals, SVD memberships and J-weighted closures are numerical (mpmath/numpy) and stay Python-only by the suite convention. Mirrors v520. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L7271)
- [(* ==== v508 round: CELEST.WP5E.GAMMA.01 -- WP5e-gamma "the sphere-axion pairing check", an honest rigid NEGATIVE result (the exchange route for pairing the three twisted sphere axions with the rigid 32 T3 residual of the v505 Atiyah-Bott ledger is killed with a certificate). Exact content mirrored: (i) the VERTEX-SPACE DIMS -- the W(D5) x W(A3)-invariant quadratics on the 8-dim glue Cartan have bidegree dims (1,0,1), total 2 = span{s5, s3}, and the invariant quartics have bidegree dims (2,0,1,0,2), total 5 = span{P1,P2,P3,T5,T3} (Weyl nullspace arithmetic, generator by generator); (ii) the TWISTED QUADRATICS + the K^(0) = -15 K^(2) collapse (the two even-sector quadratics are PARALLEL); (iii) the PRODUCT THEOREM -- (a s5 + b s3)(c s5 + d s3) = ac P1 + (ad+bc) P2 + bd P3 identically (zero T5/T3 content of every exchange image); (iv) the RANK-3 CERTIFICATE -- channels E13 = (16,-96,144,0,0), E22 = (16,32,16,0,0), E00 = 225 E22, exchange matrix rank 2 with kernel 0, rank([M | A_fix]) = 3, annihilator certificate (Phi_T5, Phi_T3, Phi_P)(A_fix) = (0, 32, 72), and the naturalness images (25/4,-11/2,97/4,0,0) / (12,-40,76,0,0) with certificate (0,0); (v) the SO(16) AB SUM (15,-18,-15,-4,40) -- T5 does not even cancel, odd classes empty; (vi) the D8 NATIVE QUARTIC 16 T8 + 12 s8^2 with Killing 28 s8 -- no T3 structure. Mirrors v508. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L4785)
- [(* v469 (i): the LR locality integers of both extension routes *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L2930)
- [experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean](../../../../experiments/lean4-carrier-rigidity/TfptCarrier/AxiomCheck.lean#L1)
- [(* ==== v519 round: WOIT.THETA.FREE.01 -- the WOIT-alpha milestone of the OS twistor bridge: the real structure exists (two-family classification, mu4 torsor, Z8/C^4/Cl(16) levels) and the exact free-RP preconditions + site-cut degeneration + continuum control. The 40-digit positive- definiteness inertia certificates of the bond-cut Grams (R3/R5) are Python-only by the suite convention (mpmath eigh at dps 40). Mirrors v519. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L7078)
- [(* ==== v515 round: CELEST.WP5E.M1.01 -- "the A3 Omega_N": milestone M1 of the CELEST.SEAM.01 back-reaction programme executed (SUCCESS on the preregistered v514 S8.1 criterion). Eight exact mirrors: (i) the residue-form identities (cover pullback 4 = |mu4|, clock phase i, weight ledger, EH anchor 2 = |Z2|), (ii) the twistor family closed (XY = Z^4 + 4 t0^2 la^2 Z^2 - t0^4 (1 - la^4)^2, O(2) reality, seam fibre = v493, self-mirror, CY-compatible clock lift), (iii) the period reduction + lockstep + clock covariance, (iv) the collision ledger (12 nodes exactly on the 8 eighth roots of unity, antipodal double nodes, clock orbits), (v) the BM anchor (2 pi i)^2 + the lens forcing N in 4Z, (vi) the source curve (four lines, K4 connectivity, forced uniform flux), (vii) the conifold nodes (Hessian det 512 t0^6, quadric rank 4, transversality), and (viii) the KILL probe + negative controls (forbidden family 0/24 vs 8/24, twelfth-root support, EH charge 2, diag(i,i) collapse). The exterior-derivative identity dK = 0 and the generic-profile sphere parametrisation stay Python-side in v515 (dict-based exterior algebra). Mirrors v515. ==== *)](../../../../verification/wolfram/tfpt_readouts_extension.wl#L6272)
