# Prüfung des Bestands gegen die Vogelperspektive (Vorbereitung H6h34)

Stand 24. November 2026. Nur lesende Prüfung, keine Rechnung im Repo, kein Commit, keine Atlas-Schreibung. Geprüft wird,
was Papiere, Ledger, Experimente und Atlas zur Architektur aus `docs/research/VOGELPERSPEKTIVE.md` (H6h33) und zu den
Lücken (A)–(G) bereits enthalten.

**Umfang.**
- Papiere: aktive Module unter `documents/` per Volltextsuche (`tfpt_1`–`tfpt_5`, `origin_theory`, `tfpt_horizon_readouts`,
  `tfpt_research_contracts`, `introduction`, `note_e8_gaussian_code`, `tfpt_safeguards`), keine Archive.
- Ledger: `verification/status_ledger.csv` (789 Zeilen), per Regex. Zitiert als „Ledger Z. n (ID)“, n ist die CSV-Zeile.
- Atlas: `semantic_explorer/research/program.json` (Revision `07de10f6082872db9a33da59`, nur gelesen) und die MCP-Suche.
- Experimente: Explorationen vom 20.–24.11., die Verträge vom 14./15.10. und gezielt ältere Verträge.

**Lesehinweis.** Statusangaben stammen aus der Quelle. Atlas-Links und KI-Zusammenfassungen sind Navigation, kein Beweis.
Punkte mit „eigene Einordnung“ sind Folgerungen dieser Prüfung. Sie stehen nicht im Bestand und sind nicht gerechnet,
außer einer kurzen sympy-Probe im Scratchpad (nicht abgelegt).

## 0. Kurzbilanz

1. **Kein direkter Widerspruch.** Kein etabliertes ([E]) oder ausgeschlossenes Ergebnis widerlegt (a)–(e). Es gibt aber
   präzise Spannungen und Bedingungen (Abschnitt 9).
2. **Großer Vorbestand zur Twistor- und Himmelssphären-Lesart.**
   - `CELEST.SEAM.01`, `WOIT.OS.TWISTOR.01` und `SEAM.EQUIV.TWISTOR.01` behandeln die Naht seit Juli 2026 als Himmelssphäre
     bzw. als Gerade im Twistorraum.
   - Deklariert sind dort auch der Twistorraum PT/Γ des A₃-ALE, die selbstduale Einschränkung als Scope-Zaun und sieben
     Kill-Tests.
   - VOGELPERSPEKTIVE §5 („Neu ist die Zuordnung“) sollte diesen Bestand zitieren.
3. **Bereits gerechnet:**
   - τ = i ⇔ J = 0 und Marken = Nullrichtungen eines Spin-2-Quadrupols (v493; Atlas, 23.11.);
   - der Lorentzkegel aus M₂(ℂ) und die Standardform von M₂(ℂ) als E₈-Karte (9.9. und 20.11.);
   - die masselose Wellengleichung aus dem Hardy-Kern auf 2×2-Matrizen (z40, z42).
4. **(A) ist im Bestand bereits als offene Aufgabe registriert:** E8.STANDARDFORM.CHART.01, Forderung D2, und die
   Atlas-Aufgabe `task:fix-nontracial-state-weight-without-parameter-20261120`.
5. **Wichtigste Spannungen:**
   - Die μ₄-Charaktere der Familien {i, −1, −i} sind nicht die von Λ²₊ = sl(2) {1, i, −i}. Die Gleichsetzung gilt nur
     mit dem Vorzeichencharakter von q = dz²/P.
   - Der Modulfluss ist auf dem Dublett eine Drehung, keine Zeittranslation.
   - In v565 ist die Nahtgerade nicht die Twistorgerade eines reellen Punkts.
   - TM1 mit μ-τ-Reflexion gibt δ = ±π/2, eingefroren ist δ = 240°.

## 1. Architektur: Vorarbeiten und Überschneidungen

### 1.1 Naht = Himmelssphäre, Möbius = Lorentz
- `documents/tfpt_3_e8_audit_bootstrap/sections/010-section-three-views-of-one-spine.tex:285-300` — Möb(Ĉ) ≅ PSL(2,ℂ) ≅
  SO⁺(3,1), die Randsphäre der Naht als Himmelssphäre und die Clock z ↦ iz als elliptische Drehung (v180). Dort [O], eine
  strukturelle Lesart. Vorläufer von z62.
- `…/tfpt_3_e8_audit_bootstrap/sections/011-section-the-celestial-and-twistor-continuum-route-label-.tex:1-30` — „The
  celestial and twistor continuum route“: Die Naht wird über die Celestial-Holographie gelesen (BHS auf ALE, Costellos
  twistorielles SDYM). Ein Programm [O]; Atlas `section:b63f8484d81bde8026ec`.
- Ledger Z. 570 (SEAM.CLOCK.RIGIDITY.01, [E]/[C]) — Theorem B spricht von der Himmelssphäre mit eindeutiger
  Spinstruktur (NS auf jedem Kreis), siehe Spannung W7.
- Ledger Z. 737 (CCC.SEAM.KINEMATICS.01, [E] Kinematik) — stereographic(Hopf(v)) = z₁/z₂, Linse S³/ℤ₄, Einstein-Lorentz-
  Brücke. Das ist die Naht als P(ℂ²), also Vorläufer von „Naht = P(V)“. Der Eintrag ist ausdrücklich kinematisch.

### 1.2 Naht im Twistorraum (Lesart T)
- Ledger Z. 556 (CELEST.SEAM.01, Open; WP1–WP5d [E]/[C]/[L]).
  - Die Naht ist der glue-äquivariante ℤ₄-Orbifoldsektor auf dem A₃-ALE ℂ²/ℤ₄.
  - Die Clock ist diag(i,1); sie dreht die Twistorsphäre um π/2.
  - Der twistorielle Bulk ist nur der selbstduale Sektor.
  - Quelle: `documents/tfpt_research_contracts/sections/010-source.tex:44-110` und `:176-200`.
- `…/011-…-route-label-.tex:322-335` (Step 15) — Auf PT′ = Tot(O(1) ⊕ O(1) → ℙ¹) schließt die Penrose-Bilanz π_*O(−2)
  exakt und äquivariant. Das Normalenbündel O(1) ⊕ O(1) ist dort gesetzt, nicht aus P1 hergeleitet.
- Ledger Z. 604 (WOIT.BETA3.DUALITY.01, [E] + [C] deklarierte Konventionen).
  - Clock-verträgliche Dualitäten PT ↔ PT* haben die Signatur (2,2).
  - Im Minkowski-Zweig trifft die Nullmenge die Nahtgerade genau im Nahtkreis |z| = 1; der OS-Schnitt ist erzwungen.
  - Im euklidischen Zweig gibt es keinen Nullpunkt.
  - Quelle: `documents/tfpt_research_contracts/sections/013-source.tex:342-390`.
  - Siehe Spannung W6.
- Ledger Z. 584 (SEAM.EQUIV.TWISTOR.01, Open [O]) — Die Twistor-Route bleibt offen: Costello–Li, globale BCOV+SDYM-Quantisierung
  und OS-Brücke. Auch in `documents/tfpt_1_architecture_e8/sections/005-source.tex:65-66` steht sie als [O].
- Atlas `premise:h6-twistor-seam-reading-20261124` und `finding:h6h33-twistor-seam-dictionary-holomorphic-field-equations-20261124`
  (conditional) — H6h33 ist eingetragen.

### 1.3 τ = i ⇔ J = 0, Marken = Nullrichtungen eines Spin-2-Objekts
- Ledger Z. 557 (CELEST.WP2.01, [E]) und `documents/tfpt_research_contracts/sections/010-source.tex:268-276`.
  - Für w² = Z⁴ + a₀ gilt I = 12a₀, J = 0 identisch, also j = 1728. Das ist τ = i ⇔ J = 0, schon v493.
  - Dieselbe Stelle gibt für die drei aufgelösten Sphären die Charaktere {i, −1, −i}.
- Atlas `finding:pauli-quadrupole-discriminant-square-torus-20261123` (conditional).
  - Für einen spurfreien reellen Quadrupol Q gilt I(F_Q) = 24 tr Q² und J(F_Q) = −576 tr Q³.
  - J = 0 ⇔ tr Q³ = 0 ⇔ Eigenwerte (a, 0, −a) ⇔ j = 1728.
  - Teil 3 von z61 rechnet das neu nach.
- Atlas `finding:tracefree-gaussian-eigenvalue-mode-square-20261123` (conditional) und
  `finding:fourmark-radial-quadrupole-shape-distribution-20261123` (conditional).
  - Der Modus des gaußschen Eigenwertmaßes liegt bei (−a, 0, a).
  - Radial allein gibt aber P(j = 1728) = 0. Erst ein sextischer Term (tr Q³)² würde das Quadrat wählen.
  - Kandidat für eine dynamische Auswahl von J = 0.
- `experiments/exploration/2026-11-23-fockspur-gso/IDEAS.md:24` (F4) — Die Idee „Träger = Sym⁴ℂ² = Weylspinor, Marken =
  Hauptnullrichtungen“ ist hier schon notiert, aber nicht gerechnet. Ihre 3+2-Lesart (Spingewichtsparität) widerspricht
  Z34 (W2).
- `experiments/exploration/2026-11-20-quaternion-frames/VERDICT.md:9,51` (A: KANDIDAT Stufe 1).
  - Der Fünferträger als Quadrupol W = Sym²₀(Im ℍ) war schon da.
  - Gegenprobe 1: Zwischen W und H⁰(O(D)) ist ein D₄-Charaktertwist (r, s) ↦ (−1, +1) nötig, siehe W2.
- `docs/research/TFPT-STAND.md:437-438` — Lean hat E₆(i) = 0. Eindeutigkeit fehlt, weil die Valenzformel fehlt.

### 1.4 Herm(2) = R^{1,3} aus dem Dublett, Standardform
- `experiments/theory-contracts/compiler-cone-object-audit/README.md:9-15,47-53` und `ORDER_PROOF.md:195-219` (9.9.2026).
  - Die maximale Ordnung M ⊂ M₂(ℚ(i)) liefert das E₈-Gitter.
  - Ihre positiven ganzzahligen Punkte bilden den Zukunftslichtkegel t ≥ |x|.
  - Ausdrücklich: „noch keine Ereignisraumzeit“.
- `experiments/exploration/2026-11-20-operational-origins/VERDICT_CAUSAL.md`:
  - `:27-45` — Herm₂(ℂ) gibt den Lorentzkegel. Reell, komplex, quaternionisch und oktonionisch geben 2+1, 3+1, 5+1 und 9+1;
    Positivität allein wählt ℂ nicht.
  - `:88-90` — „Das Lorentzkegel-Wörterbuch ist in TFPT schon bekannt und darf nicht als neu gelten.“
- `experiments/exploration/2026-11-20-quaternion-frames/VERDICT.md:67` (C) — Lokale Tomographie wählt ℂ. Das ist die
  fehlende Feldwahl zu „nur n = 2“ (z62).
- Ledger Z. 790 (E8.STANDARDFORM.CHART.01, Open), Atlas `finding:corpus-e8-chart-equals-standard-form-m2-pauli-components-20261120`
  (established_scoped).
  - Die 240 E₈-Wurzeln sind Koordinate für Koordinate die Standardform von M₂(ℂ) in Paulikomponenten von a in aΩ.
  - Familienpaare = X, Y, Z; Ankerpaar = Einskomponente × i; J = Multiplikation mit i; aΩ ↦ a*Ω ist die modulare Spiegelung.
  - Das ist dasselbe Objekt wie End(V) in z62.
- `experiments/exploration/2026-11-20-kleinstes-verbundsystem/VERDICT.md:35-40` — Zwei Qubits ergeben Cl(1,3). Drehungen
  wirken lokal, Boosts verschränken. Die Standardform hat die modularen Ladungen 0, 0, +1, −1.

### 1.5 Feldgleichungen aus Holomorphie und Hardy-Kern
- Atlas `finding:h6h21-matrix-kernel-nu-one-massless-spacetime-families-20261124` (z40, conditional) — Der Hardy-Kern ist
  der Wallach-Punkt ν = 1 von det(I − ZW†)^{−ν}. Genau dort gilt die masselose Wellengleichung (Cayley). Das ist die
  Matrixseite von Teil 4 aus z61.
- Atlas `finding:h6h22-wallach-positivity-unitarity-bound-helicity-two-20261124` (z41, conditional) — Der positive
  Helizität-2-Kern Sym⁴(X⁻¹)/det X hat Grad −6. Er ist chiral wie das Graviton in z61.
- `documents/tfpt_horizon_readouts/sections/002-section-turnaround-radius-gravitational-waves-light.tex:238-243` [C] — Ein
  einziger Lorentzkegel für EM, Gravitation und masselose Moden (v_γ = v_g = c). Das passt zu einem Kegel.

### 1.6 c₃ = 1/(8π)
- `docs/research/TFPT-STAND.md:566-567` — c₃ = 1/(|ℤ₂|·2π·χ(S²)) nach Gauß–Bonnet. Das ist dieselbe Zahl wie
  „1/(2 Blätter · 4π)“ in VOGELPERSPEKTIVE:55, also keine neue Lesart.
- Ledger Z. 684 (E8.ONEOBJECT.01) — c₃ = 1/(2πd) bei d = 4 ist ausdrücklich nur eine Umschreibung, keine P1-Herleitung.
  Dieselbe Typisierung gilt für die Himmelssphären-Lesart.
- Ledger Z. 593 (SEAM.THERMAL.KMS.01, SUCCESS) und Z. 293 (HOR.BHTHERMO.01, [E]/[C]) — β_angle = 2π = 1/(4c₃). H6 hat
  β = 4π. Die Normierung ist abzugleichen; vermutlich steckt der Faktor 2 in der Doppelüberlagerung (ungeprüft).

### 1.7 Bindende Firewalls
- Ledger Z. 739 (DIMENSION.UPLIFT.FIREWALL.01) und `documents/tfpt_research_contracts/sections/014-source.tex:183-246`.
  - F1: Kein Satz über 1+1D-Objekte schließt eine 3+1D-Aussage.
  - F2: Jede Schließung braucht eine Dimensionstypisierung.
  - Der kinematische Hopf/Kegel-Uplift (v957) schließt nichts. Das gilt auch für z61/z62.
- Ledger Z. 747 (GRAV.NONCIRCULAR.01) — Nie „KMS/modulares Spektrum ⇒ Lorentzraum ⇒ BW“.
- Ledger Z. 770 (GRAV.WEINBERG.WITTEN.01, Open) — Vor jeder Spin-2-Aussage ist die umgangene Weinberg–Witten-Prämisse
  zu nennen.
- WOIT-Kill-Tests (3), (4), (6), (7): `documents/tfpt_research_contracts/sections/013-source.tex:225-234`. Siehe (C), (G)
  und Abschnitt 9.

## 2. (A) Dublettzustand aus dem Original

- **Offene Forderung D2.** Ledger Z. 790 (E8.STANDARDFORM.CHART.01) und Atlas
  `task:fix-nontracial-state-weight-without-parameter-20261120` (open):
  - Gefordert ist ein nicht spurartiges Zustandsgewicht ohne freien Parameter.
  - Der Spurzustand hat trivialen Modulfluss.
  - Erfolgskriterium: Das Gewicht folgt aus vorhandenen Größen, und das Gitter bleibt erhalten.
- **Freie Eingaben.** `experiments/exploration/2026-11-20-operational-origins/VERDICT_CAUSAL.md:63`:
  - Für ρ = diag(p, 1−p) ist die Phase des Modulflusses (p/(1−p))^{is}.
  - p und die Umrechnung in physische Zeit sind Eingaben.
  - Ein maximal symmetrischer Zustand ρ = I/n gibt keine Zeit.
- **Weitere No-gos zum Spurzustand:**
  - `experiments/exploration/2026-11-20-foreign-origins/VERDICT.md:40` (ρ = I/2 trivial; diag(1/3, 2/3) nicht);
  - `…/2026-11-20-fibonacci-clock/VERDICT.md:27` (braidinvariante treue Zustände sind zentral);
  - `…/2026-11-21-fusion-trace-reflection/VERDICT.md:21` (bedingtes Takesaki-No-go für modulinvariante spurartige
    Teilalgebren).
- Atlas `finding:four-mark-invariant-norm-and-positive-time-nonuniqueness-20261123` (excluded_scoped) — D₄ und Positivität
  legen weder Randnorm noch Zeit fest. Eine Lösung von (A) braucht also mehr als Symmetrie.
- Atlas `finding:process-composition-and-cyclicity-do-not-select-algebra` (excluded_scoped) — Komposition und ein
  zyklischer treuer Zustand wählen die Algebra nicht. „Purifikation bleibt bedingt.“
- **Originalzustand.** Atlas `finding:h6h29-kms-polar-shift-tfd-ladder-weighted-star-20261124` und
  `experiments/exploration/2026-11-24-h6h26-source-dynamics/VERDICT.md:14-16,38-46`:
  - Volle KMS legt (1 − r)r^N fest.
  - Die kanonische Purifikation hat q = e^{−2π}.
  - Offene Prämisse: „Identifikation der beiden Blattalgebren“.
  - Das ist das Material für ein nicht spurartiges Gewicht.
- Atlas `finding:elliptic-current-state-differs-from-hardy-on-every-proper-current-interval` (established_scoped) — Der
  originale Zustand ist nicht der Hardy-Zustand (ω(J₋₁J₁) = 3,487·10⁻⁶ > 0 gegen 0). Hardy-Skalarprodukt allein reicht also
  nicht.
  - Eigene Einordnung: Das SU(2)-invariante Szegő-Produkt auf H⁰(O(1)) gibt den Spurzustand. Erst die KMS-Gewichte machen ihn
    nicht spurartig.
- **Clock-Invarianz.** Ledger Z. 261 (QGEO.SYM.01, [O]), Z. 271 (QGEO.MODULAR.01, [O]) und
  `documents/tfpt_research_contracts/sections/002-source.tex:1-20`:
  - Nach Tomita folgt [ρ, Λ_Σ] = 0 aus ω∘ρ = ω.
  - Eigene Einordnung: Ein μ₄-invarianter Dublettzustand ist diagonal in der Clockbasis {1, z}. Sein Blochvektor zeigt auf
    einen Fixpunkt 0 oder ∞, also auf die „Äonenpole“ aus Z. 737, und der Modulfluss vertauscht mit der Clock.
- Ledger Z. 728 (SEAM.STATE.DERIVATION.01, gemessen) und Z. 736 (SEAM.STATE.WIRING.SELECTOR.01, Freiheitssatz) — Auf der
  v898-Familie schließen RP und Mischung einander aus, und der 2π-KMS-Ort ist ein Punkt. Die Verdrahtung pure-I ist
  „DEPLOYMENT-CHOICE“. Diese Zustandsfreiheit muss ein neuer Ansatz respektieren.
- Ledger Z. 591 (WOIT.BETA2.OS.01, SUCCESS [C]) — Der freie OS-Quotient bei N = 8 ist thermisch, mit dim 16 = 4². Die Clock
  ist ein positiver Transferschritt mit Spektrum {1, √2 − 1}. Ein endlicher Prüfstand für „Zustand auf End(V)“.

## 3. (B) Lesart T aus P1 (O(1) ⊕ O(−1) statt O(1) ⊕ O(1))

- **T ist im Bestand nur gesetzt.** Ledger Z. 604 und Z. 556 sowie `…/011-…-route-label-.tex:322-335`:
  - CELEST und WOIT benutzen die Nahtgerade bzw. PT′ = Tot(O(1) ⊕ O(1)) als deklarierte Konvention.
  - v565: „declared conventions (seam line, …) are cited, not invented“.
  - (B) ist auch dort offen.
- **Euler-Folge.** Atlas `finding:h6h20-euler-sequence-rank-one-channel-zeta-det-selects-flat-20261124` (Z37) und
  `docs/research/TFPT-STAND.md:276-279`:
  - In 0 → O(−1) → ℂ² ⊗ O → O(1) → 0 ist das volle Dublett trivial; die Windung sitzt in der Polarisation.
  - Eigene Einordnung: Hom(O(−1), V′ ⊗ O) = O(1) ⊗ V′ = O(1) ⊕ O(1). Das gesuchte zweite O(1) entsteht, wenn das
    Partnerdublett V′ als triviales Bündel auftritt. Das ist der nächste vorhandene Baustein.
- Atlas `finding:h6h21-boundary-connection-not-weak-connection-polarisation-20261124` (Z40) — Zwei Schnitte von O(1) tragen
  keine lokale SU(2). Erst das O² der Euler-Folge trägt sie.
- **Paarzustände und Matrixkern.** Atlas `finding:h6h23-orientation-reversing-sheet-pair-states-higgs-twist-flavour-20261124`
  (Z42), mit `finding:h6h21-matrix-kernel-nu-one-…` (Z40):
  - Zwei Bosonenpaare geben ⟨W|Z⟩ = 1/det(I − W†Z), auf jedem Rang-1-Schnitt Hardy.
  - Die Lagen Z ∈ M₂(ℂ) = Hom(V, V′) erscheinen hier schon als Parameter von Paarzuständen (Dublett ⊗ Partner).
  - Nach TFPT-STAND:317 braucht der Matrixkern „eine weitere Regel“.
- **Partner.** Zweite Blatthälfte orientierungsumkehrend (8 − 8 = 0): derselbe Atlas-Befund Z42. Das stützt V′ = V̄.
  Kanonische Purifikation: `…/2026-11-24-h6h26-source-dynamics/VERDICT.md:16`, siehe (A).
- Atlas `finding:h6h18-v4-split-is-p2-colour-equals-lambda2-f-bundles-o1-6o-om1-20261124` (Z34):
  - r_* r^* O(1) = O(1) ⊕ 6O ⊕ O(−1) mit E₋ ≅ E₊^∨ und einem O(−1) ohne Nullmode.
  - Relative Selbstdualität erzwingt L = r^*O(1).
  - Das ist der Ursprung der Lücke in z61.
- **Konkurrierende Normalen-Lesarten** (siehe W6):
  - Ledger Z. 318 (TOPO.E8.01, [E]/[C]): Naht als du-Val-(−2)-Kurve.
  - `documents/origin_theory/sections/003-section-the-unique-attractor-a-gapped-transport-selects-.tex:223-236`: „seam
    normal bundle = ADE orbifold point ℂ²/Γ“, Γ = 2I, die eine offene Pfeilrichtung von SEAM.EQUIV.01.
- **Bedingung für die Partnerdoppelung.** `experiments/exploration/2026-11-20-kleinstes-verbundsystem/VERDICT.md:35-36` —
  Auf Dublett ⊗ Partner (zwei Qubits) sind Boosts verschränkend. Ohne Partner gibt es keine Boosts, also braucht (B) die
  Purifikation auch für die Lorentzgruppe als Wirkung.

## 4. (C) Nicht-selbstdualer Teil der Gravitation

- **Scope-Zaun:** `documents/tfpt_research_contracts/sections/010-source.tex:255-260` und `014-source.tex:58-66`.
  - „Self-dual is not full 4D physics“; keine Anspruch auf beide Helizitäten.
  - Kill-Test (4) von WOIT (`013-source.tex:227`): Der Penrose-Transform erreicht nur freie bzw. selbstduale Zustände.
  - Das googly-Problem ist also im Bestand als harte Wand registriert.
- Ledger Z. 740 (SEAM.BULK4D.RECON.01, Open) und `014-source.tex:256-280` — Forderung B4 heißt „both helicities (the
  anti-self-dual completion — the known hard wall of the twistor route)“. Kill-Test (2): Der Zuflussbulk ist nur
  topologisch.
- Ledger Z. 762 (GRAV.SPIN2.EMERGENCE.01, Open) — verlangt zwei Helizitäten, universelle Kopplung, Diffeo-Ward und RP.
  Z. 770 (Weinberg–Witten) gilt zusätzlich.
- Ledger Z. 786 (GRAV.FREE.CURVATURE.COVARIANT.01, [E] freies Ziel) — Ein positiv-energetisches kovariantes freies
  Krümmungsfeld mit Helizitäten ±2 existiert als Ziel. Atlas `pattern:free-curvature-versus-parent`.
- **Positivität.** Ledger Z. 80, 93 u. a. (Spin-2-Zeilen mit Positivitätsvermerk v1037,
  `experiments/theory-contracts/spin2-positivity-20261124/README.md`) — Der ganze Formfaktor ist als positiver Propagator
  ausgeschlossen [X]. Positive Spin-2-Kerne: Z41 (Wallach) und Z47 (Atlas
  `finding:h6h25-ward-structure-index-protection-single-source-spin2-channel-20261124`).
- **Kandidat „zweites Blatt = duale Chiralität“:**
  - Zweite Blatthälfte orientierungsumkehrend (Z42): Eine Orientierungsumkehr tauscht SD und ASD.
  - PT ↔ PT*-Dualität (v565): der natürliche Startpunkt für Ambitwistoren.
  - Zu Ambitwistoren oder zur Ergänzung dritter Ordnung gibt es keinen Bestand.
- **Λ² und Plebański:**
  - Ledger Z. 755 (DIMENSION.SELECTOR.4D.01, [C]): In Lorentzsignatur gilt *² = −1, die reelle SD/ASD-Teilung gibt es nur
    euklidisch. Λ²₊ ist lorentzisch komplex.
  - Atlas `finding:h6h15-plebanski-isotropy-ward-identity-fixes-spin2-not-g-20261124` (Z29): S₄ erzwingt δ_ij; Urbantke
    lorentzisch nur mit Re(a b̄) = 0.
  - Atlas `finding:h6h30-plebanski-area-rule-fixes-quartic-number-t-one-20261124`: Flächenregel t = 1.
  - Hinweis (bekannte Physik): Im chiralen Plebański-Funktional ist das spurfreie Ψ_ij auf Λ²₊ der Lagrange-Multiplikator
    und on-shell der selbstduale Weyltensor. Das verbindet (C) mit (E).
- Atlas `finding:fable-cartan-reflected-lattice-two-light-cone-coordinates-20261112` (established_scoped):
  - Die 4+4-Blockspiegelung des E₈ gibt 24 links, 24 rechts und 192 gemischte Wurzeln.
  - Das passt zu z60 (D₄ ⊕ D₄): ein möglicher Ort für zwei Chiralitäten im E₈-Gitter, aber nichtchiral gewählt.

## 5. (D) Wechselwirkung, Massen, Rolle von E₈

- **E₈ für SDYM auf Twistorraum zulässig.** `documents/tfpt_research_contracts/sections/010-source.tex:187-192`:
  - Costellos Ein-Schleifen-Axionkonsistenz (arXiv:2111.08879) lässt E₈ für SDYM auf dem Twistorraum zu (Okubo).
  - CELEST.WP1.01 (Ledger Z. 555, [E]): Das μ₄-Glue ist eine flache ℤ₄-Monodromie, die Nullmoden sind der Träger D₅ ⊕ A₃
    (60).
  - Unter T ist das die natürliche Rolle von E₈: Eichalgebra der holomorphen Twistortheorie, also Ward.
  - Ziel von WOIT.OS.TWISTOR.01 (Z. 585): wechselwirkendes A_hol = SDYM(E₈) + BCOV. Kill-Tests (4) und (5) betreffen genau
    die Wechselwirkung.
- **E₈ = ganzzahlige Struktur der Dublett-Standardform.**
  - Ledger Z. 790 und Atlas-Problem `problem:standard-form-chart-of-gaussian-e8-20261120`; Aufgabe
    `task:derive-stabilizer-lattice-from-standard-form-20261120` (open).
  - `experiments/theory-contracts/compiler-cone-object-audit/README.md:9-15`: maximale Ordnung in M₂(ℚ(i)).
  - Das ist ein zweiter Ort von E₈ direkt auf End(V) aus z62; die Ganzzahligkeit ist offen.
- **E₈-Naht starr.** z60 (`VERDICT.md:94-99`), Atlas
  `finding:h6h32-spectral-dimension-no-noncompact-four-dimensional-response-in-originals-20261124`.
  - Lokales E₈-Netz aus dem Hardy-Kozyklus: Z32, TFPT-STAND:253-258.
  - Ströme: Atlas `finding:original-e8-current-response-and-full-quarter-lifts`,
    `finding:e8-current-dictionary-all-separated-responses`.
- **Kopplungen:**
  - Atlas `finding:h6h14-rg-test-induced-equal-couplings-at-lattice-scale-fail-20261123` (Z25): Induzierte gleiche
    Kopplungen fallen gegen die Daten.
  - Ledger Z. 358/359 (PS.THRESHOLD.02 / PS.PROTON.02, [X] korrigiert mit v1036).
  - Wer Kopplungen aus Twistor-SDYM liest, muss diese Datenschranken bestehen.
- **Massen** (TFPT-STAND §3: Seed, L, K, [E] in den Annahmen) — Kein Bestand verbindet sie mit der Twistorlesart.

## 6. (E) Familiengewichte γ und TM1-Brechung

- **(λ, −λ, 0) bei J = 0 ist bekannt.** Atlas `finding:pauli-quadrupole-discriminant-square-torus-20261123`: tr Q³ = 0 gibt
  die Eigenwerte (a, 0, −a). Die prüfbare Idee aus VOGELPERSPEKTIVE:92-93 ist dort algebraisch schon enthalten; die
  Familienzuordnung fehlt.
- **Eingefrorene Werte:**
  - Ledger Z. 17 (FLAV.TH12.01, Numerical): sin²θ₁₂ = 1/3 − φ₀/2 = 0,30675.
  - Z. 361 (FLAV.TH13.01, [E]/[C]/[O]): sin²θ₁₃ = φ₀e^{−5/6} = 0,0231 mit **γ = 5/6 = tr_E Y²**. Namenskonflikt mit dem
    Clockgewicht γ aus H6h30, siehe W9.
  - Z. 476 (FLAV.PMNS.CLOSE.01, [E]/[C]/[O]): θ₂₃ = 45°, δ = 4π/3.
  - Z. 440 (THETA13.BUDGET.01): JUNO entscheidet.
  - Jede Gewichtsmechanik muss diese Zahlen tragen.
- **Flavour bei τ = ω.** Ledger Z. 475 (SEAM.EQUIV.02, [E]/[C]/[O]) und Z. 474 (E8.CMDUAL.01, [E]/[C]):
  - Der Familien- und Flavoursektor sitzt bei τ = ω (j = 0, χ_E8 = 0), dual zur Naht bei τ = i.
  - δ_PMNS = 4π/3 ist die CM-Phase bei τ = ω.
  - In Petrov-Sprache ist das der Zustand mit I = 0 (Eigenwerte λ(1, ω, ω²)), nicht J = 0. Siehe W9.
- **TM1-Hindernisse:**
  - Atlas `finding:h6h29-tm1-fifteen-integrality-obstruction-tm2-solar-sum-rule-20261124`: TM1 bei 15° mit Z₃,
    Majorana-Klein und μ-τ ist in jeder endlichen Gruppe ausgeschlossen; TM2 verletzt die Solarsummenregel.
  - Atlas `finding:h6h30-clock-operators-carry-tm1-column-frozen-values-break-tm1-20261124`: F_C und F_J tragen die
    TM1-Spalte (2/3, 1/6, 1/6). Der Winkel hängt am freien γ (γ = 1/8 gibt 14,748°).
  - Atlas `finding:h6-common-parent-clock-full-algebra-20261124` (established_scoped): F_C, F_C⁻ und F_J erzeugen M₃(ℂ); das
    eingefrorene |U_e1|² = 0,6772 ≠ 2/3.
  - Atlas `finding:h6h26-delta-6n2-residuals-only-tm2-tm1-fifteen-cp-prediction-20261124`: TM1 bei 15° mit μ-τ-Reflexion
    gibt δ = ±π/2.
- **Familienanteil der Zeit.** `experiments/theory-contracts/source-time-family-generator-20261014/README.md:9-12` — Das
  Fold-Zeitgesetz δ_k = (k+α)/8 − ½ hat einen nichtskalaren Familienanteil δ_F: {±x, ±¼, ±¼} mit x = 3/16 oder 1/16, auf
  den Marken 1+1+2. Ein vorhandener Kandidat für die Familiengewichte aus dem Original.
- **Weitere Bausteine:**
  - Atlas `finding:return-cone-needs-clifford-relation-pairings-give-lorentz-frame-positive-spin-form-has-no-spacelike-pairs-20261119`:
    Die drei Paarungen der vier Marken sind antikommutierende spurfreie Möbiusmatrizen und spannen sl(2,ℂ) auf. Damit sind
    die Familien die drei V₄-Involutionen.
  - Plebański-Hinweis aus (C): Das Weyl-Ψ auf Λ²₊ ist der Plebański-Multiplikator. Gewichte (λ, −λ, 0) brechen die
    S₄-Isotropie δ_ij aus Z29 auf D₄ (zwei invariante Formen, Z32).

## 7. (F) G und Skalen

- **Eine Einheit.** Ledger Z. 72 (VGEO.SHARPEN.01, [E]/[O]):
  - v_geo ist die einzige dimensionsbehaftete Eingabe (No-Unit-Satz v153).
  - Die Gravitation fügt keine Skala hinzu.
  - Bestätigt (F) „G bleibt eine Einheit“.
- **Absolute Skala.** Ledger Z. 367 (SCALE.OVERDET.01, [E]+[N]+[A]):
  - Zwei Routen zu M̄ stimmen auf 0,11 % überein.
  - „No derivable absolute scale is a theorem“: Die Naht ist konform mit c = 8.
- Ledger Z. 75 (GRAV.NONLINEAR.01, [E]/[C]/[O]): 8π = 1/c₃. Ledger Z. 144 und Z. 670 (ANCHOR.VGEO.01, VGEO.TORSOR.01):
  geschlossene Metrologie als ℝ₊-Torsor.
- **G im H6-Strang:**
  - Atlas `finding:h6h15-plebanski-isotropy-ward-identity-fixes-spin2-not-g-20261124`: Ward legt Spin 2 fest, nicht G.
  - Atlas `finding:h6h29-fourpoint-generating-function-depth-parameters-g-scaling-20261124`: G ist ohne Normregel frei.
  - H6h30: G = v_geo⁻² (TFPT-STAND:381).
- **Gitterabstand.** Atlas `finding:h6h9-lattice-spacing-not-fixed-edge-term-dominates-20261123` und
  `finding:h6h3-lattice-not-at-planck-scale-under-jacobson-20261123` (excluded_scoped): a/ℓ_P ist nicht festgelegt.
- Atlas `pattern:gravity-coefficient-vs-source` — Gravitationskoeffizient und lokale Quelle sind getrennte Lücken.

## 8. (G) 16 Weylfermionen je Familie aus einem Index

- **Darstellung, kein Index.** Ledger Z. 408 (CAR.SM.ANOM.01, Identity [E]) und Z. 331 (QFT4D.MATTER.01, [E]/[O]):
  - 16 = Λ^gerade ℂ⁵ = 1 + 10 + 5 ist eine anomaliefreie Generation.
  - N_fam = 3 legt b₁ = 41/10 fest.
  - Die 16 je Familie sind etabliert, aber als Darstellung.
- **Naiver Index.** Atlas `finding:h6-common-parent-family-index-obstruction-20261124` (excluded_scoped):
  - Im geraden U(5)-Fockbündel gilt c₁ = 8c₁(C); der naive Index ist 8 mod 16.
  - Ausgeschlossen ist nur diese naive Zuordnung.
- **Ältere Indexrechnung.** Atlas `finding:normal-orbit-spin-family-index-and-signed-charge-balance` (conditional, Oktober):
  - Mit dem Lift vom 22.9. und Σm = −2 ist der Produktindex −16d, gewichtet I_B = I_L = −4d.
  - Je 16 Moden in K±.
  - Andere Voraussetzungen als im Paket vom 24.11.; ein möglicher Ausweg, zu prüfen.
  - Quellen: `experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/CLAUDE18-FERMION-INDEX-REVIEW.txt`,
    `experiments/theory-contracts/source-joint-spin-family-lift-20260922/PROOF.md:36-123`.
- Atlas `finding:h6h28-gaussian-from-free-seam-z16-balance-mirror-gap-dirac-flavour-20261124` (Z53):
  - 16 Weyl mit Spin-ℤ₄, 16 ≡ 0 mod 16; ohne ν_R wären es 15.
  - Die mod-16-Bedingung ist notwendig, nicht hinreichend.
- **Naht-16 gegen Generations-16.** Atlas `finding:h6h13-seam-and-generation-sixteen-cannot-unify-direct-sum-only-20261123`
  (Z23): Sie sind nicht vereinbar, nur als direkte Summe.
- **Familienzahl aus K(D).** Atlas `finding:h6h8-mark-index-gives-three-families-and-review-corrections-20261123` (Z10,
  established_scoped):
  - ker D₊ = 3 auf K^{1/2} ⊗ O(D).
  - Λ^gerade ℂ⁵ = 16 kommt dort als Fockraum hinzu, nicht aus dem Index.
- Atlas `finding:originlab-unit-winding-fixed-s4-analytic-family-obstruction-20261114` (established_scoped) — Einheitswindung
  auf der festen analytischen Scheibe lässt keinen Familienindex 3 zu (n₀ + 3n₃ = 1). Die 3 muss aus K(D) kommen, nicht aus
  der Windung.
- `experiments/exploration/2026-11-23-fockspur-gso/VERDICT.md` (F3, Kandidat Stufe 2):
  - Kette μ₄ → 16 → 5 → 3+2.
  - Die ℤ₁₆-Anomalie mit μ₄ gibt Vielfache von 16; eine konsistente ℤ₄-Ladung gibt es nur für ungerades N, kleinstes N = 5.
  - Das ist die bisher beste Herkunft der 16, allerdings aus Anomalie und Fockraum.
- **4D-Chiralität:**
  - Ledger Z. 746 (CHIRAL4D.NOMIRROR.01, Open): 4D-Chiralität verlangt einen GW/Overlap-Operator.
  - Z. 766 (CHIRAL4D.MIRROR.DET16.01, Candidate [C]).
  - WOIT-Kill-Test (3): spiegelartige Generation (`013-source.tex:225`).
  - Atlas `finding:h6h31-chirality-from-index-without-mirror-gap-20261124`: nur auf der Naht.

## 9. Widersprüche und Spannungen (präzise)

Kein etabliertes oder ausgeschlossenes Ergebnis widerspricht (a)–(e) direkt. Die folgenden Punkte sind Bedingungen, die
eine Formulierung der Architektur erfüllen muss. W1, W4, W6, W8 und W9 wiegen am schwersten.

- **W1, zu (b) „Familien = Λ²₊ = H⁰(O(2))“: μ₄-Charaktere.**
  - Die Familien haben unter z ↦ iz die Charaktere {i, −1, −i}: tfpt_1 `005-source.tex:139-141` [E], Z10, CELEST.WP2.01.
  - Holomorphe Vektorfelder H⁰(T) = sl(2) haben {−i, 1, i}, mit invarianter Drehachse.
  - Die Identifikation K(D) ≅ T läuft über q = dz²/P, und q hat den Charakter −1 (eigene sympy-Kurzprobe).
  - Atlas `finding:four-mark-bilinear-quarter-character-20261123` bestätigt das: „Der gerade Dreierraum und die drei
    logarithmischen Formen haben verschiedene R Charaktere.“
  - Folge: Familien ≅ Λ²₊ ⊗ χ_q gilt nur mit Twist. Ohne Twist hätte eine Familie Charakter 1, im Widerspruch zum chiralen
    Familienindex von H6 (ker D₋ = 0, Charaktere {−1, i, −i}).
  - Eigene Einordnung: Im CELEST-Modell sind die Familien die drei Sphären H₂(ALE). Ihre Schnittform ist −Cartan (Z. 318
    TOPO.E8.01), sie liegen also im ASD-Teil. Die Clock wirkt auf die hyperkählerschen SD-Formen mit {1, i, −i}.
  - Außerdem tragen die Familien die A₃-Cartanform (tfpt_1 `005-source.tex:132-138`, det = 4). Ihre Abbildung auf die
    Λ²₊-Form (Killing bzw. Plebański δ_ij) ist zu prüfen; D₄ lässt zwei invariante Formen zu (Z32).
- **W2, zu (a) Träger = Spin 2.**
  - Linearisierung: Gewöhnliche Schnitte von O(4) haben Spur +1 unter r, der Spin-2-Quadrupol hat −1
    (quaternion-frames `VERDICT.md:51`). Nötig ist der D₄-Twist (r, s) ↦ (−1, +1); „eine native Gleichsetzung muss diesen
    Twist physisch liefern“.
  - Eigene Kurzprobe: Der Twist ist genau der Clock-Eigenwert −1 des Markenzustands (|2⟩ − |−2⟩)/√2, also wieder χ_q.
  - 3+2-Schnitt: Die Lesart aus fockspur-gso F4 (Parität des Spingewichts unter π) ist der Halbdrehungsschnitt. Z34/Z36
    zeigen, dass er nicht der äquivariante V₄-Schnitt des Originals ist; 4 von 5 Hyperladungen weichen ab (Atlas
    `finding:h6h18-…`). Der richtige schwache Teil ist {1 − z⁴, z²}.
  - Beobachtung (eigene Einordnung): Das Markenpolynom 1 − z⁴ ist selbst ein Basisvektor des schwachen V₄-Dubletts.
- **W3, zu (a)/(d): J = 0 reicht nicht.**
  - Ledger Z. 570 (SEAM.CLOCK.RIGIDITY.01, [E]) gibt ein Gegenbeispiel: {±1, ±(3+2√2)} ist harmonisch (j = 1728), aber das
    Deck ist nicht zentral, und es gibt keine markentreue Clock.
  - Neben J = 0 bleibt also das Ausrichtungsbit (v506, v507, v510).
  - WOIT-Kill-Test (7) (`013-source.tex:233`): Die Marken müssen inzidenzverträglich über die Raumzeit fortsetzbar sein.
  - Eigene Einordnung: J ist kubisch in Ψ und bleibt unter linearer Superposition masseloser Felder nicht erhalten. Unter T
    ist J = 0 deshalb eine Bedingung an das Feld, nicht automatisch.
  - Ledger Z. 778 (GRAV.PARENT.GAP.NO_GO.01, [E]): „the finite five-component symmetric-traceless target is only
    representation typing and is not a derived stress operator“. Das stützt „Marken gehören zu einem Zustand im Träger“.
- **W4, zu (c) und z58: Modulfluss auf dem Dublett ist eine Drehung.**
  - z62 (`ausgabe_z62.txt:10`) sagt selbst: Der Modulfluss ist eine Drehung im Ruhesystem.
  - VERDICT_CAUSAL:47-61 und Atlas `finding:h6-common-parent-modular-norm-obstruction-20261124` (excluded_scoped): Endliche
    oder beschränkte Modulflüsse können weder Boosts noch Translationen oder Dilatationen sein.
  - Die „Zeit = Modulfluss“ aus z58 ist auf End(V) also eine Drehung um die Zustandsachse. Raumzeit-Translationen müssen aus
    T kommen (z61).
  - Ledger Z. 293 (HOR.BHTHERMO.01) und Z. 593 (SEAM.THERMAL.KMS.01) lesen denselben Modulfluss als Horizont-Boost (BW).
  - GRAV.NONCIRCULAR.01 (Z. 747) verbietet „modular ⇒ Lorentz“.
  - Die beiden Rollen des Modulflusses müssen getrennt typisiert werden. Das ist eine Spannung innerhalb der Architektur.
- **W5, zu (c): Vorbestand und Feldwahl.**
  - Der Kegel ist nicht neu (VERDICT_CAUSAL:90; compiler-cone-object-audit, 9.9.).
  - „Nur n = 2“ gilt innerhalb von M_n(ℂ). Positivität allein wählt ℂ nicht (VERDICT_CAUSAL:45); ℂ kommt aus lokaler
    Tomographie (quaternion-frames C) oder aus der Holomorphie der Naht.
  - „A Lorentz-shaped state cone is not sufficient to identify physical Lorentz boosts“ (VERDICT_CAUSAL:43): x ↦ gxg† mit
    nicht unitärem g ist ein Filter, keine Zeitentwicklung.
- **W6, zu (d) T: Lage der Nahtgeraden.**
  - Im Bestand (v565, [E] mit deklarierter Konvention) hat die Nahtgerade im Minkowski-Zweig die Norm a(|z|² − 1). Die
    Nullmenge ist der Nahtkreis.
  - Eigene Einordnung (Standard-Twistorgeometrie): Unter ω = i x π ist Σ|_{L_x} = −2π† (Im x) π. Eine Gerade mit
    indefiniter Norm gehört zu einem komplexen Punkt mit raumartigem Im x, nicht zu einem reellen Punkt. Dessen Gerade läge
    ganz in PN; das ist die Himmelssphäre aus z61/z62.
  - Dichtematrizen im Zukunftskegel als Im x gehörten zu PT₊ (euklidischer Zweig, Hardy).
  - Lesart T muss festlegen, welchen reellen Typ die Nahtgerade hat. Beide Bilder legen die vier Marken (|z| = 1) auf
    Nullrichtungen.
  - Zusätzlich konkurrieren TOPO.E8.01 [C] (du-Val-(−2)-Kurve) und origin_theory `003-…:223-236` (Normalen = ℂ²/2I) als
    Lesarten der Naht. Beide sind mit O(1) ⊕ O(1) nicht identisch.
- **W7, zu (d) Naht = Himmelssphäre und E₈: Spinstruktur.**
  - `experiments/theory-contracts/seam-double-spin-structure-20261014/README.md:7-31` und `PROOF.md`.
  - Ein in die Himmelssphäre eingebetteter Nahtkreis trägt NS (v506, Theorem B). Dann liegt die Clock nicht in der
    Glueklasse A₃ ⊕ D₅, sondern in A₁ ⊕ D₆ bzw. A₁ ⊕ A₇.
  - Die Glueklasse braucht R, also die Doppelüberlagerung.
  - Die 16 Plätze müssen auf der Doppelüberlagerung liegen. Das ist mit H6 verträglich (5+3 im periodischen Sektor, Z38).
    „Naht = Himmelssphäre“ darf die Fermionen aber nicht auf einen eingebetteten Kreis legen.
- **W8, zu (b) und zum Wörterbuch „2 = Spinor bzw. Dublett“: Coleman–Mandula und Kill-Test (6).**
  - `experiments/theory-contracts/universalraum-fable-kernfragen-20260910/fable-runde3/BERICHT.md:199-200,244-245`:
    Coleman–Mandula trifft die kontinuierliche Gleichsetzung innerer Generatoren mit Raumzeit-Generatoren.
  - WOIT-Kill-Test (6) (`013-source.tex:231`): Die innere SU(2) bleibt ein Raumzeitfaktor.
  - Ledger Z. 634 (WOIT.SU2.KIN.01, [E]): Die innere su(2) trennt sich kinematisch.
  - kleinstes-verbundsystem `VERDICT.md:49`: „Dieselben drei nichtlokalen Erzeuger werden als Farbplätze und als
    Raumrichtungen gelesen.“
  - „Familien = Λ²₊“ und „Dublett = Spinor“ sind daher als Zustands- bzw. Labelräume zu lesen, nicht als Generatoren.
  - Das schwache Dublett ist nicht H⁰(O(1)): Z34 gibt {1 − z⁴, z²} ⊂ H⁰(O(4)), und Z40 sagt, zwei Schnitte von O(1) tragen
    keine lokale SU(2).
- **W9, zu (E).**
  - TM1 mit μ-τ-Reflexion gibt δ = ±π/2 (Atlas `finding:h6h26-delta-6n2-…`). Eingefroren ist δ_PMNS = 4π/3 (Ledger Z. 476,
    [C] über SEAM.EQUIV.02). TFPT-STAND:495-496 nennt δ = 240° einen zugewiesenen Texturwert.
  - Der Bestand verortet Flavour bei τ = ω (I = 0, Ledger Z. 474/475); die Idee aus (E) arbeitet bei J = 0.
  - Die eingefrorenen sin²θ₁₂ = 0,30675 und sin²θ₁₃ = 0,0231 verletzen exaktes TM1: |U_e1|² = 0,6772, Atlas
    `finding:h6-common-parent-clock-full-algebra-20261124`. Exaktes TM1 würde sin²θ₁₂ ≈ 0,318 verlangen, eigene Nachrechnung.
    Jede Gewichtsmechanik muss TM1 also in genau diesem Maß brechen.
  - Namenskonflikt: γ = 5/6 = tr_E Y² (Ledger Z. 361) gegen das freie Clockgewicht γ (= 1/8) aus H6h30.
- **W10, zu (F).** „Dimensionslose Verhältnisse brauchen eine skalentragende Antwort“ (VOGELPERSPEKTIVE:95) ist einzugrenzen.
  Der Ledger führt bereits dimensionslose Gravitations- und Kosmologieverhältnisse: ρ_Λ/M̄⁴ = (3/4π²)e^{−2α⁻¹} und
  M_scal/M̄ = c₃^{7/2} (Z. 367, Z. 75). Die Aussage gilt für die Rekonstruktion von G aus Antworten, nicht für die
  Ledger-Readouts.
- **W11, zu (C).** Der Scope-Zaun (CELEST/WOIT: nur selbstdual; B4 als harte Wand) und GRAV.SPIN2.EMERGENCE.01 (zwei
  Helizitäten) sind die bindenden Abnahmekriterien. z61 selbst respektiert sie.
- **Zu (e) z58–z62, keine Konflikte gefunden:**
  - z58: Chiralität aus dem Index bleibt auf der Naht und ist mit Firewall F1 und CHIRAL4D.NOMIRROR.01 verträglich, wenn sie
    typisiert wird.
  - z59/z60 sind verträglich mit SEAM.BULK4D.RECON.01 („no existing corpus result closes any part“).
  - z61, Teil 4, deckt sich mit WOIT-Kill-Test (4).
  - Der Atlas-Ausschluss der rohen Tetraedermarken (`finding:tetrahedral-quarter-face-phase-…`, j = 0) trifft die H6-Marken
    (j = 1728) nicht.

## 10. Wiederverwenden statt neu rechnen

1. **τ = i ⇔ J = 0 und Quadrupol:** v493 (CELEST.WP2.01) und der Atlas-Befund zum Pauli-Quadrupol. Zur Auswahl die
   Befunde zum Eigenwertmodus und zur radialen Verteilung (sextischer Term).
2. **Standardform von M₂(ℂ) = E₈-Karte:** E8.STANDARDFORM.CHART.01 mit dem Atlas-Problem und drei offenen Aufgaben. (A) ist
   dort D2, die Ganzzahligkeit D1, der Faktor 2π von c₃ D4.
3. **Twistorrahmen mit Kill-Tests:** CELEST.SEAM.01, WOIT.OS.TWISTOR.01 (α bis β₃ ausgeführt) und PT′ = Tot(O(1) ⊕ O(1))
   aus v514. Der Scope-Zaun ist zu übernehmen. Die Lage der Nahtgeraden ist mit v565 abzugleichen (W6).
4. **Paarzustände und Matrixkern:** 1/det(I − W†Z) und Hardy bei ν = 1 (Z40, Z42). Das ist der vorhandene Weg zu
   Hom(V, V′) und zur masselosen Gleichung, Material für (B).
5. **Euler-Folge, Partner, TFD:** Z37 (triviales Dublett ℂ² ⊗ O), Z42 (zweites Blatt orientierungsumkehrend) und H6h26
   (kanonische Purifikation; Blattalgebren identifizieren bleibt Prämisse). Das sind die Bausteine für (B) und (C).
6. **Träger = Sym²₀(Familien):** die exakte Folge Sym²H⁰(K(D)) → H⁰(O(D)) über q = dz²/P (Atlas, 23.11.). Sie verbindet (a)
   und (b) über den Charakter von q (W1, W2).
7. **E₈ für SDYM auf dem Twistorraum:** zulässig nach Costello (Vertrag `010-source.tex:187-192`). Das ist der vorhandene
   Anschluss für (D).
8. **Familiengewichte:** δ_F aus dem Fold-Zeitgesetz (14.10.), die TM1-Spalte der Clocks (H6h30) und die Hindernisse aus
   H6h29 und common parent.
9. **Index −16d:** aus dem Normalorbit-Lift (Oktober). Zu prüfen, ob er das 8-mod-16-Hindernis umgeht.

## 11. Grenzen dieser Prüfung

- Archivierte Papierversionen, die 152 nicht eingearbeiteten Atlas-Versuche und die meisten der 615 Theorieverträge sind
  nicht einzeln gelesen. Die Volltextsuche deckt nur die genannten Stichworte ab.
- Die Zuordnungen in W1, W2 (Charaktere), W6 (Twistorgeometrie) und W4 sind eigene Einordnungen. Sie stehen nicht im
  Bestand und sind dort nicht gerechnet.
- Atlas-Befunde sind mit ihrem Status zitiert. Ihre KI-Zusammenfassungen sind nicht geprüft.
