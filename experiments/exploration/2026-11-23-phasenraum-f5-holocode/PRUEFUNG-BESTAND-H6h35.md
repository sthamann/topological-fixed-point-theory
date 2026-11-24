# Prüfung des Bestands gegen die dritte externe Analyse (Vorbereitung H6h35)

Stand 24. November 2026. Nur lesende Prüfung: keine Rechnung im Repo, kein Commit, keine Atlas-Schreibung, kein Refresh.
Sie ergänzt `PRUEFUNG-BESTAND-H6h34.md`; was dort geprüft ist, wird nur zitiert.

**Gegenstand.** Die dritte externe Analyse vom 24.11. liegt nicht im Repo. Geprüft sind ihre sechs Behauptungen in der
Fassung des Auftrags:
1. Die kontinuierliche Zeittranslation ist P₀ = L₀ + (L₋₁ + L₁)/2 in der vorhandenen sl(2): parabolisch, Spektrum [0, ∞).
   Die TFPT-Clock L₀ ist elliptisch.
2. Die bilokale Teilchen-Loch-Auslese D_m(γ) = D₁(γ) sin(mγ)/sin γ gibt die transversale S³-Antwort nativ.
3. Die m Teilchen-Loch-Paare der Strommode m geben mit Tomita-Verdopplung die Vielfachheit m², also die S³-Vielfachheiten
   bzw. die K-Typen der minimalen SU(2,2)-Darstellung (masseloser Skalar). Vier vertauschende Translationen P_μ mit
   P₀ ≥ 0 und P² = 0 entstehen in Standardform.
4. det(I − ZW†)^{−1} ist unter den positiven holomorphen bi-U(2)-kovarianten Kernen, die Hardy- und Stromauslese treffen,
   eindeutig.
5. E₈ lässt sich nicht als freie 4D-Felder fortsetzen.
6. Drei Familien kommen aus einem relativen Toeplitz-Index; D₄ erzwingt zwei gleiche Familiengewichte.

Flavour, Massen und Teilcheninhalt sind gelöste Bausteine (TFPT-STAND §3) und werden nicht als Lücke behandelt.

**Umfang.**
- Ledger `verification/status_ledger.csv` (790 Zeilen), zitiert als „Ledger Z. n (ID)“; n ist die CSV-Zeile.
- Papiere unter `documents/` per Volltextsuche. SU(2,2), Leiterdarstellung, Wallach-Menge und minimale Darstellung kommen
  dort nicht vor; halbseitige Inklusionen und Borchers schon.
- Atlas: Revision `7114eb296d62f303591f7068` (nur gelesen), MCP-Suche und Dossiers.
- Experimente: dieser Ordner (z32, z40–z44, z58–z71), `experiments/theory-contracts/original-paper-source-reconciliation-20261102/`,
  `research/originlab/FAMILY-INDEX-20261114.md`, `docs/research/fundamental-plan/TFPT-REGIONALE-KOMPOSITION-20261112.txt` und
  der Artikel `documents/articles--2026-08-29--holomorphic_kms_extension_en`.
- Eine eigene numpy-Kurzprobe im Scratchpad (nicht abgelegt): Bilokal-Identität, Paraboliken der Marken, Gegenbeispielkern,
  Linsenvielfachheiten.

**Lesehinweis.** Statusangaben stammen aus der Quelle. Atlas-Links und KI-Zusammenfassungen sind Navigation, kein Beweis.
„Eigene Einordnung“ heißt: Folgerung dieser Prüfung, nicht im Bestand. „Bekannte Mathematik/Physik“ heißt: Standardliteratur,
hier nicht nachgerechnet.

## 0. Kurzbilanz

1. **Die Mathematik hält, ist aber fast vollständig Vorbestand.**
   - (1) steht als [E] in v438 (Ledger Z. 496) samt Lean-Kern (Z. 503) und in der Zeitentscheidung vom 5.11.
   - (3) und (4) sind z40–z44: Leiterdarstellung von U(2,2), Wallach-Punkt ν = 1, „der Kern braucht eine weitere Regel“,
     sp(4, ℝ) gegen su(2,2).
   - (5) ist für punktweise Lokalisierung schon ausgeschlossen (Atlas, 12.11.).
   - (6) ist Z10, z58, z32 und z67.
2. **Neu sind drei Punkte:**
   - Der bilokale Relativwinkel γ ist eine native Variable mit exakt der S³-Zonalform. Das korrigiert die Begründung von z71,
     nicht dessen Ergebnis.
   - Die Lesart „m Paare ⊗ Tomita-Partner = K-Typ V_j ⊗ V_j“ für die Paaralgebra aus z44.
   - Die Eindeutigkeit des Matrixkerns über die Paarauslese statt über die Wellengleichung.
3. **Last tragend sind überall dieselben drei Voraussetzungen:**
   - T′: ein Partner mit eigener U(2), im Bestand die einseitige Blattaktion (z44) bzw. V̄* (z64).
   - Bosonische Paarstatistik (z42, z44).
   - Der Vakuumzustand. Der Originalzustand ist ein Gibbszustand der elliptischen Clock und nicht P₀-invariant.
4. **Firewalls:**
   - DIMENSION.UPLIFT.FIREWALL.01 blockiert (1)–(3) als 3+1-Aussagen; (4) ist als Satz auf der Matrixkugel zulässig.
   - GRAV.NONCIRCULAR.01 greift bei (1) und (3), sobald P₀ oder SU(2,2) aus Modulfluss/KMS gewonnen oder zu BW
     weitergeführt werden.
   - GRAV.WEINBERG.WITTEN.01 trifft (1)–(4) nicht (Helizität 0), bindet aber jede Fortsetzung zur Helizität 2.
5. **Wichtigste Spannungen:** W1 (Zustand nicht P₀-invariant), W4 (fermionisches Determinantenwort), W5 (nur eine
   vertauschende Translation in der Nahtfamilie), W7 (relativer Index 3 im analytischen Einheitswindungsfall ausgeschlossen),
   W8 (D₄-Gleichheit gegen den gelösten Flavourblock).

## 1. Behauptung (1): P₀ parabolisch, L₀ elliptisch

**Vorbestand.**
- Ledger Z. 496 (SEAM.EQUIV.BW.HSMI.01, [E]/[C]), `verification/v438_seam_hsmi_borchers.py:30-34`:
  - „in the unitary lowest-weight sl(2,R) rep the line translation P_line = L0 − (L+ + L−)/2 is POSITIVE“.
  - Das ist dieselbe parabolische Richtung, nur mit dem anderen Fixpunkt.
  - Papier: `documents/tfpt_research_contracts/sections/003-source.tex:298-304`.
- Ledger Z. 503 (FORM.SEAM.BW.HSMI.01, [F]): Die Standardpaar-Relationen sind im Lean-Kern bewiesen. Borchers-Schritt und
  Kontinuums-Standardpaar sind benannte Axiome (`borchers_wiesbrock`, `seam_continuum_standard_pair`).
- PSL2-Typen und Cayley-Karte:
  - `experiments/theory-contracts/original-paper-source-reconciliation-20261102/hypotheses_20261105E/fieldint_hypotheses2/TIME-IDENTIFICATION-DECISION-20261105.md`
    §3: Drehung elliptisch, Translation parabolisch, Dilatation hyperbolisch; ∂_θ wird zu (1 + x²)∂_x/2.
  - Atlas `finding:hyp20261105e-free-quarter-flow-is-not-faithful-line-translation` (established_scoped): Der periodische
    Quarterfluss ist keine treue Linientranslation. Zitiert sind dort CLTW (arXiv:1101.2865) und Longo–Tanimoto
    (arXiv:1608.08903).
- Rotations-Gibbszustand ≠ Grundzustand:
  `documents/articles--2026-08-29--holomorphic_kms_extension_en/sections/001-section-setting-and-correction-of-terminology.tex:79-116`
  („The geometric KMS state on the real line for translations/dilations and the rotational Gibbs state on the circle are
  likewise different constructions“).
- Nur eine lichtartige Translation:
  - z59 (`VERDICT.md:365-369`) und z32 (`VERDICT.md:1571-1577`).
  - Atlas `finding:original-e8-overlap-parabolics-and-finite-region-closure-20261112` (established_scoped), Quelle
    `TFPT-REGIONALE-KOMPOSITION-20261112.txt:5-89`: Paraboliken verschiedener Endpunkte haben eine nichtverschwindende Klammer;
    sl(2, ℝ) hat höchstens eindimensionale abelsche Unteralgebren.
- Halbseitige Inklusionen für den Originalzustand:
  - z32 (`VERDICT.md:1571-1573`): Ein Gibbszustand mit Punktspektrum hat eine fastperiodische Modulgruppe, also keine echte
    halbseitige Inklusion.
  - Atlas `finding:original-four-cuts-no-local-relative-commutant-or-half-sided-interpolation-20261111` (excluded_scoped,
    quartersymmetrischer Fall).
  - Atlas `finding:e8-relative-sewing-needs-reference-transport-and-absolute-modular-time-is-two-sided-20261111`
    (excluded_scoped): Die absolute Modularzeit hat zweiseitiges Spektrum.
  - H6h22-Tabelle (`VERDICT.md:1079`): „Translationen kommen nicht aus 1D-Inklusionen, sondern aus der Matrixvariable (z40)“.
- H6h34: z70 (der periodische Modulfluss ist keine Zeittranslation), z71 (lokale Modulzeit dicht; lokale
  Translationsbeziehungen fehlen).

**Eigene Einordnung.**
- **P₀ ist keine Symmetrie des Originalzustands.**
  - Der Originalzustand ist KMS für L₀ bei β = 4π (z58) und auf dem Fockraum normal.
  - Wäre er auch P₀-invariant, dann invariant unter der von L₀ und P₀ erzeugten Gruppe. [L₀, P₀] = (L₋₁ − L₁)/2 ≠ 0, also
    erzeugen beide die ganze sl(2).
  - Eine möbiusinvariante Dichtematrix sitzt auf den möbiusinvarianten Vektoren (L₀ = 0), also höchstens auf einem Vakuum.
  - In den Quartersektoren ist das kleinste L₀ = h = 1/32 > 0 (z59). Dort gibt es gar keinen P₀-invarianten Vektor; ein
    Vakuum gibt es nur in einem Sektor mit h = 0, etwa im Vakuummodul des E₈-Netzes.
  - In P₀-Zeit wäre der Originalzustand also nicht stationär. Der Abstand zum Vakuum ist klein, aber nicht null:
    ω(J₋₁J₁) = 3,5·10⁻⁶; die niedrigste Quartermode ist zu 4,14 % besetzt.
- **P₀ wählt einen Punkt der Naht.**
  - L₀ + (e^{iφ}L₁ + e^{−iφ}L₋₁)/2 fixiert genau den Nahtpunkt −e^{−iφ}; P₀ fixiert die Marke −1.
  - Die vier Marken geben vier positive Paraboliken P^(k). Sie vertauschen nicht, und ihr μ₄-Mittel ist exakt L₀
    (Σ i^k = 0; numerisch bestätigt für h = 1/32, 1/2, 1).
  - P₀ bricht damit die Clockinvarianz (QGEO.SYM.01). Die Wahl des Fixpunkts ist die erste neu gewählte Geometrie im Sinn von
    AGENTS.md.
- **Typ:** Lichtstrahltranslation der chiralen 1D-Naht. So nennt sie auch Ledger Z. 496 („positive-energy lightlike
  translation“).

**Ergebnis.** Vorbestand. Als physikalische Zeit ist P₀ erst zulässig, wenn Zustand, Fixpunkt und 3+1-Typisierung benannt
sind.

## 2. Behauptung (2): bilokale Teilchen-Loch-Auslese

**Vorbestand.**
- z71 Teil 2 (`ausgabe_z71.txt:7-11`, `z71_lokale_zeit_4d.py:78-116`):
  - W4(τ, γ) = Σ_m (sin mγ/sin γ)(e^{−mτ} + e^{−m(β−τ)})/(1 − e^{−βm}); W4(τ, 0) ist die Stromantwort aus z58.
  - Die S³-Identität ∂²_γ|₀ = −(1/3)(m² − 1)m gilt; längs des Nahtwinkels gilt die 2D-Wellengleichung.
  - Befund: „Die 4D-Fortsetzung braucht eine Winkelvariable, die das Original nicht hat“. Atlas
    `finding:h6h34-local-modular-time-dense-4d-continuation-not-forced-origin-two-binary-20261124` (conditional).
- Stromantworten ω(J₋ₘJₘ) = m/(e^{4πm} − 1): z58 (`VERDICT.md:299`) und Atlas
  `finding:elliptic-current-state-differs-from-hardy-on-every-proper-current-interval` (established_scoped).
- Teilchen-Loch-Bausteine:
  - Atlas `finding:p1-gamma-null-word-rank-one-particle-hole-20261112` (established_scoped).
  - z63 Teil 4 (iii′): Das Teilchen/Loch-Dublett gibt den Spurzustand (KMS nur bei β = 0).
- Kriterium Spektraldimension: z60 (`VERDICT.md:225-250`).
- Einen bilokalen Relativwinkel als Variable führt der Bestand nirgends (Suche nach „bilokal“, „bilocal“, „Teilchen-Loch“
  zusammen mit sin(mγ): kein Treffer).

**Eigene Rechnung (Scratchpad) und Einordnung.**
- **Die Identität ist exakt.**
  - Für jeden eichinvarianten quasifreien Zustand mit c(x) = 1/(1 + e^{βx}) auf x ∈ ℤ + α gilt
    c(x)(1 − c(x + m)) = [c(x) − c(x + m)]/(1 − e^{−βm}).
  - Daraus folgt D_m(γ) = D₁(γ) · (sin mγ/sin γ) · (1 − e^{−β})/(1 − e^{−βm}).
  - Im Vakuum ist der letzte Faktor 1, bei β = 4π ist er 1 − O(e^{−4π}). Geprüft für α = 1/4, 3/4, 1/2 und β = ∞, 4π, 2
    (relativer Fehler 10⁻¹⁴).
- **Der Faktor ist kinematisch.**
  - sin(mγ)/sin γ kommt aus der Fensterbreite m der Paarlagen (e^{imγ} − e^{−imγ}) und hängt nicht vom Zustand ab.
  - Es ist der SU(2)-Charakter zum Spin (m − 1)/2. Auf S³ ≅ SU(2) sind die Zonalfunktionen genau diese Charaktere
    (Peter–Weyl, bekannte Mathematik). Die Übereinstimmung belegt also keine S³-Geometrie.
- **Im Original trägt das Einzelpaar eine Spinphase.**
  - In den Quartersektoren ist D₁(γ) = e^{∓iγ/2} (Vakuum, α = 1/4 bzw. 3/4). Die reine S³-Zonalform steht nur im Verhältnis
    D_m/D₁.
  - Die Summe beider Quarter ist sin(mγ)/sin(γ/2): die Zonalfunktion vom Grad 2m − 1 beim halben Winkel γ/2, also nur
    ungerade Grade (geprüft, Fehler 10⁻¹⁵).
  - Der ungetwistete Sektor α = 1/2, in dem D₁ konstant wäre, kommt im Original nicht vor (z46, z58).
- **Die Zählung ist zweidimensional.**
  - Auf Stufe m gibt es im Vakuum genau m Paare. Σ_{m≤M} m ~ M²/2 ist die Weyl-Zählung des Torus der Paarlagen S¹ × S¹
    (Spektraldimension 2). S³ braucht m².
  - Das passt zu z44: Die Paarkanäle allein schließen zu sp(4, ℝ), also 2+1.
  - Bekannte Mathematik: D_{1/2} ⊗ D_{1/2} = ⊕_{s≥1} D_s. Die m Paarzustände der Stufe m sind die Moden W^(s)_{−m}Ω
    (s = 1, …, m) der bilinearen W_{1+∞}-Ströme.
- **Typ.** γ ist der Relativwinkel zweier Nahtpunkte, keine Richtung quer zur Naht. Ihn als geodätischen S³-Winkel zu lesen,
  ist die erste neu gewählte Geometrie.

**Ergebnis.** Teilweise neu. z71 ist in der Begründung zu korrigieren: Eine Winkelvariable ist nativ vorhanden. z71s
Ergebnis bleibt: 4D ist nicht erzwungen.

## 3. Behauptung (3): m² als K-Typen von SU(2,2), vier Translationen

**Vorbestand.**
- z40 (`VERDICT.md:1230-1248`; Atlas `finding:h6h21-matrix-kernel-nu-one-massless-spacetime-families-20261124`, conditional):
  - „Der Hardy-Kern ist der erste Wallach-Punkt, die Leiterdarstellung von U(2,2) des masselosen Skalars. Das ist bekannte
    Mathematik“.
  - Vier unabhängige Hardy-Translationen reichen nicht; nötig ist ein gemeinsamer Impuls im Kegel Herm(2)₊.
- z41 (`VERDICT.md:1106-1115`): Laplace-Transformation von det(P)^{s−2} auf Herm(2)₊. Bei s → 1 sitzt das Spektralmaß auf
  m² = 0, also P² = 0 und P₀ ≥ 0 für ν = 1. Die Wallach-Menge ist die 4D-Unitaritätsschranke.
- z42 (`VERDICT.md:1012`): Zwei Bosonenpaare geben 1/det(I − W†Z).
- z44 (`VERDICT.md:914-917`; `z44_uebergangsfamilie.py:151-243`; Atlas
  `finding:h6h24-pair-algebra-sp4-su22-spinc-five-three-higgs-z-from-hopf-phase-20261124`, conditional):
  - Die vier Markenpaarkanäle schließen zu sp(4, ℝ) (2+1), erst mit der einseitigen Blattaktion L(A) = a†Aa zu su(2,2) (3+1).
  - Die Oszillatoren a, b sind bosonisch. Scope: „Dass die einseitige Blattaktion physikalisch ist, bleibt Lesart.“
- TFPT-STAND:527-529: Vier positive, kommutierende Translationen sind offen; „freier Kandidat ist der Matrixkern bei ν = 1“.
- Kegel ohne Translationen:
  - `experiments/theory-contracts/compiler-cone-object-audit/README.md` („noch keine Ereignisraumzeit“).
  - Atlas `finding:native-clock-orientation-produces-four-dimensional-state-cone-20261111` (established_scoped).
  - Atlas `finding:originlab-operator-positive-pauli-without-commuting-translations-20261113` (excluded_scoped): Positivität
    allein erzwingt keine vertauschenden Translationen.
- Standardform: Ledger Z. 790 (E8.STANDARDFORM.CHART.01), z62, z63 (ρ_m auf End(V); Modulfluss = Drehung).
- T′:
  - z64/z70: Die Signatur (1,3) entsteht nur über V̄*. Atlas `premise:h6-conjugate-partner-doublet-t-prime-20261124`.
  - z65 zeigt nur die Richtung „Twistor-Realstruktur ⇒ Tomita“. VOGELPERSPEKTIVE §4 (1) und §6 führen „T′ aus Tomita“ als
    offen.
- ℝ × S³ im Bestand: Ledger Z. 737 (CCC.SEAM.KINEMATICS.01, [E] Kinematik) mit der Lorentzbrücke (1/cos²T)(−dT² + g_S³) und
  der Linse S³/ℤ₄ als Kreuzung. Die Firewall nennt genau diesen Uplift (v957) als nicht schließend.

**Eigene Einordnung.**
- **Zählung richtig, aber nur Charakter.** V_j ⊗ V_j mit j = (m − 1)/2 hat m² Zustände (Peter–Weyl). Das ist eine
  Dimensionsübereinstimmung, kein Operatorwörterbuch (AGENTS.md; SEAM.BULK4D.RECON.01 B7 „not character-level only“;
  WOIT-Kill-Test (5) „right character but the wrong OPE“).
- **Die su(2) auf den m Paaren ist nicht benannt.**
  - Die Gewichte 2a − m + 1 stimmen.
  - Ein Heber bei festem m muss Teilchen- und Lochniveau zugleich verschieben. Er ist kein Fermionbilinear, sondern quartisch
    mit niveauabhängigen Koeffizienten (Jordan–Schwinger mit Hardy-Niveaus als Bosonquanten).
  - In W_{1+∞}-Sprache mischt er Ströme verschiedener Spins s.
- **Die Tomita-Verdopplung ist T′.**
  - Zweite SU(2) und Signatur (2,2) sind genau die einseitige Blattaktion aus z44 bzw. das Partnerdublett V̄*.
  - Offen bleibt, dass die modulare Konjugation des Originalzustands diesen Partner mit eigener U(2) liefert.
  - Die nichtkompakten Erzeuger (P_μ) verbinden beide Blätter und liegen in keiner Blattalgebra.
- **Statistik.** Leiterdarstellung und det(I − ZW†)^{−1} sind bosonisch. Fermionische Paare geben im Zweipaarsektor ein
  Determinantenwort ≠ 0 (W4). Für den Einpaarsektor ist das gleichgültig, für den Fockraum eines freien 4D-Felds nicht.
- **Zustand.** P₀ ≥ 0 und P² = 0 gelten in der Vakuumdarstellung. Der Originalzustand ist dort ein Gibbszustand des konformen
  Hamiltonoperators und nicht poincaréinvariant (wie W1).
- **Linse.** Ist die räumliche Sphäre die Linse S³/ℤ₄ des Bestands (Z. 737, CELEST), gilt m² nicht. ℤ₄ = ⟨diag(i, −i)⟩ lässt
  auf Stufe m nur die Gewichte μ ∈ 2ℤ: m = 1…6 geben 1, 0, 3, 0, 15, 0 statt 1, 4, 9, 16, 25, 36 (eigene Rechnung). Die
  Behauptung muss S³ oder S³/ℤ₄ festlegen.
- **Ertrag.** Eine freie Kinematik. SEAM.BULK4D.RECON.01 verlangt dazu B6 (wechselwirkende S-Matrix) und B7. WOIT-Kill-Test (4)
  („nur freie Zustände“) ist hier der erwartete Ausgang.

**Ergebnis.** Vorbestand bis auf die Lesart „Teilchen-Loch ⊗ Tomita = V_j ⊗ V_j“. Last tragend sind T′, bosonische Paare und
das Vakuum statt des Originalzustands. Für das Kriterium aus z60 ist das ein Kandidat für eine Antwort mit Weyl-Zählung 4, aber
nur als Lesart.

## 4. Behauptung (4): Eindeutigkeit des Matrixkerns

**Vorbestand.**
- z42 (`VERDICT.md:1012-1014`): K + ε det Z · det W̄ ist ebenso positiv und auf Rang 1 gleich; „Der Kern braucht also eine
  weitere Regel.“ Ebenso TFPT-STAND:320-321.
- z43, OriginLab `process_rules` (`VERDICT.md:1030-1031`): „Hardy plus Wellengleichung oder Kovarianz wählt ν = 1, ε = 0.
  Hardy allein lässt den det-Term offen.“
- z44 (`VERDICT.md:916-917`): Das Determinantenwort hat in K + ε det Z det W̄ das Normquadrat 4ε, im Hardy-Kern 0. Die
  fermionische Abkürzung f_i†f_{j+2}† gibt 4.
- z40/z41: Cayley-Identität (Wellengleichung genau bei ν = 1 und 0) und Wallach-Menge {0} ∪ [1, ∞).

**Einordnung.**
- **Struktur (bekannte Mathematik, Hua/Schmid).** Der Polynomraum auf M₂(ℂ) ist unter U(2) × U(2) multiplizitätsfrei. Jeder
  bi-U(2)-invariante positive holomorphe Kern ist Σ c_m s_m(ZW†) mit c_m ≥ 0, also eine Klassenfunktion F(λ₁, λ₂) der
  Eigenwerte von ZW†.
- **Hardy allein reicht nicht.** Rang-1-Schnitte legen nur F(λ, 0) fest. Gegenbeispiel (Scratchpad, wie z42):
  det(I − ZW†)^{−1} + ε det Z det W̄ ist positiv, bi-U(2)-invariant und auf Rang 1 gleich Hardy.
- **Eindeutig mit der Paarauslese.** Auf Diagonalmatrizen (Teilchen bei z₁, Loch bei z₂) ist die Vakuum-Paarauslese das
  Produkt zweier Szegő-Kerne, Π(1 − λᵢ)^{−1} = det(I − ZW†)^{−1}. Sie legt F fest, also ν = 1 und ε = 0. Der
  Strom-Zweipunkt ω(J₋ₘJₘ) allein ist eine Funktion einer Variablen und reicht dafür nicht.
- **Last tragend ist die bi-U(2)-Kovarianz.** Im Bestand ist die zweite U(2) die einseitige Blattaktion (z44, Lesart) bzw.
  der Partner in Hom(V, V′) (z42). Das ist dieselbe Voraussetzung wie T′.
- **Zustand.** Der holomorphe Kern gehört zum Hardy-/Vakuumzustand. Die Auslese des Originalzustands ist thermisch und nicht
  holomorph: Atlas `finding:elliptic-current-state-differs-from-hardy-on-every-proper-current-interval` und
  `finding:cover-green-thermal-state-is-not-original-quarter-hardy-vacuum` (beide established_scoped). „Passend zur
  Stromauslese“ kann nur die Vakuumauslese meinen.

**Ergebnis.** Als Satz auf der Matrixkugel zulässig; im Kern Vorbestand (z42, z43). Neu ist die Paarauslese als „weitere
Regel“ statt der Wellengleichung. 3+1 folgt daraus erst mit einem Wörterbuch Z ↔ Raumzeitpunkt (Cayley-Abbildung Kugel →
Röhre).

## 5. Behauptung (5): E₈ nicht als freie 4D-Felder

**Vorbestand.**
- Atlas `finding:original-e8-pure-null-ray-copy-stops-continuous-transverse-motion-20261112` (excluded_scoped), Quelle
  `TFPT-REGIONALE-KOMPOSITION-20261112.txt:91-119`: Ein exakter, zustandstreuer Einbau des ganzen E₈-Netzes auf einer
  Nullgeraden verträgt keinen stark stetigen Quertransport.
- Atlas `finding:original-e8-pointwise-transverse-profile-negative-norm-20261112` (excluded_scoped) mit
  `premise:pointwise-level-one-affine-profiles-20261114`, Quelle ebd. :121-173: Mit punktweisen Querprofilen und
  Niveau-eins-Relation ist die Zweistromnorm minus der Profilvarianz (r = cos θ gibt −1/2, ρ = 1 + cos θ gibt −17/4).
- z60 (`VERDICT.md:239-245`): Die E₈-Naht hat d = 8, kompakt und starr; keine (1,1)-Moduli.
- z40 (`VERDICT.md:1209-1214`): Niveau 1 braucht reelle (Majorana-)Eltern.
- Fortsetzung im Bestand nur wechselwirkend: SDYM(E₈) + BCOV auf dem Twistorraum (Ledger Z. 585 WOIT.OS.TWISTOR.01;
  `documents/tfpt_research_contracts/sections/010-source.tex:187-192`, Costello arXiv:2111.08879).
- Firewall-Routenordnung (Z. 739): Die Naht bleibt internes/topologisches Randmodul. Ledger Z. 773 (BULK.PRIME_COMPLETION.01):
  Die 2D-Eindeutigkeit des holomorphen (E8)₁-Netzes ist kein 4D-Satz.

**Einordnung.**
- Kein Widerspruch. Für punktweise Lokalisierung und punktweise Querprofile ist die Aussage schon ausgeschlossen.
- Nicht abgedeckt sind ausgedehnte Profile, distributionelle Grenzobjekte, operatorwertige Rückwege und ein geänderter
  Längsgrad (Scope der Atlas-Befunde). Die externe Fassung „freie 4D-Felder“ muss sich an diesem Scope messen.
- Bekannte Physik: Majorana–Weyl-Spinoren gibt es in Lorentzsignatur nur für d ≡ 2 mod 8. Die 16 reellen chiralen
  Nahtfermionen haben in 3+1 kein freies Gegenstück; das Gitter-Vertexbild braucht chirale 2D-Bosonen.
- Verträglich mit (3): Im freien 4D-Skalar ist die Lichtebenen-Stromalgebra transversal δ-normiert. Nahtströme mit Niveau 1
  können dort nur als ausgedehnte Profile vorkommen.

**Ergebnis.** Vorbestand (12.11. und z60), kein Konflikt.

## 6. Behauptung (6): drei Familien und D₄-Gewichte

**Vorbestand.**
- Drei Familien:
  - Atlas `finding:h6h8-mark-index-gives-three-families-and-review-corrections-20261123` (established_scoped):
    H⁰(K(D)) = 3 mit den Charakteren {−1, i, −i}.
  - z58 (`VERDICT.md:351`): E₋ = O(−1) + 3O, Index 3 ohne Spiegel. Atlas
    `finding:h6h31-chirality-from-index-without-mirror-gap-20261124` (conditional).
  - TFPT-STAND §3 (`docs/research/TFPT-STAND.md:71-74`): „H¹ vom Rang 3“, [E] für die Beziehung.
- Relativer Index:
  - z30 (`VERDICT.md:1644-1646`): Toeplitz-Index = −Windung.
  - z35/z36 (`VERDICT.md:1404`): Tr(P_C − U P_C U*) = 1.
  - Atlas `finding:trace-class-pair-index-fixes-relative-factor-without-primitive-ray` (established_scoped).
- Hindernis: `research/originlab/FAMILY-INDEX-20261114.md:1-90` und Atlas
  `finding:originlab-unit-winding-fixed-s4-analytic-family-obstruction-20261114` (established_scoped):
  - Für U = z^{k0}Q + z^{k3}P mit Einheitswindung auf einer analytischen Scheibe gilt n₀ + 3n₃ = 1; einen relativen
    Familienindex 3 gibt es dort nicht.
  - Dieselbe Determinante, dieselben vier Markenwerte und dieselbe Gram erlauben die relativen Indizes 0 und −3.
  - Reflexionspositivität liefert die Windungsphase nicht.
- D₄:
  - z32 Teil 5 (`VERDICT.md:1601-1603`): D₄ lässt zwei invariante Formen zu.
  - Atlas `finding:four-mark-invariant-norm-and-positive-time-nonuniqueness-20261123` (excluded_scoped): D₄-Sektoren
    1+1+1+2; D₄ und Positivität legen die Gewichte nicht fest.
  - z67 (`ausgabe_z67.txt`, Teil 1–2): gleiche Gewichte ½ auf den Familien i und −i. Atlas
    `finding:h6h34-h6-spin2-state-gives-no-flavour-weights-compatibility-20261124` (conditional).
- Gelöster Block: TFPT-STAND §3; Ledger Z. 17 (FLAV.TH12.01), Z. 361 (FLAV.TH13.01), Z. 476 (FLAV.PMNS.CLOSE.01); Atlas
  `finding:h6-common-parent-clock-full-algebra-20261124` (|U_e1|² = 0,6772).

**Einordnung.**
- „Drei“ ist Vorbestand. Zulässig ist der Index als Divisor- bzw. Residuenindex: h⁰(K(D)) − h⁰(K) = 3, also Marken − 1. Als
  Windungssymbol mit Einheitswindung auf einer analytischen Scheibe ist er bedingt ausgeschlossen. Die externe Fassung muss
  das Vergleichssymbol nennen.
- „Zwei gleiche Gewichte“ ist exakte Darstellungstheorie: Unter D₄ zerfallen die Familien in 1 ⊕ 2, die Spiegelung tauscht
  i ↔ −i. Das Verhältnis der beiden Gewichte bleibt frei.
- Keine Lücke: Gleiche Gewichte dürfen weder als Massen noch als exakte Mischung gelesen werden (W8).

**Ergebnis.** Vorbestand; verträglich nur als D₄-Gerüst unter dem gelösten Flavourblock.

## 7. Firewalls gegen (1)–(4)

| Firewall | (1) P₀ | (2) Bilokal | (3) m², P_μ | (4) Kern |
|---|---|---|---|---|
| DIMENSION.UPLIFT.FIREWALL.01 (Z. 739) | blockiert „Zeittranslation“ als 3+1 (F1, F2) | blockiert „transversale S³-Antwort“ als 3+1 | blockiert „4D-Skalar, vier Translationen“ als Schluss | Satz zulässig; blockiert die Lesart Z = Raumzeitpunkt ohne Wörterbuch |
| GRAV.NONCIRCULAR.01 (Z. 747) | blockiert, sobald P₀ aus Modulfluss/KMS gewonnen oder für BW benutzt wird | – | blockiert „KMS ⇒ SU(2,2)-Boosts ⇒ BW“ | – |
| GRAV.WEINBERG.WITTEN.01 (Z. 770) | – | – | nicht (Helizität 0); bindend für Helizität 2 | nicht; bindend für den Sym⁴-Kern aus z41 |

**Weitere bindende Zäune.**
- SEAM.BULK4D.RECON.01 (Z. 740): B6 und B7 („not character-level only“); Kill-Test „every 4D extension is generalized-free“.
- WOIT.OS.TWISTOR.01 (Z. 585): Kill-Tests (4) nur freie Zustände, (5) richtiger Charakter mit falscher OPE, (7)
  Markeninzidenz; der Scope-Zaun schließt beide Helizitäten und lokale Materie aus.
- DIMENSION.SELECTOR.4D.01 (Z. 755): „mu4 ⇒ 4D“ bleibt verboten. Das trifft (1) und (3), wenn die vier P^(k) oder die
  Tetraeder-Nullrichtungen aus z59 die Vier liefern sollen.
- E8.STANDARDFORM.CHART.01 (Z. 790): „no physics claim; exact renaming“. Das trifft (3) und (4), falls „Standardform“ diese
  Karte meint.
- CELEST.SEAM.01 (Z. 556): Die Kontinuumsexistenz des (E8)₁-Netzes ist als Eingabe unzulässig (Zirkularität).
- CHIRAL4D.NOMIRROR.01 (Z. 746) mit F1: Ein Nahtindex für (6) schließt keine 4D-Chiralität (C5).

**Zulässige Typisierung.**
- **(1)** P₀ ist die positive Lichtstrahltranslation der Möbiusdarstellung der chiralen Naht im Vakuum- bzw. Fockraum.
  - Hypothese und Schluss sind 1+1D; der Fixpunkt ist eine gewählte Marke; P₀ ist keine Symmetrie des Originalzustands.
  - Als Zeit der Raumzeit ist P₀ eine ANNAHME. Ein RESULT wäre es erst über ein bewiesenes Kontinuums-Standardpaar (Z. 503,
    benannte Axiome) in einem Zustand mit geometrischer Bogenmodularität. Der Originalzustand ist das nicht (z32).
- **(2)** Bilokale Auslese auf dem Torus der Paarlagen (1+1D). Schluss: Das Verhältnis D_m/D₁ ist der SU(2)-Charakter und
  gleich dem S³-Zonalkern; die Weyl-Zählung ist 2.
- **(3)** Unter T′ (Partner mit eigener U(2)) und bosonischer Paarstatistik trägt der verdoppelte Einpaarsektor eine unitäre
  SU(2,2)-Darstellung positiver Energie mit den K-Typen V_j ⊗ V_j.
  - Das ist freie 3+1-Kinematik (Lesart L); der Zustand ist das Vakuum.
  - SEAM.BULK4D.RECON.01 bleibt [O].
- **(4)** Satz auf der Matrixkugel: Unter bi-U(2)-Invarianz legen Hardy (Rang 1) und die Vakuum-Paarauslese (Diagonale) den
  Kern det(I − ZW†)^{−1} fest. Dimensionsneutral; die bi-U(2)-Invarianz ist Prämisse.

## 8. Widersprüche und Spannungen (präzise)

Kein etabliertes Ergebnis widerlegt die Mathematik von (1)–(6). Die Punkte betreffen Lesart und Geltung. W1, W4, W5, W7 und
W8 wiegen am schwersten.

- **W1, (1) gegen den Originalzustand.**
  - Der Originalzustand ist KMS für L₀ (z58) und nicht der Hardy-Zustand (Atlas, established_scoped).
  - Ein L₀-KMS-Zustand ist nur im Vakuumfall P₀-invariant; in den Quartersektoren gibt es keinen P₀-invarianten Vektor
    (eigene Einordnung, Abschnitt 1).
  - Der Artikel vom 29.08. trennt Rotations-Gibbszustand und geometrischen Linien-KMS-Zustand ausdrücklich.
  - „P₀ ist die Zeittranslation“ und „der Originalzustand ist stationär“ schließen sich deshalb aus.
- **W2, (1) gegen die Modulroute.**
  - Für den Originalzustand gibt es in den geprüften Scopes keine halbseitige Inklusion: z32 (fastperiodische Modulgruppe),
    Viererschnitte (excluded_scoped), zweiseitige absolute Modularzeit (excluded_scoped).
  - P₀ lässt sich daher nicht als Wiesbrock-Translation (log Δ_N − log Δ_M)/2π aus dem Originalzustand gewinnen. (1) umgeht
    das als Gruppenelement; mit GRAV.NONCIRCULAR.01 ist das nur als deklarierte Annahme verträglich.
  - Bekannte Literatur (nicht im Repo, nicht nachgerechnet): Der Modulhamiltonoperator eines Bogens im thermischen Zustand
    eines freien chiralen Fermions auf dem Kreis ist nichtlokal (Blanco–Pérez-Nadal 2019; Fries–Reyes 2019). Das passt zu z71.
- **W3, (2) gegen die Begründung in z71 (conditional).**
  - z71 sagt, die Winkelvariable fehle im Original. Die bilokale Auslese hat sie nativ, mit exakter S³-Zonalform im
    Verhältnis D_m/D₁ (Abschnitt 2).
  - z71s Ergebnis „nicht erzwungen“ bleibt: Weyl-Zählung 2, kinematisches Profil, Spinphase der Quarter. Zu korrigieren ist
    nur der Wortlaut.
- **W4, (3)/(4) gegen z44 (conditional) und `finding:p1-gamma-null-word-rank-one-particle-hole-20261112` (established_scoped).**
  - Fermionische Paare geben im Zweipaarsektor das Determinantenwort 4 bzw. 4 Σ_{i<j} s_i² s_j² statt 0.
  - Leiterdarstellung und ν = 1 verlangen bosonische Paare. Bekannte Mathematik: Fermionisch entsteht det(I + ZW†), der
    kompakte Dual.
  - Teilchen-Loch-Paare der Naht sind fermionisch. Auf dem Einpaarsektor ist das gleichgültig, für den Fockraum eines freien
    4D-Felds nicht.
- **W5, (3) gegen die regionale Komposition (established_scoped) und den Nullgeraden-Ausschluss (excluded_scoped).**
  - Innerhalb der Möbiusfamilie der Naht gibt es nur eine vertauschende Translation.
  - Punktweise Quertranslationen des reinen Nahtnetzes sind ausgeschlossen.
  - Vier vertauschende P_μ gibt es nur in einer Erweiterung, mit ausgedehnten Profilen.
- **W6, (3) gegen die Linse in Ledger Z. 737 [E].** Die Kinematik des Bestands hat S³/ℤ₄; dort gilt m² nicht (Abschnitt 3).
- **W7, (6) gegen `research/originlab/FAMILY-INDEX-20261114.md` (established_scoped, bedingter Ausschluss).**
  - Ein relativer Familienindex 3 aus einem Windungssymbol mit Einheitswindung auf einer analytischen Scheibe ist
    ausgeschlossen.
  - Aus Determinante, Markenwerten und Gram ist er nicht identifizierbar.
  - Verträglich ist nur der Residuen- bzw. Divisorindex.
- **W8, (6) gegen den gelösten Flavourblock.**
  - D₄-gleiche Gewichte (a, b, b) sind ein symmetrisches Gerüst.
  - Die eingefrorenen Werte brechen es: |U_e1|² = 0,6772 ≠ 2/3 (Atlas common-parent), sin²θ₁₂ = 1/3 − φ₀/2 (Z. 17),
    δ = 240° statt ±90° bei μ-τ-Reflexion (Z. 476; Atlas `finding:h6h26-delta-6n2-residuals-only-tm2-tm1-fifteen-cp-prediction-20261124`),
    neun verschiedene Massen (TFPT-STAND §3).
  - Als Massen- oder Mischungsaussage widerspräche (6) dem gelösten Block. Als Aussage über das D₄-Modul ist sie verträglich.
- **W9, (4) gegen den Originalzustand.** Der holomorphe Kern ist der Vakuumkern. Der Originalzustand ist nicht Hardy
  (established_scoped), und eine bloß kovarianzerhaltende Einteilchen-Identifikation von thermischer Kovarianz und
  Hardy-Projektion ist ausgeschlossen (`finding:cover-green-thermal-state-is-not-original-quarter-hardy-vacuum`,
  established_scoped).
- **W10, (1)/(3) gegen DIMENSION.SELECTOR.4D.01 (Z. 755).**
  - Kommen die vier P_μ von den vier Marken (z59, Tetraeder), greift „mu4 ⇒ 4D verboten“.
  - In der Standardform ist die Vier dim_ℝ Herm(2) = n² mit n = 2. Dann ist ℂ Eingabe: Positivität wählt ℂ nicht
    (`experiments/exploration/2026-11-20-operational-origins/VERDICT_CAUSAL.md:27-45`), lokale Tomographie schon
    (quaternion-frames C).
- **Zu (5):** kein Widerspruch.

## 9. Wiederverwenden statt neu rechnen

1. **v438 und Lean-Kern** (Ledger Z. 496, Z. 503): Positivität der parabolischen Richtung und Standardpaar-Relationen. Die
   benannten Axiome zeigen genau, was für ein RESULT fehlt.
2. **Zeitentscheidung vom 5.11.** (`TIME-IDENTIFICATION-DECISION-20261105.md`, `check_original_time_type.py`, 16/16): PSL2-Typen,
   Cayley-Karte und die Literatur CLTW/Longo–Tanimoto.
3. **Regionale Komposition vom 12.11.** mit vier Atlas-Befunden: Paraboliken mit Klammer, eine vertauschende Translation,
   Nullgeraden-Ausschluss, negative Norm punktweiser Querprofile. Prüfstein für (1), (3) und (5).
4. **Artikel `holomorphic_kms_extension`** (29.08.): Rotations-Gibbszustand gegen Grundzustand; Eindeutigkeit rotationaler
   KMS-Zustände nach Longo–Tanimoto.
5. **z40–z44:** Wallach-Kern, Laplace-Spektralmaß auf Herm(2)₊ (P² = 0), Bosonenpaare, ε-Kern, `process_rules`-Klausel,
   sp(4, ℝ) gegen su(2,2), Determinantenwort. Damit sind (3) und (4) schnell prüfbar.
6. **z58, z59, z60, z71:** Stromantworten mit Bosefaktoren, eine lichtartige Translation, Spektraldimension als Kriterium,
   W4/W2 und die S³-Identität. z71 Teil 2 lässt sich direkt um die bilokale Auslese erweitern.
7. **z63:** ρ_m, Standardform, Modulfluss = Drehung, Teilchen/Loch-Dublett = Spurzustand.
8. **OriginLab:** `reconstruction_search` (Klausel „(1,3) nur mit V̄*“) und `process_rules` (Matrixkern-Prinzipien). T′ bzw.
   die Paarauslese lassen sich dort als Klausel ergänzen, statt neu zu suchen.
9. **Für (6):** `FAMILY-INDEX-20261114.md`, Z10, z58, z32 Teil 5, z67 und der Atlas-Befund zu D₄ und Positivität.
10. **Für (5):** Ledger Z. 585 und der Vertrag `010-source.tex:187-192` (E₈ nur als SDYM auf dem Twistorraum).

## 10. Grenzen dieser Prüfung

- Die externe Analyse liegt nicht vor. Geprüft ist die Fassung des Auftrags; ihre Herleitungen sind nicht gelesen.
- Eigene Einordnungen sind W1 (Invarianzargument), die Paraboliken der Marken, die Bilokal-Identität mit Spinphase, die
  Linsenvielfachheiten und der Gegenbeispielkern. Gerechnet ist nur die numpy-Kurzprobe im Scratchpad; nichts davon ist
  abgelegt.
- Bekannte Mathematik/Physik (Peter–Weyl, Hua/Schmid, D_{1/2} ⊗ D_{1/2}, Majorana–Weyl, Literatur zu thermischen
  Modulhamiltonoperatoren) ist hier nicht nachgerechnet.
- Die meisten der 615 Theorieverträge und die 152 nicht eingearbeiteten Atlas-Versuche sind nicht einzeln gelesen. Atlas-Befunde
  sind mit ihrem Status zitiert, ihre KI-Zusammenfassungen nicht geprüft.
