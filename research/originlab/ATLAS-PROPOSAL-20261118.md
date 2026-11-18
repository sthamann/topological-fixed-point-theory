# Vorschlag für einen Atlas-Eintrag · Sitzung 17./18. November 2026

Status: Vorschlag. Im Atlas ist nichts geschrieben, kein Ledgerstatus
ist berührt. Atlasrevision beim Lesen: `c6e33e8fedf4cc2cb7bb3de1`
(lokal, erzeugt am 11. November; der lokale Stand liegt hinter
`origin/main`).

## Vor dem Eintragen zu erledigen

1. **Atlasstand abgleichen.** Der Eintrag prüft die Revision. Der
   lokale Atlas ist älter als der Remote-Stand; ohne Abgleich würde der
   Eintrag abgelehnt oder auf einen veralteten Stand geschrieben.
2. **Quellen registrieren.** Erledigt bis auf den Commit: Die
   Prüfskripte liegen jetzt in
   `experiments/theory-contracts/originlab-four-marks-session-20261118/`
   (von Git nicht ausgeschlossen), mit `MANIFEST.txt` (SHA256 und
   Ergebnis jedes Laufs) und den Ausgaben in `output/`. Alle 30 Skripte
   enden dort ohne Fehler. Der Atlas verlangt verfolgte Quellen; dafür
   fehlt noch der Commit dieses Verzeichnisses und von
   `research/originlab/THEORY-BRIEF.md`.
3. **Eigenes Problem anlegen.** Das thematisch passende Problem
   `problem:atlas-clock-defect-and-operational-reconstruction` hat 250
   Befunde und über eine Million Zeichen; ein vollständiger
   Autorendatensatz lässt sich über die Schnittstelle nicht in einem
   Stück holen, und Befunde dürfen nicht verschwinden. Vorschlag: ein
   neues Problem mit Verbindungen dorthin.
4. **Keine Promotion.** Die Befunde betreffen die freie Quelle und
   aufgezeichnete Objekte. F1 bleibt offen.

## Vorgeschlagenes Problem

- **ID:** `problem:four-marks-quartic-and-free-seam-source`
- **Frage:** Welche exakten Beziehungen bestehen zwischen dem System
  aus v117, dem aufgezeichneten Torus über den vier Marken, der freien
  elliptischen Quelle und der aufgezeichneten Markenmatrix, und welche
  Voraussetzung trennt die freie Naht vom E8-Netz?
- **Themen:** `topic:e8`, `topic:clocks`, `topic:origin`
- **Verbindungen:** zum oben genannten Problem; zu
  `finding:geometric-square-clock-is-not-theta-v117-quarter-generator`
  und
  `finding:theta-square-modular-pair-reconstructs-v117-with-distinct-geometric-q`
  (Befund 2 unten ergänzt beide); zu
  `finding:actual-elliptic-gaussian-source-determines-full-quarter-current-functional`.
- **Voraussetzungen (bestehende IDs):** `premise:typed-source-structures`,
  `premise:free-chiral-cauchy-quantization-of-canonical-elliptic-annulus`.

## Vorgeschlagene Befunde

Jeder Befund nennt Status, Aussage, Geltungsbereich und Beleg. Die
Skripte liegen in `experiments/theory-contracts/originlab-four-marks-session-20261118/`;
alle bestehen ohne Fehler.

### Belegt im genannten Rahmen

1. **`v117-a0-system-is-the-root-system-of-a-quartic`**
   Aussage: Für die vier Wurzeln X von X⁴ + zX + q = 0 mit
   q³ = 27(z⁴−1)/256 löst der Vektor (y1, y2, y3) aus dem Skript das
   unveränderte A0*-System aus v117; die polynomiale Umeichung P mit
   det P = −(9√2/40)(z⁴−1) ist symbolisch nachgewiesen. Die
   Vertauschungsmatrix der Markenschleife ist symbolisch gleich M0.
   Geltung: Das System wird wie in v117 gegeben genommen. Die
   Vertauschung ist per Fortsetzung abgelesen. Keine Herkunftsaussage.
   Beleg: `check_original_system_closed_form.py`;
   Original: `verification/v117_monodromy_weyl_a3.py`.

2. **`quarter-turn-lifts-at-the-two-poles-of-the-v117-system`**
   Aussage: Die Hebung der Vierteldrehung ist U bei z = 0 und M0⁻¹U
   (Ordnung zwei, Spektrum 1, −1, −1) bei z = ∞; derselbe Automorphismus
   wirkt auf den vier Verzweigungspunkten als Viererzyklus und auf den
   drei Zweiteilungspunkten als Vertauschung zweier.
   Geltung: wie Befund 1. Ergänzt die beiden Befunde vom 26. Oktober,
   widerspricht ihnen nicht.
   Beleg: `check_level_four_against_original.py` (O7),
   `check_original_system_closed_form.py` (C3).

3. **`recorded-torus-marks-and-mirror-census-recomputed`**
   Aussage: Auf dem aufgezeichneten Torus y² = x⁴ − 1 liegen die Marken
   bei (±1/4, ±1/4); die sechzehn antiholomorphen Abbildungen haben die
   aufgezeichneten Ordnungen und Fixkreiszahlen.
   Geltung: unabhängige Nachrechnung von `EQUATORIAL-LIFTS.txt`; nichts
   Neues außer der Lage der Marken.
   Beleg: `check_seam_mirror_census.py`.

4. **`free-kernel-between-the-four-marks`**
   Aussage: Für die aufgezeichnete Spinstruktur ist der Kern zwischen
   gegenüberliegenden Marken null, zwischen benachbarten vom Betrag
   Γ(1/4)²/(2√(2π)) mit relativer Viertelphase; die Vierpunktfunktion hat
   für die drei Paarungen die Beträge 1, 2, 1 in Einheiten des
   Betragsquadrats.
   Geltung: freie Quelle auf dem aufgezeichneten Torus.
   Beleg: `check_seam_mirror_census.py` (K3),
   `check_free_source_at_marks.py`.

5. **`recorded-mark-matrix-carries-pairing-responses-and-e8-character`**
   Aussage: Für die aufgezeichnete Markenmatrix G (v484, v485, v512,
   Diagonale null) und den Vorzeichenvektor s einer Paarung gilt
   exp(−(π/2)sᵀGs) = 2^(−8/3)|Disc|^(−1/6)|θ/η|⁸; die Summe über die drei
   Paarungen ergibt χ_E8 = j^(1/3) am Modul des Torus, mit den
   aufgezeichneten Werten 12 und 148176/25. Die Markendeterminante für
   bewegte Marken ist Π(1 − ε c3 ln R_k).
   Geltung: Gleichheit von Funktionen der Marken zwischen dem
   aufgezeichneten skalaren Nahtkern und der freien Spinorquelle. Keine
   Herleitung von SEAM.EQUIV.01.
   Beleg: `check_four_mark_response.py`,
   `check_contact_diagonal_bridge.py`,
   `check_mark_matrix_determinants.py`.

6. **`weight-of-side-pairings-is-zero-or-one-under-recorded-axioms`**
   Aussage: A1 und A3 lassen 135 = 105 + 30 Codes; der Paarcode der
   Marken gibt θ3⁸ (freie Naht), der Hamming-Code E8; A2 entspricht der
   Invarianz unter τ → τ + 1 und wählt E8. Trennende Zahl zwischen
   Quadrat und verschobener Lage: 25/16 gegen 7/4.
   Geltung: Nachrechnung und Verbindung zu SEAM.CODE.TYPEII.01; der
   dortige Rest (die physische Naht erfüllt A1 bis A3) bleibt.
   Beleg: `check_weight_b_from_record.py`.

7. **`quarter-charge-law-unique-over-offsets-and-widths`**
   Aussage: Ein gemeinsames Gaußgesetz für gerade und ungerade Summen,
   Fourier-Invarianz der Summe und selbstduale Breite lassen nur die
   Verschiebungen (±1/4, ∓1/4) und β = 4π.
   Geltung: an je drei bis vier Werten geprüft, Begründung über die
   Nullstellen von θ2; bezieht sich auf das in R066 aufgezeichnete
   Gesetz.
   Beleg: `check_recorded_elliptic_source.py` (R6).

8. **`one-sided-e8-and-diagonal-coupling-coincide-on-the-mark-circle`**
   Aussage: Die Torusspuren unterscheiden sich um 2|θ3/η|¹⁶(Im λ)²; auf
   dem Markenkreis stimmen sie überein. Am vierten Punkt −i ist die
   rechte Stressantwort 0 beziehungsweise 1; das Wort bei (1, −1, −i)
   gegen (1, i, −1) hat gleichen beziehungsweise doppelten Betrag.
   Geltung: Modellvorhersagen; für die physische Quelle liegt keine
   Antwort vor.
   Beleg: `check_reflected_seam_couplings.py`,
   `check_fourth_mark_probe.py`.

### Bedingt

9. **`seam-graph-boundary-space-of-the-smooth-torus-spinor`**
   Aussage: Ist das Feld der Spinor des glatten Torus, dann ist der
   Randraum des Inneren der Hardy-Raum der zwei Innenzellen; an den
   Marken ist nichts zu wählen, der Rahmen dreht sich um ein Achtel.
   Voraussetzung: das Feld ist dieser Spinor. Die Vollständigkeit der
   Zerlegung stützt sich auf einen klassischen Satz und ist nicht
   nachgerechnet.
   Beleg: `check_seam_graph_boundary.py`.

10. **`rational-elliptic-surface-over-the-parameter-sphere`**
    Aussage: Die Fläche mit j = 1728(1 − z⁴) hat Fasern II an den
    Marken und I4 am Pol; 82 Schnitte (10 + 32 + 40); ihr zweiter
    Torsionskubik ist die Resolvente der Quartik aus Befund 1.
    Voraussetzung: Die Gleichsetzung ihres Gitters mit dem
    aufgezeichneten E8 ist eine Hypothese; die Mordell-Weil-Struktur
    stammt aus der Literatur.
    Beleg: `check_surface_sections_and_lattice.py`.

### Ausgeschlossen im genannten Rahmen

11. **`marks-are-not-level-four-cusps`**: Eine Schleife um eine
    Stufe-4-Spitze wirkt trivial; Ordnung drei gäbe es auf Stufe 12 an
    allen sechs Spitzen. Nur diese Lesart ist ausgeschlossen.
    Beleg: `check_level_four_against_original.py` (O2).
12. **`v117-is-not-the-period-system-of-the-surface`**: unendliche
    gegen triviale Monodromie an der entarteten Faser.
    Beleg: `check_original_system_closed_form.py` (C6).
13. **`quarter-blocks-are-not-arcs-or-mark-coordinates`**: Eine Sonde
    auf einem Viertelbogen liegt zur Hälfte in jedem Block; die
    Verschiebungen ±1/4 des Bestands sind Kennzahlen der Modengruppen.
    Beleg: `check_free_source_at_marks.py` (M7),
    `check_seam_mirror_census.py` (K2).
14. **`four-square-tori-glue-is-not-the-recorded-code`**: Prüfung der
    Hypothese aus der dritten Runde; Formen 8 gegen 12, Klebung nicht
    linear über den Gaußschen Zahlen.
    Beleg: `check_gluing_and_reflection.py`.

### Hypothese (nur als Hypothese einzutragen)

15. **`sign-ensemble-form-on-the-four-marks`**: Vorzeichen an den
    Marken, Gewicht exp(−(π/2)sᵀGs), nur neutrale Belegungen. Die Form
    ist positiv vom Rang drei, ihre Referenz ist der E8-Charakter; am
    Quadrat circ(1, −2/3, 1/3, −2/3). Annahmen: Marken tragen
    Vorzeichen, nur neutrale Belegungen zählen, die nativen U- und
    Y-Einfügungen sind solche Vorzeichen. Nicht aus P1 hergeleitet.
    Beleg: `check_sign_ensemble_form.py`.

16. **`source-as-purity-of-two-opposite-arcs`** (Runde 25): Die
    gewichtete reflektierte Form auf den vier Marken ist der Tausch
    zweier Kopien der Naht auf zwei gegenüberliegenden Bögen, Referenz
    Tr ρ_A². Exakt folgt daraus: Vierteldrehung an jeder Marke,
    Verschiebungen 1/4 und 3/4, die aufgezeichnete Funktion F_α,
    b = 4π, und die Spur θ3(i)/η(i) = √2 auf dem aufgezeichneten
    Torus. Unabhängig auf einem Fermionenring bestätigt
    (Reinheitsverhältnis (1 − x)^(−1/4) je Gitterfermion). Freie Naht
    X = 1/(1 − x), E8 X = (1 − x + x²)/(1 − x); am Quadrat 2 gegen 3/2;
    Koeffizient von x² gleich (248 + 8)/256. Annahmen: die P1-Randform
    ist dieser Tausch; der Nahtzustand ist rein. Die Marken U und Y
    haben in dieser Lesart noch keinen Platz.
    Beleg: `check_two_arc_purity.py`, Suchspur `f1_two_arc_purity.py`.

17. **`purity-functional-on-the-native-tangential-operator`**
    (Runde 25, bedingt): Auf dem unveränderten h_N (96 Fasern, Moden
    −2 bis 1, Viertelverschiebung) mit gefülltem negativem Spektrum
    als Zustand: ohne Potential je Faser Reinheit 1/2 an einem Ort,
    Verhältnis 9/4 für Nachbarn, 1 für gegenüberliegende Orte. Mit
    den Kopplungen [2,2,2]: log-Verhältnis 24,39 und 24,40 für
    Nachbarn, 2,979 und 3,020 gegenüber. Kein aufgezeichneter
    Ankerwert. Annahme: der Zustand ist das gefüllte negative Spektrum
    von h_N, das der Assembler nicht als Generator ausweist.
    Beleg: `check_native_sea_purity.py`.

18. **`sign-mark-carries-the-space-of-the-e8-word`** (Runde 26,
    Wörterbuch, kein unabhängiger Beleg): An einer Marke, an der eine
    Kopie ein Vorzeichen trägt, liegen die 2⁸ Grundzustände von acht
    Nullmoden; die gerade Hälfte sind die 128 Materiewurzeln. Mit U
    als Cartan-Nullmode und Y als Paarstrom (10,6) folgt die ganze
    Tabelle von E8WORD20261116, die Dreifamilientabelle stimmt
    Eintrag für Eintrag mit `results_word.json`. Gleiche Algebra wie
    die aufgezeichnete Rechnung in anderer Schreibweise.
    Beleg: `check_marks_u_y_dictionary.py`, Abschnitt Q3.

19. **`recorded-binary-cuts-are-exchange-and-quarter-turn`** (Runde
    27, ersetzt die Fassung aus Runde 26): Auf den Zuständen mit einer
    Linie sind die aufgezeichneten Operatoren aus
    `coherent_origin.py` Matrix für Matrix swap = S (Kopientausch),
    parity = P (Vorzeichen), phase = i·S·P (Vierteldrehung). Die
    aufgezeichneten Marken Z_U = Z_Y = 0, Z_UY = i η Z_0 gelten genau
    für eine Linie in einer Kopie; η sagt, in welcher. Es bleibt keine
    Phase. Der positive Lift (1 + B)/2 ist der Projektor auf einen
    Viertelblock. Gleiche Kopien geben q = 1; der aufgezeichnete
    Vertrag nimmt q = 1 an, die aufgezeichnete Suche zählt nur
    0 < q < 1 auf. Offen: η, q und das Wörterbuch zur Originalquelle
    (das Modell hat die Rolle einer arithmetischen Kontrolle).
    Beleg: `check_shared_line_phase.py`.

21. **`flat-density-is-forced-by-the-sixty-return-values`** (Runde 28,
    belegt im genannten Rahmen; Aussage der fünfzehnten Notiz, hier
    nachgerechnet): Mit Referenz Tr ρ² und Rückkehrwert Tr(ρ X† ρ X)
    erzwingen die 60 Ankerwerte 1/4 (128 Zustände) bzw. 1/6 (96)
    ρ = 1/128 bzw. 1/96, auch mit Kohärenzen. Summe der X†X gleich
    15 bzw. 10 mal Identität; Cauchy–Schwarz; die Ströme trennen und
    verbinden alle Zustände.
    Beleg: `check_flat_sector_and_rate.py`, Abschnitt G2.

22. **`signed-and-positive-marker-differ-in-sign-for-electric-charge`**
    (Runde 28, belegt im genannten Rahmen; Aussage der Notiz): Q =
    Y + T3, schwacher Strom an der vierten Trägerkoordinate:
    signiert −1/24 und −1/36, positives Quadrat +1/48 und +1/72;
    Kreuzterm −1/32. Der signierte Wert ist ein Viertel der Differenz
    zweier positiver Auslesen.
    Beleg: Abschnitt G3.

23. **`one-charge-function-for-flat-sector-and-quarter-law`** (Runde
    28, Hypothese im Rahmen freier Fermionen auf zwei Kopien):
    F_α(Q) = Q²/2 + (α − 1/2)Q mit α aus der Markenart. PS: 1/4 und
    3/4, eindeutiger Leerzustand, R066. P und S: 0 und 1/2, je
    Fermion ein Dublett, für acht Fermionen 256 gleich gewichtete
    Grundzustände bei jeder Temperatur, gerade Hälfte 128. Nativ ist
    die Dichte ohne Potential flach, in allen 27 Konfigurationen
    nicht. Nicht gezeigt: dass P1 beide Markenarten enthält.
    Beleg: Abschnitte G5 und G6, Suchspur `twist_rule`.

25. **`recorded-half-phase-is-the-sign-of-the-chart-fermions`** (Runde
    29, belegt im genannten Rahmen): g8⁴ = (−1)^b (aufgezeichnet)
    ist gleich exp(2πi e_k·p) für jedes k und alle 240 Wurzeln: das
    Vorzeichen der sechzehn Chartfermionen um ein Bein. Materiebeine
    tragen es, der Mittelstrom nicht. Im E8-Netz selbst sieht kein
    Feld das Vorzeichen.
    Beleg: `check_sign_marks_and_complex_seam.py`, C1 und C2.

26. **`native-real-structure-maps-quarter-to-three-quarter`** (Runde
    29, belegt im genannten Rahmen): Die aufgezeichnete Realstruktur
    K antikommutiert mit der Gradierung; Konjugation mit K bildet
    den Operator mit Verschiebung 1/4 in allen 27 nativen
    Konfigurationen exakt auf den mit Verschiebung 3/4 ab.
    Beleg: Abschnitt C4.

27. **`two-copies-are-one-charged-seam-field`** (Runde 29, Hypothese
    H-C): ψ₊ = ψ1 + iψ2; Vierteldrehung = Vierteldrehung der Ladung;
    S und P = Ladungskonjugationen; die Reinheit einer reellen Kopie
    ist |⟨exp(i(π/2)Q_A)⟩|. S, P und die Vierteldrehung wirken auf
    (1, 2, −1, −2) wie die Symmetrien des Quadrats auf die vier
    Marken; der aufgezeichnete Spiegel x → i·conj(x) entspricht dem
    Kopientausch. Viertelwinkel und Markengewicht 1/2 geben sechzehn
    komplexe Fermionen, acht je Kopie. Nicht gezeigt: dass das native
    Feld mit 96 Fasern dieses Feld ist.
    Beleg: Abschnitte C3, C5 bis C7, Suchspur `charged_seam`.

28. **`native-frame-field-and-conjugate`** (Runde 30, belegt im
    genannten Rahmen): K·conj(K) = 1; die Gradierung antikommutiert
    mit K; je 48 Fasern, K bildet die eine Menge auf die andere ab;
    jeder Fourieranteil des nativen Potentials antikommutiert mit der
    Gradierung. Nativ ist die Rate zwischen den Kopien q = 1 für jede
    Funktion von h_N.
    Beleg: `check_native_quarter_count.py`, V1 und V5.

29. **`quarter-count-on-the-native-operator`** (Runde 30, bedingt):
    ⟨exp(i(π/2)Q_A)⟩ im gefüllten negativen Spektrum, alle Ortsmengen
    und alle 27 Konfigurationen: verbundene Zählung Nachbarn 12,20
    bis 12,24, gegenüber 1,47 bis 1,49 (log Betrag), Phasen innerhalb
    ±0,0027π; ohne Potential je Faser 1/√2, 3/2, 1. Kein Ankerwert.
    Beleg: V2 bis V4.

30. **`native-place-for-the-recorded-binary-cuts`** (Runde 30,
    Hypothese H-D): U = Gradierung, Y = Richtung des konstanten
    nativen Potentials, Kappe = reelle Linie mit −iYγ = η, gerader
    Transfer: Z_U = Z_Y = 0, Z_UY = iηZ_0, positiver Lift 1/2, q = 1,
    für [1,1,1] und [2,2,2]. Nicht gezeigt: dass P1 seine Marken so
    liest; der ursprüngliche Higgsmarker ist das nichtnegative
    Quadrat des Potentials und keiner der beiden Schnitte.
    Beleg: V6.

31. **`real-structure-forces-imaginary-charge-marks`** (Runde 31,
    belegt im genannten Rahmen): Hyperladung, Isospin und elektrische
    Ladung sind ungerade unter der Realstruktur, das Higgsquadrat
    (Summe der Quadrate der vier Higgsrichtungen) ist nichtnegativ,
    gerade und neutral. Auf jeder reellen Linie mit geradem Transfer
    ist die Marke des Skalars reell, jede Ladungsmarke und die
    gemischte Marke rein imaginär (Form i·η·q). Elektrische Ladung:
    Einzelmarke null, gemischte Marke nicht null. Flache Dichte:
    gemischte Hyperladungsmarke rein imaginär, etwa −1,4·10⁻⁶ i bis
    −2,4·10⁻⁶ i je nach Schritt; elektrisch null. Der Transfer ist
    gewählt; die ursprüngliche Zeitkarte ist unbekannt.
    Beleg: `check_native_charge_marks.py`.

32. **`native-operator-has-real-eigenlines-only-at-shift-zero`**
    (Runde 32, belegt im genannten Rahmen): Bei der aufgezeichneten
    Viertelverschiebung ist in allen 27 Konfigurationen keine
    Eigenlinie von h_N reell (kleinster Defekt 0,015). Ein reeller
    Modenoperator existiert nur bei Verschiebung null. Dort hat das
    konstante Potential 12 eindeutige reelle Eigenlinien, alle
    Majorana-Neutrinolinien (Werte ±0,0852 bis ±5,568). Auf ihnen
    sind die Marken (1, 0, ±1, 0), nicht (1, 0, 0, iη). Π = −iYγ
    antikommutiert mit dem Potential; η bleibt offen.
    Beleg: `check_native_real_lines.py`.

33. **Stand des Eintrags** (Runde 33): Die Sitzung ist als
    `problem:four-marks-quartic-and-free-seam-source` mit 28 Befunden
    im lokalen Atlas gespeichert, Revision
    `6ee27b08f424835e3d49af39` (Runde 33). Punkt 30 ist dort korrigiert (das
    native Potential ist −iγD_F; „Richtung des Potentials“ in Punkt 30
    war sign(D_F)). Neu eingetragen: die Kappe als Masseneigenlinie
    von D_F mit η = Vorzeichen des Werts, und der Ausschluss der
    Determinante als Ursprung der negativen gemischten Antwort.
    Beleg: `check_native_dirac_cap.py`.

34. **Stand nach Runde 34**: Der Eintrag hat jetzt 29 Befunde,
    Revision `c744bae80cc823dfe71a52cc`. Neu: Der Kandidat H-D
    liefert auf allen 27 nativen Konfigurationen vier Marken,
    gemischte Marke i·η·q mit q von 0,98842 bis 1 (bedingt, kein
    Originalproduzent). Beleg: `check_native_candidate_marks.py`.

35. **Stand nach Runde 35**: 30 Befunde, Revision
    `e3129e661d37ac382784aa8f`. Neu: Die Linie des Kandidaten ist
    das tiefste Niveau von |h_N|; bei Schritt null sind die Marken in
    allen 27 Konfigurationen exakt (1, 0, 0, iη); der Schritt folgt
    nicht aus dem Bestand. Beleg: `check_candidate_line_and_step.py`.

36. **Stand nach Runde 36**: 33 Befunde, Revision
    `e0c40348789cdfcd9d6f1691`. Neu: eine Kommutatorregel für die
    Gruppen A1 und A9 (belegt), die modulare Uhr zweier Bögen
    (belegt als Geometrie, nicht als Zeit der Quelle), der Ausschluss
    einer einzigen abelschen Vierteldrehung für Vorzeichen und
    Gewicht 1/2. Ein Satz aus Runde 29 ist zurückgezogen.
    Beleg: `check_commutator_rule_and_note16.py`.

Zu „Ausgeschlossen im genannten Rahmen“ kommt hinzu:

24. **`quarter-energies-are-no-rate-between-the-copies`** (Runde 28):
    Ein Transfer mit den Energien 1/4 und 3/4 vertauscht mit der
    Vierteldrehung, nicht mit dem Vorzeichen; er gibt
    Z_U/Z_0 = tanh(3s/4) statt null. Die aufgezeichneten Marken gelten
    nur bei s = 0.
    Beleg: Abschnitt G1.


20. **`exchange-and-sign-on-equal-copies-is-not-group-a1`** (Runde
    26): Mit U = S und Y = P auf zwei gleichen Kopien ist Z_U = Z_Y
    und die verbundene gemischte Antwort nie negativ; aufgezeichnet
    sind −1/32 und −1/72. Zwei unverbundene Kopien entgegengesetzter
    Parität geben ABAB = −1, aber Z_U = Z_UY = 0.
    Beleg: `check_marks_u_y_dictionary.py`, Abschnitt Q2.

## Nicht zum Eintragen vorgeschlagen

- Die Regel (3,1,1) der Suchspur mit dem Wert 7/4 am Quadrat: reine
  Zahlenübereinstimmung mit dem Referenzverhältnis aus Runde 17.

- Aussagen der eingefügten Notizen, die nicht nachgerechnet sind
  (Stabilität der Zerlegung, Herleitung der Determinante aus
  Kreismoden, Gegenmodelle mit verschobener Strommode).
- Reine Zahlenübereinstimmungen: Anker (1, 1, 2) gegen die Seiten des
  Dreiecks, harmonische Diagonale gegen die Gewichte am Quadrat.
- Die drei Lesarten, die in der Sitzung zurückgezogen wurden (falscher
  Torus in Runde 12 und 13, Viertelbögen, eine Messung entscheidet
  das Maß).

## Offen und unverändert

Q₀ mit Leerreferenz; ein Produzent für die ursprüngliche Naht; ob die
physische Naht A2 erfüllt; ob die Kernkonstante von P1 die
Viertelvarianz ist. Der native Aufruf liefert weiterhin keinen der
vier R113-Werte (in dieser Sitzung selbst ausgeführt).
