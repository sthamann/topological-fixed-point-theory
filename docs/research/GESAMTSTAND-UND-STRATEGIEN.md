# TFPT: Gesamtstand und Strategien zur Gesamtlösung (Stand H6h50, 25. November 2026)

Neu: Abschnitt 14 hält die Folgerunde N1–N8 und die Analysen 13–17 fest; die tiefere Welle (T1–T6, U1–U3) läuft. Abschnitt 13: G1–G3 und Analyse 12; Abschnitt 12: F1–F6 mit den Analysen 10 und 11.

Dieses Dokument fasst den ganzen Stand zusammen: gelöste Bausteine, Ursprung, Rekonstruktion und offene Tore. Daraus
entwickelt es neue Strategien. Maßgeblich bleiben Ledger und Originalquellen. Nichts hier ist befördert; Explorationsbefunde
stehen im Atlas als bedingte Befunde.

Statuskürzel:
- **E:** exakt unter genannten Voraussetzungen;
- **C:** bedingt;
- **L:** Lesart;
- **O:** offen;
- **N:** numerisch gegen Daten.

Quellen: [`TFPT-STAND.md`](TFPT-STAND.md), [`VOGELPERSPEKTIVE.md`](VOGELPERSPEKTIVE.md),
`experiments/exploration/2026-11-23-phasenraum-f5-holocode/VERDICT.md` (H6h1–H6h41), Ledger
`verification/status_ledger.csv` (TFPT.TOE.COMPLETE.01 mit den Toren T1–T8), Atlas-Problem
`problem:c3-p1-p2-shared-origin-reconciliation`.

## 1. Der Stand in einem Absatz

- **Zahlen und Teilcheninhalt:** Die liefert TFPT bereits aus zwei Eingaben: eine Generation mit 16 Zuständen, E₈, die
  Clocks, φ₀ und die Massenformel, α, die Mischungswinkel und Ω_b.
- **Eingaben:** P2 folgt aus P1. Der Kern der Suche ist also der Ursprung von P1, der Naht.
- **Ursprung (H6):** Er leitet 13 der 20 Prüfkriterien her; 7 stehen unter Lesart.
- **Rekonstruktion seit dem 24. November:** Aus dem Original allein, der quasifreien Quarterkovarianz bei β = 4π, ist die
  ganze freie und kinematische Schicht gewonnen:
  - Zustand und Uhr;
  - eine kontinuierliche positive Zeit (im Vakuumsektor);
  - die 3+1-Kinematik mit Lorentzgruppe und Orten;
  - masselose Felder.
- **Offen** ist alles Dynamische und der Ursprung selbst:
  - die wechselwirkende lokale 4D-Theorie (Yang–Mills, Kontinuum);
  - die nichtlineare Gravitation mit G;
  - der Grund, warum die Naht existiert;
  - ein Rest des Flavours.
- **Ledger:** Alle acht TOE-Tore (T1–T8) sind offen.

## 2. Gelöste Bausteine (TFPT-STAND §3, werden nicht neu geöffnet)

| Baustein | Ergebnis | Status |
|---|---|---|
| Teilchen | Λ^gerade ℂ⁵ = 1 + 10 + 5 = 16, eine Generation mit ν_R; Gruppe (SU(3) × SU(2) × U(1))/ℤ₆ | E |
| Vier Marken | Divisor vom Grad 4: drei Familienrichtungen (H¹ vom Rang 3) und fünf Plätze (h⁰ = 5) | E / C |
| E₈ | D₅ + A₃ über ℤ₄ verklebt: 240 Wurzeln, Dimension 248 | E |
| Clocks | Coxeterzahl 30, μ₄, Netz C₅/C₆/C₃₀, E₈ ↔ G₃₁ | E |
| Seed und Massen | φ₀ = 1/(6π) + 3/(256π⁴); Massenformel mit L, K für alle neun geladenen Fermionen; m̂_u/m̂_d = 55/117 | E (Transfer zu Polmassen C) |
| α | α⁻¹ = 137,03599921684 (Wurzel global eindeutig), gemessen 137,035999177(21) | E (Lesart C) |
| Mischung, Kosmologie | sin²θ₁₂ = 0,30675, CKM, Ω_b = 0,04894 | N, C |
| Prozessstruktur | GNS aus vollständigen Antworten; 68 = 5 + 63 | E |
| P2 aus P1 | fünf Plätze mit 3 + 2 aus der Viermarkennaht (P2.PARTITION.01) | E / C |

## 3. Der Ursprung von P1 (Kandidat H6): Prüfstand 13 P, 7 L

**Hergeleitet (P):**
- K1 vier Marken mit zwei Blättern; K2 τ = i; K3 μ₄ permutiert die Marken; K4 Spinstruktur (A,A);
- K5 c₃ = 1/(8π) ohne Wahl; K6 Zustand und Zeit; K7 flacher Sektor für E₈; K8 P2 über P1;
- K9 drei Familien als gleiche Kopien; K10 Familiencharaktere; K11 E₈ aus dem Prinzip;
- K14 keine freien Parameter; K17 Rahmungsanomalie.

Alles gilt unter den Explorationsprämissen A\*, C′, D, E, F.

**Unter Lesart (L), die eigentliche Front:**

| Kriterium | Was fehlt |
|---|---|
| K12 3+1 | Kinematik nativ (z79), aber Gruppe über alle Moden, ein Feld statt Spezies und der Bezug zur markierten Naht fehlen |
| K13 Gravitation | linear beide Helizitäten; nichtlinear, universelle Kopplung, G und Flächennormierung offen |
| K15 Eichdynamik | Yang–Mills in 4D nicht hergeleitet (nur Erhaltungsgruppe und Ward-Identitäten) |
| K16 Kontinuum | 4D-Regulierung und Kontinuumslimes offen |
| K18 Neutrinomischung | TM1-Spalte in den Clocks, Brechungsgewicht offen |
| K19 Existenz der Naht | warum es die Naht gibt; der Ursprung im engeren Sinn |
| K20 gemeinsamer Elternoperator | ein Operator für Naht, E₈ und Familien |

## 4. Die Rekonstruktion seit dem 24. November (z58–z81)

Leitlinie: Rekonstruktion statt Auswahl. Dynamik kommt aus den vollständigen Antworten des Originals, Chiralität darf aus
Indizes kommen.

| Schritt | Ergebnis | Status | Rechnung |
|---|---|---|---|
| Generator | Modulgenerator 4πx, CCR, gaußsche Ströme, Leiter statt Band, Chiralität aus Index | E | z58 |
| 2D-Kanäle | Tensorkanal c = 1 je Quarter; Spektraldimension der Quelle 1, der E₈-Naht 8 kompakt (keine native 4D-Antwort) | E | z59, z60 |
| Twistorgerade | π_*(η ⊗ π*O(2)) = O(1) ⊕ O(1) mit kanonischem η | E | z61, z64, z70 |
| Dublettzustand | ρ_m einzig je Windung m; m = 1 per Uhrenquotient | E / C | z63, z74 |
| Zeit | globaler Modulfluss periodisch, keine Zeit; lokale Modulzeit dicht; parabolische Wiesbrock-Zeit im Vakuumsektor | E (1D) | z70–z72 |
| Feldkern | det(I − ZW†)⁻¹ (ν = 1) durch Hardy- und Stromauslesung erzwungen | E (unter bi-U(2)) | z74 |
| 4D-Inhalt | Σm²qᵐ = konformer Skalar auf ℝ × S³ gezählt; die Gruppe ist Lesart | E / L | z73 |
| Kanal und Sphäre | zustandserhaltender Zweipfadkanal (Schwinger-Gewichte), Petz gibt ℓ(ℓ+1); thermisch exakt | E / C | z76, z80 |
| Vakuum und Orte | Araki-Woods-Vakuum = Purifikation; vertauschende P_μ, P² = 0; Orte aus gemischten Korrelationen | E | z78 |
| Quarterindex | Twistorraum mit Tomita-Partner, Lorentz = Kommutant von τ, Orte als Twistortranslationen | E | z79 |
| Markensymmetrie | auf dem Nahtdublett Ordnung 4, auf η und dem Quarterpaar 2, auf den Oszillatoren 1; auf H⁰(N) euklidisch | E | z79 |
| Uhr | R = i^(Q_a − Q_b), R² = (−1)^F: ein Spin-ℤ₄-Erzeuger | E | z81 |
| Saatfamilie | Einteilchendaten = {ρ_m auf Fermipaaren (x, −x), m ungerade} | E | z81 |
| Intervall | s² = 4(1 − Re 1/⟨X\|Y⟩) = det(Y − X) | E | z80 |
| Gravitation | linear beide Helizitäten (googly); kein Graviton aus T(z); zusammengesetzter Spin 2 ohne Pol | E / negativ | z65, z75, z80 |
| Skala | keine Skala aus dem eigenen Spektrum; ⟨T⟩ = 0 am KMS-Punkt | E | z68 |

**Korrektur im Verlauf:** z77 (H6h37) hatte im thermischen Kanal die Übergänge vertauscht. z80 korrigiert das; richtig
zugeordnet ist der Kanal bei jeder Temperatur exakt zustandserhaltend.

## 5. Was genau fehlt (Tore T1–T8 und Kriterien)

| Tor (Ledger) | Inhalt | Beitrag bisher | Offen |
|---|---|---|---|
| T1 Strukturpostulat → P1/P2 | Ursprung der Eingaben | P2 aus P1; H6: K1–K11 | Existenz der Naht (K19), Elternoperator (K20) |
| T2 Naht = (E₈)₁ | Naht als physischer Randzustand | Twistorgerade, Spinbrücke, (E₈)₁-Teilbeweise | voller Beweis im Ledger |
| T3 lokale unitäre 3+1D-Hamiltonfamilie | wechselwirkende lokale Theorie | freie Kinematik, positive Zeit 1D, Kern ν = 1 | Wechselwirkung, Lokalität, ein Feld statt Spezies |
| T4 chirales SM mit lokalem Maß | Chiralität ohne Spiegel | Index (z58), Spin-ℤ₄ mit R² = (−1)^F | 4D-Regulierung, lokales Maß |
| T5 Lorentz, Confinement, Clustering, Streuung | Infrarotverhalten | Lorentzkinematik (z79) | Dynamik |
| T6 Kopplungen und Neutrinotextur intern | Zahlen | α, Mischung, Massen (gelöst) | Neutrinozweig, Phasen, absolute Quarkamplitude |
| T7 masseloses Spin 2 mit universeller Kopplung | Gravitation | linear beide Helizitäten | nichtlinear, Universalität, G |
| T8 Anfangszustand und SK-Funktional | Kosmologie | KMS-Zustand des Originals | eindeutiger Anfangszustand, ein Funktional |

Ledgergrenzen, die jeder Weg respektieren muss:
- GRAV.PARENT.GAP.NO_GO.01: Ein gleichmäßig gegappter Elternoperator schließt einen masselosen Spin-2-Pol aus. T7 braucht
  eine kritische Folge mit lokalen Tensor-Eichbedingungen.
- GRAV.NONCIRCULAR.01: kein Schluss „modular ⇒ Bisognano–Wichmann“.
- DIMENSION.UPLIFT.FIREWALL.01: 3+1-Aussagen bleiben typisiert.

## 6. Drei Lehren der letzten Runden

1. **Zwei Faktoren.** Das Original zerfällt exakt in Ladungsteil ⊗ Oszillatorteil (Bosonisierung, z58).
   - Uhr und Marken leben im Ladungsteil. R = i^(Q_a − Q_b), R² = (−1)^F; die Marken sind eine ℤ₄-Faser über einem
     Punkt.
   - Die Raumzeitkinematik lebt im Oszillatorteil. Dort ist R trivial, und die Himmelssphäre ist rund.
   - Die Frage „I_pos für die markierte Naht“ verlangt also, Inneres und Äußeres gleichzusetzen. Wahrscheinlich ist sie
     falsch gestellt.
2. **Kinematik ja, Dynamik nein.** Alles Hergeleitete ist frei (gaußsch) und kinematisch. Nichtgaußsch und damit
   dynamisch sind im Original nur:
   - die nichtabelschen Stromantworten (z75);
   - der Ladungssektor (stark nichtgaußsch, z58 Teil 4).

   Dort muss die Wechselwirkung herkommen.
3. **Geschlossene Wege** (gescopt):
   - Graviton aus T(z) und aus zusammengesetzten freien Feldern;
   - native Spektraldimension 4 der Quelle;
   - globaler Modulfluss als Zeit;
   - T′ als markenverträgliche Identifikation;
   - gaußsche Vervollständigungen.

   Neue Wege dürfen sie nicht wiederholen.

## 7. Neue Strategien

Jede Strategie nennt Idee, Anlass, erste Rechnungen, Erfolg, Abbruch und die Tore, die sie bewegt.

### S1 Zwei-Faktoren-Architektur: Innen = Ladungsgitter, Außen = Oszillatoren

- **Idee:**
  - Innen ist alles, was R sieht: Marken, Uhr, E₈-Gitter, die 16, Flavour. Das sind genau die gelösten Bausteine.
  - Außen ist der Oszillatorteil. Er trägt Twistorraum, Lorentzgruppe, Orte und Felder.
  - Die Kopplung sind die Vertexoperatoren, also die Fermionen selbst (Ladung ⊗ Oszillator).
  - Universell ist sie, weil jedes Feld dieselbe Oszillatorhülle und denselben Zustand teilt.
- **Anlass:**
  - R ist auf der Raumzeit unsichtbar (z79).
  - R² = (−1)^F macht R zu einem Spin-ℤ₄-Erzeuger, derselben Art wie die ℤ₄ der Kette μ₄ → 16 (z81).
- **Erste Rechnungen:**
  - (a) Die Ladungen der 16 unter R: Ist R die anomalierelevante ℤ₄ (ℤ₁₆)?
  - (b) Die 4D-Lesart eines Vertexoperators. Das fermionische Quarterpaar ist ein Spinor mit g² = −1 (z79). Trägt das
    Fermion damit Helizität ±1/2 im Twistorbild?
  - (c) Die Kopplung aller Ladungssektoren an denselben Oszillator-Spannungstensor.
- **Erfolg:**
  - Materie = Ladungsgitter ⊗ Spinor;
  - Raumzeit = Oszillatoren;
  - Wechselwirkung = Vertexoperatoren mit universeller Hülle;
  - T′ und I_pos sind für die markierte Naht nicht mehr nötig.
- **Abbruch:**
  - R ist nicht die ℤ₄ der 16;
  - oder Vertexoperatoren tragen keine konsistente 4D-Helizität.
- **Bewegt:** K12, K19, K20; T1, T4.

### S2 Radiale Quantisierung: Moden sind Energiestufen, nicht Spezies

- **Idee:**
  - Die Mode m ist die m-te Energiestufe eines einzigen konformen 4D-Felds auf ℝ × S³. Die Stufe hat den K-Typ (j, j) mit
    j = (m−1)/2 und m² Zuständen.
  - L₀ ist der konforme Hamiltonoperator, elliptisch wie im Original.
  - Die Minkowskizeit ist der parabolische Erzeuger (z72), mit Vakuum und P₀ ≥ 0 (Lüscher–Mack).
  - Das Speziesproblem (eine Kopie je m, z75) verschwindet.
- **Anlass:**
  - Σm²qᵐ ist schon der Einteilchencharakter des konformen Skalars (z73).
  - Der Zweipfadkanal hat zwei native Elementarschritte (Teilchen hoch, Loch runter) mit Schwinger-CCR-Gewichten (z76,
    z80).
  - Mit dem Tomita-Partner ergeben sie 2 × 2 Heber, genau die vier P_αα̇ im kompakten Bild. Die Petz-Rückkanäle wären die
    Senker K_μ.
- **Erste Rechnungen:**
  - (a) P_αα̇ = (Schritt α) ⊗ (Partnerschritt α̇) auf ⊕_m H_m ⊗ H̄_m bauen.
  - (b) Die Kommutatoren [K, P], [P, P], [K, K] prüfen: Schließt su(2,2) mit L₀ als konformem Hamiltonoperator?
  - (c) Den Spinorindex zuordnen: Schrittart oder Quarter?
  - (d) Das thermische Original als Zustand des Felds auf S³ bei β = 4π lesen.
- **Erfolg:**
  - su(2,2) schließt mit nativen Gewichten;
  - ein Feld;
  - Vakuum und P₀ ≥ 0 hergeleitet;
  - K12 wird von Lesart zu Herleitung (freier Teil von T3, T5).
- **Abbruch:**
  - Die Kommutatoren schließen nur mit nicht nativen Gewichten (die Toeplitz-Hürde aus z73 bleibt);
  - oder der Kanal trägt die Dynamik nicht (Mehrpaaranteil (m−2)/(3m−2), z80).
- **Bewegt:** K12, K16; T3, T5 (frei).

### S3 Euklidisch zuerst: native Reflexion und Osterwalder–Schrader (Fortsetzung von WOIT.OS.TWISTOR.01)

- **Idee:**
  - Die Markensymmetrie verträgt sich auf den Nahtlagen nur mit der euklidischen Realform (z79). Das Original ist KMS,
    also von Natur aus euklidisch in der Zeit.
  - Die Lorentzform soll daher aus einer reflexionspositiven euklidischen Struktur rekonstruiert werden, statt als T′
    gesetzt zu werden.
  - T′ wäre dann ein Satz: Die Reflexion bildet den Nahtfaktor V* auf den η-Faktor ab.
- **Anlass:**
  - Im Ledger ist WOIT.OS.TWISTOR.01 offen; α- und β-Meilensteine sind frei gerechnet (v519, v522, v524, v565).
  - Neu ist der Hinweis, welche Faktoren die Reflexion vertauschen muss.
- **Erste Rechnungen:**
  - (a) Kandidaten für θ: Tomita-J, Antipodenabbildung, Teilchen-Loch-Spiegel, radiale Inversion Z → (Z†)⁻¹ auf der
    Matrixkugel.
  - (b) Reflexionspositivität der nativen Schwingerfunktionen (KMS-Korrelationen und ν = 1-Kern) für jedes θ.
  - (c) OS-Rekonstruktion: Ist der rekonstruierte Hamiltonoperator die parabolische Zeit und nicht der periodische
    Modulfluss?
- **Erfolg:**
  - eine native, reflexionspositive θ;
  - T′ als Folge;
  - positive Zeit aus Reflexionspositivität.
- **Abbruch:**
  - keine native θ ist reflexionspositiv;
  - oder die Rekonstruktion gibt nur den Modulfluss zurück.
- **Bewegt:** K12, K6; T3, T5. Mit S2 kombinierbar: Bei konformen Feldern ist die radiale Inversion die natürliche θ.

### S4 Gravitation als Deformation des Außenraums

- **Idee:**
  - Ist die Raumzeit der Twistorraum der Oszillatoren, dann ist nichtlineare Gravitation eine Deformation seiner
    komplexen Struktur. Das ist Penroses nichtlineares Graviton; googly für die zweite Helizität, ambitwistoriell für
    volle Gravitation.
  - Die native Quelle sind die einzigen nichtgaußschen Antworten: die nichtabelschen Stromantworten (z75) und der
    Ladungssektor.
- **Erste Rechnungen:**
  - (a) Die Deformationsklasse erster Ordnung (Kodaira–Spencer, H¹(PT, O(2))) aus der verbundenen Dreipunktantwort
    ⟨J^aJ^bJ^c⟩ = A f^abc und aus den Kumulanten κ₃, κ₄ des Ladungssektors.
  - (b) Ist sie ein selbstduales Graviton, und koppelt sie universell an alle Vertexoperatoren?
  - (c) Die kubische Ordnung mit dem Einstein-Vertex vergleichen.
  - (d) Die Ledgergrenze GRAV.PARENT.GAP.NO_GO.01 prüfen (kritische Folge, lokale Tensor-Eichbedingung).
- **Erfolg:** ein natives Graviton mit universeller Kopplung und Einstein-Vertex.
- **Abbruch:** Die Klasse ist reine Eichung, hat die falsche Helizität oder koppelt nicht universell.
- **Bewegt:** K13; T7.
- **Gegenhypothese S4b (Marken = Weylkrümmung):**
  - Die Marken sind die Hauptnullrichtungen eines selbstdualen Weylspinors (z61). Dann wäre die markierte Naht eine
    Twistorgerade einer gekrümmten, D₄-symmetrischen selbstdualen Einstein-Metrik (Hitchin/Tod, Painlevé VI mit vier
    Punkten).
  - **Entscheidungsexperiment gegen S1:** S1 verlangt, dass R innen wirkt (Spin-ℤ₄, Ladungen der 16). S4b verlangt, dass
    R als Raumzeitsymmetrie der Lösung wirkt.
  - z79 spricht für S1 (R trivial auf den Oszillatoren, euklidisch auf H⁰(N)). Entschieden ist es nicht.

### S5 Universalität und Flächennormierung aus dem gemeinsamen Zustand

- **Idee:**
  - Alle Sektoren teilen denselben KMS-Zustand und dieselbe Uhr. Das ist ein Kandidat für die universelle Kopplung, also
    das Äquivalenzprinzip.
  - Die Flächennormierung soll aus der Verschränkung der Purifikation Ω kommen. Die ist exakt berechenbar, je Mode endlich
    wegen q = e^(−2πm).
- **Erste Rechnungen:**
  - (a) Die Verschränkungsentropie von Ω zwischen Original und Partner, je Mode und summiert.
  - (b) Ihr Verhältnis zur TFPT-Einsteinkopplung 1/c₃ (GRAV.NONLINEAR.01) und zur Flächenregel t = 1 (z57).
  - (c) Universalität: Gleiche Uhr für Ladung und Oszillatoren heißt gleiche Kopplung an den Zeitfluss.
- **Erfolg:** eine native Entropiedichte, die 1/c₃ und die Flächennormierung trifft, ohne Bisognano–Wichmann.
- **Abbruch:**
  - Das Verhältnis ist nicht fest;
  - oder der Schnitt Original|Partner braucht die verbotene Horizontlesart (GRAV.NONCIRCULAR.01).
- **Bewegt:** K13 (Flächennormierung); T7.

### S6 Ursprung der Saat: vom Qubit zum Original

- **Idee:**
  - Die Einteilchendaten des Originals sind genau die Saatfamilie {ρ_m, m ungerade} auf den Fermipaaren (z81). Das
    Original ist ihre zweite Quantisierung.
  - Die Ursprungsfrage wird dadurch klein. Warum die Saat ρ₁, und warum alle ungeraden m?
  - Ein Kandidat für ρ₁: zwei Ja/Nein-Operationen mit Vorzeichensymmetrie des Gesetzes erzwingen M₂(ℂ) und τ = i (z80,
    bedingt). Der Zustand muss die Symmetrie brechen, sonst gibt es keine Zeit; der einzige treue KMS-Zustand je Windung
    ist ρ_m (z63).
  - „Ungerade m“ wäre die Spin-ℤ₄-Bedingung R² = −1 auf Fermionen.
- **Erste Rechnungen:**
  - (a) Folgt die ganze Familie aus ρ₁ durch eine native Regel, etwa Sym^m oder die Uhr, und die Paritätsbedingung?
  - (b) Ist die zweite Quantisierung (CAR über den Fermipaaren) die einzige Vervollständigung mit Uhr und KMS?
  - (c) Was fixiert die Lücke π der Saat (β = 4π)?
- **Erfolg:** Original = Saat ρ₁ + Spin-ℤ₄ + zweite Quantisierung; damit K19 aus einem operativen Prinzip.
- **Abbruch:** Mehrere ungleiche Vervollständigungen, oder die Familie braucht freie Wahl.
- **Bewegt:** K19, K20; T1, T8.

### S7 Szenario-Suche im Suchprogramm

- **Idee:** Die Strategien sind Wahlen in einer ganzen Theorie. OriginLab soll ganze Szenarien als Kombinationen
  aufzählen:
  - Spinorindex: Quarter, Schrittart oder η;
  - Rolle der Marken: innen, euklidische Krümmung oder Himmelssphäre;
  - Zeit: Wiesbrock, OS oder radial;
  - Moden: Spezies oder Stufen;
  - Gravitation: Deformation oder Verschränkung;
  - Ursprung: Saat oder zwei Fragen.
- **Prüfung:** Jede Wahl wird an den berechneten Klauseln aller Befunde geprüft: R-Äquivarianz, Positivität (Wallach),
  Reflexionspositivität, Schluss von su(2,2), ℤ₁₆-Anomalie, Firewalls. Erst danach kommen Grenzen.
- **Ausgabe:** Welche Szenarien überleben, welche Regeln leer sind, und wo die nächste Rechnung am meisten entscheidet.
- **Bewegt:** alle, als Steuerung.

## 8. Reihenfolge

1. **Phase 1, parallel (Architektur entscheiden):**
   - S1 (a, b): Ist R die ℤ₄ der 16, und tragen Fermionen Helizität 1/2?
   - S2 (a–c): Schließt su(2,2) aus Kanal und Partner?
   - S6 (a): Folgt die Saatfamilie aus ρ₁?
   - S7 bekommt die Universen.
2. **Phase 2:**
   - S3, mit der in Phase 1 bestätigten Struktur;
   - S4 auf dem Außenraum;
   - S5 parallel.
3. **Phase 3:** die wechselwirkende Theorie über Vertexoperatoren und nichtabelsche Ströme (K15, T3–T5), dann der
   Restflavor (K18, T6) aus den Clocks, ohne gelöste Bausteine neu zu öffnen.

**Entscheidungspunkte:**
- Scheitert S2, bleiben die Moden Spezies. Dann braucht es die zweite Quantisierung über die Moden als Regel
  (Analyse 8: Kompositionsregel für Ereignisse).
- Scheitert S1 (R ist nicht die ℤ₄ der 16), sind die Marken nicht innen. Dann wird S4b zum Hauptweg.

## 9. Was nicht getan wird

- Keine gelösten Bausteine neu öffnen (Flavour, Massen, α, 16, E₈).
- R113 ist eine Route, kein Tor.
- Keine Auswahlprinzipien, kein „modular ⇒ Bisognano–Wichmann“.
- 3+1 bleibt typisiert.
- Keine Ledgeränderung ohne Zustimmung.
- Eine vollständige Theorie ist nicht erreicht. Dieses Dokument ist ein Plan, kein Ergebnis.

## 10. Stand nach der ersten parallelen Runde (H6h42)

| Strategie | Rechnung | Ergebnis | Status |
|---|---|---|---|
| S1 Zwei Faktoren | z82 (24/24) | R ist auf dem Träger die ℤ₄ der 16 (Zentrum von Spin(10), bis auf Eichung); Spin und Statistik innen exakt; das Original erreicht nur die „2“; mit den Quarter-Oszillatoren als Lorentzgruppe haben Fermionen kein 2π-Vorzeichen | E, teilweise |
| S2 radiale Quantisierung | z83 (20/20) | su(2,2) schließt mit den Tomita-Wurzeln des nativen Hebers L^nat; ν = 1-Leiter; jede Helizität einmal; Minkowskienergie [0, ∞); Moden = Energiestufen | C (eine Identifikation) |
| S3 euklidisch, OS | z84 (21/21) | native reflexionspositive Spiegelungen; für das thermische Original gibt OS den Modulfluss zurück, positive Zeit nur im Vakuumsektor; T′ nur punktweise, R ohne Fixpunkt auf U(2) | E / C (H6h44) |
| S4 Deformation | z85 (23/23) | native Quellen geben keine universelle Grad-2-Klasse, zweite Ordnung Heisenberg statt Einstein; S4b ausgeschlossen, Marken = Hitchins vier Punkte (Symmetriedatum) | E negativ (H6h43) |
| S5 Verschränkung | z86 (29/29) | keine Flächennormierung aus Original\|Partner (Typ I, kein Horizont); eine Uhr für alle Sektoren | E negativ |
| S6 Saat | z87 (26/26) | Saat ρ₁ + Spin-ℤ₄ + CAR; Ursprung auf P_sym und P_clock eingegrenzt; Gaußprämisse aus P_clock in allen Abschneidungen | E / C |
| S7 Szenarien | OriginLab `scenarios` | von 960 Szenarien bleiben 8, alle mit Schrittart-Spinor, inneren Marken, parabolischer Zeit, Energiestufen und Saat ρ₁ | Steuerung |
| Analyse 9 | z88 (4/4) | voller 4D-Fockraum passt energietreu nicht in eine Naht; lineare Leiter aus der Hardy-Invarianz | E |

### Was sich am Gesamtbild ändert

- **Die Architektur ist eingegrenzt.**
  - Raumzeit: die radiale ν = 1-Leiter mit dem Spinorindex „Art des Elementarschritts“ (Teilchen hoch, Loch runter). Die
    Uhr R wirkt darauf trivial.
  - Innen: die Marken und die Uhr R, ein Spin-ℤ₄ mit R² = (−1)^F, auf dem Träger die ℤ₄ der 16.
  - Die Quarter-Oszillatoren tragen die Kinematik (z79), aber nicht die Statistik der Materie (z82).
- **Ein Teilchen in der Naht, viele Teilchen durch Zusammensetzung.**
  - Die Einteilchenleiter steckt in einer Naht (z83).
  - Der volle 4D-Fockraum passt energietreu nicht hinein (z88).
  - Viele Teilchen müssen also aus vielen Nähten oder Ereignissen entstehen. Das ist die Kompositionsregel, nach der
    Analyse 8 fragte.
- **Gravitation braucht einen echten räumlichen Schnitt oder eine Deformation.** Der thermische Schnitt Original|Partner
  scheidet aus (z86). Offen bleiben die Deformation (S4) und eine räumliche Kappe von S³ im rekonstruierten Feld.
- **Der Ursprung ist auf zwei Prämissen geschrumpft.** P_sym und P_clock (z87); die Leiter folgt aus der Hardy-Invarianz
  (z88).

### Nächste Schritte (aus den Rechnungen)

1. **S2:**
   - Die Identifikation „verdoppelter Diagonalschritt = L^nat“ aus einem Prinzip herleiten, etwa als einzige vollständig
     positive, J-symmetrische Fortsetzung in der Standardform.
   - Dann die Zusammensetzung vieler Nähte als zweite Quantisierung prüfen.
2. **S1:** die Verklebung Spin ×_{ℤ₂} ℤ₄. Liefern die Fermipaare einer Stufe K-Typen (j, j ± ½) mit −1 ↦ R²? Und wo liegt
   die „3“ (Paaroperatoren mit k ≡ 2)?
3. **S5 nach S2:** die Vakuumverschränkung einer räumlichen Kappe von S³ im rekonstruierten Feld, mit den nativen
   Gewichten qᵐ als Regulator. Ist S/A regulatorunabhängig?
4. **S6:**
   - die Starrheit (Gaußzustand aus P_clock) für beliebig viele Fermipaare beweisen;
   - prüfen, ob P_clock aus P_sym über den Tausch A ↔ B folgt (τ → −1/τ, Zeitneutralität).
5. **S3 und S4:** eingebaut (H6h43, H6h44).

## 11. Bilanz der Runde und die entscheidenden nächsten Rechnungen (H6h44)

**Was die Runde festgelegt hat.** Die Szenario-Suche lässt von 1152 ganzen Szenarien 12 übrig, alle mit derselben
Architektur:
- **Raumzeit:** die radiale ν = 1-Leiter aus dem nativen Heber L^nat. Der Spinorindex ist die Art des Elementarschritts.
- **Zeit:** parabolisch, im Vakuumsektor. Wiesbrock, radiale Quantisierung und OS-Halbraum geben dieselbe Zeit.
- **Innen:** Marken und Uhr R, ein Spin-ℤ₄. Die Marken sind ein Symmetriedatum (Hitchins vier Punkte), weder Krümmung noch
  Himmelssphäre.
- **Ursprung:** Saat ρ₁ + Spin-ℤ₄ + CAR mit den Prämissen P_sym und P_clock.
- **Viele Teilchen:** eine Zusammensetzung vieler Nähte. Eine Naht trägt die Einteilchenleiter.

Offen ist allein die Gravitation, mit genau zwei Wegen.

**Die entscheidenden nächsten Rechnungen**, nach Hebelwirkung:
1. **Gravitation über den radialen Kanal.**
   - Gibt die native Komposition zweier Helizität-(−2)-Sektoren in einen dritten (D = −4 ⊗ D = −4 → D = −4) auf der radialen
     Leiter die Plebański-Penrose-Klammer {h_a, h_b} aus z85?
   - Wenn ja, ist das der kubische Einstein-Vertex, nativ.
2. **Gravitation über einen räumlichen Schnitt.**
   - Die Vakuumverschränkung einer Kappe von S³ im rekonstruierten Feld, mit den nativen Gewichten qᵐ als Regulator.
   - Ist S/A regulatorunabhängig? Das wäre der Weg zu G, ohne Bisognano–Wichmann vorauszusetzen.
3. **Identifikationsprinzip und Zusammensetzung.**
   - „Verdoppelter Diagonalschritt = L^nat“ als einzige vollständig positive, J-symmetrische Fortsetzung herleiten.
   - Dann die zweite Quantisierung als Zusammensetzung vieler Nähte prüfen.
4. **Materiespin.** Die Verklebung Spin ×_{ℤ₂} ℤ₄: Tragen die Fermipaare einer Stufe K-Typen (j, j ± ½) mit −1 ↦ R²? Und
   wo liegt die „3“?
5. **P1 für die Kugel D(V, η).** Den Wallach-Parameter ν von D(V, η) aus den Antworten des Ladungssektors lesen.
6. **Starrheit.**
   - Den Gaußzustand aus P_clock für beliebig viele Fermipaare beweisen.
   - Prüfen, ob P_clock aus P_sym folgt (Tausch A ↔ B, τ → −1/τ).

## 12. Folgerunde (H6h45)

| Folgeschritt | Rechnung | Ergebnis | Status |
|---|---|---|---|
| F1 Gravitation über den Infinity Twistor | z89 (17/17) | I_π bis auf die Modulphase nativ, Einstein-Klammer exakt; native Kompositionen erzeugen sie nicht; radialer Kanal nur AdS; G frei | C / E negativ |
| F2 Räumlicher Schnitt | z90 (17/17) | native Gewichte regulieren nicht; Flächenkoeffizient lokalisierungsabhängig; induziertes 1/G = 0 bei ξ = 1/6 | E negativ |
| F3 Identifikation und Komposition | z91 (24/24) | eine Prämisse P_koh; freier 4D-Fockraum aus bosonischer Komposition (P_ident); keine native Fusion gibt Wechselwirkung | C / E negativ |
| F4 Materiespin und die „3“ | z92 (14/14) | „3“ nur algebraisch (schwaches Triplett); keine native Verklebung; Statistik im Q²/2 der radialen Energie | E negativ |
| F5 Kugel D(V, η) | z93 (11/11) | kein Kugelkern aus dem Ladungssektor; P1 ist eine Wahl; D(V, η) innen | E negativ |
| F6 Starrheit, P_clock aus P_sym | z94 (19/19) | Starrheit des Originals trägt die Nichtentartung (S(2) bewiesen, 3 und 4 Paare in Fenstern); P_clock ⇔ P_tn (Zeitneutralität des Gesetzes), unabhängig von P_sym | E / O |
| Analyse 11 | z95 (11/11) | Rechnungen bestätigt, auch die Schwelle N = 6 für identische Saaten (Korrektur H6h46) | E |

### Was die Folgerunde bedeutet

- **Kinematik ist abgeschlossen, Dynamik fehlt.**
  - Die freie 4D-Theorie steht, unter P_koh und P_ident: radiale Leiter, Vielteilchenraum aus bosonischer Komposition,
    parabolische Zeit, Infinity Twistor bis auf die Phase.
  - Jede Wechselwirkung braucht etwas, das keine native gerade oder lineare Struktur liefert: Vertexoperatoren,
    nichtabelsche Ströme oder eine Komposition mit fester Kopplung.
- **Gravitation ist auf einen Weg eingegrenzt.**
  - Die Einstein-Klammer ist da, sobald I und eine Kopplung g gegeben sind.
  - Nativ fehlen die Phase von I und g. Der radiale Kanal gibt AdS mit dem Radius der Vakuumgröße.
  - Weder ein räumlicher Schnitt noch die Verschränkung Original|Partner gibt G.
- **Materie bleibt innen.**
  - Die „3“ kommt aus dem gelösten Träger.
  - Die Statistik steckt im Q²/2 des Ladungsgitters.
  - Eine Verklebung mit dem 4D-Spin braucht eine neue Regel vom GSO-Typ.

### Nächste Schritte (aus den Rechnungen)

1. **Phase von I und Kopplung g:**
   - Fixiert die native OS-Spiegelung Θ₁ = τ (z84) die Modulphase von I?
   - Ist der P₀-graduierte kubische Vertex dann reflexionspositiv, und fixiert die Positivität g?
2. **P_koh herleiten:** Hat die natürliche unitäre Dilatation des Schritts L^nat im gemeinsamen CAR genau das
   zweidimensionale Hardy-Halbregister als Umgebung?
3. **Wechselwirkung über Vertexoperatoren:**
   - die Ladungs-Vertexfusion zweier Nähte über die gemeinsame Oszillatorhülle;
   - verschiebt sie Δ des zusammengesetzten 4D-Felds von 1 weg?
4. **GSO-artige Regel:** (−1)^{2h} = R² auf Ladung ⊗ Paarzuständen. Schließt sie unter su(2,2) und den Vertexoperatoren,
   und folgt sie aus J bei β = 4π oder aus dem ℤ₄-Kleber von E₈?
5. **Ursprung schließen:**
   - S(n) für verschiedene m mit der Nichtentartung beweisen (z94);
   - P_tn aus ⟨T⟩ = 0 (A*, E₆(τ) = 0) herleiten. Dann bliebe eine einzige Prämisse.

## 13. Folgeschritte G1–G3 und zwölfte Analyse (H6h47)

| Schritt | Rechnung | Ergebnis | Status |
|---|---|---|---|
| Analyse 12 (MacDowell–Mansouri, Plebański, Clifford, Λ–G) | z96 (18/18) | alle Rechnungen bestätigt; Clifford-Paar = Instanton/Anti-Instanton der Twistorfaserung (Ladungen durch Whitney-Summe erzwungen); Λ–G-Relation = COSMO.LAM.01 in neuen Variablen; Vorzeichen verlangt dS, nativ ist AdS | E / O |
| G1 Spin-Statistik-Verklebung | z97 (11/11) | (−1)^{2h} = R² schließt unter su(2,2), Fusion und gekleideten Vertizes (Schrittart × Ladungsschritt); nicht nativ (Prämisse P_GSO); vektorartig; Chiralität nur mit P_A3 (E₈-Kleber, A3 = Raumzeitalgebra), im Konflikt | C / E negativ |
| G2 Phase, Kopplung, MacDowell–Mansouri | z98 (11/11) | Phase von I Eichrichtung, Positivität fixiert g nicht; MM-Algebra nativ (AdS, Radius halbe Vakuumgröße); Tomita-Familie = flaches MM-Vakuum; Σ± aus Original und Partner; Plebański-Sektor nur bis Z₄ | E / C / E negativ |
| G3 Wechselwirkung über Vertexoperatoren | z99 (19/19) | Vertexroute = su(2)₁-Strom; Δ = 1 starr; marginale Deformation ohne lange Vielfachheiten; verbundene 1→2-Amplitude erst mit gewähltem Zweikörperterm | E negativ / O |
| Szenario-Suche | OriginLab (13 824) | 12 Szenarien erfüllen alles außer Graviton und Chiralität; minimale leere Regeln [graviton], [admissible, chiral_matter] | E |

### Was die Runde bedeutet

- **Die Kinematik der Gravitation ist nativ, ihre Dynamik nicht.**
  - Nativ sind:
    - die MacDowell–Mansouri-Algebra mit dem AdS-Vakuum als flacher Cartan-Verbindung;
    - τ als Punktspiegelung;
    - Σ± aus Original und Partner;
    - das flache Vierbein;
    - I bis auf die Phase.
  - Nicht nativ sind:
    - die Phase von I und die Kopplung g (also G bzw. b);
    - der Plebański-Sektor, also Gravitation oder topologisch;
    - das dS-Vorzeichen;
    - jede Krümmung F ≠ 0.
- **Spin-Statistik gibt es mit einer Prämisse, Chiralität nicht.**
  - Die Verklebung ist die Paarung „Ladungsschritt × Schrittart“.
  - Die Regel braucht P_GSO und gibt Dirac-Paare.
  - Chirale Materie verlangt die Z₄-Fassung des E₈-Klebers. Diese kollidiert in der Lesart P_A3 mit der Skalarleiter und
    mit dem gelösten Familiengitter.
- **Wechselwirkung gibt es nativ nicht.**
  - Vertexoperatoren und Ströme sind auf Niveau 1 dieselbe Struktur.
  - Alle nativen Operatoren sind einkörperig, und Δ = 1 bleibt starr.
- **Gemeinsamer Engpass.** Gravitationsdynamik (F ≠ 0 mit festem b) und Wechselwirkung (verbundene 1→2-Amplitude mit
  festem g) verlangen dasselbe Fehlende: eine nicht gaußsche Zweikörperstruktur mit festem Koeffizienten. Das passt zum
  offenen Kriterium K20, dem gemeinsamen Elternoperator.

### Nächste Schritte (aus den Rechnungen)

1. **Gemeinsame Quelle für Krümmung und Wechselwirkung:** Gesucht ist ein nativ nicht gaußscher Zustand, z. B.:
   - das Kondensat D = −4 (z89);
   - der kollektive Zustand identischer Saaten ab N = 6 (z94);
   - ein Kompositzustand (z91).

   Zu prüfen ist an seiner Tomita-Familie:
   - Ist F ≠ 0, mit festem quadratischem Koeffizienten (b, G)?
   - Gibt derselbe Zustand verbundene 1→2-Amplituden mit festem g?
2. **Zweikörperterm:** Liefert der Doppeltransfer e^{2i(φ₁−φ₂)} zusammen mit dem Tomita-Partner die fehlenden
   □-Nachkommen mit unitären Vorzeichen?
3. **Chiralität:** Gesucht ist eine Z₄-Fassung ohne Skalarkonflikt: die Skalarleiter im Sektor R = −1, gekleidet mit dem
   tiefsten R = −1-Boson aus z92. Sonst wird T4 als gescoptes Negativ festgehalten.
4. **Kleinere Prüfungen:**
   - gradierte Nahtkomposition aus einer gemeinsamen CAR;
   - native Orientierung für den Gravitationssektor;
   - Q₀ als modulare Zeit einer Kappe (dS-Vorzeichen).

## 14. Folgerunde N1–N8 und Analysen 13–17 (H6h48)

| Schritt | Rechnung | Ergebnis | Status |
|---|---|---|---|
| N1 eine Quelle für Krümmung und Wechselwirkung | z100 (10/10) | Maurer-Cartan-Lemma: keine native Tomita-Bahn ist gekrümmt; kubischer Vertex braucht ungerade Ladung und Nichtgaußizität | E negativ |
| N2 kollektiver Zustand aus sechs Saaten | z101 (13/13) | vollständig klassifiziert (6·10¹³), ganzzahlige Mehrkörperterme, keine Dynamik; J zustandsunabhängig | E negativ |
| N3 Chiralität aus E₈ ⊕ E₈ | z102 (11/11) | Klebeweg geschlossen (Coset-Satz, familienblind) | E negativ |
| N4 Orientierung und Sektor | z103 (13/13) | keine native Orientierung wählt; minimale Prämisse P_vol (ein Bit) | E negativ |
| N5 Λ > 0 aus der Kappenzeit | z104 (13/13) | Vorzeichen von Λ = Klasse der Zeit; positive Zeit heißt Λ ≤ 0; dS nur mit P_cap | C |
| N6 Zweikörperterm und induziertes G | z105 (15/15) | keine □-Nachkommen, λ nicht fixiert, G regulatorgebunden | E negativ |
| N7 gemeinsame CAR | z106 (13/13) | gradierte Komposition unter P_CAR; Spin-Statistik mit P_GSO | E / C |
| N8 Ursprung | z107 (12/12) | P_tn folgt aus A*, nicht zirkulär; Ursprung auf P_sym; S(3) auf Abspaltung reduziert | E / O |
| Analyse 13 | z108 (11/11) | Kompositionssatz, kompensierte Rotation, Rekombinationsnorm 32Δ(Δ−1), Varianzschranke | E / C |
| Analysen 14 und 15 | z109 (9/9) | Starrheit aus einer Frequenz, inverse Rekonstruktion, Achter-Konstruktion, Wilson-Dirac-Grade, modulare Eindeutigkeit | E / C |
| Analysen 16 und 17 | z116 (10/10) | Cl(3,11)-Majorana-Weyl-16, Uhr = Zentrum von Spin(10), Spur mit Uhr = ε, Nahtindex wählt einen Sektor, Austausch log 3 | E algebraisch, C nativ |

### Was die Runde bedeutet

- **Ursprung:** Er ruht nur noch auf P_sym; P_tn und P_clock folgen aus A*.
- **Jede offene Eigenschaft hat jetzt eine präzise mathematische Form:**
  - Wechselwirkung: Rekombination □φ = λφ³ über drei Nähte;
  - Gravitationsdynamik: eine Nicht-Bahn-Familie und die absolute Normierung;
  - Λ-Vorzeichen: die Klasse der Zeit;
  - Sektor: ein Bit;
  - Chiralität: der Index der vollen Algebra.
- **Neuer Kandidat für drei Lücken zugleich:** die gemeinsame reelle Clifford-Struktur mit der Uhr als Zentrum von
  Spin(10). Sie würde Chiralität, Gravitationssektor und Spin-Verklebung liefern. Offen ist, ob die native Quelle sie erzeugt.
- **Zweite Hypothese:** die Naht als Twistorgerade mit Stromalgebra (Parke–Taylor, CSW, Hodges).

### Laufende tiefe Welle

- **T1–T6 (z110–z115):**
  - Parke–Taylor aus den Nahtströmen;
  - CSW-Komposition;
  - Hodges und Normierung von G;
  - P_GSO aus Modulinvarianz;
  - Chiralität als Twistorindex;
  - Krümmung bei bewegter Algebra.
- **U1–U3 (z117–z119):**
  - erzeugt Tomita die gemeinsame Clifford-Struktur?
  - D_⊥ und die Spur mit Uhr als Wirkung, Gewicht 2c₃;
  - Transferaustausch und Petz-Prinzip statt P_koh.

## 15. Tiefe Welle T1–T6, U1–U3 und achtzehnte Analyse (H6h49)

| Schritt | Rechnung | Ergebnis | Status |
|---|---|---|---|
| T1 Parke–Taylor aus Nahtströmen | z110 (16/16) | Nenner exakt, Koeffizient durch k = 1 fest; Zähler braucht P_N4; Mehrspurterme ungedämpft | C |
| T2 CSW-Komposition | z111 (18/18) | η nativ aus I_π; es fehlt genau ein Austauschfeld (kinetischer Term, P_BF) | C |
| T3 Hodges und G | z112 (12/12) | Hodges = KLT aus nativen Klammern; Blattdiagonale braucht P_diag; G an einer Länge L | C |
| T4 P_GSO aus Modulinvarianz | z113 (13/13) | nein; pro Naht mit Modulsummation unverträglich; vier chirale Nähte geben (E₈)₁ | E negativ |
| T5 Chiralität als Twistorindex | z114 (13/13) | nur P_RJ (mit Spiegel) ist J-kovariant; Index 3 je 16 | C |
| T6 Krümmung bei bewegter Algebra | z115 (9/9) | Familie flach (kompaktes Dual U(2)) | E negativ |
| U1 Tomita und Clifford | z117 (15/15) | J_H6 = Majorana-Weyl auf dem Indexschnitt; Chiralität mit einem Bit P_glue | C |
| U2 Uhrspur als Wirkung | z118 (11/11) | MacDowell–Mansouri unter P_R; Gewichtswahl = P_vol; G nur mit 2c₃ von Hand | C / E negativ |
| U3 Austausch und Petz | z119 (13/13) | Austausch = modulares Gewicht ohne Vertex; P_koh = P_Torus + P_Petz | C |
| Analyse 18 | z120 (7/7) | Energie-Eindeutigkeit; lokales Gatter P_T − i·P_S mit erzwungenen Winkeln; ein Weylpunkt; Tr A² | E / C |

### Was die Runde bedeutet

- **Vorwärtsregel:** Das lokale Transfergatter hat keine freien Winkel. Leckfreiheit erzwingt die Vierteldrehung, also
  P_sym; die Phase ist die Wurzel der Uhr, gleichwertig die ℤ₄-Wechselphase i^(Q_A·Q_B). Ursprungsprämisse und
  Wechselwirkungsregel sind dieselbe Aussage. Offen ist, ob die spinblinde Ladungsenergie die native Nahtantwort ist.
- **Hypothese T** trägt Nenner, Farbe, Doppelkopie und Chiralität als Index. Es fehlen kinetischer Term, Zähler und eine
  Trennung der Mehrspurterme.
- **Chiralität** hat zwei bedingte Wege mit je einem Bit (P_RJ, P_glue); beide sind die ℤ₄-Regel auf Spinoren.
- **Statistik** ist der neue Engpass: P_GSO ist nicht modular nativ. Die Szenario-Suche (46 080 Szenarien) lässt 108
  Szenarien, die alles außer Graviton und modular konsistenter Statistik erfüllen.
- **Gravitation:** Kinematik, Klammern und die MacDowell–Mansouri-Spur sind nativ bis auf je ein Bit. Dynamik braucht eine
  inhomogene Familie, G eine Länge.

### Laufende Welle V1–V5 (z121–z125)

- V1: Gatter → Kette → Kontinuum, Dreinahtvertex; ist der Kreuzterm Q_A·Q_B nativ?
- V2: P_RJ aus den gekleideten Vertizes (z97).
- V3a: kinetischer Term als Tomita-Paarung zwischen D = ∓2.
- V3b: Zähler ⟨ij⟩⁴ (Gradtwist, Pati–Salam-Quartett) und Mehrspurterme.
- V5: P_sym aus kohärenter Komposition; S(n) über die Energie.

## 16. Wellen V–Z und G, neunzehnte Analyse (H6h50)

| Schritt | Rechnung | Ergebnis | Status |
|---|---|---|---|
| V1 Kette, Kontinuum, drei Nähte | z121 (14/14) | nativ trivial; unter P_null XXX-Familie, alternierende Kette integrabel, su(2)₁-Lesart | C |
| V2 P_RJ aus den Vertizes | z122 (12/12) | P_IF erzwingt P_RJ mit festem Vorzeichen; P_IF selbst Prämisse | C |
| V3a kinetischer Term | z123 (19/19) | Tomita-Paarung gibt CSW-Propagator, CSW = BCFW; bleibt P_X | C |
| V3b Zähler, Mehrspurterme | z124 (16/16) | kein nativer Gradtwist; Doppelspuranteil = konforme Supergravitation | E negativ |
| V5 P_sym aus Komposition | z125 (19/19) | ja unter P_univ | C |
| W2 Prämissen der Uhrspur | z127 (15/15) | P_R größtenteils Folge; Paarungsbit offen; Λ·G = −3/(1024πκ) | E / C |
| W3 Krümmung aus Zustandsprofilen | z128 (11/11) | eine Stabilisatorbahn, flach | E negativ |
| X1 Index gegen Chiralitätsbit | z129 (11/11) | Bit P- und C-ungerade, Index invariant; dasselbe Bit wie das Vorzeichen der Wirkung | E negativ |
| X2 P_univ aus der Standardform | z130 (17/17) | Niveauabgleich, nativ unter P_TP | C |
| Z1 P_null in zwei Faktoren | z133 (13/13) | P_null ⇔ P_sl2; P_viertel und P_ord bleiben | C |
| Z4 Yang–Mills aus Schleifen | z135 (15/15) | Plakettenstruktur und sin²θ_W = 3/8 nativ; Kopplung und Dynamik nicht | E / C |
| G1 Verschränkungsgleichgewicht | z136 (19/19) | G/R_vac² freie reine Zahl; fünf Zeilen, keine zwei gleich | E / N |
| G2 induzierte Gravitation | z137 (19/19) | Vorzeichen AdS, G·Λ = −3π/(4cL); legt P_G nicht fest | E / N |
| Analyse 19 | z144 (67/67) | P_uhr trägt P_null und P_viertel; Zopfrelation nur einkörperig; zwei isomorphe Verklebungen | E / C / L |

### Nachtrag z145: Netzprinzip für P_uhr

Geprüft wurde, ob P_uhr daraus folgt, dass Nähte Sektoren eines Netzes sind und das Gatter ihre Statistik ist. Streng
gilt das nicht: Im nativen Netz (Uhr geeicht, Sektorgruppe ℤ₁₆) ist die Statistik der Nähte ein Produkt von Uhren, das
Gatter also trivial. Das wechselwirkende Gatter verlangt eine transferierte Einheit mit Gewicht ungerade/8. Es ist die
gesteuerte Uhr zwischen den Eigenmoden der Vierteldrehung P_sym und gleichwertig eine Paarphase. Im E₈-Bild ist es die
4π-Drehung des D₅-Faktors gegen A₃, ohne freien Winkel; das ist eine Lesart. P_uhr bleibt Prämisse, Dynamik folgt nicht.

### Nachtrag z146–z153: Welle A, zwei neue Strategien, zwanzigste Analyse

- **G:** Drei weitere Wege (Relativentropie-Wirkung z146, ρ⋆ z147, gemeinsame Schleifenantwort z150) legen G nicht
  fest. Der Faktor e^(−2/α) ist im Bestand eine Leiterformel, keine Dynamik. Einziger verbleibender Hebel: ein
  nichtstörungsartiges Gewicht e^(−S) mit S = 2/α.
- **Graviton (z148):** auf eine Prämisse reduziert, P_h2 oder P_X + P_pos. Test: D(3,2) im Dubleton-Sektor.
- **Gatter:** E₈-Lesart (z149), 4D-Chern–Simons (z151) und der E₈-Code der zwanzigsten Analyse (z153) sparen keine
  Prämisse. Die Gatterlinie ist ausgeschöpft.
- **Gitterlesart (z152):** Die Naht ist die Hälfte genau einer ladungserhaltenden Nächste-Nachbar-Kette (Hüpfen t,
  Wechselpotential 2t). Ordnung, Partner, Boost und Hamiltonoperator folgen unter P_CTM. Frei, 1+1-dimensional.
- **Strategie danach:** (a) Gitterlesart in vier Dimensionen (Kähler–Dirac-Zelle, 16 Komponenten, ℤ₄-Anomalie);
  (b) geeichte Uhr auf der Kette; (c) D(3,2) im Dubleton-Sektor; (d) elliptische 4D-Chern–Simons-Klasse auf der
  Markenkurve w² = 1 − z⁴. Zurückgestellt: weitere G-Zeilen aus Normierungen, weitere E₈-Gatter.

### Was die Runde bedeutet

- **Der Kandidat** braucht neun Prämissen: A*, I_step, P_lift, P_IF, P_univ, P_uhr, P_ord, P_curv und P_G als freie
  Zahl. P_null und P_viertel sind zu P_uhr geworden („die Kompositionsphase ist eine Funktion der Uhrwerte“).
- **Wechselwirkung:** Das Gatter ist unter P_uhr fest, die alternierende Kette integrabel. Die Wechselphase ist der
  Bicharakter der Uhrgruppe ℤ₄, derselbe wie bei D₅. Eichstruktur und das Verhältnis 3/8 sind nativ, Kopplung und
  Linkdynamik fehlen.
- **Zopf- und Defektlesart (Analyse 19):** Sie trägt nur einkörperige Operatoren und eine endliche Gruppe, keine Dynamik.
  P_IF ist die Wahl zwischen zwei isomorphen Verklebungen und nicht von innen entscheidbar.
- **Gravitation:** Vier Wege geben ein Vakuum von Planckgröße mit negativem Λ. G/R_vac² ist eine freie reine Zahl; der
  Ledgerwert verlangt 9,6·10¹²² Nähte je Horizont. Zustandsprofile tragen keine Krümmung.
- **Abgebrochen und offen:** modulare Statistik (z126), I_step (z131), Uhrhebung (z132), konforme Gravitation (z134),
  Spin-2-Pol (z138), Ordnung (z139), Maßbilanz (z140), Relativentropie-Wirkung (z141), ρ⋆ (z142).
