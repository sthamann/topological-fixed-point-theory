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

### Nachtrag z154–z157: Welle B, die vier Folgerechnungen

- **Vier Dimensionen (z154):** gescopt negativ. Die 16 ist keine Kähler–Dirac-Zelle, ein hyperkubischer Halbraum
  trägt die native Leiter nicht. Verträglich ist nur die Lesart mit neuer Prämisse P_KD (Lorentz mit schwacher Spin(4)
  verdrillt); auch dann folgt die Vier nicht nativ.
- **Geeichte Uhr (z155):** P_null ist ein Gauß-Gesetz (unter P_eich + P_schnitt). Wechselwirkung und Winkel folgen
  nicht. Als Hamiltonoperator zerstört die Eichung die Naht; die Naht liegt bei g² = 0.
- **Spin 2 (z156):** D(3,2) ist im Dubleton vorhanden, unitär, mit positiver Energie; P_pos ist Folge. Offen: Auswahl
  der Kopie und Austausch zwischen Nähten (P_X in enger Form). Der Dubleton ist exakt (Di + Rac)².
- **Elliptische Klasse (z157):** Die Nahtkette liegt wörtlich auf der Markenkurve. Der Viertelwinkel ist dort eine
  Invariante der komplexen Multiplikation (δ_y = −i für jeden Kreuzungsparameter), unter P_aut. Das z120-Gatter ist das
  rationale Ende der Familie, die Naht das freie.
- **Gemeinsamer Befund:** Der Nahtzustand ist frei, das Gatter ist ein eigener Kompositionsschritt mit eigener Zahl.
  P_uhr folgt nicht aus der Kette selbst. Der Viertelwinkel hat in beiden Rechnungen dieselbe Quelle, die
  Ordnung-4-Automorphie von τ = i.
- **Strategie danach:** (a) Zwei-Naht-Antwort im Kanal D = ±4 (entscheidet P_X für den Dubleton); (b) geeichte
  Schrittkette gegen das Destri–de-Vega-Gitter aus z121; (c) Kopplung zweier lemniskatischer Ketten am Uhrpunkt und
  der Limes λ → 0; (d) Indexoperator auf K gegen d − δ (P_KD); (e) Singleton-Lesart: zweite Quantisierung des
  Singletons und Höhere-Spin-Dualität als Test für P_X. Zurückgestellt: weitere Gitter für K12, weitere G-Zeilen.

### Nachtrag z158–z161: Welle C

- **Spin-2-Austausch (z158):** Der freie Transfer überträgt den Sektor, aber als Stromüberlagerung, nicht als
  Teilchen. „Faktorisiert ⇒ P_X widerlegt“ (z156) gilt nicht. P_X (eng) wird zu P_2Q + P_sing; die Zahl der Arten
  wäre P_N, und der Ledgerwert verlangte N = 1,08·10¹²³.
- **Geeichte Schrittkette (z159):** exakt das Destri–de-Vega-Gitter aus z121. Die native Naht ist für jeden Winkel
  stationär, weil ihr modularer Erzeuger ein Cartan-Element der Flavour-Symmetrie ist. Der Winkel folgt nicht.
- **Kopplung zweier Ketten (z160):** ionisches Hubbard-Modell, nicht integrabel. In der elliptischen Familie ist die
  Masse zwischen den rationalen Enden 8e^(−πK/λ) [L].
- **Indexoperator (z161):** P_KD ist geschlossen; vier native Daten sprechen gegen die Verdrillung.
- **Gemeinsamer Befund:** wieder keine Prämisse weniger. Nach zwölf Wegen in drei Wellen ist das kein Zufall mehr.

### Nachtrag z162: Diagnose und Neuausrichtung

**Warum die bisherige Linie feststeckt.**

- **D1, Gaußsche Daten [E]:** Zustand, modulare Erzeuger, Uhr und so(3,2) sind bilinear. Bilineare Entwicklung
  erzeugt keine Kumulante. Wechselwirkung beginnt am Gatter, und der native Zustand ist gegen dessen Winkel blind.
  Aus der einzelnen Naht folgt Kinematik, keine Kopplung.
- **D2, selbstdualer Punkt [E/L]:** Sandwich (B⁺EB) und Durchgang (BEB) sind eine Familie, verschieden um
  B² = Tausch·R⁻¹. Der Gitterparameter ist Θ = tan(φ/2)/2 bzw. cot(φ/2)/2. P_viertel ist der selbstduale Punkt
  Θ = 1/2. Dort ist m·a = 0,82: eine Skala. G/R_vac² kann mit P_viertel nativ nur von der Ordnung eins sein.
- **D3, ungenutzte Zahlen:** α und die Leiter e^(−A) wurden als „gelöst, nicht neu öffnen“ behandelt und deshalb nie
  als Eingabe der Dynamik benutzt. Die offenen Zahlen des Kandidaten sind aber genau diese Größen.

**Neue Strategien, nach Rang.** Jede ändert etwas an der Fragestellung, nicht nur an der geprüften Prämisse.

1. **S-A, Skala aus Kopplung.**
   - Idee: Die Hierarchie ist die dynamisch erzeugte Vakuumskala der Kompositionskette am schwach gekoppelten Ende,
     m·a = 4e^(−(π/2)cot(ε/2)) [L]. Mit ε = πα ist das 4e^(−1/α)·1,006.
   - Anders als bisher: Die Kompositionsregel wird geändert (Durchgang statt Sandwich), und α geht als Eingabe ein.
   - Vorprüfung (z162): Form und Exponent passen unter P_durch, P_eps, P_see. Der Vorfaktor passt nicht (a²/G = 4·10²
     bis 1,6·10⁴ gegen nativ 0,2 bis 24), die zweite Ledgerskala auch nicht (Rang 9,957).
   - Entscheidende Rechnung: ε aus der geeichten Uhr (z155) mit g² = 4πα. Ist die Schrittdauer in Einheiten des
     Gauß-Gesetzes ein Viertel, ohne Wahl?
   - Erfolg: ε = g²/4 folgt. Dann wird P_G zu einem Vorfaktor, und P_viertel wird durch P_durch ersetzt.
   - Abbruch: Die Normierung ist frei. Dann ist P_eps eine Zahl auf eine Zahl gelegt (Empfindlichkeit 137) und die
     Strategie nur eine Umformulierung.
   - Bewegt: P_uhr, P_G, K15, K16.
2. **S-B, Zustand der Welt.**
   - Idee: Über dem Referenzzustand ist die Anregung der Durchgangskette ein Dirac-Teilchen mit m·a = ε; der native
     Zustand ist ein Magnongas der Dichte 3,5·10⁻⁶ je Naht. Die Hierarchie lebt dagegen im antiferromagnetischen See.
   - Entscheidende Rechnung: reduzierter Ein-Naht-Zustand des Sees gegen die nativen Niveaus π und 3π.
   - Erfolg: verträglich. Dann ist der See der Kandidat für den Tomita-Partner (P_see ersetzt P_IF-Teile).
   - Abbruch: unverträglich. Dann lebt die Welt am Referenzzustand, S-A liefert die Dirac-Kinematik (K16), aber keine
     Hierarchie.
   - Bewegt: K16, T3 (ein Feld statt Arten), P_G.
3. **S-C, laufende Kopplung als Ziel.**
   - Idee: Die gelöste α-Gleichung enthält den Einschleifen-Koeffizienten b₁ = 41/10. Jede Dynamik des Kandidaten muss
     ihn als Vakuumpolarisation liefern.
   - Entscheidende Rechnung: Ladungsantwort der geeichten Durchgangskette mit dem Inhalt 3·16 plus Higgs.
   - Abbruch: Am selbstdualen Punkt läuft nichts; wenn auch am schwachen Ende kein Logarithmus entsteht, ist die
     1+1-Kette der falsche Träger der Eichdynamik (K15 braucht dann vier Dimensionen zuerst).
   - Bewegt: K15, die Verbindung des Kandidaten zum α-Block.
4. **S-D, Zählweg.**
   - Idee: R_vac²/G = 32π³e^(2/α) ist die Dimension eines Raums von 405,4 Qubits. Als Anzahl von Nähten ist das eine
     freie Zahl (P_N), als Zustandszahl eines endlichen Objekts eine endliche Aufgabe.
   - Entscheidende Rechnung: Ist e^(2/α) im α-Block eine native Zustandssumme (die Randgleichung ist eine
     Stationaritätsbedingung von log Z)?
   - Abbruch: kein natives Objekt mit dieser Dimension. Gefahr der Zahlenspielerei ist hier am größten.
5. **S-E, No-Go festschreiben.** D1 und D2 als Regeln in OriginLab: Szenarien, die eine Kopplung aus der einzelnen
   Naht herleiten wollen, sind geschlossen; Szenarien mit P_viertel können G nicht erklären.
6. **S-F, vier Dimensionen aus dem Kollektiv.** Zurückgestellt, bis S-A und S-B entschieden sind: Die Durchgangskette
   ist ein 1+1-Lichtkegelgitter; 3+1 verlangte Durchgänge in drei unabhängigen Richtungen.

**Ehrliche Einordnung.** S-A ist kein Treffer. Sie stellt P_G von einer freien Zahl 10¹²³ auf die Frage „warum
ε = g²/4“ und einen Vorfaktor um. Ihr Wert liegt darin, dass sie die drei Diagnosen zugleich angeht und eine kleine
entscheidende Rechnung hat. Aufgegeben wird die Linie „Prämisse aus der einzelnen Naht herleiten“.

### Nachtrag z163: Gesamtdurchsicht des Bestands und fundamentale Lösungswege

**Anlass.** Die Strategien S-A bis S-F stammen nur aus den H6-Rechnungen. Sie haben den Bestand des Ledgers nicht
benutzt. Die Durchsicht von Ledger, Paper 1, Atlas und Exploration korrigiert vier Punkte:

- **b₁ = 41/10 ist hergeleitet, zweimal.** Als (4/5)·41/8 aus dem Transportbudget des Flavour-Blocks und aus den
  Hyperladungen der Träger, Σ(2/3)Y_f² + (1/3)Y_s² = 41/6, über den Satz β = a₄ (Paper 1, Abschnitt zur α-Schließung;
  v434). S-C nannte das ein offenes Ziel. Das war falsch. S-C ist gestrichen; der Wert bleibt als Prüfstein für jede
  Dynamik.
- **G ist im Ledger keine freie Zahl.** Lokal gilt G_ab + Λg_ab = (1/c₃)T_ab mit 8πG = 1/c₃ (GRAV.NONLINEAR.01,
  bedingt), Λ kommt aus der Wirkungsleiter. Die „freie Zahl G/R_vac²“ der H6-Rechnungen ist die Λ-Stufe dieser
  Leiter. Offen ist ihr Mechanismus (die Wirkungsgrammatik e^(−qα⁻¹) hat den Status einer Identität), nicht ihr Wert.
- **Das Ledger benennt die Lücke selbst** (TFPT.TOE.COMPLETE.01): „A single shared 3+1D parent satisfying T3–T8 is
  still missing.“ Gesucht ist ein gemeinsamer vierdimensionaler Elternteil, nicht die Herleitung einzelner Prämissen.
- **S-A, S-B und S-F leben auf der 1+1-Kette.** Sie fallen unter das Muster F2 unten und sind zurückgestellt.

**Was trägt (Muster der Erfolge).**

1. Jede gelöste Zahl ist ein Index, eine Anomalie oder ein Wärmekern-Koeffizient (16 = ν der Uhr, 8 = Rang,
   41/10 = a₄, 5/4 und 3/4 = quadratische Formen von D₅ und A₃), oder sie ist der stationäre Punkt einer Determinante
   an einem symmetrischen Punkt (α aus δ log det + Gegenterm + Nahtantwort = 0; τ = i; der D₄-feste polystabile Punkt).
2. Keine gelöste Zahl stammt aus einer Zeitentwicklung mit Hamiltonoperator.
3. Der Träger ist immer derselbe: die Naht als (E₈)₁ bei τ = i mit der Uhr ℤ₄, Niveau 1.

Lehre: Das fundamentale Objekt ist eine Determinante über einem Raum von Hintergründen. „Dynamik“ ist die Abhängigkeit
der Determinante vom Hintergrund.

**Was scheitert (Muster der Fehlschläge).**

- F1: Wechselwirkung aus der einzelnen Naht (zwölf Wege, D1).
- F2: Kompositionsketten in einer Dimension (sieben Wege lassen die Skala frei, D2).
- F3: vier Dimensionen aus Gitter oder Kähler–Dirac (z154, z161).
- F4: zusammengesetztes Graviton und Elternteil mit Lücke (GRAV.PARENT.GAP.NO_GO.01, GRAV.WEINBERG.WITTEN.01, z158).
- F5: Modulfluss als Zeit oder Boost (GRAV.NONCIRCULAR.01, Spannung W4).
- F6: Auswahl statt Rekonstruktion; entropische Auswahl für T8; autonome endliche Dynamik.
- F7: Twistor-String-Lesart der Naht (Mehrspurterme, konforme Supergravitation, fehlender kinetischer Term).

Lehre: Alle Fehlschläge behandeln die Naht als Baustein, der mit anderen Nähten zusammengesetzt wird. Alle Erfolge
behandeln sie als Funktional, das ausgewertet wird.

**Befund der Durchsicht: zwei Hälften, die sich nie getroffen haben.**

- **Ledger-Hälfte** (CELEST.SEAM.01, SEAM.EQUIV.TWISTOR.01, seit Juli 2026): holomorphe Theorie auf dem Twistorraum
  von ℂ²/ℤ₄, Eichgruppe E₈, Axion aus der Anomalie, Niveau 1 aus der Flussquantisierung. Geltung laut Vertrag: „the
  self-dual sector only, never full Einstein gravity“. Status: vorbereitet, nicht konstruiert.
- **H6-Hälfte** (z110–z112, z123): Die Naht auf einer Twistorgeraden liefert die Vertizes. Es fehlen der kinetische
  Teil (P_BF), der Zähler ⟨ij⟩⁴ (P_N4), und die Mehrspurterme stören.
- Jede Seite nennt als fehlend, was die andere hat. Seit H6h34 hat keine H6-Rechnung das Ledger-Programm benutzt.
- Die bekannte Struktur, die beide verbindet, ist die Twistorwirkung: holomorpher Teil plus Kopplung mal
  ∫d⁴x (log det der Naht auf der Geraden L_x) [L, Mason 2005, Boels–Mason–Skinner 2007].

**Erste Entscheidungsprüfung (z163, 22 Prüfungen, reproduziert).**

- **Anomalie [E; Bedeutung L]:** Sp_adj X⁴ = λ²(X·X)² gilt für A₁, A₂, D₄, E₆, E₇, E₈ mit λ² = h∨ + 6
  (8, 9, 12, 18, 24, 36), nicht für so(10), su(4), su(5), su(3)+su(2). Die selbstduale Hälfte verlangt im Bestand E₈
  selbst; λ² = 36 ist der Wert des Vertrags. Das ist ein vierdimensionaler Grund für E₈, unabhängig von der Naht.
- **Statistik [E; Folgerung L]:** Die Uhr-Klassen (60, 64, 60, 64) heben den quartischen Teil nur gemeinsam auf
  (4 − 8 + 12 − 8 = 0). Mit der 16 als Fermion im Bündel bleibt immer ein Rest; Spin(10) mit n Weyl-16 und m
  reellen 10 hat den Koeffizienten 4 + 2n + m. Die chirale Materie kann nicht aus dem Bündel kommen, sie muss von der
  Naht selbst kommen (ν(Uhr) = 16). Das ist eine harte Vorgabe für T4, kein Erfolg.
- **Vertex [E]:** Der zusammenhängende Teil der Stromkorrelatoren, W = log Z, ist exakt linear im Niveau und reine
  Einfachspur (su(2) bis 6 Punkte, su(3) bis 5 Punkte; nur f und k gehen ein, also auch für (E₈)₁). Der Mehrspurteil
  aus z110 gehört zu Z, nicht zu W. Der Zähler ⟨ij⟩^p folgt aus der Homogenität: p = −w_b = 4 für den Serre-Partner
  O(−4) aus z123. P_N4 entfällt; P_BF wird zum Feld b der Ledger-Hälfte.
- **Leiter [E/L, N negativ]:** Die Leiter hat zwei unabhängige Stufen; die Hubble-Stufe ist Friedmann. Die Λ-Stufe ist
  exakt ρ_Λ/M̄⁴ = q(1/α)/c₃² mit derselben Funktion q(α) = δ_top e^(−2α) wie in der α-Gleichung, also die Nahtantwort
  am dualen Punkt (e·g = 4π). Die elektroschwache Stufe ist keine Potenz davon (36,85 gegen 28,31). Ein Hinweis mit
  einer Wahl, kein Mechanismus.

**Der Kandidat für das eine Prinzip (Arbeitsname P_TW).** Ein Funktional Γ = holomorpher Teil auf dem Twistorraum
(durch Anomaliefreiheit festgelegt: E₈, Axion, Niveau 1) + Summe über die Geraden von log det der Naht im Hintergrund.
Kopplungen und Skalen sind stationäre Punkte desselben Funktionals, so wie α es im Ledger schon ist. Darin hätte jedes
offene Stück einen Ort:

- Vier Dimensionen: der Raum der Geraden (K12 wird wörtlich; der Ledger-Selektor DIMENSION.SELECTOR.4D.01 bleibt
  bedingt).
- Wechselwirkung: Abhängigkeit der Determinante vom Hintergrund. Die einzelne Naht bleibt Gaußsch, D1 ist erfüllt.
- Spin 2: Abhängigkeit derselben Determinante von der komplexen Struktur. Das ist die Forderung von
  GRAV.SPIN2.EMERGENCE.01; die universelle Kopplung wäre automatisch. Das Graviton ist ein elementares Feld der
  holomorphen Hälfte, die No-Gos F4 greifen nicht.
- Materie: die Naht (16 Majoranas), nicht das Bündel.
- Prämissen dieses Kandidaten: A* (die Naht), P_TW (Twistorwirkung als Rahmen), P_SD (die Ledger-Hälfte existiert als
  Quantentheorie; ihre eigenen Reste, der starre Rest „32 T₃“ und das Maß, bleiben), P_g (Koeffizient des Vertex),
  P_reell (Signatur und Positivität), P_inf (Einstein statt konformer Gravitation), P_mat (Materie von der Naht).
  Das sind nicht weniger Namen als zuvor. Der Unterschied: P_g, P_reell, P_inf und P_mat sind Rechenaufgaben mit
  endlichem erstem Schritt, und die neun Prämissen der Kettenlesart (I_step, P_lift, P_IF, P_univ, P_uhr, P_ord,
  P_curv, P_G) gehören zum 1+1-Bild und werden hier nicht gebraucht.

**Lösungswege, nach Rang.**

1. **W1, Kopplung aus Stationarität (T6).**
   - Idee: Die α-Gleichung ist die Stationarität einer Determinante. Der Vertex der Twistorwirkung ist dieselbe
     Determinante mit Hintergrund. Dann ist der Koeffizient des Vertex keine freie Kopplung.
   - Grundlage: α-Block (gelöst), GAUGE.DETLINE.FIXPOINT.01 (offener Vertrag, sucht genau dieses Funktional), WP3 des
     Celest-Vertrags (Axion-Koeffizient an c₃ ausgerichtet), z163 Teil C.
   - Erste Rechnung: der abelsche Sektor. Hintergrund in Hyperladungsrichtung auf einer Geraden; Vergleich der drei
     Terme von F_U(1) (Maxwell α³, Blatt −2c₃³α², Transportdeterminante mit 41) mit kinetischem Teil, Axionterm
     (λ² = 36) und log det.
   - Erfolg: Die drei Terme haben je ein Gegenstück. Dann ist P_g durch die α-Gleichung ersetzt, und
     GAUGE.DETLINE.FIXPOINT.01 hat sein Funktional.
   - Abbruch: kein Gegenstück für Blatt oder Transport. Dann gibt der Rahmen nur die Form, und der α-Block bleibt
     getrennt.
2. **W2, Gravitation aus derselben Determinante (T7).**
   - Idee: Der nicht selbstduale Gravitationsvertex ist die Determinante aus z112 (Hodges) mit dem nativen
     Unendlichkeitstwistor I_AdS = I_π − I_ω aus z123. Sein Koeffizient ist die reine Zahl Λ·G.
   - Grundlage: z112, z123, GRAV.SPIN2.EMERGENCE.01, Λ-Stufe der Leiter; Einstein-Reduktion der konformen Gravitation
     [L, Maldacena 2011, Adamo–Mason 2013].
   - Erste Rechnung: Vertex mit dem nativen I_AdS. Verschwindet der Teil, der konforme von Einstein-Gravitation
     trennt? Ist der Koeffizient Λ·G?
   - Erfolg: Einstein-Vertizes mit Koeffizient Λ·G. Dann entfallen P_X, P_h2 und P_2Q, die Klausel „graviton“ von K13
     ist gefüllt, und die freie Zahl ist die Λ-Stufe.
   - Abbruch: nur konforme Gravitation (der Geist der Weyl²-Kürzung im Ledger). Dann bleibt Gravitation selbstdual und
     der Rahmen liefert T7 nicht.
3. **W3, Signatur und Positivität (T3, T5).**
   - Idee: euklidisch zuerst. Reflexionspositivität des Vertexterms mit der nativen Spiegelung J = θ∘τ.
   - Grundlage: WOIT.OS.TWISTOR.01, z123, alte Strategie S3. Hier muss die Spannung W6 aufgelöst werden: Im Bestand ist
     die Nahtgerade nicht die Gerade eines reellen Punkts.
   - Erste Rechnung: endliche Matrixprüfung der Positivität des Kerns ⟨JJ⟩·⟨12⟩⁴ auf gespiegelten Geradenpaaren.
   - Abbruch: keine Positivität. Dann gibt der Rahmen Amplituden, aber keinen Hilbertraum.
4. **W4, Materie von der Naht (T4).**
   - Idee: Die 16 und die drei Familien sind Nullmoden der Naht im Orbifold-Hintergrund, keine Bündelkomponenten.
   - Grundlage: z163 Teil B, ν(Uhr) = 16, die drei Sphären von H₂(ALE) im Vertrag, Spannung W1 (Familiencharaktere
     {i, −1, −i} gegen {1, i, −i}).
   - Erste Rechnung: Charaktere der Nullmoden unter Uhr und Verklebung; gibt es den Twist χ_q nativ?
   - Abbruch: kein nativer Twist. Dann sind die Familien in diesem Rahmen nicht geometrisch.
5. **W5, Leiter als duale Antwort (Skalen).** Schwächster Weg, gegen Zahlenspielerei markiert. Zu prüfen ist, ob die
   α-Gleichung ein magnetisches Gegenstück hat, dessen Antwort q(1/α) ist, und ob n = 2 aus der Spin-ℤ₄-Struktur
   folgt. Die elektroschwache Stufe braucht einen anderen Mechanismus.
6. **W6, Ursprung (T1, T8).** Im Rahmen P_TW wird „warum diese Naht“ zu „warum dieser Twistorraum“: E₈ aus der
   quartischen Bedingung, Niveau 1 aus dem Fluss, ℤ₄ aus dem Orbifold. Erst sinnvoll, wenn W1 oder W2 trägt.

**Reihenfolge und Abbruchregel.** W1 und W2 zuerst, sie entscheiden, ob der Rahmen Zahlen liefert. Dann W3 und W4.
Scheitern W1 und W2 beide, ist die Twistorwirkung als Elternteil geschlossen, und es bleibt der Gitterweg
(SEAM.EQUIV.MMST.01). Nicht weiter verfolgt werden neue Kettenlesarten und jede Herleitung einer Prämisse aus der
einzelnen Naht.

**Ehrliche Einordnung.** z163 zeigt Verträglichkeit an vier Stellen, keinen Treffer aus Zahlen. Der Rahmen P_TW ist
eine Annahme. Die Ledger-Hälfte ist selbst nicht konstruiert. Kopplung, Signatur, Gravitationsteil und Materie sind
offen. Nichts ist befördert.

### Nachtrag z164: W1 und W2, erste Rechnungen

**Ergebnis in einem Satz.** Der Rahmen P_TW liefert Formen und Verhältnisse, aber keine der beiden reinen Zahlen
(Kopplung, Gravitationszahl). Keiner der beiden Wege ist entschieden; beide sind schärfer gestellt (z164, 14
Prüfungen, reproduziert).

**W1, Kopplung.**

- **Struktur [E]:** Ausintegrieren von b gibt 1/g² ∝ A²/V (A: Koeffizient des holomorphen Teils, V: Koeffizient des
  Vertex). Unter b → s·b ist nur diese Kombination invariant. Der Rahmen hat genau eine freie reine Zahl.
- **Zähler als Integral [E; Lesart L]:** ∫d⁸θ χ⁴(σ₁)χ⁴(σ₂) = ⟨12⟩⁴ mit Koeffizient eins. Liest man b als oberste
  Komponente eines Zusammenhangs, ist der Vertex genau ∫d⁴x d⁸θ log det mit dem Niveau k als Koeffizient.
  Verhältnisse von Kopplungen sind dann Niveauverhältnisse und fest; der Gesamtwert nicht.
- **P_norm ist widerlegt [N negativ]:** Der Bestand (z110) hatte unter P_norm g = k^(−1/2) = 1. Der Lauf des
  Standardmodells bis zur reduzierten Planckmasse gibt 1/α = 34,3 (U(1), GUT-Norm), 48,6 (SU(2)), 50,6 (SU(3)), also
  1/g² = 2,7 bis 4,0 statt 1.
- **Hinweis, schwach [N]:** SU(2) und SU(3) liegen dort bei 1/g² = 3,9 bis 4,0, also 1/α nahe 16π = 2/c₃ (−3,2 % und
  +0,7 %; mit Zweischleifen-Eichanteil −3,9 % und +1,6 %). Ein nativer Kandidat für die 4 ist die Quellladung
  |μ₄| = 4 aus CELEST.WP5E.M1.01. Gegenproben: Als exakte Gleichung ist das nicht haltbar. Die Kurven schneiden sich
  bei 1/α = 47 (10¹⁷ GeV), und die Vorhersagen wären α_s(M_Z) = 0,123 (+4 %) und sin²θ_W = 0,244 (+5,5 %). Eine ganze
  Zahl und eine Skala sind gewählt, U(1) passt nicht.
- **α-Gleichung [O]:** Für den kubischen Term und den Blatt-Term ist kein Gegenstück im Rahmen gerechnet. Das
  Erfolgskriterium von W1 ist nicht erfüllt, das Abbruchkriterium auch nicht.
- **Vorbestand:** Die einfache Übertragung der α-Grammatik auf SU(2) und SU(3) ist im Ledger schon ausgeschlossen
  (GAUGE.DETLINE.FIXPOINT.01, v996). Der Weg über Niveaus ist davon verschieden, löst die Aufgabe aber nicht.

**W2, Gravitation.**

- **Gewichte [E]:** Serre-Paarung O(w) mit O(−w−4); Spin 2 hat die Gewichte 2 und −6. Die reduzierte
  Hodges-Determinante ist unabhängig von der Wahl der gestrichenen Zeilen und Referenzspinoren und skaliert je Bein
  mit t⁻⁴ (n = 5, 6). Der Zähler ist ⟨ij⟩⁸, allgemein ⟨ij⟩^(4s). Eine eigene Annahme ist dafür nicht nötig.
- **Unendlichkeitstwistor [E]:** Auf L_x ist (I_π − I_ω) = (1 + det x)⟨12⟩. Die Metrik det(dx)/(1 + det x)² hat
  konstante Krümmung: euklidisch S⁴ mit R = 48, lorentzsch AdS₄ mit R = −48, also |Λ| = 12 in nativen Einheiten. Die
  Form ist nicht entartet (Pfaffsche −1). Das Abbruchkriterium „nur konforme Gravitation“ greift hier nicht.
- **Die Zahl fehlt [E]:** Flache Einstein-Gravitation ist die Kontraktion I_π − ε²I_ω mit |Λ| = 12ε². Die freie Zahl
  des Kandidaten ist genau ε². Mit der Leiter des Ledgers (Λ/M̄² = 10^−120,15) und nativen Längen wäre
  ε² = 10^−123,3 bis 10^−121,3.
- **Nicht nachgerechnet [L]:** Dass der Vertex der konformen Theorie auf Einstein-Zuständen den Einstein-Vertex mit
  Koeffizient ∝ Λ·G gibt (Maldacena 2011, Adamo–Mason 2013).

**Was daraus folgt.**

- Beide Wege enden bei derselben Art von Frage. Die zwei reinen Zahlen sind im Ledger beide an α gebunden: die
  Kopplung über die α-Gleichung, ε² über die Λ-Stufe q(1/α)/c₃². Der Engpass ist deshalb der offene Rest des α-Blocks
  selbst: „warum diese Funktion“ (EM.WARD.01, ALPHA.QUILLEN.EXACT.01), laut Ledger eine Seite von SEAM.EQUIV.01.
- **Lesart für diesen Rest [L, nicht geprüft]:** Der Rahmen hat drei Räume: die Gerade (2 reelle Dimensionen), die
  Raumzeit (4) und den Twistorraum (6). Anomaliepolynome haben dort die Grade 2, 3 und 4. Die α-Gleichung hat die
  Grade 3 (metrikunabhängig, „Chern-Niveau“), 2 (Randterm) und einen Logarithmus; ein Term vom Grad 4 fehlt, und
  genau der wird im Rahmen vom Axion aufgehoben (z163). Das ist nachträglich gelesen und zählt erst, wenn die
  Koeffizienten 1 und −2c₃³ daraus folgen.

**Nächste Rechnungen, nach Rang.**

1. **N0, die α-Funktion aus dem Rahmen.** Abstiegskette Twistorraum → Raumzeit → Gerade für den Hyperladungs-
   Hintergrund: Folgen die Grade 3 und 2 mit den Koeffizienten 1 und −2c₃³, und der Gegenterm 8b₁c₃⁶? Erfolg: P_g
   ist durch die α-Gleichung ersetzt. Abbruch: Die Koeffizienten folgen nicht; dann bleibt der Rahmen Form ohne Zahl.
2. **N1, die 4.** Trägt die Quellladung |μ₄| = 4 (v515) als Niveau in den Vertex ein (k_eff = 4)? Eine endliche
   Prüfung am Bestand des Vertrags. Wegen der Abweichung von 4 bis 6 % ist das nur mit Schwellen sinnvoll.
3. **N2, Einstein-Reduktion selbst rechnen.** Der Drei- und Vierpunkt-Vertex der konformen Theorie auf
   Einstein-Zuständen mit I_π − ε²I_ω; Koeffizient gegen Λ·G.
4. **W3 und W4** (Positivität, Materie) wie geplant.

### Nachtrag z165: einundzwanzigste externe Analyse und der Kandidat „relative graduierte Determinante“

**Anlass.** Eine externe Analyse zum Stand z163 prüft W1 bis W5 und schlägt als ein Objekt eine relative graduierte
Determinante über dem Raum der Twistorgeraden vor. Alle acht gerechneten Aussagen sind nachgerechnet und stimmen
(z165, 19 Prüfungen, reproduziert).

**Was die Analyse korrigiert.**

- **W3 war falsch gestellt [E negativ]:** Der Stromkern 1/(1 − z̄w)² ist positiv (det 7/9 an den Punkten 0 und 1/2).
  Mit dem Vertexzähler wird er (1 − z̄w)² und ist es nicht (det −7/16). Positivität gehört an den Propagator
  (Serre-Paarung mit Spiegelung, z123), nicht an den Vertex. Radiale Spiegelung der Naht, Serre-Paarung und
  vierdimensionale Zeitspiegelung sind drei verschiedene Operationen.
- **Doppelpol [E]:** Der Kern e^(−ωt)(A + Bt) hat die gespiegelte Determinante −B²(t₁ − t₂)²e^(−2ω(t₁+t₂)). Ein
  Doppelpol verletzt die Positivität; man kann der konformen Gravitation keinen Einstein-Anteil hinzufügen.
- **Hintergrund [E für A₁]:** Die aufgelöste ALE-Geometrie ist Ricci-flach und nicht konform flach (Eguchi–Hanson:
  Riem² = 384a⁸/r¹²). S⁴ und AdS₄ aus z164 sind konform flach. Verträglich sind beide nur am Orbifoldpunkt a = 0.
  Der Operator muss die Twistorgeometrie selbst deformieren können.
- **Anomalie und Kopplung:** Die Anomalie legt den nicht invarianten Teil fest, nicht den invarianten. Das deckt
  sich mit z164 (eine freie Zahl). Die Zahl verlangt die Metrik der Determinantenlinie (Bismut–Freed), den offenen
  Schritt von ALPHA.QUILLEN.EXACT.01 und SEAM.DETLINE.UNIFICATION.01.
- **W5 [E]:** Λ·G = 3/(32π³)·e^(−2/α) = 2,83508·10⁻¹²². Reiner elektromagnetischer Fluss ist nur für m ≡ 0 mod 3
  zulässig; m = 2 (α_m = 1/α) braucht Farbfluss und wird von der Eichgruppe nicht ausgezeichnet.
- **Eine Skala [E]:** Für ein positives Spektralprofil gilt f₀f₄ ≥ f₂². Das ist der strukturelle Grund hinter D2.

**Eigene Folgerechnungen.**

- **Gravitation ohne Doppelpol [E/C unter P_MM]:** Mit F = R − (Λ/3)e∧e ist εFF = Euler − (4Λ/3)(R − 2Λ). Geht die
  Gravitation wie Spin 1 als erste Ordnung plus Vertex b² ein, ist die Wirkung zweiter Ordnung, mit einfachem Pol.
  Das verlangt Λ ≠ 0, was der native Unendlichkeitstwistor liefert (z164). Die einzige Kopplung ist
  c = −3/(64πGΛ), mit der Leiter |c| = 10^119,7. Die native so(3,2)-Struktur steht im Bestand von H6.
- **Λ als relative Antwort [E/L]:** ρ_Λc₃²/M̄⁴ = φ_seam(1/α) − φ_base. Die Λ-Stufe ist die Abweichung der
  Nahtantwort von ihrem Basiswert am dualen Punkt. Am elektrischen Punkt ist dieselbe Abweichung 0,22 % des
  Basiswerts. Als absolutes Moment eines positiven Profils wäre Λ um 119 Größenordnungen ausgeschlossen.
- **Die duale Lesart ist nicht eindeutig [E/L]:** Derselbe Exponent ist 2α_m mit m = 2 oder die Maxwell-Wirkung des
  Flusses m = 4 auf der Einheitssphäre. Keine der beiden ist als reiner elektromagnetischer Fluss zulässig.
- **Native Graduierung [E; Bedeutung L]:** Sp₂₄₈(R) = Sp₂₄₈(R³) = 0 und Sp₂₄₈(R²) = −8. Im Uhr-verdrehten Sektor
  verschwindet die Zählspur (Volumenterm), die quadratische nicht: Sp(RX²) = 4p₂ − 12q₂. Die Graduierung ist die
  Uhr, nicht Supersymmetrie.
- **Der α-Block ist schon graduiert [E]:** Mit B = circ(2,−1,0,−1) (Sprungmatrix) und C = circ(0,1,2,1) (Greensche
  Matrix der Marken) gilt C = R²·(UBU), und R² vertauscht mit UBU. In H6 ist R² = (−1)^F. Daher
  det(1 − xC) = (1 − 4x)(1 + 2x)² gegen det(1 − xUBU) = (1 − 4x)(1 − 2x)²: Das Verhältnis sieht nur das ungerade
  Dublett. U ist der Charakter −1 der Vierteldrehung, derselbe wie der von q = dz²/(z⁴ − 1), der die Familien
  verdreht (W4).

**Pflichtenheft für das eine Objekt.** Γ[X] = −log(sdet D_X / sdet D_X*) mit dem uhrsymmetrischen Referenzpunkt X*.

1. Graduierung durch die Uhr (R² = (−1)^F); am Referenzpunkt kein Volumenterm.
2. Definiert auf allen zulässigen Hintergründen: Eichhintergrund, komplexe Struktur und die Auflösung weg vom
   Orbifoldpunkt.
3. Eichfelder und Gravitation in erster Ordnung mit quadratischem Vertex: einfache Pole, Gravitation mit
   so(3,2)-Krümmung.
4. Determinantenlinie mit Metrik und Zusammenhang; die α-Gleichung ist ihre Stationarität, die Λ-Stufe ihre
   relative Antwort.
5. Positivität am Propagator, mit getrennten Spiegelungen.
6. Materie als Nullmoden; Familien über den Twist mit dem Charakter von q.

**Nächste Rechnungen.**

1. **M1, Metrik der Linie.** Krümmung der Determinantenlinie der Nahtfamilie bei τ = i für den
   Hyperladungs-Hintergrund, gegen den Koeffizienten −2c₃³ der α-Gleichung. Erfolg: Der Blatt-Term folgt.
2. **M2, graduierte Kontaktentwicklung.** Gibt log der relativen Determinante von C gegen UBU mit
   x = c₃·ln 2 die Nahtantwort q(1 − q)^(−5/4)? Endliche Rechnung mit Ledgerdaten.
3. **M3, Einstein in Twistorvariablen.** Der Vertex b² mit der so(3,2)-Krümmung auf den Geraden; ersetzt N2.
4. **M4, Auflösung.** Erste Ordnung in den Auflösungsparametern: Wie ändert sich die relative Determinante?

**Ehrliche Einordnung.** Der Kandidat ist enger gefasst und an fünf Stellen mit dem Bestand verbunden. Das Objekt
selbst fehlt: Operator, Zustandsraum, Maß, Metrik der Linie und die Auswahl des zweiten Sektors. Eine volle Lösung
liegt nicht vor. Nichts ist befördert.

### Nachtrag z166: externe Arbeit „Quantisierung der Gravitation vom Nullkegel“ gegen den Bestand

**Anlass.** V. Chandrasekaran, „An algebraic quantization of gravity from the null cone“, arXiv:2611.34924
(eingereicht am 25.11.2026). Geprüft ist, ob die Arbeit dem Kandidaten aus z165 hilft (z166, 15 Prüfungen,
reproduziert).

**Einordnung der Quelle.** Die Arbeit ist drei Tage alt und nicht begutachtet. Nach eigener Angabe wurden mehrere
Rechnungen und Beweise zusammen mit einem Sprachmodell entwickelt. Sie zählt hier als ungeprüfte Literatur [Lit].
Gelesen sind Einleitung, 3.1, 4.3, 4.4, Teile von 5.1 bis 5.3, Abschnitt 6 und der Anfang von Anhang E. Nachgerechnet
sind fünf Formeln.

**Was die Arbeit tut [Lit].**

- Vakuumgravitation in vier Dimensionen mit Λ, Start von der Form erster Ordnung (Plebanski) mit Tetrade,
  Spinzusammenhang und Zweiform B.
- Die Zwangsbedingungen werden auf zwei sich schneidenden Nullkegeln eines Kausaldiamanten gelöst. Die Ecke S trägt
  die Boostladung K = A/(8πG).
- Gaußsche Zustände geben einen Hilbertraum, die Observablen eine Algebra vom Typ I. Summen über Zwischenzustände
  konvergieren absolut (Cauchy–Schwarz).
- Ein Diamant beschränkter Größe hat einen endlichdimensionalen Hilbertraum mit Dichteoperator e^(−A/(4Għ))/Z.
- Entropie A/(4Għ) + log Z, Logarithmus-Koeffizient 106/45 für Schwarzschild; Einstein-Gleichungen mit Λ bis eine
  Schleife.
- G und Λ sind Eingaben. Es gibt keine Materie, keine Eichfelder, keine Spinoren oder Twistoren.

**Nachgerechnet [E/N].**

- Paaramplitude: Σ_k 2(Δq_k)⁴e^(−2ik/s) mit (Δq_k)² = 1/(2sk) ist Li₂(e^(−2i/s))/(2s²), s = A/(4Għ).
- Entwickelt man zuerst in G, ist der Koeffizient dritter Ordnung −i·H_N, also log N. Summiert man zuerst, steht
  dort (−log s + ln 2 − 1)/s³: endlich für jedes G > 0.
- Boostgewicht e^(−2πK/ħ) = e^(−A/(4Għ)); S = ⟨s⟩ + log Z und dim = Z·Sp(ρe^s) für endliche Spektren.
- q/(2(1 − q)) − 1/2 = 106/45 für q = 257/302. Die Herkunft von 257/302 ist nicht geprüft.

**Brücken zum Bestand.**

- **Dieselbe Gravitationsform [E]:** c·εFF mit F = R − (Λ/3)e∧e und c = −3/(64πGΛ) (z165) ist die Wirkung (3.8) der
  Arbeit plus c mal Euler-Dichte. Die Arbeit quantisiert die Form, die Punkt 3 des Pflichtenhefts unter P_MM
  verlangt.
- **Positivität [Lit]:** Sie steht dort an den Zuständen: Präparation, verklebt mit ihrer Zeitspiegelung, als
  Axiom; die Zwischensummen sind durch Normen beschränkt. Das ist die Form, die Punkt 5 nach der Korrektur von W3
  verlangt.
- **Λ-Stufe als Zustandszahl [E/L]:** Mit der Λ-Stufe des Ledgers ist das Boostgewicht des kosmischen Horizonts
  s = 3π/(ΛG) = 32π⁴e^(2/α) = (3/8)/q(1/α) = 3,3243·10¹²². Gleichwertig: s·ρ_Λ/M̄⁴ = 24π² = (3/8)/c₃². Der
  Schleifenparameter der Arbeit, 1/s, ist (8/3) mal das Nahtnomen am dualen Punkt.
- **Hierarchie [E negativ]:** Die Nichtanalytizität der Arbeit ist ein Logarithmus. e^(−2/α) hat bei 0 lauter
  verschwindende Ableitungen. Der Mechanismus „erst summieren, dann entwickeln“ erzeugt die Hierarchie nicht.
- **Nahtfaktor [E/L]:** (1 − q)^(−5/4) = exp((5/4)Li₁(q)). Das ist eine Varianzsumme mit Gewichten 1/k; die
  Paaramplitude der Arbeit ist dieselbe Summe mit 1/k². Formgleichheit, keine Herleitung.
- **Offener Abgleich [O]:** Die Arbeit hat die Boostperiode 2π. Der KMS-Wert 4π des Viertelprozesses wäre in ihrer
  Normierung e^(−2s), das Gewicht von Sp(ρ²).
- **Scherung als Modul [E/L]:** Die konforme Metrik des Kegelschnitts ist |dx + τdy|²/Im τ mit τ = z² + i·e^(−z¹).
  Die Feldraummetrik der zwei Polarisationen ist die Poincaré-Metrik der oberen Halbebene. Scherungsfrei ist τ = i,
  der Punkt der Naht; sein Stabilisator ist S mit S⁴ = 1 und S² = −1, die Struktur der Uhr. Der E₈-Charakter
  j^(1/3) ist dort stationär. log|χ|² ist harmonisch, die zweiten Ableitungen sind −19,157 und +19,157: ein Sattel.
  Eine positive Normierung kann nur aus der Metrik der Determinantenlinie kommen.

**Urteil.** Die Arbeit hilft beim Gravitationsteil (Punkt 3), bei der Positivität (Punkt 5, W3) und gibt der Λ-Stufe
eine prüfbare Form. Sie hilft nicht bei den zwei freien reinen Zahlen, nicht bei der Hierarchie und nicht beim
gemeinsamen Objekt mit Materie (T3 bis T8). Bilanz unverändert 13 P / 7 L. Keine volle Lösung, nichts befördert.

**Neue Strategien aus dem Vergleich.**

- **S-I, Graviton als Modul der Naht.** Das Graviton ist die Schwankung von τ um i, je Lichtstrahl. Die Uhr ist die
  Symmetrie des scherungsfreien Punkts. Die Normierung des Gravitons käme aus derselben Linienmetrik wie der
  Blatt-Term der α-Gleichung (M1): eine Krümmung, zwei Zahlen.
- **S-II, Wechselwirkung aus der Zwangsbedingung.** Die Arbeit behält Gaußsche Zustände; die Nichtlinearität steckt
  in der Uhrgleichung (Scherung verlängert die Laufstrecke). Die Diagnose aus z162 (native Daten Gaußsch) ist dann
  kein Hindernis: Gesucht ist die Uhrgleichung der Naht, nicht ein nicht-Gaußscher Zustand.
- **S-III, Verklebungsbit aus der Zeitrichtung.** Die Klebregel der Arbeit setzt die Summe der Randladungen mit
  Außennormalen null. Die zwei Verklebungen ⟨(1,1)⟩ und ⟨(1,3)⟩ aus z144 unterscheiden sich durch ein Vorzeichen auf
  einer Seite. Eine Zeitorientierung könnte P_IF ersetzen.
- **S-IV, Λ als Zustandszahl.** log dim des kosmischen Diamanten = (3/8)/q(1/α) als Ziel einer Herleitung.

**Neue Rechnungen.**

1. **M5, Naht als Eckdaten.** Boostladung und Uhr der Arbeit gegen Uhr und KMS-Zustand der Naht; dabei der Abgleich
   2π gegen 4π. Erfolg: Die Naht ist die Ecke S mit festem Gewicht je Fläche. Scheitert der Abgleich, trägt die
   Brücke nicht.
2. **M6, Zustandszahl.** Folgt log dim = (3/8)/q(1/α) aus der relativen Nahtantwort? Verbindet M2 mit W5.
3. **M7, Graviton als Modul.** Zweite Ableitung der Nahtdeterminante mit Linienmetrik bei τ = i. Erfolg: eine
   positive Form mit festem Koeffizienten; zusammen mit M1 ein festes Verhältnis von Eich- und
   Gravitationsnormierung.
4. **M8, Uhrgleichung.** Uhroperator der Naht mit der Fokussierungsgleichung; niedrigster Vertex gegen die
   Spin-2-Form aus z164.
5. **M9, Verklebungsbit.** Klebregel mit Außennormalen in der Gittersprache von z144: Wählt sie eine der zwei
   Lagrangeschen Untergruppen?

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
