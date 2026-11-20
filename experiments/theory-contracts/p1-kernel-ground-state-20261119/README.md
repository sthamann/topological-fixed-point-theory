# P1KERN20261119 · Trägt der originale P1 Kern den Grundzustand? Kugel, Torus, Uhr und Schnitt

Datum: 19. November 2026. Auftrag: den offenen Test C_nat ausführen und selbstständig weitere
Szenarien entwickeln und prüfen. Kein Ledgerwechsel, keine Beförderung, kein Commit.

Ergebnis: 77 von 77 Prüfungen bestanden (acht Skripte). Alle positiven Aussagen sind bedingt.
Die Bedingungen stehen bei jedem Punkt.

## Was geprüft wurde

Quelle des Kerns: `documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-…tex`.
Dort steht: einseitiger Randoperator `B_Σ` (Dirac Typ), Calderón Projektor, Saatantwort
`χ_seed = 1/λ₁⁺(|B_Σ|)` und in der primitiven Wirkung der Phasenterm `(iπ/2)·Δη`.
Vergleichskern: der skalare Kugelrandkern aus `research/originlab/F1-ORIGIN-SCENARIOS-20261116.txt`, Abschnitt 5.

## 1. C_nat Test (`check_cnat_channels.py`, 12/12)

- Auf der runden Kugel S² trägt **jeder** Winkelkanal genau das Gewicht 1/(4π) = 2c₃
  (skalarer Kern bei κ = 1). Das Spektrum ist flach.
- Der originale Kern (Dirac, positive Hälfte, also die Calderón Seite) ist ebenfalls flach:
  1/(2π) = 4c₃ je Kanal, also 2c₃ je Spinorkomponente. Hier gibt es kein κ.
- Flach ist das Spektrum **nur auf S²** (geprüft S¹ bis S⁷). Das ist eine echte Auswahl der Dimension.
- Jeder Kanal verhält sich wie eine Schwingung mit h = 1. Die Kugelfläche 4π liefert diese Einheit.
- Ein positiver Dirac Kanal gibt exp(−2α) exakt. Ein skalarer Kanal gäbe exp(−α) (ausgeschlossen, z = +12042),
  beide Dirac Vorzeichen exp(−4α) (ausgeschlossen, z = −23816).
- Die Gegenpunktauslese hat den Abelwert 1/2 und gibt c₃.
- Grenze: Im skalaren Kern bleibt κ eine Eingabe. Die regularisierte Summe am selben Punkt ist nicht eindeutig.

Urteil: **bedingt bestanden**. Der originale Kern trägt genau das Gewicht des Grundzustands,
wenn die Schnittfläche eine runde S² ist und nur die positive Hälfte zählt. Eine Operatorgleichheit
mit dem vollen P1 Objekt ist damit nicht bewiesen, nur Spektrum und Gewicht.

## 2. Brücke Kugel ↔ Torus (`check_torus_bridge.py`, 13/13)

- Die Wärmespur des originalen Randoperators auf S² ist `Z⁺(β) = 2q/(1−q)²`.
- Summiert über alle Umläufe m des Zeitkreises ergibt sie exakt die Eisensteinreihe:
  `E₂(q) = 1 − 12 Σ_m Z⁺(mβ)`.
- Am quadratischen Torus (β = 2π) ist `E₂ = 3/π`. Daraus: Schwingungsenergie einschließlich
  Casimir Anteil = **−c₃**. Die Nullmode trägt +T/4 = **+c₃**. Die Summe ist null.
- Diese Balance ist erzwungen: `E₂*` hat Gewicht 2 und verschwindet am Fixpunkt der Vierteldrehung.
- Sechseck: (2/√3)·c₃, also nicht c₃.
- **Die drei Wege zu c₃ sind ein und dieselbe Zahl**: KMS Weg 1/(Iβ), Grundzustand mit h = 1,
  Toruswert E₂(i)/24. Alle sind (1/2π)/4. Sie zählen deshalb als ein Beleg, nicht als drei.
- I = 4 zerfällt in 2 × 2 (ein Gaußintegral der Nullmode, chirale Hälfte) oder gleichwertig
  in die vier reellen Richtungen von C².
- exp(−2α) ist der thermische Debye Waller Faktor einer Nullmode bei T = 1/(2π).
- (1−q)^(−5/4) ist die Zweipunktfunktion des D₅ Spinfelds (Gewicht 5/8). D₄ (−47σ) und D₆ (+51σ) fallen.
- Ob dahinter ein volles Modenprodukt steht, kann α nicht entscheiden (Verschiebung 0,03σ).

## 3. Vierteluhr, Horizont, Potential (`check_clock_horizon_ufield.py`, 14/14)

- Der originale Randoperator hat auf S² ganzzahliges Spektrum. `exp(iπB/2)` hat Ordnung 4
  mit allen vier Klassen. Der skalare Kern kennt nur zwei Klassen. S¹ und S³ geben Ordnung 8.
- |B| auf S² hat dasselbe Spektrum wie die zweidimensionale Schwingung samt Nullpunkt
  (Energie n+1, Vielfachheit n+1), mal zwei für den Spin. `exp(−iπ|B|/2)` ist die
  zweidimensionale Fouriertransformation bis auf die Phase −i.
- Horizont zuerst (G3): ein quantenmechanischer thermischer Zustand bei β = 2π gibt je Quadratur
  mindestens 2c₃, nie c₃, und die Dämpfung exp(−6,31α) (ausgeschlossen). Der Horizontweg trägt
  nur als klassische Nullmode. Das ist dieselbe Zahl wie oben.
- α als Feld (gap-ufield): U mit U′ = F hat genau ein Minimum auf (0,1), keine Barriere,
  kein zweites Vakuum. Krümmung 0,9885·3α².

## 4. Ist das Quadrat stabil? (`check_square_stability.py`, 7/7)

- τ = i ist ein **Sattel** der modulinvarianten freien Energie. Entlang der Temperatur stabil,
  entlang der Verdrehung instabil. Das echte Extremum ist das Sechseck.
- Friert die Spiegelung die Verdrehung ein (rechteckige Tori), ist das Quadrat das einzige Extremum.
- Folge: Die Spiegelung des P1 Kragens ist nötig, damit c₃ nicht zu (2/√3)c₃ rollt.

## 5. Eigenphase und Periode (`check_eta_phase_and_period.py`, 9/9)

- Dirac η von S³/Z_p direkt aus dem Spektrum: |η| = (p²−1)/(6p).
- Die Eigenphase exp(iπη/2) der primitiven Wirkung: runde Kugel 1, RP³ Ordnung 16,
  Z₄ Schnitt Ordnung 32, **Z₅ Schnitt Ordnung 5**.
- β = 2π ist die Eigenperiode des Randoperators. Der Radius fällt heraus.
- c₃ = λ₁/(I·β) mit λ₁ = 1, I = 4, β = 2π. Die Balance gilt nur bei β = 2π.

## 6. Alle ADE Schnitte (`check_ade_cuts.py`, 7/7)

- E₆, E₇, E₈: η = 167/144, 383/288, 1079/720. Phasenordnungen 576, 1152, 2880. Keine Uhr.
- A Reihe bis p = 120 und D Reihe bis n = 60: nur **Z₅ = A₄** gibt eine Uhrordnung, nämlich 5.
- Allgemein: Ordnung p genau dann, wenn 24 | p²−1. Kleinster Fall p = 5, danach 7, 11, 13.

## 7. Schnitt gegen Klebedaten von E₈ (`check_cut_glue_match.py`, 9/9) · nur Muster

- Z₄ Schnitt: die zwei Spinstrukturen geben |η| = 5/8 und 3/8. Das sind genau die Gewichte der
  Spinfelder von D₅ und A₃ (Normen 5/4 und 3/4), also die Klebedaten der Karte D₅ ⊕ A₃ in E₈.
- Die Summe η₊ + |η₋| ist p/4. Gewicht 1 (ein Klebevektor der Norm 2) gibt es nur bei p = 4.
  Als Teilung D_n ⊕ D_(8−n) wählt das Paar n = 5.
- Z₅ Schnitt: η = 4/5 = Norm der A₄ Klebeklasse. Gleichheit η = (p−1)/p gilt nur bei p = 5.
- Ehrliche Grenze: Z₄ trifft auf der Stufe „η = Norm/2“, Z₅ auf der Stufe „η = Norm“.
  Eine gemeinsame Regel gibt es nicht. Vielfache von 1/8 treten auch bei p = 2 auf.
  Das ist ein Muster, keine Herleitung.

## 8. Erzwingt das Original die runde S² und die positive Hälfte? (`check_original_forcing.py`, 6/6)

Gelesen im Original (Randkernquelle, Abschnitte 001 und 002; Papier Abschnitt 013):
- Die Schnittfläche Σ ist dort nur „glatt, kompakt, orientiert“. Rund steht nirgends.
- Die **positive Hälfte steht im Original**: Polarisation = Calderón Projektor des APS Randproblems,
  Saatantwort aus dem ersten positiven Eigenwert λ₁⁺.
- Die markierte Geometrie ist eine **flache** Kugel mit vier Marken (Defizit π je Marke, Gauß Bonnet),
  gemeinsamer Divisor μ₄ in P¹. Das ist ein Kissen T²/Z₂.

Gerechnet:
- Die Doppelüberlagerung von P¹ über μ₄ hat Doppelverhältnis 2 und j = 1728: sie ist **exakt der
  quadratische Torus τ = i**. Der Torus muss also nicht über die Eigenperiode einer runden Kugel kommen.
- Liegen die vier Marken auf einem Kreis (reelles Doppelverhältnis), ist der Torus rechteckig, und die
  freie Energie hat ihr einziges Extremum bei λ = 1/2, also bei μ₄.
- Verschiebt man die Marken, ändert sich die Zahl (λ = 1/4: −0,0413 statt −c₃). Die Lage der Marken trägt alles.
- „Rund“ ist im Mittel unnötig: Auf jedem zweidimensionalen Schnitt gibt das Weylgesetz 1/(2π) = 4c₃
  je Spektraleinheit (am flachen Torus numerisch 0,1577 gegen 0,1592 bei K = 400). Die runde S² hat es exakt je Kanal.
- Am flachen Kissen ist das Spektrum nicht ganzzahlig. Die Vierteluhr wirkt dort geometrisch
  (Multiplikation mit i auf Z[i]), nicht als exp(iπB/2).

Urteil: **Positive Hälfte: ja, steht im Original. Runde S²: nein, nicht erzwungen und nicht nötig.**
Nötig ist stattdessen „zweidimensionaler Schnitt“ plus „Marken bei μ₄“. Die Lage der Marken ist laut
Atlasdiagnose vom 16.11. durch P1, P2 und E₈ nicht festgelegt. Neuer Kandidat für ihre Auswahl:
Spiegelung (Marken auf einem Kreis) plus Extremum der freien Energie. Das ist eine Hypothese, kein Beweis.

## Arbeitsbild

Nach Abschnitt 8 kürzer: flaches Kissen mit vier Marken bei μ₄ → Doppelüberlagerung ist der quadratische
Torus → Nullmode trägt T/4 = c₃. Die runde S² ist nur das exakt lösbare Vergleichsmodell. Die Fünferuhr sitzt nicht im
E₈ Schnitt, sondern im A₄ Schnitt (Z₅), passend zu „Ordnung 5 braucht Rang 4“ vom selben Tag.

## Nicht belegt

- Dass die vier Marken bei μ₄ liegen (die runde S² ist nicht nötig, die positive Hälfte steht im Original).
- Warum genau eine chirale Nullmode zählt.
- δ_top = 48c₃⁴ und φ_tree = (4/3)c₃ sind in diesem Bild nicht erklärt (keine Mechanik gefunden).
- Die Uhrenliste {4, 5, 8, 12, 15, 20, 24, 30} wurde vor dem ADE Lauf festgelegt, stammt aber aus
  dem Vertrag LANDKARTE20261119 desselben Tages.
- Der Bezug zu S = A/4 ist nur eine Rechenlesart.
- Der gewichtete Originalschnitt (R113) bleibt blockiert.

## Wiederholen

```
cd experiments/theory-contracts/p1-kernel-ground-state-20261119
for s in cnat_channels torus_bridge clock_horizon_ufield square_stability eta_phase_and_period ade_cuts cut_glue_match original_forcing; do python3 -B check_$s.py; done
```
