# TFPT · Ergebnisse der Versuche vom 18. November 2026

Grundlage ist das Versuchspaket „TFPT · Experimente aus der Unterhaltung vom 18. November 2026“ (F1 bis F10, KI1 bis KI6, S1 bis S3) und der anschließende Teil der Unterhaltung über Quantencomputer, Zehn‑Linien‑Geometrie, Reflexionsdynamik und Entfernung aus Korrelationen. Jeder Versuch hat ein eigenes Skript in diesem Ordner, die Rohwerte liegen unter `out/`.

Alle Läufe sind Nebenversuche. Es wurde nichts in den Atlas, das Ledger oder den Fundamentalplan eingetragen, und keine Aussage der Theorie ändert dadurch ihren Status.

## Kurzfazit

1. **Die Mathematik hält.** Alles, was die Quelldokumente über Quartikcode, Projektor, Zehn‑Linien‑Geometrie, Rückkehrpulse und E₈‑Abschluss behaupten, ließ sich unabhängig nachrechnen, meist auf Maschinengenauigkeit.
2. **Für KI liegt der Nutzen in bekannter Codierungstheorie.** Der Rang‑5‑Projektor ist auf klassischen Zuständen exakt die XOR‑Paritätsprüfung von vier Symbolen. Ein Netz mit diesem Baustein lernt passende Regelwelten mit deutlich weniger Daten. Ein Netz mit der gewöhnlichen, gleich informierten Paritätsnachricht ist aber in jedem Vergleich mindestens gleich gut und meist besser. Ein eigener TFPT‑Vorteil gegenüber diesem Vergleich ist nicht belegt.
3. **Die Schutzenergie H_mom taugt nicht als Suchenergie.** Sie kann gültige und ungültige Viererzellen nicht trennen. Brauchbar sind nur die Teile, die gewöhnlichen Nebenbedingungen entsprechen.
4. **Als Quantencode ist der Fünferraum echt besonders, aber klein.** Er erkennt jeden Fehler auf einem Register und stellt ein bekannt verlorenes Register exakt wieder her. Zufällige Rang‑5‑Räume im selben symmetrischen Raum können das nicht. Der Abstand ist 2, wie beim Lehrbuchcode [[4,2,2]].
5. **Auf zwei Qubits ist die Grundtatsache messbar**: Vierter Moment 1/32 statt 1/35. Das ist die bekannte Aussage, dass die Cliffordgruppe ein 3‑Design und kein 4‑Design ist. Eine Messung prüft Quantenmechanik und Gerät, keine TFPT‑Vorhersage.
6. **Physik**: Drei Raumdimensionen, vier Marken, E₈ und die Ladungen folgen in den Rechnungen jeweils aus einer zusätzlich gesetzten Regel, die benannt ist. Massen, Spin 2, Kosmologie und die Zeitkarte sind **nicht gerechnet**, weil die gewichtete Originalquelle (R113) nicht ausgewertet vorliegt.
7. **Signaturtests**: Für S1 bis S3 gibt es Empfindlichkeitsrechnungen. Ein echter Test braucht bei S1 die Auger‑Daten, bei S2 und S3 einen vorab festgelegten Effekt, den die Unterhaltung nicht liefert.

## Übersicht

| Versuch | Skript | Stand | Kernaussage |
|---|---|---|---|
| F1 Viermarkendeterminante, Dirac, Sektorgewicht | `f1_marks_dirac_weight.py` | gerechnet (eigene Rekonstruktion der Formel) | Markenseite legt Produkt und Sektorverhältnisse fest, drei Daten bleiben offen |
| F2 E₈ durch Abschluss | `f2_f3_f6_lattice_algebra.py` | gerechnet | 112 Wurzeln plus ein Spinorstrom schließen zu 240, mit Cartan 248 |
| F3 Vier Marken aus einer Regel | `f2_f3_f6_lattice_algebra.py` | gerechnet, in einer Klasse | nur m = 4 erlaubt die gerade unimodulare Verklebung |
| F4 Drei Raumdimensionen | `f4_spectral_dimension.py` | gerechnet | 3 folgt aus der Nachbarschaftsregel, nicht aus den Zellen |
| F5 Zeit, Positivität, Entropiepfeil | `f5_readers_time_positivity.py` | Teil gerechnet | Leser mit 1 bis 3 Registern sehen nichts, Zeitkarte ungerechnet |
| F6 Materie, Familien | `f2_f3_f6_lattice_algebra.py` | Teil gerechnet | Ladungen und Anomalien stimmen, Chiralität netto null |
| F7 Massen und Wechselwirkungen | – | **ungerechnet** | braucht Full96 mit offenen Bosonendpunkten |
| F8 Spin 2, Einstein | – | **ungerechnet** | braucht Zeitkarte und 3+1D‑Spektrum |
| F9 Kosmologische Naht | – | **ungerechnet** | braucht gewichtete P1‑Naht |
| F10 Naturkonstanten | `f10_alpha_reproduction.py` | Reproduktion | α⁻¹ liegt 1,9 σ über CODATA 2022 |
| Q1 Zwei‑Qubit‑Messung | `q1_two_qubit.py` | simuliert | 1/32 gegen 1/35 mit einigen tausend Schüssen trennbar |
| Q2 Quartikcode als Fehlerkorrektur | `q2_code_error_correction.py` | gerechnet | Abstand 2, ein verlorenes Register exakt zurück |
| Q3 Rückkehrpulse, Pulsreflexion | `q3_q4_q5_dynamics_geometry.py` | gerechnet | Formeln exakt, 1 % Kopplungsfehler kostet etwa 0,15 % |
| Q4 Zehn gleichwinklige Linien | `q3_q4_q5_dynamics_geometry.py` | gerechnet | optimal (Welch‑Schranke), frei optimierte Rahmen erreichen dasselbe |
| Q5 Entfernung aus Verschränkung | `q3_q4_q5_dynamics_geometry.py` | gerechnet | zurückgelesen wird der eingesetzte Kopplungsgraph |
| KI1 Weltmodell, Planung | `e1_three_way.py`, `ki1_generic_control.py`, `ki1_planning.py`, `e7_world_model_consistency.py` | gerechnet, verkleinert | kein eigener TFPT‑Gewinn gegen die gleich informierte gewöhnliche Regelnachricht |
| KI2 A₃‑Zellen | Vorlauf `a3_pilot.py`, `e3_invariants.py` | gerechnet | A₃‑Projektion ist die gewöhnliche K₄‑Graphprojektion |
| KI3 Lokale Selbstorganisation | `e2_local_vs_global.py`, `ki3_async_local.py` | gerechnet | lokale Regeln genügen mit Regelwissen, lange Schlussketten sind empfindlich gegen Verspätung |
| KI4 Quartikcode, Schutzenergie | `e0_structure.py`, `e4_energy_reasoning.py`, `ki4b_energy_candidates.py`, `ki4c_ki7_projector_reflection.py`, `q2_code_error_correction.py` | gerechnet | Projektion hilft nur bei passenden Daten, H_mom hilft nicht |
| KI5 Adaptive Rechenzeit | `e5_adaptive_compute.py` | gerechnet | Stopp mit Konsistenzprüfung spart 46 bis 82 % Runden ohne falschen Stopp |
| KI6 E₈‑Symmetrie | `ki6_e8_equivariance.py`, `e6_e8_representations.py` | gerechnet | hilft nur bei Aufgaben über das Gitter selbst, sonst reicht allgemeine Geometrie |
| KI7 Reflexionsdynamik | `ki4c_ki7_projector_reflection.py` | gerechnet | hält Information über 300 Schritte, bekannte Architekturklasse |
| S1 Kubisches Himmelsmuster | `s1_s2_s3_signatures.py` | Empfindlichkeit | ab etwa 20 000 Ereignissen ist ein Muster von 25 % Amplitude auffindbar, echte Daten nicht geladen |
| S2 Botschaft in Konstanten | `s1_s2_s3_signatures.py` | Empfindlichkeit | mindestens 10 unabhängige stabile Bits nötig, vorhanden ist 1 |
| S3 Getrennter Eingriff | `s1_s2_s3_signatures.py` | Empfindlichkeit | Blockzahl wächst mit 1/δ² |

## Physik

### F1 · Viermarkendeterminante, Dirac und Sektorgewicht

Die Unterhaltung nennt keine Formel. Gerechnet wurde der naheliegende Kandidat: die Marken e₁, e₂, e₃ und ∞ der Weierstraß‑Kurve, deren Abstände d_ij² den achten Potenzen der drei geraden Thetafunktionen entsprechen. Geprüft mit 95 Stellen an zwei reellen und zwei komplexen Punkten.

- Vandermonde‑Quadrat der Marken gleich π¹² η²⁴ Π(θ_ν/η)⁸, das Produkt der drei Dirac‑Determinanten ist 2, das Doppelverhältnis ist λ(τ). Abweichung jeweils unter 10⁻⁹⁴.
- Der E₈‑Charakter lässt sich vollständig durch die Marken schreiben: χ = 2^(5/3) Σd² / (Πd²)^(1/3). Koeffizienten 1, 248, 4124, 34752, mit Spinorsektor 248, ohne 120.
- Die Suche über Spinorgewichte w zeigt: nur w = 1 ist unter S invariant (Defekt 3·10⁻⁹⁵, sonst 0,25 bis 0,54).

**Erkenntnis.** Die Markenseite legt das Produkt und die Verhältnisse der Sektoren fest. Nicht festgelegt werden durch die Determinantenrelation: das Gesamtgewicht ½, der Zweig der dritten Wurzel (c/24 = 1/3) und die Wahl S₊ gegen S₋. Genau diese drei Daten müsste der Originalbeleg mit Gewicht und Normierung liefern.

**Grenze.** Die Formel ist meine Rekonstruktion. Ob sie der in der Unterhaltung gemeinte Viermarkenbeleg ist, entscheidet erst die gebundene Originalquelle.

### F2 · E₈ durch Abschluss

Die 112 Wurzeln von D₈ plus **ein** Spinorstrom der Norm 2 schließen unter Addition in zwei Schritten zu 240 Wurzeln, mit den 8 Cartan‑Strömen 248. Ohne Kandidat bleibt es bei 112. Der konjugierte Spinor liefert ebenfalls 240. Beide Spinoren zusammen scheitern, weil ihr Skalarprodukt 3/2 nicht ganzzahlig ist. Über die Naht D₅ ⊕ A₃ (52 Wurzeln) und einen Klebevektor entsteht dasselbe in drei Schritten. Die Jacobi‑Identität wurde an 19 742 Wurzeltripeln ohne Verletzung geprüft.

**Erkenntnis.** Der Abschluss ist zwingend, sobald der eine Spinorstrom zugelassen ist. Die Zulassung selbst ist die Eingabe.

### F3 · Vier Marken aus einer Regel

In der Klasse D_(m+1) ⊕ A_(m−1) für m = 2 bis 40 hat nur m = 4 passende Diskriminantengruppen und eine gerade unimodulare Verklebung. Für Worte aus zwei Involutionen gilt: (AB)² = +1 ergibt die Kleinsche Vierergruppe mit vier eindimensionalen Charakteren, (AB)² = −1 ergibt acht Elemente mit einer zweidimensionalen Darstellung.

**Grenze.** Die Aussage gilt in dieser einen Gitterklasse. Die Gewichtung der Übergänge ist nicht hergeleitet. Das geordnete ABAB‑Wort am gewichteten Zweifensterkern ist nicht gerechnet.

### F4 · Drei große Raumdimensionen

Gemessen wurde die Spektraldimension über den Wärmekern, geeicht an Z², Z³ und Z⁴.

| Aufbau | Spektraldimension |
|---|---|
| Differenzgraph von 3 Marken (A₂) | 2,01 |
| Differenzgraph von 4 Marken (A₃) | 3,45 → 3,11 mit wachsender Größe, wie die Eichung Z³ |
| Differenzgraph von 5 Marken (A₄) | 5,35 → 4,43 |
| Viererzellen im Gitter verklebt | 2,02 |
| Viererzellen zufällig verklebt | 5,6 bis 6,2, kein fester Wert |
| A₃ mit zufällig umverdrahteten Kanten | 6,6 bis 9,2, fällt dann ab |

**Erkenntnis.** Vier Marken geben die Dimension 3 nur über die Lesart „Nachbar heißt eine Markendifferenz“. Dieselben Zellen anders verklebt geben 2 oder gar keine Dimension. Die Drei steckt in der Nachbarschaftsregel. Das Fenster 2,8 bis 3,2 wird ab 8 000 Knoten erreicht.

### F5 · Zeit, Quantenmechanik und Entropiepfeil

Gerechnet am Quartikcode selbst, ohne die R113‑Fenster.

- **Eingeschränkte Leser.** Wer 1, 2 oder 3 der vier Register sieht, kann den Codezustand nicht vom symmetrischen Zustand unterscheiden (relative Entropie 0). Erst alle vier Register zeigen log₂ 7 = 2,807 Bit.
- **Globale Kontrolle.** Ein Paarpuls lässt die globale Entropie bei 0, das Leck aus dem Code folgt exakt (8/9) sin² t.
- **Positivität.** Für die Zwei‑Zellen‑Kopplung mit dem dokumentierten Vorzeichen H = −K ist der Grundzustand eindeutig, die Lücke 2,5 und die reflektierte Gramform positiv (kleinster Eigenwert 0,2). Mit umgekehrtem Vorzeichen ist der Grundzustand 15‑fach entartet und die Gramform hat einen negativen Modus (−0,41).

**Erkenntnis.** Ein lokaler Informationsverlust entsteht hier allein durch die Einschränkung des Lesers, global bleibt alles umkehrbar. Das Vorzeichen der Kopplung entscheidet über Positivität.

**Ungerechnet.** Die R113‑Fenster, die Zeitkarte T(s+t) = T(s)T(t), der Kausalkegel und die Lorentzsignatur. Übergabepunkt: die tatsächliche gewichtete reflektierte Originalform Q₀ mit zwei getrennt markierbaren Fenstern und absoluter Leerreferenz.

### F6 · Materie und drei chirale Familien

Unter D₅ ⊕ A₃ zerfallen die 240 Wurzeln in 64 + 64 + 40 + 12 + 60. Das sind vier 16 und vier 16‑quer, die Nettochiralität ist null. Die Ladungen sind 0, ±1/3, ±2/3, ±1, alle Anomaliesummen verschwinden, sin² θ_W = 3/8 an der Vereinigungsskala.

**Erkenntnis.** Die Algebra liefert die richtigen Ladungen, aber von sich aus weder Chiralität noch die Zahl drei. Beides muss aus der Auswahl der Quelle kommen.

**Ungerechnet.** Pole und Residuen bei wachsendem Volumen, Spektralfluss, Spiegelentkopplung und die typabhängige Rückkehr an der frühen Higgsmarke. Übergabepunkt: die gemeinsame gewichtete R113‑Quelle bei derselben Zeit und äußerer Auswahl.

### F7, F8, F9 · Massen, Spin 2, kosmologische Naht

Diese drei Auswertungen bleiben **unberechnet**. Sie lassen sich nicht ersatzweise an einem Spielzeug nachstellen, ohne das Ergebnis vorwegzunehmen.

| Versuch | Fehlende Eingabe (Übergabepunkt laut Atlas, Stand 16. November 2026) |
|---|---|
| F7 | Die zwei offenen Bosonendpunkte in I_A(w_A s_A) und die gebundene Referenzkontraktion N₁ᴿ. Die Full96‑Fermionform mit Majorana‑Anteil und die 15 gemeinsamen Radialproben liegen vor. |
| F8 | Die akzeptierte Zeitkarte aus F5 und ein gemeinsames physisches 3+1D‑Spektrum mit thermodynamischem Limes. Diese Aufgabe steht im Atlas als blockiert. |
| F9 | Ein aus derselben Quelle berechneter Hintergrund auf beiden Seiten der gewichteten P1‑Naht und die Übergangskarte für Störungen. Ohne sie gibt es keine Signatur, die man in Planck‑Daten einspeisen könnte. |

### F10 · Naturkonstanten

Die vorhandene Rechnung `verification/v3_em_alpha.py` läuft durch (14 Prüfungen bestanden) und liefert α⁻¹ = 137,035 999 216 8.

| Vergleich | Wert | Abstand |
|---|---|---|
| CODATA 2022 | 137,035 999 177(21) | +1,9 σ |
| Rubidium‑Rückstoß | 137,035 999 206(11) | +1,0 σ |
| Elektron g−2 | 137,035 999 166(15) | +3,4 σ |
| Cäsium‑Rückstoß | 137,035 999 046(27) | +6,3 σ |

**Erkenntnis.** Der Wert liegt nahe an der Rubidium‑Messung. Rubidium und Cäsium widersprechen einander selbst um 5,5 σ. Ob TFPT hier richtig liegt, entscheidet die nächste unabhängige Messung.

**Grenze.** Das ist eine Reproduktion, kein Vorhersagetest. Ein eingefrorenes Bündel mehrerer Größen mit gemeinsamer Kovarianz wurde nicht gebildet. Die Werte für CODATA 2022 und Elektron g−2 sind gegen die NIST‑Tabelle und die Literatur geprüft. Die Werte für Rubidium und Cäsium sind aus dem Gedächtnis eingetragen und noch gegen die Originalarbeiten zu prüfen.

## Quantenversuche

### Q1 · Was sich auf einem echten Zwei‑Qubit‑System messen lässt

Ein TFPT‑Register ist C⁴, also zwei Qubits. Die 60 Quellstrahlen sind die 60 Stabilisatorzustände zweier Qubits.

- Die Rahmenpotentiale stimmen für t = 1, 2, 3 mit dem Zufallswert überein und weichen bei t = 4 ab: **1/32 statt 1/35**.
- Messvorschrift: Prüfzustand präparieren, in den 15 Stabilisatorbasen messen, den Mittelwert von p⁴ schätzen. Stabilisator‑ und Bellzustand geben 0,03125, der magische Zustand T⊗T gibt 0,02852.
- Simulation mit 5 000 Schüssen je Basis: Trennung 7,7 σ ohne Rauschen, 6,1 σ bei 3 % Depolarisation und 2 % Auslesefehler. Rauschen verschiebt allerdings die Absolutwerte (Stabilisator 0,0277), ein Gerät muss also geeicht sein.
- Die Reflexion I − 2|Singulett⟩⟨Singulett| ist das SWAP‑Gatter. Als Orakel‑ und Diffusionspaar findet sie einen von vier Einträgen in einem Schritt mit Wahrscheinlichkeit 1.

**Erkenntnis.** Der Versuch ist heute auf jedem Zwei‑Qubit‑Gerät machbar. Er bestätigt eine bekannte Eigenschaft der Cliffordgruppe. Eine TFPT‑eigene Quantenvorhersage bräuchte eine Zahl, die die Standardtheorie nicht liefert.

**Nicht gemacht.** Es wurde keine Hardware angesteuert.

### Q2 · Quartikcode als Fehlerkorrektur

| Code | Verletzung bei 1 Register | Rückgewinnung nach Verlust eines Registers |
|---|---|---|
| TFPT‑Quartikcode | 5·10⁻¹⁷ | 1,000 (auch im schlechtesten Fall) |
| zufälliger Rang 5 im selben Sym⁴ | ≥ 0,29 | 0,63 |
| zufälliger Rang 5 im vollen Raum | 0,19 bis 0,23 | 0,82 |

Bei zwei Registern liegt die Verletzung bei 0,8, der Abstand ist also 2. Mit Nachauswahl im verrauschten Speicher: bei 3 % Fehler je Register werden 89 % akzeptiert, die Güte ist 0,9996 gegen 0,885 ungeschützt.

**Erkenntnis.** Das ist der eine Versuch, in dem der TFPT‑Raum gegen den strengen Vergleich (Zufallsraum im selben symmetrischen Raum) klar gewinnt. In absoluten Zahlen ist es ein kleiner Erkennungscode: 8 Qubits für 2,32 logische Qubits bei Abstand 2.

### Q3 · Rückkehrpulse und Pulsreflexion

Zwei Zellen mit Kopplung g. Für m = 2, 3, 5 Zyklen stimmen Gesamtzeit (m−1)π/Δ, Rückkehr in den Code (Leck 10⁻¹⁵) und größtes Leck unterwegs (0,75, 0,556, 0,36) mit den Formeln überein. Der aktive Block bekommt exakt ein Minuszeichen (21 Eigenphasen +1, 4 Eigenphasen −1). Ein Kopplungsfehler von 1 % hinterlässt ein Leck von 0,12 bis 0,17 %, einer von 5 % etwa 3 bis 4 %. Die Pulsreflexion hat das Spektrum 0⁷⁰, 2⁵ und bildet J₀ auf W ab.

**Erkenntnis.** Die Pulse sind ein exaktes, einstellbares Zwei‑Zellen‑Gatter. Sie sind nicht von sich aus fehlertolerant, der Fehler wächst quadratisch mit der Fehlkalibrierung.

**Grenze.** Die Pulsreflexion wurde für allgemeine Isometrien mit J₀†W = sI geprüft, der tatsächliche Encoder W_ev wurde nicht nachgebaut.

### Q4 · Zehn gleichwinklige Linien

C² = 9I gilt, alle Überlappungen haben Betrag 1/6, es gibt sechs Vorzeichenrahmen und alle sind Petersen‑Graphen. Die Kohärenz 1/3 erreicht die Welch‑Schranke, ist also bestmöglich. Als Überlagerungsspeicher (10 Einträge in 5 Dimensionen) bei Rauschen 0,2: TFPT‑Linien 98,7 %, frei optimierter Rahmen 98,4 %, Zufallsrahmen 90,7 %.

**Erkenntnis.** Die Geometrie ist optimal. Eine numerische Optimierung findet aber gleich gute Rahmen (beste Kohärenz 0,3334). Der Vorteil liegt darin, die Lösung geschlossen zu kennen, nicht in besserer Leistung.

### Q5 · Entfernung aus Verschränkung

Logische Zellen auf Ring und Linie mit der Kopplung −K zwischen Nachbarn. Die gegenseitige Information fällt auf dem Sechserring mit dem Graphabstand (1,60, 0,58, 0,51 Bit). Auf der offenen Linie sind die beiden Enden etwas stärker korreliert als Abstand 4, und der Siebenerring hat einen entarteten Grundzustand.

**Erkenntnis.** Aus Korrelationen lässt sich ein Abstand zurücklesen. Es ist der Abstand des Kopplungsgraphen, den man selbst eingesetzt hat. Über eine Raumdimension sagen diese Größen nichts.

## KI‑Versuche

Gemeinsamer Aufbau: Eine Welt ist eine Belegung von Zellen mit vier Symbolen. Gültig ist sie, wenn jede aktive Viererzelle in der Gruppe null ergibt. Die Gruppe ist entweder XOR (passt zum Quartikcode) oder Z₄ (gleiche Form, falsche Gruppe, als Gegenprobe). Verglichen werden stets dieselben Netze mit einem einzigen ausgetauschten Baustein.

Vorab die Strukturprüfung (`e0_structure.py`): Auf klassischen Zuständen ist Q_P exakt die XOR‑Prüfung (Diagonale = ¼·[x₁⊕x₂⊕x₃⊕x₄ = 0]), der Träger von P sind die 64 XOR‑Null‑Tupel, und der fünfte Codevektor ist „alle vier verschieden“. Die Schutzenergie H_mom hat dagegen auf gültigen und ungültigen Tupeln dieselben Werte (xxyy gültig 0,875, xxxy ungültig 0,875).

### KI1 · Weltmodell mit neuen Regeln, Gegenfakten und Planung

**(a) Vervollständigen verdeckter Welten** (`e1_three_way.py`, 4×4‑Welten, fünf Regelfamilien, 3 000 Trainingsschritte, 3 Startwerte). Anteil vollständig richtig gelöster Welten, XOR:

| Trainingswelten | Aufgabe | Transformer | lokales Netz | mit TFPT‑Projektor | mit Zufallsprojektor |
|---|---|---|---|---|---|
| 30 | bekannte Regelkombination | 0,08 | 0,52 | **0,92** | 0,58 |
| 30 | neue Regelkombination | 0,07 | 0,29 | **0,87** | 0,52 |
| 30 | mehr verdeckte Zellen | 0,00 | 0,01 | **0,27** | 0,02 |
| 300 | bekannte Regelkombination | 0,08 | 0,80 | **1,00** | 0,90 |
| 300 | neue Regelkombination | 0,07 | 0,58 | **0,99** | 0,85 |
| 300 | mehr verdeckte Zellen | 0,00 | 0,06 | **0,73** | 0,22 |
| 3 000 | bekannte Regelkombination | 0,08 | 0,42 ± 0,29 | 1,00 | 1,00 |
| 3 000 | neue Regelkombination | 0,07 | 0,34 ± 0,20 | 1,00 | 1,00 |
| 3 000 | mehr verdeckte Zellen | 0,00 | 0,03 | 0,91 | 0,92 |

In der Z₄‑Welt ist der Vorsprung kleiner (300 Welten, mehr verdeckte Zellen: 0,52 gegen 0,33 für den Zufallsprojektor und 0,09 für das lokale Netz), und bei 3 000 Welten ist der Zufallsprojektor auf neuen Regelkombinationen besser (0,99 gegen 0,92). Der Transformer erreicht auch mit 30 000 Welten und fünffachem Training nur 0,24.

**(b) Gleich informierte Kontrolle** (`ki1_generic_control.py`): dasselbe Netz mit der gewöhnlichen exakten Faktornachricht „vier Symbole ergeben null“.

| Welt | Trainingswelten | Aufgabe | TFPT‑Projektor | gewöhnliche Nachricht, richtige Gruppe | gewöhnliche Nachricht, falsche Gruppe |
|---|---|---|---|---|---|
| XOR | 30 | bekannte Regelkombination | 0,92 | **0,96** | 0,80 |
| XOR | 30 | neue Regelkombination | 0,87 | **0,93** | 0,62 |
| XOR | 30 | mehr verdeckte Zellen | 0,27 | **0,43** | 0,18 |
| XOR | 300 | mehr verdeckte Zellen | 0,73 | **0,84** | 0,49 |
| XOR | 3 000 | mehr verdeckte Zellen | 0,91 | **0,98** | 0,86 |
| Z₄ | 30 | neue Regelkombination | 0,63 | **0,92** | 0,63 |
| Z₄ | 300 | mehr verdeckte Zellen | 0,52 | **0,92** | 0,54 |
| Z₄ | 3 000 | mehr verdeckte Zellen | 0,83 | **0,99** | 0,84 |

Die gewöhnliche Nachricht mit der richtigen Regel ist in jeder Zeile besser als der TFPT‑Projektor, auch in der XOR‑Welt, zu der er passt. In der Z₄‑Welt verhält sich der TFPT‑Projektor wie die gewöhnliche Nachricht mit der falschen (XOR‑)Regel: 0,52 gegen 0,54. Das bestätigt die Strukturprüfung, der Projektor **ist** die XOR‑Regel.

**(c) Planung mit Eingriffen** (`ki1_planning.py`). Sichtbare Zellen sind Schalter, verdeckte sind Lampen (im Mittel 5,3). Eine Aktion setzt einen Schalter um, die Welt antwortet mit dem neuen Lampenzustand. Ziel: eine Lampe auf ein vorgegebenes Symbol bringen. Jedes Modell sagt für 8 Kandidaten die Folge voraus, derselbe Planer wählt, die exakte Welt urteilt. 1 000 Trainingswelten, 1 000 gepaarte Aufgaben je Gruppe, 3 Startwerte. Anteil erfolgreicher Pläne:

| Gruppe | Zufall | Transformer | lokales Netz | gewöhnliche Regelnachricht | TFPT‑Projektor | Zufallsprojektor | TFPT minus gewöhnlich (95 %‑Intervall) |
|---|---|---|---|---|---|---|---|
| XOR, neue Welten | 0,20 | 0,21 | 0,78 | 1,00 | 1,00 | 0,98 | 0,000 (0,000 bis 0,000) |
| XOR, neue Welten, zwei Schritte | 0,07 | 0,08 | 0,67 | 0,96 | 0,97 | 0,94 | +0,009 (0,000 bis +0,017) |
| XOR, neue Regelkombination | 0,19 | 0,20 | 0,53 | 0,92 | 0,93 | 0,80 | +0,008 (−0,002 bis +0,018) |
| XOR, neue Regelkombination, zwei Schritte | 0,06 | 0,08 | 0,31 | 0,82 | 0,84 | 0,66 | +0,016 (+0,002 bis +0,030) |
| Z₄, neue Welten | 0,20 | 0,20 | 0,97 | 1,00 | 1,00 | 1,00 | −0,001 (−0,002 bis 0,000) |
| Z₄, neue Welten, zwei Schritte | 0,07 | 0,09 | 0,95 | 0,96 | 0,96 | 0,97 | 0,000 (−0,006 bis +0,007) |
| Z₄, neue Regelkombination | 0,19 | 0,19 | 0,75 | 0,93 | 0,88 | 0,87 | −0,054 (−0,071 bis −0,038) |
| Z₄, neue Regelkombination, zwei Schritte | 0,06 | 0,06 | 0,61 | 0,84 | 0,77 | 0,76 | −0,067 (−0,086 bis −0,049) |

Regelverletzungen in den vorhergesagten Zuständen (XOR, neue Regelkombination): lokales Netz 95 %, Zufallsprojektor 38 %, TFPT 0,9 %, gewöhnliche Regelnachricht 0,4 %.

**(d) Lange Rollouts** (`e7_world_model_consistency.py`): Ein Weltmodell mit Fehlerrate 0,3 % je Schritt ist nach 50 Schritten nie mehr exakt. Mit Konsistenzreparatur nach jedem Schritt bleiben 93 % exakt, nach 200 Schritten 77 %. Die Reparatur mit der falschen Gruppe hilft nicht.

**Entscheidung.** Die im Paket gesetzte Schwelle war 5 Prozentpunkte Vorsprung vor dem besten gleich informierten Vergleich. In der XOR‑Welt liegt der Vorsprung des TFPT‑Projektors bei 0 bis 1,6 Punkten, in der Z₄‑Welt verliert er 5 bis 7 Punkte. **Ein eigener TFPT‑Gewinn ist damit nicht belegt.** Belegt ist der große Gewinn durch eingebautes Regelwissen überhaupt: 40 bis 52 Punkte vor dem lokalen Netz bei neuen Regelkombinationen. Bei genügend Daten lernt ein Netz mit beliebigem Zufallsprojektor dasselbe.

**Grenze.** Verkleinerte Fassung des Protokolls: 16 Zellen statt Graphen mit 8 bis 48 Knoten, Horizont 2 statt bis 64, keine Lagerwelt, keine kalibrierten Zustandsmengen, 3 statt 5 Startwerte, eine Datenstufe.

### KI2 · A₃‑Zellen und Lernen durch Invarianten

Der vorhandene Vorlauf `a3_pilot.py` wurde neu gerechnet und stimmt Ziffer für Ziffer mit dem abgelegten Ergebnis überein.

- Bei 16 Trainingspaaren und verschobenem Zyklusrauschen: rohe Merkmale 8,5·10⁻⁴, A₃‑Projektion 8,1·10⁻⁵, **allgemeine K₄‑Graphprojektion 8,1·10⁻⁵**, zufälliger Rang 3 0,86. A₃ und Graphprojektion sind dieselbe Abbildung.
- Gegenfall, in dem die Zyklen das Signal tragen: roh 4,6·10⁻⁴, A₃ 0,74. Der Projektor entfernt dann genau die gesuchte Information.

Der Zusatzversuch `e3_invariants.py` mit den fünf Quartik‑Invarianten zeigt beide Seiten: Passt das Ziel zur Symmetrie, erreichen die fünf Merkmale R² = 0,92 bei 1 000 Beispielen, rohe Merkmale liegen unter null. Passt es nicht, erreichen die fünf Merkmale R² ≈ 0.

**Entscheidung.** Ein eigener TFPT‑Beitrag müsste aus einer zusätzlichen typgetreuen Operation kommen. Die liegt nicht vor. Der Versuch mit S₄‑Augmentation und TFPT‑Phasenwirkung im Hauptversuch ist nicht gebaut.

### KI3 · Lokale Selbstorganisation ohne globale Attention

`e2_local_vs_global.py`: Gitter aus überlappenden 2×2‑Zellen, trainiert auf 6×6, geprüft auf 6×6, 12×12 und 16×16. „Kette“ heißt, dass die Schlussfolgerung das ganze Gitter durchqueren muss. Anteil gelöst, XOR:

| Modell | 6×6 zufällig | 6×6 Kette | 12×12 zufällig | 12×12 Kette | 16×16 zufällig | 16×16 Kette |
|---|---|---|---|---|---|---|
| Transformer | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 |
| lokales Netz | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 |
| mit Zufallsprojektor | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 | 0,00 |
| mit TFPT‑Projektor, 12 Runden | 0,94 | 1,00 | 0,66 | 0,06 | 0,60 | 0,00 |
| mit TFPT‑Projektor, 4n Runden | 0,97 | 1,00 | 0,77 | 0,61 | 0,70 | 0,30 |

Die mittlere Schlusstiefe der Kettenrätsel ist 6,8, 15,9 und 22,1 Runden. Mit 12 Runden kann die Ferninformation auf 12×12 gar nicht ankommen, das erklärt die 0,06.

Gleich informierte Kontrolle auf dem Gitter (4n Runden):

| Welt | Modell | 6×6 zufällig | 6×6 Kette | 12×12 zufällig | 12×12 Kette | 16×16 zufällig | 16×16 Kette |
|---|---|---|---|---|---|---|---|
| XOR | TFPT‑Projektor | 0,97 | 1,00 | 0,77 | 0,61 | 0,70 | 0,30 |
| XOR | gewöhnliche Nachricht, richtige Gruppe | 0,98 | 1,00 | **0,80** | **0,89** | **0,75** | **0,44** |
| XOR | gewöhnliche Nachricht, falsche Gruppe | 0,67 | 0,92 | 0,15 | 0,14 | 0,04 | 0,01 |
| Z₄ | TFPT‑Projektor | 0,86 | 0,95 | 0,43 | 0,73 | 0,20 | 0,55 |
| Z₄ | gewöhnliche Nachricht, richtige Gruppe | **0,98** | **1,00** | **0,83** | **0,95** | **0,79** | **0,60** |

Auch hier ist die gewöhnliche, gleich informierte Nachricht mindestens gleich gut und auf großen Gittern besser.

Verspätete Updates und Zustandsverlust (`ki3_async_local.py`, dieselben Netze, kein Nachtraining, 4n Runden):

| Welt | Modell | Aufgabe | synchron | 20 % der Zellen setzen je Runde aus | dasselbe bei gleicher Zahl von Updates | halber innerer Zustand einmal gelöscht |
|---|---|---|---|---|---|---|
| XOR | TFPT‑Projektor | 6×6 zufällig | 0,97 | 0,92 | 0,93 | 0,95 |
| XOR | TFPT‑Projektor | 12×12 zufällig | 0,77 | 0,71 | 0,72 | 0,73 |
| XOR | TFPT‑Projektor | 12×12 Kette | 0,61 | 0,35 | 0,35 | 0,33 |
| XOR | gewöhnliche Nachricht | 6×6 zufällig | 0,98 | 0,93 | 0,94 | 0,97 |
| XOR | gewöhnliche Nachricht | 12×12 zufällig | 0,80 | 0,68 | 0,69 | 0,80 |
| XOR | gewöhnliche Nachricht | 12×12 Kette | 0,89 | 0,48 | 0,50 | 0,46 |
| Z₄ | gewöhnliche Nachricht | 12×12 Kette | 0,95 | 0,58 | 0,59 | 0,63 |

Bei kurzen Schlussketten kosten verspätete Updates und ein einmaliger Zustandsverlust nur 3 bis 12 Punkte. Bei langen Ketten halbiert sich die Lösungsrate, und mehr Runden gleichen das nicht aus. Die Netze wurden nur synchron trainiert.

Kommunikation: Auf 12×12 liest das lokale Modell in 4n Runden rund 70 000 Nachbarwerte, globale Attention mit vier Schichten rund 83 000. Auf 6×6 ist die globale Variante sparsamer (5 200 gegen 7 200), auf 16×16 die lokale (173 000 gegen 262 000). Der Vorteil der lokalen Variante liegt also in der Skalierung, nicht in kleinen Welten.

**Entscheidung.** Lokale Regeln genügen, wenn die Zelloperation die Regel kennt, und sie übertragen auf größere Gitter, solange genug Runden laufen. Ohne Regelwissen lernt hier weder das lokale noch das globale Modell die Aufgabe im gegebenen Budget.

**Grenze.** Keine öffentlichen Maze‑, Sudoku‑ oder ARC‑Aufteilungen, kein Vergleich mit CRAM.

### KI4 · Quartikcode, Schutzenergie und Fixpunktsuche

**(a) Projektion.** `ki4c_ki7_projector_reflection.py`: 256 Eingangsgrößen, 40 Trainingsbeispiele, Ridge‑Regression, echter Projektor gegen Zufall, PCA, PLS und alle 256 Größen. Testfehler geteilt durch Varianz (1,0 = nicht besser als der Mittelwert):

| Daten | TFPT Rang 5 | Zufall Rang 5 | PCA Rang 5 | PLS Rang 5 | alle 256 |
|---|---|---|---|---|---|
| Signal in den TFPT‑Raum gelegt (Funktionskontrolle) | **0,15** | 1,12 | 1,01 | 0,92 | 0,91 |
| Signal in zufälligem Fünferraum | 1,08 | 1,09 | 1,00 | 0,90 | 0,89 |
| Signal in zufälligem Fünferraum von Sym⁴ | 0,98 | 1,06 | 1,01 | 0,91 | 0,90 |
| Handschrift‑Ziffern | 0,95 | 0,96 | 0,77 | 1,14 | **0,59** |

Der Gewinn aus dem Mini‑Test der Unterhaltung entsteht also nur, wenn die Daten passend gemacht werden. Verlustkorrektur mit maskierten Registern steht in Q2: dort gewinnt der echte Projektor gegen den Zufallsraum im selben Sym⁴.

**(b) Energie.** `e4_energy_reasoning.py` (ohne jeden gelernten Parameter, exakt gleichverteilter Start) und `ki4b_energy_candidates.py` (acht Zufallsstarts je Rätsel, gleiche Iteration):

| Energie | 4×4, 1–5 verdeckt | 4×4, 6–7 verdeckt | Gitter 6×6 | Gitter 24×24 Kette | niedrige Energie verbirgt falschen Zustand (4×4, 6–7 verdeckt) |
|---|---|---|---|---|---|
| gewöhnliche Nebenbedingung 1 − Pr[XOR = 0] | 1,00 | 1,00 | 1,00 | 1,00 | 0 % |
| generisch gelernte Energie | 0,99 | 0,97 | 0,39 | – | 0 % |
| I − P | 0,87 | 0,35 | 0,08 | – | 41 % |
| H_mom | 0,11 | 0,00 | 0,00 | – | 95 % |

Die ersten vier Spalten sind Lösungsraten (bei der gelernten Energie: bester von acht Kandidaten, Gitter mit Zufallsstart). Rätsel, die sich nicht durch lokale Weitergabe lösen lassen, löst keine der Energien. Bei Zufallsstart statt gleichverteiltem Start bricht auch die gewöhnliche Nebenbedingung auf Gittern ein (6×6 Kette: 4 %). Für Sudoku 4×4 löst die Energie „alle verschieden“ (der fünfte Codevektor) 100 %, H_mom 0 %. Die parameterfreie Dynamik aus P‑Kontraktionen allein löst auf Gittern kein einziges Rätsel, weil der symmetrische Anteil von P unterbestimmte Zellen zu gleichen Symbolen zieht.

**Entscheidung.** H_mom gewinnt nicht gegen I − P, sondern verliert deutlich. I − P verliert gegen die gewöhnliche Nebenbedingung. Nach der Regel des Pakets ist damit der Code der Kandidat und die Energie nicht. Niedrige Energie bei falscher Lösung zeigt, dass H_mom als Kostenfunktion am Ziel vorbeigeht.

### KI5 · Adaptive Rechenzeit

`e5_adaptive_compute.py`, auf den Netzen aus KI3, Höchstbudget 6n Runden. Netz mit TFPT‑Projektor, XOR:

| Aufgabe | feste Tiefe: gelöst / Runden | Stopp bei Stillstand: gelöst / Runden / falsche Stopps | Stopp mit Konsistenzprüfung: gelöst / Runden / falsche Stopps |
|---|---|---|---|
| 6×6 Kette | 1,00 / 36 | 1,00 / 8,7 / 0 % | 1,00 / 6,7 / 0 % |
| 12×12 zufällig | 0,77 / 72 | 0,75 / 14,5 / 23 % | 0,77 / 22,8 / 0 % |
| 12×12 Kette | 0,61 / 72 | 0,60 / 22,8 / 38 % | 0,61 / 38,7 / 0 % |
| 16×16 Kette | 0,30 / 96 | 0,29 / 41,5 / 60 % | 0,30 / 75,2 / 0 % |

Die Rundenzahl bis zur Lösung entspricht der Schlusstiefe (Verhältnis 0,98 bis 1,26).

**Entscheidung.** Der Stopp mit expliziter Konsistenzprüfung spart 46 bis 82 % der Runden ohne Verlust und ohne einen einzigen falschen Stopp, auf dem schwersten Block 22 %. Der Stopp allein nach Stillstand hält in 23 bis 60 % der Fälle bei einem falschen Ergebnis an und scheidet als Gültigkeitszertifikat aus. Die Konsistenzprüfung ist die gewöhnliche Regelprüfung. Für XOR ist sie identisch mit dem Stabilisatortest 1 − 4⟨Q_P⟩.

**Grenze.** Ein trainierter generischer Stopp und die gezielt gebauten Plateaus und Zweierzyklen fehlen. Angegeben sind Mittelwerte der Runden, keine p50‑ und p95‑Werte.

### KI6 · Wann hilft die E₈‑Symmetrie?

`ki6_e8_equivariance.py`: Eingaben sind Tupel von E₈‑Wurzeln. Fünf Lesarten: rohe Koordinaten, roh mit Weyl‑Augmentation, Gram‑Merkmale (alle Skalarprodukte, also allgemeine Geometrie), Normalform in der E₈‑Kammer, Normalform in der Untergruppe D₈. Trefferquote auf neuen Tupeln:

| Aufgabe | Beispiele | roh | roh + Augmentation | Gram | E₈‑Kammer | D₈‑Kammer |
|---|---|---|---|---|---|---|
| A: ist a + b − c eine Wurzel? | 30 | 0,51 | 1,00 | 0,98 | 0,98 | 0,86 |
| A, zurückgehaltene Orbitfamilien | 30 | 0,54 | 1,00 | 0,97 | 0,73 | 0,72 |
| A, zurückgehaltene Orbitfamilien | 3 000 | 1,00 | 1,00 | 1,00 | 0,67 | 0,53 |
| B: welche Sorte orthogonaler Vierergruppe? | 30 | 0,49 | 0,50 | 0,50 | **1,00** | 0,86 |
| B | 3 000 | 0,51 | 0,54 | 0,50 | **1,00** | 1,00 |

Aufgabe C bricht die Symmetrie durch eine feste Hintergrundrichtung h. Bleibt h verborgen, liegen Gram und E₈‑Kammer bei 0,50. Wird h mitgegeben, erreicht Gram 0,98 bei 30 Beispielen, die E₈‑Kammer 0,80, die Untergruppe, die h festhält (D₇), 0,93. Rohe Koordinaten erreichen 0,86. Die Invarianz unter allen acht einfachen Spiegelungen ist für Gram und E₈‑Kammer exakt (Abweichung 0), für rohe Netze 0,2.

Ergänzend `e6_e8_representations.py`: Als Quantisierer hat E₈ den kleinsten Fehler (zweites Moment 0,0717 gegen 0,0833 für das Würfelgitter) und hält 1 000 verrauschte Schritte mit Nachrasten zu 95 % aus (D₈ 82 %, Würfelgitter 6 %). In einem Autoencoder für Ziffern ist der Unterschied zu D₈ dagegen vernachlässigbar.

**Entscheidung.** E₈‑Äquivarianz bringt genau dort etwas, wo die Aufgabe eine Gittereigenschaft abfragt, die allgemeine Geometrie nicht sehen kann (Aufgabe B). Wo Skalarprodukte genügen, ist allgemeine Geometrie gleich gut und überträgt besser auf neue Orbitfamilien. Bei gebrochener Symmetrie ist die echte Untergruppe die bessere Wahl. Ein Transfer auf die Weltmodelle aus KI1 wurde nicht versucht, weil es dort keine geprüfte wirkungstreue E₈‑Darstellung gibt.

### KI7 · Reflexionsdynamik für lange Rechenketten

`ki4c_ki7_projector_reflection.py`: Ein Rechenschritt aus Spiegelungen U = I − 2dd† und einer Verschiebung, gegen ein tanh‑Netz und ein GRU. Aufgabe: ein Symbol nach T leeren Schritten wieder nennen (Zufall 12,5 %).

| Schritte | Reflexion + Verschiebung | tanh‑Netz | GRU |
|---|---|---|---|
| 20 | 100 % | 100 % | 100 % |
| 100 | 100 % | 63 % | 13 % |
| 300 | 100 % | 25 % | 13 % |

Der Gradient bleibt über 1 000 Schritte exakt erhalten, bei den Vergleichsnetzen verschwindet er.

**Entscheidung.** Die Eigenschaft ist echt und nützlich. Sie gehört zur bekannten Klasse orthogonaler und unitärer rekurrenter Netze, der TFPT‑Pulsoperator ist ein Vertreter davon. Die Vergleichsnetze hatten nur 600 Trainingsschritte und keine besondere Initialisierung.

## Signaturtests

### S1 · Kubisches Gitter am Himmel

Gerechnet wurde nur der phänomenologische Zweig: ein Himmelsmuster 1 + κ·K₄ mit der kubischen Harmonischen K₄ = n_x⁴ + n_y⁴ + n_z⁴ − 3/5, über dem bekannten Dipol (6,5 %, Rektaszension 100°, Deklination −24°), gesehen durch die Exposition eines Detektors auf −35,2° Breite. Die Gitterorientierung wird über 200 Lagen abgesucht, die Nullverteilung enthält diese Suche (entspricht etwa 48 unabhängigen Versuchen, Schwelle 4,1 σ für ein globales p von 0,001).

| Ereignisse | kleinstes κ mit 90 % Nachweiswahrscheinlichkeit | mit Aufteilung in Such‑ und Bestätigungshälfte |
|---|---|---|
| 2 000 | 0,79 | nicht erreichbar |
| 20 000 | 0,25 | 0,32 |
| 100 000 | 0,10 | 0,13 |

Die Grenze fällt wie 1/√N. Fehlt der Dipol im Nullmodell, steigt die Teststatistik ohne jedes kubische Signal nur leicht (Median 2,0 bis 2,6 σ, unter der Schwelle).

**Erkenntnis.** Mit den rund 50 000 Ereignissen oberhalb 8 EeV, die Auger insgesamt berichtet, wäre ein kubisches Muster ab etwa 15 % Amplitude auffindbar. Der öffentliche Datensatz ist kleiner und entsprechend unempfindlicher.

**Nicht gemacht.** Die echten Auger‑Ereignisse wurden nicht heruntergeladen, das braucht deine Zustimmung. Die Umrechnung von κ in eine Gitterweite fehlt, weil dafür eine feste Emissions‑ und Ausbreitungsannahme nötig ist. Natürliche Diskretheit und externe Implementierung wären in dieser Messung ohnehin nicht unterscheidbar.

### S2 · Botschaft in Naturkonstanten

- Für ein globales p unter 0,001 braucht ein vorab festgelegter Bitstring mindestens **10** unabhängige stabile Bits. Mit 16 ausprobierten Dekodern sind es 14, mit 1 000 Dekodern 20. Jedes tolerierte falsche Bit kostet weitere 4 bis 5.
- Ist die Botschaft echt und liegen die Residuen bei 2 σ, hat ein Test mit 12 Bits 70 % Nachweiswahrscheinlichkeit, mit 20 Bits 97 %. Bei 1 σ reichen auch 40 Bits kaum (40 %, mit zwei tolerierten Fehlern 80 %).
- Korrelierte Messungen entwerten Bits stark: zehn gleiche Vorzeichen treten bei Korrelation 0,3 in 3,6 % der Fälle zufällig auf, statt in 0,1 %.

**Erkenntnis.** Aus diesem Paket liegt genau ein reproduziertes Residuum vor (α). Sollstring und Dekoder gibt es nicht. Der eigentliche Botschaftstest ist damit **nicht ausführbar**, und die Zahl der dafür nötigen unabhängigen Konstanten ist weit größer als das, was TFPT derzeit eingefroren vorhersagt.

### S3 · Kausal getrennter Eingriff

Benötigte Blöcke je Gruppe für 90 % Nachweiswahrscheinlichkeit bei p < 0,001:

| Verschiebung δ in Standardabweichungen | 1 | 0,3 | 0,1 | 0,03 | 0,01 |
|---|---|---|---|---|---|
| Blöcke je Gruppe | 39 | 425 | 3 823 | 42 473 | 382 250 |

Simulation mit 400 Blöcken: Der Randomisierungstest hält seine Fehlalarmrate ein (4,2 % bei Soll 4,5 %), auch unter starker langsamer Drift (4,4 %). Eine Verschiebung von 0,5 Standardabweichungen wird ohne Drift in 96 % der Läufe gefunden, unter starker Drift nur in 1 %. Der Vergleich benachbarter Blöcke hebt das auf 12 %. Verspätet gezogene Placebo‑Bits schlagen nie an.

**Erkenntnis.** Der Test ist statistisch sauber, aber Drift ist der eigentliche Gegner. Ein Protokoll muss die Blockdauer deutlich kürzer wählen als die Driftzeit des Interferometers.

**Grenze.** Das ist reine Versuchsplanung. Die Unterhaltung nennt weder Zeitfenster noch Effektgröße, und TFPT sagt eine solche Abhängigkeit nicht vorher. Über die Natur wird hier nichts geprüft.

## Vorbehalte

- **Rechenumgebung.** Der Rechner war während der Läufe durch fremde Trainingsläufe stark ausgelastet. Das betrifft nur die Laufzeit.
- **Kleine Welten.** Die KI‑Welten haben 16 bis 256 Zellen. Aus dem Paket fehlen: die Lagerwelt mit Schlüsseln und Türen, Graphen mit 24 bis 48 Knoten, Horizonte 16 bis 64, öffentliche Maze‑, Sudoku‑ und ARC‑Aufteilungen, ein gelernter generischer Stopp, eine gelernte generische Energie.
- **Zahl der Wiederholungen.** Drei Trainingsstartwerte je Modell statt fünf.
- **Transformer‑Vergleich.** Der Transformer löst in allen Läufen höchstens 24 % der Welten, auch mit zehnfacher Datenmenge. Er ist hier ein schwacher Vergleich. Die tragende Aussage stützt sich auf die Vergleiche zwischen den lokalen Netzen.
- **Schwellen.** Wo ein Grenzwert erst nach Sicht der Ergebnisse gewählt wurde, ist die Aussage explorativ.

## Nachrechnen

```bash
cd experiments/tfpt-ai-architectures-20261118
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python torch numpy scikit-learn sympy mpmath scipy
.venv/bin/python f1_marks_dirac_weight.py     # jedes Skript läuft einzeln und schreibt nach out/
```

Die Reihenfolge für die Netzläufe ist `e1_three_way.py`, `e2_local_vs_global.py`, dann `e5_adaptive_compute.py`, `ki1_generic_control.py`, `ki3_async_local.py` (sie nutzen die gespeicherten Netze aus E2). `ki1_planning.py` und alle übrigen Skripte sind unabhängig.
