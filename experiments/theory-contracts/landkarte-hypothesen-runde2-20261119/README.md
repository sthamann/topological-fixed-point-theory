# HYPCHECKB20261119 · Prüfung der zehn Hypothesen der zweiten Runde

Datum: 19. November 2026. Kein Ledgerwechsel, keine Beförderung, kein Atlaseintrag, kein Commit.
`physical_solution = false`. F1 bleibt offen.

## Frage

Die zweite Runde der Landkarte (`tfpt_landkarte/content-hypothesen2.js`, Nummern 11 bis 20) entstand aus den
Urteilen der ersten Runde und aus einer parallelen Prüfung in einer anderen Sitzung. Hier wird jede an ihrem
eigenen Test gerechnet. Wo die parallele Sitzung Daten benutzt, die nicht im Repository liegen, steht das als Grenze.

## Ergebnis in einem Satz

Keine der zehn löst F1. Vier fallen (11 als Zeit, 14, 17, 19), eine ist richtig, steht aber schon im Original (13),
drei tragen exakt in ihrem Kern und bleiben ohne Bindung an das Original (12, 15, 16), eine trägt nur halb (20),
eine ist ohne die fremden Daten nicht entscheidbar (18).

## Urteile

| Nr | Hypothese | Urteil | Prüfungen |
|---|---|---|---|
| 11 | Zeit zwischen Viererschnitt und Fünferschnitt | Geometrie exakt, als Zeit gefallen | 12/12 |
| 12 | Marken als Ausnahmestellen der Familie der Stufe 4 | Trägt als Mathematik, Bindung an P1 offen | 16/16 |
| 13 | Die Vier aus Torus und Fouriertransformation | Richtig, steht aber schon im Original | 11/11 |
| 14 | Raten als Anteile der Spinsektoren | Gefallen | 10/10 |
| 15 | Viertelphase als Abstand der Seiten | Kern exakt, Bindung an v724 gefallen | 11/11 |
| 16 | Sechs Paare geben quadratisch 3 + 1 | Form exakt, für fermionische Kanäle gefallen | 13/13 |
| 17 | CP Verletzung als ungerader Teil | Gefallen (2,3 σ bei der Phase) | 10/10 |
| 18 | TFPT als Bruchteil einer OPH Zelle | Nicht entscheidbar, Daten fehlen | 10/10 |
| 19 | Zeitpfeil als Umlauf ohne Rückkehr | Gefallen | 9/9 |
| 20 | Normierung aus dem Verkleben mit Drehung | Halb: c festgelegt, aber auf null | 9/9 |

### 11 · Nachbarknoten
- Exakt: v = 4ω₅ − 5ω₆ hat Norm 20. U(t) = exp(2πi t v) ist ein Kreis in E₈; U(1/4) ist der Erzeuger der
  Viererklebung (Zentralisator D₅ ⊕ A₃), U(1/5) der der Fünferklebung (A₄ ⊕ A₄), allgemein A₄ ⊕ A₃. Takt 20.
- Abbruchregel greift: v liegt ganz in D₅ und ist die bekannte Ladung X von SO(10) → SU(5) × U(1)
  (−5, −1, 3 auf dem Spinor). Auf den fünf Plätzen wirkt U nur als Zahl; weder Vierteluhr noch Fünfertakt entstehen.
- Als Zeit gefallen: periodisch, ganzzahliges Spektrum.

### 12 · Familie der Stufe 4
- Exakt: y² = x(x + t²)(x + 1) hat Fasern I₄, I₄ (bei 0, ∞) und I₂, I₂ (bei ±1), Doppelverhältnis −1, Punkt der Ordnung 4.
- Kern: Die Teilung durch den Punkt der Ordnung 2 gibt die Kurve bei m(t) = (1 + t)/(1 − t). m dreht die vier
  Stellen im Kreis. Die Vierteldrehung der Marken ist also eine Halbierung der Kurve (Isogenie vom Grad 2),
  bis auf die Verdrehung mit −1 (isomorph erst mit i). Die Halbdrehung t → −1/t ist eine echte Symmetrie.
- Exakt: Die Fixpunkte ±i von m tragen den quadratischen Torus (j = 1728). Das sind die zwei Pole neben den Marken.
- Exakt: A₃ ⊕ A₃ ⊕ A₁ ⊕ A₁ liegt in D₅ ⊕ A₃ und hat in E₈ Quotient Z₄ × Z₂, dieselbe Gruppe wie die Torsion der Fläche.
- Abbruchregel greift nicht. Nicht gezeigt: dass P1 die Basis ist; keine Abbildung von den fünf Trägerfunktionen
  auf die Schnitte. Spannung: im Original sind die vier Marken gleichartig.

### 13 · Die Vier
- Exakt: Geschlecht (n − 2)/2; flache Kugeln mit Kegelpunkten gibt es genau vier, (2,2,2,2), (2,3,6), (2,4,4), (3,3,3).
  Nur (2,2,2,2) hat eine freie Lage.
- Das Original leitet die Vier schon so her (flache Z₂ Klasse, Gauss Bonnet, v216; läuft durch) und nennt die
  quadratische Lage ausdrücklich zusätzlich. Flach, Torus und Fouriertransformation sind dieselbe Bedingung.
- Zirkel beim Teil mit Schieber und Uhr: die kommen im Original aus der Klebegruppe, also von P2.

### 14 · Thetaanteile
- Exakt: 1/6, 2/3, 1/6 am Quadrat. Am Rechteck 2:1: 0,00045 · 0,51471 · 0,48485. Am Sechseck ist die Summe null.
- Abbruchregel greift: Die Vielfachheiten 6 und 6 folgen nicht, und die Raten des Originals haben keinen Formparameter.
- Nachgerechnet (parallele Aussage): gemischte Ladungsantwort senkrechter Wurzeln ist in E₈ immer 1/3, in D₁₆ nur
  1/9 oder 7/9. Erklärung: 1/3 ist die Zahl jeder gleichmäßigen Verteilung (Σx²y² : Σx⁴ = 1 : 3).

### 15 · Viertelphase
- Satz, für alle vier Z₄ Formen geprüft: Sind zwei Seiten verklebbar, ist das Verhältnis ihrer Gaußphasen gleich dem
  Kommutator von Schieber und Uhr der ersten Seite. Für D₅ und A₃: i. Seitentausch gibt −i.
- v724 ist reell. Der umlaufabhängige Teil der Dreierspur ist bei allen vier Tripeln exakt null.
- Ein Ordnungssignal Tr(M PₐP_bP_c) − Tr(M PₐP_cP_b) gibt es bei allen vier Tripeln, aber mit verschiedenen
  Beträgen (0,082 bis 0,157). Abbruchregel greift für diese Bindung.
- Grenze: Die Größe der parallelen Sitzung liegt nicht im Repository; geprüft ist die hier definierte.

### 16 · Paarform
- Exakt: Σ_{i<j} x_i x_j hat Signatur (1, n − 1), also 3 + 1 genau bei vier Größen. Ganzzahlig U ⊕ A₂(−1),
  Determinante −3. Die vier Achsen sind lichtartig. Der Viererzyklus ist eine Drehspiegelung der Ordnung 4.
- Nachgerechnet (parallele Aussage): Kanten = A₁ + E + T₂, antisymmetrische Paare = T₁ + T₂, linearer Rang höchstens 3.
- Quadratisch mit S₄ allein: Rang 10 erreichbar. Mit voller Familiensymmetrie SU(4): Rang null (Abbruchregel greift
  für die Lesart mit Ladungen). Die sechs fermionischen Kanäle tragen die Keilform mit Signatur (3,3).
- Nicht gezeigt: gleiche Gewichte im Original; aus einer Form allein folgt kein Ort.

### 17 · CP
- Nachgerechnet: J₊ = 3,32702 · 10⁻⁵, J₋ = −2,78956 · 10⁻⁵ aus der Originalformel (v88, v561). Bei der Umkehr kippt
  nur π/3, der Zusatz 3λ² nicht.
- Ungerader Teil 3,0583 · 10⁻⁵: −0,16 σ gegen J = 3,08(14) · 10⁻⁵ (voll: +1,76 σ).
- Zweittest Phase gegen γ = 65,7(3,0)°: voll +0,98 σ, ungerade Lesart −2,27 σ. Zusammen schlechter (5,2 gegen 4,1).
- Der gerade Teil ist messbar (|V_td| unterscheidet die Orientierungen). Das Original erklärt den Überschuss bei J
  schon über die Skala von s₂₃ (dann +0,34 σ).

### 18 · OPH Zelle
- Rekonstruktion: Q = n · Σr²/(Σr)² gibt für die Raten 1/3, 1/3, 2/3 genau 9/8. Das ist eine Lesart, keine
  hinterlegte Definition. Die OPH Werte sind exakt 8 − 96/N.
- Ersatzrechnung am Wurzelgraphen von E₈ (Spektrum 0, 28, 48, 58, 60): Q = 1,014, nicht 7,6.
- Korrektur: Die Spiegelungen von D₅ und A₃ geben fünf Klassen (40, 12, 64, 64, 60), nicht vier. Zusammengefasst
  Q = 1,008, nicht 9/8.
- Nicht entscheidbar ohne das OPH Netz und die Definition des Quotienten.

### 19 · Zeitpfeil
- Gefallen: τ und τ + 4n sind derselbe Torus mit demselben Spektrum.
- Gefallen: Auf Thetafunktionen der Stufe N hat T Ordnung 2N; die Zählung wird wieder zur Uhr.
- Abbruchregel greift: Das flache Kissen hat eine feste Lücke 4π² je Flächeneinheit.
- Bleibt: Ein Zeitpfeil braucht einen wachsenden Zustandsraum.

### 20 · Normierung
- Mit festem Umfang bricht exp(−c · Fläche) die Invarianz unter der Vierteldrehung für jedes c ≠ 0. Genau ein c bleibt.
- Der Wert ist null, nicht c₃. c₃ ist die Balance am Fixpunkt (+1/2 und −1/2, geteilt durch 4π).
- Mit fester Fläche ist jedes c erlaubt (das ist der parallele Befund). Ohne absolute Länge hat c keinen Platz.

## Neue exakte Befunde

1. Ein Kreis in E₈ läuft durch die Klebeerzeuger des Viererschnitts und des Fünferschnitts; er ist die Ladung X von SO(10), Takt 20.
2. Auf der Familie der Stufe 4 ist die Vierteldrehung der vier Ausnahmestellen eine Halbierung der Kurve; ihre Fixpunkte tragen den quadratischen Torus.
3. A₃ ⊕ A₃ ⊕ A₁ ⊕ A₁ hat in E₈ den Quotienten Z₄ × Z₂.
4. Verhältnis der Gaußphasen zweier verklebbarer Z₄ Seiten = Kommutator von Schieber und Uhr.
5. Die Paarsumme auf vier Größen ist die Lorentzform U ⊕ A₂(−1).
6. Flache Kugeln: nur der Fall mit lauter Halbdrehungen hat eine freie Markenlage.

## Nachgerechnete Aussagen der parallelen Sitzung

- Gitterprobe: 0 oder 1/3 gegen 1/9 oder 7/9. Stimmt.
- Sechs Paarübergänge: linearer Rang höchstens 3. Stimmt (als Darstellungen von S₄).
- J₊ und J₋. Stimmen.
- Ordnungssignal in v724: vorhanden in der hier definierten Größe.
- Nicht nachrechenbar: OPH Quotienten, Energiekriterium des Vakuumwechsels.

## Wiederholen

```
cd experiments/theory-contracts/landkarte-hypothesen-runde2-20261119
for s in h11_nachbar h12_beauville h13_genus h14_theta h15_phase h16_paarform h17_cp h18_zelle h19_pfeil h20_norm; do python3 -B check_$s.py; done
```
