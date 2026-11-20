# HYPCHECK20261119 · Prüfung der zehn Hypothesen der Landkarte

Datum: 19. November 2026. Kein Ledgerwechsel, keine Beförderung, kein Commit. Fünf Befunde nachträglich im Atlas eingetragen (siehe contract_index.json).
`physical_solution = false`. F1 bleibt offen.

## Frage

Die Landkarte (`tfpt_landkarte/content-hypothesen.js`) nennt zehn Hypothesen, wie TFPT fundamental gelöst
werden könnte. Jede hatte einen Test und eine Abbruchregel. Hier wird jede so weit gerechnet, wie es ohne
neue Annahmen geht. Wo die Bindung an das Original am gesperrten Schnitt R113 hängt, steht das als Grenze.

## Ergebnis in einem Satz

Keine der zehn löst F1. Fünf fallen in ihrer wörtlichen Form (4, 5, 6, 8, 10), bei einer fällt die genannte
Formel und eine exakte tritt an ihre Stelle (7), vier tragen exakt in einer engeren Fassung (1, 2, 3, 9).
Sechs exakte Befunde bleiben, die vorher nicht im Korpus standen.

## Urteile

| Nr | Hypothese | Urteil | Prüfungen |
|---|---|---|---|
| 1 | Marken als Hauptrichtungen eines Krümmungsspinors | Enge Fassung exakt, starke Fassung ausgeschlossen | 10/10 |
| 2 | Vier Qubits, Marken als GHZ Zustand | Auswahl der Lage exakt, Bindung an das Original nicht gezeigt | 7/7 |
| 3 | Ein Flussquant auf der Naht | Rechnung exakt, Gleichsetzung mit der Windung des Originals nicht gezeigt | 7/7 |
| 4 | Klebung als Vierteldrehung plus halbe Verschiebung | Wörtlich gefallen, exakter Ersatz gefunden | 11/11 |
| 5 | P1 und P2 als Faser und Umlauf einer Funktion | Wörtlich gefallen, exakter Ersatz gefunden, Zeitlesart gefallen | 14/14 |
| 6 | Verklebung ist Verschränkung | Statik exakt, Dynamik gefallen, eigene Abbruchregel greift | 11/11 |
| 7 | Schnitt als Knoten des E₈ Diagramms | Genannte Formel gefallen, exakte gemeinsame Regel gefunden | 13/13 |
| 8 | Wachstum heißt ein Qubit anfügen | Gefallen, eigene Abbruchregel greift | 13/13 |
| 9 | Quelle dort, wo alle Schranken scharf sind | Teilweise: Randpunkte ja, Eindeutigkeit nein | 12/12 |
| 10 | Frobenius Elemente als Ereignisprotokoll | Gefallen, eigene Abbruchregel greift | 12/12 |

### 1 · Krümmungsspinor
- Exakt: Die Eichalgebra des Standardmodells (Dimension 12) ist der Zentralisator der Halbdrehung exp(iπJz)
  der Hauptdrehgruppe SL(2) auf C⁵. Andere Drehwinkel geben Dimension 8, 6, 4, 4. Y = 1/12 − (5/12)R.
- Exakt: z⁴ − 1 hat I = −1, J = 0. Möbius Stabilisator: allgemein 4, Quadrat 8, Tetraeder 12.
- Ausgeschlossen: dieselbe Fünf zugleich als Spin 2 der Raumzeit mit vertauschender Eichwirkung. Die
  Darstellung ist irreduzibel, nur Vielfache der Eins vertauschen mit der ganzen Gruppe.
- Offen: ob der originale Kern außer μ₄ Möbiuskovarianz trägt.

### 2 · GHZ
- Exakt: 3 Qubits: 1080 Stabilisatorzustände, 18 symmetrisch, 12 mit drei verschiedenen Sternen, alle
  gleichseitig. 4 Qubits: 36720, 18 symmetrisch, 12 mit vier verschiedenen Sternen, alle Quadrate.
- Grenze: Das Prinzip wählt das regelmäßige n Eck, nicht die Vier.
- Typgrenze: Die Fünf des Originals ist Fix_P(Sym⁴ C⁴), die Majorana Lesart braucht Sym⁴ C². Als
  Darstellungen von A₅ sind beide gleich (Charaktere 5, 1, −1, 0, 0). Die Gleichsetzung verlangt die Wahl
  einer Ikosaedergruppe in S₆, die im Original nicht steht.

### 3 · Flussquant
- Exakt: Ein Flussquant durch die Einheitszelle gibt ℓ² = 1/(2π) und Varianz 1/(4π) = 2c₃ je Quadratur. Das
  ist die Konvention h = 1, jetzt als Folge der Quantisierung.
- Exakt: Das zurückgezogene Trägerbündel hat auf dem Torus Grad 8 mit acht Schnitten = 5 gerade + 3
  ungerade; die ungeraden verschwinden an den Marken.
- Grenze: Fluss 1 auf dem Torus ist kein Rückzug von der Kugel (halber Fluss). Quadrat und Rechteck 2:1
  werden nicht getrennt.
- Nicht gezeigt: Die Einheitswindung des Originals ist die Chernzahl des Endlinienbündels über C/K, nicht
  der Fluss durch den Markentorus.

### 4 · Affine Vierteldrehung
- Exakt: z → iz + ½ ist ein Viererzyklus auf den Marken und steigt zur Kugel ab.
- Gefallen: Diese Markendrehung ist das Coxeter Element von A₃ und wirkt auf der Klebegruppe trivial. Die
  Klebegruppe selbst wirkt als Skalar i. Der Klebegenerator ist der Kommutator von Schieber und Uhr.
- Ersatz, exakt: Der Coxeter Viererzyklus ist genau dann eine Möbiussymmetrie, wenn die Marken harmonisch
  liegen. Marken plus Pole haben 24 Möbiussymmetrien; die Bolzakurve hat 48 = GL(2,3) mit Zentrum der
  Ordnung 2 und ist damit nicht Picards G48 (Zentrum μ₄).
- Offen: Pole im Original (R113 gesperrt).

### 5 · Faser und Umlauf
- Exakt: Faser von z⁴ ist A₃ mit Umlauf der Ordnung 4. A₂ ⊗ A₄ gibt E₈, Umlauf der Ordnung 30 = c₂ ⊗ c₃ ⊗ c₅.
- Gefallen: D₅ ⊕ A₃ ist kein Teildiagramm von E₈, entsteht also durch keine Deformation der Singularität.
- Ersatz, exakt: die elliptische Fläche mit Punkt der Ordnung 4, y² + xy − ty = x³ − tx², Δ = t⁴(16t + 1),
  Fasern I₄, I₁, I₁*. Die Torsion der Ordnung 4 ist der Klebevektor. Drei Fasern haben j = 1728.
- Gefallen: Zeit als Umlauf. Alle drei Umläufe (T⁴, T, −T) haben unendliche Ordnung.
- 3/4 + 1/2 = 5/4 bleibt Arithmetik ohne Abbildung.

### 6 · Verschränkung
- Exakt: Sektorgewichte, konjugierte Fouriermatrizen, c = 5 und 3, Klebevektor maximal verschränkt und fest
  unter S ⊗ S und T ⊗ T.
- Nicht eindeutig: Der feste Raum ist zweidimensional (diagonal und antidiagonal).
- Gefallen: Dynamik. Die Wirkung ist ein Produkt, unitär, endliche Gruppe (projektiv 192), ohne Raten.
- Abbruchregel greift: reduzierter Zustand 1/4 auf vier Sektoren, die Referenz ist I/5 auf fünf Plätzen.

### 7 · Knoten des E₈ Diagramms
- Gefallen: k(p − k)/(2p) als Randinvariante.
- Exakt, p = 2 bis 12, alle k: η_k(p) = (p² − 1)/(6p) − k(p − k)/p. Der zweite Term ist die Klebenorm der
  Klasse k von A_{p−1}. Die zwei Muster aus P1KERN sind zwei Fälle dieser einen Formel.
- Exakt: η verschwindet bei Einheitswindung k = 1 nur für p = 5.
- Exakt: Knotentafel des erweiterten Diagramms (Marke 4: 5/4 + 3/4; Marke 5: 4/5 + 6/5; Marke 6 dreiteilig).
- Grenze: 5/8 der D₅ Seite ist nur bei p = 4 eine Linsenraumgröße. Nicht gezeigt: P1 sitzt an einem solchen Schnitt.

### 8 · Turm
- Exakt: Hamming Code → E₈, φ = 1 + J verdoppelt Normen, Stufe 16 ohne Wurzeln mit 4320 kürzesten Vektoren.
- Abbruchregel greift: J invariante Unterräume haben gerade Dimension, D₅ ⊕ A₃ (5 und 3) kann nicht erhalten bleiben.
- Gefallen: Raten. Alle singulären Werte sind √2.
- Nicht entscheidbar hier: Ort als Tensorfaktor.

### 9 · Schärfe
- Exakt: Die Sprungraten der S₃ Familie sind nichtnegativ genau für 2/9 ≤ q ≤ 1/4 (Reproduktion, zusätzlich
  am Matrixlogarithmus auf zwei Registern gegengeprüft).
- Neu: Bei q = 2/9 springen nur die drei Spiegelungen, und q = (2/3)(1/3). Bei q = 1/4 fällt die Spiegelung
  r59 weg. Das Prinzip lässt zwei Punkte übrig, nicht einen.
- Exakt: γ = m wird im Viererring erreicht.
- Nicht verbunden: keine Abbildung von γ = m auf q. Für κ und das Seitenverhältnis gibt es keine Schranke.

### 10 · Frobenius
- Exakt: |a_p| = 2a mit p = a² + b² (a ungerade), 0 für p ≡ 3 mod 4. Werte 2, 6, 10, 22 an 5, 13, 41, 137.
- Vergleich mit 1000 Zufallskurven (Festlegung vorher, seed 20261119): mindestens drei Treffer aus vier
  erreicht etwa die Hälfte. Auf 166 Primzahlen hat die Kurve weniger Treffer als fast alle Zufallskurven.
- Bleibt: Frobenius in (Z/5)* = Z/4 als Etikett je Primzahl. Ohne Reihenfolge, Rate oder Verknüpfung.

## Neue exakte Befunde

1. Eichalgebra = Zentralisator der Halbdrehung der Hauptdrehgruppe auf C⁵.
2. Gemeinsame Schnittregel η_k(p) = (p² − 1)/(6p) − k(p − k)/p; Nullstelle bei k = 1 nur für p = 5.
3. Elliptische Fläche mit Punkt der Ordnung 4: Fasern I₄, I₁, I₁*, Torsion = Klebevektor.
4. Coxeter Viererzyklus der Marken ist Möbiussymmetrie genau in harmonischer Lage; Klebegenerator = Kommutator.
5. Acht Thetafunktionen auf dem Torus = 5 gerade + 3 ungerade.
6. q = 2/9 ist der einzige Punkt der Familie mit reinen Spiegelungssprüngen.

Alle sechs sind Mathematik auf vorhandenen Objekten. Keiner bindet P1 an P2, keiner liefert Dynamik.

## Grenzen

- Die Prüfungen arbeiten mit den im Atlas registrierten Objekten, nicht mit dem gewichteten Originalschnitt R113.
- „Gefallen“ gilt für die hier formulierte Fassung, nicht für jede denkbare Abwandlung.
- Die Typgrenze bei 2 und die Flussgleichsetzung bei 3 sind offene Fragen an das Original, keine Ausschlüsse.

## Wiederholen

```
cd experiments/theory-contracts/landkarte-hypothesen-20261119
for s in h01_weyl h02_ghz h03_fluss h04_affin h05_milnor h06_bell h07_kac h08_turm h09_scharf h10_frobenius; do python3 -B check_$s.py; done
```
