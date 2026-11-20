# KERNSIX20261120: Sechs Folgeschritte am originalen Randkern

Datum: 20. November 2026. Theorievertrag mit begrenzten Urteilen. Nichts wird befördert, kein
Ledger und kein Paper wird geändert. Ausgangspunkt ist P1BIND20261120
(`../p1-kernel-marker-cut-mirror-20261120`).

Reproduktion (69 Prüfungen, alle bestanden):

```bash
cd experiments/theory-contracts/kernel-followups-20261120
for s in 1_state 2_dictionary 3_time 4_marks 5_cubic 6_gravity; do ../../../.venv/bin/python -B check_$s.py; done
```

## Kurzfassung

| Schritt | Frage | Urteil |
|---|---|---|
| 1 | Liefert die Wirkung des Originals den Fockzustand von P_C? | bedingt ja: genau für das Profil f = −log\|x\| |
| 2 | Wörterbuch vom Kern zum E₈ Wort | Zustand, Form und Marker folgen; die Aufteilung 3 + 2 kommt nicht aus dem Kern |
| 3 | Zeit mit Länge | bedingt ja: Kragen der Breite ℓ gibt Temperatur 1/(2ℓ); ℓ ist frei |
| 4 | Lage der vier Marken | bedingt ja: auf der Spiegelfamilie ist das Quadrat das einzige Minimum |
| 5 | Kubik der Alpha Gleichung | nicht hergeleitet; Normalform und alle Koeffizienten zugeordnet |
| 6 | Raumzeit und Gravitation | statische Raumzeit R × Σ kinematisch; Gravitation nicht hergeleitet |

## Schritt 1: Randzustand aus der Wirkung (`check_1_state.py`, 12/12)

Quelle: `documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-…tex`,
Zeilen 302–318 (primitive Wirkung), Zeilen 319–343 (Osterwalder–Schrader Rekonstruktion),
`002-…tex` Zeilen 90–98 (Reflexionspositivität, ein elliptischer Kragenerzeuger) und 172–181
(Determinantenlinien der polarisierten Blöcke), `003-…tex` Zeile 12 (D_rel als „fermionic
comparison operator").

- Für den Kragenoperator D = γ_n(∂_n + B) laufen die Randwerte der Lösungen für jedes Innere
  gegen den positiven Spektralprojektor von B, mit der Rate 2 λ₁. λ₁ ist die Saatskala des
  Originals.
- Das fermionische Gaußintegral desselben Kragens gibt einen reinen Zustand, dessen
  Zweipunktfunktion dieser Projektor ist. Das gilt auch für einen Innenzustand, der nicht
  gaußsch ist.
- Die primitive Wirkung S = tr f(D_rel/χ) − tr f(D_ref/χ) + (iπ/2) Δη ist genau für
  f = −log|x| der negative Logarithmus der relativen Fermiondeterminante. Der Phasenterm
  (iπ/2) Δη ist schon für sich die Phase dieser Determinante. Für dieses Profil fällt χ heraus.
- Für ein anderes Profil (geprüft: Gauß) ist der Betragsteil keine Determinante.
- Die Determinantenlinie eines polarisierten Blocks ist seine Fock Vakuumlinie.

Urteil: Die Voraussetzung (a) aus P1BIND ist ersetzt durch eine schärfere: Das deklarierte Profil
der primitiven Wirkung ist f = −log|x|, die Wirkung also die eines Fermionfelds. Dann ist der
Randzustand der Fockzustand von P_C. Das Original lässt f offen („any declared spectral test
profile"). Für die fermionische Lesart sprechen im Original der Phasenterm, die
Determinantenlinien, der Name von D_rel und die Rekonstruktion aus Reflexionspositivität. Ein
Satz, der f festlegt, steht in den gelesenen Abschnitten nicht.

## Schritt 2: Wörterbuch zum E₈ Wort (`check_2_dictionary.py`, 16/16)

Quelle: `documents/archive--tfpt-45--source_extracts--02_carrier_source/sections/001-…tex`,
Zeilen 170–203 (Y = 1/12 + (5/12) ε_car, X = 6Y = −2 P₃ + 3 P₂).

- Die Hyperladung des Originals ist die spurfreie affine Funktion der Trägerpolarisation. 6Y ist
  der Markervektor (−2, −2, −2, 3, 3) des aufgezeichneten Worts.
- Die von G und ι_C erzeugte Kernalgebra ist eine einzige Kopie der 2 × 2 Matrizen. ι_C
  verschwindet auf jeder Normalhälfte, und auf jedem Block, der unter G und ι_C invariant ist,
  hat ι_C die Spur null. Eine Aufteilung 3 + 2 kann dort nicht liegen.
- Im Fockraum ist der Kernzustand auf einer Normalhälfte die normierte Spur. Ein Sprung von der
  einen Hälfte in die andere hat im Kernzustand das Gewicht 1/4 und wirkt bis auf ein Vorzeichen
  wie der Strom a_i⁺ a_j innerhalb einer Hälfte.
- Acht Nullmoden, flacher Zustand, Hyperladung des Originals: Markergewicht tr_E(Y²)/4 = 5/24,
  jeder der 60 Ströme zwischen Träger und Familie 1/4, Isospin 1/8 (Verhältnis 5/3), gemischt
  verbunden −1/32 und −1/72, drei Familien komprimiert 1/6. Das sind die aufgezeichneten Zahlen.
- Wörtliche Lesart (Trägermoden sind Eigenmoden von ι_C und werden vom Kernzustand besetzt):
  scharfe Ladung 1, Gewicht 1, keine Streuung. Das widerspricht der Aufzeichnung.

Urteil: Zustand (flach), Form (Spiegelform Tr(F*G)/2^d) und Marker (Hyperladung als Funktion
von ε_car) des E₈ Worts folgen aus Kern und Original. Die Zahl 1/4 ist das Gewicht eines
Sprungs zwischen den Hälften. Ausgeschlossen ist, dass die Trägermoden Eigenmoden von ι_C sind.
Im Wörterbuch wirkt ε_car auf den Moden einer Hälfte und vertauscht mit der Kernalgebra; die
Aufteilung 3 + 2 liefert der Trägerschritt, nicht der Kern. Nicht hergeleitet: warum genau die
Ströme vom Typ (10,6) der Schnitt sind, die Paarströme a a und a⁺ a⁺ als Kernoperatoren, und die
Auswahl von drei Familien.

## Schritt 3: Zeit aus dem Kragen (`check_3_time.py`, 10/10)

- Ein Kragen der Breite ℓ mit zwei Flächen trägt einen reinen Zustand. Auf eine Fläche
  eingeschränkt ist er C = (1 + exp(−2ℓB))⁻¹, der thermische Zustand von B bei inverser
  Temperatur 2ℓ. Im Fockraum ist die reduzierte Dichtematrix exp(∓2ℓ dΓ(B))/Z.
- Breite null gibt C = 1/2, den flachen Kernzustand mit Modularoperator 1 aus P1BIND.
- Nullmoden von B bleiben für jede Breite exakt flach. Dort lebt das E₈ Wort.
- Erster Hauptsatz dS = 2ℓ dE für jede kleine Änderung des Zustands.
- Der modulare Fluss exp(−i s 2ℓ B) ist die Dirac Entwicklung auf R × Σ. Auf dem Gitterring
  bleibt das Gewicht außerhalb des Lichtkegels unter 10⁻⁶.

Urteil: Zeit mit Länge entsteht bedingt aus der Abweichung, die P1BIND gefordert hat: Der
Kernzustand ist der Grenzfall Breite null, eine endliche Breite gibt die Temperatur 1/(2ℓ) und
den Fluss von B. Voraussetzungen: Schritt 1 und zwei Flächen des Kragens als die zwei Seiten.
Offen: Das Originaldatum legt keine Breite fest; in Saateinheiten ist T χ = 1/(2ℓ λ₁) frei.

## Schritt 4: Lage der vier Marken (`check_4_marks.py`, 10/10)

Modell: Vier Marken auf dem Nahtkreis sind die Verzweigungspunkte eines rechteckigen Torus
τ = i t. Wirkung S_f(t) = Σ f(λ/λ₁) über das Dirac Spektrum, geteilt durch den kleinsten
Eigenwert wie im Original.

- Antiperiodisch auf beiden Zyklen: S(t) = S(1/t) für jedes Profil; das Quadrat t = 1 ist das
  einzige Minimum, mit positiver Krümmung (vier Profile). Führender Term: Fläche mal λ₁² =
  π²(t + 1/t).
- Periodisch auf beiden Zyklen: Maximum mit Knick. Antiperiodisch auf einem Zyklus: kein
  stationärer Punkt am Quadrat.
- Verdrillung τ = x + i senkt die Wirkung nach beiden Seiten: Unter allen Tori ist das Quadrat
  ein Sattel.
- Die Vierteldrehung bildet die Impulse nur am Quadrat auf sich ab.
- Marken bei den vierten Einheitswurzeln (Doppelverhältnis 1/2) entsprechen genau t = 1.

Urteil: Unter drei Voraussetzungen sitzen die Marken bei μ₄: (i) die Marken liegen auf dem
Spiegelkreis (rechteckige Familie), (ii) Spinstruktur antiperiodisch auf beiden Zyklen, (iii) der
Modul wird durch das Minimum der mit λ₁ normierten Spektralwirkung des elliptischen Überlagerers
gewählt. (iii) ist ein Modell der „modularen freien Energie" des Originals, nicht das Original
selbst; die relative Wirkung mit D_ref und der η Term sind nicht nachgebildet. Ohne (i) ist das
Quadrat kein Minimum.

## Schritt 5: Kubik (`check_5_cubic.py`, 12/12)

Quelle: `documents/archive--tfpt-45--source_extracts--03_em_flavor_source/sections/001-…tex`,
Zeilen 124–215.

- Wurzel der eingefrorenen Gleichung: 1/α = 137,0359992168407125…
- Normalform in der inversen Kopplung w = c₃/α = 1/(2e²): w + κ w³ = 1/(2c₃²) = 32π² mit
  κ = (b₁/2π) log(1/φ). Das ist der stationäre Punkt des konvexen Quartiks
  w²/2 + κ w⁴/4 − w/(2c₃²); es gibt genau eine reelle Wurzel.
- Der logarithmische Term trägt 98,3 % der Bilanz. Ohne ihn wäre α = 2c₃³, 58-mal zu klein.
  Die führende Näherung w = (32π²/κ)^(1/3) gibt 1/α = 137,83.
- Koeffizienten: 4c₃ = 1/(2π). b₁ = (4/3) N_fam + (1/6) tr_{E₊}(Y²)/tr_E(Y²) = 4 + 1/10 mit der
  Hyperladung des Originals; die Summe von Y² über eine Familie ist 16 mal das Markergewicht 5/24
  aus Schritt 2. Ohne den Beitrag des Blocks E₊ (b₁ = 4) läge 1/α bei 138,16.
- Eine Laufrelation in einer Schleife wäre linear in w. Die Gleichung ist keine solche.

Urteil: Die Kubik ist nicht hergeleitet. Neu sind die Normalform und die Zuordnung: Alle
Koeffizienten außer den Potenzen (w¹, w³) und dem Quellterm 1/(2c₃²) sind Kern- und
Trägerdaten. Das Original nennt den quartischen Teil einen lokalen Gegenterm. Offen ist ein
Funktional der inversen Kopplung mit genau diesem quadratischen und quartischen Teil.

## Schritt 6: Raumzeit und Gravitation (`check_6_gravity.py`, 9/9)

- Die Entropie des Flächenzustands aus Schritt 3 ist endlich und proportional zur Fläche der
  Naht: S = (9ζ(3)/2) c₃ A/ℓ², geprüft auf dem flachen Torus und der runden Sphäre.
- Als S = A/(4G) gelesen: G = ℓ²/(18 ζ(3) c₃) = 1,1616 ℓ². Die Breite ist frei, also folgt kein
  Wert.
- Breite null: log 2 je Mode, die Modenzahl 4c₃ A L² wächst mit der Abschneideskala L.
- Der modulare Erzeuger vertauscht mit den Translationen der Naht: ein Gleichgewicht bei einer
  Temperatur. Er ist kein Boost, also gibt es keinen Horizont und keine Unruh Lesart.

Urteil: Kinematisch folgt aus Schritt 3 eine statische Raumzeit R × Σ mit Lichtkegel, bedingt
wie dort. Gravitation ist nicht hergeleitet. Von den Zutaten einer Herleitung nach Jacobson
liegen vor: der erste Hauptsatz (exakt) und ein Flächengesetz ohne Abschneideskala. Es fehlen:
eine von Ort zu Ort veränderliche Breite (Lapse), eine Horizontlesart des Flusses und ein Wert
der Breite.

## Was sich an F1 ändert

F1 bleibt in Bearbeitung. Die eine Voraussetzung aus P1BIND ist auf eine Aussage über das Profil
der Originalwirkung zurückgeführt. Neu offen und benannt: die Breite ℓ des Kragens (Schritte 3
und 6), die Wahl des Funktionals für den Modul (Schritt 4), das Funktional der Kubik (Schritt 5),
die Herkunft der Aufteilung 3 + 2 auf einer Normalhälfte (Schritt 2).
