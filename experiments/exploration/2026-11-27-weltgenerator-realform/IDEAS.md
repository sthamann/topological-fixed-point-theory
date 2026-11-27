# Gesamtszenarien aus den Grundlagen: Weltgenerator und Realform (27. November 2026)

Modus: Exploration. Auftrag: zurück zu den Grundlagen und zu dem, was funktioniert (`docs/research/TFPT-STAND.md`,
Abschnitte 2 bis 5), ganze Szenarien dafür, wie alles zusammenhängt, und eine Suche, die vollständige Lösungen
enthalten und erkennen kann. Keine Prämisse des Kandidaten H6 ist Startpunkt.

## Befund vor den Ideen

Der Atlas enthält fast jede einzelne Lesart bereits: Reinheit zweier Bögen (18.11.), Jones-Index 4 in
c₃ (v813), G₃₁ und Gauß-Gitter, Zwei-Qubit-Clifford, Twistorgerade, Himmelskugel und Weyl-Spinor (Z61). Es fehlt kein
lokaler Baustein. Es fehlte eine Suche, die ganze Welten erzeugt und an allen Zielen zugleich misst, und es fehlte
die Frage, welche Lesart welche Rolle aus ihrem eigenen Grundobjekt liefert. Nicht im Bestand gefunden (Atlas- und
Repositoriumssuche): die Realform E₈(−24), die geschlossene 2+1-Welt mit Kegelteilchen und die Aussage, dass
„dieselbe Vier“ nur für zwei Blätter gilt.

## Diversitätsprüfung

Vorhandene Ordner unter `experiments/exploration/`: Fibonacci-Uhr (Anyonen), fremde Ursprünge (Graphen, Automaten,
Polytope), kleinstes Verbundsystem (zwei Qubits), operative Ursprünge, Quaternionrahmen, Vertauschung und Auslese,
Fusionsspur, Fockspur und GSO, Phasenraum F₅ und Holocode (H6), Quelldynamik, unabhängige Ursprungsprinzipien
(Wechselwirkung, Quadrikenbüschel, kausale Signale, Singularitäten). Neu als Gebiete: Realformen von Lie-Algebren
(Ansatz 1), Gravitation in 2+1 Dimensionen (Ansatz 3), kritischer Index der Operatoralgebren zusammen mit kritischer
Statistik (Ansatz 4).

## Ansätze

Jeder Ansatz nennt ein Grundobjekt, was daraus folgen soll, die Messgröße, die freien Parameter und die erste freie
Wahl. In der Suche heißen sie Lesarten (`READINGS` in `research/originlab/whole_solution_search.py`).

1. **Realform (Lie-Theorie): Raumzeit ist die nichtkompakte Form des Familienfaktors.**
   Objekt: das geschlossene Gitter D₅ ⊕ A₃ in der reellen Form mit kompakter Eichgruppe und nichtkompaktem
   Familienfaktor. Folgen soll: Signatur 3+1 (so(4,2) ist die konforme Algebra von R^{3,1}), Kausalstruktur, die
   Leiter mit positiver kontinuierlicher Zeit, die Uhr μ₄ als Zentrum von SU(2,2) mit Quadrat (−1)^F. Messgröße:
   Dimension und Signatur der Raumzeit. Freie Parameter: keine. Erste freie Wahl: dass die physische Realform nicht
   die kompakte ist, und dass die 4 von A₃ das Nahtdublett mit seinem Tomita-Partner ist.
2. **Spin 2 auf der Himmelskugel (Spinorgeometrie).** Objekt: ein Spin-2-Zustand auf der Richtungskugel eines
   Beobachters. Folgen soll: d = 4 (nur S² ist eine Riemannsche Fläche), zwei Blätter (π₁ der Drehgruppe), vier
   Marken (Grad 2s), Träger 5 = 2s + 1. Messgröße: 3+1, Teilcheninhalt. Erste freie Wahl: setzt eine lorentzsche
   Mannigfaltigkeit mit Gravitationsfeld voraus; τ = i (J = 0) bleibt unerklärt. Bestand: Z61, Z62.
3. **Geschlossene 2+1-Welt (Gravitation).** Objekt: geschlossene Raumscheibe der 2+1-Gravitation mit gleichen
   Kegelteilchen. Folgen soll: flache Sphäre (die Zwangsbedingung ist Gauß–Bonnet), c₃ = 1/(8π) als Masse je
   Defizitwinkel, Marken als Teilchen der Spannung 1/8. Messgröße: Flächennormierung G. Erste freie Wahl: warum 2+1,
   warum gleiche Teilchen, die Einheit G = 1.
4. **Kritischer Index (Operatoralgebren, kritische Statistik).** Objekt: die kleinste Inklusion mit unendlicher
   Tiefe (Jones-Index 4, Schleifengewicht 2, Vier-Zustands-Potts). Folgen soll: die Vier und die Minimalität aus
   einem Satz, zwei Blätter, der Übergang von endlich zu Kontinuum. Messgröße: c₃. Erste freie Wahl: warum die Welt
   genau am kritischen Index sitzt. Bestand: v813.
5. **Reinheit zweier Gebiete (Quanteninformation).** Objekt: Tr ρ² eines chiralen Vakuums auf zwei Gebieten.
   Folgen soll: zwei Blätter, vier Marken, Genus 1, Vierteltwiste, μ₄ als fermionischer Tausch, (−,−)-Spinstruktur.
   Messgröße: c₃, Zustand und Zeit. Erste freie Wahl: warum genau zwei Gebiete, welcher globale Zustand.
   Bestand: Befunde vom 18.11.
6. **Kleinstes anomaliefreies System (Topologie).** Objekt: das kleinste chirale System ohne ℤ₄- und ohne modulare
   Anomalie. Folgen soll: 16, N = 5, E₈. Erste freie Wahl: warum ℤ₄ und warum das kleinste. Bestand: F3, R102.
7. **Irreduzible Isometrie (Operatortheorie).** Objekt: einseitiger Shift. Folgen soll: Einheitswindung. Erste freie
   Wahl: Irreduzibilität; keine kontinuierliche Zeit. Bestand: Z30.
8. **Diamant des Beobachters (lokale Quantenphysik).** Objekt: das Vakuum auf dem Kausaldiamanten mit
   Verschränkungsgleichgewicht. Folgen soll: Zustand, Zeit, Reflexion, Einsteingleichung. Erste freie Wahl: setzt
   Raumzeit und Vakuum voraus (Prämisse D).
9. **Familien als selbstduale Zweiformen** und 10. **hintergrundfreie Naht** (H6, Prämisse A\*) stehen als Lesarten
   in der Suche, damit der Bestand mit den neuen Ansätzen im selben Maß gemessen wird.

## Die Gesamthypothese, die die Ansätze verbindet

Ein Grundelement: die ℤ₂ der Spinoren (eine Drehung um 2π ist −1). Daraus folgt in der Kette, die w2 nachrechnet:

- Zwei Blätter. Nur dafür fallen Index k², Markenzahl 2k/(k−1) und Liftordnung 2k zusammen, alle gleich 4.
  „Dieselbe Vier überall“ ist damit eine Eigenschaft von k = 2 und kein Zufall, der einzeln zu erklären wäre.
- Vier Marken mit Defizit π und Spannung 1/8; c₃ = 1/(8π) ist dieselbe Zahl wie die Einsteinkonstante je Defizit.
- Träger 5 mit 3+2 (einzige Zahl mit genau einer positiven Teilung), Familien 3; D_N ⊕ A_f verklebt bis N = 9
  nur bei (5, 3), also E₈ (darüber auch (13, 3), dort ist N − f ≠ 2).
- E₈ mit kompakter Spin(10) hat genau eine Realform mit nichtkompaktem Familienfaktor, E₈(−24), und der Faktor ist
  SU(2,2): 3+1 mit Lorentzsignatur.
- 3+1 hat die Drehgruppe SO(3) mit π₁ = ℤ₂. Die Kette schließt sich; d = 4 ist ihr einziger Fixpunkt.

## Begründung der Auswahl zum Rechnen

Gerechnet wird, was die Ansätze verbindet oder trennt, und der Weltgenerator selbst:

- w1: Keime (k, n, N, f, d) zu Welten ausrechnen, Regelsuche mit Gegenprinzipien, Schluss, Zahlen der
  Nachbarwelten, Verbünde der Lesarten, Hebel, Kontrolle mit eingepflanzter Lösung.
- w2: die exakten Brücken B1 bis B8, zuerst der Realformsatz (Ansatz 1), dann Kegel (3, 5), Index (4), McKay-Schnitt,
  Schleife, Teilung und Verklebung (6).

## Basislinien

- Regelsuche: sieben Gegenprinzipien, die die aufgezeichnete Welt nicht erfüllt (ℤ₃-Blätter, A-Verklebung, Ebene als
  Raum, Kreuzprodukt in R⁷, S⁶ als Himmel, gerader Träger, ungerade Marken). Ohne sie könnte keine Regel eine andere
  Welt wählen.
- Schluss: über alle 16 464 Keime gerechnet, unabhängig vom Klauselpool.
- Zahlen: alle 336 Keime (k, N, f) mit denselben Gesetzformen. Die verallgemeinerten Formen sind Konventionen
  dieses Generators (im Kopf von Abschnitt D genannt), für (2, 5, 3) geben sie die Ledgerwerte.
- Verbünde: eingepflanzte vollständige Lesart, mit drei Gegenproben.
- Leiter: so(d,2) für d = 3 bis 7 und die kompakte so(6) als Kontrolle.

## Nachtrag zum zweiten Durchgang (27. November)

Der Abgleich im ersten Durchgang lief nur gegen `TFPT-STAND.md` und die Atlas-Suche. Die Twistorbrücke des Ledgers
(WOIT.*, 9 Behauptungen, dazu CELEST mit 23) hatte er nicht gefunden. Ansatz 1 ist deshalb neu gefasst: Neu sind die
Realform E₈(−24), die Gleichsetzung des Familienfaktors mit der Twistoralgebra und die ℤ₄-Regel; die Signatur (2,2)
aus Uhr und reeller Struktur ist Bestand. Als elfte Lesart steht die Twistorbrücke selbst in der Suche. Zwei Rollen
kamen aus den eingefrorenen Vorhersagen dazu (Chiralität ohne Spiegel, paritätsungerade Antwort des Weltzustands).

## Nachtrag zum dritten Durchgang (27. November)

Der Flavortest (w5) hat die Aussage `family_is_twistor` in der Fassung des zweiten Durchgangs verworfen: Als
Cartan von sl(4) gelesen hat die Familienuhr die Coxeter-Klasse und ist nicht die Twistoruhr ρ₄. Der Bestand liest
die Familien aber selbst anders, als H⁰(K(D)) = Möbius-Algebra der Nahtgeraden (h6h33). In dieser Lesart ist die
Familienuhr −Ad(ρ₄), und die Familien sind die so(2,1) der Nahtgeraden in so(4,2). Ansatz 1 steht deshalb in
schwächerer Form: Die Signatur der Realform kommt aus dem Quadrat der Uhr; welche der beiden Lesarten des
Familienraums der A₃-Faktor von E₈ ist, ist offen. Die Transportmatrizen Q, K, R, L folgen aus keinem der Ansätze.
