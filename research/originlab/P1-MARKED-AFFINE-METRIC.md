# Markierte Verklebung und reflektierte Form auf affinem Grad eins

Der Test verwendet das feste affine E₈ Vakuummodul auf Niveau eins.
Die Cartanwirkung, die Wurzeloperatoren und ihre Phasen stammen aus der
vorhandenen Farbprüfung. Im ursprünglichen R113 Aufruf bleibt das
gewichtete Randdatum samt absoluter Leerreferenz offen.

## Feinere Graduierung

Für Klebeklasse g modulo 4 und Farbe c modulo 3 ist der markierte
Grad d = 9g + 8c modulo 12. Er hat d(p₂) = 1. Auf jedem Gittervektor
ist g modulo 2 die Spinorparität. Damit ist d modulo 6 die gemeinsame
Spinor und Farbgraduierung in der Normierung an p₂. Der geordnete
Weg mit p₂, γ und −p₁ hat Schritte 1, 2 und 9 modulo 12.

Die Zählung erfasst alle Gitterladungen bis Normquadrat 6 und fügt
acht neutrale Heisenberg Oszillatoren hinzu. Für die sechs groben Klassen
ergeben sich auf affinem Grad 1 die Dimensionen (54, 48, 33, 32, 33, 48)
und auf Grad 2 die Dimensionen (762, 720, 657, 608, 657, 720).
Im feineren Grad 0 sind die Wurzelkomponenten A₃, A₂, A₁ und A₁.
Eine zusätzliche Cartanrichtung ist abelsch. Die 24 Wurzeln im
feineren Grad 6 tragen Y = ±1/2 und schwache Gewichte ±1.
Die grobe Klasse 0 fasst die beiden feinen Klassen zusammen.
Ihre Dimensionen auf affinem Grad 1 sind 30 und 24, auf Grad 2
sind es 382 und 380. Die feinere Zählung benötigt eine getrennte
Auflösung der Klebeklasse in der Quelle.

## Vier Gewichtungen und drei unabhängige Beziehungen

Sei G die reflektierte Form vor Zeiteinfügung und äußerer Auswahl
in einem als affiner Grad eins identifizierten Raum. Schreibe
G(v,w) = ⟨v,Dw⟩ relativ zur festen kompakten affinen Form.
Wahren die neutralen Nullmoden ihre gegebenen Adjungierten, kommutiert
D mit ihrer Wirkung. Die vier Summanden unter E₆ und A₂ sind die beiden
Adjungierten und die beiden geladenen irreduziblen Module mit
Dimensionen 78, 8, 81 und 81. Auf den geladenen Summanden sind die
Cartangewichte eindimensional und die durch neutrale Wurzeln verbundenen
Orbits jeweils vollständig. Schurs Lemma gibt je Summand einen freien
reellen Skalar.

Wahrt auch der markierte geladene Nullmodus J^{p₂}_0 sein festes
Adjungiertes J^{-p₂}_0, kommutiert D mit dieser zusätzlichen Wirkung.
Jeder ihrer nichtverschwindenden Blöcke erzwingt gleiche Skalare auf
seinen beiden Summanden. Die exakte Wurzelklammer liefert fünf solche
gerichteten Blöcke. Die Gleichungen haben Rang drei und lassen den
gemeinsamen Skalar übrig. Die positive neutrale Gewichtung (1, 1, 1, 2)
verletzt diese geladene Adjungierung.

Die Niveau eins Beziehung zwischen J^{p₂}_{-1} und J^{-p₂}_1 bindet
den gemeinsamen Skalar an die Norm des auf demselben Träger gebundenen
Vakuums. Die absolute ursprüngliche Referenz und der physische
Zeitschritt müssen aus derselben Randpaarung kommen. Der Test
bestimmt eine relative Form unter einer überprüfbaren Hypothese.
Die zeitversetzte Form Q₀(t), die gemischten R113 Formen und der
interagierende volle Träger bleiben eigenständige Quellenauswertungen.

## Ausführung

`python3 -B research/originlab/main.py --out /tmp/originlab-marked-source --campaign colour_ward --campaign r113 --campaign r113_synthesis --campaign f1`

Der Farbbericht enthält die Klassenzählung und die Wurzelzeugen unter
`marked_source`. Der R113 Bericht führt die ursprünglichen vier Marken
separat.
