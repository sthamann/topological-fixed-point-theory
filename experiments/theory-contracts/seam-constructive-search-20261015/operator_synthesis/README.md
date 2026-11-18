# Automatische Operatorsynthese

Hier wird nach einer **Bauvorschrift für eine Abbildung** gesucht. Die frühere
Suche nach zehn Phasenexponenten bleibt abgeschlossen. Ihre Zeit- und
Zustandsannahmen werden nicht als neue Suchvariablen ausgegeben.

## Zwei ausführbare Suchverfahren

**1. Operatorprogramme mit lernenden Gegenbeispielen (`clifford_search.py`).**
Der Generator bildet Wörter aus den zehn Majorana-Generatoren des gegebenen
Fünf-Moden-CAR-Trägers. Wortlänge, beteiligte Generatoren und linearer bzw.
antilinearer Typ sind variabel. Die CAR-Regeln entfernen äquivalente Wörter,
bevor sie Rechenzeit verbrauchen. Ein Prüfer untersucht alle Feldgeneratoren,
Ladung, die korrekte Viertelwirkung und eine vollständige reelle Basis aller
hermiteschen Einteilchenoperatoren. Eine verletzte Beziehung wird automatisch
zur zusätzlichen Bedingung. Gespeicherte Bewertungen werden bei geänderten
Bedingungen oder geändertem Programm verworfen.

Der erste Lauf ist eine **Wiedergewinnungskontrolle eines bekannten Ergebnisses**:
2048 selbst erzeugte Wort-/Typ-Klassen lassen zunächst zwei Kandidaten übrig.
Eine automatisch gefundene imaginäre Mischkomponente unterscheidet sie und
liefert genau `gamma_0 gamma_2 gamma_4 gamma_6 gamma_8 K`. Die Zielmatrix und
dieses Wort werden nicht vorgegeben. Die zugrunde liegende CAR-Darstellung,
ihre Besetzungsbasis und die geforderten Transportgleichungen sind Eingaben.

**2. Inverse gemeinsame Operatorsynthese (`intertwiner_search.py`).**
Dieser Löser braucht kein vorgegebenes Wort. Aus beliebig vielen zusammen
gehörenden Operatorpaaren sucht er eine gemeinsame Matrix `S` mit
`S A_j = B_j S` bzw. `S conjugate(A_j) = B_j S`. Er löst die Gleichungen exakt
über den gaußschen rationalen Zahlen. Ein Gegenbeispiel ergänzt die ganze
verletzte Operatorgleichung. Adjungierte werden ausdrücklich mit geprüft.
Ein invertierbarer Treffer kann durch seine positive Gram-Matrix polar
normalisiert werden; bei skalarer Gram-Matrix wird die unitäre Matrix direkt
ausgegeben. Die berechnete Dimension des gesamten Intertwinerraums enthält
alle Gleichungen, auch wenn schon ein früher Kandidat erfolgreich war.

Ein Nullraum der Dimension null liefert einen konkreten, durch Entfernen
einzelner Bedingungen minimalen Konfliktsatz. Ein erfolgloser begrenzter
Koeffizientenversuch heißt dagegen `full_rank_search_inconclusive`, nicht
Unmöglichkeitsbeweis. Die aktuelle generische Implementierung ist auf gleiche
Dimensionen bis 16 begrenzt; der spezialisierte CAR-Löser arbeitet auf 32
Dimensionen. Neue unendliche Felddarstellungen sind dadurch nicht abgedeckt.

## Ausführen

Im Repository:

```sh
python3 experiments/theory-contracts/seam-constructive-search-20261015/operator_synthesis/clifford_search.py --run-dir experiments/theory-contracts/seam-constructive-search-20261015/operator_synthesis/runs/recovery-seed7 --seed 7
python3 experiments/theory-contracts/seam-constructive-search-20261015/operator_synthesis/intertwiner_search.py experiments/theory-contracts/seam-constructive-search-20261015/operator_synthesis/example_operators.json --output /tmp/tfpt-inverse-result.json
```

Das JSON enthält Operatorpaare und erklärte Voraussetzungen, **keine Zielabbildung**.
Einträge sind ganze Zahlen, Brüche als Strings oder `[Realteil, Imaginärteil]`.
Beliebige Python-Ausdrücke und Fließkommazahlen werden nicht ausgeführt bzw.
nicht als exakte Eingaben akzeptiert. Das beigefügte Beispiel ist ein kleiner
synthetischer Funktionstest, kein zusätzliches TFPT-Modell.

## Drei weitere, sachlich unterschiedliche Suchrichtungen

1. **Vollständige positive Wortkerne vervollständigen.** Statt einen Hamiltonoperator
   zu raten, unbekannte gemeinsame Korrelationswerte suchen. Positivität,
   Adjungierung, die ursprünglichen Operatorrelationen und Reflexionsbedingungen
   werden gekoppelte Constraints. Eine konsistente, fortsetzungsstabile Lösung
   erlaubt die vorhandene GNS-Rekonstruktion. Markierung, Zustand und Zeit sollen
   aus demselben Kern kommen. Eine endliche positive Momentmatrix allein genügt
   nicht; erforderlich sind Erweiterbarkeit, die tatsächliche Zeit-/Feldwirkung
   und später der kontrollierte Grenzübergang. **Noch nicht implementiert.**
2. **Nichtkommutative Gleichungen statt Zahlen regressieren.** Die Wortgrammatik um
   dünne Summen, Kommutatoren und lokale Kopplungsterme erweitern. Auf H200 viele
   Strukturen vorsortieren; nur exakt normalgeordnete Identitäten oder belastbare
   Fehlergrenzen zählen anschließend. Die jetzige inverse Matrixsuche kann
   Abbildungen auch außerhalb einzelner Clifford-Wörter finden, ist aber keine
   symbolische Summen- oder Kontinuumssuche. **Diese Erweiterung ist offen.**
3. **Unterscheidungsversuche aktiv synthetisieren.** Überlebende Konstruktionen
   nach dem kleinsten Feldwort suchen lassen, bei dem ihre tatsächlichen Antworten
   auseinandergehen. Das priorisiert die nächste Rechnung. In der CAR-Kontrolle
   ist diese Funktion auf die feste vollständige Einteilchenbasis begrenzt;
   eine freie Suche nach längeren Versuchswörtern ist noch zu ergänzen.

Für die fundamentale Quellenauswahl ist Richtung 1 am unmittelbarsten.
Die beiden jetzt laufenden Löser dienen der automatischen Konstruktion und
Prüfung ihrer Operatorabbildungen. Sie wählen die positive Quelle nicht selbst.
Die methodischen Grundlagen sind CEGIS
([Solar-Lezama](https://people.csail.mit.edu/asolar/papers/thesis.pdf)) und für
die geplante Kernsuche nichtkommutative Momenthierarchien
([Pironio, Navascués, Acín](https://arxiv.org/abs/0903.4368)). Deren allgemeine
Konvergenzaussagen übertragen sich nur unter den jeweiligen Voraussetzungen,
nicht automatisch auf eine unbeschränkte TFPT-Feldtheorie.

## Beweisgrenze und Ressourcen

`lean/OperatorTransport.lean` beweist zwölf allgemeine Transportlemmata über
Summen, Produkte, Potenzen, Adjungierte und explizite Konjugation. Der
Lean-Receipt belegt diese Lemmas; konkrete Matrizen prüft der unabhängige exakte
Python-Pfad. Die Analytik eines Logarithmus oder eines unbeschränkten Operators
ist damit nicht formalisiert. Kein Lean-Satz erklärt einen numerischen Treffer
automatisch zur physischen Quelle.

Die kleine Kalibrierung ist auf der CPU bereits unter einer Sekunde fertig.
Auf den zwei verfügbaren H200 wurde deshalb **kein Scheinlanglauf** gestartet.
Sie werden erst für große Wort-/Kopplungsfamilien bzw. eine geeignete numerische
Vorsuche gebraucht; der abschließende exakte Prüfer bleibt davon getrennt.

**Verdict:** funktionierende Suchwerkzeuge und bekannte Quelle als Kontrollfall.
Keine neue P1/P2-Herleitung, kein ausgewählter ursprünglicher Zustand und kein
physisches Ereignisgesetz. Keine Änderung an Ledger, Verifikationen oder Paper.
