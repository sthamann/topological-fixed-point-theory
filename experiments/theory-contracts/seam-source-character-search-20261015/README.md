# Automatische Suche nach der Markierung aus derselben Quelle

Die [Fortsetzung zur ursprünglichen Operatorabbildung](source_reconstruction/README.md)
prüft jetzt Calderón-Projektor und Nahtvertreter gemeinsam. Sie enthält eine
bedingte Feld-/Zustands-/Zeitrekonstruktion, die vollständige Klassifikation
stetiger skalarer Hardy-Symbole und den Originalquellenaudit der behaupteten
Minimalität. Die physische Auswahl aus P1/P2 ist dadurch nicht geschlossen.

Der Lauf ist abgeschlossen: Aus der gegebenen Zeitwirkung und der vollständigen
Quellladung wurde der gesamte Raum der Vierercharaktere bestimmt. Von 16
Charakteren haben zwölf genau Ordnung vier; vier davon erzeugen D₅+A₃.
Der Rang-5-Projektor wird aus ihren neutralen Wurzeln berechnet, statt als
Zielprojektor in die Kandidaten eingebaut zu werden.

Das ist eine vollständige Klassifikation **dieser bedingten Suchfamilie**.
Es ist keine vollständige physikalische TFPT-Lösung. Die vorhandene chirale
Quelle, ihre Faltung und ihre E₈-Erweiterung sind Voraussetzungen. Der Lauf
prüft die tatsächlichen Feldphasen und übernimmt den bekannten Zustandskonflikt,
damit ein interner Treffer nicht als physischer Quellennachweis erscheint.

## Was die Suche konkret macht

1. Aktuellen Atlas-Kontext und die verwendeten Originalquellen festhalten.
2. Aus den beiden Quellenoperatoren einen zweidimensionalen Charakterraum
   bilden. Ein bewiesener Gittersatz reduziert den kontinuierlichen Raum
   vollständig auf 16 Fälle; ein willkürlicher Rasterabbruch entfällt.
3. Alle 240 E₈-Wurzeln exakt auswerten, neutrale Komponenten anhand ihrer
   Wurzelgraphen und Ränge identifizieren und den Projektor rekonstruieren.
   Ein zweiter Lauf findet dieselben vier Kandidaten direkt in den ursprünglichen
   Quellenkoordinaten, ohne die Zielabbildung M oder einen Zielprojektor zu lesen.
   Erst danach wird die Übereinstimmung der Abbildungen geprüft.
4. Jeden Treffer auf die ursprünglichen acht komplexen Felder zurückziehen.
   Dort ergeben sich vier Phasen mit je zwei Feldern; ein elementarer
   Fünf-plus-drei-Träger folgt daraus noch nicht.
5. Zustand, Zeit und Viertelphase beibehalten. Die schon bekannte fehlende
   Gleichheit der unveränderten vollständigen Quellen wird weitergeführt.
6. Alle Kandidaten, ihre Voraussetzungen und offenen Beweisverpflichtungen
   speichern. Quellenänderungen machen eine Wiederverwendung ungültig.

Die ergänzende Gruppe hat 64 Elemente. Für diesen vollständig reduzierten
Suchraum sind H200 oder genetische Heuristiken unnötig. Die frühere
262.144-Kandidaten-Suche wird nicht erneut gestartet.

## Belege und Verwendung

- [Mathematische Ableitung](PROOF.txt), einschließlich Vollständigkeit und
  genauer Grenze zwischen Operatorzusammensetzung und physischer Auswahl.
- [Suchergebnis](results.json): alle Kandidaten, Projektoren, Feldphasen,
  Herkunft und verbleibende Verpflichtungen.
- [Lean-Zertifikat](lean/SourceCharacter.lean) und
  [ausgeführte Prüfung](lean-receipt.json): sieben überprüfte Sätze über
  Kompensation, Blatttransport und den vollständigen Gitterschnitt.
- `search.py`: ausführbare Suche und Wiederprüfung. Im Experimentordner:

  ```sh
  python3 search.py verify --output results.json
  python3 -m unittest test_search.py
  ```

  `python3 search.py search --output results.json` berechnet die abgeschlossene
  Familie erneut und holt einen frischen Atlas-Kontext. Das ist eine
  Reproduktion, keine neue Suchklasse. `test_search.py` enthält die Gegenfälle
  und Herkunftsprüfungen.

Die physische Viertelphase bleibt erhalten. Der zusammengesetzte interne
Operator `W = exp(4πiK) exp(-iπαQ/2)` ist kein Ersatz für die ursprüngliche
Zeitentwicklung. Dass gerade W die physische Nahtmarkierung ist, muss durch
die ursprüngliche Feldabbildung folgen. Diese Identifikation wird nicht aus
einem passenden Ergebnis rückwärts als Voraussetzung eingefügt.

## Fortsetzung ohne Wiederholung derselben Sackgasse

Die nächste Konstruktion benötigt den tatsächlichen geladenen Nahtoperator
und seinen gemeinsamen Definitionsbereich mit dem ursprünglichen Collaroperator.
Ein präziser Test ist die orientierte Shift-Kovarianz auf diesem Bereich,
zusammen mit dem Transport des Calderón-Zustands und der konkreten Feldwirkung.
Sie darf nicht lediglich aus dem topologischen Windungsindex angenommen werden.
Ein neuer Suchlauf muss an dieser Quelle etwas begründet ergänzen; eine andere
Phase innerhalb der vollständig klassifizierten Familie liefert keine neue
physikalische Lösung.

Experiment-only; keine Änderung des Ledgers, keine Paper-Promotion und kein
Abschluss von T1–T8. Formal überprüfte Algebra ersetzt keinen Herkunftsbeweis.
