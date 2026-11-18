# Nahtoperator, Feldwirkung, Zustand und Zeit

Diese Fortsetzung prüft die Herkunft des gemeinsamen Quellenoperators. Die
abgeschlossene E₈-Phasensuche wird nicht wiederholt. Das Ergebnis ist eine
bedingte Rekonstruktion samt exakt prüfbarer Herkunftsgrenze; eine vollständige
P1/P2-Herleitung wurde nicht gefunden.

## Der entscheidende Unterschied

Eine Windungszahl beschreibt eine Klasse von Nahtwirkungen. Um daraus die
tatsächliche Feldwirkung zu bestimmen, muss die Naht mit den vom ursprünglichen
Randkern erlaubten Zuständen zusammenpassen. Für den Randprojektor C und den
Nahtoperator U prüft der Rekonstruktionsweg nach Wold zuerst `(1-C) U C = 0`.
Dieser Test ist keine allgemeine physikalische Zulässigkeitsbedingung.

Auch eine glatte, streng orientierungserhaltende Naht mit Windung eins kann
diesen Test verfehlen. Der exakte Gegenfall im geprüften skalaren Randmodell
liefert die Leckage `15/4096`. Er zeigt die Reichweite des Windungskriteriums;
er ist keine Widerlegung sämtlicher TFPT-Quellen. Solche nichtisometrischen
Nahtwirkungen können trotzdem gültige Zustände in der gegebenen CAR-Darstellung
erzeugen. Der folgende Energietest schließt sie ausdrücklich ein.

Wenn der erste Test und die weiteren ausgewiesenen Erzeugungsbedingungen
gelten, lässt sich die Modenstruktur aus dem Operatorpaar rekonstruieren.
Wird zusätzlich die tatsächliche Zeitwirkung identifiziert, folgt die gemeinsame
Feld-, Zustands- und Zeitabbildung. Ein Windungsschritt verschiebt die Feldmode
um eins, erhöht die relative Ladung um eins und verändert die Energie um
`Q + alpha`. Bei `alpha=1/4` besitzt der erste erzeugte Zustand Ladung eins
und Energie `1/4`. Die Viertelphase bleibt erhalten. Ein solcher geladener
Windungsschritt lässt den Vakuumzustand nicht unverändert.

Der Rang eins dieser Modenerzeugung ist nicht der fünfteilige P2-Träger.
Die Zuordnung zum tatsächlichen Träger und zur E₈-Erweiterung muss ebenfalls
aus der ursprünglichen Quelle folgen.

## Ein stärkeres Auswahlresultat aus derselben Quelle

Bei unverändertem Zeitoperator und unverändertem Referenzvakuum hat jede skalare
Nahtwirkung `u = exp(i phi)` mit Windung eins die erzeugte Energie

`E(u) = alpha + Mittelwert[(phi' - 1)^2] / 2`.

Damit ist die gleichmäßige Windung `u = exp(i chi) z` das eindeutige Minimum,
bis auf eine Gesamtphase. Das gilt analytisch für die gesamte angegebene
H¹-Klasse, nicht nur für die 55 gerechneten Kandidaten. Ein verformter
Blaschke-Faktor mit Parameter `1/4` und sein nichtisometrischer Gegenfall
haben beide den exakten Energieüberschuss `1/15`. Der Energietest liefert
also eine positive Auswahl im selben Quellenmodell, ohne dessen Zeit passend
umzubauen. Er vergleicht die erzeugten Zustände, nicht nur die Operatornamen.

Die Eindeutigkeit bezieht sich auf skalare Multiplikationsoperatoren. Für
allgemeine unitäre Feldabbildungen können zusätzliche Vakuumsymmetrien dieselbe
Zustandsenergie besitzen. Und die Regel, dass die ursprüngliche P1-Naht dieses
Minimum wählen muss, ist bisher selbst nicht hergeleitet. Der Satz ersetzt
ebenso wenig die Herkunft von Feldalgebra, Vakuum und Zeitoperator.

## Was jetzt automatisch geprüft wird

- Nahtkandidaten als exakte rationale Funktionen, einschließlich Kürzung
  äquivalenter Präsentationen, Windung, Feldraum-Erhaltung und Zeitwirkung.
- Rekonstruktion des Operators und Vergleich mit dem vorgegebenen Zeitoperator.
  Eine nachträglich passend geänderte Zeit wird getrennt ausgewiesen.
- Tatsächliche CAR-Feldaktionen mit fermionischen Vorzeichen, Ladungsänderung,
  Zustandstransport, Energie und Viertelphase auf dem endlichen Anregungskern.
- Exakte Energie aller Kandidaten bei derselben ursprünglichen Zeit, samt
  unabhängiger Kontrolle über ihre Fourierkoeffizienten.
- Reproduktion der Ergebnisse und Bindung an die verwendeten Quellen.

Die analytische Klassifikation gilt für alle stetigen skalaren unimodularen
Symbole im Hardy-Modell: Eine windende Isometrie vom Grad eins ist ein einzelner
Blaschke-Faktor. Alle diese Operatorpaare gehören zu einer Rekonstruktionsklasse.
Der Computerlauf mit 55 Symbolen ist eine exakte Diagnose dieser Aussage,
kein Beweis durch endlich viele Beispiele und keine Suche über alle P1-Quellen.

Die Originalargumentation zur minimalen Quelle wurde ebenfalls geprüft:
Ein minimaler vierteiliger Defektwert beweist keine eindeutige Operatorstruktur.
Die dafür erforderliche Klassifikation wird im zitierten Beweis nicht geliefert.
Der Forschungsatlas führt diese Stelle jetzt als konkrete Beweisverpflichtung.

Auch das spätere Masterfunktional beseitigt diese Lücke nicht: In
`_archive/tfpt-45/source_extracts/02_carrier_source.tex:3020–3063` wird die
eindeutige diskrete Quelle bereits fixiert. Eine Barriere weist ihr den Wert
null und jeder anderen Quelle unendlich zu. Die anschließende Optimierung
betrifft die kontinuierlichen Größen auf diesem gewählten Zweig. Sie liefert
keinen unabhängigen Beweis seiner Auswahl und ist nicht identisch mit dem
oben hergeleiteten Energievergleich erzeugter Zustände. Insbesondere ist das
Alpha dieser elektromagnetischen Optimierung nicht automatisch der hier
festgehaltene Viertelphasenparameter.

## Ausführen

Im Ordner dieser Datei:

```sh
python3 search.py search --output results.json
python3 search.py verify --output results.json
python3 -m unittest test_search.py
```

Der Lauf endet nach der angegebenen Diagnose. Er läuft nicht unbegrenzt im
Hintergrund weiter und benötigt keine GPU. Neue Darstellungen derselben
klassifizierten Familie liefern keine weitere Herkunftsherleitung.

## Belege

- [PROOF.txt](PROOF.txt): Originalstellen, allgemeine Rekonstruktion,
  vollständige skalare Klassifikation und gemeinsame Feld-/Zustands-/Zeitwirkung.
- [results.json](results.json): tatsächlicher Diagnoselauf und Herkunftspflichten.
- [Lean](lean/SourceReconstruction.lean) und [Prüfbeleg](lean-receipt.json):
  Kompressionsidentität und zweiseitige Modenrekurrenz, jeweils mit Voraussetzungen.
- [Atlas-Eingabe](atlas_snapshot.json): verwendeter überprüfter Forschungsstand.
- [Unabhängige Prüfung](REVIEW.txt).

Experiment-only. Kein neuer Ledger-Beweis, keine physische P1/P2-Auswahl,
kein Ereignisgesetz und kein T1–T8-Gesamtabschluss. Die mathematische
Rekonstruktion darf ihre eigenen Eingangsvoraussetzungen nicht beweisen.
