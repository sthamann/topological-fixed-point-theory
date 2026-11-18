# Welche Architektur TFPT beim Lösen tatsächlich hilft

24. Oktober 2026. Auf dem Mac ausgeführte Vergleiche.

**Der stärkste getestete Kern ist eine exakte Suche über mögliche Konstruktionen, die gezielt unterscheidende Versuche auswählt und bewiesene Schritte wiederverwendet. Die hier getesteten lokalen Sprachmodelle liefern bislang keinen zusätzlichen vollständigen Lösungserfolg. Deshalb würde ich unsere Hardware zuerst um diesen überprüfbaren Forschungskern herum einsetzen.**

## Was sich gegenüber dem ersten Prototyp geändert hat

Der erste Versuch ließ ein Modell acht Binärwerte raten oder eine fertige Elimination aufrufen. Das prüfte vor allem den technischen Aufrufpfad. Diesmal mussten vollständige mathematische Bedingungen für vier zusammenhängende Konstruktionen gefunden werden:

1. Vier ursprüngliche E8 Transportwege mit einer räumlichen Reflexion.
2. Dieselben Wege mit einer anderen Reflexionsmarkierung.
3. Eine zusätzliche zyklische Rotation.
4. Rotation und Reflexion mit ihrer gemeinsamen Verträglichkeitsbedingung.

Die Aufgaben umfassen 40 beziehungsweise 48 Binärvariablen. Das Modell bekommt die ursprünglichen Operatorregeln und die verlangten Beziehungen, aber keine fertig hergeleitete Gleichungsmatrix. Der Prüfer erkennt sowohl fehlende Bedingungen als auch zusätzliche falsche Einschränkungen. Sein Gegenbeispiel ist eine konkrete Konstruktion, die Modellformel und ursprüngliche Operatoren unterschiedlich bewerten.

Alle vier Aufgaben gehören zum gleichen mathematischen Nahtprogramm. Vier Varianten sind kein Beleg für eine Übertragung auf unabhängige Theorien. Gegenüber dem ersten Bitversuch ist die Herleitungsaufgabe jedoch substanziell erweitert.

## Was ausprobiert wurde

| Ansatz | Ergebnis |
|---|---|
| Qwen 3.6 35B, vollständige Gleichungen direkt ausgeben | 0 von 4 vollständigen Aufgaben korrekt |
| Dasselbe Modell mit Denkmodus und größerem Tokenlimit | 0 von 4; drei Antworten erschöpften das Limit ohne Ausgabe, eine war mathematisch zu stark |
| Direkte Ausgabe mit neutralen Wiederholungsparametern | 0 von 4 |
| Bis zu zwei zusätzliche Selbstkorrekturen ohne Prüferfeedback | 0 von 4 |
| Dieselben ersten Vorschläge, bis zu zwei konkrete Gegenbeispiele zur Reparatur | 0 von 4 |
| Einzelne Beziehungen herleiten, prüfen, speichern und zusammensetzen | Kein vollständiger Abschluss mit den getesteten Qwen Modellen |
| Benannte Variablen und eine elementare Normalformbibliothek | In der ausgeführten Reflexionsaufgabe kein vollständiger Abschluss |
| Exakte Hypothesensuche mit zufälligen Beobachtungen | Alle 4 korrekt, durchschnittlich 72,55 Entdeckungsbeobachtungen |
| Exakte Hypothesensuche in fester Basisreihenfolge | Alle 4 korrekt, 129 Entdeckungsbeobachtungen |
| Exakte Hypothesensuche mit gezielter Auswahl | Alle 4 korrekt, 43 Entdeckungsbeobachtungen |

Die drei Zeilen mit neutraler Decoderkonfiguration verwenden dieselben Anfangsvorschläge. Damit ist Selbstkorrektur mit Gegenbeispielkorrektur direkt vergleichbar. Es wurden maximal drei Versuche pro vollständiger Aufgabe zugelassen. Längere Denkläufe mit Tokenlimits 4096 und 8192 wurden separat gespeichert. Die Antwortbudgets der verschiedenen Ansätze sind unterschiedlich; daraus wird keine Aussage über Recheneffizienz der Sprachmodelle abgeleitet.

Bei der Zerlegung wurden zusätzlich das vorhandene `qwen2.5-coder:32b` und verschiedene Eingabeformen ausprobiert. Einzelne triviale Beziehungen bestanden, eine vollständige Aufgabe bestand nicht. Ein weiterer Versuch mit kurzen Denkläufen wurde nach wiederholter Tokenerschöpfung beendet. Seine unvollständigen Aufgaben zählen nicht als fertig ausgeführter Vergleich.

Die Modellkennungen und Inhaltskennungen stammen aus dem lokalen Ollama Dienst. Diese Läufe prüfen die vorhandenen quantisierten Installationen und die verwendeten Aufrufparameter. Sie erlauben keine Aussage, dass alle größeren oder anders trainierten Modelle dieselben Grenzen haben.

## Was die erfolgreiche Variante tatsächlich lernt

Sie erhält eine ausdrücklich begrenzte Sprache möglicher Operatoren: XOR Summen von `I`, `A`, `K`, `AK`. Diese Sprache folgt aus den vorhandenen Clock und Reflexionsregeln. Verschiedene Schreibweisen reduzieren sich in dieser Quelle auf acht verschiedene Operatoren.

Die Suche hält diese Möglichkeiten offen. Sie wählt einen Eingabewert, bei dem die noch möglichen Operatoren möglichst unterschiedliche Antworten liefern. Anschließend wertet sie den ursprünglichen Gruppenoperator auf genau dieser Eingabe aus und verwirft die unpassenden Kandidaten. Sie bekommt keine ganze Zielmatrix als Antwort.

Damit wurden 17 unterschiedliche Wortbeziehungen gelernt, gespeichert und für 30 benötigte Beziehungen in den vier Konstruktionen wiederverwendet. Aus den gespeicherten Beziehungen werden vollständige Gleichungen zusammengesetzt. Diese werden erneut an den ursprünglichen Gruppenbedingungen geprüft und anschließend vom tatsächlichen Klassifikationsprogramm verwendet.

Die 43 Beobachtungen zählen nur die Entdeckung. Danach kommen in allen Varianten dieselben 344 Basisprüfungen für die 17 verschiedenen Zertifikate sowie die Prüfung der zusammengesetzten Systeme hinzu. Der Faktor 3 gegenüber 129 Basisbeobachtungen ist daher eine Einsparung bei der Entdeckung und keine gemessene Verdreifachung der gesamten Laufgeschwindigkeit. Die 20 Zufallswiederholungen benötigten 65 bis 78 Entdeckungsbeobachtungen.

Ein neuer Python Prozess hat die gespeicherten Beziehungen geladen, erneut geprüft und damit die nachgelagerte Klassifikation reproduziert. Es waren keine neuen Entdeckungsbeobachtungen nötig. Das ist Speicherung und Wiederverwendung eines gelernten symbolischen Modells. In diesem Vergleich wurden keine neuronalen Gewichte trainiert. Der frühere kleine LoRA Versuch bleibt separat dokumentiert; er verbesserte das Format, aber keine mathematische Testantwort.

## Mathematischer Fortschritt im vorhandenen Programm

Die exakte Auswertung liefert für die ausdrücklich festgelegten zusätzlichen Markierungen:

| Zusätzliche Struktur | Klassen mit Markierung modulo Charakterkonjugation | Zugrunde liegende Klassen im bisherigen Gamma Moduli |
|---|---:|---:|
| Reflexion0 | 65536 | 2176 |
| Reflexion1 | 65536 | 2560 |
| Rotation | 4096 | 256 |
| Gemeinsame Rotation und Reflexion0 | 16384 | 64 |

Die beiden Spalten zählen verschiedene Objekte. Rechts wird vergessen, welche zusätzlichen Lifts gewählt wurden. Die 64 bedeuten, dass genau 64 der zuvor 526336 Klassen eine gemeinsame Markierung der hier verlangten Art tragen. Das ist kein Beweis, dass die Natur eine dieser Klassen auswählt.

Zusätzlich wurde der reine Mapping Class Quotient des ursprünglichen geordneten Vierpunktmodells bestimmt. Die beiden relevanten Dehn Twists wirken hier durch kommutierende Involutionen. Aus 526336 Gamma Klassen werden 157696 Klassen unter dieser zusätzlichen Äquivalenz:4096 Orbits mit einer ursprünglichen Klasse, 46080 mit zwei und 107520 mit vier. Alle 1048576 Punkte des vorgelagerten Charakterquotienten wurden vollständig erfasst. Normalformen sind gespeichert.

Die topologische Erzeugung durch zwei Dehn Twists ist Standardwissen; sie wird beispielsweise in Abschnitt2 des [Manuskripts von Ferihe Atalan](https://archive.mpim-bonn.mpg.de/5138/1/mpim-preprint_2024-26.pdf) angegeben. Die Wirkung auf den hier verwendeten ursprünglichen E8 Lifts und die genannten Orbitzahlen wurden lokal hergeleitet und gerechnet.

Diese mathematischen Ergebnisse stammen aus der ausgearbeiteten Algebra und ihren exakten Programmen. Sie sind kein den getesteten lokalen Sprachmodellen zugeschriebener Forschungserfolg. Die Gleichungen und Beweise stehen in [MATH-SEAM-EQUIVARIANCE.md](THEOREMS.md).

Die 168 linearen Koeffizientenspalten wurden über einen zweiten Weg mit den ursprünglichen Kozykeloperatoren geprüft. Die Ränge wurden zusätzlich mit einer anderen Implementierung über GF(2) kontrolliert. 27 Tests bestehen, darunter Speicherung, tatsächliche Wiederverwendung und Ablehnung nachträglich verfälschter Lernresultate.

## Wie ich daraus das Gesamtsystem bauen würde

**Ein gemeinsamer Forschungsstand hält fest, was die vollständige Lösung verlangt.** Ausgangspunkt sind der aktuelle Atlas und die Originalquellen. Jede neue Herleitung nennt ihre Eingaben und das konkrete noch fehlende Resultat. Die aktuellen mathematischen Rekonstruktionen der lokalen E8 Quelle werden verwendet. Die weiterhin fehlende Identifikation der ursprünglichen P1 Quelle mit dieser Realisierung wird dadurch nicht vorausgesetzt. Auch Dimension, chirales wechselwirkendes Kontinuum, physikalische Dynamik, Zustandsauswahl und gemeinsame Gravitation behalten ihre eigenen Verpflichtungen.

**Ein starker Vorschlagsgeber sucht die nächste tragende Konstruktion.** Er soll geeignete Darstellungen, Abbildungen, Zerlegungen, Invarianten oder zusätzliche notwendige Bedingungen vorschlagen. Bekannte lineare Elimination und standardisierte Gruppenprodukte werden direkt gerechnet. Die getesteten lokalen Modelle haben die anspruchsvollere Rolle noch nicht bestanden. Vor einer großen Anpassung muss ein stärkeres Basismodell auf denselben Aufgaben und anschließend auf anderen vorhandenen Satzfamilien überzeugen.

**Der ausführbare Kern hält mehrere Kandidaten und ihre Unterschiede offen.** Für jede Kandidatenklasse gibt es einen passenden Prüfer. Bei endlicher Operatoralgebra funktioniert dieser Teil bereits. Bei analytischen Fragen brauchen wir entsprechend belegte Abschätzungen, Operatorbereiche und formale Beweisschritte. Eine neue Geometrie oder ein frei gewählter Zustand ist ein zusätzlicher Eingang und wird nicht durch eine gute Modellantwort zur Herleitung.

**Die nächste Frage wird nach ihrem Entscheidungseffekt gewählt.** Der lokale Vergleich zeigt dieses Prinzip in einer begrenzten Klasse. Bei offenen Fragen ist die Kandidatenmenge zunächst unvollständig. Dort muss das System auch erkennen, dass seine Sprache erweitert werden muss. Es darf aus einer leeren Kandidatenliste keinen allgemeinen Unmöglichkeitssatz machen. Der Beitrag eines trainierten Planers wäre, solche Strukturschritte und unterscheidenden Versuche schneller zu finden, wenn die vollständige exakte Auswahl zu teuer wird.

**Nur vollständig überprüfte Schritte werden wiederverwendet und zum Training zugelassen.** Gespeichert werden Frage, gültige Voraussetzungen, Vorschlag, konkreter Fehlschlag, Reparatur und abschließendes Zertifikat. Die bisherigen 17 Beziehungen sind ein Anfang für diesen Datenfluss, kein ausreichender Trainingskorpus für einen allgemeinen Forscher. Testdaten müssen ganze mathematische Familien zurückhalten. Weitere Zufallsinstanzen aus derselben Familie messen diese Übertragung nicht.

Für einen eigenen neuronalen Planer würde ich zuerst die Vorhersage erfolgreicher nächster Schritte und die Rangfolge von Versuchsvorschlägen trainieren. Später folgt das Erzeugen ganzer geprüfter Programme. Freie physikalische Erzählungen und richtige Endzahlen ohne Herleitung sind keine Trainingsziele. Der Erfolg wird daran gemessen, ob die Gesamtpipeline zusätzliche offene Aussagen korrekt entscheidet und ob sie ihre eigenen Resultate im nächsten Problem tatsächlich verwendet.

## Einsatz unserer Hardware

Die folgende Verteilung ist ein Ausbauvorschlag. Ausgeführt wurden die Versuche auf dem Mac; die anderen Rechner wurden nicht umkonfiguriert oder gestartet.

| Rechner | Vorgesehene Aufgabe |
|---|---|
| Zwei H200 | Zunächst stärkere Modelle gegeneinander testen; danach Anpassung eines Planers an geprüften Herleitungsverläufen |
| RTX PRO6000 | Unabhängiger Modellvergleich, Kritik an Vorschlägen und Suche nach Gegenfällen; exakte Algebra bleibt bei geeigneten Rechenwerkzeugen |
| Fünf DGX Sparks | Unterschiedliche bereits begründete Kandidaten oder Beweisschritte parallel bearbeiten und ihre Resultate an denselben Prüfer zurückgeben |
| Mac | Gemeinsamer Forschungsstand, direkte Entwicklungsversuche, exakte kleine Rechnungen und Einsicht in die gespeicherten Schritte |

Konkrete Modellgrößen, Verteilung von Gewichten und Durchsatz müssen auf den tatsächlichen Geräten geprüft werden. Die Zahl der Rechner wird kein Ersatz für eine fehlende tragende Abbildung. Mehr parallele Vorschläge helfen erst, wenn sich erkennen lässt, welcher Vorschlag das Gesamtproblem wirklich weiterbringt.

## Ausführen und Belege

```bash
.venv/bin/python -m tfpt_local.active_search --out results/active-search-new
.venv/bin/python -m tfpt_local.active_search --restore results/active-search-final/active-20261024.json
.venv/bin/python -m tfpt_local.independent_check
.venv/bin/python -m tfpt_local.mapping_class
.venv/bin/python -m pytest -q
```

Die wichtigsten Datensätze:

* `results/flat-control/`: Vergleich gleicher Anfangsvorschläge mit Selbstkorrektur und exakten Gegenbeispielen.
* `results/research-v2/` und `results/research-v3/`: ursprüngliche Gesamtversuche und Denkbudgets.
* `results/decomposed/`, `results/decomposed-coder/`, `results/named*/`: Zerlegung, Modell und Eingabevarianten einschließlich ehrlicher Fehlschläge.
* `results/active-search-final/summary.json`: drei Suchstrategien und 20 Zufallswiederholungen.
* `results/active-search-final/active-20261024.json`: gespeicherte geprüfte Beziehungen, zusammengesetzte Programme und tatsächliche Klassifikation.
* `results/active-search-final/restore.json`: erneute Verwendung in einem frischen Prozess.
* `results/independent-check.json`: zweiter Prüfweg.
* `results/pure-mapping-class.json` und `results/pure-mapping-class-normal-forms.npz`: vollständiger endlicher Quotient und seine Normalformen.

Der aktuelle überprüfte Fortschritt liegt in exakter Konstruktion, gezielter Experimentauswahl, Wiederverwendung und den genannten neuen mathematischen Klassifikationen. Ein neuronaler Lernvorteil für die vollständige TFPT Herleitung ist durch diese Versuche noch nicht nachgewiesen.
