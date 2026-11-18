# Tatsächlicher Programmstand am 6. November 2026

Der vorhandene Suchmotor kann jetzt räumlich veränderliche Higgsfelder mit dem ganzen ursprünglichen internen 96er Operator zusammensetzen und prüfen. Daneben gibt es einen getrennten Integrationsweg für tatsächlich gelieferte ursprüngliche Feldatome und ihre Gewichte. Der echte Integrationsstart wurde ausgeführt. Er hält an der fehlenden ursprünglichen Faserberechnung an, statt unbekannte Gewichte durch Fitness oder Einsen zu ersetzen.

## Zwei verschiedene Aufrufe

`runtime/launch.py --mode diagnostic_search` verwendet den vorhandenen `gpu_search.py` unverändert. Die ganzzahligen Gene wählen deklarierte rationale Fourierkoeffizienten aller vier realen Higgskoordinaten. Der Operator enthält alle 96 internen Komponenten, das vollständige gelieferte Majoranaelement, beide ursprünglichen Vierteldomänen und die Mischung verschiedener räumlicher Fouriermoden. Die geprüften Bedingungen sind Selbstadjungiertheit, Majoranasymmetrie, Fourieradjunktion, Graduierung, Realstruktur und die exakte Eingabestruktur. Numerische Restprüfungen erhalten ihre Toleranz. `feasible` bedeutet hier ausschließlich, dass diese Bedingungen pro Konfiguration erfüllt sind. Es bezeichnet keinen physikalischen Gesamtabschluss.

`runtime/launch.py --mode integrate` ruft den echten Beitragsauswerter auf. Dessen Konfigurationseinheit enthält die ganze ursprüngliche Historie mit beiden bosonischen Endpunkten oder ein vollständiges ursprüngliches Faserpaket mit solchen Endpunkten. Die vier unnormierten Matrixbeiträge W0, WU, WY und WUY bleiben bis zur Aggregation in der ursprünglichen Wortbasis. Die ursprüngliche Atommasse steht separat. Der gegenwärtig implementierte Plan enumeriert jedes gelieferte ursprüngliche Atom genau einmal. Verschiedene Atome mit gleicher sichtbarer Konfiguration werden beide gezählt. Weder die genetische Auswahl noch das Archiv bestimmen ein Integrationsmaß. Eine andere Vorschlagsregel oder Quadratur wird nicht stillschweigend angenommen.

Nach vollständiger Atomabdeckung folgen die gemeinsame Gramquotientenkarte, die einzige äußere physikalische Projektion, die Randpaarung und der absolute Referenzfaktor. Erst danach werden ZU/Z0, ZY/Z0 und ZUY/Z0 sowie deren verbundene Kombination gebildet. Der räumliche Einteilchenoperator h_N wird dabei niemals als gemeinsamer dynamischer OS Generator ausgegeben.

Vor dem Integrationsstart müssen echte Beiträge und eine explizite nichtleere Liste implementierter Ensembleprüfungen vorliegen. Noch nicht gemessene Prüfergebnisse werden erst nach der Aggregation ermittelt. Fehlende Implementierungen bleiben offen. Gegenwärtig existiert als ausführbarer Callback die Atomabdeckung. Die physikalischen OS, Ward, Schwinger Dyson, vollständigen Quellen und Grenzwertprüfungen sind noch nicht implementiert.

## Tatsächlich ausgeführte Wege

* `assembly-check.json`: sechs Prüfungen der ursprünglichen full96 Assemblierung, einschließlich 384er Modenmatrix und erhaltener Majoranadaten.
* `e2e-contribution-v2.log`: zwei Generationen des vorhandenen Suchmotors mit insgesamt vier verschiedenen räumlichen Konfigurationen. Alle implementierten Bedingungen pro Konfiguration bestanden. Die ursprünglichen Quellenbeiträge bleiben mangels Producer unausgewertet.
* `integrate-arithmetic-v2.log` und `runs/integrate-arithmetic-v2/`: vollständige numerische Ausführung des Beitrags und Aggregationswegs mit dem bereits vorhandenen ausdrücklich arithmetischen Kontrollpaket. Beide Atome gleicher sichtbarer Konfiguration wurden verwendet. Komplexe Gesamtwerte und verbundene Antwort wurden berechnet. Dieser Lauf liefert keinen ursprünglichen physikalischen Wert.
* `full-integration-start-v2.log` und `runs/full-integration-start-v2/`: tatsächlicher Integrationsstart mit dem vorhandenen ursprünglichen Diagnoseeingang. Exitcode 2, weil `original_blocked_fiber_packet` fehlt. W0, WU, WY, WUY und Originalmasse bleiben null. Die fünf geforderten physikalischen Ensembleprüfungen sind als fehlende Implementierungen ausgewiesen.
* `ensemble-registration-regression.json`: fehlende, leere, doppelte oder unvollständig implementierte Anforderungslisten ergeben keinen bestandenen leeren Prüfungsumfang.

Die früheren Dateien mit `e2e-v1` und `full-start-attempt-v1` dokumentieren den überholten ersten API Stand. Maßgeblich ist die getrennte Beitrags und Aggregationsschnittstelle v2.

## Start nach Installation im ursprünglichen TFPT Repository

Der installierte Paketpfad ist `experiments/theory-contracts/original-paper-source-reconciliation-20261102/fieldsearch_20261106P`. Die folgenden Befehle werden im TFPT Repository ausgeführt. Das Programm findet dessen Wurzel anhand seiner tatsächlichen Quelldateien über die übergeordneten Verzeichnisse. Außerhalb des Repositorys wird `TFPT_REPO_ROOT` ausdrücklich auf das tatsächliche TFPT Repository gesetzt. Ein ungültiger gesetzter Pfad wird zurückgewiesen. Absolute frühere Quellenpfade bleiben als Eingabe akzeptiert.

```sh
PYTHONDONTWRITEBYTECODE=1 experiments/tfpt-discovery/.venv/bin/python experiments/theory-contracts/original-paper-source-reconciliation-20261102/fieldsearch_20261106P/runtime/launch.py --domain experiments/theory-contracts/original-paper-source-reconciliation-20261102/fieldsearch_20261106P/domain/native-diagnostic-v1.json --out work/fieldsearch_20261106P/diagnostic --mode diagnostic_search --generations 2 --population 4 --batch 2
```

```sh
PYTHONDONTWRITEBYTECODE=1 experiments/tfpt-discovery/.venv/bin/python experiments/theory-contracts/original-paper-source-reconciliation-20261102/fieldsearch_20261106P/runtime/launch.py --domain experiments/theory-contracts/original-paper-source-reconciliation-20261102/fieldsearch_20261106P/domain/native-diagnostic-v1.json --out work/fieldsearch_20261106P/original-integration --mode integrate
```

```sh
PYTHONDONTWRITEBYTECODE=1 experiments/tfpt-discovery/.venv/bin/python experiments/theory-contracts/original-paper-source-reconciliation-20261102/fieldsearch_20261106P/runtime/launch.py --domain experiments/theory-contracts/original-paper-source-reconciliation-20261102/fieldsearch_20261106P/domain/arithmetic-control-domain-v2.json --out work/fieldsearch_20261106P/arithmetic-control --mode integrate
```

Der ursprüngliche Integrationsaufruf hält mit dem aktuellen Eingang reproduzierbar an der benannten fehlenden Operation an. Für eine tatsächliche vollständige Suche muss derselbe ursprüngliche Producer zu jeder zugelassenen Historienkonfiguration die vier markierten unnormierten Faserpakete einschließlich tatsächlicher Originalmasse liefern. Der bestehende räumliche Higgsoperator allein berechnet diese Pakete nicht. Fehlende ursprüngliche Kopplungen, körperliche Selektoren, Referenzfasern, Zeitdaten oder Domänen werden nicht gewählt.

Diese portable Kopie enthält die datierten früheren Protokolle unverändert. Deren absolute Pfade und Hashes bezeichnen den damaligen Lauf. `STATUS-original-snapshot.json` ist dessen erhaltenes Manifest. Die neue Pfadumstellung mit einmaligem Kontrolllauf und Originalstart steht in `PORTABILITY.json` und `portable_runs/`. Die arithmetischen Quellenformeln und der Quellenadapter bleiben unverändert.

## Herkunft und Wiederaufnahme

Der Diagnoseeingang bewahrt den vorhandenen reproduzierbaren v252 Zufallsdatenfall ausdrücklich als Diagnose. Er legt keine physikalischen Yukawawerte oder Ränge fest. Die räumliche Beispielvariation entspricht dem vorhandenen ursprünglichen variablen Hintergrund. Die Extraktion der ursprünglichen Higgsfunktion vermeidet die Nebeneffekte ihres Prüfskripts. Kanonische Quellen wurden nicht verändert.

Der vorhandene Suchmotor speichert sein Archiv und seine Wiederaufnahmedaten unverändert. Zur Quellenidentität gehören die Domänendatei, der neue Assembler und Startpfad, ursprüngliche v252 und Higgsfunktion, der Quellenadapter sowie die Bytes eines externen Quellenpakets und seiner angegebenen lokalen Quellen. Eine andere Quellenmasse darf dadurch keine alten Archivwerte wiederverwenden. Der Integrationsweg verwendet die Atomfolge und kein nach Eindeutigkeit gefiltertes Sucharchiv.
