# Neue Q Runtime mit tatsächlich angeschlossenen Prüfungen

Die kopierfähige Einheit besteht aus `runtime/`, `domain/` und `checks/`. Die daneben laufenden Forschungsdateien unter `source/` und `recovery/` gehören nicht zu dieser Programmkopie. Die beiden P Versionen sind unverändert.

`runtime/original_source_adapter.py` ist eine byteidentische Kopie des geprüften P Quellenconsumers. `runtime/source_adapter.py` delegiert dessen echte Beiträge und gemeinsame Aggregation und ergänzt die fünf mathematischen Callbacks aus `checks/original_checks.py`. Der tatsächliche Feldplugin lädt diesen Adapter und nimmt alle beteiligten Dateien in seine Quellenidentität auf. Eine primitive Beitragsberechnung wird durch diese Verdrahtung nicht vorgetäuscht.

Der native Start ruft den tatsächlichen Beitragsauswerter auf. Fehlt seine ursprüngliche Faserberechnung, endet er mit Exitcode 2. Die fünf registrierten Callbackfunktionen werden in diesem Fall tatsächlich aufgerufen und geben für ihre fehlenden ursprünglichen Eingaben `MISSING_INPUT` zurück. Damit unterscheiden sich eine nicht implementierte Funktion und eine implementierte Berechnung mit fehlenden Originaldaten.

Die Wahrheit der Ensembleprüfungen ist keine Vorbedingung der Integration. Nach der Aggregation bedeutet Exitcode 0, dass die vier Z Werte vorliegen und sämtliche ausdrücklich verlangten Prüfungen bestanden sind. Fehlen dann erforderliche Prüfdaten oder scheitert eine Prüfung, sind die numerisch berechneten Z weiterhin gespeichert, aber der Startpfad liefert Exitcode 2 und `SOURCE_AGGREGATED_REQUIRED_CHECKS_NOT_PASSED`. Auch ein bestandener arithmetischer Kontrolllauf bleibt `whole_physical_search_completed=false`.

## Tatsächlich ausgeführt

* Bestehender kleiner Checker: 23 von 23 Kontrollen bestanden. Darunter der echte Pfad `field_problem → ADAPTER → registry`, ein erfolgreicher an aktuelle Beiträge und Ergebnisse gebundener Callback und die Zurückweisung eines veränderten aktuellen ZU. Die vorhandene arithmetische Kontrollquelle wird im Paket byteidentisch mitgeführt. Nur deren ausdrücklich arithmetische Zahlen werden für diese exakte Prüfung rational restauriert.
* Nativer ursprünglicher Start: Exitcode 2, tatsächlicher W Aufruf, `original_blocked_fiber_packet` fehlt, alle fünf mathematischen Prüfungen `MISSING_INPUT`.
* Vorhandener arithmetischer Integrationskontrollfall mit der ausdrücklich allein verlangten Atomabdeckung: Exitcode 0 nach gemeinsamer Integration beider Atome.
* Derselbe unveränderte arithmetische Fall mit den fünf physikalischen Prüfungen als Anforderungen: Die vier Z Werte werden berechnet; die fehlenden ursprünglichen Prüfeingaben bleiben offen. Exitcode 2 nach Aggregation. Es werden keine passenden physikalischen Eingaben erfunden.

Die beiden zusätzlichen Dateien `domain/arithmetic-required-five-*.json` ändern ausschließlich die verlangten Prüfungs IDs, den Dateiverweis und die zugehörige SHA Bindung. Sie enthalten keine neuen physikalischen Parameter oder Gewichte.

## Start aus diesem Arbeitsordner

```sh
TFPT_REPO_ROOT=/Users/stefanhamann/Projekte/tfpt-theoryv4 PYTHONDONTWRITEBYTECODE=1 /Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/tfpt-discovery/.venv/bin/python work/contribution20261106Q/runtime/launch.py --domain work/contribution20261106Q/domain/native-diagnostic-v1.json --out work/contribution20261106Q/runs/next-native-start --mode integrate
```

Nach Installation im tatsächlichen TFPT Repository erkennt der vorhandene Pfadhelfer dessen Wurzel über seine übergeordneten Verzeichnisse. Dann kann `TFPT_REPO_ROOT` entfallen; Paket und Ausgabeordner werden relativ zum Repository angegeben.

Die Einzelbefunde und Dateihashes stehen in `runtime-wiring.json`. `runtime-wiring.diff` enthält ausschließlich die Änderung am vorhandenen Aufrufpfad. Die ganze mathematische Prüfung und ihre Grenzen stehen in `checks/ORIGINAL-CHECKS-AND-SCOPE.md`.
