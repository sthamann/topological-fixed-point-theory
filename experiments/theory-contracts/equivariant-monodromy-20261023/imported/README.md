# TFPT: Monodromie und zentrale Erweiterung

Forschungsstand dieser mathematischen Ausarbeitung: 23. Oktober 2026.

## Dateien

* `REPORT.md`: Definitionen, vollständige analytische Beweise, endliche Klassifikation, alle 17 Vertreter.
* `01_monodromy.py` bis `05_cocycle.py`: exakte reproduzierbare Rechnungen.
* `expected_results.json`: vollständiges erwartetes Klassifikationszertifikat einschließlich der Clockmatrizen.
* `run_checks.py`: Neuberechnung in einem frischen temporären Verzeichnis; prüft das gesamte Ergebnis gegen das Zertifikat.
* `sources_manifest.json`: genau bezeichnete Originaldateien und Leseprovenienz.
* `verification_summary.json` und `verification.log`: tatsächlicher Ausführungsnachweis dieses Durchlaufs.

## Ausführung

Benötigt werden Python 3, NumPy und SymPy. Getestet mit NumPy 2.3.5 und SymPy 1.14.0.

```bash
python run_checks.py --optimized
```

Das führt alle fünf Skripte normal und unter Python `-OO` aus. Prüffehler lösen RuntimeError aus, unabhängig von Assertions. Ergebnisse werden in temporären Verzeichnissen neu erzeugt. Das mitgelieferte Ergebniszertifikat wird nicht als Recheninput verwendet.

Die Originaldateien wurden über den verbundenen GitHub-Lesezugriff eingesehen. Sie werden nicht als mitgelieferte, vollständig erneut ausgeführte Originalsuite ausgegeben. Die mathematischen Anfangsdefinitionen sind im Bericht angegeben; alle erforderlichen Daten werden vom Prüfer daraus neu erzeugt.

Keine Fernspeicherung, keine Änderung des Repositoriums und kein Netzaufruf durch die Prüfer.
