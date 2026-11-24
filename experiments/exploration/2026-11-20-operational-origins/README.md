# TFPT: operative Ursprungsprüfungen vom 20. November 2026

**Status: zwei bedingte Kandidaten auf Stufe 1, keine vollständige Lösung.**

Dieses Explorationspaket enthält sieben Ursprungsansätze, drei Hauptrechnungen,
vollständig ausgezählte endliche Basislinien innerhalb der angegebenen
Modellklassen und mehrere gezielte Gegenprüfungen. Es baut
auf keinem gesetzten TFPT Randkern, Träger oder Clocknetz als Startobjekt auf.

## Lesen

* `IDEAS.md`: alle sieben Prinzipien, Messgrößen, freien Eingaben und Wildcard.
* `VERDICT.md`: Ergebnisse, Stufen, Basislinien, präzise Ausschlüsse und offene Brücken.
* `AUDIT.md`: Abgleich mit dem aktuellen Repository und allen zwölf Katalogblättern.
* `QEC_REPORT.md`, `VERDICT_CAUSAL.md`, `QEC_BRIDGE_AUDIT.md`: ergänzende Herleitungen.

## Reproduzieren

Python mindestens 3.10 und SymPy genügen; kein Netzwerk während der Rechnung.
Aus einem Checkout dieses Repositorys:

```bash
python experiments/exploration/2026-11-20-operational-origins/run_all.py
```

Bei einer isolierten Kopie des Ordners den Pfad zum bestehenden Ursprungssucher
angeben:

```bash
python run_all.py --originlab-root /absolute/path/to/research/originlab
```

Der Adapter importiert `genesis_rules.search` aus dem vorhandenen Suchprogramm.
Die anderen Skripte lassen sich unabhängig davon einzeln ausführen. Jede Ausgabe
ist JSON; `run_all.py` aktualisiert nur die zugehörigen Ausgabedateien hier.

## Hauptrechnungen und Kontrollen

| Datei | Inhalt |
|---|---|
| `qec_origin.py` | Minimaler perfekter Fehlercode, Alphabetselektion, vollständige Korrektur und Zugangstests, zyklische Basislinie |
| `syndrome_spread.py` | Alle GL4(F2) Abbildungen, Syndromsymmetrien, perfekter Kern, Orbits, Intertwiner und Charaktergegenprobe |
| `causal_origin.py` | SO(d) Klammerselektor, exakter Lorentzkegel, endliche Boostschranke |
| `sector_origin.py` | Zwei Sektoren als Fixpunktalgebra, minimale Blockgrößen, ablatierte Voraussetzungen und Kanalgegenbeispiel |
| `qec_bridge_audit.py` | Unabhängiger phasengenauer Codeprojektor und Rang aller Dreiplatzreduktionen |
| `qec_lift_check.py` | Physische symplektische Lifts und deren logisches Bild |
| `arithmetic_baseline.py` | Vollständige Formelzählung bei festgelegter Grammatik und Komplexität |
| `reversible_control.py` | Alle reversiblen binären Gatter bis drei Bits |
| `originlab_adapter.py` | Vorhandener Genesis Sucher auf neue operative Klassen angewandt |

Alle empirischen Zahlkonstanten bleiben außerhalb der Suche. Die resultierenden
Ganzzahlen sind nicht als unabhängige statistische Bestätigungen zu zählen.
Insbesondere werden 16 Syndrome nicht mit 16 Fermionen und drei logische
Pauli Richtungen nicht mit drei Familien identifiziert.

Keine Änderungen an Papieren, Atlas, Ledger oder dem laufenden Fundamentalplan.
