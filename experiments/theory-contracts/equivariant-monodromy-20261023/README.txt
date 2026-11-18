MATH.SEAM.EQUIVARIANT.02

THEOREMS.txt enthält den vollständigen Liftsatz, die explizite Auswertung
an G31, die Gruppenabbildungen und die Beweise der Äquivalenzaussagen.
INTEGRAL.tsv gibt alle17 integralen Kohomologiegruppen an.
CLOCKS.tsv gibt jeden der62 Clockvertreter mit seinem Obstruktionspaar an.
results.json enthält außerdem ihre exakten Matrizen und Gegenbeispiele
gegen die Vollständigkeit schwächerer Invarianten.

Reproduzieren aus diesem Verzeichnis, Python3 mit NumPy und SymPy:

  python3 -B -OO imported/01_monodromy.py
  python3 -B -OO imported/02_clocks.py
  python3 -B -OO imported/03_bundles.py
  python3 -B -OO imported/04_integral.py
  python3 -B -OO imported/05_cocycle.py
  python3 -B -OO checker.py

Die Reihenfolge ist nötig, weil jeder Prüfer auf den zuvor erzeugten
vollständigen Daten arbeitet. Das erwartete Zertifikat ist kein
Enumerationseingang. Der letzte Prüfer vergleicht die neu berechneten
Originaldaten vollständig mit diesem Zertifikat.

Zusätzlich wurde das ursprüngliche run_checks.py --optimized ausgeführt:
Es rechnet alle fünf importierten Prüfer in frischen temporären
Verzeichnissen normal und unter -OO und vergleicht sämtliche Ergebnisse.
Der tatsächliche Lauf steht in fresh-original-verification.log.
Der zusätzliche Prüfer verwendet aktive Fehlerprüfungen ohne assert.

Die Gruppen und die Äquivalenzbegriffe sind vor jeder Zahl festgelegt.
Insbesondere sind strikte Lifts, zentral kovarianter Transport und
projektive Orbifoldrelationen verschiedene gezählte Objektklassen.
Alle verwenden dieselbe vorhandene Erweiterung E über G31.
Das nichtzentrale Vorzeichenmodell der Ordnung1024 wird über konkrete
Gitterabbildungen verglichen und nicht mit E identifiziert.

Quellen und Downloadprovenienz: source_manifest.json und THEOREMS.txt§0.
Die Rechenprogramme haben keinen Netzbedarf. Ältere Forschungsdateien
wurden nicht geändert. Dieser Contract bleibt im Bereich experiments.
