MATH.SEAM.G31.LATTICE.03

THEOREMS.txt beweist die Nichtspaltung der ursprünglichen charakterwertigen
Gittererweiterung über G31 durch ein siebenzeiliges Originalzertifikat.
Auf allen22 echten Untergruppen des16er Zeugen existieren dagegen Sektionen.
Die62 ursprünglichen Clockmonodromien besitzen insgesamt11356 Liftklassen;
die17 ungeclockten Monodromien besitzen835456 Liftklassen.

Dieselbe ursprüngliche Gittererweiterung wird durchgehend benutzt. Die
zentrale Determinantenerweiterung aus MATH.SEAM.EQUIVARIANT.02 ist ein
anderer Vertrag. Ihre früheren Ausschlüsse werden nicht zurückgenommen.

Dateien:
  checker.py: Originalkokykel, kurzes Widerspruchszertifikat, alle Untergruppen,
              Vergleich zum endlichen Quotienten, vollständiges binäres G31 Bild.
  quarter_lifts.py: sämtliche11356 Normalformen und Automorphismenordnungen.
  audit_lifts.py: unabhängige Burnsidezählung aller62 Fälle sowie alle17 freien
                 Monodromiequotienten.
  CLOCK_LIFTS.tsv, ORDINARY_LIFTS.tsv: vollständige kompakte Ergebnistabellen.
  results.json, quarter-lifts.json, burnside-audit.json: exakte Zertifikate.
  probe.py und probe-results.json: reproduzierbare Suche des kurzen Zeugen.
  REVIEW.txt, claude-review.txt: genaue Reichweite der zusätzlichen Beweisprüfung.
  verification_summary.json: identische normale und optimierte Ausführung.

Ausführung im vorhandenen TFPT Repository mit Python3, NumPy und SymPy:
  python3 -B -OO probe.py
  python3 -B -OO checker.py
  python3 -B -OO quarter_lifts.py
  python3 -B -OO audit_lifts.py

Benötigt werden die bereits neu erzeugten Originaldaten unter
../equivariant-monodromy-20261023/imported/work/. Wenn diese noch fehlen,
vorher dort die mitgelieferten ursprünglichen Prüfer01 bis04 der Reihe
nach ausführen. Die Kokykelfunktionen werden durch AST direkt aus
../seam-joint-ward-state-20261014/operator_transport_check.py gelesen;
das übrige alte Modul wird nicht ausgeführt.

Normale Ausführung und -OO wurden auf identische JSON Ergebnisse geprüft.
Die globale Nichtspaltung ist durch sieben direkt im Originalkokykel
nachgerechnete notwendige Gleichungen belegt. Die endlichen Zählungen
ergänzen die angegebenen algebraischen Beweise.

Keine neuen Operatoren, keine andere Koeffizientengruppe, keine Änderung
von Papieren oder Ledger. Alle Dateien bleiben unter experiments.
