# Zwei H200 für die Suche nach einer gemeinsamen TFPT-Lösung

Stand: 15. Oktober 2026. Konkreter Vorschlag auf Grundlage vorhandener Rechner,
Programme und neuer Quellensätze. Dafür wurde kein neuer GPU-Lauf gestartet.

## Was vorhanden ist

Auf `fieldai` waren bei der Live-Abfrage zwei **H200 NVL mit jeweils
143771 MiB Speicher** erreichbar und ohne aktive Computeprozesse. Die
Momentaufnahme reserviert die Karten nicht. Eine zusätzliche RTX PRO 6000
bleibt für diesen Vorschlag unberücksichtigt.

Der vorhandene Rechenpfad
`/home/aime/tfpt-gateb-gpu-20261005/source_mixing_20261008/`
enthält einen QMC-Backend mit stabilisierter komplexer Algebra und
Sampler-Kontrollen. Der gelesene Bericht vom 8. Oktober meldet sieben bestandene
Kontrollen sowie CPU/GPU-Abweichungen von rund 9,33e−15 für die Phase und
1,14e−13 für das Loggewicht im kleinen Kontrollsystem. Diese alten Kontrollen
wurden jetzt nicht neu ausgeführt. Sie belegen weder Gleichgewicht noch ein
physikalisches Ergebnis der früheren 128er-Produktion. Der ausdrückliche Stopp
vom 9. Oktober bleibt erhalten; der Vorschlag startet diesen Lauf nicht neu.

## Der entscheidende neue Unterschied

Die bisherige Suche darf zwei Informationsfragen nicht verwechseln:

1. **Welches quadratische Feld liegt vor?** Vier tatsächliche Clockbilder von
   P5 erzeugen über ihre Doppelpol-OPE alle 36 reell-symmetrischen Richtungen.
   Eine exakte rationale Inverse rekonstruiert die vollständige Matrix aus
   36 gemischten Antworten; bei vorgegebener Spur sind 35 Richtungen unabhängig.
   Vier Generatorfelder bedeuten nicht vier Messwerte.
2. **Welche geladene Zeitwirkung hat die Quelle?** Im vorgegebenen E8-Fockraum
   bleiben alle ausschließlich aus Heisenberg-/Jordan-Feldern erzeugten
   Vakuumantworten blind für den Twist δ, wenn jede einzelne Insertion alle
   acht Cartan-Ladungen erhält. Auf dem erreichten Sektor p=0 ist Kδ=N.
   Das gilt bei beliebiger Wortlänge. Mehr solche Momente oder GPU-Stichproben
   können die fehlende Information nicht liefern.

Geladene Zwischenzustände unterscheiden dagegen Kandidaten. Bei normierten
Primärzuständen und bekanntem gemeinsamem Zeitmaß gilt

    Cδ(p,t) = exp[−t(|p|²/2 + δ·p)],
    yi = −log[Cδ(bi,t)/Cδ(−bi,t)]/(2t),
    δ = (Bᵀ)⁻¹ y.

Acht unabhängige ursprüngliche Basisladungen, jeweils mit Plus-/Minus-Antwort,
reichen für ein uneingeschränktes δ. Ein engerer begründeter Kandidatenraum
kann weniger Probes benötigen. Diese globale Ladungspräparation ist noch kein
aus der ursprünglichen Naht hergeleitetes lokales physisches Instrument.
Die Antworten bestimmen δ, nicht automatisch q, die Vierteltransition oder
die gesamte wechselwirkende Quelle.

## Der Anschluss ist jetzt am tatsächlichen Quelloperator geprüft

Die originale Round20-Gitter-/Skalarquelle besitzt bereits onsite Ladungsfelder
T und balancierte Transporte W. Ein Start im total neutralen Sektor C0 macht
nicht sämtliche Ladungsproben blind: T führt vorübergehend in einen geladenen
Sektor. Diese Propagation muss im bestehenden Solver erhalten bleiben.

Die unveränderte Quelle behandelt entgegengesetzte Ladungen bei einem
ladungskonjugationsinvarianten Start jedoch gleich. Selbst beliebige feste
Peierlsphasen, einschließlich Viertelphasen und Flux, beheben diesen konkreten
Unterschied nicht: Die Kombination A aus Ladungskonjugation und komplexer
Konjugation bleibt erhalten. Bei A-invariantem Start sind die positiven
euklidischen Plus-/Minus-Rückkehrantworten weiterhin gleich.

**Ein großer Scan fester Hoppingphasen entfällt deshalb für dieses Ziel.**
Die Originaloperatorrechnung und ein unabhängiger Review haben dies bestätigt.
Der Test schließt andere Starts, andere aus der Naht begründete Feldabbildungen
oder A-brechende Quellenterme nicht aus. Diese müssen aber zunächst tatsächlich
aus der ursprünglichen Naht stammen.

Ein zusätzliches uniformes η·Q würde innerhalb C0 verschwinden und dennoch in
geladenen Zwischenzuständen das Verhältnis exp(−2tη·p) erzeugen. Das ist ein
präzises Antwortkriterium, keine Herleitung von η. Die Gleichsetzung η=Bᵀδ
passt nur diesen ungeraden Faktor; sie identifiziert nicht die ganze Quelle.

## So arbeitet die automatische Suche

**Atlas → Kandidat → exakter Filter → GPU-Antworten → Zertifikat → Atlas.**

Der Atlas liefert kompakte Resultate, Voraussetzungen und frühere Ausschlüsse.
Ein KI-Agent formuliert daraus eine konkrete Operator-/Quellenkonstruktion und
benennt jede zusätzliche Annahme. Er darf nicht durch Parameteranpassung eine
fehlende Herkunftsherleitung ersetzen. Vor einer großen Rechnung entscheidet
`assess_research_plan`, ob der Ansatz in derselben Klasse bereits erledigt,
ausgeschlossen oder für das Ziel unzureichend ist.

Dann prüft eine exakte Rechnung Ladungen, Cocycle/Kleinzeichen, Viertelphase,
Statistik, Zustandswirkung und Zeitoperator. Nur überlebende Kandidaten kommen
auf die GPUs. H200 A berechnet größere Quellen- und Mehrzeitantworten; H200 B
verwendet unabhängige Starts und passende Gegenkandidaten. Die Rollen können
vertauscht werden. Ein lokales KI-Modell kann Kandidaten vorschlagen, ersetzt
aber keinen dieser Prüfungen.

Ein numerischer Treffer wird anschließend durch Symbolik, rationale Rechnung,
passende Lean-Sätze oder kontrollierte Fehlergrenzen abgesichert. Der Atlas
speichert den konkreten Satz mit seinem Geltungsbereich. Ein Widerspruch stoppt
die betroffene Kandidatenklasse. Ein unentschiedener Lauf bleibt unentschieden.
Eine Änderung der Voraussetzungen erlaubt einen neuen Ansatz, muss aber im
Vergleich ausdrücklich sichtbar sein.

## Wo welche Rechenleistung sinnvoll ist

| Aufgabe | Geeigneter Weg | Aussage eines Erfolgs |
|---|---|---|
| Gitter-, Ladungs-, Cocycle- und Clockverträglichkeit | Exakte CPU/Symbolik | Das tatsächliche Operatorwörterbuch besteht die benannten Bedingungen |
| 36-dimensionale neutrale und achtdimensionale geladene Inverse | Bereits gelöste kleine exakte Rechnung | Rekonstruktion im jeweils angegebenen Antwortvertrag |
| Große Fock-/Gitterantworten und geordnete Zeitfolgen | GPU-Algebra, Krylov oder validierter QMC-Pfad | Numerische Antworten derselben deklarierten Quelle |
| Spiegelphase des vorhandenen Theta1-Modells | QMC nach erneuter Mischungsdiagnose | Bedingter Fortschritt bei einer chiralen Realisierung |
| Auswahl von Quelle, Instrument und Raumzeit aus P1/P2 | Herleitung plus tatsächliche Operatorprüfung | Fundamentaler Fortschritt, den ein guter Zahlenfit allein nicht liefert |

Die erste große Rechnung benötigt somit ein bereits bestandenes Operatorwörterbuch
und einen unterscheidenden Antwortvertrag. Erfolg bedeutet startunabhängige
Antworten mit kontrollierten Fehlern und erhaltenen Strukturbedingungen.
Dauerhafte Startabhängigkeit, dominierende Phasenfehler oder ein falsches
Wörterbuch führen zum Solverwechsel beziehungsweise Ausschluss des Ansatzes.
Ein längerer identischer QMC-Lauf behebt ungeklärte Mischung nicht von selbst.

## Der konkrete gemeinsame Suchauftrag

Vorhandene Hamiltonzeit, thermischer Zustand, interner Carry und der globale
W4-Lift sollen wiederverwendet werden. Die native Gitter-/Skalarquelle der
Runden 20–26 und der separat stipulierte chirale E8-Fockraum sind jedoch zwei
Quellverträge. Gesucht ist ihre aus den ursprünglichen Prinzipien begründete
Operator-/Zustandsverbindung mit tatsächlich verfügbaren geladenen Feldern,
Zeitwirkung und Instrumenten. Die neuen Antwortkriterien machen Kandidaten
prüfbar, ohne diese Verbindung bereits vorauszusetzen.

Die passenden gespeicherten Fragen lauten:

- `problem:common-source-state-event-selection`
- `problem:clock-linked-neutral-response-reconstruction`
- `problem:clock-source-covariant-transport`
- `problem:charged-twist-response-identifiability`
- `problem:native-charged-seam-selection`

Die Originalbelege liegen unter
`experiments/theory-contracts/clock-jordan-response-reconstruction-20261015/`
`experiments/theory-contracts/clock-source-twist-transport-20261015/` und
`experiments/theory-contracts/native-charged-conjugation-20261015/`.
Der erste Checker bestand 46 Prüfgruppen, der J-Lift 123 und sein
Identifizierbarkeits-Addendum 23. Die All-Charge-/Allwort-Aussagen beruhen auf
analytischen Beweisen; die endlichen Checks prüfen deren Implementierung.

Der native Gegenladungschecker bestand 24, die feste-Phasen-Erweiterung 16
exakte Prüfgruppen; beide wurden unabhängig wiederholt.
