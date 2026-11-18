# Konstruktionssuche für die gemeinsame Nahtwirkung

Diese Suche konstruiert exakte Operatorrelationen. Sie verwendet keine Sprachmodelle,
keine API-Schlüssel und keinen frei an die Zielwerte angepassten Hamiltonoperator.
Das Forschungsziel bleibt eine ursprüngliche Quelle, aus der Markierung, Zustand
und Zeit gemeinsam folgen. Der erste implementierte Suchraum bearbeitet den
geladenen Defektanschluss der vorhandenen achtkanaligen CAR-Quelle. Diese Quelle
und ihre Zeitwirkung sind Voraussetzungen; ihre Auswahl aus P1/P2 wird durch
einen Treffer nicht bewiesen.

Aktueller Forschungsstand: [SYNTHESIS.md](SYNTHESIS.md). Die Fortsetzung trennt algebraische Treffer, positive Zustandsfortsetzung und tatsächliche Quellenverträglichkeit. Der erste GPU-Lauf ist abgeschlossen.

## Tatsächlicher Ablauf

1. `atlas_inputs.py` fragt den gemeinsamen Atlas über dessen MCP ab. Die vollständigen
   relevanten Problemkarten, bisherigen Ergebnisse und präzisen Originalstellen
   werden in `atlas_packet.json` festgehalten. Veränderte Originalquellen stoppen
   die Übernahme in ein Zertifikat.
2. `domain_spec.json` benennt feste Quelldaten, erlaubte Änderungen und offene
   Beweisverpflichtungen. Gesucht werden nur acht gegenseitige Phasen und ein
   Permutationsvorzeichen. Eine zusätzliche Phasennormierung wird als Eichwahl
   entfernt. Weder Zustand noch Energie werden passend optimiert.
3. Zwei genetische Populationen starten ohne vorgegebenen Treffer auf GPU 1 und
   GPU 2 des Rechners `fieldai`, den beiden H200. Mutation, Kreuzung, Auswahl und
   Austausch verwenden die tatsächlich ausgewerteten Constraints. Ein gemeinsames
   Archiv verhindert doppelte Bewertungen innerhalb eines Laufs.
4. Nach 20 Generationen prüft jede Karte einen disjunkten Teil des vollständigen
   Raums von 262.144 kanonischen Genotypen. Das trennt genetisches Finden von einer
   vollständigen endlichen Klassifikation. Danach ist dieser Suchraum beendet;
   die Karten werden nicht mit Wiederholungen beschäftigt.
5. Der lokale Sammler übernimmt Treffer automatisch in `proof_check.py`. Kandidaten
   dürfen ausschließlich kleine ganze Zahlen liefern. Ein festes Lean-Theorem
   prüft daraus universelle algebraische Identitäten. Kandidaten können keine
   eigenen Axiome, Theoremtexte oder ausführbaren Programme einschleusen.
6. Ergebnisse, tatsächliche CUDA-Nutzung, Quellversionen, Kandidaten und
   Lean-Ausgaben bleiben zusammen im Laufverzeichnis. Der Atlas erhält einen
   ausdrücklich bedingten Forschungseintrag; kein Treffer wird automatisch in
   den Theorie-Ledger oder als physikalische Lösung übernommen.

## Starten und nachsehen

Im Verzeichnis dieses Experiments:

```sh
python3 atlas_inputs.py --verify
python3 manage.py start --generations 20 --population 128 --batch 128 --hours 6
python3 manage.py status
python3 manage.py stop
```

Ein unveränderter bereits vollständig klassifizierter Suchraum wird wiederverwendet. Nur `--repeat` fordert eine ausdrückliche Reproduktion derselben Klasse an.

Die sechs Stunden sind eine Obergrenze. Der endliche Lauf endet früher, sobald
sein Suchraum vollständig klassifiziert ist. `current_run.json` nennt den
zuletzt gestarteten Lauf. Die Remote-Ausführung verwendet einen unveränderlichen
Code-Snapshot unter `/home/aime/tfpt-constructive-search-20261015/code/`.
GPU 0 wird nicht verwendet. Der Mac führt Lean aus; die H200 prüfen die
Kandidaten. Während des Sammelns verhindert ein an den Sammler gebundener
Ruhezustandschutz den automatischen Schlaf des Mac.

Ein abgebrochener GPU-Prozess kann mit genau seinem gespeicherten Aufruf aus
`launch.json` wiederaufgenommen werden: Die Engine lädt den passenden Checkpoint
und lehnt geänderte Eingaben ab. Den Sammler startet man erneut mit
`python3 manage.py collect --run /absoluter/pfad/zum/run.json`.

## Was ein grünes Ergebnis bedeutet

`feasible` aus der GPU ist ein vorläufiger Constraint-Treffer. `ok` in einem
Lean-Beleg bedeutet ausschließlich, dass dessen ausdrücklich ausgewiesenes
Theorem vom Kernel akzeptiert wurde. Die ursprünglichen Produktionszertifikate deckten die Zustandsfortsetzung noch nicht ab. Die Fortsetzung in `continuation/state_expectation/` beweist nun die allgemeine endliche Wortpositivität in Lean und die konkrete Quellenübertragung analytisch. Sie bestimmt keinen physischen Zustand auf einer größeren Observablenalgebra.

Vier Phasenpräsentationen reduzieren sich durch die vorhandene Parität auf zwei globale markierte P-Lifts. Der originale Einzelwindungsschritt lässt beide zu; die zusätzlich stipulierte E8-Stromgruppe wählt bedingt das P-Zeichen +1. Der Zustandsvergleich zeigt, weshalb gleiche Verschiebungsantworten verschiedene physische Grunddichten verbergen. Der vorhandene CAR/QWZ-Quellenausschluss gilt unabhängig von allen vier Phasenbelegungen. Details, Voraussetzungen und der konkrete nächste Suchvertrag stehen in `SYNTHESIS.md`.

## Dateien und überprüfbare Ergebnisse

- `DESIGN.md`, `domain_spec.json`: mathematische Suchklasse, Voraussetzungen und Grenzen.
- `problem.py`: Konstruktion, Mutation, Kreuzung und exakte ganzzahlige GPU-Constraints.
- `gpu_search.py`: zwei Populationen, Archiv, Checkpoints und endliche Klassifikation.
- `proof_check.py`, `lean/FixedCertificate.lean`: exakte Gegenprüfung und Lean-Zertifikat.
- `test_problem.py`: unabhängige Normalform-, Zeit-, Mutanten- und CPU/CUDA-Prüfungen.
- `test_gpu_search.py`, `proof_tests.py`: Suchbetrieb und Beweispfad einschließlich Negativfällen.
- `runs/*/snapshot.json`: GPU-Status und Klassifikation des konkreten Laufs.
- `runs/*/proofs/*/receipt.json`: Ergebnis, Umfang, Axiome und Versionen jedes Beweises.
- `source_compatibility.py`, `runs/*/proofs/*/assessment.json`: separate Bewertung der unveränderten Quellenverträglichkeit; `status` zeigt sie neben den Algebra-Zertifikaten.
- `results.json`: kompakter, tatsächlich geprüfter Stand dieses Experiments.

Die erste Version automatisiert die Konstruktion und Zertifizierung dieser
ausdrücklich definierten Algebra. Sie ist kein unbegrenzter automatischer
Erfinder neuer physikalischer Axiome. Erweiterungen müssen eine neue, aus
Originalquellen begründete Suchvariable und einen entscheidenden Test benennen.
