# TFPT Wissensatlas

Ein lokaler Atlas über TFPT-Claims, einzelne Resultate, gemeinsam erforderliche Voraussetzungen, Prüfprogramme, Lean-Sätze, Herkunft der Eingaben und offene physikalische Anschlüsse. Der ursprüngliche Ledger bleibt unverändert und maßgeblich.

## Papers und einzelne Abschnitte

Die Papier-Ebene übernimmt die exakten `document:`-/`section:`-IDs und physischen
Fundstellen des zentralen Registers. Vollständiger Text und quellengebundene
KI-Beschreibungen erschließen Abschnitte; ihre Einordnung als geprüftes
Forschungsergebnis ist eine getrennte Aufgabe. `reviewed:false` bleibt sichtbar.
Aktuelle Inventar-, Abschluss-, Fehler- und Quellenfrischezahlen liefert
`paper_coverage`; vollständiges Inventar bedeutet nicht abgeschlossene KI-Auswertung.

Praktischer MCP-Weg: `paper_coverage` → `search_papers`/`list_papers` →
`paper_outline(document_id)` → `paper_section(section_id)` → `paper_evidence`
und `read_source` → `research_state` zur konkreten Forschungsfrage. Elternkontext,
Voraussetzungen, Gegenbefunde und ursprüngliche Ledger-/Lean-Geltung erhalten.
Die [gemeinsame Werkzeugreferenz](skills/tfpt-knowledge/references/papers-navigation.md)
beschreibt Argumente, Pagination, Frische und Datenschutz; der
[Navigationsleitfaden](../docs/SEMANTIC_NAVIGATION.md) trennt die Abdeckungsebenen.
Listen/Suche beginnen öffentlich. Private Inhalte erfordern expliziten Scope
oder eine bekannte genaue ID und dürfen nicht in öffentliche Artefakte fließen.

Die Papierdaten werden aus dem Register aufgebaut und anhand des vollständigen
Abschnittstexts ausgewertet. Vom Repository-Verzeichnis aus:

```sh
python3 semantic_explorer/papers.py --offline       # Inventar + gültigen Cache laden, ohne API
python3 semantic_explorer/papers.py --include-local --embeddings
python3 semantic_explorer/papers.py --check         # nur lesen: Quellen, Kontext und Abdeckung
```

Der Auswertungslauf setzt fehlende oder ungültige Einträge fort; unveränderte,
gültige Ergebnisse werden wiederverwendet. `--section section:ID` begrenzt einen
gezielten Lauf. Ohne expliziten Scope-Wechsel bleibt eine zuvor gewählte lokale
Ergänzung erhalten; `--public-only` schließt private Eingaben ausdrücklich aus.
`refresh.py` erneuert die Papieransicht offline, erhält passende Vektoren und
überschreibt keinen laufenden Import. Ein API-Lauf nutzt dieselbe vorhandene
Schlüsselkonfiguration wie der Atlas; Schlüssel gehören weder in Quellen noch
in Ausgaben. Neu erzeugte Beschreibungen verwenden robustes Klartext-Mathe und
gegebenenfalls ae/oe/ue; Originalzitate behalten ihre exakten UTF-8-Zeichen.

## Forschungsstand statt Wiederholung

Die Startseite führt durch dauerhaft gespeicherte Forschungsfragen. Jede Frage trennt belegte Teilresultate, zusätzliche Voraussetzungen, eng begrenzte Ausschlüsse und die nächsten Aufgaben. Ein historischer Gesamtvertrag mit `PARTIAL` öffnet dadurch nicht alle seine bereits erledigten Teilaufgaben erneut. Der Bestand unter `research/program.json` wird in denselben Atlas eingebunden, nicht als weiterer Theoriegraph geführt.

1. `research_state({"query":"Gedächtnis Zukunftsquotient GNS","lang":"de","max_chars":12000})` findet passende Fragen ohne Kenntnis von Dateinamen. `problem_index` bleibt die Navigation, `problems` enthält nur vollständige Pakete; bei Platzmangel die genaue Frage und das ausgewiesene Budget nachladen.
2. `assess_research_plan` prüft eine gefundene Strategie-ID mit ausdrücklich angegebenen Voraussetzung-IDs. Fehlende Voraussetzungen geben `scope_not_established`, veränderte Quellen `needs_review`. Eine Textähnlichkeit liefert nur Kandidaten. Ein Wiederholungsversuch braucht eine tatsächlich andere, genau benannte Voraussetzung.
3. `research_state({"query":"problem:future-quotient-and-gns-reconstruction","editable":true,"max_chars":60000})` liefert den vollständigen bearbeitbaren Datensatz und seine Revision. Nach Originalquellenprüfung aktualisiert `record_research_update` Befund, Aufgabe und Voraussetzungen gemeinsam. Pflichtargumente sind `expected_revision`, `problem` und `reason`; neue Voraussetzungen kommen in `premises`. Konkurrierende Änderungen führen zu einem Revisionskonflikt. Frühere Datensätze bleiben unter `research/history/` erhalten.

Quellen sind mit Pfad, SHA-256 und Zeilenbereich verbunden. Änderungen entziehen den betroffenen Ergebnissen ihre wirksame Wiederverwendung, ohne ihren historischen Inhalt zu löschen. Ein konstruktiver Lösungsnachweis und eine abgeschlossene Prüfung mit negativem Ergebnis bleiben verschiedene Aufgabenarten. Alle Voraussetzungen eines Belegs müssen bei seiner Verwendung erhalten bleiben.

Die erste geprüfte Entscheidungsschicht umfasst sieben zentrale Forschungsfragen; sie behauptet keine wissenschaftliche Vollprüfung des gesamten Korpus. `pending_attempts` zählt thematisch verknüpfte, noch nicht in diesen Entscheidungsstand eingearbeitete Versuche — keine neuen wissenschaftlichen offenen Fragen. `topic_coverage` macht fehlende Einordnungen sichtbar; der vorhandene Quellenbestand bleibt über `#discover`, die Suche und die bisherigen MCP-Werkzeuge zugänglich. Die Prüfung eines Datensatzes ist keine mathematische Wahrheitsprüfung und verändert den Originalledger nicht.

## Öffnen

```sh
./semantic_explorer/run.sh serve
```

Im Browser: <http://127.0.0.1:8731>. Der Atlas aktualisiert seine Daten beim nächsten Aufruf, wenn `data/knowledge.json` oder der ergänzende Kontext `data/context.json` atomar erneuert wurde. Die Statusanzeige unterscheidet vollständige Inventarisierung, KI-Auswertung, offene Extraktionen und Fehler.

## Daten aktualisieren

Python 3 und Rust/Cargo sind erforderlich. Der Extraktor verwendet nur die Python-Standardbibliothek. Er liest einen vorhandenen Schlüssel aus `OPENAI_API_KEY` oder einer expliziten Datei, deren Inhalt nicht ausgegeben oder in die Wissensdatei kopiert wird.

```sh
python3 semantic_explorer/extract.py --key-file /path/to/project/.env --workers 8 --embeddings
```

Der API-Aufruf ist kostenpflichtig. Standardmodell: `gpt-5.4-mini`, mit striktem Structured-Outputs-Schema; für die semantische Suche `text-embedding-3-small` mit 256 Dimensionen. Modell mit `--model` ändern. `--claims CLAIM.ID ...` wählt gezielte Aktualisierungen, `--offline` erstellt die Inventarisierung ohne API-Aufruf. Hashbasierter Cache verhindert unveränderte Wiederholungen; fehlgeschlagene Einträge werden bei einem neuen Lauf erneut versucht. Eine laufende Extraktion ist durch eine PID-Lockdatei gegen doppelte Writer geschützt. Vor dem Entfernen einer alten Lockdatei prüfen, ob deren Prozess noch lebt.

`data/knowledge.json` enthält den vollständigen ursprünglichen Ledger, statisch inventarisierte Scripts und Lean-Deklarationen sowie strukturierte KI-Annotationen und genaue Quellenverweise. Die mehrwertigen Strukturen passen in JSON ohne künstliche CSV-Zellenkonventionen. Der Graph kann später in andere Formate exportiert werden; die kanonische Datei bleibt JSON.

## Herkunft und Beweiskraft

Jeder erkannte syntaktische Eingabekandidat muss in einer KI-Einordnung vorkommen: tragende Eingabe, technische Einstellung, berechneter Zwischenwert oder ungeklärt. Nicht vollständig belegte Herkunft und ausgelassene Quellenbereiche bleiben sichtbar. Die statische Prüfung umfasst Zuweisungen, numerische Literale, Funktionsstandardwerte, Umgebungszugriffe und direkte lokale Importe. Sie ist keine vollständige dynamische Datenflussanalyse; indirekte Imports, Metaprogrammierung, unvollständige Quellbereiche und Laufzeitkonfiguration können zusätzliche Eingaben enthalten.

Die Annotationen unterscheiden gesetzte Axiome, externe Daten, fest eingetragene Werte, freie Parameter und berechnete Größen. Ein importierter Wert ist nicht automatisch aus TFPT hergeleitet. Rechenidentität, Zahlenübereinstimmung und kommutierende Matrizen begründen keine stärkere physikalische Identifikation. Gemeinsam benötigte Voraussetzungen laufen über einen Herleitungsknoten; sie werden nicht in mehrere scheinbare Einzelbeweise zerlegt.

Der Extraktor führt die Forschungs-Scripts nicht aus. Lean-Deklarationen sind statisch erfasst; formale Build- oder Kernel-Validierung ist damit nicht behauptet. KI-Auswertung ist stets `reviewed:false`. Der Gravitationsvertrag `seam-gravity-source-bridge-20261013` liefert einen gezielten Quellenfall für gesetzte Flächennormierung und den fehlenden modular-geometrischen Anschluss.

## Orientierung

Der Quellenabgleich vom 15.10.2026 ist auf der Startseite und unter
`/website-audit.html` erreichbar. `research/website-audit-20261015.json` hält
Websiteaussagen, Originalstellen, PDF-Seiten und den tatsächlich umgesetzten
Änderungsstand fest. Das Original-PDF und die öffentliche Website bleiben
dabei unverändert. Sieben positive, bedingungsgebundene Muster ergänzen die
vorher überwiegend auf offene Brücken gerichtete Übersicht. Die ursprünglichen
DET16- und Dynamikartikel sind direkt im Korpus inventarisiert.

`classification_reviews.json` korrigiert sechs zu breite Ausschlussanzeigen:
Ein ausgeschlossener Teilansatz oder eine erfolgreiche Must-fail-Kontrolle
hebt das daneben berichtete positive Resultat nicht auf. Jede Zuordnung ist
an den exakten Ledgerzeileninhalt gebunden. Bei einer Änderung wird sie beim
erneuten Kontextaufbau nicht angewandt und als veraltet markiert. Der Ledger
selbst bleibt maßgeblich und unverändert. Die Alpha-Entscheidungsschicht enthält
jetzt zusätzlich den vorhandenen globalen Eindeutigkeitsbeweis als erledigte
Teilaufgabe; Auswahl der Eingaben und physischer Transfer bleiben getrennt.

`orientation.json` behält die Seitenverweise auf die historische 25-seitige deutsche Kernessenz und verbindet aktuelle Ledger-Aussagen mit TFPT-Themen, Grundstrukturen und offenen Zielen T1–T8 sowie Universalraum. Die Karte unterscheidet vorhandene mathematische Strukturen, bedingte physikalische Anschlüsse und gesuchte Herkunftsherleitungen. Gemeinsame Seeds und Konsistenzprüfungen sind als korrelierte Zusammenhänge eingeordnet, nicht als unabhängige empirische Bestätigungen.

## Themen, Lösungsketten und bisherige Versuche

Themenansichten zeigen getrennt: im Ledger belegte Resultate mit Geltungsbereich, bedingte Anschlüsse, offene Fragen und gemischte Aussagen mit belegten Teilen bei offener Gesamtaussage. Der Stand der KI-Auswertung ist davon unabhängig. Ein Rechenbeleg schließt weder sein übergeordnetes Thema noch T1–T8 automatisch.

`build_context.py` übernimmt den Status aus dem unveränderten Ledger und inventarisiert bestehende Forschungsversuche, Gegenbefunde und Vertragsquellen. Auch neuere `contract_index.json` werden berücksichtigt. Der Altindex liefert zusätzliche Verweise, keine erneute Prüfung der Originalarbeit. Ein `PASS` bezeichnet den jeweiligen Vertrag; reine Indexeinträge zählen nicht als unabhängige Experimente. Thematische Zuordnungen können sich überschneiden.

```sh
python3 semantic_explorer/build_context.py
```

Die Lösungskette zeigt Voraussetzungen, Herleitungen, Ergebnisse und zugeordnete Prüfprogramme. Besitz- und Prüfbeziehungen unterscheiden sich von logischen Folgerungen. Querverweise sind separat zuschaltbar. Begrenzte Ausschnitte sind entsprechend markiert; eine dargestellte KI-Verknüpfung ist keine fachliche Abnahme.

Fehlgeschlagene und noch offene Extraktionen lassen sich mit `extract.py --resume-errors` nachbearbeiten. Erfolgreiche, quellenaktuelle Annotationen werden erhalten. `complete_import.py --key-file /path/to/project/.env` wartet auf einen vorhandenen Import und führt höchstens zwei solche Nachbearbeitungen durch. Der Status bleibt bei verbleibenden Fehlern ausdrücklich unvollständig.

## MCP und Skill

```sh
./semantic_explorer/run.sh mcp
```

Der Rust-Service bietet `research_state`, `assess_research_plan`, `record_research_update`, `graph_status`, `theory_overview`, `research_brief`, `research_context`, `compare_aspects`, `list_topics`, `search_knowledge`, `get_dossier`, `trace_reasoning`, `audit_inputs`, `read_source`, `topic_overview`, `solution_chain` und `list_attempts`. Kommunikation erfolgt als MCP über stdio. Der mitgelieferte Skill liegt unter `skills/tfpt-knowledge/SKILL.md`. Der Server ist zusätzlich als `tfpt-knowledge` in Codex und in `.cursor/mcp.json` registriert; neue MCP-Konfigurationen werden nach dem entsprechenden Client-Reload verfügbar.

Wenn ein einzelner geprüfter Problemdatensatz größer als das Antwortbudget ist,
liefert `research_state` eine ausdrücklich als `packet_state: "partial"` markierte
Seite. Die enthaltenen Befunde bleiben jeweils vollständig mit Geltungsbereich,
Premissen und Quellenverweisen. `budget.continuation` enthält den unveränderten
Aufruf samt `detail_offset` und `expected_revision` für die nächste Seite. Eine
zwischenzeitliche Atlasänderung verwirft die Fortsetzung, statt Seiten aus zwei
Revisionen zu mischen. Eine solche Seite belegt nur die ausgegebenen Datensätze,
nicht den Abschluss des gesamten Problems.

Hinzu kommen `list_papers`, `search_papers`, `paper_outline`, `paper_section`,
`paper_evidence` und `paper_coverage`. Ein bereits laufender MCP-Client muss neu
verbunden werden; die tatsächlich angebotenen Deskriptoren sind maßgeblich.

Der HTTP-Service bindet standardmäßig nur an `127.0.0.1`. Quellen können nur über im Graphen registrierte Referenzen gelesen werden, müssen innerhalb des Repositorys liegen und werden gegen ihre Hashes geprüft. Schlüssel-/Umgebungsdateien und Traversal-Pfade sind ausgeschlossen. Ein geänderter Quellenhash zeigt einen veralteten Auszug an und verlangt Neuauswertung; er widerlegt keine mathematische Aussage.

Kompatibilitätsleser können den aktuellen, zusammengeführten Atlas-Snapshot über `./semantic_explorer/run.sh export` als einzelnes JSON-Dokument auf stdout lesen. Das Kommando schreibt keinen Graphen zurück auf die Festplatte; es umfasst dieselben aktiven Overlays wie der Dienst.

## Prüfen

```sh
cd semantic_explorer
PYTHONPATH=.. python3 -m unittest discover -v -p 'test_*.py'
cargo test --manifest-path server/Cargo.toml
```

Die Tests prüfen Quellen- und Verweisgrenzen, MCP-Antworten, Eingabeformen, fehlende Auditabdeckung, gemeinsam benötigte Voraussetzungen und die Unterscheidung zwischen Querverweis und Beweis.

## Ein Atlas, eine Aktualisierungskette

`data/knowledge.json` ist der kanonische, vom Dienst geladene Wissensgraph. Die Quelldatei-Inventare und Projekt-Overlays werden in derselben Aktualisierungskette erzeugt; `source_inventory.py` liest Ledgers und freigegebene Originalartefakte nur vorübergehend im Speicher und schreibt keinen zweiten Graphen. Es führt keine Forschungsprogramme aus und liest den alten `verification/theory_graph.json` nicht.

Der empfohlene lokale Ablauf ist:

```sh
python3 semantic_explorer/refresh.py
python3 semantic_explorer/refresh.py --check
./semantic_explorer/run.sh export
```

Die Aktualisierung inventarisiert offline und erhält den Quellen-Cache, baut anschließend die Atlas-Quelldaten und den Forschungskontext und aktualisiert die strukturelle Qualitätsprüfung. `--check` liest nur und prüft registrierte Quellen sowie Ledger-Zeilen gegen den letzten Import. `run.sh export` gibt den aktiven, zusammengeführten Atlas-Snapshot als ein JSON-Dokument auf stdout aus; es erzeugt und speichert keinen weiteren Graphen. KI-Annotationen, Modellprüfungen und Rohquellen sind Navigations- und Reviewmaterial. Sie gelten nicht als mathematischer Beweis oder unabhängige Validierung.


## Gesamtbild und kompakte Forschungsarbeit

`/#map` zeigt die quellengebundene Karte P1/P2 → Naht/Cover → D5/A3/E8 → Clocks/Anker/Readouts sowie die offenen lokalen, dynamischen und gravitativen Anschlüsse. Ein Knoten trägt seinen Geltungsbereich und Quellenauszüge. Die Karte ist eine Orientierung, kein weiterer Beweis.

Themen beginnen mit Rolle im Gesamtbild, vorhandenem Kern, Verbindungen und der nächsten fehlenden Brücke. Der Themenkern wird durch transparente Navigationsregeln aus Claim-ID und Hauptsatz ausgewählt; verwandte Themen sind zuschaltbar. Diese Fokussierung ändert keine Graphbeziehungen und keine Ledger-Verdicts. „Status nicht automatisch zugeordnet“ bezeichnet eine Klassifikationslücke, keine offene physikalische Frage.

Ergänzendes Korpuspaket nach `research_state`: `research_context(id, max_chars=12000, offset=0)` erstellt ein begrenztes, deterministisches Forschungspaket aus Status und Geltungsbereich, expliziten gemeinsamen Seed-/Objekt-IDs, gespeicherten Graphverbindungen, früheren Versuchen, scoped Ausschlüssen, Eingaberisiken und Originalquellen. Gleiche Namen oder Symbole werden nicht zu einem gemeinsamen Objekt erklärt; Verknüpfungen und Arbeitsmuster sind keine Beweise. `/api/context?id=topic:gravity&max_chars=12000&offset=0` liefert dieselbe Ansicht für die Oberfläche. `research_brief(id, max_chars=12000)` bleibt der kompakte Belegsatz. `compare_aspects(ids, max_chars=12000)` vergleicht 2–6 Aspekte über explizite gemeinsame Eingaben und Verknüpfungen. Die aufrufende KI bildet die Schlussfolgerung mit benannten Voraussetzungen; der Server erfindet keine Folgerung. Auslassungen, Gesamtzahlen und Fortsetzung bleiben sichtbar. Erst danach gezielt Dossiers und Originalquellen öffnen. Der Skill fordert vor einer Wiederholung eines früheren Versuchs die Benennung einer tatsächlich geänderten Voraussetzung.

`search_knowledge` liefert kompakte Treffer statt vollständiger Vertragsdateien. `research_brief` bewahrt Status, fehlende Prüfdeckung und Seitenposition auch bei kleinen Antwortgrenzen; bei 1.000 Zeichen passen unter Umständen keine Versuchskarten. Dann mit größerer Grenze und unverändertem Offset fortsetzen. Die Antwortgrenze bezieht sich auf den JSON-Inhalt, nicht auf den zusätzlichen MCP-Protokollrahmen. Ausschlüsse innerhalb eines `PARTIAL`-Vertrags bleiben mit ihrem begrenzten Geltungsbereich sichtbar.

Die Suche meldet ihren tatsächlichen Modus (`hybrid` oder `lexical`). Semantische Vektoren sind nur für einen Teil der Einträge verfügbar; exakte Text-Hashes verhindern die Wiederverwendung eines Vektors für einen geänderten Text. Vorhandene Cache-Vektoren bleiben während laufender Import-Checkpoints verfügbar. `graph_status` weist diese Abdeckung aus. Ein passender Treffer ist eine Navigationshilfe; seine Voraussetzungen und Quellen müssen weiterhin gelesen werden.

Ein praktischer Forschungsweg ist: Thema eingrenzen, zwei Aspekte vergleichen, bisherige Versuche und Ausschlüsse lesen, die entscheidende Originalstelle öffnen und erst dann eine neue Schlussfolgerung formulieren. Die Karte enthält dafür auch ausdrücklich als Arbeitsinferenz markierte Verbindungen. Diese ersetzen keinen Herkunftsnachweis und keinen Beweis.

Freie Querfragen mit erkannten Fachbegriffen liefern in `research_state`
zuerst eine kompakte Auswahl von Mustern und Befunden aus mehreren Problemen.
Die Karten behalten Beleg, Voraussetzungen und Gegenbefunde. Mit den
Kennungen lassen sich vollständige Befunde und ihre Belege öffnen.
`theory_overview` zählt alle erfassten Befunde nach Status und Themen.
Mehrfach zugeordnete Themen dürfen sich überlappen. Mit
`search_knowledge(type="theory_pattern", query="")` lässt sich der
Musterkatalog öffnen. Die Suche trennt redaktionell erfasste Verbindungen
von Kandidaten, die sich dieselbe Originalstelle teilen.

`compare_aspects` akzeptiert auch Kennungen der Theorieübersicht
(`map:` und `pattern:`). Übereinstimmende Textabschnitte werden bei
der Suche nach ihrem Span Hash gebündelt und zeigen die weiteren
Fundstellen. Bei fachlichen Fragen bekommen Originalabschnitte Vorrang
vor Changelog und Verifikationsindex. Quellen, die in ein anderes
Repository umgezogen sind, liefern ihren fest referenzierten Link samt
Revision, Hash und Zeilenbereich.

## Korpus, Sprache und Qualität

`build_corpus.py` inventarisiert TeX-Dokumente mit Includes, Forschungsordner und anerkannte Ergebnis-/Vertragsartefakte. `data/corpus.json` dokumentiert Abdeckung, nicht auflösbare Referenzen und Größenbegrenzungen. `corpus_aspects.json` ergänzt quellenbezogene Themenaspekte; eine bloße Erwähnung bleibt ein Relevanzlink.

Die Oberfläche unterstützt Deutsch und Englisch. Themenorientierung und Gesamtkarte sind zweisprachig. `translate.py --key-file /path/to/project/.env --follow` ergänzt englische Claim-Zusammenfassungen anhand des genauen Ausgangstexts. Originalquellen und noch nicht übersetzte Detailannotation bleiben als Original sichtbar. Die Aktualisierung erhält die aktuelle Ansicht; neue Daten werden in der Statusleiste angezeigt.

`audit_extractions.py --follow` prüft Verweise, lokale Ableitungen und Abdeckung gegen die Quellenpakete. Das ist keine allgemeine semantische Zertifizierung. Quellenstichproben fanden tatsächlich falsche Parameterabhängigkeiten und Selbstverweise; diese Fälle wurden erneut ausgewertet. Der Extraktor verwendet für beobachtete ignorierte Aufrufargumente standardmäßig `gpt-5.6-terra` (`--risk-model`), sonst `gpt-5.4-mini`. Auch das stärkere Modell bleibt fehlbar. `review_semantics.py --run-review --claim EM.FP.01 --key-file /path/to/project/.env --out data/EM-review.json` erstellt einen getrennten zweiten KI-Prüfvorschlag mit quellengeprüften wörtlichen Belegen. Ein Vorschlag verändert weder Ledger noch wissenschaftlichen Status.

## Gemeinsamer Zugang für lokale Agenten

Codex und Cursor verwenden den MCP-Dienst `tfpt-knowledge` und denselben
kanonischen Skill. Claude Code hat sowohl einen Projekteinstieg (`CLAUDE.md`,
`.mcp.json`) als auch einen globalen TFPT-Skill; Claude Desktop verwendet
denselben lokalen Startbefehl. Nach Konfigurationsänderungen müssen bereits
laufende Clients ihre MCP-Verbindung neu laden. Neue Sitzungen beginnen über
die TFPT-Regeln mit `graph_status`, `theory_overview` und `research_context`.

Der alte generierte TFPT-Graph, sein Primer, die drei Query-/Build-CLIs und
der zugehörige Cursor-Hook sind entfernt. Ursprüngliche Ledger, Beweise,
Forschungsprogramme und der spezialisierte RH-Quellenkatalog bleiben erhalten
und werden vom Atlas erschlossen. Der Codegraph dient weiterhin ausschließlich
der Codesuche. Corpus-, Kontext- und Übersetzungsdateien sind interne Teile
desselben Atlas und keine konkurrierenden Theoriegraphen.

ChatGPT im Web erhält durch lokale Dateien oder Skills allein keinen Zugriff
auf einen localhost-Dienst. Dafür ist eine separat erreichbare und
authentifizierte Connector-Anbindung erforderlich. Die lokale Codex-Anbindung
ist davon unabhängig.

### Kontextsparende Detailabfragen

`get_dossier(id, max_chars=24000)` liefert über MCP standardmäßig eine kompakte Ansicht mit Ergebnissen, Voraussetzungen, offenen Punkten, Quellenpfaden und ausgewiesenen Auslassungen. `expanded=true` ist ausdrücklich zu setzen, wenn die vollständige Ansicht gewünscht ist; auch diese muss in `max_chars` passen. Die Browser-Detailansicht bleibt davon unabhängig. `research_context` und `research_brief` akzeptieren neben Themen und Claims auch explizite Ziel-IDs (etwa `goal:universalraum`) und einzelne Versuche (`contract:…`); Versuche liefern ihren eigenen Belegsatz, Ziele verwenden ihre direkt verknüpften Claims und Versuche. Themen kennzeichnen zusätzliche relevante Claims getrennt von ihren Kern-Claims.

## Neue Kernessenz und gezielte Lösungssuche (15.10.2026)

Die neue Fassung liegt als `TFPT-Kernessenz-Deutsch-2026-10-15-v2.pdf` (20 Seiten) und editierbare TeX-Datei im Repository. Der Korpus liest die neue TeX-Quelle; alte PDF-Seitenverweise bleiben historische Belege.

`research_state` kennt jetzt auch die vollständige Clock-Jordan-Rekonstruktion, den kovarianten J-Lift, die exakte Twistblindheit neutraler Vakuumwörter und die achtdimensionale geladene Antwortinverse. Hamiltonzeit, Carry und globaler gerader W4-Lift sind als vorhandene, jeweils bedingte Teilresultate erfasst. Die verschiedenen Quellverträge werden nicht stillschweigend gleichgesetzt.

Die [H200-Suchstrategie](research/H200-Suchstrategie-2026-10-15.md) nutzt diese Unterscheidung: exakte Operatorfilter zuerst, GPUs für große passende Antworten und unabhängige Gegenprüfungen, anschließend Zertifizierung und Rückführung in denselben Atlas. Es wurde dafür kein neuer GPU-Lauf gestartet.
