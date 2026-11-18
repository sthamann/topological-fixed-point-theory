# SOLUTIONSEARCH20261106R4: laufende Originalquellenforschung

Am 6. November 2026 wurde die tatsächliche automatische Originalquellenforschung gestartet. Herleitung und unabhängige Gegenprüfung laufen beide mit OpenAI `gpt-6.1-sol` über den bestehenden Codex Zugang. Die alte endliche Ausgabekampagne bleibt historisch erhalten und angehalten. Ihr `NO_COMPLETION_ROUTE` gilt für ihre registrierte endliche Klasse, nicht für diesen neuen Forschungsablauf.

Der neue Ablauf liest Originaldefinitionen, Atlasbefunde mit Voraussetzungen und Ausschlüssen, bisherige Experimente, Verifikationen und das vollständige eigene Rundenjournal. Die vollständigen Bestände sind als lesbare lokale Dateien erreichbar. Sie werden nach Bedarf gelesen und nicht als angeblich vollständig verarbeiteter Kurzprompt ausgegeben. Quellen und Stand werden pro Runde protokolliert. Jeder Ansatz muss eine entscheidende Frage, neue Information, zusätzliche Annahmen und offene Nachweise nennen. Der zweite, getrennte Sol Aufruf prüft den Ansatz und seine tatsächlichen Rechnungen; daraus entsteht automatisch die nächste Aufgabe.

Das geprüfte Claude Feedback ist enthalten: ursprüngliche Feldgewichte und derselbe ganze Caller, keine flache Ersatzdichte, Maßableitungen und Randfluss erhalten, Q5 mit echten Eingaben statt bloßen Funktionsnamen, keine Gleichsetzung endlicher Algebra oder Modellzustimmung mit ursprünglicher Physik. Aufrufe können Quellcode lesen und tatsächliche lokale Rechnungen im eigenen Ausgabeverzeichnis durchführen. Originalquellen und Prüfer bleiben schreibgeschützt.

Beobachteter Stand bei Einrichtung: drei vollständige Herleitungs-/Kritikpaare gespeichert, vierte Herleitung aktiv; alle sieben begonnenen Modellaufrufe mit festem gpt-6.1-sol. Berechnungen und Kritik sind unter ihren jeweiligen Bedingungen dokumentiert. Dies schließt weder die ursprüngliche gemeinsame Feldintegration noch T1..T8 ab. Der Host erkennt keine fundamentale Lösung aufgrund von Modellzustimmung an. F1, Originalowner und bestehende Abschlüsse bleiben erhalten. Eine wissenschaftliche Übernahme erfordert Quellenprüfung der konkreten Resultate.

## Aktiver Lauf und Wiederaufnahme

Laufverzeichnis: `/Users/stefanhamann/Documents/Codex/2026-11-05/erk/work/solution-search20261106R/runs/sol-original-source-research`.

Start mit `research_runtime/launch.py --repo TFPT_REPO --out RUN_DIRECTORY --rounds 0`. Null bedeutet fortlaufend, ohne festes Rundenende. Derselbe Aufruf übernimmt gespeicherte Rollen und Runden. Eine Dateisperre verhindert doppelte Controller. Einen laufenden gesunden Controller nicht für Routineänderungen neu starten. Der tatsächliche aktuelle Stand steht in `worker.json`; abgeschlossene Paare in `research-rounds.jsonl`, einzelne Rollen vor der nächsten Rolle in `rounds/*/round-state.json`. Eingaben, Befehle, Ausgaben und Aufrufbelege werden gespeichert.

Der lokale Monitor unter http://127.0.0.1:8768/ wählt standardmäßig `sol-original-source-research`. Historische Läufe sind ausdrücklich historisch gekennzeichnet. Auswahl oder Einfrieren der Ansicht verändert die Forschung nicht. Ein Aufrufproblem wird sichtbar als Fehler gespeichert; es gibt keinen Qwen/Ollama Ersatz. Ollama wurde auf Nutzeranweisung beendet.

Keine neue Veröffentlichung von Paper oder öffentlicher Webseite und keine physikalische Ledgerpromotion durch diesen Betriebsstand.
