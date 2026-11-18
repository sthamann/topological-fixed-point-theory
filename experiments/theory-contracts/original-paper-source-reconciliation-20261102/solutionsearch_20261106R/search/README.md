# Typisierte Suche nach Quellenprogrammen

Dieser Kern durchsucht Programme aus einer ausdrücklich gelieferten Operationsregistry. Ein Programm wird ausgeführt und anschließend vom tatsächlich eingesetzten Zielprüfer bewertet. Fehlende ursprüngliche Eingaben sind Befunde einzelner Kandidaten. Sie stoppen andere Kandidaten nicht.

Die Sprache ist endlich begrenzt durch Registry, Eingänge und die erklärte maximale AST Größe. Sie ist kein Anspruch, alle möglichen ursprünglichen Quellenprogramme zu umfassen. Genetische Fitness und physikalische Masse werden nicht gleichgesetzt.

## Schnittstelle

`Operation` oder kompatible Dictionaries enthalten:

```python
{
    'id': 'native.constructor',
    'input_types': ('OriginalKernelProgram',),
    'output_type': 'WholeSourceCandidate',
    'execute': trusted_python_callable,
    'version': actual_native_dependency_fingerprint,
}
```

`InputBinding(type_name, value, fingerprint)` enthält tatsächliche ursprüngliche Werte oder `Unavailable(reason)`. Ein unverfügbarer Eingang erhält seinen eigenen tatsächlichen Fingerprint, der sich beim Vorliegen neuer Originaldaten ändern muss. Originaldaten und relevante Kontexte gehören in den angegebenen Eingabefingerprint. Native fehlende Eingaben können `MissingInput` auslösen. Weitere konkrete native Ausnahmeklassen werden ausdrücklich mit `missing_input_exceptions=(NativeInputMissing,)` angegeben.

Ein AST hat ausschließlich diese Form:

```json
{"op":"apply","id":"native.constructor","args":[{"op":"input","id":"given_original_input","args":[]}]}
```

Nur IDs der gelieferten Registry und der gelieferten Eingänge werden akzeptiert. JSON enthält keine Pythonprogramme, Imports oder Shellbefehle. Die Registry selbst ist vertrauenswürdiger vorhandener Pythoncode, kein vom Modell importiertes Modul.

```python
engine = Synthesizer(
    operations=registry_operations(),
    inputs=original_bindings,
    validator=actual_original_verifier,
    required_obligation_ids=original_obligation_ids,
    target_scope=actual_original_verifier_scope,
    input_fingerprint=actual_original_input_fingerprint,
    out_dir=run_directory,
    target_type='WholeSourceCandidate',
    verifier_dependencies=original_verifier_and_proof_files,
    missing_input_exceptions=(NativeInputMissing,),
    context=actual_original_context,
)
```

Die Operationsliste des Source Agenten unter `source_recovery/native_constructors.py` passt unmittelbar zu diesem Vertrag. Sie enthält tatsächliche originale Kernel, gerade Grassmannpolynome, höhere Vertices, Sewing und markierte Quellenjets. Ihre Herkunfts und Anwendungsprüfung bleibt Aufgabe des injizierten ursprünglichen Verifiers.

Der Zielprüfer erhält `(candidate, context)` und gibt ein Dictionary mit `obligations`, optional `scope` und `evidence` zurück. Jede verlangte Obligation enthält `id` und einen Status aus `PASS`, `FAIL`, `MISSING_INPUT`, `NOT_PROVED`, `ERROR`. Ein fehlendes oder fehlerhaftes Dictionary bestätigt keine Obligation. Kandidateneigene Erfolgslabels werden nicht ausgewertet.

`SCOPED_HIT` entsteht ausschließlich, wenn sämtliche ausdrücklich verlangten Obligationen des Zielprüfers `PASS` sind. `FAIL` ergibt einen gescheiterten Kandidaten. Ein ausführbares, aber unbelegtes Programm bleibt `PARTIAL`. Der Zielscope stammt vom vertrauenswürdigen Prüfer beziehungsweise seinem aufrufenden Programm und wird nicht aus dem Kandidaten übernommen. Für den F1 Auftrag gibt Root die fünf vorhandenen Q Prüfungen, die ursprüngliche Quellenidentität und den ursprünglichen markierten Zeitwortabgleich vor. Ein voller TFPT Zielprüfer muss seine weitergehenden ursprünglichen Anforderungen selbst angeben.

## Zwei Wege durch denselben Prüfer

```python
engine.run(max_cost=8, max_asts=20000)
```

Die vollständige begrenzte Enumeration läuft fair nach `1 + Summe(Argumentkosten)`. Der Status nennt die letzte vollständig behandelte Kostenstufe und unterscheidet vollständige Abdeckung der deklarierten endlichen Sprache von einem ausgeschöpften Suchbudget. Keine dieser Angaben behauptet vollständige Abdeckung aller physikalischen Quellen.

```python
record = engine.evaluate_program(proposed_ast, max_cost=64)
```

Ein gezielter Vorschlag wird sofort über dieselben Typen, tatsächlichen Konstruktoren und den gleichen Zielprüfer ausgewertet. Er wartet nicht, bis die Enumeration alle Programme niedrigerer Kosten gebaut hat. Das ist für längere Quellenprogramme entscheidend. Eine solche Auswertung behauptet keine vollständige Enumeration bis Kosten 64. Ungültige Vorschläge werden lokal als `INVALID_AST` gespeichert.

`engine.enqueue(ast)` stellt zusätzliche typgeprüfte ASTs in `queue.jsonl`. Die Enumeration übernimmt diese innerhalb ihrer eigenen Kostenbegrenzung. Längere gezielte Modellvorschläge können stattdessen unmittelbar `evaluate_program` verwenden. `engine.publish_status()` schreibt den aktuellen Stand auch unabhängig vom Enumerationslauf.

## Speicherung und Wiederaufnahme

`results.jsonl` speichert jeden ausgewerteten AST zusammen mit AST Hash, Suchfingerprint, Ausführung und tatsächlichem Verifierbefund. Der Suchfingerprint umfasst die tatsächlichen Eingangsfingerprints und vollständigen Eingangswerte, Operationssignaturen und Versionen, Operationscode und gebundene Zustände, den Kerncode, Verifiercode und tatsächlich gebundene Closure, Default und Instanzdaten samt referenzierten Prüffunktionen, den tatsächlichen Kontext, typprüfende Callbacks und deren Zustand, ausdrücklich angegebene Originalabhängigkeiten, Scope und verlangte Obligationen. Sympywerte und NumPyarrays werden vollständig erfasst. Für unbekannte opaque Objekte gibt es keinen bloßen repr Fallback. Diese gebundenen semantischen Eingaben bleiben während einer Instanz unverändert; für geänderte Eingaben oder Prüferkonfiguration wird eine neue Instanz aufgebaut, die den geänderten Suchfingerprint automatisch erkennt.

Derselbe AST mit demselben Fingerprint wird nicht erneut geprüft. Ändert sich ein tatsächlicher Eingang oder relevanter Code, werden alte Ergebnisse nicht als Ergebnisse der geänderten Quelle wiederverwendet. Erfolgreiche Zwischenwerte können bei Bedarf durch ihre reinen Konstruktoren rekonstruiert werden; solche Abhängigkeitswiederholungen werden getrennt gezählt. Gespeicherte Darstellungstexte werden niemals als Python zurückgeladen.

Nach jeder direkten Programmauswertung und nach jedem Enumerationslauf entstehen:

* `status.json` mit Zielscope, Fingerprint, Verifiercodehash und tatsächlichen Suchzahlen.
* `best-unresolved.json` mit dem noch offenen Kandidaten mit den wenigsten offenen verlangten Obligationen, anschließend nach Kosten geordnet. Dies ist eine Softwareübersicht und kein physikalisches Gewicht.
* `scoped-hits.json` mit den vom eingesetzten Zielprüfer bestätigten Treffern.
* Rejectionsdateien für unbekannte oder falsch typisierte Vorschläge.

Die Dateien werden von einem lokalen Writer geführt. Der Kern startet keinen Dienst und wählt kein zusätzliches physikalisches Maß.

## Tatsächliche Methodenkontrollen

`method-control-run-v3.log` und `method-controls.json` enthalten 19 bestandene Kontrollen auf dem bereits vorhandenen ursprünglichen v330 Transfer. Dieser Fall prüft die Programmsuche und behauptet keinen F1 oder vollständigen TFPT Erfolg.

Die Suche findet das korrekte begrenzte Programm mit dem wirklichen endlichen Prüfer, verwirft ein korrumpiertes Programm und lässt einen ausführbaren symbolisch unbelegten Kandidaten offen. Eine fehlende primitive Operation stoppt die anderen Programme nicht. Geprüft sind außerdem Kostenfairness, Queue Deduplikation, Wiederaufnahme, geänderte Eingangsfingerprints, abgelehnte unregistrierte Befehle und die Wirkungslosigkeit eigener Erfolgslabels.

Der entscheidende Direktpfad wurde ebenfalls ausgeführt: Nach einem ausgeschöpften Enumerationsbudget wird ein höheres Programm unmittelbar gebaut und vom gleichen Prüfer bewertet. Der Status behauptet weiterhin nur die tatsächliche kleine Enumerationsabdeckung. Ein fehlendes Faultdictionary bleibt offen. Die zusätzliche Regression verändert bei gleichem Prüfercode tatsächlich seine geschlossene Referenz, während der deklarierte äußere Fingerprint unverändert bleibt: Der neue Zustand erhält einen neuen tatsächlichen Suchfingerprint und übernimmt den alten Treffer nicht. Gebundene Prüferdaten, Kontext und typprüfende Closures werden ebenfalls unterschieden.

Reproduktionsbefehl mit einem frischen Kontrollordner:

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/tfpt-discovery/.venv/bin/python work/solution-search20261106R/search/check_synthesis.py --out work/solution-search20261106R/search/another-method-control
```

P und Q werden von diesem Kern nicht verändert. Die tatsächliche ursprüngliche Gesamtlösung hängt von den gefundenen Konstruktionen und den wirklich bestandenen Originalobligationen ab.
