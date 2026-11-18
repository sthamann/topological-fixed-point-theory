**Die automatische Suche kann jetzt mit gemischten Zuständen weiterarbeiten.**

Im bisherigen Programm gab es einen konkreten Abbruch: Es konnte nach solchen
Zuständen suchen, sie anschließend aber nur in einer unpassenden reinen
Darstellung abschließend prüfen. Dieser Pfad ist erweitert und ausgeführt.

Der Kontrolllauf funktioniert wie ein Rätsel mit einem fehlenden Feld:

1. Die bereits bekannte Quelle liefert 51 Antworten.
2. Die Maschine findet zwei mathematisch zulässige Ergänzungen.
3. Sie erkennt selbst die kleinste neue Frage, bei der diese Ergänzungen
   unterschiedliche Antworten geben.
4. Die vorhandene Quelle beantwortet genau diese Frage.
5. Danach sind alle 16 Einträge des gesuchten Operators in dieser Darstellung
   bestimmt. Seine tatsächliche Wirkung stimmt mit der alten Quelle überein.

Die Maschine kann damit Kandidaten erzeugen, Unterschiede finden und eine neue
geprüfte Antwort unmittelbar weiterverwenden. Sie braucht dafür keine reine
Zustandswahl und keine vorausgesetzte freie Wickregel. Ein weiterer Kontrollfall
prüft die neutralen Transfers zwischen zwei Moden. Insgesamt bestehen 37 passende
Tests, einschließlich absichtlich falscher Daten und der bisherigen Regressionen.

**Was damit noch nicht gelöst ist:** Der Versuch rekonstruiert eine schon
gegebene Quelle. Für die ursprünglich gesuchte P1 Quelle fehlt weiterhin die
Herleitung ihrer lokalen Einsetzung und der daraus folgenden gemeinsamen
Mehrfeldantwort. Die gezielt gelesenen Originalstellen definieren diese Antwort
über ein bereits eingeführtes Erzeugungsfunktional. Sie liefern an dieser Stelle
noch keinen zusätzlichen auswertbaren Zahlenwert.

Der nächste wissenschaftliche Test steht bereits in den alten Unterlagen:
Welche verbundene Vierfeldantwort liefert die tatsächliche P1 Quelle zwischen
den Randstücken? Der automatische Sucher kann eine solche Antwort jetzt prüfen
und weiterverwenden. Er darf sie nicht aus einem bevorzugten Kandidaten erfinden.

Die vollständige fundamentale TFPT Lösung ist damit weiterhin nicht bewiesen.
Diese Arbeit repariert den automatischen Suchweg und hält die noch benötigte
ursprüngliche Herleitung sichtbar.

[Prüfergebnisse](RESULTS.json) und [vollständige Argumentation](PROOF.txt).
