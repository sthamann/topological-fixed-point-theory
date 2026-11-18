# Ist die genetische Suche hier sinnvoll?

Im aktuellen Feldprogramm bringt ein längerer genetischer Lauf noch keine physikalisch gerichtete Suche. Ich habe den gesamten vorhandenen Raum mit 27 Konfigurationen ausgewertet. Alle 27 erhalten dieselbe Bewertung. Der Algorithmus kann also verschiedene Felder erzeugen, aber anhand der gegenwärtigen Zielfunktion keines davon als physikalisch besser erkennen.

Der Test verwendet genau den bestehenden Feldplugin und seine unveränderte Bewertung. Vier Ergebnisse stammen aus dem bereits ausgeführten Lauf, die übrigen 23 wurden durch denselben Programmeinstieg ergänzt. Die 23 zusätzlichen Auswertungen brauchten lokal rund 0,063 Sekunden. Jede Konfiguration enthält die vollständige 384×384 Matrix aus 96 internen Komponenten und vier räumlichen Fouriermoden. Es wurde kein zusätzlicher genetischer Mehrgenerationenlauf gestartet.

| Tatsächlicher Befund | Ergebnis |
| --- | --- |
| Vollständig ausgewertete deklarierte Konfigurationen | 27 von 27 |
| Unterschiedliche Zielfunktionssignaturen | 1 |
| Signatur für jede Konfiguration | `[0, 1]` |
| Verletzungen der implementierten Bedingungen pro Konfiguration | 0 |
| Größter numerischer Rest dieser Operatorbedingungen | 0 |
| Verfügbare ursprüngliche Beiträge W0, WU, WY, WUY | jeweils 0 von 27 |
| Verfügbare ursprüngliche Feldmassen | 0 von 27 |

Die erste Zahl der Bewertung zählt harte Verletzungen. Sie ist überall null. Die zweite Zahl zeigt, dass die ursprünglichen Quellenbeiträge fehlen. Sie ist überall eins und liefert keinen physikalischen Vergleich zwischen den Feldern. Auch `feasible=true` bedeutet hier ausschließlich, dass die implementierten Bedingungen pro Konfiguration erfüllt sind.

Die Felder sind tatsächlich verschieden: Die diagnostische Norm der Mischung verschiedener räumlicher Moden reicht von null bis etwa 3,9406. Diese Norm ist kein aus der ursprünglichen gemeinsamen Quelle hergeleitetes Auswahlkriterium. Sie wurde daher nicht nachträglich zur Fitness gemacht. Unterschiedliche Gene und Diversitätsmechanismen können den Raum weiter erkunden, liefern bei der einheitlichen Zielfunktion aber keinen physikalisch begründeten Fortschritt.

Für die jetzigen 27 Konfigurationen ist die vollständige Aufzählung bereits die abgeschlossene sinnvolle Suche innerhalb dieses ausdrücklich kleinen Diagnosefensters. Für einen größeren Raum wird die genetische Gegenfallsuche dann nützlich, wenn die ursprüngliche Quelle eine berechenbare Antwort oder eine konkrete Verletzung liefert, anhand derer Kandidaten verglichen werden können. Der erste fehlende Programmeingang bleibt `original_blocked_fiber_packet` mit den tatsächlichen markierten Beiträgen und der ursprünglichen Feldmasse.

Die gemeinsame Quellenintegration behält ihren eigenen Auswertungsweg. Genetische Auswahl, bevorzugte Kandidaten und ein Archiv eindeutiger Felder dürfen dabei keine ursprünglichen Gewichte oder die Häufigkeit gezogener Konfigurationen ersetzen. Dieser Test berechnet keine ganze Feldquelle und beweist keinen physikalischen Gesamtabschluss.

Daten und alle 27 Einzelbewertungen stehen in `results.json`, die beobachtete Ausführung in `run.log`. `assess.py` verwendet die vorhandene Plugin API. Die eingefrorenen Dateien unter runtime, domain und source wurden gegen STATUS.json geprüft und sind unverändert.
