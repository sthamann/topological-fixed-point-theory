# Kohärente F1 Suche

Der Anhang legt eine kohärente Zweikopienquelle, eine reine Fourierreferenz,
die rohe relative Auslese auf dem K4 Block, geordnete Marker, einen
additiven Generator und einen kohärenten Viererdecoder vor. Die
[Gegenprüfung mit Fable 5.2]([private review])
priorisiert eine ausführbare gemeinsame Quelle, Energiefalsifikationen,
einen wiederaufnehmbaren Suchlauf und unabhängige formale Abnahme.

| Schritt | Implementierung und Kriterium |
|---|---|
| 1. Suchfamilie | Reine Fourierquelle, Mixer, rohe Differenzauslese, affine und nichtlineare Generatoren, beide Markerorientierungen, äußerer und dazwischengeschalteter Abschluss sowie Viererdecoder mit Energiebilanz prüfen. Der gemischte Viererzyklus, die ausgelassenen Zwischenrichtungen und die nichtlineare Clockalias dienen als Gegenproben. |
| 2. Fall auswerten | Die vier geordneten Antworten und drei reflektierten Gramblöcke aus einem erzeugten Transfer ausrechnen. Die unabhängig gewählte Referenz bleibt ein eigener Eingang. Die Ergebnisse durch den vorhandenen R113 Verbraucher mit exakten Werten kontrahieren und mit der direkten Formel vergleichen. |
| 3. Formal prüfen | Die aus der Mixerkomposition folgende Rekursion des Generators und die geordnete Markerwirkung für alle Besetzungsgrade in Lean prüfen. G0 bis G7 ausschließlich anhand gemeinsamer, physikalisch formulierter Aussagen annehmen. Die bestehende gepinnte Beweisannahme verwenden. |
| 4. Suche ausführen | Rationale Kandidaten nach der Summe ihrer reduzierten Nenner erschließen, die Summe pro Runde vergrößern und die vorhandenen Checkpoints und autonomen Jobs verwenden. Beweisbudget separat behandeln. Endliche Gleichheiten und Falsifikationen gelten für den jeweiligen Modelltyp. |
| 5. Originalbindung | Den Originalaufruf separat ausführen. Für die Annahme sind insbesondere der physisch gebundene Quellenvertrag, die ursprüngliche gemeinsame gewichtete Randform, der vollständige Cap, die absolute Referenz, die Zeitkarte sowie die gemeinsamen lokalen und Kontinuumsaussagen erforderlich. Fehlende Originalwerte bleiben im Protokoll sichtbar. |

Fables Produktinvarianzargument für zahlendiagonale Zustände setzt
eine von null verschiedene Vakuummasse voraus. Seine Phase aus
\(AB=iR\) wird als geordnete Auslese behandelt. Eine Identifikation
mit einer Phase des physikalischen Generators wäre eine weitere
Annahme. Der Decoder wird direkt gegen den additiven Zahlgenerator
geprüft. Führt ein Vakuumreservoir die Energiedifferenz ab, trägt es
unterschiedliche Zahlaufzeichnungen. Gewöhnliches Weglassen dieses
Reservoirs löscht die logische Kohärenz.

Auf dem isolierten logischen Block ist auch ein Ausgangstransfer mit
neuer Schrittvariable (q^2) exakt verträglich. Die Wahl der Skala muss
mit der ursprünglichen physikalischen Zeitkarte verbunden werden.
Für die unmittelbare Markeridentifikation gibt es eine weitere
Gegenprobe. Der originale skalare Higgsmarker (U) ist im entsprechenden
Fasertyp nichtnegativ. Der Oszillatormarker (B) hat negative Eigenwerte.
Die positive Projektion ((I+B)/2) ändert sowohl den Einzelwert als auch
die geordnete Antwort. Eine abgeleitete Originalzuordnung bleibt
erforderlich.
