# Feldsuche und genetischer Algorithmus

Der Suchmotor wurde um räumlich veränderliche Higgsfelder und einen getrennten Integrationsweg erweitert. Die kleine Diagnoseklasse mit 27 Konfigurationen ist vollständig ausgewertet. Ein genetischer Algorithmus hat darin derzeit kein physikalisches Auswahlsignal: Alle Konfigurationen erhalten dieselbe Bewertung. **Die vollständige physikalische Suche läuft noch nicht**, weil die ursprünglichen gemeinsamen Quellenbeiträge noch nicht berechnet werden.

## Was tatsächlich funktioniert

Die Feldassemblierung verwendet den vorhandenen ganzen internen Operator als 96 mal 96 Matrix und erhält die gegebenen Yukawa und Majoranadaten. Mit den vier ursprünglichen räumlichen Fouriermoden entsteht eine 384 mal 384 Matrix. Die Moden bleiben gekoppelt. Diese innere Matrixdimension bezeichnet weder eine Raumzeitdimension noch eine Zahl unabhängiger Teilchensorten.

Der unveränderte bestehende genetische Suchmotor wurde mit dem neuen Feldplugin ausgeführt. Zwei Generationen verarbeiteten vier verschiedene räumliche Konfigurationen. Anschließend wurde der vollständige kleine Raum mit 27 Konfigurationen geprüft, wobei die vier vorhandenen Bewertungen wiederverwendet wurden. Die Eingabedaten dieses Laufs sind die ausdrücklich als Diagnose gekennzeichneten ursprünglichen v252 Daten. Sie legen keine physikalischen Yukawawerte oder Ränge fest.

Der Integrationsweg ist davon getrennt. Er zählt jedes gelieferte ursprüngliche Integrationsatom mit seiner Masse und erhält komplexe Beiträge. Verschiedene Atome mit derselben sichtbaren Feldkonfiguration werden beide verwendet. Die vier unnormierten markierten Matrixbeiträge werden gemeinsam summiert. Erst danach folgen die äußere Projektion, die Randpaarung, die Referenznormierung und die verbundene Antwort. Genetische Fitness und das Sucharchiv sind kein Integrationsmaß.

Dieser numerische Beitragspfad bestand 18 gezielte Kontrollen und die unabhängige Gegenlese. Ein vollständig ausgeführter arithmetischer Kontrolllauf prüfte insbesondere zwei verschiedene Atome mit derselben sichtbaren Konfiguration, komplexe Gewichte und eine nichtorthogonale Gramform. Seine Zahlen sind Softwarekontrollen und keine physikalischen TFPT Ergebnisse.

## Was der GA Test entschieden hat

| Beobachtung | Ergebnis |
|---|---:|
| Verschiedene Feldkonfigurationen im deklarierten Raum | 27 |
| Erfüllte implementierte Strukturbedingungen | 27 von 27 |
| Größter numerischer Rest dieser Bedingungen | 0 |
| Verschiedene Zielfunktionssignaturen | 1 |
| Gemeinsame Signatur | `[0, 1]` |
| Verfügbare ursprüngliche Quellenbeiträge | 0 von 27 |

Die erste Komponente der Signatur zählt die bekannten Strukturverletzungen. Die zweite zeigt den fehlenden Quellenbeitrag an. Da beide für alle Kandidaten gleich sind, kann die bestehende Auswahl keinen davon bevorzugen. Die unterschiedlichen Matrizen und ihre unterschiedliche Modenmischung sind reale Unterschiede, aber noch kein ursprüngliches physikalisches Bewertungskriterium.

Für diese kleine Klasse ist die vollständige Aufzählung deshalb die passende Methode. Ein längerer genetischer Lauf mit derselben Bewertung würde den fehlenden Quellenaufruf nicht erzeugen: Die Gene verändern hier Feldkoeffizienten und keine Programmkörper.

## Wo der tatsächliche vollständige Start anhält

Der Integrationsstart wurde auch von der Hauptinstanz tatsächlich aufgerufen. Er endet mit Exitcode 2 und `ORIGINAL_SOURCE_INTEGRATION_NOT_STARTED`.

Die noch auszuführende ursprüngliche Operation erhält im neuen Programm diesen Schnittstellennamen:

```text
original_blocked_fiber_packet(V_n, x, endpoint_words, ordered_times)
```

Sie muss für dieselbe ursprüngliche Historie die tatsächliche Maßmasse sowie die vier vollständigen unnormierten Quellenbeiträge liefern. Dieser Programmname steht für die bereits verlangte ursprüngliche Integration und ist kein zusätzliches physikalisches Postulat. Ein tatsächlich ausgewerteter identischer Randquellenaufruf könnte das fundamentale Ziel auch direkt erfüllen. Die vorhandene Fermionrechnung berechnet eine Antwort bei vorgegebenem Hintergrund. Das bestimmt noch nicht das gemeinsame Gewicht dieses Hintergrunds und aller übrigen Felder. Ebenso berechnet der räumliche Einteilchenoperator keinen gemeinsamen dynamischen Hamiltonoperator.

Als nicht implementiert sind außerdem die fünf zugehörigen Gesamtprüfungen ausgewiesen: gemeinsame markierte Quelle, ursprüngliche Reflexionspositivität, Quanten Ward und Schwinger Dyson Kontakte, ganzer physischer Selektor samt Referenz und Eta sowie die ursprüngliche Grenzfolge. Das Programm verlangt keine bereits bestandenen Messungen vor Beginn einer Integration. Es verlangt, dass die Beiträge und Prüfungen überhaupt ausführbar implementiert sind.

Die Originalquellen wurden an dieser Stelle erneut rückwärts verfolgt. Das proper time Profilmaß, die C6 Positivität, die Labelprojektoren und die vorhandenen Verkettungsregeln behalten ihre Ergebnisse. Sie liefern in den geprüften ursprünglichen Aufrufen den fehlenden gemeinsamen numerischen Faserbeitrag nicht bereits als ausgeführte Rechnung. Der direkte Randquellenweg bleibt zulässig. Der vorhandene lineare quasifreie Kontaktweg allein erzeugt jedoch skalare Higgskontakte; ein anderer passender vollständiger Randquellenwert wurde hier nicht gewonnen.

## Wo ein genetisches Verfahren sinnvoll werden kann

Für einen größeren Kandidatenraum ist GA als äußere Suche nach Gegenfällen oder nach besonders aussagekräftigen Konfigurationen interessant, sobald dafür eine tatsächliche ursprüngliche Antwort berechnet wird. Die bereits festgelegten Kopplungen bleiben dabei erhalten. Eine neue Zielfunktion darf die fehlende Physik nicht durch passend optimierte Zahlen ersetzen.

Für die gewichtete Integration ist ein populationsbasiertes Verfahren mit korrekter Gewichtung die passendere Anschlussidee. Population Annealing verwendet Populationen und Vermehrung, berechnet diese aber aus den tatsächlichen statistischen Gewichten. Die Primärarbeiten behandeln insbesondere positive Gibbsverteilungen und Fehlerkontrolle in Spinmodellen. Ihre Übertragung auf die ursprüngliche TFPT Quelle wäre erst mit deren tatsächlichem Gewicht und passender Behandlung ihrer Phasen begründet. Sie ist hier noch nicht implementiert. Siehe [Machta 2010, Version 2](https://arxiv.org/abs/1006.0252v2) und [Wang, Machta und Katzgraber 2015, Version 2](https://arxiv.org/abs/1508.05647v2).

Es gibt auch eine frühe Arbeit über eine Kombination genetischer Schritte mit Metropolis für zweidimensionale reine SU(2) Gittertheorie. Sie behandelt ausdrücklich die korrekte Gleichgewichtsgewichtung und ihre Grenzen. Das ist ein Beispiel für die mögliche Verbindung, kein Nachweis einer vollständigen TFPT Simulation. Siehe [Yamaguchi 1998, Version 2](https://arxiv.org/abs/hep-lat/9808001v2).

Die aktuelle Entscheidung lautet deshalb: den genetischen Suchmodus erhalten, die kleine erklärte Klasse als ausgewertet sichern und den nächsten Entwicklungsaufwand auf die tatsächliche ursprüngliche Beitragsberechnung richten. Ein großer physikalischer Lauf wurde nicht als gestartet ausgegeben.
