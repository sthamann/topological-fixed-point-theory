# OriginLab 2.0 · Planabnahme

Stand 14. November 2026. TFPT Fundamentalplan 361, Atlasrevision
`32bb9da024c4595c2a8380f7`. Maßstab sind alle 13 Abschnitte
der OriginLab 2.0 Spezifikation.
Die Abnahme bezieht sich auf den kumulativen Patch und den beschriebenen
Prüflauf. Status des kanonischen Forschungsplans bleibt F1 offen.

Ergänzung vom selben Tag. [Quellenabgleich](SOURCE-BRIDGE-2026-11-14.md),
Atlasrevision `dace8ddec7ef89c319e6bfed`. Im zweiten geprüften Lauf
standen cvc5, egglog und FLINT zur Verfügung. Alle 59 Tests bestanden.
Der native R113 Aufruf meldet weiterhin den fehlenden Quellenproduzenten.

| Plan | Ausgeführt oder implementiert | Offene Abnahme |
| --- | --- | --- |
| §1 Stand und §2 Infrastruktur | Revisionsgebundene Eingabe, Quellhashes und getrennte Ausführung für Original und Modelle. Direkter R113 Aufruf und Originalprüfer R005, R006, R069, R104. | Der Atlas enthält auch wissenschaftliche Aussagen ohne ausführbare Eingabeverträge. Die neuen Änderungen liegen im Patch. |
| §3 Architektur | Originalquellenleser, endliche Regelsuche, exakte Verifikation, Gegenbeispiele, Beobachtungsplanung, lokale Queue und revisionsgebundene Review Übergabe. | Keine durchgängig ausführbare physikalische Beweiskette. |
| §4 Werkzeuge | SymPy und Z3 eingesetzt. Lean prüft den R102 Rangschluss (n mid 4). Pawl bestätigt B01 bis B09 als technische Abnahme. cvc5, egglog und FLINT sind im bestehenden Code angebunden. | cvc5, egglog und python-flint fehlen in dieser Laufumgebung. Aesop, quimb, cotengra, OQuPy, CVXPY/SCS, LeanDojo, PySR und Ray wurden nicht angebunden. Allgemeine Gleichheitssättigung und zertifizierte physikalische Integrale sind offen. |
| §5 R113 | Originalverbraucher ausgeführt. Gleichzeitige Kontraktion der vier Marken bei gelieferten Formen. Neue Gauß Referenz berechnet vier gemeinsame Ableitungen inklusive direktem Kontakt und singulärer Matrix. Der rechteckige Löser prüft gelieferte Cartan, Farb, Zeit und Stromoperatoren mit Adjungierten. | Ursprünglicher Produzent für gemeinsame gewichtete Randform und absolute Referenz fehlt. Endliche Berezin und Schur Routinen besitzen nur Teilinputs. Vollständige Wechselwirkungen, unabhängige numerische Integration und Fehlergrenzen fehlen. |
| §6 Sprache und Suche | Typisierte endliche Prozessgrammatik, Vorwärtssuche, exakte Matrixsemantik, Gegenbeispielreplay, Mutation, deklarierte Modellvorschläge, exakte inverse Teilaufgaben und Minimalitätskriterien. | Keine allgemeine präphysikalische Grammatik mit emergenten `Event`, `State`, `Clock`, `Source` und `Geometry`. Kein vollständiger Quotient beliebiger gradierter Äquivalenzen und keine automatische Theoriebildung über alle Kandidatenklassen. |
| §7 Aktive Suche | R006 findet neue Wörter der Tiefe fünf. R104 prüft kohärenten Speicher und Phasentiefe acht. Komplexe R113 Kontrolle trennt (+i/4) und (-i/4). Die Viertelclock Planung entscheidet die 30 endlichen Varianten in höchstens drei Abfragen. | Zulässigkeit der geplanten Wörter auf der physikalischen Quelle und Auswertung mit derselben Referenz. Keine allgemeinen Suchgarantien. |
| §8 G0 bis G7 | Der F1 Bericht zeigt belegte Eingaben, eng begrenzte Ausschlüsse, offene Gatter und die Blockade von G3 und G4. Ersatzwerte bleiben ausgeklammert. | Die vollständigen physikalischen Prüfer für Ursprung, Dynamik, Emergenz, Kontinuum und empirischen Vergleich sind offen. Die Ausführung prüft keine unbekannten physikalischen Beweise automatisch. |
| §9 Benchmarks | B01 bis B09 gemeinsam ausgeführt. R005 rekonstruiert eine Abbildung aus 2.048 Klassen, R006 prüft je 2.245 Momente und entdeckt das kürzeste Trennwort, R069 führt 53 Originalkontrollen aus. B07 bestätigt die erwartete Blockade, B09 hält F1 für eine Modellrechnung offen. | Die übergreifende Blindtest Sammlung sowie die vollständigen Manipulationstests an Projektor, Referenz, Integrationsgrenzen und Spinlift fehlen. B07 braucht nach einem echten Quellenanschluss einen neuen positiven Referenzfall. |
| §10 Betrieb | Lokale Queue mit Worker Wiederanlauf, Manifesten, Negativreplay und Review Übergabe. Pawl Ausführungsgate für die gemeinsame Benchmark Suite. | Keine verteilte Ausführung, kein vollautomatisches Atlas Review und keine vollständige Adaption der Pawl Forschungssteuerung. |
| §11 Forschungspfade A bis D | A liefert einen reproduzierten Originalblocker und genaue Anschlussbedingungen. B prüft endliche Algebra und Glues. C prüft positive Transfer und Fortsetzungsmodelle. D prüft gelieferte regionale Algebraaktionen, gemeinsame Spektren und endliche Bewegung. | A braucht den Originalproduzenten. B braucht eine primitive Regel, die P1/P2 und Ladung erzeugt. C braucht eine Herleitung der Übergangsamplituden. D braucht drei tatsächliche unabhängige Raumrichtungen und eine kontrollierte (3+1) Rekonstruktion. |
| §12 Umsetzung und §13 Metriken | Konsolidierte Suite und pro Lauf überprüfbare Belege, Kandidatenzahl, Gegenbeispiele, markierte Quellwerte und offene Ableitungen. | Keine automatische vollständige Messung aller verbleibenden primitiven Freiheitsgrade, keine vollständige Beweiskette und kein positiver R113 Lauf. |

## Die früheste physikalische Lücke

Der kanonische Plan priorisiert echte regionale Trägerbewegung. R113 und
diese regionale Frage können bis zur gemeinsamen Quellbindung unabhängig
untersucht werden. Der vorhandene Seitentausch erhält das Produkt, bewegt
aber keine Trägerregion. Benötigt werden zwei zusätzliche unabhängige
Originalschnittbewegungen neben der Kreisrichtung, mit demselben
positiven Zeittransport und Zustandskegel.

Für R113 fehlen weiterhin das gemeinsame Feldmaß und die gewichtete
Randform mit den geladenen Zwischenwegen vor dem äußeren Selektor, der
Kappe und derselben absoluten Leerreferenz. Der nächste überprüfbare
Anschluss ist eine originale Faser oder eine äquivalente integrierte
Randform. Die kleinste positive Abnahme berechnet (Z_0,Z_U,Z_Y,Z_{UY})
samt direkten Kontakten aus einem gemeinsamen Aufruf und prüft ein
zurückgehaltenes geordnetes Zeitrechteck. Ein negativer Test muss die
getestete Klasse und ihre Voraussetzungen benennen.

Für F1 braucht es darüber hinaus die gemeinsame Herkunft von P1/P2 und
der primitiven Ladung, eine kausale lokale (3+1) Theorie, chirale
Materie, einen wechselwirkenden universellen masselosen Spin 2 Sektor
sowie kontrollierte Grenzübergänge und unabhängige empirische Tests.
Keiner der neun bestandenen Benchmarks schließt diese Beweispflichten.
