# IDEAS: operative Ursprungsprinzipien für TFPT

Datum: 20. November 2026. Modus: Exploration.
Geprüfter Ausgangsstand: `a1273187b82488d34d2a474470152a41d42d4055`.

## Ergebnis vorweg

Keine vollständige Lösung. Zwei präzise bedingte Kandidaten erreichen Stufe 1:
minimale perfekte Quantenfehlerkorrektur und eine isotrope Korrespondenz zwischen
Observablen und Dynamik. Ein dritter Ansatz erzeugt die algebraische Form von P2,
verlegt die entscheidende Sektorasymmetrie jedoch in ein neues Postulat. Die
übrigen Ideen werden durch kleine Auswahltests eingeschränkt.

Kein Ansatz beginnt mit P1, P2, vier Marken, D5+A3, E8, einem TFPT Clocknetz oder
einer TFPT Zahlkonstante. Bekannte TFPT Strukturen sind Vergleichsziele. Die
wohlbekannte Mathematik der Quantenfehlerkorrektur wird nicht als neue Mathematik
beansprucht. Neu ist hier ihre ausdrücklich geprüfte Eignung als Ursprung und
die Prüfung konkreter Abbildungsbrücken.

## Übersicht

| ID | Gebiet | Grundprinzip | Erhoffte TFPT Struktur | Entscheidende Messgrößen | Erste freie Struktur und freie Parameter |
|---|---|---|---|---|---|
| A | Quanteninformation | Ein unbekannter Zustand bleibt unter einem beliebigen lokalen Fehler rekonstruierbar, bei minimalem Träger und perfekter Packung | Fünf Plätze, Zugangsschwelle drei gegenüber zwei, eventuell Materie und E8 über Fehlersymmetrien | Mindestlänge, Alphabetgröße, Distanz, Syndrome, Gruppenwirkung, Intertwiner | Tensorprodukt Quantenmechanik; ein logisches System; ein lokaler Fehler; Nichtentartung und perfekte Packung. Kein Zahlenfit, Zustand und Dynamik frei. |
| B | Informationsgeometrie und Jordan Algebra | Jede Observablendifferenz bezeichnet einen reversiblen Generator, dessen alternierende bilineare Wirkung vollständig isotrop ist | Komplexer elementarer Quantenträger und eine Signatur 3+1 | Dimension rotationsäquivarianten Klammerraums, Determinantenform, Boostrelation | Nichtklassischer homogener selbstdualer Kegel vom Rang zwei, volle SO(d) Symmetrie, besondere Dynamikkorrespondenz. Skala kappa, Orientierung und Zustand offen. |
| C | Renormierung und Operatoralgebren | Die nach Informationsverlust verbleibenden Observablen besitzen den kleinsten treuen Träger mit zwei nichtisomorphen nichtkommutativen einfachen Sektoren | Markierte algebraische Zerlegung 3+2, interne Eichalgebra | Minimale Blockgrößen, Fixraumrang, Kanalpositivität, Ablation jeder Sektorforderung | Zwei Sektoren und deren Nichtisomorphie sind gesetzt; bedingte Erwartung oder geeigneter bistochastischer Kanal; Zeit und physische Ladungsnormierung frei. |
| D | Reversible Berechnung und kausale Umschreibungen | Die kleinste reversible lokale Regel mit echter Nichtlinearität erzeugt Geschichte durch ihre erlaubten Anwendungen | Zeitordnung, Wechselwirkung, später kausale Dimension und Materiesektoren | Alle lokalen Permutationen, Anteil nichtaffiner Regeln, Regelkonfluenz, Wachstumsdimension | Binäres Alphabet, reversible lokale Gatter und eine Zusammensetzungsregel. Unter 38976 nichtaffinen Dreibitgattern ist noch keine Regel ausgewählt. |
| E | Kausalmengen und kritisches Wachstum | Eine nur von Vorfahren abhängige Ereignisgeburt ist unabhängig von der Reihenfolge unverbundener Geburten und liegt am kritischen Verzweigungspunkt | Zeitrichtung und Raumzeit ohne eingebettetes Gitter | Pfadunabhängigkeit, Verstöße gegen kausale Kovarianz, Verzweigungsrate, Intervallstatistik und Spektraldimension | Partielle Ordnung, Markov Wachstum und lokaler Gewichtstyp. Selbst das kleine kritische Modell p0=p2=a, p1=1-2a enthält a. |
| F | Selbstreferenz und algorithmische Information | Das Naturgesetz ist die kürzeste Beschreibung, die ihre eigene Auswertung konsistent beschreibt | Gemeinsame Auswahl von Regel, Zustand und endlichen Strukturen | Beschreibungslängen und Stabilität ihres Rankings bei gleichwertigen Universalinterpretern | Universalinterpreter, Kodierung, Maß und Bedeutung von Konsistenz. Ohne diese gibt es kein wohldefiniertes Ranking. |
| G | Kategorien und topologische Zustandssummen | Alle lokalen Zerlegungen desselben Prozesses liefern dieselbe Amplitude; Umtriangulieren ist eine exakte Redundanz | Positiver Randkern und diskreter Weg zum Feldintegral | Pentagon und Pachner Identitäten, Positivität, Sektoren, lokale propagierende Moden | Kategorie, Dimension der Verklebung, Objektalphabet und Gewichtung. Diese Daten werden durch Konsistenz allein nicht ausgewählt. |

## A. Information muss sich mit minimalem Aufwand gegen einen lokalen Fehler erhalten

Ein logisches q Niveau System wird isometrisch in n gleichartige q Niveau
Subsysteme kodiert. Kein vorhandenes Fünferobjekt wird eingesetzt. Der Schutz
gegen einen beliebigen einzelnen lokalen Fehler verlangt Distanz d mindestens
drei. Die nicht assistierte Quantum Singleton Schranke lautet

    n - 1 >= 2(d - 1) >= 4.

Somit n mindestens fünf. Die explizite Konstruktion in `qec_origin.py` erreicht
dieses Minimum. Nichtentartete perfekte Packung lautet

    q**n = q * (1 + n*(q**2 - 1)).

Am unabhängig bestimmten Minimum faktorisiert die Gleichung zu

    (q**2 - 1)*(q**2 - 4) = 0.

Für ganze q mindestens zwei ist q gleich zwei erzwungen. Fünf Plätze und das
binäre Alphabet sind in dieser Modellklasse Ausgaben. Die Quantenmechanik, das
Fehlermodell und perfekte Packung sind Eingaben.

### Test und Abbruchkriterien

Der Code muss alle Fehler mit Gewicht höchstens eins exakt korrigieren, jedes
Paar muss den logischen Zustand vollständig verbergen und jedes Dreierkomplement
muss Rekonstruktion erlauben. Anschließend muss eine vorgeschlagene P2 oder
Materieabbildung die natürliche Symmetrie respektieren. Eine Dimensionsgleichheit
allein reicht nicht.

Ergebnis: Alle Korrekturtests bestehen, aber eine intrinsische Markierung 3+2
scheitert an der transitiven C5 Symmetrie. Die natürliche A5 Syndromdarstellung
hat Orbits 1+15; gerade Teilmengen der fünf Plätze haben 1+10+5. Ein linearer
Intertwiner über F2 existiert nicht. Auch phasenkorrigierte monomiale Lifts über C
scheitern an der Charakterspur eines Elements der Ordnung drei. Siehe VERDICT.

### Was erhalten bleibt

Fünf Plätze sind echt selektiert. Die vier Syndrombits ergeben eine Spreadgeometrie
mit fünf Geraden. Ihr kanonischer perfekter Gruppenkern ist die alternierende
Gruppe A5. Das ist ein konkreter möglicher Anschluss an projektive Symmetrien,
aber noch keine physische binäre Ikosaederwirkung und kein E8 Ursprung.

Quellen: [Laflamme et al.](https://arxiv.org/abs/quant-ph/9602019),
[Cleve, Gottesman, Lo](https://arxiv.org/abs/quant-ph/9901025),
[Grassl, Huber, Winter](https://arxiv.org/abs/2010.07902).

## B. Die Dynamik eines elementaren Systems muss in seinen eigenen Observablen leben

Starte mit einem noch nicht dimensionierten Rang zwei Spin Factor
V = R plus R**d. Fordere eine nichtverschwindende alternierende bilineare
Abbildung von zwei Observablendifferenzen zu einer Observablendifferenz, die
unter allen orientierungserhaltenden Rotationen äquivariant ist:

    Hom_SO(d)(Lambda**2 R**d, R**d).

Dieser Raum hat genau bei d gleich drei Dimension eins und sonst Dimension null.
Der Überlebende ist die Kreuzproduktklammer bis auf einen Skalenfaktor. Der
symbolische Test enthält d von zwei bis acht; ein Paritätsbeweis in VERDICT deckt
alle d ab. Das daraus entstehende Herm2(C) besitzt die Determinantenform
t**2-x**2-y**2-z**2.

Das ist eine bedingte Auswahl eines 3+1 Antwortkegels. Eine physische Ereignismenge,
ein lokales Observablennetz, Zustandsauswahl und ein Kontinuum sind damit nicht
gegeben. Nichtunitäre Kongruenzen H -> B H B* liefern Lorentztransformationen
dieses Kegels, aber im Allgemeinen keine deterministischen spurtreuen Kanäle.

### Negativkontrollen

Unter O(d) statt SO(d) verschwindet die Klammer auch bei d gleich drei.
Orientierung wurde also nicht abgeleitet. Ohne besondere Dynamikkorrespondenz
existiert der positive Lorentzkegel für jedes d. Gewöhnliche lokale Tomographie
allein selektiert die gewünschte Algebra nicht.

Für endliche hermitesche K und P erzwingt i[K,P]=cP mit reellem c ungleich null
P gleich null. Sogar der relative Fehler in voller Hilbert Schmidt Norm ist
mindestens eins. Eine endliche Clock kann deshalb nicht einfach als exakter
unitärer Horizontboost umbenannt werden. Sektorweise Grenzwerte bleiben offen.

Quellen und Abgrenzung: [Müller und Masanes](https://arxiv.org/abs/1206.0630)
setzen Raum mit unbekannter Dimension voraus; ihr Satz erzeugt keine Raumzeit
aus dem Nichts. [Barnum und Wilce](https://arxiv.org/abs/1202.4513) setzen unter
anderem bereits einen komplexen Qubit voraus. Die hier betrachtete zusätzliche
Klammerforderung darf nicht als Folgerung bloßer Tomographie dargestellt werden.

## C. Die Physik liegt in den Observablen, die eine Informationsreduktion überleben

Eine exakte bedingte Erwartung E mit E**2=E bildet M_N(C) auf eine Unteralgebra
der dauerhaft zugänglichen Observablen ab. Als Minimalhypothese verlangt man
zwei nichtisomorphe nichtkommutative einfache Sektoren. Die kleinsten zulässigen
Matrixgrößen sind dann zwei und drei. Der treue minimale Träger hat Dimension
fünf, seine zentralen Projektoren markieren die beiden Summanden.

Dies ist ein sehr kurzer algebraischer Weg zu einer P2 ähnlichen Struktur. Die
Ablation zeigt aber, wo sie herkommt:

| Zusatzbedingungen bei zwei Summanden | Kleinster Träger |
|---|---|
| Keine | C plus C |
| Beide nichtkommutativ | M2 plus M2 |
| Nichtisomorph | C plus M2 |
| Beide nichtkommutativ und nichtisomorph | M2 plus M3 |

Die Spur null Bedingung im Zentrum liefert lediglich ein Ladungsverhältnis;
die physische Hyperladungsnormierung ist nicht erklärt. Der ausdrücklich
konstruierte Kanal E(X)=PXP+QXQ ist positiv, spurtreu, unital und idempotent.
Sein Fixraum hat Dimension 13. Die Krausprojektoren sind eine Realisierung der
ausgewählten Algebra, keine Herleitung der Auswahl aus Dynamik.

Ein beliebiger unitaler CP Projektor genügt nicht: Sein Bild muss unter dem
ursprünglichen Produkt keine Unteralgebra sein. `sector_origin.py` enthält ein
explizites Gegenbeispiel. Auch die gleichmäßige Mittelung über alle zehn 3+2
Markierungen erzeugt kein markiertes Y: Mittelwert null, zweites Moment I/6.

Urteil: OFFEN trotz exaktem algebraischem Teiltreffer. Die Nichtisomorphie braucht
eine unabhängige physikalische Begründung. Bisher wird die Asymmetrie postuliert.

## D. Reversible Nichtlinearität als kleinster Prozessbaustein

Geometrie entsteht erst aus Abhängigkeiten der Regelanwendungen. Vorher gibt es
nur Zustände und reversible lokale Gatter. Die vollständige minimale Prüfung
zählt alle Permutationen von einem, zwei und drei Bits. Bei zwei Bits sind alle
24 reversiblen Gatter affin über F2. Bei drei Bits sind von 40320 Gattern genau
1344 affin und 38976 nichtaffin.

Die erste echte Nichtlinearität braucht in dieser Klasse drei Bits. Diese Drei
ist eine Gatterarität, keine Familienzahl und keine Raumdimension. Der Ansatz
braucht zusätzlich eine Regel für überlappende Anwendungen, kausale Unabhängigkeit
und Quantisierung. Eine Toffoliregel oder ein bestimmter Hypergraph wird hier
nicht passend ausgewählt. OFFEN, Stufe 0.

## E. Kritisches kausales Wachstum als Erzeuger von Zeit und Kontinuum

Geburtsschritte hängen ausschließlich vom bereits vorhandenen Vorfahrenraum ab.
Kommutierende Geburten müssen in jeder Reihenfolge dieselbe Gesamtamplitude
ergeben. Zusätzlich soll weder Aussterben noch explosives Wachstum bevorzugt
werden. Das klingt wie ein selektierendes Minimalprinzip, bestimmt aber nicht
einmal die einfachste Verzweigungsverteilung eindeutig:

    p0 = p2 = a, p1 = 1-2*a, 0 <= a <= 1/2.

Für jedes a ist die Summe eins und der mittlere Nachwuchs eins. Kritikalität
allein beseitigt den freien Parameter nicht. Eine Quantenfassung müsste zudem
Interferenz, eine positive Zustandsform und ein Kontinuum mit stabiler
Intervallstatistik liefern. Keine Raumdimension wird hier vorausgesetzt oder
behauptet. OFFEN, Stufe 0.

[Rideout und Sorkin](https://arxiv.org/abs/gr-qc/9904062) zeigen bereits, dass
Kovarianz und bestimmte Kausalitätsforderungen eine Familie klassischer
Wachstumsdynamiken erlauben. Der neue Anspruch wäre eine belegte Auswahl aus
dieser Familie, nicht die bloße Existenz kausaler Geburt.

## F. Der kürzeste selbstkonsistente Naturcode

Eine Beschreibung soll ihre eigene Auswertung oder Grobreduktion reproduzieren.
Der kürzeste solche Fixpunkt soll die Welt auswählen. Die rein sprachfreie
Version scheitert bereits vor einer TFPT Rechnung: Für jede berechenbare
Kandidatentheorie T und einen universellen Präfixinterpreter U konstruiert man

    U_T(0) = T,
    U_T(1p) = U(p).

U_T ist universell und präfixfrei, und T hat die eindeutige Einbitbeschreibung.
Somit kann jede berechenbare Theorie durch eine legitime Wahl des Interpreters
zum kürzesten Kandidaten werden. Die Kosten des Interpreters mitzuzählen braucht
wieder eine Referenzkodierung. Das ist kein Ausschluss aller Selbstreferenzideen;
es widerlegt nur eine eindeutige Auswahl durch unbestimmte Beschreibungslänge.

VERWORFEN in dieser präzisen Form, Stufe 0. Ein physisch definierter
Informationsverlust oder eine rekonstruierbare universelle Kostenfunktion
würde die Voraussetzungen verändern.

## G. Konsistenz aller Zerlegungen als Naturgesetz

Eine Amplitude wird aus elementaren Prozessstücken zusammengesetzt. Jede
zulässige lokale Umzerlegung muss dieselbe Antwort liefern. Die resultierenden
Pentagon und Pachner Gleichungen sind starke algebraische Bedingungen, die
potenziell einen positiven Randkern liefern können. Sie wählen jedoch die
Dimension, die Kategorie und die Gewichte nicht allein aus. Bereits die triviale
Kategorie und verschiedene Gruppenkategorien liefern unterschiedliche
konsistente Antworten.

Ein rein topologischer endlicher Zustandssummentyp garantiert auch keine lokalen
propagierenden Gravitationsfreiheiten. Eine kritische Erweiterung mit lokalen
Observablen wäre ein neuer Ansatz und müsste gesondert geprüft werden. OFFEN,
Stufe 0. Diese konkrete Form wurde im durchsuchten TFPT Katalog nicht gefunden.

Primärquelle: [Barrett und Westbury, Invariants of piecewise linear 3 manifolds](https://www.math.ru.nl/~mueger/TQFT/BW2.pdf).

## Ausgesonderte Wiederholungen

Die Auswahl von E8 als kleinster nichttrivialer bosonischer invertierbarer
chiraler Phase ist bereits in `F1-E8-ANOMALY-INFLOW-ORIGIN-20261116.md` vorhanden.
Der spezielle Frobenius 1+i über F2 gehört ebenfalls zum Bestand. Beide sind
deshalb keine neuen Ansätze in dieser Sitzung. Der bekannte Lorentzkegel selbst
ist ebenso kein neuer Treffer; neu geprüft ist die zusätzliche Klammerauswahl.

## WILDCARD: Raumzeit als Reihenfolge von Wiederherstellbarkeit

Vielleicht sind die elementaren Ereignisse keine Veränderungen eines Zustands
auf einer schon vorhandenen Bühne. Ein Ereignis könnte vielmehr bedeuten, dass
eine zuvor nicht rekonstruierbare Unterscheidung erstmals rekonstruierbar wird.
Definiere A vor B, wenn B die Information von A wiederherstellen kann, die
Umkehrung aber unmöglich ist. Nach Quotientierung wechselseitiger Rekonstruktion
entsteht eine partielle Informationsordnung.

Eine Raumzeit wäre dann ein großer konsistenter Zusammenhang solcher
Wiederherstellungsbeziehungen. P1 wäre ein Übergang zwischen unvereinbaren
Rekonstruktionsregimen, P2 ein minimaler persistenter Zugangstyp. Der Ansatz
berührt Fehlerschutz, modulare Theorie und holographische Rekonstruktion, setzt
aber keines der TFPT Randobjekte voraus.

Die abwegige Möglichkeit: Die sichtbare Zeitrichtung wäre nicht die Richtung,
in der etwas gerechnet wird, sondern die Richtung, in der Zugriff verloren
geht oder wieder gewonnen wird. Noch fehlt die Auswahl der Kanalmenge und des
gemeinsamen Zustands. Stufe 0, OFFEN. Ein erster entscheidender Test wäre eine
nichttriviale halbseitige modulare Inklusion einer aus dieser Ordnung erzeugten
unendlichen Algebra, ohne eine geometrische Zeit von außen einzusetzen.
