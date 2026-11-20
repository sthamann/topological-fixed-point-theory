# VERDICT: exakt geprüfte Ursprungsansätze

20. November 2026. Ausgangsstand `a1273187b82488d34d2a474470152a41d42d4055`.

## Wissenschaftlicher Status

**Keine vollständige TFPT Lösung. Kein Ansatz erreicht Stufe 2.**
Ein Vorsprung vor sämtlichen bisherigen Wegen ist nicht bewiesen. Neu sind
konkrete bedingte Strukturtests und nachgerechnete Hürden, die bisherige
Dimensionsähnlichkeiten durch überprüfbare Abbildungskriterien ersetzen.

Stufe 0 bedeutet hier: keine der vom Auftrag definierten Stufen erreicht.
Stufe 1 ist ein mathematischer Teiltreffer unter offen genannten Voraussetzungen,
keine physikalische Bestätigung.

| ID | Urteil | Stufe | Ergebnis | Stärkster Einwand | Erste freie Wahl |
|---|---|---:|---|---|---|
| A | KANDIDAT für die Fünferstruktur | 1 | n=5 und q=2, exakte Korrektur, Zugangsschwelle, kanonischer A5 Kern | Keine markierte 3+2 Auswahl; natürliche Syndromdarstellung ist nicht die Materiedarstellung | Tensorprodukt Quantenmechanik, ein logisches System, ein lokaler Fehler, perfekte Packung |
| B | KANDIDAT für einen 3+1 Antwortkegel | 1 | Nichttriviale SO(d) Klammer nur bei d=3 | Besonderes Dynamikpostulat; Antwortkegel ist keine Raumzeit; Orientierung vorausgesetzt | Rang zwei HSD Kegel und volle SO(d) Kovarianz der Klammer |
| C | OFFEN | 1, nur algebraisch | M2 plus M3 und Trägerdimension fünf nach zwei Sektorpostulaten | Nichtisomorphie setzt die entscheidende Asymmetrie; physische P2 Brücke fehlt | Genau zwei ungleiche nichtkommutative einfache Sektoren |
| D | OFFEN | 0 | Erste nichtaffine reversible binäre Gatter haben Arity drei | Kein Gatter und kein globaler kausaler Prozess ausgewählt | Binäre reversible lokale Berechnung |
| E | OFFEN | 0 | Kritisches Wachstum lässt freie Gewichte übrig | Weder eindeutige Amplitude noch Quantisierung und Dimension | Markov Wachstum partieller Ordnungen |
| F | VERWORFEN in der sprachfreien Form | 0 | Universalinterpreter kann jeden Kandidaten zum kürzesten machen | Keine referenzunabhängige Auswahl durch Beschreibungslänge allein | Interpreter und Kodierung |
| G | OFFEN | 0 | Umzerlegungskonsistenz ist kompatibel mit vielen Zustandssummen | Kategorie, Dimension und lokale Dynamik werden nicht selektiert | Typ der Verklebung und Kategorie |

## 1. Positive Resultate aus der Ausführung

### A. Minimaler perfekter Quantenspeicher

`qec_origin.py` verwendet ausschließlich exakte Brüche und symbolische
Faktorisierung. Die Singleton Schranke ist ein analytischer Eingangssatz und
wird nicht aus einem endlichen Suchfenster erschlossen.

    n - 1 >= 2(d - 1), d >= 3  =>  n >= 5
    q**4 - 5*q**2 + 4 = (q-2)*(q-1)*(q+1)*(q+2)
    q ganz, q>=2 => q=2.

Tatsächlich ausgeführte Ergebnisse:

| Test | Ergebnis |
|---|---:|
| Stabilisatorrang | 4 |
| Distanz | 3 |
| Logische Dimension | 2 |
| Fehleralternativen | 16 |
| Orthogonalitätsbedingungen der 32 Fehlerbilder | 1024 von 1024 |
| Paare ohne logische Information | 10 von 10 |
| Dreiergruppen mit Rekonstruktion | 10 von 10 |
| Maximale Dreiplatzmischung des sechsteiligen Choi Zustands | 20 von 20 |
| Dimensionen der physischen Paritätssektoren | 16 und 16 |
| Geschützter Code im geraden Paritätssektor | 1 |

Die 16 Syndrome, der 16 dimensionale Paritätssektor und eine 16 dimensionale
chirale Materiedarstellung sind verschiedene Objekte. Insbesondere ist der
geschützte Quantenzustand zweidimensional.

### Syndromsymmetrie

Die drei Einzelplatzfehler bilden je eine projektive Gerade in F2**4. Die fünf
Geraden zerlegen alle 15 Nichtnullvektoren. `syndrome_spread.py` prüft alle
20160 invertierbaren linearen Abbildungen. Die Erhaltungsgruppe hat Ordnung 360,
Bild S5 auf den fünf Geraden und geradenweisen Kern C3. Die direkt berechnete
abgeleitete Reihe hat die Ordnungen

    360 -> 180 -> 60 -> 60.

Der terminale perfekte Kern wirkt treu als alle geraden Permutationen, also als
die alternierende Gruppe A5. Das ist nicht die Lie Algebra A5. Reine physische
Platzpermutationen erhalten den Code dagegen nur in einer Gruppe der Ordnung
zehn, D5. Ein unabhängiger Test berücksichtigt dafür den phasengenauen
32 mal 32 Codeprojektor.

Eine universelle zentrale Erweiterung von A5 und eine anschließende McKay
Konstruktion sind mögliche mathematische Zusatzoperationen. Sie sind in dieser
Rechnung keine physisch abgeleitete Wirkung und kein berechnetes E8 Strommodell.

### B. Auswahl der Dimension eines Antwortkegels

    dim Hom_SO(d)(Lambda**2 R**d, R**d) = 1 für d=3, sonst 0.

Die ausgeführte Tabelle:

| d | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Freie alternierende Koeffizienten vor Kovarianz | 2 | 9 | 24 | 50 | 90 | 147 | 224 |
| Invariante Dimension | 0 | 1 | 0 | 0 | 0 | 0 | 0 |

Beweis über das Suchfenster hinaus: Ein Koeffizient c_ij^k verändert bei einer
Pi Rotation in der ab Ebene sein Vorzeichen entsprechend der Anzahl der
auftretenden a und b Indizes. Alle Indexhäufigkeiten müssen deshalb dieselbe
Parität besitzen. Insgesamt sind es drei Indizes. Alle gerade ist unmöglich;
alle ungerade verlangt d ungerade und d höchstens drei. d=1 scheitert an der
Alternierung. Bei d=3 bleibt der Epsilon Tensor, und die infinitesimale
Rotationskovarianz reduziert seine Koeffizienten auf eine Skala. Jacobi und
Pauli Identitäten sind exakt geprüft.

Für H=tI+x sigma_x+y sigma_y+z sigma_z ergibt sich exakt

    det H = t**2 - x**2 - y**2 - z**2.

Die Transformation mit B=diag(2,1/2) erhält die Determinante exakt. Die daraus
folgenden Zahlen 17/8 und 15/8 sind beliebige rationale Prüfdaten eines Boosts,
keine physikalischen Zahlentreffer.

## 2. Exakte Hürden für die stärksten Ideen

### H1. Zugangsschwelle ist keine markierte 3+2 Zerlegung

Schon C5 wirkt transitiv auf den fünf physischen Plätzen. Eine kanonische
Teilmenge müsste invariant sein. Die einzigen invarianten Teilmengen haben
Größe null und fünf. Eine Markierung von zwei oder drei Plätzen ist deshalb
aus diesem unmarkierten Code allein unmöglich.

Das schließt eine spontane, bedingte oder umgebungsabhängige Markierung nicht
aus. Diese würde aber eine eigene Zustands und Dynamikregel verlangen.

Die gleichmäßige klassische Mittelung über alle zehn möglichen Markierungen
repariert das nicht:

    mean(Y) = 0, mean(Y**2) = I/6.

Die Aussage betrifft diese lineare Mittelung, nicht kohärente konditionierte
Zweige oder eine globale Zustandssumme.

### H2. Die natürliche 16er Syndromdarstellung ist nicht Lambda gerade C5

Die kanonische A5 Wirkung besitzt auf Syndromvektoren die Orbits 1+15. Auf den
geraden Teilmengen der fünf Geraden besitzt dieselbe Gruppenwirkung die Orbits
1+10+5. Die linearen Intertwinergleichungen über F2 haben Rang 16 bei 16
Unbekannten. Ihr Lösungsraum ist null.

Auch über C lässt sich die naive Identifikation nicht durch Phasen retten:

| Elementordnung | Fixierte Syndromlabels | Spur auf Lambda gerade C5 |
|---|---:|---:|
| 1 | 16 | 16 |
| 2 | 4 | 0 |
| 3 | 1 | 4 |
| 5 | 1 | 1 |

Ein monomialer unitärer Lift einer Ordnung drei Syndrompermutation hat höchstens
einen nichtverschwindenden Diagonaleintrag, also Betrag der Spur höchstens eins.
Die Exterioralgebra hat Spur vier. Ein Intertwiner ist auch mit solchen Phasen
unmöglich. Eine andere nichtmonomiale Darstellung oder ein größerer Raum fällt
nicht unter diesen Ausschluss.

### H3. Der perfekte Syndromkern erzeugt keine icosahedrale logische Qubitwirkung

`qec_lift_check.py` hebt alle 360 Syndromautomorphismen in physische lokale
symplektische Abbildungen mit Platzpermutation an. Alle erhalten den
Stabilisatorraum. Das volle logische symplektische Bild hat Ordnung sechs.
Alle 60 Elemente des perfekten Kerns wirken auf logischen Pauli Labels als
Identität. Phasen und Pauli Korrekturen bleiben dabei getrennte Daten.

Eine treue projektive A5 Wirkung lässt sich innerhalb der logischen
Einqubit Cliffordgruppe nicht erhalten. Der Schluss betrifft diesen Lift,
nicht jede denkbare Quantenwirkung von A5.

### H4. 63 ist hier keine erreichbare Gedächtnisdimension

Für jedes Dreiergebiet A lautet die kodierte Operatorabbildung

    F_A(X) = Tr_Ac(V X V*), X in M2.

Sie ist injektiv, da aus A rekonstruiert werden kann. Ihr Rang ist daher vier;
die hermiteschen spurfreien Richtungen haben Dimension drei. Für alle zehn
Dreiergebiete wurde dies exakt berechnet. Bei reinem logischem Eingang lautet
das Spektrum 1/4 vierfach, null vierfach.

Die Zahl 63=4**3-1 zählt die volle spurfreie Dreiqubit Algebra. Sie zählt nicht
die erreichbare kodierte Operatorfamilie. Deshalb ist 68=5+63 nicht reproduziert.
Zusätzliche Fehleraufzeichnungen oder Dynamik könnten den Raum erweitern und
müssten separat gerechnet werden.

### H5. Eine endliche hermitesche Clock ist kein exakter nichttrivialer Boost

Für endliche hermitesche K und P und reelles c gilt durch Zyklizität der Spur

    Tr(P i[K,P]) = 0,
    ||i[K,P]-cP||_HS**2 = ||i[K,P]||_HS**2 + c**2 ||P||_HS**2.

Für c ungleich null ist i[K,P]=cP somit nur bei P=0 möglich. Für P ungleich null
ist der relative Fehler, normiert durch |c| ||P||_HS, mindestens eins. Der
analytische Beweis gilt in jeder endlichen Dimension; zusätzlich wurde die
Identität symbolisch für allgemeine hermitesche 2 mal 2 Matrizen ausgeführt.

Ausgeschlossen ist eine globale Identifikation in dieser Norm. Nicht
ausgeschlossen sind ausgewählte Sektoren, andere Konvergenzbegriffe und
unbeschränkte Operatoren im unendlichen Grenzwert.

### H6. Symmetrischer endlicher Zustand erzeugt keine nichttriviale Modularzeit

Ein unter allen Unitaries invarianter normierter Zustand ist rho=I/n. Dann ist
rho**(is) A rho**(-is)=A. Eine nichttriviale endliche Modularzeit verlangt mehr
Struktur als maximale Symmetrie. Eine Inklusion mit einem nichttrivialen
Zustand in unendlicher Dimension ist eine andere Voraussetzung.

## 3. Basislinien

### Konstruktionen gleicher syntaktischer Komplexität

Der QEC Zensus verwendet alle Pauli Wörter derselben Länge als zyklische
Stabilisatorsaaten. Dies ist keine vollständige Klassifikation aller Codes.

| Plätze | Alle Wörter | Kommutierend und Rang n minus eins | Korrigierend |
|---|---:|---:|---:|
| 1 | 4 | 1 | 0 |
| 2 | 16 | 3 | 0 |
| 3 | 64 | 9 | 0 |
| 4 | 256 | 48 | 0 |
| 5 | 1024 | 75 | 30 |
| 6 | 4096 | 216 | 0 |

Bei n=5 sind 30/1024 gleich 15/512 erfolgreich; bedingt auf kommutierende Saaten
korrekten Rangs sind es 30/75 gleich 2/5. Die zwei Stabilisatorräume bilden in
dieser Klasse eine Äquivalenzklasse unter Platzpermutationen. Die Null bei n=6
ist kein allgemeiner Ausschluss von Sechsercodes. Fünf, 16 und Zugangsschwelle
sind korrelierte Konsequenzen und zählen als ein Paket.

Für den Klammeransatz besitzt bei gleicher Definition unter allen d nur d=3
eine nichttriviale Lösung. Das ist eine strukturelle Eindeutigkeit in der
genannten Klasse, kein statistischer p Wert. Werden O(d) statt SO(d), andere
Darstellungen oder andere Dynamikpostulate erlaubt, ändert sich die Klasse.

Für C enthält die vollständige Zweiblockliste mit Trägergröße höchstens zwölf
36 Paare. 25 haben zwei nichtkommutative Summanden, 20 zusätzlich nichtisomorphe.
Mit beiden Bedingungen und minimalem Träger gewinnt genau (2,3). Alle vier
Teilmengen der zwei Bedingungen wurden mit dem vorhandenen
`genesis_rules.search` abgeglichen. Die einzige zweiklauselige Regel produziert
(2,3), die drei anderen produzieren (1,1), (1,2) beziehungsweise (2,2). Ein
Vorsprung über frei wählbare Klauselvokabulare ist damit nicht belegt.

### Wörtliche Zahlenspielkontrolle mit Ausdrücken

Es wurden keine gemessenen Dezimalzahlen angepasst. Zusätzlich zählt
`arithmetic_baseline.py` jede zulässige Formel gleicher Komplexität aus den
Atomen 1, q=2, n=5 und den Operationen plus, minus, mal, geteilt und Potenz.
Gezählt werden geordnete Syntaxbäume mit genau null, einer oder zwei Operationen;
Kommutationsduplikate und Identitäten bleiben enthalten. Potenzen haben ganze
Exponenten von null bis acht; 0**0 und Division durch null sind ausgeschlossen.
Nach jedem Knoten sind Betrag des Zählers und Nenner auf eine Million begrenzt.
Diese ganze Grammatik und Ressourcenbegrenzung ist offengelegt, nicht universell.

| Operationen | Zulässige Formeln | Ergebnis 3 | Ergebnis 4 | Ergebnis 5 | Ergebnis 16 |
|---|---:|---:|---:|---:|---:|
| 0 | 3 | 0 | 0 | 1 | 0 |
| 1 | 45 | 3 | 4 | 4 | 0 |
| 2 | 1294 | 54 | 70 | 89 | 10 |

Beispielsweise treffen zehn Ausdrücke mit zwei Operationen die 16 exakt. Eine
bloße Rechnung auf 16 wäre somit schwach. Der QEC Treffer besteht in den
Korrekturidentitäten und der Minimalität, nicht in dieser Integerarithmetik.
Die Atome n und q sind hier bereits aus dem bedingten Modell bekannt; diese
Tabelle ist keine unabhängige Wahrscheinlichkeit seiner Entdeckung. Die
Gruppenordnungen werden über vollständige Gruppenwirkungen geprüft, nicht
über passende Zahlenausdrücke.

## 4. Welche TFPT Ergebnisse werden tatsächlich getragen?

| Ziel | Ergebnis dieser Sitzung |
|---|---|
| P1 mit Orientierung, Reflexionspositivität, Einheitswindung, c3 | Von keinem Ansatz vollständig erzeugt |
| P2 als physischer markierter Träger | Nicht erzeugt; A liefert unmarkierte fünf Plätze, C eine bedingt gewählte Algebra |
| Chirale Standardmodellgeneration | Natürliche A Syndromidentifikation widerlegt; andere Abbildung offen |
| Drei Familien und vier Marken | Nicht erzeugt; drei logische Pauli Richtungen sind keine Familien |
| Physisches E8 Strommodell | Nicht erzeugt; A5 eröffnet nur eine zusätzliche mathematische Anschlussfrage |
| Seed, Alpha, Massen und Mischung | Keine neue gemeinsame Herleitung; bestehende Rechenketten nicht erneut getestet |
| Zustand und Zeit | Nicht ausgewählt; endliche symmetrische Modularzeit trivial |
| 3+1 Raumzeit | B wählt nur eine Antwortsignatur; Ereignisse, Metrik und lokale Felder fehlen |
| Kontinuum und Gravitation | Nicht erzeugt; endlicher Boosttest präzisiert eine notwendige Hürde |

Die vorhandenen TFPT Resultate bleiben als Zielkatalog bestehen. Das bedeutet
nicht, dass ein neuer Ursprungsansatz sie bereits trägt. Hier wurden keine
Ledger oder Papierbehauptungen befördert.

## 5. Konkreter nächster Lösungsversuch

Die stärkste noch offene gemeinsame Frage lautet: Kann ein eigenständiger
Prozess Schutz, nichtkommutative Dynamik und eine asymmetrische Zugangsstruktur
gleichzeitig auswählen? Für einen belegten Fortschritt über diese Sitzung
hinaus muss eine Konstruktion mindestens Folgendes erfüllen:

1. Sie bestimmt einen Kanal beziehungsweise eine Amplitudenregel und ihren
   Zustand, ohne eine 3+2 Projektorenwahl in die Eingabe zu schreiben.
2. Ihr Fehler oder Gedächtnisraum besitzt einen überprüfbaren Intertwiner zur
   Materiedarstellung. Der hier widerlegte monomiale Syndromweg reicht nicht.
3. Eine wachsende Algebra hat einen angegebenen lokalen Observablensektor und
   Zustand. Auf diesem muss ein nichttrivialer modularer Boostgrenzwert
   nachgewiesen werden. Eine globale endliche Hilbert Schmidt Näherung reicht
   nach H5 nicht.
4. Erst an diesem einen Prozess werden Einheitswindung und absolute
   Randnormierung c3 geprüft. Frei normierte U1 Richtungen gelten nicht als
   Herleitung von c3.

Dies ist eine präzise offene Aufgabe, keine behauptete Lösung. Ein Rang,
Orbit oder Charaktervergleich, der scheitert, soll den jeweiligen Zweig sofort
stoppen, bevor weitere TFPT Konstanten darauf aufgebaut werden.

## 6. Nachvollziehbarkeit und Umfang der Prüfung

Alle Skripte wurden ausgeführt. `run_all.py` reproduziert die JSON Ausgaben und
nutzt für Klauselablationen den vorhandenen Ursprungssucher im Repository.
Unabhängige Gegenrechnungen prüfen den phasengenauen Code, die logische Wirkung
und die erreichbare Operatorfamilie. Es wurden keine Lean Beweise angelegt oder
behauptet. Die vorhandene Qodo Prüfung war in dieser Laufzeit nicht verfügbar;
das ist keine abgeschlossene externe Qodo Begutachtung.

Der aktuelle Prüfkatalog wurde auf allen zwölf Blättern durchsucht. Atlas war
nicht live erreichbar, weil in dieser Sitzung keine entsprechenden Werkzeuge
vorlagen. AUDIT.md nennt Quellen, Suchbegriffe, Vorläufer und Ausschlussgrenzen.
Die Aussage, dass ein passender Ursprungstest nicht gefunden wurde, ist eine
begrenzte Neuheitsaussage.
