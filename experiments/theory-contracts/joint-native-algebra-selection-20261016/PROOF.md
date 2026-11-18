# Gemeinsamer nativer Operatorursprung: präzise Vereinheitlichung

16. Oktober 2026. Forschungs-ID `UR.NATIVE.JOINT_OPERATOR.01`.

## Auftrag und Aussagegrenze

Dieser Text verwendet den vorhandenen TFPT-Kern als mathematische Eingabe. Er ersetzt weder P1/P2 noch die vorhandenen Verifikationen. Gesucht werden tatsächliche gemeinsame Darstellungen, keine passenden Dimensionszahlen. Die Herkunft der Eingaben aus einer physischen Randwirkung wird durch eine algebraische Rekonstruktion nicht automatisch bewiesen.

Die mathematischen Eingaben sind das Hamming-/E8-Ereignissystem mit seiner gewählten gaußschen Struktur J, seine native Viererdarstellung, der symmetrische vierte Tensorgrad und die bereits vorhandene Einregister-Erasure-Konstruktion. Die Wahl von J und der physische Einsatz des symmetrischen vierten Grades bleiben erklärte Voraussetzungen. Unter diesen Voraussetzungen müssen die folgenden drei Träger nicht unabhängig gewählt werden.

## 1. Ein gemeinsames Darstellungsobjekt

Sei V=C^4 die tatsächliche native Darstellung. Die 60 hermiteschen Reflexionen r_j des vorhandenen Hamming/E8/J-Systems erzeugen G=G31. Die vorhandenen Gruppenresultate liefern |G|=46080, eine normale Pauli-Untergruppe P von Ordnung 64, Z(G)=mu4 und G/P=Sp(4,2)=S6. Hier werden die konkreten Matrizen benutzt, nicht lediglich eine abstrakte isomorphe Gruppe. P ist innerhalb dieses G nicht frei auszuwählen: P=O2(G) ist seine größte normale2-Untergruppe und daher charakteristisch. Der Pauli-Unterbau kann aus dem vorhandenen Ereignisverband zurückgewonnen werden.

Der Hamming-Construction-A-Anschluss an das E8-Gitter und die anschließenden60 J-Bahnen sind vorhandene mathematische Ergebnisse. Daraus folgt hier noch keine J-erhaltende Identifikation mit jedem separat markierten D5/A3-Glue- oder Stromwörterbuch; diese muss ihre tatsächliche Abbildung bewahren.

Der gemeinsame Ursprung ist das markierte Darstellungsobjekt (G,V,{r_j},J) zusammen mit Tensorprodukt, Symmetrisierung und dem vorhandenen Erasure-Intertwiner. Die bloße associative Matrixalgebra End(V)=M4 reicht dafür nicht: Sie besitzt keine unitale fünfdimensionale Darstellung. Das Gruppen-/Tensorproduktdatum muss erhalten bleiben.

### 1.1 Der Fünfer ist ein bestimmter Fixraum

Es gilt

    V5 = Fix_P(Sym^4 V).

Der zentrale Operator J=iI erzwingt für positive homogene P-Invarianten einen durch vier teilbaren Grad. Grad vier ist daher der erste mögliche positive Grad. Dass darin genau fünf Richtungen verbleiben, lässt sich ohne numerische Suche zeigen. Schreibe die Koordinaten von V als x_00,x_01,x_10,x_11. Die beiden Z-Paulis erzwingen gleiche Parität aller vier Exponenten. Im Grad vier bleiben somit nur die Muster (4,0,0,0), (2,2,0,0) und (1,1,1,1). Die X-Paulis wirken durch Translation der vier binären Koordinaten. Die resultierenden Bahnsummen sind genau

    sum_a x_a^4,
    x_00^2 x_01^2 + x_10^2 x_11^2,
    x_00^2 x_10^2 + x_01^2 x_11^2,
    x_00^2 x_11^2 + x_01^2 x_10^2,
    x_00 x_01 x_10 x_11.

Nach den gewöhnlichen Tensor-Normierungen sind dies exakt die fünf Spalten des vorhandenen `_code_basis`-Isometrieoperators V: C^5 -> (C^4)^{tensor4}, nicht eine nachträglich gewählte Einbettung. Normalität von P macht diesen Fixraum G-invariant. Die Originalereignisse wirken darauf als die Standarddarstellung von S6; die 60 Reflexionen fallen auf 15 Transpositionen mit je vier Urbildern zusammen:

    r_j^{tensor4} V = V rho5(r_j).

Wichtige Typgrenze: Im VOLLEN vierten Tensorprodukt ist der P-Fixraum 16-dimensional, nicht fünf-dimensional. Die Charakterrechnung ergibt (1/64) sum_p tr(p)^4 = 16. Erst Sym^4 reduziert ihn auf fünf. Die Symmetrisierung darf nicht stillschweigend weggelassen oder schon als physisch ausgewählt bezeichnet werden.

### 1.2 Derselbe Code bestimmt den Erasure-Zwanziger

Der vorhandene Code korrigiert den Verlust eines der vier Register. Der tatsächliche Erasure-Intertwiner W bildet V20=V5 tensor conjugate(V) isometrisch in V^{tensor3} ab und erfüllt

    r_j^{tensor3} W = W [rho5(r_j) tensor conjugate(r_j)].

Somit gilt für JEDES Gruppenwort g, mit derselben Reihenfolge,

    rho20(g) = rho5(g) tensor conjugate(rho4(g)).

Die dekodierte Erasure-Präparation und das logische Readout sind

    E(rho)=rho tensor I4/4,
    R(X)=Tr_4 X.

Native unitäre Ereignisse erfüllen beide Intertwining-Relationen: Ad_rho20(g) E = E Ad_rho5(g) und R Ad_rho20(g)=Ad_rho5(g) R. Die erste gilt für alle logischen Zustände, die zweite für alle gemeinsamen Eingaben X. Gemeinsame klassische Ereignisgewichte erben beide Relationen einschließlich jeder festgehaltenen Ereignisfolge.

Dies vereinheitlicht die native Ereigniswirkung, den quartischen Code und dessen Erasure-Fortsetzung durch tatsächliche Tensor- und Intertwiner-Abbildungen. Die gemeinsame Ursache ist mathematisch konkret; sie ist innerhalb des vorhandenen Kerns keine neue unabhängige physische Bestätigung. Abschnitt3 erweitert diese Operatorverbindung um den tatsächlichen Galois-Träger und seine eindeutige Darstellungsbahn.

## 2. Gemeinsame Algebra und Eindeutigkeit

V4 ist irreduzibel, denn die 16 projektiven Paulis spannen M4. Auf V5 ist P trivial und S6 wirkt standard-irreduzibel. Auf V20 wirkt P als I5 tensor conjugate(P), daher ist sein Kommutant M5 tensor I4. Die zusätzliche S6-Wirkung reduziert ihn auf C I20. Folglich ist auch V20 irreduzibel. Ihre unterschiedlichen Dimensionen machen die drei Darstellungen paarweise nichtisomorph.

Für H=V4 directsum V5 directsum V20 folgt mittels Burnside/Artin-Wedderburn

    A_joint = image(C[G] -> End(H)) = M4 directsum M5 directsum M20,
    A_joint' = C directsum C directsum C.

Die Dimension der Algebra ist 16+25+400=441. Der 29-dimensionale Summenraum ist Buchhaltung der drei Darstellungen, keine hergeleitete Raumzeit- oder Teilchenzahl. A_joint ist der kleinste Quotient der Gruppenalgebra, der genau diese drei Wortdarstellungen erhält: Sein Kern ist der Schnitt ihrer drei Gruppenalgebra-Kerne. Das ist eine präzise universelle Minimalität für diese vorgegebenen Ansichten, keine globale Minimalität über alle denkbaren TFPT-Ursprünge.

Dieser Quotient allein ist außerdem nicht tensorabgeschlossen: Andere Tensorgrade können weitere Irreduzible enthalten. Die volle gemeinsame Ursache bleibt das Gruppen-/Darstellungsobjekt mit seinen Tensoroperationen.

### 2.0 Gemeinsame Zustandsnormierung

Der berechnete Kommutant beantwortet auch eine Zustandsfrage. Jede unter der gesamten nativen Gruppe invariante Dichtematrix auf H hat notwendig die Form

    rho = w4 I4/4 directsum w5 I5/5 directsum w20 I20/20,
    wd>=0, w4+w5+w20=1.

Treue verlangt wd>0, wählt die Gewichte aber nicht aus. Innerhalb jedes einzelnen irreduziblen Blocks ist der gruppeninvariante Zustand eindeutig; die relative Gewichtung der drei Ansichten ist es unter dieser Bedingung nicht. Diese Aussage betrifft Gruppeninvarianz, nicht automatisch alle stationären Zustände eines gegebenen anisotropen Prozesses. Insbesondere liefert ein bloßer Dimensionsquotient5/4 noch keine physische relative Suszeptibilität oder Wirkungsnormierung.

### 2.1 Normalisatoren: zwei unterschiedliche Fragen

Der unitäre Normalisator der ASSOCIATIVEN Algebra ist U4 x U5 x U20. Er ist viel größer als der Normalisator der ursprünglichen Ereignisgruppe. Da P=O2(G) charakteristisch ist, normalisiert jeder U in N_U4(G) auch P. Dessen unitärer Normalisator ist die Zwei-Qubit-Cliffordgruppe bis auf Phasen. Projektiv ist dies genau G, also

    N_U4(G)=U(1) G.

Auf der Dreierdirektsumme können die ungleich großen einfachen Blöcke nicht vertauscht werden. Die treue Viererdarstellung legt den induzierten Gruppenautomorphismus fest; Schurs Lemma bestimmt die anderen Blöcke bis auf Phasen. Daher

    N_U(H)(rho(G)) = U(1)^3 rho(G).

Dies ist eine neue Folgerung aus den vorhandenen Gruppenfakten, kein hier behaupteter bereits formal bewiesener Lean-Satz. Antiunitäre Normalisatoren werden nicht klassifiziert.

### 2.2 Was daraus zu 3+2 folgt

Irreduzibilität von V5 verbietet einen einzelnen nichttrivialen G-invarianten Projektor. Das verbietet KEINE kanonische Familie von Projektoren, die G untereinander vertauscht. Eine solche Familie existiert tatsächlich: Für die 15 vollständigen Paarungen der sechs S6-Punkte sei t=(ab)(cd)(ef). Auf dem sechsdimensionalen Permutationsraum hat t die Eigenwerte +1 dreimal und -1 dreimal. Nach Entfernen des konstanten Vektors bleiben auf V5 die Vielfachheiten +1 zweimal, -1 dreimal. Somit ist

    Q_t=(I5+rho5(t))/2

ein Rang-zwei-Projektor mit einem Rang-drei-Komplement. Dies ist die vorhandene 15er-Familie von nativen 2+3-Antworträumen; eine weitergehende Zuordnung zu den Pauli-Kontexten benötigt die tatsächlichen Originaloperatoren und folgt nicht aus der gemeinsamen Anzahl15. Sie ist ein positiver Strukturzusammenhang, keine Auswahl eines bestimmten schwachen bzw. farbigen Sektors. Die Gruppe permutiert diese15 Möglichkeiten transitiv. Der entscheidende Vergleich mit dem tatsächlich vorhandenen P2=(6Y5+2I)/5 wurde jedoch zusätzlich ausgeführt: KEINER der15Projektoren Q_t ist P2. Das ist exakt und strukturell bewiesen. P2 kommutiert mit Farb-S3_(234). Wäre P2=Q_t, so müsste die zugehörige feste-punktfreie Involution t im Zentralisator S3_(015) dieser Farbgruppe liegen. Dort gibt es keine feste-punktfreie Involution. Die native15er-Antwortfamilie reproduziert die bestehende physikalische Markierung daher nicht. Ihre gleiche Rangaufteilung ist keine zulässige Operatoridentifikation. Ein beliebiger U5-Basiswechsel würde dieses Problem nur durch eine neue, die ursprünglichen Gruppenwirkungen verändernde Einbettung überdecken. Für die bestehende Markierung benötigt man die tatsächliche zusätzliche Quellwirkung, nicht die Auswahl eines der15Antwortprojektoren.

## 3. Clockrelationen richtig typisieren

Die zentrale Drehung hat auf den drei Trägern die Wirkungen

    rho4(J)=iI4, rho5(J)=I5, rho20(J)=-iI20.

Der Compilerclock c=J sigma besitzt auf V4 Ordnung zwölf mit c^4=sigma und c^9=J. Auf V5 ist sein Bild nur rho5(sigma), von Ordnung drei. Er kann daher nicht mit einem quartischen Clock der Ordnung vier identifiziert werden.

Der vorhandene markierte RR-Clock R5 besitzt das charakteristische Polynom (z-1)^2(z+1)(z^2+1). Er ist ein kohärentes Algebraelement und kein einzelnes natives Gruppenwort. Der exakte Test: R5 kommutiert mit der schwachen S2- und der Farb-S3-Wirkung. Deren Zentralisator in S6 hat nur die beiden Elemente1 und die schwache Transposition; er enthält kein Element der Ordnung vier. Die treue Standarddarstellung kann R5 daher nicht als Gruppenmatrix enthalten. Der ergänzende Vergleich mit allen720 Matrizen bestätigt dies. Gleiche Spur mit einer S6-Vierzyklusmatrix hätte keine Identifikation begründet. Eine Wirkung eines solchen kohärenten Wortpolynoms auf V20 muss eigens aus demselben Polynom geprüft werden: Aus R5 allein folgt nicht automatisch R5 tensor I4.

Der Galoisclock C5 lebt zunächst auf Q(zeta5), nicht auf dem nativen Vierer. Seine Relation G C5 G^{-1}=C5^2 findet als F20=C5 semidirect C4 eine Gruppenstruktur im Normalisator eines Fünfzyklus von S6. Auf V5 lautet die Restriktion 1 directsum V_cyc,4. Dieser Anschluss wurde nun an allen sechs nach Fixierung des vorhandenen Punktes verbleibenden C5-Untergruppen geprüft. Für jede wurden tatsächliche native G31-Matrizen C,G mit C^5=G^4=I und GCG^{-1}=C^2 konstruiert. Ein möglicher zentraler Liftfehler wurde dabei explizit ausgeschlossen. Die simultanen Gleichungen C_native A=A C_original und G_native A=A G_original besitzen jeweils einen eindimensionalen Lösungsraum mit Rang(A)=4. Es existiert damit für jeden Kandidaten eine tatsächliche invertierbare Abbildung zum ursprünglichen v419-Galois-Träger, eindeutig bis auf einen Skalar; sie ist keine bloße Dimensionsübereinstimmung. Auch die Metrik ist geprüft: In der ursprünglichen Potenzbasis ist die kanonische zyklotomische Spurmetrik T=5I4−ones(4,4), mit Diagonale4 und Nebendiagonale−1. Für alle sechs Abbildungen gilt exakt A* A=kT mit k=2 bzw.2/5. A/sqrt(k) ist daher eine Isometrie dieser ursprünglichen Spurmetrik in die native euklidische Metrik. Dies legt keine physische Wirkungs- oder Kopplungsnormierung fest. Die sechs Kandidaten sind zunächst in festgehaltenen Koordinaten benannt, nicht sechs nachgewiesene physische Alternativen. Die Gleichwertigkeit wurde auf der logischen Markierungsebene geprüft: S2_(01) x C3_(234) wirkt auf den sechs Untergruppen regulär, S2_(01) x S3_(234) transitiv. Diese Stabilisatoren erhalten die vorhandenen logischen Markierungen. Es handelt sich daher um EINE ausgezeichnete Bahn bzw. Darstellungsklasse, nicht sechs unabhängige physische Wahlmöglichkeiten. Ein wörtlicher Vertreter bleibt koordinatenabhängig. Daraus folgt noch keine Behauptung, dass jede verbindende native Viererliftmatrix mit dem ursprünglichen P1-/Source-Clock kommutiert; dies ist die verbleibende anders getypte Clockfrage.

## 4. Dynamische Auswahl: Originalmenü und Algebraabschluss unterscheiden

### 4.0 Eindeutigkeit innerhalb des ursprünglichen 60-Ereignismenüs

Für einen Poissonprozess ausschließlich aus den tatsächlichen 60 Reflexionen sei

    L_d = sum_{j=1}^{60} lambda_j (Ad_rho_d(r_j)-id), lambda_j>=0.

Die gemeinsame lineare Abbildung der 60 Raten auf die VOLLSTÄNDIGEN Vierer- und Fünfergeneratoren besitzt Rang60. Einzeln betragen die Ränge55 auf V4 und15 auf V5; die gemeinsame Nullität ist null. Die Rechnung verwendet ganzzahlige native Pauli-Adjunktmatrizen und direkt aus denselben Strahlen rekonstruierte Quartikmatrizen: Für die unnormierte Codewortbasis B wird (2r_j)^{tensor4}B=B A_j auf allen256Zeilen exakt geprüft, mit A_j^2=256I und rho5(r_j)=A_j/16. Ein numerisches S6-Matching wird für diesen Rangnachweis nicht benutzt. Ein60x60Minor der881x60-Matrix hat modulo1009 die Determinante72 und zertifiziert den vollen Rang. Folglich bestimmen L4 UND L5 sämtliche lambda_j eindeutig und damit auch L20. Es gibt in diesem ursprünglichen primitiven Menü keinen Zwanziger-Ratentwin bei unveränderten vollständigen Vierer-/Fünfergeneratoren.

Diese Aussage ist stärker als ein Vergleich weniger Populationen oder Transferwerte, aber enger als eine Klassifikation aller Quellen: Die vollständigen Generatoren und das primitive Menü sind Voraussetzungen, keine aus P1 hergeleiteten Daten. Ein einzelner Kanal bei nur einer Zeit kann seinen Generator nicht ohne Weiteres eindeutig bestimmen. Eine zusätzliche Nulloperation/skalare Phase wäre zudem im Generator unsichtbar und gehört nicht zu den hier gezählten60 Reflexionen.


Der ältere S3-Quellenprozess besitzt eine zusätzliche Vorzeichendarstellung im Zwanziger, die Vierer und Fünfer nicht sehen. Sein kohärenter Sprung z=1-2e_sign verändert die Erasure-Fortsetzung. Der neue exakte Erasure-Test liefert sogar einen expliziten Eingang mit Spurentfernung1/4 nach diesem zusätzlichen Sprung. Die natürliche logische Fortsetzung R5 tensor I4 kommutiert mit ihm; sie ist ohne einen ursprünglichen gemeinsamen Clocklift kein Beleg für die Wirkung des tatsächlichen Source-Clocks. Der alte Befund ist daher KEIN Gegenbeispiel gegen die STÄRKERE Forderung, E samt allen logischen Zuständen und Records unverändert zu erhalten. Ebenso darf volle G31-Kovarianz nicht nur dem Twin auferlegt werden, wenn der anisotrope Ausgangsprozess selbst schon nicht vollständig farb-S3-kovariant ist.

Die strengere verbleibende Frage lässt sich an einer klar definierten Klasse vollständig beantworten. Betrachtet wird eine unsichtbare Restoperation Phi auf V20, die (i) vollständig positiv und spurerhaltend ist, (ii) das VOLLSTÄNDIGE logische Readout für jede gemeinsame Eingabe erhält, Tr_4 Phi(X)=Tr_4 X, und (iii) unter der nativen G31-Produktdarstellung kovariant ist. Dies ist stärker als bloß Phi(E rho)=E rho. Es ist keine Klassifikation sämtlicher physischer Quellen.

### 4.1 Faktorwirkung wird hergeleitet

Im Heisenbergbild fixiert T=Phi* punktweise M5 tensor I4. Jede dortige unitäre Matrix U erfüllt T(U)=U und T(U*U)=T(U)*T(U)=I, ebenso für UU*. Daher gehört sie zum multiplikativen Bereich von T. Es folgt für alle A,B

    T((A tensor I)Y)=(A tensor I)T(Y),
    T(Y(A tensor I))=T(Y)(A tensor I).

Insbesondere kommutiert T(I tensor B) mit M5 tensor I, ist also I tensor D*(B). Damit ist notwendig

    Phi = id5 tensor D

mit einem CPTP-Gaugekanal D. Die Faktorwirkung wurde nicht als passende Kopplung eingesetzt, sondern folgt hier aus der ausdrücklich starken vollständigen Nichtstörung des logischen Readouts.

### 4.2 Native Symmetrie klassifiziert D vollständig

Pauli-Konjugationen zerlegen die spurlosen Vierermatrizen in 15 verschiedene eindimensionale Charakterräume. Der Cliffordquotient permutiert diese 15 Pauli-Richtungen transitiv. G31-Kovarianz erzwingt daher auf ihnen einen einzigen reellen Eigenwert q und auf I4 die Identität. Genau die Kanäle

    D_q(X)=qX+(1-q) tr(X) I4/4

bleiben. Die Eigenwerte der UNNORMIERTEN Choi-Matrix sind (1+15q)/4 einmal und (1-q)/4 fünfzehnmal (bei normierter Choi-Matrix jeweils Nenner16), deshalb

    -1/15 <= q <= 1.

Dies ist die vollständige diskrete Klassifikation unter (i)-(iii). Sie betrifft die unsichtbare Restoperation im erweiterten Algebra-Menü, nicht neue Ratenfreiheiten innerhalb des zuvor vollständig bestimmten60-Reflexionsprozesses. Für eine stetige Halbgruppe ab der Identität gilt q(t+s)=q(t)q(s), q(0)=1. Stetigkeit verhindert eine Nullstelle oder einen Vorzeichenwechsel bei endlicher Zeit; die CPTP-Schranken geben exakt

    q(t)=exp(-gamma t), gamma>=0.

Native Symmetrie und die vollständige logische Nichtstörung bestimmen somit die FORM dieser Restdynamik, aber nicht gamma oder die physische Zeiteinheit.

### 4.3 Gemeinsame normierte Algebra-Realisierung ohne frei eingesetzten Matrixprojektor

Aus den zentralen J-Charakteren folgt das ursprüngliche Wortpolynom

    e20=(1/4) sum_{k=0}^3 i^k J^k,

das auf V4,V5,V20 genau (0,0,I20) wirkt. Für die 16 projektiven Paulis p_m, p_0=1, setze im dargestellten Quellenalgebraabschluss

    u_m=(1-e20)+e20 p_m,
    K_m=sqrt(w_m) u_m,
    w_0=(1+15q)/16, w_{m>0}=(1-q)/16.

Da e20 zentral ist, ist jedes u_m unitär. Es gilt K_m* K_m=w_m I in ALLEN drei Darstellungen und sum_m K_m* K_m=I. Auf Vierer und Fünfer ist der grobe Kanal die Identität; auf dem Zwanziger genau id5 tensor D_q. Die gemeinsamen Recordgewichte sind sektor- und zustandsunabhängig. Für jeden einzelnen Record m gilt

    K_m,20 E(rho) K_m,20* = E(K_m,5 rho K_m,5*) = w_m E(rho).

Die gleiche Gleichung bleibt für jede gemischte Folge solcher Records und nativer Ereignisse wahr. Der Kanal ist G31-kovariant; das markierte Instrument ist äquivariant, wenn die Pauli-Recordklassen mittransformiert werden. Die einzelnen festgehaltenen Labels sind nicht alle invariant. Diese Äquivarianz mit16 projektiven Labels gilt sektorweise bzw. auf A_joint, nicht auf der größeren Algebra End(H29) mit zusätzlichen Kohärenzen zwischen den drei Blöcken: Pauli-Konjugationsphasen werden nach dem Sektorgating zu zentralen Blockphasen und sind auf solchen Zusatzkohärenzen sichtbar. Eine Erhöhung auf64 phasenmarkierte Labels wäre eine andere Recordverfeinerung und wird hier nicht stillschweigend vorgenommen.

Der Restkanal kommutiert mit jedem nativen Produkt-Adjunktkanal und damit mit jedem aus diesen Ereignissen gebildeten Random-Unitary-Grundprozess. Er kann deshalb mit dem vorhandenen Prozess komponiert werden: Vierer-/Fünferdynamik und die Erasure-Relation bleiben erhalten, obwohl die Gauge-Korrelation sich ändert. Dafür muss der Grundprozess selbst nicht vollständig G31-kovariant sein. Die zusätzlichen Records sind jedoch eine Erweiterung des ursprünglichen Ereignisalphabetes. Vergisst man sie, bleiben die alten marginalen Vierer-/Fünferdaten gleich; ein vollständig etikettiertes ursprüngliches Recordgesetz wird dadurch nicht wörtlich zum selben Gesetz.

Die kontinuierliche Realisierung lautet

    L=(gamma/16) sum_{m=1}^{15}(Ad_u_m-id).

Sie ist auf Vierer und Fünfer null und auf dem Zwanziger id5 tensor gamma(D_0-id). Für ein spurloses Pauli Q gilt die zweizeitige tracial normierte Korrelation

    C_Q(t)=(1/20) tr[(I5 tensor Q) Phi_t(I5 tensor Q)]=exp(-gamma t).

Damit sind die verbleibenden Prozesse innerhalb dieser Klasse mathematisch unterscheidbar. Ein physischer Zugang zu Q und die Identifikation dieses Tests mit einer realen Feldantwort sind hier NICHT hergeleitet. I20/20 ist stationär; für gamma>0 ist die Fixalgebra M5 tensor I4, daher wird kein eindeutiger globaler Zustand ausgewählt.

WICHTIGE MENÜGRENZE: u_m sind kohärente Elemente der vorhandenen Gruppenalgebra, im Allgemeinen keine einzelnen ursprünglichen Gruppenwörter oder Reflexionen. Die Verwendung von e20 als Dynamik-Gate ist trotz seiner algebraischen Kanonizität eine zusätzliche Auswahl. Ein ursprünglicher Vertrag, der nur klassische Gewichte über native Reflexionen erlaubt, muss separat geprüft werden. Diese Konstruktion behauptet keine Erlaubnis solcher Eingriffe aus P1. Ebenso folgt aus G31-Kovarianz keine Kovarianz unter beliebigen kohärenten Wortpolynomen oder unter einem noch nicht identifizierten ursprünglichen P1-Clock.

Negativkontrolle: e20 allein oder z20=1-2e20 wirkt auf V20 skalar. Ad_z20 ist dort die Identität und dämpft keine innere Korrelation. Erst die nichtskalare Pauli-Wirkung liefert die oben berechnete Rate.

## 5. Physikalische Konsequenzen und Evidenzklassen

**A — notwendige mathematische Konsequenzen der erklärten Kerndaten:** der konkrete symmetrische Quartikfixraum, die Ereignis-/Erasure-Intertwiner, die gemeinsame Blockalgebra, ihre Kommutanten/Normalisatoren und die oben klar eingegrenzte Kanalklassifikation. Sie teilen Eingaben und zählen nicht als unabhängige physische Treffer. Das Hamming/E8/J-Konstrukt, der Vierer und seine Quartik sind keine drei unabhängigen Experimente.

**B — unabhängige prospektive physische Vorhersagen:** Diese Runde liefert keine neue Messbestätigung und noch keine parameterfreie physische Zahl. Die Korrelation C_Q ist ein mathematischer Diskriminator mit freiem gamma. Der bestehende eingefrorene Vorhersagekatalog muss nach tatsächlichem Festlegungsdatum, späterer Datenverfügbarkeit und gemeinsam verwendeten Seeds bewertet werden. Ein Katalogeintrag allein beweist weder Prospektivität noch Unabhängigkeit.

**C — zusätzliche physische Identifikationen:** V20 hier ist der dekodierte Erasure-Träger, nicht ohne Brückensatz der geladene E8-Stromsektor. Die äußere gerade Algebra Lambda_even C5 hat Dimension16; dies identifiziert sie allein weder mit einem physischen Spinor noch mit dem vollständigen Pauli-Fixraum im vierten Tensorgrad. Der physische 3+2-Schnitt, lokale Fermionen, Higgswirkung und zulässige Randinterventionen sind durch die obigen endlichen Darstellungsresultate nicht ausgewählt.

Die Original-Alphaquelle setzt unter anderem B_gamma=5/4, eine konkrete Abhängigkeit q(alpha), eine Logarithmusdeterminante und einen lokalen Closureterm ein. Der vorhandene positive-Root-Satz beweist die eindeutige Lösung der so gegebenen Funktion, nicht die physische Auswahl ihrer Koeffizienten. Die Dimensionen fünf und vier in einer gemeinsamen Darstellung begründen keine Suszeptibilitätsnormierung 5/4. Dafür müsste dieselbe tatsächliche Quellantwort samt Maß, Skala und Feldabbildung diese Koeffizienten liefern.

Bei Flavor ist die bereits abgeleitete Ratio55/117 innerhalb ihres angegebenen Stratum kein neu offenes Problem. Die absolute Amplitude/der physische Uf*-Anker benötigt weiterhin einen aus derselben Quelle bestimmten Readout. Die gemeinsame finite Darstellung ersetzt diese fehlende Normierung nicht.

## 6. Was an unabhängigen Entscheidungen entfällt

Bei festgehaltenem nativen Vierer, Pauli-Wirkung und symmetrischem Grad vier entfallen eine freie Wahl der fünf Codebasisrichtungen, eine freie Wahl ihrer Ereigniswirkung und eine freie Wahl der Erasure-Sektorwirkung. Die vorhandene physikalische3+2-Markierung folgt dagegen nicht aus der nativen15er-Antwortfamilie; deren erste direkte Identifikation ist exakt ausgeschlossen. Basiswechsel ändern nur Koordinaten. Die gemeinsame Algebra muss nicht angepasst werden, sondern folgt aus diesen Wirkungen. Innerhalb der starken Readout-/Kovarianzklasse entfallen eine unabhängig angenommene Faktorwirkung und 15 getrennte Pauli-Dämpfungsfaktoren; genau ein q bzw. eine Rate gamma bleibt.

Nicht entfallen die gaußsche Ausgangswahl J, der physische Grund für den symmetrischen vierten Grad, die Auswahl eines physischen 3+2-Schnitts, das primitive Ereignismenü, der Zustands-/Zeitskalenanker, die weiterhin nicht identifizierten ursprünglichen P1-/markierten kohärenten Clocklifts und die physische Feld-/Kopplungsnormierung. Die kleinste nächste Herkunftsprüfung betrifft diese tatsächlichen Quellenoperatoren; sie ist kein Anlass, den vorhandenen gemeinsamen Kern wieder zu verwerfen.

## Reproduktion und Quellen

Die zusätzliche exakte Prüfung der Galois-Spurmetrik steht in `galois_metric_check.py` und `galois_metric_results.json`.

Originale: `tfpt_explorer/process.py` (insbesondere `_code_basis`, `_invariants`), `tfpt_explorer/quartic_selection.py`, `verification/v634_st31_structure.py`, `verification/v419_seam_galois_carrier.py`, der bestehende Vertrag `compiler-quartic-carry-20260919` sowie die unabhängige Quellenprüfung vom8.Oktober. Die genaue Dateifassung steht in `source_manifest.json`. `check_hidden_gauge.py` und `hidden_gauge_results.json` unterscheiden exakte symbolische Rechnungen, exakte ganzzahlige Ränge und ergänzende Matrixresiduen.

Die Clifford-/Code-Verknüpfung ist als mathematisches Umfeld bekannt; vgl. [Nebe/Rains/Sloane](https://arxiv.org/abs/math/0001038). Die verwendete Normalisator-Konvention entspricht [Tolar](https://arxiv.org/abs/1810.10259). Diese Primärquellen ersetzen nicht die Prüfung der konkreten TFPT-Operatoren. Für diese Runde wird keine historische Priorität der allgemeinen Cliffordstruktur beansprucht.
