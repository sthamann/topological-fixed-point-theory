# Vier Marken, ursprüngliche Quelle und Messvergleich

13. November 2026. OriginLab prüft den Vorschlag im dauerhaften Suchlauf.
Ein gebundenes Originalpaket kann dieselben Operatorbedingungen prüfen.

## Was die vier Marken auswählen

Die orientierten Ränder der vier Tetraederflächen erzeugen das ganzzahlige
Zyklusgitter des Kantengraphen. Ihre Gramform ist $4I-J$. Drei Ränder
bilden über den Spannbaum eine ganzzahlige Basis mit Determinante 16.
Halbierung ergibt vier Gewichte mit Normquadrat $3/4$. Ihre Differenzen
erzeugen $A_3$. Der Index des Wurzelgitters im Gewichtsgitter ist vier.
Die vier Dreiecksrichtungen und ihre Gegenorientierungen sind die
kürzesten ganzzahligen Zyklen. Ihre Klassen haben Ordnung vier.
Die orientierte Randabbildung liefert das Vorzeichen der
Markenpermutation. Beim Vierzyklus ergibt sich das Spektrum $1,i,-i$.
Die räumliche Realisierung braucht die ursprünglichen Operatoren.

Über $\mathbb F_2$ sind die Graphzyklen und die kleinen Umläufe der
vierfach punktierten Kugel bereits als $S_4$ Darstellungen äquivalent.
Beide sind der Quotient des vierdimensionalen Permutationsraums nach
der Summe der vier Generatoren. Beide haben Fixraumdimension null.
Der Dualraum hat eine feste Richtung. Die geometrische Obstruktion
bleibt bestehen. Die Schnittpaarung der kleinen Umläufe ist null.
Ein isotroper Unterraum kann in einem symplektischen Viererraum
höchstens zwei Dimensionen haben.

Für den vollen Flächenraum gibt es genau eine nichtentartete
$S_4$ invariante alternierende Paarung. Der getrennte Ansatz mit festem
Anker hat unter $S_4$, $A_4$ und $C_4$ keine solche Paarung.
Unter $C_3$ erhält er genau eine. Der Flächenraum ist eine nicht
gespaltene Erweiterung des Graphzyklusraums.

Die volle Paarung hat **zwei** $S_4$ invariante quadratische
Verfeinerungen. Beide sind ungerade. Die geschlossene Oberfläche
hat in beiden $q=0$. Unmarkiert vertauscht eine symplektische
Transformation über der unveränderten Randabbildung die beiden
Verfeinerungen. Bei festgehaltener einzelner Flächenoperation trennt
deren Quadrat $+I$ oder $-I$ die Varianten. OriginLab ermittelt
die Abfrage $q=0$ für $+I$ und $q=1$ für $-I$ und prüft beide
Antworten. Die ursprüngliche Phase der Operation gehört zum Test.
Die vier Generatoren mit
Quadrat $-I$ benötigen daher eine zusätzliche Auswahl. Der Ankershift
von v796 verwendet eine ausgewählte Einzelfläche und liefert die
gerade Partnerform. Die äußere $S_6$ Brücke bewahrt die Unterscheidung
zwischen geometrischer Deckwirkung und Familienclock.

Die Cliffordkontrolle gibt sechzehn Operatorwörter, sechs mit Quadrat
$+I$ und zehn mit Quadrat $-I$. Die volle Markensymmetrie erlaubt
unterschiedliche positive Transferoperatoren. Sie lässt im Viererblock
die skalare Richtung und die Summe der vier hermiteschen
Cliffordgeneratoren zu. Mit Erhaltung der Cliffordparität bleibt nur
ein skalarer Transfer.

[Prüfer](selection.py),
[Integralflüsse und Schnittgitter](https://www.numdam.org/item/BSMF_1997__125_2_167_0/),
[Kurventest v784](../../verification/v784_curve_code_2torsion.py),
[Äußere Brücke v796](../../verification/v796_duad_syntheme_bridge.py).

## Wo Quantenoperation und primitive Ladung zusammenpassen

Eine homogene Operation der primitiven Viererladung kann kein
nichtverschwindendes neutrales Quadrat haben. Ihre doppelte Ladung
wäre zwei modulo vier. Der ausführbare Test prüft Quadratphasen,
alle sechs Antikommutatoren und ihre Wirkung auf dem zyklischen
Raum. Dazu prüft er die drei Cartan Kommutatoren jedes geladenen
Worts auf demselben Raum. Zwei verschiedene Rückkehroperatoren
können dieselbe einzelne äußere Spur haben. Der Prüfer vergleicht
deshalb auch die Operatoren auf dem reflektierten Quotienten.

Ein endlicher Träger mit drei kommutierenden integralen
$A_3$ Cartanladungen und der Wirkung der geraden
Markenpermutationen hat eine scharfe
Untergrenze. Ein nichtnulliges geordnetes Wort aus allen vier
primitiven Gewichten durchläuft vier verschiedene Klassen
modulo Wurzelgitter. Liegt sein besuchtes Gewicht der Klasse null
bei Gewicht null, erzwingen die Permutationen vier Gewichtsorbits
der Größen $1,4,6,4$. Liegt das Gewicht der Klasse null außerhalb
des Ursprungs, haben alle vier Orbits mindestens vier Elemente.
Für jedes Gewicht außerhalb des Ursprungs sind höchstens drei
der zwölf geraden Permutationen im Stabilisator. Die Schranke
gilt daher auch für beide Erweiterungen auf $S_4$.
Damit braucht der Träger mindestens 15 Dimensionen. Das Modell
der echten Teilmengen erreicht die Schranke, wobei ein hinzugefügtes
viertes Element zur leeren Teilmenge zurückführt.

Mit vier komplexen Cliffordoperationen auf demselben Träger muss
seine Dimension ein Vielfaches von vier sein. Vier CAR Moden
realisieren die gemeinsame Algebra auf 16 Zuständen. Die geladenen
Erzeuger $a_i^\dagger$ tragen die vier Gewichte. Aus ihnen und ihren
Adjungierten entstehen $\Gamma_i=a_i^\dagger-a_i$. Mit einer
kohärenten Summe von Vakuum und voller Besetzung liefert das
bezüglich der Cartanladungen neutrale Wort den Wert $1/2$.
Reine Permutation der vier Erzeuger macht das Wort alternierend.
Ein unter dieser Wirkung invarianter Zustand gibt dafür null.

Die orientierte Graphwirkung gibt eine zweite Möglichkeit. Eine
ungerade Permutation kehrt jede Schleife um. Auf demselben
16er Träger kann sie Erzeugung und Vernichtung vertauschen.
Der Prüfer verifiziert alle 24 Matrizen und ihre Gruppenrelationen.
Das Wort $W=a_1^\dagger\cdots a_4^\dagger$ geht unter ungerader
Permutation in $-W^\dagger$ über. Das hermitesche Wort
$i(W-W^\dagger)$ bleibt invariant. Zwei invariante Dichtematrizen
aus den kohärenten Strahlen $|0\rangle\pm i|F\rangle$ geben
dafür exakt $+1$ und $-1$. Damit ermöglicht die orientierte
Symmetrie eine Antwort, legt ihren Wert aber nicht fest.
Das Modell setzt Ladungskonjugation voraus. Eine zusätzliche
Mittelung über die gesamte Teilchenzahl löscht beide Antworten.
P1 muss die ursprüngliche Orientierung, die Zulässigkeit des
Worts und die Zustandsauswahl gemeinsam liefern.

Die vier Higgsrichtungen des gepinnten v252 Quarkblocks liefern
keinen solchen festen Viererträger. Der Kommutator von
$Y_uY_u^\dagger$ und $Y_dY_d^\dagger$ ist invertierbar.
Eine gemeinsame Familiengerade wäre ein Nullvektor dieses
Kommutators. Das reproduziert den vorhandenen Originalausschluss.
v775 bleibt für die direkte Identifikation der endlichen
Codeklassen mit Materie maßgeblich.

[Native Quarkprobe](../../experiments/theory-contracts/clifford-native-quark-source-bridge-20261110/PROOF.txt),
[Code und Materie v775](../../verification/v775_gaussian_class_d5_purity.py).

## Eine vollständig rechenbare Spektralquelle im kleinsten Vergleich

Das Metrologie Original definiert
$e^{-\operatorname{Tr} f_*(D_V^2/\Lambda^2)}d\nu_V$ und nennt
die BRST reduzierte Konfigurationsbasis. Der überprüfte native
Aufruf liefert keinen ausführbaren vollständigen $D_V$,
keine Auswertung des Basismaßes und keine vier physikalischen Marken.
Die Behauptung, ein unabhängiges neues Postulat sei zwingend,
braucht einen Eindeutigkeits oder Unabhängigkeitsnachweis.
Offen ist, welche Prämisse die ursprüngliche gewichtete
Integration bestimmt.

Die ausführbare Variante benutzt $x,y\in\{0,7/20\}$ mit je
Basismasse $1/2$ für zwei Halbintervalle. $D_F(x)$ ist der v252
Diracoperator mit Higgs $H=(0,x)$, allen vier Yukawa Matrizen
und dem festen Majoranablock des Diagnosefixtures. Die Halbintervalle
werden im 192er Diracoperator durch $p\Gamma_F$ mit $p=1/4$ gekoppelt.
Das Profil ist $f(s)=\frac14 e^{-s}$ bei $\Lambda=1$.
Die Fermionantwort ist die gewöhnliche Exterior Spur des
geordneten Transfers mit Halbzeit $2/5$. Für jedes Halbintervall
gilt $X_x=\exp[-s(p\Gamma_F-i\Gamma_FD_F(x))]$.
Der Faktor ist $\det(1+X_xX_y)$. Die Leerreferenz
verwendet die gleichen Regeln mit $x=y=0$. Das Higgswort ist $x^2$.
Die Yukawamarke ist die Ableitung des späten Fermionbeins bei
festgehaltenem Spektralgewicht. Die vier Historien werden vollständig
mit Intervallarithmetik summiert.

Der kleine Vergleich hat eine positive reflektierte Matrix der
Größe zwei und $C_{UY}$ rund $-6.49992038\cdot 10^{-5}$.
Der zweite Vergleich verwendet nur nichtnegative Radien
$0,7/20,7/10,21/20,7/5$, Halbzeit $1/200$ und gleiche Basismassen.
Mit $g=0$ ist der Fermionkern auf diesen fünf Profilen positiv
definit. Bei $g=1/100$ liefert der gemeinsame Kern für den
Koeffizientenvektor $(2408,-9082,17224,-14745,4096)/4096$
die zertifiziert negative reflektierte Norm rund $-1.73324\cdot 10^{-4}$.
Bei $g=1/10$ ist dieselbe Fünfpunktmatrix wieder positiv definit.
Damit ist auch für ein positives Ausgangsmaß und einen punktweise
positiven Faktor ein eigenständiger Nachweis der Reflexionspositivität
nötig. Schon die fünf radialen Testfunktionen widerlegen die
allgemeine Behauptung für diese Spektralvariante.

Bei unbeschränktem Lebesgue Basismaß und beschränktem Profil
gilt für jede endliche Diracdimension
$e^{-\operatorname{Tr} f}\ge e^{-d\sup f}>0$.
Das Gewicht hat damit unendliche Masse. Das gilt hier auch mit
dem gewöhnlichen Fermionfaktor. $X_xX_y$ ist ähnlich zu
$X_x^{1/2}X_yX_x^{1/2}>0$, daher ist $\det(1+X_xX_y)\ge 1$.
Ein endliches Basismaß, eine begründete Beschränkung der Domäne
oder eine koerzive Ergänzung ist für diese Variante erforderlich.
Die Variante verwendet ein deklariertes endliches Basismaß.
Die vier physikalischen R113 Werte bleiben unberechnet.

[Originaler Caller Trace](../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/triangulation_20261106O/caller/ACTUAL-PRIMITIVE-REGULATOR-CALLER-TRACE.txt),
[Originales Basismaß](../../documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex),
[Ausführbarer Spektralvergleich](spectral_candidate.py),
[Spektralwirkung](https://arxiv.org/abs/hep-th/9606001),
[Endliche Dirac Ensembles](https://arxiv.org/abs/1510.01377).

## Alpha und Leptonen

Die fünf verschachtelten Alpha Werte und ihre Abweichungen sind
reproduziert. Die volle Gleichung ergibt $\alpha^{-1}=137.0359992168407$.
Die Abstände betragen rund $1.90\sigma$ zu CODATA, $0.99\sigma$
zu Rubidium, $3.39\sigma$ zur aus dem magnetischen Elektronenmoment
mit SM berechneten Zahl und $6.33\sigma$ zu Cäsium, jeweils
in Einheiten der Unsicherheit des Vergleichswerts. Die Experimente
sind untereinander in Spannung. Der genaue Quillen Anschluss bleibt
laut Ledger offen. Die Ablationen zeigen die Empfindlichkeit der
Gleichung. Die historische Wahl ihrer Koeffizienten folgt daraus nicht.

Die eingefrorenen Leptonwerte sind Quellenmassen. Als direkte
Polmassenwerte liegen sie um $2.198\%$ beziehungsweise $0.215\%$
zu hoch. Die Rechnung bei $M_Z$ aus den veröffentlichten
$\overline{\mathrm{MS}}$ Yukawas vergrößert beide Abstände.
Das widerlegt in diesen beiden Konventionen die unmittelbare
Identifikation. Eine Aussage über jede Skala braucht eine
spezifizierte Abbildung und einen kontrollierten Verlauf.
Die 3 Prozent in v307 sind ein Vergleichsfenster. Die benötigten
relativen Transferfaktoren sind etwa $0.97850$ und $0.99786$.
Ein gemeinsames Nachstellen des einen Seeds reicht nicht.
Die beiden Quellenformeln erzwingen
$(m_\mu/m_\tau)^2/(m_e/m_\mu)=16/21$, gemessen sind rund $0.73106$.
Ein physischer Vergleich braucht eine flavorabhängige
Quellenabbildung. Ihre Herkunft aus R113 bleibt offen.

[Eingefrorenes Register](../../verification/predictions_frozen.json),
[Quillen Ziel](../../verification/v382_alpha_quillen_exact.py),
[CODATA](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf),
[Rubidium](https://doi.org/10.1038/s41586-020-2964-7),
[Elektron](https://doi.org/10.1103/PhysRevLett.130.071801),
[Cäsium](https://doi.org/10.1126/science.aap7706),
[PDG 2026](https://pdg.lbl.gov/2026/listings/particle_properties.html),
[Laufende Leptonwerte](https://doi.org/10.1103/fdcc-ycph).
