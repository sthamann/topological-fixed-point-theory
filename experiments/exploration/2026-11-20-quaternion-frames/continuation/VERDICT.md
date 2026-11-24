# Quaternionenquelle, nativer Träger und kohärentes Gedächtnis

Fortsetzung vom 20. November 2026. Ausgangspunkt ist die Quaternionenkonstruktion in `../VERDICT.md`, Commit `b209e747198332afcc207f917fa41f1ea9202597`. Die Rechnung produziert unter ausdrücklichen Modellannahmen den aufgezeichneten sichtbaren P2 Träger samt Projektor und reproduziert den gesamten nativen 68 dimensionalen Einzelwurzelbaustein bei gegebener Vertexregel. Alle 204 Operatorbilder stimmen exakt mit der nativen Fockrealisierung überein. Im selben endlichen Modell entsteht eine geometrische Vierteldrehung als Holonomie. Der native P1 Randkern, die physische Ausführung der Interventionen und der volle Elternprozess sind weiterhin offen. Keine neue Herleitung von Alpha, Flavor, wechselwirkendem 3+1 Kontinuum oder Gravitation.

## 1. Twist und tatsächlicher Fünferträger

Der komplexifizierte Quadrupolraum ist die Spin 2 Darstellung. Restriktion einer spurfreien Quadrik auf den projektiven Nullkegel identifiziert ihn mit $H^0(\mathbb P^1,K^{-2})$. In der durch das harmonische Quadrat bestimmten Karte $D=\{z^4=-1\}$ ist $H^0(K^{-2}(-D))$ eine intrinsische Linie. Ihr Schnitt $(1+z^4)(\partial_z)^2$ transformiert unter $r:z\mapsto iz$ mit Minus und unter $s:z\mapsto1/z$ mit Plus. Division durch den Schnitt liefert den gewöhnlichen Raum $H^0(\mathcal O(D))$ tensoriell mit dieser Linie. Die Linie trägt das Determinantenzeichen der drei Familienrichtungen. Der in der ersten Sitzung gefundene Twist entsteht damit aus Geometrie und Familienantwort. Das Programm prüft die konkrete invertierbare Nullkegelabbildung, die Generatoren und alle vier Zeichenalternativen.

Ein zweiter, eigenständiger Schritt erreicht die tatsächlichen nativen Quartiken. Die vier projizierten Marken sind ein Torsor der Kleingruppe $V_4$. Auf den Marken wirken Translationen und Multiplikation mit den reellen Charakteren. Fordert man für jede Verschiebung einen hermiteschen involutiven Repräsentanten als scharfe binäre Probe, braucht ein antikommutierendes Produkt eine Phase $\pm i$. Der erzeugte endliche Heisenberg Operatorverband enthält daher das Zentrum $i$. Invarianz setzt Grad durch vier teilbar voraus. Bei Grad vier verlangt Charakterinvarianz, dass die XOR Ladung der belegten Marken null ist. Die Translationsorbits sind

- die Summe aller $x_i^4$,
- die drei Paarungssummen $x_i^2x_j^2+x_k^2x_l^2$,
- das Produkt $x_0x_1x_2x_3$.

Dies sind exakt die fünf nativen Quartiken mit roher Fockmetrik $(96,8,8,8,1)$. Das Zentrum $i$ wirkt auf den vier positiven Einteilchenrichtungen und fixiert die Quartiken. Das stimmt mit dem nativen sichtbaren Cartan $J$ überein. Die geometrische Vierteldrehung auf dem Quadrupol wirkt mit Multiplizitäten $(1,1,2,1)$ und hat eine andere Aufgabe. Ein Skalarzeichen kann sie nicht mit dem sichtbaren nativen $J=I$ identifizieren. Der Quartikraum und der geometrische Quadrupolraum haben auch bei natürlicher geometrischer Wirkung unterschiedliche Charaktere. [Native Quartiken](https://github.com/sthamann/topological-fixed-point-theory/blob/a1273187b82488d34d2a474470152a41d42d4055/experiments/theory-contracts/joint-native-algebra-selection-20261016/check_mixed_return.py#L50)

### Der native P2 Projektor aus zwei Beobachtungen

Die 15 Ereignisse des Quaternioneninstruments sind die Kanten eines vollständigen Graphen auf sechs Punkten. Verschiedene erlaubte aufeinanderfolgende Kanten $va,vb$ treffen sich in $v$. Die zugehörigen Transpositionen liefern als geordnetes Produkt den Dreierzyklus auf $\{v,a,b\}$. Die fehlende Kante $ab$ gehört zu genau einer der fünf bereits abgeleiteten Paarungen. Diese Paarung verbindet $v$ mit genau einem Punkt $q$ außerhalb des bewegten Tripels. Damit sind der native Typ des Dreierzyklus, der ausgezeichnete feste Punkt $q$ und die beiden übrigen festen Plätze durch dieselbe Ereignisfolge gegeben. Die exakte Zählung liefert alle 60 Kombinationen aus bewegtem Tripel und Anker zweimal, mit umgekehrter Zyklusrichtung. Als unabhängige, nur auf der Markovmatrix beruhende Kontrolle erzeugt der erste Austritt aus dem Tripel ebenfalls alle Kombinationen. Die Austrittswahrscheinlichkeit ist an jedem Innenschritt $1/2$, jeder neue Punkt hat Wahrscheinlichkeit $1/6$ pro Schritt. Das zeigt die Bedingtheit der Regelwahl und liefert keine absolut ausgezeichnete Geschichte.

Schreibt man $n_i$ für die sechs zentrierten Simplexvektoren, haben die fünf Differenzen $D_q=(n_j-n_q)_{j\ne q}$ Gram $I+\mathbf1\mathbf1^T$. Das positive Sternereignis $F_q=\sum_{e\ni q}E_e$ hat Spektrum $(1,1/6,1/6,1/6,1/6)$. Die gemeinsame positive Konditionierung $F_q^{-1/2}E_eF_q^{-1/2}$ erzeugt die orthonormalen Polarrichtungen. Sie ist implementierbar als positive Vorfilterung $F_q^{1/2}$ und anschließende projektive Messung. Im aufgezeichneten nativen Simplexrahmen ergibt die Markierung $q=5$ und bewegtes Tripel $\{2,3,4\}$ genau den vorhandenen Projektor $P_2=F_\sigma-c_qc_q^\dagger$. Die bloße Teilung der Quartikbasis hat dagegen quadrierten Projektorabstand $(21-6\sqrt6)/25$. Die alte positive Polarformel und der native Simplex sind Vergleichsdaten. Neu ist die Auswahl von Zyklus und Anker aus der Quaternionenpaarung. Aus dem Rang eins Messeffekt folgt die logische Reflexion algebraisch über $I-6E_e$. Ihre kohärente Ausführung ist eine andere Intervention als das Registrieren eines Ergebnisses. Ein gemeinsamer physischer Operatorlift zwischen beiden ist weiterhin erforderlich. [Nativer Projektor und Clock](https://github.com/sthamann/topological-fixed-point-theory/blob/a1273187b82488d34d2a474470152a41d42d4055/experiments/theory-contracts/joint-native-algebra-selection-20261016/MARKED-CLOCK-LIFT.md#L51)

## 2. Ein natives Gedächtnissignal aus der Geometrie

Für die normierte native Richtung $f_4=x_0x_1x_2x_3$ liefern die Originalmatrizen exakt

$$
\langle L\rangle=32,\qquad
\left(\langle LJ^kL\rangle\right)_{k=0}^3=(6144,0,-2048,0).
$$

Unabhängig davon hat die axiale geometrische Vierteldrehung auf der Zweiebene Eigenwerte $(1,-1)$ und auf der Dreierseite $(-1,i,-i)$. Ihr Übergangsraum $\operatorname{Hom}(W_2,W_3)$ hat vier Charaktergewichte $(1,2,1,2)$ und exakte Spurfolge $(6,0,-2,0)$. Durch $32^2$ dividiert ergeben die Originalantworten dieselbe Folge. Mit normierter Spur sind beide Folgen $(1,0,-1/3,0)$. Unter den drei ungerichteten Typen eigener Block und Übergangsblock trifft genau der Übergangsblock. Die Folge ist ein positiver zyklischer Wortkern in der Phase, da ihre vier Fouriergewichte positiv sind. Die Spurwahl auf dem Übergangsraum und ihre physische Zuordnung zum ausgezeichneten $f_4$ bleiben Annahmen. Vom einzelnen übereinstimmenden Zyklus folgt keine Abbildung des ganzen Antwortsystems.

Das Minus bei $k=2$ schließt **jede** direkte Semantik dieser nativen Worte als Folge positiver Zustandsabbildungen mit positivem Effekt aus. Das gilt auch für größere Systeme. Ein kohärenter Vergleich von $Lf_4$ und $J^kLf_4$ über zwei Pfade liefert gültige Plusportwahrscheinlichkeiten $(1,1/2,1/3,1/2)$. Der negative Halbdrehungswert bezeichnet relatives Interferenzgewicht. Das ist die erste genaue Messvorschrift für einen aus beiden Modellen vergleichbaren Teilkern. [Originalmatrizen R023](https://github.com/sthamann/topological-fixed-point-theory/blob/a1273187b82488d34d2a474470152a41d42d4055/docs/research/fundamental-plan/evidence/RECURSIVE-WORD-SOURCE-PROOF.txt#L140)

### Nativer Cartanoperator und Quellkoeffizienten aus der Quaternionensaat

Ein geordneter nichttrivialer Flagwechsel orientiert über sein quartisches Winkelmoment die erste Achse $i$. Das Zeichen der nächsten Achse wird so gewählt, dass ihr Skalarprodukt mit $i$ positiv ist. Ihr größerer absoluter Überlapp mit einer der beiden übrigen Achsen wählt $j$, mit positivem Überlapp als Zeichenkonvention. $k=ij$ schließt den Quaternionenrahmen. Die Nichtentartung ist für alle 120 Wechsel exakt geprüft. Diese gerichtete Auswahl und die folgende Auswahl der Wurzel als $\beta=-j$ sind deklarierte Prozessregeln.

Der Rahmen hebt sich zu $Q_8$ in der binären Ikosaedergruppe. Linksmultiplikation mit $i$ ist Ordnung vier und liefert den echten Cartanoperator $J$, statt die Quadratdrehung mit ihm zu identifizieren. Auf dem Ikosianergitter wird durch vier orthogonale $J$ Paare die Originalkarte auf allen 240 Wurzeln exakt erhalten. Der zu $1+i$ senkrechte Schnitt im rationalen Quaternionenraum hat 12 $A_3$ Wurzeln. Die orthogonalen Wurzeln im vollen Gitter ergeben 40 $D_5$ Wurzeln. Eine gewählte geladene Chiralität liefert vier Fasern zu je 16 Wurzeln. Maximales Skalarprodukt von $\beta$ mit ihren tetraedrischen Gewichten wählt genau eine Faser aus. Das spurfreie zweite Moment der 48 übrigen Wurzeln ergibt $R=4P_{A_3}-12P_w$, wobei $P_w$ auf die normierte Richtung der entfernten Marke projiziert. Die tangentiale Antwort ist

$$v=-R\beta+\frac{\langle\beta,R\beta\rangle}{2}\beta=4(1-i).$$

Die beiden aktiven gaußschen Koeffizienten sind damit $b=-i$ und $c=4-4i$ auf verschiedenen Paaren. Die Stärke ist eine Folge des normierten Wurzelmoments in dieser Klasse. Sie wurde nicht aus den Zielamplituden eingesetzt. Die Verwendung dieser Schnitt und Momentregel, die Festlegung der geladenen Chiralität und die Übertragung der nativen Vertexoperation gehören weiter zu den physisch offenen Prämissen.

### Vollständiger Vergleich des nativen Wortkerns

Ein unabhängiges Programm übernimmt die beiden geometrisch berechneten Koeffizienten. Es extrahiert die Modenkoeffizienten gemäß der nativen Vertexregel und rechnet die Wortabschlüsse exakt neu, ohne die gespeicherte 68er Basis zu verwenden. Die zwei restlichen Paare sind Zuschauer. Daraus entstehen ein 50er und drei 6er Beiträge mit Gesamtcharakter $(23,14,17,14)$. Nullsetzen der Wurzelverschiebung liefert stattdessen $(10,4,4,2)$. Der kleine $f_4$ Bereich ergibt $32,6144,0$ bereits alleine. In einem getrennten Vergleich werden anschließend alle 68 gespeicherten Originalpolynome unter $L,J,J^\dagger$ geprüft. Alle 204 Bilder sind exakt identisch. Damit stimmen in dieser gemeinsamen Fockkarte sämtliche Wörter des gewählten Einzelwurzelbausteins überein. Das Ergebnis umfasst weder den ganzen geladenen Feldintegralprozess noch eine Herleitung seiner Vertexregel aus der Flagwahrscheinlichkeit.

Der Charakter unterscheidet weitere naheliegende Modelle. Eine achtstufige Gedächtnisdichte mit vollem spurfreien 63er Raum und reversibler vollständig positiver Uhr kann die native Halbdrehungsspur nicht tragen. Für die Restriktion einer gewöhnlichen endlichen $SO(3)$ Darstellung auf eine Vierteldrehung gilt $d+2\operatorname{tr}r-3\operatorname{tr}r^2=0\pmod8$. Das folgt für jeden ganzzahligen Spin durch Periodizität nach Spin vier. Der native Ausdruck ist $44=4\pmod8$. Eine rein geometrische Vierteldrehung auf beliebig vielen gewöhnlichen Rotationstensoren reicht also nicht. Die beiden Tests lassen zusätzliche spinorielle oder nicht voll rotationsäquivariante Sektoren zu.

## 3. Rand, Orientierung und Uhr

Auf jedem der 120 erlaubten gerichteten Übergänge liefert maximale Überlappung der konkreten vier Quadratpunkte eine eindeutige $D_4$ Zuordnung. In den 32 geschlossenen Dreieckswegen pro Ausgangsereignis treten genau zwölf Identitäten, 16 Spiegelungen und je zwei Vierteldrehungen jeder Richtung auf. Die geometrische Vierteldrehung ist somit eine Holonomie derselben Geschichten. Das Orientierungsvorzeichen wird erst durch gerichteten Transport lesbar.

Die Transfers auf 60 Zuständen für Ecken und Seiten sind strikt positiv. Der orientierte 120er Lift hat genau drei negative Eigenwerte $(4-2\sqrt5)/12$. Dies ist ein Ausschluss des **gleichzeitigen** exakten stochastischen Lifts derselben markierten Einzelschritte, globaler Vierteldrehungskovarianz und Positivität über den Zeitschnitt. Denn wenn der Punkttransport $k\mapsto\pi(k)$ feststeht und die globale Vierteldrehung $(k,s)\mapsto(k+s,s)$ lautet, erzwingt Kommutation auf jedem erlaubten Ziel $\pi(k+s)=\pi(k)+s'$. Damit ist der deterministische Orientierungstransport der einzige kovariante stochastische Lift. Er hat die negative Signatur. Zwei konstruktive Auswege haben unterschiedliche Antworten.

- Zwei Schritte zu blocken macht den Transfer positiv und ändert die Einzelschrittprobe.
- Die drei negativen trivialen Uhrmoden abzutrennen gibt einen 117 dimensionalen positiven Hilbertraumtransfer mit derselben globalen Uhr. Die 60 gerade orientierten Antworten bleiben dabei erhalten. Die vorzeichenempfindlichen Antworten ändern sich und das reduzierte Modell ist kein stochastischer Lift aller 120 Ergebnisse. Sein Transfer kommutiert mit der Uhr. Bei einem davon getrennten nativen $L$ muss gerade die Vermischung der Uhrcharaktere hergeleitet werden.

Zufällige Orientierungswechsel bei Bewegungen geben eine dritte Kontrolle. Die kleinste Flippwahrscheinlichkeit zur Positivität ist $(5-2\sqrt5)/10$. Sie lässt das unorientierte Gesetz fest und bricht die globale Vierteldrehungskovarianz. Die so bestimmten Transfers sind Modellrekonstruktionen. Sie liefern den nativen Dirac Randoperator und seine Normierung $c_3$ nicht.

## 4. Warum selbst alle Flaggeschichten keine physische Zeit festlegen

Die ursprüngliche skalare Flagmatrix hat Eigenwerte $1,1/2,1/6$ mit Multiplizitäten $1,5,9$. Sie liefert ein diskretes reflektionspositives Modell mit gleichverteiltem Zustand und Energieabständen $\log2$ und $\log6$. Die Abstände sind inkommensurabel. Ihre reelle Zeitentwicklung kann bei keinem gemeinsamen Maßstab eine projektiv endliche Vierteldrehung sein. Dieselbe Matrix kann auch keine positive Zeit eines kontinuierlichen Markovprozesses auf genau diesen 15 Zuständen sein. Disjunkte Flaggen haben einmal Übergang null und zweimal Übergang $1/36$, was der Positivität einer eingebetteten zeitlich homogenen Markovhalbgruppe widerspricht.

Zwei unendliche Pfadalgebren reproduzieren **jedes** unmarkierte Flagwort. Belässt man einen gewichteten Selbstpfeil mit Gewicht $1/3$ und Nichtselbstpfeile $1/12$, erzeugen Rücklaufgewichte $\log3$ und $2\log12$ eine dichte additive Gruppe. Spaltet man jeden Selbstpfeil in vier gleichwahrscheinliche Markenpfeile auf, haben alle zwölf Mikropfeile Gewicht $1/12$. Dann ist die Gruppe $(\log12)\mathbb Z$. Die konservative Graph KMS Konstruktion liefert entsprechend Faktortyp $\mathrm{III}_1$ beziehungsweise $\mathrm{III}_{1/12}$. Die Gleichheit aller groben Worte folgt durch Summieren über $4^k$ Mikroweisen bei $k$ Wiederholungen, und die Primzahlbewertungen machen die dichte Gruppe exakt. Ein detailliertes physisches Markenexperiment oder die native kohärente Quellenantwort muss die Eltern unterscheiden. [Thomsen, Satz 3.3](https://arxiv.org/abs/1805.05637)

## Entscheidung

**Kandidat mit bedingter P2 Konstruktion und exakt reproduziertem R023 Einzelwurzelbaustein.** Die neue Auswahl stimmt mit dem tatsächlichen P2 Projektor unter der nativen Simplexkarte überein. Die Geometrie erzwingt das passende Charakterzeichen. Unter der gewählten Schnitt und Momentregel trifft der Quaternionenweg auch die nativen Koeffizienten und den gesamten geordneten R023 Wortkern. Gemeinsame Nenner wie $6,16,68$ wurden nur dann als Treffer gewertet, wenn eine Operatorantwort oder eine Darstellung zusätzlich übereinstimmt.

**Die vollständige Lösung bleibt offen.** Der nächste entscheidende Test ist die physische Herleitung der hier noch benutzten Vertexoperation mit ihrer Fockadjungierten aus der gemeinsamen markierten Geschichte. Die R023 Prüfung ist dafür ein fertiger Referenztest. Die Auswahl einer einzelnen physisch ausgezeichneten Geschichte und der Übergang zur vollständigen geladenen Quelle sind zusätzliche Pflichten. Der unendliche Kandidat braucht zusätzlich die konkrete native Calderón Antwort mit orientierter Einheitswindung und absolutem $c_3$, bevor von P1, physischem Zustand und Zeit gesprochen werden kann. [Definition und Normierung des nativen P1](https://github.com/sthamann/topological-fixed-point-theory/blob/a1273187b82488d34d2a474470152a41d42d4055/documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex#L48)

## Ausführung

Python 3 und SymPy. Vom Hauptverzeichnis des Repositorys aus, mit diesem Ordner als `p`:

```sh
p=experiments/exploration/2026-11-20-quaternion-frames/continuation
python3 "$p/seam_character.py"
python3 "$p/native_carrier.py"
python3 "$p/memory_tests.py" --native docs/research/fundamental-plan/evidence
python3 "$p/native_blocks.py" --native docs/research/fundamental-plan/evidence/RECURSIVE-WORD-SOURCE.json
python3 "$p/flag_time.py"
python3 "$p/flag_holonomy.py"
```
