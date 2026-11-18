# Ursprüngliche Quelle und regionale Bewegung

## Gewichteter Beitrag auf einer festen Geschichte

Der R032 Vergleich bei $p=1/4$, $t=4/5$, $x=1/2$ verwendet die 42
geladenen Paare des nativen v252 Charts und die spiegelgerade Folge
$(+x,-x,-x,+x)$. Für jedes Paar ergibt die Verklebung $s=e^{-tE}$,
wobei $E=\sqrt{p^2+x^2y^2}$. Das Verhältnis des leeren Exteriorsektors
gegen S ist das Produkt der $e^{-t(E-p)}$. Die Spur enthält zusätzlich
das Produkt der $(2+2\cosh g)/(2+2\cosh(tp))$, mit
$g=2\operatorname{arsinh}((p/E)\sinh(tE/2))$.

| Logarithmischer Quotient | Ursprungsgeschichte |
| --- | ---: |
| Verklebter Skalar gegen S | −19,73994045338434 |
| Unprojizierte geladene Spur gegen S | −19,71988306021885 |
| Verklebter Skalar gegen nachträgliche polare Wahl | −19,54420754979753 |
| Skalar und unprojizierte Spur jeweils gegen T | 0 |

OriginLab berechnet die Intervalle mit Arb aus den originalen
Yukawamatrizen und prüft die veröffentlichten Rundungen der Spur und des
polaren Vergleichs. Bei T stimmen die bedingten Matrizen und Teilskalare
überein. Für das vollständige R113 Integral fehlen die operatorwertigen
bosonischen Enddaten samt Historiengewichten, der physische Selektor,
die Restfaktoren und die gemeinsame absolute Leerreferenz. Die absolute
Antwort $Z_0$, $F(0)$ und $C_{UY}$ bleiben offen.
[Originaler R113 Vergleich](../../experiments/theory-contracts/chat-triangulation-20261108/R113-ORIGINAL-UND-REFERENZ-20261109.txt)

## Regionale Generatoren und neutraler Rückweg

In $x=\tan(\theta/2)$ lauten die vier Viertelkonjugate im gemeinsamen
Koordinatenfeld $x^2,(x-1)^2/2,1,(x+1)^2/2$. Sie haben linearen
Rang drei. Die Lieklammer $[a\partial_x,b\partial_x]=(ab'-ba')\partial_x$
hat als lineare Abbildung auf dem Keilquadrat in der Basis
$(1\wedge x,1\wedge x^2,x\wedge x^2)$ Determinante 2. Ein
kommutierendes Paar hat also Keilprodukt null. Im gesamten bekannten
Möbius Spann gibt es höchstens eine unabhängige kommutierende
Richtung. Eine Erweiterung durch räumlich ausgedehnte oder gemischte
Randprofile bleibt zu prüfen.
[Originaler Trägertransport](../../docs/research/fundamental-plan/TFPT-TRAEGERTRANSPORT-PICARD-KUBIK-20261111.txt)

Für die volle Wurzelladung $\gamma$ und den longitudinalen Grad eins
gilt auf demselben invarianten Bereich $F^*F=I+E^*E$.
Der fallende Operator $F$ muss injektiv in den neutralen Grad null
abbilden. Sei dort der Grundgrad nichtnegativ. Dann gilt auf dem neutralen
Grad null $F=0$ und $FE=I$. Damit ist $F$ auf dem geladenen Grad eins
ein Isomorphismus zum neutralen Grad null. Für zwei zusätzliche
orthogonale geladene Richtungen werden zwei zusätzliche neutrale
Rückwege benötigt. Der bekannte Grad zwei enthält einen neutralen
Weylrückweg vom Rang acht, dessen geprüfte Ableitung auf der einzelnen
Wurzellinie bleibt. Das liefert einen Kontrollwert für neue regionale
Einbettungen.
[Regionale Komposition](../../docs/research/fundamental-plan/TFPT-REGIONALE-KOMPOSITION-20261112.txt)

## Prüfungen aus den parallelen Sitzungen

Die [drei Prüfläufe](../../experiments/theory-contracts/f1-parallel-session-checks-20261114/README.md)
liefern Kontrollfälle für endliche Verklebung, Uhrholonomie und einen
chiralen Index. Der gezielte nächste Originaltest liest den neutralen
Grad null im operatorwertigen R113 Leerfermionblock. Eine neue
Querbewegung muss eine ganze Trägeralgebra im selben positiven Zustand
und unter demselben Zeittransport verschieben. Der geordnete geladene
Übergang und seine Referenz müssen dabei erhalten bleiben.

Atlas `problem:p0-common-geometric-and-charged-process-20261112`,
Revision `dace8ddec7ef89c319e6bfed`. F1 bleibt offen.

## Punktweise Querprofile mit mehrdimensionalem Grundgrad

Betrachte den untwisted affinen Strom mit gemeinsamer longitudinaler
Graduierung $L_0\geq 0$, einem gemeinsamen invarianten Produktbereich
und unitärer, integrierbarer Wirkung der konstanten $E_8$ Ströme auf
Niveau eins. Die Profile bilden eine unitale kommutative $*$ Algebra.
Ihre Stromklammer benutzt das punktweise Produkt und eine lineare,
zentrale $*$ Aktion $C$ mit $C(1)=I$. Für die volle Wurzel $\gamma$ gilt
in der Konvention der [regionalen Komposition](../../docs/research/fundamental-plan/TFPT-REGIONALE-KOMPOSITION-20261112.txt)

$$[F_1(a^*),E_{-1}(b)]=C(a^*b)-H_0(a^*b),\qquad E_{-1}(a)^*=F_1(a^*).$$

Sei $M=\ker L_0$. Positive Moden töten $M$. Die konstanten Nullmoden
können auf $M$ nur dominante $E_8$ Gewichte der affinen Stufe eins
tragen. Die exakten höchsten Wurzelkoeffizienten aus
[v844](../../verification/v844_message_doily_rank.py) sind
$(2,2,3,3,4,4,5,6)$ bis auf Reihenfolge. Die Niveaugrenze
$\langle\lambda,\theta^\vee\rangle\leq1$ erzwingt daher $\lambda=0$.
Alle profilierten Nullmoden töten $M$ ebenfalls, denn
$[X_0(1),Y_0(a)]=[X,Y]_0(a)$ und $[\mathfrak e_8,\mathfrak e_8]=\mathfrak e_8$.

Die konstante Niveau eins Relation macht $E_{-1}(1)$ zu einer Isometrie
von $M$ auf den **gesamten** Grad eins mit voller Ladung $\gamma$.
Für jeden dortigen Vektor $v$ gilt
$\|F_1(1)v\|^2=\|v\|^2+\|E_{-1}(1)v\|^2$.
$v-E_{-1}(1)F_1(1)v$ liegt im Kern von $F_1(1)$ und ist daher null.
So folgt direkt

$$E_{-1}(a)|_M=E_{-1}(1)C(a)|_M,\qquad C(a^*b)|_M=C(a)^*C(b)|_M.$$

Auf dem von $M$ durch die konstanten Ströme erzeugten Bereich erzeugen
die Profile keine weitere Stromdarstellung. Die profilierten Nullmoden
kommutiert man durch die konstanten negativen Moden. Auf $M$ verschwinden
sie. Durch Induktion entsteht für $X_0(a)$ auf jedem konstant erzeugten
Wort dieselbe Wirkung wie $C(a)X_0(1)$. Für beliebiges $n$ folgt
$X_n(a)=C(a)X_n(1)$ aus $[Y_0(a),Z_n(1)]=[Y,Z]_n(a)$ und der
Perfektheit. Zentrale $C$ kommutieren hierbei mit der Graduierung und
allen benutzten Strommoden. Bei beschränkter zentraler Wirkung oder
endlichem $M$ ist die Wirkung von $C$ auf dem Grundgrad kommutativ.
Im endlichen Fall zerfällt sie in gemeinsame Charakterräume. Jede
stetige unitäre Bewegung, die diese Algebra als Familie erhält, fixiert
deren endlich viele Minimalprojektoren, weil sie sie nur permutieren kann.
Unendliche zentrale Algebren können stetig normalisiert werden. Daraus
allein folgt weder ein regionaler Transport noch eine Kopplung.

Für eine gewichtete R113 Kandidatenform seien $G$ der neutrale Gram,
$K_{ab}$ der Gram der geladenen Bilder und $0$ das Einheitsprofil.
Auf dem Nullquotienten misst OriginLab
$K_{ab}-K_{a0}G^+K_{0b}$ vor äußerer Auswahl. Der gemeinsame positive
Rest ist die Gramform der geladenen Anteile außerhalb des gemessenen
Einheitsrückwegs. Sein Rang gibt eine Untergrenze für zusätzliche
neutrale Rückwege, wenn die Niveau eins Voraussetzungen gelten. Bei
vollständigem Grundgrad schließt jeder positive Rest die Kandidatenklasse
aus. Der Test prüft außerdem die angegebenen punktweisen Produktrelationen.
Zwei Auswertungssektoren für $s^2=1$ bestehen mit Wirkung $\mathrm{diag}(1,-1)$.
Ein skalares Mittel mit $C(s)=0$ und $C(s^2)=1$ erfüllt den gemessenen
Ward Produktteil, besitzt aber einen geladenen Rest vom Rang eins.
Dieser Rest kann bei nur einem neutralen Grundzustand nicht zurückkehren.
[Code](affine_profiles.py) und [R113 Verbraucher](source.py).

Für den Anschluss an räumliche Dynamik ist nun eine geordnete Antwort
zwischen zwei neutralen Profilsektoren auf dem ursprünglichen
R113 Leerfermionblock ein direkter Test. Bei kommutierendem bosonischem
Transfer ist $P_i C_t^{00} P_j=0$ für verschiedene Profilsektoren.
Ein nichtverschwindender Wert auf derselben gewichteten Form wäre eine
Mischung durch den ursprünglichen Zeittransport. Alternativ muss eine
regionale Einfügung die punktweise Stromfortsetzung verlassen. In
beiden Fällen bleibt die Bewegung ganzer lokaler Algebren eigens zu prüfen.
Weder der native Aufruf noch die endlichen Kontrollen liefern diesen Wert.
