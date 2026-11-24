# Welchen zweiten Zeitschritt erlaubt der markierte Cliffordkanal?

**Ergebnis.** Wird der erste Schritt exakt als markierter Kanal gemessen, muss ein stationärer euklidischer Prozess mit vollständig positiven Zwischenkarten beim zweiten Schritt messbar vom zweimal frisch ausgeführten Kanal abweichen. Die garantierte Abweichung beträgt mindestens **0,0621126793** in der unten definierten Operatornorm. Schon vollständige Positivität zur halben Schrittzeit erzwingt mindestens **0,0319096754**. Die beiden Schranken betreffen alle Dimensionen eines möglichen Gedächtnisses.

## Voraussetzungen und Normierung

Auf dem Operatorraum $\mathcal K=M_4(\mathbb C)$ gilt das Skalarprodukt $\langle X,Y\rangle=\operatorname{Tr}(X^\dagger Y)/4$. Die fünf hermiteschen Cliffordgeneratoren erfüllen $\Gamma_i\Gamma_j+\Gamma_j\Gamma_i=2\delta_{ij}I$. Die orthogonalen Projektoren $P_0,P_5,P_{10}$ bezeichnen die skalare Richtung, die fünf Generatoren und die zehn Richtungen $i\Gamma_j\Gamma_k$ mit $j<k$. Die Beobachtung nach einer Zeiteinheit sei

$$
K_1=\Phi=\frac{7\,\mathrm{id}+\sum_j\operatorname{Ad}_{\Gamma_j}}{12}
=P_0+aP_5+bP_{10},\qquad a=\tfrac13,\quad b=\tfrac23.
$$

Ein stationärer euklidischer Prozess hat hier die präzise Bedeutung $K_t=V^\dagger e^{-tH}V$ für dieselbe Isometrie $V$ zu allen Zeiten und einen selbstadjungierten Operator $H\geq0$. Die Erweiterung ist ein Hilbertraum. Ein quantenphysischer Generator auf einer erweiterten Observablenalgebra wäre eine zusätzliche Aussage. $K_1(I)=I$ und $e^{-H}\leq1$ implizieren $e^{-H}VI=VI$. Daher gilt $K_t(I)=I$ für jedes $t$, und Selbstadjungiertheit von $K_t$ auf $\mathcal K$ liefert Spurerhaltung. Vollständige Positivität einer Zwischenkarte auf $M_4$ wird gesondert vorausgesetzt. Alle Fehler werden in der von $\mathcal K$ induzierten Operatornorm gemessen.

## Exakter Abschluss scheitert an der halben Zeit

Mit $T=e^{-H}$ und $P=VV^\dagger$ gilt

$$
K_2-K_1^2=V^\dagger T(1-P)TV\geq0.
$$

Falls $K_2=\Phi^2$, folgt $(1-P)TV=0$. Das Bild von $V$ ist invariant und wegen $T=T^\dagger$ reduzierend. Der Funktionalkalkül liefert $K_{1/2}=\sqrt{\Phi}$. Für eine über $\mathrm{Spin}(5)$ gemittelte, unital wirkende Karte $P_0+uP_5+vP_{10}$ lauten die normierten Choiwerte pro Pauliwort

$$
p_0=(1+5u+10v)/16,\quad p_5=(1-3u+2v)/16,\quad p_{10}=(1+u-2v)/16.
$$

Die Vielfachheiten betragen $1,5,10$. Bei der positiven Quadratwurzel ist $16p_{10}=1+\sqrt{1/3}-2\sqrt{2/3}<0$. Damit ist exakter Abschluss mit einer vollständig positiven Halbzeitkarte ausgeschlossen. Dieser Schluss setzt keine Symmetrie des unbekannten $H$ oder der Isometrie voraus.

## Quantitative Schranke

Für jeden Sektor vom Rang $d\in\{5,10\}$ sei $m_d(t)=\operatorname{Tr}_{\mathcal K}(P_dK_t)/d$. Die Spektralzerlegung von $e^{-H}$ ergibt $m_d(t)=\int_0^1 x^t\,d\nu_d(x)$ für Wahrscheinlichkeitsmaße $\nu_d$. Die Mittelung der Karte über Spin(5) erhält vollständige Positivität und diese Sektormittelwerte. Die drei Darstellungen der Ränge $1,5,10$ treten jeweils einmal auf, sodass ihr gemittelter Operator $P_0+m_5(t)P_5+m_{10}(t)P_{10}$ ist.

**Nur Halbzeitpositivität.** Jensen liefert $m_5(1/2)\leq\sqrt a$. Die Choiungleichung $p_{10}\geq0$ erzwingt $m_{10}(1/2)\leq(1+\sqrt a)/2$. Hölder ergibt $m_{10}(1)\leq m_{10}(1/2)^{2/3}m_{10}(2)^{1/3}$. Daraus folgt

$$
\|K_2-\Phi^2\|\geq m_{10}(2)-b^2
\geq\frac{b^3}{[(1+\sqrt a)/2]^2}-b^2
=\frac{28-16\sqrt3}{9}.
$$

Die Konstante ist unter **dieser** Voraussetzung ein scharfes Infimum. Für den Operator $e^{-H/2}$ wähle auf $P_5$ deterministisch $\sqrt a$. Auf jeder $P_{10}$ Richtung verwende eine Zweipunktverteilung mit Mittel $m=(1+\sqrt a)/2$ und zweitem Moment $b$. Ein unterer Punkt $\epsilon>0$ liefert den oberen Punkt $h=(b-m\epsilon)/(m-\epsilon)$ mit Gewicht $(m-\epsilon)/(h-\epsilon)$. Für hinreichend kleines $\epsilon$ liegen beide Punkte in $(0,1]$. Diagonale Operatoren und der mit den Quadratwurzeln der Wahrscheinlichkeiten belegte Einbettungsvektor realisieren diese Maße. Die Halbzeitkarte liegt exakt auf dem positiven Choirand. Für $\epsilon\downarrow0$ konvergiert ihr Zweizeitfehler zur Schranke. Die Konstruktion sichert vollständige Positivität nur zur geprüften Halbzeit und nach dem ersten Schritt.

**Vollständige Positivität zu allen Zeiten.** Setze $f(t)=1+m_5(t)-2m_{10}(t)$. Durch die Choiungleichung gilt $f(t)\geq0$, und $f(1)=1+a-2b=0$. Deshalb gilt $f'(1)=0$. Die Ableitung existiert, da $x\log x$ bei null stetig fortsetzbar ist. Mit der Konvexität dieser Funktion erhält man $m_{10}'(1)=m_5'(1)/2\geq a\log(a)/2$. Unter dem Wahrscheinlichkeitsmaß $d\widetilde\nu=x\,d\nu_{10}/b$ ist $\int\log(x)\,d\widetilde\nu=m_{10}'(1)/b$. Jensen liefert

$$
\boxed{\displaystyle\|K_2-\Phi^2\|\geq m_{10}(2)-b^2\geq b\exp\!\left(\frac{a\log a}{2b}\right)-b^2=\frac23\,3^{-1/4}-\frac49.}
$$

Dafür genügt vollständige Positivität in einer offenen Umgebung des ersten Schritts. Die optimale Allzeitkonstante ist offen.

## Ausführbare Kontrollen und Gegenmodelle

| Prozess | Sektorwerte von $K_2$ auf $P_5,P_{10}$ | Zweizeitfehler | Aussage |
|---|---|---:|---|
| Frische, unabhängige Ausführung von $\Phi$ | $1/9,4/9$ | $0$ | Exakte diskrete Komposition. Die lineare Interpolation der frischen periodischen Schritte verletzt euklidische Positivität. |
| Zwei Energien $0,\log3$ | $1/9,5/9$ | $1/9$ | Stationäre euklidische Kompression und alle Zwischenkarten vollständig positiv. |
| Eine dauerhaft gespeicherte Wahl aus fünf Messrichtungen | $2/9,11/18$ | $1/6$ | Explizite Erweiterung mit fünf klassischen Gedächtniswerten und reversibler Dephasierung in jedem Zweig. |

Für die zweite Zeile setze $K_t=P_0+3^{-t}P_5+(1+3^{-t})P_{10}/2$. Der Energiesektor $P_{10}$ sieht je zur Hälfte $0$ und $\log3$. Die normierten Choiwerte sind $(3+5x)/8,(1-x)/8,0$ mit $x=3^{-t}$. Für die dritte Zeile verwende die fünf nichtselektiven Messkarten $D_j=(\mathrm{id}+\operatorname{Ad}_{\Gamma_j})/2$. Bei festem $j$ ist $D_j+6^{-t}(\mathrm{id}-D_j)$ eine dephasierende Halbgruppe. Ein gleichverteiltes, zeitlich festes $j$ liefert die angegebene Kompression auf der direkten Summe. Beide Konstruktionen ergeben denselben ersten Schritt. In der ersten Zeile liefert die lineare Interpolation auf dem Fünfersektor bei $t=1/4$ den negativen Spiegeltest $K_{2t}-K_t^2=-1/36$ auf dieser Richtung.

Reproduktion im Ordner dieser Datei mit `python3 continuous_time.py` (SymPy). Das Skript kontrolliert die exakten Choiwerte, beide Schranken, die zwei kontinuierlichen Konstruktionen und den negativen Interpolationstest.

## Physische Bedeutung und Quellenabgleich

Ein Experiment braucht den markierten ersten Schritt **und** eine zweite Messung mit derselben Zustandspräparation, derselben Zeitdefinition und derselben Einbettung des beobachteten Raums. Der relevante Vergleichswert ist die mittlere Rückkehr der zehn Bivektorrichtungen. Bei der unabhängigen Wiederholung beträgt sie $4/9$. Unter allen durchgehend vollständig positiven stationären euklidischen Modellen muss sie mindestens $(2/3)3^{-1/4}$ betragen. Der vorausgesetzte erste Schritt fixiert die Zeiteinheit. Der aktuelle physische TFPT Prozess hat diese drei Identifikationen und die gemessene Zweizeitantwort noch nicht gemeinsam geliefert.

Im Atlas ist die Kompressionsidentität samt Reduktion bei verschwindendem Defekt bereits unter `finding:original-os-finite-jet-is-compression-with-single-time-closure-test-20261112` dokumentiert. Die Eigenwertzuordnung des markierten Kanals und der nichtpositive Zwischenwert stehen im Quellenvergleich vom 21. September. Ein weiterer Befund schließt auf demselben $M_4$ alle homogenen GKSL Logarithmen sämtlicher positiver ganzzahliger Potenzen aus. Die hier abgeleiteten quantitativen Schranken gelten auch für beliebig große stationäre euklidische Erweiterungen unter den genannten Positivitätsannahmen. In der gezielten Suche fand sich hierfür kein früherer Eintrag.

[Markierter Kanal und zusätzlich angenommene Zuordnung der Orbits](../../../theory-contracts/marked-seam-cp-lift-20260921/SOURCE_PROVENANCE.md) · [Kompression des ursprünglichen OS Transfers](../../../theory-contracts/p1-context-source-interface-20261112/PROOF.txt) · [Ausschluss homogener GKSL Logarithmen](../../../theory-contracts/source-continuous-marked-time-search-20260922/PROOF.md)
