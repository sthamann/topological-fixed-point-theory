# Auswahl des Fano Falls aus den regulären dicken sphärischen Rangzweigebäuden

**Satz.** Wählt man von einer festen Startkammer eine gleichverteilte Zielkammer und verlangt gleiche rohe Wahrscheinlichkeiten für sämtliche von der Identität verschiedenen Weylendpunkte nach Mischung innerer Beschriftungswechsel, dann gilt genau $m=3,\ q_s=q_t=2$. Die zwei nichttrivialen Drehungen und drei Spiegelungen entstehen nach der Auswahl. Die drei Möglichkeiten, welche Spiegelung als längster Endpunkt erscheint, erhalten jeweils Gewicht $1/3$.

## Voraussetzungen

Das endliche sphärische Rangzweigebäude hat zunächst einen beliebigen Diedertyp $I_2(m)$, $m\geq2$. Eine Tafel vom Typ $i$ enthält $q_i+1$ Kammern. Regularität und Dicke bedeuten, dass $q_i$ nur vom Typ abhängt und $q_i>1$ gilt. Die beiden Parameter dürfen verschieden sein, soweit das Gebäude dies erlaubt. Für den Endpunkt $w=\delta(c,d)$ gilt die natürliche Kammerzählung. Ein Beschriftungswechsel ersetzt $w$ durch $gwg^{-1}$ für ein beliebiges Element der Weylgruppe $W$. Dessen Verteilung ist beliebig, aber unabhängig von der Zielkammer. Mit „aktiv“ ist jedes $w\ne1$ gemeint. Die Gleichheit betrifft die tatsächlichen Wahrscheinlichkeiten nach dem Beschriftungswechsel. Eine gemeinsame Konditionierung auf Aktivität ändert das Kriterium nicht.

Gesetzt sind die Gebäudeklasse, die uniforme Ziehung, die Identitätsausnahme, die zugelassenen Beschriftungswechsel und das Kriterium gleicher Sichtbarkeit. Weder $m$ noch Panelgrößen, Zahl der Endpunkte, eine Fanoebene oder eine bestimmte Wurzelfamilie sind Eingaben. Die Anwendung dieses Auslesekriteriums in der Physik bleibt eine Auswahlannahme.

## Beweis

Für eine reduzierte Darstellung $w=s_1\cdots s_\ell$ enthält die Weylsphäre von $w$ genau $q_w=\prod_i q_{s_i}$ Kammern. Man zählt also Endkammern anhand eines fest gewählten reduzierten Worttyps und addiert keine Darstellungen desselben Weylelements. Insbesondere erzwingt die längste Zopfrelation für ungerades $m$ die Gleichheit $q_s=q_t$. Nach einer inneren Umbenennung bleiben die Summen über jede Konjugationsklasse erhalten. Daher ist gleiche aktive Rohsichtbarkeit nur möglich, wenn alle aktiven Klassen dieselbe mittlere Zahl von Kammern pro Element haben. Umgekehrt erreicht die Gleichverteilung über alle inneren Umbenennungen genau diese Klassenmittel. [J. Parkinson, Buildings and Hecke algebras, Proposition 2.1](https://extern.maths.usyd.edu.au/u/jamesp/1.pdf)

**Bei $m\geq4$** setze $r=st$ und $Q=q_sq_t>1$. Die Klasse einer Drehung $r^k$ enthält nur $r^k$ und $r^{-k}$, gegebenenfalls dasselbe Element. Die Zählungen für $r,r^{-1}$ sind $Q$. Für $r^2,r^{-2}$ sind die alternierenden Wörter der Länge vier reduziert und liefern $Q^2$. Bei $m=4$ ist $r^2$ das zentrale längste Element, dessen beide kürzeste Wörter denselben Wert $Q^2$ geben. Die beiden Klassen sind für jedes $m\geq4$ verschieden. Ihre pro Element erhaltenen Massen schließen Gleichheit aus.

**Bei $m=2$** ist $W$ abelsch. Die drei aktiven Endpunkte haben Gewichte $q_s,q_t,q_sq_t$. Das letzte Gewicht ist wegen Dicke strikt größer als jedes der ersten beiden.

**Bei $m=3$** gilt $q_s=q_t=q$. Beide Drehungen haben Rohgewicht $q^2$. Unter den drei Spiegelungen haben die beiden einfachen Gewicht $q$ und die längste Gewicht $q^3$. Sei $p_i$ die Gesamtwahrscheinlichkeit aller inneren Umbenennungen, unter denen die längste auf Spiegelung $i$ fällt. Das Rohgewicht an $i$ beträgt $q+(q^3-q)p_i$. Gleichheit mit $q^2$ erzwingt

$$
p_i=\frac{q^2-q}{q^3-q}=\frac1{q+1}.
$$

Wegen $\sum_i p_i=1$ folgt $q=2$ und $p_i=1/3$. Umgekehrt genügen genau diese drei Gesamtgewichte. Die einzelnen Wahrscheinlichkeiten der sechs inneren Umbenennungen lassen sich innerhalb der jeweiligen Zweierfaser frei aufteilen.

**Existenz.** Die eindimensionalen und zweidimensionalen Unterräume von $\mathbb F_2^3$ bilden ein projektives Inzidenzgebäude vom Typ $I_2(3)$. Sieben Punkte liegen auf je drei Geraden, also gibt es 21 inzidente Flaggen. Nach der beschriebenen Mischung trägt jeder der fünf aktiven Endpunkte genau $4/21$, bedingt auf Aktivität genau $1/5$. Die Startflagge hat Wahrscheinlichkeit $1/21$. Die drei Spiegelungen und die zwei nichttrivialen Drehungen bilden eine markierte $3+2$ Zerlegung. Diese Zählung liefert keine Amplitude für den physikalischen Rand.

## Kontrollen

| $m$ | $(q_s,q_t)$ | Kammern insgesamt | Verschiedene aktive Klassenmittel vor Division durch die Kammerzahl |
|---:|---|---:|---|
| 2 | $(2,2)$ | 9 | 2, 4 |
| 3 | $(2,2)$ | 21 | 4 |
| 3 | $(3,3)$ | 52 | 9, 11 |
| 4 | $(2,4)$ | 135 | 8, 10, 17, 64 |
| 6 | $(2,2)$ | 189 | 4, 14, 16, 64 |

Das Skript prüft zusätzlich sämtliche formalen Parametersysteme mit $2\leq m\leq12$ und $2\leq q_s,q_t\leq6$, mit der notwendigen Paritätsbedingung. Genau **1 von 175** besteht. Die Vergleichsmenge enthält auch formale Typen, deren Realisierbarkeit als endliches Gebäude nicht beansprucht wird. Der Beweis über $r$ und $r^2$ erfasst alle höheren $m$ unabhängig von dieser endlichen Kontrolle. Zielabhängige innere Umbenennungen erhalten ebenfalls die Klassensummen und somit die Auswahl von $m,q$. Die Identifikation der drei Gewichte $p_i$ verwendet die erklärte Unabhängigkeit.

Reproduktion mit der Python-Standardbibliothek über **python3 dihedral.py**. Eine zweite Darstellung der verwendeten Sphärenformel steht bei [Abramenko, Parkinson und Van Maldeghem, Gleichung (1.4)](https://www.maths.usyd.edu.au/u/jamesp/13.pdf).
