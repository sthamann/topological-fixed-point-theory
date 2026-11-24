# Zwei gemessene Flaggen bestimmen ein Teilgitter

**Ergebnis.** Ein geordnetes Paar gegenüberliegender Fano Flaggen bestimmt genau ein Apartment und daraus einen von sieben eingebetteten Familienräumen. Im Modell kann damit ein einzelnes protokolliertes Rückkehrereignis das Teilgitter A₃⊕D₅ auswählen. Die Mittelung über die Ereignismarke lässt die Auswahl verschwinden.

## Die Marke aus dem Rückkehrereignis

Eine Flagge besteht aus einem Punkt p und einer inzidenten Linie ℓ. Die beiden Flaggen (p,ℓ) und (q,m) sind gegenüberliegend, wenn p nicht auf m und q nicht auf ℓ liegt. Sei r der einzige Schnittpunkt der Linien. Dann sind p,r,q drei nichtkollineare Punkte. Sie bilden ein Apartment, also ein Dreieck einschließlich seiner sechs inzidenten Eckflaggen. Jedes Apartment, das die beiden Gegenflaggen enthält, muss p und q als Ecken und ℓ und m als Seiten enthalten. Seine dritte Ecke ist daher r.

Alle Punkte sind von null verschiedene Vektoren in F₂³. Die Vektorsumme

$$
u=p+r+q
$$

ist von null verschieden und unabhängig von der Reihenfolge der drei Ecken. Unter jedem linearen Koordinatenwechsel g wird sie zu g(u). Die Beziehung benutzt die gemessenen Flaggen und die projektive Inzidenz. Sie benötigt kein vorab gewähltes Apartment. Die Addition ist intrinsisch rekonstruierbar, denn für zwei verschiedene Fano Punkte ist ihre Summe der dritte Punkt ihrer Linie.

Es gibt 21 Flaggen und acht Gegenflaggen für jede. Die 168 geordneten Gegenpaare ergeben 28 Apartments mit je sechs Paaren. Genau vier Apartments haben dieselbe Marke u. Damit treten alle sieben u jeweils 24 Mal auf.

## Welche Gleichverteilung der Rücklauf liefert

Im verwendeten Rückkehrprotokoll zieht man zunächst beide Flaggen unabhängig und gleichverteilt. Zu jeder elementaren Panelbewegung gibt es zwei gleichgewichtete Fortsetzungen. Für gegenüberliegende Flaggen erreicht unter einer festgelegten kürzesten Panelreihenfolge genau eine der acht dreischrittigen Galerien die gewünschte Gegenflagge. Dies gilt für beide reduzierten Wörter aba und bab. Die bedingte Rückkehrchance eines solchen Paares ist deshalb 1/8. Das Protokoll zählt den längsten relativen Typ einmal.

| Ereignis | Wahrscheinlichkeit |
|---|---:|
| Längster Typ und erfolgreiche Rückkehr | 1/21 |
| Bestimmtes geordnetes Gegenpaar, bedingt auf diesen erfolgreichen Typ | 1/168 |
| Bestimmtes Apartment unter derselben Bedingung | 1/28 |
| Bestimmte Marke u unter derselben Bedingung | 1/7 |

Hält man bereits die erste Flagge fest, bleiben genau vier mögliche u übrig, die vier Punkte außerhalb ihrer Linie. Jedes dieser u hat dann Wahrscheinlichkeit 1/4. Die Gleichverteilung über alle sieben setzt die Mittelung über die erste Flagge voraus.

## Der bedingte Familienraum

Jedes u paart die sechs alten Punkte x≠u zu drei Paaren {x,x+u}. Die Koordinate 0 ist die zusätzliche Prüfkoordinate des Fano Codes. Sie bildet mit u das Referenzpaar. Die Auslesung summiert über jedes der drei aktiven Paare. Ihr unsichtbarer Raum hat die orthonormale Basis

$$
\frac{e_x-e_{x+u}}{\sqrt{2}},
\qquad \text{je ein Vertreter }x\text{ pro aktivem Paar}.
$$

Für festes u ist dies der Familienraum F(u). Das in GITTER.md hergeleitete E₈ Gitter schneidet F(u) in einem primitiven D₃≅A₃ mit zwölf Wurzeln und sein orthogonales Komplement in einem primitiven D₅ mit vierzig Wurzeln. Der Quotient des Gesamtgitters durch die direkte Summe ist zyklisch von Ordnung vier. Vier verschiedene Apartments mit gleichem u führen zu denselben Teilgittern, denn die Paarung hängt allein von u ab. Verschiedene u geben verschiedene Familienräume. Der Punkt u ist nämlich der einzige alte Punkt, an dem jede Familienrichtung verschwindet. Die zwölf Wurzeln spannen den jeweiligen Raum auf. Das Ereignisverfahren erreicht somit **sieben der 7 560 eingebetteten A₃ Wurzelsysteme**. Diese Vergleichszahl zählt Untergitter im fest beschrifteten E₈, keine physikalisch verschiedenen Welten.

## Was ohne die Marke messbar bleibt

Sei P(u) der orthogonale Projektor auf F(u). Jedes ungeordnete Paar verschiedener alter Punkte x,y kommt genau bei u=x+y vor. Die Summe aller Kantenprojektoren ist deshalb der halbe Graphlaplaceoperator des vollständigen Graphen auf sieben Punkten. Auf den alten Punktkoordinaten folgt mit der Einsmatrix J₇

$$
\frac{1}{7}\sum_{u\ne0}P(u)
=\frac12\left(I_7-\frac{J_7}{7}\right).
$$

Auf der neuen Prüfkoordinate ist der Mittelwert null. Das gesamte Spektrum hat zwei Nullwerte und sechsmal den Wert 1/2. Die alten Richtungen mit Komponentensumme null bleiben gleichwertig. Ein fester dreidimensionaler Familienraum ist aus diesem gemittelten Projektor nicht bestimmbar. Die konditionierte Zuordnung ist unter allen 168 Elementen von GL(3,2) kovariant; das Mittel ist invariant.

**Physische Reichweite.** Für eine Auswahl pro Ereignis reichen die bereits definierten Flaggen und ihre erfolgreiche Rückkehr des längsten Typs. Zusätzlich interpretiert die Konstruktion die Differenzen innerhalb derselben Komplementklasse als Familienrichtungen. Ein über viele Ereignisse persistenter Sektor erfordert eine Dynamik des ausgelesenen u oder einen Zustandsnachweis, der diese Marke behält. Aus der hier berechneten Einereignisverteilung ergibt sich weder eine Persistenzzeit noch die Auswahl eines kosmischen Vakuums.

**Reproduktion.** `python3 apartment.py` enumeriert Flaggen, Galerien, Apartments, die 240 Konstruktion A Wurzeln und alle linearen Symmetrien exakt mit der Python Standardbibliothek. Die Gesamtzahl der 7 560 A₃ Einbettungen wurde in `fano_lattice.py` separat über alle einfachen A₃ Wurzeltripel gezählt.
