# Wann sehen alle aktiven Operationen gleich wahrscheinlich aus?

**Klassifikationssatz.** Unter allen endlichen regulären dicken sphärischen Gebäuden mit nichtkommutierender Coxetergruppe besitzt genau der Typ $A_2$ mit Panelparameter $q=2$ die nachstehende Gleichsichtbarkeit. Das zugehörige Gebäude ist bis auf Isomorphie die Fanoebene. Die Klassifikation umfasst beliebigen Rang, Produktgruppen und unterschiedliche reguläre Panelparameter. Sie setzt keinerlei TFPT Zahlen ein.

## Welche Wahrscheinlichkeit wird verglichen?

Sei $(W,S)$ das endliche Coxetersystem eines Gebäudes. Ein Panel vom Typ $s$ bietet neben der aktuellen Kammer genau $q_s$ Fortsetzungen. Regularität bedeutet, dass $q_s$ von der Kammer unabhängig ist. Dicke verlangt $q_s\geq2$. Bei gleichverteilten Kammern $F,G$ hat die relative Position $w=\delta(F,G)$ das unnormierte Gewicht $q_w=\prod_j q_{s_j}$ für eine reduzierte Darstellung $w=s_1\cdots s_k$. Das Produkt zählt die Endkammern im Abstand $w$ von festem $F$.

Wir wechseln zusätzlich die Bezeichnung der Gruppenoperationen durch eine gleichverteilte **innere** Automorphie. Das Gewicht pro einzelnem beschrifteten Endpunkt wird dann

$$
A(w)=\frac1{|W|}\sum_{g\in W}q_{gwg^{-1}}
=\frac1{|[w]|}\sum_{v\in[w]}q_v.
$$

Die Auswahlregel lautet $A(w)=k$ für **alle** $w\ne e$. Sie vergleicht einzelne aktive Operationen über sämtliche Konjugationsklassen hinweg. Die gesamte Masse einer Klasse beträgt dagegen $|[w]|k$. Endlichkeit, reguläre dicke Panels, freie gleichverteilte Kammerwahl, innere Präsentationswechsel und diese starke Gleichsichtbarkeit sind Modellannahmen. Weshalb eine physische Naht sie realisieren sollte, bleibt offen.

Auch eine beliebige, von Start und Ziel unabhängige Verteilung innerer Umbenennungen kann genau unter derselben Bedingung sämtliche aktiven Gewichte ausgleichen. Jede Umbenennung erhält die Summe jeder Konjugationsklasse. Gelingt der Ausgleich, sind deren pro Element gemittelte Werte daher schon gleich. Umgekehrt erreicht die uniforme innere Mittelung genau die Klassenmittel. Der Produktbeweis unten darf diese uniforme Mittelung verwenden, auch wenn die prüfweise zugelassene ursprüngliche Verteilung Korrelationen zwischen Komponenten enthält. Im einzig verbleibenden Fall müssen die drei effektiven Beschreibungen gleiche Gewichte tragen, die Verteilung innerhalb jeder Zweierfaser bleibt frei.

## Gemeinsamer Panelparameter

Zunächst sei $q_s=q\geq2$. Schreibe $r=|S|$, $N=|W|$ und $P(t)=\sum_{w\in W}t^{\ell(w)}$. Das Vorzeichen $\varepsilon(w)=(-1)^{\ell(w)}$ ist ein nichttrivialer Charakter und damit eine Klassenfunktion mit Summe null. Gleichsichtbarkeit liefert unmittelbar

$$
P(q)=1+(N-1)k,\qquad P(-q)=1-k.
$$

Wegen $k>1$ muss $P(-q)<0$ sein. Die Gleichungen erzwingen außerdem $P(q)/|P(-q)|=N-1+N/(k-1)>N-1$. Die Gradformel von Solomon mit den Invariantengraden $d_i$ ergibt

$$
P(t)=\prod_{i=1}^{r}(1+t+\dots+t^{d_i-1}),\qquad
\frac{P(q)}{|P(-q)|}=\left(\frac{q+1}{q-1}\right)^r
\prod_{d_i\text{ ungerade}}\frac{q^{d_i}-1}{q^{d_i}+1}\leq3^r.
$$

Somit ist $N-1<3^r$ notwendig. Der Faktor $3$ benutzt die ganzzahlige Mindestverzweigung zwei.

**Ab Rang vier.** Ein zusammenhängendes Coxeterdiagramm lässt eine Kette zusammenhängender Standardparabolgruppen $W_1\subset\cdots\subset W_r=W$ mit aufsteigenden Rängen zu. Im irreduziblen reellen Spiegelungsraum der Dimension $j$ fixiert $W_{j-1}$ einen von null verschiedenen Vektor $v$, da seine $j-1$ einfachen Spiegel gemeinsam eine Gerade festhalten. Der Orbit $W_jv$ spannt den Raum. Sein Schwerpunkt ist invariant und folglich null. Mit höchstens $j$ Orbitpunkten wäre dies eine nichttriviale lineare Abhängigkeit, die das Aufspannen verhindert. Also $[W_j:W_{j-1}]\geq|W_jv|\geq j+1$. Mit $|W_1|=2$ folgt $N\geq(r+1)!$. Schon $5!-1=119>3^4$ und das Verhältnis wächst für höhere Ränge. Damit ist jeder irreduzible Fall ab Rang vier ausgeschlossen.

**Rang drei.** Die irreduziblen Typen $A_3,B_3,H_3$ haben Ordnung $24,48,120$. Die letzten beiden verletzen $N-1<27$. Die Grade $2,3,4$ von $A_3$ ergeben $P(-q)>0$, ebenfalls unmöglich.

**Rang zwei.** Bei $I_2(m)$ ist $N=2m$. Für $m\geq5$ verletzt $2m-1\geq9$ die strikte Schranke. Bei $I_2(4)$ ergeben die Grade $2,4$ ein positives $P(-q)$. Es bleibt $I_2(3)=A_2$. Die drei Spiegelungen haben Längen $1,1,3$, die beiden Dreizyklen jeweils Länge $2$. Es gilt

$$
A(\text{Spiegelung})=(2q+q^3)/3,\qquad
A(\text{Dreizykel})=q^2.
$$

Gleichheit ist äquivalent zu $q(q-1)(q-2)=0$. Unter den dicken Parametern gewinnt genau $q=2$. Hier sind $P(2)=21$ und $k=4$. Ein dickes sphärisches Gebäude vom Typ $A_2$ ist eine projektive Ebene. Die projektive Ebene der Ordnung zwei ist eindeutig und besitzt 21 inzidente Punkt-Geraden-Paare als Kammern.

## Ungleiche Verzweigung, Produktgruppen und Randfälle

Bei ungeradem $m_{st}$ geben die beiden alternierenden reduzierten Ausdrücke für das längste Element von $\langle s,t\rangle$ identische Kammerzahlen. Da die Buchstabenhäufigkeiten um eins abweichen, folgt $q_s=q_t$. Nach der Klassifikation endlicher irreduzibler Coxeterdiagramme sind unterschiedliche Parameter damit nur in $B_r$, $F_4$ und $I_2(2j)$ möglich. In allen diesen Fällen ist das längste Element $w_0$ zentral. Für jeden anderen Endpunkt gilt $\ell(w^{-1}w_0)=\ell(w_0)-\ell(w)>0$. Die reduzierte Fortsetzung liefert $q_{w_0}=q_w q_{w^{-1}w_0}>q_w$. Die Einerklasse von $w_0$ besitzt daher ein strikt größeres Mittel als jede andere Klasse. Das schließt diese weiteren Fälle aus.

Bei einem echten Produkt $W'\times W''$ faktorisiert das Gewicht und die innere Mittelung wirkt komponentenweise. Zwei nichttriviale Endpunkte $(u,e),(e,v)$ müssten beide das Gewicht $k>1$ tragen. Für $(u,v)$ wäre es $k^2$. Daher scheitert jeder Typ mit mindestens zwei nichttrivialen Faktoren, auch mit unterschiedlichen Dicken. $A_1$ erfüllt die Bedingung für jedes $q$ trivial bei einem einzigen aktiven Endpunkt. Der geforderte nichtkommutierende Austausch schließt ihn aus. Beim leeren Coxetersystem gibt es keine aktiven Endpunkte. Ohne Dicke, etwa beim dünnen Apartment mit allen $q_s=1$, werden beliebig viele weitere Systeme gleich sichtbar.

## Messgröße und Gegenproben

| Vergleich | Ergebnis |
|---|---|
| $A_2,q=2$, vor innerer Mittelung | Gewichte nach Längen $1,2,2,4,4,8$. |
| $A_2,q=2$, nach innerer Mittelung | Gewichte $1,4,4,4,4,4$. Leerschritt $1/21$. Jede einzelne aktive Operation hat bedingt auf Aktivität Wahrscheinlichkeit $1/5$. |
| $A_2,q=3$, bloße Konjugationsinvarianz | Drei Spiegelungen besitzen einzeln Gewicht $11$, die zwei Dreizyklen jeweils $9$. Die Verteilung ist vollständig konjugationsinvariant und verletzt die geforderte Gleichsichtbarkeit. |
| Separate Bedingung erfolgreicher Rückgalerie im Typ $A_2$ | Unabhängig vom regulären $q$ kompensiert $q^{-\ell(w)}$ das Ausgangsgewicht $q^{\ell(w)}$. Alle sechs relativen Positionen, einschließlich der leeren, erhalten danach je $1/6$. |

Die letzte Zeile beschreibt eine alternative Auslesung desselben Kammerpaares. Der Rücktest verwendet die zur tatsächlichen Weylposition passende reduzierte Galerie. Ein Wechsel der Präsentation muss die zugehörigen Generatoren mittransportieren. Eine nur am bereits gemittelten Label nachträglich angesetzte, unkorrelierte Rückkehrwahrscheinlichkeit ist eine andere Operation.

Der Unterschied ist am Quantenkanal messbar, sobald der leere Zweig die Identität und ein aktiver Zweig die vorgegebene nichtselektive Cliffordmessung $D_j=(\mathrm{id}+\operatorname{Ad}_{\Gamma_j})/2$ ausführt. Aus den **rohen** gemittelten Gewichten folgt $(11\,\mathrm{id}+2\sum_j\operatorname{Ad}_{\Gamma_j})/21$ mit Sektoreigenwerten $1,5/21,13/21$. Erst die **Rückselektion** ergibt $(7\,\mathrm{id}+\sum_j\operatorname{Ad}_{\Gamma_j})/12$ und die markierten Sektorwerte $1,1/3,2/3$. Die Klassifikation entscheidet die endliche Geometrie in der angegebenen Klasse. Die Wahl der Messoperation und ihre Identifikation als physische P1 Auslesung bleiben zusätzliche Schritte.

Reproduktion im Ordner dieser Datei mit `python3 coxeter_sichtbarkeit.py` (nur Standardbibliothek). Das Programm berechnet reduzierte Längen, Kammergewichte und vollständige wirkliche Konjugationsklassen für die kleinen Kontrollfälle. Die dihedrale Kontrolle deckt 112 Kombinationen von $m=3,\ldots,18$ und $q=2,\ldots,8$ ab. Die allgemeine Aussage folgt aus dem obigen Beweis.

## Literatur

[Abramenko, Parkinson und Van Maldeghem, Kammerzählung in regulären Gebäuden](https://www.sciencedirect.com/science/article/pii/S0021869317301308) · [Davis, The Geometry and Topology of Coxeter Groups](https://people.math.osu.edu/davis.12/second_edition.pdf) · [Bredon, Solomons Wachstumsformel und Invariantengrade, Kapitel 5](https://homeweb.unifr.ch/kellerha/pub/NBredon-Thesis.pdf) · [Spencer, Konjugation einfacher Spiegel und zentrale längste Elemente, Propositionen 7.1.20 und 7.1.21](https://webhomes.maths.ed.ac.uk/~igordon/Spencerthesis.pdf) · [Projective and affine planes, Eindeutigkeit der Fanoebene](https://maths.qmul.ac.uk/~lsoicher/designtheory.org/library/encyc/topics/pap.pdf)
