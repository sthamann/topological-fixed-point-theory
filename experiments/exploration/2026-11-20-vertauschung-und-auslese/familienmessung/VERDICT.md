# Lässt sich die Familiengeometrie aus ihrer Messgeschichte zurückgewinnen?

**KANDIDAT, Stufe 1.** Aus drei aufeinanderfolgenden Marken lässt sich in diesem Prozess die Fano Ebene rekonstruieren. Der bisherige Vierertransfer trägt in einer gewählten Spinorauslesung die Geometrie ebenfalls bis zur dritten Messzeit weiter. Die ersten beiden unmarkierten Transfers enthalten bereits Gedächtnis, lassen die Inzidenztabelle aber unbestimmt. Das Modell liefert eine kontinuierliche, reversible Dynamik für eine physikalische Zusatzannahme über Wartezeiten. Die ursprüngliche Zeitspiegelung und eine orientierte Windung aus dem Prozess bleiben offen.

## Welches Experiment wurde verändert?

Das [vorherige Experiment](../fortsetzung/APARTMENT.md) wählt zwei gegenüberliegende Flaggen und misst einen von sieben Familienräumen. Jetzt wird beim nächsten Ereignis die Zielflagge als neue Ausgangsflagge übernommen. Dies ist eine neue, überprüfbare Annahme. Alternativen sind unabhängiges Neuziehen und dauerndes Festhalten derselben Marke. Die bedingte Gegenflagge ist gleichverteilt unter acht Möglichkeiten. Der optionale erfolgreiche Rücktest hat für jeden solchen Übergang dieselbe Wahrscheinlichkeit 1/8. Wer den Rückweg tatsächlich läuft, muss die Zielflagge separat speichern, um das hier untersuchte Protokoll auszuführen.

Für einen Übergang von (p,L) nach (q,M) ist r der Schnittpunkt beider Geraden. Als Marke wird u=p+r+q im binären Dreierraum gelesen. Die Rechnung verwendet die symmetrischen 21 mal 21 Übergangsmatrizen A_u mit Einträgen 0 oder 1/8. Ihre Summe T ist stochastisch. Anfangszustand ist die stationäre Gleichverteilung. Für ein Wort ist seine Wahrscheinlichkeit die Summe aller Einträge des entsprechenden Matrizenprodukts geteilt durch 21.

## Die erste zeitlich beobachtbare Inzidenz

Nach einer beobachteten Marke ist die Wiederholungswahrscheinlichkeit 1/4. Unabhängige Marken liefern 1/7, eingefrorene Marken eins. Für zwei verschiedene letzte Marken a,b gilt die vollständige bedingte Verteilung:

| Nächste Marke | Wahrscheinlichkeit pro Marke |
|---|---:|
| b | 1/4 |
| a | 1/12 |
| Die dritte Marke der Geraden durch a,b | 1/6 |
| Jede der vier übrigen Marken | 1/8 |

Folglich findet man die dritte Marke anhand ihrer Häufigkeit, auch wenn keinerlei Geradentabelle zum Auswerten übergeben wurde. Der Abstand der beiden Alternativen unter den fünf neuen Namen beträgt 1/24. Alle 42 geordneten Paare wurden geprüft. Auf denselben sieben Namen gibt es 30 beschriftete Fano Ebenen. Alle 30 haben identische Einzel- und Paarstatistiken. Genau eine erklärt die vollständigen Dreierdaten. Dieser Vergleich ist modellintern. Die sieben bekannten Namen und ausreichend genaue Messwahrscheinlichkeiten werden vorausgesetzt.

Die lineare Dimension dieser Ausgabe ist exakt 21. Sei h_w=A_w 1. Für zwei Wörter u,v liefert das Beobachtungswort aus dem umgekehrten u und dem angehängten v die Gramform h_u^T h_v/21. Die symmetrischen Übergangsmatrizen machen diese Formel unmittelbar nachprüfbar. Die Spanne einschließlich Wortlänge null, eins und zwei hat die Ränge 1, 7 und 21. Bereits Wahrscheinlichkeiten bis Länge vier bezeugen den letzten Rang. Längere Wörter können ihn wegen der 21 Flaggen nicht erhöhen. Beim Test werden die A_u mit acht multipliziert, was die Spalten nur einzeln skaliert.

Jede stationäre Realisierung mit einem festen Instrument auf einer komplexen Matrixalgebra M_d und ohne weitere Erinnerung faktorisiert dieselben Wortwahrscheinlichkeiten durch einen Operatorraum der Dimension höchstens d². Damit ist d mindestens fünf. Die Realisierbarkeit mit d=5 wird dadurch offen gelassen. Die Zahl 21 ist ein exaktes Minimum für lineare Wortrealisierungen und für die hier gegebene klassische Markovrealisierung mit emissionsbeschrifteten Kanten. Der ursprüngliche 68-dimensionale L/J Prozess hat andere Eingriffe und andere Antwortdaten.

## Was aus dem tatsächlichen Vierertransfer folgt

Der ursprüngliche Septembervertrag setzt bedingt Φ=P₀+P₅/3+2P₁₀/3. Wir lesen die markierten Wörter aus dessen ursprünglicher Cliffordquelle ein und bilden die Karte in der Matrixeinheitenbasis. D ist die Projektion X↦tr(X)I/4. Schon der vollständige Kanal ohne Namen der fünf Achsen legt die lineare Abbildung S=3 id+4D−6Φ fest. Auf den drei Sektoren hat sie die Werte +1,+1,−1.

Eine lineare Rechnung für alle Matrixeinheiten E sucht ein B mit S(E)B=BEᵀ. Der Lösungsraum hat Dimension eins. Der im Originalkoordinatensystem gefundene reelle Erzeuger hat obere rechte Einheitsmatrix und untere linke negative Einheitsmatrix als Zweierblöcke, die übrigen Blöcke sind null. Damit ist S(X)=BXᵀB†, Bᵀ=−B und B·konj(B)=−I. Die umgekehrte Produktregel S(XY)=S(Y)S(X) ist für alle 256 Paare von Matrixeinheiten geprüft. Die quaternionische Struktur Θ=B·komplexe Konjugation ist also bis auf Phase aus diesem bedingten Kanal rekonstruierbar. Der physische Selbstverknüpfungsstrom und seine absolute Normierung ergeben sich aus dieser Operation allein weiterhin nicht. [Ursprünglicher Vertrag](../../../theory-contracts/marked-seam-cp-lift-20260921/SOURCE_PROVENANCE.md)

## Sieben mögliche Vierertransfers in einem gemeinsamen Messraum

Der gemittelte Familienprojektor des Vorgängerversuchs hat als Träger den sechsdimensionalen reellen Raum V der Nullsummenvektoren auf den sieben Marken. Für jede Marke u liegt die zugehörige dreidimensionale Differenzebene in der fünf-dimensionalen Tangentialebene senkrecht zur Simplexrichtung n_u=(e_u−1/7)/||e_u−1/7||. Ihr Orthogonalkomplement in dieser Tangentialebene hat Dimension zwei. Diese Tangentialebene ist ein neues Ausleseobjekt und darf nicht ohne weiteren Nachweis mit dem fünf-dimensionalen D₅ Gitterraum identifiziert werden.

Zusatzannahme ist ein gemeinsames chirales Spinregister für V und eine feste Identifikation seines Matrixraums mit dem nativen Laborregister. Die etablierte Isomorphie Spin(6)=SU(4) identifiziert die fünfzehn reellen spurlosen hermiteschen Observablen bis auf eine gemeinsame Normierung mit Λ²V. Eine ausgewählte Marke ergibt darin den Fünferraum n_u∧V. Unter ihrem Spin(5) Stabilisator spalten die übrigen zehn Bivektoren ab. Man kann einen benannten Fünferraum mit dem des zuvor ausgelesenen Originalkanals identifizieren. Die Ausrichtung des gesamten Siebenerrahmens gegenüber der ursprünglichen Quelle ist dabei eine Kalibrierung. [Darstellungstheoretische Grundlage](https://mathweb.ucsd.edu/~lni/math251B/LIEGROUPS-Donaldson-2007.pdf)

Sei P_u der Rangfünfprojektor auf n_u∧V. Auf dem spurlosen Sektor hat der einzelne Kanal die Matrix (2I−P_u)/3. Er ist in einem gemeinsamen Register eine unitäre Konjugation des nativen Kanals, also vollständig positiv und spurerhaltend. Jede zweifache Überlappung verschiedener P_u hat Spur 10/9. Die sieben Projektoren sind linear unabhängig, ihr Mittel ist I/3. Bei bekannten sieben Kanalmatrizen kann vollständige Eintakttomographie somit die sieben Mischgewichte bestimmen. Eine ausschließlich vergessene uniforme Marke ergibt den Skalar 5/9. Die tatsächliche, fest markierte Septembereintaktkarte besitzt dagegen die zwei verschiedenen Eigenwerte 1/3 und 2/3. Eine Identifikation der ursprünglichen Eintaktkarte mit dem uniform vergessenen Prozess scheitert exakt.

| Gesetz der Marken | Erster Takt, spurlos | Zweiter Takt, spurlos |
|---|---:|---:|
| Frisch und unabhängig | 5/9 | 25/81 |
| Jede Zielflagge wird weitergereicht | 5/9 | 101/324 |
| Anfangsmarke uniform, danach fest | 5/9 | 1/3 |

Im dritten Takt der weitergereichten Flagge spaltet der skalare Transfer auf. Die Eigenwerte sind 24427/139968 mit Vielfachheit sieben und 195479/1119744 mit Vielfachheit acht. Der Abstand ist genau 7/124416. Die Rechnung propagiert für jeden Zielzustand die volle geordnete Matrix auf Λ²V. Ersetzt man bei festgehaltenen Messachsen die Inzidenztabelle durch alle 30 Beschriftungen, entstehen **30 verschiedene dritte Transfermatrizen mit denselben Eigenwerten**. Die Lage der Spektralräume liefert die Information. Die Spurüberlappungen der Rang sieben Projektoren sind 7, 13/3, 3 und 7/3. Die kleinste quadrierte Hilbert-Schmidt-Distanz zweier dritter Transfermatrizen ist 49/2902376448. Reine Spektroskopie ohne Orientierungsinformation trennt keine der 30 Tabellen. Das Labor müsste die vollständige Prozessmatrix mit entsprechender Genauigkeit sehen. Der Originalprüfer enthält keine solche Serie siebenfach beschrifteter bedingter Matrizen.

## Von drei Ereignissen zu einer Windungsparität

Die im Originalregister bestimmten fünf hermiteschen Γ_j ergeben zusammen mit B sechs schiefsymmetrische Formen B₀=B und B_j=iΓ_j B. Sie erfüllen B_i B_j†+B_j B_i†=2δ_ij I. Die quadratische Pfaffsche Form von B(n)=Σ n_i B_i ist auf diesem reellen Raum genau Pf(B)·Σ n_i². Das Skript kontrolliert alle Diagonal- und Mischkoeffizienten, außerdem die exakte Überlappung 10/9 anhand zweier direkt gebildeter nativer 4×4 Matrizen. Die Einheitsnormalen geben also unitäre schiefsymmetrische Formen. Der Transfer kann das Vorzeichen von B(n) nicht sehen. Bei festem reellem Normalraum ist sein Parameter eine Linie in RP⁵.

Zwischen zwei verschiedenen Simplexlinien nehmen wir den kürzesten projektiven Großkreisbogen. Diese Wahl der Zwischenbewegung ist eine neue Hypothese. Die kanonisch zentrierten Einheitsvektoren haben Skalarprodukt −1/6. Zur stetigen Hebung des kürzesten Bogens auf die Sphäre muss jeweils die nächste Normale mit positivem Skalarprodukt gewählt werden. Jeder Schritt zwischen verschieden benannten Marken wechselt deshalb in der kanonischen Kalibrierung das Vorzeichen. Nach einem geschlossenen Dreieck a,b,c,a endet die gehobene Normale bei −n_a. Das Vorzeichenprodukt ist unabhängig von der Wahl der sieben Vertreter, da jeder lokale Vorzeichenwechsel zweimal auftritt. Der Beweis gilt für alle 35 Dreiecke und alle 64 verschiedenen Vorzeichengrammatrizen. Die Fano Inzidenz ist für diese Parität unerheblich. Bei weitergereichten Flaggen treten im stationären Vierwortprotokoll geschlossene Dreiecke auf einer Geraden mit Wahrscheinlichkeit 3/128 und die übrigen mit 7/128 auf. Unter geschlossenen Dreiecken liegen damit 3/10 auf einer Geraden. Bei unabhängigen uniformen Marken wären es 7/35=1/5.

Auf der reellen Hebung ist die Pfaffsche Form konstant, ihr Endpunkt als Matrix jedoch das Negative des Anfangs. Schließt man auch die Matrix durch eine stetige Phase exp(iφ), muss Δφ=(2k+1)π sein. Da Pf(exp(iφ)B(n))=exp(2iφ)Pf(B(n)), hat **jede solche geschlossene komplexe Hebung ungerade Pfaffsche Windung**. Die beiden Hebungen mit minimaler absoluter Phasenänderung haben Windungen +1 und −1. Die Tatsache, dass das beobachtete Dreieck die Parität erzwingt, ist stärker als das bloße Angeben eines Kreises mit Windung eins. Die Auswahl der physischen Phase, die Orientierung und die TFPT Interpretation als Selbstverknüpfungsstrom bleiben offen. Erst zusätzliches Dividieren durch vier komplexe Registerdimensionen und eine Periode 2π gibt den Betrag 1/(8π). Aus dem diskreten Protokoll allein folgt diese Normierung nicht.

## Der Familienprojektor als messbare Rückkehrwirkung

Eine stärkere Folgerung entsteht, wenn die projektiven Bögen auch als paralleler Transport der Registerachsen ausgeführt werden. Dies ist eine zusätzliche Prozessentscheidung. Auf jedem Bogen verwende man die kürzeste eigentliche Drehung in der von beiden Normalen aufgespannten Ebene, mit Identität auf dem orthogonalen Komplement. In den sieben Punktkoordinaten ist die vollständige Rückkehr um das Dreieck u,b,c,u genau das Produkt zweier orthogonaler Spiegelungen. Die erste wechselt das Vorzeichen der Normalen n_u, die zweite das der Differenz e_b−e_c. Der Code multipliziert die drei rationalen Rodriguesmatrizen. Für jeden der 15 Paartypen eines festen u wird die geschlossene Formel geprüft. Die Punkttransitivität überträgt das Resultat auf alle u.

Auf dem chiralen Spinregister ist diese Halbdrehung bis auf Phase eine Cliffordkonjugation mit der Achse d=(e_b−e_c)/√2. Das folgt, weil der Spinlift das Produkt der beiden Cliffordvektoren der Drehebene ist. Eine Rückkehr und ihre Umkehr haben dieselbe adjungierte Wirkung. Auf den fünf ursprünglichen Richtungen des lokalen Cliffordraums wirkt die beobachtbare Konjugation als 2P_d−I. Deshalb ist die Achsenmatrix P_d selbst aus einer linearen Antwort zugänglich, ohne eine Projektionsmessung auf diese Richtung vorauszusetzen.

Bedingt auf das erhaltene u und eine geschlossene Dreiecksgeschichte beträgt das Gewicht der drei kollinearen Paare zusammen 3/10, das der zwölf übrigen 7/10. Innerhalb jeder Bahn sind die Gewichte gleich. Die 15 Differenzrichtungen haben uniformes Mittel P_T/5 auf dem fünfdimensionalen Tangentialraum T. Die drei kollinearen Paare sind disjunkte Komplementpaare und ihre Achsenprojektoren summieren zum früheren Familienprojektor P_F. Daher folgt aus den gemessenen Dreieckshäufigkeiten allein

$$
M_u=\frac{7}{40}P_T+\frac{1}{24}P_F,
\qquad P_F=24\left(M_u-\frac{7}{40}P_T\right).
$$

Die Dreierfamilie ist der Eigenraum zum größeren Eigenwert 13/60; der restliche Zweierraum trägt 7/40. Bei uniformem Ziehen der fünfzehn Paare fällt diese Aufspaltung vollständig weg. Damit entfällt im angegebenen Prozess das zuvor separate Postulat, die Unterschiede komplementärer Punkte als Familienraum auszulesen. Seine Funktion übernehmen das Weiterreichen der Flagge und der minimale Drehtransport. Beide stehen weiterhin zur physischen Prüfung. Der Vergleich benutzt den tatsächlich früher konstruierten Projektor aus apartment.py in sämtlichen Einträgen. Im zugehörigen E₈ Gitter trägt dieser Raum zwölf A₃ Wurzeln und sein Orthogonalraum die vierzig D₅ Wurzeln. Die beiden Eigenräume des fünf-dimensionalen Tangentialraums sind nicht identisch mit der vollständigen Verteilung der acht Gitterdimensionen.

## Für den ursprünglichen Dreierleser reichen zwei Kontraste

Der Septemberprüfer liefert schon im ursprünglichen Dreierleser die Achsen F (Familiensumme) und A (Anker). Die originalen Pauliwörter erfüllen F=±Γ_i und A=±iFΓ_j, j verschieden von i. trine.py prüft das aus den Originalgeneratoren. In Bivektorkoordinaten vertreten sie bis auf Vorzeichen n₀∧x und x∧y mit orthonormalen n₀,x,y. Die drei ursprünglichen Trineeffekte spannen I,F,A. Deshalb sind die beiden normierten Rückkorrelationen FF und AA aus den präparierten Achsen oder aus der vollständigen ursprünglichen trinen Korrelationsmatrix rekonstruierbar.

Bei einem beliebigen Einheitsnormalen n seien a=n·n₀, b=n·x und c=n·y. Das Quadrat der Projektion eines einfachen Einheitsbivektors auf n∧V ist die quadrierte Länge seiner Kontraktion mit n. Daher liefert Φ_n die Rückantworten FF=(2−a²−b²)/3 und AA=(2−b²−c²)/3. In der Kombination Δ=(FF−1/3)+(2/3−AA) verschwindet b², und es bleibt Δ=(1−a²+c²)/3. Der Originaltransfer hat Δ=0. Das wählt die Linie ±n₀ sogar aus allen reellen Normalen aus, sobald die Spinorauslesung vorausgesetzt wird. Jede andere der sechs Simplexmarken trägt mindestens 35/108 zu Δ bei. Für ein Gemisch mit Gesamtgewicht w außerhalb der Referenz gilt Δ≥35w/108. Die Rechnung ist unabhängig von der tangentialen Kalibrierung der anderen Marken.

Damit erzwingt die tatsächliche feste Eintaktantwort innerhalb dieser Modellklasse die Beibehaltung der Referenzmarke im Laborregister. Eine gleichzeitige Behauptung, dieselbe feste Auslesung messe zu jedem Takt unverändert den Originaltransfer und ihr ungesehenes Wandern erzeuge die obige Dreieckswindung, ist in dieser Form widerlegt. Zwischen beidem vermitteln könnte ein registrierter oder mitbewegter Messrahmen. Welcher Rahmen im ursprünglichen Randfunktional physisch weitergegeben wird, ist offen. Die beiden Einzelantworten messen keine absolute Benennung der sieben Fano Punkte.

## Physische Zeit als Wartezeit zwischen Ereignissen

Für den markierten Übergang und den Spintransport zusammen verwende man als Observablenalgebra die direkten Summanden M₄ für jede der 21 Flaggen. Ein Ereignisschritt J zieht eine Gegenflagge mit Gewicht 1/8 und wendet den zu ihrer Marke gehörenden Kanal an. Der Übergang ist in beiden Richtungen gleich wahrscheinlich, die Marke ist bei Umkehr gleich, und jeder Kanal ist bezüglich der Spur selbstadjungiert. Daher ist J ein selbstadjungierter vollständig positiver, unitaler, spurerhaltender Kontraktionsoperator im normierten Hilbert-Schmidt-Raum der direkten Summe. Die uniforme Flaggenverteilung mit Spin I/4 ist sein eindeutiger stationärer Zustand. Auf Flaggenfunktionen gilt die Irreduzibilität von T, und auf den spurlosen Komponenten ist der Normfaktor höchstens 2/3, per Dreiecksungleichung und stochastischer Mittelung.

Nimmt man eine gedächtnislose kontinuierliche Wartezeit, zeitliche Homogenität und ausschließliche Sprünge zu Gegenflaggen an, erzwingt die Inzidenzsymmetrie dieselbe Rate auf allen 168 erlaubten gerichteten Übergängen. Das Skript berechnet den vollständigen Orbit. Nach Wahl der Zeiteinheit ist die Rate pro Ausgangsflagge eins. Ohne zusätzliche Zwischenentwicklung ist der gemeinsame kontinuierliche Transfer

$$
K(t)=\exp[t(J-I)]=e^{-t}\sum_{k\ge0}\frac{t^k}{k!}J^k .
$$

Die Summanden sind physische Kanäle. Der Generator H=I−J ist positiv. Die euklidische Reflexionsform K(s+t) ist deshalb für alle Zeiten positiv. An einer vorgegebenen Ereigniszahl wird weiterhin J hoch k gemessen. An einer vorgegebenen Dauer wird über diese Zahlen gemittelt. Die frühere Schranke für K(1)=Φ und K(2)=Φ² wird dadurch respektiert. Die stationäre Spur hat trivialen modularen Fluss. Diese Konstruktion begründet weder die absolute Sekunde noch die ursprüngliche Nahtzeit.

Der Vorzeichenleser nimmt +1 auf einer und −1 auf einer zweiten Marke an, sonst null. Mit v=(A₁−A₂)1 sind seine stationären Korrelationen bei k dazwischenliegenden unbeobachteten Sprüngen vᵀTᵏv/21. Für k=0,1,2 folgen 1/28, −1/168, 1/224. Ein positives euklidisches e^(−H) kann als feste Einfügung desselben Spaltes zwischen denselben Lesevektoren den negativen Wert nicht liefern. Für physisch kontinuierliche Zeit macht das Poissonmittel daraus dagegen die positive Korrelationsdichte

$$
e^{-t}\left[\frac{\cosh(t/\sqrt8)}{28}-\frac{\sqrt8}{168}\sinh(t/\sqrt8)\right].
$$

Die Gleichung benutzt T²v=v/8. Ihre logarithmische Anfangssteigung ist −7/6 in Einheiten der Gesamtsprungrate. Diese Vorhersage betrifft die Dichte zweier zeitlich getrennter markierter Sprünge ohne den Selbstkontakt bei t=0.

## Was die 64 Zyklen für P1 bedeuten

Der zusammenhängende Gegenflaggengraph hat 21 Knoten, 84 ungerichtete Kanten und ganzzahligen Zyklenrang 64. Die Inzidenzgruppe ist transitiv auf gerichteten Kanten einschließlich ihrer Umkehrung. Eine invariante reelle alternierende Kantenfunktion ist daher null. Jede invariante Kohomologieklasse hätte durch Mittelung einen invarianten Repräsentanten. Es gibt also keine von allen Inzidenzsymmetrien erhaltene nichttriviale reelle oder ganzzahlige Windung. Für ganze Koeffizienten nutzt man die Torsionsfreiheit der ersten Graphkohomologie. Auch bei Fixierung einer einzelnen Familienmarke ist die invariante Dimension null. Beim Signum des Markenstabilisators entstehen vier mögliche Richtungen. Bei Erhaltung eines ganzen ungerichteten Apartments gibt es sieben invariante Richtungen. Die Charakterrechnung steht in zyklen.py. Eine ausgewählte orientierte Linie aus diesen Alternativen müsste physisch begründet werden.

Der Rang 64 kann daher nicht mittels einer kanonischen invarianten Windung in eins plus 63 zerlegt werden. Eine Gleichsetzung mit den ursprünglichen 63 Gedächtnisrichtungen wäre unbelegt. Der geprüfte Atlasbefund dazu lautet finding:original-return-and-native-j-have-minimal68-dimensional-all-word-memory. Der benachbarte Vorversuch contract:joint-native-algebra-selection-20261016 behandelt bereits gemeinsame native Operatoralgebren. Hier werden eine andere Ereignisverkettung, eine direkt prüfbare Dreiwortstatistik und ihre Spinorauslesung untersucht.

## VERDICT und WILDCARD

Die geänderte erste freie Wahl ist das Weiterreichen der Zielflagge. Für den kontinuierlichen Prozess kommen gedächtnislose Wartezeiten und für die Quantenauslesung die Spinquantisierung sowie die Kalibrierung hinzu. Stufe 1 gilt für die expliziten endlichen Quanten- und Gedächtnisstrukturen. Die Originaldaten für eine auf dieselbe Marke konditionierte Dreiereignisantwort fehlen. P1, physische Identität der fünf Plätze, unbegrenzte Raumzeitlokalität, Gravitation und der vollständige Rückgewinn aller TFPT Ergebnisse sind offen.

Die ungewöhnlichste Möglichkeit ist eine tomographisch lesbare Familiengeometrie, obwohl die erste Messung auf jedem spurlosen Operator nur denselben Skalar zeigt. Für dieses Modell ist das eine exakte, widerlegbare Aussage. Bei vollständig kalibrierten Achsen müsste die dritte Prozessmatrix eine der 30 berechneten Orientierungen besitzen. Ohne solche Daten gibt es keine empirische Bestätigung.
