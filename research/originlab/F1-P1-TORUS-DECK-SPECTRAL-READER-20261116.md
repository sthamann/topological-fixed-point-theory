# F1 · Torus Double Cover als spektraler P1 Lesekopf: neuer Kandidat und harte Gegenproben

**16. November 2026. Unabhängige mathematische Ursprungshypothese. Kein physischer F1 Beweis.** Die Originalquelle enthält weiterhin `whole_source_bundle=null`; die alte Fourierkreisvariable ist ausdrücklich räumlich und nicht OS Zeit. Diese Notiz ersetzt weder `Q_0` noch `P_adm` noch die absolute R113 Leerreferenz. Sie setzt an einer anderen Stelle an als die schon dokumentierten elliptischen Determinanten und die Klein Viersektoralgebra.

## 1. Nicht P1 als Punkt, sondern P1 als Quotientenoperation

Beginne **bedingt** mit einem kompakten flachen zweidimensionalen Torus
`T_y=R²/Lambda_y`, Rechteckaspekt `y>0` und Gesamtfläche `A=1`.
Der *ungerichtete* geometrische Deckoperator ist

`Jf(z)=f(-z)`, `J²=I`, `[J,Delta_y]=0`.

Seine vier Fixpunkte `2z=0 mod Lambda_y` sind exakt
`{0,1/2,tau/2,(1+tau)/2}`, die natürliche Torus Zweitorsion `(Z2)^2`.
Der Quotient `T_y/{z~ -z}` ist topologisch die Riemannsche Sphäre
`P¹` mit **vier Verzweigungswerten**, durch Riemann Hurwitz
`2g(T)-2=-4+4=0`. Die induzierte Metrik ist die *konische*
Pillowcase Metrik mit vier Kegelwinkeln `pi`, NICHT die glatte runde
Standardkugel. Der Divisor und die geometrische Metrik dürfen nicht
verwechselt werden.

Ein möglicher TFPT P1 Lesekopf wäre die **Operation** aus
Deckprojektion `P_+=(I+J)/2`, Quotientenidentifikation und anschließend
aus der Quelle gewähltem Quanteninstrument. P1 müsste in diesem Szenario
weder eine ursprüngliche Singularität noch ein stationärer Fixpunkt sein.
Ob die physischen TFPT P1 Operatoren diese Bedeutung haben, ist offen.

## 2. Der Faktor 1/(8pi) aus demselben spektralen Doppeldecker

Für den **skalaren** Laplaceoperator auf dem flachen Einheitstorus
mit Fouriermoden `(m,n) in Z²` gilt

`lambda_(m,n)=4*pi²*(y*m²+n²/y)`,
`Z_y(t)=Tr(exp(-t Delta_y))=sum_(m,n)exp(-t*lambda_(m,n))`.

Die Inversion vertauscht `(m,n)` mit `(-m,-n)`. Nur `(0,0)`
bleibt fest. Daher sind die beiden spektralen Paritätsspuren für
**JEDES t>0, JEDEN Torusaspekt und JEDE symmetrische Fourierabschneidung** exakt

`Z_+(t)=(Z_y(t)+1)/2`,
`Z_-(t)=(Z_y(t)-1)/2`,
`Z_+(t)-Z_-(t)=1`.

Poissons Summenformel liefert
`Z_y(t)=1/(4*pi*t)*sum_(p,q in Z²) exp(-(p²/y+y*q²)/(4t))`.
Somit ist die dimensionslose normierte führende Wärmeresidue exakt

`lim_(t->0+) t*(Z_+(t)-1/2)=1/(8*pi)`.

Es ist ein echter **mathematischer** Herkunftsweg für dieselbe Zahl
aus einer Projektionsoperation auf einer Fläche mit vier Zweigpunkten.
Er ist KEIN Nachweis, dass der ursprüngliche TFPT Koeffizient
`c3` physikalisch genau diese skalare spektrale Residue ist.
Die Wahl effektiver spektraler Dimension zwei, der Fläche eins,
des physikalischen Operators und des Referenzmaßes sind Zusatzannahmen.
Bei beliebiger Fläche A ist der Rohkoeffizient `A/(8*pi)`.
Die Halbierung des Weyl Koeffizienten wäre auch bei anderen
zweidimensionalen Involutionen möglich; der bloße Wert ist
**kein** eindeutiger Torus oder Viermarkenselektor.

Eine **zusätzliche** scharfe Signatur ist die *unveränderte
Paritätsspurdifferenz* `Z_+-Z_-=1` über alle spektralen Skalen.
Bei hyperelliptischen Doppeldeckern der Sphäre mit `r=2g+2`
einfachen Verzweigungspunkten gibt Donnellys Fixpunktformel
`lim_(t->0) Tr(J exp(-t Delta))=r/4=(g+1)/2`, während für die
skalare konstante Nullmode `lim_(t->infty)=1` gilt. Die Forderung
**skalenkonstanter Spur 1** erzwingt damit `g=1,r=4` innerhalb
dieser Klassen als **notwendige**, nicht hinreichende Bedingung.
Sie ist selbst ein noch unbegründetes physikalisches Prinzip.

Numerische Kontrollen, direkte Fourierreihe mit `|m|,|n|<=40`:

| Aspekt y | t | `t*(Z_+-1/2)` | `Z_+-Z_-` |
|---|---:|---:|---:|
| 1 | 0.005 | 0.03978873577297383 | 1 |
| 1.7 | 0.005 | 0.03978873577298724 | 1 |

Der Grenzwert ist `1/(8*pi)=0.03978873577297383...`.
Die angegebenen Fourierkontrollen gelten bei den genannten endlichen t.
Bei festem Cutoff und anschließend t gegen null würde die numerische
Residue dagegen gegen null gehen. Die Residuenfunktion verwendet deshalb
eine stabile Poisson/Jacobi Auswertung; die volle spektrale Grenzaussage
wird nicht aus einer fest abgeschnittenen Summe gefolgert.
**Die Quadratlage y=1 wird durch diese Quelle NICHT ausgewählt.**
Das frühere elliptische Determinantenmodell kann als separater
Modulus Selektor geprüft werden, darf aber nicht still in das
Spektralmodell hineindefiniert werden.

## 3. Die zwingende Klein Gegenprobe

Auf derselben rechteckigen Torusmetrik kann stattdessen die
**freie orientierungsumkehrende** Involution
`K(x,u)=(x+1/2,-u)` wirken. Der Quotient ist die Kleinsche
Flasche, nicht die verzweigte Sphäre. Die genaue
equivariante Wärmespur ist

`Tr(K exp(-t Delta))=sum_(m in Z)(-1)^m exp(-4*pi²*t*y*m²)`
`=1/sqrt(4*pi*t*y)*sum_(k in Z)exp(-(k+1/2)²/(4*t*y))`.

Sie geht für `t->0` exponentiell gegen **null**, während die
Torusinversion **exakt 1** hat. Für `y=1,t=0.005`
ergibt sich `2.9734390e-5`, für `y=1.7,t=0.005`
etwa `0.00392079`. Beide freien und verzweigten Versionen
haben den führenden halben zweidimensionalen Weyl Koeffizienten,
aber NUR der verzweigte Torus dieses Modells trägt die
skalenunabhängige unpaarige scalar Nullmode Signatur.

Das ist der präzisere Unterschied für eine spätere
Original `Q0` Prüfung, nicht die bloße Übereinstimmung mit `1/(8pi)`.

## 4. Die Spin Doppelüberlagerung ist eine zweite, andere Überlagerung

Die Deckinversion `z->-z` dreht den orientierten Tangentialraum
um `pi`. Auf einem chiralen Spin(2) Faserraum wird diese
Drehung als `U_spin=exp(i*pi/2)=+i` beziehungsweise mit
entgegengesetztem Lift als `-i` dargestellt. Der
globale spinorische Operator ist schematisch `U=iJ` und erfüllt

`U²=-I`, `U⁴=I`.

Daher können **Viertelphasen** aus einer geometrischen
Halbdrehung durch das Spin Double Cover kommen, ohne eine
vorgegebene `Z4` Pfadphase zu behaupten. Die Existenz der
mathematischen Spinhebung ist jedoch KEINE Auswahl ihres
Vorzeichens oder des physischen TFPT q0. Die skalare Projektion
`(I+J)/2` ist nicht automatisch dieselbe Operation wie das
spinorische / fermionische Quellenintegral.

Die vier flachen Spinstrukturen auf dem Torus haben
Fourier Halbcharakteristiken `(alpha,beta) in {0,1/2}²`.
Nur `(0,0)` besitzt eine unpaarige Fourier Nullmode und
`Tr(J e^(-tH_spin))=1` in der elementaren spektralen
Inversion. Die drei übrigen besitzen dort `0`. Dieser
Unterschied kann den Spinsektor diagnostizieren; eine
ungerade Spinstruktur kann aber zugleich ungesättigte
fermionische Nullmoden und einen verschwindenden Pfaffian
in der **unmarkierten** Quelle erzeugen. Originales `Z0`,
zulässige Nullmodeneinfügungen und die absolute Referenz
entscheiden, ob das akzeptabel ist.

## 5. Zwei übereinandergelegte Tori sind meist KEIN Torus mehr

Zwei unabhängige einfach verzweigte **elliptische** Doppelüberlagerungen
derselben `P¹` mit jeweils vier Zweigwerten `B1,B2` ergeben
für `k=|B1 intersect B2|<4` nach **Normalisierung des Faserprodukts**
einen zusammenhängenden Vierfachdeck mit Gruppe `(Z2)²`. Bei gemeinsamen
Verzweigungswerten ist das rohe Faserprodukt nodal; die folgende Formel
nennt den geometrischen Genus der Normalisierung. Riemann Hurwitz ergibt

`2g(C)-2=4*(-2)+2*|B1 union B2|`, folglich `g(C)=5-k`.

| Gemeinsame Zweigwerte k | Genus des verbundenen Faserprodukts |
|---:|---:|
| 0 | 5 |
| 1 | 4 |
| 2 | 3 |
| 3 | 2 |
| 4 | Normalisierung: zwei getrennte Kopien desselben Torus |

Eine **unverzweigte** doppelte Überlagerung eines Torus hat
dagegen wieder Genus 1. Die orientierte Doppelüberlagerung einer
Kleinschen Flasche ist ein Torus ohne Zweigpunkte. Die Wörter
„doppelter Torus“ bezeichnen also drei mathematisch verschiedene
Ursprungsszenarien. Ihre Vorhersagen für Verzweigung und Spur
unterscheiden sich messbar. Keine Wahl davon ist schon
physikalisch erzwungen.

## 6. P1 als holografischer Lesekopf oder als dynamischer Transformationspunkt

Auf Dichtematrizen ist die legitime **unbedingte** Deckauslese
`E(rho)=(rho+J rho J)/2`: ein CP, spurtreuer,
idempotenter Kanal, der Paritätskohärenz verwirft. Er ist
auf dem Gesamtzustandsraum **nicht injektiv**:
die beiden orthogonalen Zustände
`|+>=(|0>+|1>)/sqrt(2)` und
`|->=(|0>-|1>)/sqrt(2)` liefern unter
`J=diag(1,-1)` beide `I/2`.
**Bloßes Projizieren ist daher keine holografische Vollrekonstruktion.**
Dafür braucht die ORIGINALQUELLE eine identifizierte
Codeunteralgebra, ein konsistentes bulk-to-boundary Instrument
und einen überprüfbaren Rekonstruktionskanal. Eine reversible
Stinespring Dilatation mit gespeichertem Ausleseergebnis wäre
ein eigener, physikalisch zu prüfender Quellentyp.

Eine ehrlich dynamische vollständig positive Alternative ist die
aus derselben Deckaktion gebildete GKLS Form

`L_gamma(rho)=gamma/2*(J rho J-rho)`,
`T_t(rho)=E(rho)+exp(-gamma*t)*(rho-E(rho))`, `gamma>=0`.

Dies ist eine CPTP Halbgruppe. Die ideale P1 Auslese `E`
wird als `t->infty` Grenze erreicht, nicht als endlicher
Zeitpunkt einer beschränkten invertierbaren Halbgruppe.
Ein konkreter reiner Pluszustand wird auf eine positive Matrix
mit Eigenwerten `(1±exp(-gamma*t))/2` abgebildet. Das ist
eine mathematische **Dynamik** aus der Projektionsalgebra,
aber `gamma`, physische OS Zeit, Raumlokalität und
Rückkopplung an Materie sind NICHT hergeleitet.
Die dazu parallele `exp(-t Delta)` Diffusion besitzt
ein natürliches positives mathematisches `Delta`, ist
aber keine automatisch aus TFPT erzeugte Lorentzzeit.

Mit `U=iJ` ist `Ad_U=Ad_J`: diese **inkohärente**
Deckhalbgruppe ist blind für das Vorzeichen der Viertelphase.
Die originalen geordneten R113 Phasenantworten verlangen
einen *kohärenten* gemischten Quellenmarker, nicht bloß
Paritätstwirling.

## 7. Was, wenn die Struktur VOR P1 größer als E8 ist?

Eine stimmige konkurrierende Möglichkeit ist ein
`E8_1 tensor E8_1` Parent (`c=16`) oder ein
anderes rangsechzehn even unimodulares Gitter wie `D16+`.
Beim symmetrischen Zusammenfassen zweier *unabhängiger*
affiner Level eins Ströme gilt für die Diagonalströme
`J=J1+J2` zwingend **Level zwei**;
bei drei Kopien Level drei. Eine vom Parent
durch die Austauschgruppe erzwungene diagonale Auswahl
liefert deshalb NICHT automatisch die native
`E8_1` Stromalgebra. Level eins verlangt zusätzliche
Selektion, Projektion oder ein anderes explizites Interface.
`E8²` und `D16+` haben zudem beide 480 Wurzeln plus
16 Cartanrichtungen, aber unterschiedliche Lie Algebren:
bloße `496` Zählung kann die Parent Struktur nicht auswählen.
Ein `c=24` Parent oder Niemeier/Leech Szenario ist eine
weitere Klasse und muss dieselben Level, Anomalie,
Nullmoden und Einquellenbedingungen erfüllen.

## 8. Die notwendige Originalquellenprobe und Stoppregeln

Originale Eigentümer:
- [Native Originaldomain](../../experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/domain/native-diagnostic-v1.json): `whole_source_bundle=null`; kreisförmige Koordinate ist räumlich, NICHT OS Zeit.
- [R113 Original und Referenz](../../experiments/theory-contracts/chat-triangulation-20261108/R113-ORIGINAL-UND-REFERENZ-20261109.txt): `Q_j[b]=s_j[b] Gamma(R_j[b])`, gesamte äußere `P_adm` beidseitig, absolute leere Referenz bleibt erhalten.
- [Direkter gemeinsamer Quellenkern](F1-DIRECT-SINGLE-ORIGINAL-SOURCE-20261116.md).
- [Elliptischer Determinant als separater Quadrat Selektor](P1-ELLIPTIC-DETERMINANT-SOURCE-20261115.md).
- [Klein Viersektoralgebra und unterbestimmter Spin](F1-P1-KLEIN-INFINITE-TWISTED-ORIGIN-20261116.md).
- [F1 Quantenmarginal No Go](F1-P1-PROJECTIVE-CONSISTENCY-AND-COHERENT-GLUING-20261116.md).
- [E8 chiral anomaly und minimaler Doppelcover](F1-E8-ANOMALY-INFLOW-ORIGIN-20261116.md).

**Akzeptanzbedingung eines echten Erfolgs:** Eine EINZIGE ursprüngliche
gewichtete, OS/CP positive Quelle mit tatsächlichem `P_adm`
muss selbst den physikalischen Deckoperator, Spinstruktur,
zulässige Geschichte, vollständige `Q0/Z0, ZU, ZY, ZUY`
und beide realen Blockstufen liefern. Aus ihr müssen die
physikalische gleiche GNS/Energieform mit Higgs, Familien,
Clock und Stress, 3+1D chirale Lokalität, das korrekte
native E8 Niveau und universeller masseloser Spin 2 folgen.
Der bloße exakte Topologie, Spektral und Channel Test
weder ersetzt noch erfüllt diese Anforderungen.

**Konkrete Ausschlüsse:** (i) Der halbe Weyl Koeffizient
allein ist nicht toruseindeutig. (ii) Der Grenzwert ohne
absolute Originalreferenz ist nicht automatisch `c3`.
(iii) Torusinversion wählt das Quadrat nicht. (iv) Eine
Paritätsauslese ist kein verlustloses allgemeines
Hologramm. (v) Ihre dephasierende Uhr bestimmt keine
physische Zeitskala. (vi) Zwei gleichberechtigte E8
Parentkopien liefern diagonal Level 2 statt Level 1.

Reproduzierbare Standardbibliothek:
`p1_torus_deck_spectral_reader.py`,
`test_p1_torus_deck_spectral_reader.py`.
Lokaler kontrollierter Lauf am 16. November:
**14/14 Tests bestanden**, einschließlich zweier zusätzlicher
Kontrollen des kleinen Zeitgrenzwerts. Fourier und Poissonform
unabhängig gegengeprüft. Status F1 ausdrücklich **OFFEN**.

Externes Fundament:
- H. Donnelly, asymptotische Fixpunktformel für equivariante Wärmeoperatoren: https://gdz.sub.uni-goettingen.de/download/pdf/PPN235181684_0224/PPN235181684_0224.pdf
- Spin(2) Doppelcover (RMP 88, 035001): https://link.aps.org/accepted/10.1103/RevModPhys.88.035001
- Elliptischer Vierzweigdecker und Riemann Hurwitz: https://www.mdpi.com/2227-7390/14/8/1287
- Holografische Quantenrekonstruktion und Recovery: https://arxiv.org/abs/1704.05839
