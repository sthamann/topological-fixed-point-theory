# TFPT. Welche Holonomie den Rückkehrkern tatsächlich erzeugt

12. November 2026. Exakte Abnahme der vorgeschlagenen Holonomie und des Anschlusses an R113.

Die Viertelholonomie ist auf dem ursprünglichen $E_8$ Gitter reproduzierbar. Ihre zwei gleichenergetischen Rückwege gehören zum zweidimensionalen Fixraum des $E_6$ Faktors. Die vom $A_2$ Faktor gedrehte Ebene ist ein anderer Raum. Auf ihr verschwindet derselbe Rückkehrkern. Der physikalische R113 Selektor muss daher mit dem richtigen Raum verglichen werden, bevor Energie oder Gewicht übertragen werden.

## 1. Vollständige Klassifikation für den gepinnten Coxeteroperator

Der vorhandene ganzzahlige Seifert Intertwiner trägt die globale Monodromie in die Weylgruppe des ursprünglichen $E_8$ Gitters. Für den gepinnten Coxeteroperator $C$ habe ich alle 2.240 verschiedenen $A_2$ Coxeterfaktoren $B=s_rs_s$ aus Wurzelpaaren mit $(r,s)=-1$ gebildet. Aus $C=AB$ folgt jeweils eindeutig $A=CB^{-1}$. Die ersten acht Potenzspuren prüfen das ganze charakteristische Polynom. Die 20 Treffer haben den $E_6$ Coxetercharakter. Ihr bewegter Rang sechs trägt jeweils 72 Wurzeln, der feste Rang zwei sechs Wurzeln. Die 72 bewegten Wurzeln bilden einen $E_6$ Wurzelraum. Der punktweise Stabilisator des Fixraums ist seine Weylgruppe. Keine dieser Wurzeln hat unter $A$ eine Bahn der Länge drei. In der rationalen Spektralzerlegung $\Phi_{12}\Phi_3$ hat daher jede Wurzel eine primitive Komponente. Das liefert einen regulären Eigenvektor zur primitiven zwölften Einheitswurzel und damit die $E_6$ Coxeterklasse. [Brillon, Ramazashvili, Schechtman und Varchenko zur Seifert und Coxeter Korrespondenz](https://armj.math.stonybrook.edu/pdf-Springer-final/017-0065.pdf) · [Springers Satz über reguläre Elemente](https://mathweb.ucsd.edu/~nwallach/Springer-reg.pdf).

Für $W=ABA^{-1}B^{-1}$ ergibt sich

| Anzahl | Bahnen unter Konjugation durch $C$ | Polynom von $W$ | Rückkehr auf $\operatorname{Fix}(A)$ |
| ---: | ---: | --- | --- |
| 5 | eine Bahn der Länge 5 | $(x-1)^4(x^2+1)^2$ | Rang zwei, Energie 1 |
| 15 | eine Bahn der Länge 15 | $(x-1)^6(x^2+x+1)$ | Rang eins, Energie $5/3$ |

Das ist eine Klassifikation der geordneten Weylfaktoren für den festgehaltenen Coxeteroperator. Sie wählt keine physikalische P1 Deformation oder ausgezeichnete Vanishing Paths. Zwei Coxeteroperatoren auf orthogonalen $E_6$ und $A_2$ Untergittern würden kommutieren. Deren Kommutator wäre die Identität. Im untersuchten Produkt ist die gedrehte $A_2$ Ebene gegenüber dem festen $A_2$ des $E_6$ Faktors schräg. Im Viertelfall gilt exakt $P P_{A_2,\mathrm{Pfad}}P=P/3$. Das unterscheidet zwei verschiedene $A_2$ Rollen auf demselben Gitter.

## 2. Der genaue endliche Kandidat und seine Auswahl

Setze $H=2I-W-W^{-1}$, $P$ als orthogonalen Projektor auf $\operatorname{Fix}(A)$ und $Q=I-P$. Für alle fünf Viertelfälle ergibt die Originalrechnung

$$
PHP=P,\qquad
B_{\mathrm R}^\dagger B_{\mathrm R}=P,\qquad
(QHQ)B_{\mathrm R}=B_{\mathrm R},\qquad B_{\mathrm R}=QHP.
$$

Die mittlere Identität liefert zwei Zwischenrichtungen. Die letzte legt deren gemeinsame Energie fest. In der durch $B_{\mathrm R}$ transportierten Basis ist der aktive Generator die Blockmatrix mit vier Identitätsblöcken. Der projizierte Transfer ist $(1+e^{-2t})P/2$ und der Gedächtniskern $e^{-t}P$. Auf der vom $A_2$ Pfad gedrehten Ebene gilt dagegen $PHP=2P$ und $B_{\mathrm R}=0$. Die beiden Räume sind deshalb durch dieselbe Rückkehrprobe eindeutig unterscheidbar.

Bei den fünf Viertelholonomien hat jedes $H$ eine Bahn der Länge sechs unter der ursprünglichen Zwölferclock $c=J\sigma$. Genau eine der fünf Orientierungen hat im gepinnten gemeinsamen Rahmen das Clockmittel $I_8$. In den verwendeten $D_8^+$ Koordinaten lauten ihre $A_2$ Wurzeln $-e_2-e_6$ und $e_6-e_7$. Das rekonstruiert den im eingefügten Text angegebenen Fall aus den Originalmatrizen. Rang zwei sondert die Fünferbahn aus, wenn der Selektor der $E_6$ Fixraum ist. Das isotrope Clockmittel sondert einen Vertreter aus, falls diese Eigenschaft von der physischen Quelle verlangt wird.

Die direkte Summe über die sechs Orientierungen erhält die Zwölferclock. Der blockweise Rahmenwechsel $\operatorname{diag}(I,c,\ldots,c^5)$ zieht den Hamiltonoperator auf sechs gleiche Kopien zurück. Die Clock wird zur zyklischen Verschiebung mit Minuszeichen am letzten Übergang. Damit bleibt die Antiperiodizität erhalten. Der Transport erzeugt allein keine Kopplung zwischen den sechs Hamiltonblöcken. Genuin unterscheidbar sind relative Einfügungen, für die nicht alle Felder mit demselben Rahmen konjugiert werden. Ein streng invarianter Vakuumvektor ist in dieser ungeraden Registerdarstellung wegen $\mathbf C^6=-I$ ausgeschlossen.

## 3. Der kleinste Originaltest für R113

Der ursprüngliche R113 Vertrag enthält bereits den ganzen komplementären Rückkehrterm. Für ein identifiziertes gemeinsames selbstadjungiertes $H$ und festes $P$ definiere $A_0=PHP$, $C_Q=QHQ$, $G=B_{\mathrm R}^\dagger B_{\mathrm R}$, $K=B_{\mathrm R}^\dagger C_Q B_{\mathrm R}$ und $L=B_{\mathrm R}^\dagger C_Q^2 B_{\mathrm R}$. Mit $\Pi$ dem Projektor auf $\operatorname{Ran}B_{\mathrm R}$ gilt am endlichen Regulator

$$
D=L-KG^{-1}K=(C_QB_{\mathrm R})^\dagger(1-\Pi)(C_QB_{\mathrm R}).
$$

Der Ausdruck gilt auf dem Träger von $G$, bei fehlender Invertierbarkeit mit seiner Pseudoinversen. Ein positiver Eigenwert von $D$ zeigt einen zusätzlichen Rückweg. $D=0$ schließt den erreichten Zwischenraum. Eine einzige aktive Energie folgt erst zusätzlich aus $K=\lambda G$. Tatsächlich ist $D=0$ für jeden Projektor automatisch, wenn $H$ ein quadratisches Polynom mit konstanten Koeffizienten erfüllt. Die Projektion der quadratischen Gleichung auf den $QP$ Block gibt $C_Q B_{\mathrm R}=B_{\mathrm R}(\nu I-A_0)$, falls $H^2=\nu H+\kappa I$. Daher bestehen sowohl der goldene Kandidat als auch beide Holonomieklassen diese Nullprobe. Bei $D=0$ und $K=\lambda G$ unterscheidet der direkte Block $A_0=\lambda I$ die Viertelholonomie vom goldenen Kandidaten mit $A_0=2\lambda I$ in der festgelegten Nullenergiekonvention. Die Differenz $A_0-\lambda I$ ist auch bei einem gemeinsamen Energieoffset invariant. Die Viertelholonomie verlangt zusätzlich $G=\lambda^2 I$ bei gemeinsamer Zeitskalierung.

Wenn der Originalaufruf nur den komprimierten Transfer $S(t)=Pe^{-tH}P$ ausgibt, genügen vier Ableitungen am Ursprung. Auf dem in der ursprünglichen Nullzeit Gramform normierten $P$ Raum sei $M_j=(-1)^jS^{(j)}(0)$. Dann gilt ohne Kommutationsannahme

$$
\begin{aligned}
 A_0&=M_1,\\
 G&=M_2-A_0^2,\\
 K&=M_3-A_0^3-A_0G-GA_0,\\
 L&=M_4-A_0^4-A_0^2G-A_0GA_0-GA_0^2-G^2-A_0K-KA_0 .
\end{aligned}
$$

Das folgt aus den Blockprodukten von $H$ und wurde mit nichtkommutierenden symbolischen Blöcken exakt ausgeführt. Für den unbeschränkten physikalischen Generator sind der tatsächliche Zeitparameter und die nötigen gemeinsamen Ableitungsdomänen Teil des Tests. Eine normierte Randantwort mit frei zeitabhängigem Nenner darf nicht als $S(t)$ eingesetzt werden. [Feshbach und Schur Methode](https://arxiv.org/abs/2105.02058).

Falls eine Originalrechnung die drei Momente als positive Mischung gemeinsamer markierter Kanäle liefert, gibt es einen zusätzlichen exakten Test. Setze $A_b=G_b^{-1}K_b$ und $A=G^{-1}K$ für die gemischten Momente. Dann ist $D$ die Summe des gemittelten Einzelfaserdefekts und des positiven Terms $\int (A_b-A)^\dagger G_b(A_b-A)\,d\mu(b)$. Selbst wenn jeder Einzelfaserdefekt verschwindet, kann die gemeinsame Energieordnung weitere Rückwege erzeugen. Diese Formel setzt die positive Mischung auf derselben Endbasis voraus. Das gekoppelte operatorwertige R113 Modell wurde hier nicht durch eine solche Mischung ersetzt.

## 4. Physischer Stand

Die beiden historischen lokalen Oszillatorpartitionen $(2)$ und $(1,1)$ sind belegte Beiträge der geordneten Kreisquelle. Ihre gleiche Energie und eine mögliche Gramnormierung legen weder den ursprünglichen Selektor noch die physisch gewichtete gemeinsame Zeitantwort fest. Die Originaldatei TFPT-FAMILIENSTROM-PIN-GANZWORT-UND-TEILWOERTER-20261110.txt hält auch fest, dass die Summe bei anderen Ordnungen auslöschende Beiträge enthält und ein mittlerer Ramond Filter die betreffenden Moden löscht.

Der unveränderte lokale Aufruf wurde erneut ausgeführt. Der direkte Weg liefert original_full_source_operation_not_supplied beim Eintrag Q_0[a,b]=I_orig((Theta F_a) theta_T F_b). Der gestufte Weg endet mit Exit 2, vor der Integration, bei null Originalatomen. Physikalische $G,K,L,D$ wurden deshalb nicht ausgewertet. Die vorhandene tangentiale Einteilchenmatrix und die gegebenen konformen Zwischenworte wurden nicht als Ersatz für diesen Aufruf verwendet.

Der aktuelle Originalbefund TFPT-R113-GELADENE-PAAR-RUECKKEHR-20261112 bestimmt bereits den passenden Feldaufruf. Der Mittelstrom $J^\gamma$ ist das gerade D8 Fermionpaar $(1,6)$. In der konformen Längszeit gilt $f e^{-tL_0} e\Omega=e^{-t}\Omega$. Das ist ein Rückweg auf der einzigen Vakuumlinie. Es ist keine positive Rang-zwei-Matrix auf dem Holonomieprojektor. Ein direktes Spinorpaar im nackten R113 Fockraum verfehlt sogar das neutrale Wiederholungswort mit $7/4$ anstelle von null. Für die tatsächliche gewichtete Paarantwort muss der passende Minor von $\wedge^2 R_t[b]$ innerhalb derselben bosonischen Integration und vor der äußeren Auswahl mitgeführt werden. Der Leerfermionblock allein liefert diesen Eintrag nicht.

Bei der in der letzten Cartanherleitung geprüften Identifikation $J=D|_{\mathfrak h'}$ entspricht die sechste Clockpotenz der ursprünglichen Spinorparität $D^2$. Das Register hat diese Parität $-1$ auf jedem Vektor. Der ursprüngliche D8 Paarzustand und sein neutraler Vakuumrückweg haben Parität $+1$. Jede lineare Abbildung des Registers auf diese Zustände, die dieselbe Halbphase erhält, verschwindet. Diese genaue Identifikation braucht einen geraden Darstellungsanteil oder eine aus der Quelle hergeleitete paarige Auslesung. Eine Umbenennung der beiden Zwischenzustände genügt nicht.

Das lokale $E_8$ Netz und die frühere geladene Cartanverkettung bleiben die Ladungs und Stressprobe. Eine holonomische Weylwirkung auf dem internen Zyklusraum trägt für sich keine zusätzlichen Regionalalgebren. Der Faktor $4/9$, sein neunfacher Antwortsektor und die beiden zusätzlichen lokalen Raumrichtungen bleiben herzuleiten. Die neue Rechnung bestimmt erstmals, welchen endlichen Selektor und welche vier Originalableitungen der Holonomievorschlag tatsächlich verlangt.

## Reproduktion

- check_tfpt_holonomy_r113_20261112.py mit dem TFPT Repository als erstem Argument. Die native Coxetermatrix und Zwölferclock werden über den bestehenden Cartanchecker aus dem originalen Triality Zertifikat gelesen.
- R113 Original: experiments/theory-contracts/original-local-source-variable-background-20261102/PROJECTED-JOINT-SOURCE-GENERATOR-20261103.txt, Abschnitt 4.
- Nativer Wiederholungslauf: docs/research/fundamental-plan/p1-original-marken-und-metrik-20261112/replay_tfpt_native_p1.py.
