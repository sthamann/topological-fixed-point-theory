# Fünf binäre Messungen, ein diskreter Transfer und sein Pfaffsches Signal

**Ergebnis.** Auf dem vierdimensionalen Spinorregister reproduziert ein konkretes Messprotokoll den bereits im TFPT Compiler bedingt konstruierten Transfer in allen 256 Matrixeinträgen. Derselbe Satz binärer Projektoren liefert einen Pfaffschen Kreis mit Windung eins. Die Randnormierung der TFPT verlangt darüber hinaus eine Identifikation dieses Stroms mit der Selbstverknotung der ursprünglichen Naht. Der bedingte Zweizeitversuch im letzten Abschnitt prüft eine weitere unabhängige Lücke, die Wahl des physischen Messinstruments.

## Ein Schritt auf der vollen Operatoralgebra

Die fünf hermiteschen Matrizen erfüllen
\[
\{\Gamma_i,\Gamma_j\}=2\delta_{ij}I_4.
\]
In der hier verglichenen Darstellung sind sie \(I\otimes X,I\otimes Z,X\otimes Y,Z\otimes Y,Y\otimes Y\). Das Ergebnis hängt von den Cliffordrelationen ab und bleibt bei Umnummerierung, Vorzeichenwechsel und simultärem Basiswechsel erhalten.

Ein Schritt verwendet die Identität und fünf binäre Messungen mit den Projektoren \(P_{j,\pm}=(I\pm\Gamma_j)/2\). Für einen unbeobachteten Ausgang ergibt die Lüdersregel
\[
L_j(A)=P_{j,+}AP_{j,+}+P_{j,-}AP_{j,-}
      =\tfrac12(A+\Gamma_j A\Gamma_j).
\]
Sind die sechs Möglichkeiten gleich wahrscheinlich, folgt
\[
\Phi=\frac{\mathrm{id}+\sum_jL_j}{6}
     =\frac7{12}\mathrm{id}+\frac1{12}\sum_j\operatorname{Ad}_{\Gamma_j}.
\]
Die Wahl der Möglichkeiten, ihrer Gewichte und der Zustandsaktualisierung gehört zum Protokoll. Die Gleichheit der Endkarten beweist für sich allein keine dieser Wahlen.

| Operatorraum | Dimension | Fünf Messungen ohne Leerschritt | Sechs Möglichkeiten |
| --- | ---: | ---: | ---: |
| \(I\) | 1 | 1 | 1 |
| \(\Gamma_i\) | 5 | \(1/5\) | \(1/3\) |
| \(i\Gamma_i\Gamma_j,\ i<j\) | 10 | \(3/5\) | \(2/3\) |

Für einen festen Einzelgenerator kommutiert genau einer der fünf Generatoren mit ihm. Für ein Zweierprodukt antikommutieren genau die beiden beteiligten Generatoren mit dem Produkt. Damit erhält man die Tabelle auf einer vollständigen orthogonalen Basis. Der direkte Vergleich der beiden Superoperatoren ergibt 256 von 256 identische Einträge. Innerhalb der fünf Messungen ohne Leerschritt maximiert die Gleichverteilung die Relaxationslücke eindeutig. Bei Gewichten \(w_i\) sind die beiden nichttrivialen Eigenwertfamilien \(w_i\) und \(1-w_i-w_j\). Der Mittelwert aller zehn Paarsummen beträgt \(2/5\), sodass der größte nichttriviale Eigenwert mindestens \(3/5\) erreicht. Gleichheit erfordert zehn gleiche Paarsummen. Der optimale Prozess ohne Leerschritt verfehlt damit den publizierten Endpunkt.

Der einzige stationäre Dichteoperator des Sechsfachprotokolls ist \(I_4/4\). Die endliche Transferfolge besitzt einen positiven reflektierten Kern. Mit \(\langle A,B\rangle=\operatorname{Tr}(A^\dagger B)/4\) gilt für \(t_i\in\mathbb N_0\)
\[
\sum_{ij}\langle A_i,\Phi^{t_i+t_j} A_j\rangle
=\left\|\sum_i\Phi^{t_i}A_i\right\|^2\geq0.
\]
Dies liefert eine mathematische Zeiteinheit je Lesevorgang. Eine physikalische Dauer oder ein wechselwirkendes Feldintegral folgt daraus bislang nicht. Für genau diesen Endpunkt und seine positiven ganzzahligen Potenzen schließt der bestehende Nachweis einen zeitunabhängigen GKSL Generator auf demselben \(M_4\) aus.

## Die ursprüngliche Auslesung auf einem eingebetteten \(M_2\)

Der Quellenchecker rekonstruiert aus den ursprünglichen Wörtern \(A=iW_{\mathrm{A\_BIT}}\) und \(F=W_{\mathrm{FSIG}}\) zwei antikommutierende hermitesche Involutionen. Sie erzeugen im Viererregister eine Kopie von \(M_2\) mit Multiplizität zwei. Es gilt \(\Phi(A)=2A/3\), \(\Phi(F)=F/3\) und \(\Phi(iAF)=iAF/3\). Die letzte Gleichung folgt aus dem dort geprüften Transfer des Besetzungsoperators für \(c=(F+iA)/2\).

Die drei positiven Effekte der ursprünglichen Trineauslesung sind
\[
E_1=\tfrac13(I+\tfrac{\sqrt3}{2} A+\tfrac12 F),\quad
E_2=\tfrac13(I-\tfrac{\sqrt3}{2} A+\tfrac12 F),\quad
E_3=\tfrac13(I-F).
\]
Sie summieren sich zur Identität und haben jeweils das Spektrum \((0,0,2/3,2/3)\). Die Beziehung \(\Phi(E_i)=\sum_j B_{ij}E_j\) ergibt die vollständige bekannte Übergangsmatrix
\[
B=\tfrac1{18}\begin{pmatrix}13&1&4\\1&13&4\\4&4&10\end{pmatrix}.
\]
Sechs Anwendungen liefern \(B^6=T_{\mathrm{v221}}\). Bereits der ursprüngliche Checker enthält allerdings eine Familie positiver Auslesungen mit genau derselben Übergangsmatrix. In \(E_{1,2}\) kann man die Koeffizienten von \(\pm A,F\) durch \(u,v\) und in \(E_3\) den von \(-F\) durch \(2v\) ersetzen, solange \(u^2+v^2\leq1\) und \(|2v|\leq1\).

Der ursprüngliche Befund ist eine exakte, bedingte Brücke zwischen dem markierten Compiler und dem diskreten Nahttransfer. Die physische Auswahl des Instruments und der gemeinsame unitäre Feldprozess sind offen. Der multiplikative Bereich dieser Karte besteht ausschließlich aus den Vielfachen der Identität. [Ursprünglicher Vergleich](../../../../experiments/theory-contracts/marked-seam-cp-lift-20260921/README.md) · [Quellenchecker mit den Wörtern und Effekten](../../../../experiments/theory-contracts/marked-seam-cp-lift-20260921/checker.py) · [Ausschluss eines homogenen Generators](../../../../experiments/theory-contracts/source-continuous-marked-time-search-20260922/PROOF.md)

## Warum Determinante und Pfaffsche Form verschiedene Windungen sehen

Für die fünf konkreten Matrizen hat die antiunitäre Abbildung \(\Theta=DK\), \(D=(XZ)\otimes I\), \(K\) komplexe Konjugation, die Eigenschaften \(\Theta^2=-I\) und \(\Theta\Gamma_j\Theta^{-1}=\Gamma_j\). Jeder der zehn Projektoren \(P_{j,\pm}\) ist invariant und hat komplexen Rang zwei. Eine unter \(\Theta\) invariante komplexe Gerade ist wegen der orthogonalen Kramerspartner unmöglich.

Der Zählkreis eines binären Ausgangs ist \(u_P(z)=I+(z-1)P\) mit \(z=e^{i\theta}\). Weil \(PD=DP^{\mathsf T}\), ist \(u_P(z)D\) schiefsymmetrisch. In der komplexen Viererdarstellung folgt für sämtliche zehn Vorzeichenwahlen exakt
\[
\det u_P(z)=z^2,\qquad
\frac{\operatorname{Pf}(u_P(z)D)}{\operatorname{Pf}(D)}=z .
\]
Die zweite Identität folgt auch daraus, dass die linke Seite ein Polynom ist, dessen Quadrat \(z^2\) ist und dessen Wert bei \(z=1\) eins beträgt. Der skalare Pfaffsche Kreis ist unabhängig vom gewählten Platz und Vorzeichen. Allgemeiner hat jeder stetige punktweise zeitreflektierte unitäre Kreis mit \(\Theta u(\theta)\Theta^{-1}=u(-\theta)\) gerade Determinantenwindung. Seine Werte bei \(0\) und \(\pi\) sind quaternionisch unitär, also von Determinante eins. Die beiden Halbkreise tragen dann gleiche ganzzahlige Windung.

Mit der **komplexen Spurnormierung** \(1/4\) des stationären Viererregisters ergibt die logarithmische Pfaffsche Ableitung pro Winkel
\[
\frac{1}{4}\frac{1}{2\pi i}\frac{d}{d\theta}
\log\frac{\operatorname{Pf}(u_P(e^{i\theta})D)}{\operatorname{Pf}(D)}
=\frac{1}{8\pi}.
\]
Dieselbe Formel mit der Determinante liefert \(1/(4\pi)\). Die Wahl des Pfaffschen Stroms und seiner Normierung muss am ursprünglichen Randfunktional gerechtfertigt werden. Die Gleichheit der Zahl liefert keine solche Identifikation. Der Winkel zählt ein Ereignis. Die diskrete Zeit zählt Anwendungen von \(\Phi\). [Ursprüngliche Definition von P1](../../../../documents/tfpt_1_architecture_e8/sections/001-source.tex)

## Ein konditionierter Versuch trennt zwei sonst identische Prozesse

Schreibe die Krausoperatoren als \(K_0=\sqrt{7/12}\,I\) und \(K_j=\Gamma_j/\sqrt{12}\). Sie sind orthogonal und linear unabhängig. Der Choi Rang beträgt sechs. Jede binäre Zerlegung des Kanals in zwei vollständig positive Ausgangskarten entspricht deshalb eindeutig einem Effekt \(0\leq H\leq I_6\) auf seinem minimalen Umgebungsraum. Das sind 36 reelle Freiheitsgrade im Inneren. Fordert man für jeden Eingang exakt Wahrscheinlichkeit \(1/2\) für den ersten Ausgang, bleiben lokal 20. Die 16 unabhängigen Bedingungen stammen von den Produkten \(K_a^\dagger K_b\), die den gesamten hermiteschen Operatorraum aufspannen.

Schon die unter Vorzeichenwechsel und Platzpermutation verträgliche Unterfamilie \(H=\operatorname{diag}(a,b,b,b,b,b)\) lässt eine kontinuierliche Wahl. Für \(7a+5b=6\), \(0\leq a,b\leq1\), beträgt die erste Ausgangswahrscheinlichkeit stets \(1/2\). Die erste Teilkarte und die komplementäre Teilkarte summieren sich für alle Werte zum selben \(\Phi\). Auf \(\Gamma_1\) wirkt die konditionierte erste Teilkarte mit Faktor \((7a-3b)/6\).

| Parameter \((a,b)\) | Erste Ausgangswahrscheinlichkeit | Erwartungswert von \(\Gamma_1\) nach Auswahl des ersten Ausgangs |
| --- | ---: | ---: |
| \((6/7,0)\) | \(1/2\) | \(+1\) |
| \((1/7,1)\) | \(1/2\) | \(-1/3\) |

Man präpariert dazu einen Zustand im \(+1\) Eigenraum von \(\Gamma_1\), führt einen instrumentierten Schritt aus und misst nach Auswahl des ersten Ausgangs erneut \(\Gamma_1\). Die unbedingten Transferdaten stimmen in beiden Fällen auf jedem Operator und nach beliebig vielen Schritten überein. Ein ursprünglicher, physisch identifizierter Zweizeitdatensatz mit genau dieser Präparation und dieser Ausgangsmarkierung könnte die beiden Prozesse unterscheiden. Die Werte sind Vorhersagen der angegebenen Instrumente.

Reproduktion mit Python 3 und der Standardbibliothek über [signed_clifford.py](signed_clifford.py). Die endliche Rechnung prüft den vollständigen Operatorvergleich, die zehn Pfaffschen Polynome, die Originaldarstellung von \(\Theta\) und beide konditionierten Werte.
