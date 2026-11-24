# H6h27: Verbundene Vierpunktantwort physischer Randquellen

**Ausgangspunkt:** H6h27, main b5cc3d407, insbesondere Z50 und Z51. Dieser Bericht berechnet den bisher offenen Beitrag der quartischen Potentialfreiheit zur verbundenen Antwort der physischen Quellen. Die Aussage gilt auf Baumordnung um den gemeinsamen flachen Hintergrund. Sie ist keine Behauptung über ein bereits konstruiertes Quantenfunktional.

**Ergebnis:** Der masslose Pol verschwindet aus jedem relativen Kantenprofil. Dennoch bleibt die nicht amputierte Vierpunktantwort generischer TT Randquellen im Infraroten endlich und im Allgemeinen von der Quartikkopplung abhängig. Reine masslose Eigenzustände sehen diesen Beitrag nicht. Für den unendlichen flach gewichteten Stern ergibt sich sogar eine geschlossene erzeugende Funktion aller tiefenabhängigen Quartikkopplungen.

## 1. Wirkung, Normierung und Vorzeichen

Es gibt eine dynamische Metrik je Knoten eines endlichen Baums. Eine Kante e verbindet Elternknoten v mit Kindknoten w. Verwende

\[
g_v=\eta+\frac{2h_v}{M_v},\qquad
S_{\rm EH}=\sum_v\frac{M_v^2}{2}\int\sqrt{-g_v}\,R[g_v],
\qquad
S_J=\sum_v\int \frac{h_v^{\mu\nu}J_{v,\mu\nu}}{M_v}.
\tag{1}
\]

Die physische Quellenvariable ist somit \(\phi_v=h_v/M_v\). Die Faktoren \(1/M_v\) sind Teil der Quellenregel, nicht nachträglich wählbare Residuen.

Auf jeder Kante sei

\[
S_{\rm pot}=-\sum_e k_e\int\sqrt{-g_v}\,V_{\tau_e}(X_e),
\quad X_e=\sqrt{g_v^{-1}g_w},
\]

\[
V_{\tau_e}(X)=\frac{-4+e_1(X)+e_3(X)-4e_4(X)}2
+\frac{\tau_e}{2}\det(I-X).
\tag{2}
\]

Die Planckmassen \(M_v\) und Steifigkeiten \(k_e\) seien positiv. Die verwendete Familie hat denselben flachen Hintergrund, dieselbe Fierz und Pauli Masse und dieselben kubischen Vertizes für alle \(\tau_e\). Die übliche Bimetrik ist ein zusätzlicher physikalischer Wirkungsansatz; ihr Ursprung ist nicht aus dem Hardyraum allein abgeleitet. Grundlage der verwendeten Potentialklasse: Hassan und Rosen, <https://arxiv.org/abs/1109.3515>.

Definiere den linearen relativen Kantenwert

\[
D_{e,\mu\nu}
=\frac{h_{w,\mu\nu}}{M_w}-\frac{h_{v,\mu\nu}}{M_v}.
\tag{3}
\]

Bei Lorentzsignatur \((-+++)\) gilt

\[
X_e-I=\eta^{-1}D_e+O(h^2),\qquad
\det(I-X_e)=-\det(D_{e,\mu\nu})+O(h^5).
\]

Der Faktor \(\sqrt{-g_v}=1+O(h)\) ändert die vierte Ordnung nicht. Deshalb lautet die Änderung der **lorentzschen Lagrangedichte**

\[
\boxed{\Delta\mathcal L_L^{[4]}
=+\frac12\sum_e k_e\tau_e\det(D_{e,\mu\nu}).}
\tag{4}
\]

Nach konsistenter Fortsetzung auch der Tensorindizes lautet die entsprechende **euklidische Wirkung**

\[
\boxed{\Delta S_E^{[4]}
=+\frac12\sum_e k_e\tau_e\int d^4y\,\det D_e(y).}
\tag{5}
\]

Für den verbundenen euklidischen Funktional verwenden wir

\[
W_E[J]=\log\int Dh\,\exp(-S_E[h]+S_J[h,J]).
\tag{6}
\]

Die folgende Quellenantwort erhält daher ein Minuszeichen aus \(-\Delta S_E\). Dieses Minuszeichen darf nicht ungeprüft als Vorzeichen einer lorentzschen Streuamplitude übernommen werden. Dort gehören Tensorfortsetzung, Feynmanfaktoren und gegebenenfalls LSZ Amputation dazu. Die volle euklidische Gravitation wird hier nicht als bereits konvergentes positives Maß vorausgesetzt; die Rechnung ist eine klassische beziehungsweise perturbative Baumrechnung.

## 2. Exakte lineare Profile mit physischen Quellen

Sei \(I_{ev}=+1\) am Elternknoten, \(I_{ew}=-1\) am Kindknoten und sonst null. Ferner

\[
\mathsf M=\operatorname{diag}(M_v),\quad
K=\operatorname{diag}(k_e),\quad
Q=I^T K I,\quad
L=\mathsf M^{-1}Q\mathsf M^{-1}.
\tag{7}
\]

Für jede externe Quelle wird angenommen

\[
p^\mu J_{\mu\nu}(p)=0,\qquad
\delta^{\mu\nu}J_{\mu\nu}(p)=0.
\tag{8}
\]

Der TT Projektor wirkt dann identisch auf sie. Die skalare kanonische Kovarianz und die physische Quellenkovarianz sind

\[
G(x)=(x\,1+L)^{-1},\qquad
C^\phi(x)=\mathsf M^{-1}G(x)\mathsf M^{-1}
=(x\mathsf M^2+Q)^{-1},\quad x=p_E^2>0.
\tag{9}
\]

Eine physische Quelle am Knoten b erzeugt an der Kante e das Profil

\[
\boxed{
R_{eb}(x)=
\frac{G_{wb}(x)}{M_wM_b}
-\frac{G_{vb}(x)}{M_vM_b}
=C^\phi_{wb}(x)-C^\phi_{vb}(x).
}
\tag{10}
\]

Damit ist

\[
D_e^{(1)}(p)=\sum_bR_{eb}(p^2)J_b(p).
\tag{11}
\]

Diese Formel enthält beide Normierungen: \(1/M_b\) an der Quelle und \(1/M_v\) beziehungsweise \(1/M_w\) in der relativen Metrik. Die Identifikation eines internen Stroms mit genau dieser Tensorquelle bleibt eine vorgelagerte physikalische Aufgabe.

## 3. Der verbundene Vierpunktbeitrag

Mit \([4]\) bezeichnen wir den homogenen Anteil vierter Ordnung in J. Auf Baumordnung gilt exakt

\[
\boxed{
\Delta W_E^{[4]}[J]
=-\frac12\sum_e k_e\tau_e
\int d^4y\,
\det\!\left[
\int\frac{d^4p}{(2\pi)^4}e^{ipy}
\sum_bR_{eb}(p^2)J_b(p)
\right].
}
\tag{12}
\]

Hier bedeutet \(\Delta W=W_\tau-W_0\). Der Ausdruck ist nicht die gesamte Vierpunktfunktion.

**Beweis:** Der neue Vertex beginnt quartisch. Der lineare klassische Lösungsanteil \(h^{(1)}=G\mathsf M^{-1}J\) genügt daher in diesem Vertex. Höhere Lösungsanteile erzeugen mindestens fünf Quellen. Die Änderung der klassischen Lösung selbst trägt in vierter Ordnung nicht zusätzlich bei, weil die Variation der ursprünglichen Wirkung minus Quellen am stationären Punkt verschwindet. Alle Austauschdiagramme aus zwei kubischen Vertizes bleiben unverändert, ebenso alle übrigen quartischen Kontakte. Bei vier externen Quellen enthält der neue Beitrag genau einen quartischen Vertex mit vier linearen Quellenbeinen. Das gilt auf dieser Ordnung auch bei einer endlichen Änderung von \(\tau_e\).

Um jede Fakultät eindeutig festzulegen, definiere

\[
\mathcal V_{\det}(T_1,T_2,T_3,T_4)
=[z_1z_2z_3z_4]\det\!\left(\sum_{i=1}^4z_iT_i\right).
\tag{13}
\]

Dieser Ausdruck ist symmetrisch und multilinear. Für vier gleiche Argumente ist er \(24\det T\). Äquivalent ist er die doppelte Levi und Civita Kontraktion der vier Tensoren, ohne zusätzliche Fakultät.

Für Quellen \(J_{b_i}(p_i)=s_iT_i\) lautet die vierfache Ableitung von W:

\[
\boxed{
\Delta\Gamma^{(4)}_{J,E}(p_i;b_i,T_i)
=-\frac12\sum_e k_e\tau_e
\left[\prod_{i=1}^4R_{eb_i}(p_i^2)\right]
\mathcal V_{\det}(T_1,T_2,T_3,T_4).
}
\tag{14}
\]

Der Faktor \((2\pi)^4\delta^{(4)}(\sum_i p_i)\) ist in (14) unterdrückt. In der üblichen Darstellung \(W^{[4]}=(1/4!)\Gamma^{(4)}[J,J,J,J]\) ist keine weitere Fakultät einzusetzen.

Für nicht TT Quellen muss (11) durch die vollständigen Tensorpropagatoren mit den jeweiligen Spuranteilen ersetzt werden. Die skalare Formel allein zertifiziert solche Quellen nicht.

## 4. Kontrolle an einer einzelnen Kante

Für zwei Metriken g und f gilt

\[
m^2=k\left(M_g^{-2}+M_f^{-2}\right),\qquad
M_{\rm eff}^2=M_g^2+M_f^2.
\]

Direkte Inversion von (9) ergibt

\[
\boxed{
R_{eg}(x)=-\frac1{M_g^2(x+m^2)},\qquad
R_{ef}(x)=+\frac1{M_f^2(x+m^2)}.
}
\tag{15}
\]

Obwohl jede einzelne Metrik einen masslosen Pol besitzt, enthalten diese relativen Profile nur den massiven Propagator. Vier g Quellen liefern

\[
\boxed{
\Delta\Gamma^{(4)}_{J,E}
=-\frac{k\tau}{2M_g^8}
\frac{\mathcal V_{\det}(T_1,T_2,T_3,T_4)}
{\prod_i(p_i^2+m^2)}.
}
\tag{16}
\]

Für eine f Quelle ändert sich der zugehörige Faktor von \(-1/M_g^2\) zu \(+1/M_f^2\). Das ist eine physische Quellenabhängigkeit und kein negatives Spektralgewicht.

## 5. Warum der masslose Pol exakt verschwindet

Auf jedem endlichen zusammenhängenden Baum besitzt L genau die Nullmode

\[
\eta_v=\frac{M_v}{M_{\rm eff}},
\qquad M_{\rm eff}^2=\sum_vM_v^2.
\]

Die universelle masslose Kopplung ist \(\eta_v/M_v=1/M_{\rm eff}\). Daher

\[
C^\phi_{vb}(x)=\frac1{M_{\rm eff}^2x}
+C^\phi_{{\rm reg},vb}(x).
\tag{17}
\]

In der Differenz (10) kürzt sich der erste Term für jede Quelle b exakt. Da der endliche Baum eine positive Lücke oberhalb seiner einzigen Nullmode hat, ist \(R_{eb}(x)\) bei null analytisch.

Es gibt eine nützliche geschlossene Formel für diesen Grenzwert. Sei \(S_e\) der vollständige Teilbaum unterhalb der Kante e und

\[
M_{S_e}^2=\sum_{v\in S_e}M_v^2.
\]

Für beliebige verteilte physische Quellen folgt

\[
\boxed{
k_eD_e^{(1)}(0)
=\sum_{v\in S_e}J_v
-\frac{M_{S_e}^2}{M_{\rm eff}^2}\sum_vJ_v.
}
\tag{18}
\]

**Beweis:** Die physische lineare Gleichung ist \((x\mathsf M^2+Q)\phi=J\). Setze \(\phi=(\sum J)/(M_{\rm eff}^2x)\,1+\phi_0+O(x)\). Dann

\[
Q\phi_0=J-\mathsf M^2\,1\,\frac{\sum J}{M_{\rm eff}^2}.
\]

Summiert man diese Gleichung über den Teilbaum \(S_e\), verschwinden alle inneren Kantenbeiträge. Übrig bleibt genau \(k_e(\phi_{0,w}-\phi_{0,v})\), also (18).

Für eine Quelle nur am Wurzelknoten c gilt insbesondere

\[
\boxed{R_{ec}(0)=-\frac{M_{S_e}^2}{k_eM_{\rm eff}^2}.}
\tag{19}
\]

Nullimpuls bedeutet in (18) und (19) immer den regulären Kantenlimes. Die singuläre volle Kovarianz \(C^\phi(0)\) wird nicht invertiert.

## 6. Drei verschiedene Aussagen über kleine Impulse

Halte die Quellentensoren und ihre Richtungen fest und skaliere \(p_i=\epsilon\widehat p_i\). Für generische nichtverschwindende Tensorstruktur gilt:

| Größe | Führende Ordnung des neuen Beitrags |
|---|---|
| Nicht amputierte verbundene Randantwort | \(O(\epsilon^0)+O(\epsilon^2)\) |
| Mit vier inversen Randkovarianzen amputierter Vertex | \(O(\prod_i p_i^2)=O(\epsilon^8)\) |
| Exakt auf masslose Eigenzustände projizierter Kontakt | Identisch null |

Die zweite Zeile folgt aus \(C^\phi_{cc}(p^2)^{-1}=M_{\rm eff}^2p^2+O(p^4)\). Sie beschreibt einen Vertex der Randvariable, keine reine masslose Eigenmode bei endlichem Impuls.

Die dritte Zeile ist stärker: Für eine reine Nullmode gilt \(h_v=\eta_vH_0\), also \(D_e=0\) auf jeder Kante. Wegen der Multilinearität verschwindet (14) schon dann, wenn **ein einziges** der vier externen Beine eine reine masslose Eigenmode ist. Für LSZ mit der Randvariable ist der masslose Residuenfaktor \(Z_\phi=1/M_{\rm eff}^2\); die vier Faktoren \(p_i^2/\sqrt{Z_\phi}\) lassen den endlichen Beitrag ebenfalls verschwinden.

### 6.1 Ein tatsächlicher TT Zeuge

Der Nachweis der Sichtbarkeit darf nicht nur beliebige Matrizen verwenden. Wähle euklidische Impulse

\[
p_1=\epsilon(1,0,0,0),\quad
p_2=\epsilon(0,1,0,0),\quad
p_3=-p_1,\quad p_4=-p_2,
\]

und

\[
T_1=T_3=A=\operatorname{diag}(0,1,1,-2),\qquad
T_2=T_4=B=\operatorname{diag}(1,0,1,-2).
\tag{20}
\]

Jeder Tensor ist spurfrei und transversal zu seinem eigenen Impuls; die Impulssumme ist null. Trotzdem ist

\[
\det(aA+bB)=-2ab(a+b)^2,\qquad
\boxed{\mathcal V_{\det}(A,B,A,B)=-16.}
\tag{21}
\]

Damit ist die nichtverschwindende Antwort für echte TT Quellen ausdrücklich belegt. Sind hingegen alle Impulse kollinear, besitzen sämtliche TT Tensoren einen gemeinsamen Nullvektor. Dann verschwindet die gesamte Determinantenpolarisierung. Diese spezielle Kinematik eignet sich nicht zur Bestimmung von \(\tau\).

Der endliche konstante Grenzwert ist ein analytischer Kontakt im effektiven Funktional der Quellen. Er ist kein zusätzlicher massloser Langreichweitenanteil. Bei räumlich getrennten Quellen beschreibt der vollständige Kern massive Ausbreitung, nicht einen neu entstandenen \(1/p^2\) Pol. Bei einem mikroskopischen Vergleich müssen Quellenkonvention und lokale Kontaktgegenstücke im selben Schema festgelegt werden.

## 7. Endlicher Stern, beliebige Tiefenwerte

Setze \(\lambda=\mu^2\), \(r=q^2\), und

\[
M_{j,n}=\frac{M_c}{2}q^{n+1},\qquad
k_{j,n}=\lambda w_nM_{j,n}^2,
\quad n=0,\ldots,m-1.
\tag{22}
\]

Der ursprüngliche Stern hat \(w_n=1\); der quadratische Elternoperator der zwei Strommoden führt zum gesondert angenommenen Gewicht \(w_n=n+1\).

Schreibe \(S_\ell(r)=\sum_{a=0}^{\ell-1}r^a\). Dann

\[
M_{\rm eff}^2=M_c^2S_{m+1},\qquad
M_{S_{j,n}}^2=M_{j,n}^2S_{m-n}.
\]

Für zentrale Quellen liefert (19)

\[
\boxed{
R_{(j,n),c}(0)
=-\frac{S_{m-n}(r)}
{\lambda w_nM_c^2S_{m+1}(r)}.
}
\tag{23}
\]

Definiere den Armmittewert \(\overline\tau_n=\frac14\sum_{j=0}^3\tau_{j,n}\). Die zentrale Antwort hängt nur von diesem Mittelwert ab, weil ihre linearen Profile auf allen Armen identisch sind. Der skalare Vorfaktor von \(\mathcal V_{\det}\) lautet

\[
\boxed{
F_m(0)=
-\frac1{2\lambda^3M_c^6S_{m+1}^4}
\sum_{n=0}^{m-1}
\frac{\overline\tau_n\,r^{n+1}S_{m-n}^4}{w_n^3}.
}
\tag{24}
\]

Eine ungebrochene Vierteldrehung fordert gleiche \(\tau\) auf den vier Armen bei derselben Tiefe. Sie identifiziert aber keine verschiedenen Tiefen. Auch die Gleichheit des Potentials beim Umkehren einer einzelnen Kante legt keine Beziehung zwischen verschiedenen Kanten fest. Ein einziges gemeinsames \(\tau\) für den gesamten Stern ist daher eine zusätzliche Homogenitätsregel.

Als einfache Kontrolle hat ein Stern mit nur einer Kante je Arm das exakte Profil

\[
R_{(j,0),c}(x)
=-\frac1{M_c^2[x+\lambda(1+r)]}.
\tag{25}
\]

Sein Infrarotkontakt bei gemeinsamem \(\tau\) ist
\(-\tau r/[2\lambda^3M_c^6(1+r)^4]\), in Übereinstimmung mit (24).

### 7.1 Unendlicher Infrarotlimes

Für \(0<r<1\), beschränkte \(\tau_n\) und beispielsweise \(w_n\ge1\) ist die Summe absolut dominiert. Für jede feste Tiefe geht \(S_{m-n}/S_{m+1}\) gegen eins. Daher

\[
\boxed{
F_\infty(0)=
-\frac1{2\lambda^3M_c^6}
\sum_{n=0}^{\infty}\frac{\overline\tau_n r^{n+1}}{w_n^3}.
}
\tag{26}
\]

Für ein gemeinsames \(\tau\) sind die zwei besonderen Fälle

\[
F_{\rm flat}(0)=-\frac{\tau}{2\lambda^3M_c^6}\frac r{1-r},
\qquad
F_{\rm weighted}(0)=-\frac{\tau}{2\lambda^3M_c^6}\operatorname{Li}_3(r).
\tag{27}
\]

Diese Werte betreffen den einzelnen klassischen Quartikkontakt. Sie beweisen keine Kontrolle über Quantenkorrekturen der tiefen Knoten.

## 8. Geschlossene Vierpunktfunktion des unendlichen flachen Sterns

Für \(w_n=1\) lässt sich die Summe sogar bei allgemeinen euklidischen äußeren Impulsen exakt ausführen. Definiere

\[
x=\frac{p^2}{\lambda},\qquad
z(x)=
\frac2{x+1+r+\sqrt{(x+1+r)^2-4r}}.
\tag{28}
\]

Für \(x>0\) gilt \(0<z(x)<1\), und \(z(0)=1\) ist der Grenzwert. Die abfallende Lösung der radialen Rekursion ist \(qz(x)\). Aus

\[
rz^2-(x+1+r)z+1=0,\qquad
x=\frac{(1-z)(1-rz)}z
\]

folgt für die rohe Randresolvente

\[
C_{\rm raw}(x)=\frac z{1-z}.
\tag{29}
\]

Die kanonische Zentrumsantwort ist \(\lambda^{-1}C_{\rm raw}(p^2/\lambda)\); die physische Quellenkovarianz enthält zusätzlich \(M_c^{-2}\). Auf einem einzelnen Arm ist der kanonische Fourierfaktor \(1/2\) zu berücksichtigen. Er kürzt sich gegen denselben Faktor in \(M_{j,n}\). Damit erhält man

\[
\phi_{j,n}^{(1)}(p)
=\frac{C_{\rm raw}(x)}{\lambda M_c^2}z(x)^{n+1}J_c(p),
\]

\[
\boxed{
R_{(j,n),c}(p^2)
=\frac{C_{\rm raw}(x)}{\lambda M_c^2}[z(x)-1]z(x)^n
=-\frac{z(x)^{n+1}}{\lambda M_c^2}.
}
\tag{30}
\]

Der letzte Schritt verwendet \(C_{\rm raw}(z-1)=-z\). Die Aufhebung des masslosen Pols ist hier somit bei jedem Impuls ausdrücklich sichtbar.

Setze

\[
P=\prod_{i=1}^4z(p_i^2/\lambda).
\]

Für tiefenabhängige Quartikkopplungen lautet (14)

\[
\boxed{
\Delta\Gamma^{(4)}_{J,E}
=F(P)\,\mathcal V_{\det}(T_1,T_2,T_3,T_4),\qquad
F(P)=-\frac1{2\lambda^3M_c^6}
\sum_{n=0}^{\infty}\overline\tau_n(rP)^{n+1}.
}
\tag{31}
\]

Für beschränkte \(\overline\tau_n\) konvergiert diese Reihe absolut für die gesamte euklidische Region \(0\le P\le1\). Bei einem gemeinsamen \(\tau\) wird sie geometrisch:

\[
\boxed{
F(P)=-\frac{\tau}{2\lambda^3M_c^6}
\frac{rP}{1-rP}.
}
\tag{32}
\]

Dies ist die vollständige Änderung dieses verbundenen Baumkontakts, mit vier äußeren linearen Profilen, nicht bloß eine Entwicklung bei null.

Aus \(z(x)=1-x/(1-r)+O(x^2)\) folgt

\[
F(P)=
-\frac{\tau r}{2\lambda^3M_c^6(1-r)}
+\frac{\tau r}{2\lambda^3M_c^6(1-r)^3}
\sum_i\frac{p_i^2}{\lambda}
+O(p^4).
\tag{33}
\]

Bei festem \(\lambda,M_c,\tau\) beginnt der kleine q Beitrag deshalb mit \(q^2\).

### 8.1 Eine kontrollierte Aussage zur Tiefe

Beschränkt man nur die Quartikkontakte auf die ersten N Tiefen und lässt den freien unendlichen Elternoperator unverändert, ist für ein gemeinsames \(\tau\) der ausgelassene Rest exakt

\[
F_{\rm tail}(P)
=-\frac{\tau}{2\lambda^3M_c^6}
\frac{(rP)^{N+1}}{1-rP}.
\tag{34}
\]

Für \(|\tau_n|\le T\) gilt allgemeiner

\[
|F_{\rm tail}(P)|
\le\frac{T}{2\lambda^3M_c^6}\frac{r^{N+1}}{1-r}.
\tag{35}
\]

Dies ist keine Gleichsetzung mit einem endlichen physischen Stern: Dessen freier Randoperator und Endbedingung ändern sich ebenfalls, weshalb er mit (9) und (14) gerechnet werden muss. Ebenso ersetzt die Schranke keinen Nachweis, dass sämtliche Schleifen und Wechselwirkungen der stark gekoppelten tiefen Region klein sind.

## 9. Wie die Kopplung tatsächlich bestimmt werden könnte

Die Rechnung liefert einen konkreten Matchingkanal. Sie liefert ohne unabhängige mikroskopische Vierpunktfunktion noch keinen Wert.

Zuerst müssen dieselbe physische Quellenzuordnung, dieselben Normierungen und dieselben zwei und drei Punktdaten verwendet werden. Dann ist die vollständige Vierpunktfunktion bei \(\tau=0\), einschließlich ihrer Austauschdiagramme und sonstigen Kontakte, abzuziehen. Nur der Rest ist mit (14) beziehungsweise (31) zu vergleichen.

### 9.1 Test der Annahme eines gemeinsamen \(\tau\)

Für zentrale Quellen mit \(\mathcal V_{\det}\ne0\) definiere den skalaren Rest

\[
F_{\rm micro}(p_i,T_i)
=\frac{\Gamma_{\rm micro}^{(4)}
-\Gamma_{\tau=0}^{(4)}}{\mathcal V_{\det}(T_i)}.
\]

Falls die Familie mit gemeinsamem \(\tau\) zutrifft, muss

\[
\boxed{
\tau=
-2\lambda^3M_c^6\,\frac{1-rP}{rP}\,
F_{\rm micro}(p_i,T_i)
}
\tag{36}
\]

bei sämtlichen zulässigen Impulsen und nichtverschwindenden Tensorprojektionen denselben Wert liefern. Unterschiedliche Werte oder zusätzliche unabhängige Tensorstrukturen widerlegen diese eingeschränkte Wechselwirkungsfamilie.

### 9.2 Bestimmung der tiefenabhängigen Werte

Ohne Homogenitätsregel ist \(F(P)\) eine erzeugende Funktion. Kennt man sie unabhängig auf einem offenen Intervall und erfüllt sie die bewiesene analytische Form, sind ihre Koeffizienten eindeutig:

\[
\boxed{
\overline\tau_n
=-\frac{2\lambda^3M_c^6}{r^{n+1}}\,
[P^{n+1}]F(P).
}
\tag{37}
\]

Die TT Konfiguration (20) mit variablem \(\epsilon>0\) liefert \(P=z(\epsilon^2/\lambda)^4\) und durchläuft ein offenes Intervall. Damit ist diese Rekonstruktion nicht nur eine formale Wahl unabhängiger äußerer Impulse.

Das ist eine Eindeutigkeitsaussage bei exakten Daten. Die praktische Bestimmung sehr tiefer Koeffizienten ist wegen der kleinen Faktoren \(r^{n+1}\) schlecht konditioniert. Der bloße Infrarotwert bestimmt ohnehin nur die gewichtete Summe (26). Beispielsweise verändern \(\delta\tau_0=r\) und \(\delta\tau_1=-1\) den flachen Infrarotwert überhaupt nicht, liefern bei endlichen Impulsen aber einen nichtverschwindenden Rest proportional zu \(r^2P(1-P)\).

Ohne Viertelsymmetrie können rein zentrale Quellen nur die Armmittewerte bestimmen. Um Unterschiede zwischen Armen zu erkennen, benötigt man zusätzliche nichtzentrale Quellenprofile.

### 9.3 Warum der vorhandene Gaußzustand keinen Wert setzt

Ein Gaußzustand eines einzelnen Strommodus legt nicht automatisch die verbundene gravitative Vierpunktfunktion fest. Der metrische Quellenoperator und die nichtlineare Prozessregel sind zusätzliche Informationen.

Würde man stattdessen einen vollständig gaußschen Prozess mit linearem physischem Quellenabbild fordern, verschwänden alle verbundenen höheren Kumulanten. Das wäre bei endlichem G bereits mit den nichtverschwindenden Einstein Vertizes auf dritter Ordnung unvereinbar. Es würde nicht lediglich \(\tau=0\) auswählen. Außerdem blieben bei \(\tau=0\) die anderen Kontakte und Austauschdiagramme der Vierpunktfunktion bestehen.

Auch Positivität des Zweipunkts oder universelle masslose Quellenresiduen bestimmen (36) nicht. Bedingungen aus höheren Amplituden können den zulässigen Bereich einschränken; bei einem masslosen Graviton müssen die Voraussetzungen gewöhnlicher Vorwärtspositivität eigens geprüft werden. Der masslose Kanal verhindert eine automatische Übernahme einfacher Positivitätsschranken. Primärquelle: <https://arxiv.org/abs/2007.12667>. Die unabhängige eikonale Analyse einer verbleibenden Familie von Bimetriken findet sich in <https://arxiv.org/abs/1712.10020>; sie ersetzt ebenfalls keine Bestimmung dieser gesamten Vierpunktfunktion aus TFPT.

## 10. Reproduzierbare Prüfung und Reichweite

Die neuen Dateien sind research/originlab/process_boundary_fourpoint.py und research/originlab/test_process_boundary_fourpoint.py. Der Helfer boundary_fourpoint_audit() liefert direkt JSON verträgliche Aussagen und 14 exakte algebraische Kontrollen.

Die fokussierte Testsuite prüft 14 Gruppen: physische Faktoren auf einer einzelnen Kante, allgemeine äußere Impulse, Vergleich des vollständigen endlichen Massenoperators mit den vorhandenen Sternmodulen, die Teilbaumidentität auch für nichtzentrale Quellen, exakte Polaufhebung, reine Nullmoden, die Fakultät des Determinantenvertex, tatsächlich zulässige TT Tensoren, die achte weiche Ordnung der Randamputation, unabhängige Kantenkoeffizienten, Armmittewerte, die geschlossene unendliche Summe sowie die fehlende Bestimmbarkeit der Tiefenwerte aus dem Infrarotwert allein. Alle 14 Testgruppen bestanden.

Die bewiesene Aussage ist ein expliziter klassischer Vierpunktkanal unter den angegebenen metrischen Eingaben. Noch nicht bewiesen sind das gemeinsame mikroskopische Funktional mit Naht, chiraler Materie, Eichfeldern und Higgs, eine Auswahl der gesamten Folge \(\tau_n\), die Werte von \(M_c\) und \(\lambda\), oder die vollständige Quantenkontrolle des unendlichen Sterns.
