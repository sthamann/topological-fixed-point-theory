# Zwei schaltbare Anregungen in derselben thermischen Quelle

## Ergebnis

Die Rechnung liefert einen sauberen negativen Befund für die einfachste Spin 2 Vermutung. Zwei unabhängig schaltbare Stromanregungen existieren in derselben thermischen Quelle. Ihr gemeinsamer Zustand liegt auf chiralem Niveau zwei, die OS Form ist positiv und beide Einfügungsordnungen stimmen exakt überein. Trotzdem entsteht kein neuer gemeinsamer Kanal. Die normierte Vierfachantwort faktorisiert vollständig:

\[
R_{11}=R_{10}R_{01},\qquad
\log\frac{R_{11}R_{00}}{R_{10}R_{01}}=0.
\]

Das Niveau zwei kommt hier nur durch Addition zweier Gewicht eins Anregungen zustande. Es liefert weder den Sugawara Stresstensor noch einen vierdimensionalen Spin 2 Zustand.

Die weiterführende Paarleiterhypothese ist mathematisch tragfähig, sobald man auf (n_A=n_B) projiziert. Genau diese Projektion ist aber kein Zustandsgesetz der thermischen Quarterquelle. Der unabhängige KMS Produktzustand liegt nur mit Wahrscheinlichkeit (	anh(2\pi)) im neutralen Differenzzahlsektor. Seine Dominanz kommt fast vollständig vom Vakuum und erzeugt keine Tomita Paarung.

## Tatsächliche Quelle

Ich verwende genau die in Z58 rekonstruierte quasifreie thermische Quelle, keine neu gewählte Oszillatorquelle:

\[
C_\alpha(r)=\frac{1}{1+e^{4\pi(r+\alpha)}},
\qquad \alpha\in\left\{\frac14,\frac34\right\}.
\]

Aus ihr folgt auf dem Einteilchenraum

\[
h_\alpha(r)=\log\frac{1-C_\alpha(r)}{C_\alpha(r)}
=4\pi(r+\alpha).
\]

Das Jacobi Tripelprodukt in Z58 zerlegt dieselbe Quelle in den Ladungsanteil und die thermischen Stromoszillatoren. In der dort festgelegten See Normalordnung lautet der volle Generator

\[
K=4\pi\sum_{\alpha=1/4,3/4}
\left[
\frac{Q_\alpha^2}{2}
+\left(\alpha-\frac12\right)Q_\alpha
+\sum_{m\geq1}mN_{\alpha,m}
\right].
\]

Die Einfügungen sind

\[
A=J^{(1/4)}_{-1},\qquad
B=J^{(3/4)}_{-1}.
\]

Sie sind gerade, ladungsneutral und gehören zur in Z58 ausdrücklich aufgebauten Stromalgebra. Weil sie in verschiedenen Quarterblöcken der vollen Quelle liegen, sind sie unabhängig schaltbar und kommutieren. Die gemeinsame Domäne ist der endliche Teilchenraumkern aus endlichen Stromwörtern. Beide Einfügungen sitzen bei euklidischer Zeit \(\tau=\beta/4\), mit \(\beta=4\pi\), im positiven Halbkreis \(0<\tau<\beta/2\). Die thermische Spiegelung ist \(\theta(\tau)=\beta-\tau\).

Der verwendete Kernel ist der periodische Stromkern aus Z84,

\[
S_J(\Delta\vartheta,\Delta\tau)
=\sum_{m\ne0}e^{im\Delta\vartheta-m\Delta\tau}
\frac{m}{1-e^{-\beta m}},
\qquad 0<\Delta\tau<\beta.
\]

Für die Fouriermode \(m=1\) und den Spiegelabstand \(\beta/2\) ergibt er

\[
g=\frac{1}{2\sinh(\beta/2)}
=0{,}0018674492441428359842.
\]

Quellen: [Z58, Generator, Stromwörter und Faktorisierung](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/exploration/2026-11-23-phasenraum-f5-holocode/z58_rekonstruktion.py:63), [Z84, wirklicher thermischer Stromkern und OS Gram](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/exploration/2026-11-23-phasenraum-f5-holocode/z84_euklidisch_os.py:418).

## Vier Antworten und beide Reihenfolgen

Die OS Antwort eines Wortes \(W\) ist

\[
R(W)=Z_0^{-1}\operatorname{Tr}
\left(e^{-K/2}W^\dagger e^{-K/2}W\right).
\]

| Schalter | Wort | normierte Antwort | absolute Zweimodenantwort bei Energie 1 |
|---|---:|---:|---:|
| 00 | \(1\) | 1 | 1,0000069747211972578 |
| 10 | \(A\) | 0,0018674492441428359842 | 0,0018674622690806639102 |
| 01 | \(B\) | 0,0018674492441428359842 | 0,0018674622690806639102 |
| 11 | \(AB\) | 0,0000034873666794496494375 | 0,0000034873910028599512054 |
| 11, umgekehrt | \(BA\) | 0,0000034873666794496494375 | 0,0000034873910028599512054 |

Die beiden 11 Antworten sind exakt gleich. Das ist keine numerische Näherung, denn \(A\) und \(B\) sind gerade Operatoren in verschiedenen Blöcken. Der verbundene logarithmische Anteil ist exakt null.

## OS Gram Matrix

In der Wortbasis \((1,A,B,AB,BA)\) ist

\[
G_{ij}=R(W_i,W_j)=
\begin{pmatrix}
1&0&0&0&0\\
0&g&0&0&0\\
0&0&g&0&0\\
0&0&0&g^2&g^2\\
0&0&0&g^2&g^2
\end{pmatrix}.
\]

Ihre Eigenwerte sind

\[
0,\quad 2g^2,\quad g,\quad g,\quad 1.
\]

Die Matrix ist positiv semidefinit. Der einzige Nullvektor ist \(AB-BA\). Die OS Rekonstruktion identifiziert also die beiden Schreibweisen desselben Zustands, erzeugt aus der Reihenfolge aber keinen zusätzlichen Freiheitsgrad.

## Zeit und die zwei KMS Ordnungen

Die unsymmetrischen Antworten lauten

\[
\omega(J_{-1}J_1)=\frac{e^{-4\pi}}{1-e^{-4\pi}}
=0{,}0000034873545178081166,
\]

\[
\omega(J_1J_{-1})=\frac{1}{1-e^{-4\pi}}
=1{,}0000034873545178081.
\]

Ihr Verhältnis ist \(e^{4\pi}=286751{,}3131366533\). Damit rekonstruiert dieselbe Antwort die modulare Frequenz \(4\pi\). Der OS Wert ist genau das geometrische Mittel der beiden KMS Ordnungen.

Das ist thermische modulare Zeit. Für die volle thermische Quelle ist der rekonstruierte Liouville Generator zweiseitig. Eine positive physische Zeit entsteht nach Z84 erst im Vakuumgrenzfall. Die Rechnung liefert daher keinen neuen Hamiltonoperator für eine vierdimensionale Raumzeit.

## Absolutes und normiertes \(Z_0\)

Eine absolute Zahl nenne ich nur für den tatsächlich gerechneten Zweimodenfaktor. Gewählt sind die Energie \(\varepsilon=1\) je Modenquant und

\[
K_{AB}=4\pi(N_A+N_B),\qquad
q=e^{-4\pi}.
\]

Dann ist

\[
Z_{0,AB}=\operatorname{Tr}e^{-K_{AB}}
=\frac1{(1-q)^2}
=1{,}0000069747211972578.
\]

Ladungssektoren, höhere Strommoden und alle anderen Quellenfaktoren sind in dieser absoluten Zahl ausdrücklich nicht enthalten. Sie kürzen sich in den normierten Zweimodenantworten heraus. Eine Verschiebung \(K\mapsto K+c\) verändert auch diesen absoluten Wert, während die normierten Antworten und die OS Gram Matrix invariant bleiben.

## Polarshift der Paarleiter

Für die bosonisierten Moden schreibe

\[
B^\dagger=a^\dagger b^\dagger.
\]

Seine Polarzerlegung ist

\[
B^\dagger=S\sqrt{(N_A+1)(N_B+1)},\qquad
S|n_A,n_B\rangle=|n_A+1,n_B+1\rangle.
\]

Auf dem vollen Produktfockraum gilt

\[
S^\dagger S=1,\qquad
SS^\dagger=1-P_\partial,
\]

mit

\[
P_\partial=\mathbf1_{\{n_A=0\ \mathrm{oder}\ n_B=0\}}.
\]

Der Defektraum ist unendlichdimensional. Die Differenzzahl

\[
D=N_A-N_B
\]

wird von \(S\) erhalten, und jeder Sektor \(D=d\), \(d\in\mathbb Z\), trägt einen eigenen unilateralen Shift. Der volle Shift ist deshalb kein einzelner Fredholmoperator mit Index minus eins. Erst auf dem neutralen Sektor

\[
\mathcal H_0=\operatorname{span}\{|n,n\rangle:n\ge0\}
\]

gilt

\[
\ker S_0=0,\qquad
\ker S_0^\dagger=\mathbb C|0,0\rangle,qquad
\operatorname{ind}S_0=-1.
\]

Hier ist die Hardy beziehungsweise Szegő Struktur tatsächlich vorhanden. Sie folgt jedoch erst nach der Wahl von \(P_0=\mathbf1_{D=0}\).

## Wie viel des thermischen Zustands ist neutral?

Der unabhängige Produktzustand der beiden Moden hat

\[
p(n_A,n_B)=(1-q)^2q^{n_A+n_B}.
\]

Daraus folgen exakt

\[
\omega(P_0)
=(1-q)^2\sum_{n\ge0}q^{2n}
=\frac{1-q}{1+q}
=\tanh(2\pi)
=0{,}9999930253396106106,
\]

und

\[
\omega(1-P_0)
=\frac{2q}{1+q}
=0{,}000006974660389389395.
\]

Die große Neutralwahrscheinlichkeit ist kein dynamisches Matching. Das Vakuum allein hat schon das Gewicht

\[
p_{00}=(1-q)^2=0{,}9999930253274491387.
\]

Das Gewicht angeregter neutraler Paare beträgt nur

\[
\omega(P_0)-p_{00}
=1{,}2161471886681454\cdot10^{-11}.
\]

Damit ist (P_0) als Spektralprojektion der rekonstruierten Oszillatordifferenz mathematisch verfügbar. Es ist jedoch keine native Zustandsbedingung der Quarterquelle. Insbesondere ist die Z76 Uhr (Q_{1/4}-Q_{3/4}\) eine Differenz der Ladungsnullmoden. Sie ist nicht dieselbe Größe wie (N_A-N_B) der beiden hier ausgewählten Stromoszillatoren.

Für den Polardefekt erhält man

\[
\omega(P_\partial)=1-q^2
=0{,}9999999999878384433,
\]

während der Innenraum \(n_A,n_B\ge1\) nur das Gewicht \(q^2=1{,}2161556709409308\cdot10^{-11}\) trägt. Der Grundzustand induziert daher keine Konditionierung auf \(D=0\). Er macht lediglich sowohl Neutralität als auch den Polarrand numerisch fast sicher.

## Drei verschiedene Präparationen

Die folgenden drei Konstruktionen dürfen nicht gleichgesetzt werden.

### Unabhängige kohärente Quellen

Das Produkt zweier Exponentialvektoren hat vor der Projektion die Koeffizienten

\[
\frac{\alpha^{n_A}\beta^{n_B}}{\sqrt{n_A!n_B!}}.
\]

Nach \(P_0\) bleiben

\[
P_0|\alpha\rangle|\beta\rangle
=\sum_{n\ge0}\frac{(\alpha\beta)^n}{n!}|n,n\rangle.
\]

Der Gramkern ist deshalb

\[
K_{\mathrm{coh}}(t,u)
=\sum_{n\ge0}\frac{(\bar t u)^n}{(n!)^2}
=I_0\!\left(2\sqrt{\bar t u}\right),
\qquad t=\alpha\beta.
\]

Das ist ein Besselkern und kein Szegőkern.

### Nachträglich konditionierter thermischer Produktzustand

Nach der Konditionierung auf \(D=0\) sind die Wahrscheinlichkeiten proportional zu

\[
q^{2n}.
\]

Das Verhältnis aufeinanderfolgender Paarwahrscheinlichkeiten ist also \(q^2=e^{-8\pi}\). Der Zustand bleibt eine diagonale Mischung.

### Echte Tomita Paarung

Die thermische Standarddarstellung eines einzelnen Modus hat dagegen den reinen Vektor

\[
|\Omega_{\mathrm{TFD}}\rangle
=\sqrt{1-q}\sum_{n\ge0}q^{n/2}|n\rangle_L|n\rangle_R.
\]

Er liegt exakt in \(D=0\), besitzt die Hardy Koeffizienten \(z^n\) mit \(z=\sqrt q=e^{-2\pi}\) und den Szegőkern

\[
K_{\mathrm H}(z,w)=\frac1{1-\bar z w}.
\]

Die beiden Besetzungszahlen sind hier perfekt korreliert. Im unabhängigen Quarterprodukt ist ihre Kovarianz null. Tomita links und rechts sind eine Darstellung desselben thermischen Modus. Die Quarter \(1/4\) und \(3/4\) sind zwei physische Blöcke der Z58 Quelle. Ihre direkte Identifikation mit Tomita links und rechts wäre eine zusätzliche Abbildung.

## Folgt die P1 Skala?

Nein. Der Polarshift, seine Differenzzahlsektoren und der Index minus eins bleiben bei jeder positiven Modenenergie gleich. Die thermischen Gewichte ändern sich dagegen mit \(q=e^{-\beta\varepsilon}\). Die Paarleiter bestimmt daher weder \(\beta\) noch \(\varepsilon\).

Nach Einsetzen des bereits aus Z58 übernommenen Werts \(\beta=4\pi\) gilt numerisch

\[
\frac1{2\beta}=\frac1{8\pi}=c_3.
\]

Der Faktor \(1/(2\beta)\) wird von der Paarleiter aber nicht erzeugt. Ohne eine ursprüngliche Gleichung, die genau diese Umrechnung mit der P1 Randnormierung verbindet, ist die Gleichheit nur eine Rückschreibung des schon eingesetzten \(\beta\).

## Woher das Spin 2 Wörterbuch wirklich kommt

Jede der beiden Stromanregungen hat chirales Gewicht eins. Das gemeinsame Wort \(AB\) hat deshalb Gewicht zwei. Diese Addition ist noch kein Spin 2 Wörterbuch.

Der chirale Stresstensor entsteht in Z58 erst durch die Sugawara Identität. Für zwei abelsche Ströme ist sein Niveau zwei Zustand proportional zur diagonalen Kombination

\[
L_{-2}|0\rangle
\propto
\frac12\left[
(J^{(1/4)}_{-1})^2
+(J^{(3/4)}_{-1})^2
\right]|0\rangle.
\]

Das hier gerechnete Kreuzwort

\[
J^{(1/4)}_{-1}J^{(3/4)}_{-1}|0\rangle
\]

ist davon linear unabhängig. Sein Label Gewicht zwei stammt allein aus \([L_0,J_{-1}]=J_{-1}\). Ein vierdimensionales Spin 2 Teilchen bräuchte zusätzlich eine hergeleitete Lorentzwirkung, einen Tensorindex, die richtige Charakterfolge und eine physische positive Zeit. Z75 hat bereits gezeigt, dass der Virasoro Charakter mit dem vierdimensionalen Spin 2 Charakter auf keiner Stufe übereinstimmt. Die neue Rechnung schließt nun die einfachere Hoffnung, dass zwei unabhängige Quarteranregungen durch gemeinsames Einschalten automatisch den fehlenden Zustand erzeugen.

## Kontrollen

Der direkte endliche Spurtest mit zwei, vier und acht Besetzungsstufen konvergiert auf den analytischen Kernelwert. Eine zweite endliche Summe bestätigt sowohl \(\omega(P_0)\) als auch \(\omega(P_\partial)\). Der Besselgram, der Hardykern, das Tomita Folgengewicht und eine geänderte Modenenergie dienen als getrennte Kontrollen.

Als gezielter Gegenfall bleibt jede Einzelantwort unverändert, während die gemeinsame Antwort künstlich mit \(1+\varepsilon\), \(\varepsilon=0{,}1\), multipliziert wird. Dann ist

\[
\log\frac{R_{11}R_{00}}{R_{10}R_{01}}
=\log(1{,}1)
=0{,}09531017980432486.
\]

Der Test erkennt eine echte gemeinsame Kopplung also sofort. Ihr Fehlen in der Z58 Quelle ist ein Resultat der Quasifreiheit und Blockfaktorisierung, kein Mangel des Detektors.

Ausführung:

```text
PYTHONDONTWRITEBYTECODE=1 python3 thermische_quelle.py
Alle Prüfungen bestanden: True (18/18)
```

Reproduktionsdatei: [thermische_quelle.py](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/exploration/2026-11-26-unabhaengige-ursprungsprinzipien/thermische_quelle.py:1)

## Stärkster Schluss und nächste Hypothese

Der stärkste Schluss lautet: Die originale quasifreie Z58 Quelle trägt zwei unabhängige Quarterstromschalter und einen positiven gemeinsamen Niveau zwei Zustand, aber keine gemeinsame Dynamik zwischen ihnen. Die Paarleiter besitzt auf dem neutralen Differenzzahlsektor eine echte Hardy Struktur mit Index minus eins. Die thermische Quarterquelle wählt diesen Sektor jedoch nicht aus. Positivität, gleiche Quelle, Niveau zwei und vertauschbare Reihenfolge reichen zusammen nicht für Spin 2. Das ist ein Ergebnis über diese beiden Quartermoden, keine Widerlegung aller möglichen Randeinfügungen oder der OS Rekonstruktion.

Die erste echte Lücke liegt vor jeder Gravitationserklärung. Z58 setzt die quasifreie blockdiagonale Quarterquelle voraus. Die engere ursprüngliche Annulusquelle besitzt zusätzlich eine gemeinsame Zweirandkohärenz. Die gezielte nächste Hypothese lautet deshalb:

> Die zwei physisch schaltbaren Anregungen sind Zweirandströme des vollen Annulusprojektors und nicht die zwei faktorisierten Quarterströme. Seine feste Kreuzkovarianz erzeugt ohne angepassten Parameter eine nichtverschwindende verbundene OS Vierpunktantwort auf Niveau zwei.

### Gezielte Prüfung der schon vorhandenen Randabbildungen

Der Atlas und die Originaldateien enthalten mehr als nur zwei abstrakte Seiten:

1. Der geschlossene Annulusgraph besitzt den festen Zweirandprojektor

   \[
   P_A=\begin{pmatrix}1-C_A&X_A\\X_A&C_A\end{pmatrix},
   \qquad X_A^2=C_A(1-C_A).
   \]

2. Der alte signierte Stromgraph enthält tatsächliche E8 Strombilder auf \(R_{\mathrm{oben}}\), NS und \(R_{\mathrm{unten}}\) in beiden Quartern. Sein Geltungsbereich ist ausdrücklich der alte signierte Quellgraph, nicht der volle physische P1 Rand.

3. Z76 identifiziert die Quarter Ladungsdifferenz mit der modularen Vierteluhr und rekonstruiert einen zustandserhaltenden Zweipfadkanal. Z76 liefert keine Abbildung der beiden bosonisierten Quartermoden auf die beiden geometrischen Annulusränder.

4. Ein lokales Weylpaar mit dem richtigen Level eins Kokykel ist konstruiert. Sein tatsächliches P1 Bild bleibt bedingt.

Gelesene Originalstellen: [voller Annulusprojektor](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/native-arc-charge-source-interface-20261026/PROJECTOR-STATE-REVIEW.txt:1), [beide alten Randflächen und Quarterstrombilder](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/chat-triangulation-20261108/TFPT-P1-SPINWEG-UND-NAHT-20261110.txt:63), [lokaler P1 Weyltest](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/original-admissible-source-closure-20261025/P1-W8-LOCAL-MAP-AND-TWOPOINT-TEST.txt:27) und [Z76 Quarteruhr und Zweipfadkanal](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/exploration/2026-11-23-phasenraum-f5-holocode/z76_zweipfadkanal.py:1).

Damit fehlt nicht bloß eine Formel für den Weyloperator. Es fehlt die gemeinsame ursprüngliche Abbildung, die gleichzeitig

\[
(J^{(1/4)}_{-1},J^{(3/4)}_{-1})
\longmapsto
(R_{\mathrm{oben}},R_{\mathrm{unten}})
\]

mit Zustand, Adjunktion, Kreuzkovarianz, Reflexion, Zeit und Support verträglich macht. Genau diese Abbildung würde entscheiden, ob \(D=N_A-N_B\) eine native Randbedingung ist oder nur die nachträgliche Projektion dieser Rechnung.

Der kleinste Test muss deshalb nicht den Annulusprojektor neu erfinden. Er muss die vorhandenen Randstrombilder durch den vorhandenen lokalen Weyltest in denselben positiven Zweirandzustand heben und prüfen, ob dessen Kreuzblock exakt das Tomita Matching \(n_A=n_B\) liefert. Danach wird dieselbe Vierertabelle mit dem vollen Projektorkern gerechnet. Bleibt der verbundene Anteil null oder scheitert die gemeinsame Abbildung, ist diese Quellenroute beendet. Ist er ungleich null, folgt noch kein Graviton. Dann muss sein Niveau zwei Eigenraum unter einer unabhängig hergeleiteten Lorentzwirkung den vierdimensionalen Spin 2 Charakter tragen.

Atlasbezug, nur gelesen: `finding:h6h31-reconstruction-modular-generator-gaussian-currents-ladder-20261124`, `finding:h6h44-native-reflection-positivity-os-gives-clock-for-thermal-original-tprime-pointwise-20261125`, `finding:h6h39-nonabelian-multipoint-no-scalar-4d-no-graviton-from-virasoro-20261124`, `finding:actual-two-face-projector-reduces-free-source-test-to-cross-coherence`, `finding:tfpt-signed-triplet-two-faces-20261110-same-cut-old-two-modes`, `finding:local-current-phase-and-conditional-cyclic-p1-w8-image`. Der Atlas wurde nicht verändert.
