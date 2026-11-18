# Gemeinsamer Kreis-, Mark- und Sewing-Test des bedingten Filters

## Prüfhypothese

Für diesen Test werden die vier Markpunkte explizit als \(z_k=e^{i\theta_k}\), \(\theta_k=k\pi/2\), identifiziert und die Sewing-Zeit gleich dem positiven Parameter \(\alpha\) gesetzt. Das ist eine zusätzliche Suchhypothese. Der gefundene Modenfilter hat \(\kappa_ns_n^2=1/n\) für \(n\ge2\).

## Kreisantwort und gemeinsamer Referenzabzug

Für \(0<\theta<2\pi\) und \(0<r<1\) gilt die absolut konvergente Abel-Summe

\[
\sum_{n\ge1}\frac{r^n\cos(n\theta)}n
=-\frac12\log(1-2r\cos\theta+r^2).
\]

Der Grenzwert \(r\uparrow1\) ist \(-\log(2|\sin(\theta/2)|)\); zieht man den fehlenden \(n=1\)-Term \(\cos\theta\) ab, folgt

\[
K_R(\theta)=\sum_{n\ge2}\frac{\cos(n\theta)}n
=-\log\!\left(2\left|\sin\frac\theta2\right|\right)-\cos\theta.
\]

Am Diagonalpunkt ist die Reihe divergent. Verwende für den Matrixvergleich denselben harmonischen Counter \(H_N=\sum_{n=1}^N1/n\), aber nur eintragsweise auf der Diagonale:

\[
K_{\rm full,N}[i,j]=\sum_{n=1}^N\frac{\cos(n(\theta_i-\theta_j))}n-H_N\delta_{ij},
\quad
K_{R,N}[i,j]=\sum_{n=2}^N\frac{\cos(n(\theta_i-\theta_j))}n-H_N\delta_{ij}.
\]

Der Abel-Grenzwert bestimmt die Off-Diagonale; der Counter setzt die volle Diagonale auf 0 und die gefilterte Diagonale auf \(-1\), denn der \(n=1\)-Term wurde entfernt, ohne den Counter zu ändern. Auf allen vier Punkten gilt exakt

\[
K_{\rm full}=-\frac{\log2}{2}C,
\qquad K_R=K_{\rm full}-C_{\cos},
\]

mit

\[
C=\operatorname{circ}(0,1,2,1),\qquad
C_{\cos,ij}=\cos(\theta_i-\theta_j).
\]

Also hat \(K_R\) bei den vier Punkten Diagonale \(-1\), Nachbarwert \(-\log2/2\) und Gegenpunktwert \(1-\log2\). Es wird nicht als positive endliche Kovarianz ausgegeben: die Eigenwerte der gemeinsamen-counter Matrix sind \(-\log4\), \(\log2-2\) (zweifach) und 0, somit nichtpositiv.

Sei \(T\) der zyklische Vierteldrehungsoperator und \(E_+=(I+T^2)/2\). \(C_{\cos}\) trägt nur die unter der Halbdrehung ungeraden Fouriermoden, daher

\[
E_+C_{\cos}E_+=0,
\qquad E_+K_RE_+=E_+K_{\rm full}E_+.
\]

Mit der Originalquelle \(G_{\rm orbit}=-c_3\log2\,C\) und \(K_{\rm full}=-(\log2/2)C\) folgt der exakte gerade-Sektor-Vergleich

\[
2c_3 E_+K_RE_+=E_+G_{\rm orbit}E_+.
\]

Das ist ein echter, aber bedingter Kreis-/Mark-Treffer. Der harmonische Diagonal-Counter ist eine für diesen Vergleich gesetzte Referenz und nicht aus P1 hergeleitet. Wenn man stattdessen den gefilterten Kern separat auf Diagonale 0 renormiert, zieht man \(H_N-1\) statt \(H_N\) ab und addiert dadurch \(I\): \(K_R^{(\mathrm{diag0})}=K_R+I\). Dann verschiebt sich die gerade Projektion um \(E_+\), und der obige Match gilt nicht mehr. Diese Referenzwahl bleibt deshalb explizit.

## Sewing-Zeit: zwei verschiedene Mark-Charaktere

Ein positiver Sewing-Operator gewichtet einen Zylinder-Modengrad \(n\) mit \(L_0=n\) und damit mit \(e^{-\alpha n}\). Eine positive Summe mit mindestens zwei verschiedenen Graden kann nicht für alle \(\alpha>0\) eine feste Konstante mal \(e^{-2\alpha}\) sein.

**Skalarer konstanter Markport.** Der direkte Port \(p_0=(1,1,1,1)/2\) mittelt die vier Marken ohne Zusatzphase und selektiert \(n=0\pmod4\). Für die gefilterten Gewichte \(1/n\) ergibt sich

\[
4\sum_{m\ge1}\frac{e^{-4m\alpha}}{4m}
=\sum_{m\ge1}\frac{e^{-4m\alpha}}m
=-\log(1-e^{-4\alpha}).
\]

Die ersten Terme sind \(e^{-4\alpha}+\tfrac12e^{-8\alpha}+\cdots\). Es gibt gar keinen Grad-2-Anteil; bei Cutoff \(N\ge8\) tragen bereits die positiven Grade 4 und 8 bei.

**Gewicht-2-Feld auf der Ebene.** Schreibe \(R(z)=\sum_nR_nz^{-n-2}\). Der Erzeugungskoeffizient \(R_{-n}\) erscheint als \(z^{n-2}\). Der konstante Viermarkenmittelwert ist deshalb

\[
\frac14\sum_{k=0}^3 i^{k(n-2)}=
\begin{cases}1,&n=2\pmod4,\\0,&\text{sonst.}\end{cases}
\]

Denn \(i^{k(n-2)}=i^{-2k}i^{kn}\): der Faktor \(i^{-2k}=(-1)^k\) ist durch das Feldgewicht 2 festgelegt. Er ist die Twist-/Charakteränderung gegenüber dem skalaren Greenport, keine frei wählbare neue Projektion. Mit den gefilterten Gewichten entsteht

\[
\sum_{m\ge0}\frac{e^{-(4m+2)\alpha}}{4m+2}
=\frac12e^{-2\alpha}+\frac16e^{-6\alpha}+\frac1{10}e^{-10\alpha}+\cdots.
\]

Für \(\alpha>0\) lässt sich die Reihe auch als \(\tfrac12\operatorname{artanh}(e^{-2\alpha})=\tfrac14\log\!\frac{1+e^{-2\alpha}}{1-e^{-2\alpha}}\) schreiben. Der positive Grad-6-Term bleibt im vollen Turm erhalten. Auch dieser Port ist daher keine feste Konstante mal \(e^{-2\alpha}\). Ohne Filter sind die positiven Feldnormgewichte \(\kappa_{4m+2}=8(4m+2)((4m+2)^2-1)\); der ungefilterte R-Normturm enthält ebenso mehrere Grade.

Die skalare \(P_0\)-Auslesung auf \(0\pmod4\) und die Gewicht-2-Feldauslesung auf \(2\pmod4\) sind unterschiedliche Abbildungen. Man darf nicht dieselbe unverschobene \(P_0\)-Regel zugleich auf den skalaren Greenkern und das Feld \(R(z)\) anwenden, ohne die Gewicht-2-Phase mitzunehmen.

## Geltungsbereich

Der positive Kontrollpfad trennt den literal-Tower-Test von einer gemeinsamen Multiport-Hypothese. Der einzelne Quellenzustand \(R_{-2}\Omega\) hat \(L_0=2\) und Normquadrat 48; wählt man explizit seinen Rang-eins-Projektor \(P_{L_0=2}\), dann gilt

\[
\left\|e^{-\alpha L_0/2}R_{-2}\Omega\right\|^2=48e^{-2\alpha}.
\]

Daneben besteht der separate gefilterte skalare Kreis-/Markport den bedingten even-Sektor-Match mit dem gemeinsamen Counter. Diese zwei Ports sind logisch getrennt: Der positive Ein-Zustand-Kontrollpfad rettet nicht die pure-Graduierung des vollständigen Towers unter derselben festen Projektion.

Der negative Sewing-Befund gilt nur für die wörtliche Auswertung des vollständigen positiven Towers mit fester Markprojektion und \(\alpha L_0\)-Sewing. Er schließt keinen eigens definierten Kontaktport, keine andere Quelle oder Referenzkontraktion und keine eigenständige Feld-zu-Skalar-Abbildung aus. Eine Rand-Portabbildung muss kein Operatorhomomorphismus sein; eine solche Forderung wird hier nicht als allgemeines Gate verwendet. Ebenso ist \(\alpha\) hier nur die getestete positive Sewing-Zeit, nicht die daraus physisch hergeleitete TFPT-Kopplung.

**Verdict:** Der inverse Casimirfilter hat einen exakten bedingten geraden Kreis-/Mark-Match, wenn beide Kerne denselben harmonischen Diagonal-Counter verwenden. Für die wörtliche Ein-Port-Auswertung des vollen positiven Sewing-Turms entstehen unter beiden durch die Gewichte festgelegten Charakteren mehrere Grade und damit kein reiner Grad-2-Kern. Ein separater Rang-eins- \(L_0=2\)-Port liefert hingegen erwartungsgemäß \(48e^{-2\alpha}\), sodass kein allgemeiner Multiport-No-Go folgt. Auswahl von \(P_{L_0=2}\), Filter, Kopplung und Counter aus P1 bleibt offen.

Das reproduzierbare Skript [source_filter_joint_gate.py](source_filter_joint_gate.py) verwendet SymPy für exakte Matrix- und Logarithmusidentitäten. Seine Ausgabe liegt in [results/ports.json](results/ports.json):

```sh
python3 source_filter_joint_gate.py --outputpath results/ports.json
```
