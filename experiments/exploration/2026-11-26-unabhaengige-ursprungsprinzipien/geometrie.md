# Unabhängiger Ursprung aus einem reellen Quadrikenbüschel

## Ergebnis

Der stärkste neue Ansatz ist ein reelles Grad 4 del Pezzo Büschel. Er beginnt bei einer glatten komplexen Fläche als kleinstem vollständigen Schnitt zweier Quadriken. Er setzt keinen TFPT Randkern, keine vier Marken, keine Fünferzahl, keine Teilung 3+2, kein $D_5$, keine 16 und kein $E_8$ ein.

Aus diesem einen Grundobjekt folgen zunächst drei bekannte TFPT Strukturen gemeinsam:

1. Das Büschel besitzt fünf singuläre Quadriken.
2. Das zum kanonischen Divisor orthogonale Picardgitter ist $D_5$.
3. Die Fläche besitzt 16 Geraden. Sie tragen die Geometrie der 16 Halbspinorgewichte von $D_5$.

Eine reelle Form mit Signatur $(4,1)$ schärft das Ergebnis. Hat der zugehörige selbstadjungierte Operator tatsächlich ein nichtreelles Eigenwertpaar, dann erlaubt der Witt Index kein zweites Paar. Von den fünf Spektralmarken sind daher genau drei reell und zwei komplex konjugiert. Das liefert eine operatorinvariante Zerlegung

\[
\mathbb R^5=\mathbb R^3\oplus\mathbb R^2.
\]

Nach Komplexifizierung ist dies eine Darstellung $\mathbb C^3\oplus\mathbb C^2$, wobei die reelle Struktur die beiden Spektralbasen im Zweierblock vertauscht. Damit erzeugt die Geometrie erstmals den Strukturverbund $5$, $3+2$, $D_5$ und $16$ ohne eine TFPT Saat.

Das Urteil ist **Kandidat Stufe 1**. Die Rechnung erzeugt einen echten P2 Strukturkandidaten, aber noch nicht P2 als physische Schnittstelle. Weder die Signatur $(4,1)$ noch das Auftreten des nichtreellen Paars noch die Eichwirkung folgen bisher aus dem Grad 4 del Pezzo Prinzip.

## Warum gerade dieser Ausgangspunkt

Zwei unabhängige Quadriken in $\mathbb P^4$ schneiden sich in einer komplexen Fläche. Durch Adjunktion gilt

\[
K_X=\mathcal O_X(-5+2+2)=\mathcal O_X(-1),
\]

also ist $X$ eine del Pezzo Fläche vom Grad $H^2=4$. Die erste gewählte Struktur ist damit nicht die Zahl fünf, sondern ein glatter Flächenschnitt aus genau zwei quadratischen Gleichungen. Die fünf kommt erst aus der Determinante des $5\times5$ Büschels.

Für zwei symmetrische Matrizen $J,A$ ist die Diskriminante

\[
d(\lambda)=\det(A-\lambda J)
\]

eine Quintik. Bei glattem Basisschnitt hat sie fünf einfache Nullstellen. Jede Nullstelle bezeichnet eine singuläre Quadrik vom Rang 4. Die Primärliteratur beschreibt dieselbe Struktur als Grad 5 Spektralüberlagerung. Ihr voller geometrischer Symmetrierahmen ist $W(D_5)$, mit der natürlichen Wirkung auf den fünf singulären Quadriken.

## Der scharfe reelle Test

Der Test verwendet

\[
J=\operatorname{diag}(1,1,1,1,-1).
\]

Damit hat $J$ Signatur $(4,1)$. Für eine reelle symmetrische Matrix $A$ ist

\[
T=J^{-1}A=JA
\]

$J$ selbstadjungiert, denn $T^T J=JT=A$.

### Beispiel mit $3+2$

\[
A_{3+2}=
\begin{pmatrix}
1&0&0&0&0\\
0&2&0&0&0\\
0&0&3&0&0\\
0&0&0&0&1\\
0&0&0&1&0
\end{pmatrix}.
\]

Auch $A_{3+2}$ hat Signatur $(4,1)$. Exakt ergibt sich

\[
\det(A_{3+2}-\lambda J)
=(\lambda-1)(\lambda-2)(\lambda-3)(\lambda^2+1).
\]

Die letzten beiden reellen Koordinaten tragen

\[
T_{\rm Paar}=\begin{pmatrix}0&1\\-1&0\end{pmatrix},
\qquad T_{\rm Paar}^2=-I.
\]

Auf dieser Ebene hat $J$ Signatur $(1,1)$. Sie ist die neutrale Ebene des konjugierten Paars $+i,-i$. Die übrigen drei Eigenlinien besitzen die reellen Eigenwerte $1,2,3$. Das ist die verlangte $3+2$ Zerlegung als tatsächliche Operatorwirkung.

Alle fünf Nullstellen sind einfach. Für jede Nullstelle hat $A_{3+2}-\lambda J$ Rang 4, und ihr eindimensionaler Kern ist nicht isotrop für die komplex bilineare Form $x^T Jx$. Der Vertex keiner singulären Quadrik liegt daher im gemeinsamen Basisschnitt. Der komplexe del Pezzo Schnitt ist glatt.

### Gleich starke Gegenkontrolle ohne $3+2$

\[
A_{5+0}=\operatorname{diag}(1,2,3,4,-5)
\]

hat dieselbe Signatur $(4,1)$ wie $J$ und $A_{3+2}$. Trotzdem gilt

\[
\det(A_{5+0}-\lambda J)=
(\lambda-1)(\lambda-2)(\lambda-3)(\lambda-4)(\lambda-5).
\]

Hier sind alle fünf Marken reell. Auch dieser komplexe Basisschnitt ist glatt. Die Signatur erlaubt das konjugierte Paar, sie wählt es nicht.

### Warum höchstens ein Paar möglich ist

Seien $Tv=\lambda v$ und $Tw=\mu w$. Aus der $J$ Selbstadjungiertheit folgt für die hermitesche Paarung

\[
(\lambda-\bar\mu)\,v^*Jw=0.
\]

Für einen nichtreellen Eigenwert ist deshalb $v^*Jv=0$. Die reelle Ebene aus Realteil und Imaginärteil von $v$ ist nicht entartet und hat je eine positive und eine negative Richtung. Jedes konjugierte Paar verbraucht damit eine negative Richtung. Bei Signatur $(4,1)$ gibt es nur eine. Deshalb kann ein reguläres Büschel höchstens zwei nichtreelle Eigenwerte besitzen. Da die Quintik reell ist, erzwingt ein tatsächlich vorhandenes Paar genau drei weitere reelle Wurzeln.

Nakatsukasa und Noferini beweisen die entsprechende scharfe Inertiaschranke allgemein: Hat die zweite Form nur eine negative Richtung, sind bei $n=5$ mindestens $n-2=3$ verallgemeinerte Eigenwerte reell. Die beiden exakten Matrizen oben zeigen beide erlaubten Fälle.

## Welche P2 Darstellung entsteht

Die Rechnung erzeugt zunächst die Spektraldarstellung der fünf singulären Quadriken. Im Paarzweig zerfällt sie in den von den drei reellen Eigenlinien getragenen Dreierblock und die reelle Ebene des konjugierten Paars. Nach Komplexifizierung lautet der Blocktyp

\[
V_5=V_3\oplus V_2,
\qquad \dim_{\mathbb C}(V_3,V_2)=(3,2).
\]

Die komplexe del Pezzo Geometrie liefert außerdem das $D_5$ Wurzelgitter und 16 Geraden. Im Standardmodell des $D_5$ Wurzelsystems gibt es 40 Wurzeln und 16 Gewichte in einer Halbspinordarstellung. Der Test reproduziert beide Zähler exakt. Unter dem $S_5$ Anteil von $W(D_5)$ zerfallen die 16 Halbspinorgewichte nach gerader Minuszahl in Orbits der Größen $1,10,5$.

Das ist noch keine Darstellung der physischen Eichgruppe. Aus der Blockzerlegung folgen weder lokale $SU(3)$ und $SU(2)$ Verbindungen noch die Hyperladungen

\[
(-1/3,-1/3,-1/3,1/2,1/2).
\]

Auch die Identifikation der 16 Geraden mit $\Lambda^{\rm gerade}\mathbb C^5$ als physischem Fermionraum braucht eine äquivariante Operatorabbildung. Der Kandidat liefert den geometrischen Träger, nicht seine Eichwirkung.

## Negativkontrolle ohne reelle Paarstruktur

Über den komplexen Zahlen besitzt ein unmarkiertes Grad 4 del Pezzo Büschel nur eine ungeordnete Fünfermenge. Die volle $S_5$ Wirkung hat auf den zehn möglichen $3+2$ Partitionen einen einzigen Orbit und keine feste Partition. Eine gewählte Zweiermenge hat den Stabilisator $S_2\times S_3$ von Ordnung 12 und damit Index 10.

Eine faktorisierte binäre Quintik wie

\[
(u^3-2v^3)(u^2-3v^2)
\]

trägt zwar sofort eine Teilung in Grad 3 und Grad 2. Diese Teilung wurde dann mit der Faktorisierung eingesetzt. Der reelle Paarzweig ist stärker, weil die Teilung aus der Spektralart des $J$ selbstadjungierten Operators folgt. Die Gegenkontrolle mit $A_{5+0}$ zeigt zugleich, dass diese Spektralart eine zusätzliche, noch nicht hergeleitete Phase ist.

## Abgrenzung zur vorhandenen Singularitätsroute

Der Atlas enthält die Brieskorn Kette

\[
x^2+y^3+z^5=0
\longrightarrow \mu=8
\longrightarrow E_8
\longrightarrow h=30
\]

bereits als `TOPO.BRIESKORN.01`. Sie erzeugt unter gesetzten Exponenten das Milnorgitter und die Coxeteruhr, aber keinen neuen Ursprung der Exponenten. Der del Pezzo Test wiederholt diese Route nicht.

Zusätzlich wurde ein 4D $N=2$ Wörterbuch geprüft, weil es eine standardisierte physikalische Abbildung von einer isolierten Calabi Yau Dreifachsingularität auf Coulomboperatoren, Massenparameter und das BPS Ladungsgitter bereitstellt. Nur unter der ausdrücklich zusätzlichen Identifikation

\[
5=r_{\rm Coulomb},\qquad 3=f_{\rm Flavor}
\]

müsste die Milnorzahl $\mu=2r+f=13$ sein. Die suspendierte $E_8$ Singularität

\[
x_0^2+x_1^2+x_2^3+x_3^5=0
\]

hat jedoch $\mu=8$ und das exakte Deformationsspektrum

\[
\left\{\frac{15}{8},\frac32,\frac54,\frac98,
\frac78,\frac34,\frac12,\frac18\right\}.
\]

Damit ist $(r,f)=(4,0)$. Die einzige diagonale Brieskorn Pham Form mit $\mu=13$ ist bis auf Permutation $(2,2,2,14)$; sie liefert $(r,f)=(6,1)$. Das verwirft genau dieses Wörterbuch. Es verwirft keine allgemeine Singularitätsgeometrie, denn TFPT Trägerplätze sind formal keine Coulomboperatoren und drei Generationen sind formal kein $N=2$ Flavourrang.

## Stärkster Einwand und nächster entscheidender Test

Der stärkste Einwand sitzt jetzt an einer klaren Stelle. Grad 4 del Pezzo erzeugt $5$, $D_5$ und $16$. Die reelle Signatur $(4,1)$ begrenzt die Spektraltypen. Erst das tatsächliche nichtreelle Paar erzeugt $3+2$. Beide Büschel im Test haben dieselben Signaturen und sind komplex glatt, aber nur eines besitzt das Paar. Die Paarphase ist daher unberechnet, nicht durch Grad 4 ausgewählt.

Der kleinste nächste Test braucht keine weitere Geometriewelt. Gesucht ist ein intrinsisches Funktional oder eine Stabilitätsbedingung auf reellen Grad 4 del Pezzo Büscheln, deren Minima genau im einfachen Einpaarzweig liegen. Erfolg wäre eine begründete Auswahl des offenen Spektralgebiets mit drei reellen Wurzeln und einem konjugierten Paar. Abbruchkriterium ist ein gleichwertiges glattes Minimum im Fünfreellzweig. Erst danach lohnt die Prüfung, ob dieselbe Auswahl eine natürliche hermitesche Wirkung auf $V_3\oplus V_2$ erzeugt und ob deren Stabilisator die physische $3+2$ Eichwirkung trägt.

## Primärquellen

* Brendan Hassett, Andrew Kresch und Yuri Tschinkel, [On the Moduli of Degree 4 Del Pezzo Surfaces](https://math.nyu.edu/~tschinke/papers/yuri/13noheighttwo/noheighttwo16.pdf). Die Arbeit beschreibt den Schnitt zweier Quadriken, die binäre Quintik, 16 Geraden, das primitive $D_5$ Picardgitter sowie die Grad 5 Spektralüberlagerung und ihre $W(D_5)$ Wirkung.
* B. È. Kunyavskii, A. N. Skorobogatov und M. A. Tsfasman, [Del Pezzo Surfaces of Degree Four](https://www.numdam.org/item/MSMF_1989_2_37__1_0/). Abschnitt 3 beweist die fünf singulären Quadriken, den Rang 4 jedes singulären Mitglieds, das $D_5$ Gitter und die Konfiguration der 16 Geraden.
* Yuji Nakatsukasa und Vanni Noferini, [Inertia Laws and Localization of Real Eigenvalues for Generalized Indefinite Eigenvalue Problems](https://arxiv.org/abs/1711.00495). Die scharfe Inertiaschranke gibt bei nur einer negativen Richtung mindestens drei reelle Eigenwerte und höchstens ein nichtreelles Paar.
* Dan Xie und Shing Tung Yau, [4d N=2 SCFT and Singularity Theory Part I: Classification](https://arxiv.org/abs/1510.01324). Gleichungen 2.8 und 2.16 geben die Skalendimensionen der Miniversaldeformation und $\mu=2r+f$; der Schnitt der verschwindenden Zyklen ist der BPS Quiver Kandidat.

## Ausführung

Ausgeführt wurde ausschließlich

```text
python3 -B experiments/exploration/2026-11-26-unabhaengige-ursprungsprinzipien/geometrietest.py
```

Der Lauf bestand alle exakten Assertions. Atlas, Ledger, Plan, Papers und bestehende Experimente wurden nicht verändert.
