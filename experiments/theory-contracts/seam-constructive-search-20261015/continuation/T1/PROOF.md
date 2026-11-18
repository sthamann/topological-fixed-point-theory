# Fortsetzung: der echte Einzelwindungsschritt \(T_1\)

## Ergebnis

Der ursprüngliche Einzelwindungsschritt lässt sich exakt auf die gekreuzte Defektdarstellung anheben. Er wählt mit den bisher festen Originaldaten aber **weder** die Orientierung \(\chi=\pm i\) **noch** den (P)-Liftzweig \(\kappa=\pm1\) aus.

Eine scheinbare Auswahl \(\kappa=+1\) entsteht nur unter der zusätzlichen Forderung

\[
 \widehat C X\widehat C^{-1}=\nu X,
\]

also wenn der adjoinierte Halbsektorgenerator (X) ein skalares Eigenfeld des Kanalzyklus sein soll. Die ursprüngliche ganzzahlige CAR-Quelle enthält (X) nicht und begründet diese Forderung nicht.

Ohne sie existiert die positive Familie

\[
 \widehat C_r|n,c\rangle
 =\gamma(n)\zeta_8^{rc}|Cn,c\rangle,
 \qquad r\in\mathbb Z/8,
\]

mit Defektgrad \(c\in\{0,1\}\). Gerade (r) geben \(\kappa=+1\), ungerade (r) geben \(\kappa=-1\). Alle besitzen dieselbe ursprüngliche Integer-CAR-Wirkung, dieselbe Ladungswirkung und dieselbe (K_\alpha)-Wardidentität.

Das genaue Verdict ist

```text
CONDITIONAL_T1_LIFT_FAMILY_NO_KAPPA_SELECTION
```

## 1. Feste Originaldaten

Aus `seam-joint-ward-state-20261014/PROOF.txt`, Zeilen 752–861:

\[
 T_1(p_0,\ldots,p_7)=(p_7+1,p_0,\ldots,p_6),
\]

\[
 T_1^4(p)=Pp+a,\qquad T_1^8(p)=p+q,
\]

mit

\[
 q=(1^8),\qquad a=(1,1,1,1,0,0,0,0).
\]

Die feste Zeit ist

\[
 K_\alpha(p,N)=N+\frac{|p|^2}{2}+\delta_\alpha\cdot p,
 \qquad (\delta_\alpha)_j=\frac{j+\alpha}{8}-\frac12,
\]

für \(\alpha=1/4,3/4\). (W_1) ist ungerade und trägt Ladung eins.

Aus Zeilen 1003–1087 derselben Datei:

\[
 \epsilon_{\mathbb Z}(n,m)=(-1)^{\sum_{i>j}n_i m_j},
\]

\[
 P_{\rm CAR}|n\rangle=\eta(n)|Pn\rangle,
 \qquad \eta(n)=(-1)^{A(n)B(n)},
\]

und

\[
 W_4=Z_aE_aP_{\rm CAR}.
\]

Diese Originalgrößen werden nicht verändert.

## 2. Der source-kanonische Achtzyklus auf der Integer-CAR-Algebra

Setze

\[
 C(n_0,n_1,\ldots,n_7)=(n_7,n_0,\ldots,n_6),
 \qquad u=e_0.
\]

Dann ist (T_1(p)=Cp+u\). Der geordnete Cocycle benötigt unter (C) die Korrektur

\[
 \gamma(n)=(-1)^{n_7(n_0+\cdots+n_6)}.
\]

Definiere

\[
 \widehat C E_n\widehat C^{-1}=\gamma(n)E_{Cn},
 \qquad \widehat C|0\rangle=|0\rangle.
\]

Für alle (n,m\in\mathbb Z^8) gilt exakt

\[
 \gamma(n)\gamma(m)\epsilon(Cn,Cm)
 =\epsilon(n,m)\gamma(n+m).
\]

Außerdem ist \(\gamma(e_j)=1\) für alle Elementarladungen. Daher führt \(\widehat C\) die tatsächliche phasenfreie Kanalrotation

\[
 \widehat C E_{e_j}\widehat C^{-1}=E_{e_{j+1\bmod8}}
\]

aus. Ein zusätzlicher Ganzgittercharakter würde bereits diese feste Einteilchenwirkung ändern.

Vier Schritte reproduzieren die vorhandene (P)-Phase:

\[
 \prod_{r=0}^{3}\gamma(C^rn)=\eta(n),
 \qquad \widehat C^4=P_{\rm CAR}
\]

auf der ursprünglichen ganzzahligen Algebra.

## 3. Natürlich normierter Fock-Implementer

Setze

\[
 W_1=-Z_uE_u\widehat C.
\]

In der festen Konvention (Z_u|p\rangle=(-1)^{u\cdot p}|p\rangle) kompensiert das Minus (Z_u|u\rangle=-|u\rangle). Damit

\[
 W_1|0\rangle=|u\rangle.
\]

Der Originaltext bestimmt den Implementer zunächst nur bis auf eine zentrale Phase. Diese Normierung wählt den positiv beschrifteten Einteilchenzustand; alle Konjugationsaussagen sind von der globalen Phase unabhängig.

Auf Ladungszuständen gilt

\[
 W_1|n\rangle
 =(-1)^{n_7(n_0+\cdots+n_6+1)}|Cn+u\rangle.
\]

Viermaliges Anwenden liefert einschließlich aller geordneten Kleinzeichen

\[
 W_1^4=Z_aE_aP_{\rm CAR}=W_4.
\]

Die Einteilchen-Bündelwirkung ist

\[
 B_1(x)=e^{ix/8}C.
\]

Mit \(D_\alpha=\operatorname{diag}(e^{2\pi i(k+\alpha)/8})\) gilt

\[
 D_\alpha C=e^{\pi i/4}CD_\alpha,
\]

also

\[
 B_1(x+2\pi)D_\alpha=D_\alpha B_1(x),
 \qquad B_1(x)^4=e^{ix/2}P=B_4(x).
\]

Die Viertelphase bleibt somit schon beim Einzelwindungsschritt erhalten.

## 4. Allgemeiner sektor-diagonaler Lift

Die gekreuzte Erweiterung besitzt

\[
 XE_nX^{-1}=\chi(n)E_n,
 \qquad X^2=i^zE_q,
\]

mit \(\chi(e_j)=+i\) oder \(-i\) einheitlich für alle Kanäle.

Da (Cs=s\), erhält ein sektor-diagonaler Lift die absolute Ladung. Er darf auf dem zusätzlichen Defektgrad trotzdem eine Phase tragen:

\[
 \widehat C_r|n,c\rangle
 =\gamma(n)\zeta_8^{rc}|Cn,c\rangle,
 \qquad \zeta_8=e^{\pi i/4}.
\]

Auf (c=0\) ist dies für jedes (r) genau der ursprüngliche source-kanonische \(\widehat C\). Auch (T_1(p)=Cp+u\) und die Zeitwirkung sind unabhängig von (r\).

Die Konjugationswirkung auf (X\) hat jedoch zwei Eingangsgrad-Phasen:

\[
 \operatorname{Ad}(\widehat C_r)(X):
 \quad \zeta_8^r\ \text{auf }c=0,
 \qquad \zeta_8^{4-r}\ \text{auf }c=1.
\]

Der Unterschied um den Viererterm folgt aus dem Normalform-Carry (X^2=i^zE_q) und der kombinierten Cocycle-Identität

\[
 \frac{\epsilon(n,q)\gamma(n+q)}
 {\gamma(n)\epsilon(Cn,q)}=-1.
\]

Der Quotient \(\gamma(n+q)/\gamma(n)\) allein ist nicht konstant; erst die vollständige geordnete Normalform liefert das feste Minuszeichen.

Nur bei (r=2\) oder (r=6\) stimmen beide Phasen überein. Ausschließlich in dieser zusätzlichen Unterklasse ist (X) ein skalares Eigenfeld:

\[
 \widehat C X\widehat C^{-1}=\nu X,
 \qquad \nu=+i\ \text{oder}\ -i.
\]

Dann ergibt das Anwenden auf (X^2=i^zE_q) korrekt \(\nu^2=\gamma(q)=-1\). Diese Rechnung beweist aber nur eine Aussage **unter** Skalar-Kovarianz; sie beweist deren Herkunft nicht.

Die gemischte Integerrelation fordert lediglich

\[
 \chi(Cn)=\chi(n),
\]

was für beide uniformen Orientierungen \(\chi=\pm i\) gilt.

## 5. Beide \(\kappa\)-Zweige besitzen einen vierten Ursprung

Für die Familie \(\widehat C_r\) folgt direkt

\[
 \widehat C_r^4|n,c\rangle
 =\eta(n)(-1)^{rc}|Pn,c\rangle,
 \qquad \widehat C_r^8=1.
\]

Damit

\[
 \boxed{\kappa=(-1)^r}.
\]

Gerade (r) liefern \(\kappa=+1\), ungerade (r) liefern \(\kappa=-1\). Ein expliziter Gegenfall gegen eine automatische Auswahl ist

\[
 \widehat C_-|n,c\rangle
 =\gamma(n)\zeta_8^c|Cn,c\rangle.
\]

Er erfüllt

\[
 \widehat C_-^4|n,c\rangle
 =\eta(n)(-1)^c|Pn,c\rangle,
 \qquad \widehat C_-^8=1,
\]

und ist auf (c=0\) exakt der ursprüngliche Integer-CAR-Zyklus. Seine (X\)-Konjugation ist gradingabhängig: \(\zeta_8\) auf Eingangsgrad null und \(-\zeta_8^{-1}=\zeta_8^3\) auf Eingangsgrad eins.

Für den vollständigen Winding-Operator

\[
 W_{1,r}=-Z_uE_u\widehat C_r
\]

kommen die festen Wirkungen

\[
 Z_uXZ_u^{-1}=e^{i\pi u\cdot s}X=iX,
 \qquad E_uXE_u^{-1}=\chi(u)^{-1}X
\]

hinzu. Sie verschieben beide gradingabhängigen Exponenten gleich und ändern weder ihre Differenz noch \(\kappa=(-1)^r\). Somit gilt für beide Zweige

\[
 W_{1,r}^4=W_{4,\kappa}
\]

auf der vollen zweigradigen Darstellung.

Die korrekte Grenze ist damit scharf: \(\kappa=+1\) folgt nur aus der zusätzlich gewählten scalar eigenoperator covariance für (X\). Die Originalquelle gibt keine Wirkung auf (X\) vor.

## 6. Exakte Zeitidentität auf beiden Graden

Für jedes (p\in\mathbb Z^8\cup(\mathbb Z^8+s)\) gilt durch direktes Einsetzen

\[
 \boxed{
 K_\alpha(Cp+e_0,N)-K_\alpha(p,N)
 =\frac{Q(p)+\alpha}{8}}
\]

für \(\alpha=1/4,3/4\). Verwendet werden nur

\[
 \delta_{j+1}-\delta_j=\frac18,
 \qquad \delta_0-\delta_7=-\frac78,
 \qquad \frac{|Cp+e_0|^2-|p|^2}{2}=p_7+\frac12.
\]

Summation entlang (p,T_1p,T_1^2p,T_1^3p) ergibt

\[
 K_\alpha(T_1^4p,N)-K_\alpha(p,N)
 =\frac{Q(p)}2+\frac34+\frac\alpha2,
\]

also die vorhandene (W_4)-Wardidentität. \(\chi\), (r) und \(\kappa\) treten in dieser Energie nicht auf. Die Originalzeit unterscheidet die Liftzweige daher nicht.

## 7. Erste zusätzliche Wahl

Die erste Zusatzwahl ist nicht ein neuer Hamiltonoperator oder Zustand. Sie lautet:

> Soll die Zykluswirkung auf dem adjoinierten (X\) skalar sein, oder darf sie vom Defektgrad abhängen?

Die Skalarwahl reduziert auf (r=2,6\) und damit \(\kappa=+1\). Ohne sie bleiben alle (r\in\mathbb Z/8\) und beide \(\kappa\)-Zweige möglich. Da (X\) nicht zur ursprünglichen Integer-CAR-Quelle gehört, ist diese Wahl dort nicht beantwortet.

Ein Fortschritt zur Auswahl müsste einen bereits quellenbegründeten gemischten Einzelwindungsmarker liefern. Die bekannten Ladungs-, Paritäts-, Zeit- und Viertelidentitäten bleiben für den gradingabhängigen Gegenfall identisch.

## 8. Reproduzierbare Prüfung

```bash
python3 experiments/theory-contracts/seam-constructive-search-20261015/continuation/T1/check_t1_lift.py
```

Der Checker schreibt `result.json` und prüft:

- 65.536 Cocycle-Paare für den exakten Integer-C-Lift,
- alle 256 Ladungsparitäten für \(\widehat C^4=P_{\rm CAR}\) und die vierte Winding-Potenz,
- beide \(\chi\)-Orientierungen und alle acht sektor-diagonalen r-Lifts,
- explizit die skalare Unterklasse r=2,6 und den \(\kappa=-1\)-Gegenfall r=1,
- die vollständigen affinen Koeffizienten der Zeitidentität auf Ganz- und Halbgitter für beide Viertelwerte,
- die Bündelkovarianz von B1 und dessen vierte Potenz B4.

Die Paritätsreduktion ist vollständig, weil die Cocycle- und Liftphasen quadratische Ganzzahlpolynome modulo zwei sind. Die rationale Zeitidentität ist oben zusätzlich universell hergeleitet.
