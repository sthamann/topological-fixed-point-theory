# Fortsetzung: Stromgruppenlift des echten Achtzyklus

## Ergebnis

Der Stromgruppenlift ist **in der bereits stipulierten positiv definiten E8-Fockdarstellung mit gemeinsamem Cocycle und dem dort definierten globalen Ladungstransport X** stark genug, die diskrete Familie aus [`../T1/PROOF.md`](../T1/PROOF.md) auf diesem Raum zu verengen:

\[
 r\in\{2,6\},\qquad \nu=\pm i,\qquad \kappa=+1.
\]

Die skalare Kovarianz des globalen X wird unten direkt aus seiner Kommutation mit den festen SU(8)-Wurzelströmen hergeleitet. Dafür wird X ausdrücklich nicht mit dem lokalen Feld Y(s,z) gleichgesetzt. Der ergänzende Zustandsbefund ist, dass der Halbzustand

\[
 s=\tfrac12(1,1,1,1,1,1,1,1)
\]

ein (SU(8))-Singulett ist. Die beiden Vorzeichen (\pm i) bleiben zwei Hebungen desselben ganzzahligen Endpunkts. Weder die bisherige Viertelmarkierung noch (K_\alpha) unterscheidet sie. Eine vollständige gemeinsame Darstellung mit den odd Z8-Vektorfeldern wird daraus nicht behauptet.

Hier bezeichnet \(\kappa\) das **Vorzeichen** in \(\widehat P X\widehat P^{-1}=\kappa X\). Der Suchcode speichert dagegen einen Exponenten \(\kappa_{\rm exp}\in\mathbb Z/4\), mit \(\kappa=i^{\kappa_{\rm exp}}\). Das Ergebnis \(\kappa=+1\) bedeutet daher \(\kappa_{\rm exp}=0\); das frühere Gegenbeispiel \(\kappa=-1\) hat Exponent 2.

Das genaue Verdict lautet

```text
CONDITIONAL_CURRENT_GROUP_NARROWING_R_2_6_KAPPA_PLUS
```

Es gilt nicht schon für die rohe Integer-CAR-Endpunktwirkung. Dort bleiben die gradingabhängigen Gegenbeispiele mit (\kappa=-1) aus der vorigen Fortsetzung gültig. Für die vollständigen vier Suchpräsentationen auf Z8 vereinigt mit Z8+s ist eine zusätzliche kompatible Sektordarstellung erforderlich.

## 1. Was die Originalquellen tatsächlich voraussetzen

Die ursprüngliche Nahtrechnung konstruiert den positiven Einzelwindungsschritt auf den ganzzahligen CAR-Moden und seinen eingeschränkten Fock-Implementer nur bis auf eine zentrale Phase; sie führt keinen Halbsektorgenerator (X) ein ([`seam-joint-ward-state-20261014/PROOF.txt`, Zeilen 752–766](../../../seam-joint-ward-state-20261014/PROOF.txt)). Auf Ladungen ist

\[
 T_1(p)=Cp+e_0,
 \qquad C(p_0,\ldots,p_7)=(p_7,p_0,\ldots,p_6),
\]

mit (T_1^4(p)=Pp+a) und (T_1^8(p)=p+q) (ebenda, Zeilen 791–811).

Die tatsächlich passende positive Quelle ist §20 desselben Vertrags: Dort ist die unitäre integrable E8-Level-1-Darstellung auf \(L_s=D_8\cup(D_8+s)\) ausdrücklich **stipuliert**, einschließlich echter unitärer SU(2)-Stromaktionen auf der ganzen Darstellung ([Zeilen 934–959](../../../seam-joint-ward-state-20261014/PROOF.txt)). §21 verwendet den geordneten ursprünglichen CAR-Cocycle, seine explizite Fortsetzung auf \(L_s\) und denselben globalen Transport \(X=E_s\), mit \(X^2=-E_q\) und \(XE_dX^{-1}=(-1)^{s\cdot d}E_d\) für \(d\in D_8\) ([Zeilen 1003–1087](../../../seam-joint-ward-state-20261014/PROOF.txt)). Die Quelle begrenzt das Resultat ausdrücklich auf die gerade D8/E8-Darstellung; odd Z8-Vektorfelder sind dort nicht zugleich lokal eingeschlossen.

Die ältere [`source-sector-reconstruction-20260920/PROOF.txt`](../../../source-sector-reconstruction-20260920/PROOF.txt) behandelt dagegen \(\Gamma=\mathbb Z^{10}\) mit Signatur (9,1) und einen anderen Halbvektor. Ihre volle Vertex-/Ramondprämisse ist **kein Transportbeweis** in den positiv definiten Achtkanalraum \(\mathbb Z^8+s\). Sie begründet hier weder das SU(8)-Singulett noch die fehlenden odd Intertwiner. Eine Erweiterung auf alle odd Z8-Sektoren muss deren gemeinsame Cocycle-, Adjunktions- und Oszillatordarstellung zusätzlich liefern.

Die verwendete Vertexfaktorisierung und die Exponentiation von Gewichts-eins-Nullmoden sind für positiv definite gerade Gitter in Dong–Nagatomo, [*Automorphism Groups and Twisted Modules for Lattice Vertex Operator Algebras*, §2.2, §§2.3–2.4](https://arxiv.org/pdf/math/9808088), beschrieben. Diese Primärquelle liefert die Darstellungsformeln; sie identifiziert keine andere physische Quelle mit dem hier stipulierten E8-Lift.

## 2. Der konkrete Generatorweg zum Achtzyklus

Für die sieben einfachen (A_7)-Wurzeln

\[
 \alpha_i=e_i-e_{i+1},\qquad i=0,\ldots,6,
\]

setze in der Vektordarstellung den tatsächlichen Einparameterweg

\[
 n_i(t)=\exp\!\left[t\frac\pi2
 (E_{i+1,i}-E_{i,i+1})\right].
\]

Bei (t=1) ist der nichttriviale Zweierblock

\[
 \begin{pmatrix}0&-1\\1&0\end{pmatrix},
\]

also (\det n_i(t)=1). Mit

\[
 R(t)=n_0(t)n_1(t)\cdots n_6(t),
 \qquad R=R(1)
\]

gilt auf der Standardbasis exakt

\[
 Re_j=e_{j+1}\quad(0\le j<7),
 \qquad Re_7=-e_0,
 \qquad R^8=-I.
\]

Der fehlende Endpunktsign ist der U(1)-Torusfaktor. Die beiden konkreten
Wege zum selben Endpunkt sind

\[
 Z_\pm(t)=\exp(\pm i\pi t H_{e_0}),
 \qquad
 Z_\pm(1)=Z_{e_0}=\operatorname{diag}(-1,1,\ldots,1).
\]

Damit sind

\[
 h_\pm(t)=Z_\pm(t)R(t)
\]

zwei ausdrücklich angegebene Stromgruppenwege mit (h_\pm(0)=I) und demselben
Endpunkt, dem tatsächlichen phasenfreien Kanalzyklus,

\[
 \boxed{C=Z_{e_0}R},
 \qquad C^4=P,
 \qquad C^8=I.
\]

Dies ist eine konkrete Stromgruppenzerlegung desselben (C), nicht die Wahl
eines unbestimmten Matrixlogarithmus und nicht bloß eine Permutation von
Ladungslabels. Der Test betrifft den linearen Zyklusfaktor im bereits festen
affinen Winding (T_1(p)=Cp+e_0); der Ladungscarry (e_0) wird nicht ersetzt.

## 3. Warum der Ramond-Halbzustand skalar transformiert

Für jede einfache Wurzel gilt

\[
 \alpha_i\cdot s=0.
\]

Im angenommenen Lattice-Vertexmodul folgt daher

\[
 H_{\alpha_i,0}|s\rangle=0.
\]

Auch die beiden Wurzelnullmoden annihilieren den Zustand. Denn im OPE von (J_{\pm\alpha_i}(z)) mit (V_s(0)) ist der führende Potenzfaktor

\[
 z^{\pm\alpha_i\cdot s}=z^0.
\]

Für ein Gewicht-eins-Stromfeld ist der Nullmodenkoeffizient der (z^{-1})-Term; ein solcher Term fehlt hier. Deshalb annihilieren die drei (\mathfrak{su}(2))-Generatoren jeder einfachen Wurzel (|s\rangle). Der gesamte durch die (n_i)-Wege erzeugte (SU(8))-Faktor (R) fixiert (|s\rangle) exakt.

Dieser Zustandsbefund allein ist noch keine Gleichung für den globalen unitären X. Die erforderliche Brücke folgt direkt aus den festen Wurzelströmen. Für \(\alpha=e_i-e_j\) gilt auf demselben Raum

\[
 J_\alpha(z)=\mathcal E^-_\alpha(z)\mathcal E^+_\alpha(z)
 E_\alpha z^{\alpha\cdot H_0},
\]

wobei Oszillatorfaktoren, Gruppenalgebra-/Kleinoperator \(E_\alpha\) und Cocycle genau diejenigen der festgelegten E8-Fockdarstellung sind. Der globale X ist hier der **Gruppenalgebra-Ladungstransport** \(E_s\), tensoriert mit der Identität auf den Heisenberg-Oszillatoren. Das vollständige Feld \(Y(s,z)\) enthält zusätzliche Oszillator- und Nullmodenfaktoren und ist nicht X.

Für jede dieser Wurzeln ist \(Q(\alpha)=\alpha\cdot s=0\). Daher liefern die gegebene Cocyclerelation und die Ladungswirkung

\[
 [X,E_\alpha]=0,\qquad
 [X,z^{\alpha\cdot H_0}]=0,\qquad
 [X,\mathcal E^\pm_\alpha(z)]=0.
\]

Somit kommutiert X mit jedem Koeffizienten \(J_{\alpha,m}\), insbesondere mit den kompakten Nullmodengeneratoren der sieben SU(2)-Wege. Dies ist eine Rechnung mit demselben Wurzelcocycle, keine Identifikation gleich großer Zustandsräume. Auf der dichten gemeinsamen Domäne der endlichen Ladungs-/Oszillatorzustände sind alle Faktoren definiert, und die Domäne wird durch X und die Strommoden erhalten. Die stipulierte unitäre Integrabilität liefert die echten kompakten Gruppenaktionen. Ihre Nullmoden erhalten die endlichen konformen Gewichtsräume; die dortigen Vektoren sind analytisch. Da X ein beschränkter unitärer Ladungstransport ist und mit den Generatoren auf diesem Kern kommutiert, gilt die Kommutation auch mit ihren unitären Exponentialen:

\[
 [X,R(t)]=0.
\]

Der feste Torusoperator erfüllt unabhängig davon

\[
 Z_\pm(t)X Z_\pm(t)^{-1}
 =e^{\pm i\pi t e_0\cdot s}X
 =e^{\pm i\pi t/2}X.
\]

Für die so gebildeten Stromgruppenhebungen \(\widehat C_\pm=Z_\pm(1)\widehat R\) folgt daher als **globale Operatorgleichung auf dem gegebenen E8-Raum**

\[
 \widehat C_\pm X\widehat C_\pm^{-1}=\pm iX.
\]

Die gradingabhängigen Operatoren

\[
 \widehat C_r|n,c\rangle
 =\gamma(n)\zeta_8^{rc}|Cn,c\rangle
\]

mit ungeradem (r) aus der T1-Fortsetzung sind gültige Darstellungen der diskreten gekreuzten Algebra. Ihre Konjugationsphase auf X hängt jedoch vom Eingangsgrad ab. Sie stimmen daher nicht mit dieser aus denselben Wurzelströmen und demselben Torus gebildeten zusammenhängenden Stromgruppenhebung überein. Das ist der zusätzliche Auswahlgehalt der stipulierten E8-Stromdarstellung. Ohne diese Identifikation der Stromgruppe bleiben sie gültige Gegenmodelle.

Für beide Suchcharaktere \(\chi=\pm i\) ist die Generatorrelation auf den SU(8)-Wurzeln dieselbe, weil \(Q(\alpha)=0\). Das beweist die gemeinsame SU(8)-Kommutation, aber nicht die Vereinbarkeit jedes odd \(E_{e_j}\)-Intertwiners mit einem einzelnen fest gewählten Feldcocycle \(Y(s,z)\). Dafür ist die volle kompatible Sektordarstellung erforderlich. Auch die Permutationsinvarianz formaler Oszillatorfaktoren würde ohne den nachgewiesenen gemeinsamen Ladungs-/Kleinfaktor diese Lücke nicht schließen.

## 4. Die zwei echten Spinhebungen

Es gibt zwei gerichtete Toruswege mit demselben Endpunkt auf der ganzzahligen Quelle:

\[
 Z_+(t)=\exp(+i\pi tH_{e_0}),
 \qquad
 Z_-(t)=\exp(-i\pi tH_{e_0}).
\]

Für jedes (n\in\mathbb Z^8) stimmen ihre Endpunkte überein,

\[
 e^{+i\pi n_0}=e^{-i\pi n_0}=(-1)^{n_0},
\]

und beide liefern zusammen mit (R) dieselbe Matrix (C). Auf dem Halbzustand unterscheiden sie sich jedoch:

\[
 Z_+(1)|s\rangle=+i|s\rangle,
 \qquad
 Z_-(1)|s\rangle=-i|s\rangle.
\]

Äquivalent haben die Determinantenwege die Windung (+1) beziehungsweise (-1), und ihre Halbdeterminante endet bei (+i) beziehungsweise (-i). Dies liefert

\[
 r=2\quad\text{oder}\quad r=6\pmod 8.
\]

In beiden Fällen

\[
 \nu^2=-1,\qquad \nu^4=+1.
\]

Deshalb ist die vierte Potenz auf dem Defektgrad trivial und der zugehörige (P)-Lift hat

\[
 \boxed{\kappa=+1}.
\]

Die beiden Wege unterscheiden sich um einen vollen U(1)-Determinantenloop. Der bisherige Integer-CAR-Endpunkt sieht diesen Loop nicht. Einen der beiden Wege als (den) Logarithmus von (Z_{e_0}) zu erklären wäre genau die zusätzliche, bisher unbegründete Wahl.

Der im Original positiv orientierte Operator (e^{i\theta}) ändert daran nichts
automatisch. Für nichtganzzahliges (t) bildet (e^{it\theta}) den fest
getwisteten Einteilchenraum in einen Raum mit verschobener Randphase ab; er
ist daher ohne eine zusätzliche Trivialisierung kein Weg in der festen
globalen (U(8))-Faser. Er fixiert den affinen Ladungscarry und die Richtung
der Wardverschiebung, aber noch keine Quadratwurzel des Determinanten auf dem
hinzugenommenen Ramond-Sektor.

## 5. Was Viertelphase und Zeit unterscheiden

Beide Hebungen besitzen dieselbe Einteilchen-Endpunktwirkung. Daher bleibt

\[
 D_\alpha C=\zeta_8CD_\alpha
\]

für (\alpha=1/4,3/4) unverändert. Ebenso gilt auf Ganz- und Halbgitter

\[
 K_\alpha(Cp+e_0)-K_\alpha(p)
 =\frac{Q(p)+\alpha}{8}.
\]

Weder (r) noch das Vorzeichen des Determinantenwegs kommt in einer dieser Identitäten vor. Die vorhandene Viertelmarkierung und die feste Zeit wählen ( +i ) gegen ( -i ) daher nicht aus. Die beiden Braidingpräsentationen (\chi=\pm i) werden durch diese SU(8)-Wurzelrechnung nicht unterschieden. Im gegebenen globalen Operatorrahmen sind sie durch das vorhandene Paritätsdressing \(X_+=\Pi X_-\) äquivalent, solange X nicht selbst als zusätzliche konkrete Operation markiert wird; siehe [`../branch_equivalence/PROOF.md`](../branch_equivalence/PROOF.md). Eine lokale Feldidentifikation folgt daraus nicht.

Auch der räumliche Bündelweg (B_1(x)=e^{ix/8}C) schließt diese Lücke nicht.
Seine Determinante windet bei (x:0\to2\pi) einmal positiv, wodurch eine
stetige Halbdeterminante ihr Vorzeichen wechselt. Das bestimmt die relative
Monodromie entlang des Kreises, aber nicht ihren Anfangswert bei (x=0):
(+i)e^{ix/2} und (-i)e^{ix/2} erfüllen dieselbe Vektor-Bündelkovarianz. Die
Werte (\alpha=1/4,3/4) treten in (D_\alpha) ein, wählen aber keinen dieser
beiden Anfangswerte.

Die erste verbleibende Zusatzwahl ist somit genau:

> Welche orientierte Hebung des U(1)-Determinantenwegs wirkt auf dem Ramond-Halbzustand?

Das ist eine Ramond-/Spinliftkonvention, kein neuer Hamiltonoperator und kein neuer Zustand. Eine Quellenwahl dafür müsste eine Wirkung des **Wegs**, nicht nur seines ganzzahligen Endpunkts, festlegen.

## 6. Geltungsgrenze

Der neue Schluss ist eine Verengung innerhalb des ausdrücklich gegebenen positiv definiten E8-Stromgruppenlifts:

```text
diskreter Integer-CAR-Endpunkt: r=0,...,7 und kappa=±1
stipulierter E8-Lift mit gemeinsamem X/Stromcocycle: r=2 oder 6 und kappa=+1
noch offen: Orientierung +i gegen -i und kompatible volle odd-Z8-Sektordarstellung
```

Er leitet weder diese Vertexkonstruktion aus der endlichen Quelle noch ihren Zustand, Energie oder physische Auswahl aus P1/P2 her. Er überträgt insbesondere nicht die anders signierte Γ=(9,1)-Quelle auf diesen Raum und behauptet keine lokale Vollalgebra mit odd CAR-Feldern. Der frühere even-sector Vertrag erlaubt ohne die zusätzliche Wahl dieser konkreten Stromgruppenhebung weiterhin beide (P)-Lifts ([`seam-joint-ward-state-20261014/PROOF.txt`, Zeilen 1034–1054](../../../seam-joint-ward-state-20261014/PROOF.txt)).

## 7. Reproduzierbare exakte Prüfung

```bash
python3 experiments/theory-contracts/seam-constructive-search-20261015/continuation/current_group/check_current_group.py
```

Der Rechner prüft ausschließlich exakte Identitäten:

- die sieben (SU(8))-Titsfaktoren und (C=Z_{e_0}R),
- (R^8=-I), (C^4=P) und (C^8=I),
- die Singulettbedingung (\alpha_i\cdot s=0) für alle sieben einfachen Wurzeln,
- die exakten Charakter-, Nullmoden- und Toruskoeffizienten der direkten globalen X/Wurzelstrom-Kommutation,
- die zwei gerichteten Halbdeterminantenphasen (r=2,6) und (\kappa=+1),
- die Gleichheit beider Torusendpunkte auf allen ganzzahligen Testladungen,
- die unveränderte rationale Einzelwindungs-Wardidentität auf Ganz- und Halbgitter sowie die Viertelbündelidentität.

Die vollständige Feldfaktorisierung, der gemeinsame Kern und die Integration sind die oben offen ausgewiesenen analytischen Schritte innerhalb der stipulierten E8-Darstellung. Der endliche Rechner prüft ihre diskreten Koeffizienten und ersetzt diese Darstellungsprämisse nicht.
