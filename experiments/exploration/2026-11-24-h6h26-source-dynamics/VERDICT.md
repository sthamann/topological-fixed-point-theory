# TFPT H6h27: Zustand, metrischer Prozess und Vierpunktantwort

Stand: 24. November 2026. Fortsetzung des Forschungsstands H6h26, anschließend mit H6h27 auf Commit b5cc3d4070a1d869cafc3c414cf0cb8b2d4f1401 zusammengeführt.

## Ergebnis und Geltungsbereich

Die verbundene Vierpunktantwort ist für die quartische Gegenfamilie jetzt explizit berechnet. Sie kann den freien Parameter messen, wählt ihn aber nicht ohne eine unabhängig bestimmte mikroskopische Vierpunktfunktion. Außerdem lässt die Kantenregel zunächst einen Parameter je Kante frei. Die Vierteldrehung setzt die vier Arme bei gleicher Tiefe gleich; die Gleichheit verschiedener Tiefen ist eine weitere Regel.

Die Zustandsfrage lässt sich teilweise schließen: Auf einem festgelegten regulären Stromoszillator erzwingt die volle KMS Bedingung den thermischen Gaußzustand. Der dazu konstruierte quadratische Zweimodenoperator ist jedoch nicht der bisherige ungewichtete Hardy Operator. Er besitzt denselben geometrischen Grundzustand, aber ein anderes, diskretes Spektrum. Beide Operatoren lassen sich in einen positiven linearen metrischen Stern einbetten.

Dies sind exakte Folgerungen innerhalb ausdrücklich benannter mathematischer und physikalischer Voraussetzungen. Eine gemeinsame wechselwirkende Quantentheorie von Gravitation, Naht, chiraler Materie, Eichfeldern und Higgs ist damit nicht konstruiert. Es erfolgt keine Beförderung im Ledger und keine neue Behauptung, die gesamte Theorie sei eindeutig.

| Frage | Ergebnis dieser Fortsetzung | Verbleibende Voraussetzung |
| --- | --- | --- |
| Gaußzustand des Strommodus | Eindeutig aus voller KMS Bedingung | Reguläre CCR Algebra, physische Zeit, thermische Normierung |
| Wert von q | q = exp(−2π) in der kanonischen Purifikation | βω = 4π und Identifikation der beiden Blattalgebren |
| Positiver Präparationsoperator | Quadratische Klasse wählt A†(N+1)A bis auf Skala | Zwei CCR Moden, relative Ladungserhaltung, Blatttausch |
| Metrischer Stern | Eine universell gekoppelte Nullmode; gewichtete Variante vollständig berechnet | Metrische Wirkung und minimale physische Quellenkopplung |
| Quartische Wechselwirkung | Vierpunktantwort explizit, Parameter pro Tiefenorbit | Mikroskopischer Quellenabgleich oder neue Homogenitätsregel |
| G und massive Skala | Freie Normierungen bleiben nachweisbar | Physische Skalenregel |
| TM1 bei 15 Grad | Direkte endliche Gruppenroute im präzisen Restklassenansatz ausgeschlossen | Dynamische Ausrichtung in einer kleineren Restsymmetrie |
| Naht, chirale 4D Materie, gemeinsame Prozessregel | Nicht geschlossen | Eigene nichtperturbative Konstruktion |

Die vollständigen Beweise stehen in den begleitenden Dateien unter proofs. Die nachfolgenden Abschnitte erklären die entscheidenden Übergänge und ihre Grenzen.

## 1. Die volle KMS Bedingung liefert mehr als ein Besetzungsmoment

Für einen regulären Oszillator mit [a,a†] = 1 und festgelegter Dynamik α_t(a) = exp(−it)a macht Stationarität die normale Zustandsdichte diagonal im Zahloperator N. Mit analytischen Matrixeinheiten liefert die KMS Bedingung

\[
p_{n+1}=e^{-\beta}p_n.
\]

Damit folgt eindeutig

\[
\rho_\beta=(1-r)r^N,\qquad r=e^{-\beta}.
\]

Die kanonische Purifikation bezüglich einer festgelegten konjugierten rechten Algebra ist

\[
|\eta_q\rangle=\sqrt{1-q^2}\sum_{n\ge0}q^n|n,n\rangle,
\qquad q=e^{-\beta/2}.
\]

Für die vorgegebene thermische Normierung β = 4π folgt q = exp(−2π). Die Gaußform dieses einzelnen Strommodus ist damit keine zusätzliche unabhängige Prämisse. Die Existenz und physische Bedeutung der Stromalgebra sowie ihrer Zeitentwicklung sind weiterhin Voraussetzungen. Insbesondere ist damit weder der vollständige E8 Zustand noch die physische Raumzeit hergeleitet.

Ein exaktes Gegenbeispiel verhindert den schwächeren Schluss nur aus der mittleren Besetzung:

\[
\rho_{\rm alt}=\rho_\beta+
\epsilon(P_0-2P_1+P_2),\qquad
\epsilon=(1-r)r/4.
\]

Diese Dichte ist positiv, treu, normiert und hat denselben Wert von ⟨N⟩. Ihre volle KMS Relation und ihre Wick Identitäten scheitern. Eine Übereinstimmung des Zweipunktmoments allein ist daher keine Zustandsrekonstruktion.

Der Beweis benutzt den normalen regulären Oszillatorfaktor. Die Literatur behandelt die entsprechenden quasifreien Zustände und ihre KMS Voraussetzungen; die hier benutzte Matrixeinheitenrechnung macht die benötigten Annahmen ausdrücklich sichtbar. Siehe [Dereziński und Gérard, Kapitel 17](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/AC95B25B70BCC64D668945FE64D09239/9781009290821c17_423-474.pdf/quasifree_states.pdf) sowie [Cottrell et al. zur Thermofeldkonstruktion](https://arxiv.org/abs/1811.11528).

## 2. Der Stromoperator ist nicht der Hardy Shift

Die Aussage von Z48 über Modenetiketten bleibt richtig: Die Kommutatorwirkung von J₋₁ verschiebt einen Fermionindex um eins und damit den entsprechenden Viertelindex um vier. Auf dem Stromfockraum gilt jedoch

\[
J_{-1}|n\rangle=\sqrt{n+1}|n+1\rangle.
\]

Schon die Normenfolge 1, 1, 2 für Ω, J₋₁Ω, J₋₁²Ω widerspricht einer Isometrie. Der Shift ist erst die Polarphase

\[
V=J_{-1}(N+1)^{-1/2}.
\]

Die Identifikation von V mit S⁴ in einem Restklassensektor ist eine gesonderte, zulässige Abbildung. Sie identifiziert nicht die gesamte Stromalgebra mit einer einzigen Hardy Leiter. Die vollständigen Gradmultiplizitäten sind ebenfalls verschieden. Die zugrunde liegende Stromalgebra ist in [Alexandrov und Zabrodin](https://arxiv.org/abs/1212.6049) beschrieben.

### Ein positiver quadratischer Operator für den Thermofeldzustand

Setze

\[
d_L=\frac{a_L-q a_R^\dagger}{\sqrt{1-r}},
\qquad
d_R=\frac{a_R-q a_L^\dagger}{\sqrt{1-r}}.
\]

Beide vernichten den kanonischen Thermofeldzustand. Unter den Voraussetzungen

1. genau diese zwei CCR Moden,
2. ein positiver quadratischer Operator mit diesem eindeutigen Grundzustand,
3. Erhaltung von N_L−N_R,
4. Austausch der beiden Moden

ist der Operator bis auf Ω > 0 festgelegt:

\[
H_q=\Omega(d_L^\dagger d_L+d_R^\dagger d_R).
\]

Im Paarsektor |n,n⟩ wird daraus

\[
K_\lambda=\lambda A_q^\dagger(N+1)A_q,
\quad
(A_qf)_n=qf_n-f_{n+1},
\quad
\lambda=\frac{2\Omega}{1-r}.
\]

Die modularen Aussagen sind ebenfalls präzise: Die modulare Konjugation dieses bosonischen Standardraums hat J² = +1. Sie darf nicht ohne Weiteres mit einer fermionischen Blattoperation mit Quadrat −1 identifiziert werden.

Die zusätzliche Gewichtung N+1 ist entscheidend. Der frühere Operator λA†A hat denselben Grundzustand, ist aber in diesen mikroskopischen CCR Variablen nicht der oben bestimmte quadratische Operator. Andere positive Gewichte W(N) erhalten ebenfalls den Grundzustand. Die Zustandswahl allein beseitigt diese dynamische Freiheit nicht.

## 3. Das vollständige Spektrum und der gewichtete metrische Stern

Für den quadratischen Operator gilt mit Δ = λ(1−r)

\[
\operatorname{spec}K_\lambda=\{\Delta n:n=0,1,2,\ldots\}.
\]

Sein Randmaß ist positiv und diskret:

\[
C_{\rm raw}(x)
=(1-r)\sum_{n\ge0}\frac{r^n}{x+\Delta n}.
\]

Nach Normierung des masselosen Residuums auf eins:

\[
C(x)=\frac1x+\sum_{n\ge1}\frac{r^n}{x+\Delta n}.
\]

Damit ist die skalare Källén–Lehmann Positivität explizit. Für den inversen Kern H = 1/C folgen

\[
H(x)=x-a x^2+b x^3+\cdots,
\]

\[
a=\frac{-\log(1-r)}{\Delta},\qquad
b=\frac{\log^2(1-r)+\operatorname{Li}_2(r)}{\Delta^2}.
\]

Ein willkürlich gesetzter Wert von a kann weiterhin durch die freie Skala Δ angepasst werden. Daraus folgt kein unabhängiger Wert von q.

### Einbettung in alle vier Arme

Die Zentralmetrik besitzt die Norm M_c. In Tiefe n des Arms j sei

\[
M_{j,n}=\frac{M_c}{2}q^{n+1}.
\]

Die kanonische Differenzmatrix B ist dieselbe wie in H6h27. Mit Kantensteifigkeit

\[
k_{j,n}=\lambda(n+1)M_{j,n}^2
\]

lautet der lineare Massenoperator

\[
L=\lambda B^\dagger W B,\qquad W_{j,n}=n+1.
\]

Die diskrete Fouriertransformation der vier Arme gibt den vollständigen Vierteltransfer mit diesem Gewicht. Es wurde nicht nur der zentrale Sektor betrachtet: Neben der symmetrischen Leiter liegen drei identische Dirichletkompressionen vor. Deren Eigenwerte liegen zwischen je zwei benachbarten symmetrischen Eigenwerten. Deshalb gibt es im unendlichen gewichteten Modell nur endlich viele Moden unter jeder festen endlichen Massenschwelle. Beim alten beschränkten Bandoperator war diese Eigenschaft nicht gegeben.

Für endliche physische Arme ist der letzte Diagonaleintrag λD, nicht die bloße obere Ecke der unendlichen Matrix. Die positiven endlichen Graphen konvergieren im starken Resolventensinn zum selbstadjungierten unendlichen Operator. Daraus folgt kein Beweis der Existenz einer vollständigen unendlichen wechselwirkenden Quantentheorie.

### Universelle Nullmode

Für jede positive Gewichtung bleibt

\[
\eta_v=\frac{M_v}{M_{\rm eff}},
\qquad
M_{\rm eff}^2=\sum_vM_v^2.
\]

Mit der ausdrücklich angenommenen minimalen metrischen Quellenkopplung h_v:T_v/M_v koppelt die Nullmode an alle T_v gleich:

\[
\frac{\eta_v}{M_v}=\frac1{M_{\rm eff}}.
\]

Im unendlichen Stern gilt

\[
M_{\rm eff}^2=\frac{M_c^2}{1-r},
\qquad
G=\frac{1-r}{8\pi M_c^2}.
\]

Dies beweist Universalität innerhalb des metrischen Modells. Es beweist nicht, dass ein bestimmter mikroskopischer Strom aus TFPT bereits der physische Energieimpulstensor ist. Die gewöhnlichen HR Paarpotentiale besitzen die bekannte klassische Zwangsstruktur; siehe [Hassan und Rosen](https://arxiv.org/abs/1109.3515) und für die Graphfrage [Hinterbichler und Rosen](https://arxiv.org/abs/1203.5783).

## 4. Was die Quellenargumente tatsächlich leisten

Bei festem q ist der Vierteltransfer unter μ₄ invariant, nicht unter der gesamten kontinuierlichen Kreisgruppe. Eine Viertelmittelung lässt auch höhere Moden mit Index durch vier teilbar stehen. Selbst eine zusätzliche kontinuierliche Mittelung lässt neutrale zusammengesetzte Operatoren wie J₋ₙJₙ bestehen.

Die Ladung J₀ ist deshalb nicht allein aufgrund ihrer Invarianz als universeller Energieimpulstensor ausgewiesen. Der metrische Ansatz liefert eine gravitative Ward Struktur, wenn die Materiewirkung bereits diffeomorphismenkovariant und minimal an die jeweiligen Metriken gekoppelt ist. Die fehlende Aufgabe ist der Abgleich mit der ursprünglichen mikroskopischen Quelle.

Die MICROSCOPE Abschätzung aus Z48 ist unter ihrem eigenen masselosen Quellenansatz verwendbar. Bei gleichem normierten Monopol folgt

\[
\eta_{AB}
=\frac{q(h_A-h_B)}
{1+q(h_A+h_B)/2}.
\]

Ein gemeinsamer höherer Anteil verletzt das Äquivalenzprinzip nicht. Eine grobe Schranke |η| ≲ 10⁻¹⁴ entspräche bei kleinem Nennerkorrekturterm |Δh| ≲ 5,35 × 10⁻¹². Das ist weder ein exakt ausgewertetes Konfidenzintervall noch ein Beweis, dass sämtliche höheren Quellenmoden verschwinden. Das tatsächliche Experiment berichtet getrennte statistische und systematische Unsicherheiten; siehe [MICROSCOPE, Endergebnis](https://arxiv.org/abs/2209.15487).

Auch die Zweimetrikenform ist nur ein Ersatz mit abgestimmten Momenten. Ein einzelner massiver Pol kann nicht das ganze Marchenko–Pastur Band reproduzieren. Bei gleichem massivem Gesamtgewicht und gleicher erster Inverskorrektur gilt exakt b_Band−b_Pol = a². Die volle Spektralfunktion unterscheidet beide bereits auf Zweipunktebene.

## 5. Die quartische Freiheit ist jetzt als Quellenantwort berechnet

### Exakte Klassifikation des Paarpotentials

In der fünfparametrigen HR Klasse sind die zwei Flachheitsbedingungen und die normierte Fierz–Pauli Masse drei unabhängige lineare Bedingungen. Die vollständige verbleibende Familie lautet

\[
\beta=\beta_V+\alpha(-2,1,0,-1,2)
+\gamma(1,-1,1,-1,1),
\]

mit β_V = (−2, 1/2, 0, 1/2, −2). Eine ungerichtete Paarpotentialdichte verlangt β_n = β_{4−n} und setzt α = 0.

In der bisherigen Notation ist γ = t/2:

\[
V_t(X)=V(X)+\frac t2\det(I-X).
\]

Die vollständige Kubik ist damit festgelegt. Der Determinantenterm beginnt in vierter Ordnung. Diese Aussagen beruhen nun auf exakten Koeffizientenidentitäten für die ganze Klasse und nicht auf zwei zufälligen Taylorproben.

Noch stärker: Die vollständige homogene metrische Quartik bestimmt innerhalb dieser festen fünfparametrigen HR Klasse alle fünf β Werte. Deterministische Monomkoeffizienten haben Rang fünf. Das bedeutet, eine unabhängig bekannte volle Quartik könnte das Paarpotential festlegen. Es bedeutet nicht, dass diese Quartik bereits bekannt ist. Feldumdefinitionen, andere Wirkungen und zusätzliche höhere Ableitungen sind nicht Teil dieses Satzes.

### Ein Parameter pro Kante

Bei D Armlagen gibt es zunächst 4D unabhängige Parameter t_e. Ein Schnitt einer einzelnen Baumkante liefert einen Zeugen ihrer Unabhängigkeit: Auf jeder der beiden entstehenden Komponenten werden die Metriken gleich gewählt; dann verschwindet die Determinantendeformation auf allen anderen Kanten.

μ₄ reduziert die Freiheit auf genau D Parameter t_n. Nur eine zusätzlich verlangte Homogenität über alle Tiefen macht daraus eine einzige Zahl t.

### Verbundene Vierpunktfunktion

Definiere den kanonischen Propagator G(p²) = (p²I+L)⁻¹ und die physische Kantenantwort auf eine Quelle am Knoten b:

\[
R_{e,b}(p^2)
=\frac{G_{c(e),b}(p^2)}{M_{c(e)}M_b}
-\frac{G_{p(e),b}(p^2)}{M_{p(e)}M_b}.
\]

Für konservierte spurfreie Quellen und

\[
W_E[J]=\log\int e^{-S_E+\sum_v h_v:J_v/M_v}
\]

lautet die Differenz der verbundenen Vierpunktantworten auf Baumebene exakt

\[
\Delta W_E^{(4)}
=-\frac12\sum_e k_e t_e
\prod_{i=1}^4R_{e,b_i}(p_i^2)\,
\mathcal D(T_1,T_2,T_3,T_4).
\]

Dabei ist 𝒟 der Koeffizient von z₁z₂z₃z₄ in det(Σ_i z_iT_i). Es gibt keinen zusätzlich versteckten Faktor 4!. Die Impulserhaltung ist vorausgesetzt.

Der Grund für die überschaubare Formel: Quadratische und kubische Vertizes hängen nicht von t ab. Daher heben sich alle Austauschdiagramme beim Vergleich mit t = 0 auf. Die Formel ist die Differenz der vollständigen Baumantworten, nicht die gesamte Vierpunktfunktion.

Für reine masselose Eigenmoden verschwindet jeder R Unterschied. Bereits ein einziges reines masseloses externes Bein setzt diesen Kontakt auf null. Generische Randquellen enthalten dagegen massive Anteile und sehen ihn. Ein expliziter konservierter spurfreier Zeuge mit nichtkollinearen Impulsen hat 𝒟 = −16.

### Geschlossene Antwort des unendlichen flachen Sterns

Setze r = q², x_i = p_i²/λ und

\[
z_i=\frac{2}
{x_i+1+r+\sqrt{(x_i+1+r)^2-4r}},
\qquad P=\prod_{i=1}^4z_i.
\]

Für vier zentrale Quellen ist die Kantenantwort der Tiefe n

\[
R_n(p^2)=-\frac{z(p^2/\lambda)^{n+1}}{\lambda M_c^2}.
\]

Bei gemeinsamem t folgt

\[
\Delta W_E^{(4)}
=-\frac{t}{2\lambda^3M_c^6}
\frac{rP}{1-rP}\,\mathcal D.
\]

Bei verschiedenen, beschränkten Tiefenparametern lautet der entsprechende Faktor

\[
F(P)=-\frac1{2\lambda^3M_c^6}
\sum_{n\ge0}t_n(rP)^{n+1}.
\]

Das ist eine erzeugende Funktion der freien Wechselwirkungen. Wenn die mikroskopische Vierpunktfunktion unabhängig und mit gleicher Quellennormierung gegeben wäre, würden ihre analytischen Koeffizienten die t_n bestimmen. Der Grenzwert kleiner Impulse allein liefert nur eine gewichtete Summe.

### Unterscheidung der beiden Dynamiken im Infraroten

Der endliche Baum erfüllt die allgemeine Flussidentität

\[
R_{e,b}(0)
=\frac{\mathbf1_{b\in S_e}-M_{S_e}^2/M_{\rm eff}^2}{k_e}.
\]

Für vier zentrale Quellen und konstantes t ergibt sich daraus:

\[
F_{\rm flach}(0)
=-\frac{t}{2\lambda^3M_c^6}\frac r{1-r},
\]

\[
F_{\rm gewichtet}(0)
=-\frac{t}{2\lambda^3M_c^6}\operatorname{Li}_3(r).
\]

Die Vierpunktantwort unterscheidet also auch den Hardy Operator vom quadratischen CCR Operator. Die unamputierte Antwort bleibt im Allgemeinen endlich. Vier inverse masselose Randpropagatoren bringen dagegen den Faktor ∏p_i²; die reine masselose S Matrix erhält daraus keinen t Beitrag.

Im Zweimetrikenfall ist t auch in einer physikalischen massiven Helizität 0 Amplitude sichtbar:

\[
\Delta\mathcal M_{0000}(s,0)
=\frac{t\,s(s-4m^2)}
{9M_{\rm rel}^2m^2}.
\]

Diese Amplitudendifferenz allein unterliegt keinem isolierten Positivitätsgebot. Der masselose Vorwärtspol und die übrige Amplitude müssen bei einer Dispersionsanalyse behandelt werden. Die bekannten eikonalen Bedingungen lassen ebenfalls eine quartische Familie übrig; siehe [Bonifacio et al.](https://arxiv.org/abs/1712.10020).

### Weshalb ein Gaußzustand t nicht automatisch auf null setzt

Ein Gaußzustand eines Strommodus bestimmt nicht die gravitative Vierpunktfunktion. Ein vollständig gaußscher Prozess mit linearem physischem Quellenabbild hätte zwar verschwindende höhere verbundene Kumulanten. Dann wäre aber schon die wechselwirkende Einstein Struktur bei endlichem G zu erklären; das Problem würde nicht bloß durch t = 0 gelöst.

Für ein Matching werden die gesamte Tensorstruktur, die Impulsabhängigkeit, die übrigen Kontakte, Austauschdiagramme und lokale Quellenterme benötigt. Erst danach ist eine Gleichung für t oder die Folge t_n eine physische Bestimmung.

## 6. Endliche Tiefe, Kopplung und G

Für d = n+1 gilt M_d = M_cq^d/2. Mit ρ₁ = q/2, ρ_d = q für d ≥ 2 und w_d = 1 beziehungsweise d sind die lokalen Größen des isolierten Metrikpaares

\[
M_{{\rm rel},d}=\frac{M_d}{\sqrt{1+\rho_d^2}},
\qquad
m_{{\rm edge},d}^2=\lambda w_d(1+\rho_d^2).
\]

Die übliche lokale Abschätzung Λ₃ = (m_edge²M_rel)^{1/3} ist keine Berechnung aller kollektiven longitudinalen Moden. Die Abhängigkeit der tatsächlichen Skala von der Graphstruktur wird beispielsweise in [Scargill und Noller](https://arxiv.org/abs/1511.02877) untersucht.

Mit dem deklarierten Beispiel M_c = 2,435 × 10¹⁸ GeV und √λ = 1 TeV folgen:

| Geforderter lokaler Abstand | Flacher Stern | Gewichteter Stern |
| --- | ---: | ---: |
| Λ₃/m_edge > 1 | 5 Lagen | 5 Lagen |
| Λ₃/m_edge > 10 | 4 Lagen | 4 Lagen |
| Λ₃/m_edge > 100 | 3 Lagen | 3 Lagen |

Auf der fünften Lage sind die Verhältnisse nur etwa 3,024 und 2,312. Die fünf Lagen sind daher keine Aussage über fünf großzügig kontrollierte Energiehierarchien.

Der räumliche normierte Nullmodenschwanz nach D Lagen ist q^{2(D+1)}. Nach fünf Lagen ist das q¹² ≈ 1,80 × 10⁻³³. Dies unterdrückt nicht automatisch sämtliche Quanteneffekte, denn die physische Nullmodenkopplung η_v/M_v bleibt universell. Die frühere Beschriftung über n_max war um eine Lage verschoben.

### Ein konstruktiver endlicher Abschluss

Mit auf 0,…,D abgeschnittenen Oszillatoren liefert

\[
\mathcal H_D=\frac{\lambda}{2}
\left[
(a_D-qb_D^\dagger)^\dagger(a_D-qb_D^\dagger)
+(b_D-qa_D^\dagger)^\dagger(b_D-qa_D^\dagger)
\right]
\]

genau den physischen endlichen gewichteten Paarsektor samt Enddiagonale λD. Der Grundzustand ist die abgeschnittene geometrische Purifikation. Ihre Randänderung ist ausdrücklich

\[
[a_D,a_D^\dagger]
=I-(D+1)|D\rangle\langle D|.
\]

Das ist eine positive endliche Konstruktion für jede gewählte Tiefe. Sie wählt D nicht aus und erfüllt am Endpunkt nicht mehr die unendliche CCR Algebra des ursprünglichen Eindeutigkeitssatzes.

### Warum die vorliegenden Strukturdaten G nicht bestimmen

Für jedes c > 0 erhält die Transformation

\[
M_v\mapsto cM_v,\qquad k_e\mapsto c^2k_e
\]

den kanonischen Massenoperator, q, die Nullmodenrichtung und alle normierten Spektralformen. Die physischen Quellenkopplungen ändern sich um 1/c und G um 1/c². Ohne zusätzliche physische Normregel kann G deshalb nicht aus genau diesen Daten folgen.

Wenn die Einstein Terme bereits in der Wirkung vorausgesetzt sind, folgen ihre Feldgleichungen durch Variation. Jacobsons Annahme wäre dann eine zusätzliche mikroskopische Begründungsroute, keine noch benötigte algebraische Ableitung derselben bereits eingesetzten Wirkung.

## 7. Flavour: Eine ganze Suchrichtung ist jetzt ausgeschlossen

Für den genau benannten Zielansatz mit geladener Z₃,

\[
T=\operatorname{diag}(1,\omega,\omega^2),
\]

voller Majorana Klein Restsymmetrie und μτ Reflexion gilt für die Spalten u_i der Mischungsmatrix

\[
S_i=2u_i u_i^\dagger-I,
\qquad
\operatorname{tr}(TS_i)=3|U_{ei}|^2-1.
\]

TM1 bei 15 Grad ergibt die Spuren

\[
1,\qquad-\frac{2-\sqrt3}{4},
\qquad-\frac{2+\sqrt3}{4}.
\]

Die beiden letzten Zahlen erfüllen z²+z+1/16 = 0 und sind keine algebraischen ganzen Zahlen. Eine endliche unitäre Gruppe besitzt dagegen für jedes Gruppenelement eine Spur, die Summe von Einheitswurzeln und deshalb algebraisch ganz ist.

Damit ist diese exakte direkte Realisierung in jeder endlichen Gruppe ausgeschlossen, unabhängig von ihrer Größe. Weitere Scans größerer Gruppen würden dieses Hindernis nicht lösen. Eine semidirekte Route mit nur einer verbleibenden Z₂ und CP ist nicht betroffen, enthält aber einen freien Winkel; vgl. [Ding, King und Neder](https://arxiv.org/abs/1409.8005).

Ein stabiles Vakuumpotential kann jeden dieser Winkel tragen:

\[
V(x,y)=\lambda_f(x^2+y^2-v^2)^2+\kappa(bx-ay)^2.
\]

Seine Minima liegen bei ±v(a,b)/√(a²+b²), mit positivem Hesseoperator. Der entsprechende reale Massenteil liefert tan(2θ) = b/a. Für 15 Grad muss also weiterhin b/a = 1/√3 dynamisch ausgewählt werden. Stabilität allein wählt diesen Wert nicht.

Auch der eingeschränkte rationale thermische Projektoransatz wurde geprüft. In der dokumentierten Basis liegen die Einträge in Q(q) und der Nebendiagonalfaktor in √2 Q(q). Für einen nicht entarteten Block würde bei q = exp(−2π) die exakte 15 Grad Bedingung ein Quadrat gleich 24 erzwingen. Der Transzendenzbeweis und die rationale Funktionenkörperrechnung schließen eine nicht entartete, eindeutig ausgewählte 15 Grad Richtung in genau dieser Klasse aus. Der entartete Block A = C, B = 0 wählt überhaupt keine Richtung. Es ist kein Ausschluss beliebiger Thetafunktionen, anderer Koeffizientenkörper oder nichtlinearer Vakuummechanismen.

### Vergleich mit dem benannten Datensatz

TM1 bei 15 Grad mit μτ Reflexion liefert

\[
\sin^2\theta_{13}=\frac{2-\sqrt3}{12}\simeq0{,}0223291,
\quad
\sin^2\theta_{12}\simeq0{,}318107,
\quad
\sin^2\theta_{23}=\frac12,
\quad
\delta=\pm90^\circ.
\]

Für NuFIT 6.1, normale Ordnung mit IC24 und SK, liegen die drei Winkel und δ = 270 Grad in den angegebenen marginalen 3σ Bereichen. δ = 90 Grad liegt außerhalb des dort angegebenen Bereichs 125 bis 365 Grad. Dies ist ein Vergleich mit eindimensionalen Bereichen, kein gemeinsamer Gütefit und keine bereits dynamisch ausgewählte CP Phase. Quelle: [offizielle NuFIT Tabelle](https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf).

Für die direkte TM2 Klasse hilft eine größere Gruppe bezüglich des Reaktorwinkels zwar bei einzelnen Werten. Ihre exakte Solarsummenregel sin²θ₁₂ = 1/[3(1−sin²θ₁₃)] > 1/3 liegt aber bereits oberhalb der 3σ Obergrenze 0,3295 dieses Datensatzes. Die strukturelle TM2 Aussage der direkten Δ(6n²) Klasse ist in [King, Neder und Stuart](https://arxiv.org/abs/1305.3200) hergeleitet.

## 8. OriginLab und Reproduzierbarkeit

Die bestehende Lane wurde erweitert. Der Suchmotor bleibt genesis_rules.search. Die neue Untersuchung wird durch process_source_dynamics gesammelt und von process_rules.run aufgerufen.

Neu oder präzisiert sind:

- exakte Zertifikate statt endlicher Stichproben für Matrixkerne, Spektralkerne, Transferindex und HR Potentialklassen;
- getrennte Angaben für tatsächliche Pole und die historische Kostenfunktion;
- getrennte Aussagen über alle passenden Kandidaten und nur den kleinsten Kandidaten;
- korrekte Unterscheidung der mikroskopischen Vakuumsubtraktion von der Zentrierung eines 4D Feldes;
- ein Universum verschiedener positiver Elternoperatoren mit demselben Thermofeldzustand;
- alle vier gewichteten metrischen Sektoren, Quellenfaktoren, Quartik, Tiefenkontrolle und die Flavour Gegenbeweise;
- eine explizite Kennzeichnung, dass ein gemeinsames t über alle Tiefen nicht aus μ₄ folgt.

Die Schalter transfer, state und metric beschränken weiterhin die entsprechenden Kandidatenklassen. Nach den Zertifikatskorrekturen sind ihre aktuellen JSON Ausgaben keine byteidentischen Wiederholungen alter Versionen. Die historischen Ergebnisdateien und die unveränderten Originalskripte bleiben als Quellen erhalten.

Insbesondere war die Auswahl ohne die Skalenkonvention zuvor zu stark als eindeutig beschrieben: Das korrigierte Spiegelgewicht gilt auch für das Transfergedächtnis bei q = 1/√2. Im endlichen alten Universum bleiben deshalb zwei verschiedene ausgewählte Familien, solange ihre unterschiedlichen thermischen oder Spiegelannahmen nicht zusätzlich festgelegt werden. Keine davon ist dadurch eine eindeutig hergeleitete physische Prozessregel.

Die reproduzierbare Befehlsfolge und der tatsächliche Prüfstand werden in VALIDATION.md festgehalten. Der Begleitbericht enthält die vollständigen Einzelbeweise und die maschinenlesbaren Resultate.

## 9. Die noch fehlende Schließung ist jetzt als Gleichung formulierbar

Die nächste Auswahl darf nicht aus einer weiteren passend gewählten Klausel bestehen. Benötigt wird eine physisch definierte mikroskopische Quellenfunktion mit denselben Normen und Ward Identitäten wie die metrische Antwort:

\[
W_{\rm micro}^{(2)}=W_{\rm metric}^{(2)},\qquad
W_{\rm micro}^{(3)}=W_{\rm metric}^{(3)},\qquad
W_{\rm micro}^{(4)}=W_{\rm metric}^{(4)}.
\]

Die Zweipunktfunktion kann Skalen und Spektralgewicht prüfen. Die Dreipunktfunktion muss die gravitative Selbstkopplung und die Quellenidentifikation tragen. Die volle Vierpunktfunktion kann dann die verbleibenden t_n bestimmen. Der berechnete Ausdruck F(P) sagt genau, welche Information dafür fehlt und wie sie in einer benannten Operatorfamilie ausgewertet würde.

Unabhängig davon bleiben die physische Naht und ihre Kondensation, ein lokales chirales Fermionmaß in 3+1 Dimensionen, die gemeinsame Eichwirkung und die Flavour Ausrichtung offen. Der neue Stand liefert somit konkrete positive Operatoren, eine nichtlineare metrische Vergleichswirkung und ihre messbaren Unterschiede. Er liefert weiterhin keine vollständig bestimmte Theorie von allem.
