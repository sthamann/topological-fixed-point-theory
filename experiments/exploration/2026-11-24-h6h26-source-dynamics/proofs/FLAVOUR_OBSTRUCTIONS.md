# H6h26: Flavour Auswahl, zwei exakte Grenzen und die verbleibende Dynamik

Stand: 24. November 2026. Originalstand: `0172a74bdea0d3f751372813598dfe3f777a5930`.

## Ergebnis

Der TM1 Punkt mit zusätzlichem Mischungswinkel `vartheta=pi/12` und mu-tau Reflexion besitzt zwei präzise Auswahlgrenzen:

1. Er kann nicht aus einem geladenen `Z3` und dem vollständigen Majorana Kleinrest als Untergruppen **irgendeiner endlichen Flavourgruppe** folgen. Der Beweis benutzt eine nicht ganzalgebraische Gruppenspur und hängt nicht an einer endlichen Suchgrenze.
2. Eine endliche rationale Gewichtung der angegebenen kanonischen Projektoren durch `q=exp(-2 pi)` kann diesen Winkel ebenfalls nicht nichtentartet auswählen. Der Satz gilt nur im unten genannten Koeffizientenkörper und in der angegebenen realen Basis.

Der semidirekte Weg mit nur einem neutralen `Z2` und CP bleibt zulässig. Eine stabile lokale Potentialfamilie zeigt den noch freien dynamischen Richtungsparameter ausdrücklich. Ein größerer Gruppenname oder die bloße Kenntnis des thermischen `q` bestimmt diesen Parameter nicht.

Die neuen Dateien `research/originlab/process_flavour_obstructions.py` und `test_process_flavour_obstructions.py` enthalten die exakten Zertifikate, negative Kontrollen und einen separat gekennzeichneten Intervallvergleich. Elf Testgruppen bestehen. Bestehende Lane und Statusdateien wurden für diese Arbeit nicht geändert.

## 1. Was Z49 tatsächlich untersucht

Originalquellen sind [z49_flavourgruppen.py](https://github.com/sthamann/topological-fixed-point-theory/blob/0172a74bdea0d3f751372813598dfe3f777a5930/experiments/exploration/2026-11-23-phasenraum-f5-holocode/z49_flavourgruppen.py) und die Funktion `flavour_groups` in [process_rules.py](https://github.com/sthamann/topological-fixed-point-theory/blob/0172a74bdea0d3f751372813598dfe3f777a5930/research/originlab/process_rules.py).

Die untersuchten Werte sind

\[
n\in\{2,3,4,5,6,7,8,10,12,16\}.
\]

Das ist kein vollständiger Scan aller ganzen Zahlen bis 16. Bei ungeradem `n` kann `Delta(6 n²)` allerdings schon aus seiner Gruppenordnung keinen Kleinrest enthalten: Der Zweieranteil seiner Ordnung ist nur zwei. Unter den ausgelassenen geraden Werten bis 16 bleibt damit `n=14`; auch dort trifft die analytische direkte Klasse keines der hier betrachteten Reaktorfenster.

Z49 verlangt einen vollständigen Kleinrest im neutralen Sektor und einen zyklischen geladenen Rest mit drei verschiedenen Eigenwerten. Die Monomialgruppen werden mit ganzzahligen Phasenexponenten erzeugt; die gemeinsame Eigenbasis und deren Betragsquadrate werden anschließend numerisch berechnet und gerundet. Das Prädikat `mixes` bedeutet dort, dass kein Eintrag verschwindet. `theta13_window` prüft irgendeinen Eintrag der Betragsmatrix; erst `fits_data` prüft eine gemeinsame Zeilen und Spaltenzuordnung aller drei Winkel.

Die ursprüngliche Lane benutzt die gerundeten Bereiche

`sin²(theta12) in [0.27,0.35]`, `sin²(theta13) in [0.020,0.024]`, `sin²(theta23) in [0.43,0.60]`,

als ungefähre NuFIT 6.0 Bereiche. Diese sind von der bereits geprüften Originaltabelle [NuFIT 6.1, normale Ordnung mit SK](https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf) zu unterscheiden. Deren Reaktorbereich lautet `[0.02064,0.02418]`, der Solarbereich `[0.2893,0.3295]`. Alle diese Vergleiche sind marginale Intervallvergleiche und keine gemeinsame Likelihood.

### Größere n treffen das Reaktorfenster, beheben aber TM2 nicht

Für die direkte Klasse mit geladenem `Z3` und vollständigem neutralem Kleinrest existiert eine analytische Beschreibung für beliebige gerade `n`: [King, Neder und Stuart, 2013, Gleichung 20](https://arxiv.org/abs/1305.3200). Der zweite Vektor ist trimaximal, die übrigen enthalten Sinus und Kosinus von `pi gamma/n` und ihren Verschiebungen um `pi/6`. Zeilen und Spaltenzuordnungen sind zu berücksichtigen.

Die erweiterte Prüfung dieser analytischen Klasse ergibt:

| Vergleich | Erstes gerades n | Gruppenordnung | Ein Reaktorwert im Fenster |
| --- | ---: | ---: | --- |
| Ursprüngliches Z49 Fenster | 18 | 1944 | `(2/3) sin²(pi/18) = 0.0201024597...` |
| NuFIT 6.1, NO mit SK | 22 | 2904 | `(2/3) sin²(2 pi/33) = 0.0238773557...` |

Der Helper verwendet dafür Intervallarithmetik bei den exakt rationalen Winkeln. Dies ist eine Prüfung der analytisch klassifizierten direkten Familie, kein neuer vollständiger Gruppen oder Modellscan. Die größeren Reaktortreffer dürfen nicht als vollständige Fits bezeichnet werden: Für sämtliche exakt erhaltenen TM2 Muster gilt

\[
\sin^2\theta_{12}=\frac1{3(1-\sin^2\theta_{13})}>\frac13>0.3295.
\]

Damit scheitert diese direkte Klasse unabhängig von `n` bereits am aktuellen marginalen Solarbereich. Eine semidirekte `Z2+CP` Konstruktion ist davon nicht erfasst. Diese Unterscheidung ist auch in [Ding, King und Neder, 2014](https://arxiv.org/abs/1409.8005) ausdrücklich vorgenommen: Weniger als der volle Kleinrest lässt einen reellen Mischungsparameter übrig.

## 2. Delta(96) reproduziert das alte Muster mit einem anderen neutralen Rest

In der Standarddarstellung mit `eta=i` seien

\[
a=\begin{pmatrix}0&1&0\\0&0&1\\1&0&0\end{pmatrix},\qquad
b=-\begin{pmatrix}0&0&1\\0&1&0\\1&0&0\end{pmatrix},\qquad
c=\operatorname{diag}(i,-i,1).
\]

Dann sind

\[
K_1=c^2=\operatorname{diag}(-1,-1,1),\qquad
K_2=abc=\begin{pmatrix}0&i&0\\-i&0&0\\0&0&-1\end{pmatrix}
\]

kommutierende Involutionen. Beispielsweise erhält dieser Kleinrest die symmetrische Majorana Matrix

\[
M=\begin{pmatrix}-2&i&0\\i&2&0\\0&0&2\end{pmatrix}.
\]

Ihre Takagi Massen sind die drei verschiedenen Zahlen `1,2,3`. Zusammen mit dem geladenen Dreierzyklus ergibt ihre Eigenbasis, bis auf Permutationen, das Muster

\[
|U|^2=\begin{pmatrix}
a_0&1/3&c_0\\1/3&1/3&1/3\\c_0&1/3&a_0
\end{pmatrix},\qquad
a_0=\frac{2-\sqrt3}{6},\quad c_0=\frac{2+\sqrt3}{6}.
\]

Der Helper prüft die Majorana Invarianz direkt. Die volle Vierteluhr mit Charakteren `(-1,i,-i)` in dieser Eigenbasis erhält dieselbe nichtentartete Masse dagegen nicht. Z49 zeigt somit eine neue zulässige Herkunft derselben Betragsmatrix. Es widerlegt nicht die frühere Paarentartung bei einer tatsächlich ungebrochenen vollen neutralen `Z4`.

## 3. Ein gruppengrößenunabhängiger Ausschluss des vollen Kleinansatzes bei 15 Grad

Wähle im geladenen Massenrahmen

\[
T=\operatorname{diag}(1,\omega,\omega^2),\qquad \omega=e^{2\pi i/3}.
\]

Für eine normierte PMNS Spalte `u_i` ist der zugehörige Majorana Generator mit Determinante eins

\[
S_i=2u_i u_i^\dagger-I.
\]

Mu-tau Reflexion verlangt `|u_(mu i)|²=|u_(tau i)|²`. Daher folgt, ohne eine Matrixelementapproximation,

\[
\operatorname{tr}(TS_i)
=2\sum_{\alpha=e,\mu,\tau}T_{\alpha\alpha}|u_{\alpha i}|^2
=3|u_{ei}|^2-1.
\]

Für TM1 bei `vartheta=pi/12` gilt

\[
x=|U_{e3}|^2=\frac{2-\sqrt3}{12},\qquad
|U_{e2}|^2=\frac13-x,\qquad |U_{e1}|^2=\frac23.
\]

Damit sind die drei Spuren

\[
\operatorname{tr}(TS_1)=1,\qquad
\operatorname{tr}(TS_2)=-\frac{2-\sqrt3}{4},\qquad
\operatorname{tr}(TS_3)=-\frac{2+\sqrt3}{4}.
\]

Die letzten beiden haben das irreduzible primitive Minimalpolynom

\[
16z^2+16z+1,
\]

also das monische Minimalpolynom `z²+z+1/16`. Ihre Körpernorm in `Q(sqrt(3))` ist `1/16`; sie sind keine algebraischen ganzen Zahlen.

Wäre `TS_i` ein Element einer endlichen Gruppe in einer gewöhnlichen unitären Darstellung, hätten seine Eigenwerte endliche Ordnung. Sie wären Wurzeln der Einheit, und ihre Summe wäre eine algebraische ganze Zahl. Der Widerspruch beweist den Ausschluss für jede endliche Gruppengröße. Der allgemeine Einsatz von Wurzeln der Einheit bei direkten Restklassifikationen ist in [Fonseca und Grimus, 2014](https://arxiv.org/abs/1405.3678) dokumentiert; die obige kurze Rechnung ist der hier benötigte spezielle Nachweis.

### Warum der semidirekte Weg weiterhin offen bleibt

Im semidirekten Modell muss nur `S1` aus der vorgegebenen endlichen Flavourgruppe stammen. Dessen Spurzertifikat ist unauffällig. Die beiden weiteren Involutionen sind zwar algebraische Symmetrien einer gegebenen nichtentarteten Majorana Masse, müssen aber keine Elemente der fundamentalen endlichen Flavourgruppe sein. Ihre gemeinsame Erzeugung mit `T` darf dann unendlich sein, ohne den semidirekten Ansatz zu verletzen.

Die Aussage sperrt deshalb einen **vollen direkten Kleinursprung** des exakten 15 Grad Punktes mit dem genannten geladenen `Z3` und mu-tau Reflexion. Sie sperrt weder `Z2+CP` noch allgemeine dynamische Flavourmodelle.

## 4. Was die transzendente Quarterzahl im rationalen Projektoransatz leisten kann

Die relevante reale Ebene des neutralen Restes `(03)` besitzt die orthonormale Basis

\[
u=(1,-1,-1,1)/2,\qquad w=(0,1,-1,0)/\sqrt2.
\]

Sie liegt im dreidimensionalen Summennullraum der vier Marken. In der verwendeten TM1 Zuordnung misst `vartheta` die Rotation in dieser Ebene relativ zur tribimaximalen Basis.

Sei nun ein symmetrischer Operator aus endlich vielen kanonischen rationalen Projektoren und rationalen Gewichtsfunktionen konstruiert. Präzise angenommen werden Matrixeinträge in `Q(t)` vor der Kompression; etwa bleiben normierte Rang eins Projektoren auf Vektoren mit rationalen Komponenten in diesem Körper, weil sich die Normquadrate rational ausdrücken. In der oben fixierten Basis hat die Kompression die Form

\[
M_+(t)=\begin{pmatrix}A(t)&\sqrt2 B(t)\\\sqrt2 B(t)&C(t)\end{pmatrix},
\qquad A,B,C\in\mathbb Q(t).
\]

Für eine nichtentartete Eigenrichtung gilt

\[
\tan(2\vartheta)=\frac{2\sqrt2 B}{A-C}.
\]

Die notwendige Bedingung für `vartheta=pi/12` ist folglich

\[
(A-C)^2=24B^2.
\]

### Transzendenz und Einbettung

Die Quellzahl `q=exp(-2 pi)` ist nach Gelfond–Schneider transzendent: Sie ist der Wert von `(-1)^(2i)` auf dem Zweig `log(-1)=i pi`; Basis und Exponent sind algebraisch, der Exponent ist irrational. Siehe [Gelfond, Originalarbeit 1934](https://www.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=4924&wshow=paper).

Die Auswertung `Q(t) -> R` bei `t=q` ist daher eine injektive Körperabbildung. Kein nichtverschwindendes rationales Polynom kann an `q` verschwinden, insbesondere kein Nenner. Gälte die Winkelgleichung bei `q`, müsste sie schon eine Identität in `Q(t)` sein.

Aber `24` ist kein Quadrat in `Q(t)`. Äquivalent ist `X²-24` über diesem rein transzendenten Erweiterungskörper weiterhin irreduzibel. Für `B != 0` würde die Gleichung verlangen, dass `(A-C)/B` eine Quadratwurzel von 24 in `Q(t)` ist. Also bleiben nur `B=0` und `A=C`: ein skalarer Operator, der überhaupt keine Richtung auswählt.

Die gleiche Argumentation gilt bei endlich vielen rationalen Exponenten der Boltzmannzahl, wenn ein gemeinsamer Wurzelparameter `t=exp(-r pi)` mit rationalem `r != 0` benutzt wird. Dieser Parameter ist ebenfalls transzendent.

### Grenzen dieses Satzes

Die Aussage setzt die genannte Basis und den angegebenen Koeffizientenkörper voraus. Sie gilt nicht automatisch bei zusätzlichen algebraischen Koeffizienten, anderen Insertionsoperatoren, unendlichen Thetafunktionen oder einer nichtlinearen Vakuumgleichung. Insbesondere wird **keine** algebraische Unabhängigkeit von `pi` und `exp(pi)` vorausgesetzt oder behauptet. Gewichte mit unabhängig verbleibenden `pi` Koeffizienten können außerhalb von `Q(q)` liegen.

Der Helper prüft rationalen Körperbesitz ausdrücklich und weist beispielsweise ein neu eingesetztes `1/sqrt(24)` oder `sin(t)` zurück, statt den Satz dort fälschlich anzuwenden. Der bereits bekannte Dreieckenkandidat mit `(A,B,C)=(7/3,-1/3,2/3)` hat den exakten Rest `(A-C)²-24B²=1/9`, nicht null.

Die bekannte Zahl `q` kann in diesem eingeschränkten linearen Projektoransatz also Stärken und Gewichte verändern, aber den exakten 15 Grad Winkel nicht nichtentartet erzeugen. Das ist kein Ausschluss beliebiger S4 Modelle.

## 5. Die erlaubte nichtlineare Z2+CP Dynamik besitzt einen sichtbaren Freiheitsgrad

Im TM1 Rahmen hat eine passende Majorana Masse den Block

\[
M_\nu^{\rm TBM}=
\begin{pmatrix}m_1&0&0\\0&b_0&i d\\0&i d&c_0\end{pmatrix},
\qquad b_0,c_0,d\in\mathbb R.
\]

Nach der festen CP Rephasierung ist der untere Block eine beliebige reelle symmetrische Matrix. Schreibe ihn als

\[
M'=mI+g\begin{pmatrix}x&y\\y&-x\end{pmatrix}.
\]

Dann gilt `tan(2 vartheta)=y/x`. Die Masse lässt sich zurück in den genannten CP Rahmen transportieren; bei `m>|g|sqrt(x²+y²)>0` sind ihre beiden Takagi Massen verschieden.

Wenn nur der erhaltene neutrale `Z2` und diese CP Transformation vorgegeben sind, sind `x,y` zwei reelle CP gerade Massenspurionen. Beide Komponenten sind durch die Restsymmetrie zulässig. Schon die stabile Potentialfamilie

\[
V_{a,b}(x,y)=\lambda(x^2+y^2-v^2)^2+\kappa(bx-ay)^2,
\qquad \lambda,\kappa,v>0,\quad(a,b)\ne(0,0),
\]

wählt isolierte Minima

\[
(x,y)=\pm v\frac{(a,b)}{\sqrt{a^2+b^2}},\qquad
\tan(2\vartheta)=\frac ba.
\]

Der radiale Hessianeigenwert ist `8 lambda v²`; der angulare ist `2 kappa(a²+b²)`. Beide sind strikt positiv. Das verbleibende Problem ist somit nicht bloß eine flache Richtung eines unvollständigen Potentials: Verschiedene erlaubte stabile Potentiale wählen verschiedene Winkel. Für 15 Grad braucht man ausdrücklich `b/a=1/sqrt(3)`.

Dies ist eine Gegenfamilie effektiver Potentiale unter den angegebenen Restannahmen, kein Anspruch auf eine vollständige S4 invariante ultraviolette Wirkung. Eine Herkunft der beiden Richtungskoeffizienten aus dem gemeinsamen Quellfunktional könnte die Freiheit schließen. Die bloße Tatsache, dass `q` feststeht, tut dies nicht. Beispielsweise verändern positive Faktoren `lambda(q)` und `kappa(q)` nur Radial und Angularmassen, nicht die noch freie Richtung `b/a`.

## 6. Prüfung einer naheliegenden minimalen skalaren Winkelselektion

Eine echte S4 Doppeltdarstellung faktorisiert durch `S4/V4 = S3`. In passenden realen Doppeltdaten `z=x+i y` sind die grundlegenden Winkelinvarianten der dihedralen Wirkung `I2=|z|²` und `I3=Re(z³)=x(x²-3y²)`. Ein renormierbares skalares Potential enthält im Allgemeinen den freien kubischen Koeffizienten von `I3`; es erzwingt nicht die gewünschte Halbwinkelrelation.

Ein zusätzliches Vorzeichen `z -> -z` würde den kubischen Term verbieten. Dann könnte ein positiver Term `I3²` als erster Winkelterm die Richtungen `arg(z)=pi/6+k pi/3` auswählen. Eine anschließend postulierte Identifikation von `arg(z)/2` mit dem Mischungswinkel enthielte auch 15 Grad. Das ist jedoch noch keine TFPT Herleitung:

- Das zusätzliche Vorzeichen, die skalare Darstellung und die Kopplung an den neutrinen Massenblock müssten aus dem Prozess folgen.
- Die sechs Vakuumrichtungen sind nicht ohne weitere relative Ausrichtung eine einzige PMNS Vorhersage.
- Vor allem wirkt `V4` auf jeder S4 Doppeltdarstellung trivial. Enthält ein vollständig S4 invariantes Majorana Modell nur solche Doppelts und Singuletts als brechende Vakuumfelder, bleibt dieser volle Kleinrest erhalten. Die Doppeltdynamik allein kann daher nicht den oben ausgeschlossenen TM1 Punkt erzeugen. Eine neue triplette oder anderweitige Kopplungsstruktur wäre erforderlich.

Damit wird kein `cos(12 vartheta)` nach Zielwert als Erklärung eingesetzt. Die minimale tatsächlich nachgewiesene Aussage bleibt: Unter `Z2+CP` existiert eine stabile Familie mit freiem Richtungsverhältnis. Unter einer vollständigen endlichen Kleinherkunft ist der genaue 15 Grad Punkt ausgeschlossen. Ein nichtlinearer Auswahlmechanismus außerhalb dieser Grenzen bleibt eine echte nächste Aufgabe; seine erste zusätzliche Struktur ist nun genau benannt.

## Reproduktion

Vom Verzeichnis `research/originlab` aus:

```bash
python -m unittest test_process_flavour_obstructions
```

Das JSON fähige Ergebnis liefert `flavour_obstruction_audit()`. Es markiert `physical_solution=False`, lässt den semidirekten Weg ausdrücklich offen und trennt exakte Algebra von den marginalen Intervallvergleichen.
