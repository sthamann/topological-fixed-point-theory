# TFPT · Der APS Rand entscheidet über die Higgsdämpfung

10. November 2026. Auswertung der vorhandenen Spin4 × v252 Produktrealisierung für ein flaches Produkt mit drei tangentialen Richtungen und normalabhängigem Higgsfeld. Verglichen werden die beiden bereits unterschiedenen APS Vorschriften. Die Identifikation mit der ganzen physischen Quelle bleibt eine eigene Voraussetzung.

## 1. Geprüfter Originalpfeil und tatsächliche Operatoren

Die lokale, noch nicht integrierte v252 Wärmeantwort aus R031/R041 enthält für `H=r(s) e_2`, `u=r²` und flache Metrik

```text
a |r'|² + b r⁴ + d_M r² − (a/3) u'',
a = Tr96(Y²)/(8π²), b = Tr96(Y⁴)/(8π²), d_M = Tr96(Y² M²)/(2π²).
```

`Y=L(e_2)` und `M=D_M` sind unverändert. Die feste Skalarkrümmung null dient hier genau dem schon benannten Produktanschluss. Auf dem günstigsten Profil `r(s)=R/(1+kRs)` mit `k=√(b/a)` ist `u''=6 k² r⁴`. Deshalb heben sich die ersten beiden Terme mit der **Volumendivergenz** punktweise auf. Damit ist der Randbeitrag keine optionale Ergänzung der zuvor hergeleiteten kubischen Schranke.

In inward Eigenlänge `s` am linken Ende ist für jeden tangentialen Impuls `p∈R³` oder auf dem flachen Torus

```text
D_p = N (∂_s+B_p(s)),      N=−iγ_n,
B_p(s)=P(p)+J(M+r(s)Y),   J=iγ_n γ_5,
N* = −N, J*=J, J²=1, {N,B_p}=0,
B_p(s)²=|p|²+(M+r(s)Y)²,   D_p²=−∂_s²+B_p(s)²−B_p'(s).
```

Tensorprodukte mit dem 96er Raum sind unterdrückt. `P(p)` ist die masselose tangentiale Matrix, `P(p)²=|p|²` und `{P(p),J}=0`. Keine Kommutation von `M` und `Y` wird benutzt. Das alte feste APS setzt links `P_+(P(p))ψ=0`, rechts komplementär. Die mitbewegte Vorschrift setzt stattdessen die Projektion des tatsächlichen `B_p` am jeweiligen Rand ein. Der alte eindimensionale Aufruf verwendet die feste Vorschrift mit `h0=pγ_F`. Ihre Verwendung in der 4D Produktrealisierung ist genau die hier zu prüfende Fortsetzung. Die geometrische Originaldefinition allein identifiziert diese Fortsetzung nicht mit dem physischen `P_adm`.

[Volle lokale Koeffizienten R031/R041](/Users/stefanhamann/Projekte/tfpt-theoryv4/docs/research/fundamental-plan/RESULTS.md:463) · [Fester originaler APS Aufruf](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/ORIGINAL-ENERGY-SOURCE-PROOF.txt:96) · [Geometrischer Operator und gebundener Vergleich](/Users/stefanhamann/Projekte/tfpt-theoryv4/documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/001-subsection-from-bridge-to-geometric-completion.tex:30)

## 2. Ausgeführte APS Rechnung mit vollständigem normalem Jet

Friert man den Randwert ein und schreibt `B(s)=B_0+sF+O(s²)`, bleiben bei Variation der Normalableitung die Randprojektoren konstant. Die zum Quadrat gehörige Dirichlet Gruppe verschwindet am Rand. Die andere Gruppe hat die repulsive Robin Bedingung `∂_s ψ=h ψ`. Für die mitbewegte Vorschrift ist `h=E=√(|p|²+μ²)` auf einem Paar aus den Eigenwerten `±E`. Hier ist `μ` ein Eigenwert der unveränderten inneren Matrix `Φ_0=M+r(0)Y`. Jedes Paar hat zwei Spin Kopien pro innerer Eigenlinie. Bei der festen freien Vorschrift ist `h=|p|` und der Massenblock hat in jeder Gruppe verschwindende Diagonalkompression.

Für den einzelnen inneren Eigenwert sei `c` der Diagonaleintrag von `{Φ_0, Φ'_0}`. Duhamel variiert den Wärmetrace durch `−t Tr(δ(D²) e^(−tD²))`. Im mitbewegten Fall hat `δ(D²)` auf den beiden Gruppen im linearen Normaljet die Einträge `c(s−1/(2E))` (Dirichlet) und `c(s+1/(2E))` (Robin). Das folgt aus `P_± F P_± = ±P_±{B_0,F}P_±/(2E)` auf dem jeweiligen Eigenblock. Im festen Fall ist `F` außerdiagonal bezüglich der freien Gruppen, und beide diagonalen Einträge sind `c s`. Die Beiträge zu höheren Taylorjets sind für den hier betrachteten Randkoeffizienten zu klein.

Für `x=s/√t`, `q=h√t` lautet der Unterschied zum Halbwert des ungebundenen Diagonalkerns für Dirichlet `−t^(−1/2) e^(−x²)/(2√π)`, für Robin sein Gegenstück minus `t^(−1/2) q F_q(x)`, wobei

```text
F_q(x)= exp(2qx+q²) erfc(x+q).
I_0(q)=∫ F_q(x) dx = [1−exp(q²)erfc(q)]/(2q).
q I_1(q)=1/(2√π)−I_0(q)/2,   I_1(q)=∫ x F_q(x) dx.
```

Alle x Integrale gehen von null bis unendlich. Die zwei Identitäten folgen aus `F_q'=2qF_q−2exp(−x²)/√π`. Nach den x Integralen bleiben für den gemeinsamen Paartrace, vor dem Faktor `−t c` und dem Normalmaß `√t`,

```text
J_bewegt(q)=1/(4q)−1/(2√π),
J_fest(q)=−q I_1(q).
```

Der dreidimensionale tangentiale Hochimpulsfaktor liefert das radiale Maß `q² exp(−q²) dq`. Exakt gilt

```text
∫ q² exp(−q²) J_bewegt(q) dq = 0,
∫ q² exp(−q²) J_fest(q) dq   = −1/16.
```

Für den ersten Wert heben sich die beiden Beiträge `1/8` auf. Für den zweiten wird `∫ q erfc(q) dq=1/4` verwendet. Alle q Integrale gehen von null bis unendlich. Die auf dem Torus verwendete reskalierte Impulssumme hat denselben Grenzwert. Die höchstens `1/q` Singularität ist in drei tangentialen Dimensionen integrierbar. Ein fester Vierteltwist liefert keinen p Nullmodus. Das ist eine Dimensionsabhängigkeit dieses Koeffizienten und kein Selektor für die Raumzeitdimension.

Der angegebene Duhamel Wert ist der **Netto** Normaljet nach Abzug der rein inneren, ableitungsfreien ersten Variation. Beim Vergleich mit der kanonischen lokalen Bulkform muss deren `−a Δu/3` in die Trennung von Volumen und Rand eingesetzt werden. Er darf nicht ein zweites Mal auf den Netto Wert addiert werden. Mit beiden Spin Kopien und `∑c=A u'(0)`, `A=Tr96 Y²`, ergibt sich links

```text
Netto a4 Normaljet, bewegt = 0,
Netto a4 Normaljet, fest   = (a/2) u'(0).
```

Am rechten Ende gilt dieselbe Formel in der dortigen inward Koordinate, also mit `−u'(L)`. Im 4D Produkt mit normalabhängiger glatter Masse gilt die Rechnung für jeden einzelnen Profilwert, nicht nur für infinitesimal kleine Feldamplituden. Für die Begründung der Restordnung friert man in der Randzone den Wert bei `s=0` ein. Für jedes feste glatte Profil sind `δ(B²)−δB'` beschränkt und die genannte APS Bedingung bleibt beim Einfrieren dieselbe. Bei mindestens zwei Duhamel Einsätzen ist die Randkorrektur `O(t^(1/2))`: Der gebundene und der freie Kern unterscheiden sich nur innerhalb der Gauß Randzone von Dicke `O(√t)`. Dirichlet und Robin mit `h≥0` sind durch freie und Neumann Bildkerne nach oben beschränkt. Nach Integration der tangentialen Impulse hat ein n-facher Einsatz einen Randfaktor `O(t^(n−3/2))`. Die im linearen Einsatz ausgelassenen `O(s²)` Massenterme und `O(s)` Ableitungsterme tragen ebenfalls höchstens `O(√t)`; im bewegten Fall ist die diagonale Ableitungswirkung zusätzlich `O(1/|p|)`, im festen Fall null. Die inneren Duhamel Beiträge bleiben die vollen lokalen Terme aus R031/R041. Die Rechnung fordert keine in R gleichmäßige Wärmeentwicklung.

Ein konstanter Randwert bringt hier keinen weiteren Higgsbeitrag der Ordnung `t⁰` von der Randzone. Im bewegten Paar ist die integrierte Randspur `−[e^(−tE²)−erfc(E√t)]/2`. Im festen Paar ersetzt man den Robin Parameter durch `|p|` und behält `e^(−tμ²)` separat. Nach tangentialer Skalierung hat die bewegte Formel die Ordnungen `t^(−3/2)`, `t^(−1/2)` und einen Rest `O(t^(1/2)|log t|)`. Das folgt aus `erfc(√x)=(1/π)∫_1^∞e^(−x v)/(v√(v−1)) dv` und der Dreier Gaußspur. Im festen Fall liefert das separate `e^(−tμ²)` dieselbe Aussage ohne diesen logarithmischen Rest. Für den kompakten flachen Torus kann die Gaußsumme mit Poisson in demselben Integral ersetzt werden. Die nichtnullen Dualterme und der Teil `tv≥1` sind zusammen `O(√t)`. Die gegenüberliegende Fläche beeinflusst diese lokalen Randkoeffizienten bei fester Länge nur exponentiell klein.

## 3. Entscheidung für die geometrische Residuenwirkung

Für die gerade S Referenz im vorhandenen R032 Vergleich ist das Higgsprofil im Referenzbein null; die reelle geometrische Higgsdifferenz erhält den ganzen vorstehenden Term. Unter der ausdrücklichen gemeinsamen Produkt, Spur und Domainidentifikation ergibt sich auf dem Profil aus Abschnitt 1 daher für den kinetischen und quartischen Teil einschließlich beider Randenden

```text
bewegt:  +(2/3) √(ab) [R³−r(L)³],
fest:    −(1/3) √(ab) [R³−r(L)³].
```

Denn `u'(0)−u'(L)=−2k[R³−r(L)³]`. Die H abhängigen quadratischen Volumenterme einschließlich des Majoranabeitrags wachsen auf derselben Familie nur wie `O(R)`, da `∫r²=(R−r(L))/k`. Damit kann die **direkte Fortsetzung der festen freien APS Domäne in den lokalen 4D Residuenaufruf** für die benannte S Referenz und zugelassene Normalprofile die erforderliche Dämpfung nicht erzeugen. Wenn freie Normaljets bei festem Higgsrandwert zugelassen sind, ist die lokale Wirkung mit dem festen Netto Randterm sogar bei festem R nach unten unbeschränkt. Der in der vorigen Rechnung angegebene Zusatz `√δ f(s/δ)` liefert beide Vorzeichen des Normaljets bei beschränkter Volumenenergie. In der bewegten Vorschrift fällt dieser Randjet aus der reellen Wirkung, sodass auf dem festen Kragen die kubische Untergrenze aus den positiven Gradienten und Quartiktermen erhalten bleibt. Die Bewegung des Projektors für wechselnde Randwerte gehört allerdings zum ganzen Quellenaufruf.

In der T Referenz bei einer spiegelgeraden ganzen Geschichte sind Original und Referenzoperator gleich. Sind auch ihre Domains dieselben, kürzt sich die relative geometrische Antwort vollständig. Die S und T Aussagen wählen keine physische Referenz. Sie unterscheiden die beiden im Original tatsächlich vorhandenen Vergleichsverträge. Die komplexe Eta Phase wird in der reellen Dämpfung nicht als Ersatz für einen Betrag verwendet.

[Vorherige Schranke und Jetprofil](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/recursive-mirrored-source-20261103/TFPT-R113-RANDDAEMPFUNG-UND-NULLPROBE-20261110.md:35) · [Die gebundenen S und T Referenzkarten](/Users/stefanhamann/Projekte/tfpt-theoryv4/docs/research/fundamental-plan/RESULTS.md:343) · [Warum der bewegte Aufruf keine feste Randbedingung ersetzt](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/APS-DOMAIN-AND-STATE.txt:105) · [Dowker, Gilkey und Kirsten zur APS Wärmeasymptotik](https://arxiv.org/abs/hep-th/0004020)

## 4. Die volle positive Spektralwirkung liefert auf derselben variablen Profilfamilie keine kubische Dämpfung

Auf dem flachen kompakten Produkt mit der tatsächlichen unveränderten inneren Matrix lässt sich das ganze Spektrum **ohne Modenschnitt** nach oben in der Wärmespur kontrollieren. Das entscheidet den direkten Übertrag von Abschnitt 3 auf den positiven Laplaceprofilaufruf.

Schreibe `x=s+ε`, `ε=1/(kR)`, also `r=1/(kx)`. Für jede tangentiale Mode ist die quadratische Form von `D_p²` die Norm `|| (∂+B_p)ψ ||²` auf der jeweiligen APS Formdomain. Nach partieller Integration besteht sie aus `|p|²`, der Form von `−∂²+Φ²−JΦ'` und den beiden Randbeiträgen. Bei bewegter APS sind diese Randbeiträge `〈ψ,|B_p|ψ〉`, und `|B_p|≥|Φ|`. Daher dominieren sie die entsprechenden Normalbeiträge von `|| (∂+JΦ)ψ ||²`. Bei fester freier APS ist die Diagonalkompression von `JΦ` auf den zugelassenen Randkomponenten null, und die tatsächlichen Beiträge sind `|p| ||ψ||²`. In beiden Fällen folgt exakt `D_p² ≥ |p|²`.

Beschränkt man die Formdomain auf Dirichlet an beiden Enden, hat die Einschränkung pro p höchstens endliche Kodimension `2d` mit `d=384`. Im Inneren ist der führende singuläre Term

```text
( (Y/k)² + J(Y/k) ) / x².
```

Seine kleinste mögliche Eigenzahl ist `−1/4+δ` mit `δ=min_j (|y_j|/k−1/2)²`, wobei `y_j` die tatsächlichen Eigenwerte von Y sind. Die exakte Reduktion von `det(4A Y²−B I)` für die vorhandenen dyadischen Matrixeingaben modulo `65537` ergibt **1763**; `i` wird auf 256 abgebildet. Daher ist `δ>0` exakt, ohne Annahme einer numerischen Schwelle. Die numerische Distanz zum kritischen Verhältnis beträgt rund 0,06475. Der Majoranaterm `r{M,Y}` wird durch `−C/x` abgeschätzt. Mit `δ_0=min(δ,1/8)` gilt `C/x≤(δ_0/2)/x²+C²/(2δ_0)`. Die Dirichlet Hardy Ungleichung am linken Ende `∫|ψ|²/x²≤4∫|ψ'|²` liefert damit für die Dirichlet Form

```text
D²_D,p ≥ |p|² + 2δ_0 (−∂²_D) − C²/(2δ_0).
```

Minmax und die endliche Kodimension geben für die APS Form, unabhängig von R,

```text
Tr exp(−t D²_APS)
 ≤ Σ_p exp(−t |p|²) · [2d+d exp(t C²/(2δ_0)) Σ_{n≥1} exp(−2tδ_0 (πn/L)²)].
```

Für `0<t≤1` ist das `O(t^(−2))` mit R unabhängiger Konstante. Für `t≥1` genügt die Positivität und die Spur bei `t=1`. Ein positives Laplaceprofil mit `∫(1+s^(−2))dν(s)<∞` hat damit für feste Skala auf dieser ganzen variablen Familie gleichmäßig beschränkte Spektralspur. Der Faktor `exp(−Tr f(D²/χ²))` ist nach unten durch eine positive, von R unabhängige Konstante beschränkt. Eine feste freie S Referenz ändert dies nur durch ihren festen Faktor. Damit ist der aus der lokalen Residuenwirkung gewonnene kubische Abfall **kein globaler Großfeldterm des vollen positiven Profilgewichts**, auch wenn die Randwahl die Residuenwirkung stabilisiert. Dieses Resultat geht über den früheren konstanten Higgsvergleich hinaus. Es betrifft die normale Profilfamilie auf dem festen flachen Produkt und die genannte positive Profilklasse. Es entscheidet weder einen gemeinsamen χ und R Grenzübergang noch die Stütze oder Dichte der tatsächlichen Konfigurationsbasis.

Falls die normale Variable in der tatsächlichen Feldkarte zugleich die Variable des R113 Anregungsintegrals ist, wächst dessen bekannte positive große Feldantwort auf genau der Profilfamilie im Logarithmus weiterhin wie `R log R`, da `∫r² log r ds=[r log r−r]_{r(L)}^R/k` und `∫r²=O(R)`. Die Identifikation dieser Variablen ist hier nicht benutzt, um die Wärmeschranke zu beweisen. Ohne die tatsächliche Konfigurationsbasis folgt aus einer Familie von Feldgeschichten kein divergierendes Feldintegral.

[Älterer Vergleich von vollem Profil und lokaler Quartik](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105K/fullprofile/WHOLE-PROFILE-TO-QUARTIC-TRANSFER-AND-DOMAIN-20261105.txt:13) · [Ursprüngliche Argumente, Profil und Konfigurationsbasis](/Users/stefanhamann/Projekte/tfpt-theoryv4/documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex:136)

## 5. Die stabilere Randwahl ist im selben reinen Fockzustand keine Rahmenumbenennung

Verglichen werden bei fester Feldidentifikation die besetzten negativen Projektionen `P_0` der freien tangentialen Matrix und `P_Φ` der vollständigen tangentialen Matrix am gleichen Rand. Für jeden Impuls p und jeden tatsächlichen inneren Eigenwert μ gilt im Spin Viererblock

```text
Tr_spin |P_Φ−P_0|² = 2[1−|p|/√(|p|²+μ²)].
```

Das folgt aus `(P_Φ−P_0)²=(1/2)(1−|p|/E) I_4` und der Antikommutation der Massenrichtung mit dem freien Symbol. Daher ist die Summe pro p `Tr96(Φ²)/|p|²+O(|p|^−4)`. Auf dem drei dimensionalen Produkt mit Impulsabstand eins und festem Vierteltwist liefert der kugelförmige Impulsschnitt bis `Λ` den führenden Wert `4π Λ Tr96(Φ²)`; der Rest ist `o(Λ)`. Die feste Modenprobe bis 64 ergibt 0,950635444 des führenden Terms. Die Divergenz folgt analytisch aus der Modenzählung. Für zwei verschiedene volle Randmassen `Φ,Ψ` liefert die Entwicklung `P_Φ−P_Ψ=−J(Φ−Ψ)/(2|p|)+O(|p|^−2)` den entsprechend führenden Koeffizienten `Tr96[(Φ−Ψ)²]`. Dies erfordert keine Kommutation der beiden Randmassen. Die Aussage vergleicht beide im selben Rahmen. Wird eine physische Eichabbildung tatsächlich mitgeführt, ist die Differenz im transportierten Rahmen zu bilden.

Folglich sind die beiden reinen quasifreien Darstellungen, wenn ihre Kovarianzen gerade diese Projektionen sind, auf der fest identifizierten vollen CAR Algebra nicht quasiäquivalent. Der Hilbert Schmidt Test für Projektionen ist genau der Test der Differenz ihrer Wurzeln. Eine direkte Summe der gleichartigen Flächenbeiträge kann diese Normquadrate nicht aufheben. Diese Aussage betrifft die reinen Projektorvakuumdarstellungen. Sie identifiziert keine integrierte oder algebraisch neu transportierte TFPT Referenz mit einer davon. Die früher geprüfte zweidimensionale Pin Quelle hatte eine logarithmische Divergenz. Hier wird die lineare Divergenz auf dem ausdrücklich vierdimensionalen Produkt für die volle v252 Randmatrix ausgewertet.

Auf endlichem Schnitt ist außerdem das Quadrat der Vakuumüberlappung für gewöhnliche CAR pro p das Produkt `∏_μ [(1+|p|/E_μ)/2]²`. Beide negativen Spinrichtungen pro μ tragen denselben Hauptwinkel. Der logarithmische führende Abfall pro p ist `−Tr96 Φ²/(2|p|²)`. Gleiche Spur von `Φ²` kann diese skalare führende Divergenz im Quotienten aufheben. Die Darstellungen zu zwei verschiedenen Massen bleiben im festen Rahmen trotzdem getrennt. Ein absoluter Referenzskalar kann daher den fehlenden Feld und Zustandsübergang nicht ersetzen.

[Vorhandener fester und bewegter APS Vergleich](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/APS-DOMAIN-AND-STATE.txt:105) · [Araki zur Quasiäquivalenz der CAR Zustände](https://ems.press/journals/prims/articles/2511)

## 6. Zweiter Lösungsversuch über den vorhandenen Exterior Randfaktor

Der ausgeführte alte feste Cap enthält pro besetztem Paar exakt `C_l=a_p+(1−a_p)e^(−2lE)` mit `a_p=(1+|p|/E)/2`. Im konstanten tangentialen Produkt mit den tatsächlichen 96 Massen gilt die gleiche Formel auf den zwei besetzten Spinrichtungen jeder inneren Eigenlinie. Der Zeitterm gegenüber `a_p` ist im logarithmischen Quotienten für festes `l>0` durch `2 e^(−2l|p|)` beschränkt. Er verändert keinen kubischen Großfeldkoeffizienten. Die angegebene Winkelrechnung lässt sich als reine Überlappung auf demselben einzelnen Rand lesen. Eine Identifikation dieser Überlappung mit der integrierten Randdichte einer variablen Geschichte ist nicht angenommen.

Auf dem ursprünglichen Kreis ist die gewöhnliche Determinante der Korrektur wohldefiniert. Auf dem hier untersuchten 4D Produkt ist der führende logarithmische Winkelterm pro innerem Eigenwert und p `−μ²/(2|p|²)`, mit der linearen Divergenz aus Abschnitt 5. Sowohl die minimale Subtraktion genau dieses Terms als auch die Carleman Fredholm Determinante `det₂` sind in drei tangentialen Richtungen endlich. Die erste fügt pro Eigenwert und p `μ²/(2|p|²)` hinzu. Die zweite fügt `2(1−C_l)` hinzu und behält nur `log C_l+1−C_l` pro Spinrichtung. Das sind mathematische Auswertungen **zweier verschiedener** Fortsetzungen des Randfaktors. Die bereits bekannte zentrierte CAR Quelle verwendet `det₂` zusammen mit einer bestimmten Spurkorrektur. Ihre Gleichheit mit einer der beiden hier geprüften Cap Subtraktionen müsste aus dem ursprünglichen Feldaufruf folgen.

Die führenden Großfeldkoeffizienten sind für das unveränderte `M+RY` exakt. Setze `y=|p|/|μ|` und `q=y/√(1+y²)`. Die beiden radialen dimensionslosen Integrale sind

```text
I_min = ∫_0^∞ y² [2 log((1+q)/2)+1/(2y²)] dy = π/3−4/9 > 0,
I_2   = ∫_0^∞ y² [2 log((1+q)/2)+1−q] dy     = π/3−10/9 < 0.
```

Für den ersten Wert differenziert man vor der Skalierung das Radialintegral nach `μ²`. Bei `y=sinh u` bleibt `|μ| [1/3+(π/2−1)]`. Integration in `μ²` gibt den angegebenen Wert. Die Differenz der beiden Ableitungen reduziert sich auf `(|μ|/2)∫[y³/(1+y²)^(3/2)−1]dy=−|μ|`. Damit folgt der zweite Wert. Die Diskretisierung auf dem flachen Torus gibt dieselben führenden Riemannintegrale, weil die Integranden als radiale Funktionen am Ursprung höchstens wie `1/y²` und am unendlichen Ende wie `1/y⁴` wachsen beziehungsweise abfallen. Die Beiträge außerhalb kompakter radialer Intervalle lassen sich mit diesen Potenzschranken gleichmäßig kontrollieren. Kleine oder null werdende Eigenwerte tragen im normierten kubischen Grenzwert null. `R^-1(M+RY)→Y` genügt für den Grenzkoeffizienten; kein Austausch nichtkommutierender Massenmatrizen ist nötig.

Pro Einheit Randvolumen ist der Koeffizient im Loggewicht folglich `I·Tr96|Y|³/(2π²)`. Beim tatsächlichen v252 Input ergeben sich **+32,1347686635** bei minimalem Abzug und **−3,4074440087** bei `det₂`. Der minimale Abzug verstärkt. `det₂` würde als einzelner Faktor dämpfen. Selbst wenn man den separaten `det₂` Faktor ohne Doppelerfassung zum festen lokalen Residuenvertrag addieren dürfte, reicht er dort nicht. Aus `Tr|Y|³≤√(AB)` folgt allgemein

```text
(10/9−π/3) Tr|Y|³/(2π²) < √(AB)/(24π²) = √(ab)/3.
```

Für den ausgeführten Input bleibt in dieser bedingten Kombination ein positives kubisches Loggewicht von **+1,3314451927**. Für den vollen positiven Profilaufruf wäre die nackte `det₂` Dämpfung prinzipiell groß genug gegen die bekannte subkubische Verstärkung, aber gerade ihre finite Spurkorrektur ist physikalisch zu begründen.

Auch der Kompositionsvertrag entscheidet auf derselben Größenordnung. Für die bloße Produktidentität zweier gleicher Cap Faktoren gilt `log det₂(C_l²)−2 log det₂(C_l)=−Tr(C_l−I)²`. Der kubische Koeffizient dieses Spurterms pro Fläche ist `(π/4−2/3)Tr|Y|³/(2π²)`, im tatsächlichen Input **+6,3299701615**. Das ist der Betrag der logarithmischen Produktkorrektur. Das Integral folgt durch Integration von `y²(1−q)²/2`. Die wirkliche unprojizierte Verklebung ist separat: `C_(2l)−C_l²=a_p(1−a_p)(1−e^(−2lE))²`. Sie enthält den zwischenzeitlichen Rückweg. Damit ist die frei gewählte, separat renormierte Cap Multiplikation auch mit einer formal endlichen Determinante kein Ersatz für die Originalverklebung.

[Ursprünglicher fester Exterior Cap](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/FIXED-CAP-DETERMINANT.txt:9) · [Zentrierte CAR Quelle samt Spurkorrektur](/Users/stefanhamann/Projekte/tfpt-theoryv4/docs/research/fundamental-plan/RESULTS.md:1346)

## 7. Tatsächlich ausgeführter direkter Randzugang

Als Alternative zur Fortsetzung über das Volumengewicht wurde der schon vorhandene ganze R032 Übergang der reflektierten v252 Historie erneut an drei Schnittpositionen verklebt. Die relative leere Fermionspur gegen S beträgt `6,088338720414977·10^−8`; ihr logarithmischer metrischer erster Jet ist `−5,010567639853645` und der zweite `2,76155191476368`. Die direkte unabhängige Differenzrechnung bestätigt diese Werte. Gegen T sind alle drei null, weil dieser bedingte Referenzaufruf den ganzen gleichen Verlauf mitführt. Eine vorzeitige Ladungskompression verliert einen echten Rückweg mit Operatornorm `0,24511333989136047`. Diese Zahlen betreffen eine vorgegebene Geschichte und ihre ursprüngliche gebundene Fermionreferenz.

Der alternative **direkte** unnormierte Endpunktaufruf und der Faseraufruf wurden anschließend auf demselben zusammengesetzten Originaldatensatz `original-variable-H-quarter-full96-diagnostic-v1`, Gene `[2,1,1]`, ausgeführt. Die fünf kinematischen Residuen der tatsächlichen Matrix sind null. Der direkte Aufruf meldet `original_full_source_operation_not_supplied` mit `evaluate_original_full_blocked_endpoint_source`. Der Faseraufruf meldet `original_fiber_operation_not_supplied` mit `original_blocked_fiber_packet`. Alle vier `Z` und vier `W` sind null im Datensinn `None`. Das sind keine verschwindenden physikalischen Antworten. Beide vorhandenen Rechenpfade sind Verbraucher der fehlenden ursprünglichen markierten Geometrie und Bosonenoperation. Die bereits funktionierende Verklebung der gegebenen Geschichte liefert die absoluten Gewichte der Geschichten nicht. Der bekannte direkte Rand oder IR Weg ist damit in seinem vorhandenen Originalaufruf versucht, ohne eine neue unabhängige Raumzeit oder Feldbasis einzuführen.

[Unveränderter ganzer Historienaufruf](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/chat-triangulation-20261108/check_original_whole_two_strips_20261109.py) · [Direkter Endpunktaufruf und Faserverbraucher](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime/original_source_adapter.py:279)

## 8. Verifikation und Folgerung für die gemeinsame Quelle

Der beigefügte Checker verwendet den unveränderten v252 Builder mit den vorhandenen R032 Eingaben. 28 direkte Kontrollen bestanden. Dazu gehören das volle tangentiale 384er Spektrum, die Polarisation der festen Randgruppe, zwei getrennte Robin Quadraturen, der Grenzwert des gesamten variablen Jetbeitrags, die ganze Projektordifferenz samt Modensumme, der exakte modulare Nachweis und die Quadraturen für beide Cap Abzüge samt Produktkorrektur. Die Koeffizienten in dieser Diagnosenormierung sind `a=4.63844220424`, `b=43.5736027032`, mit Netto Kubik `+9.47777840282` beziehungsweise `−4.73888920141`. Die allgemeine Wärmeschranke ist analytisch und wurde nicht aus einem Modenschnitt extrapoliert.

Für die gemeinsame Quelle ist der Pfad nun entschieden, soweit diese konkreten Konstruktionen reichen. Die alte feste APS Fortsetzung liefert in diesem 4D Produkt nicht die benötigte Dämpfung aus dem lokalen Residuum. Der bewegte Fall liefert sie nur im lokalen Residuenvertrag. Im vollen Profil muss die ursprüngliche Konfigurationsbasis oder ein anderer tatsächlich vorhandener Quellanteil die Masse kontrollieren. Die historische Originaldefinition enthält die Konfigurationsbasis als Eingabe ihrer Finitevolumenformel; die untersuchte Native Laufzeit liest fertige Beitragsformen und ihre Maßmassen. Die vorliegenden Werte liefern keine solche Maßmasse. Der ursprüngliche IR oder Randweg bleibt ein möglicher Zugang, wenn er selbst die ganzen markierten Beiträge oder einen quellenidentischen unnormierten Übergang liefert. Der ganze unmarkierte Transfer kann die eingesetzten Feldmarken nicht ersetzen. Die bestehende Faltung und der ursprüngliche sichtbare Zustand sind gemeinsam mit den physischen Higgs, Eich und Metrikmarken zu rekonstruieren. Beim Transport einer bewegten Randwahl muss zusätzlich die algebraische Referenz oder die tatsächliche relative Zustandskonstruktion aus dem Original folgen. Auch dafür werden dieselben Originalzeiten, die wirkliche Auswahl, die gebundene Referenz und die beiden geordneten Normen benötigt. F1 und die physische Gesamtherleitung bleiben offen.

[Native Beitragsfunktion](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime/original_source_adapter.py:366) · [Ursprüngliche R113 leere Randmatrix und gemeinsame Marken](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/chat-triangulation-20261108/R113-ORIGINAL-UND-REFERENZ-20261109.txt:7)
