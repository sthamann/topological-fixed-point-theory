# TFPT · Das Spektralgewicht an der Spiegelnaht

10. November 2026. Der flache Spin4 Produktzweig mit unverändertem v252 Higgs und Majorana, symmetrischem Kragen und APS Projektoren wird weiter ausgewertet. Die Frage ist, ob die bereits bekannte relative Spektralwirkung zusammen mit dem gewöhnlichen gepaarten Fermioncap die benötigte positive Verklebung selbst liefert. Geprüft werden die beiden üblichen Argumentkarten des Exponentialprofils, die gebundenen R032 Referenzen und ein ausdrücklich feldunabhängiger Referenzvergleich. Die physischen Identifikationen des gesamten `S_adm`, der ursprünglichen Zeit und des `P_adm` bleiben eigene Voraussetzungen.

## 1. Eine notwendige Form am selben Schnitt

Sei `r=x u+y v` mit `u` auf der positiven Hälfte und `v(s)=u(−s)` auf der negativen Hälfte. Beide Variationen bleiben von den äußeren Rändern und dem Schnitt getrennt. Für einen positiven ganzen Kern `K(x,y)` müssen die vier Einträge für `x,y∈{0,ε}` eine positiv semidefinite Matrix bilden. Bei positivem Basiswert folgt `∂x∂y log|K|(0,0)≥0`. Jeder Faktor, der getrennt auf nur einer Hälfte wirkt, fällt in diesem Test heraus. Für den reellen geometrischen Faktor `K=e^(−S)` steht hier das Negative der gemischten Ableitung von `S`. Eine Phase ändert diese notwendige Betragsbedingung nicht. Zur allgemeinen Reflexionspositivität siehe [Jaffe und Janssens](https://arxiv.org/abs/1607.07126).

## 2. Unveränderte Masse und gemischte Wärmespur

Auf der Normalgeraden ist das Quadrat je Tangentialimpuls `D_p²=−∂²+|p|²+(M+rY)²−J_m Y r'`. Hier gilt `J_m²=1`, `Tr_spin J_m=0`. Die beiden Insertionsmatrizen sind `{M,Y}u−J_m Y u'` und ihr gespiegeltes Gegenstück. Die zweite lokale Variation `2Y²uv` verschwindet bei disjunkten Trägern. Duhamel und zweimaliges partielles Integrieren liefern die gemischte Wärmespur als Integral gegen `u(s)v(t)` mit Kernel `4T e^(−T|p|²) Ψ_T(|s−t|)`. Für den üblichen nichtlokalen Wärmeformalismus siehe [Iochum und Vassilevich](https://doi.org/10.1063/1.4774116).

Der Originaloperator `M` trägt nur sechs sterile Richtungen. `Y` hat zwischen diesen sechs exakt den Nullblock. Deshalb treten nur Paare ohne Majoranamasse oder Paare aus einer Nullrichtung und einer Masse `m` auf. `w_m` bezeichnet die Summe von `|Y_ij|²` mit beiden Einsatzordnungen. Mit der Faltung `L_m(T,d)=∫₀ᵀ g_s(d) g_(T−s)(d)e^(−m²s) ds`, `g_s(d)=exp(−d²/(4s))/√(4πs)`, gilt für `d>0`

```text
Ψ_T(d) = Σ_m w_m (m² − ∂d²) L_m(T,d).
Laplace_T L_m = exp[−d(√z+√(z+m²))]/(4√z√(z+m²)).
```

Setze `l_d(T)=d exp(−d²/(4T))/(2√π T^(3/2))`, `G_d=l_d*(e^(−m²·)l_d)` und `F_m(T)=e^(−m²T/2) I₀(m²T/2)`. Dann gilt für positive `d,T`

```text
(∂d²−m²)L_m = G_d + (F'_m * G_d)/2.
```

Dies folgt aus `Laplace_T (∂d²−m²)L_m = (1+√z/√(z+m²)) Laplace_T G_d/2`. `G_d≤l_(2d)` und seine Masse ist `e^(−d|m|)`. Bei `d→0+` nähert sich `G_d` der Delta Masse am Ursprung. Deshalb geht `(m²−∂d²)L_m` bei festem `T>0` gegen `−F'_m(T)/2`, positiv für jede echte Masse. Der masselose Beitrag geht gegen null. Die tatsächlich sechs schweren Richtungen und ihre Yukawakopplungen geben `Ψ_1(0+)=6.99796918634`; unabhängige Quadratur am getrennten Abstand `0.01` ergibt `5.04874597889`. Ein äußeres positives Profilmaß mit den gewöhnlichen endlichen inversen Wärmemomenten lässt den Grenzübergang zu. Dafür genügen `G_d≤C/T`, `|F'_m|≤m²/2` und die positive Tangentialspur. Bei feldunabhängigem Referenzanteil verletzt der geometrische Faktor allein die notwendige Betragsbedingung an genügend nahen, getrennten Varianten, sofern das Profil positive Masse bei `T>0` hat.

Für `T=1` ist das Ergebnis auch am endlichen symmetrischen Kragen `[-8,8]` mit dem an den äußeren Flächen am gemeinsamen Hintergrund eingefrorenen massiven APS Projektor nachweisbar. Der mittlere Eigenwert von `M²` liegt exakt in `(6/5,3/2)`, sein Gewicht ist exakt größer als `40` und `Tr Y²<367`. Das liefern rationale Eigenwertintervalle und der Spektralprojektor als Polynom in `M²`. Für `0<d≤0.01` gilt

```text
Ψ_1(d) ≥ −367 l_(2d)(1) + 40 · (3/10)e^(−3/2) · (e^(−d√(3/2))−erf(d)) > 0.49.
```

Die Abschätzung nutzt `−F'_m(1)≥(m²/2)e^(−m²)` und `∫₀¹G_d≥e^(−d|m|)−erf(d)`. Die beiden kompakten glatten Variationen können im positiven beziehungsweise negativen Intervall zwischen den absoluten Abständen `0.002` und `0.005` gewählt werden, mit Integral eins und totaler Variation höchstens `4000`. Bei `H=0` zerfällt der tatsächliche Randwärmekern in skalare Dirichlet/Robin Kanäle mit Robinparameter `√(|p|²+m²)≥0`. Deren Kerne liegen zwischen dem Dirichlet und dem Neumann Kern des Intervalls. Auf diesen Trägern ist die Differenz dieser Kerne für `0<s≤1` durch `100 s^(−1/2)exp(−60/s)` nach oben beschränkt. Dies folgt aus der Gaußbildsumme, weil die erste Randspiegelung mindestens `15.99` entfernt ist. Duhamel mit den begrenzten L1 Normen der Insertionspotentiale erhält deshalb die strikte positive gemischte Ableitung. Die Fehlergrenze ist pro `p` kleiner als `10^(-8)e^(−|p|²)`, während die untere Schranke des Volumenwerts `4·0.49 e^(−|p|²)` ist. Der Vergleich summiert über alle Tangentialmoden.

## 3. Der bekannte Fermioncap und die Prüfung beider Exponentialprofile

Am selben konstanten Hintergrund nimmt der gewöhnliche feste Randvektor des bewegten APS Zweigs die negativen Eigenrichtungen von `h_p` ein. In erster Ordnung bewegt eine Feldvariation `u(s)` den normierten Exterior Randvektor senkrecht zu sich selbst mit Matrixelementen `Y_tan,+− ∫₀^∞u(s)e^(−(E_+−E_−)s)ds`. `Y_tan=−iΓ_F Y`, und `E_−<0<E_+` sind die tatsächlichen vorzeichenbehafteten Einteilchenenergien. Dies ist die Duhamel Ableitung des ganzen Halbtransfers mit der am äußeren Rand unveränderten Präparation. Der gemischte Logarithmus der reflektierten Cap Paarung ist die Gramform dieser Tangenten. Für das Spiegelpaar ist er doppelt so groß. Bei `|p|=1/4` treten exakt vier verschiedene Raten auf, nämlich `0.5, 0.78442609, 1.42379103, 4.42134315`. Ihr gemischter Gesamtwert direkt am Schnitt je quadrierter Profilfläche ist `335.46507075`. Separate Seefaktoren beeinflussen die gemischte Logarithmusform nicht.

Für einen endlichen Tangentialblock eliminiert `u=∏_j(∂_s−λ_j)g` alle diese Raten exakt. Die geometrische Form auf den verbleibenden Variationen hat bei einem einzelnen Wärmeprofil den Kern `∏_j(∂_d+λ_j)^2 Ψ_T(d)`. Er fällt bei großen Abständen wie ein negatives Polynom mal `exp(−d²/T)`. Das folgt direkt durch Ableiten der Gaußintegrale und den Sattel bei zwei gleichen Wärmestrecken. Der negative Kern kann damit nicht überall eine negativ semidefinite Form liefern. Für zwei verschobene schmale Profile ist das Verhältnis der beiden Diagonalen zum Quadrat des Nebeneintrags im Limes `exp(−2a²/T)<1`, bei Verschiebung `a`. Am tatsächlichen Viertelmodus und `T=1` ist das Verhältnis bei Abständen `10,10.75,11.5` bereits `0.30580304`. Glatte kompakte Profile und ein ausreichend weiter APS Kragen erhalten das strikte Ergebnis. Ein positives Wärmeprofil mit endlichem oberem Wärmeträger hat im selben Argument eine überexponentiell fallende nichtverschwindende Form. Eine positive Hankelform müsste durch ihre Zweierminoren logarithmisch konvex sein und kann dann nicht so fallen. Diese Fortsetzung bleibt auf endliche Tangentialblöcke begrenzt.

Für die zweite v259 Lesart, den geraden ganzen Operatorfaktor `exp(−a|D|)` mit `a>0`, ist der Wärmeträger unbeschränkt. Die Antwort lässt sich hier vollständig als Frequenzdichte ausrechnen. Schreibe `λ_±=√(|p|²+m²)±|p|` und `Δ=(λ²−m²)²−4λ²|p|²`. Die Laplacedichte des gemischten geometrischen **Loggewichts** je `p,m` lautet

```text
ρ_(p,m)(λ) = (2a/π) w_m (λ²−m²) cos(a√Δ/(2λ))/√Δ  · 1_(λ>λ_+).
```

Zur Herleitung setze `A=|p|²+m²u`, `b=u(1−u)`. Das Integral von `T e^(−T|p|²)L_m` gegen das Poisson Unterordnungsmaß ist

```text
W(d) = (a/(8π)) ∫₀¹ exp[−√(A/b)√(d²+a² b)] du/√(A b).
```

Der inverse Laplacetransformierte von `exp(−k√(d²+c²))` ist `δ_k−kc J₁(c√(λ²−k²))/√(λ²−k²)` oberhalb `k`. Das folgt durch Ableiten in `k` aus der `J₀` Identität für `exp(−k√(d²+c²))/√(d²+c²)`. Mit `k=√(A/b)` und `c=a√b` liegen die Integrationsgrenzen an den beiden Nullstellen von `λ²b−A`. Das Besselintegral wird durch `∫₀^π J₁(z sinθ)dθ=2(1−cos z)/z` exakt auf den Kosinus reduziert. Die Identität folgt auch termweise aus der Potenzreihe von `J₁`. Anwendung von `4(∂d²−m²)` liefert `ρ`. Direkte Integrale über `u` und unabhängige Rücktransformationen über `λ` stimmen bis `2.33e−15` überein.

Bei einem endlichen Tangentialblock haben alle Summanden für hinreichend große `λ` dieselbe führende Kosinusphase `aλ/2`. Die summierte Dichte hat offene negative Intervalle. Bei `a=1`, dem ganzen Tangentialtorus mit Impulsen `(n+1/4,j,k)` und `λ=8` ist sie auch **ohne Modenabschneidung** strikt negativ. Nur `260` Impulse können an diesem Frequenzpunkt beitragen. Rationale Eigenwertintervalle und gerichtete Intervallrechnung geben `[-24426.655592,−24426.655146]` und schließen alle Schwellen an dieser Stelle aus. Die Gegenviertelkarte hat dieselben Beträge. Der bekannte Cap kann diese offene negative Dichte mit seinem reinen Punktspektrum nicht aufheben. Die Laplaceform kann durch glatte Halbvariationen auf diesem Intervall getestet werden. Man nimmt eine feste schmale glatte Funktion und bildet Linearkombinationen ihrer positiven Verschiebungen. Ihre Laplacetransformen sind das gemeinsame positive Profil mal Polynome in `e^(−hλ)`. Polynome approximieren auf `[0,1]` stetige lokalisierte Funktionen und die exponentiell gewichtete signierte Dichte und Cap Masse sind endlich. Das liefert eine endliche glatte Variation mit negativem Gesamtjet. Für einen einzelnen ursprünglichen Viertelmodus ist zusätzlich eine explizite Ableitungsprobe mit vier ausgeschalteten Cap Raten und `32` weiteren Ableitungen berechnet. Ihr skaliertes geometrisches Ergebnis ist `−1154.89382552`.

Für jede solche feste kompakte Halbvariation folgt die Übertragung auf einen genügend weiten endlichen symmetrischen APS Kragen durch Konvergenz der Dirichlet und Neumann Wärmekernschranken. Der Tangentialtorus mit festem Viertel hat eine positive erste Energie. Die Duhamel Integrale sind nach den skalar dominierenden Neumann Kernen mit den tatsächlichen kompakten Insertionspotentialen gleichmäßig kontrolliert. Das Poissonmaß fällt am Nullende exponentiell und am großen Ende sichert die tangentiale Energielücke die Integrabilität. Der gemischte Cap Beitrag ist für Variationen in positivem Abstand zur Naht über alle Moden absolut summierbar. Diese Aussage betrifft seinen gemischten Logarithmus auch dann, wenn ein ganzer skalarer Cap Grenzwert noch renormiert werden muss.

## 4. Mitbewegte R032 Referenzen statt eingefrorenem Nenner

Der vorige Test mit feldunabhängiger Referenz setzt nicht geometrische Deckspiegelung und OS Halbierung gleich. Wenn man beides dagegen im gegebenen R032 Chart identifiziert und unabhängige Halbprofile zulässt, bewegt sich auch das Referenzbein. Für `H=xu+yv` ist im S Zweig der Referenzhiggs `(x−y)(u−v)/2` ohne Majorana, im T Zweig `(x+y)(u+v)/2` mit Majorana. Sei `Q_M` die vollständige quadratische Variation der ganzen Spektralspur am konstanten Majoranahintergrund. Setze `L_M=Q_M(u,u)=Q_M(v,v)`, `X_M=Q_M(u,v)` und verwende Nullindices für die freie Referenz. Dann sind die relativen gemischten Wirkungen exakt

```text
X_rel,S = X_M + (L_0−X_0)/2,
X_rel,T = (X_M−L_M)/2.
```

In diesen Formeln bleibt die lokale zweite Potentialvariation erhalten. Für Wärmezeit `T`, normalem Wellenvektor `k` und ein Null/Masse Paar ist der vollständige Multiplikator bis auf den gemeinsamen positiven Vorfaktor

```text
T(m²+k²) ∫₀¹ exp[−T(m²u+k²u(1−u))]du − 1 − exp(−Tm²).
```

Bei `k=0` ist er exakt `−2e^(−Tm²)`. Für die freie Referenz ist er bei großem `Tk²` gleich `4/(Tk²)+O((Tk²)^(-2))` und positiv. Beide Vorzeichen bleiben für das ganze Poissonprofil erhalten, bei hohen `k` durch dominierte Konvergenz nach Multiplikation mit `k²`. Die Zahlkontrolle des Viertelmodus gibt im freien Bein pro Gewicht `−1.96270735` bei `k=0` und `+0.02015529` bei `k=16`. Die ganze Majoranaspur gibt bei `k=0` `−597.15569633`. Für breit eingehüllte langsame Profile im T Zweig oder hochfrequente Profile im S Zweig kann der Abstand beider Hälften von der Naht groß gewählt werden. Die Kreuzwerte `X_M,X_0` verschwinden dann im Limes, genauso die Cap Tangenten durch die tangentiale Lücke. Der verbleibende relative gemischte Wirkungswert ist in beiden Versuchen strikt positiv. Auch die beiden tatsächlich mitbewegten Referenzformen geben in dieser Identifikation keinen positiven Gesamtanschluss über beliebig unabhängig verbundene Halbprofile. Beschränkt die wirkliche Konfigurationsbasis die gesamten Geschichten auf deckgerade Felder oder hat ihre OS Zeit eine andere Geometrie, ist diese Anwendung gesondert auszuwerten.

## 5. Dieselbe gemischte Antwort im chiralen Majoranaoperator

Die in der vorherigen Rechnung hergeleitete Produktrealität trägt `p` auf `−p`, vertauscht die Gesamtchiralität und erhält die entsprechende äußere APS Domäne. Auf der positiven Gesamtchiralität ist die schiefe Majoranaform des gegebenen Produktoperators für ein inverses Fourierpaar blockweise `(0,A;−Aᵀ,0)`. Ihre Pfaffianlinie hat bei festem Rahmen den Betrag der Determinante von `A`. Das `A` Bein wird bis auf feste Isometrien durch den chiralen Operator `N(∂+B_p)` gegeben. Auf dieser Chiralität ist `B_p` die direkte Summe der beiden bereits explizit bewiesenen R032 Sektoren. Die Rechnung benötigt für den Betrag weder eine freie Pfaffianphase noch eine Identifikation dieses Produktes mit dem physischen `P_adm`.

Am konstanten Hintergrund mit äußeren APS Präparationen ist das Greenbein von `L=∂+B_0` exakt `e^(−(s−t)B_0)P_+` für `s>t`, beziehungsweise `−e^(−(s−t)B_0)P_−` für `s<t`. Beide erfüllen die tatsächlichen einseitigen APS Bedingungen und der Sprung ist die Identität. Für disjunkte Halbvariationen `V_u,V_v` ist die gemischte Ableitung der relativ gebundenen chiralen Determinante `−Re Tr(L⁻¹V_u L⁻¹V_v)`. Pro Tangentialmodus sind beide `L⁻¹V` Hilbert Schmidt, und diese Kreuzspur folgt unmittelbar aus der zweiten Ableitung von `det₂(I+L⁻¹δB)` am Ursprung. Lokal additive Normierungen und die lineare Spurkompensation ändern sie nicht. Ein weiterer möglicher nichtlokaler Faktor der Originalquelle müsste dagegen explizit mitgerechnet werden.

In den beiden Greenbeinen setzen sich die Exponenten zu `exp[−(E_+−E_−)(s−t)]` zusammen, mit positiven Quadraten `|Y_tan,+−|²`. Der gemischte Beitrag der chiralen relativen Majoranaform ist damit derselbe diskrete Typ wie der Cap Beitrag, bis auf die tatsächlich verwendete Paarmultiplikation. Am Originalviertel ergeben sich mit beiden Cliffordsektoren `275.2241709559` bei Abstand `0.3` und `195.5747078024` bei `0.9`. Die unabhängige Differenzrechnung mit den zwei vollständigen Einfügungen in der Cap Matrix stimmt bis `2.7e−7` überein. Diese Ableitung bestimmt die getrennte chirale Kreuzantwort am tatsächlich geprüften Produkt. Ein Originalterm über denselben Schnitt muss weiterhin die im ganzen Profil nachgewiesene negative kontinuierliche Dichte ausgleichen oder die tatsächliche Konfigurationsdomäne anders bestimmen.

## 6. Was im ursprünglichen ersten Block zu berechnen bleibt

In der feldunabhängigen Poissonrechnung muss jede separat additive gemischte Logantwort der restlichen Quelle auf einem offenen Frequenzband eine positive kontinuierliche Komponente in mindestens der Größe des dortigen Fehlbetrags liefern, sofern sie selbst als translationsinvariante Spektralform darstellbar ist. Ein ursprünglicher gekoppelter Aufruf kann den verbleibenden Term auch ohne solche Faktorisierung liefern. Getrennte Halbgewichte, feste Referenznenner, eine Eta Phase und zusätzliche reine diskrete Cap Raten können die gezeigte Verletzung in dieser Kontraktion nicht reparieren. Ebenso ändert der in QFT007 enthaltene diagonale IR Regulator `r I` diesen gemischten Jet nicht. Die primitive Definition unterscheidet den Spektralsamen ausdrücklich vom späteren Ghost, Boson und Selektoraufbau. Hier liegt der erste konkrete Test dieses späteren Aufbaus. Er erfordert die tatsächliche physische Halbkonfiguration, ihre Randbasis und die unnormierte gemeinsame Integration.

Originalquellen sind [der primitive Spektralsame](/Users/stefanhamann/Projekte/tfpt-theoryv4/documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-section-primitive-core-and-carrier-interface-handoff.tex:302), [die beiden bewegten Referenzen](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/ORIGINAL-HIGGS-BOUND-REFERENCE-PROOF.txt:48), [die beiden Argumentkarten des Profils](/Users/stefanhamann/Projekte/tfpt-theoryv4/docs/research/fundamental-plan/evidence/CURRENT-ACTION-ARGUMENT-PROOF.txt:37) und [der gemeinsame regulierte Originalaufruf](/Users/stefanhamann/Projekte/tfpt-theoryv4/documents/archive--tfpt-45--source_extracts--04_qft_source/sections/007-subsection-renormalized-physical-and-scheme-readout-laye.tex:3). Der Checker `check_r113_majorana_heat_reflection_20261110.py` benutzt unverändert die v252 Diagnosematrix. Er unterscheidet exakte rationale Aussagen, Intervallgrenzen und unabhängige Quadraturen. Der vollständige physische Quellenwert und seine Selektor und Maßidentität sind unberechnet.
