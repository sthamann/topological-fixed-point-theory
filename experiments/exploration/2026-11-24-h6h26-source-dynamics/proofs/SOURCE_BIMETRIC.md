# H6h26: physikalische Quellen, Transferband und eine konkrete Vierpunktprüfung

**Ausgangspunkt:** H6h26 auf `main 0172a74bd`. Geprüft wurden `experiments/exploration/2026-11-23-phasenraum-f5-holocode/z48_zuordnung_quelle.py`, insbesondere Teile 3 und 4, außerdem die zugehörigen Definitionen in `research/originlab/process_rules.py` und die Aussagen in `docs/research/TFPT-STAND.md`. Dieser Bericht ist eigenständig; er setzt keine historischen Berichte voraus.

**Neue Ergebnisse:** Eine vollständige Bimetrik Quellenantwort, ein strenger Momententest zwischen dem Transferband und einem einzelnen massiven Pol, eine für alle euklidischen Impulse gültige Fehlerschranke der Einpolapproximation und ein expliziter physikalischer Vierpunktvertex, der die verbleibende Potentialfreiheit misst.

Der Wert `q=exp(−2π)` wird in diesen Rechnungen als festgelegt behandelt. Ob die geometrische Nome gerade dem Amplitudenparameter des physikalischen Transfers entspricht, ist eine vorgelagerte Identifikation; diese Rechnung ersetzt sie nicht.

## 1. Universelle gravitative Quellen benötigen mehr als einen invarianten Modus

Eine lineare Fouriermode `J_n` transformiert unter einer vollen Kreisrotation mit einem Faktor `exp(inθ)`. Die lineare Gruppenmittelung erhält nur `J_0`. Daraus folgen jedoch zwei wichtige Einschränkungen:

* Unter der diskreten Vierteldrehung allein bleiben sämtliche Moden `n≡0 mod4` invariant.
* Selbst unter voller Kreisrotation sind zusammengesetzte Dichten wie `J_(−n)J_n` invariant. Rotationsinvarianz allein entfernt sie nicht.

Die Aussage über einen linearen invarianten Strommodus ist daher kein allgemeiner Satz über alle zulässigen Quellenoperatoren. Zusätzlich legt sie nicht fest, dass der betreffende Operator der physikalische Energie und Impulstensor ist.

### 1.1 Was der tatsächliche Z48 Code prüft

Teil 3 erzeugt einige Dichten `1+A cos(kθ)` und prüft numerisch, dass ihr voller Kreismittelwert eins ist. Dieser Test bestätigt eine Fourieridentität. Er prüft weder einen Noethererhaltungssatz aus einer Wirkung noch einen gravitativen Tensorvertex noch eine Wardidentität. Die anschließende Aussage, Ward verlange die eindeutige Quelle `|0⟩`, wird im Code als positiver Befund gesetzt.

Für den tatsächlich verwendeten Transfer ist außerdem die Symmetrie kleiner als die des numerischen Mittelungsverfahrens. Mit `R_θ|n⟩=exp(inθ)|n⟩` gilt exakt

\[
R_\theta A_q R_\theta^\dagger
=e^{-i\theta}A_{q e^{4i\theta}},\qquad
A_q=S^\dagger-qS^3.
\]

Bei festem reellem `q>0` ist `L_q=A_q^†A_q` daher unter `μ₄` invariant, nicht unter jeder Kreisrotation. Eine vollständige kontinuierliche Symmetrie würde eine zusätzliche Transformation des Transferparameters oder eine andere Eichkonstruktion benötigen. Unter `μ₄` ist schon

\[
v_h=|0\rangle+h|4\rangle
\]

für jedes `h` invariant. Die behauptete eindeutige Quellenauswahl folgt also nicht aus der vorhandenen Vierteldrehung.

### 1.2 Der stärkere physikalische Test

Falls `J_0` eine gewöhnliche innere Ladung ist, darf sie nicht ohne weitere Konstruktion mit gravitativer Masse identifiziert werden. Photonen, neutrale Körper und Bindungsenergie sind besonders deutliche Kontrollen. Gleiche Ladung für zwei Zustände mit verschiedenen Energien würde für sich genommen keine Gleichheit von gravitativer und träger Masse liefern.

Die konstruktive nächste Prüfung ist eine aus demselben Prozess abgeleitete metrische Quelle

\[
T^{\mu\nu}_A=\frac2{\sqrt{|g|}}\frac{\delta S_A}{\delta g_{\mu\nu}}
\]

für jeden Materiesektor und ein entsprechender Vertex `Γ_A^(μν)(p',p)`. Beispielsweise lautet die skalare gravitative Wardidentität bei `k=p'−p`, bis auf die gemeinsame Kopplungsnorm,

\[
k_\mu\Gamma_A^{\mu\nu}(p',p)
=p^\nu D_A^{-1}(p')-p'^\nu D_A^{-1}(p).
\]

Der universelle Softgravitonfaktor enthält die Viererimpulse, nicht nur innere Ladungen. Weinbergs ursprünglicher Satz verbindet gerade Lorentzinvarianz und die Entkopplung unphysikalischer Gravitonpolarisationen mit der Universalität gravitativer Kopplungen. Quelle: <https://doi.org/10.1103/PhysRev.135.B1049>.

Für den neuen TFPT Claim sind deshalb drei Dinge zu zeigen: der Operator wird tatsächlich auf `T_(μν)` abgebildet, seine Normierung folgt der Energie einschließlich Bindungsenergie, und die Kopplung erfüllt dieselbe Wardidentität in allen miteinander wechselwirkenden Materiesektoren. Ein gemeinsamer Nullmodus liefert dafür einen möglichen Kanal; die Quellenabbildung muss ihn physikalisch ausfüllen.

## 2. Vollständige lineare Quellenantwort der Bimetrik

Seien

\[
\alpha=\frac{M_f}{M_g},\qquad M_{\rm eff}^2=M_g^2+M_f^2,
\]

und verwende `g=η+2h_g/M_g`, `f=η+2h_f/M_f`. Die normierten kanonischen Eigenvektoren um den gemeinsamen flachen Hintergrund sind

\[
u_0=\frac1{M_{\rm eff}}\binom{M_g}{M_f},\qquad
u_m=\frac1{M_{\rm eff}}\binom{-M_f}{M_g}.
\]

Minimale Quellen koppeln als `h_gT_g/M_g+h_fT_f/M_f`. Im Eigenzustandsraum ergibt dies

\[
\frac{H_0}{M_{\rm eff}}(T_g+T_f)
+\frac{H_m}{M_{\rm eff}}\left(-\alpha T_g+\alpha^{-1}T_f\right).
\tag{1}
\]

Für konservierte Quellen definiere

\[
\mathcal P_0(T,T')=T_{\mu\nu}T'^{\mu\nu}-\frac12TT',\qquad
\mathcal P_m(T,T')=T_{\mu\nu}T'^{\mu\nu}-\frac13TT'.
\]

Die euklidische Austauschantwort lautet dann, mit gemeinsamen konventionellen Amplitudenfaktoren unterdrückt,

\[
\mathcal A(p)=\frac1{M_{\rm eff}^2}
\left[
\frac{\mathcal P_0(T_g+T_f,T'_g+T'_f)}{p^2}
+\frac{\mathcal P_m(-\alpha T_g+\alpha^{-1}T_f,
                    -\alpha T'_g+\alpha^{-1}T'_f)}{p^2+m^2}
\right].
\tag{2}
\]

Im TT Anteil sind die beiden Spektralresiduenmatrizen

\[
\frac1{M_{\rm eff}^2}
\begin{pmatrix}1&1\\1&1\end{pmatrix},\qquad
\frac1{M_{\rm eff}^2}
\begin{pmatrix}\alpha^2&-1\\-1&\alpha^{-2}\end{pmatrix}.
\tag{3}
\]

Beide sind positive Matrizen vom Rang eins. Ein negatives Kreuzresiduum im zweiten Ausdruck ist damit vereinbar; es ist kein negativer Eigenwert des Spektralmaßes.

Für ausschließlich an `g` gekoppelte sichtbare Materie ist der normierte TT Propagator

\[
C_g(p^2)=\frac1{p^2}+\frac{\alpha^2}{p^2+m^2}.
\tag{4}
\]

Diese minimale gemeinsame metrische Kopplung ist in der klassischen Testkörpergrenze universell. Ein massiver Yukawaanteil verletzt die schwache Äquivalenz nicht automatisch. Alle Testkörper erfahren dann dieselbe zusätzliche Beschleunigung pro träger Masse.

Die Quellenkonvention und die Tatsache, dass eine gemeinsame Materiemetrik beide Eigenzustände sieht, stimmen mit der primären Bimetrikliteratur überein: <https://arxiv.org/abs/1712.10020>, Abschnitt 4; <https://arxiv.org/abs/2101.08795>, Abschnitt 2. Eine direkte Kopplung desselben Materiesektors an beide Metriken ist eine andere nichtlineare Theorie mit eigener Constraintprüfung.

## 3. MICROSCOPE: massloser Quellenfehler und massive Zusatzkraft unterscheiden

Das Endresultat lautet

\[
\eta(\mathrm{Ti},\mathrm{Pt})=
[-1.5\pm2.3\ \mathrm{stat}\pm1.5\ \mathrm{syst}]\times10^{-15}.
\]

Die statistische Unsicherheit wird mit einer Standardabweichung angegeben. Gemessen wurde der differentielle freie Fall von Titan und Platinlegierungen im Feld der Erde. Primärquelle: <https://arxiv.org/abs/2209.15487>.

### 3.1 Die tatsächliche Z48 Annahme betrifft den masslosen Pol

Z48 verwendet ausdrücklich den Überlapp

\[
\kappa_A=\sqrt Z\left(\rho_{0,A}+q\rho_{4,A}+\cdots\right),
\qquad Z=1-q^2.
\]

Die zusätzliche Harmonische verändert damit im angegebenen Gedankenmodell den **masslosen** Pol. Sie wird an dieser Stelle nicht als bloßer massiver Yukawaanteil behandelt.

Setzt man zusätzlich `ρ_(0,A)/M_A` für alle Testkörper gleich, beschränkt sich auf die erste Zusatzharmonische und schreibt `h_A=ρ_(4,A)/ρ_(0,A)`, so ist die Beschleunigung pro träger Masse proportional zu `1+q h_A`. Der gemeinsame Erdfaktor und `√Z` kürzen sich aus dem Eötvösparameter heraus. Exakt folgt

\[
\boxed{\eta_{AB}^{(0)}=
\frac{q(h_A-h_B)}{1+q(h_A+h_B)/2}.}
\tag{5}
\]

Die Näherung `η≈qΔh` ist unter genau diesen Quellen und Normierungsannahmen korrekt. **Da dieser Ansatz den masslosen Pol verändert, ist seine Grenze unabhängig von der massiven Reichweite `1/μ`.** Eine Verlegung der massiven Linie an die Planckskala würde diesen angenommenen masslosen Quellenfehler nicht abschirmen.

Die grobe Hilfsgrenze `|η|<10^(−14)` ergibt bei kleiner mittlerer Harmonischer `|Δh|≲10^(−14)/q≈5.35×10^(−12)`. Diese Größenordnung ist eine bedingte Übersetzung der Messung, kein separat bestimmtes Konfidenzintervall und keine Herleitung der zugrunde gelegten Quellenabbildung.

Aus (5) folgt auch nicht, dass die Quelle exakt `|0⟩` sein muss. Derselbe nichtverschwindende Wert `h_A=h_B` ergibt universellen freien Fall. Hinreichend kleine Unterschiede bleiben experimentell erlaubt. Selbst eine theoretisch geforderte strikte Universalität legt zunächst das Verhältnis von effektiver Ladung und träger Energie fest; sie bestimmt nicht eindeutig den mikroskopischen Quellenvektor.

Nichtuniverselle Kopplungen eines reinen masslosen Spin 2 Feldes müssen außerdem erst mit seinen gravitativen Wardidentitäten verträglich gemacht werden. (5) ist eine phänomenologische Diagnose, keine bereits konsistente vollständige Theorie solcher Kopplungen.

### 3.2 Ein massiver Effekt benötigt hingegen eine Reichweitenrechnung

Wenn eine zusätzliche Quelle stattdessen an die massive Spin 2 Mode koppelt, ist eine andere Rechnung notwendig. Sei `χ_A` ihre Ladung pro träger Energie und die universelle masselose Ladung eins. Aus (2) folgt im linearen, ungescreenten, nichtrelativistischen Bereich für punktförmige Körper

\[
V_{EA}(R)=-\frac{G_N M_E M_A}{R}
\left[1+\frac43\chi_E\chi_A e^{-mR}\right].
\tag{6}
\]

Der Faktor `4/3` folgt aus dem Unterschied der Spurterme in (2). In Einheiten `c=ℏ=1` ist `1/m` die Reichweite. Für zwei Testkörper gilt

\[
\eta_{AB}^{(m)}(R)=
\frac{\frac43\chi_E(\chi_A-\chi_B)(1+mR)e^{-mR}}
{1+\frac23\chi_E(\chi_A+\chi_B)(1+mR)e^{-mR}}.
\tag{7}
\]

Hier sind die Erdkopplung, die physikalische Masse, die Quellengeometrie und gegebenenfalls Abschirmung wesentlich. Für eine homogene Erdkugel mit Radius `R_E` multipliziert außerhalb der Kugel zusätzlich

\[
F_E(z)=\frac{3(z\cosh z-\sinh z)}{z^3},\qquad z=mR_E
\]

den Yukawaanteil. Eine reale Verteilung verlangt die entsprechende räumliche Integration. Die Anwendung einer linearen Bimetrikformel verlangt außerdem einen Test des Vainshtein Regimes; eine Primäranalyse dazu ist <https://arxiv.org/abs/2101.08795>.

Die gemeinsame minimale Kopplung aller sichtbaren Materie an `g` setzt `χ_A=χ_B=−α`. Dann ist (7) exakt null, obwohl eine universelle zusätzliche Kraft bestehen kann. Falls `χ_A−χ_B=qΔh` angesetzt wird und die Erdkopplung selbst von Ordnung `q` ist, wäre dieser massive Effekt von Ordnung `q²Δh`, zusätzlich zur Reichweitenfunktion. Er ist vom tatsächlichen masslosen Z48 Gedankenmodell (5) verschieden.

MICROSCOPE besitzt eine separate Analyse kurzreichweitiger Yukawakräfte innerhalb des Sensors. Sie nutzt andere Messungen und liefert nicht einfach die Erdfallgrenze von `10^(−15)` für beliebige Massen: <https://arxiv.org/abs/2102.00022>.

**Konsequenz für H6h26:** Die Näherung in Z48 darf unter ihrer masslosen Quellenannahme bestehen bleiben. Die Schlussfolgerungen „Ward ist hergeleitet“ und „nur `|0⟩` ist möglich“ folgen daraus nicht. Für eine tatsächliche massive Bimetriksignatur ist stattdessen (7) mit der richtigen Geometrie und dem passenden nichtlinearen Regime anzuwenden.

## 4. Bimetrikgewicht und vollständiges Transferband sind verschiedene Daten

Setze

\[
q=e^{-2\pi},\qquad r=q^2=e^{-4\pi},\qquad Z=1-r.
\]

Der normierte Hardy Propagator hat massives Gesamtgewicht

\[
w=\frac rZ=\frac1{e^{4\pi}-1}
\simeq3.4873545\times10^{-6}.
\tag{8}
\]

Die Gleichsetzung `α²=w` in (4) ist deshalb ein mathematisch konsistentes Matching des gesamten massiven Gewichts. Sie ist zugleich ein Matching des Verhältnisses zwischen UV und IR Residuum: Für den Hardy Kern ist dieses Verhältnis `1/Z`, für die Bimetrik `1+α²`.

Sie ist jedoch keine Gleichheit der vollen Spektralantwort. Der Hardy Kern hat ein kontinuierliches Maß auf

\[
[\mu^2(1-q)^2,\ \mu^2(1+q)^2],
\]

während (4) genau ein massives Atom hat. Die Eindeutigkeit der Stieltjestransformation verbietet für `q>0` ihre Identität. Der volle Banddurchmesser beträgt relativ zu `μ²` etwa `0.00746977`; das gesamte Bandgewicht ist zugleich sehr klein.

### 4.1 Der erste nichttriviale Momententest

Normiere beide Antworten auf masseloses Residuum eins. Schreibe `H=x−ax²+bx³−cx⁴+…`. Für das Transferband gilt

\[
a_{\rm MP}=\frac r{\mu^2Z^2},\quad
b_{\rm MP}=\frac{r(1+r)}{\mu^4Z^4},\quad
c_{\rm MP}=\frac{r(1+3r+r^2)}{\mu^6Z^6}.
\tag{9}
\]

Ein einzelner massiver Pol mit Gesamtgewicht `w` hat

\[
a_{\rm bi}=\frac w{m^2},\quad
b_{\rm bi}=\frac{w(1+w)}{m^4},\quad
c_{\rm bi}=\frac{w(1+w)^2}{m^6}.
\]

Fordert man außer `w=r/Z` auch `a_bi=a_MP`, ist seine Masse zwingend

\[
m^2=\mu^2Z.
\tag{10}
\]

Danach bleiben dennoch die exakten Unterschiede

\[
\boxed{b_{\rm MP}-b_{\rm bi}=a^2},\qquad
c_{\rm MP}-c_{\rm bi}=\frac{r^2(3+r)}{\mu^6Z^6}.
\tag{11}
\]

Die relativen Unterschiede sind bei der kleinen Nome klein, aber nicht null. Insbesondere ist der relative Unterschied von `b` gleich `r/(1+r)≈3.48733×10^(−6)`. Keine Wahl der einzelnen Polmasse kann gleichzeitig `w,a,b` treffen.

### 4.2 Ein invariant formulierter Spektraltest

Für eine positive massive Maßkomponente `ρ` seien

\[
w=\int\rho(ds),\qquad a=\int\frac{\rho(ds)}s,
\qquad m_{-2}=\int\frac{\rho(ds)}{s^2}=b-a^2.
\]

Die Cauchy und Schwarz Ungleichung liefert

\[
\Delta_{\rm KL}=w(b-a^2)-a^2\ge0.
\tag{12}
\]

Gleichheit gilt genau für ein einzelnes massives Atom. Das MP Band besitzt dagegen

\[
\Delta_{\rm KL}^{\rm MP}=\frac{r^3}{\mu^4Z^5}>0.
\tag{13}
\]

Dieser Test erfasst den Unterschied unabhängig von der Benennung der Hilfskanäle. Die Stromnorm liefert höchstens `w`; für eine exakte Bimetrikableitung müsste der Ursprung zusätzlich die Sättigung (12) ergeben. Das ursprüngliche nichtdegenerierte Transferband erfüllt sie gerade nicht.

### 4.3 Kontrollierte Einpolapproximation mit festem Fehlerzeichen

Es lässt sich sogar ein Vergleich für alle euklidischen `x>0` beweisen. Setze `u=1/s` und verwende das Wahrscheinlichkeitsmaß `ρ/w`. Dann ist der massive Antwortanteil

\[
C_{\rm massiv}(x)=w\,\mathbb E\frac u{1+xu}.
\]

Die Funktion `φ(u)=u/(1+xu)` ist strikt konkav. Bei gleichem `w` und `a=w E[u]` folgt deshalb

\[
C_{\rm MP,massiv}(x)
<\frac w{x+m^2}=C_{\rm bi,massiv}(x),\qquad x>0,
\tag{14}
\]

mit der Matchmasse (10). Der exakte Jensenrest ist

\[
C_{\rm bi,massiv}-C_{\rm MP,massiv}
=\frac{wx}{(1+x\bar u)^2}
\mathbb E\left[\frac{(u-\bar u)^2}{1+xu}\right]>0.
\tag{15}
\]

Für kleine `x` beginnt er mit `a²x`. Somit ist die Bimetrik eine mögliche kontrollierte Einpolapproximation. Ihre exakte Identifikation mit dem unveränderten Hardy Band wäre dagegen falsch. Alternativ müsste ein nichtlokaler Selbstenergiebeitrag das integrierte Band erhalten; dann handelt es sich nicht mehr um die bloße lokale Zwei Metriken Theorie.

## 5. Ein neuer konkreter Vierpunkttest für die Potentialfreiheit

Die bisherige Potentialfamilie lautet

\[
V_\tau(X)=\frac{-4+e_1(X)+e_3(X)-4e_4(X)}2
+\frac\tau2\det(I-X).
\tag{16}
\]

Sie hat dieselben Zweier und Dreiervertizes für jedes `τ`. Damit ist die Differenz ihrer Vierpunktamplituden rein durch den neuen Kontaktvertex bestimmt: sämtliche Austauschdiagramme bleiben gleich.

### 5.1 Kanonischer relativer Tensor

Definiere

\[
M_{\rm rel}^{-2}=M_g^{-2}+M_f^{-2},\qquad
m^2=\frac{k_{\rm edge}}{M_{\rm rel}^2}.
\]

Für das kanonische massive Feld `H` gilt um den flachen Hintergrund

\[
X-I=\frac1{M_{\rm rel}}\eta^{-1}H+O(H^2,H_0H).
\]

Bei Signatur `(-+++)` folgt aus dem Potentialbeitrag `−k_edge V_tau` der quartische Unterschied

\[
\Delta\mathcal L_4
=\frac{\tau m^2}{2M_{\rm rel}^2}\det(H_{\mu\nu}).
\tag{17}
\]

Das Vorzeichen berücksichtigt `det(η^(−1))=−1`. Alle folgenden Amplituden verwenden diese explizite Konvention.

### 5.2 Physikalische massive Helizität 0 statt beliebiger Testmatrizen

Wähle die Vorwärtskinematik

\[
p_\pm=(E,0,0,\pm k),\qquad E^2=k^2+m^2,
\qquad s=4E^2.
\]

Mit den transversalen Vektoren `e_x,e_y` und den kovarianten longitudinalen Vektoren

\[
\ell_\pm=(-k/m,0,0,\pm E/m)
\]

seien die normierten Helizität 0 Tensoren

\[
\epsilon_\pm=\frac1{\sqrt6}
\left(e_xe_x^T+e_ye_y^T-2\ell_\pm\ell_\pm^T\right).
\]

Sie sind transversal, spurfrei und haben Norm eins. Direktes Ausrechnen ergibt

\[
\det(A\epsilon_++B\epsilon_-)
=\frac{4E^2k^2}{9m^4}AB(A+B)^2.
\tag{18}
\]

Für vier äußere Beine setze `A=z₁+z₃`, `B=z₂+z₄`. Der Koeffizient von `z₁z₂z₃z₄` ist `32E²k²/(9m⁴)`. Damit ist die Änderung der vollständigen Baumamplitude exakt

\[
\boxed{\Delta\mathcal M_{0000}(s,0)
=\frac{\tau\,s(s-4m^2)}{9M_{\rm rel}^2m^2}.}
\tag{19}
\]

Insbesondere

\[
\frac12\partial_s^2\Delta\mathcal M_{0000}
=\frac\tau{9M_{\rm rel}^2m^2}.
\tag{20}
\]

Hier verschwindet der quartische Effekt nicht wie bei vier speziellen rein räumlichen Tensoren in derselben Ruhestellung. Er ist eine konkrete messbare Potentialabhängigkeit einer physikalischen Helizitätsamplitude. Die Rechnung bestimmt die Amplitudendifferenz, nicht die noch fehlende vollständige Baselineamplitude.

## 6. Was unabhängige Konsistenzbedingungen schon auswählen

Die Primäranalyse von Bonifacio, Hinterbichler, Joyce und Rosen, <https://arxiv.org/abs/1712.10020>, verwendet eikonale Zeitverzögerungen für eine massive und eine masselose Spin 2 Mode. Unter ihren EFT und Skalentrennungsannahmen wird die kubische Potentialkombination `c₃=1/4` gefordert. Genau dies entspricht der palindromischen Familie `β_n=β_(4−n)`. Der quartische Parameter bleibt frei.

In deren Konventionen entspricht (16)

\[
c_3=\frac14,\qquad d_5=\frac{\tau-4}{96}.
\tag{21}
\]

Die vollständige Familie besitzt daher einen unabhängigen positiven kubischen Konsistenzhinweis. Dieser Test entscheidet `τ` ausdrücklich nicht.

Für eine skalare Materiequelle mit unserer massiven Ladung `χ` und universeller masseloser Kopplung übersetzen sich die beiden führenden eikonalen Bedingungen zu

\[
-\alpha\le\chi\le\alpha^{-1}.
\tag{22}
\]

Die minimale `g` Quelle liegt am unteren, die minimale `f` Quelle am oberen Rand. Dies ist ein zusätzlicher konkreter Test für einen behaupteten Quellenprozess. Er ist keine hinreichende Bedingung für eine vollständige UV Theorie.

### 6.1 Keine naive Vorwärtspositivität bei einem masselosen Graviton

Gewöhnliche Vorwärtspositivität würde in einer geeigneten Theorie mit Masselücke einen polsubtrahierten zweiten Ableitungskoeffizienten begrenzen. In der Bimetrik bleibt aber ein masseloser gravitativer `t` Kanalpol, typischerweise von der Form `s²/t`. Deshalb darf man nicht einfach die bekannte massive Positivitätsformel einsetzen oder nach einem beliebigen Polabzug `f>0` behaupten.

Die ausdrückliche Analyse dieses Problems ist <https://arxiv.org/abs/2007.12667>. Gravitationsabhängige Positivitätsschranken mit zusätzlichen Annahmen über Reggeverhalten werden beispielsweise in <https://arxiv.org/abs/2007.15009> hergeleitet. Die Arbeit <https://arxiv.org/abs/1910.11799> liefert nützliche Methoden für wechselwirkende massive Spin 2 EFTs, ersetzt aber nicht die Behandlung des masselosen Pols der hier vorliegenden Theorie.

Der Unterschied (19) ist dennoch wohldefiniert: Die identischen Austauschpole heben sich zwischen zwei Mitgliedern derselben `τ` Familie auf. Daraus folgt keine Positivität dieser Differenz. Es liefert vielmehr eine exakt bekannte lineare Abhängigkeit, die in einen vollständig spezifizierten gravitativen Dispersions oder Matchingtest eingesetzt werden kann.

## 7. Der nächste stärkste konkrete Test

**Erstens: die physikalische Zweipunktantwort identifizieren.** Aus derselben mikroskopischen Quellenregel ist der komplette konservierte Tensorquellenprozess herzuleiten, einschließlich Wardidentität und Massennormierung. Wenn er das vollständige MP Band besitzt, kann die lokale Bimetrik nur eine Approximation sein. Die zwei Werte `w` und `a` reichen nicht; (12) oder schon `b` entscheidet diese Frage exakt.

**Zweitens: den experimentellen Quellenprozess berechnen.** Zunächst ist zu entscheiden, ob eine Artenabhängigkeit das masslose Residuum wie in (5) verändert oder einen massiven Zusatzkanal wie in (7) betrifft. Im ersten Fall ist die Energie und Quellenidentifikation entscheidend; im zweiten zusätzlich die Masse, Erdkopplung, Geometrie und Abschirmung. Gemeinsame minimale Materiekopplung erfüllt WEP in der klassischen Testkörpergrenze, auch wenn universelle Yukawakorrekturen bestehen.

**Drittens: den quartischen Parameter an einer physikalischen Amplitude testen.** Aus dem gemeinsamen Ursprung ist ein Vierquellenprozess in der kanonischen massiven Helizität 0 zu bestimmen. Nach expliziter Identifikation derselben Zweier und Dreiervertizes ist seine Differenz zur Baseline frei von den gemeinsamen Austauschpolen. Der Koeffizient in (20) würde dann unmittelbar `τ` bestimmen. Ohne einen solchen Ursprungseingang können eikonale Bedingungen die kubische Familie eingrenzen; eine weitere gravitative Positivitätsanalyse kann Parameterbereiche testen, bestimmt aber nicht automatisch einen einzigen Wert.

Diese Testfolge ist stärker als eine weitere reine Zählung: Sie entscheidet die Operatoridentität, liefert eine wirkliche Vergleichsgröße für MICROSCOPE und macht die verbliebene nichtlineare Freiheit durch eine explizite Amplitude prüfbar.

## 8. Reproduzierbarkeit

Die neuen Module `research/originlab/process_source_bimetric.py` und `research/originlab/test_process_source_bimetric.py` gehören zu diesem Nachweis. Zwölf Testgruppen bestehen. Das direkt JSON serialisierbare `source_bimetric_audit()` enthält zwanzig exakte Auditidentitäten. Dazu gehören die vollständigen Quellenresiduenmatrizen, der Unterschied zwischen tatsächlicher Vierteldrehung und zusätzlich unterstellter Kreisrotation, universelle nichtverschwindende Harmonische als Gegenkontrolle, die zwingende Einpolmatchmasse, die abweichenden Spektralmomente, der positive Jensenrest, die physikalischen Helizitätstensoren und ihre Normierung, die Kombinatorik der vier äußeren Beine, die Amplitude (19) sowie die Parameterabbildung (21).

Öffentliche Rechenfunktionen sind `bimetric_source_data`, `massless_harmonic_eotvos`, `massive_yukawa_eotvos`, `transfer_band_vs_bimetric`, `quartic_forward_data` und `eikonal_family_data`. Binäre Fließkommazahlen werden nicht als exakte Eingaben akzeptiert. Die numerischen Nomewerte dienen der Größenordnung, nicht der Entscheidung der analytischen Aussagen.

Es wurde weder eine vollständige Basisamplitude noch ein kompletter gravitativer Loop noch eine neue experimentelle Datenanalyse vorgetäuscht. Die hier gewonnenen Beziehungen sind exakte algebraische Konsequenzen der jeweils ausdrücklich genannten Quellen und Potentialkonventionen.
