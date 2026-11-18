# Ganze neutrale Quelle des bestehenden R048 Radialintegrals

5. November2026. Liveauftrag Plan142, Atlas d1ff953daecc63e28671b479. Keine kanonischen oder öffentlichen Änderungen.

Der bereits vorhandene R048 Radialcaller ist jetzt für die ganze neutrale Quelle u und die komplexe skalare Hilfsquelle geschlossen ausgewertet. Sein Determinantenquotient, die absolute ursprüngliche Baseline und sämtliche neutralen Momente bleiben erhalten. Die erste Gleichsetzung dieses integrierten Subcallers mit der ganzen ursprünglichen Higgsquelle scheitert dagegen an den ausdrücklich vorhandenen ursprünglichen kubischen und quartischen Vertices. Dafür wird unten der erste nichtverschwindende ursprüngliche Schwinger Dyson Rest ausgerechnet.

## Genau der vorhandene Caller

COMMON-CURVATURE-RESPONSE-PROOF.txt:117–174 gibt im benannten flachen gebrochenen stationären Zweig μ0²>0 die kanonische Radialkarte

Hc=(0,v+η/√2), v²=μ0²/(2λ), m_rad²=4λv², ξ=1/6,

K=K_R=−Δ+m_rad²+ξR, l=J−gR, g=√2ξv.

Alle Yukawas und MR bleiben in den ursprünglichen Koeffizienten λ,μ0². Metrik, übrige Felder und der tatsächliche ursprüngliche Regulator sind die konditionalen Argumente dieser gegebenen Subintegration. K ist der ganze originale regulierte Radialoperator, keine Einmodennäherung. Die ursprüngliche Existenzbedingung dieses quadratischen Integrals ist K>0. Wo der Regulator einen zusätzlichen originalen IR Beitrag enthält, gehört er unverändert zu K.

## Tatsächlich integrierte ganze neutrale Quelle

Setze u=|Hc|²=v²+√2vη+η²/2. Eine neutrale Testquelle q wird mit dem ursprünglichen Volumenpaar ∫dvol_g q u gekoppelt. Sie ist der Quellparameter des bereits vorhandenen invarianten Feldpolynoms, keine neue physische Kopplung. Regulatorweise kann sie als zusammengefasste ursprüngliche lineare H Quellenableitung definiert werden. Ihr Support und alle übrigen Quellrestriktionen werden beibehalten. Für reelle genügend kleine q gelte K_q=K−M_q>0; M_q ist der tatsächliche Multiplikationsoperator mit q. Dies reicht für die analytische Auswertung nahe q=0. Die Quelljets selbst bleiben dieselben formalen Koeffizienten der ursprünglichen Subquelle.

Die exakte quadratische Quellenidentität ist

−½(η,Kη)+(l,η)+(q,u)
=−½(η,K_qη)+(l_q,η)+v²(q,1),

l_q=l+√2v q.

Mit G=K⁻¹, G_q=K_q⁻¹ und der tatsächlichen ursprünglichen nullsourcigen q Baseline Z_quad[R,J] folgt durch vollständige Gaußintegration

\[
Z_u[q;R,J]=Z_{\rm quad}[R,J]\,
e^{v^2(q,1)}\det(I-GM_q)^{-1/2}
\exp\!\left\{\frac12(l_q,G_q l_q)-\frac12(l,G l)\right\}.
\]

Der Ausdruck enthält die ganze Nichtkommutativität von K und M_q. Sie werden nicht gemeinsam diagonalisiert. Auf K_q>0 ist der positive Determinantenzweig durch q=0 gebunden; det(I−GM_q) stimmt mit der selbstadjungierten Darstellung det(I−K^(-1/2)M_qK^(-1/2)) überein. Kein neues Normierungsmaß wird gewählt. Die ursprünglichen 2π und Zellfaktoren, S_base, Referenz und alle gegebenen absoluten Faktoren liegen unverändert in Z_quad[R,J]. Die Formel setzt diese Baseline nicht auf1. Der Quotient bei q=0 ist exakt1.

Damit sind sämtliche neutralen Quellmomente dieses ganzen vorhandenen Subcallers durch Ableitung bestimmt, einschließlich ihrer ursprünglichen räumlichen und zeitlichen Moden. Das ist eine tatsächliche Feldintegration der bestehenden Quadratik; die verbleibenden ursprünglichen anderen Felder sind hier konditionale Argumente.

Für b(x)=(Gl)(x) und a(x)=v+b(x)/√2 ergibt sich insbesondere

\[
\overline u(x)=a(x)^2+\frac12G(x,x),
\quad
C_u(x,y)=2a(x)a(y)G(x,y)+\frac12G(x,y)^2.
\]

G(x,x) ist der durch denselben ursprünglichen Regulator bestimmte Kontakt. Er wird nicht verworfen oder mit einer neuen Renormalisierung gleichgesetzt. Auf R=J=0 gilt b=0 und daher ⟨u⟩=v²+G(x,x)/2. Der klassische Wert v² ist bereits in dieser tatsächlich vorhandenen Subquelle nicht das ganze Quantenmoment.

Die vorhandene Krümmungsantwort ist konsistent enthalten: δ_RK=ξδ_R, δ_Rl=−gδ_R. Einschließlich des S_base Beitrags ξv²R liefert die vollständige Ableitung der oben integrierten Quadratik genau −ξ⟨u(x)⟩ und alle übrigen tatsächlichen S_base Ableitungen. Ein externes R Quellenbein wird dabei nicht zur unabhängigen Integration über die skalare Krümmung umgedeutet; R bleibt die Krümmung der gegebenen konditionalen Metrik.

## Ganze skalare Hilfsquelle anschließen

Für denjenigen ursprünglichen lokalen kanonischen Zweig, dessen Nichtminimalität F(H)=F0+c_u u lautet, ist μ=F0+c_u u+2αR. Im oben gegebenen R048 Zweig ist c_u=ξ. α ist der tatsächliche R² Koeffizient desselben Callers und wird nicht aus einer anders normalisierten M Schreibweise übernommen. Mit kovarianter Hilfsquelle θ und der bereits exakt bewiesenen komplexen Hilfsintegration gilt

\[
Z_{u,\rho}[q,\theta;R,J]
=e^{(\theta,F_0+2\alpha R)-\alpha(\theta,\theta)}
 Z_u[q+c_u\theta;R,J].
\]

Das ist die ganze neutrale gemeinsame Quelle der bestehenden quadratischen H Subintegration und ihres normiert eingesetzten komplexen Hilfsbeins. Auf diesem festen Hintergrund folgen

⟨ρ(x)⟩=F0+c_u ū(x)+2αR(x),

Cov(u(x),ρ(y))=c_u C_u(x,y),

Cov(ρ(x),ρ(y))=c_u²C_u(x,y)−2αδ_g(x,y).

Der letzte Kontakt ist der tatsächliche komplexe Hilfskontakt. Das ganze physische μ Wort besitzt c_u²C_u ohne diesen Hilfskontakt. Bei Zellquellen θ_iρ_i lautet der Kontakt entsprechend −2αδ_ij/V_i. V_i bezeichnet das originale Zellvolumen und wird nicht mit dem Vakuumparameter v verwechselt. Die gebundene Referenz und jede tatsächlich mitgeführte Referenzableitung behalten die ursprüngliche Vorschrift. Dies ist R048s benannter Referenzbranch; eine andere H abhängige R032 Referenz darf nicht nachträglich eingefroren werden.

## Kleinster entscheidender Test gegen den ganzen Originalcaller

Dieselbe Originalstelle nennt ausdrücklich die noch vorhandenen Radialvertices

V_int(η)=√2λvη³+λη⁴/4.

QFT007:119–145 schreibt die ursprüngliche ganze Quelle ausdrücklich als Gaußbasis mit sämtlichen höheren S_adm Ableitungen. DYNAMICAL-HIGGS-SOURCE-PROOF.txt:34–78 erhält zudem die ganze konditionale H Faser, ihre tatsächliche affine Yukawaquelle J−j und gegebenenfalls S_extra(δ_J,ψ). Keines dieser Stücke ist durch K bestimmt oder entfallen.

Schon ohne weitere Felder lautet der erste volle Radialkontakt auf R=J=0

K⟨η(x)⟩+3√2λv⟨η(x)²⟩+λ⟨η(x)³⟩=0.

Der Versuch, die soeben integrierte ursprüngliche Quadratik unverändert als ganze radiale Quelle einzusetzen, gibt hierin den exakt bestimmten Rest

\[
3\sqrt2\lambda v\,G(x,x)\ne0
\]

auf dem angegebenen positiven gebrochenen regulierten Zweig. Eine weitere Einfügung η(y) liefert als zusätzlichen Rest 3λG(x,x)G(x,y). Dies ist der entscheidende Gegenfall innerhalb des bereits vorhandenen Originalsubcallers und seiner Originalvertices, keine neue freie Hilfsfamilie. Zusätzliche ursprüngliche Felder und ihre Ströme müssen im ganzen Caller mitintegriert werden und können diese Antworten ändern. Der Test behauptet deshalb keinen globalen Ausschluss der Theorie; er verbietet die automatische Ersetzung ihrer ganzen Quelle durch die vorhandene Quadratik.

Die erste noch nicht ausgeführte Erweiterung ist die ursprüngliche Gaußintegration mit dem UNGEKÜRZTEN Faktor

exp[−∫(√2λvη³+λη⁴/4)−S_extra(η,ψ)]

innerhalb derselben Integration und mit der tatsächlichen affinen Quelle J−j. Danach bleiben die ursprünglichen übrigen Feldintegrationen auszuführen. Keine Integration wird mit einer unendlichen Reihe vertauscht.

Die häufige Schreibweise

exp[−∫(√2λv δ_J³+λ δ_J⁴/4)−S_extra(δ_J,ψ)] Z_quad[R,J−j]

ist hier NUR eine formale Diagrammoperation mit ursprünglicher Reihenfolge und Paarung. Eine konvergente Exponentiation dieses Differentialoperators und ihre Vertauschung mit der Integration wurden nicht nachgewiesen. Auch bei positiver Quartik folgt aus der Existenz des ursprünglichen Integrals keine Konvergenz seiner Gaußstörungsreihe. Eine analytische Bedeutung dieser Schreibweise darf hier ausschließlich durch das gerade bezeichnete ungekürzte ursprüngliche Reweighting definiert werden. Sie ist keine bereits ausgeführte volle Integration. Die neutrale q Quelle der Quadratik ist oben exakt bekannt; gerade die tatsächliche Auswertung der ursprünglichen höheren Vertices ist noch nötig. Die unabhängige geometrische Integration einschließlich der ganzen räumlichen Stresswörter bleibt ebenfalls erhalten.

λ ist in R048 bereits durch λ=8π²B/(dA²) gebunden und auf dem benannten Zweig positiv: d>0, A=Tr DY²>0 und B=Tr DY⁴>0 bei nichtverschwindender ursprünglicher Yukawafamilie. Es wird weder neu angepasst noch ein quartischer Gegenfaktor ergänzt. Bei degenerierten Y, A=0 oder außerhalb μ0²>0 gilt diese gebrochene Radialkarte nicht. Der von R048 exakt integrierte Caller bleibt ausdrücklich die QUADRATIK bei festem Hintergrund und ursprünglichem endlichem Regulator. Er ist weder die ganze ursprüngliche H Faser noch die vollständig integrierte Γ_ren. Die Typgrenze folgt zudem ausdrücklich aus QFT007:119–145.

Eine bare zusätzliche räumliche Metrikableitung ersetzt kein integriertes Stresswort: δ_g logZ enthält die erwartete ganze Aktions und Maßantwort, ihre gemischte u Ableitung deren volle Kovarianz und echten Kontakte. Der erste Test dieser Identifikation liefert folglich keine unabhängig ausgerechneten Werte für I_orig[R] oder I_orig[uR]. Es wurde keine weitere unintegrierte Stressvertexrunde gerechnet und kein neuer homogener Eliminator aufgebaut.

Die alte O1 Kovarianzklassifikation liefert keinen vergessenen ursprünglichen Wert für I_orig[u]: ORIGINAL-O1-HIGGS-LAW-PROOF.txt:1–20 benennt die Kovarianz cP und c ausdrücklich als zusätzliche, nicht aus dem Index hergeleitete Testannahmen. Die dokumentierte Kosmologierückreaktion v361 erhält ihre vorhandenen Alpha, Λ und Einsteinbezüge, berechnet in ihrem tatsächlichen Aufruf aber weder diesen neutralen H Quellenjet noch I_orig[uR] oder den gemischten UY Wert. Keine dieser vorhandenen Aussagen wird zurückgenommen oder als neue offene globale Quellenkarte deklariert.

## Ausführung und Belege

12 von12 exakten Kontrollen bestehen. Geprüft sind die ganze quadratische Quellenverschiebung, neutrale erste und zweite Quellkoeffizienten, ihre Hilfsquellenverkettung und Kontakte sowie beide expliziten ursprünglichen Schwinger Dyson Reste. Die zwei markierten Quellenplätze im Prüfer extrahieren Koeffizienten des ganzen vorhandenen G; sie sind keine physische Zweimodentrunkierung. Kein ursprünglicher primitiver Zahlenwert und keine physische Promotion wird behauptet.

Originaldateien, Zeilen und vollständige SHA256 stehen in check_original_radial_composites.json. Tragend sind COMMON-CURVATURE-RESPONSE-PROOF.txt:117–174, SHA256 2b05629937b9d46883d4c456c5bb30acecdc59e7876767eb2207d33026683df2; DYNAMICAL-HIGGS-SOURCE-PROOF.txt:34–78,88–116, SHA256 eda5bb664b1fb20370e59b39dbd5e3cb04d79864a8d1d0a0eb6c791cdaf98f49; Source45 QFT007:1–28,119–145, SHA256 9c1c587cbd61fc4034649d778d16f1ee0e8c1e1189885babb3d8148a6fe72d7c. Die bereits geschlossenen Wick, Referenz, BRST, U und Zeitkarten bleiben unverändert verwendbar.
