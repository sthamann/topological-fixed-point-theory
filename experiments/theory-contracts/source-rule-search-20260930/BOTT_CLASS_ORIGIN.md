# TFPT: gemeinsamer Klassenweg aus den Originaldaten

Stand 30.09.2026. Forschungs-ID UR.SOURCE.RULE_SEARCH.01. Verdict PARTIAL.
Die fundamentale Gesamtlösung ist nicht gefunden. Der folgende Klassenweg ist eine explizite mathematische Konstruktion aus benannten TFPT-Daten; seine Identifikation mit dem wirklichen Quellenoperator bleibt offen.

## Ergebnis, das den Suchweg verändert

Unter den ursprünglichen P1-Typen gibt es einen präzisen Weg von der Nahtwindung zur fundamentalen E8-Klasse:

\[
[u_\Sigma]\in K^1(S^1)
\quad\boxtimes\quad
([S^\pm]-[1])\in\widetilde K^0(S^2)
\ \longrightarrow\
\widetilde K^1(S^3)
\ \longrightarrow\
|c_2^{\rm basic}(E_8)|=1.
\]

Dabei ist S² die im P1-Hardening bereits verwendete kompaktifizierte orientierte Normalfläche. Eine runde Metrik oder volle Spin(3)-Symmetrie wird für diese Aussage nicht gebraucht. Die Sphäre hat eine eindeutige Spinstruktur; ihre beiden chiralen Spinorlinien haben erste Chernklassen ±1. Das äußere Produkt mit der Nahtklasse ist daher ein Bott-Generator.

Die tatsächlichen Standardmodell-Ladungen und die beiden nativen Clockmatrizen bestimmen gleichzeitig einen einzigen semisimplen Singlet-Kanal:
\[
\mathfrak{su}(2)_\lambda,\qquad
\lambda=(-\tfrac12,\ldots,-\tfrac12).
\]
Seine Einbettung in E8 hat Index eins. Dieser Weg vergleicht keine bloßen Dimensionszahlen; er verwendet Klassenprodukt, Clutching, tatsächliche Gruppenwirkung und die normierte Liealgebraform.

## Expliziter Test des Klassenprodukts

Setze
\[
P(n)=\tfrac12(I+n\cdot\sigma),\quad
u(z,n)=I+(z-1)P(n),\quad
U(z,n)=u(z,n)u(z,n_0)^{-1}.
\]
Es gilt det U=1. An z=1 und n=n₀ ist U=I; die Karte geht deshalb auf dem Smash-Produkt S¹∧S²≅S³ in SU(2)≅S³ über.

Mit z=exp(iθ), n₀=e₃, a=cos(θ/2), b=sin(θ/2) lautet U=q₀I+iq·σ:
\[
q_0=a^2+b^2n_z,\qquad
q=ab(n-e_3)+b^2(n\times e_3).
\]
Der Checker berechnet die normierte geometrische Gradzahl direkt. Für n=(sinβ cosφ,sinβ sinφ,cosβ) ist der Volumen-Jacobian
\[
J=\tfrac12(1-\cos\beta)\sin\beta\sin^2(\theta/2),
\qquad
\int J\,d\theta\,d\beta\,d\phi=2\pi^2.
\]
Der geometrische Grad ist somit exakt +1 in der angegebenen Orientierung. Der Betrag der geklammerten E8-Grundklasse ist eins. Vorzeichen von Spinor-Chiralität, ungeradem Cherncharakter und Clutching-Konvention werden dadurch nicht stillschweigend gleichgesetzt.

Der Klassenweg verwendet Standard-Bott-Theorie; explizite Clifford-Repräsentanten und verbindende Abbildungen sind in Schulz-Baldes–Stoiber, [The generators of the K-groups of the sphere](https://arxiv.org/abs/2403.18944), beschrieben. Die hier angegebene Karte und ihr Grad wurden zusätzlich symbolisch gerechnet.

## Auswahl des nativen SU(2)-Kanals

Aus der vorhandenen Chevalley-Algebra werden sämtliche Wurzelräume ausgewählt, die mit den Standardmodell-Wurzelgeneratoren und Hyperladung kommutieren. Das ergibt 20 Wurzeln: zwölf A₃-Wurzeln, vier C-Felder und ihre vier Gegenfelder. Sie bilden mit vier Cartanrichtungen ein A₄-System. Alle 120 nichtverschwindenden wurzelwertigen Klammerausgaben bleiben im System.

Die zusätzliche Hyperladungs-Cartanrichtung gehört ebenfalls zum vollen SM-Kommutanten. Dieser lautet auf Liealgebraebene su(5)⊕u(1)Y; die Aussage betrifft keine globale Gruppenkonstruktion oder Eichsymmetrieherleitung.

Die native Clockwirkung auf dem zugehörigen Fundamental zerfällt in zwei triviale Linien und einen irreduziblen Dreierblock:
\[
M_5=\operatorname{diag}(1,1,M),\qquad
U_5=\operatorname{diag}(1,1,U).
\]
Der exakt berechnete gemeinsame Kommutant in gl(5) ist Mat₂⊕ℂI₃. Nach Spur- und Realitätsbedingung bleibt in su(5) die Algebra su(2)⊕u(1); zusammen mit Hyperladung bleiben zwei abelsche Richtungen neben genau einem semisimplen su(2).

Die beiden Wurzeln dieses su(2) sind ±λ. Es gilt κE8(Hλ,Hλ)=2=κSU2(H,H). Das ist die Index-eins-Einbettung. Die Cartanrichtungen sind weiterhin nicht als ganze Quellenwahl eindeutig.

## Verbindung mit den bereits untersuchten Antworten

Im schon gewählten nativen Level-eins-Strommodul liefert dieselbe Hλ-Richtung zwei Antworten:
\[
J=H_\lambda/\sqrt2,\quad
\langle J_nJ_{-n}\rangle=n,\quad
X_n=J_n/n,\quad
\langle X_nX_{-n}\rangle=1/n;
\]
\[
\sum_{\beta\in64}(\lambda\cdot\beta)\,\beta\beta^{T}=-R,
\qquad \operatorname{tr}R=0,\qquad \|R\|^2=48.
\]
Hier sind die 64 Wurzeln das deklarierte chirale Quellenpaket. Das ganze E8-Wurzelsystem enthält 240 Wurzeln. Die Gleichung ist eine native, ladungsgewichtete Paarantwort; sie leitet noch keine elektromagnetische Wechselwirkung, Vakuumauswahl oder Nahtzeit ab. Das Vorzeichen folgt aus der Orientierung von λ.

Der zuvor geprüfte direkte interne Wurzelloop ist auf den getesteten endlichen Fasern kein P1-Einheits-Clutching. Die Determinantenwindungen sind 0 auf dem vollständigen Paket und −12 auf den 48 zurückbehaltenen Wurzeln. Diese endlichen Spuren entscheiden keinen separat definierten Fredholm-Spektralfluss.

## Entscheidende Operatorprüfung

Eine topologische Grundklasse ist noch keine physikalische Determinantenlinie. Im E8-Adjungierten ergibt die tatsächliche Wurzelwirkung
\[
\operatorname{Tr}_{248}H_\lambda^2
=2\cdot2^2+112\cdot1^2=120.
\]
Im SU(2)-Doppelten ist die Vergleichsspur 2. Eine gewöhnliche adjungierte E8-Diracdeterminante trägt daher den Faktor 60 gegenüber der Grundform. Man darf die Unit-Bott-Klasse nicht zu einer Unit-Anomalie dieser endlichen Operatorfamilie erklären. Die Darstellungsspuren und Dynkinindex-Normierung behandelt [Baek–Neher–Zainoulline, Abschnitt 4](https://ems.press/content/serial-article-files/26143); der konkrete Faktor wurde hier aus allen 240 nativen Wurzeln gerechnet.

Auch die Normaloperatoren sind zu unterscheiden. Auf CP¹ ist K=O(−2), die Spinorlinie O(−1). Der ungetwistete kompakte Spin-Diracoperator hat Index null. Ein Twist mit O(k) liefert Index k; der kanonische komplexe SpinC-Dolbeaultoperator hat Index eins. Der normale Thom-/Bott-Symboltyp ist wiederum ein anderer Gegenstand. Ein Odd-Indexpairing der obigen U-Karte mit einem passenden dreidimensionalen Produkt-Diracoperator kann Betrag eins liefern. Welches dieser Pairings der tatsächliche Nahtoperator realisiert, ist aus den gelesenen ursprünglichen Daten nicht hergestellt. Zum ungetwisteten sphärischen Operator siehe [Abrikosov, Dirac operator on the Riemann sphere](https://arxiv.org/abs/hep-th/0212134).

Die Viertelrotation z→iz bewirkt
\[
U(iz,n)=u(i,n)\,U(z,n)\,u(i,n_0)^{-1}.
\]
Das ist keine Konjugation mit einer konstanten nativen Clockmatrix. Die Rotation erhält auch nicht den gewählten kollabierten Basisschnitt des Smash-Produkts. Eine geeignete Bundle-/Transportwirkung muss geprüft werden; eine physikalische Clock-Identifikation ist nicht bewiesen.

## Korrektur eines früheren geometrischen Anschlusses

Im archivierten Randtheorem wird Rτ=2P_C−I ohne eigenständige Definition der transportierten Rτ-Wirkung behauptet. Als Gleichheit mit einem lokalen geometrischen Spinlift scheitert sie schon am Dirac-Hauptsymbol: C(−ξ)=−C(ξ), während ein lokaler Lift am selben Randpunkt nicht von ξ abhängt.

Der exakte Produktkollar D=iσ₁∂x+σ₂(−i∂θ) mit antiperiodischen Moden und Dirac-erhaltendem Lift R=σ₂ liefert RP_CR⁻¹=I−P_C. Die passende Normal-Clifford-Beziehung komplementiert die Polarisation. Dieser Standardbefund ist in [Kirk–Lesch, Abschnitt 2](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/kirklesch2.pdf) und [Guillarmou–Moroianu–Park, Proposition 12](https://arxiv.org/abs/1009.3179) formuliert.

Das widerlegt keine gesamte TFPT-Konstruktion. Es entscheidet die wörtliche lokale Identifikation; ein gesonderter nichtlokaler Transport müsste definiert und hergeleitet werden. Der endliche 3+2-Träger kann unter zwei komplementierenden Involutionen nicht gleichzeitig invariant sein: diese geben eine Bijektion zwischen beiden Polarisationsteilen. Eine entsprechende Kompression beziehungsweise relative Indexabbildung ist deshalb explizit nötig, falls diese Invarianz vorausgesetzt werden soll.

## Präzise verbleibende Herleitung

Es fehlt die Operatorabbildung, welche den wirklichen P1-Response mit dem hier konstruierten Bott-/E8-Grundklassenweg identifiziert und zugleich die native Clockwirkung, positive Quantisierung, Zustandswahl und Kopplung realisiert. Eine passende topologische Karte ersetzt diese Operatorgleichheit nicht.

Die Originalquellen liefern hierfür unterschiedliche, klar getrennte Eingaben:

- P1-K¹-Typ und Spektralfluss: 03_em_flavor_source.tex Zeilen 11–25 und 02_carrier_source.tex Zeilen 2862–2873.
- Kompaktifizierte normale S²: tfpt_1_architecture_e8.tex Zeilen 1639–1655 und origin_theory.tex Zeilen 1187–1194.
- Allgemeiner Kollargenerator und Zustand als Eingabe: 01_boundary_kernel_source.tex Zeilen 69–90 und 480–523.
- Diskrete Defektminimalität: dort Zeilen 543–559. Sie minimiert keine hier ausgewählte normale Operatorfamilie.
- Nativer Current-Rahmen: current_block_geometry.py, actual_current_source, mit Chevalley-Algebra v498 und nativen U/M-Matrizen.

P1 bleibt im Ledger ein deklarierter Ursprungstyp. Die obige Arbeit leitet aus diesen Typen einen Klassenweg und einen ausgewählten nativen semisimplen Kanal ab; sie beweist nicht den tieferen Ursprung von P1 oder einen physikalisch vollständigen 3+1D-Compiler. Es gibt daraus noch keine Primzahldynamik, keinen RH-Beweis und keine überlegene AI-Architektur.

Alle vier beigefügten Checker liefen erfolgreich mit exakter Algebra. Physikalische Gates wurden nicht als geschlossen markiert; Core, Ledger und Originalpapiere wurden nicht verändert.

