# Ein Stromoszillator: Zustandsauswahl und korrekter positiver Elternoperator

## Ergebnis

Die native Vorgabe einer regulären Stromalgebra mit voller KMS Dynamik bei
\(\beta\omega=4\pi\) legt die kanonische Thermofeldamplitude tatsächlich fest:

\[
q=e^{-2\pi}.
\]

Für diesen einzelnen Oszillator ist die Gaußform des Zustands eine Folgerung der
KMS Bedingung. Die bloße mittlere Besetzung genügt dafür nicht. Die Identifikation
der beiden physikalischen Blattalgebren bleibt eine gesonderte Voraussetzung.

Der passende quadratische Elternoperator auf dem Paarsektor ist

\[
K_\lambda=\lambda A_q^\dagger(N+1)A_q,
\qquad (A_q f)_n=qf_n-f_{n+1}.
\]

Er hat dieselbe geometrische Nullmode wie der bisherige flache Hardy Operator
\(\lambda A_q^\dagger A_q\), aber eine andere Dynamik. Sein vollständiges Spektrum
und sein positives Randmaß sind exakt berechenbar. Innerhalb einer benannten
quadratischen Klasse wählen relative Ladungserhaltung und Blatttausch diesen
Elternoperator bis auf seine Skala eindeutig.

Die Zertifikate stehen in
[process_gaussian_parent.py](../../../../research/originlab/process_gaussian_parent.py).
Die begleitenden 13 exakten Tests bestanden in 1,068 Sekunden. Sie kontrollieren
Faktoren, Normen, Rekursionen und Gegenbeispiele; die allgemeinen Aussagen folgen
aus den nachstehenden Beweisen, nicht aus einem endlichen Scan.

## 1. Was Z48 tatsächlich identifiziert

Geprüft wurden die unveränderten Snapshots von Commit
`0172a74bd`, insbesondere `z48_zuordnung_quelle.py` und `process_rules.py`.
Z48 überprüft die Umnummerierung

\[
k=4(r+\alpha),\qquad r\mapsto r+1\Longrightarrow k\mapsto k+4.
\]

Die Identität der Modenetiketten ist richtig. Sie identifiziert jedoch nicht den
Stromoperator auf dem Zustandsraum mit einem unilateralen Shift.

In einer CAR Konvention mit
\(J_{-1}=\sum_r c_{r+1}^\dagger c_r\) gilt

\[
[J_{-1},c_r^\dagger]=c_{r+1}^\dagger.
\]

Das ist eine Wirkung durch Kommutator auf Feldoperatoren. Auf Zuständen bleibt
der zusätzliche Term

\[
J_{-1}c_r^\dagger\Omega
=c_{r+1}^\dagger\Omega+c_r^\dagger J_{-1}\Omega.
\]

Der normierte Stromoszillator erfüllt dagegen

\[
[J_1,J_{-1}]=1,\quad
|n\rangle=\frac{J_{-1}^n}{\sqrt{n!}}\Omega,\quad
J_{-1}|n\rangle=\sqrt{n+1}|n+1\rangle.
\]

Schon \(\|J_{-1}^2\Omega\|^2=2\) widerspricht der Isometrienorm
\(\|S^8\Omega\|^2=1\). Ein exakter CAR Zeuge mit vier Fermionmoden und zwei
anfangs besetzten tiefsten Moden bestätigt die Normfolge \(1,1,2\).
Die ersten beiden Schritte benötigen keine abgeschnittenen Zustände.

Der korrekte Shift entsteht erst als Polarphase:

\[
V=J_{-1}(N+1)^{-1/2},\quad V^\dagger V=1,
\quad VV^\dagger=1-|0\rangle\langle0|.
\]

Eine zusätzlich angegebene Isometrie \(U|n\rangle=e_{4n+j}\) identifiziert
\(UVU^\dagger\) mit \(S^4\) auf genau einem Restklassensektor.
Sie identifiziert nicht \(J_{-1}\) selbst mit \(S^4\), und sie ersetzt nicht
die vollständige Stromalgebra: deren Vakuummodul hat die Gradmultiplizitäten
\(1,1,2,3,5,7,\ldots\), die isolierte \(J_{-1}\) Leiter jeweils nur eins.
Die zugrunde liegende CAR und Stromalgebra findet sich bei
[Alexandrov und Zabrodin, Abschnitt 2.9](https://arxiv.org/abs/1212.6049).

## 2. Volle KMS wählt den Zustand des einzelnen Oszillators

Voraussetzungen sind eine reguläre Darstellung eines CCR Paars
\(a=J_1\), \(a^\dagger=J_{-1}\), die festgelegte Dynamik
\(\alpha_t(a)=e^{-it}a\) und die volle KMS Bedingung bei \(\beta>0\).
Auf dem regulären endlichen Oszillatorfaktor liegt eine normale Zustandsdichte
vor. Stationarität macht sie diagonal im nicht entarteten Zahloperator:

\[
\rho=\sum_{n\ge0}p_n|n\rangle\langle n|.
\]

Die beschränkten Matrixeinheiten
\(E_{n,n+1}=|n\rangle\langle n+1|\) sind analytische Elemente der Dynamik.
Die KMS Relation mit \(A=E_{n,n+1}\), \(B=E_{n+1,n}\) ergibt

\[
\omega(A\alpha_{i\beta}(B))=\omega(BA)
\quad\Longrightarrow\quad e^{-\beta}p_n=p_{n+1}.
\]

Normierung liefert eindeutig

\[
\rho_\beta=(1-e^{-\beta})e^{-\beta N},\qquad
\langle N\rangle=\frac1{e^\beta-1}.
\]

Der thermische Oszillatorzustand ist damit gaußsch; eine separate Gaußannahme
für diesen Modus ist überflüssig. Dies ist keine Aussage über alle Ladungssektoren
oder den gesamten E8 Zustand. Insbesondere ersetzt die Hypothese einer
physikalischen KMS Zeit keine Herleitung dieser Zeit aus den Marken.

Für die in H6h26 vorgegebene Fermibesetzung
\(c_r=(1+e^{\beta(r+\alpha)})^{-1}\) gilt exakt

\[
c_{r+1}(1-c_r)=\frac{c_r-c_{r+1}}{e^\beta-1}.
\]

Die Summe über alle ganzen \(r\) teleskopiert wegen
\(c_{-\infty}=1\), \(c_{+\infty}=0\) zu \((e^\beta-1)^{-1}\), unabhängig
von \(\alpha\). Das bestätigt den Strommittelwert ohne endlichen Abschneider.
Die Übertragung der **vollen** KMS Struktur auf den regulären Stromoszillator
ist stärker als diese einzelne Summenidentität.

### Exaktes Gegenbeispiel zur Auswahl nur über den Mittelwert

Setze \(r=e^{-\beta}\), \(Z=1-r\), \(\epsilon=Zr/4\) und

\[
\rho_{\rm alt}=\rho_\beta+
\epsilon(P_0-2P_1+P_2).
\]

Dieser Zustand ist positiv, normiert und treu. Sein Mittelwert \(\langle N\rangle\)
bleibt unverändert. Damit bleibt auch die thermische Zweipunktrelation zwischen
\(\langle a^\dagger a\rangle\) und \(\langle aa^\dagger\rangle\) unverändert.
Die volle KMS Bedingung scheitert aber bereits an

\[
p'_1-rp'_0=-(2+r)\epsilon\ne0,
\]

und \(\langle a^{\dagger2}a^2\rangle\) weicht um \(2\epsilon\) vom Wick Wert ab.

## 3. Kanonische Purifikation und ihre Identifikationsvoraussetzung

Bei festgelegter konjugierter rechter Algebra ist die kanonische Purifikation

\[
|\eta_q\rangle=\sqrt{1-q^2}\sum_{n\ge0}q^n|n,n\rangle,
\qquad q=e^{-\beta/2}.
\]

Mit dem tatsächlich vorgegebenen \(\beta=4\pi\) ist somit
\(q=e^{-2\pi}\) korrekt. Hier liegt kein Fehler um den Faktor zwei vor.
Ohne Festlegung der rechten Algebra beziehungsweise der kanonischen positiven
Standardform ist die Purifikation nicht eindeutig: bereits
\(q^n\mapsto q^ne^{in\theta}\) bewahrt beide reduzierten thermischen Zustände
und die Zweimoden Gaußform. Siehe die Definition und Phasenfamilie in
[Cottrell et al., Gleichungen 1.1 und 1.2](https://arxiv.org/abs/1811.11528).

Auch die Zuordnung \(|n\rangle\mapsto|n,n\rangle\) muss die Observablen
mitführen. Als Zustand **eines** Oszillators ist \(\sum q^n|n\rangle\) für
\(q>0\) nicht gaußsch. Die notwendige Gleichung eines reinen Gaußzustands
\((ua+va^\dagger-\alpha)|\psi\rangle=0\) erzwingt aus ihren ersten drei
Komponenten
\(uq(\sqrt3+\sqrt2-3)=0\), im Widerspruch zu
\(|u|^2-|v|^2=1\). Die Zweimoden Lesart ist daher eine substantielle
Algebraidentifikation, keine bloße Umbenennung von Koeffizienten.

## 4. Ein bedingt eindeutiger quadratischer Elternoperator

Definiere \(Z=1-q^2\) und

\[
d_L=\frac{a-qb^\dagger}{\sqrt Z},\qquad
d_R=\frac{b-qa^\dagger}{\sqrt Z}.
\]

Diese Operatoren erfüllen die beiden unabhängigen CCR und vernichten
\(|\eta_q\rangle\). Der positive Elternhamiltonian ist

\[
H_q=\Omega(d_L^\dagger d_L+d_R^\dagger d_R)
=\frac{\Omega}{Z}
\big[(1+q^2)(N_L+N_R)+2q^2-2q(ab+a^\dagger b^\dagger)\big].
\]

Die Standardkonstruktion und ihre Lücke werden in
[Cottrell et al., Abschnitt 3.1](https://arxiv.org/abs/1811.11528) angegeben.
Hier werden anschließend ihr Paarsektor und dessen Randmaß ausgewertet.

**Eindeutigkeit innerhalb der benannten Klasse.** Verlange einen positiven
quadratischen Elternoperator auf genau diesen zwei CCR Moden mit dem angegebenen
einzigen Grundzustand, Erhaltung von \(Q=N_L-N_R\) und Blatttausch.
Im \(d\) Vakuum erzwingt die Grundzustandseigenschaft das Verschwinden linearer
und anomaler Erzeugungsterme. Es bleibt \(d^\dagger Gd\) mit \(G>0\).
Die \(Q\) Symmetrie macht \(G\) diagonal; Blatttausch setzt die beiden Einträge
gleich. Also bleibt nur \(\Omega>0\) sowie eine additive Energiekonstante frei.

In der kanonischen KMS Standardform können die Symmetrien genauer
formuliert werden. Dort sind
\(\Delta_{\rm mod}=\rho_\beta\otimes\rho_\beta^{-1}=e^{-\beta Q}\)
und die modulare Konjugation \(J_{\rm mod}\), welche die beiden Faktoren vertauscht
und komplex konjugiert, bereits gegeben. Die Forderungen
\([H,\Delta_{\rm mod}]=0\) und \(J_{\rm mod} H J_{\rm mod}=H\) erzwingen innerhalb derselben
quadratischen Klasse erneut \(G=\Omega I\).
Zusätzlich verlangt wird damit, dass der Präparationsoperator die
modularen Symmetrien respektiert; deren konkrete Form braucht nicht
gesondert geraten zu werden.

Dies ist eine Auswahlregel für einen **Präparationshamiltonian**. Sie setzt ihn
nicht mit dem ursprünglichen einseitigen KMS Generator gleich. Der Zustand
\(|\eta_q\rangle\) ist für jedes \(q\) unter Blatttausch invariant. Ein kleines
\(q\) beweist daher keine spontane Brechung dieser Symmetrie durch den Zustand.

## 5. Paarsektor, Eigenwerte und Normierung

Der Sektor \(Q=0\) besitzt die Basis \(|n,n\rangle\). Dort sind

\[
ab|n,n\rangle=n|n-1,n-1\rangle,\qquad
a^\dagger b^\dagger|n,n\rangle=(n+1)|n+1,n+1\rangle.
\]

Mit \(\lambda=2\Omega/Z\) wird die Einschränkung deshalb

\[
(K_\lambda f)_n=\lambda\big([(1+q^2)n+q^2]f_n
-q(n+1)f_{n+1}-qn f_{n-1}\big)
=\lambda A_q^\dagger(N+1)A_qf.
\]

Die mathematische Diagonalisierung ist vollständig. Für
\(\vartheta=\operatorname{artanh}q\) ist
\(U=\exp[\vartheta(a^\dagger b^\dagger-ab)]\) unitär,
\(d_L=UaU^\dagger\), \(d_R=UbU^\dagger\), und \([U,Q]=0\).
Damit sind \(U|m,m\rangle\) eine vollständige Eigenbasis des Paarsektors mit

\[
E_m=2\Omega m=\Delta m,
\qquad\Delta=\lambda(1-q^2).
\]

Die Disentanglung

\[
U=e^{qa^\dagger b^\dagger}\,
Z^{(N_L+N_R+1)/2}\,e^{-qab}
\]

ergibt \(\langle0,0|U|m,m\rangle=\sqrt Z(-q)^m\).
Das Randmaß ist folglich

\[
d\mu_0(s)=Z\sum_{m=0}^\infty q^{2m}\delta_{\Delta m}(ds),
\qquad \mu_0([0,\infty))=1.
\]

Alternativ liefern die Meixner Polynome
\(F_n(m)=q^n\,{}_2F_1(-n,-m;1;1-q^{-2})\) die normierte Spektraltransformation.
Ihre Orthogonalitätsnorm mit Gewicht \(Zq^{2m}\) ist eins; die Standardnormen
stehen in [DLMF, Tabelle 18.19.1](https://dlmf.nist.gov/18.19).
Wegen \(\sum_n[\lambda q(n+1)]^{-1}=\infty\) ist zudem kein freier
selbstadjungierter Randparameter am unendlichen Ende vorhanden. Zum verwendeten
Jacobi Kriterium siehe
[Kostenko und Teschl, Beweis von Theorem 2.1](https://mat.univie.ac.at/~gerald/ftp/articles/DispLag.pdf).

## 6. Vollständige Randresolvente und Niederenergiekoeffizienten

Für \(x>0\) gilt exakt

\[
C_{\rm raw}(x)=Z\sum_{m=0}^\infty\frac{q^{2m}}{x+\Delta m}
=\frac Z\Delta\Phi(q^2,1,x/\Delta),
\]

wobei \(\Phi\) die Lerch Funktion ist. Es gelten

\[
\frac Zx\le C_{\rm raw}(x)\le\frac1x,
\qquad\|(x+K_\lambda)^{-1}\|=\frac1x.
\]

Nach den ersten \(M\) Summanden ist der positive Rest durch
\(q^{2M}/(x+\Delta M)\) beschränkt. Das ist eine Fehlergarantie für numerische
Auswertung; die Positivität folgt bereits aus dem expliziten Maß.

Normiere \(C=C_{\rm raw}/Z\), damit das masselose Residuum eins ist. Für
\(r=q^2\), \(\ell=-\log(1-r)\) und \(H=C^{-1}\) folgt

\[
H=x-ax^2+bx^3-cx^4+O(x^5),
\]

\[
a=\frac\ell\Delta,\qquad
b=\frac{\ell^2+\operatorname{Li}_2(r)}{\Delta^2},\qquad
c=\frac{\ell^3+2\ell\operatorname{Li}_2(r)+\operatorname{Li}_3(r)}{\Delta^3}.
\]

Das massive Gesamtgewicht bleibt \(r/(1-r)\). Die UV Steigung von \(H\) ist
\(1-r\). Die Skala \(\lambda\), ihre physikalische Dimension und ihre Verbindung
mit \(G\) werden durch diesen Satz nicht festgelegt.

## 7. Gleicher Zustand ist noch nicht derselbe Prozess

Der flache Hardy Elternoperator hat für \(q>0\) ein beschränktes kontinuierliches
Band, während \(K_\lambda\) eine unbeschränkte arithmetische Punktleiter besitzt.
Beide sind positiv und haben denselben Nullvektor sowie dieselben masselosen und
massiven Gesamtgewichte. Sogar ihre ersten beiden Randmassenmomente stimmen:

\[
\mu_1=\lambda r,\qquad \mu_2=\lambda^2r(1+r).
\]

Erst beim dritten Moment zeigt sich der Unterschied:

\[
\mu_3^{\rm Gauß}=\lambda^3r(1+4r+r^2),\qquad
\mu_3^{\rm Hardy}=\lambda^3r(1+3r+r^2).
\]

Damit liegt ein exaktes Gegenmodell zur Identifikation allein über Nullmode,
Besetzung oder Gesamtresiduum vor. Gleichzeitig gibt es eine konkrete neue Route:
Die volle KMS Bedingung wählt den Einmodenzustand; die benannte quadratische Klasse
mit zwei Symmetrien wählt den gewichteten Elternoperator. Bei einer Erweiterung
zu vier metrischen Armen ist dieser Satz zunächst auf den symmetrischen Sektor
anzuwenden. Die drei anderen Sektoren sind die Dirichletkompression auf \(n\ge1\),
nicht automatisch dieselbe Zweioszillatorleiter.

Eine vollständige Theorie folgt daraus nicht. Die noch zu beweisende Brücke ist
die gemeinsame physikalische Zuordnung der Stromzeit, der beiden Blattalgebren,
des Präparationsoperators und des vierdimensionalen Feldprozesses.

Die Tomita Konjugation hat hier \(J_{\rm mod}^2=+1\). Sie wird nicht mit einem fermionischen Spinor oder Deckoperator vom Quadrat \(-1\) identifiziert. Ein mit K kommutierender antiunitärer Operator vom Quadrat \(-1\) auf demselben Hilbertraum würde die einfache Nullmode durch Kramers Entartung ausschließen. Der Übergang zwischen diesen Darstellungen ist ein eigener Nachweis.

## 8. Unabhängige Prüfung der Erweiterung auf vier metrische Arme

Geprüft wurden die neuen Dateien `process_weighted_star.py` und
`test_process_weighted_star.py`. Die folgenden Aussagen gelten für
\(0<q<1\) und \(\lambda>0\). Sie benötigen keine Annahme über einen
unendlichen wechselwirkenden Quantengravitationsgrenzfall.

### Vollständige Fourieridentität

Im endlichen Stern ordnet \(U\) die native Vertexbasis den vier Armcharakteren
zu; \(V\) ist die entsprechende Transformation auf den Kanten. Bereits
\(BU=VA\) ist eine exakte Identität des gesamten rechteckigen Operators.
Die neue Tiefengewichtung \(W\) hat an jeder Tiefe vier gleiche Einträge und
kommutiert deshalb mit \(V\). Folglich gilt auf allen Sektoren gemeinsam

\[
U^\dagger B^\dagger WB U=A^\dagger WA.
\]

Der zentrale symmetrische Block ist \(K_\lambda\); die drei anderen Blöcke
sind seine Dirichletkompression auf \(n\ge1\). Diese Aussage identifiziert
keinen dieser drei Blöcke mit einem zweiten Paaroszillator.

### Interlacing und endliche Anzahl von Moden unter fester Energie

Schreibe \(D\) für die Dirichletkompression. Die Entkopplung von \(e_0\)
ändert den selbstadjungierten Jacobi Operator nur durch eine beschränkte
Störung endlichen Ranges. Daher besitzt auch \(D\) kompakte Resolvente.
Für \(q>0\) ist \(e_0\) zyklisch und alle Gewichte
\(Zq^{2n}\) seines Spektralmaßes sind strikt positiv.

Auf jedem Intervall \((-\Delta(n+1),-\Delta n)\) gilt

\[
C'_{\rm raw}(x)
=-Z\sum_{m\ge0}\frac{q^{2m}}{(x+\Delta m)^2}<0.
\]

Die Randwerte sind links \(+\infty\), rechts \(-\infty\).
Die Schurkomplementformel identifiziert die jeweils genau eine Nullstelle
von \(C_{\rm raw}\) mit einem Eigenwert von \(D\), mit umgekehrtem Vorzeichen.
Damit hat \(D\) genau einen einfachen Eigenwert

\[
\Delta n<d_n<\Delta(n+1),\qquad n=0,1,\ldots.
\]

Insbesondere sind die drei Zusatzsektoren strikt positiv. Für
\(N(E)\), die Gesamtzahl der Eigenwerte der Masse zum Quadrat bis einschließlich
\(E\), inklusive Multiplizitäten, folgt

\[
N(E)\le4\left\lfloor\frac E\Delta\right\rfloor+4.
\]

An den Leiterpunkten \(E=k\Delta\) gilt sogar exakt \(N(E)=4k+1\).
Anders als beim flachen Elternoperator liegen somit nicht unendlich viele
massive Moden in einem beschränkten Massenband. Das allein kontrolliert
noch keine Summe virtueller Moden in einer wechselwirkenden Theorie.

### Physischer Endpunkt und starker Resolventengrenzfall

Für einen endlichen radialen Graphen mit \(m\) Kanten ist
\(K^{(m)}=\lambda A_m^\dagger W_m A_m\), und der letzte Diagonaleintrag ist
\(\lambda m\). Die Hauptkompression des unendlichen Operators hätte dort
stattdessen \(\lambda[m+q^2(m+1)]\). Die beiden endlichen Randbedingungen
sind also verschieden; der physische endliche Graph behält die exakte
geometrische Nullmode.

Trotzdem konvergieren die physischen Graphen stark im Resolventensinn.
Erweitere \(K^{(m)}\) auf \(\ell^2(\mathbb N_0)\) durch einen Nulloperator auf
dem orthogonalen Rest und nenne ihn \(\widehat K^{(m)}\).
Für jedes endlich getragene \(\psi\) stimmt
\(\widehat K^{(m)}\psi=K_\lambda\psi\) ab hinreichend großem \(m\) exakt.
Der Carleman Test zeigt, dass die endlich getragenen Vektoren einen
Operatorkern von \(K_\lambda\) bilden.

Für \(x>0\) und \(g=(x+K_\lambda)\psi\) gilt dann

\[
[(x+\widehat K^{(m)})^{-1}-(x+K_\lambda)^{-1}]g
=(x+\widehat K^{(m)})^{-1}(K_\lambda-\widehat K^{(m)})\psi
=0
\]

ab hinreichend großem \(m\). Die Menge solcher \(g\) ist dicht;
beide Resolventen haben Norm höchstens \(1/x\). Damit folgt die starke
Konvergenz auf dem gesamten Hilbertraum. Derselbe Beweis gilt für die
Dirichletblöcke und nach der festen Armtransformation für den ganzen Stern.

Das ist kein Beweis einer Konvergenz der Hamiltonoperatoren in Operatornorm
und kein Beweis einer unendlichen nichtlinearen Quantentheorie. Für die
gewählte Nullfortsetzung kann sogar die Resolventenkonvergenz nicht in
Operatornorm gelten: ihr unendlicher Nullraum behält die Resolventennorm
\(1/x\), während die Grenzresolvente kompakt ist. Der künstliche Nullraum
ist nur ein Hilfsmittel der gemeinsamen Einbettung, kein zusätzlicher
physischer Modeninhalt des endlichen Graphen.
