# QTheta, vorzeichenempfindliche Higgsnorm und Rückkehr

11. November 2026 · TFPTQTHETABACKGROUND20261111

Die vorangegangene Verfeinerung bei konstantem Hintergrund ließ das Vorzeichen der zusätzlichen Higgswelle in allen geschlossenen groben Spektralformen verschwinden. Wenn die drei ursprünglichen Hintergrundharmonischen des Q Assemblers mitlaufen, trägt bereits eine Rückkehr zwischen zwei groben Plätzen das Vorzeichen. Zusätzlich kann die schon definierte eichinvariante Higgsnorm dieselbe Information tragen. Damit lässt sich die Bedingung für eine vollständige Integration der einzelnen untersuchten Feldrichtung ohne einen zusätzlichen geladenen Feldoperator formulieren.

## Geordnete Rückkehr am unveränderten Assembler

Der [Originalassembler](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime/assemble.py:70) setzt aus vier reellen Higgskoordinaten die vorhandene v252 Matrix und die räumlichen Fourierblöcke zusammen. Im vorangegangenen Test standen die drei ursprünglichen räumlichen Vorschlagskoeffizienten auf null. Hier ist der tatsächlich im Domainpaket definierte Hintergrund mit den Genen [2,2,2] eingeschaltet. Im ursprünglichen Feldchart lautet er mit derselben reellen Richtung für den Testkoeffizienten a
```text
H_1(x) = b cos x + i c sin x + a cos 8x,
H_2(x) = 1 + i d cos 2x,
b = 13/100, c = 9/100, d = 7/100.
```
Die ursprünglichen vier Plätze sind −2,−1,0,1. Der Test benutzt wie zuvor eine temporäre Verfeinerung auf die 14 Plätze −7,…,6 und a=−u,0,+u mit u=13/100. Das Originalpaket, die Viertelphase, die v252 Kopplungen und die Majoranamatrix bleiben gleich. Die neu zugelassene Amplitude und die feinere räumliche Abschneidung sind die Testannahmen.

Für jeden dieser Werte gilt `P h(a) P = P h(0) P` im unveränderten groben räumlichen Block. Ein ungerader Beitrag in a benötigt mindestens einen Sprung um acht. Mit insgesamt höchstens drei Schritten ist der kleinste verbleibende Abstand mindestens vier, da die weiteren Schritte höchstens zwei beitragen. Der grobe Raum hat maximal Abstand drei. Die Momente bis zur dritten Potenz sind deshalb gerade in a.

Im vierten Moment gibt es für den Übergang −2→0 genau drei zulässige Anordnungen der Schritte +8,−2,−2,−2. Ein vierter Weg verließe das Intervall nach −8. Sei B der feine Block pro Einheit a für +8 und D der ursprüngliche Block für −2. Dann gilt exakt
```text
½[(h(a)^4)_(0,−2) − (h(−a)^4)_(0,−2)]
   = a (D³ B + D² B D + D B D²).
```
Der Koeffizient hat Frobeniusnorm 0,0145192056102. Der Vergleich mit der unabhängigen Multiplikation des ganzen feinen Operators weicht um 3,8·10⁻¹⁸ ab. Die zweite Potenz enthält auf dem groben Raum weiterhin den positiven geraden Hin und Rückweg des direkten achten Übergangs. Seine zusätzliche Spur bei |a|=0,13 beträgt 3,09470025650. Der vierte Moment ist ein räumliches Operatormoment im gegebenen Chart.

## Farbzentrum, Cauchy Schnitt und Halbwendung

Der v252 Farbzentrumsoperator wurde durch die reale Bimodulwirkung berechnet. Der Koeffizient des vierten Moments kommutiert damit und mit der endlichen Graduierung. Seine Norm auf dem 24 dimensionalen farbneutralen Einteilchenblock des konkreten Übergangs −2→0 beträgt 0,0140632641712. Die Farbprüfung ersetzt den vollständigen Fock Selektor nicht.

Für die räumliche Halbwendung x→x+π zusammen mit der elektromagnetischen Drehung exp(iπQ_em) gilt auf dem ganzen betrachteten feinen Operator
```text
G h(a) G* = h(−a),   G_n = (−1)^n exp(iπQ_em).
```
Die elektromagnetische Drehung ändert das Vorzeichen der oberen Higgskoordinate. Die Halbwendung dreht zuvor die beiden ersten Harmonischen um und erhält die zweite und achte Harmonische. Daher hält die gemeinsame Wirkung den ganzen ursprünglichen Hintergrund fest und kehrt nur a um. Die Gleichung wurde am gesamten v252 Feldoperator verifiziert. G erhält den groben Raum als Unterraum, wirkt darin jedoch nicht als Identität. Alle unmarkierten groben und vollen Spektralspuren sind aus diesem Grund gerade in a. Der ungerade Matrixanteil liegt zwischen entgegengesetzten G Vorzeichen im farbneutralen Träger. Für einen Selektor, dessen ganzer Bildraum tatsächlich in einem einzigen G Eigenraum liegt, verschwindet die beidseitig ausgewählte ungerade Einfügung. Eine positive Rückkehr über den entgegengesetzten Raum kann zugleich erhalten bleiben, falls dort keine zusätzliche Zwischenprojektion eingesetzt wird.

Die vorhandene feste [Zweirandform des Produktkragens](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/chat-triangulation-20261108/TFPT-SCHWACHE-GREEN-RANDANTWORT-20261109.txt:47) lässt sich auf derselben ursprünglichen groben Hintergrundmatrix auswerten. Mit Kragenlänge s=1/8, A=(1+exp(−2sh))⁻¹ und V=(sqrt(A),sqrt(1−A))ᵀ liest sie denselben Operator X auf beiden Flächen als
```text
V* diag(X,X) V = sqrt(A) X sqrt(A) + sqrt(1−A) X sqrt(1−A).
```
In der Eigenbasis von h multipliziert dies jedes Matrixelement mit einem strikt positiven Faktor. Die Ablesung ist im endlichen Chart injektiv. Für den gesamten ungeraden vierten Moment pro Einheit a im farbneutralen Träger beträgt die Norm vor der Cauchy Auswahl 0,172345344722 und danach 0,167719444725. Bei entgegengesetzten Flächenzeichen ergibt sich 0,0744065667011. Der kleinste gleichseitige Multiplikator beträgt 0,781716661137. Der Unterschied der beiden Randleser zeigt, warum die wirkliche Randmarkierung mitgeführt werden muss. Die Rechnung benutzt die vorhandene Produktform an der festgehaltenen, gegebenen räumlichen Hintergrundmatrix. Sie identifiziert die physische Randpaarung, das physische QΘ und den Gapoperator nicht mit diesem Cauchy Chart.

## Originale Higgsnorm als direktes Integral derselben drei Momente

Der bereits definierte [ursprüngliche U Quelleneinsatz](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/original-paper-source-reconciliation-20261102/shorttime_20261105L/original_marks/ORIGINAL-U-Y-AND-SOURCE-SUBTRACTION.txt:6) benutzt die verschmierte eichinvariante Norm U(x)=H(x)*H(x). Schreibe U_k für den reellen Cosinuskoeffizienten mit Normierung 1/π auf [0,2π]. Für den hier untersuchten Hintergrund gilt unmittelbar
```text
U_7 = b a,       U_16 = a²/2.
```
Der siebte Anteil entsteht aus 2ab cos x cos 8x. Der sechzehnte stammt aus a² cos²8x. Die native Feldfaltung gibt für die drei Amplituden die Paare (−0,0169; 0,00845), (0; 0), (+0,0169; 0,00845). Dieselben Koeffizienten lassen sich durch lineare Auslesung der Higgsnorm an antipodalen Orten bestimmen. Die Norm wird auf der gegebenen Feldrichtung zu einem eindeutigen Leser für a, obwohl die ursprüngliche grobe Diracmatrix diese Richtung vollständig verfehlt.

Für eine tatsächliche gemeinsame Randfunktion L[S] auf dieser konditionalen Faser reichen daher die bereits typisierten Einfügungen 1, U_7 und U_16 am selben Zeitort. Mit derselben Ordnung, Randpaarung, dem vollen äußeren P_adm und derselben Leerreferenz lauten die Momente
```text
m0 = L[1],       m1 = L[U_7]/b,       m2 = 2 L[U_16].
w+ = ½(m2/u² + m1/u),
w− = ½(m2/u² − m1/u),
w0 = m0 − m2/u².
```
Die Gleichungen gelten komponentenweise auch für Matrixformen, sofern alle drei Eingaben durch dieselbe lineare Originalkontraktion auf der gegebenen dreipunktigen Faser dargestellt sind. Weder die äußere Auswahl noch eine innere QΘ Einfügung muss dazu mit a vertauschen. Die Reihenfolge muss in allen drei Eingaben dieselbe bleiben. Ein anderes Feldensemble mit weiteren gleichzeitig variablen Moden erzeugt weitere Beiträge zu U_7 und U_16. Die obige Inversion ist dann nur nach tatsächlicher Konditionierung auf die untersuchte Richtung anwendbar. Wenn die physische volle Randfunktion unter G invariant ist, folgt auf dieser Faser m1=0 und damit w+=w−. Die Invarianz ist an tatsächlicher Zeit, Rand und Leerreferenz zu prüfen.

Der offene Matrixübergang und seine direkte adjungierte Rückkehr bleiben als unabhängiger Vergleich verfügbar. Die beiden U Einsätze liefern einen Weg über die vorhandene neutrale Quellfamilie und vermeiden eine zusätzliche Annahme, dass ein nackter geladener Einteilchenoperator zur admissiblen Algebra gehört.

## Tatsächlicher Lauf und erste noch fehlende Originaldaten

Die [Definition des ursprünglichen Selektors](/Users/stefanhamann/Projekte/tfpt-theoryv4/documents/archive--tfpt-45--source_extracts--04_qft_source/sections/002-section-closure-theorems-from-the-admissibility-complex.tex:78) verlangt Nahtpolarisation, Farbzentrum, QΘ und den Gap von D_adm im selben Raum. Der installierte Q Adapter liefert auf der unveränderten Domain kein ursprüngliches Feldpaket. Der native Integrationslauf endete mit Exitcode 2 und null Originalatomen. Sein tatsächlicher Konsument hat vier feste Felder Z0, ZU, ZY und ZUY für mitgelieferte Formen. Er erzeugt die beiden verschmierten U Formen, den ganzen Selektor, die Zeitverschiebung und die Referenzmasse nicht selbst. Die [euklidische Gitterkonstruktion von Osterwalder und Seiler](https://doi.org/10.1016/0003-4916(78)90039-8) liefert hierfür das allgemeine Beispiel einer tatsächlich konstruierten gemeinsamen reflektierten Quelle, aber keine Werte für das TFPT Paket.

Der kleinste wirkliche Originalaufruf verwendet den leeren Wert und die beiden verschmierten U Werte mit gleichem Rand, tatsächlicher Zeit, demselben vollen Selektor und gebundener Referenz. Für den Vergleich mit dem ursprünglichen Auftrag sind danach auch der offene Eintrag und die vollständige Rückkehr am selben Schnitt zu kontrahieren. Erst diese Originalformen entscheiden die physische Verteilung. Z₀, F(0), C_UY und H_tr,pos bleiben unberechnet.

[Reproduzierbarer Prüflauf](check_qtheta_background_return_20261111.py)

## Der offene Übergang und die vollständige Rückkehr im vorhandenen Kragen

Im feinen farbneutralen Träger sei R der räumliche Viermodenschnitt und B=(I−R)(h(a)−h(0))R der offene Übergang in die zwei zusätzlichen Partnerplätze. Der ursprüngliche nichtkonstante Hintergrund wird in der festen Zweirandform mit 14 Plätzen als h(0) eingesetzt. Die Kragenlänge ist s=1/8. Schreibe wieder V für ihre Cauchy Einbettung und P_C=VV*. Auf den beiden Flächen trägt derselbe Eintrag B den Operator O=diag(B,B). Dann sind
```text
offener Eintrag            V* O V,
vollständige Rückkehr       V* O* O V,
Rückkehr mit Zwischenauswahl (V* O V)* (V* O V).
```
Ihre Differenz ist exakt `V* O* (I−P_C) O V`. Sie ist eine positive Gramform des verworfenen Zwischenanteils. Dieser Vergleich behält die beiden Randflächen in derselben Ausrichtung wie die gleichseitige Probe oben. Auf dem farbneutralen Einteilchenraum der Verfeinerung und bei a=0,13 ergibt der offene Eintrag die Frobeniusnorm 1,09817167160. Die Spur der vollständigen Rückkehr beträgt 1,41772687578. Eine zusätzliche Cauchy Projektion in der Mitte reduziert sie auf 1,20598102030. Der Verlust 0,21174585549 wurde unabhängig über die orthogonale zweite Cauchy Einbettung berechnet. Auch wenn auf beiden äußeren Rändern dieselbe Auswahl feststeht, darf eine weitere innere Auswahl nicht mit der ursprünglichen offenen Geschichte verwechselt werden.

Diese Werte sind Matrixkontraktionen am vorhandenen Kragen und farbneutralen Einteilchenraum. Die leere Matrixkontraktion in der normierten Cauchy Einbettung ist V*V=I. Die leere physische Masse verlangt dagegen weiterhin den tatsächlichen unnormierten Randwert samt Feldintegration und Leerreferenz. QΘ und der physische Gapoperator wurden nicht aus der vorhandenen Kragenform abgeleitet.
