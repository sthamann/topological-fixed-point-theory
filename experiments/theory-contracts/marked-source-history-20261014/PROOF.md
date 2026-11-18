# Gemeinsame Markierung, Zustand und Ereignisse: exakte Fortsetzung des vorhandenen Ansatzes

14. Oktober 2026 · **PARTIAL** · UR.SOURCE.COHERENT_MARK_HISTORY.01

## Ergebnis und Auftrag

Der Auftrag ist eine Herleitung von Markierung, Zustand und Ereignisgesetz aus der ursprünglichen TFPT-Quelle. Diese vollständige Herleitung ist hier **nicht erbracht**. Insbesondere wird der vorhandene bedingte Nahtkanal nicht nachträglich zur ursprünglichen P1/P2-Quelle erklärt.

Der positive Fortschritt ist ein gemeinsamer kohärenter Ausdruck für Zustandsentwicklung und Markentransport auf den tatsächlichen Compileroperatoren. Ein klassisches Vollprotokoll ist dafür nicht nötig. Der weitere Versuch, daraus allein das vollständige Ereignisgesetz zurückzugewinnen, ist jedoch nachweislich nicht eindeutig: Zwei explizite stationäre, markiert-kovariante Prozesse liefern alle gleichen ungestörten Transfers, unterscheiden sich aber bei einer zeitlich aufgelösten Abfrage. Damit ist eine konkrete vermeintliche Abkürzung geschlossen; keine Unmöglichkeit der TFPT wird behauptet.

## 1. Welche Quelle tatsächlich verwendet wird

`marked-seam-cp-lift-20260921` liefert unter seiner ausdrücklich zusätzlichen orbit-skalaren CPTP-Annahme

\[
\Phi=\Pi_0+\tfrac13\Pi_5+\tfrac23\Pi_{10}
=\tfrac7{12}\operatorname{id}+\tfrac1{12}\sum_{j=1}^5\operatorname{Ad}_{g_j}.
\]

Die tatsächlichen fünf hermiteschen Clifford-Wörter auf C4 haben binäre Labels
0111, 1011, 1101, 1110, 1111. Die 16 ursprünglichen Pauli-Wörter W_v sind orthogonal für Tr(A*B)/4 und multiplizieren projektiv nach v+w in F2^4. Die Projektoren Pi betreffen Operatorräume, keine Teilchensorten oder Raumzeitdimensionen. Der frühere Contract enthält bereits den exakten positiven Anschluss an die ursprüngliche Drei-Zustands-Matrix B und B^6.

Der Kanal hat genau den normierten stationären Zustand rho*=I4/4: Seine beiden nichtskalaren Eigenwerte sind kleiner als eins. Das ist die eindeutige stationäre Dichte dieses endlichen Kanals, kein abgeleitetes physisches Vakuum und keine Auswahl einer anfänglichen Rang-3-Marke.

Auch seine algebraische Fünfermarkierung kann aus demselben Kanal zurückgelesen werden:

\[
\Pi_5=\tfrac92(\Phi-\mathrm{id})(\Phi-\tfrac23\mathrm{id}).
\]

Dieser spektrale Projektor rekonstruiert die schon vorausgesetzte markierte Operatorstruktur. Er beweist nicht deren Herkunft aus der rohen Naht und identifiziert sie nicht mit dem acht-dimensionalen E8-Stromprojektor P5.

## 2. Ein gemeinsamer kohärenter Ausdruck statt klassischem Vollrecord

Zunächst sei die **zusätzliche** Prozessannahme unabhängiger Wiederholungen festgehalten. Setze g0=I, p0=7/12 und pj=1/12. Für eine Geschichte h=(a1,...,an) gilt

\[
g_{a_n}\cdots g_{a_1}=\omega(h)W_{v(h)},\qquad |\omega(h)|=1.
\]

Die Phasen werden erhalten. Im kohärenten Lift

\[
V_n\psi=\sum_h\sqrt{p_h}\,\omega(h)W_{v(h)}\psi\otimes|h\rangle
\]

können die disjunkten Geschichten mit gleichem Endwort zusammengefasst werden:

\[
|e_v^{(n)}\rangle=q_n(v)^{-1/2}\sum_{h:v(h)=v}\sqrt{p_h}\omega(h)|h\rangle,
\quad V_n\psi=\sum_v\sqrt{q_n(v)}W_v\psi\otimes|e_v^{(n)}\rangle.
\]

Nur v mit q_n(v)>0 werden verwendet. Die e_v sind orthonormal. Aus der Fouriertransformation auf der ursprünglichen Pauli-Labelgruppe bzw. der Kanaldiagonalisierung folgt **pro Wort**

\[
q_n(0)=\frac{1+5(1/3)^n+10(2/3)^n}{16},
\]
\[
q_n(O_5)=\frac{1-3(1/3)^n+2(2/3)^n}{16},\qquad
q_n(O_{10})=\frac{1+(1/3)^n-2(2/3)^n}{16}.
\]

Für n=2 sind dies 3/8, 7/72 und 1/72. Alle 16 Gewichte sind dann positiv. Faltung mit p0>0 erhält strikte Positivität, daher ist der Choi-Rang für jedes n>=2 genau 16; bei n=1 ist er 6. Der Rang bezieht sich auf den abgeschlossenen Endpunktkanal, nicht auf den Speicher eines Prozesses mit offenen Zwischenschritten.

Für eine bereits vorgegebene Markierung P setzt derselbe Ausdruck

\[
P_v=W_vPW_v^*,\qquad
\widehat P_n=\sum_vP_v\otimes|e_v^{(n)}\rangle\langle e_v^{(n)}|.
\]

Dann gilt exakt

\[
\widehat P_n V_nP=V_nP.
\]

Zustand und mitgeführte Markierung bleiben zusammen unterstützt. Diese Rechnung benötigt weder Messung noch Dephasierung noch die Wahl eines klassischen Pointerregisters. Sie löst die endliche Transportfrage für ein gegebenes P. Sie wählt weder P noch V_n physisch aus. Die Formeln für den Endpunktkanal gelten auch für jede andere Vervollständigung mit denselben Phi^n; der angegebene spezielle Geschichtslift verwendet unabhängige Schritte.

Der unabhängige Nachtest zeigt zudem: Der minimale Ein-Schritt-Umgebungszustand am Fixpunkt ist diag(7/12,(1/12)I5). Sein Spektrum wählt keine fünf Pointerachsen. Drei verschiedene Instrumente mit sechs, zwei oder einem Ausgang ergeben denselben Kanal und dieselbe Markkovarianz.

## 3. Weshalb ein Endwort kein vollständiger Ablauf ist

Wähle zwei tatsächliche antikommutierende gi,gk und A=(I+i gk)/sqrt(2). Die Geschichten (0,i) und (i,0) enden ohne Einfügung beide bei gi. Nach Einfügung von A zwischen den Schritten ergeben sich jedoch gi A und A gi. Es gilt

\[
\operatorname{Tr}[(g_iA)^*(Ag_i)]/4=0.
\]

Die beiden unitären Operatoren sind nicht nur durch eine Phase verschieden; ihre Konjugationsabbildungen unterscheiden sich. Eine Verdichtung, die nur das Endwort behält, kann diesen offenen Zeitschlitz nicht beantworten. Das ist ein Operatorvergleich, kein Vergleich gleicher Dimensionen oder Spektren.

## 4. Zwei vollständige zeitliche Prozesse mit allen gleichen ungestörten Transfers

Dieser Abschnitt konstruiert **Gegenbeispiele zu einem Schluss aus reduzierten Daten**, keine vorgeschlagenen physikalischen TFPT-Modelle. Insbesondere wird nicht behauptet, beide seien vollständige P1/P2-Realisierungen oder erfüllten eine noch nicht gegebene ursprüngliche Mehrzeitwirkung.

Prozess A wählt die sechs vorhandenen Einzelschrittoperatoren unabhängig mit p.

Für Prozess B sei h ein geordnetes Tripel verschiedener Labels aus {1,...,5}. Delta_h hat +1 auf dem orientierten Dreizyklus, -1 auf dem umgekehrten Dreizyklus und sonst null. Definiere

\[
J_h(a,b)=p_ap_b+\epsilon\Delta_h(a,b),\qquad\epsilon=1/288.
\]

Alle Einträge sind positiv: Geänderte Einträge sind 1/144 +/- 1/288. Die Zeilen- und Spaltensummen von Delta verschwinden, somit haben beide Einzelereignisse genau die Randverteilung p. Weil Delta antisymmetrisch ist und v_a+v_b=v_b+v_a gilt, verändert es keine Zweischritt-Endwortwahrscheinlichkeit. Die möglichen Clifford-Phasen fallen in der unkontrollierten Konjugation heraus; sie werden hier nicht bei kohärenten Summen gleichgesetzt.

Wähle h einmal gleichverteilt aus den 60 geordneten Tripeln. Bedingt auf h werden unabhängige Ereignispaare mit J_h auf einer beidseitig unendlichen ganzzahligen Zeitachse erzeugt. Die beiden möglichen Paaranordnungen werden mit Wahrscheinlichkeit 1/2 gemischt.

Diese Mischung definiert einen konsistenten stationären Prozess: Eine Einerschiebung tauscht die beiden Paaranordnungen, eine Zweierschiebung verschiebt unabhängig gleichverteilte Paare. Die uniforme Menge der Tripel ist unter allen S5-Permutationen invariant. Die tatsächliche signierte Clifford-Wirkung ändert an den Zufallsunitär-Kanälen nur die Permutation, da skalare Vorzeichen aus der Konjugation herausfallen. Prozess B besitzt daher dieselbe markierte S5-Kovarianz wie A. Beide haben den stationären Systemzustand I4/4.

**Alle ungestörten Intervalle stimmen überein.** Jedes zusammenhängende Intervall zerfällt in ganze Paare und höchstens zwei einzelne Randereignisse aus verschiedenen Paaren. Bedingt auf h sind diese Stücke unabhängig. Ganze Paare liefern genau Phi^2, einzelne Randereignisse Phi; die Pauli-Konjugationen kommutieren. Folglich ist die Gesamtwirkung für jedes n und jede unkorrelierte eingesetzte Anfangsdichte exakt Phi^n, unabhängig von h und der Paaranordnung. Gemeint sind unbedingte Intervallkanäle bei frischem Systemeingang, keine bedingten Kanäle nach Auslesen der Vergangenheit.

Trotzdem unterscheiden sich die vollständigen Prozesse. Bei vier zeitlich getrennten Bell-Proben wird jeweils ein normierter maximal verschränkter C4-C4-Zustand eingesetzt. Die sechs Choi-Vektoren (ga tensor I)|Omega> sind orthonormal, daher lassen sich die sechs Fehlerlabels unterscheiden, ohne Zugriff auf eine vorausgesetzte physische Umwelt zu verlangen. Das sind mathematisch erlaubte Operationen des hier untersuchten M4-Prozesses; ihre physische TFPT-Realisierung wird nicht vorausgesetzt.

Für die beobachtete Labelsequenz (1,2,1,2) gilt

\[
\Pr_A=\frac1{12^4}=\frac1{20736},\qquad
\Pr_B=\frac1{20736}+\frac1{552960}=\frac{83}{1658880}.
\]

Herleitung: E Delta_12=0 und E Delta_12^2=3/10, weil genau 3 von 10 ungerichteten Dreiermengen die Labels 1 und 2 enthalten. Bei paarbündigem Beginn ist die Wahrscheinlichkeit E[J_h(1,2)^2]=(p1p2)^2+epsilon^2(3/10); beim anderen Beginn ist sie p1p2 E[J_h(2,1)]=(p1p2)^2. Mitteln der beiden Anordnungen ergibt die Differenz 3 epsilon^2/20=1/552960. Bis einschließlich drei Schritten kann höchstens ein vollständiges Paar auftreten, dessen Störung linear wegmittelt; dort stimmen sogar alle zeitaufgelösten Wortwahrscheinlichkeiten überein.

Dies beweist streng: Selbst **alle** Phi^n, ihre stationäre Dichte und die markierte Symmetrie wählen das Mehrzeitgesetz nicht eindeutig. Es ist kein Beweis gegen die Existenz eines quellenselektierten Gesetzes. Eine ausdrückliche Markov-Annahme würde Prozess B ausschließen, wäre aber zusätzliche Information, die ihrerseits aus der Quelle begründet werden muss.

### 4a. Der Unterschied bleibt auch bei ergodischer, mischender Dynamik bestehen

Das obige globale h ist im Allgemeinen nicht ergodisch. Um diese mögliche Zusatzforderung direkt zu prüfen, wird derselbe Gegenbeweis ohne konserviertes verborgenes h und ohne starre Paarperiode konstruiert. Dies ist eine Robustheitsprüfung des Schlusses, kein neuer physikalischer Ansatz.

Die Blocklänge beträgt unabhängig 1 oder 2 mit je Wahrscheinlichkeit 1/2. Ein Einerblock emittiert p, ein Zweierblock J_h. Zwischen Blöcken ändert sich h mit der Übergangsmatrix M=(I+U)/2, wobei U alle 60 Frames gleich gewichtet. Ein endlicher verborgener Markovprozess realisiert dies direkt:

- R_h bedeutet Beginn eines Blocks. Mit Wahrscheinlichkeit 1/2 emittiert er a~p als Einerblock und geht nach dem M-Übergang in R_h'; mit Wahrscheinlichkeit 1/2 emittiert er a~p als erste Hälfte und geht nach T_(h,a).
- T_(h,a) emittiert b mit Wahrscheinlichkeit J_h(a,b)/p_a und geht anschließend nach dem M-Übergang nach R_h'.

Die stationären Massen sind mu(R_h)=1/90 und mu(T_(h,a))=p_a/180. Sie summieren sich zu eins und bleiben unter der Übergangsmatrix erhalten. Alle Paarwahrscheinlichkeiten und alle Frame-Resets sind positiv; von jedem Zustand ist jeder andere erreichbar, und R_h besitzt eine positive Selbstschleife. Der endliche Prozess ist daher irreduzibel, aperiodisch und mischend. Seine stationäre Ausgabe ist entsprechend ergodisch und mischend. Die uniforme Frameverteilung und M bewahren die markierte S5-Kovarianz ensembleweit.

Der Intervallbeweis bleibt exakt: Volle Blöcke haben unabhängig von h den Endpunktkanal Phi oder Phi^2, abgeschnittene Randblöcke haben Einzelmarginale p. Bedingt auf die Frames und die Blockzerlegung faktorisieren ihre Ausgaben. Daher ergeben sämtliche ungestörten Intervalle weiterhin Phi^n für jede frisch eingesetzte Anfangsdichte. Bis drei zeitaufgelösten Nutzungen tritt höchstens ein vollständiger Zweierblock auf; dessen Störung mittelt weg.

Für die Viererfolge (1,2,1,2) kann nur die Anordnung zweier vollständiger Zweierblöcke einen quadratischen Unterschied beitragen. Die stationäre Wahrscheinlichkeit, am ersten Platz eines Zweierblocks zu beginnen, beträgt 1/3; der nächste Block hat mit Wahrscheinlichkeit 1/2 Länge zwei. Insgesamt ist die betreffende Anordnung daher 1/6 wahrscheinlich. Die Framekorrelation zwischen benachbarten Blöcken multipliziert E Delta_12^2=3/10 mit 1/2. Somit gilt jetzt

\[
\Pr_{B,\mathrm{mix}}-\Pr_A
=\tfrac16\cdot\tfrac12\cdot\tfrac3{10}\epsilon^2
=\frac1{3317760},\qquad
\Pr_{B,\mathrm{mix}}=\frac{161}{3317760}.
\]

Ergodizität oder Abklingen zeitlicher Korrelationen allein beseitigen die gezeigte Unterbestimmtheit folglich ebenfalls nicht. Exakte Markov-Unabhängigkeit oder ein ausreichend bestimmtes ursprüngliches Mehrzeitfunktional wären stärkere Angaben. Auch die ergodische Konstruktion ist kein zweiter vollständiger TFPT-Ursprung.

## 5. Anschluss an die ursprüngliche Quelle und Forschungsentscheidung

Der frühere Vorschlag, aus dem kompatiblen Nahtkanal plus gespeicherter Marke schon den gemeinsamen Ursprung zu gewinnen, war zu stark. Die Markierungs-Transportfrage lässt sich kohärent lösen; die vollständige Quellenwahl kann nicht durch mehr Potenzen desselben Transfers entschieden werden.

Die relevanten Originalquellen liefern positive, jedoch unterschiedlich vorausgesetzte Teile:

- `compiler-process-reconstruction-20260918`: Die tatsächlichen 60 E8-Operationsstrahlen besitzen eine vollständige komplexe Komposition. Ihre uniforme ungelesene Ausführung ist bereits als vollständig depolarisierend berechnet und trägt die nichttrivialen Nahtmoden nicht. Dieses gescheiterte Ereignisgesetz wird hier nicht wiederholt.
- `compiler-correlated-event-clock-20260919`: Ein ausgeführter und gespeicherter 60-Ereignisprozess erlaubt den korrekten mitgeführten Rahmen. Die Quelle benennt selbst, dass Messverfeinerung und klassisches Protokoll nicht alle ursprünglichen Interferenzen erhalten.
- `source-register-response-20260920`: Für den ausdrücklich angenommenen gemeinsamen QWZ-Fock-Lift existieren ein eindeutiger endlicher Gesamtgrundzustand (r=2,N=3,E=-4) und exakte geladene Mehrzeitantworten. Dies ist ein positiver Zustandsselektionsbeweis im angegebenen Modell, keine noch einmal offene numerische Frage. Dieser Lift ist bislang nicht als dieselbe primitive Quelle des markierten M4-Kanals nachgewiesen. Seine Zustandsauswahl darf daher nicht in den M4-Beweis eingesetzt werden.
- `moving-history-memory-round25`, `history-determinant-solver-round25` und `history-sum-control-round26`: Im ursprünglichen **deklarierten** g=0-Ladungs-/Skalarmodell liegt bereits eine vollständige regulierte euklidische Geschichtsdarstellung vor. Ihr positiver Transfer Q_delta, die exakte skalare Determinante, die ursprünglichen vorzeichenbehafteten E8-Hops und die absolut konvergente neutrale Summe sind vorhandene Ergebnisse. Sie sind keine neue Konstruktion dieser Sitzung. Die beiden Round25-Originalchecker wurden unverändert erneut ausgeführt.

Gerade dieser letzte positive Zusammenhang ist für die Fortsetzung wesentlich. Die Originalrechnung liefert mit G=M^-1

\[
\partial_i\partial_k\log\det M=-\delta^2G_{ik}G_{ki},\qquad
\partial_i\partial_k\log Z_\phi=+\tfrac12\delta^2G_{ik}G_{ki}.
\]

Der wiederholte Originaltest enthält ausdrücklich eine **nichtverschwindende zeitübergreifende Antwort im tatsächlichen Modell**. Ein Produkt isolierter Ladungsschrittgewichte würde diese bereits vorhandene Antwort entfernen. Der ursprüngliche skalare Zustand oder sein exakt erhaltenes Gedächtnis muss im gemeinsamen Objekt bleiben. Die unintegrierte größere Quelle kann durchaus eine Transferkomposition besitzen; ausgeschlossen ist hier das stillschweigende Weglassen ihres induzierten Gedächtnisses.

Der Atlas zeigte diese drei älteren Contracts zunächst nur mit PASS, generischem Scope und leerem `what_we_have`. Das ist keine Abwesenheit der mathematischen Ergebnisse. Die obigen positiven Aussagen wurden deshalb an ihren Originalbeweisen geprüft und als explizite, begrenzte Contract-Metadaten nachgetragen. Die absolute Konvergenz der Round26-Summe wurde hier am Beweis gelesen, nicht als vollständig ausgeführte Partitionssumme ausgegeben.

Auch dieses vorhandene Geschichtsgesetz hat benannte Eingaben: ein gewähltes räumliches Gitter, a,m,u,nu,J, einen euklidischen Zeitregulator und g=0. Weder ihre Auswahl aus P1/P2 noch die Identifikation dieses Quellraums mit dem markierten M4-Prozess ist damit bewiesen. Es wäre falsch, zu behaupten, im Atlas existiere überhaupt kein Mehrzeitfunktional; ebenso falsch wäre es, dieses vorhandene bedingte Funktional ohne die passende Feld- und Operationsabbildung zum vollständigen TFPT-Ursprung zu erklären.

Der konkrete Anschluss muss eine **gemeinsame geordnete Quellenwirkung** liefern, aus der der markierte Kanal, seine zulässigen Eingriffe und die Markierungsauslese gemeinsam folgen. Ein bloßer Endpunktvergleich reicht nach Abschnitt 4 nachweislich nicht. Für den Markov-Kandidaten lautet der erste zu prüfende geordnete Ausdruck Phi o Ad_A o Phi; der explizite gi/ gk-Test in Abschnitt 3 verhindert ein falsches Zusammenfassen von Geschichten. Der Vier-Schritt-Test prüft darüber hinaus versteckte zeitliche Korrelationen, die sogar alle Ein-, Zwei- und Drei-Schritt-Wortdaten dieser beiden Kandidaten überleben.

Dies ist die sachliche Grenze des hier verfolgten Ansatzes. Kein vollständiger T1–T8-Vertrag wird geschlossen. Der entscheidende Unterschied zu einer allgemeinen Lückenliste ist der ausgeführte Nachweis, welche Daten die gemeinsame Herleitung **nicht** bestimmen, sowie der konkrete geordnete Ausdruck, den eine tatsächliche Ursprungskonstruktion zusätzlich liefern muss.

## 6. Validierung, Herkunft und Literatur

`endpoint_compression_check.py` lädt die Original-Clifforddefinitionen, prüft die komplexen Produktphasen, exakte Faltung für n=1,...,8, Choi-Ränge und den Zwischenoperator. Der allgemeine Rang- und Transportbeweis steht oben; acht Stichproben ersetzen ihn nicht.

`history_identifiability_check.py` prüft mit rationaler Arithmetik alle 60 Paarverteilungen, ihre Positivität, Marginalen, identischen Endwortverteilungen, sämtliche Wortwahrscheinlichkeiten bis Länge drei und den exakten Vier-Schritt-Unterschied. Allgemeine Konsistenz, Stationarität und die Aussage für alle Intervalle folgen aus der beschriebenen beidseitigen Blockkonstruktion, nicht aus endlichen Zeitstichproben.

`instrument_selection_check.py` prüft den ursprünglichen markierten minimalen Lift und die verschiedenen Instrumentverfeinerungen. Die Endpunkt- und Instrumentrechnungen wurden unabhängig bearbeitet; die stationäre Mehrzeitkonstruktion zusätzlich analytisch gegengeprüft.

Diese Anwendung beansprucht keine neue allgemeine Quantentheorie. Die Unterscheidung zwischen reduzierten dynamischen Abbildungen und operationalen Mehrzeitprozessen ist in der Primärliteratur etabliert: [Pollock et al., Operational Markov condition for quantum processes, 2018](https://arxiv.org/abs/1801.09811). Der spezifische Beitrag hier sind die exakten Folgerungen und Gegenprozesse auf den tatsächlichen TFPT-Cliffordwörtern, ihre Markentransport-Gleichung und der berechnete Vier-Schritt-Zeuge.
