# Audit des eingefügten Prozess-/Quellentextes

**Gegenstand.** Geprüft wurden nur Abschnitt 2 (behaupteter Intertwiner der vier Fünferräume) und Abschnitt 5 (nativer Kanal, Lüders-Folge, Lindblad-Ausschluss, Paar-Hamiltonian). Der eingefügte Text ist Sekundärmaterial, keine Quelle. Repo unverändert.

## Kurzurteil

1. **Der Fünfer-Intertwiner ist in der angegebenen Form unmöglich.** Eine isometrische Abbildung zwischen zwei fünfdimensionalen Räumen ist unitär/surjektiv. Aus `T C_v = P_v T` folgt daher unitäre Äquivalenz von `C_v` und `P_v`, insbesondere identische Eigenwertmultiplizitäten. `4+1` und `3+2` widersprechen dem. Die konkrete Spurprobe lautet `tr(C_v)=3`, `tr(P_v)=1` (bis zu gemeinsamem Vorzeichen wären es `-3` und `-1`); Spuren können unter Intertwining nicht wechseln. Der Text kann höchstens zwei verschiedene natürliche Wirkungen beschreiben, aber keine Äquivalenz unter genau diesem Generator.
2. **Die Lüders-Matrix ist exakt, aber die Beschreibung muss den Kanal zwischen den Messungen nennen und `J_3` definieren.** Mit `J_3=11^T/3`, dem nativen Kanal `Phi` und den originalen Effekten gilt tatsächlich die autonome konditionierte Übergangsmatrix `M=(B+J_3)/2`. Ohne den dazwischenliegenden `Phi`-Schritt ist die Lüders-Folge eine andere Matrix, `3G`, mit Spektrum `{1,1/2,1/2}`. Das physische Instrument ist weiterhin nicht aus P1/P2 ausgewählt.
3. **Der Rangtest liefert tatsächlich einen allgemeinen endlichen homogenen Lindblad-Ausschluss für den Viererkanal.** Das ist stärker als der alte Vertrag, der nur die markierte orbit-skalare Halbgruppenklasse geprüft und entsprechend eng formuliert hatte. Für einen endlichdimensionalen zeit-homogenen QMS ist der Choi-Träger für jedes `t>0` von gleicher Dimension. Hier springt er von `rank J(Phi)=6` auf `rank J(Phi^2)=16`; daher kann `Phi` auf demselben `M_4` nicht `exp(tau L)` für einen GKLS-Generator sein. Der Text darf dieses neue Argument nicht als bereits im alten Contract bewiesen ausgeben. Die zusätzlich genannten Ränge für Fünfer- und geladenen Zwanzigerprozess wurden in diesem gezielten Audit nicht an Originalquellen bestätigt.
4. **Die thermische Choi-Präparation stimmt nur nach expliziter Einschränkung auf `W_6=P_1\oplus P_5`.** Im vollen 16-dimensionalen Raum ist der Gibbszustand bei endlichem `beta` vollrangig und enthält zwingend den Zehnersektor. Der Originalvertrag arbeitet ausdrücklich im `1+5`-Bad und nennt den zusätzlichen `4g`-Zehnersektor des vollen Raums. Die Einschränkung auf `W_6` ist daher eine zusätzliche Quell-/Sektorauswahl, keine Folge des vollen Paar-Hamiltonians allein.

## 1. Fünferräume und der behauptete Intertwiner

Für den im Text selbst angegebenen Quartikcode ist die `4+1`-Zerlegung direkt sichtbar. Unter `D=diag(-1,1,1,1)` auf jedem der vier Register sind `q_0` und die drei `aabb/ccdd`-Zustände gerade; `q_4`, das jedes Symbol genau einmal enthält, ist ungerade. Auf `H^0(P^1,O(4))` hat `(u,v)->(-u,v)` auf der Monombasis `u^k v^(4-k)` dagegen drei gerade (`k=0,2,4`) und zwei ungerade (`k=1,3`) Richtungen.

Angenommen, `T` sei die behauptete Isometrie zwischen beiden fünfdimensionalen Räumen und `T C_v=P_v T`. Dann ist `T` invertierbar und `P_v=T C_v T^{-1}`. Das verlangt dasselbe charakteristische Polynom. Tatsächlich wären es

- Quartikcode: `(x-1)^4(x+1)`,
- Riemann–Roch-Wirkung: `(x-1)^3(x+1)^2`.

Damit ist der Satz als geschrieben widerlegt. Man könnte einen Intertwiner erst nach Änderung des Operators, der Darstellung oder des Zielraums suchen; keine dieser Änderungen ist im Text angegeben.

Der Atlas bestätigt nur die schwächere Vorgeschichte: `claim:ARCH.RRCAR.01` führt `h^0=5` und `rank H_1=3` als exakte Riemann–Roch-Arithmetik, die physische Identifikation des Carriers mit `H^0(P^1,O(mu4))` aber nur als `[C]`. Es gibt dort keinen Beleg für den neuen vierfachen Operator-Intertwiner. Registrierte Quellen sind unter anderem `source:301753fcdb13e7f57ae9` (`verification/v189_riemann_roch_carrier.py`, Zeilen 1–86) und `source:7a302f0bab969f351261` (Ledgerzeile 262). Der direkte Atlas-Aufruf liefert `source:301...` vollständig. Zwei frühere Leerantworten kamen aus dem lokalen Helper, der JSON-RPC-Fehlerfelder verwirft; sie sind kein belegter Server-Retrievalfehler. Die spektrale Widerlegung ist davon unabhängig.

**Folge:** Auch die Wiederholung dieses `T C_v=P_vT` als T1-Schritt in Abschnitt 8 des eingefügten Textes ist nicht tragfähig.

## 2. Lüders-Folge: was `M` genau bedeutet

Der Originalvertrag `contract:marked-seam-cp-lift-20260921` belegt den bedingten Kanal

`Phi = (7/12) id + (1/12) sum_a Ad(Gamma_a)`

und die Heisenberg-Wirkung `Phi(E_i)=sum_j B_ij E_j`; siehe `source:4fc6776f318cd172d39d`, Zeilen 3–20, sowie `source:f279d6fa4d14536a2311`, Zeilen 84–100. Er sagt zugleich ausdrücklich, dass noch kein physisches Messinstrument ausgewählt ist (`source:4fc...`, Zeile 20; `source:f279...`, Zeilen 139–144).

Die im eingefügten Text undefinierte Größe muss sein

`J_3 = (1/3) 1 1^T`,

also der Rang-eins-Projektor auf die uniforme Dreierichtung; genau diese Definition steht im Checker (`source:f279...`, Zeile 96). Mit `E_i^2=(2/3)E_i` gilt `sqrt(E_i)=sqrt(3/2)E_i`. Der entscheidende exakte Identitätstest ergab für alle `i,j`

`sqrt(E_i) Phi(E_j) sqrt(E_i) = M_ij E_i`,

mit

`M=(B+J_3)/2 = (1/36)[[19,7,10],[7,19,10],[10,10,16]]`.

Deshalb ist die beobachtete Folge wirklich Markovsch: Nach einem gespeicherten Ausgang `i`, anschließendem `Phi`-Schritt und nächster Lüders-Messung ist die Verteilung des nächsten Ausgangs die `i`-te Zeile von `M`, unabhängig vom verborgenen Zustand innerhalb des Lüders-Unterraums. Ein exakter Scratch-Check bestätigte die Operatoridentität und auch die dritte Messung für alle Vorgeschichten.

Ohne `Phi` zwischen zwei Messungen ergibt sich dagegen aus dem originalen Gramoperator `G_ij=Tr(E_iE_j)/4`

`3G = [[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]]`,

nicht `M`. Die Formulierung „wiederholtes Messen mit Lüders-Zweigen liefert M“ ist daher unvollständig; korrekt ist „wiederholte Folge *Lüders-Messung – nativer Phi-Schritt – Lüders-Messung*“. Die Schranke `lambda_max(E_i)=2/3 < B_ii=13/18` zeigt korrekt, dass `B` selbst keine solche konditionierte Born-Übergangsmatrix der `E_i` sein kann.

Die Algebra bestimmt damit eine kompatible konditionierte Instrumentfolge. Sie leitet weiterhin weder deren physische Auswahl noch ein klassisches Pointerregister aus P1/P2 her.

## 3. Choi-Rang und allgemeiner homogener Lindblad-Ausschluss

Der alte Contract beweist nur: Der in derselben markierten orbit-skalaren Klasse erzwungene Halbzeitschritt hätte einen negativen Choiwert (`source:f279...`, Zeilen 75–78), und beschränkt den Schluss ausdrücklich auf diese Klasse (`source:4fc...`, Zeile 24; `source:f279...`, Zeile 144).

Der neue Ranggedanke ist dennoch mathematisch stärker und korrekt. Für einen endlichdimensionalen normstetigen QMS `T_t=exp(tL)` wende die Standard-Trajektorien-/Supportdarstellung auf den maximimal verschränkten Vektor im System plus untätiger Ancilla an. Der Choioperator ist genau der daraus entstehende Zustand. Sein Träger hat für jedes `t>0` die Form `P_t S(Omega)`, wobei `P_t` invertierbar ist und der durch Lindblad-Wörter und iterierte Kommutatoren erzeugte Raum `S(Omega)` nicht von `t` abhängt. Also ist `rank J(T_t)` für alle positiven Zeiten konstant. Eine passende Primärquelle ist Fagnola–Girotti, Lemma C.2, *Irreducibility of Quantum Markov Semigroups, uniqueness of invariant states and related properties* (2026), deren analytischer Trajektorienbeweis genau diese zeitunabhängige Supportspanne liefert; die Choi-Fassung folgt durch die untätige Ancilla.

Beim vorhandenen Viererkanal sind die sechs Krauswörter `I,Gamma_1,...,Gamma_5` Hilbert–Schmidt-orthogonal und haben positive Gewichte, also `rank J(Phi)=6`. Nach zwei Schritten sind alle 16 Clifford-Endwörter positiv belegt:

- Identität: `3/8`,
- jedes der fünf `Gamma_a`: `7/72`,
- jedes der zehn Bivektoren: `1/72`.

Daher `rank J(Phi^2)=16`. Wäre `Phi=T_tau`, wäre `Phi^2=T_{2tau}`; die konstante positive-Zeit-Rangregel wird verletzt. Das schließt **jeden** zeit-homogenen GKLS-Generator auf demselben `M_4` aus.

Der Ausschluss betrifft nicht:

- zeitabhängige oder nicht-Markovsche Generatoren,
- größere System-/Umwelträume und anschließendes Ausspuren,
- diskrete Collision-Modelle mit frischen Zellen,
- konditionierte Instrumente/Records,
- TFPT insgesamt.

Für die behaupteten weiteren Rangfolgen `15->25` und `20->150->400` fehlen in diesem gezielten Quellenpaket registrierte Originalstellen. Sie sollten separat mit den jeweiligen tatsächlichen Kanälen belegt werden; der allgemeine Satz allein erzeugt diese Zahlen nicht.

## 4. Paar-Hamiltonian und Gibbszustand

Aus `A=sum_a Gamma_a tensor Gamma_a` mit Eigenwerten `5,-3,1` auf den Rängen `1,5,10` folgt für `H_pair=g(5I-A)`:

- `0` auf `P_1`,
- `8g` auf `P_5`,
- `4g` auf `P_10`.

Bei `beta g=(log 7)/8` ist der **volle** Gibbszustand

`rho_beta = [P_1 + 7^(-1) P_5 + 7^(-1/2) P_10] / [1+5/7+10/sqrt(7)]`.

Er hat Rang 16 und ist nicht `tau_Phi=(7P_1+P_5)/12`. Das folgt schon allgemein: Der Gibbszustand eines endlichen Hamiltonoperators ist bei endlicher Temperatur strikt positiv und kann keinen Zehnersektor exakt auslöschen.

Auf dem ausdrücklich komprimierten Raum `W_6=P_1+P_5` ist die Aussage dagegen exakt:

`rho_(beta,W6) = [P_1+(1/7)P_5]/[1+5/7] = (7P_1+P_5)/12`.

Dies entspricht dem Originalvertrag. `source:a045cc82fb62819a74e7`, Zeile 31, sagt ausdrücklich: im `1+5`-Raum liegen die Energien bei `0,8g`, der volle `1+5+10`-Raum enthält zusätzlich `4g`. Zeilen 32–36 behandeln ausschließlich dieses `1+5`-Bad und nennen es eine Approximationsvorschrift, keine mikroskopische Kühlungsherleitung. Zeilen 38–53 geben die positive bedingte Paar-an-Pulskonstruktion an, lassen aber reine Grundpräparation, `kappa/g`, Pulszeit und vollständige Quellenwahl als Eingaben offen.

Die korrekte Formulierung lautet daher: **„Der auf `W_6` restringierte Gibbszustand von `H_pair|W6` ist bei `beta g=log(7)/8` gleich dem Choi-Zustand.“** Nicht korrekt ist eine thermische Präparation dieses rangsechs Choi-Zustands durch den vollen 16-dimensionalen Hamiltonoperator ohne zusätzliche Sektorauswahl.

## Quellen- und Atlasbefund

- `contract:marked-seam-cp-lift-20260921` — `PARTIAL`, bedingter diskreter CP-/Readout-Anschluss; Originalquellen `source:4fc6776f318cd172d39d`, `source:f279d6fa4d14536a2311`.
- `contract:joint-singlet-pulse-source-20261013` — `PARTIAL`; Originalbeweis `source:a045cc82fb62819a74e7`, Checker `source:b5bdec945d5196948e81`.
- `claim:ARCH.RRCAR.01` — exakte `h^0=5`, `rank H_1=3`, aber Carrieridentifikation `[C]`; `source:301753fcdb13e7f57ae9` lässt sich direkt vollständig lesen.
- Parserbefund: `get_dossier(expanded=true,max_chars=14000)` wird korrekt mit JSON-RPC `-32602` abgewiesen, weil das vollständige Dossier 38.037 Zeichen benötigt. Der lokale Helper liest nur `result.content` und schreibt bei einem top-level `error` fälschlich eine leere Ausgabe. Der kompakte direkte Dossieraufruf funktioniert. Die früheren Leerantworten sind daher nicht als Atlas-Serverfehler zu führen.

**Gesamtbewertung des eingefügten Textes:** Abschnitt 2 enthält mit dem `4+1`/`3+2`-Intertwiner einen direkten algebraischen Widerspruch. Abschnitt 5 mischt drei verschiedene Evidenzstufen: eine exakt rekonstruierbare konditionierte Lüders-Folge, einen neuen und tatsächlich allgemeinen Choi-Rang-Ausschluss, sowie eine nur im gewählten `W_6` richtige thermische Präparation. Die Überschrift „Vollständige Lösungen“ und die spätere Verwendung als gemeinsamer fundamentaler Elternprozess sind durch die Originalverträge nicht gedeckt.


## Ergänzung zur Quellenpriorität, 14. Oktober 2026

Der Vergleich mit dem engeren Vertrag vom 21. September war korrekt, die Einordnung als neue Verschärfung im Projekt war unvollständig: `source-continuous-marked-time-search-20260922` beweist bereits den Ausschluss sämtlicher zeitlich homogener GKLS-Generatoren auf demselben M4 für alle positiven ganzzahligen Potenzen des angegebenen Kanals. Der hier geprüfte Choi-Rangsprung ist ein alternativer Beweis für Phi, keine erstmalige Schließung dieser Projektfrage. Sein mathematischer Geltungsbereich bleibt unverändert. Größere Umgebungen, GNS-unitäre Zeit und zeitabhängige Generatoren sind dadurch nicht ausgeschlossen.
