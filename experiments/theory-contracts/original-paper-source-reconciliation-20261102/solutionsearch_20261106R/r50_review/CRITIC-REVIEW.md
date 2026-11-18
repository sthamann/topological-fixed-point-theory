# Unabhängige Kritik R50

Der Ausschluss von HYP-R50-NONCOMPACT-JOINT-PERRON-LAW besteht die unabhängige Prüfung. Die Hypothese ersetzt tatsächlich die kompakte Basis aus R47 durch eine neue primitive Präparationsregel auf ganzem C². Ihr verlangter normalisierbarer unterster Eigenvektor existiert nicht. Das schließt diese partielle Ursprungshypothese aus, liefert aber kein ursprüngliches nu_Vn und keine vier ganzen Beiträge.

## Tatsächliche Ausführung

Der gelesene Proposer wurde als repeat_proposer.py unverändert in dieses Scratch kopiert und direkt mit dem gelieferten Python und -B ausgeführt. Exit0, RESULT_SHA256=b5e64ed1c63f7ea824d4b168fb1c15601ccca0df1fe497c8aa5342e9513c479b. repeat-comparison.json bestätigt bytegleiche Ergebnisse. Der deklarierte Test ist kein Beitragsproduzent: unnormalized_values ist null. Ein Vergleich von Z0,ZU,ZY,ZUY ist daher nicht möglich und wird nicht behauptet.

independent_check.py baut die 96 Matrizen aus den ursprünglichen Feldlabels, q(H), vier getrennten Yukawamatrizen und Majorana neu auf, ohne den Proposer zu importieren. Direkter Lauf mit demselben Python und -B: Exit0, RESULT_SHA256=bdbfe6c45b5cf23646c8b1950d5ac46bcaeb3c4ad48b6b836784a92725c14cee. 72 Matrixproben kontrollieren Gap und Energieobergrenze, auch mit verschiedenen komplexen invertierbaren Yukawamatrizen und einem nichtskalaren Majoranaanker. Das sind numerische Gegenproben; der allgemeine Satz folgt aus dem folgenden analytischen Argument.

## Angriff auf die Variationsschranke

Die ursprüngliche Karte ist q(H)diag(Yup,Ydown), q(H)*q(H)=|H|²I. Daher besitzt ihr Yukawaanteil kleinsten Singularwert |H|y_min und Norm höchstens |H|y_max. Der feste Majoranablock verschlechtert die untere Schranke höchstens um M. Auf dem ganzen vierdimensionalen Ball um H_R=(0,R) mit Radius L ist der Vorzeichengap mindestens g_R=(R-L)y_min-M>0. Es wird weder ein einzelner Strahl als Volumenmenge behandelt noch ein unbekannter Jacobian eingesetzt.

Die gamma-ungerade selbstadjungierte Phi führt zu h²=p²+Phi² und genau 48 negativen Eigenwerten. Für den vollen Λ48 Sektor gilt deshalb Esea+min Spec dGamma_48(h)=0. Das matrixwertige Potential W ist positiv. Der Heatanteil T=Tr exp(-h²) ist strikt positiv, und G+c_ref=-kappa Delta+W+T ist als Friedrichsform nichtnegativ.

Ein Eigenvektor an der Untergrenze müsste verschwindende gewöhnliche Ableitungsenergie besitzen. Er wäre ein konstanter Exteriorvektor und damit auf C² kein nichtnull L² Vektor. Dieses Argument allein genügte noch nicht: Ein Grundzustand könnte oberhalb der behaupteten Untergrenze liegen. Gerade diese Lücke schließt die neue vierdimensionale Variationsfolge.

Für den am Mittelpunkt eingefrorenen negativen Projektor Q_R und seinen Slatervektor ist F_R(H)=Esea(H)+Tr(Q_R h(H)) am Mittelpunkt samt erster Ableitung null. Die Vorzeichengruppen sind getrennt, auch bei Entartungen innerhalb der Gruppen. Ihre zweite Ableitungsformel ist 2 sum_(i<0,j>0)|B_ji|²/(lambda_j-lambda_i). Mit Tr(B²)<=96 y_max²|delta|² gilt die im Proposer verwendete, ausreichend großzügige Schranke F_R<=48 y_max²|delta|²/g_R. Der konstante Slatervektor erzeugt keine unterschlagene Berryableitung.

Unabhängige exakte Radialintegration ergibt integral_0^1 r³(1-r²)^4 dr=1/60. Mit |S³|=2pi² ist phi=sqrt(30)/pi (1-|x|²)_+² tatsächlich normiert; seine kinetische Energie ist 16 und sein zweites Moment 2/7. Damit gilt

0 <= <Psi_R,(G+c_ref)Psi_R>
 <= 16kappa/L²+(96/7)y_max²L²/g_R+96exp(-p²-g_R²).

Für L=R^(1/4) verschwindet die rechte Seite. Somit ist inf Spec G=-c_ref wirklich erreicht als Grenzwert von Rayleighquotienten, aber kein Eigenwert. Eine Energieverschiebung ändert diese Nichtexistenz nicht. Die normierten Testvektoren sind ein Variationswerkzeug, keine Normalisierung ursprünglicher Quellenfasern.

Die Kontrollreferenz einschließlich Majorana ist 90exp(-1/16)+6exp(-17/16)=86.62072016867464. Der wiederholte Test liefert die berichteten oberen Schranken 5.41978021978022 bei R=16 und 0.004520102257757547 bei R=16777216. Die analytische Abschätzung beweist die Grenze; die Stichproben und numerischen Heatnullen tun dies nicht.

## Quellenbezug, Neuheit und Grenze

Gelesene Originale: verification/v252_full_finite_triple.py:61-133; DYNAMICAL-HIGGS-SOURCE-PROOF.txt:14-55; check_dynamical_higgs_source.py:20-50; FIXED-CAP-DETERMINANT.txt:16-33,82-103; PROJECTED-JOINT-SOURCE-GENERATOR-20261103.txt:15-95; QFT004:137-155; Metrology005:136-178; EXACT-CUTOFF-HIGGS-MEASURE-PROOF.txt:1-138. Die genauen langen Pfade stehen im Proposerbundle und im Abschluss dieser Kritik.

Die interne Karte ist ein echter ursprünglicher Consumer. Sie als gegeben zu verwenden setzt das jetzt ausgeschlossene primitive Grundzustandsgesetz nicht voraus. premises_presupposing_target bleibt leer. Gewöhnliche Kinetik, Hilbertpaarung d4H, der gewählte gemeinsame Generator und die Forderung einer normalen untersten Präparation sind dagegen neue Ursprungshypothesen. Sie folgen nicht aus P1 und werden nicht mit QFT004s tatsächlichem Perronabschluss identifiziert.

R055 behandelt Integrabilität eines beschränkten Skalargewichts, nicht das Spektruminfimum eines matrixwertigen Schrödingeroperators. R37 schloss zusätzliche exakte Nullenergie bei einer vorgewählten Capbasis aus; R50 gibt diese Basis auf und bestimmt das wirkliche Infimum ohne vorgegebene Grundenergie. R40s Halbtraceargument enthält keine Bosonkinetik. Keiner dieser Sätze impliziert unter denselben Prämissen die gesamte jetzige Ausschlussaussage. already_implied_by bleibt leer. Die einschlägigen Atlasbefunde und die Journalreviews R37/R40/R41/R43/R49 wurden zum Scopevergleich herangezogen.

Der entscheidungsrelevante Fortschritt ist candidate_excluded, nicht input_derived oder candidate_computed. Die neue Regel für den fehlenden primitiven Slot scheitert an ihrer eigenen notwendigen Existenzbedingung. Alle offenen ganzen Faktoren bleiben offen: Eta, Orientierung, Phasen, Restfelder, tatsächlicher äußerer Selector, absolute Referenz, ursprüngliche Zeit und Marken. Ein anderer ursprünglicher Domainanschluss könnte die Testpakete entfernen und ist nicht ausgeschlossen. Eine einzelne Frequenz bei endlichem internem Rang ist keine ganze Kreisquelle und kein physischer 3+1D Abschluss. S/Z/OS/Gamma/state/time/reference, c3/P1/P2, E8, Clocks, Flavor/Yukawa/alpha und die bekannten exakten Teilintegrationen bleiben erhalten. T1..T8 werden nicht geschlossen.

## Nächste entscheidende Aufgabe

Keine weitere Grundzustandsoptimierung auf diesem Generator. Die vorhandene feste APS Randvektorregel aus R43 zu einem vollständigen endlichen Quellenproduzenten ergänzen, statt erneut ihre bereits geprüfte Capnormierung zu studieren: primitive gemeinsame Rang/IR Regel, Randvektor genau einmal, ganzer operatorwertiger Transfer, alle übrigen Skalar/Referenz/Eta/Orientierungs/Restfaktoren und äußerer Selector als Originaleingaben oder ausdrücklich neue Gesetze festlegen. Zuerst prüfen, ob die tatsächlich vorgeschlagene Randpräparation im gewählten endlichen Feldraum integrierbar ist; ein fester endlicher Cap allein dämpft die unbeschränkte Higgsrichtung nicht. Bei Verwendung der R43 gemeinsamen Rang/IR Regel deren bekannten Normhüllensatz wiederverwenden. Keine statische Bosonzeit aus R44 und kein skalarer Capstreifen aus R45 übernehmen. Danach die vier ganzen unnormierten Beiträge im selben ausführbaren Caller berechnen und den ersten neuen Skalar oder Zeitanschluss angreifen. Die ursprüngliche Auswahl der Randregel bleibt eine Hypothese; R50 leitet sie nicht her.

## Zugriff

Korrektur, Integrandaudit, Operatornotizen, beide aktuellen Claude Reviews und eingebetteter Gesamtkontext wurden gelesen. Die direkten globalen Skilldateien und der direkte globale Kontextpfad waren durch die Sandbox gesperrt; der mitgelieferte Gesamtkontext war zugänglich. Keine externen MCP, Agenten, Modelle oder Netze. Ausschließlich dieses Scratch wurde beschrieben. Keine neue globale Nichtfundaussage.
