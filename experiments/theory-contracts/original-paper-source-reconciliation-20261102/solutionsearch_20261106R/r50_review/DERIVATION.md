# R50: Quellenmaß aus gemeinsamem Grundzustand auf ganzem C²

## Konkrete Frage und erste neue Annahme

HYP-R50-NONCOMPACT-JOINT-PERRON-LAW ersetzt HYP-R47-COMPACT-C2-PRIMITIVE vollständig: Kein Radius7/4, keine gerundeten Zellen und kein vorgegebenes normiertes Kugelgesetz. Als primitive neue Regel wird die Konfigurationsmarginale eines normalisierbaren untersten Eigenvektors eines gemeinsamen Generators gewählt. Der primitive Raum ist L²(C²,d⁴H) tensor Λ⁴⁸(C⁹⁶). Die gewöhnliche Koordinatenkinetik und diese Grundzustandsregel sind neue Annahmen, keine Ableitung aus P1. Die Randpräparation und der Pushforward wären Psi_min und nu(A)=C_abs integral_A ||Psi_min(H)||² d⁴H. C_abs bleibt offen. Kein freies omega für QFT004 wird ergänzt.

Die konkrete Generatorhypothese lautet

G = -κ Δ_H + Esea(H) I + dGamma_48(h(H)) + T(H) - c_ref,
T(H)=Tr_96 exp(-h(H)²), Esea(H)=Tr_96 |h(H)|/2,
h(H)=p gamma_F-i gamma_F Phi(H), c_ref=T(0).

κ>0 ist offen, die Kontrolle benutzt κ=.3; p=1/4 ist eine deklarierte einzelne ursprüngliche Viertelfrequenz. Phi(H) enthält die ganze96 q(H) Karte und den festen Majoranablock. Die Rechenkontrolle verwendet Yf=I3, MR=I3 ausdrücklich als Kontrollbelegung. Die ursprünglichen physikalischen Flavorwerte werden nicht verändert. Eine einzelne Frequenz ist keine ganze Kreisquelle und keine 3+1D Physik.

Begründung dieser Hypothese: Sie versucht erstmals hier, die kompakte primitive Basis durch eine gemeinsame stationäre Präparation auf dem unverkleinerten C² zu ersetzen. Kinetik, ganze Exteriorquelle und der bisher benutzte Seeanker sind gemeinsam enthalten. Das wäre ein konkreter Maßproduzent, falls sein unterster Eigenvektor existierte. Es wird keine Nullenergie vorgeschrieben und keine Basisprojektivität verlangt. Der Test entscheidet genau die benötigte Existenz dieser Präparation.

## Analytischer entscheidender Test

Für die ganze96 Kontrollkarte gilt für alle H:

Phi*=Phi, {gamma_F,Phi}=0, h²=p²+Phi², Tr h=0.

Weil p ungleich0 ist, hat h stets48 negative und48 positive Eigenwerte. Auf dem vollständigen Λ48 Sektor ist daher

min Spec dGamma_48(h)=sum_(lambda<0) lambda=-Esea(H).

Somit W(H)=Esea(H)I+dGamma_48(h(H)) ist positiv und besitzt punktweise einen Nullvektor. Das ist eine ganze Slateridentität, kein Ersatz durch zwei oder vier Fermionmoden. T(H)>0 und folglich G>=-c_ref als Friedrichsoperator der geschlossenen semibeschränkten Form. Ein Eigenvektor zu -c_ref müsste insbesondere integral|grad Psi|²=0 erfüllen. Ein konstanter Vektor auf ganzem C² ist nur bei Psi=0 in L². Deshalb existiert am unteren Rand kein nichtnull Eigenvektor. Noch zu zeigen ist, dass dieser Rand tatsächlich das Spektruminfimum ist. Dies fehlt im alten Nullenergietest R37 und wird jetzt ohne adiabatische Näherung erledigt.

### Vierdimensionale Testpakete statt eines bloßen Strahls

Für Yf=I3 besitzt der Yukawateil L(H) exakt L(H)²=|H|² I96. Der erhaltene Majoranablock hat Operatornorm1. Auf jeder Kugel |H-H_R|<=L um H_R=(0,R) gilt daher mit g=R-L-1>0

min |Spec h(H)|>=sqrt(p²+g²).

Diese Schranke gilt auf einem ganzen vierdimensionalen Volumen, unabhängig von der Richtung des Störvektors. Es wird kein unbekannter Restgrad geraten.

Sei v_R der normierte vollständig besetzte negative Slatervektor von h(H_R), Q_R dessen Einteilchenprojektor. v_R wird auf dem ganzen Testpaket konstant gehalten. Es gibt daher keine übergangene Berryverbindung oder abgeleitete Fermionkinetik. Für δ=H-H_R ist

F_R(H)=Esea(H)+<v_R,dGamma_48(h(H))v_R>
      =Esea(H)+Tr(Q_R h(H)).

F_R(H_R)=0 und DF_R(H_R)=0. Die zweite Ableitungsformel der Summe der negativen Eigenwerte liefert für B=Dh[δ]

D² Esea[δ,δ]=2 sum_(i negativ,j positiv) |B_ji|²/(lambda_j-lambda_i)
            <= Tr(B²)/g =96 |δ|²/g.

Entartungen innerhalb derselben Vorzeichengruppe ändern diese Projektorformel nicht; die Gruppen sind durch den genannten Gap getrennt. Taylor entlang des geraden Segments im selben Testball ergibt exakt

0<=F_R(H)<=48 |δ|²/g.

Wähle phi(x)=sqrt(30)/pi * (1-|x|²)_+² auf R4. Sie ist normiert und in H1, und die im Produzenten symbolisch integrierten Konstanten sind

integral |grad phi|²=16,
integral |x|² |phi|²=2/7.

Mit Psi_R(H)=L^-2 phi((H-H_R)/L) v_R folgt deshalb

0 <= <Psi_R,(G+c_ref)Psi_R>
   <= 16κ/L² + (96/7)L²/(R-L-1)
      +96 exp[-p²-(R-L-1)²].                         (R50.1)

Setzt man L=R^(1/4), geht die ganze rechte Seite gegen0. Die normierten Pakete verlassen jede feste kompakte Menge. Ihre vierdimensionale Unterstützung ist Teil des Beweises. Damit ist inf Spec G=-c_ref und der Rand kein Eigenwert. Man braucht keine Behauptung, dass alle Punkte des übrigen Spektrums absolut stetig seien.

Für eine Folge von Dirichletkugeln mit gegen -c_ref minimierenden Grundenergien folgt sogar Verlust jeder lokalen Masse: Auf |H|<=B gilt T(H)>=96 exp[-(|p|+B+1)²]=m_B>0. Jede normierte minimierende Folge erfüllt integral_(|H|<=B)||Psi||² <=(<G>+c_ref)/m_B ->0. Ein bloßes Verschieben des Energieankers behebt diese fehlende L² Präparation nicht.

### Keine Gleichsetzung der vorgegebenen Yukawas nötig

Der analytische Ausschluss gilt allgemeiner für vier beliebige vorgegebene invertierbare Yukawamatrizen, ohne sie gleichzusetzen. Setze y_min=min_f sigma_min(Y_f)>0, y_max=max_f ||Y_f|| und M=||D_M||. Aus der ursprünglichen gemeinsamen q(H) Karte folgen sigma_min(L(H))=|H| y_min und ||Dh[δ]||<=y_max |δ|. Ersetze daher g durch g_R=(R-L)y_min-M und die rechte Seite von (R50.1) durch

16κ/L² + (96/7)y_max² L²/g_R +96 exp[-p²-g_R²].

Für L=R^(1/4) geht diese ebenfalls gegen0. Der untere Rand ist wieder -c_ref für die zu diesen tatsächlichen Matrizen gehörige Referenz. SingularY wird nicht mit ausgeschlossen. Die numerische Identitätsbelegung ist nur eine ausführbare Kontrolle des ganzen Constructors; der Satz benötigt keine Gleichsetzung verschiedener physischer Yukawas. Die Kontrollreferenz lautet exakt c_ref=90 exp(-1/16)+6 exp(-17/16), einschließlich des nichtnull Majoranaankers.

Die normierten Slatervektoren und Testpakete sind ausschließlich Vektoren eines Variationsbeweises. Sie normieren weder ursprüngliche Quellenfasern noch den ganzen Caller um.

## Beobachtete Rechnung

Direkter Lauf des gelieferten Python mit -B, Exit0. test_perron_source.py extrahiert ausschließlich die benannten Konstruktoren und Konstanten per AST aus verification/v252_full_finite_triple.py. Es importiert weder dessen Laufzeithooks noch lokale Modelle und schreibt nur hier.

RESULT_SHA256=b5e64ed1c63f7ea824d4b168fb1c15601ccca0df1fe497c8aa5342e9513c479b

c_ref=86.62072016867464, also Spektruminfimum -86.62072016867464 für die deklarierte Kontrollbelegung. Fünfzehn Matrixdatensätze prüfen die ganze96 Karte, einschließlich komplexer H Richtungen und Majorana. Die positiven Rayleighobergrenzen über dem Infimum aus (R50.1) fallen für R=16,256,4096,65536,1048576,16777216 von 5.41978021978022 auf 1.1742174160500853, .28975759376420007, .07233532857426309, .01808077864611043 und .004520102257757547. Die Formel beweist die Grenze; diese sechs Zahlen allein beweisen keinen Grenzwert. Numerisch unterlaufende Heatwerte werden als0 dargestellt; die strikt positiven exakten Exponentialterme und ihre Exponenten bleiben angegeben. Kein numerischer Nullwert wird zum exakten Nullbeweis.

Diese Zahlen sind Variationsschranken, keine vier unnormierten Quellenbeiträge. Da der verlangte primitive Grundvektor nicht existiert, entsteht aus dieser Regel kein ausführbarer ganzer Quellenproduzent. Z0,ZU,ZY,ZUY wurden daher nicht erfunden. Der Produzent ist ein ausführbarer Ausschlusstest einer partiellen Quellenhypothese, nicht candidate_computed.

## Angriff, Quellenfaktoren und enger Scope

Ausgeschlossen ist genau HYP-R50-NONCOMPACT-JOINT-PERRON-LAW mit gewöhnlicher C² Kinetik, der deklarierten endlichen ganzen96 Karte, bewegtem Seeanker, positivem Heatprofil und normaler Grundzustandspräparation vor einem ursprünglichen äußeren Selector. Noch offene ursprüngliche Faktoren werden nicht still als eins gesetzt. Der Satz gilt für den benannten Teilproduzenten und erklärt, warum er den primitiven Input nicht liefern kann. Er behauptet nicht, dass der ursprüngliche Gesamtcaller nur aus diesem Generator besteht.

Die hypothetische partielle Quellenbeschreibung candidate-partial.json erhält offene Slots für tatsächliche Phasen, Eta, Orientierung, Restfelder, absolute Referenz, ursprünglichen Selector und Zeitmarken. Solche Faktoren können einen anderen Generator, andere Domain oder einen anderen Raum auswählen. Eine tatsächliche physische Einschränkung könnte die Testpakete entfernen; das wurde nicht vorausgesetzt oder ausgeschlossen. Es wird keine innere P Einfügung vorgenommen und kein ursprünglicher äußerer Quotient geändert. QFT004s Perronabschluss des tatsächlichen Transfers bleibt unberührt.

Gegenüber R055 liegt ein anderer mathematischer Gegenstand vor: R055 betrifft Integrabilität eines skalaren beschränkten Gewichts auf flachem Volumen. Ein Schrödingeroperator mit beschränktem Potential kann durchaus einen gebundenen Grundzustand besitzen. Hier werden zusätzlich das nichtnegative ganze matrixwertige Potential, dessen exakte Seekomplementaufhebung und eine vierdimensionale Variationsfolge bewiesen. Daher folgt der gesamte R50 Satz nicht aus R055. Gegenüber R37 wird die Capbasis aufgegeben, keine Energie0 gefordert und das tatsächliche Spektruminfimum samt Nichterreichen ermittelt. R40s HS Trace Divergenz enthält keine Kinetik und impliziert diesen Satz ebenfalls nicht. already_implied_by bleibt für den gesamten neuen Satz leer.

Der erste zusätzliche Input bleibt ein Ursprungsgesetz für nu_Vn. R50 setzt es nicht als bereits gegeben voraus, sondern versucht es durch eine konkrete neue gemeinsame Präparationsregel zu erzeugen und widerlegt genau diese Regel. Es wird keine globale Quellenabwesenheit behauptet und keine neue Nichtfundsuche als Ergebnis gezählt.

## Gelesene Originalstellen

* verification/v252_full_finite_triple.py:61-133: Index, ganze96 Karte, Majorana, gamma_F und J.
* experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/DYNAMICAL-HIGGS-SOURCE-PROOF.txt:14-30,34-55,129-143: q(H), tatsächlicher gemeinsamer Consumer, ursprüngliche Integration und KO Einschränkung.
* experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/check_dynamical_higgs_source.py:20-50: tatsächlicher voller q(H) Matrixconstructor.
* experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/FIXED-CAP-DETERMINANT.txt:16-33,82-103: h, ganze besetzte See und Unterscheidung von Cap, Bulk und ursprünglichem Skalar.
* experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/EXACT-CUTOFF-HIGGS-MEASURE-PROOF.txt:19-100: verwendeter, nicht erneut behaupteter R055 Scope.
* experiments/theory-contracts/original-local-source-variable-background-20261102/PROJECTED-JOINT-SOURCE-GENERATOR-20261103.txt:15-35,43-95,98-125: ganzer Caller, tatsächlicher Integrator, ursprüngliche Zeit und ausschließlich äußerer Selector.
* documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex:136-178: tatsächliches nu_Vn Symbol und bedingte Koerzivität.
* research-000049/critic-001/scratch/CRITIC-REVIEW.md:49-57: keine Perronprojektion auf unverändert vorgewählter Kugel als Herkunftsherleitung.

Zugriff: Die direkten globalen Skilldateien waren durch die Sandbox gesperrt. Der gelieferte vollständige globale Kontext wurde lokal aus full-knowledge-context.json ausgelesen und gelesen. atlas-selection.json protokolliert die einschlägigen Befunde; R37, R40 und R41 wurden aus dem vollständigen Journal verglichen. Korrektur, Integrandaudit, beide aktuellen Claude Reviews und Operatornotizen wurden zuerst gelesen. Keine externen MCP, Netze, Modelle oder Agenten. Kein unabhängiger Kritiker wird hier behauptet.

## Nächster entscheidender Schritt

Den bereits ursprünglichen festen APS Rand statt der nicht vorhandenen bewegten Grundzustandspräparation als Quelle verwenden: die ganze feste Capamplitude mit Majorana und gebundener Referenz aus R41/R43 als tatsächliche Randbedingung eines gemeinsamen Bosontransfers einsetzen. Zuerst den durch diesen Rand erzeugten unnormierten Konfigurationspushforward berechnen und prüfen, ob die fehlende zusätzliche Wirkung oder das Randgewicht die vorstehende Fluchtrichtung wirklich unterdrückt, ohne den Cap nochmals als unabhängigen Bulkfaktor zu zählen. Die Wahl dieser Randregel bleibt eine benannte Ursprungshypothese, bis sie aus dem ursprünglichen Caller folgt. Keine erneute reine Capmomentenstudie und keine Kugelprojektion. Bei tragfähigem gemeinsamen Produzenten alle vier Beiträge mit sämtlichen explizit gesetzten oder original ausgewerteten Faktoren ausrechnen.

Primitive Quelle, ursprünglicher Regulator/Domainanschluss, physische Zeitmarken, absolute Masse und T1..T8 bleiben offen. S/Z/OS/Gamma/state/time/reference, c3/P1/P2, E8, verschränkte Clocks, Flavor/Yukawa/alpha und bekannte Teilintegrationen bleiben in ihrem exakten Scope erhalten.
