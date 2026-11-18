# Gewichte: tatsächliche vorhandene Regeln und neue exakte Auswertung

5. November 2026. Read-only kanonische Arbeit. Ausgang ist der tatsächliche Source45 QFT004 Aufruf, nicht ein neues Modell. Gegenstand sind drei verschiedene Objekte: trennende diagonale Einfügungen, dynamische Transfers und die daraus abgeleiteten Zustands beziehungsweise Spektralgewichte.

## 1. Trennende Diagonaloperatoren verlangen keine besonderen Amplitudenzahlen

Source45 QFT004:282–335 benötigt eine Diagonalalgebra, die die retained Basiswörter trennt. Die anschließende Burnside Argumentation benutzt ihre gemeinsamen Eigenräume und einen tatsächlich strikt positiven P_V. Sie benutzt keinen vorgeschriebenen Zahlenwert d_r(w). Dieses Symbol steht nicht im dortigen Original.

Auf dem endlichen tatsächlich gegebenen retained Wortbereich W seien L_1,...,L_k die dortigen kommutierenden Count und Labeloperatoren. Sei l_j(w) der wirkliche Wert des jten Labels. Für seine endliche Wertemenge A_j definiere

`P_(j,a) = product_(b in A_j, b!=a) (L_j-b I)/(a-b)`.

Dann gilt durch Einsetzen auf jedem Basiswort exakt `P_(j,a)|w> = 1_(l_j(w)=a)|w>`. Also ist

`E_w = product_j P_(j,l_j(w))`

der Projektor auf alle Wörter mit demselben vollständigen Labeltuple wie w. Wenn die ursprünglichen Tuple die Wörter tatsächlich trennen, ist E_w genau der Rang1 Projektor `|w><w|`. Daraus folgt ohne freie Gewichtwahl die ganze Diagonalalgebra. Umgekehrt genügt jede bereits trennende Familie D_V für dieselbe Spektralprojektorkonstruktion; ihre projektorwertige Umkodierung erhält ihre erzeugte Diagonalalgebra.

Diese Konstruktion beweist nicht zusätzlich die Trennung eines unvollständig spezifizierten Worttuples. Gleiche Counts bei verschieden geordneten Wörtern ergeben einen gemeinsamen Blockprojektor. Der ursprüngliche retained Wort und Labelvertrag ist daher beizubehalten. Die Konstruktion bestimmt auch nicht den dynamischen T_V: Projektoren anstelle tatsächlicher Transfer oder Krausamplituden einzusetzen wäre ein neuer Quellenaufruf.

Folge für die vorherige Formulierung: Die nicht vorhandene konkrete numerische d_r Schreibweise darf für den Burnside Schritt nicht als zusätzlich fehlende physische Dynamikpflicht ausgegeben werden. Tatsächlich gebraucht wird dort die ursprüngliche Trennung; für Quellenwerte braucht man zusätzlich den tatsächlich gegebenen Transfer.

## 2. Tatsächliche v330 Gewichte für sämtliche diagonalen Einfügungen

Der unveränderte Originalcaller `verification/v330_qgamb_admissible_measure.py:symmetric_transfer:42–50` liefert

`T = P0 + r P2 + s P3`, mit `r=64/729`, `s=1/729`,

`u0=(1,1,1)/sqrt(3)`, `u2=(1,-1,0)/sqrt(2)`, `u3=(1,1,-2)/sqrt(6)` und `Pj=uj uj^T`.

Die tatsächliche ganze Matrix ist

```
T = [[1651/4374, 1267/4374, 728/2187],
     [1267/4374, 1651/4374, 728/2187],
     [ 728/2187,  728/2187, 731/2187]].
```

Für jedes ursprüngliche diagonale Testobservable `D=diag(a,b,c)` ergibt Spektralzerlegung exakt

`<u0,D T^n D u0> = A + B r^n + C s^n`,

`A=(a+b+c)^2/9`, `B=(a-b)^2/6`, `C=(a+b-2c)^2/18`.

Beweis: `A_j=<u0,D Pj D u0>=|<uj,D u0>|²`. Einsetzen der drei ausdrücklich gelieferten uj ergibt die angezeigten Formeln. Der verbundene Correlator ist `B r^n+C s^n`; das konstante A ist gerade `<D>²`. Dies sind abgeleitete Spektralgewichte der tatsächlichen Quelle und Einfügung, keine freien Modellparameter.

Die stationären Koordinatengewichte sind genau `pi_i=|u0_i|²=1/3`. Da T positive Einträge hat und jede Zeilensumme1 ist, besitzt jede tatsächlich eingefügte Folge der ursprünglichen Koordinatenprojektoren die ganze geordnete Antwort

`<u0,E_i0 T^n1 E_i1 ... T^nk E_ik u0> = (1/3) product_j (T^nj)[i_(j-1),i_j]`.

Alle Projektoren und Transfers bleiben auf demselben ursprünglichen Dreizustandsraum. Dies ist eine vollständige explizite Gewichtregel für diesen eigenen Caller, keine Herleitung seines Anschlusses an die primitive S_adm Quelle oder an ein 3+1D Materie/Geometriefeld.

Ausgeführt: `derive_v330_weights.py`, **15/15 PASS**. Allgemeine symbolische Gewichtsidentität, Orthogonalprojektoren, ganze ursprüngliche Matrix, Perronvektor und Raten, symbolische Momente n=0..7 sowie tatsächlicher unveränderter v330 correlator und ganze zweistufige Labelwörter geprüft. Allgemeine n und Wortzahl folgen aus der obigen Spektralzerlegung beziehungsweise Projektorproduktformel, nicht aus der endlichen Anzahl Tests.

## 3. Was v171 und v175 tatsächlich schon liefern

`v171_os_moment_cluster.py:55–64` setzt A=B=C=1 als exemplarische beliebige nichtnegative Gewichte für den allgemeinen Hankel/OS Momentensatz. Das ist kein Producer von dynamischen Amplituden. Sobald eine wirkliche Quelle samt Observable vorliegt, sind A,B,C wie in Abschnitt2 berechenbar. Ausgeführt: v171 **6 bestanden,0 fehlgeschlagen**.

`v175_net_existence_full_cone.py:94–121` nimmt ausdrücklich eine Einteilchenkontraktion mit acht Eigenwerten1, den beiden Raten r,s und sechs weiteren Eigenwerten `r*0.1^k`. Der CAR Funktor produziert daraus für jede Besetzung S das tatsächliche Exteriorgewicht `product_(j in S) t_j`. Die ganze Quelle ist im gegebenen Kontraktionszweig bestimmt. Die Prüfroutine leitet nicht die sechs eingegebenen weiteren Einteilchenwerte aus dem primitiven Quellenaufruf her. Ausgeführt: v175 **5 bestanden,0 fehlgeschlagen**.

`v330_qgamb_admissible_measure.py` zusätzlich unverändert ausgeführt: **6 bestanden,0 fehlgeschlagen**, eigener ausgewiesener Umfang `admissible-sector OS measure built; metric sector open`. Bestehende Resultate werden dadurch bestätigt, kein wissenschaftlicher Geltungsbereich erweitert.

## 4. KMS: eine tatsächliche Gewichtsregel bei geliefertem gemeinsamen Generator

`v258_dirac_covariance_induction.py:65–76` liefert die tatsächlichen inversen Operationen `H=log((1-C)C^-1)` und `C=(1+exp(H/mu))^-1`. Für eine tatsächlich quasifreie finite CAR Quelle sind die Besetzungsgewichte im Eigenmodusrahmen daher `product_j c_j^n_j (1-c_j)^(1-n_j)`. `v813_p1_index_kms.py:264–273`, gaussian_rho, setzt genau dieses Produkt als ganze Dichtematrix um.

Bei der vollständigen endlichen Operatoralgebra und tatsächlich gegebenem Zeitgenerator H fixiert KMS die Gibbsgewichte: in einer Eigenbasis folgen aus den Matrixeinheiten die Verhältnisse `rho_j/rho_i=exp[-beta(E_j-E_i)]`; die Normierung liefert `rho=exp(-beta H)/Tr exp(-beta H)`. Der vorhandene physische reduzierende Selektor bleibt im tatsächlichen Trace erhalten. Auf einer neutralen Teilalgebra ohne verbindende Matrixeinheiten ist diese volle-Faktor Folgerung nicht automatisch zulässig.

Diese Rechnung ist kein zusätzlicher Quellenselektor. Wird H zuerst aus derselben noch unausgewerteten omega definiert, produziert die Rückformel keinen neuen Quellenwert. Bei dynamischen Bosonen und ganzen nichtfaserweisen Selektoren ist auch das Produkt der bedingten Fermiongewichte nicht schon der ganze gemeinsame Feldkernel.

`v259_modular_cutoff_kappa.py:54–63` berechnet für sein deklariertes KMS Profil `f(u)=exp(-u)` die drei Profilmomente f0=f2=f4=1. Das sind die Gewichte der analytischen Spektralprofilrealisierung. Ihre Identifikation mit einer Konfigurationsdichte der gesamten wechselwirkenden Quelle folgt nicht aus den drei Momenten.

## 5. Quillen und Determinanten: tatsächliche Faktoren mit tatsächlicher Einfügungsart

`ELLIPTIC-QUILLEN-RESPONSE.txt:29–89` im vorhandenen original-admissible-source-closure-20261025 Contract wertet die tatsächliche freie elliptische Verbindung aus: `N(z)=exp(-2pi u²)|theta3(z|i)/theta3(0|i)|²`. Domänen, Nullstelle und gebundener Referenzpunkt bleiben enthalten. Die Einfügung ist eine flache U1 Verbindung im angegebenen elliptischen Quellenzweig. Große Eichtransformationen haben Normquote1; die alte Occupation/Windungszählantwort1-q ist eine andere Einfügung (:125–152). N kann deshalb nicht bloß in die Stelle eines unbekannten allgemeinen retained Wortgewichts kopiert werden.

`RESIDUE-REFERENCE-AND-FULL-INTEGRATION.txt:19–56` im bestehenden Craig Contract gibt zusätzlich den tatsächlichen geometrischen Faktor `R_geo=exp(6chi0² E_D-Gamma_D4)` mit bereits gebundener Referenz und Eta. Dessen ausgeschriebene lokale Koeffizienten sind vorhanden. Der gleiche Originalreview erhält den exakten Unterschied von Res/FP zu einer ganzen relativen Determinante; eine Verdoppelung der Referenz oder ein Weglassen des vollen Wärmespurrests ist unzulässig.

Der gemeinsame vollständige Quellenwert enthält solche Faktoren an ihrem tatsächlichen Platz. Keiner der gelesenen Producer identifiziert allein dadurch eine ganze primitive gemeinsame Feldintegration mit einem ausgewerteten T_V. Das ist ein begrenzter Input/Consumer Befund; keine Behauptung physischer Mehrdeutigkeit oder einer notwendigen neuen Zustandswahl.

## Konsequenz für die Fortsetzung

Die trennenden Diagonaloperatoren können unter dem ursprünglichen vollständigen Labelvertrag konstruiert werden. Im v330 Zweig sind sogar alle Spektral, stationären und ganzen geordneten Labelgewichte explizit ausgewertet. KMS/CAR und Quillen liefern ebenfalls konkrete vorhandene Gewichtregeln in ihren tatsächlichen Quellenklassen. Diese Ergebnisse wiederverwenden und die drei Gewichttypen getrennt führen.

Der zusätzliche Schritt zur verlangten gemeinsamen ursprünglichen Dynamik ist eine tatsächliche Gleichheit des entsprechenden Quellenaufrufs beziehungsweise seiner ganzen Feld/Selektoreinfügungen. Die Behauptung einer weiteren unabhängig zu erfindenden numerischen Diagonalgewichtsfunktion gehört nicht vor diese bereits vorhandenen relativen Wege.
