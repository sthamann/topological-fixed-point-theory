# Ganze nichtlineare ursprüngliche R104 Zahlquelle

5. November 2026. Tatsächliche Auswertung der vorhandenen R104 Quelle. Kein neues Gewicht, Spektrum, Zustand, Oszillator oder physischer Zeitmaßstab wird gewählt. H, N und die Gaugeimplementierung sind die bereits angegebenen Originalobjekte.

## Ergebnis

Für die ganze ursprüngliche R104 Klasse ist das erzeugende Funktional sämtlicher zeitgeordneter Gesamtzahlantworten konstruktiv auf einen sechzehndimensionalen Charakter und eine zweidimensionale zeitgeordnete Matrixentwicklung reduziert. Die Reduktion gilt für beliebige stückweise konstante reelle oder komplexe Quellenprofile, ohne Taylortrunkierung. Sie ergänzt den schon ausgewerteten ganzen Mirrorwortreader um eine vollständige nichtlineare gaugeinvariante Quellenfamilie.

Die allgemeine ganze Gaußquellenverkettung bleibt die ursprüngliche R104 Summe mit unverändertem Resttrace und tatsächlichem Referenznenner. Bei tatsächlich gegebener kanonischer Mirroraktion r(g) kann für jede Zelle direkt der unten berechnete volle Faktor Z_r[h] eingesetzt werden. Eine zuvor unbelegte v1013 Gaugeeinbettung oder primitive S_adm Identität folgt daraus nicht.

## Originale und feste Operatoren

`DET16-WHOLE-SOURCE-REDUCTION.txt:12–43,135–178` und der ursprüngliche `check_det16_source.py` geben auf derselben65536 dimensionalen Fockfaser

`H=2(I-P_phi)`, `phi=(|0>+exp(i vartheta)|F>)/sqrt(2)`, `vartheta=v1002.PHI=pi/7`, `N=sum_(j=0..15) n_j`.

Die ursprüngliche Gaugeaktion ist Gamma(r), r unitär und det r=1, soweit diese tatsächliche Darstellung im Quellenaufruf gegeben ist. N kommutiert mit jeder kanonischen Gamma(r). Die Wärmezeiten sind dieselben ursprünglichen Zeiten; beta wird durch die Intervalle bestimmt. h ist ausschließlich der äußere Quellparameter, keine neue physische Kopplung oder Kalibrierung.

## Ganze Quellenformel und Beweis

Der Raum K=span{|0>,|F>} und sein Komplement sind unter H, N und Gamma(r) invariant. Direktes Einsetzen des tatsächlichen P_phi ergibt

```
H2 = [[1, -exp(-i vartheta)],
      [-exp(i vartheta), 1]],
N2 = diag(0,16),
Gamma(r)|K = diag(1,det r).
```

Auf dem Komplement ist H=2I, während N nach der tatsächlichen Teilchenzahl diagonal ist. Für ein Profil aus Intervallen delta_t_j mit Quellwerten h_j sei beta=sum delta_t_j und chi=sum h_j delta_t_j. In der dargestellten Operatorreihenfolge folgt daher exakt

`Z_r[h] = exp(-2 beta)[det(I+exp(chi)r)-1-exp(16chi)det r]`

`          + Tr2[diag(1,det r) product_j exp(-delta_t_j(H2-h_j N2))]`.

Die Profil LISTE in `corner_source` und `native_full_source` ist ausdrücklich die ursprüngliche **Operatorproduktreihenfolge von links nach rechts**. Wird eine Quelle dagegen chronologisch von früh nach spät eingegeben, ist für die übliche Zeitordnung mit späteren Zeiten links die Liste umzukehren. Der geordnete Produktwert und der native Vergleich verwenden dieselbe explizite Reihenfolge; eine chronologische Eingabeliste wird nicht stillschweigend als zeitgeordnetes Produkt ausgegeben.

Der Charakter `Tr Gamma(exp(chi)r)=det(I+exp(chi)r)` summiert alle ursprünglichen Besetzungsgrade. Die beiden subtrahierten Eckbeiträge werden genau einmal durch die wirkliche wechselwirkende zweidimensionale Entwicklung ersetzt. Kein charakteristischer Faktor wird dividiert. Deshalb gelten dieselben Formeln auch bei det(I+r)=0.

Die Aussage betrifft die ganze Profilquelle, nicht eine abgeschnittene Anzahl ihrer Ableitungen. Für allgemeine beschränkte zeitabhängige Quellen erhält man den ursprünglichen geordneten ODE Propagator durch die üblichen normkontrollierten Profilapproximationen auf dieser endlichen Faser; die stückweise konstante Formel ist schon vollständig ausführbar. Eine neue unendlichdimensionale Grenzaussage wird hier nicht vorausgesetzt.

## Tatsächliche gesamte Zentrumquelle

Für den ursprünglichen globalen Zentrumprojektor `P4=(1/4)sum_(q=1,i,-1,-i)q^N` liefert der Binomialcharakter

`(1/4)sum_q(1+q exp(chi))^16 = 1+1820exp(4chi)+12870exp(8chi)+1820exp(12chi)+exp(16chi)`.

Folglich ist die ganze ursprüngliche Zentrumquelle

`Z4[h]=exp(-2beta)[1820exp(4chi)+12870exp(8chi)+1820exp(12chi)]`

`       +Tr2 product_j exp(-delta_t_j(H2-h_j N2))`.

Bei h=0 reduziert sie auf das ursprüngliche `Z4=1+16511 exp(-2beta)`; bei beta0 auf den tatsächlichen Rang16512. Die globale Projektion wird nicht zwischen die Einfügungen gesetzt. Die Zahlquelle selbst ist gaugeinvariant innerhalb der kanonischen Klasse; einzelne geladene CAR Buchstaben erhalten dadurch keine automatische physische Lokalitätsabbildung.

## Konstantes Profil, vollständig nichtlinear

Für konstantes h ist der Mittelwert der beiden Eigenwerte von H2-hN2 gleich1−8h. Sein spurloser Teil B erfüllt exakt `B²=(1+64h²)I`. Daher

`Tr2 exp[-beta(H2-hN2)] = 2 exp[-beta(1-8h)] cosh(beta sqrt(1+64h²))`.

Dieser Ausdruck ist eine ganze Funktion von h; die scheinbare Wurzelverzweigung verschwindet wegen der geraden Coshfunktion. Beim tatsächlich getesteten komplexen h=i/8 ist die Diskriminante0 und der Ausdruck bleibt regulär. Kein Tayloransatz und keine Division durch die verschwindende Diskriminante werden benutzt.

## Zusammenhang mit den ursprünglichen Zahlwörtern

Die ganze Kickquelle `Tr[Gamma(r) E_s exp(xN) E_t exp(yN)]` besitzt dieselbe Komplement plus Ecke Zerlegung. Ihr xy Koeffizient ist exakt das bereits unabhängig hergeleitete geordnete N Zweizeitwort. Im Eckblock trägt dieser Koeffizient **auch** den Beitrag256 a_s a_t der beiden ursprünglichen Koordinatenzustände. Dieser Anteil wird in der ganzen Charakterformel durch den korrespondierenden Eckabzug genau korrekt verarbeitet.

Für P4 gilt deshalb weiterhin

`W_N,4(s,t)=1115136 a_s a_t+128(b_s a_t+a_s b_t)+64 b_s b_t`,

`C_N,4(s,t)=[64(a_s+a_t)+58240a_s a_t]/[1+16511a_s a_t]`,

mit a_x=exp(-2x), b_x=1-a_x. Der echte Gleichzeitgrenzfall ist `Tr(P4 exp(-beta H)N²)=1115136a_beta+128b_beta`. Er wird nicht mit einem über die Zeit integrierten Konstantquellenkoeffizienten verwechselt.

## Tatsächlich ausgeführte unabhängige native Prüfung

Die Eingabe der nativen Gegenrechnung ist die **volle ursprüngliche65536er sparse Matrix** H=2(I-P_phi) und der vollständige65536er Diagonaloperator N. Sie benutzt weder die angezeigte H2 noch die Charakterformel für die Traceauswertung. Ihre tatsächliche sparse Graphstruktur wird algorithmisch bestimmt:65535 Komponenten, davon eine nichttriviale Komponente der Größe2. Die Exponentialmatrizen werden aus diesen tatsächlich gefundenen Komponenten zusammengesetzt und auf dem ganzen Originalraum multipliziert. Die erwartete H2 ist nur ein davon getrennter Vergleich.

Zusätzlich wurde die volle sparse Exponentialwirkung ohne jede Komponentenzerschneidung für das nichtkonstante komplexe Profil auf beiden nichttrivialen Basisbeinen geprüft; sie stimmt mit dem nativen ganzen Matrixprodukt überein.

Der erste Versuch mit generischem sparse Pade expm auf der gesamten Matrix wurde im spsolve Schritt kontrolliert unterbrochen: dieser Pfad arbeitete65536 dünne rechte Seiten ab. Es wurden daraus keine fertigen Endwerte behauptet. Die komponentenbasierte native Rechnung behebt nur diese numerische Darstellungsschwäche derselben Matrix; sie ändert keinen Quelleninput. Eine anfänglich unvollständige Eck N2 Hilfsformel wurde um256a_s a_t berichtigt. Die tatsächlichen ganzen Profilvergleiche bestanden bereits; nach der Korrektur sind alle folgenden Endkontrollen bestanden.

**Endstand:69/69 Kontrollen bestanden.42 native ganze Profiltracevergleiche**, aus sieben Profilen und je sechs tatsächlichen Randeinsätzen: die vier ursprünglichen Zentren, P4 und ein nichtzentraler ursprünglicher Spin10 Cartaneinsatz. Zusätzlich zwei native N2 beziehungsweise Gleichzeitkoeffizienten und die volle sparse Exponentialwirkung geprüft. Nach diesem bestandenen Lauf erfolgte kein weiterer Testlauf.

Profile: h0, h.1, h−.2, h.1+.07i, h=i/8 sowie zwei nichtkonstante dreistufige reelle beziehungsweise komplexe Profile. Alle Quellparameter sind diagnostische Einfügungen, keine physische Temperatur oder h Kalibrierung.

- Maximale absolute Abweichung der42 Profiltraces:1.6528599314e−9.
- Maximal skalierte Abweichung der Profiltraces:7.2051407266e−13, mit Skala max(1,|native|).
- Die separaten großen N2 Grenztraces haben einen absoluten Summationsrest7.8289303929e−9 bei Traces um1.5e5. Das kompakte JSON Feld `max_absolute_error` enthält dieses tatsächlich größere Maximum einschließlich dieser zusätzlichen Kontrollen.
- Native nichtkonstante reelle Zentrumquelle bei beta1, Intervallen(.2,.3,.5) und h(.1,−.2,.3): `Z4=5507.576023081181`, separate Reduktion `5507.576023081229`.
- Native nichtkonstante komplexe Zentrumquelle bei demselben Zeitprofil und h(.1+.07i,−.05+.02i,.12−.09i): `3704.472739001493−773.3106026048499i`, separate Reduktion `3704.4727390011976−773.3106026048607i`.

Diese numerischen Restgrößen werden nicht als exakte Nullen ausgegeben. Die allgemeine endliche Aussage beruht auf der tatsächlichen invarianten Raumzerlegung und dem ganzen Quellenoperator; die native Rechnung ist eine unabhängige begrenzte Gegenprüfung.

## Ablage und Geltungsbereich

`number_source_original_det16.py`, `.json`, `.log` im selben Ordner. JSON enthält die kompakt verwendbaren Felder `checked`, `total`, `native_comparisons`, `max_absolute_error`, `scope`, zusätzlich den getrennten Profil und Grenzcheck Fehlerumfang sowie sämtliche tatsächlich ausgeführten Profile und Werte.

Die ursprüngliche gesamte R104 Zahlquellenfamilie ist ausgewertet und kann als ganzer h abhängiger Mirrorfaktor in die bereits vorhandene tatsächliche Gaußsumme eingesetzt werden. Der übrige Restgitterprozess, seine tatsächliche physische Gaugekarte und die Identifikation mit der primitiven S_adm Quelle sind dadurch nicht still berechnet. Alle bisherigen ursprünglichen Quellen, Referenz und Feldtransporte sowie gegebene Zustands und Zeitdaten bleiben erhalten.
