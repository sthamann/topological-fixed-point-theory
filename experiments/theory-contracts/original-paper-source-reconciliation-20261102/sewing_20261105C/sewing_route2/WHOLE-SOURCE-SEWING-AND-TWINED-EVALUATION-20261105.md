# Tatsächliche Quellenverkettung und ausgewertete R104 Einfügungen

5. November 2026. Alle kanonischen Dateien bleiben unverändert. Es wird kein neues Gewicht, Spektrum, Feldmodell, Zustand oder Zeitparameter gewählt.

## Bereits geschlossene ursprüngliche Verkettungen

1. `FULL-SOURCE-CORRECTION.txt:61–90` im ursprünglichen lokalen Quellencontract gibt auf derselben tatsächlichen Basis `I_N` die geordnete Identität `Z_rel[F]=omega_R(w_R^-1 w_A F)`, alternativ `omega_A(F)/omega_A(w_A^-1 w_R)`. Die Rechnung schließt Einfügungen und den Nullquellfaktor ein. Sie verlangt keine zweite freie Zustandswahl. Nullstellen behalten die ursprüngliche Integration, nicht eine ungültige Quotientenfortsetzung.
2. `REFERENCE-SOURCE-COMPLETION.txt:75–108` trägt den gesamten tatsächlich gegebenen S/Z Aufruf durch den vollständigen festen Chart U=(E,NE): Bosonfelder und ihr Gewicht, Berezinmaß, Ref, Eta, Domänen und Einfügungen bleiben erhalten. `Z_rel^U[J_new]=Z_rel[J]` gilt regulatorweise, ebenso jedes ursprüngliche Mehrzeitwort. Diese semantische Feld und Referenzabbildung ist geschlossen. Ein zusätzliches Gibbs Referenzbild ist optional.
3. `DET16-WHOLE-SOURCE-REDUCTION.txt:135–178` gibt die wirkliche Partials Spur des im Original separierten Mirrorfaktors. Mit tatsächlicher Produktgaugeaktion ist die ganze Quelle

`Z_phys[J] = |G|^-|Vint| sum_g Tr_R[U_R(g) exp(-beta H_R) ordered_rest_source] product_x W_m,x[g_x; ordered_mirror_source]`.

Das ursprüngliche Quellenverhältnis behält seinen tatsächlichen Referenznenner. Der Resttrace wird nicht verworfen oder durch ein neu gewähltes Gewicht ersetzt. Die ursprüngliche r(g) Gaugekarte muss wirklich geliefert sein. Die bekannte Spin10 Zentrumseinfügung ist kein Nachweis der v1013 Gaugeeinbettung.

## Jetzt tatsächlich ausgewertete Erweiterung des R104 Readers

Quelle ist unverändert `H_m=2(I-P_phi)`, `phi=(|0>+exp(i vartheta)|F>)/sqrt(2)`, `vartheta=v1002.PHI=pi/7`, auf den ursprünglichen16 CAR Moden. Für jedes ursprüngliche Intervall bleibt `E_t=a_t I+b_t P_phi`, `a_t=exp(-2t)`, `b_t=1-a_t`.

Der vorhandene zyklische Rang1 Reader wird um die tatsächliche Randoperation `U_q=q^N`, q=1,i,−1,−i erweitert. Sie wird am ursprünglichen Trace-Rand eingefügt, ohne zusätzliche physische Projektion zwischen die CAR Buchstaben. Jeder originale Rang1 Zwischenblock wird mit dem echten Zwischenbesetzungsfaktor ausgewertet. Der Null-P-Term behält den vollständigen Fockcharakter. Eine Division durch `(1+q)` oder einen Determinanten erfolgt nirgends.

Erste unabhängige Ausführung: 128 ursprüngliche native sparse Tracevergleiche, 93/93 zusätzliche Kontrollen bestanden. Nach Ergänzung des gaugeinvarianten Gesamtzahlwortes: **138 native sparse Tracevergleiche,105/105 Kontrollen bestanden**. Größte absolute Abweichung2.3283e−10 bei großen unnormierten Zahlworttraces, relative Abweichung in der Größenordnung1e−15. Die Annahmeschranke des Vergleichs ist `2e-10*max(1,abs(native))`; bei kleinen determinantischen Antworten gelten zusätzliche absolute Kontrollen. Keine gleichmäßige numerische Fehlergarantie für beliebige Wortlängen wird behauptet.

Getestet: ursprüngliche Buchstabenwörter Länge0,2,4,8,16,32 bei diagnostischen beta0,1,2pi; alle vier originalen Zentrumseinfügungen; ungleichmäßige und gleiche16er Zeitintervalle; echte globale Zentrumprojektion P4. Die allgemeine Aussage folgt aus der ursprünglichen Rang1 Zerlegung und der vollständigen CAR Bitwirkung, nicht aus dieser endlichen Testanzahl.

### Ganze tatsächliche16er Zeitantwort

Für das ursprüngliche Wort `c_0 ... c_15`, Intervalle t_0,...,t_16 und beta=sum t_j gilt

`W_q = exp(i vartheta)/2 sum_(j=0..15) (-1)^j q^(16-j) b_(delta_j) exp[-2(beta-delta_j)]`,

mit `delta_0=t_0+t_16`, `delta_j=t_j` für j>0.

Beweis: Genau ein P Schnitt kann beitragen. Jeder nichtleere echte Teilblock von1 bis15 Annihilatoren hat verschwindendes phi Matrixelement. Beim Schnittj entsteht die wirkliche zyklische CAR Permutationssignatur `(-1)^(j(16-j))=(-1)^j`; die Zentrumsoperation sieht N=16−j. Der Koeffizient b und alle übrigen ursprünglichen Wärmefaktoren bleiben erhalten. Zwei unabhängige Reviewer haben diese Randphasenrechnung hergeleitet; der native sparse Vergleich ist eine andere Rechnung.

Bei Gleichzeitigkeit bleibt nur delta_0=beta und die bekannte gleiche Zeitkohärenz `b_beta exp(i vartheta)/2`. Bei gleichen zyklischen Zeitintervallen verschwindet die ungetwinte Antwort wegen der Alternationssumme; die q=−1 Antwort tut dies nicht. Die Mehrzeitquelle enthält also tatsächlich zusätzliche Randphaseninformation.

### Tatsächlicher Nullstellensektor und globale Zentrumprojektion

`Z_q[0]=a_beta(1+q)^16+b_beta`. Bei beta0 und q=−1 ist dieser Wert0. Das tatsächlich markierte Wort `product_(j=0..15) n_j` besitzt dort trotzdem unnormierten Koeffizienten1. Eine Normierung durch den verschwindenden Zentrumscharakter wäre falsch. Der Reader erhält den Koeffizienten ohne solche Division.

Der tatsächliche globale Zentrumprojektor `P4=(1/4)sum_q U_q` hat Rang16512 und

`Z4=16512 a_beta+b_beta=1+16511 exp(-2beta)`.

Seine ganzen markierten Koeffizienten sind genau das Mittel der vier berechneten W_q. P4 bleibt am ursprünglichen Rand. Dies ist die kanonische Zentrumprojektion, kein voller Spin10, v1013 oder P1 Gaußnachweis. Einzelne geladene CAR Buchstaben werden dadurch nicht automatisch zu lokalen physischen Observablen; geordnete Koordinatenquellen und tatsächlich gaugeinvariante Beine werden getrennt.

### Ein tatsächlich gaugeinvariantes ganzes Wort

Der ursprüngliche Gesamtzahloperator `N=sum n_j` kommutiert mit jeder kanonischen Gamma(r). Es gilt `<phi|N|phi>=8`, `<phi|N²|phi>=128`. Für s+t=beta folgt für die ganze Zentrumquelle exakt

`C_N,4(s,t) = [64(exp(-2s)+exp(-2t))+58240 exp(-2beta)]/[1+16511 exp(-2beta)]`.

Dies ist der verbundene thermische Zweizeitkoeffizient; `<N>4=8`. Die unnormierte Quelle lautet

`W_N,4 = 1115136 a_s a_t+128(b_s a_t+a_s b_t)+64 b_s b_t`.

Diese unabhängig erhaltenen Formeln wurden an den tatsächlichen nativen65536er Diagonaloperatoren und dem originalen P_phi geprüft. Diagnostisch s=.3,t=.7,beta1 ergibt `C_N,4=3.548539037434813`; bei beta2pi,s=.3beta,t=.7beta ergibt sich `1.5963310797881707`. Das ist keine Auswahl einer physischen Temperatur und keine primitive Gesamttheorievorhersage.

## Ganze kanonische Gamma(r) Randoperation konstruktiv geschlossen

Für jedes tatsächlich gegebene r der ursprünglichen kanonischen Klasse gilt in der originalen Exteriorbasis

`<T|Gamma(r)|S> = det r[T,S]` bei gleicher Besetzungszahl, sonst0.

Jeder ursprüngliche CAR Monomblock B ist auf Bitzuständen eine partielle vorzeichenbehaftete Bijektion. Daher wird jeder Rang1 Block `<phi|B_L Gamma(r) B_R|phi>` durch vier tatsächliche Endpunktpaare x,y in{0,F} berechnet: Setze `B_R|y>=a_R|S>` und `B_L*|x>=a_L|T>`, dann ist der Summand

`conjugate(phi_x*a_L) phi_y a_R det r[T,S]`.

Das ist eine konkrete endliche Rechnung ohne Feldexpansion oder inversen Kovarianzkern. Im Null-P-Term liefert die ursprüngliche CAR Bitwirkung `B|S>=sign(S)|S'>` die volle Spur

`Tr(Gamma(r)B)=sum_S sign(S) det r[S,S']`.

Ungültige CAR Ausgangsbesetzungen fallen weg. Es gibt höchstens2^16 originale Bitzustände; bei n benutzten Moden nur2^(16−n) mögliche freie Anfangsbits. Der tatsächliche Umfang ist fix16, keine Suche über neue Modelle. Der ganze zyklische Reader bleibt unverändert; lediglich die beiden tatsächlichen Blockauswerter werden übergeben. Er enthält weiterhin jeden inneren Grad, den ursprünglichen Zeitparameter und den absoluten Nullquellfaktor.

Konstruktive Aussage für die ganze ursprüngliche kanonische R104 Klasse: Der bisherige Faktor `a det(I+r(g))+b` lässt sich für jedes markierte geordnete Mirrorwort durch diesen vollständig berechenbaren W_m(g;word) ersetzen. Er kann direkt in die ursprüngliche Gaußsumme eingesetzt werden, sobald deren tatsächliches r(g) gegeben ist. `det(I+r)=0` ist kein Problem; kein Schritt dividiert dadurch. Ein nichtkanonischer Implementer oder eine andere Feld/Gaußkarte behält ihre ursprünglichen Voraussetzungen.

### Unabhängige native Gegenprüfung jenseits des Zentrums

Verwendet wurde der tatsächliche positive Weyl Spin10 Generator G_(0,2) von v1002 und dessen ursprünglicher diagnostischer Winkel .2718281828459045. r=exp(i angle G) ist nichtzentral, unitär und hat Determinante1. Die Gegenrechnung benutzt ausdrücklich **keine Exterior-Minoren**: Der ursprüngliche native second_quantize erzeugt dGamma(G), und sparse Exponentialwirkung `exp(i angle dGamma(G))` wird auf den wirklichen Bitbasisvektoren ausgeführt.

Native gesamte ursprüngliche Wärmewortmatrizen wurden für density32, offdiagonal30 und Omega16 verwendet. Ergebnisse native versus Minor-Renewal:

| Wort | Unnormierter nativer Koeffizient | Absolute Abweichung |
|---|---|---|
| density32 | .22892760363074674 + .0000000000000000000051 i | 1.67e−16 |
| offdiagonal30 | .01833739347023202 | 3.47e−18 |
| Omega16 | .008237638259284204 + .003967037504589112 i | 1.91e−17 |

**20/20 Kontrollen bestanden**, einschließlich des singulären r=−I, und der Rückkehr zu allen ursprünglichen Zentrumwortantworten. Die allgemeine Klassenkonstruktion beruht auf den gerade bewiesenen Minor und Rang1 Identitäten. Diese drei nichtzentralen nativen Fälle sind begrenzte unabhängige numerische Gegenprüfungen.

## Was dadurch verbunden ist und was getrennt bleibt

Die tatsächliche ganze R104 Gauge-Wortverkettung ist jetzt auch mit markierten Mirrorbeinen konstruktiv auswertbar. Die alten ursprünglichen whole U und Rückgewichtungsabbildungen bleiben geschlossen. Ihre Anwendung verlangt keinen neuen Refzustand, keine neue Gewichtsfunktion und keine erneute Auswahl von omega.

Die Quelle der R104 Tabelle ist aber nicht allein durch diese Verkettung mit der primitiven spektralen S_adm Quelle identifiziert. Die direkte vollständige U Karte trägt dieselbe gegebene ursprüngliche Quelle; sie macht aus einer bedingten quasifreien Cliffordquantisierung nicht den nichtquasifreien DET16 Parent. R104 benennt selbst seine tatsächliche Feld/Zeiteinfügungs und Gaugeantezedenz. Die behauptete volle v1013 Mirror Gaugeeinbettung ist in den gelesenen Originaltable/Provenienzstellen nicht ausgeführt. Die DET4/Z6 Schattenprovenienz ist ausdrücklich keine volle Spin10 oder dynamische Linkquelle. Kein solcher Ersatz wird hier gewählt.

## Tatsächliche Artefakte

- `twined_original_det16.py`, `.json`, `.log`: ganzer Zentrumwortreader mit nativer Gegenrechnung, P4 und gaugeinvariantem N Zweizeitwort.
- `gamma_original_det16.py`, `.json`, `.log`: ganze kanonische Gamma(r) Einfügung mit Exterior-Minoren und unabhängiger echter sparse Exponentialwirkung.
- Kanonischer unveränderter Reader: `experiments/theory-contracts/original-local-source-variable-background-20261102/check_det16_source.py:74–130`.
- Kanonische native Generatoren: `verification/v1002_spin10_mirror_projector.py:114–152,164–201`.

Kein gemeinsames primitives S_adm Integral, vollständiges physisches Spektrum oder 3+1D Gesamtabschluss wird aus diesen Auswertungen behauptet. Die ganze vorliegende Quellenverkettung hat aber jetzt tatsächliche markierte Werte und einen auswertbaren allgemeinen Mirrorfaktor.
