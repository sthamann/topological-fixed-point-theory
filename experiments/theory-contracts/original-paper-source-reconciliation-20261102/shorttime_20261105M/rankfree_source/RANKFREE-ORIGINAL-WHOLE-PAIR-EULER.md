# Rangfreie Auswertung desselben ursprünglichen Gegenmoduspaares

Stand 2026-11-05. Dies ergänzt die bereits geprüfte ganze96 Produktformel und ihre physische384 Gegenmoduspaarung, ohne einen physikalischen Yukawarank zu wählen. Die Dateien und der eingefrorene frühere Bericht in shorttime20261105K_pfaffian bleiben unverändert. Das Ergebnis betrifft den globalen euklidischen Eulerquellenmarker, nicht das lokale R113/HYPF Wort Y_ham(s).

## 1. Ursprung und tatsächliche komplexe Basisänderung

Aus der ursprünglichen v252 LR Konvention ist nach dem vorhandenen schwachen Orbittransport H=(0,sqrt(u)) der ganze Neutrinoblock

    Xi = [[0, sqrt(u) bar(Ynu)],
          [sqrt(u) Ynu^dagger, MR]],     MR^T=MR.

Der bereits erhaltene ganze96 Faktor lautet

    Q(p,u) = det3(z+u Au)^6 det3(z+u Ad)^6 det3(z+u Ae)^2 Qnu(p,u),
    Af=Yf^dagger Yf, z=p^2>0,
    Qnu=det6(z+Xi^dagger Xi).

Die tatsächliche384 Gegenmoduspaarung hat Pf M_pair=Q^2. Insbesondere ist ein einzelner Q Faktor keine zusätzlich postulierte unabhängige physische Pfaffianquelle.

Sei Ynu=L D R^dagger eine Singulärwertzerlegung mit D=diag(y1,...,yr,0,...,0), yi>0, k=3-r. Die UNITÄRE Kongruenz W=diag(L,bar(R)) gibt

    W^T Xi W = [[0, sqrt(u) D], [sqrt(u) D, M]],
    M=R^dagger MR bar(R),    M^T=M.

Die komplexe Konjugation ist durch die originale Xi Konvention gebunden. M=R^T MR R ohne gleichzeitigen Konventionswechsel wäre hier falsch. Die Determinante z+Xi^dagger Xi bleibt unter der angegebenen unitären Kongruenz unverändert. Diese Basisänderung ist keine neue Quelle, Yukawawahl oder neue Familienprojektion.

Zerlege nur zur Rechnung die wirkliche ganze Matrix nach ihrem Rang:

    M=[[M11,M10],[M10^T,M00]],
    J=diag(y1^-2,...,yr^-2), G=u D^2+z I3,
    S0=Ik+M00^dagger M00/z.

Leere Produkte und Determinanten sind1; leere Spuren sind0. Bei r=0 ist J leer und die spätere Rückkehrspur0. Bei r=3 ist die Nullfaser leer.

## 2. Ganze Neutrinodeterminante und ihre beiden höchsten Koeffizienten

Der ursprüngliche ganze6 Block ergibt über Xi Xi^dagger und den exakten Schurkomplementausdruck

    Qnu=det(G) det(G+z M G^-1 M^dagger).                    (1)

Da z>0 ist G invertierbar für jeden Rang und jedes u>=0. Die Formel enthält keine Pseudoinverse einer physisch singulären Masse.

Eine ausdrücklich endliche, polynomiale Form von (1) ist

    Qnu = sum_(ell=0)^3 z^ell
          sum_(|I|=|Jset|=ell) |det M[I,Jset]|^2
          product_(i notin I) gi product_(j notin Jset) gj,
    gi=u yi^2+z (i aktiv), gi=z (i null).                 (2)

Dies folgt durch zwei Cauchy-Binet Anwendungen. Die sechs ursprünglichen Neutrinoplätze und alle Majoranaeinträge bleiben erhalten. Es sind20 Minorenprodukte, keine neu gewählte kleine physische Matrix. (2) ist eine direkt auswertbare ganze Neutrinopolynomregel für alle Ränge.

Setze Delta=product_(i=1)^r yi^2=pdet(Ynu^dagger Ynu). Dann

    Qnu=a_nu u^(2r)+b_nu u^(2r-1)+O(u^(2r-2)),
    a_nu=z^k Delta^2 det(z Ik+M00^dagger M00),
    b_nu/a_nu=2z Tr(J)+2 T_M,
    T_M=Tr[J M10 S0^-1 M10^dagger].                      (3)

Bei r=0 ist Qnu unabhängig von u, a_nu=z^3 det(z I3+MR^dagger MR) und b_nu=0. Die negative Potenz in der allgemeinen Schreibweise ist dann kein zusätzliches Polynomglied. T_M>=0, weil z>0 und S0 positiv definit sind.

Zur unabhängigen Herleitung des zweiten Koeffizienten:

    G^-1=P_null/z+J/u+O(u^-2),
    R0=z Ik+M00 M00^dagger.

Das zweite Schurkomplement aus (1) hat in der Nullfaser den ersten Korrekturblock

    z M10^T J bar(M10)
       -M00 M10^dagger J M10 M00^dagger.

Die Aktivfaserdeterminante trägt Tr[J(z I+M10 M10^dagger)] bei; det(G) liefert ein weiteres z Tr J. Mit

    I-M00^dagger R0^-1 M00
       =z(z I+M00^dagger M00)^-1

werden der Aktivfaserterm mit M10 M10^dagger und der negative Rückkehrterm zusammen zu T_M. Der verbleibende positive Nullfaserterm ist ebenfalls T_M: hierfür wird tatsächlich M00^T=M00 eingesetzt und die Spur transponiert. Ohne diese ursprüngliche komplexe Symmetrie wäre die Gleichsetzung nicht begründet. M11 fällt in dieser Ordnung aus; es wurde nicht aus der Quelle gelöscht.

Für r=2 ist M00 ein komplexer Skalar und

    a_nu=z Delta^2(z+|M00|^2),
    T_M=z/(z+|M00|^2) ||D_active^-1 M10||_F^2.

Damit kann Majorana bei singulärem Ynu bereits im führenden und im nächsten Koeffizienten auftreten. Die alte Vollrangbehauptung, Majorana beginne erst eine Ordnung später, bleibt genau auf ihren Vollrangzweig beschränkt.

## 3. Alle ursprünglichen Familien und Farben, ohne Rangannahme

Für jedes f=u,d,e,nu sei rf=rank Yf und A_f^+ die Moore-Penrose Inverse auf den tatsächlichen positiven Singulärwertplätzen. Die geladenen Faktoren erfüllen exakt

    det3(z+u Af)^sf
      =a_f u^(sf rf) [1+sf z Tr(A_f^+)/u+O(u^-2)],
    a_f=z^(sf(3-rf)) pdet(Af)^sf,
    (su,sd,se)=(6,6,2).

Definiere

    kQ=6ru+6rd+2re+2rnu,
    Sigma_plus=3Tr(Au^+)+3Tr(Ad^+)+Tr(Ae^+)+Tr(Anu^+),
    S=z Sigma_plus+T_M,
    aQ=a_u a_d a_e a_nu.

Dann lautet die tatsächliche ganze96 Polynomantwort

    Q=aQ u^kQ+2aQ S u^(kQ-1)+O(u^(kQ-2)),

und die ganze384 Gegenmoduspfaffianantwort

    Pf M_pair=aQ^2 u^(2kQ)
              +4aQ^2 S u^(2kQ-1)+O(u^(2kQ-2)).          (4)

Für kQ=0 sind alle Yukawas0, S=0 und die ganze Quelle ist H unabhängig. Weil z>0 ist aQ>0 für jeden tatsächlichen Rang. Die neue Formel braucht weder einen geratenen Rang noch eine Division durch det Ynu. Die Spuren auf positiven Singularplätzen und die auf der ursprünglichen Nullfaser berechnete Rückkehrspur sind gemeinsam durch den Koeffizientenquotienten aus (1) basisunabhängig.

## 4. Vollständig integrierter globaler Eulerkontakt in derselben IR Klasse

Behalte den ursprünglichen isotropen IR Faktor exp(-r_IR u/2) in den vier reellen Koordinaten der gegebenen Higgsdublette. Der gesamte gehaltene Rest w(H) sei beschränkt und habe unter H=y/sqrt(r_IR) einen dominierten Grad0 Strahlengrenzwert w_infinity(y/|y|). Sein führender Winkelintegralwert sei nichtzero. Richtungsabhängige und komplexe w sind damit weiterhin erlaubt; alle gegebenen Referenzfaktoren bleiben im gleichen Caller.

Der Compositequellenparameter g multipliziert ausschließlich den originalen H linearen Fermionvertex bei festgehaltenen physikalischen Yukawaparametern, Bosonkoeffizienten, Majorana, Referenz und Eta. Es ist die ursprüngliche Halbgewichtpaarung, so dass d_g P=-I_fermion[S_Y^E] gilt. Diese Ableitung wird nicht als Variation sämtlicher Y abhängigen Terme umgedeutet.

(4) hat führenden H Homogenitätsgrad4kQ. Dieser konstante Euleranteil fällt vor jedem Grenzübergang im verbundenen Ausdruck aus. Für kQ>=1 ist das einzige radiale Momentverhältnis

    E_gamma4[u^(2kQ-1)]/E_gamma4[u^(2kQ)]
       =1/[2(2kQ+1)].

Daher ist B/A=2S/(2kQ+1) und die tatsächliche radiale Integration ergibt

    lim_(r_IR downarrow0) Cov_rIR(U,S_Y^E)
       = -4 B/A
       = -8 [p^2 Sigma_plus+T_M]/(2kQ+1).               (5)

Für kQ=0 gilt derselbe Ausdruck als0. Es ist keine Konvergenzrate des gehaltenen Restgewichts nötig, da die führende Eulerkonstante bereits algebraisch verschwindet.

Alle Ränge3 geben kQ=48, T_M=0 und exakt den früheren Wert -8p^2 Sigma/97 zurück. Sind die geladenen Ränge3 und der tatsächliche Neutrinorang2, dann ist kQ=46 und

    C=-8 p^2/93 * [Sigma_plus
                  +||D_active^-1 M10||_F^2/(p^2+|M00|^2)].

Die Berechnung gilt für p^2>0. Der Fall p=0 wird nicht durch Vertauschen der beiden Grenzwerte eingeschmuggelt; die singuläre reine Massenfaser hat ihren eigenen bereits bekannten Polynomgrad und kann verschwinden.

Der Quellengeltungsbereich ist unverändert: freie gewöhnliche Spin4 Produktpaarung, konstantes ganzes H, tatsächliches Gegenmoduspaar, globaler euklidischer Eulerkontakt und der benannte begrenzte Strahlengrenzwert. Die ursprüngliche ungekürzte positive R048 Quartik exp(-lambda u^2) hat Strahlengrenzwert0; (5) wird nicht auf diesen anderen Caller übertragen. Nichtkonstante H Felder koppeln Impulse und erlauben kein Produkt unabhängiger Hintergrundmoden. Das lokale R113 Y_ham(s), der ganze Modenturm, die gekoppelte Geometrie und die ganze primitive markierte OS Historienquelle sind hier nicht integriert. Die neue Rechnung beseitigt die Rangvoraussetzung der vorhandenen Modenantwort; sie behauptet keine gesamte primitive Lösung.

## 5. Prüfung und Originalbezug

check_rankfree_original_pfaffian.py prüft die beiden höchsten Koeffizienten aus der endlichen20-Minoren-Polynomregel für alle vier möglichen tatsächlichen Neutrinoränge symbolisch. Sämtliche sechs komplexen symmetrischen Majoranaeinträge und alle positiven Singulärwerte bleiben Symbole. Zusätzlich prüft es die notwendige komplexe Schurrückkehridentität, die Zählung der ganzen Familien/Farben und die exakte radiale Eulerintegration. Keine numerischen physikalischen Yukawas, Majoranaeinträge, Impulse oder Ränge werden gewählt. Ergebnis:24/24 exakte symbolische Kontrollen bestanden.

Die erste ganze96/384 Paarungsidentität und die v252 Konvention wurden zuvor in shorttime20261105K_pfaffian gegen die tatsächlich ursprünglichen Operatorblöcke geprüft (64/64 Kontrollen). Dieser Zusatz verändert weder den originalen Operator noch die damaligen Dateien. Die neue Rechnung beginnt an dessen ganzer6 Neutrinodeterminante, nicht an einer neuen physikalischen Projektion.

Originalquellen mit aktuellen Hashes stehen in provenance.json: verification/v252_full_finite_triple.py:80–128; experiments/theory-contracts/craig-original-source-bridge-20261103/ACTUAL-PRODUCT-PAIRING.txt:48–105; tatsächlicher actual_product_pairing.py; experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/DYNAMICAL-HIGGS-SOURCE-PROOF.txt:34–78,121–144. Der Rangbelegcheck aus K bleibt erhalten: es wird weder tatsächlicher Vollrang noch tatsächlicher Rang2 des physischen Ynu behauptet.

## 6. Strikt negativer verbundener Eulerkontakt ohne Strahlengrenzwertannahme

Dieselbe bereits ausgerechnete ganze Polynomquelle liefert eine zusätzliche qualitative Antwort ohne IR Grenzwert. Für z=p^2>0 hat sie

    Q(u)=sum_(j=0)^kQ a_j u^j,     a_j>=0,
    Q(0)=z^45 det3(z I3+MR^dagger MR)>0.                 (6)

Die nichtnegativen Koeffizienten folgen direkt aus der ganzen20-Minorenregel (2): jeder Faktor gi besitzt nichtnegative Koeffizienten und jedes Minorenbetragsquadrat ist nichtnegativ. Die geladenen Faktoren besitzen ebenfalls nichtnegative Koeffizienten, weil Af positiv semidefinit ist. Ist irgendeine ursprüngliche Yukawamatrix nichtzero, dann kQ>0 und zugleich a_kQ>0. Damit besitzt Q mindestens zwei unterschiedliche nichtverschwindende Grade. Kein physikalischer Rang ist dafür zu wählen.

Lasse jetzt dnu(H) den tatsächlich gegebenen restlichen H Integrationscaller einschließlich aller gehaltenen Quellen, Koordinatenfaktoren und Referenzbeine bezeichnen. Er wird nicht neu gewählt. Für diesen Satz gelte: dnu ist ein positives Maß, alle betrachteten Faktoren sind beim fermionischen Compositeparameter g fest, und

    0<Z=int Q(u)^2 dnu(H)<infinity,
    int u Q(u)^2 dnu(H)<infinity.

Nichtpositive oder komplexe ursprüngliche Restgewichte erfüllen die Satzvoraussetzung nicht. Winkelabhängige positive Gewichte sind erlaubt, da keine Winkelmittelung oder radiale Ersatzdichte gewählt wird.

Die tatsächliche affine Halbgewichtquellenableitung ist punktweise

    I_fermion[S_Y^E](H)=-d_g[Q(g^2 u)^2]_(g=1),
    y(u):=I_fermion[S_Y^E](H)/Q(u)^2
         =-4u Q'(u)/Q(u).                               (7)

Die Division ist hier für alle u>=0 legitim, denn (6) garantiert Q(u)>0. Es wird kein singulärer allgemeiner Feldpfaffian normalisiert. Definiere lediglich zur Auswertung der bereits gegebenen Quelle

    pi_j(u)=a_j u^j/Q(u).

Für u>0 ist dies eine tatsächliche positive endliche Gradverteilung und

    y(u)=-4 E_pi(u)[j],
    d_(log u)y(u)=-4 Var_pi(u)(j)<=0.                   (8)

Bei mindestens einer nichtzero Yukawa ist die Varianz für alle u>0 strikt positiv. Ferner y(0)=0 und -4kQ<=y(u)<=0. Daher ist y auf [0,infinity) strikt fallend; die gegebenen Momentbedingungen reichen auch für die markierten Integrale aus.

Nach der tatsächlichen Berezinintegration ist mit dmu=Z^-1 Q(u)^2 dnu(H)

    Cov(U,S_Y^E)=Cov_mu(u,y(u)).                         (9)

Für zwei unabhängige Kopien U,U' dieser EINEN gegebenen positiven Quelle gilt exakt

    2 Cov_mu(U,y(U))
      =E[(U-U')(y(U)-y(U'))].                           (10)

Jeder Integrand in (10) ist nichtpositiv und bei U!=U' strikt negativ. Folglich

    Cov(U,S_Y^E)<0                                      (11)

wenn irgendeine ursprüngliche Yukawa nichtzero und U unter dem tatsächlichen Integrationscaller nicht konstant ist. Bei allen Yukawas0 oder konstantem U ist die Kovarianz0. Es braucht weder ein beschränktes Restgewicht noch einen Grad0 Strahlengrenzwert, weder eine IR Entfernung noch eine kleine Kopplung oder Gaußnäherung.

Insbesondere darf eine tatsächlich vorhandene positive ungekürzte quartische H Gewichtung vollständig im selben Caller dnu stehen. Wo dieser Caller wirklich den oben geprüften konstanten ganzen H Produktpfaffian verwendet, gilt (11) auch für eine solche konfinerende Gewichtung. Der Zahlenwert aus (5) wird dadurch nicht auf die Quartik übertragen. Die ursprünglichen R048 radialen/globale Domain und Jacobiandaten bleiben erforderlich und müssen als tatsächliche positive Bestandteile von dnu vorliegen. Keine lokale radiale Gaugekarte wird globalisiert.

Die Aussage beweist damit eine strikte Kopplungsantwort des schon vorhandenen ganzen konstanten Produktmodenpaars für jede passende ursprüngliche positive H Gewichtung. Die Geltungsgrenzen aus Abschnitt4 bleiben erhalten: S_Y^E ist die globale affine euklidische Fermionmarke; kein lokales Y_ham(s), kein variabler H Modenturm, keine gekoppelte primitive OS Feldhistorie und kein komplexer/signed Referenzcaller werden ohne ihre passende ursprüngliche Identifikation eingeschlossen. Die24 symbolischen Kontrollen oben betreffen die Rangkoeffizienten. (6)–(11) sind die hier dargestellte direkte algebraische und integrale Folgerung, kein zusätzlich behaupteter Leanbeweis.

## 7. Eine gemeinsame Higgsintegration für die gesamte tatsächlich vorgegebene finite Modenfamilie

Sei J die tatsächlich vorgegebene kompatible finite Menge von m unterschiedlichen ursprünglichen Gegenmoduspaaren (p_j,-p_j), jeweils p_j^2>0, ohne Nullmodus und ohne doppelte Zählung eines Paars. Jede Paarung muss im ursprünglichen Regulator samt seinen tatsächlichen gegebenen beziehungsweise dualen Domains dieselbe oben geprüfte freie Produktform besitzen; Impulse oder Moden werden hier nicht neu gewählt. Alle Paare verwenden EIN gemeinsames konstantes ganzes H und dieselben vorgegebenen Y und MR.

Dann ist die tatsächliche gemeinsame Quelle P_total(u)=product_(j in J) Q(p_j,u)^2, nicht eine Produktintegration mit voneinander unabhängigen H_j. Mit S_j=p_j^2 Sigma_plus+T_M(p_j) und A_total=product_j aQ(p_j)^2 bleibt ihr absoluter führender Faktor erhalten: P_total=A_total u^(2kQ m)[1+4(sum_j S_j)/u+O(u^-2)].

Die ursprüngliche einmalige vierreelle Gaußintegration hat deshalb das Momentverhältnis 1/[2(2kQ m+1)]. Für denselben dominierten Grad0 Grenzwert des gemeinsam gehaltenen Restgewichts mit nichtverschwindendem führenden Winkelintegral folgt B/A=2 sum_j S_j/(2kQ m+1). Der Eulermarker g skaliert den H linearen Fermionvertex in ALLEN tatsächlichen Paaren gemeinsam, so dass derselbe bereits algebraisch verbundene Ausdruck -4B/A gilt.

Damit ist der ganze gemeinsame finite Familienwert exakt C_total=-8 sum_(j in J)[p_j^2 Sigma_plus+T_M(p_j)]/(2kQ m+1). Er ist wegen des gemeinsamen Nenners ausdrücklich nicht die Summe separat normalisierter Einpaarantworten. Bei allen Yukawas0 gilt0; das leere J gibt ebenfalls0.

Auch P_total besitzt nichtnegative Polynomkoeffizienten und einen positiven konstanten Koeffizienten. Für m>=1 und irgendein Y!=0 ist der gemeinsame Score -2u P_total'/P_total streng fallend. Deshalb gilt der strikte Vorzeichensatz aus Abschnitt6 unverändert für jede passende tatsächliche positive gemeinsame H Gewichtung mit nichtkonstantem U und endlichem ersten U Moment, ohne IR Grenzwertannahme.

Dies integriert die EIN gemeinsame H Variable der gesamten tatsächlich kompatibel vorgegebenen finiten konstanten Familie. Nichtkonstantes H, dadurch gekoppelte Impulse, der unendliche Cut Grenzwert sowie die ganze gekoppelte physische primitive OS Quelle werden nicht durch diese Faktorisation ersetzt. Ein tatsächlich vorhandener Nullmodus müsste als sein eigener ursprünglicher Quellenfaktor mitgeführt werden und ist hier ausdrücklich nicht Teil von J.
