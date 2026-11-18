# Ungekürztes ursprüngliches Radialintegral mit ausführbarer Restgrenze

5. November2026. Fortsetzung des vorhandenen R048 Callers, keine neue physische Feldfamilie. Die ursprüngliche Metrik und die übrigen Felder bleiben konditionale Argumente. Der unveränderte endliche Regulator liefert die endlichen Feldkoordinaten und ihre wirkliche Gaußpräzision. Es wird kein physischer Wert für λ,v,K oder die Koordinatenzahl N erfunden.

## Tatsächliche ursprüngliche Integrationsoperation

COMMON-CURVATURE-RESPONSE-PROOF.txt:117–174 enthält K_R, λ=8π²B/(dA²)>0, den gebrochenen Zweig v>0 und V_int(η)=√2λvη³+λη⁴/4. Die neutrale Quelle q verschiebt die tatsächlich vorhandene Quadratik zu K_q=K_R−M_q und l_q=J−gR+√2v q. Ihre ganze Gaußintegration Z_quad,u[q] wurde in SHORTTIME20261105I ausgewertet, einschließlich Determinantenquotient und ursprünglichem Z0. Setze nun ausschließlich diese vorhandenen Vertices wieder ein:

Z_rad,u[q]=Z_quad,u[q] E_(G_q,b_q)[exp(−∑_i V_i V_int(η_i))],

G_q=K_q⁻¹, b_q=G_q l_q.

Hier bedeutet V_i>0 das ursprüngliche regulierte Zellvolumen. Es ist nicht das Vakuum v. Die tatsächliche endliche Präzisionsmatrix in den gewählten Feldkoordinaten enthält diese Zellfaktoren. Ein nackter nicht selbstadjungierter Differenzialoperator darf nicht als symmetrische Gaußpräzision eingesetzt werden. Für andere ursprüngliche Basisfunktionen sind ihre tatsächlichen Auswertungsmatrizen und Integrationsdaten zu transportieren. Der vorliegende ausführbare Code erwartet die entsprechende vorhandene Zellkoordinatendarstellung.

Diese Gleichung definiert das ungekürzte endliche Reweighting der tatsächlich gegebenen R048 Quadratik mit ihren ausdrücklich bezeichneten kubischen und quartischen Vertices auf deren Gaußintegrationsbasis. Es wird keine Kopplungs Taylorreihe benutzt oder mit der Integration vertauscht. Die tatsächliche Baseline und nur die eta unabhängigen ursprünglichen Referenzfaktoren bleiben mit Z_quad,u[q] außen gebunden. H abhängige Referenzjets behalten ihre eigene ursprüngliche Innenintegration oder Subtraktion; nur tatsächlich polynomiale Beine sind durch die unten beschriebenen P Einfügungen erledigt, übrige gehören weiterhin zum ursprünglichen S_extra und Domainvertrag. Dies ist das vorhandene radiale Selbstwechselwirkungsstück; S_extra und die übrigen gemeinsamen Felder sind dadurch nicht verschwunden.

Die R048 Radialkarte ist eine lokale schwache Gaugeorbitkarte. Die Integration η∈ℝ^N der gegebenen Gaußquadratik mit den literalen V3/V4 ist dadurch nicht automatisch ein globaler Pushforward der ursprünglichen Higgsintegration. Die ganze tatsächliche Koordinatendomain, Winkel, Ghost und Jacobiandaten bleiben Teil des ursprünglichen Callers. Soweit sie tatsächlich als gelieferte polynomiale Einfügungen auf derselben Gaußbasis dargestellt sind, werden sie unten unverändert mitintegriert. Andernfalls gehören sie zum noch nicht ausgewerteten S_extra beziehungsweise zur ursprünglichen Domain und werden nicht als erledigt gezählt. Die globale Polynomuntergrenze gilt auch auf einer engeren ursprünglichen Domain. Die hier bewiesene Mittelpunkt Lipschitzgrenze für f alleine gilt jedoch nicht automatisch für einen zusätzlich eingefügten unstetigen Domainindikator. Eine solche Einschränkung muss ihre tatsächliche Domainintegration beziehungsweise passende zusätzliche Fehlerkontrolle behalten. Keine unabhängige neue Maßwahl wird verlangt oder vorgenommen.

## Exakte globale Schranke am ganzen ursprünglichen K

Schreibe a=√2v, β_i=λV_i und x_i=η_i+3a. Dann ist

W_i(η_i)=β_i V_int(η_i)/λ+27β_i a⁴/4
=β_i/4 x_i²[(x_i−4a)²+2a²]≥c_i x_i²,

c_i=β_i a²/2>0.

Der erste Ausdruck bedeutet β_i(η_i⁴/4+aη_i³)+27β_i a⁴/4; der zweite folgt aus η_i−a=x_i−4a. Setze C=∑27β_i a⁴/4=∑27λV_i v⁴. Das ungekürzte Originalgewicht ist exakt

exp(−∑V_i V_int)=exp(C) f,

f=exp(−∑W_i), 0<f≤1.

Diese Originalpolynomfaktorisierung und der daraus folgende exponentielle Bound sind Lean kernelgeprüft. C kann mit verfeinertem Regulator groß werden; keine uniforme UV oder Kontinuumsschranke wird daraus behauptet.

Mit der tatsächlichen Gaußzerlegung η=b_q+Lz, LLᵀ=G_q, z standardnormal in ℝ^N, integriert man f über den GANZEN ursprünglichen K_q. Es wird keine diagonale oder freie Feldnäherung eingesetzt.

Für k≥1 ist

M_k(c)=sup_(x∈ℝ)|x|^k exp(−c x²)=(k/(2ec))^(k/2), M_0=1.

Weil W'_i=β_i(x_i³−6a x_i²+9a²x_i), besitzt f die explizite globale Lipschitzkonstante

D_i=β_i[M_3(c_i)+6|a|M_2(c_i)+9a²M_1(c_i)],

L_f≤sqrt(λ_max(G_q)) ||D||₂.

Das benutzt die tatsächliche größte Kovarianzeigenzahl, keine willkürlich gesetzte Kopplung oder Metrik. Die ganze Quelle q kann in jedem tatsächlich ursprünglichen positiven K_q Bereich eingesetzt werden.

## Ausführbare deterministische Quadratur und Rest

Zerlege den Würfel [−T,T]^N in Würfel mit Kantenlänge h. An Mittelpunkt z_c wird f(b_q+Lz_c) exakt nach dem ungekürzten Originalpolynom ausgewertet. Jede Zelle erhält ihre tatsächliche Gaußwahrscheinlichkeit

p_c=∏_j[Φ(z_c,j+h/2)−Φ(z_c,j−h/2)].

Q=∑_c p_c f(b_q+Lz_c).

Die ganze deterministische analytische Restgrenze lautet

|E f−Q|≤L_f h sqrt(N)/2+2N exp(−T²/2).

Der erste Term folgt aus der Mittelpunktdistanz und der globalen Lipschitzkonstante; der zweite ist die Gaußtail Uniongrenze. Für eine gewünschte absolute Schranke ε kann beispielsweise T=sqrt(2 log(4N/ε)) und h≤ε/(L_f sqrtN) gewählt werden. Die wirkliche Gitterkantenlänge wird nach unten angepasst; ihre tatsächliche Punktzahl wird vor Ausführung ausgegeben. Dieser elementare Integrator ist exponentiell in N und behauptet keine schnelle Ganztheorieberechnung.

Der entsprechende ursprüngliche unnormierte Integralrest ist Z_quad,u[q] exp(C) mal dieser Rest, mit seinem tatsächlichen komplexen Baselinefaktor. Die Baseline wird nicht neu normiert.

Für q_n Quadraturen mit ε_n=ε0*2^(−n) gilt |Q_(n+1)−Q_n|≤ε_n+ε_(n+1). Diese tatsächlich konstruierbare Folge erfüllt die summierbare Defektvoraussetzung des vorhandenen Lean Satzes PositiveDescentMaster.exists_limit/tail_bound. Für f und u_i f sind auch die exakten Zellenapproximanten nichtnegativ; der positive endliche Vektorkern master_core kann dann mit diesen tatsächlich produzierten Zahlen angewandt werden. Komplexe oder signierte polynomiale Beine benutzen den Cauchysatz auf Real und Imaginärteilen, ohne künstliche Positivitätsbehauptung. Der Nutzen ist eine konvergente Auswertung des gegebenen endlichen Integrals; kein neuer Zustand wird aus einer frei gewählten Folge postuliert. In Lean formalisiert und neu kernelgeprüft ist hier die Polynomschranke. Die Gaußtail, Lipschitz und Quadraturgrenzen wurden analytisch hergeleitet und unabhängig vom Root geprüft; sie sind nicht als zusätzlich Lean bewiesene Gaußintegration ausgegeben.

## Neutrale Momente und ursprüngliche polynomiale Einfügungen

Die gleiche Quadratur erfasst die tatsächlichen numerators. Im ursprünglichen Radialchart ist

u_i=(x_i−2a)²/2.

Für jedes tatsächlich ursprüngliche endliche Polynom P(η), beispielsweise eine neutrale u Einfügung oder ein originaler Berezinquellenkoeffizient, wird P(x−3a) als ∑_p c_p x^p ausgeschrieben. Dann gilt global

B_P=∑_p |c_p|∏_i M_(p_i)(c_i)≥|P f|.

Definiere das tatsächliche Ableitungspolynom

Q_(P,j)=∂_(x_j)P−P β_j(x_j³−6a x_j²+9a²x_j).

Die gleiche monomiale Schranke gibt D_(P,j)=B_(Q_(P,j)) und

L_P≤sqrt(λ_maxG_q) ||D_P||₂.

Für die Gaußzellenquadratur Q_P dieses ganzen ursprünglichen numerators ist daher

|E(Pf)−Q_P|≤L_P h sqrtN/2+B_P 2N exp(−T²/2).

Der Integrator berechnet diese Schranken direkt aus den gelieferten ORIGINALpolynomkoeffizienten. Die Beträge dienen ausschließlich der Fehlerschranke. Das komplexe beziehungsweise signierte Originalpolynom wird unverändert in der Quadratur ausgewertet; kein Betrag eines Pfaffians wird als Gewicht eingesetzt.

Mindestens Z und alle angeforderten u_i numerators werden gemeinsam ausgewertet. Mehrere ganze Quellenpolynome können in denselben Gaußzellen mitlaufen. Ein tatsächlich fehlender ursprünglicher Polynomkoeffizient wird nicht ergänzt. Nichtpolynomiales S_extra und übrige dynamische Felder sind nicht von dieser polynomialen Einfügungsregel erledigt.

## Nenner und adaptive tatsächliche Genauigkeit

Der positive reine Radialnenner I_Z=E f liegt in [Q_Z−e_Z,Q_Z+e_Z]∩[0,1]. Ein Quotient wird erst gebildet, wenn dieser Bereich strikt von0 getrennt ist. Für ein tatsächlich ursprüngliches komplexes Nennerpolynom muss stattdessen sein Fehlerkreis |I_Z−Q_Z|≤e_Z die0 ausschließen. Dann gilt

|I_P/I_Z−Q_P/Q_Z|
≤e_P/(|Q_Z|−e_Z)+|Q_P| e_Z/(|Q_Z|(|Q_Z|−e_Z)).

Der ausführbare Pfad kann h und T bei unveränderten physikalischen Daten adaptiv verfeinern, bis diese Grenze die angegebene absolute und relative Toleranz erfüllt. Er behält den ursprünglichen Z0 Faktor und meldet eine nicht getrennte Null. Für spätere gemeinsame Feldintegration sollen ursprüngliche polynomiale Grassmann beziehungsweise Pfaffianbeine UNNORMIERT bleiben. Die Positivität der radialen Gaußfaser rechtfertigt keine Division durch einen verschwindenden oder komplexen konditionalen Pfaffiannenner. Der ganze physische Z0 Quotient entsteht erst nach den tatsächlich zugehörigen ganzen Integrationen.

Eine weitere mathematische Kontrollschranke ist E exp(−V_int,total)≥exp(−E V_int,total). Der Gaußerwartungswert auf der rechten Seite ist aus dem tatsächlichen b_q,G_q exakt polynomial bestimmt. Sie beweist Existenz und Nichtnull des reinen radialen positiven Stücks, nicht die Nichtnull eines anders gewichteten ganzen Pfaffiancallers.

## Implementierung, Nachweis und noch fehlende ursprüngliche Eingabe

integrate_original_radial_tilt.py implementiert die unveränderte ganze Kovarianz, deterministische normalgewichtete Zellen, sämtliche polynomialen Einfügungsgrenzen, positive oder komplexe Nennerprüfung und adaptive Toleranz. Er hat keine Defaultwerte für physisches K, λ, v, V_i oder N. Er verlangt die tatsächliche originale Präzisionsmatrix, den ursprünglichen Gaußmittelwert, die ursprünglichen Zellvolumina sowie λ und v. Die Rechenzahlbegrenzung und Toleranzen sind reine numerische Algorithmenparameter. Eine nicht ausführbare Punktzahl wird gemeldet, bevor eine Integration beginnt.

Der Code benutzt aktuell gewöhnliche Floats, Cholesky/eigvalsh und ndtr. Seine analytische Tail und Quadraturfehlergrenze ist von Rundungsfehlern getrennt. Ausgegebene Endpunkte sind ausdrücklich VOR einer Rundungseinschließung; rounding_enclosed=false. Ein vollständig gerichtetes numerisches Intervall ist damit nicht zertifiziert. Die tragende Polynomschranke ist kernelgeprüft. Root hat den ausführbaren Integrator zusätzlich in sechs mathematischen Methodenkontrollen unabhängig geprüft: Z und u gegen45stellige direkte Quadratur, korrelierte ganze Kovarianz gegen Hermite48/80, erhaltener exakter Nullnenner und adaptive u Toleranz bei unveränderten Kontrollparametern. Alle sechs bestehen, einschließlich der Zurückweisung einer angeforderten, nicht gelieferten UY Marke. Der adaptive Lauf verwendete vier Verfeinerungen und295532 finale Zellen; die Abweichung zur unabhängigen Quadratur8.63e−7 liegt unter der analytischen Grenze vor Rundung0.003019. Diese Kontrollparameter sind ausdrücklich keine tatsächlichen nativen TFPT Felddaten und keine zusätzlichen physischen Quellwahlen. Ergebnisse und Geltungsbereich stehen in work/shorttime20261105J_gravity/integrator-math-controls.json und ROOT-NUMERICAL-METHOD-CONTROLS.txt. Ohne tatsächliche native Originaldaten wurde keine physische Integralzahl oder primitive Gesamtintegration berechnet.

In der tatsächlich gelesenen R048 Quelle ist K_R als Operator mit gegebenem Hintergrund und Regulator formuliert; eine konkrete physische endliche K_R Präzisionsmatrix und Zellliste sind dort nicht numerisch eingesetzt. Die erste benötigte Einfügung für eine tatsächliche native Ausführung ist jetzt genau diese vorhandene Originalregulatordarstellung samt den ganzen tatsächlichen Koordinaten, Domain, Ghost und Jacobiandaten, nicht die Auswahl eines neuen Feldmaßes. Das quartische Reweighting des gegebenen endlichen Gaußcallers ist konstruktiv und mit expliziter numerischer Restformel bestimmt. Seine Gleichsetzung mit der globalen ursprünglichen Higgsgaugeintegration folgt nicht allein aus der lokalen Radialkarte. Ganzer ursprünglicher UY Wert, gekoppeltes R Integral, zusätzliche Domainantworten, S_extra, UV/Kontinuum und vollständiger physischer Kernel sind dadurch nicht behauptet.

Originalbeleg: COMMON-CURVATURE-RESPONSE-PROOF.txt:117–174, SHA2562b05629937b9d46883d4c456c5bb30acecdc59e7876767eb2207d33026683df2. Ganze H/Fermionquelltypen und S_extra: DYNAMICAL-HIGGS-SOURCE-PROOF.txt:34–78,88–144, SHA256eda5bb664b1fb20370e59b39dbd5e3cb04d79864a8d1d0a0eb6c791cdaf98f49. Originalradial und Hilfsquellenidentität: WHOLE-NEUTRAL-R048-SOURCE-AND-FIRST-FULL-OBSTRUCTION.md aus SHORTTIME20261105I. Lean Beleg und Axiome: SourceAxiomAudit.lean/.log; gesamte aktuelle Leansuche: LEAN-ORIGINAL-AUDIT-AND-CONSTRUCTIVE-LEVER.md.

Der unabhängige Review fand einen behobenen Implementierungsfehler: nicht gelieferte angeforderte normalize_insertions Namen werden jetzt vor Planung und Quadratur abgewiesen. Insbesondere darf eine fehlende UY Einfügung keine leere erfolgreiche Toleranzprüfung erzeugen. Die zusätzliche gezielte Regression und alle übrigen Methodenkontrollen wurden nach der Korrektur ausgeführt; dies erzeugt keinen neuen physikalischen Quellenwert.
