# TFPT F1: Quellenphase, originale Wardidentitäten und verbleibende Herkunftsfrage

**Stand: 23. November 2026. Fortsetzung unter Einbeziehung des Anhangs „Eingefügter Text(20261123-070045).txt“.**

Geprüft wurde der erneut abgefragte Originalstand von [sthamann/topological-fixed-point-theory, Commit 5b3f26d7befd7550127fe5d60fa94c899ed81ed3](https://github.com/sthamann/topological-fixed-point-theory/commit/5b3f26d7befd7550127fe5d60fa94c899ed81ed3). Frühere Ergebnisse dieser Untersuchung werden als Vorbefunde behandelt. Keine neue Änderung am Originalrepository oder Atlas wurde vorgenommen.

## Ergebnis

Die vollständige F1 Herkunft ist weiterhin offen. Diese Fortsetzung schließt aber eine konkrete Frage des Anhangs positiv und liefert zwei stärkere Auswahlmechanismen:

1. **Der tatsächliche ursprüngliche Tensor W mit 480 Kopplungen erfüllt sämtliche nichtabelschen Standardmodell Wardidentitäten exakt.** Damit ist die lokale Verträglichkeit dieses W mit der Grassmannverbindung bewiesen.
2. **Eine ausdrücklich definierte positive fundamentale Determinantenquelle wählt global 3+2 über ein offenes Parameterintervall.** Die vollständige logarithmische Wirkung genügt; ein quartischer Abbruch ist nicht erforderlich.
3. **Auch die vollständige E8 Determinantenquelle besitzt eine lokal stabile Phase in der gewünschten Richtung.** Ein ganzes Intervall ist exakt zertifiziert. Das geht über den schon im Oktober dokumentierten Sattelbefund der isolierten achten Ordnung hinaus.
4. **Die Polarisation bestimmt einen lokalen Diracoperator unter genau benannten Zusatzdaten eindeutig.** Der pseudodifferentielle Gegenbeweis des Anhangs widerlegt diesen lokalen Starrheitssatz nicht.

Der wichtigste verbleibende Schritt ist nun präziser: P1 muss das tatsächlich verwendete Quellenmaß, die Darstellung und den Kopplungskanal bestimmen und dieselbe Quelle muss ihre lokalen Felder, Zeit und höheren Antworten hervorbringen. Die folgenden Resultate belegen Teile dieser Verbindung, setzen andere Teile aber noch voraus.

## 1. Der Originaltensor W besteht den geforderten Test

Der Anhang schlägt vor, die Hyperladungsprüfung um die vollständigen nichtabelschen Wardgleichungen zu ergänzen. Dafür wurde das ursprüngliche Archiv `experiments/theory-contracts/rr-residue-memory-20260921/native_tensor.npz` benutzt, keine passend neu konstruierte Ersatzkopplung.

Sein SHA256 lautet:

```
3f00a0892157a691e4ef26920c3702772c7bf76a250fad98ff04a079bb64b763
```

W ist eine 60 mal 2016 Matrix mit 480 Einträgen ±1. Der unabhängige Aufbau aus den ursprünglichen Clifford und Exteriorräumen reproduziert alle Einträge; außerdem gilt WW†=8I₆₀. Für jeden Generator wurde die gesamte Identität

\[
W\,d\Lambda^2\rho_{\mathcal F}(T)
=\rho_{\mathcal B}(T)W
\]

an allen 120960 Matrixpositionen geprüft.

| Familie | Zahl der Generatoren | Ergebnis |
|---|---:|---|
| SU(3) Farbe | 8 | exakt null Residuum |
| SU(2) schwach | 3 | exakt null Residuum |
| Hyperladung | 1 | exakt null Residuum |
| U(5) mit ursprünglichem Determinanttwist | 25 | exakt null Residuum |
| Spin(10) | 45 | exakt null Residuum |
| SU(4) Familie | 15 | exakt null Residuum |

Diese Familien überlappen; die Zahlen werden nicht zu einer Zahl unabhängiger Symmetrien addiert. Die U(5) Zielwirkung wurde unabhängig aus der natürlichen Darstellung samt Determinanttwist vorgegeben. Für die zusätzlich getesteten Spin(10) Generatoren wurde die induzierte Zehnerwirkung bestimmt und ihre Leckagefreiheit, Faktorisierung und Erhaltung der symmetrischen Paarung gesondert geprüft.

Eine negative Kontrolle verdoppelt nur einen vorhandenen Koeffizienten. Alle Ladungslabels bleiben erhalten, und der Hyperladungstest besteht weiterhin. Sechs nichtdiagonale Farbidentitäten scheitern. Die Prüfung unterscheidet also wirklich vollständige Kovarianz von bloßer Ladungserhaltung.

### Was damit geschlossen ist

Für jede Verbindung A in der getesteten Algebra folgt

\[
\boxed{dW+A_{\mathcal B}W-WA_{\Lambda^2\mathcal F}=0.}
\]

Insbesondere ist die kanonische Grassmannverbindung lokal mit dem Originaltensor W kompatibel. Die kleinere Symmetrie der quadratischen Compilerabbildung μ ist hierfür kein Hindernis: μ und W sind unterschiedliche Abbildungen zwischen unterschiedlichen Räumen.

W†W liefert entsprechend einen invariant kontrahierten internen quartischen Kernel. Eine physische verbundene Vierpunktfunktion benötigt darüber hinaus passende Propagatoren, einen Zustand und Dynamik. Diese folgen nicht aus der Tensoridentität allein.

Auch der früher geprüfte sichtbare Fünferkanal, hier zur Unterscheidung \(\mathcal E\) genannt, ist ein anderes Objekt. Sein kontinuierlicher Kovarianzkern bleibt skalar. Der positive W Befund repariert seine fehlende Kovarianz nicht automatisch. Der gemeinsame Quellenanschluss von W, \(\mathcal E\) und lokalen Feldern bleibt zu bestimmen.

Belege und vollständige Rechnung: `gauge_geometry/REPORT.md`, `check_native_W.py`, `native_W_results.json` und die unveränderten Originalquellen.

## 2. Die Auswahl von 3+2 lässt sich ohne quartischen Abbruch beweisen

Für eine spurfreie hermitesche Matrix Σ auf C⁵ sei pₖ=trΣᵏ. Die Rechnung des Anhangs ist im wesentlichen richtig:

\[
V_0(\Sigma)=\operatorname{tr}(\Sigma^2-I)^2
\]

hat als globale Minima

\[
\Sigma=\pm\frac1{\sqrt7}
U\operatorname{diag}(2,2,2,-3,-3)U^\dagger,
\qquad V_{\min}=5/7.
\]

Die vollständige stationäre Liste enthält zusätzlich den im Anhang fehlenden Typ diag(−1,0,0,0,1) mit V=3. Er ist ein Sattel und verändert das globale Minimum nicht. Die Minimalmenge besteht aus **zwei** getrennten Grassmannorbits. Die Symmetrie Σ→−Σ ist noch keine Identifikation der beiden Komponenten als Eichäquivalenz.

Die scharfe Schranke lautet

\[
\boxed{\frac7{30}\leq\frac{p_4}{p_2^2}\leq\frac{13}{20}.}
\]

Die untere Gleichheit kennzeichnet 3+2, die obere 4+1. Für die übliche quartische Familie −ap₂+bp₄+cp₂² ergibt sich bei a,b>0 und c+7b/30>0 die robuste 3+2 Auswahl. Diese SU(5) Potentialfamilie und ihre Symmetriebrechung sind bekannte Physik; siehe [Pogosian und Vachaspati, Domain Walls in SU(5), Gleichungen 2 bis 6](https://arxiv.org/pdf/hep-ph/0007045). Der neue Anschluss besteht in der vollständigen positiven Quelle und ihrem konkreten Verhältnis zum ursprünglichen Spiegelungsoperator.

### Ein vollständiger positiver Quellenzeuge

Betrachte das endliche Maß

\[
\boxed{
d\mu_{\rho_5}(\Sigma)=Z^{-1}
e^{-\rho_5 p_2}\det(I+\Sigma^2)\,d\Sigma,
\qquad\rho_5>0.
}
\]

Es ist positiv und integrierbar. Seine vollständige logarithmische Wirkung ist

\[
\Gamma_{\rho_5}(\Sigma)
=\rho_5p_2-\operatorname{tr}\log(I+\Sigma^2).
\]

**Für 1/3<ρ₅<1 bestehen sämtliche globalen Minima genau aus den beiden 3+2 Orbits**

\[
\Sigma_{\min}\sim\pm v\operatorname{diag}(2,2,2,-3,-3),
\]

\[
\boxed{
v^2=\frac{6-13\rho_5+\sqrt{25\rho_5^2-12\rho_5+36}}{72\rho_5}.
}
\]

Der Beweis klassifiziert alle stationären Eigenwerttypen der vollständigen logarithmischen Wirkung. Dreiwertige Konkurrenten und 4+1 haben im behaupteten Intervall negative Hessianrichtungen; der Ursprung ist ein Maximum. Koerzivität garantiert ein globales Minimum. Damit bleibt ausschließlich 3+2. Der vollständige Beweis und die exakten Hesseformen stehen im Teilbericht `quartic_source/REPORT.md`.

Dies behebt ein reales Problem des bloßen quartischen Ansatzes: Ohne konfinerenden Zusatz fällt die negative logarithmische Determinantenwirkung bei großen Feldern unbegrenzt. Ein Minimum des abgeschnittenen Polynoms würde dann keinen stabilen Zustand der vollständigen Quelle belegen. Das hier verwendete positive quadratische Maß beseitigt diese Instabilität ausdrücklich.

### Konkreter Anschluss an die originale Spiegelung

Im gelesenen P1 Original gilt für die involutive geometrische Spiegelung τ

\[
D_{\mathrm{rel}}=\tfrac12(D-\tau D\tau),
\qquad\{D_{\mathrm{rel}},\tau\}=0.
\]

Die neu vorgeschlagene Kopplung

\[
H_\Sigma=D_{\mathrm{rel}}\otimes I_5+g\tau\otimes\Sigma
\]

besitzt auf einem konstanten internen Feld und einem paarungserhaltend abgeschnittenen Bereich die exakte Eigenschaft

\[
H_\Sigma^2=D_{\mathrm{rel}}^2\otimes I_5+g^2I\otimes\Sigma^2.
\]

Das ursprüngliche negative logarithmische Profil kann damit das positive gepaarte Determinantengewicht liefern. **Die Kopplung τ⊗Σ, die fundamentale Darstellung und das Matrixmaß sind jedoch noch neue Voraussetzungen.** Die originale Spiegelung allein wählt diesen Kanal nicht aus. Eine kommutierende Kopplung I⊗Σ gibt andere Vorzeichen; eine feldabhängige Referenzsubtraktion oder der im Original zusätzlich geschriebene positive logarithmische Regulatorterm kann die Gesamtwirkung verändern.

Die spektrale Paarung von euklidischen Fermionen und das Determinantenargument sind Standardwerkzeuge; siehe [Kaplan, Chiral Symmetry and Lattice Fermions, Abschnitt 1.8](https://arxiv.org/pdf/0912.2560) sowie [Dunne, Functional Determinants in Quantum Field Theory](https://saalburg.aei.mpg.de/wp-content/uploads/sites/25/2017/03/dunne.pdf). Die Anwendung auf diese konkrete TFPT Kopplung ist hier ein bedingter Vorschlag.

Ein endliches invariantes Maß hat weiterhin ⟨Σ⟩=0. Die Minimumorbits sind Aussagen über die Wirkung. Eine Konzentrationsgrenze e^(−NΓ), die Auswahl einer Orientierung und ein physischer gebrochener Zustand erfordern zusätzliche Schritte.

## 3. Die Darstellung entscheidet: 5, 16 und E8 sind nicht austauschbar

Für denselben positiven Determinantenkoeffizienten wurden die vollständigen Darstellungsspuren berechnet:

| Darstellung | Tr_R Σ² | Tr_R Σ⁴ | Richtungspräferenz der quartischen Näherung |
|---|---|---|---|
| 5 | p₂ | p₄ | 3+2 |
| 10 | 3p₂ | 3p₂²−3p₄ | 4+1 |
| 16=Λ^gerade5 | 4p₂ | 3p₂²−2p₄ | 4+1 |
| adj SU(5) | 10p₂ | 6p₂²+10p₄ | 3+2 |
| adj E8 | 60p₂ | 36p₂² | keine Formauswahl |

Für E8 ist auch Tr₂₄₈Σ⁶=30p₂³ rein radial. Insbesondere darf ein Erfolg im fundamentalen Fünferraum nicht ohne neue Rechnung der tatsächlichen Sechzehnerdarstellung zugeschrieben werden. Ebenso ist die vollständige E8 Antwort keine Summe ausschließlich günstig wirkender Fünferbeiträge: Die anderen Darstellungen heben deren quartische Formabhängigkeit exakt auf.

Das Fehlen primitiver E8 Invarianten vierten und sechsten Grades ist bekannt; siehe [Cederwall und Palmkvist, The octic E8 invariant](https://arxiv.org/pdf/hep-th/0702024). Die obigen Identitäten wurden unabhängig mit den ursprünglichen Wurzeln als vollständige Polynome geprüft.

### Warum die vollständige E8 Quelle trotzdem interessant bleibt

Der ursprüngliche Versuch einer Auswahl durch den isolierten Term achten Grades wurde bereits am 13. Oktober 2026 archiviert. Die Richtung

\[
h=\frac1{\sqrt{30}}(2,2,2,-3,-3,0,0,0)
\]

ist dort ein Sattel. Diese bekannte Sackgasse wird nicht als neue Untersuchung ausgegeben.

Die vollständige Wirkung

\[
\boxed{
\Gamma_{\rho_8}(H)=\rho_8\|H\|^2
-\sum_{\alpha\in\Delta(E_8)}
\log\bigl(1+(\alpha\cdot H)^2\bigr)
}
\]

verhält sich bei endlicher Kopplung anders. Setze H=√q h und

\[
\rho_8(q)=\sum_\alpha
\frac{(\alpha\cdot h)^2}{1+q(\alpha\cdot h)^2}.
\]

Dann ist H radial stationär. **Für jedes 10≤q≤15 ist der entsprechende E8 Orbit ein striktes lokales Minimum modulo Konjugation.** Die drei normalen Hessianfunktionen sind auf diesem ganzen Intervall durch positive Bernsteinzertifikate exakt kontrolliert. Auch die radiale Krümmung −4qρ₈′(q) ist positiv.

Der dazugehörige Parameterbereich ist

\[
\frac{52660}{5643}\leq\rho_8\leq\frac{3075}{247}
\quad\approx[9.33191565,12.44939271].
\]

Bei q=10 sind die kartesischen normalen Krümmungen ungefähr 3.181899, 0.408990 und 2.659621 mit Multiplizitäten 8,3,24; radial ergibt sich 34.081245. Die 212 übrigen Richtungen sind die Nullrichtungen des E8 Konjugationsorbits. Insgesamt sind somit alle 36 normalen Richtungen positiv. Eine unabhängige Gegenrechnung mit selbst erzeugten Wurzeln bestätigt sowohl die vollständige Cartanhessian als auch den Übertrag auf die adjungierte Algebra.

Die zugehörige Grafik und sämtliche rationalen Funktionen liegen im Teilordner `representation_source`. Die Hauptrechnung und die unabhängige Gegenprüfung unterscheiden ausdrücklich Hessianwerte auf der Einheitsrichtungssphäre von kartesischen Werten; letztere enthalten den Faktor 1/q.

Eine ergänzende Suche mit je 41 Starts findet bei q=10 und q=15 keinen tieferen Endwert. Das ist Suchinformation, **kein globaler Beweis**. Bei q=20 bleibt die Richtung noch lokal stabil, aber die Suche findet bereits eine bessere Konkurrenzrichtung. Gerade dieses Beispiel zeigt die notwendige Trennung zwischen lokaler Stabilität und globaler Auswahl.

### Die E8 Symmetrie ist dabei größer als das Standardmodell

Für das feste Y enthält E8 28 neutrale Wurzeln und acht Cartanrichtungen. Der Zentralisator ist auf Liealgebraebene

\[
\mathfrak{su}(5)\oplus\mathfrak{su}(3)
\oplus\mathfrak{su}(2)\oplus\mathfrak u(1),
\]

also 36 dimensional. Der sichtbare SU(5) Teil liefert die vertrauten zwölf neutralen Standardmodellgeneratoren. Der vollständige E8 Zentralisator enthält zusätzliche Symmetrien. Außerdem bezeichnet 3+2 die Richtung relativ zur gewählten sichtbaren Einbettung. Diese Markierung und die weitere Reduktion werden durch die neue Potentialrechnung nicht hergeleitet.

Der neue positive Befund lautet daher präzise: **Ein vollständiges logarithmisches Quellenmodell kann die gewünschte Richtung über einen Bereich stabilisieren, obwohl seine niedrigen Taylorordnungen das nicht leisten.** Die Identifikation dieses Modells mit der primitiven Quelle, seine globale Auswahl und die physische Interpretation bleiben offen.

## 4. Grassmanngeometrie: richtiger Symmetrieraum, noch eingeschränkte Felder

Die 66 Krümmungspaare der zwölf reellen Tangentialrichtungen von Gr₃(C⁵) spannen exakt

\[
\mathfrak s(\mathfrak u(3)\oplus\mathfrak u(2)).
\]

Dieser positive Befund des Anhangs ist bestätigt. Durch den nun geprüften W ist zudem ein kompatibler ursprünglicher interner Kopplungstensor vorhanden.

Trotzdem entstehen aus einer Projektorbewegung allein nicht alle benötigten unabhängigen Eichfeldkonfigurationen. Zwei präzise Rechnungen lokalisieren die Grenze:

**Erstens:** Um einen konstanten Projektor beginnt die Berryverbindung bei zweiter Ordnung in der Störung X:

\[
A=\tfrac12(X^\dagger dX-dX^\dagger X)+O(X^4),
\qquad F=dX^\dagger\wedge dX+O(X^4).
\]

Eine daraus gebildete Wirkung tr(F∧*F) beginnt deshalb quartisch. Der gewöhnliche quadratische kinetische Term unabhängiger Yang und Mills Felder folgt hier nicht.

**Zweitens:** Für die komplementären Berrykrümmungen im flachen festen C⁵ gilt punktweise

\[
\boxed{
\operatorname{tr}(F_3\wedge F_3)
+\operatorname{tr}(F_2\wedge F_2)=0.
}
\]

Eine reine Farbverbindung mit nichtverschwindender entsprechender Vierform und F₂=0 ist damit bereits lokal ausgeschlossen. Volle Holonomie auf dem Parameterraum ist stärker als bloße Labels, aber schwächer als die Freiheit beliebiger physischer Eichfelder.

### Was sich konstruktiv reparieren lässt

Auf S⁴ erzwingt ein triviales E₃⊕E₂ die Relation k₃+k₂=0. Eine unabhängig gewählte Verbindung auf denselben Bündeln ändert deren Topologie nicht. Für beliebige vorgegebene k₃,k₂ genügt hingegen ein Rangzweikomplement K₂ mit c₂(K₂)=−(k₃+k₂):

\[
E_3\oplus E_2\oplus K_2\cong\underline{\mathbb C}^{7}.
\]

C⁶ genügt bei nichtverschwindender Summe nicht, weil ein Linienkomplement auf S⁴ keine zweite Chernklasse tragen kann. **Sieben ist hier die minimale topologische Hülle auf S⁴**, keine abgeleitete neue Teilchenzahl und keine universelle minimale Dimension für die Darstellung jeder Verbindung.

Die Darstellung einer gegebenen Verbindung als Berryverbindung in genügend großem Raum ist möglich. Das ist das klassische Resultat über universelle Verbindungen von [Narasimhan und Ramanan](https://dml.cz/bitstream/handle/10338.dmlcz/700905/Toposym_01-1961-1_68.pdf), mit lokalen Dimensionsfragen etwa bei [Schlafly](https://msp.org/pjm/1982/98-1/pjm-v98-n1-p15-s.pdf). Im Prüfpaket ist zusätzlich eine ausdrückliche lokale Konstruktion A=U†dU enthalten. Ihr erster frei eingesetzter Input ist A selbst. Sie repräsentiert eine Verbindung und leitet deren Dynamik nicht her. Die topologische Klassifikation auf S⁴ folgt aus der Clutchingbeschreibung; siehe [Hatcher, Vector Bundles and K Theory, Abschnitt 1.2](https://pi.math.cornell.edu/~hatcher/VBKT/VB.pdf).

## 5. Die Zeitmehrdeutigkeit besitzt eine konkrete lokale Grenze

Die Familie des Anhangs

\[
B_\varepsilon=B_0+\varepsilon B_0(1+B_0^2)^{-1}
\]

behält Polarisation, Hauptsymbol, geometrische Chiralität und die vorhandene Vierteldrehung, verändert aber Lücke und Zeitentwicklung. Das ist ein gültiger Gegenbeweis gegen Eindeutigkeit in der größeren pseudodifferentiellen Operatorenklasse.

Der Zusatz verletzt jedoch die Bedingung erster Ordnung

\[
[[B,a],b]=0\quad(a,b\in C^\infty(T^2)).
\]

Ein konkreter Fourierzeuge wurde exakt berechnet. Diese Bedingung steht unter anderem im spektralen Rekonstruktionsansatz von [Connes](https://arxiv.org/pdf/0810.2088). Die weiteren Voraussetzungen dieses Rekonstruktionstheorems werden hier nicht als für TFPT bewiesen angenommen.

**Eigener lokaler Starrheitssatz:** Fixiert seien der flache antiperiodische Torus, die Darstellung der Funktionen, das Dirac Hauptsymbol, der Cliffordrahmen und die geometrische Chiralität Γ=σ₃⊗I. Dann erzwingt für einen lokalen Operator D=D₀+V(x) die wörtliche Gleichheit signD=signD₀ bereits **D=D₀**.

Der Beweis hat zwei Schritte. Die Kommutatorgleichung mit der Polarisation zwingt V zunächst in die Form I₂⊗Φ(x). Die Fouriergleichung zwingt Φ zur Konstanz. Die feste geometrische Chiralität setzt schließlich Φ=0. Ohne den letzten Schritt bleiben kleine konstante innere Potentiale innerhalb der Lücke möglich.

Der lokale Verbindungsansatz des Anhangs ist damit keine stärkere Gegenprobe mit derselben vollständigen Polarisation: Er verändert deren positive und negative Unterräume bereits in erster Ordnung. Er zeigt eine Mehrdeutigkeit bei weniger festgehaltenen Daten.

Das ist ein konstruktiver Fortschritt, aber keine Ableitung von Metrik, Hauptsymbol, Spinstruktur oder Raumzeit. Bei einer anderen Graduierung mit innerem Faktor bleiben wiederum lokale Gegenbeispiele möglich. Die genauen Voraussetzungen und Gegenbeispiele stehen in `spectral_clock/README.md`.

## 6. A2 Dreieck und 8/24 Rückkehr verlangen unterschiedliche Nachweise

Die vollständige Wurzelzählung bestätigt das Dreieck

\[
e_0-e_2,\quad e_2-e_4,\quad e_4-e_0
\]

mit Ladungen 0,−5/6,+5/6. In E8 existieren genau sechs Wurzeln jeder Ladung ±5/6; alle zwölf sind die verbindenden SU(5) Wurzeln. Das orthogonale Wurzelsystem des Dreiecks ist E6. Diese Struktur ist real vorhanden.

Die tatsächlichen Quellenmatrizen liefern jedoch verschiedene Rückkehrbegriffe:

| Konkret geprüfter Operator | Lineare Ordnung | Ordnung auf Strahlen |
|---|---:|---:|
| Nativer c=J₀σ | 12 | 6 |
| Ursprüngliche Faltung C(±i), C⁸=±iI | 32 | 8 |
| Bedingte Erweiterung C(±i)⊗P₃ | 96 | 24 |

Beim nativen c gilt c⁸=σ². Abtastung bei 0,8,16,24 Schritten besucht also das Dreieck, obwohl die erste vollständige Rückkehr unter c schon bei zwölf liegt. Die konkrete Matrix J₀ dreht vier reelle Ebenen gleichzeitig. Ihr Spin(8) Lift hat daher Ordnung vier, während eine einzelne geometrische Spin Vierteldrehung Ordnung acht hat. Das ist noch keine Identifikation des ursprünglichen 16 Majorana Focklifts oder des Raumzeit Spinlifts.

Der abstrakte Satz bleibt uneingeschränkt richtig: Haben unabhängige Operatoren echte lineare Ordnungen acht und drei, hat ihr Tensorprodukt Ordnung 24. Wegen der Teilerfremdheit verkürzen zentrale Phasen diese Ordnung nicht. Offen ist die Identifizierung genau solcher Faktoren in derselben ursprünglichen Quelle.

Außerdem erhält die native σ das festgehaltene Y auf 198 von 240 Wurzeln nicht; c verändert es auf 224. Ein aktiver physischer Übergang benötigt damit Ladungsüberträge oder mitbewegte Marken. Ein Wurzelorbit allein ist noch kein ladungserhaltender physischer Takt.

## 7. Was vom Schlussbild des Anhangs bestehen bleibt

Die Folge M₂ → quadratische Fünferkonstruktion → 3+2 → Grassmanngeometrie → 16 → Spin(10)/E8 verbindet mathematisch sinnvolle Objekte. Diese Untersuchung liefert nun einen echten lokalen Anschluss zwischen der Grassmannverbindung und dem ursprünglichen W sowie kontrollierte Potentialmechanismen. Sie identifiziert aber noch nicht alle Pfeile mit ein und derselben primitiven Quelle.

Insbesondere liefert der Lorentzkegel in Herm₂(C) noch keine lokale vierdimensionale Raumzeit mit Feldgleichungen. Die Verknüpfung von Skalen und Zeitentwicklung durch denselben Spektraloperator liefert nicht automatisch eine weitere unabhängige Koordinate. In dieser Runde wurde kein fehlender Raumzeit oder Gravitationsbeweis ersetzt.

Die folgende Übersicht nennt für jede ernsthaft geprüfte Fortsetzung den ersten noch eingesetzten Input. Alle Einträge sind **bedingte Anschlüsse**, keine unabhängigen Fundamentalursprünge:

| Ansatz | Erreichter Fortschritt | Erster noch freier Input |
|---|---|---|
| Reflektierter fundamentaler Fermionkern | positive vollständige Quelle, globales 3+2 Minimum | fundamentale Kopplung τ⊗Σ und Matrixmaß |
| Vollständige E8 Determinante | exaktes stabiles Intervall, 36 positive Normalrichtungen | E8 Kopplung, Maß und ρ₈ |
| Originales W plus Grassmannverbindung | sämtliche Wardidentitäten und ∇W=0 | physische Verbindung, Propagatoren und Zustand |
| Größerer universeller Verbindungsraum | konstruktive Darstellung beliebiger gegebener lokaler Verbindungen | die Verbindung A selbst |
| Topologische Ergänzung auf S⁴ | beliebige k₃,k₂ in trivialem C⁷ | Komplementbündel und physische Zuordnung |
| Lokale spektrale Starrheit | eindeutiger D bei identischer Polarisation | feste Metrik, Hauptsymbol und geometrische Chiralität |
| Native A2 und Uhranalyse | exakte Ladungen, tatsächliche Ordnungen und Phasen | gemeinsamer Quelloperator samt Markentransport |

## 8. Der nächste entscheidende Ursprungstest

Die engste offene Frage ist kein weiteres Zählen von Symmetriegeneratoren. Der Originalprozess muss selbst bestimmen, welches Feld in welcher Darstellung wie an D_rel koppelt und welches Quellenmaß verwendet wird. Dafür sind vier überprüfbare Aussagen nötig:

1. **Darstellung und Kopplung aus P1:** Eine konkrete Originaloperation muss den verwendeten internen Operator liefern und zwischen dem kommutierenden und antikommutierenden Kanal entscheiden. Die Existenz eines Fünferraums allein reicht nicht.
2. **Gesamte relative Wirkung:** Primitive Determinante, Referenz und Regulatorterm müssen gemeinsam ausgewertet werden. Erst daraus dürfen ρ und die stabilisierenden Antworten gelesen werden.
3. **Gemeinsame lokale Dynamik:** Derselbe Anschluss muss den positiven W Befund mit einer kovarianten Dynamik verbinden. Ein eigenständiger Yang und Mills Term wäre eine mögliche bedingte Vervollständigung, aber eine neue Annahme, solange er nicht aus der Quelle folgt.
4. **Physischer Rahmen:** Der ursprüngliche Zustand, lokale Raumzeit, korrekte Zeitrekonstruktion und die gemeinsame Transportregel für alle Marken müssen festgelegt werden, ohne die später gewünschten Symmetrien als Auswahlkriterium einzusetzen.

Die neuen Resultate schaffen eine belastbare Ausgangslage für diese Tests. Sie ersetzen sie nicht. Insbesondere gelten **globales Minimum des fundamentalen Modells**, **lokales Minimum des vollständigen E8 Modells** und **exakte Kovarianz des ursprünglichen W** mit jeweils unterschiedlichen Voraussetzungen.

## 9. Reproduktion

Das Prüfpaket enthält die fünf Hauptprüfer dieser Runde und eine unabhängige E8 Gegenprüfung. Der ausführende Runner arbeitet in einer temporären Kopie; die archivierten Resultate und Quellhashes bleiben erhalten. Er benötigt keinen Netzwerkzugang.

Die drei abgegrenzten symbolischen Hauptblöcke enthalten 80 Prüfungen zur quartischen und vollständigen fundamentalen Quelle, 28 zur darstellungsabhängigen E8 Quelle und 797 zur Spektral und Uhranalyse. Hinzu kommen 13 unabhängige E8 Prüfungen sowie die gesonderten vollständigen W und Geometrieprüfer. Diese heterogenen Aussagen werden nicht zu einer vermeintlichen Zahl physischer Beweise addiert.

Die W Identitäten und algebraischen Zertifikate sind exakt. Die lokale universelle Verbindung wird zusätzlich numerisch geprüft; die konkurrierenden E8 Minima werden ausdrücklich nur explorativ gesucht. Das Paket weist diese Ebenen getrennt aus.

Enthalten sind außerdem die Originaleingaben, der neue Anhang, die verwendeten P1 und Kernquellen, der bereits bekannte oktische Befund, die detaillierten Teilberichte, die exakten Ergebnisdateien, eine Stabilitätsgrafik und ein SHA256 Manifest. Die vollständigen Teilbeweise stehen in den Fachberichten; dieser Bericht hält die gemeinsame Schlussfolgerung und ihre Voraussetzungen fest.

**Abschließender Status:** Der ursprüngliche nichtabelsche W Anschluss ist positiv geprüft. Zwei vollständige Quellenmodelle besitzen kontrollierte 3+2 Phasen mit unterschiedlichen Stärken. Die gemeinsame Herkunft dieser Quellen, ihrer lokalen Dynamik und der Raumzeit aus P1 ist weiterhin die offene F1 Aufgabe.
