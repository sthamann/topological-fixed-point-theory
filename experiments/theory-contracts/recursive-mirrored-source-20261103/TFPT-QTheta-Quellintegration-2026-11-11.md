# Der ursprüngliche QTheta Aufruf und der erste verborgene Feldanteil

11. November 2026, MEZ

**Der installierte QTheta Aufruf integriert keine physische Feldverteilung.** Die grobe Matrix, die ihm bisher als möglicher Einstieg zur Verfügung steht, kann diese Feldverteilung auch nicht bestimmen. Ich habe das am unveränderten Quellassembler mit einer verschachtelten räumlichen Verfeinerung geprüft. Zwei Felder liefern im installierten Block exakt dieselbe Matrix. Im feineren Block tragen sie unterschiedliche spektrale Gewichte und einen unterschiedlichen Rückweg durch zwei der groben Fourierplätze. Aus dem ersten offenen Übergang in diese zusätzliche Mode lässt sich unter einer genau angegebenen Voraussetzung die Verteilung des verborgenen Freiheitsgrads zurückgewinnen.

## Randpaarung, Zeit und Leerreferenz im ursprünglichen Aufruf

Die [QFT Definition des Selektors](/Users/stefanhamann/Projekte/tfpt-theoryv4/documents/archive--tfpt-45--source_extracts--04_qft_source/sections/002-section-closure-theorems-from-the-admissibility-complex.tex:78) gibt den gemeinsamen äußeren Projektor P aus Nahtzeichen, Farbzentrum, Determinanteninvolution und Spektralauswahl. Es gilt QΘ P = P. Der [ursprüngliche OS Abschnitt](/Users/stefanhamann/Projekte/tfpt-theoryv4/documents/archive--tfpt-45--source_extracts--04_qft_source/sections/008-section-relative-reflection-positivity-with-fermions.tex:235) gibt die reflektierte Gramform und aus der Zeitverschiebung auf deren Nullquotient eine positive Kontraktionshalbgruppe. Diese Konstruktion gibt die physische Zeit an, sobald die tatsächliche gemeinsame Gramform und die verschobenen Matrixelemente vorhanden sind. Die räumliche Viertelkoordinate des Assemblers ist in dessen ursprünglicher Domain ausdrücklich als räumlich markiert.

Die ganze Paarung auf diesem Raum lässt sich im endlichen ursprünglichen Wortchart so schreiben. Sei G der reflektierte Gram der Wörter, V eine Isometrie für seine Quotientennorm, C die tatsächliche unnormierte Randform im Quotienten, P der dortige physische Projektor und Q[S] die aus allen Feldern integrierte Form mit innerer Einfügung S. Für den vom installierten Aufruf verlangten Matrixkontraktionsfall ist der Wert
```text
B(S) = a Tr(C P V*Q[S]V P)/r,   r = Z_ref[0,0,0].
```
Dabei ist a die im Paket separat geführte absolute Quelle. C, P und der Referenznenner bleiben für den Vergleich S=I und S=QΘ fest. Die Lage von S bleibt zwischen den ursprünglichen Halbworten. Wenn S Bosonenrandwerte mischt, wird zuerst über beide offenen Seiten mit dem ganzen Operatorkern kontrahiert. S wird nicht als Zahl pro einzelner Geschichte behandelt. Falls derselbe ursprüngliche Zeitschnitt durch tatsächlich ausgewertete Formen Q₀ und Q_t auf einem reduzierenden Wortchart gegeben ist, gilt V*Q₀V=I und der positive Transfer ist V*Q_tV. Sein Generator ist auf seinem tatsächlichen Definitionsbereich −log(V*Q_tV)/t. Für mehrere Zeiten muss die Halbgruppenverklebung gelten. Eine einzelne räumliche Einteilchenmatrix liefert dafür keine der beiden Formen. Das entspricht der Voraussetzung eines tatsächlich ausgewerteten reflektierten Transfers bei [Osterwalder und Seiler](https://www.sciencedirect.com/science/article/pii/0003491678900398/pdf?md5=183a8aa4c17f212b3e6acff2a50f8f1a&pid=1-s2.0-0003491678900398-main.pdf).

Für das leere Wort, dessen beide Hälften in der rekonstruierten admissiblen Zeit auf Ran P bleiben, folgt B_leer(QΘ)=B_leer(I). Das ist eine exakte Identität des vollständigen Integrals, auch wenn sein Zahlenwert fehlt. Die [Referenzformel im Paper](/Users/stefanhamann/Projekte/tfpt-theoryv4/documents/archive--tfpt-45--source_extracts--04_qft_source/sections/008-section-relative-reflection-positivity-with-fermions.tex:143) fixiert r als separat auszuwertenden unmarkierten Referenzaufruf. Die Gleichheit der zwei leeren Antworten setzt weder r noch B_leer gleich eins. Die in der Wirkung enthaltene Operatorsubtraktion mit D_ref ist von diesem äußeren Nenner zu unterscheiden. Eine APS oder Pfaffianphase muss mit ihrer Randlinie verklebt werden. Die allgemeine Verklebungsregel dafür gibt [Dai und Freed](https://arxiv.org/abs/hep-th/9405012). Der erneut geprüfte v173 Test benutzt eine konstruierte endliche Matrix und liefert keine solchen Randdaten für den ursprünglichen Aufruf.

## Die gewichtete Feldintegration beim Verfeinern

Schreibe y für die feinere Konfiguration und x=π(y) für die sichtbare Konfiguration. Wenn die tatsächlich gegebene Grundverteilung im Regulator disintegriert ist als ν_m(dy)=ν_n(dx) κ_x(dy), dann ergibt die Integration über die verborgenen Felder exakt
```text
M_n^S(x) = ∫ κ_x(dy) w_m(y) K_m^S(y).
```
Hier ist K_m^S die volle gebundene Matrixform nach den für diese Faser zulässigen Kontraktionen. Für einen nichtdiagonalen Schnitt hat diese Formel dieselbe Bedeutung erst auf dem gemeinsamen Raum beider Historien mit beiden Randargumenten. Das folgt durch Disintegration des vollständigen Integrals, nicht aus der Kompression eines Diracoperators. Wenn die Referenzen r_n und r_m separat vorgegeben sind, ist Gleichheit der relativen Antworten für alle beschränkten groben Testfunktionen und alle Matrixelemente äquivalent dazu, dass die entsprechenden integrierbaren Dichten nach Division durch ihren jeweiligen Referenznenner fast überall gleich sind. Für die absolute Antwort muss die ursprüngliche Relation der beiden Referenzen erhalten bleiben.

Der [Paperbeweis zur projektiven geometrischen Quelle](/Users/stefanhamann/Projekte/tfpt-theoryv4/documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/005-section-absolute-dimensionless-metrology-from-the-bounda.tex:137) verweist an dieser Stelle auf die Kompatibilität der Diracoperatoren. Das reicht allein für die obige Dichtegleichung nicht aus. Der folgende tatsächliche Assemblertest zeigt, dass selbst ein identischer grober Operator einen feinen Beitrag zur spektralen Wirkung nicht bestimmt. Damit ist die fehlende Gleichheit konkret und prüfbar. Der Test schließt eine Herleitung mit der noch zu bestimmenden ursprünglichen Grundverteilung und vollständiger Fasergewichtung nicht aus.

## Ausgeführter Test im Originalassembler

Der [installierte Assembler](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime/assemble.py:79) baut h_N aus dem unveränderten v252 Yukawa und Majorana Träger, der Viertelphase α=1/4 und den vier räumlichen Fourierplätzen −2,−1,0,1. Die Originaldomain blieb unverändert. Für die Probe habe ich eine temporäre Kopie mit 14 Plätzen −7,…,6 und einem zusätzlichen reellen Higgskoeffizienten a cos(8x) in der vorhandenen Richtung 0 verwendet. Alle drei ursprünglichen nichtkonstanten Vorschlagskoeffizienten wurden für den Vergleich auf null gesetzt, der ursprüngliche konstante Hintergrund blieb erhalten. Die 14 Plätze, die zusätzliche Richtung und a=0,13 sind Testparameter. Der Lauf erzeugt daraus keine physische Zeit oder Feldverteilung.

Schreibe A=−iγ_F V₀. Dann ist A=A* und die zusätzliche Mode verbindet n mit n±8 durch den Block aA/2. Das liefert sechs ungeordnete Paare in der Verfeinerung, davon genau zwei mit einem Ende im ursprünglichen Block: (−2,6) und (1,−7). Auf den groben Plätzen ist P h(a)P=P h(0)P für alle a. Dagegen gilt exakt
```text
Δ Tr h² = (14−8) a² Tr(A²)/2,
Δ P h² P = (a²/4) A² auf den Plätzen −2 und 1,
            0 auf den anderen beiden Plätzen.
```
Die feine Wärmespur ändert sich für hinreichend kleines s um −s ΔTr h²+O(s²). Das gilt ebenso für die Differenz der Wärmespuren gegen denselben festgehaltenen Vergleichsoperator. Ob der in der vollständigen Wirkung tatsächlich verwendete Referenzoperator mit dem verborgenen Feld mitläuft, wird dabei nicht festgelegt.

Der vollständige Rückweg eines Paars ist ebenfalls bestimmbar. Für die beiden ursprünglichen Diagonalblöcke h_x und h_y liefert das Schurkomplement am groben Platz x
```text
[(h(a)−z)⁻¹]ₓₓ = [h_x−z−(a²/4) A(h_y−z)⁻¹A]⁻¹.
```
Bei z=i hat der eingeschobene Rückweg für a≠0 einen strikt positiven Imaginärteil auf der Bildrichtung von A. Der Rückweg verschwindet daher nicht. Am tatsächlich benutzten v252 Diagnosesatz ergaben sich ΔTr h²=18,56820153897, ΔTr P h²P=3,09470025650 und ΔTr exp(−0,001 h²)=−0,01744880296. Die Änderung des groben Resolventenblocks bei z=i hatte Frobeniusnorm 0,00411936742. Direkte Paarinversion und Schurformel stimmen bis 1,2·10⁻¹⁵ überein. Diese Diagnosewerte betreffen räumliche Spektralfunktionen. Sie sind keine Messung der physischen Zeit oder der ganzen relativen Wirkung.

## Welche markierten Werte die verborgene Faser bestimmen

Für das tatsächlich getestete Paar mit a und −a gibt es eine exakte Diagonalunitäre, die auf allen groben Zuständen die Identität ist und die beiden feinen Operatoren ineinander überführt. Denn jeder durch die achte Harmonische verbundene Paarblock enthält höchstens einen groben Platz. Das Vorzeichen auf dem anderen Platz wechselt den Übergang. Deshalb gilt für jede Borel Funktion f
```text
P f(h(a))P = P f(h(−a))P .
```
Selbst der vollständige komprimierte spektrale Leser kann die Gewichte der beiden Vorzeichen nicht unterscheiden. Der offene Übergang zum Partnerplatz unterscheidet sie. Im Test ist der erste mit A kontrahierte Übergang ±23,80538658843. Welche Seite physisch gemeinsam mit QΘ und dem ursprünglichen Eichtransport zugänglich ist, muss aus der originalen Randabbildung folgen. Auf eine Vorzeichensymmetrie der vollständigen markierten Quelle wurde dabei nicht geschlossen.

Für genau diesen einzelnen verborgenen Fourierkoeffizienten gibt es einen vollständigen bedingten Rückleser. Sei eine endliche skalare Faserverteilung σ auf einem vorgegebenen kompakten Amplitudenintervall gegeben. Dann ist
```text
Tr[(h(a)^(2k))ₓₓ] = 4^(−k) Tr(A^(2k)) a^(2k) + niedrigere gerade Potenzen.
Tr[A (h(a)^(2k+1))ᵧₓ] = 2^(−(2k+1)) Tr(A^(2k+2)) a^(2k+1)
                            + niedrigere ungerade Potenzen.
```
In einer Produktentwicklung trägt nur die durchgehend zwischen x und y wechselnde Folge die jeweils höchste Potenz. Ihre Koeffizienten sind für k≥1 im geraden und k≥0 im ungeraden Fall strikt positiv. Der integrierte nullte Wert gibt die absolute Masse von σ. Aus den integrierten geraden und ungeraden Koeffizienten gewinnt man die Momente von σ der Reihe nach, da der jeweils führende Koeffizient ungleich null ist. Auf dem kompakten Intervall bestimmen alle Momente die endliche Verteilung. Für die drei tatsächlich im Test benutzten Amplituden −u,0,+u genügen schon m₀,m₁,m₂: Die beiden nichtnull Gewichte sind (m₂/u²±m₁/u)/2, das Gewicht bei null ist m₀−m₂/u². Dazu braucht es den unnormierten nullten Wert, einen Rückweg und einen offenen Übergang. Diese Inversion setzt eine gemeinsame skalare Faserverteilung auf genau diesem Koeffizienten voraus. Sie weist der physischen TFPT Quelle keine solche Reduktion zu.

## Tatsächlicher Quellenstatus

Der Quellenaufruf wurde auf der ursprünglichen, unveränderten Domain separat vor dem Fourierexperiment ausgeführt. Sein [Endpunkt und Faseradapter](/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/runtime/original_source_adapter.py:280) meldet `evaluate_original_full_blocked_endpoint_source` und `original_blocked_fiber_packet`. Alle vier Werte Z0, ZU, ZY, ZUY und die originale Maßmasse wurden als null im Sinn von fehlender Ausgabe zurückgegeben. Die belegten Bundles im Repository sind als arithmetische Kontrollen registriert. Die ursprünglichen Verifikationen v173, v252, v259 und v330 liefen mit 4/4, 11/11, 6/6 und 6/6 bestandenen Checks. Der [neue reproduzierbare Test](check_qtheta_refinement_20261111.py) lief ebenfalls erfolgreich.

Die ganze physische Integration benötigt die ursprüngliche unnormierte Übergangsform über die beiden Ränder samt wirklicher Zeitverschiebung, die QΘ Wirkung auf den offenen Halbworten und den zugehörigen Referenzwert. Die hier ausgeführte Rechnung bestimmt den ersten konkreten Teil der Feldinformation, die bei der Verwendung der groben Matrix verloren geht. Für eine Fortsetzung am Original ist der offene Übergang mit seiner adjungierten Rückkehr auf derselben tatsächlichen reflektierten Randpaarung auszuwerten, zusammen mit der unnormierten leeren Masse. Z₀, F(0), C_UY und der physische Zeitgenerator bleiben unberechnet.
