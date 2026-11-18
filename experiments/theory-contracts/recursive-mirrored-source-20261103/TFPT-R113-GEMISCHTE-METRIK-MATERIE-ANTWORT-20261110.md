# R113: gemeinsame metrische und materielle Antwort des ursprünglichen Fermionübergangs

Der unveränderte reflektierte R032 Verlauf trägt im selben Fermionübergang sowohl eine metrische als auch eine materielle Antwort. Die gemischte Antwort bleibt nach Einbeziehung des ursprünglichen Fermionsees und derselben freien S Referenz ungleich null. Der berechnete Übergang betrifft die gegebene Geschichte. Die Gewichtung verschiedener Geschichten bleibt eine eigene, noch unberechnete Operation.

## Originaleingabe und Rechnung

Übernommen sind die sechs Feldwerte, die sechs Zeitlängen, die 96 inneren Richtungen, die v252 Yukawakarte, der unveränderte Majoranaterm und die S Referenz aus `../chat-triangulation-20261108/check_original_whole_two_strips_20261109.py`. Der bereits vorhandene Logradius `σ` ändert `p` zu `exp(−σ)p`. `h` ändert ausschließlich den ursprünglichen Higgsanteil `Y(t)` zu `(1+h)Y(t)`. Majorana und die sechs Zeitlängen bleiben gleich. Gemessen wird die Summe der drei rechtshändigen geladenen Leptonenbesetzungen in der gegebenen Familienbasis. Das ist eine konkrete Fermionmarke, keine Ersatzdefinition für die ursprünglichen `q/r` Wörter oder für `Y_em`.

Das ganze einteilchige Übergangsprodukt sei `R`. Die Fermionbesetzung ist `n = Tr[P_R (I−(I+R)⁻¹)]`. Ein oberer Dreiecksblock trägt die Ableitungen in `σ` und `h` mit beiden Einsetzungsordnungen durch sämtliche sechs Übergänge. Für `G=(I+R)⁻¹` folgt der gemischte Jet der Besetzung aus `G R_σh G − G R_σ G R_h G − G R_h G R_σ G`. Keine interne Ladungsprojektion wird vor dem Endpunkt eingeführt. Der unnormierte S relative Wert ist `N = exp(L)n`, `L=log(Z_A/Z_S)` auf dieser Geschichte, mit der Energie des ursprünglichen gefüllten Fermionsees. Die S Referenz hat eine eigene metrische Antwort und bleibt vom Higgsparameter unabhängig. Die Rechnung differenziert auch den Seefaktor.

| Am Originalpunkt | Wert | metrisch `∂σ` | Higgs `∂h` | gemischt `∂σ∂h` |
|---|---:|---:|---:|---:|
| `n` | 1.647798481396 | −0.146783642927 | −0.003892652596 | +0.003917913585 |
| `L` | −16.614305487553 | −5.010567639854 | −17.687494973609 | −2.231558818935 |
| `N` | 1.003235529772e−7 | −5.116146334305e−7 | −1.774709336895e−6 | +8.826729377567e−6 |

Die unabhängige Rechnung diagonalisiert an den verschobenen Eingabepunkten jede Hermitesche Matrix neu und differenziert den ganzen Übergang numerisch. Sie ergibt für die beiden gemischten Besetzungsableitungen `0.003917864644` und `8.826749009134e−6`. Für alle fünf möglichen Zerlegungen des sechsfachen Produkts stimmt der normalisierte Jet mit dem ganzen Produkt bis auf `6.7e−16` überein. Der ausführbare Beleg ist `check_r113_mixed_metric_matter_20261110.py` neben dieser Datei.

## Was ein ganzer physischer Übergang zusätzlich braucht

Die beiden Antworten kommen hier aus demselben Operator und teilen dieselben Zwischenzustände. Der S relative Wert ist unnormiert bezüglich des gegebenen Fermionübergangs, aber er enthält noch kein ursprüngliches Gewicht der Geschichte. Eine gemeinsame Geschichte mal einem markenunabhängigen Skalar hat dieselbe normierte Besetzung. In einem Ensemble mit gewichtetem unmarkiertem Integral `Z=∫w` und Besetzung `⟨n⟩=∫wn/Z` tritt bereits in der ersten Ableitung der Zusatz `Cov(n,∂σ log w)` auf, sofern die positiven Dichten und ihre Ableitungen dies erlauben. Ein fester Verlauf und sein nichtverschwindender gemischter Jet entscheiden weder diesen Zusatz noch das Vorzeichen des ganzen gemischten Jets. Das ist auch der Grund, weshalb eine Geschichte ohne die tatsächliche Auswahl und Maßmasse die absolute Antwort nicht festlegt.

Die übertragbare Rechnung ist die markierte ganze Übergangskontraktion vor der Normierung, einschließlich der relativen Leeramplitude und beider Einsetzungsordnungen. Die physische Anwendung braucht dieselbe ursprüngliche Feldbasis, tatsächliche Auswahl, Zeit, Referenz und gemeinsame Gewichte, oder einen aus dem Original bewiesenen unnormierten Randcaller, der genau diese Beiträge liefert. Auf dieser Quelle kann dieselbe Art von Rechnung zusammen mit der `m=2` Norm und der `m=−1` Kontaktprobe ausgewertet werden. Der vorliegende gemischte Wert ist weder der physische `C_UY` noch eine Ableitung des Einstein Terms oder ein gemeinsames physisches Spektrum.
