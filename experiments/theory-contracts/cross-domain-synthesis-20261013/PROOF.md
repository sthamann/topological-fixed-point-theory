# Kompatibilität des vorhandenen Nahtkanals mit einer festen oder mitgeführten Markierung

13. Oktober 2026 · Endlicher algebraischer Folgesatz · Keine physische T1–T8-Schließung

## Ausgangspunkt

Der bestehende Vertrag `marked-seam-cp-lift-20260921` konstruiert auf M4(C)

\[
\Phi=\Pi_0+\tfrac13\Pi_5+\tfrac23\Pi_{10}
=\tfrac7{12}\mathrm{id}+\tfrac1{12}\sum_{j=1}^5\mathrm{Ad}_{g_j}.
\]

Die Pi sind orthogonale Hilbert-Schmidt-Projektoren auf Operatorräume mit Dimensionen 1, 5, 10. Dies sind keine Projektoren auf 1, 5, 10 physische Teilchensorten. Die fünf tatsächlichen Wörter sind hermitesche Clifford-Generatoren auf C4 mit Determinante eins. Der Kanal ist unital und spurtreu. Seine vorhandenen positiven Trine-Effekte realisieren die ursprüngliche Matrix B und B⁶; das ist ein früheres Ergebnis, kein Ergebnis des vorliegenden Nachtests.

## Exakter universeller Rang-3-Test

Die Rechnung mit den tatsächlichen fünf Wörtern bestätigt die Tensoridentität

\[
\Pi_{\mathrm{sym}}\left(\sum_j g_j\otimes g_j-I_{16}\right)
\Pi_{\mathrm{sym}}=0.
\]

Für jeden normierten komplexen Vektor psi liegt psi⊗psi im symmetrischen Tensorraum. Daher gilt mit q=|psi><psi|

\[
\sum_j\langle\psi,g_j\psi\rangle^2=1,
\qquad \operatorname{Tr}(q\Phi(q))=\tfrac7{12}+\tfrac1{12}=\tfrac23.
\]

Jeder Rang-3-Projektor ist P=I−q für ein solches q. Unitalität liefert

\[
\operatorname{Tr}((I-P)\Phi(P))
=1-\operatorname{Tr}(q\Phi(q))=\tfrac13.
\]

Für den **gleichmäßig gemischten Zustand rho=P/3** lautet die Wahrscheinlichkeit außerhalb desselben festen Unterraums deshalb exakt **1/9**. Diese Zahl wird nicht für jeden beliebigen Zustand in P behauptet. Sie genügt, die Invarianz des Unterraums auszuschließen, denn ein invariantes P müsste jeden in P unterstützten Zustand dort halten.

Für q haben die drei orthogonalen Operatoranteile die quadrierten Normen 1/4, 1/4 und 1/2: der skalare Anteil folgt aus Tr(q)=1, der Fünferanteil aus der Tensoridentität, der Rest aus Tr(q²)=1. Mit den vorhandenen Kanaleigenwerten ergibt sich somit

\[
L_n=\operatorname{Tr}((I-P)\Phi^n(P/3))
=\frac14-\frac1{12}\left(\frac13\right)^n
-\frac16\left(\frac23\right)^n.
\]

Insbesondere L0=0, L1=1/9 und L6=343/1458. Bei sehr vielen Schritten strebt L_n gegen 1/4. Die Rechnung enthält einen direkten Matrixabgleich für n=1 und n=6 zusätzlich zur allgemeinen Herleitung.

## Warum eine gleichzeitige Carrierrotation den Test nicht beseitigt

Der vorhandene Gruppenvertrag `compiler-quartic-holonomy-20260919` liefert

\[
\widetilde G\longrightarrow(\mathrm{Spin}(10)\times SU(4))/\mathbb Z_4,
\qquad h\mapsto[\Lambda_{\rm even}R_4(h),h].
\]

Das zugehörige geladene Darstellungsstück ist (16,4). Der Originaltest v783 weist die 16 Paulioperatoren in G31 nach; sein unabhängiger Lauf bestand mit 32 Prüfungen. Die fünf verwendeten g_j haben Determinante eins und liegen damit in G_tilde. Ihre korrelierten Lifts sind U_j=S_j⊗g_j. Eine zusätzliche zentrale Phase fällt aus der Konjugation heraus.

Wendet man **bedingt dieselben Krausgewichte** auf diese vorhandenen Gruppenlifts an, bleibt der Carrier von rho_full=(I16/16)⊗(P/3) in jedem Zweig maximal gemischt. Die Familienranddichte entwickelt sich genau mit Phi. Deshalb gilt dieselbe Leakageformel gegenüber dem festen 48-dimensionalen Projektor I16⊗P. Dies setzt weder den physisch wirkenden Kanal mit I16⊗Phi gleich noch beweist es, dass das ursprüngliche TFPT-Quellengesetz diese Krausausführung auswählt.

Die Schlussfolgerung ist begrenzt: Die unprotokollierte Ausführung dieses bestimmten Kanals kann keine geschlossene feste Drei-von-vier-Auswahl darstellen. Eine tatsächliche Identifikation dieses Projektors mit der physisch ausgewählten W48-Struktur ist weiterhin ein eigener Beweisauftrag. Die Rechnung ersetzt keine Abbildung von M4-Operatoren auf quadratische E8-Stromoperatoren.

## Der positive Anschluss: gespeicherter Zweig und bewegte Markierung

Für die vorhandenen random-unitären Zweige ist

\[
\mathcal I_j(\rho)=p_jU_j\rho U_j^\dagger,
\qquad \widehat P_j=U_j\widehat P U_j^\dagger.
\]

Ist rho anfangs in P_hat unterstützt, ist I_j(rho) in P_hat_j unterstützt. Das folgt unmittelbar durch Einsetzen. Wird der Zweig j als klassischer Record behalten, lässt sich auch die Markierung P_hat_j eindeutig mitführen. Über mehrere Schritte gilt dies für die geordnete Folge der U_j. Vergisst man die Zweige und liest anschließend gegen den ursprünglichen festen Projektor, entsteht die oben berechnete Leakage.

Dies zeigt einen algebraisch verfügbaren Ausweg. Es wählt **kein physisches Instrument** aus: Verschiedene Instrumente können denselben mittleren Kanal besitzen. Auch ist die Speicherung des Records keine aus P1/P2 bereits bewiesene Dynamik. Sie muss samt physischem Ausleseprotokoll aus derselben Quelle kommen. Ein mitbewegter Projektor bezeichnet zudem eine andere Auslese als der festgehaltene Projektor: Physisch definierte Ladungen dürfen nicht allein durch Umbenennen des Sektors als erhalten gelten. Ob dieser Transport eine zulässige Symmetrie, ein physischer Übergang oder nur eine Änderung der Auslese ist, muss die gemeinsame Feldabbildung entscheiden.

## Zusammenhang mit der übrigen TFPT-Synthese

Der aktuelle Wurzeltest bestätigt separat: Der bestehende Ladungscharakter bestimmt P5, und ein transportierter Charakter bestimmt den durch C/J transportierten P5. Zusammen geben beide Ergebnisse einen konkreten Grund, eine dynamische Markierung samt Record zu untersuchen. Die beiden Tests liegen auf unterschiedlichen Operatorräumen und werden hier nicht durch eine erfundene Abbildung identifiziert.

Der entscheidende nächste Quellenbeweis wäre eine explizite gemeinsame Ausführung, deren Records und Marks sowohl den vorhandenen Nahttransfer als auch die geladenen W48-Antworten liefern. Anschließend bleiben die gemeinsame Feldzeit, wechselwirkende 3+1D-Rekonstruktion und metrische Antwort zu beweisen. Der vorliegende Folgesatz schließt keinen dieser physischen Gesamtverträge.

## Reproduktions- und Neuheitsscope

Der Checker benutzt die tatsächlichen fünf Projektwörter, prüft ihre Clifford-Beziehungen, Determinanten und die Tensoridentität. Der allgemeine Rang-3-Satz folgt aus dieser Identität und der obigen Rechnung, nicht aus zufälligen Zustandsstichproben. Die mittransportierte Markierung und der korrelierte Carrierlift werden algebraisch begründet.

Der neue Beitrag dieser Sitzung ist der quantitative gemeinsame Kompatibilitätstest samt seiner Einordnung in die bestehende Mark-/Record-Frage. CP-Dilatation, Kraus-Instrumente, Clifford-Algebra und die vorhandenen Gruppenlifts sind bekannte bzw. bereits im Projekt enthaltene Strukturen. Eine weltweite mathematische Neuheit wird nicht behauptet.
