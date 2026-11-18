# Exakte ursprüngliche skalare Hilfsintegration mit Quellen

5. November 2026. Ergänzung zu SHORTTIME20261105G. Atlas beim Einstieg aa7a2fae187f8d99d93468a1. Nur Scratchdateien geändert.

Die komplexe skalare Hilfsintegration lässt sich einschließlich ihrer Quellen und neutralen Kontakte exakt ausführen. Hierfür ist innerhalb des vorhandenen euklidischen Quellenaufrufs keine neue Wickkarte oder Maßwahl erforderlich. Das Ergebnis eliminiert das eingeführte Hilfsbein und gibt eine bestimmte ursprüngliche zusammengesetzte Einfügung zurück. Die noch nicht ausgewerteten ursprünglichen Higgs und Geometriefelder werden dadurch nicht eliminiert.

## Ausgeführte Quellenrechnung

Auf der bereits gegebenen homogenen Residuenklasse seien α>0 und v>0 das tatsächliche regulierte Zellvolumen. Mit ursprünglichem F und R setze μ=F+2αR, ρ=F+iξ. Für eine Hilfsquelle θ, die als θρ in den Exponenten eingeht, gilt exakt

\[
\sqrt{\frac{v}{4\pi\alpha}}\int_{\mathbb R}d\xi\,
 e^{-vFR-v\xi^2/(4\alpha)-iv\xi R+\theta(F+i\xi)}
=e^{-v(FR+\alpha R^2)}
 e^{\theta\mu-\alpha\theta^2/v}.
\]

Dies ist eine tatsächlich ausgeführte reelle Gauß Fourierintegration auf der komplexen ρ Kontur. Keine reelle positive ρ Achse wird gewählt. Der absolute Vorfaktor bleibt erhalten. Die normalisierte Antwort des Hilfsbeins lautet daher

\[
\langle\rho\rangle_{\mathrm{aux}|b}=\mu,\qquad
\langle\rho^2\rangle_{\mathrm{aux}|b}=\mu^2-\frac{2\alpha}{v},\qquad
\langle\rho\rho\rangle_{\mathrm{aux}|b,c}=-\frac{2\alpha}{v}.
\]

Insbesondere ist die neutrale quadratische Einfügung nicht allein μ². Der kontaktkorrigierte Ausdruck ρ²+2α/v bildet auf μ² ab. Für verschiedene Zellen ist

\[
\langle\rho_i\rho_j\rangle_{\mathrm{aux}|b}
=\mu_i\mu_j-\delta_{ij}\frac{2\alpha_i}{v_i}.
\]

Das negative Kontaktvorzeichen ist eine komplexe Hilfsfeldantwort, keine negative Varianz eines neu gewählten physischen Wahrscheinlichkeitsmaßes. Mit einer als v kρ geschriebenen skalaren Quelle ist θ=v k und der Quellfaktor entsprechend exp[v kμ−vαk²]. Im Kontinuum bleibt der durch den ursprünglichen Regulator induzierte diagonale Kontakt, formal −2αδ_g(x,y), erhalten.

Der tatsächlich gedruckte ursprüngliche Koeffizient wird ohne Fit eingesetzt:

\[
\mu=\frac{\chi_0^2}{4\pi^2}+\frac{u}{96\pi^2}+2\alpha R,
\qquad u=\Phi^\dagger\Phi.
\]

α bleibt der tatsächliche Koeffizient desselben ursprünglichen Zweigs. Die früher gefundene rohe und kanonische M Konvention wird nicht still gleichgesetzt.

## Ganze Quelle und bereits geschlossene Karten

Bezeichne die ursprüngliche unnormierte Integration einschließlich aller übrigen Felder, Regulatoren, Selektoren, Eta und Referenzdaten mit I_orig. Die ursprünglichen externen Quellen J werden nachstehend explizit eingesetzt. Die Quellenidentität nach zuerst ausgeführter Hilfsintegration lautet

\[
Z_{\mathrm{aux}}[J,\theta]
=I_{\mathrm{orig}}\!\left[
 e^{\langle J,b\rangle}
 e^{\sum_i[\theta_i\mu_i(b)-\alpha_i\theta_i^2/v_i(b)]}
\right].
\]

Bei dynamischer Geometrie bleibt v_i(b) innerhalb dieser ursprünglichen Integration. Der Kontaktfaktor wird nicht als feldunabhängiger äußerer Skalar herausgezogen. Insbesondere

\[
\left.\partial_{\theta_i}Z_{\mathrm{aux}}\right|_0=I_{\mathrm{orig}}[\mu_i],
\quad
\left.\partial_{\theta_i}\partial_{\theta_j}Z_{\mathrm{aux}}\right|_0
=I_{\mathrm{orig}}[\mu_i\mu_j-2\alpha_i\delta_{ij}/v_i].
\]

Alle bestehenden Quellen J und ihre tatsächlichen Kontakte bleiben erhalten. Die Hilfsquelle θ ist ein Diagnosebein für die explizit angegebene ursprüngliche zusammengesetzte Einfügung, kein neuer freier physischer Input.

Für jeden tatsächlich gelieferten ursprünglichen Caller kann dieselbe Identität als normierte Eins eingesetzt werden: der Hilfsintegrand wird zusätzlich mit exp[+v(FR+αR²)] multipliziert. So ist die Erweiterung auch dann wohldefiniert, wenn nicht behauptet wird, dass die lokale Residuenwirkung die ganze S_adm ist. Eine nützliche Ersetzung des FR+αR² Faktors setzt dessen tatsächliches Auftreten in genau diesem Caller voraus; sämtliche Restterme bleiben unverändert. Γ_res wird hier nicht mit S_adm, W oder Γ_ren gleichgesetzt.

Die nullsourcige Hilfsintegration ist punktweise exakt die Eins. Daher werden weder die ursprüngliche Referenz noch Z[0] neu normalisiert oder ein zweiter Referenzfaktor abgezogen. Derselbe vorhandene volle U Chart trägt die ursprüngliche Quelle weiter. Der ursprüngliche BRST Quotient und seine Wardregel bleiben dieselben; die direkt substituierte physische zusammengesetzte Quelle erhält die bereits vorausgesetzte Wardgültigkeit. Keine zusätzliche globale Wick, Gaugefermion oder Referenzpflicht entsteht.

Die ausführbare Integrationsreihenfolge ist ausdrücklich Hilfsbein zuerst. Ein Vertauschen mit der ursprünglichen Geometrie oder Higgsintegration wird daraus nicht abgeleitet. Ebenso wird die euklidische komplexe Karte nicht mit einer unabhängigen positiven Lorentz Scalaronquantisierung oder einer H5 Wärmequelle identifiziert.

## Erste noch tatsächlich benötigte ursprüngliche Einfügung

Bereits der erste neue Hilfsquellenwert verlangt genau

\[
I_{\mathrm{orig}}[\mu_i]
=\frac{\chi_0^2}{4\pi^2}Z_{\mathrm{orig}}[0]
 +\frac{I_{\mathrm{orig}}[u_i]}{96\pi^2}
 +2\alpha I_{\mathrm{orig}}[R_i]
\]

im hier gegebenen konstanten χ0 Zweig und bei konstantem α. Bei variablem α bleibt dessen Produkt innerhalb der Integration. Die zweite Antwort verlangt die ganze gemeinsame Einfügung I_orig[μ_iμ_j], zusätzlich zum jetzt bestimmten lokalen Kontakt. Eine neue Wardregel oder klassische stationäre Gleichung setzt diese ursprünglichen Momente nicht ein. Die hier ausgeführte Gaußoperation produziert somit den genauen Quellkontakt und die vollständige Einfügungsabbildung, aber noch keinen fehlenden primitiven gemischten UY Wert oder den ganzen physisch integrierten Kurzzeitkernel.

Die nächste Rechnung muss an diesen tatsächlichen ursprünglichen Einfügungen beziehungsweise dem bereits festgelegten P1 Wortjet ansetzen. Es gibt keinen Anlass, zuerst eine neue positive Hilfsfelddichte, H5 Wärme oder zweite Wickkarte zu bauen.

## Originalbelege und ausgeführte Kontrolle

Originalstellen mit SHA256 sind in check_auxiliary_sources.json beigefügt:

* CRAIG20261103F, RESIDUE-REFERENCE-AND-FULL-INTEGRATION.txt:20–60: ursprünglicher Caller und lokale Koeffizienten;143–173: ganze Integration und Typgrenze. SHA256 353f446053335bc35b539b3dfa660b5fa68514be43a5668fcebd2298d79bb7f1.
* REFERENCE-SOURCE-COMPLETION.txt:82–140: vollständiger Quellen und Referenztransport. SHA256 5ddc9c63c5f20f42befba6ae32fc57fdecb61c4b70dd9ab2d7b26fe667380cbb.
* FULL-FUNCTOR-REVIEW.txt:65–125,130–162: ursprüngliche Bosonoperation und ihre tatsächlich mitgeführten Quellen. SHA256 81b0c0588236e02ded5d3893abcd804674ccc0b90d14fd589c977f9a606cc94e.
* WEIGHT-FULL.txt:45–95: bereits vorhandener physischer Quotient und Wardgültigkeit in ihrem ursprünglichen Geltungsbereich. SHA256 e6776276275e6ceb9db60f60c0f7d9a7f85383c00dd2558a18f84faae2ec1438.
* Source45 QFT007:1–28: ursprünglicher ganzer Quellenaufruf mit Regulator. SHA256 9c1c587cbd61fc4034649d778d16f1ee0e8c1e1189885babb3d8148a6fe72d7c.

check_auxiliary_sources.py hat15 von15 exakten Kontrollen bestanden. Die direkte reelle Gauß Fourierintegration und die unabhängig über R Ableitungen integrierten Einfügungen bis Grad4 stimmen überein. Geprüft sind absoluter Faktor, ganze Quellabhängigkeit, erster und zweiter Wert, neutraler Kontakt und seine Korrektur sowie die originale F Substitution. Dies ist eine ausgeführte skalare Hilfsintegration; keine Simulation des gemeinsamen ursprünglichen Feldmaßes und keine physische Statuspromotion.
