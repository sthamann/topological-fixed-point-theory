# R104 und R113 · Ergebnisse und nächster Originaltest

13. November 2026.

## R113 · Orientierung fehlt selbst bei vollständiger radialer Information

Der native Aufruf für die vollständige Quelle meldet
`original_full_source_operation_not_supplied`. Der Aufruf einer Faser
meldet `original_blocked_fiber_packet`. Die vier Originalmarken, der
geladene markierte Rückweg und der Leerwert derselben Kappe sind
unberechnet. Die Registrierung der Operation liefert keinen Quellenwert.

In der vorhandenen R032 Diagnose mit fester Verbindung und
Farbzentrum Z3 lässt sich die Unterbestimmtheit genauer zeigen.
Die vier reellen Historien laufen außen mit $h=\pm 7/20,\pm 2/5$
und innen mit $h=-11/20$, bei $p=1/4$, $T=4/5$.
Innerhalb jedes Vorzeichenpaars ist die gesamte Bahn von $h^2$ gleich.
Der geordnete Yukawamark $y_i=K_{Y,i}/K_{0,i}$ bei festem Seefaktor unterscheidet
die Vorzeichen. Alle vier bedingten Werte wurden unabhängig mit dem
ursprünglichen dichten Fréchet Aufruf verglichen.

Setze $\delta_a=y_{a+}-y_{a-}$, $\delta_b=y_{b+}-y_{b-}$
und $\epsilon=1/1024$. Die effektiven Massen der beiden
Vergleichsmaße sind

$$
q_\pm=(1/4,1/4,1/4,1/4)\ \pm\epsilon(\delta_b,-\delta_b,-\delta_a,\delta_a).
$$

Die Grundmassen sind $q_i/K_{0,i}$. Alle vier $K_{0,i}$ und alle
Massen sind rigoros positiv eingeschlossen. Beide radialen Historien
behalten jeweils Masse $1/2$. Die Änderung in $Z_Y$ ist exakt
proportional zu $\delta_b\delta_a-\delta_a\delta_b=0$.
Damit sind die ganze **gewichtete** radiale Historienverteilung und
$Z_0,Z_U,Z_Y$ identisch bei einer gemeinsamen skalaren Referenz.

Die Differenz der verbundenen Antworten ist

$$
C_+-C_-=2\epsilon(u_a-u_b)\delta_a\delta_b
\quad\in\quad [-0.008364,-0.008363].
$$

Die direkte Kontraktion gibt ungefähr $-0.0034126125$ und $0.0049511810$.
Die Vergleichsmaße sind frei gewählte Zeugen in der bedingten R032 Klasse.
Der Ausschluss betrifft die Ermittlung der gemischten Antwort aus den
genannten niedrigeren Daten. Die physikalische Frage ist die ursprüngliche
Verteilung der geladenen Orientierung innerhalb jeder radialen Historie.
Ein gemeinsamer skalarer Referenzfaktor kürzt sich in $C_{UY}$. Für
absolute Werte braucht man seine tatsächliche Bindung.

[Originale R032 Zeitfenster](../../experiments/theory-contracts/chat-triangulation-20261108/TFPT-R113-GETRENNTE-ZEITFENSTER-20261109.txt),
[Originaler nativer Aufruf](../../experiments/theory-contracts/chat-triangulation-20261108/check_original_qtheta_word_20261109.py),
[Intervallrechnung](r113.py).

## R104 · Kohärenz braucht eine geeignete Referenz

Der getrennte Originalfaktor hat 16 Moden und
$H=2(I-|\phi_\theta\rangle\langle\phi_\theta|)$. Die Petz Probe
gibt ein Bit bedingte gegenseitige Information und Spurabstand $1/2$.
Der kohärente Wortraum wird bei CAR Tiefe 8 vollständig. Die endliche
Wärmeentwicklung besitzt Operator Schmidt Rang 4 bei einem einzelnen
Modus am Schnitt und 5 bei mindestens zwei Moden auf jeder Seite.
Mit Zahlselektion modulo vier ist
$Z_4(\beta)=1+16511\exp(-2\beta)$. Die geordneten Fermionwörter
werden mit ihren Vorzeichen und allen eingeschobenen Zeitintervallen
absolut berechnet. Eine dichte kleine Rechnung prüft auch einen
nichtkommutierenden Zahlsektorfilter.

Jede Phase von R104 hat dieselben reinen Zahlwortantworten bei
gleichzeitigem Rahmentransport. Eine festgehaltene Determinantenreferenz
unterscheidet die Phase. Zwei kohärente Bereiche liefern unter gemeinsamer
Ladungsrotation eine relative Antwort proportional zu
$\exp(i(\theta_B-\theta_A))$. Für unabhängig drehbare regionale Rahmen
braucht der gemeinsame Operator eine Verbindung mit entgegengesetzten
Endpunktladungen. Der Vergleich arbeitet mit der jeweils erklärten
Fermionzahl. Eine Gleichsetzung mit Hyperladung oder primitiver
R102 Ladung bräuchte eine originale Abbildung.

[Originaler R104 Quellenvertrag](../../experiments/theory-contracts/original-local-source-variable-background-20261102/DET16-WHOLE-SOURCE-REDUCTION.txt),
[R104 Wortauswertung](coherence.py),
[Relationale Prüfer](relational.py).

## Zwei Wege und der gemeinsame Ladungskern

In einer abelschen Vergleichsverbindung sieht Ladung $q$ einen
Wegunterschied $\Phi$ als $e^{iq\Phi}$. Sind die gelieferten
ganzzahligen Ladungen $q_i$, ist der kleinste unsichtbare positive
Winkel $2\pi/\gcd(q_i)$. Jeder gemeinsame unsichtbare Winkel erfüllt
$\gcd(q_i)\Phi\in 2\pi\mathbb Z$, weil der größte gemeinsame Teiler
eine ganzzahlige Linearkombination der $q_i$ ist. Umgekehrt ist jedes
$q_i$ ein Vielfaches dieses Teilers.
Bei $\Phi=\pi/8$ sieht Ladung 16 denselben Wert wie bei Null. Ladung 4
oder eine primitive Ladung sehen einen Unterschied. Eine einzelne
Determinantenantwort kann daher die ursprüngliche Verbindung nicht
eindeutig bestimmen. Für einen vorhandenen Originalkandidaten ist die
Frage, welche zusätzlichen physikalisch zugelassenen Ladungen denselben
Pfad sehen und wie sie auf demselben Quellenkern wirken.

Graphschleifen und berandete Flächen sind getrennt prüfbar. Nichttriviale
Holonomie auf einer Schleife ohne ausgefüllte Fläche allein identifiziert
noch keine lokale Feldstärke. Der physikalische P1 Nachweis braucht
zwei zulässige regionale Wege, ihre gemeinsame absolute Kontraktion,
den Reflexionskern und den tatsächlichen regionalen Transport.

[Relationale Kontrollen](relational.py),
[Endpunktkontraktion](source.py),
[Regionale Originalfrage](../../docs/research/fundamental-plan/TFPT-REGIONALE-KOMPOSITION-20261112.txt).

## Der erste konkrete Originalkandidat

Das ursprüngliche E8 Stromnetz liefert einen geladenen Mittelstrom
$\gamma$ mit $Y=-1/3$ und nichttrivialer Farbtrialität. Der native
schwache Aufsteiger im R032 Vergleich hat $Y=0$ und $Q=+1$.
Die beiden Übergänge können im tatsächlichen Ladungswörterbuch
verschiedene Quellwerte tragen. Die überprüften nativen Q Generatoren
sind bezüglich $\exp(4\pi iY)$ blockdiagonal. Produkte davon erzeugen
das farbwechselnde Mittelbein ebenfalls nicht. Der nächste Originaltest
muss daher das ladungsneutrale ganze E8 Wort mit seinen geladenen
Zwischenzuständen vor der äußeren Auswahl am R113 Kern binden.
Der ausgeführte Originaltest bestätigt $-i/4$ auf dem gegebenen
E8 Vakuum. Seine absolute gewichtete R113 Antwort ist unberechnet.

[Originaler E8 Cartantest](../../docs/research/fundamental-plan/check_tfpt_cartan_process_20261112.py),
[Nativer P1 Ladungsvergleich](../../docs/research/fundamental-plan/p1-original-marken-und-metrik-20261112/TFPT-P1-ORIGINAL-MARKEN-UND-METRIK-20261112.txt).
