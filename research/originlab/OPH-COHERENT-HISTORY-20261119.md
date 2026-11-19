# OPH Randwege und TFPT · 19. November 2026

Forschungsstand OPHCOHERENT20261119. Die folgenden Resultate betreffen deklarierte endliche Modelle. Die ursprüngliche gewichtete TFPT Quelle, die native Zeitkarte sowie 3+1 dimensionale Lokalität und wechselwirkende Gravitation bleiben offen.

## Gruppe und Verklebung

Die affine McKay Konstruktion des OPH Spindubletts liefert E₈. Die Löschung des vierdimensionalen spinoriellen Knotens liefert D₅⊕A₃. Die gelöschte Wurzel projiziert mit quadrierten Normen 5/4 und 3/4 auf die beiden Faktoren und erzeugt die zyklische Verklebung der Ordnung vier. Die Wurzelklassen zählen 52, 64, 60, 64. Beide Galoisgraphen liefern dieselben Aussagen. Die OPH Außenalgebra auf dem deklarierten 3+2 Träger stimmt mit dem TFPT Materieträger einschließlich Hyperladung überein, bis auf den neutralen Λ⁰ Zustand, den TFPT mitführt.

Der antipodale Quotient des regulären OPH Ikosaeders ist eine projektive Ebene. Seine flache Orientierungscharakteristik hebt sich nicht zu einer primitiven C₄ Holonomie. Auf dem Möbiusband nach Entfernen einer Scheibe ist eine solche Abbildung der Fundamentalgruppe möglich. Zwei antipodale Paare von Flächenzentren liefern unter den 45 Kombinationen keine harmonische Quadratmarkierung. Ihre quadrierten Richtungskosinus sind 1/9 oder 5/9. Das Ergebnis gilt für die genannte Wahl regulärer Flächenzentren.

## Spintransport und gemeinsame Phasen

Der Vergleich des TFPT Paars D,T mit den 60 OPH Portrotationen ergibt eine Antikommutation der unmodifizierten Spinhebungen gegenüber der erforderlichen Kommutation von D und T². In einem Modell mit vier neutralen Zusatzrichtungen korrigiert ein Vorzeichenwechsel auf zwei Richtungen die gemeinsame Paritätszerlegung von (7,3,3,3) zu (6,4,4,2). Eine Version mit extern vorgesehener Folge von 64 Kantenblöcken erzeugt diese Korrektur im Sektor mit einer Last. In der deklarierten algebraischen Spinhebung stimmt dann die gemeinsame Wirkung auf allen 248 Generatoren überein. Die beiden D₈ Halbspinoren bleiben gemeinsam isospektral. Der ausgerichtete Achtladungsleser trennt sie mit +1/2 und −1/2. Die homogene OPH Zellentwicklung allein erhält die Rotationssymmetrie und erzeugt die richtungsabhängige Korrektur nicht.

SL(2,F₅) enthält fünf binär tetraedrische Untergruppen. Beide OPH Dubletts beschränken sich dort auf den TFPT Picard Spinor. Im 128 dimensionalen Halbspinor der deklarierten Porthebung gilt 8·4₋⊕16·6₋. Ein eindeutiges Dublett entsteht im Tensorprodukt des spinoriellen McKay Quartetts mit einem ausgewählten OPH Tripel. Der gemittelte Projektor auf die stationären Komponenten der Flächendrehungen hat auf den Sektoren 2, 2′, 4 und 6 die Eigenwerte 0, 0, 1/2, 1/3. Auf dem gewählten Tensorprodukt isoliert sein Kern das Dublett.

## Geordnete Suche und Ladung

Ein adressiertes Instrument auf fünf OPH Ports liefert den positiven Filter diag(5/8,1). Das normierte inverse Instrument hat mit ihm die zustandsunabhängige Rückkehrwahrscheinlichkeit 25/64. Die Erklärung verwendet das OPH Kantenmodell bei κ=8, die kohärente Unterraumabfrage und eine vorgegebene Steuerung.

Eine kleinere Grammatik verwendet ausschließlich eine lokale Phasendrehung und die Projektion auf den nächsten gerichteten Flächenrand. Alle 60 Randkarten liefern dieselben Matrizen. Reine Wörter erzeugen projektiv eine diskrete Gruppe über gaußschen ganzen Zahlen. Kohärentes Zusammenführen zweier gleich langer Rückwege ergibt beim ersten Kandidaten der geprüften Klasse nach sechs Schritten und gemeinsamer Phasenkorrektur diag(63,65)/128 und diag(65,63)/128. Die beiden Zweige werden nach Normierung invers. Dreieckswege und lokale Phasen liefern obere und untere unipotente Matrizen mit zwei reell unabhängigen komplexen Richtungen. Konjugation durch den Filter erzeugt beliebig kleine Schritte. Der projektive Abschluss ist deshalb PSL(2,C).

Die Differenz der beiden Ausgaben bei festem Basisinput liegt exakt zwischen −1/64 und +1/64. Eine vorherige klassische Wegaufzeichnung ergibt null. Der Extremwert fixiert im betrachteten Zweiroutenmodell die vollständige relative Kohärenz. Von den 16 C₄ Phasenbelegungen zerfallen acht auf dem ausgewählten Rückkehrtyp projektiv in einen neutralen Spinortransfer und eine Ladungsphase. Die übrigen acht erzeugen entgegengesetzte Boosts in den geraden und ungeraden Graden. In dieser Klasse kann man danach die zwei Lorentzwirkungen unabhängig steuern. Ein deklarierter Austausch vollständiger lokaler Register koppelt Ladung und Spinor über zwei Regionen bei gemeinsamer Lorentzwirkung. Der positive endliche Kontaktcheck ergibt 9/100 für zwei gekoppelte Standortvariablen. Diese Variablen sind keine originalen R113 Einfügungen.

## Ausführung und Quellen

- `python3 -B research/originlab/tfpt_oph_spintransport_20261119.py`
- `python3 -B research/originlab/tfpt_oph_ursprungsregel_20261119.py`
- `python3 -B research/originlab/main.py --finite-controls-only --campaign f1_oph_coherent_histories --out /tmp/tfpt-oph-f1`
- [OPH Nahtoperator](https://github.com/FloatingPragma/observer-patch-holography/blob/main/model/polynomial_sv/SPEC.md)
- [OPH Portrotationen](https://github.com/FloatingPragma/observer-patch-holography/blob/main/Lean/Screen/A5PortAction.lean)
- [OPH McKay Rechendaten](https://github.com/FloatingPragma/observer-patch-holography/blob/f52fa4dab609845b20a80b43ea6bd9cc11a19f56/code/a5_closure/receipts/external/oph_sl2f5_mckay_e8_4ae2148a.json)
- [TFPT Operatorvergleich](../../docs/research/fundamental-plan/TFPT-Gesamtkonsolidierung-2026-11-12.md)

Nächster Anschluss in dieser Modellklasse: Eine gemeinsame positive Entwicklung muss Vakuummischung und geladene Rückkehr auf derselben ursprünglichen Zeit und derselben Referenz tragen. Die endlichen OPH Ergebnisse bestimmen diesen Anschluss noch nicht.
