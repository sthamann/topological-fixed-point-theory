# IDEAS — neue Fibonacci Ursprungsspur

**Gebiete:** unitäre Fusionskategorien, Zopfgruppen, Topologie, Quanteninformation, Operatoralgebren und 4D Topologie. Stand 20. November 2026.

## Ursprung, nicht TFPT als Input

Annahme: die einfachste unitäre **nichtpunktförmige** Fusionskategorie von Rang zwei. Sie besitzt die Fibonacci Regel tau ⊗ tau = 1 ⊕ tau. Die Wahl von Rang zwei und Nichtinvertibilität ist das erste freie Prinzip; die Natur muss es noch auswählen.

Der erste nichttriviale Assoziativitätstest ist das Pentagon für vier Faktoren. Deren fünf Klammerungen bilden einen fünfeckigen Kohärenzgraphen. Vier tau Ereignisse erzeugen zwei mögliche Gesamtladungssektoren mit genau 2 beziehungsweise 3 Fusionsgeschichten. Ihre direkte Summe hat fünf Dimensionen und eine **kanonische Endladungsaufteilung 3+2**. Ein geschlossenes System mit Gesamtvakuum besitzt aber nur zwei Dimensionen; die 3+2 Summe verlangt einen physisch erlaubten Randladungssektor.

## Prüfbare Strukturergebnisse

1. **Pentagon:** F Matrix in jedem Sektor; die Wege mit zwei beziehungsweise drei Reassoziationen sind exakt identisch. Die erste arity vier liefert die fünf Klammerungen.
2. **Braid-Symmetrie:** Auf dem Fusionraum 2+3 erzeugen B1,B2,B3 unitäre Zopfdarstellungen. Numerisch ist die Lie Abschlussspur 3+8=11, die Dimension von su2 ⊕ su3; das wird von bekannten universellen Jones Dichtheitssätzen unterstützt. Die natürliche sektorverträgliche spurfreie Lie Algebra ist su3 ⊕ su2 ⊕ u1, wobei die kontinuierliche u1 **nicht** durch Braidwörter generiert wird.
3. **Exakte μ4 Kandidatenstruktur:** Für C=B1 B2 B3 und J=−C^5 ergeben die Matrizen im Feld Q(zeta_10): J^4=I, J²≠I, det(J)=1. Auf dem zweidimensionalen Sektor gilt das charakteristische Polynom x²−1, auf dem dreidimensionalen x³+x²+x+1. Gesamt-Eigenwerte sind 1,−1,−1,i,−i. Das Rechenprogramm exact_core.py belegt dies ohne Toleranz.
4. **Drei geometrische Richtungen aus derselben Vier:** Vier unterschiedliche Ereignispunkte auf S² liefern H1 Rang 3 und eine natürliche ganzzahlige reduzierte S4 Darstellung A3 = {x∈Z⁴ | sum x=0}. Der dortige Viererzyklus hat dasselbe charakteristische Polynom x³+x²+x+1 wie der entsprechende 3D Braid Uhrsektor. **Aber es gibt keine Äquivalenz der vollständigen B4 Darstellungen**: im geometrischen Modul sind Transpositionen involutiv, im Fibonacci Modul nicht.
5. **Fünf geometrische Plätze:** Vier Punkte auf CP1 bestimmen einen Divisor D vom Grad 4 mit h0(O(D))=5; der hyperelliptische Doppelcover hat Genus 1; eine *zusätzlich gewählte* zyklische Vierfachüberlagerung hat Genus 3 und C4 Deckgruppe. Es fehlt eine kanonische physische Identifikation der beiden fünfdimensionalen Objekte.
6. **Hyperladung:** Spurfreiheit auf C³ ⊕ C² bestimmt den relativen Generator proportional zu diag(−1/3,−1/3,−1/3,+1/2,+1/2). Die absolute Skalierung ist frei. Unter dem **zusätzlichen** CAR Fermionisierungsaxiom ergeben fünf ungerade Einteilchenmoden exakt 16 Zustände in der geraden äußeren Algebra und die bekannten Standardmodell Multiplets; die intrinsische Fibonacci Statistik ist aber anyonisch, nicht gewöhnlich fermionisch.
7. **E8 ist nur bedingt:** Bei *zusätzlich* verfügbarer CAR/Spin10 Struktur (D5) und einem integralen A3 mit passender normierter Glue Auswahl kann wie bisher E8 entstehen. Keines dieser physikalischen Brückenglieder folgt schon aus der Fibonacci Fusionsregel. Ihre chirale Ladung c=14/5 mod8 ist nicht die des E8_1 Vakuums (c=8).
8. **3+1D ist nicht emergent:** Walker–Wang und Crane–Yetter liefern zwar 3+1 dimensionale topologische Realisierungen aus einer modularen Fusionskategorie, aber sie setzen den 3D Raum bzw. das 4D Bordismusobjekt voraus, und sie erzeugen keinen Einstein Graviton.

## Baselines und exakte Gegenmodelle

- In der Klassifikation der **beiden Rang-zwei-Fusionsringtypen** ergibt der punktförmige Fall tau²=1 nur (1,0), der nichtpunktförmige Fibonacci Fall tau²=1+tau genau (2,3). Nach Bedingung auf Nichtinvertibilität ist die Eindeutigkeit Teil der Definition, kein p Wert.
- Für formale positive Fusionsringe X²=1+nX, n=0..6, trifft nur n=1 die 3+2. Die Ringe n>=2 sind keine gültigen Rang-zwei-Fusionskategorien, also nur algebraische Kontrollen.
- Unter arity 2..7 ergeben nur vier Faktoren genau fünf Klammerungen, aber diese Stelligkeit wurde durch den Pentagontest gewählt.
- Unter 38 **nachträglich** betrachteten Vorzeichen und Potenzen ±C^k, k=1..19, ergeben exakt (5,−1) und (15,−1) zugleich determinant1, primitive Ordnung4 und den A3 Charakter im dreidimensionalen Sektor. Diese Basislinie darf nicht als vorregistrierte Wahrscheinlichkeit gelesen werden.
- Die gleiche Kategorie erlaubt verschiedene Kreuzverhältnisse. λ=1/2 liefert j=1728, λ=1/3 dagegen j=21952/9. Keine Selbstwahl der Quadratstruktur oder von P1/c3.
- Jeder vollständig Braid-invariante treue Zustand auf der beobachtbaren Algebra M2 ⊕ M3 ist sektoral skalar und zentral. Daher ist sein Tomita Modularfluss identisch, also **keine physische Zeit** aus der symmetrischen endlichen Algebra allein.
- Der einfache Walker–Wang/Crane–Yetter Lift ist topologisch und kann keine dynamische Lorentzmetrik mit Einstein Wirkung erklären.

## Physischer Gate zum nächsten Fortschritt

Konstruiere **einen einzigen ursprünglichen lokalen Rand zu Bulk Prozess**, der einen Ladungsreservoir Rand physisch erzeugt (damit 2+3 wirklich gemeinsam zugelassen sind), danach eine OS positive *nichtzentrale* Zustandszeit ohne freien Gewichtsfaktor auswählt, die anyonische Oberfläche lokal und anomalfrei in 3+1 fermionische Materie überführt und dabei zugleich einen Lorentz Kontinuumsgrenzwert mit echter Gravitationsdynamik und c3=1/(8π) liefert. Ohne diese vier Beweise bleibt es **Stufe 1, KANDIDAT**, nicht F1.

## Literatur

- https://www.damtp.cam.ac.uk/user/tong/qhe/four.pdf
- https://web.math.ucsb.edu/~zhenghwa/data/research/pub/twoeigenvalue-02.pdf
- https://people.tamu.edu/~rowell/DRWv10.pdf
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.085134
- https://msp.org/gt/2018/22-5/gt-v22-n5-p05-s.pdf
