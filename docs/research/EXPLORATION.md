# Explorationsmodus

Entscheidung des Nutzers vom 20. November 2026. Dieser Modus gilt, wenn nach neuen
Ideen, einem Neuansatz, Kreativität oder einer breiten Suche gefragt wird. Er soll
verhindern, dass jede Sitzung wieder am selben Objekt weiterbaut.

## Warum es diesen Modus gibt

Bis zum 20. November las jede Sitzung zuerst rund 40.000 Zeilen Plan, Ergebnisse
und Fortschrittsstand, alle zur selben Aufgabe F1. Danach wurde der „Nächste
Schritt“ des letzten Vertrags fortgesetzt. Auch ausdrückliche Neuansätze landeten
so binnen Stunden wieder bei denselben Objekten. Der Bestand ist wertvoll, aber er
gehört an das Ende einer Idee, nicht an ihren Anfang.

## Ablauf

1. **Vielfalt prüfen.** `ls experiments/exploration/` und die Gebiete der letzten
   fünf Ordner notieren. Mindestens zwei der neuen Ansätze müssen aus Gebieten
   kommen, die dort nicht vorkommen.
2. **Ideenphase ohne Bestand.** Nicht PLAN.md, RESULTS.md, progress.json, Atlas
   oder Verträge lesen. Mindestens fünf Ansätze aus mindestens drei verschiedenen
   Gebieten (Liste unten) in `IDEAS.md` aufschreiben. Je Ansatz: Grundobjekt,
   Regel oder Prinzip, welche Zielgröße herauskommen könnte, freie Parameter.
3. **Sperrvokabular.** In der Ideenphase nicht als Ausgangspunkt verwenden: Naht,
   Seam, Marken, P1, Randkern, τ = i, Viermarkensphäre, D5 + A3, E8 Verklebung,
   Träger 3 + 2, Clocknetz, R113, F1, Universalraum. Wenn ein Ansatz von selbst
   dorthin führt, ist das ein Befund und wird so notiert, nicht der Startpunkt.
4. **Schnelltest.** Für die zwei bis drei stärksten Ansätze die kleinste
   entscheidende Rechnung bauen (Skript im Ordner, exakt wo möglich). Frage:
   Kommt eine Zielgröße heraus, ohne dass sie hineingesteckt wurde?
5. **Zahlenspielkontrolle.** Jeder Treffer braucht eine Basislinie: wie viele
   Ausdrücke gleicher Komplexität aus denselben Bausteinen treffen dieselbe Größe
   mit derselben Genauigkeit? Ein Treffer, den auch zufällige Ausdrücke erreichen,
   zählt nicht.
6. **Abgleich mit dem Bestand.** Erst jetzt Atlas (`search_knowledge`,
   `list_attempts`) und Prüfkatalog `docs/catalog/TFPT_Pruefkatalog_2026-11-20.xlsx`
   befragen: Wurde das schon versucht, unter welchen genauen Voraussetzungen ist
   es ausgeschlossen? Ein Ausschluss trifft nur, wenn seine Voraussetzungen
   wirklich übereinstimmen. Gibt es Anschluss an bekannte TFPT Werte?
7. **Urteil** in `VERDICT.md`: **Kandidat**, **verworfen** oder **offen**, mit
   Begründung, freien Parametern, getroffenen Zielen und Basislinie.

Ein Explorationsergebnis wird nicht automatisch in PLAN, progress oder Atlas
übernommen. Das entscheidet eine Fortsetzungssitzung oder der Nutzer.

## Kandidat

Ein Ansatz ist Kandidat, wenn er mindestens eine Zielgröße der Liste ohne
TFPT Eingaben und ohne angepasste Konstante reproduziert und die Basislinie
deutlich schlägt. Ein Kandidat ist ein positives Ergebnis und wird so berichtet.
Er ist kein Beweis und keine Ledgeränderung.

## Zielgrößen (gemessen)

| Größe | Wert | Quelle |
|---|---|---|
| Feinstrukturkonstante α⁻¹ | 137,035999177(21) | CODATA 2022 |
| Fermionfamilien | 3 | |
| Raumzeit | 3 + 1, Lorentzsignatur | |
| Eichgruppe | SU(3) × SU(2) × U(1), 8 + 3 + 1 Bosonen | |
| Fermionen je Familie | 15 (16 mit rechtshändigem Neutrino) | |
| Hyperladungen | 1/6, 2/3, −1/3, −1/2, −1 (Q = T₃ + Y) | |
| sin²θ_W(M_Z), MS-bar | ≈ 0,2313 | PDG 2024 |
| α_s(M_Z) | 0,1180 | PDG 2024 |
| Koide Q der geladenen Leptonen | 0,666661 (2/3) | PDG Massen |
| m_μ / m_e, m_τ / m_μ | 206,768, 16,817 | PDG 2024 |
| m_p / m_e | 1836,1527 | CODATA 2022 |
| \|V_us\| (Cabibbo) | 0,2243 | PDG 2024 |
| sin²θ₁₂, sin²θ₂₃, sin²θ₁₃ (PMNS) | 0,307, ≈ 0,47–0,57, 0,0220 | NuFIT 6.0 |
| Ω_Λ, Ω_b h² | 0,685, 0,0224 | Planck 2018 |

Werte vor Gebrauch an der Quelle prüfen; die Tabelle ist eine Arbeitsliste.

TFPT Anschlusswerte, **nur für Schritt 6**: c₃ = 1/(8π), φ₀ = 1/(6π) + 48 c₃⁴,
E8 mit 248/240, α Wurzel 137,03599921684.

## Gebiete für die Ideenphase

Zelluläre Automaten und Umschreibsysteme · Kausalmengen · Tensor- und
Spinnetzwerke · Quantenfehlerkorrektur allgemein · Informationsgeometrie und
Fishermetrik · Zufallsmatrizen · Perkolation und kritische Phänomene ·
Renormierung als dynamisches System · Variations- und Extremalprinzipien ·
Spieltheorie und Evolution · Thermodynamik der Berechnung · Graphspektren und
Netzwerke · Matroide und extremale Kombinatorik · Knoten und niedrigdimensionale
Topologie · Kategorien und Operaden · Komplexitätstheorie · Lernende Systeme ·
Nichtkommutative Geometrie · Zahlentheorie (nur mit neuem Objekt, nicht E8).

## Suchprogramm

`research/originlab/genesis_rules.py --foreign` sucht Regeln in Universen ohne
TFPT Vokabular (Graphen, Automaten, Gruppen, Polytope, Codes) und bewertet nach
den Zielzahlen dieser Liste mit Zufallsbasislinie. Erweitern statt ersetzen.
