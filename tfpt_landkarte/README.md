# TFPT Landkarte

Eine lokale Lese und Denkhilfe für TFPT: Herleitung, Dynamikkandidaten, Muster,
Ursprung von P1 und P2, Universalraum, Origin, OPH und der vollständige Atlasstand.

## Starten

```sh
python3 -m http.server 8790 --bind 127.0.0.1 --directory tfpt_landkarte
```

Dann http://127.0.0.1:8790 öffnen. Die Seite läuft auch direkt aus `index.html`.

## Daten

- `data/atlas.js` wird erzeugt. Nach Änderungen im Atlas neu bauen:
  `python3 -B tfpt_landkarte/build_data.py`. Das Skript liest nur
  `semantic_explorer/research/program.json`, `semantic_explorer/theory_map.json`
  und `website/public/verification/predictions_frozen.json`.
- `content-kern.js` und `content-weit.js` sind von Hand geschriebene Karten mit
  Status, Zahlen, Stichworten, Verbindungen und Quellenpfaden (Stand 19. November 2026).

## Geltung

Die Karten vereinfachen. Status und die Punkte der Dynamiktabelle sind eine
Einordnung des Korpus und keine Ledgereinträge. Maßgeblich bleiben
`verification/status_ledger.csv`, die Originalquellen und der lebende Atlas.
Notizen im Ideenraum liegen nur im Browser (localStorage).
