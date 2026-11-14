# F1 Kontrollrechnungen aus den Sitzungen zur Verklebung und Chiralität

| Prüfskript | Ergebnis und Geltungsbereich |
| --- | --- |
| `tfpt_k4_z4_sewing_audit.py` | 41 bestandene Prüfungen. Viermarkengruppe, E8 Verklebung, Z4 Gaußprojektion und Rückkehr im endlichen Modell. |
| `tfpt_f1_sewing_rigidity_probe.py` | 51 bestandene Prüfungen. Keine ausgezeichnete primitive Z4 Charakterwahl unter voller S4 Symmetrie. Zwei konjugierte isotrope Klebungen, freie elliptische Markenkonfiguration und ordnungsabhängige Rückkehr. Gleiche Sektorgeschwindigkeit unter der Voraussetzung einer gemeinsamen chiralen Translation. |
| `tfpt_f1_six_dimensional_index_audit.py` | 7 bestandene Prüfungen. Die explizite Flusswahl auf dem Sechsertorus erreicht Indexbetrag drei und behält elf Spiegelpaare. Die geprüfte Anomaliebedingung und der einfache Orbifoldweg liefern keinen vollständigen Parent. |

Jedes Skript läuft einzeln mit Python, NumPy und SymPy. Zu jedem Skript
liegt die Ausgabe des erneuten Laufs mit gleichem Basisnamen als Text vor.
Die Modelle prüfen mögliche Auswahlregeln und liefern Gegenproben.
Der gemeinsame gewichtete P1 Aufruf mit absoluter Referenz bleibt offen.

Sitzungen vom 13. November 2026. Die ersten beiden Dateien stammen aus
„Dynamische Verklebung TFPT“, die dritte aus der Untersuchung des
sechsdimensionalen chiralen Ursprungs.
