#!/usr/bin/env python3
"""Z21: Ledgerzensus fuer die TOE-Lueckenkarte (nur lesend).

Liest verification/status_ledger.csv und schreibt nichts ausser auf stdout.
Zwei Klassifikationen:

1. Offizielle Anzeigeklasse, exakt nachgebaut aus
   verification/v649_discipline_audit.py (row_class): AXIOM, X, O, C, E, OTHER.
   Dort zaehlen [N] (numerisch) und [I]/[L]/[F]/[M] als E, [P]/[B]/[R] als C,
   [A] als O, und jedes Statusfeld mit "contract" als O.
2. Feinklasse "schwaechstes Glied" fuer diese Karte: Ein zusammengesetzter
   Status ist so stark wie sein schwaechster Teil. Reihenfolge
   O (offen/[A]) < C (bedingt) < N (numerisch) < E (exakt). Zusaetzlich
   getrennt: AXIOM, REGEL (Firewall/Decision), und ob ein Kill-Test [X]
   angehaengt ist. Diese Feinklasse ersetzt keinen Ledgerstatus; sie
   befoerdert und degradiert nichts.

Aufruf: /Users/stefanhamann/Projekte/tfpt-theoryv4/.venv/bin/python z21_ledger_zensus.py
"""
import csv
import re
import sys
from collections import Counter
from pathlib import Path

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LEDGER = ROOT / "verification" / "status_ledger.csv"

# ---- 1. offizielle Anzeigeklasse (Kopie von v649.row_class, unveraendert) ----
E_TAGS = {"E", "I", "L", "F", "N", "M"}
C_TAGS = {"C", "P", "B", "R"}
O_TAGS = {"O", "A"}
X_TAGS = {"X"}
E_WORDS = ("identity", "formal", "lattice", "numerical", "measured",
           "exact", "red-team")
C_WORDS = ("physical", "branch policy", "geometric premise", "conjecture",
           "reduction")


def tags_of(status):
    tags = set()
    for grp in re.findall(r"\[([A-Za-z/+, ]+)\]", status):
        for t in re.split(r"[/+, ]+", grp):
            t = t.strip().upper()
            if t:
                tags.add(t)
    return tags


def row_class(status):
    low = status.lower()
    tags = tags_of(status)
    if low.startswith("axiom"):
        return "AXIOM"
    if tags & X_TAGS:
        return "X"
    if tags & O_TAGS or low.startswith("open") or "contract" in low:
        return "O"
    if tags & C_TAGS:
        return "C"
    if tags & E_TAGS:
        return "E"
    if any(low.startswith(w) for w in E_WORDS):
        return "E"
    if any(low.startswith(w) for w in C_WORDS):
        return "C"
    return "OTHER"


# ---- 2. Feinklasse: schwaechstes Glied ----
FINE_E = {"E", "I", "L", "F", "M", "E-NUM", "E-FLOAT", "PROOF", "IDENTITY",
          "E*"}
FINE_N = {"N"}
FINE_C = {"C", "P", "B", "R"}
FINE_O = {"O", "A"}


def fine_tags(status):
    tags = set()
    for grp in re.findall(r"\[([^\]]+)\]", status):
        for t in re.split(r"[/+,| ]+", grp):
            t = t.strip().upper()
            if t:
                tags.add(t)
    return tags


def fine_class(status):
    low = status.lower().strip()
    tags = fine_tags(status)
    if low.startswith("axiom"):
        return "AXIOM"
    if low.startswith("firewall") or low.startswith("decision (registered"):
        return "REGEL"
    if tags & FINE_O or low.startswith("open"):
        return "O"
    if tags & FINE_C:
        return "C"
    if tags & FINE_N:
        return "N"
    if tags & FINE_E:
        return "E"
    # Wortfallback in derselben Reihenfolge
    if "contract" in low and "closed" not in low:
        return "O"
    if any(w in low for w in ("physical", "conjecture", "candidate",
                              "branch policy", "geometric premise",
                              "conditional", "near-miss", "input")):
        return "C"
    if low.startswith("numerical") or "numerical" in low[:30]:
        return "N"
    if any(low.startswith(w) for w in ("identity", "formal", "lattice",
                                       "exact", "proven", "closed finite",
                                       "red-team", "measured", "conceptual")):
        return "E"
    return "SONST"


def has_kill(status):
    return "X" in fine_tags(status) or "kill" in status.lower()


# ---- kuratierte tragende Eintraege je Schicht (IDs aus dem Ledger) ----
SCHICHTEN = [
    ("S0 Ursprung der Eingaben (P1, P2, Einheit) -- TOE T1, F1 G0-G2", [
        "AX.P1.01", "AX.P2.01", "P2.PARTITION.01", "TFPT.IRREDUCIBLE.01",
        "QGEO.REALIZE.01", "QGEO.PILLOW.01", "DIMENSION.SELECTOR.4D.01",
        "ANCHOR.VGEO.01", "METROLOGY.CLOSURE.01", "VGEO.SHARPEN.01"]),
    ("S1 Naht = (E8)_1 als physischer Randzustand -- TOE T2, F1 G1/G3", [
        "SEAM.EQUIV.01", "SEAM.EQUIV.MMST.01", "SEAM.EQUIV.TWISTOR.01",
        "SEAM.SIMPLECURRENT.GENERATOR.01", "SEAM.DETLINE.UNIFICATION.01",
        "SEAM.STATE.RPMIXING.01", "SEAM.EQUIV.UNIQ.01",
        "ALPHA.QUILLEN.EXACT.01"]),
    ("S2 2D-Naht -> lokale wechselwirkende 3+1D-Theorie -- TOE T3/T5, F1 G4", [
        "DIMENSION.UPLIFT.FIREWALL.01", "SEAM.BULK4D.RECON.01",
        "BULK.PRIME_COMPLETION.01", "WOIT.OS.TWISTOR.01",
        "QFT4D.OS.RECON.01", "QFT4D.LATTICE.FUNDAMENTAL.01",
        "TFPT4D.LATTICE.ACTION.01", "DYN.UNITARY.DILATION.01",
        "DYN.MARKOV.EMBED.01", "QFT4D.EG.ALLORDER.01",
        "CCC.SEAM.KINEMATICS.01", "CCC.SEAM.CROSSOVER.01"]),
    ("S3 Chirale Eichtheorie und Fermioninhalt -- TOE T4, F1 G5", [
        "CHIRAL4D.NOMIRROR.01", "CHIRAL4D.MIRROR.DFILTER.BLOCK.01",
        "CHIRAL4D.MIRROR.FACTORIZED.BOUND.01", "CAR.SM.ANOM.01",
        "CONTRACT.QFT4D.DIRAC.01", "PS.NCG.FLUCT.01"]),
    ("S4 Kopplungen, Higgs/EW-Brechung, Yukawa, Flavor, starkes CP -- TOE T6", [
        "GAUGE.DETLINE.FIXPOINT.01", "HIGGS.FREESEAM.01", "SM.RG.01",
        "QCD.LAMBDA.01", "FR.MPME.01", "FR.QCD.BUDGET.01", "FLAV.CP.01",
        "FLAV.PMNS.03", "FLAV.PMNS.CLOSE.01", "FLAV.TH13.01",
        "FLAV.TH13.PRESSURE.01", "FLAV.NUSCALE.01", "FLAV.NUSCALE.05",
        "FLAV.NUSCALE.06", "FLAV.NU.TEXTURE.MECHANISM.01",
        "FLAV.NU.PARENT.INTERTWINER.01", "FR.KOIDE.03", "FR.POLE.SOLVE.01",
        "FLAV.STRONGCP.02"]),
    ("S5 Gravitation jenseits Jacobson, Quantengravitation -- TOE T7, F1 G6", [
        "GRAV.NONCIRCULAR.01", "GRAV.NONLINEAR.01",
        "GRAV.ENTROPY.EQUILIBRIUM.01", "GRAV.SPIN2.EMERGENCE.01",
        "GRAV.WEINBERG.WITTEN.01", "GRAV.PARENT.GAP.NO_GO.01",
        "GATE.METRIC.01", "QG.AMB.01", "GRAVITY.COMPLETE.01"]),
    ("S6 Anfangszustand und Kosmologie -- TOE T8", [
        "FTRANSFER.SK.RHO0.01", "FTRANSFER.GENERATING.01",
        "LAMBDA.BRANCH.01", "GRAV.H0.01", "COSMO.NSTAR.01", "GATE.NSTAR.01",
        "AS.RECONCILE.01", "FR.DM.01", "FR.DM.02", "FR.RELIC.SOLVE.01",
        "FR.ETAB.01", "FR.ETAB.04", "FR.BOLTZMANN.SOLVE.01", "COSMO.OMB.01"]),
    ("S7 Datenpruefung und Statistik -- TOE T8/Validated, F1 G7", [
        "PRED.JOINTLIKELIHOOD.01", "OBS.TRANSDUCTION.01", "PRED.LIVE.01",
        "KILL.BOARD.01", "PS.PROTON.02"]),
    ("Dach", ["TFPT.TOE.COMPLETE.01", "RESIDUAL.CERTIFICATION.01"]),
]

# Wo die mechanische Klasse den Statustext falsch liest (von Hand gelesen,
# keine Statusaenderung; der Ledgertext bleibt massgeblich).
LESUNG = {
    "CCC.SEAM.KINEMATICS.01": "Text: Kinematik [E], CCC-Lesart [C]; das [O] "
                              "im Status gehoert zu WOIT.OS.TWISTOR.01",
    "QFT4D.LATTICE.FUNDAMENTAL.01": "Text: typisierte Entscheidung; '[E]' "
                                    "steht dort verneint ('no [E] claim')",
    "GRAV.NONCIRCULAR.01": "Disziplinregel (Firewall), kein Ergebnis",
    "DIMENSION.UPLIFT.FIREWALL.01": "Disziplinregel (Firewall), kein Ergebnis",
}

FOKUS = ["DIMENSION.UPLIFT.FIREWALL.01", "SEAM.BULK4D.RECON.01",
         "WOIT.OS.TWISTOR.01", "QFT4D.OS.RECON.01", "CHIRAL4D.NOMIRROR.01",
         "GRAV.NONCIRCULAR.01", "QG.AMB.01", "GATE.METRIC.01",
         "CCC.SEAM.CROSSOVER.01"]


def short(s, n):
    s = " ".join(s.split())
    return s if len(s) <= n else s[: n - 3] + "..."


def main():
    with open(LEDGER, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    by_id = {r["claim_id"]: r for r in rows}
    active = [r for r in rows
              if (r.get("active") or "").strip().lower() == "true"]

    checks = []

    def check(name, ok):
        checks.append((name, ok))

    print("Z21 Ledgerzensus (nur lesend): verification/status_ledger.csv")
    print("=" * 78)
    print(f"Zeilen gesamt {len(rows)}, aktiv {len(active)}, "
          f"inaktiv {len(rows) - len(active)}")
    check("Ledger gelesen, mindestens 700 Zeilen", len(rows) >= 700)

    # 1. offizielle Anzeigeklasse
    print("\n1. Offizielle Anzeigeklasse (Nachbau v649.row_class)")
    for label, sel in (("alle", rows), ("aktiv", active)):
        c = Counter(row_class(r["status"]) for r in sel)
        print(f"   {label:6s}: " + ", ".join(
            f"{k} {c.get(k, 0)}" for k in
            ("E", "C", "O", "X", "AXIOM", "OTHER")))
    off_act = Counter(row_class(r["status"]) for r in active)

    # 2. Feinklasse
    print("\n2. Feinklasse 'schwaechstes Glied' (nur fuer diese Karte; "
          "keine Statusaenderung)")
    fine = Counter(fine_class(r["status"]) for r in active)
    for k in ("E", "N", "C", "O", "AXIOM", "REGEL", "SONST"):
        print(f"   {k:6s} {fine.get(k, 0)}")
    kill = sum(1 for r in active if has_kill(r["status"]))
    print(f"   davon mit Kill-Test oder [X] im Statusfeld: {kill}")
    check("Feinklassen summieren sich auf die aktiven Zeilen",
          sum(fine.values()) == len(active))

    # 2b. offene Zeilen: rein offen vs gemischt
    o_rows = [r for r in active if row_class(r["status"]) == "O"]
    pure = [r for r in o_rows
            if r["status"].lower().startswith("open")
            or fine_tags(r["status"]) <= FINE_O]
    print(f"\n   Anzeigeklasse O aktiv: {len(o_rows)}; davon Status beginnt "
          f"mit 'Open' oder traegt nur [O]/[A]: {len(pure)}; "
          f"gemischt (exakte Teile + offener Rest): {len(o_rows) - len(pure)}")
    print("   Rein offene aktive Zeilen:")
    for r in pure:
        print(f"     {r['claim_id']:40s} {short(r['status'], 70)}")
    check("offizielle O-Zahl gleich der Zahl offener aktiver Zeilen",
          off_act.get("O", 0) == len(o_rows))

    # 3. TOE-Tore aus TFPT.TOE.COMPLETE.01
    print("\n3. Die acht TOE-Tore (TFPT.TOE.COMPLETE.01) und ihre Kinder")
    toe = by_id["TFPT.TOE.COMPLETE.01"]["claim"]
    head = toe.split("SEPARATION")[0]
    gates = []
    for seg in re.split(r"\s(?=T\d )", head):
        if not re.match(r"T\d ", seg):
            continue
        groups = [m for m in re.finditer(r"\(([^()]*)\)", seg)
                  if re.search(r"[A-Z0-9]+\.[A-Z0-9_.]+\.\d\d", m.group(1))]
        if not groups:
            continue
        last = groups[-1]
        gates.append((seg[:2], seg[3:last.start()].strip(), last.group(1)))
    check("TOE-Vertrag nennt genau acht Tore", len(gates) == 8)
    for g, title, deps in gates:
        print(f"   {g} {short(title, 66)}")
        ids = []
        for d in [d.strip() for d in deps.split(",") if d.strip()]:
            m = re.match(r"^(.*\.)(\d\d)-\.(\d\d)$", d)
            if m:
                ids += [f"{m.group(1)}{k:02d}" for k in
                        range(int(m.group(2)), int(m.group(3)) + 1)]
            else:
                ids.append(d)
        for cid in ids:
            r = by_id.get(cid)
            if r is None:
                print(f"       {cid:38s} FEHLT IM LEDGER")
                continue
            print(f"       {cid:38s} {row_class(r['status']):5s} "
                  f"{fine_class(r['status']):5s} {short(r['status'], 52)}")
    print("   Status des Vertrags:",
          short(by_id["TFPT.TOE.COMPLETE.01"]["status"], 120))

    # 4. Fokusliste aus dem Auftrag
    print("\n4. Fokuszeilen des Auftrags")
    for cid in FOKUS:
        r = by_id.get(cid)
        ok = r is not None and r["active"] == "true"
        check(f"Fokuszeile {cid} vorhanden und aktiv", ok)
        if r:
            print(f"   {cid:32s} {row_class(r['status']):5s} "
                  f"{fine_class(r['status']):5s} {short(r['status'], 80)}")

    # 5. tragende Zeilen je Schicht
    print("\n5. Tragende Zeilen je Schicht (Anzeige | Fein | Kill | Status)")
    missing = []
    for title, ids in SCHICHTEN:
        print(f"\n   {title}")
        cnt = Counter()
        for cid in ids:
            r = by_id.get(cid)
            if r is None:
                missing.append(cid)
                print(f"     {cid:36s} FEHLT")
                continue
            fc = fine_class(r["status"])
            cnt[fc] += 1
            print(f"     {cid:36s} {row_class(r['status']):5s} {fc:5s} "
                  f"{'X' if has_kill(r['status']) else '-'}  "
                  f"{short(r['status'], 64)}")
            if cid in LESUNG:
                print(f"     {'':36s} Lesung: {LESUNG[cid]}")
        print("     Bilanz (fein): " + ", ".join(
            f"{k} {v}" for k, v in sorted(cnt.items())))
    check("alle kuratierten IDs im Ledger", not missing)

    # 6. Atlas-Vergleich (nur Hinweis)
    print("\n6. Hinweis: Das Atlas-Inventar zaehlt mit eigenem Vokabular "
          "(theory_overview: open 26, mixed 546, recorded_established 160, "
          "conditional 34 ...). Es ist Navigation; massgeblich bleibt das "
          "Ledger.")

    print("\nPruefungen")
    for name, ok in checks:
        print(f"   {'bestanden' if ok else 'GESCHEITERT'}: {name}")
    nfail = sum(1 for _, ok in checks if not ok)
    print(f"\n{len(checks) - nfail}/{len(checks)} bestanden. "
          "Nichts geschrieben, nichts befoerdert.")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
